use crate::{
    alignment::{Alignment, AlignmentSequence, Strand, TargetGroup},
    alphabet::{GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, UTF8_TO_DIGITAL_NUCLEOTIDE},
    formats::{fasta, AlignmentFormat, FormatCheck, SeekableReader},
    sequence_store::SequenceStore,
    uleb::{self, Cigar, CigarIterator, CigarSegment},
    util::VecMap,
};
use anyhow::{anyhow, Context};
use std::{collections::HashMap, f32, fmt::Display, slice, sync::Arc};
use std::{
    io::{self, SeekFrom},
    string::FromUtf8Error,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BPAFError {
    #[error("io error while reading file.")]
    IOError(#[from] io::Error),
    #[error("file of size {0} is too small to be a valid BPAF")]
    FileToSmall(u64),
    #[error("invalid BPAF header {0:?}")]
    InvalidHeader(BPAFHeader),
    #[error("unsupported BPAF header {0:?}")]
    UnsupportedHeader(BPAFHeader),
    #[error("invalid BPAF footer {0:?}")]
    InvalidFooter(BPAFFooter),
    #[error("unable to decode invalid integer.")]
    InvalidULEB(u64),
    #[error("invalid field, expected size: {0}, actual size: {1}")]
    InvalidField(u64, u64),
    #[error("Reserved bits were set for a record.")]
    InvalidRecordFlags,
    #[error("invalid cigar string, expected length is {0}, but computed length is {1}")]
    InvalidCigar(u64, u64),
    #[error("failed to decode string entry due to: {0}")]
    UTF8Error(#[from] FromUtf8Error),
    #[error("found invalid nucleotide code in sequence: {0}")]
    InvalidNucleotide(char),
}

#[derive(Debug)]
pub struct BPAFHeader {
    pub magic: [u8; 5],
    pub version: u8,
    pub flags: u8,
}

#[derive(Debug)]
pub struct BPAFFooter {
    pub table_offset: u64,
    pub magic: [u8; 5],
    pub footer_offset: u64,
}
#[allow(unused)]
struct ULEBEntry {
    size: u64,
    data: Vec<u8>,
}

struct ULEBEntryIterator<'a, R: SeekableReader> {
    reader: &'a mut R,
    bytes_taken: u64,
    bytes_allowed: Option<u64>,
}

impl<'a, R: SeekableReader> ULEBEntryIterator<'a, R> {
    fn new(reader: &'a mut R, bytes_allowed: Option<u64>) -> Self {
        Self {
            reader,
            bytes_taken: 0,
            bytes_allowed,
        }
    }
}

fn io_decode_next_uleb(
    bytes: &mut impl SeekableReader,
    max_to_take: Option<u64>,
    byte_counter: &mut u64,
) -> Option<Result<u64, BPAFError>> {
    let mut result_int: u64 = 0;
    let mut shift: u8 = 0;
    let mut bytes_read: u64 = 0;

    let max_bytes = 10.min(max_to_take.unwrap_or(10));

    for _ in 0..max_bytes {
        let mut b = 0u8;
        let success = bytes.read_exact(slice::from_mut(&mut b));

        match success {
            Ok(()) => {
                bytes_read += 1;
                *byte_counter += 1;
                let b_checked = b;
                result_int |= ((b_checked & 0b0111_1111) as u64) << shift;
                if (b_checked & 0b1000_0000) == 0 {
                    return Some(Ok(result_int));
                }
                shift += 7;
            }
            Err(err) => {
                if !matches!(err.kind(), io::ErrorKind::UnexpectedEof) {
                    return Some(Err(err.into()));
                }
                break;
            }
        }
    }

    if bytes_read == 0 {
        return None;
    }
    Some(Err(BPAFError::InvalidULEB(result_int)))
}

impl<R: SeekableReader> Iterator for ULEBEntryIterator<'_, R> {
    type Item = Result<ULEBEntry, BPAFError>;

    fn next(&mut self) -> Option<Self::Item> {
        let size = io_decode_next_uleb(
            self.reader,
            self.bytes_allowed
                .map(|v| v.saturating_sub(self.bytes_taken)),
            &mut self.bytes_taken,
        )?;

        match size {
            Ok(size) => {
                //let size = size as usize;
                let mut data = vec![0u8; size as usize];
                let success = self.reader.read_exact(&mut data);

                match success {
                    Ok(_) => {
                        self.bytes_taken += size;
                        let data_len = data.len() as u64;
                        if data_len != size {
                            Some(Err(BPAFError::InvalidField(size, data_len)))
                        } else {
                            Some(Ok(ULEBEntry {
                                size: size,
                                data: data,
                            }))
                        }
                    }
                    Err(err) => Some(Err(err.into())),
                }
            }
            Err(err) => Some(Err(err)),
        }
    }
}

mod bpaf_record_flags {
    pub const ORIENT_C: u8 = 0x01;
    pub const HAS_SCORE: u8 = 0x02;
    pub const HAS_EVAL: u8 = 0x04;
    pub const HAS_DIV: u8 = 0x08;
    pub const HAS_BITSCORE: u8 = 0x10;
    /// Bits 5-7; must be zero in version 0.
    pub const RESERVED: u8 = 0xE0;
}

mod bpaf_feature_flags {
    pub const INCLUDES_SEQUENCES: u8 = 0x01;
    // TODO: Possible proposal later...
    //pub const INCLUDES_MATRICIES: u8 = 0x02;
}

#[allow(unused)]
#[derive(Debug)]
pub struct BPAFRecord {
    pub query_id: u64,
    pub target_id: u64,
    pub query_start: u64,
    pub query_length: u64,
    pub target_start: u64,
    pub target_length: u64,
    pub matrix_id: u64,
    pub strand: Strand,
    pub score: Option<f32>,
    pub e_value: Option<f64>,
    pub divergence: Option<u16>,
    pub bit_score: Option<f32>,
    pub cigar: Cigar,
}

impl Display for BPAFRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            self.query_id,
            self.target_id,
            self.query_start,
            self.query_length,
            self.target_start,
            self.target_length,
            self.matrix_id,
            self.strand.to_string(),
            self.score.map(|v| v.to_string()).unwrap_or("".to_string()),
            self.e_value
                .map(|v| v.to_string())
                .unwrap_or("".to_string()),
            self.divergence
                .map(|v| v.to_string())
                .unwrap_or("".to_string()),
            self.bit_score
                .map(|v| v.to_string())
                .unwrap_or("".to_string()),
        )
    }
}

macro_rules! read_le_from_array {
    ($array:ident, $offset:ident, $int_type:ty) => {
        || -> Result<$int_type, BPAFError> {
            const SIZE: usize = size_of::<$int_type>();
            let end = $offset + SIZE;
            let arr_bytes: [u8; SIZE] = $array[$offset..end]
                .try_into()
                .map_err(|_e| BPAFError::InvalidField($array.len() as u64, end as u64))?;
            $offset += SIZE;
            Ok(<$int_type>::from_le_bytes(arr_bytes))
        }()
    };
}

fn parse_ulebs<const NUM: usize>(
    data: &mut impl Iterator<Item = u8>,
) -> Result<([u64; NUM], usize), BPAFError> {
    let mut first_ulebs = [0u64; NUM];
    let mut byte_offset = 0;

    for i in 0..first_ulebs.len() {
        let (val, bytes_taken) =
            uleb::decode_next_uleb(data).map_err(|e| BPAFError::InvalidULEB(e.0))?;
        first_ulebs[i] = val;
        byte_offset += bytes_taken;
    }

    Ok((first_ulebs, byte_offset))
}

// Note: BPAF as specified in repeat masker identifies target as consensus sequences and queries as whole genomes.
// In aurora we define them the exact opposite. This BPAF reader uses aurora definitions.
fn parse_bpaf_record(entry: ULEBEntry) -> Result<BPAFRecord, BPAFError> {
    let data = entry.data;

    let (
        [target_id, query_id, target_start, target_length, query_start, query_length, matrix_id],
        mut byte_offset,
    ) = parse_ulebs(&mut data.iter().copied())?;

    let flags = read_le_from_array!(data, byte_offset, u8)?;

    let strand = if (flags & bpaf_record_flags::ORIENT_C) == 0 {
        Strand::Forward
    } else {
        Strand::Reverse
    };

    if (flags & bpaf_record_flags::RESERVED) != 0 {
        return Err(BPAFError::InvalidRecordFlags);
    }

    let score = ((flags & bpaf_record_flags::HAS_SCORE) != 0)
        .then(|| read_le_from_array!(data, byte_offset, f32))
        .transpose()?;
    let e_value = ((flags & bpaf_record_flags::HAS_EVAL) != 0)
        .then(|| read_le_from_array!(data, byte_offset, f64))
        .transpose()?;
    let divergence = ((flags & bpaf_record_flags::HAS_DIV) != 0)
        .then(|| read_le_from_array!(data, byte_offset, u16))
        .transpose()?;
    let bit_score = ((flags & bpaf_record_flags::HAS_BITSCORE) != 0)
        .then(|| read_le_from_array!(data, byte_offset, f32))
        .transpose()?;

    let (cigar_pair_count, bytes_taken) =
        uleb::decode_next_uleb(&mut data[byte_offset..].iter().copied())
            .map_err(|e| BPAFError::InvalidULEB(e.0))?;
    byte_offset += bytes_taken;

    let cigar_length = cigar_pair_count * 2 + 1;

    let mut target_count = 0;
    let mut query_count = 0;

    let cigar_iter: Result<Cigar, CigarSegment> =
        CigarIterator::new(data[byte_offset..].iter().copied())
            .take(cigar_length as usize)
            .map(|r| {
                if let Ok(v) = r {
                    match v {
                        CigarSegment::Aligned(c) => {
                            target_count += c;
                            query_count += c;
                        }
                        CigarSegment::QueryGap(c) => {
                            target_count += c;
                        }
                        CigarSegment::TargetGap(c) => {
                            query_count += c;
                        }
                    }
                }
                r
            })
            .collect();
    let cigar = cigar_iter.map_err(|e| BPAFError::InvalidULEB(i64::from(e) as u64))?;

    if target_count != target_length {
        return Err(BPAFError::InvalidCigar(target_length, target_count));
    }
    if query_count != query_length {
        return Err(BPAFError::InvalidCigar(query_length, query_count));
    }

    Ok(BPAFRecord {
        query_id,
        target_id,
        query_start,
        query_length,
        target_start,
        target_length,
        matrix_id,
        strand,
        score,
        e_value,
        divergence,
        bit_score,
        cigar,
    })
}

fn parse_string_table_record(entry: ULEBEntry) -> Result<String, BPAFError> {
    String::from_utf8(entry.data).map_err(|e| e.into())
}

pub struct SequenceEntry {
    pub sequence_id: u64,
    pub start: u64,
    pub remaining: u64,
    pub sequence: Vec<u8>,
}

fn parse_sequence_record(entry: ULEBEntry) -> Result<SequenceEntry, BPAFError> {
    let ([sequence_id, start, remaining], bytes_read) =
        parse_ulebs(&mut entry.data.iter().copied())?;
    let sequence: Result<Vec<u8>, BPAFError> = entry.data[bytes_read..]
        .iter()
        .map(|byte| {
            UTF8_TO_DIGITAL_NUCLEOTIDE
                .get(byte)
                .copied()
                .filter(|&v| !matches!(v, GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL))
                .ok_or_else(|| BPAFError::InvalidNucleotide(*byte as char))
        })
        .collect();

    Ok(SequenceEntry {
        sequence_id,
        start,
        remaining,
        sequence: sequence?,
    })
}

const BPAF_MAGIC: &'static [u8] = b"BPAF\x01";

struct BPAFReader<R: SeekableReader> {
    reader: R,
    header: BPAFHeader,
    footer: BPAFFooter,
    file_offset: u64,
    table_offsets: Vec<u64>,
}

impl<R: SeekableReader> BPAFReader<R> {
    pub fn check_header(reader: &mut R) -> Result<BPAFHeader, BPAFError> {
        let mut header_bytes = [0u8; 7];
        reader.read_exact(&mut header_bytes)?;

        let header = BPAFHeader {
            magic: header_bytes[..5].try_into().unwrap(),
            version: header_bytes[5],
            flags: header_bytes[6],
        };

        if header.magic != BPAF_MAGIC {
            Err(BPAFError::InvalidHeader(header))
        } else if header.version != 0 {
            Err(BPAFError::UnsupportedHeader(header))
        } else {
            Ok(header)
        }
    }

    fn check_footer(reader: &mut R) -> Result<BPAFFooter, BPAFError> {
        let orig_offset = reader.stream_position()?;
        let mut footer_bytes = [0u8; 13];
        let into_file = reader.seek(io::SeekFrom::End(-13))?;
        if into_file < 10 {
            return Err(BPAFError::FileToSmall(into_file + 13));
        }
        reader.read_exact(&mut footer_bytes)?;
        reader.seek(io::SeekFrom::Start(orig_offset))?;

        let footer = BPAFFooter {
            table_offset: u64::from_le_bytes(footer_bytes[..8].try_into().unwrap()),
            magic: footer_bytes[8..].try_into().unwrap(),
            footer_offset: into_file,
        };

        if footer.magic != BPAF_MAGIC {
            Err(BPAFError::InvalidFooter(footer))
        } else {
            Ok(footer)
        }
    }

    pub fn new(mut reader: R) -> Result<Self, BPAFError> {
        let file_offset = reader.stream_position()?;
        let header = Self::check_header(&mut reader)?;
        let footer = Self::check_footer(&mut reader)?;

        let first_table_offset = footer.table_offset;

        Ok(Self {
            reader: reader,
            header,
            footer,
            file_offset,
            table_offsets: vec![first_table_offset],
        })
    }

    fn skip_table(&mut self, offset: &mut u64) -> Result<(), BPAFError> {
        let entries = io_decode_next_uleb(
            &mut self.reader,
            Some(self.footer.footer_offset - *offset),
            offset,
        )
        .unwrap_or(Err(BPAFError::InvalidField(1, 0)))?;

        for _ in 0..entries {
            let entry_size = io_decode_next_uleb(
                &mut self.reader,
                Some(self.footer.footer_offset - *offset),
                offset,
            )
            .unwrap_or(Err(BPAFError::InvalidULEB(0)))?;

            if *offset + entry_size > self.footer.footer_offset {
                return Err(BPAFError::InvalidField(
                    entry_size,
                    self.footer.footer_offset - *offset,
                ));
            }
            self.reader.seek(SeekFrom::Current(entry_size as i64))?;
            *offset += entry_size;
        }

        Ok(())
    }

    fn load_table_offsets(&mut self, up_to: usize) -> Result<u64, BPAFError> {
        let mut offset = *self
            .table_offsets
            .last()
            .unwrap_or(&self.footer.table_offset);
        let start_idx = self.table_offsets.len();

        self.reader.seek(io::SeekFrom::Start(offset))?;

        for _ in start_idx..=up_to {
            self.skip_table(&mut offset)?;
            self.table_offsets.push(offset);
        }

        Ok(offset)
    }

    fn seek_to_table(&mut self, table_idx: usize) -> Result<(u64, u64), BPAFError> {
        let mut offset = if let Some(&val) = self.table_offsets.get(table_idx) {
            self.reader.seek(SeekFrom::Start(val))?;
            val
        } else {
            self.load_table_offsets(table_idx)?
        };

        let entries = io_decode_next_uleb(
            &mut self.reader,
            Some(self.footer.footer_offset - offset),
            &mut offset,
        )
        .unwrap_or(Err(BPAFError::InvalidField(1, 0)))?;

        Ok((entries, offset))
    }

    fn read_table(
        &mut self,
        offset: usize,
    ) -> impl Iterator<Item = Result<ULEBEntry, BPAFError>> + '_ {
        let res = self.seek_to_table(offset);
        let (entries, bytes_allowed) = res
            .as_ref()
            .map(|v| (v.0, self.footer.footer_offset - v.1))
            .unwrap_or((0, 0));

        // TODO: Add way of recording next table offset if iterator is fully consumed...
        res.err()
            .map(|v| Err(v.into()))
            .into_iter()
            .chain(
                ULEBEntryIterator::new(&mut self.reader, Some(bytes_allowed))
                    .take(entries as usize),
            )
            .take_while(|v| v.is_ok())
    }

    pub fn read_target_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(0)
            .map(|r| r.map(parse_string_table_record).flatten())
    }

    pub fn read_query_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(1)
            .map(|r| r.map(parse_string_table_record).flatten())
    }

    pub fn read_matrix_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(2)
            .map(|r| r.map(parse_string_table_record).flatten())
    }

    pub fn read_target_sequences(
        &mut self,
    ) -> Option<impl Iterator<Item = Result<SequenceEntry, BPAFError>> + '_> {
        (self.header.flags & bpaf_feature_flags::INCLUDES_SEQUENCES != 0).then(|| {
            self.read_table(3)
                .map(|r| r.map(parse_sequence_record).flatten())
        })
    }

    pub fn read_query_sequences(
        &mut self,
    ) -> Option<impl Iterator<Item = Result<SequenceEntry, BPAFError>> + '_> {
        (self.header.flags & bpaf_feature_flags::INCLUDES_SEQUENCES != 0).then(|| {
            self.read_table(4)
                .map(|r| r.map(parse_sequence_record).flatten())
        })
    }

    pub fn read_records(&mut self) -> impl Iterator<Item = Result<BPAFRecord, BPAFError>> + '_ {
        let seek_err: Option<Result<BPAFRecord, BPAFError>> = self
            .reader
            .seek(io::SeekFrom::Start(self.file_offset + 7))
            .err()
            .map(|v| Err(v.into()));
        let bytes_taken = self.footer.table_offset - 7;

        seek_err
            .into_iter()
            .chain(
                ULEBEntryIterator::new(&mut self.reader, Some(bytes_taken))
                    .map(|r| r.map(parse_bpaf_record).flatten()),
            )
            .take_while(|r| r.is_ok())
    }
}

impl FromIterator<SequenceEntry> for SequenceStore {
    fn from_iter<T: IntoIterator<Item = SequenceEntry>>(iter: T) -> Self {
        let mut seq_store = SequenceStore::new();

        for entry in iter {
            seq_store.add_sequence(
                entry.sequence_id as usize,
                entry.start as usize,
                &entry.sequence,
            );
        }

        seq_store
    }
}

pub struct BPAFFormat {}

fn map_entry_id(
    new_map: &VecMap<String>,
    old_map: &Vec<String>,
    value: Result<SequenceEntry, BPAFError>,
) -> Result<SequenceEntry, anyhow::Error> {
    let unwrapped_value = value?;
    let seq_name = old_map
        .get(unwrapped_value.sequence_id as usize)
        .with_context(|| {
            format!(
                "Sequence id: {} is out of bounds!",
                unwrapped_value.sequence_id
            )
        })?;
    let new_seq_id = new_map
        .key(seq_name)
        .with_context(|| format!("Can't find sequence '{}' in provided sequences!", seq_name))?;

    Ok(SequenceEntry {
        sequence_id: new_seq_id as u64,
        start: unwrapped_value.start,
        remaining: unwrapped_value.remaining,
        sequence: unwrapped_value.sequence,
    })
}

impl AlignmentFormat for BPAFFormat {
    fn format_check(
        primary_reader: &mut impl SeekableReader,
        _secondary_reader: Option<&mut impl SeekableReader>,
    ) -> anyhow::Result<FormatCheck> {
        let header = BPAFReader::check_header(primary_reader);
        match header {
            Ok(_) => Ok(FormatCheck::Valid),
            Err(e @ BPAFError::InvalidHeader(_)) => Ok(FormatCheck::Invalid(e.to_string())),
            Err(e @ BPAFError::UnsupportedHeader(_)) => Ok(FormatCheck::Invalid(e.to_string())),
            Err(v) => Err(v.into()),
        }
    }

    fn read(
        primary_reader: &mut impl SeekableReader,
        secondary_reader: Option<&mut impl SeekableReader>,
        substitution_matrices: VecMap<crate::substitution_matrix::SubstitutionMatrix>,
    ) -> anyhow::Result<crate::alignment::AlignmentData> {
        let mut reader = BPAFReader::new(primary_reader)?;

        let queries_old = Result::<Vec<String>, BPAFError>::from_iter(reader.read_query_names())?;
        let targets_old = Result::<Vec<String>, BPAFError>::from_iter(reader.read_target_names())?;
        let sub_matrix_names_old: Vec<String> =
            Result::<Vec<String>, BPAFError>::from_iter(reader.read_matrix_names())?;

        let mut queries_new = VecMap::new();
        let mut targets_new = VecMap::new();
        let mut target_groups = Vec::new();

        // Iterate records to extract names actually used from the BPAF, we only keep those...
        for (idx, record_or_error) in reader.read_records().enumerate() {
            let error = |e| anyhow!("Record {}: {}", idx + 1, e);
            let record = record_or_error?;

            if record.target_id as usize >= targets_old.len() {
                return Err(error("Target id out of bounds."));
            }
            if record.query_id as usize >= queries_old.len() {
                return Err(error("Query id out of bounds."));
            }
            if record.matrix_id as usize >= sub_matrix_names_old.len() {
                return Err(error("Matrix id out of bounds."));
            }

            let target_id = targets_new.insert(targets_old[record.target_id as usize].clone());
            let query_id = queries_new.insert(queries_old[record.query_id as usize].clone());
            let matrix_id = substitution_matrices
                .key_by_name(&sub_matrix_names_old[record.matrix_id as usize])
                .with_context(|| {
                    format!(
                        "No matrix with name: {}",
                        sub_matrix_names_old[record.matrix_id as usize]
                    )
                })?;

            let target_start = record.target_start as usize;
            let target_end = (record.target_start + record.target_length - 1) as usize;

            let target_group = match target_groups.get_mut(target_id) {
                Some(group) => group,
                None => {
                    target_groups.push(TargetGroup {
                        target_id,
                        target_start,
                        target_end,
                        alignments: vec![],
                        tandem_repeats: vec![],
                    });
                    target_groups.last_mut().unwrap()
                }
            };

            let (query_start, query_end) = match record.strand {
                Strand::Forward => (
                    record.query_start as usize,
                    (record.query_start + record.query_length - 1) as usize,
                ),
                Strand::Reverse => (
                    (record.query_start + record.query_length - 1) as usize,
                    record.query_start as usize,
                ),
                _ => return Err(error("Invalid strand value!")),
            };

            let sequence = AlignmentSequence {
                query_seq: Arc::default(),
                target_seq: Arc::default(),
                cigar: record.cigar,
            };

            target_group.alignments.push(Alignment {
                sequence,
                query_id,
                target_start,
                target_end,
                query_start,
                query_end,
                strand: record.strand,
                // Set later...
                id: 0,
                substitution_matrix_id: matrix_id,
            });
            target_group.target_start = target_group.target_start.min(target_start);
            target_group.target_end = target_group.target_end.max(target_end);
        }

        let mut query_lengths: HashMap<usize, usize> = HashMap::new();

        let mut query_store = match reader.read_query_sequences() {
            Some(iter) => Result::<SequenceStore, anyhow::Error>::from_iter(iter.map(|r| {
                let r_new = map_entry_id(&queries_new, &queries_old, r);
                if let Ok(v) = r_new.as_ref() {
                    query_lengths.insert(
                        v.sequence_id as usize,
                        (v.start.saturating_sub(1)) as usize
                            + v.sequence.len()
                            + v.remaining as usize,
                    );
                }
                r_new
            }))?,
            None => SequenceStore::new(),
        };

        let mut target_store = match reader.read_target_sequences() {
            Some(iter) => Result::<SequenceStore, anyhow::Error>::from_iter(
                iter.map(|v| map_entry_id(&targets_new, &targets_old, v)),
            )?,
            None => SequenceStore::new(),
        };

        match secondary_reader {
            Some(second_reader) => {
                if target_store.sequence_count() > 0 || query_store.sequence_count() > 0 {
                    eprintln!("BPAF already contains sequence data, extending with FASTA data...");
                }

                fasta::parse_fasta_file(
                    second_reader,
                    &targets_new,
                    &queries_new,
                    &mut target_store,
                    &mut query_store,
                    Some(&mut query_lengths),
                )?;
            }
            None => {
                if target_store.sequence_count() == 0 || query_store.sequence_count() == 0 {
                    return Err(anyhow!(
                        "BPAF contains no sequences, please provide a fasta file!"
                    ));
                }
            }
        }

        let query_sequences = query_store.into_index();
        let target_sequences = target_store.into_index();

        // Get alignment references for every sequence...
        for group in target_groups.iter_mut() {
            let target_id = group.target_id;
            for al in group.alignments.iter_mut() {
                let (q_start, q_end) = al.ordered_query_range();
                let t_start = al.target_start;
                let t_end = al.target_end;

                al.sequence.target_seq = target_sequences
                    .find(target_id, t_start, t_end)
                    .ok_or(anyhow!("Unable to find needed target sequence!"))?;
                al.sequence.query_seq = query_sequences
                    .find(al.query_id, q_start, q_end)
                    .ok_or(anyhow!("Unable to find needed query sequence!"))?;
            }
        }

        Ok(crate::alignment::AlignmentData {
            target_groups,
            target_name_map: targets_new,
            query_name_map: queries_new,
            substitution_matrices,
            target_sequences,
            query_sequences,
            query_lengths,
        })
    }

    fn name() -> &'static str {
        "BPAF"
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;
    use std::io::Cursor;

    const BPAF_FILE: &'static [u8] = include_bytes!("../../fixtures/test/human-1mb.bpaf");

    #[test]
    fn minimal_empty_file_bytes() -> anyhow::Result<()> {
        #[rustfmt::skip]
        let expected: Vec<u8> = vec![
            b'B', b'P', b'A', b'F', 0x01,   // header MAGIC
            0x00,                           // version
            0x00,                           // extensions
            0x00,                           // qids count = 0
            0x00,                           // tids count = 0
            0x00,                           // mats count = 0
            0x07, 0, 0, 0, 0, 0, 0, 0,      // tables_start = 7
            b'B', b'P', b'A', b'F', 0x01,   // footer MAGIC
        ];

        let mut r = BPAFReader::new(Cursor::new(expected))?;
        assert_eq!(r.read_records().collect_vec().len(), 0);
        assert_eq!(r.read_query_names().collect_vec().len(), 0);
        assert_eq!(r.read_target_names().collect_vec().len(), 0);
        assert_eq!(r.read_matrix_names().collect_vec().len(), 0);
        assert!(r.read_query_sequences().is_none());
        assert!(r.read_target_sequences().is_none());

        Ok(())
    }

    #[test]
    fn test_table_reading() -> anyhow::Result<()> {
        let mut reader = BPAFReader::new(Cursor::new(BPAF_FILE))?;

        let query_names =
            Result::<Vec<_>, BPAFError>::from_iter(reader.read_query_names().take(10))?;
        assert_eq!(
            query_names,
            [
                "DF0000001",
                "DF0000002",
                "DF0000003",
                "DF0000004",
                "DF0000005",
                "DF0000007",
                "DF0000016",
                "DF0000017",
                "DF0000023",
                "DF0000024"
            ]
        );

        let target_names = Result::<Vec<_>, BPAFError>::from_iter(reader.read_target_names())?;
        assert_eq!(target_names, ["Human"]);

        let matrix_names = Result::<Vec<_>, BPAFError>::from_iter(reader.read_matrix_names())?;
        assert_eq!(
            matrix_names,
            [
                "20p41g.matrix",
                "20p49g.matrix",
                "20p51g.matrix",
                "20p53g.matrix",
                "20p45g.matrix",
                "20p43g.matrix",
                "20p39g.matrix"
            ]
        );

        assert_eq!(
            Result::<Vec<_>, BPAFError>::from_iter(reader.read_records().take(10))?
                .iter()
                .map(|v| v.to_string())
                .collect_vec(),
            [
                "0,0,14,135,1802,130,0,+,493,0.00000000000000037196267569848027,,80.22068",
                "0,0,22,166,2496,160,0,-,420,0.0000000000010478045375488144,,68.76076",
                "0,0,183,74,3970,74,0,-,283,0.0000031226739722775843,,47.253777",
                "0,0,18,103,5332,103,0,-,338,0.000000007859169247470169,,55.887966",
                "0,0,108,152,17834,146,0,-,588,0.000000000000000000012051981341223158,,95.134285",
                "0,0,84,61,30825,59,0,+,239,0.00037484061777854835,,40.346428",
                "0,0,65,122,33206,106,0,+,282,0.0000034816404241435072,,47.096794",
                "0,0,29,193,62693,182,1,+,441,0.000000000000007411720003251736,,75.904106",
                "0,0,113,149,82975,134,1,-,460,0.000000000000000833618593949468,,79.05645",
                "0,0,4,222,90191,211,1,+,584,0.0000000000000000000005343299388178203,,99.62968"
            ]
        );

        Ok(())
    }
}

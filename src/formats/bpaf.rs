use crate::{
    alignment::Strand,
    formats::SeekableReader,
    uleb::{self, Cigar, CigarIterator, CigarSegment},
};
use std::{f32, slice};
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
    #[error("invalid cigar string, expected length is {0}, but computed length is {1}")]
    InvalidCigar(u64, u64),
    #[error("failed to decode string entry due to: {0}")]
    UTF8Error(#[from] FromUtf8Error),
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

pub struct BPAFRecord {
    query_id: u64,
    target_id: u64,
    query_start: u64,
    query_length: u64,
    target_start: u64,
    target_length: u64,
    matrix_id: u64,
    strand: Strand,
    score: Option<f32>,
    e_value: Option<f64>,
    divergence: Option<u16>,
    bit_score: Option<f32>,
    cigar: Cigar,
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

fn parse_bpaf_record(entry: ULEBEntry) -> Result<BPAFRecord, BPAFError> {
    let mut byte_offset = 0;
    let mut first_ulebs = [0u64, 7];
    let data = entry.data;

    for i in 0..first_ulebs.len() {
        let (val, bytes_taken) = uleb::decode_next_uleb(&mut data[byte_offset..].iter().copied())
            .map_err(|e| BPAFError::InvalidULEB(e.0))?;
        first_ulebs[i] = val;
        byte_offset += bytes_taken;
    }

    let flags = read_le_from_array!(data, byte_offset, u8)?;

    let strand = if (flags & bpaf_record_flags::ORIENT_C) == 0 {
        Strand::Forward
    } else {
        Strand::Reverse
    };

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

    if target_count != first_ulebs[5] {
        return Err(BPAFError::InvalidCigar(first_ulebs[5], target_count));
    }
    if query_count != first_ulebs[3] {
        return Err(BPAFError::InvalidCigar(first_ulebs[3], query_count));
    }

    Ok(BPAFRecord {
        query_id: first_ulebs[0],
        target_id: first_ulebs[1],
        query_start: first_ulebs[2],
        query_length: first_ulebs[3],
        target_start: first_ulebs[4],
        target_length: first_ulebs[5],
        matrix_id: first_ulebs[6],
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

struct BPAFReader<R: SeekableReader> {
    reader: R,
    header: BPAFHeader,
    footer: BPAFFooter,
    file_offset: u64,
    table_offsets: Vec<u64>,
}

impl<R: SeekableReader> BPAFReader<R> {
    const BPAF_MAGIC: &'static [u8] = b"BPAF\x01";

    pub fn check_header(reader: &mut R) -> Result<BPAFHeader, BPAFError> {
        let mut header_bytes = [0u8; 7];
        reader.read_exact(&mut header_bytes)?;

        let header = BPAFHeader {
            magic: header_bytes[..5].try_into().unwrap(),
            version: header_bytes[5],
            flags: header_bytes[6],
        };

        if header.magic != Self::BPAF_MAGIC {
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

        if footer.magic != Self::BPAF_MAGIC {
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
            self.skip_table(&mut offset);
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

    pub fn read_query_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(0)
            .map(|r| r.map(parse_string_table_record).flatten())
    }

    pub fn read_target_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(1)
            .map(|r| r.map(parse_string_table_record).flatten())
    }

    pub fn read_matrix_names(&mut self) -> impl Iterator<Item = Result<String, BPAFError>> + '_ {
        self.read_table(2)
            .map(|r| r.map(parse_string_table_record).flatten())
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

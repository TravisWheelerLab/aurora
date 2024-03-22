use anyhow::Result;
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::{fmt, hash};

use serde::{ser::SerializeStruct, Serialize, Serializer};

use crate::alphabet::{
    ALIGNMENT_ALPHABET_STR, DASH_UTF8, FORWARD_SLASH_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
    NUCLEOTIDE_ALPHABET_UTF8, PLUS_UTF8, UTF8_TO_DIGITAL_NUCLEOTIDE,
};
use crate::pipeline::StrSliceExt;
use crate::substitution_matrix::SubstitutionMatrix;

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    #[default]
    Unset,
}

impl Strand {
    pub fn from_str(str: &str) -> Self {
        match str {
            "+" => Self::Forward,
            "-" => Self::Reverse,
            str => panic!("unknown strand str: {}", str),
        }
    }
}

impl Serialize for Strand {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unset => write!(f, "?"),
        }
    }
}

#[derive(Default, Debug, Eq)]
pub struct Alignment {
    pub target_seq: Vec<u8>,
    pub query_seq: Vec<u8>,
    pub target_start: usize,
    pub target_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
    pub id: usize,
    pub query_id: usize,
    pub substitution_matrix_id: usize,
}

impl Alignment {
    #[allow(dead_code)]
    pub fn from_str(str: &str) -> Self {
        let tokens: Vec<&str> = str.split('\n').collect();

        let target = tokens[0];
        let query = tokens[1];
        assert_eq!(target.len(), query.len());

        let target_seq = target.to_digital_nucleotides();
        let query_seq = query.to_digital_nucleotides();

        let target_len = target_seq
            .iter()
            .filter(|&&b| b != GAP_OPEN_DIGITAL && b != GAP_EXTEND_DIGITAL)
            .count();

        let query_len = query_seq
            .iter()
            .filter(|&&b| b != GAP_OPEN_DIGITAL && b != GAP_EXTEND_DIGITAL)
            .count();

        Self {
            target_seq,
            query_seq,
            target_start: 1,
            target_end: target_len,
            query_start: 1,
            query_end: query_len,
            strand: Strand::Forward,
            id: 0,
            query_id: 0,
            substitution_matrix_id: 0,
        }
    }
}

impl PartialEq for Alignment {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl hash::Hash for Alignment {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: T|{}-{} Q|{}-{}",
            self.query_id, self.target_start, self.target_end, self.query_start, self.query_end
        )
    }
}

impl Serialize for Alignment {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Alignment", 5)?;
        state.serialize_field("query", &self.query_id)?;
        state.serialize_field("queryStart", &self.query_start)?;
        state.serialize_field("queryEnd", &self.query_end)?;
        state.serialize_field("targetStart", &self.target_start)?;
        state.serialize_field("targetEnd", &self.target_end)?;
        state.serialize_field("row", &self.query_id)?;
        state.serialize_field("strand", &self.strand.to_string())?;
        state.end()
    }
}

pub fn caf_str_to_digital_nucleotides(caf_str: &str) -> (Vec<u8>, Vec<u8>) {
    //  Robert's notes on the CAF format:
    //
    //      *** THE GENOME IS THE QUERY HERE ***
    //      *** THE SUBJECT IS THE MODEL/TE  ***
    //      Yet another Compressed Alignment Format (yaCAF or just CAF).
    //      This format was developed for the use case where sequence databases
    //      may not be available for either the query or the subject of an
    //      alignment and where it's still desirable to communicate the alignment
    //      in a semi-succinct fashion.
    //
    //      Three basic inline string operators are provided: "/" for substitutions,
    //      "+" for insertions (relative to the query) and "-" for deletions.
    //
    //      For example the exact alignment:
    //
    //        Query: AATTGG
    //        Subj : AATTGG
    //
    //      would not need any of these operators and would be encoded using
    //      the single string "AATTGG".
    //
    //      Substitutions are encocoded as query_base/subj_base.  For example:
    //        Query: AAGAA
    //                 |
    //        Subj : AACAA
    //
    //      would be encoded as: "AAG/CAA"
    //
    //      Finally gaps are encoded by surrounding the deleted sequence or the
    //      inserted sequence (relative to the query) by either "+" or "-".  For
    //      instance the following alignment:
    //
    //        Query: AAGCTA--A
    //        Subj : AA--TAGGA
    //
    //      would be encoded as: "AA-GC-TA+GG+A"
    #[derive(Copy, Clone, Debug)]
    pub enum CafState {
        Match,
        TargetGap,
        QueryGap,
        Mutation,
    }

    let mut prev_state = CafState::Match;
    let caf_str_bytes = caf_str.as_bytes();

    let mut target_bytes_digital: Vec<u8> = vec![];
    let mut query_bytes_digital: Vec<u8> = vec![];

    let mut ali_idx = 0usize;
    for &utf8_byte in caf_str_bytes {
        let new_state = match utf8_byte {
            b if NUCLEOTIDE_ALPHABET_UTF8.contains(&b) => match prev_state {
                CafState::Mutation => CafState::Match,
                other => other,
            },
            DASH_UTF8 => match prev_state {
                CafState::TargetGap => CafState::Match,
                _ => CafState::TargetGap,
            },
            PLUS_UTF8 => match prev_state {
                CafState::QueryGap => CafState::Match,
                _ => CafState::QueryGap,
            },
            FORWARD_SLASH_UTF8 => CafState::Mutation,
            unknown => panic!(
                "unknown byte: {}",
                String::from_utf8(vec![unknown]).unwrap()
            ),
        };

        let digital_byte = match UTF8_TO_DIGITAL_NUCLEOTIDE.get(&utf8_byte) {
            Some(byte) => *byte,
            None => 255,
        };

        match (prev_state, new_state) {
            (CafState::Match, CafState::Match) => {
                // AA-GC-TA
                // ^^    ^^
                // this will treat mutation starts as a match
                // position, but we will retroactively fix it
                // when we pop in out of the mutation state later
                target_bytes_digital.push(digital_byte);
                query_bytes_digital.push(digital_byte);
                ali_idx += 1;
            }
            (CafState::Mutation, CafState::Match) => {
                // A A G / C A A
                //         ^
                // we need to back up and fix the target sequence
                // because it will have been erroneously set as if
                // it were a match position a couple states back
                query_bytes_digital[ali_idx - 1] = digital_byte;
            }
            (CafState::TargetGap, CafState::TargetGap) => {
                // AA-GC-TA
                //    ^^
                target_bytes_digital.push(digital_byte);
                match query_bytes_digital[ali_idx - 1] {
                    GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                        query_bytes_digital.push(GAP_EXTEND_DIGITAL)
                    }
                    _ => query_bytes_digital.push(GAP_OPEN_DIGITAL),
                }
                ali_idx += 1;
            }
            (CafState::QueryGap, CafState::QueryGap) => {
                // AA+GC+TA
                //    ^^
                match target_bytes_digital[ali_idx - 1] {
                    GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                        target_bytes_digital.push(GAP_EXTEND_DIGITAL)
                    }
                    _ => target_bytes_digital.push(GAP_OPEN_DIGITAL),
                }
                query_bytes_digital.push(digital_byte);
                ali_idx += 1;
            }
            // ----
            // AA+GC+TA
            //      ^
            (CafState::QueryGap, CafState::Match) |
            // AA-GC-TA
            //      ^
            (CafState::TargetGap, CafState::Match) |
            // AA-GC-TA
            //   ^
            (CafState::Match, CafState::TargetGap) |
            // AA+GC+TA
            //   ^
            (CafState::Match, CafState::QueryGap) |
            // AAG/CAA
            //    ^
            (CafState::Match, CafState::Mutation) => {
                // valid transitions that have no effect
            }
            // ----
            (prev, new) => panic!("invalid CAF state transition: {:?} -> {:?}", prev, new),
        }

        prev_state = new_state;
    }
    (target_bytes_digital, query_bytes_digital)
}

/// A simple Vec-based map that facilitates mapping
/// between usize keys and type <T> values
pub struct VecMap<T: std::cmp::PartialEq> {
    values: Vec<T>,
}

impl<T: std::cmp::PartialEq> VecMap<T> {
    pub fn new() -> Self {
        Self { values: vec![] }
    }

    pub fn from(values: Vec<T>) -> Self {
        Self { values }
    }

    /// Inserts the value and returns the key. If the
    /// value was already in the VecMap, return the key.
    pub fn insert(&mut self, value: T) -> usize {
        match self.contains(&value) {
            true => self.key(&value),
            false => {
                self.values.push(value);
                self.values.len() - 1
            }
        }
    }

    pub fn get(&self, key: usize) -> &T {
        debug_assert!(key < self.values.len(), "invalid key: {key}");
        &self.values[key]
    }

    pub fn contains(&self, value: &T) -> bool {
        self.values.contains(value)
    }

    pub fn size(&self) -> usize {
        self.values.len()
    }

    /// Get the key associated with the value.
    /// This panics if the value is not in the VecMap.
    pub fn key(&self, value: &T) -> usize {
        self.values
            .iter()
            .enumerate()
            .find(|(_, n)| *n == value)
            .expect("key not found")
            .0
    }
}

impl<T: std::cmp::PartialEq + fmt::Debug> fmt::Debug for VecMap<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (0..self.size()).for_each(|key| {
            writeln!(f, "{key}: {:?}", self.values[key]).unwrap();
        });
        Ok(())
    }
}

#[derive(Debug, Default)]
pub struct TandemRepeat {
    pub id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub consensus_pattern: String,
    pub period: usize,
    pub scores: Vec<f64>,
}

#[derive(Deserialize)]
#[serde(rename_all = "PascalCase")]
struct UltraJson {
    pub repeats: Vec<UltraRecord>,
}

#[derive(Deserialize)]
#[serde(rename_all = "PascalCase")]
struct UltraRecord {
    pub sequence_name: String,
    pub start: usize,
    pub length: usize,
    pub consensus: String,
    pub period: usize,
    pub position_score_deltas: Vec<f64>,
}

/// A group of alignments that share
/// the same target sequence
pub struct TargetGroup {
    pub target_id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub alignments: Vec<Alignment>,
    pub tandem_repeats: Vec<TandemRepeat>,
}

/// This holds all input alignments along
/// with the maps for target names, query
/// names, and substitution matrices.
pub struct AlignmentData {
    pub target_groups: Vec<TargetGroup>,
    pub target_name_map: VecMap<String>,
    pub query_name_map: VecMap<String>,
    pub query_lengths: HashMap<usize, usize>,
    pub substitution_matrices: VecMap<SubstitutionMatrix>,
}

impl AlignmentData {
    pub fn from_caf_and_ultra_and_matrices<C: Read, U: Read, M: Read>(
        caf: C,
        ultra: Option<U>,
        matrices: M,
    ) -> Result<Self> {
        // Robert's notes on CAF:
        //   0: score - bit, raw, complexity adjusted or evalue.
        //   1: Percent Substitution - Percent of mismatched non-gap characters in the alignment
        //   2: Percent Deletion - Percent of deletion characters in alignment
        //   3: Percent Insertion - Percent of insertion characters in alignment
        //      *** THE GENOME IS THE QUERY HERE
        //   4: Query Sequence ID
        //   5: Query Start - 1-based, fully closed
        //   6: Query End - 1-based, fully closed
        //   7: Query Remaining - Remaining length of query sequence
        //      *** THE SUBJECT IS THE MODEL/TE
        //   8: Subject Sequence ID - Subject sequence is generally the TE family model for our use cases
        //   9: Subject Classification - [optional] The Dfam/RepeatMasker classification for the TE family
        //  10: Subject Start - 1 based, fully closed
        //  11: Subject End - 1 based, fully closed
        //  12: Subject Remaining - Remaining length of subject sequence
        //  13: Orientation - 0=plus_strand, 1=negative_strand
        //  14: Overlap - [optional] Overlapping annotations from RepeatMasker are flagged using this field
        //  15: Linkage_ID - [optional] RepeatMasker linkage id
        //  16: CAF encoded alignment string
        //  17: Matrix - [optional] The matrix used in scoring the alignment encoded as ##p##g.matrix or simply ##p##g
        //
        //  example record:
        //    0: 199
        //    1: 12.12
        //    2: 0.00
        //    3: 0.00
        //    4: 6
        //    5: 18172245
        //    6: 18172277
        //    7: 58603
        //    8: DF0000023
        //    9: <empty>
        //   10: 1
        //   11: 33
        //   12: 2673
        //   13: 1
        //   14: <empty>
        //   15: <empty>
        //   16: ACCT/GGA/TCT/CGTGGCCT/CGGGGGTTGGGGACCCCTG
        //   17: 14p41g.matrix

        let substitution_matrices = VecMap::from(SubstitutionMatrix::parse(matrices));

        let mut target_groups: Vec<TargetGroup> = vec![];
        let mut target_name_map: VecMap<String> = VecMap::new();
        let mut query_name_map: VecMap<String> = VecMap::from(vec!["skip".into()]);
        let mut query_lengths: HashMap<usize, usize> = HashMap::new();
        query_lengths.insert(0, 0);

        let caf_lines = BufReader::new(caf).lines();

        caf_lines
            .map(|l| l.expect("failed to read line"))
            .filter(|l| !l.is_empty())
            .enumerate()
            .for_each(|(line_num, line)| {
                let tokens: Vec<&str> = line.split(',').collect();

                let target_name = tokens[4].to_string();
                let target_start =
                    str::parse::<usize>(tokens[5]).expect("failed to parse target start");
                let target_end =
                    str::parse::<usize>(tokens[6]).expect("failed to parse target end");
                let query_name = tokens[8].to_string();
                let query_start =
                    str::parse::<usize>(tokens[10]).expect("failed to parse query start");
                let query_end = str::parse::<usize>(tokens[11]).expect("failed to parse query end");
                let query_remaining =
                    str::parse::<usize>(tokens[12]).expect("failed to parse query remaining");

                let strand = match tokens[13] {
                    "0" => Strand::Forward,
                    "1" => Strand::Reverse,
                    _ => {
                        panic!()
                    }
                };

                let (query_start, query_end) = match strand {
                    Strand::Forward => (query_start, query_end),
                    Strand::Reverse => (query_end, query_start),
                    _ => unreachable!(),
                };

                let (target_seq, query_seq) = caf_str_to_digital_nucleotides(tokens[16]);
                let substitution_matrix_name = tokens[17].to_string();

                let target_id = target_name_map.insert(target_name);
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

                let query_id = query_name_map.insert(query_name);
                match strand {
                    Strand::Forward => {
                        query_lengths.insert(query_id, query_end + query_remaining);
                    }
                    Strand::Reverse => {
                        query_lengths.insert(query_id, query_start + query_remaining);
                    }
                    Strand::Unset => panic!(),
                }
                let substitution_matrix_id = substitution_matrices
                    .values
                    .iter()
                    .enumerate()
                    .find(|(_, m)| m.name == substitution_matrix_name)
                    .expect("unknown substitution matrix")
                    .0;

                target_group.alignments.push(Alignment {
                    target_seq,
                    query_seq,
                    target_start,
                    target_end,
                    query_start,
                    query_end,
                    strand,
                    id: line_num + 1,
                    query_id,
                    substitution_matrix_id,
                });
            });

        if let Some(buf) = ultra {
            let buf_reader = BufReader::new(buf);
            let ultra_json: UltraJson = serde_json::from_reader(buf_reader)?;
            ultra_json
                .repeats
                .into_iter()
                .enumerate()
                .for_each(|(idx, r)| {
                    // TODO: fix the VecMap API to better handle this
                    if !target_name_map.contains(&r.sequence_name) {
                        target_name_map.insert(r.sequence_name.clone());
                    }

                    let target_id = target_name_map.key(&r.sequence_name);

                    if let Some(group) = target_groups.get_mut(target_id) {
                        group.tandem_repeats.push(TandemRepeat {
                            // TODO: figure out if ultra uses 0- or 1-based indexing
                            id: idx + 1,
                            target_start: r.start,
                            target_end: r.start + r.length - 1,
                            consensus_pattern: r.consensus,
                            period: r.period,
                            scores: r.position_score_deltas,
                        })
                    }
                });
        }

        target_groups
            .iter_mut()
            .for_each(|g| g.alignments.sort_by(|a, b| a.id.cmp(&b.id)));

        target_groups.iter().for_each(|g| {
            debug_assert!(g
                .alignments
                .iter()
                .zip(g.alignments.iter().skip(1))
                .all(|(a, b)| a.target_start <= b.target_start));
        });

        Ok(Self {
            target_groups,
            target_name_map,
            query_name_map,
            query_lengths,
            substitution_matrices,
        })
    }
}

impl Alignment {
    #[allow(dead_code)]
    pub fn print(&self) {
        self.target_seq
            .iter()
            .for_each(|&byte| print!("{}", ALIGNMENT_ALPHABET_STR[byte as usize]));
        println!();

        self.target_seq
            .iter()
            .zip(self.query_seq.iter())
            .for_each(|(&t, &q)| if t == q { print!("|") } else { print!(" ") });
        println!();

        self.query_seq
            .iter()
            .for_each(|&byte| print!("{}", ALIGNMENT_ALPHABET_STR[byte as usize]));
        println!();
    }
}

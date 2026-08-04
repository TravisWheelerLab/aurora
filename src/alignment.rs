use anyhow::Result;
use itertools::Itertools;
use std::collections::HashMap;
use std::sync::Arc;
use std::{fmt, hash};

use serde::{ser::SerializeStruct, Serialize, Serializer};

use crate::alignment::CigarSegment::{Aligned, QueryGap, TargetGap};
use crate::alphabet::{
    NucleotideAlignmentType, NucleotideByteUtils, A_DIGITAL, C_DIGITAL, GAP_EXTEND_DIGITAL,
    GAP_OPEN_DIGITAL, G_DIGITAL, T_DIGITAL,
};
use crate::sequence_store::SequenceIndex;
use crate::substitution_matrix::SubstitutionMatrix;
use crate::util::VecMap;

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

pub struct ULEBS(Vec<u8>);

impl ULEBS {
    fn iter(&self) -> ULEBIterator<'_> {
        return ULEBIterator {
            ints: &self.0,
            offset: 0,
        };
    }
}

impl<T: Iterator<Item = u64>> From<T> for ULEBS {
    fn from(values: T) -> Self {
        let mut data: Vec<u8> = Vec::new();

        for value in values {
            let mut value = value;
            while value > 0x80 {
                data.push((value & 0x7F) as u8 | 0x80);
                value = value >> 7;
            }
            data.push((value & 0x7F) as u8);
        }

        ULEBS(data)
    }
}

pub struct ULEBIterator<'a> {
    ints: &'a [u8],
    offset: usize,
}

fn decode_next_uleb(arr: &[u8], mut offset: usize) -> Result<(u64, usize), (u64, usize)> {
    let mut result_int: u64 = 0;
    let mut shift: u8 = 0;

    let move_by = (arr.len() - offset).min(9);

    for _ in 0..move_by {
        let b = arr[offset];
        offset += 1;
        result_int |= (b as u64 & 0b01111111) << shift;
        if (b & 0x80) == 0 {
            return Result::Ok((result_int, offset));
        }
        shift += 7;
    }

    Result::Err((result_int, offset))
}

impl Iterator for ULEBIterator<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.offset >= self.ints.len() {
            return None;
        }

        let (val, next_offset) = decode_next_uleb(self.ints, self.offset).unwrap();
        self.offset = next_offset;
        Some(val)
    }
}

pub struct Cigar(ULEBS);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CigarSegment {
    Aligned(u64),
    TargetGap(u64),
    QueryGap(u64),
}

impl CigarSegment {
    pub fn count(&self) -> &u64 {
        let (Aligned(count) | TargetGap(count) | QueryGap(count)) = self;
        return count;
    }

    pub fn count_mut(&mut self) -> &mut u64 {
        let (Aligned(count) | TargetGap(count) | QueryGap(count)) = self;
        return count;
    }
}

impl From<CigarSegment> for i64 {
    fn from(value: CigarSegment) -> Self {
        match value {
            Aligned(val) | QueryGap(val) => val as i64,
            TargetGap(val) => -(val as i64),
        }
    }
}

impl Cigar {
    pub fn iter(&self) -> CigarIterator<'_> {
        CigarIterator(ULEBIterator {
            ints: &self.0 .0,
            offset: 0,
        })
    }
}

impl FromIterator<i64> for Cigar {
    fn from_iter<T: IntoIterator<Item = i64>>(iter: T) -> Self {
        let cigar = Cigar(
            iter.into_iter()
                .enumerate()
                .map(|(i, v)| {
                    if i & 1 == 0 {
                        v as u64
                    } else {
                        zig_zag_encode(v)
                    }
                })
                .into(),
        );
        assert!(cigar.0 .0.len() > 0 && cigar.0 .0.len() % 2 == 1);
        cigar
    }
}

impl FromIterator<CigarSegment> for Cigar {
    fn from_iter<T: IntoIterator<Item = CigarSegment>>(iter: T) -> Self {
        iter.into_iter().map(|v| i64::from(v)).collect()
    }
}

fn zig_zag_decode(num: u64) -> i64 {
    (num >> 1) as i64 ^ -((num & 1) as i64)
}

fn zig_zag_encode(num: i64) -> u64 {
    if num >= 0 {
        (num as u64) << 1
    } else {
        ((-num as u64) << 1) - 1
    }
}

pub struct CigarIterator<'a>(ULEBIterator<'a>);

impl Iterator for CigarIterator<'_> {
    type Item = CigarSegment;

    fn next(&mut self) -> Option<Self::Item> {
        if self.0.offset & 1 == 0 {
            self.0.next().map(|v| Aligned(v))
        } else {
            self.0.next().map(|v| {
                let z = zig_zag_decode(v);
                if z >= 0 {
                    QueryGap(v)
                } else {
                    TargetGap(v)
                }
            })
        }
    }
}

pub struct AlignmentSequence {
    pub target_seq: Arc<(usize, Vec<u8>)>,
    pub query_seq: Arc<(usize, Vec<u8>)>,
    pub cigar: Cigar,
}

impl Default for AlignmentSequence {
    fn default() -> Self {
        Self {
            target_seq: Arc::default(),
            query_seq: Arc::default(),
            cigar: Cigar(ULEBS(Vec::default())),
        }
    }
}

impl PartialEq for AlignmentSequence {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl Eq for AlignmentSequence {}

#[derive(Default, Eq)]
pub struct Alignment {
    pub sequence: AlignmentSequence,
    pub target_start: usize,
    pub target_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
    pub id: usize,
    pub query_id: usize,
    pub substitution_matrix_id: usize,
}

struct AlignmentIterator<'a, I: Iterator<Item = CigarSegment>, F: Fn(CigarSegment) -> bool> {
    cigar_seq: I,
    seq: &'a [u8],
    reverse: bool,
    gap_check: F,
    is_gap: bool,
    cigar_remaining: usize,
    cigar_steps: usize,
    offset: usize,
}

impl<'a, I: Iterator<Item = CigarSegment>, F: Fn(CigarSegment) -> bool>
    AlignmentIterator<'a, I, F>
{
    fn new(cigar: I, seq: &'a [u8], reverse: bool, gap_check: F) -> Self {
        Self {
            cigar_seq: cigar,
            seq,
            reverse,
            gap_check: gap_check,
            is_gap: false,
            cigar_remaining: 0,
            cigar_steps: 0,
            offset: 0,
        }
    }
}

impl<I: Iterator<Item = CigarSegment>, F: Fn(CigarSegment) -> bool> Iterator
    for AlignmentIterator<'_, I, F>
{
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.offset >= self.seq.len() {
            return None;
        }

        if self.cigar_remaining == 0 {
            if let Some(next_val) = self.cigar_seq.next() {
                let gap: i64 = next_val.into();
                self.cigar_remaining = gap.abs() as usize;
                self.is_gap = (self.gap_check)(next_val);
                self.cigar_steps = 0;
            }
        }
        if self.cigar_remaining == 0 {
            return None;
        }

        let char = if self.is_gap {
            if self.cigar_steps == 0 {
                GAP_OPEN_DIGITAL
            } else {
                GAP_EXTEND_DIGITAL
            }
        } else {
            let c = if self.reverse {
                self.seq[self.seq.len() - (self.offset + 1)]
            } else {
                self.seq[self.offset]
            };
            self.offset += 1;
            c
        };

        self.cigar_remaining -= 1;
        self.cigar_steps += 1;

        Some(char)
    }
}

impl Alignment {
    pub fn target_aligned_sequence(&self) -> impl Iterator<Item = u8> + '_ {
        let offset = self.sequence.target_seq.0;

        AlignmentIterator::new(
            self.sequence.cigar.iter(),
            &self.sequence.target_seq.1[self.target_start - offset..=self.target_end - offset],
            false,
            |v| matches!(v, TargetGap(_)),
        )
    }

    pub fn query_aligned_sequence(&self) -> impl Iterator<Item = u8> + '_ {
        let is_rev = matches!(self.strand, Strand::Reverse);
        let (start, end) = if is_rev {
            (self.query_end, self.query_start)
        } else {
            (self.query_start, self.query_end)
        };
        let offset = self.sequence.query_seq.0;

        AlignmentIterator::new(
            self.sequence.cigar.iter(),
            &self.sequence.target_seq.1[start - offset..=end - offset],
            is_rev,
            |v| matches!(v, QueryGap(_)),
        )
    }

    /// Compute the kimura80 score for a slice of the consensus sequence.
    pub fn kimura80(&self, query_start: usize, query_end: usize) -> f64 {
        let is_forward = match self.strand {
            Strand::Forward => true,
            Strand::Reverse => false,
            Strand::Unset => panic!("Strand is not set!"),
        };

        let mut aligned_positions: u64 = 0;

        // Count the CpG weighted transitions and transversions...
        let mut transitions10x: u64 = 0;
        let mut transversions: u64 = 0;

        let mut query_offset: usize = self.query_start;
        let mut prior_pair = (GAP_OPEN_DIGITAL, GAP_EXTEND_DIGITAL);

        let query_iter = self
            .query_aligned_sequence()
            .zip(self.target_aligned_sequence())
            .filter_map(|(q, t)| {
                let old_query_offset = query_offset;
                let old_prior_pair = prior_pair;
                if matches!(q, A_DIGITAL | C_DIGITAL | T_DIGITAL | G_DIGITAL) {
                    prior_pair = (q, t);
                }
                if matches!(q, GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL) {
                    return None;
                }

                if is_forward {
                    query_offset += 1;
                } else {
                    query_offset -= 1;
                }

                let past_start = if is_forward {
                    old_query_offset >= query_start
                } else {
                    old_query_offset <= query_start
                };

                if past_start {
                    Some((old_query_offset, old_prior_pair.0, old_prior_pair.1, q, t))
                } else {
                    None
                }
            })
            .take_while(|&val| {
                if is_forward {
                    val.0 <= query_end
                } else {
                    val.0 >= query_end
                }
            });

        for (_i, q_p, t_p, q_c, t_c) in query_iter {
            let is_cpg_group = q_p == C_DIGITAL && q_c == G_DIGITAL;
            let current_state = NucleotideAlignmentType::from_pair(q_c, t_c);
            let prior_state = NucleotideAlignmentType::from_pair(q_p, t_p);

            aligned_positions += matches!(
                current_state,
                NucleotideAlignmentType::Match
                    | NucleotideAlignmentType::Transition
                    | NucleotideAlignmentType::Transversion
            ) as u64;

            if is_cpg_group {
                match current_state {
                    NucleotideAlignmentType::Transversion => {
                        if matches!(prior_state, NucleotideAlignmentType::Transition) {
                            // Correct prior value so it's 1/10th as expected...
                            transitions10x -= 9;
                        }
                        transversions += 1;
                    }
                    NucleotideAlignmentType::Transition => {
                        match prior_state {
                            // Don't add anything, count double as a single transition...
                            NucleotideAlignmentType::Transition => {}
                            // Add 1/10th for anything else...
                            _ => {
                                transitions10x += 1;
                            }
                        }
                    }
                    _ => {
                        if matches!(prior_state, NucleotideAlignmentType::Transition) {
                            // Correct prior value so it's 1/10th as expected...
                            transitions10x -= 9;
                        }
                    }
                }
            } else {
                match current_state {
                    NucleotideAlignmentType::Transversion => {
                        transversions += 1;
                    }
                    NucleotideAlignmentType::Transition => {
                        transitions10x += 10;
                    }
                    _ => {}
                }
            }
        }

        let p = (transitions10x as f64) / ((10 * aligned_positions) as f64);
        let q = (transversions as f64) / (aligned_positions as f64);

        (-50.0 * ((1.0 - 2.0 * p - q) * (1.0 - 2.0 * q).sqrt()).ln()).abs()
    }

    pub fn ordered_query_range(&self) -> (usize, usize) {
        if matches!(self.strand, Strand::Reverse) {
            (self.query_end, self.query_start)
        } else {
            (self.query_start, self.query_end)
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

impl fmt::Debug for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: T|{}-{} Q|{}-{}",
            self.query_id, self.target_start, self.target_end, self.query_start, self.query_end
        )
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "#{}", self.id)?;
        writeln!(f, "#T {}..={}", self.target_start, self.target_end)?;
        writeln!(f, "#Q {}..={}", self.query_start, self.query_end)?;
        let mid_line: String = self
            .sequence
            .cigar
            .iter()
            .map(|seg| match seg {
                Aligned(v) => "|".repeat(v as usize),
                TargetGap(v) | QueryGap(v) => " ".repeat(v as usize),
            })
            .collect();

        writeln!(
            f,
            "{}",
            self.target_aligned_sequence()
                .collect_vec()
                .to_utf8_string()
        )?;
        writeln!(f, "{mid_line}")?;
        writeln!(
            f,
            "{}",
            self.query_aligned_sequence().collect_vec().to_utf8_string()
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

#[derive(Debug, Default)]
pub struct TandemRepeat {
    pub id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub consensus_pattern: String,
    pub period: usize,
    pub scores: Vec<f64>,
}

/// A group of alignments that share
/// the same target sequence
#[allow(dead_code)]
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
    pub target_sequences: SequenceIndex,
    pub query_sequences: SequenceIndex,
}

impl AlignmentData {
    #[allow(dead_code)]
    pub fn allocation_size(&self) -> usize {
        self.target_groups
            .iter()
            .flat_map(|g| &g.alignments)
            .map(|a| a.sequence.cigar.0 .0.capacity() + std::mem::size_of::<Alignment>())
            .sum::<usize>()
            + self.substitution_matrices.capacity() * std::mem::size_of::<SubstitutionMatrix>()
            + self.query_lengths.capacity() * std::mem::size_of::<usize>()
            + self.query_name_map.capacity() * std::mem::size_of::<String>()
            + self
                .query_name_map
                .values()
                .map(|s| s.capacity())
                .sum::<usize>()
            + self.target_name_map.capacity() * std::mem::size_of::<String>()
            + self
                .target_name_map
                .values()
                .map(|s| s.capacity())
                .sum::<usize>()
            + self.target_sequences.allocation_size()
            + self.query_sequences.allocation_size()
    }
}

impl Alignment {
    #[allow(dead_code)]
    pub fn print(&self) {
        println!("{}", self);
    }
}

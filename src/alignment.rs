use anyhow::Result;
use itertools::Itertools;
use std::collections::HashMap;
use std::mem;
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
use crate::util::{StrSliceExt, VecMap};

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
    #[allow(dead_code)]
    fn iter(&self) -> ULEBIterator<'_> {
        return ULEBIterator {
            ints: &self.0,
            offset: 0,
        };
    }
}

impl FromIterator<u64> for ULEBS {
    fn from_iter<T: IntoIterator<Item = u64>>(iter: T) -> Self {
        let mut data: Vec<u8> = Vec::new();

        for value in iter {
            let mut value = value;
            while value >= 0b1000_0000 {
                data.push((value & 0b0111_1111) as u8 | 0b1000_0000);
                value = value >> 7;
            }
            data.push((value & 0b0111_1111) as u8);
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

    let move_by = (arr.len() - offset).min(10);

    for _ in 0..move_by {
        let b = arr[offset];
        offset += 1;
        result_int |= ((b & 0b0111_1111) as u64) << shift;
        if (b & 0b1000_0000) == 0 {
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
        CigarIterator {
            inner_iter: ULEBIterator {
                ints: &self.0 .0,
                offset: 0,
            },
            return_count: 0,
        }
    }
}

impl FromIterator<i64> for Cigar {
    fn from_iter<T: IntoIterator<Item = i64>>(iter: T) -> Self {
        let mut count = 0;

        let cigar = Cigar(
            iter.into_iter()
                .enumerate()
                .map(|(i, v)| {
                    count += 1;
                    if i % 2 == 0 {
                        v.abs() as u64
                    } else {
                        zig_zag_encode(v)
                    }
                })
                .collect(),
        );
        assert!(count > 0 && count % 2 == 1);
        cigar
    }
}

impl FromIterator<CigarSegment> for Cigar {
    fn from_iter<T: IntoIterator<Item = CigarSegment>>(iter: T) -> Self {
        iter.into_iter().map(|v| i64::from(v)).collect()
    }
}

#[inline]
fn zig_zag_decode(num: u64) -> i64 {
    let is_neg = num & 1;
    let msk = !(is_neg.wrapping_sub(1));
    ((num >> 1) ^ msk) as i64
}

#[inline]
fn zig_zag_encode(num: i64) -> u64 {
    let is_neg = (num < 0) as u64;
    let msk = !(is_neg.wrapping_sub(1));
    (((num as u64) ^ msk) << 1) + is_neg
}

pub struct CigarIterator<'a> {
    inner_iter: ULEBIterator<'a>,
    return_count: usize,
}

impl Iterator for CigarIterator<'_> {
    type Item = CigarSegment;

    fn next(&mut self) -> Option<Self::Item> {
        let next_val = self.inner_iter.next().map(|v| {
            if self.return_count % 2 == 0 {
                Aligned(v)
            } else {
                let z = zig_zag_decode(v);
                let z_abs = z.abs() as u64;
                if z >= 0 {
                    QueryGap(z_abs)
                } else {
                    TargetGap(z_abs)
                }
            }
        });
        self.return_count += 1;
        next_val
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

pub struct RawSequences {
    pub target_seq: Vec<u8>,
    pub query_seq: Vec<u8>,
    pub cigar: Cigar,
}

pub fn digital_nucleotides_to_original_sequences(
    target_gapped_seq: Vec<u8>,
    query_gapped_seq: Vec<u8>,
    strand: Strand,
) -> RawSequences {
    let mut target_seq = Vec::new();
    let mut query_seq = Vec::new();
    let mut cigar = Vec::new();

    for (target_val, query_val) in target_gapped_seq.iter().zip(query_gapped_seq.iter()) {
        let next_val = match (*target_val, *query_val) {
            (GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL) => {
                panic!("Gaps in both sequences!")
            }
            (GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL, q_val) => {
                query_seq.push(q_val);
                CigarSegment::TargetGap(1)
            }
            (t_val, GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL) => {
                target_seq.push(t_val);
                CigarSegment::QueryGap(1)
            }
            (t_val, q_val) => {
                target_seq.push(t_val);
                query_seq.push(q_val);
                CigarSegment::Aligned(1)
            }
        };

        if let Some(prior_val) = cigar.last_mut() {
            if mem::discriminant(prior_val) == mem::discriminant(&next_val) {
                let count = prior_val.count_mut();
                *count += *next_val.count();
                continue;
            }
        }
        cigar.push(next_val);
    }

    if matches!(strand, Strand::Reverse) {
        query_seq.reverse();
    }

    query_seq.shrink_to_fit();
    target_seq.shrink_to_fit();
    println!("Pre-Save: {:?}", cigar);

    RawSequences {
        target_seq,
        query_seq,
        cigar: cigar.into_iter().collect(),
    }
}

impl Alignment {
    #[allow(dead_code)]
    pub fn from_str(str: &str) -> Self {
        Self::from_str_with_target_offset(str, 1)
    }

    #[allow(dead_code)]
    pub fn from_str_with_target_offset(str: &str, target_start: usize) -> Self {
        let tokens: Vec<&str> = str.split('\n').collect();

        let target = tokens[0];
        let query = tokens[1];
        assert_eq!(target.len(), query.len());

        let seqs = digital_nucleotides_to_original_sequences(
            target.to_digital_nucleotides(),
            query.to_digital_nucleotides(),
            Strand::Forward,
        );

        println!("{}", str);
        println!(
            "{:?}, {:?}, {:?}",
            seqs.target_seq.to_debug_utf8_string(),
            seqs.query_seq.to_debug_utf8_string(),
            seqs.cigar.iter().collect_vec()
        );

        let target_len = seqs.target_seq.len();
        let query_len = seqs.query_seq.len();

        Self {
            sequence: AlignmentSequence {
                target_seq: Arc::new((target_start, seqs.target_seq)),
                query_seq: Arc::new((1, seqs.query_seq)),
                cigar: seqs.cigar,
            },
            target_start: target_start,
            target_end: target_start + target_len - 1,
            query_start: 1,
            query_end: query_len,
            strand: Strand::Forward,
            id: 0,
            query_id: 0,
            substitution_matrix_id: 0,
        }
    }

    pub fn target_aligned_sequence(&self) -> impl Iterator<Item = u8> + '_ {
        let offset = self.sequence.target_seq.0;

        println!(
            "{} {} {}, {:?}",
            self.target_start,
            self.target_end,
            self.sequence.target_seq.0,
            self.sequence.target_seq.1
        );

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
            &self.sequence.query_seq.1[start - offset..=end - offset],
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::Xoshiro256PlusPlus, RngExt, SeedableRng};

    #[test]
    fn check_zig_zag_encoding() {
        for i in -200..=200 {
            assert_eq!(i, zig_zag_decode(zig_zag_encode(i)));
            if i >= 0 {
                assert_eq!((i.abs() as u64 * 2), zig_zag_encode(i));
            } else {
                assert_eq!((i.abs() as u64 * 2) - 1, zig_zag_encode(i));
            }
        }
    }

    #[test]
    fn test_uleb_encoding() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        let vals = (0..1000)
            .map(|_| rng.random_range(..=u64::MAX))
            .collect_vec();

        let ulebs: ULEBS = vals.clone().into_iter().collect();
        let vals_decoded = ulebs.iter().collect_vec();

        assert_eq!(vals, vals_decoded);
    }

    #[test]
    fn test_cigar_encoding() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        let vals = (0..1001)
            .map(|idx| {
                let v = rng.random_range(i64::MIN..=i64::MAX);
                if idx % 2 == 0 {
                    v.abs()
                } else {
                    v
                }
            })
            .collect_vec();

        let cigar: Cigar = vals.clone().into_iter().collect();

        let vals_decoded: Vec<i64> = cigar.iter().map(|v| v.into()).collect_vec();

        assert_eq!(vals, vals_decoded);
    }
}

use anyhow::Result;
use std::io::{BufRead, BufReader, Read};

use crate::alphabet::{
    ALIGNMENT_ALPHABET_STR, DASH_UTF8, FORWARD_SLASH_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
    NUCLEOTIDE_ALPHABET_UTF8, PLUS_UTF8, UTF8_TO_DIGITAL_NUCLEOTIDE,
};

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    #[default]
    Unset,
}

impl Strand {
    pub fn to_string(&self) -> String {
        match self {
            Strand::Forward => "+".to_string(),
            Strand::Reverse => "-".to_string(),
            Strand::Unset => panic!("can't convert Strand::Unset to String"),
        }
    }
}

#[derive(Default, Debug)]
pub struct Alignment {
    pub target_name: String,
    pub query_name: String,
    pub query_id: usize,
    pub target_seq: Vec<u8>,
    pub query_seq: Vec<u8>,
    pub target_start: usize,
    pub target_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
    pub score: usize,
    pub substitution_matrix: String,
    pub gap_init: f64,
    pub gap_extend: f64,
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
            // rust note: this if statement is called a match guard
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

impl Alignment {
    pub fn from_caf<R: Read>(caf: R) -> Result<(Vec<Self>, Vec<String>)> {
        let mut alignments = vec![];
        let mut query_names: Vec<String> = vec!["skip".to_string()];

        //  Robert's notes on the CAF record format:
        //      0: score - bit, raw, complexity adjusted or evalue.
        //      1: Percent Substitution - Percent of mismatched non-gap characters in the alignment
        //      2: Percent Deletion - Percent of deletion characters in alignment
        //      3: Percent Insertion - Percent of insertion characters in alignment
        //         *** THE GENOME IS THE QUERY HERE
        //      4: Query Sequence ID
        //      5: Query Start - 1-based, fully closed
        //      6: Query End - 1-based, fully closed
        //      7: Query Remaining - Remaining length of query sequence
        //         *** THE SUBJECT IS THE MODEL/TE
        //      8: Subject Sequence ID - Subject sequence is generally the TE family model for our use cases
        //      9: Subject Classification - [optional] The Dfam/RepeatMasker classification for the TE family
        //     10: Subject Start - 1 based, fully closed
        //     11: Subject End - 1 based, fully closed
        //     12: Subject Remaining - Remaining length of subject sequence
        //     13: Orientation - 0=plus_strand, 1=negative_strand
        //     14: Overlap - [optional] Overlapping annotations from RepeatMasker are flagged using this field
        //     15: Linkage_ID - [optional] RepeatMasker linkage id
        //     16: CAF encoded alignment string
        //     17: Matrix - [optional] The matrix used in scoring the alignment encoded as ##p##g.matrix or simply ##p##g
        //
        //     GAP_PARAMS = {
        //         14: {"open": -35, "ext": -6},
        //         18: {"open": -33, "ext": -5},
        //         20: {"open": -30, "ext": -5},
        //         25: {"open": -27, "ext": -5}
        //     }
        //
        // 0 199
        // 1 12.12
        // 2 0.00
        // 3 0.00
        // 4 6
        // 5 18172245
        // 6 18172277
        // 7 58603
        // 8 DF0000023
        // 9 <empty>
        // 10 1
        // 11 33
        // 12 2673
        // 13 1
        // 14 <empty>
        // 15 <empty>
        // 16 ACCT/GGA/TCT/CGTGGCCT/CGGGGGTTGGGGACCCCTG
        // 17 14p41g.matrix

        let caf_lines = BufReader::new(caf).lines();

        for line_result in caf_lines {
            let line = match line_result {
                Ok(line) => line,
                Err(_) => break,
            };

            if line.is_empty() {
                continue;
            }

            let tokens: Vec<&str> = line.split(',').collect();
            let score = str::parse::<usize>(tokens[0]).expect("failed to parse score");
            let target_name = tokens[4].to_string();
            let target_start =
                str::parse::<usize>(tokens[5]).expect("failed to parse target start");
            let target_end = str::parse::<usize>(tokens[6]).expect("failed to parse butt end");
            let query_name = tokens[8].to_string();
            let query_start = str::parse::<usize>(tokens[10]).expect("failed to parse query start");
            let query_end = str::parse::<usize>(tokens[11]).expect("failed to parse query end");
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
            let substitution_matrix = tokens[17].to_string();

            let query_id = match query_names
                .iter()
                .enumerate()
                .find(|&(_, n)| *n == query_name)
            {
                Some((id, _)) => id,
                None => {
                    query_names.push(query_name.clone());
                    query_names.len() - 1
                }
            };

            alignments.push(Alignment {
                target_name,
                query_name,
                target_seq,
                query_seq,
                query_id,
                target_start,
                target_end,
                query_start,
                query_end,
                strand,
                score,
                substitution_matrix,
                gap_init: -30.0,
                gap_extend: -5.0,
            });
        }

        Ok((alignments, query_names))
    }

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

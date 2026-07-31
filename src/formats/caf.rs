use crate::alignment::{Alignment, AlignmentData, Strand, TargetGroup};
use crate::formats::{AlignmentFormat, FormatCheck, SeekableReader};
use crate::util::{read_non_empty_lines, VecMap};
use anyhow::{anyhow, Context};
use std::collections::HashMap;

struct CAFFormat {}

impl AlignmentFormat for CAFFormat {
    fn format_check(
        primary_reader: impl SeekableReader,
        secondary_reader: Option<impl SeekableReader>,
    ) -> anyhow::Result<FormatCheck> {
        if secondary_reader.is_some() {
            return Ok(FormatCheck::Invalid(
                "CAF format doesn't require a second file!".to_string(),
            ));
        } else {
            Ok(FormatCheck::Valid)
        }
    }

    fn read(
        primary_reader: impl SeekableReader,
        secondary_reader: Option<impl SeekableReader>,
        substitution_matrices: VecMap<crate::substitution_matrix::SubstitutionMatrix>,
    ) -> anyhow::Result<crate::alignment::AlignmentData> {
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

        let mut target_groups: Vec<TargetGroup> = vec![];
        let mut target_name_map: VecMap<String> = VecMap::new();
        let mut query_name_map: VecMap<String> = VecMap::from(vec!["skip".into()]);
        let mut query_lengths: HashMap<usize, usize> = HashMap::new();
        query_lengths.insert(0, 0);

        read_non_empty_lines(primary_reader).try_for_each(|line_unchecked| {
            let (line_num, line) = line_unchecked?;

            let error_msg =
                |msg, col| move || format!("{} at line '{}', column '{}'", msg, line_num, col);

            let error_msg_str = |msg: String, col: usize| {
                move || format!("{} at line '{}', column '{}'", msg, line_num, col)
            };

            let tokens: Vec<&str> = line.split(',').collect();

            if tokens.len() < 18 {
                return Err(anyhow!(
                    "line {} does not have at least 18 columns!",
                    line_num
                ));
            }

            let target_name = tokens[4].to_string();
            let target_start = str::parse::<usize>(tokens[5])
                .with_context(error_msg("failed to parse target start", 5))?;
            let target_end = str::parse::<usize>(tokens[6])
                .with_context(error_msg("failed to parse target end", 6))?;
            let query_name = tokens[8].to_string();
            let query_start = str::parse::<usize>(tokens[10])
                .with_context(error_msg("failed to parse query start", 10))?;
            let query_end = str::parse::<usize>(tokens[11])
                .with_context(error_msg("failed to parse query end", 11))?;
            let query_remaining = str::parse::<usize>(tokens[12])
                .with_context(error_msg("failed to parse query remaining", 12))?;

            let strand = match tokens[13] {
                "0" => Ok(Strand::Forward),
                "1" => Ok(Strand::Reverse),
                v => Err(anyhow!(error_msg_str(
                    format!("invalid strand value: '{}'", v),
                    13
                )())),
            }?;

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
                .values()
                .enumerate()
                .find(|(_, m)| m.name == substitution_matrix_name)
                .with_context(error_msg_str(
                    format!("unknown substitution matrix '{}'", substitution_matrix_name),
                    17,
                ))?
                .0;

            target_group.alignments.push(Alignment {
                target_seq,
                query_seq,
                target_start,
                target_end,
                query_start,
                query_end,
                strand,
                id: 0, // We fix this later, we don't know if these are sorted yet...
                query_id,
                substitution_matrix_id,
            });

            Ok(())
        })?;

        Ok(AlignmentData {
            target_groups,
            target_name_map,
            query_name_map,
            query_lengths,
            substitution_matrices,
        })
    }

    fn name() -> &'static str {
        "CAF"
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
    target_bytes_digital.shrink_to_fit();
    query_bytes_digital.shrink_to_fit();
    (target_bytes_digital, query_bytes_digital)
}

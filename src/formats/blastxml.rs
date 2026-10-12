use crate::alignment::{
    digital_nucleotides_to_original_sequences, Alignment, AlignmentData, AlignmentSequence,
    RawSequences, Strand, TargetGroup,
};
use crate::alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE;
use crate::formats::{AlignmentFormat, FormatCheck, SeekableReader};
use crate::sequence_store::SequenceStore;
use crate::util::VecMap;
use anyhow::Context;
use quick_xml::de::from_reader;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

// When validating or checking if the file is BLAST XML, we'll only look at the top-most tags to
// save on memory and speed.

#[derive(Debug, Serialize, Deserialize)]
enum BlastOutputValidator {
    #[serde(rename = "BlastOutput")]
    BlastOutput {
        #[serde(rename = "BlastOutput_param")]
        params: EmptyStruct,
        #[serde(rename = "BlastOutput_iterations")]
        iterations: EmptyStruct,
    },
}

#[derive(Debug, Serialize, Deserialize)]
struct EmptyStruct {}

// BLAST XML FORMAT:

// NOTE: Query mapping to target is intentional, in blast the query is actually what aurora calls the target, or the genome.
// In blast, the actual query sequences are referred to as subject or hit sequences.

#[derive(Debug, Serialize, Deserialize)]
enum BlastOutput {
    #[serde(rename = "BlastOutput")]
    BlastOutputTag {
        #[serde(rename = "BlastOutput_param")]
        params: BlastOutputParam,
        #[serde(rename = "BlastOutput_iterations")]
        iterations: BlastOutputIterations,
    },
}

#[derive(Debug, Serialize, Deserialize)]
struct BlastOutputParam {
    #[serde(rename = "Parameters")]
    value: BlastParameters,
}

#[derive(Debug, Serialize, Deserialize)]
struct BlastParameters {
    #[serde(rename = "Parameters_matrix")]
    matrix: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct BlastOutputIterations {
    #[serde(rename = "Iteration")]
    value: Vec<Iteration>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Iteration {
    #[serde(rename = "Iteration_query-def")]
    target_name: String,
    #[serde(rename = "Iteration_query-len")]
    target_length: String,
    #[serde(rename = "Iteration_hits")]
    hits: IterationHits,
}

#[derive(Debug, Serialize, Deserialize)]
struct IterationHits {
    #[serde(rename = "Hit")]
    value: Vec<Hit>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Hit {
    #[serde(rename = "Hit_def")]
    query_name: String,
    #[serde(rename = "Hit_len")]
    query_length: usize,
    #[serde(rename = "Hit_hsps")]
    high_scoring_pairs: HighScoringPairs,
}

#[derive(Debug, Serialize, Deserialize)]
struct HighScoringPairs {
    #[serde(rename = "Hsp")]
    value: Vec<HighScoringPair>,
}

#[derive(Debug, Serialize, Deserialize)]
struct HighScoringPair {
    #[serde(rename = "Hsp_query-from")]
    target_start: usize,
    #[serde(rename = "Hsp_query-to")]
    target_end: usize,
    #[serde(rename = "Hsp_hit-from")]
    query_start: usize,
    #[serde(rename = "Hsp_hit-to")]
    query_end: usize,
    #[serde(rename = "Hsp_qseq")]
    target_sequence: String,
    #[serde(rename = "Hsp_hseq")]
    query_sequence: String,
}

pub struct BlastXMLFormat {}

fn parse_gapped_sequence_string(seq: &str) -> Result<Vec<u8>, anyhow::Error> {
    seq.as_bytes()
        .iter()
        .map(|v| {
            UTF8_TO_DIGITAL_NUCLEOTIDE
                .get(v)
                .copied()
                .ok_or_else(|| anyhow::anyhow!("Invalid nucleotide character: '{}'", *v as char))
        })
        .collect()
}

impl AlignmentFormat for BlastXMLFormat {
    fn format_check(
        primary_reader: &mut impl SeekableReader,
        secondary_reader: Option<&mut impl SeekableReader>,
    ) -> anyhow::Result<FormatCheck> {
        if secondary_reader.is_some() {
            return Ok(FormatCheck::Invalid(
                "BLAST_XML format doesn't require a second file!".to_string(),
            ));
        }

        let blast_output: Result<BlastOutputValidator, _> = from_reader(primary_reader);

        match blast_output {
            Ok(_) => Ok(FormatCheck::Valid),
            Err(err) => Ok(FormatCheck::Invalid(err.to_string())),
        }
    }

    fn read(
        primary_reader: &mut impl SeekableReader,
        _secondary_reader: Option<&mut impl SeekableReader>,
        substitution_matrices: crate::util::VecMap<crate::substitution_matrix::SubstitutionMatrix>,
    ) -> anyhow::Result<crate::alignment::AlignmentData> {
        let BlastOutput::BlastOutputTag { params, iterations } = from_reader(primary_reader)?;

        let mut target_groups: Vec<TargetGroup> = vec![];
        let mut target_name_map: VecMap<String> = VecMap::new();
        let mut query_name_map: VecMap<String> = vec!["skip".to_string()].into_iter().collect();
        let mut query_lengths: HashMap<usize, usize> = HashMap::new();
        let mut target_store: SequenceStore = SequenceStore::new();
        let mut query_store: SequenceStore = SequenceStore::new();
        query_lengths.insert(0, 0);

        let substitution_matrix_id = substitution_matrices
            .key_by_name(&params.value.matrix)
            .with_context(|| {
                format!(
                    "Substitution matrix '{}' not provided!",
                    params.value.matrix
                )
            })?;

        for iteration in iterations.value.iter() {
            for hit in iteration.hits.value.iter() {
                for hsp in hit.high_scoring_pairs.value.iter() {
                    let target_id = target_name_map.insert(iteration.target_name.clone());
                    let target_group = match target_groups.get_mut(target_id) {
                        Some(group) => group,
                        None => {
                            target_groups.push(TargetGroup {
                                target_id,
                                target_start: hsp.target_start,
                                target_end: hsp.target_end,
                                alignments: vec![],
                                tandem_repeats: vec![],
                            });
                            target_groups.last_mut().unwrap()
                        }
                    };

                    let strand = if hsp.query_start <= hsp.query_end {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    };

                    let RawSequences {
                        target_seq,
                        query_seq,
                        cigar,
                    } = digital_nucleotides_to_original_sequences(
                        parse_gapped_sequence_string(&hsp.target_sequence)?,
                        parse_gapped_sequence_string(&hsp.query_sequence)?,
                        strand,
                    );

                    let query_id = query_name_map.insert(hit.query_name.clone());
                    query_lengths.insert(query_id, hit.query_length);

                    if hsp.target_start > hsp.target_end {
                        return Err(anyhow::anyhow!(
                            "Alignments file has a backward target sequence. Start: {}, End: {}.",
                            hsp.target_start,
                            hsp.target_end
                        ));
                    }

                    target_store.add_sequence(target_id, hsp.target_start, &target_seq);
                    query_store.add_sequence(
                        query_id,
                        hsp.query_start.min(hsp.query_end),
                        &query_seq,
                    );

                    target_group.target_start = target_group.target_start.min(hsp.target_start);
                    target_group.target_end = target_group.target_end.max(hsp.target_end);
                    target_group.alignments.push(Alignment {
                        sequence: AlignmentSequence {
                            target_seq: Arc::default(),
                            query_seq: Arc::default(),
                            cigar,
                        }, // Run a second loop to resolve this...
                        target_start: hsp.target_start,
                        target_end: hsp.target_end,
                        query_start: hsp.query_start,
                        query_end: hsp.query_end,
                        strand,
                        id: 0, // We fix this later, we don't know if these are sorted yet...
                        query_id,
                        substitution_matrix_id,
                    });
                }
            }
        }

        let target_sequences = target_store.into_index();
        let query_sequences = query_store.into_index();

        // Get alignment references for every sequence...
        for group in target_groups.iter_mut() {
            let target_id = group.target_id;
            for al in group.alignments.iter_mut() {
                let (q_start, q_end) = al.ordered_query_range();
                let t_start = al.target_start;
                let t_end = al.target_end;

                al.sequence.target_seq = target_sequences
                    .find(target_id, t_start, t_end)
                    .ok_or(anyhow::anyhow!("Unable to find needed target sequence!"))?;
                al.sequence.query_seq = query_sequences
                    .find(al.query_id, q_start, q_end)
                    .ok_or(anyhow::anyhow!("Unable to find needed query sequence!"))?;
            }
        }

        Ok(AlignmentData {
            target_groups,
            target_name_map,
            query_name_map,
            query_lengths,
            substitution_matrices,
            target_sequences,
            query_sequences,
        })
    }

    fn name() -> &'static str {
        "BLAST_XML"
    }
}

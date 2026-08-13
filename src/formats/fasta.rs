use std::collections::HashMap;
use std::io::BufRead;

use crate::util::read_non_empty_lines;
use crate::{
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE, formats::bpaf::Timer, sequence_store::SequenceStore,
    util::VecMap,
};

pub fn parse_fasta_file(
    reader: &mut impl BufRead,
    target_names: &VecMap<String>,
    query_names: &VecMap<String>,
    target_sequences: &mut SequenceStore,
    query_sequences: &mut SequenceStore,
    mut query_lengths: Option<&mut HashMap<usize, usize>>,
) -> anyhow::Result<()> {
    let mut timer = Timer::new();

    let mut current_sequence_name = None;
    let mut current_sequence = Vec::new();

    let mut try_add_prior =
        |prior_sequence_name: Option<&String>,
         prior_sequence: &[u8],
         q_length_map: &mut Option<&mut HashMap<usize, usize>>| {
            if let Some(prior_name) = prior_sequence_name {
                if let Some(target_id) = target_names.key(prior_name) {
                    target_sequences.add_sequence(target_id, 1, prior_sequence);
                } else if let Some(query_id) = query_names.key(prior_name) {
                    query_sequences.add_sequence(query_id, 1, prior_sequence);
                    if let Some(map) = q_length_map {
                        map.insert(query_id, prior_sequence.len());
                    }
                }
            }
        };

    read_non_empty_lines(reader).try_for_each(|line_unchecked| {
        let (_, line) = line_unchecked?;
        let line_clean = line.trim();

        if let Some(new_full_name) = line_clean.strip_prefix(">") {
            timer.segment(format!("Read sequence {current_sequence_name:?}").as_ref());
            let new_name = new_full_name.split_whitespace().next().unwrap();
            timer.segment(format!("Header for sequence {new_name}").as_ref());

            try_add_prior(
                current_sequence_name.as_ref(),
                &current_sequence,
                &mut query_lengths,
            );
            timer.segment(format!("Add sequence {current_sequence_name:?} to index").as_ref());

            current_sequence.clear();
            current_sequence_name = Some(new_name.to_string());
        } else {
            if current_sequence_name.is_none() {
                return Err(anyhow::anyhow!("FASTA invalid, no name before sequence!"));
            }

            current_sequence.extend(
                line.into_bytes()
                    .iter()
                    .filter(|&&v| v != b'-' && v != b'+')
                    .filter_map(|v| UTF8_TO_DIGITAL_NUCLEOTIDE.get(v)),
            );
        }

        try_add_prior(
            current_sequence_name.as_ref(),
            &current_sequence,
            &mut query_lengths,
        );

        Ok(())
    })?;

    Ok(())
}

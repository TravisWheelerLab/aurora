use std::collections::HashMap;
use std::io::BufRead;

use itertools::Itertools;

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
         prior_sequence: &mut Vec<u8>,
         q_length_map: &mut Option<&mut HashMap<usize, usize>>| {
            let seq = prior_sequence
                .drain(..)
                .filter(|&v| v != b'-' && v != b'+')
                .filter_map(|v| UTF8_TO_DIGITAL_NUCLEOTIDE.get(&v))
                .copied()
                .collect_vec();

            if let Some(prior_name) = prior_sequence_name {
                eprintln!("Sequence '{}' length: {}", prior_name, seq.len());
                if let Some(target_id) = target_names.key(prior_name) {
                    target_sequences.add_sequence(target_id, 1, &seq);
                } else if let Some(query_id) = query_names.key(prior_name) {
                    query_sequences.add_sequence(query_id, 1, &seq);
                    if let Some(map) = q_length_map {
                        map.insert(query_id, seq.len());
                    }
                }
            }
        };

    timer.segment("Start reading FASTA");

    while reader.read_until(b'>', &mut current_sequence)? != 0 {
        timer.segment("Read till next '>'");
        let idx = if current_sequence_name.is_none() {
            0
        } else {
            if current_sequence.last() == Some(&b'>') {
                current_sequence
                    .iter()
                    .rposition(|&v| v == b'\n')
                    .unwrap_or(0)
            } else {
                current_sequence.len()
            }
        };
        let line_start = String::from_utf8(current_sequence.split_off(idx))?;
        let bracket_no_ws = line_start.trim_start();

        println!("{:?}", bracket_no_ws);

        let next_name = match bracket_no_ws {
            ">" => {
                let mut next_line = String::new();
                if reader.read_line(&mut next_line)? == 0 {
                    return Err(anyhow::anyhow!(
                        "No name after '>' in FASTA at end of file."
                    ));
                }

                let nxt_name = next_line.split_whitespace().next().unwrap();

                if nxt_name.len() == 0 {
                    return Err(anyhow::anyhow!(
                        "No name after '>' in FASTA at end of file."
                    ));
                }

                Some(nxt_name.to_string())
            }
            "" => None,
            _ => {
                return Err(anyhow::anyhow!("Invalid characters before '>' in FASTA!"));
            }
        };
        timer.segment(format!("Read next sequence name '{next_name:?}'").as_str());

        try_add_prior(
            current_sequence_name.as_ref(),
            &mut current_sequence,
            &mut query_lengths,
        );

        timer.segment(format!("Added prior sequence '{current_sequence_name:?}'",).as_str());

        current_sequence_name = next_name;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::alphabet::NucleotideByteUtils;
    use crate::sequence_store::SequenceStore;
    use crate::util::VecMap;
    use std::io::Cursor;

    #[test]
    fn can_parse_basic_fasta() -> anyhow::Result<()> {
        let fasta_strings = [
            b">seq\nATGCNATG-?+AT".to_vec(),
            b"\n\n   >seq\nATGCNATG-?+AT\n\n".to_vec(),
            b"\n\n   >  seq\n\n\nATG\nCN\nATG-?+AT\n\n".to_vec(),
        ];

        for fasta_str in fasta_strings {
            let mut file = Cursor::new(fasta_str);
            let targets: VecMap<String> = VecMap::new();
            let queries: VecMap<String> = ["seq"].iter().map(|v| v.to_string()).collect();
            let mut t_seqs = SequenceStore::new();
            let mut q_seqs = SequenceStore::new();
            let mut q_lengths = HashMap::new();

            parse_fasta_file(
                &mut file,
                &targets,
                &queries,
                &mut t_seqs,
                &mut q_seqs,
                Some(&mut q_lengths),
            )?;

            assert_eq!(q_seqs.sequence_count(), 1);
            assert_eq!(t_seqs.sequence_count(), 0);
            assert_eq!(q_lengths.get(&0), Some(&10));
            assert_eq!(
                q_seqs
                    .into_index()
                    .find(0, 1, 10)
                    .unwrap()
                    .1
                    .to_debug_utf8_string(),
                "ATGCNATGAT"
            );
        }

        Ok(())
    }
}

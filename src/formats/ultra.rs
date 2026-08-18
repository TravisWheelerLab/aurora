use crate::alignment::{AlignmentData, TandemRepeat};
use anyhow::{Context, Result};
use serde::Deserialize;
use std::io::{BufReader, Read};

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

pub fn load_ultra_file(alignment_data: &mut AlignmentData, ultra_file: impl Read) -> Result<()> {
    let buf_reader = BufReader::new(ultra_file);
    let ultra_json: UltraJson =
        serde_json::from_reader(buf_reader).context("Failed to load provided ultra file.")?;

    ultra_json.repeats.into_iter().for_each(|r| {
        let target_id = alignment_data
            .target_name_map
            .insert(r.sequence_name.clone());

        if let Some(group) = alignment_data.target_groups.get_mut(target_id) {
            group.tandem_repeats.push(TandemRepeat {
                // TODO: figure out if ultra uses 0- or 1-based indexing
                id: 0, // The id is set durring normalization...
                target_start: r.start,
                target_end: r.start + r.length - 1,
                consensus_pattern: r.consensus,
                period: r.period,
                scores: r.position_score_deltas,
            });

            group.target_start = group.target_start.min(r.start);
            group.target_end = group.target_end.max(r.start + r.length - 1);
        }
    });

    Ok(())
}

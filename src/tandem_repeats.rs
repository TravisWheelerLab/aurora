use std::io::{BufReader, Read};

use anyhow::Result;
use itertools::Itertools;
use serde::Deserialize;

use crate::alignment::VecMap;

#[derive(Debug)]
pub struct TandemRepeat {
    pub target_id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub consensus_pattern: String,
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
    pub position_score_deltas: Vec<f64>,
}

impl TandemRepeat {
    pub fn from_ultra_json<R: Read>(
        data: R,
        target_name_map: &VecMap<String>,
    ) -> Result<Vec<Self>> {
        let buf_reader = BufReader::new(data);
        let ultra_json: UltraJson = serde_json::from_reader(buf_reader).unwrap();

        let repeats = ultra_json
            .repeats
            .into_iter()
            .map(|r| TandemRepeat {
                target_id: target_name_map.key(&r.sequence_name),
                target_start: r.start,
                target_end: r.start + r.length,
                consensus_pattern: r.consensus,
                scores: r.position_score_deltas,
            })
            .collect_vec();

        repeats.iter().for_each(|r| println!("{:?}", r));

        Ok(repeats)
    }
}

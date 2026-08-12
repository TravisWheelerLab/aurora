use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

pub struct SequenceStore {
    sequences: HashMap<usize, BTreeMap<usize, Vec<u8>>>,
}

pub struct SequenceIndex {
    sequences: HashMap<usize, Vec<Arc<(usize, Vec<u8>)>>>,
}

impl SequenceStore {
    pub fn new() -> Self {
        Self {
            sequences: HashMap::new(),
        }
    }

    pub fn add_sequence(&mut self, sequence_id: usize, start: usize, sequence: &[u8]) {
        if sequence.len() == 0 {
            return;
        }
        let end = start + sequence.len() - 1;
        let entry = self
            .sequences
            .entry(sequence_id)
            .or_insert_with(|| BTreeMap::new());

        let mut overlaps: Vec<(usize, usize)> = entry
            .range(..=start)
            .next_back()
            .map(|v| (*v.0, v.0 + v.1.len() - 1))
            .filter(|v| v.1 >= start)
            .iter()
            .copied()
            .collect();

        overlaps.extend(
            entry
                .range((start + 1)..=end)
                .map(|v| (*v.0, v.0 + v.1.len() - 1)),
        );

        let new_start = overlaps.first().map(|v| v.0).unwrap_or(start).min(start);
        let new_end = overlaps.last().map(|v| v.1).unwrap_or(end).max(end);

        let mut new_seq = vec![0_u8; new_end - new_start + 1];

        new_seq[start - new_start..=end - new_start].copy_from_slice(sequence);

        for (start, end) in overlaps.iter() {
            if let Some(seq) = entry.remove(start) {
                new_seq[start - new_start..=end - new_start].copy_from_slice(&seq);
            }
        }

        entry.insert(new_start, new_seq);
    }

    pub fn sequence_count(&self) -> usize {
        self.sequences.len()
    }

    pub fn into_index(self) -> SequenceIndex {
        SequenceIndex {
            sequences: self
                .sequences
                .into_iter()
                .map(|(id, tree)| (id, tree.into_iter().map(|v| Arc::new(v)).collect()))
                .collect(),
        }
    }
}

impl SequenceIndex {
    pub fn find(
        &self,
        sequence_id: usize,
        start: usize,
        end: usize,
    ) -> Option<Arc<(usize, Vec<u8>)>> {
        assert!(start <= end, "Invalid inclusive range!");

        if let Some(sorted_seq_list) = self.sequences.get(&sequence_id) {
            if let Some(selected_seq) = sorted_seq_list
                [..sorted_seq_list.partition_point(|v| v.0 <= start)]
                .iter()
                .next_back()
            {
                let ent_start = selected_seq.0;
                let seq = &selected_seq.1;
                if end <= (ent_start + seq.len() - 1) {
                    return Some(selected_seq.clone());
                }
            }
        }

        None
    }

    pub fn allocation_size(&self) -> usize {
        self.sequences.iter().map(|v| 16 + v.1.capacity()).sum()
    }
}

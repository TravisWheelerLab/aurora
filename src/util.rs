use std::{
    collections::HashMap,
    io::{self},
};

use itertools::Itertools;
use thiserror::Error;

use crate::{
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    balanced_tree::{AVLIndexSet, SetInsert},
};

/// A simple Vec-based map that facilitates mapping
/// between usize keys and type <T> values
/// use avl index tree for O(n log n) value based lookup.
/// key lookup is constant.
pub struct VecMap<T: std::cmp::Ord> {
    values: Vec<T>,
    tree: AVLIndexSet<usize>,
}

impl VecMap<String> {
    // Maintains all seqeunces from vector in-order instead
    // appending an index to the end of duplicate entries.
    pub fn from_vec_raw(values: Vec<String>) -> Self {
        let mut new_self = VecMap::new();
        let mut counts: HashMap<&String, usize> = HashMap::new();

        for val in values.iter() {
            let entry = counts.entry(val).or_insert(0);
            new_self.insert(format!(
                "{}{}",
                val,
                if *entry == 0 {
                    "".to_string()
                } else {
                    entry.to_string()
                }
            ));
            *entry += 1;
        }
        new_self
    }
}

impl<I: std::cmp::Ord> FromIterator<I> for VecMap<I> {
    fn from_iter<T: IntoIterator<Item = I>>(iter: T) -> Self {
        let mut new_self = VecMap::new();
        for val in iter {
            new_self.insert(val);
        }
        new_self
    }
}

impl<T: std::cmp::Ord> VecMap<T> {
    pub fn new() -> Self {
        Self {
            values: vec![],
            tree: AVLIndexSet::new(),
        }
    }

    pub fn values(&self) -> std::slice::Iter<'_, T> {
        self.values.iter()
    }

    /// Inserts the value and returns the key. If the
    /// value was already in the VecMap, return the key.
    pub fn insert(&mut self, value: T) -> usize {
        match self.tree.add(|idx| self.values[idx].cmp(&value)).unwrap() {
            SetInsert::Found(idx) => idx,
            SetInsert::New(new_idx) => {
                self.values.push(value);
                new_idx
            }
        }
    }

    pub fn get(&self, key: usize) -> &T {
        debug_assert!(key < self.values.len(), "invalid key: {key}");
        &self.values[key]
    }

    #[allow(dead_code)]
    pub fn contains(&self, value: &T) -> bool {
        self.key(value).is_some()
    }

    pub fn size(&self) -> usize {
        self.values.len()
    }

    pub fn capacity(&self) -> usize {
        self.values.capacity()
    }

    /// Get the key associated with the value.
    pub fn key(&self, value: &T) -> Option<usize> {
        self.tree.search(|idx| self.values[idx].cmp(value))
    }
}

impl<T: std::cmp::Ord + std::fmt::Debug> std::fmt::Debug for VecMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        (0..self.size()).for_each(|key| {
            writeln!(f, "{key}: {:?}", self.values[key]).unwrap();
        });
        Ok(())
    }
}

#[derive(Error, Debug)]
#[error("unknown byte: {0}")]
pub struct InvalidByte(u8);

pub trait StrSliceExt {
    fn try_to_digital_nucleotides(self) -> Result<Vec<u8>, InvalidByte>;
}

impl StrSliceExt for &str {
    fn try_to_digital_nucleotides(self) -> Result<Vec<u8>, InvalidByte> {
        self.as_bytes()
            .iter()
            .map(|byte| {
                UTF8_TO_DIGITAL_NUCLEOTIDE
                    .get(byte)
                    .copied()
                    .ok_or_else(|| InvalidByte(*byte))
            })
            .collect()
    }
}

pub fn read_non_empty_lines(
    readable: impl io::BufRead,
) -> impl Iterator<Item = io::Result<(usize, String)>> {
    readable
        .lines()
        .enumerate()
        .map(|(lineno, v)| match v {
            Ok(s) => Ok((lineno + 1, s)),
            Err(e) => Err(e),
        })
        .filter_ok(|(_lineno, s)| !s.trim().is_empty())
}

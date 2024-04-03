use crate::{
    alignment::{Alignment, AlignmentData, TandemRepeat, TargetGroup},
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    chunks::ProximityGroup,
    substitution_matrix::SubstitutionMatrix,
};

/// This is a silly function that can extract a type and return it's size in bytes.
///
/// The alternative would be to use `std::mem::size_of_val(&val)`, but that won't
/// always work in every case. For example, say we have a `Vec<T>` of length `0`.
/// We would not be able to get the size of T because we can't actually produce
/// a reference to an element in the vector, since it's empty.
///
/// Instead, this function expects a typed closure which is never actually invoked,
/// but we are able to extract the type we want and call `std::mem::size_of()` on it.
fn get_size_of_return_type<F, T, U>(_f: F) -> usize
where
    F: FnOnce(T) -> U,
{
    std::mem::size_of::<U>()
}

pub trait AllocationSize: Sized {
    fn heap_size(&self) -> usize;
    fn total_size(&self) -> usize {
        self.heap_size() + std::mem::size_of::<Self>()
    }
}

impl<T> AllocationSize for Vec<T>
where
    T: AllocationSize,
{
    fn heap_size(&self) -> usize {
        self.iter().map(|x| x.total_size()).sum()
    }
}

impl AllocationSize for usize {
    fn heap_size(&self) -> usize {
        0
    }
}

impl AllocationSize for u8 {
    fn heap_size(&self) -> usize {
        0
    }
}

impl AllocationSize for f32 {
    fn heap_size(&self) -> usize {
        0
    }
}

impl AllocationSize for f64 {
    fn heap_size(&self) -> usize {
        0
    }
}

impl<K, V> AllocationSize for std::collections::HashMap<K, V>
where
    K: AllocationSize,
    V: AllocationSize,
{
    fn heap_size(&self) -> usize {
        self.values().map(|v| v.total_size()).sum::<usize>()
            + self.keys().map(|k| k.total_size()).sum::<usize>()
    }
}

impl AllocationSize for String {
    fn heap_size(&self) -> usize {
        self.capacity() * 8
    }
}

impl AllocationSize for Alignment {
    fn heap_size(&self) -> usize {
        let size = get_size_of_return_type(|x: Self| x.query_seq[0]);
        self.query_seq.len() * size + self.target_seq.len() * size
    }
}

impl AllocationSize for TandemRepeat {
    fn heap_size(&self) -> usize {
        let size = get_size_of_return_type(|x: Self| x.scores[0]);
        self.scores.len() * size + self.consensus_pattern.heap_size()
    }
}

impl AllocationSize for TargetGroup {
    fn heap_size(&self) -> usize {
        self.alignments.heap_size() + self.tandem_repeats.heap_size()
    }
}

impl<T> AllocationSize for VecMap<T>
where
    T: AllocationSize + PartialEq,
{
    fn heap_size(&self) -> usize {
        self.values().map(|v| v.total_size()).sum()
    }
}

impl AllocationSize for SubstitutionMatrix {
    fn heap_size(&self) -> usize {
        0
    }
}

impl AllocationSize for AlignmentData {
    fn heap_size(&self) -> usize {
        self.target_groups.heap_size()
            + self.target_name_map.heap_size()
            + self.query_name_map.heap_size()
            + self.query_lengths.heap_size()
            + self.substitution_matrices.heap_size()
    }
}

impl<'a> AllocationSize for ProximityGroup<'a> {
    fn heap_size(&self) -> usize {
        0
    }
}

/// A simple Vec-based map that facilitates mapping
/// between usize keys and type <T> values
pub struct VecMap<T: std::cmp::PartialEq> {
    values: Vec<T>,
}

impl<T: std::cmp::PartialEq> VecMap<T> {
    pub fn new() -> Self {
        Self { values: vec![] }
    }

    pub fn from(values: Vec<T>) -> Self {
        Self { values }
    }

    pub fn values(&self) -> std::slice::Iter<T> {
        self.values.iter()
    }

    /// Inserts the value and returns the key. If the
    /// value was already in the VecMap, return the key.
    pub fn insert(&mut self, value: T) -> usize {
        match self.contains(&value) {
            true => self.key(&value),
            false => {
                self.values.push(value);
                self.values.len() - 1
            }
        }
    }

    pub fn get(&self, key: usize) -> &T {
        debug_assert!(key < self.values.len(), "invalid key: {key}");
        &self.values[key]
    }

    pub fn contains(&self, value: &T) -> bool {
        self.values.contains(value)
    }

    pub fn size(&self) -> usize {
        self.values.len()
    }

    /// Get the key associated with the value.
    /// This panics if the value is not in the VecMap.
    pub fn key(&self, value: &T) -> usize {
        self.values
            .iter()
            .enumerate()
            .find(|(_, n)| *n == value)
            .expect("key not found")
            .0
    }
}

impl<T: std::cmp::PartialEq + std::fmt::Debug> std::fmt::Debug for VecMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        (0..self.size()).for_each(|key| {
            writeln!(f, "{key}: {:?}", self.values[key]).unwrap();
        });
        Ok(())
    }
}

pub trait StrSliceExt {
    fn to_digital_nucleotides(self) -> Vec<u8>;
}

impl StrSliceExt for &str {
    fn to_digital_nucleotides(self) -> Vec<u8> {
        self.as_bytes()
            .iter()
            .map(|byte| {
                *UTF8_TO_DIGITAL_NUCLEOTIDE
                    .get(byte)
                    .unwrap_or_else(|| panic!("unknown byte: {byte}"))
            })
            .collect()
    }
}

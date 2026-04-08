use crate::alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE;

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

    pub fn values(&self) -> std::slice::Iter<'_, T> {
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

    pub fn capacity(&self) -> usize {
        self.values.capacity()
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

#[derive(Default)]
pub struct ULEBS(Vec<u8>);

impl ULEBS {
    #[allow(dead_code)]
    fn iter(&self) -> impl Iterator<Item = u64> + '_ {
        return ULEBIterator {
            ints: self.0.iter().copied(),
        }
        .map(|v| v.expect("Invalid ULEB found!"));
    }
}

impl FromIterator<u64> for ULEBS {
    fn from_iter<T: IntoIterator<Item = u64>>(iter: T) -> Self {
        let mut data: Vec<u8> = Vec::new();

        for value in iter {
            let mut value = value;
            while value >= 0b1000_0000 {
                data.push((value & 0b0111_1111) as u8 | 0b1000_0000);
                value = value >> 7;
            }
            data.push((value & 0b0111_1111) as u8);
        }

        ULEBS(data)
    }
}

pub struct ULEBIterator<T: Iterator<Item = u8>> {
    ints: T,
}

fn decode_next_uleb(bytes: &mut impl Iterator<Item = u8>) -> Result<(u64, usize), (u64, usize)> {
    let mut result_int: u64 = 0;
    let mut shift: u8 = 0;
    let mut bytes_taken = 0;

    for _ in 0..10 {
        if let Some(b) = bytes.next() {
            bytes_taken += 1;
            result_int |= ((b & 0b0111_1111) as u64) << shift;
            if (b & 0b1000_0000) == 0 {
                return Result::Ok((result_int, bytes_taken));
            }
            shift += 7;
        } else {
            break;
        }
    }

    Result::Err((result_int, bytes_taken))
}

impl<T: Iterator<Item = u8>> Iterator for ULEBIterator<T> {
    type Item = Result<u64, u64>;

    fn next(&mut self) -> Option<Self::Item> {
        match decode_next_uleb(&mut self.ints) {
            Result::Ok((val, _)) => Some(Ok(val)),
            Result::Err((val, bytes_taken)) => {
                if bytes_taken == 0 {
                    None
                } else {
                    Some(Err(val))
                }
            }
        }
    }
}

#[derive(Default)]
pub struct Cigar(ULEBS);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CigarSegment {
    Aligned(u64),
    TargetGap(u64),
    QueryGap(u64),
}

impl CigarSegment {
    pub fn count(&self) -> &u64 {
        let (Self::Aligned(count) | Self::TargetGap(count) | Self::QueryGap(count)) = self;
        return count;
    }

    pub fn count_mut(&mut self) -> &mut u64 {
        let (Self::Aligned(count) | Self::TargetGap(count) | Self::QueryGap(count)) = self;
        return count;
    }
}

impl From<CigarSegment> for i64 {
    fn from(value: CigarSegment) -> Self {
        match value {
            CigarSegment::Aligned(val) | CigarSegment::QueryGap(val) => val as i64,
            CigarSegment::TargetGap(val) => -(val as i64),
        }
    }
}

impl Cigar {
    pub fn iter(&self) -> impl Iterator<Item = CigarSegment> + '_ {
        CigarIterator {
            inner_iter: ULEBIterator {
                ints: self.0 .0.iter().copied(),
            },
            return_count: 0,
        }
        .map(|v| v.expect("Invalid ULEB found in Cigar string!"))
    }

    pub fn capacity(&self) -> usize {
        self.0 .0.capacity()
    }
}

impl FromIterator<i64> for Cigar {
    fn from_iter<T: IntoIterator<Item = i64>>(iter: T) -> Self {
        let mut count = 0;

        let cigar = Cigar(
            iter.into_iter()
                .enumerate()
                .map(|(i, v)| {
                    count += 1;
                    if i % 2 == 0 {
                        v.abs() as u64
                    } else {
                        zig_zag_encode(v)
                    }
                })
                .collect(),
        );
        assert!(count > 0 && count % 2 == 1);
        cigar
    }
}

impl FromIterator<CigarSegment> for Cigar {
    fn from_iter<T: IntoIterator<Item = CigarSegment>>(iter: T) -> Self {
        iter.into_iter()
            .enumerate()
            .map(|(i, v)| {
                let valid_offset =
                    matches!(v, CigarSegment::QueryGap(_) | CigarSegment::TargetGap(_)) as usize;
                debug_assert_eq!(i % 2, valid_offset);
                i64::from(v)
            })
            .collect()
    }
}

#[inline]
pub fn zig_zag_decode(num: u64) -> i64 {
    let is_neg = num & 1;
    let msk = !(is_neg.wrapping_sub(1));
    ((num >> 1) ^ msk) as i64
}

#[inline]
pub fn zig_zag_encode(num: i64) -> u64 {
    let is_neg = (num < 0) as u64;
    let msk = !(is_neg.wrapping_sub(1));
    (((num as u64) ^ msk) << 1) + is_neg
}

pub struct CigarIterator<T: Iterator<Item = Result<u64, u64>>> {
    inner_iter: T,
    return_count: usize,
}

impl<T: Iterator<Item = Result<u64, u64>>> Iterator for CigarIterator<T> {
    type Item = Result<CigarSegment, CigarSegment>;

    fn next(&mut self) -> Option<Self::Item> {
        let into_val = |v: u64| -> CigarSegment {
            if self.return_count % 2 == 0 {
                CigarSegment::Aligned(v)
            } else {
                let z = zig_zag_decode(v);
                let z_abs = z.abs() as u64;
                if z >= 0 {
                    CigarSegment::QueryGap(z_abs)
                } else {
                    CigarSegment::TargetGap(z_abs)
                }
            }
        };

        let next_val = self
            .inner_iter
            .next()
            .map(|v| v.map(into_val).map_err(into_val));
        self.return_count += 1;
        next_val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use rand::{rngs::Xoshiro256PlusPlus, RngExt, SeedableRng};

    #[test]
    fn check_zig_zag_encoding() {
        for i in -200..=200 {
            assert_eq!(i, zig_zag_decode(zig_zag_encode(i)));
            if i >= 0 {
                assert_eq!((i.abs() as u64 * 2), zig_zag_encode(i));
            } else {
                assert_eq!((i.abs() as u64 * 2) - 1, zig_zag_encode(i));
            }
        }
    }

    #[test]
    fn test_uleb_encoding() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        let vals = (0..1000)
            .map(|_| rng.random_range(..=u64::MAX))
            .collect_vec();

        let ulebs: ULEBS = vals.clone().into_iter().collect();
        let vals_decoded = ulebs.iter().collect_vec();

        assert_eq!(vals, vals_decoded);
    }

    #[test]
    fn test_cigar_encoding() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        let vals = (0..1001)
            .map(|idx| {
                let v = rng.random_range(i64::MIN..=i64::MAX);
                if idx % 2 == 0 {
                    v.abs()
                } else {
                    v
                }
            })
            .collect_vec();

        let cigar: Cigar = vals.clone().into_iter().collect();

        let vals_decoded: Vec<i64> = cigar.iter().map(|v| v.into()).collect_vec();

        assert_eq!(vals, vals_decoded);
    }
}

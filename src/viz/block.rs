use serde::{Serialize, Serializer};

use crate::{alignment::Strand, results::Annotation};

use super::BedRecord;

///
///
///
///
pub struct Block {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub query_length: Option<i32>,
}

impl std::fmt::Display for Block {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.query_length {
            Some(l) => write!(f, "{},{},{},{}", self.id, self.start, self.end, l),
            None => write!(f, "{},{},{}", self.id, self.start, self.end),
        }
    }
}

impl Serialize for Block {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BlockGroup {
    pub id: String,
    pub visual_start: usize,
    pub visual_end: usize,
    pub align_start: usize,
    pub align_end: usize,
    pub strand: Strand,
    pub query: String,
    pub target: String,
    pub left: Block,
    pub right: Block,
    pub aligned: Vec<Block>,
    pub inner: Vec<Block>,
}

impl BlockGroup {
    pub fn from_joined_annotations(joins: &mut [&Annotation]) -> Self {
        joins.sort_by_key(|a| a.target_start);

        let first = joins.first().unwrap();
        let last = joins.last().unwrap();

        let mut id_cnt = 0usize;
        let mut id_fn = || {
            let id = format!("{}-{}-{}", first.region_id, first.join_id, id_cnt);
            id_cnt += 1;
            id
        };

        let mut aligned = vec![];
        let mut inner = vec![];
        (0..joins.len() - 1)
            .map(|idx| (joins[idx], joins[idx + 1]))
            .for_each(|(a, b)| {
                aligned.push(Block {
                    id: id_fn(),
                    start: a.target_start,
                    end: a.target_end,
                    query_length: None,
                });
                inner.push(Block {
                    id: id_fn(),
                    start: a.target_end,
                    end: b.target_start,
                    query_length: Some(b.query_start as i32 - a.query_end as i32 + 1),
                });
            });

        aligned.push(Block {
            id: id_fn(),
            start: last.target_start,
            end: last.target_end,
            query_length: None,
        });

        let align_start = first.target_start;
        let align_end = last.target_end;

        // TODO: off by one?
        let visual_start = first.target_start - first.query_start;
        // TODO: need model length information to get this
        let visual_end = last.target_end + 0;

        let left = Block {
            id: id_fn(),
            start: visual_start,
            end: first.target_start,
            query_length: None,
        };

        let right = Block {
            id: id_fn(),
            start: last.target_end,
            // TODO: need model length information to get this
            end: last.target_end,
            query_length: None,
        };

        Self {
            id: format!("{}-{}", first.region_id, first.join_id),
            visual_start,
            visual_end,
            align_start,
            align_end,
            strand: first.strand,
            query: first.query_name.clone(),
            target: first.target_name.clone(),
            left,
            right,
            aligned,
            inner,
        }
    }

    pub fn from_bed_record(bed: &BedRecord) -> Self {
        let mut id_cnt = 0usize;
        let mut id_fn = || {
            let id = format!("bed-{}-{}", bed.id, id_cnt);
            id_cnt += 1;
            id
        };
        let mut aligned = vec![];
        let mut inner = vec![];

        (1..bed.block_count - 1)
            .map(|b| {
                (
                    bed.block_sizes[b],
                    bed.block_starts[b],
                    bed.block_starts[b + 1],
                )
            })
            .for_each(|(size, start, next_start)| {
                // aligned blocks have a positive start
                if start >= 0 {
                    // start is an offset from chrom_start
                    // +1 (since starts are 0-based)
                    let block_start = bed.chrom_start + start as usize + 1;
                    aligned.push(Block {
                        id: id_fn(),
                        start: block_start,
                        // -1 to include <size> positions in the length
                        end: block_start + size as usize - 1,
                        query_length: None,
                    });
                }
                // unaligned blocks have -1 start
                else {
                    let last_aligned = aligned.last().expect("aligned group is empty");
                    inner.push(Block {
                        id: id_fn(),
                        // this unaligned block starts 1 position
                        // after the previous aligned block
                        start: last_aligned.end + 1,
                        // it also ends 1 before the start of the next aligned block
                        // don't need to -1 (since starts are 0-based)
                        end: bed.chrom_start + next_start as usize,
                        // the size field encodes the length of the model
                        // that is projected into the inner unaligned block
                        query_length: Some(size),
                    });
                }
            });

        let left = Block {
            id: id_fn(),
            // +1 because starts are 0-based
            start: bed.chrom_start + 1,
            // don't need to -1 (since starts are 0-based)
            end: bed.thick_start,
            query_length: None,
        };

        let right = Block {
            id: id_fn(),
            // +1 since the unaligned starts 1 after
            // the thick_end (which is already base-1)
            start: bed.thick_end + 1,
            // the last unaligned block ends at chrom_end
            // don't need to +1 (since ends are 1-based)
            end: bed.chrom_end,
            query_length: None,
        };

        Self {
            id: format!("bed-{}", bed.id),
            visual_start: bed.chrom_start,
            visual_end: bed.chrom_end,
            align_start: bed.thick_start,
            align_end: bed.thick_end,
            strand: bed.strand,
            query: bed.name.clone(),
            target: bed.chrom.clone(),
            left,
            right,
            aligned,
            inner,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{alignment::Strand, viz::bed::BedRecord};

    use super::BlockGroup;

    #[test]
    fn test_block_group_from_bed() {
        let bed = BedRecord {
            score: 0,
            chrom: "".to_string(),
            name: "L1MC5a#LINE/L1".to_string(),
            strand: Strand::Reverse,
            reserved: 0,
            id: 0,
            description: "".to_string(),
            chrom_start: 9_122,
            chrom_end: 17_126,
            thick_start: 11_504,
            thick_end: 11_675,
            block_count: 3,
            block_sizes: vec![2_382, 171, 5_451],
            block_starts: vec![-1, 2_382, -1],
        };

        // correct under 1-based indexing:
        //    unaligned:   9,123..11,504    11,504 -  9,123 + 1 = 2,382
        //    aligned:    11,505..11,675    11,675 - 11,505 + 1 =   171
        //    unaligned:  11,676..17,126    17,126 - 11,676 + 1 = 5,451
        let block_group = BlockGroup::from_bed_record(&bed);
        let aligned = block_group.aligned.first().unwrap();

        assert_eq!(block_group.left.start, 9_123);
        assert_eq!(block_group.left.end, 11_504);

        assert_eq!(aligned.start, 11_505);
        assert_eq!(aligned.end, 11_675);

        assert_eq!(block_group.right.start, 11_676);
        assert_eq!(block_group.right.end, 17_126);

        let bed = BedRecord {
            score: 0,
            chrom: "".to_string(),
            name: "L1P4a#LINE/L1".to_string(),
            strand: Strand::Reverse,
            reserved: 0,
            id: 0,
            description: "".to_string(),
            chrom_start: 18_184_829,
            chrom_end: 18_191_691,
            thick_start: 18_188_369,
            thick_end: 18_189_819,
            block_count: 7,
            block_sizes: vec![3_540, 76, 0, 959, 0, 14, 1_872],
            block_starts: vec![-1, 3_540, -1, 3_724, -1, 4_976, -1],
        };

        let block_group = BlockGroup::from_bed_record(&bed);

        // left:  18_184_830..18_188_369
        assert_eq!(block_group.left.start, 18_184_830);
        assert_eq!(block_group.left.end, 18_188_369);

        // start: 18,184,830 + 3,540     = 18,188,370
        // end:   18,188,370 +    76 - 1 = 18,188,445
        // len:   18,188,445 - 18,188,370 + 1 = 76
        let aligned_1 = &block_group.aligned[0];
        assert_eq!(aligned_1.start, 18_188_370);
        assert_eq!(aligned_1.end, 18_188_445);

        // start: 18,188,445 + 1         = 18,188,446
        // end:   18,184,830 + 3,724 - 1 = 18,188,554
        // query len: 0
        let inner_1 = &block_group.inner[0];
        assert_eq!(inner_1.start, 18_188_446);
        assert_eq!(inner_1.end, 18_188_553);
        assert_eq!(inner_1.query_length.unwrap(), 0);

        // start: 18,184,830 + 3,724     = 18,188,554
        // end:   18,188,554 +   959 - 1 = 18,189,512
        // len:   18,189,512 - 18,188,554 + 1 = 959
        let aligned_2 = &block_group.aligned[1];
        assert_eq!(aligned_2.start, 18_188_554);
        assert_eq!(aligned_2.end, 18_189_512);

        // start: 18,189,512 + 1         = 18,189,513
        // end:   18,184,830 + 4,976 - 1 = 18,189,805
        // query len: 0
        let inner_2 = &block_group.inner[1];
        assert_eq!(inner_2.start, 18_189_513);
        assert_eq!(inner_2.end, 18_189_805);
        assert_eq!(inner_2.query_length.unwrap(), 0);

        // start: 18,184,830 + 4,976     = 18,189,806
        // end:   18,189,806 +    14 - 1 = 18,189,819
        // len:   18,189,819 - 18,189,806 + 1 = 14
        let aligned_3 = &block_group.aligned[2];
        assert_eq!(aligned_3.start, 18_189_806);
        assert_eq!(aligned_3.end, 18_189_819);

        assert_eq!(block_group.right.start, 18_189_820);
        assert_eq!(block_group.right.end, 18_191_691);
    }
}

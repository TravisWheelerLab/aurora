use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use itertools::Itertools;
use serde::{Serialize, Serializer};

use crate::{
    alignment::{Alignment, Strand},
    alphabet::{NucleotideByteUtils, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, PAD_DIGITAL},
    matrix::MatrixDef,
    results::Annotation,
};

///
///
///
///

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AuroraSodaData {
    target_start: usize,
    target_end: usize,
    target_seq: String,
    aurora_ann: Vec<AuroraAnnotationGroup>,
    reference_ann: Vec<String>,
    // alignments: ,
}

impl AuroraSodaData {
    pub fn new(
        matrix_def: &MatrixDef,
        alignments: &[Alignment],
        annotations: &[Annotation],
        reference_bed_path: impl AsRef<Path>,
        names_path: impl AsRef<Path>,
    ) -> Self {
        let target_start = matrix_def.target_start;
        let target_length = matrix_def.num_cols;

        let mut target_seq_digital_bytes = vec![PAD_DIGITAL; target_length];
        alignments
            .iter()
            .flat_map(|a| {
                a.target_seq
                    .iter()
                    .filter(|&&c| c != GAP_OPEN_DIGITAL && c != GAP_EXTEND_DIGITAL)
                    .enumerate()
                    .map(|(idx, c)| (idx + a.target_start - target_start, c))
            })
            .for_each(|(col_idx, &char)| target_seq_digital_bytes[col_idx] = char);

        let target_seq = target_seq_digital_bytes.into_utf8_string();
        let unique_join_ids: Vec<usize> = annotations.iter().map(|a| a.join_id).unique().collect();

        let aurora_annotations: Vec<AuroraAnnotationGroup> = unique_join_ids
            .iter()
            .map(|&id| {
                AuroraAnnotationGroup::from_joined_annotations(
                    &mut annotations
                        .iter()
                        .filter(|&a| a.join_id == id)
                        .collect::<Vec<&Annotation>>(),
                )
            })
            .collect();

        let reference_bed_file = File::open(reference_bed_path).expect("failed to open");

        Self {
            target_start: aurora_annotations
                .iter()
                .map(|a| a.visual_start)
                .min()
                .unwrap(),
            target_end: aurora_annotations
                .iter()
                .map(|a| a.visual_end)
                .max()
                .unwrap(),
            target_seq,
            aurora_ann: aurora_annotations,
            reference_ann: vec![],
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AuroraAnnotationGroup {
    id: String,
    visual_start: usize,
    visual_end: usize,
    align_start: usize,
    align_end: usize,
    strand: Strand,
    label: String,
    left: Ann,
    right: Ann,
    aligned: Vec<Ann>,
    inner: Vec<Ann>,
}

struct Ann {
    id: String,
    start: usize,
    end: usize,
    query_length: Option<i32>,
}

impl std::fmt::Display for Ann {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.query_length {
            Some(l) => write!(f, "{},{},{},{}", self.id, self.start, self.end, l),
            None => write!(f, "{},{},{}", self.id, self.start, self.end),
        }
    }
}

impl Serialize for Ann {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl AuroraAnnotationGroup {
    fn from_joined_annotations(joins: &mut [&Annotation]) -> Self {
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
                aligned.push(Ann {
                    id: id_fn(),
                    start: a.target_start,
                    end: a.target_end,
                    query_length: None,
                });
                inner.push(Ann {
                    id: id_fn(),
                    start: a.target_end,
                    end: b.target_start,
                    query_length: Some(b.query_start as i32 - a.query_end as i32 + 1),
                });
            });

        aligned.push(Ann {
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

        let left = Ann {
            id: id_fn(),
            start: visual_start,
            end: first.target_start,
            query_length: None,
        };

        let right = Ann {
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
            label: first.query_name.clone(),
            left,
            right,
            aligned,
            inner,
        }
    }

    fn from_bed_record(bed: &BedRecord) -> Self {
        let mut id_cnt = 0usize;
        let mut id_fn = || {
            let id = format!("bed-{}-{}", bed.id, id_cnt);
            id_cnt += 1;
            id
        };

        //   5                       <-- count
        //   2939,235,11,196,11      <-- sizes
        //   -1,2940,-1,3389,-1      <-- starts
        //
        // example visual:
        //     2939, -1           234, 2940     11, -1     196, 3389     11, -1
        //   |-----------------[             ]/   --   \[             ]---------|
        //
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
                    aligned.push(Ann {
                        id: id_fn(),
                        // start is an offset from thick_start
                        start: bed.thick_start + start as usize,
                        end: bed.thick_start + start as usize + size,
                        query_length: None,
                    });
                }
                // unaligned blocks have -1 start
                else {
                    let last_aligned = aligned.last().expect("aligned group is empty");
                    inner.push(Ann {
                        id: id_fn(),
                        // this unaligned block starts 1 position
                        // after the previous aligned block
                        start: last_aligned.end + 1,
                        // it also ends 1 before the start
                        // of the next aligned block
                        end: bed.thick_start + next_start as usize - 1,
                        // the size field encodes the length of the model
                        // that is projected into the inner unaligned block
                        query_length: Some(size as i32),
                    });
                }
            });

        let left = Ann {
            id: id_fn(),
            start: bed.chrom_start,
            end: bed.thick_start,
            query_length: None,
        };

        let last_aligned = aligned.last().unwrap();
        let last_size = bed.block_sizes.last().unwrap();
        let right = Ann {
            id: id_fn(),
            start: last_aligned.end + 1,
            // the last unaligned block
            // ends at thick_end
            end: bed.thick_end,
            query_length: Some(*last_size as i32),
        };

        Self {
            id: bed.id.to_string(),
            visual_start: bed.chrom_start,
            visual_end: bed.chrom_end,
            align_start: bed.thick_start,
            align_end: bed.thick_end,
            strand: bed.strand,
            label: bed.name.clone(),
            left,
            right,
            aligned,
            inner,
        }
    }
}

pub struct BedRecord {
    chrom: String,
    chrom_start: usize,
    chrom_end: usize,
    name: String,
    score: usize,
    strand: Strand,
    thick_start: usize,
    thick_end: usize,
    reserved: usize,
    block_count: usize,
    block_sizes: Vec<usize>,
    block_starts: Vec<i32>,
    id: usize,
    description: String,
}

impl BedRecord {
    fn from_str(record_str: &str) -> Self {
        // repetetive element annotation schema:
        //   0 string  chrom;          "Reference sequence chromosome or scaffold"
        //   1 uint    chromStart;     "Start position of visualization on chromosome"
        //   2 uint    chromEnd;       "End position of visualation on chromosome"
        //   3 string  name;           "Name repeat, including the type/subtype suffix"
        //   4 uint    score;          "Divergence score"
        //   5 char[1] strand;         "+ or - for strand"
        //   6 uint    thickStart;     "Start position of aligned sequence on chromosome"
        //   7 uint    thickEnd;       "End position of aligned sequence on chromosome"
        //   8 uint    reserved;       "Reserved"
        //   9 uint    blockCount;     "Count of sequence blocks"
        //   10 lstring blockSizes;     "A comma-separated list of the block sizes(+/-)"
        //   11 lstring blockStarts;    "A comma-separated list of the block starts(+/-)"
        //   12 uint    id;             "A unique identifier for the joined annotations in this record"
        //   13 lstring description;    "A comma separated list of technical annotation descriptions"
        let tokens: Vec<&str> = record_str.split_whitespace().collect();

        Self {
            chrom: tokens[0].to_string(),
            chrom_start: tokens[1].parse::<usize>().unwrap(),
            chrom_end: tokens[2].parse::<usize>().unwrap(),
            name: tokens[3].to_string(),
            score: tokens[4].parse::<usize>().unwrap(),
            strand: Strand::from_str(tokens[5]),
            thick_start: tokens[6].parse::<usize>().unwrap(),
            thick_end: tokens[7].parse::<usize>().unwrap(),
            reserved: 0,
            block_count: tokens[9].parse::<usize>().unwrap(),
            block_sizes: tokens[10]
                .split_whitespace()
                .collect::<Vec<&str>>()
                .iter()
                .map(|d| d.parse::<usize>().unwrap())
                .collect(),
            block_starts: tokens[11]
                .split_whitespace()
                .collect::<Vec<&str>>()
                .iter()
                .map(|d| d.parse::<i32>().unwrap())
                .collect(),
            id: tokens[12].parse::<usize>().unwrap(),
            description: tokens[13].to_string(),
        }
    }

    fn from_joined_annotations(joins: &mut [&Annotation]) -> Self {
        joins.sort_by_key(|a| a.target_start);

        let first = joins.first().unwrap();
        let last = joins.last().unwrap();

        // TODO: probably need to map from target_name to "chr<num>"
        let chrom = first.target_name.clone();
        let chrom_start = first.target_start - first.query_start;
        // TODO: going to need extra information for this
        let chrom_end = last.target_end;
        let name = first.query_name.clone();
        let score = 1.0;
        let strand = first.strand;
        let thick_start = first.target_start;
        let thick_end = last.target_end;

        // we have:
        //   a block for every join
        //   a block in between each join
        //   two flanking blocks (unaligned projections to the left & right)
        //
        //   (num joins) + (num_joins - 1) + (2)
        //   (num joins * 2) + 1
        let block_count = joins.len() * 2 + 1;

        //   5                       <-- count
        //   2939,235,11,196,11      <-- sizes
        //   -1,2940,-1,3389,-1      <-- starts
        //
        // example visual:
        //     2939, -1           234, 2940     11, -1     196, 3389     11, -1
        //   |-----------------[             ]/   --   \[             ]---------|
        //

        // this creates the left flanking unaligned block
        let mut block_sizes = vec![first.target_start - chrom_start + 1];
        let mut block_starts: Vec<i32> = vec![-1];

        (0..joins.len() - 1)
            .map(|idx| (joins[idx], joins[idx + 1]))
            .for_each(|(a, b)| {
                // create the aligned block
                // *note: the target coordinates of the
                //        annotation are used here
                block_sizes.push(a.target_end - a.target_start + 1);
                block_starts.push((a.target_start - chrom_start) as i32);

                // create the unaligned block
                // *note: the query coordinates are used here to
                //        indicate the length of the query that
                //        is projected into the unaligned block
                block_sizes.push(b.query_start - a.query_end + 1);
                block_starts.push(-1);
            });

        // create the last aligned block
        block_sizes.push(last.target_end - last.target_start + 1);
        block_starts.push((last.target_start - chrom_start) as i32);

        // creat the right flanking unaligned block
        block_sizes.push(0);
        block_starts.push(-1);

        Self {
            chrom,
            chrom_start,
            chrom_end,
            name,
            score: todo!(),
            strand,
            thick_start,
            thick_end,
            block_count,
            block_sizes,
            block_starts,
            id: todo!(),
            description: todo!(),
            reserved: todo!(),
        }
    }
}

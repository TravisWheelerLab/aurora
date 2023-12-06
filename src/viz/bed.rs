use crate::{alignment::Strand, annotation::Annotation};

use super::BlockGroup;

pub struct BedRecord {
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    pub name: String,
    #[allow(dead_code)]
    pub score: usize,
    pub strand: Strand,
    pub thick_start: usize,
    pub thick_end: usize,
    #[allow(dead_code)]
    pub reserved: usize,
    pub block_count: usize,
    pub block_sizes: Vec<i32>,
    pub block_starts: Vec<i32>,
    pub id: usize,
    #[allow(dead_code)]
    pub description: String,
}

impl std::fmt::Display for BedRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {} {} {} {} {} {} {} {}",
            self.chrom_start,
            self.chrom_end,
            self.strand,
            self.thick_start,
            self.thick_end,
            self.block_count,
            self.block_sizes
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<String>>()
                .join(","),
            self.block_starts
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<String>>()
                .join(","),
            self.name,
        )
    }
}

impl PartialEq for BedRecord {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom
            && self.chrom_start == other.chrom_start
            && self.chrom_end == other.chrom_end
            && self.strand == other.strand
            && self.thick_start == other.thick_start
            && self.thick_end == other.thick_end
            && self.block_count == other.block_count
            && self.block_sizes == other.block_sizes
            && self.block_starts == other.block_starts
    }
}

impl BedRecord {
    pub fn from_tokens(tokens: &[&str]) -> Self {
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
                .split(',')
                .collect::<Vec<&str>>()
                .iter()
                .map(|d| d.parse::<i32>().unwrap())
                .collect(),
            block_starts: tokens[11]
                .split(',')
                .collect::<Vec<&str>>()
                .iter()
                .map(|d| d.parse::<i32>().unwrap())
                .collect(),
            id: tokens[12].parse::<usize>().unwrap(),
            description: tokens[13].to_string(),
        }
    }

    pub fn from_str(record_str: &str) -> Self {
        let tokens: Vec<&str> = record_str.split_whitespace().collect();
        Self::from_tokens(&tokens)
    }

    pub fn from_joined_annotations(joins: &mut [&Annotation]) -> Self {
        joins.sort_by_key(|a| a.target_start);

        let first = joins.first().unwrap();
        let last = joins.last().unwrap();

        // TODO: probably need to map from target_name to "chr<num>"
        let chrom = first.target_name.clone();
        let chrom_start = first.target_start - first.query_start;
        // TODO: going to need extra information for this
        let chrom_end = last.target_end;
        let name = first.query_name.clone();
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

        // this creates the left flanking unaligned block
        let mut block_sizes: Vec<i32> = vec![(first.target_start - chrom_start + 1) as i32];
        let mut block_starts: Vec<i32> = vec![-1];

        (0..joins.len() - 1)
            .map(|idx| (joins[idx], joins[idx + 1]))
            .for_each(|(a, b)| {
                // create the aligned block
                // *note: the target coordinates of the
                //        annotation are used here
                block_sizes.push((a.target_end - a.target_start + 1) as i32);
                block_starts.push((a.target_start - chrom_start) as i32);

                // create the unaligned block
                // *note: the query coordinates are used here to
                //        indicate the length of the query that
                //        is projected into the unaligned block
                block_sizes.push(b.query_start as i32 - a.query_end as i32 + 1);
                block_starts.push(-1);
            });

        // create the last aligned block
        block_sizes.push((last.target_end - last.target_start + 1) as i32);
        block_starts.push((last.target_start - chrom_start) as i32);

        // creat the right flanking unaligned block
        block_sizes.push(0);
        block_starts.push(-1);

        Self {
            chrom,
            chrom_start,
            chrom_end,
            name,
            score: 0,
            strand,
            thick_start,
            thick_end,
            block_count,
            block_sizes,
            block_starts,
            id: first.join_id,
            description: "".to_string(),
            reserved: 0,
        }
    }

    pub fn from_block_group(group: &BlockGroup) -> Self {
        let mut block_sizes: Vec<i32> = vec![];
        let mut block_starts: Vec<i32> = vec![];

        // +1 for subtracting across the interval
        block_sizes.push((group.left.end - group.left.start + 1) as i32);
        block_starts.push(-1);

        (0..group.inner.len())
            .map(|idx| (&group.aligned[idx], &group.inner[idx]))
            .for_each(|(aligned, inner)| {
                // aligned
                block_sizes.push((aligned.end - aligned.start) as i32);
                block_starts.push((aligned.start - group.visual_start) as i32);

                // inner
                block_sizes.push(inner.query_length.expect("inner block has no query_length"));
                block_starts.push(-1);
            });

        let last_aligned = group.aligned.last().unwrap();

        // +1 for subtracting across the interval
        block_sizes.push((last_aligned.end - last_aligned.start + 1) as i32);
        block_starts.push((last_aligned.start - group.visual_start - 1) as i32);

        // +1 for subtracting across the interval
        block_sizes.push((group.right.end + 1 - group.right.start) as i32);
        block_starts.push(-1);

        Self {
            chrom: group.target.clone(),
            chrom_start: group.visual_start,
            chrom_end: group.visual_end,
            name: group.query.clone(),
            score: 0,
            strand: group.strand,
            thick_start: group.align_start,
            thick_end: group.align_end,
            block_count: 2 + group.aligned.len() + group.inner.len(),
            block_sizes,
            block_starts,
            id: 0,
            description: "".to_string(),
            reserved: 0,
        }
    }
}

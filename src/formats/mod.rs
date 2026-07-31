mod caf;
mod ultra;

use anyhow::Result;
use std::io::{BufRead, Seek};

use crate::{alignment::AlignmentData, substitution_matrix::SubstitutionMatrix, util::VecMap};

pub enum FormatCheck {
    Valid,
    Invalid(String),
}

trait SeekableReader: BufRead + Seek {}
impl<T: BufRead + Seek> SeekableReader for T {}

pub trait AlignmentFormat {
    fn format_check(
        primary_reader: impl SeekableReader,
        secondary_reader: Option<impl SeekableReader>,
    ) -> Result<FormatCheck>;

    fn read(
        primary_reader: impl SeekableReader,
        secondary_reader: Option<impl SeekableReader>,
        substitution_matrices: VecMap<SubstitutionMatrix>,
    ) -> Result<AlignmentData>;

    fn name() -> &'static str;
}

// let substitution_matrices = VecMap::from(SubstitutionMatrix::parse(matrices)?);

fn assign_alignment_ids(alignment_data: &mut AlignmentData) {
    // Sort all alignment entries...
    let mut ali_id = 0;

    alignment_data.target_groups.iter_mut().for_each(|g| {
        g.alignments
            .sort_by(|a, b| a.target_start.cmp(&b.target_start));

        g.alignments.iter_mut().for_each(|v| {
            v.id = ali_id;
            ali_id += 1;
        });
    });
}

mod caf;
mod ultra;

use crate::{alignment::AlignmentData, substitution_matrix::SubstitutionMatrix, util::VecMap};
use anyhow::{self, Context};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

pub enum FormatCheck {
    Valid,
    Invalid(String),
}

pub trait SeekableReader: BufRead + Seek {}
impl<T: BufRead + Seek> SeekableReader for T {}

pub trait AlignmentFormat {
    fn format_check(
        primary_reader: &mut impl SeekableReader,
        secondary_reader: Option<&mut impl SeekableReader>,
    ) -> anyhow::Result<FormatCheck>;

    fn read(
        primary_reader: &mut impl SeekableReader,
        secondary_reader: Option<&mut impl SeekableReader>,
        substitution_matrices: VecMap<SubstitutionMatrix>,
    ) -> anyhow::Result<AlignmentData>;

    fn name() -> &'static str;
}

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

macro_rules! _try_formats_helper {
    ($vec:ident, $reader1:ident, $reader2:ident, $sub_matrix:ident, [$format:ty, $($formats:ty),+]) => {
        _try_formats_helper!($vec, $reader1, $reader2, $sub_matrix, [$format]);
        _try_formats_helper!($vec, $reader1, $reader2, $sub_matrix, [$($formats),+]);
    };
    ($vec:ident, $reader1:ident, $reader2:ident, $sub_matrix:ident, [$format:ty]) => {
        let format_indicator = <$format>::format_check(&mut $reader1, $reader2.as_mut())?;

        $reader1.seek(SeekFrom::Start(0))?;
        if let Some(v) = $reader2.as_mut() {
            v.seek(SeekFrom::Start(0))?;
        }

        match format_indicator {
            FormatCheck::Valid => {
                return Ok(<$format>::read(&mut $reader1, $reader2.as_mut(), $sub_matrix)?);
            },
            FormatCheck::Invalid(err_msg) => {
                $vec.push(format!("{} not recognized due to: '{}'", <$format>::name(), err_msg));
            }
        }
    };
}

macro_rules! try_formats {
    ($reader1:ident, $reader2:ident, $sub_matrix:ident, [$($formats:ty),+]) => {
    || -> anyhow::Result<AlignmentData> {
        let mut error_vec = Vec::new();

        _try_formats_helper!(error_vec, $reader1, $reader2, $sub_matrix, [$($formats),+]);

        Err(anyhow::anyhow!(format!("Unrecognized format. Details:\n{}", error_vec.join("\n"))))
    }()
    };
}

pub fn load_alignments(
    primary_file: &impl AsRef<Path>,
    secondary_file: Option<&impl AsRef<Path>>,
    matrices: &impl AsRef<Path>,
    ultra_file: Option<&impl AsRef<Path>>,
) -> anyhow::Result<AlignmentData> {
    let mut primary_reader = BufReader::new(File::open(primary_file.as_ref()).context(format!(
        "failed to open alignments file: '{}'",
        primary_file.as_ref().to_str().unwrap_or("?")
    ))?);
    let mut secondary_reader = if let Some(v) = secondary_file {
        Some(BufReader::new(File::open(v.as_ref()).context(format!(
            "failed to open supplementary alignments file: '{}'",
            v.as_ref().to_str().unwrap_or("?")
        ))?))
    } else {
        None
    };

    let substitution_matrices = VecMap::from(SubstitutionMatrix::parse(
        File::open(matrices.as_ref()).context(format!(
            "failed to load substitution matricies file: '{}'",
            matrices.as_ref().to_str().unwrap_or("?")
        ))?,
    )?);

    let mut alignment_data = try_formats!(
        primary_reader,
        secondary_reader,
        substitution_matrices,
        [caf::CAFFormat]
    )?;

    if let Some(ultra_file) = ultra_file {
        ultra::load_ultra_file(
            &mut alignment_data,
            File::open(ultra_file).context(format!(
                "failed to open ultra file: '{}'",
                ultra_file.as_ref().to_str().unwrap_or("?")
            ))?,
        )?;
    }
    assign_alignment_ids(&mut alignment_data);

    Ok(alignment_data)
}

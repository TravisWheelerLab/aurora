mod bpaf;
mod caf;
mod fasta;
mod ultra;

use crate::alignment::Strand;
use crate::{alignment::AlignmentData, substitution_matrix::SubstitutionMatrix, util::VecMap};
use anyhow::{self, Context};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

pub enum FormatCheck {
    Valid,
    Invalid(String),
}

pub trait SeekableReader: BufRead + Read + Seek {}
impl<T: BufRead + Seek + Read> SeekableReader for T {}

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

fn normalize_alignment_data(alignment_data: &mut AlignmentData) {
    // Sort all entries...
    let mut ali_id = 0;
    let mut tr_id = 1;

    alignment_data.target_groups.iter_mut().for_each(|g| {
        for al in g.alignments.iter() {
            if matches!(al.strand, Strand::Reverse) {
                eprintln!("{}", al);
            }
        }

        g.alignments
            .sort_by(|a, b| a.target_start.cmp(&b.target_start));

        g.alignments.iter_mut().for_each(|v| {
            v.id = ali_id;
            ali_id += 1;
        });

        g.tandem_repeats
            .sort_by(|a, b| a.target_start.cmp(&b.target_start));

        g.tandem_repeats.iter_mut().for_each(|v| {
            v.id = tr_id;
            tr_id += 1;
        });

        g.target_start = g
            .alignments
            .iter()
            .map(|v| v.target_start)
            .chain(g.tandem_repeats.iter().map(|v| v.target_start))
            .min()
            .expect("Empty target group!");

        g.target_end = g
            .alignments
            .iter()
            .map(|v| v.target_end)
            .chain(g.tandem_repeats.iter().map(|v| v.target_end))
            .max()
            .expect("Empty target group!")
    });

    panic!("Uh oh...")
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

fn validate_alignment_data(alignment_data: &AlignmentData) -> bool {
    for t_grp in alignment_data.target_groups.iter() {
        if !t_grp
            .alignments
            .is_sorted_by(|a, b| a.target_start <= b.target_start)
        {
            eprintln!("Alignments are not sorted!");
            return false;
        }

        if !t_grp
            .tandem_repeats
            .is_sorted_by(|a, b| a.target_start <= b.target_start)
        {
            eprintln!("Tandem repeats are not sorted!");
            return false;
        }

        for ali in t_grp.alignments.iter() {
            if ali.target_start < t_grp.target_start || ali.target_end > t_grp.target_end {
                eprintln!(
                    "Alignment not within target group bounds! Alignment: {}, Target Range: ({}, {})",
                    ali, t_grp.target_start, t_grp.target_end
                );
                return false;
            }
        }

        for repeat in t_grp.tandem_repeats.iter() {
            if repeat.target_start < t_grp.target_start || repeat.target_end > t_grp.target_end {
                eprintln!(
                    "Repeat not within target group bounds! Repeat: {:?}, Target Range: ({}, {})",
                    repeat, t_grp.target_start, t_grp.target_end
                );
                return false;
            }
        }
    }

    true
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

    let substitution_matrices: VecMap<SubstitutionMatrix> =
        SubstitutionMatrix::parse(File::open(matrices.as_ref()).context(format!(
            "failed to load substitution matricies file: '{}'",
            matrices.as_ref().to_str().unwrap_or("?")
        ))?)?
        .into_iter()
        .collect();

    let mut alignment_data = try_formats!(
        primary_reader,
        secondary_reader,
        substitution_matrices,
        [bpaf::BPAFFormat, caf::CAFFormat]
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
    normalize_alignment_data(&mut alignment_data);

    debug_assert!(validate_alignment_data(&alignment_data));
    // Check for valid skip state...
    assert_eq!(alignment_data.query_name_map.get(0), "skip");

    Ok(alignment_data)
}

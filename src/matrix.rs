use std::collections::HashMap;

use crate::{
    alignment::Strand,
    alphabet::{GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL},
    chunks::ProximityGroup,
    collapse::AssemblyGroup,
};

/// This is an auxiliary data structure that
/// describes the layout of a column major
/// (in memory) sparse matrix.
#[derive(Clone)]
pub struct MatrixDef {
    /// The logical number of rows in the matrix (also the number of alignments)
    pub num_rows: usize,
    /// The number of columns in the matrix (also the length of the target sequence)
    pub num_cols: usize,
    /// The logical start position of the target sequence
    pub target_start: usize,
    /// The total number of active cells in the matrix
    pub num_cells: usize,

    // -- sparse data --
    /// For each column, a list of the rows that are active in that column
    pub active_rows_by_col: Vec<Vec<usize>>,
    /// For each column, a list of the query model consensus positions
    pub consensus_positions_by_col: Vec<Vec<usize>>,
    /// For each column, a list of the indices that map to the alignment
    /// that corresponds to each active row
    pub ali_ids_by_col: Vec<Vec<usize>>,

    // -- row data --
    /// For each row, the start and end column indices
    pub col_range_by_logical_row: Vec<(usize, usize)>,
    /// For each row, the numeric ID of the query that corresponds to that row
    pub query_id_by_logical_row: Vec<usize>,
    /// For each row, the strand of the alignment that corresponds to that row
    pub strand_by_logical_row: Vec<Strand>,
}

impl MatrixDef {
    pub fn from_assembly_group(group: &AssemblyGroup) -> Self {
        let num_rows = group.assemblies.len() + 1;
        let num_cols = group.target_end - group.target_start + 1;

        let mut num_cells = num_cols;

        let mut active_rows_by_col = vec![vec![0usize]; num_cols];
        let mut consensus_positions_by_col = vec![vec![0usize]; num_cols];
        let mut ali_ids_by_col = vec![vec![0usize]; num_cols];

        let mut col_range_by_logical_row = vec![(0usize, num_cols); num_rows];
        let mut query_id_by_logical_row = vec![0usize; num_rows];
        let mut strand_by_logical_row = vec![Strand::default(); num_rows];

        group
            .assemblies
            .iter()
            .enumerate()
            .for_each(|(assembly_idx, assembly)| {
                assembly
                    .alignments
                    .iter()
                    .zip(assembly.alignments.iter().skip(1))
                    .for_each(|(a, b)| {
                        debug_assert!(a.target_start <= b.target_start);
                    });

                let row_idx = assembly_idx + 1;

                query_id_by_logical_row[row_idx] = assembly.query_id;
                strand_by_logical_row[row_idx] = assembly.strand;

                // find the start/end of the assembly
                // relative to the start of the group
                //
                // we have something like:
                //
                //     (         [    -------------     ]           )
                //     ^         ^         ^            ^           ^
                //  chrom      group    assembly      group       chrom
                //  start      start                   end         end
                //
                let assembly_col_start_in_matrix = assembly.target_start - group.target_start;
                let assembly_col_end_in_matrix = assembly.target_end - group.target_start;
                let assembly_len = assembly_col_end_in_matrix - assembly_col_start_in_matrix + 1;

                num_cells += assembly_col_end_in_matrix - assembly_col_start_in_matrix + 1;

                // set the entire logical row defined
                // by the assembly to be active
                //
                // we have something like:
                //
                //  --------------************-------------***********---------
                //        ^            ^           ^           ^          ^
                //    alignment      skip      alignment      skip    alignment
                //
                col_range_by_logical_row[row_idx] =
                    (assembly_col_start_in_matrix, assembly_col_end_in_matrix);
                (assembly_col_start_in_matrix..=assembly_col_end_in_matrix).for_each(|col_idx| {
                    active_rows_by_col[col_idx].push(row_idx);
                });

                // for each alignment, figure out the consensus positions that
                // map to each target position/column index in the group
                let mut assembly_consensus_map: HashMap<usize, Vec<usize>> = HashMap::new();
                assembly.alignments.iter().for_each(|ali| {
                    let mut consensus_positions = vec![0usize; assembly_len];

                    let mut consensus_position = ali.query_start as i32;

                    // the column index relative to the assembly
                    let mut assembly_col_idx = ali.target_start - assembly.target_start;

                    let addend = match assembly.strand {
                        Strand::Forward => 1,
                        Strand::Reverse => -1,
                        Strand::Unset => {
                            panic!("Alignment.strand is unset in call to MatrixDef::new()")
                        }
                    };
                    for (&target_char, &query_char) in
                        ali.target_seq.iter().zip(ali.query_seq.iter())
                    {
                        consensus_positions[assembly_col_idx] = consensus_position as usize;
                        match target_char {
                            GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                                // if the target character is a gap,
                                // we advance the consensus position
                                // and leave the column index alone
                                consensus_position += addend;
                            }
                            _ => {
                                // if the target character is not a gap,
                                // we advance the column index regardless
                                // of what the query character is
                                assembly_col_idx += 1;
                                match query_char {
                                    GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                                        // if the query character is a gap, we use
                                        // the previous query consensus position
                                    }
                                    _ => {
                                        // if the query character isn't a gap,
                                        // we advance the consensus position
                                        consensus_position += addend;
                                    }
                                }
                            }
                        }
                    }
                    debug_assert_eq!(
                        consensus_positions[ali.target_start - assembly.target_start],
                        ali.query_start
                    );

                    debug_assert_eq!(
                        consensus_positions[ali.target_end - assembly.target_start],
                        ali.query_end
                    );

                    assembly_consensus_map.insert(ali.id, consensus_positions);
                });

                let mut last_matrix_col_idx = assembly_col_start_in_matrix;
                let mut last_consensus_position = 0usize;

                assembly.alignment_ranges.iter().for_each(|range| {
                    let ali = assembly
                        .alignments
                        .iter()
                        .find(|a| a.id == range.ali_id)
                        .expect("failed to match alignment with range");

                    let consensus_positions = assembly_consensus_map
                        .get(&ali.id)
                        .expect("failed to find consensus positions");

                    let new_matrix_col_idx =
                        range.assembly_col_start + assembly_col_start_in_matrix;

                    // place the last consensus position of the previous alignment
                    // at every empty position leading up to this alignment
                    let spaces = new_matrix_col_idx - last_matrix_col_idx;
                    ((last_matrix_col_idx + 1)..(last_matrix_col_idx + spaces)).for_each(
                        |matrix_col_idx| {
                            consensus_positions_by_col[matrix_col_idx]
                                .push(last_consensus_position);
                            ali_ids_by_col[matrix_col_idx].push(0);
                        },
                    );

                    (range.assembly_col_start..=range.assembly_col_end)
                        .map(|assembly_col_idx| {
                            (
                                assembly_col_idx,
                                assembly_col_idx + assembly_col_start_in_matrix,
                            )
                        })
                        .for_each(|(assembly_col_idx, matrix_col_idx)| {
                            let consensus_position = consensus_positions[assembly_col_idx];
                            consensus_positions_by_col[matrix_col_idx].push(consensus_position);
                            ali_ids_by_col[matrix_col_idx].push(ali.id);
                        });

                    last_matrix_col_idx = range.assembly_col_end + assembly_col_start_in_matrix;
                    last_consensus_position = consensus_positions[range.assembly_col_end];
                });
            });

        Self {
            num_rows,
            num_cols,
            target_start: group.target_start,
            num_cells,
            active_rows_by_col,
            consensus_positions_by_col,
            ali_ids_by_col,
            col_range_by_logical_row,
            query_id_by_logical_row,
            strand_by_logical_row,
        }
    }

    pub fn new(group: &ProximityGroup) -> Self {
        let target_start = group.target_start;
        let target_end = group.target_end;

        // +1 for the skip state
        let num_rows = group.alignments.len() + 1;
        // +1 for subtracting across the interval
        let num_cols = target_end - target_start + 1;

        // add cells for the skip state, which has a cell in every column
        let mut num_cells = num_cols;

        let mut active_rows_by_col = vec![vec![0usize]; num_cols];
        let mut consensus_positions_by_col = vec![vec![0usize]; num_cols];
        let mut ali_ids_by_col = vec![vec![0usize]; num_cols];

        // NOTE: we are setting the first row (skip state) to span the entire matrix
        let mut col_range_by_logical_row = vec![(0usize, num_cols); num_rows];
        let mut query_id_by_logical_row = vec![0usize; num_rows];
        let mut strand_by_logical_row = vec![Strand::default(); num_rows];
        // we'll say the skip state is the forward strand for convenience
        strand_by_logical_row[0] = Strand::Forward;

        group
            .alignments
            .iter()
            .enumerate()
            .for_each(|(ali_idx, ali)| {
                // the index of the row that the alignment maps to
                // it's +1 because of the skip state
                let row_idx = ali_idx + 1;

                query_id_by_logical_row[row_idx] = ali.query_id;
                strand_by_logical_row[row_idx] = ali.strand;

                // the start/end of the alignment in terms of the columns of the matrix
                let col_start = ali.target_start - target_start;
                let col_end = ali.target_end - target_start;
                num_cells += col_end - col_start + 1;

                col_range_by_logical_row[row_idx] = (col_start, col_end);
                (col_start..=col_end).for_each(|col_idx| active_rows_by_col[col_idx].push(row_idx));

                // ali_col_idx is the column index that corresponds to the
                // current alignment position; we need to manually keep
                // track of this because of gaps in the target sequence
                let mut ali_col_idx = col_start;

                // ali_consensus_position is the consensus position that
                // corresponds to the current alignment position we are
                // processing; we need to manually keep track of this
                // because of gaps in the query sequence
                let mut ali_consensus_position = ali.query_start;

                match ali.strand {
                    Strand::Forward => {
                        debug_assert!(ali.target_start < ali.target_end);
                        debug_assert!(ali.query_start < ali.query_end);
                    }
                    Strand::Reverse => {
                        debug_assert!(ali.target_start < ali.target_end);
                        debug_assert!(ali.query_start > ali.query_end);
                    }
                    Strand::Unset => {
                        panic!("Alignment.strand is unset in call to MatrixDef::new()")
                    }
                }

                for (&target_char, &query_char) in ali.target_seq.iter().zip(ali.query_seq.iter()) {
                    match target_char {
                        GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                            // if the target character is a gap, we advance the
                            // consensus position and leave the column index alone
                            match ali.strand {
                                Strand::Forward => ali_consensus_position += 1,
                                Strand::Reverse => ali_consensus_position -= 1,
                                Strand::Unset => {
                                    panic!("Alignment.strand is unset in call to MatrixDef::new()")
                                }
                            }
                        }
                        _ => {
                            match query_char {
                                GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                                    // if the query character is a gap, we use
                                    // the previous query consensus position
                                }
                                _ => {
                                    // if the query character isn't a gap,
                                    // we advance the consensus position
                                    match ali.strand {
                                        Strand::Forward => ali_consensus_position += 1,
                                        Strand::Reverse => ali_consensus_position -= 1,
                                        Strand::Unset => {
                                            panic!(
                                            "Alignment.strand is unset in call to MatrixDef::new()"
                                        )
                                        }
                                    }
                                }
                            }

                            debug_assert!(ali_col_idx >= col_start && ali_col_idx <= col_end);

                            consensus_positions_by_col[ali_col_idx].push(ali_consensus_position);
                            ali_ids_by_col[ali_col_idx].push(ali.id);
                            ali_col_idx += 1;
                        }
                    }
                }
            });

        Self {
            num_rows,
            num_cols,
            target_start,
            num_cells,
            active_rows_by_col,
            consensus_positions_by_col,
            ali_ids_by_col,
            col_range_by_logical_row,
            query_id_by_logical_row,
            strand_by_logical_row,
        }
    }
}

pub struct Matrix<'a, T>
where
    T: Clone + Copy + Default + std::fmt::Display,
{
    pub def: &'a MatrixDef,
    pub data: Vec<Vec<T>>,
}

impl<'a, T> Matrix<'a, T>
where
    T: Clone + Copy + Default + std::fmt::Display,
{
    /// Create a new Matrix from a MatrixDef.
    pub fn new(def: &'a MatrixDef) -> Self {
        let mut data = vec![];

        def.active_rows_by_col
            .iter()
            .for_each(|v| data.push(vec![T::default(); v.len()]));

        Self { def, data }
    }

    /// Copy the data from other into self.
    /// This will panic if the MatrixDefs of self and other are incompatible.
    pub fn copy_fill(&mut self, other: &Self) {
        let mut value_map: HashMap<usize, T> = HashMap::new();
        (0..self.num_cols()).for_each(|col| {
            other
                .col_slice(col)
                .iter()
                .enumerate()
                .for_each(|(sparse_row, value)| {
                    let ali_id = other.ali_id_sparse(sparse_row, col);
                    value_map.insert(ali_id, *value);
                });
            (0..self.col_length(col)).for_each(|sparse_row| {
                let ali_id = self.ali_id_sparse(sparse_row, col);

                let value = value_map.get(&ali_id).expect("failed to find target value");
                self.set_sparse(sparse_row, col, *value)
            })
        });
    }

    // -----------
    // data access
    // -----------
    pub fn get(&self, row: usize, col: usize) -> T {
        let sparse_row_index = self.logical_to_sparse_row_idx(row, col);
        self.data[col][sparse_row_index]
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> &mut T {
        let sparse_row_index = self.logical_to_sparse_row_idx(row, col);
        &mut self.data[col][sparse_row_index]
    }

    pub fn set(&mut self, row: usize, col: usize, value: T) {
        let sparse_row_index = self.logical_to_sparse_row_idx(row, col);
        self.data[col][sparse_row_index] = value;
    }

    pub fn get_skip(&self, col: usize) -> T {
        self.data[col][0]
    }

    pub fn set_skip(&mut self, col: usize, value: T) {
        self.data[col][0] = value;
    }

    pub fn set_sparse(&mut self, sparse_row: usize, col: usize, value: T) {
        self.data[col][sparse_row] = value;
    }

    pub fn get_sparse(&self, sparse_row: usize, col: usize) -> T {
        self.data[col][sparse_row]
    }

    pub fn num_cols(&self) -> usize {
        self.def.num_cols
    }

    pub fn num_rows(&self) -> usize {
        self.def.num_rows
    }

    pub fn target_start(&self) -> usize {
        self.def.target_start
    }

    pub fn col_slice(&self, col: usize) -> &[T] {
        &self.data[col]
    }

    pub fn col_slice_mut(&mut self, col: usize) -> &mut [T] {
        &mut self.data[col]
    }

    pub fn col_length(&self, col: usize) -> usize {
        self.data[col].len()
    }

    pub fn col_range_by_row(&self, row: usize) -> (usize, usize) {
        self.def.col_range_by_logical_row[row]
    }

    pub fn initial_active_cols(&self) -> Vec<usize> {
        (0..self.num_cols())
            .filter(|&col_idx| self.def.ali_ids_by_col[col_idx].iter().any(|&id| id != 0))
            .collect()
    }

    pub fn query_id_of_row(&self, row: usize) -> usize {
        self.def.query_id_by_logical_row[row]
    }

    pub fn query_id_of_cell_sparse(&self, sparse_row: usize, col: usize) -> usize {
        self.def.query_id_by_logical_row[self.def.active_rows_by_col[col][sparse_row]]
    }

    pub fn ali_id(&self, row: usize, col: usize) -> usize {
        let sparse_row_index = self.logical_to_sparse_row_idx(row, col);
        self.def.ali_ids_by_col[col][sparse_row_index]
    }

    pub fn ali_id_sparse(&self, sparse_row: usize, col: usize) -> usize {
        self.def.ali_ids_by_col[col][sparse_row]
    }

    pub fn strand_of_row(&self, row: usize) -> Strand {
        self.def.strand_by_logical_row[row]
    }

    pub fn strand_of_cell_sparse(&self, sparse_row: usize, col: usize) -> Strand {
        self.def.strand_by_logical_row[self.sparse_to_logical_row_idx(sparse_row, col)]
    }

    pub fn consensus_position_sparse(&self, sparse_row: usize, col: usize) -> usize {
        self.def.consensus_positions_by_col[col][sparse_row]
    }

    pub fn logical_to_sparse_row_idx(&self, row: usize, col: usize) -> usize {
        self.def.active_rows_by_col[col]
            .binary_search(&row)
            .unwrap_or_else(|_| panic!("invalid logical row: {row} in column: {col}"))
    }

    pub fn sparse_to_logical_row_idx(&self, sparse_row: usize, col: usize) -> usize {
        self.def.active_rows_by_col[col][sparse_row]
    }

    pub fn contains_cell(&self, row: usize, col: usize) -> bool {
        for &row_idx in self.def.active_rows_by_col[col].iter() {
            if row_idx == row {
                return true;
            };
        }
        false
    }

    #[allow(dead_code)]
    pub fn print(&self) {
        (0..self.num_rows()).for_each(|row_idx| {
            (0..self.num_cols()).for_each(|col_idx| {
                //
                if self.contains_cell(row_idx, col_idx) {
                    print!("{:8.3} ", self.get(row_idx, col_idx));
                } else {
                    print!("{:>8.3} ", "x");
                }
            });
            println!();
        });
    }
}

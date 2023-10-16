use crate::{
    alignment::Strand,
    alphabet::{GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL},
    chunks::ProximityGroup,
};

/// This is an auxiliary data structure that
/// describes the layout of a column major
/// (in memory) sparse matrix.
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

    // -- row data --
    /// For each row, the start and end column indices
    pub col_range_by_row: Vec<(usize, usize)>,
    /// For each row, the numeric ID of the alignment that corresponds to that row
    pub id_by_row: Vec<usize>,
    /// For each row, the strand of the alignment that corresponds to that row
    pub strand_by_row: Vec<Strand>,
}

impl MatrixDef {
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

        // NOTE: we are setting the first row (skip state) to span the entire matrix
        let mut col_range_by_row = vec![(0, num_cols); num_rows];
        let mut id_by_row = vec![0usize; num_rows];
        let mut strand_by_row = vec![Strand::default(); num_rows];
        // we'll say the skip state is the forward strand for convenience
        strand_by_row[0] = Strand::Forward;

        group
            .alignments
            .iter()
            .enumerate()
            .for_each(|(ali_idx, ali)| {
                // the index of the row that the alignment maps to
                // it's +1 because of the skip state
                let row_idx = ali_idx + 1;

                id_by_row[row_idx] = ali.query_id;
                strand_by_row[row_idx] = ali.strand;

                // the start/end of the alignment in terms of the columns of the matrix
                let col_start = ali.target_start - target_start;
                let col_end = ali.target_end - target_start;
                num_cells += col_end - col_start + 1;

                col_range_by_row[row_idx] = (col_start, col_end);
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
            col_range_by_row,
            id_by_row,
            strand_by_row,
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
    pub fn new(def: &'a MatrixDef) -> Self {
        let mut data = vec![];

        def.active_rows_by_col
            .iter()
            .for_each(|v| data.push(vec![T::default(); v.len()]));

        Self { def, data }
    }

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
        self.def.col_range_by_row[row]
    }

    pub fn initial_active_cols(&self) -> Vec<usize> {
        (0..self.num_cols())
            .filter(|&col_idx| self.def.active_rows_by_col[col_idx].len() > 1)
            .collect()
    }

    pub fn id_of_row(&self, row: usize) -> usize {
        self.def.id_by_row[row]
    }

    pub fn id_of_cell_sparse(&self, sparse_row: usize, col: usize) -> usize {
        self.def.id_by_row[self.def.active_rows_by_col[col][sparse_row]]
    }

    pub fn strand_of_row(&self, row: usize) -> Strand {
        self.def.strand_by_row[row]
    }

    pub fn strand_of_cell_sparse(&self, sparse_row: usize, col: usize) -> Strand {
        self.def.strand_by_row[self.sparse_to_logical_row_idx(sparse_row, col)]
    }

    pub fn consensus_position_sparse(&self, sparse_row: usize, col: usize) -> usize {
        self.def.consensus_positions_by_col[col][sparse_row]
    }

    pub fn logical_to_sparse_row_idx(&self, row: usize, col: usize) -> usize {
        for (sparse_idx, &logical_idx) in self.def.active_rows_by_col[col].iter().enumerate() {
            if logical_idx == row {
                return sparse_idx;
            }
        }
        panic!("invalid logical row: {row} in column: {col}");
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

pub enum RepresentativeType {
    DEFAULT,
    SMALLEST,
    LARGEST,
}

pub struct UnionFind {
    representative_mode: RepresentativeType,
    parents: Vec<usize>,
    sizes: Vec<usize>,
    representative: Vec<usize>,
}

/// Optimal implementation of union-find with consistent representatives for each group...
/// The representative for a group is always the smallest element in the group.
impl UnionFind {
    pub fn new(size: usize, mode: RepresentativeType) -> Self {
        Self {
            representative_mode: mode,
            parents: (0..size).collect(),
            sizes: vec![0; size],
            representative: (0..size).collect(),
        }
    }

    fn _root(&self, mut node: usize) -> usize {
        while self.parents[node] != node {
            node = self.parents[node];
        }
        node
    }

    fn _collapse_path(&mut self, mut node: usize, root: usize) {
        while self.parents[node] != root {
            let new_node = self.parents[node];
            self.parents[node] = root;
            node = new_node;
        }
    }

    pub fn union(&mut self, node1: usize, node2: usize) {
        let mut root1 = self._find_root(node1);
        let mut root2 = self._find_root(node2);

        if root1 == root2 {
            return;
        }

        // If 1st root is smaller, swap so we merge smaller into larger tree...
        if self.sizes[root1] < self.sizes[root2] {
            let tmp = root1;
            root1 = root2;
            root2 = tmp;
        }

        self.parents[root2] = root1;
        self.sizes[root1] += self.sizes[root2];
        match self.representative_mode {
            RepresentativeType::DEFAULT => self.representative[root1] = root1,
            RepresentativeType::SMALLEST => {
                self.representative[root1] =
                    self.representative[root1].min(self.representative[root2])
            }
            RepresentativeType::LARGEST => {
                self.representative[root1] =
                    self.representative[root1].max(self.representative[root2])
            }
        };
    }

    fn _find_root(&mut self, node: usize) -> usize {
        let root = self._root(node);
        self._collapse_path(node, root);

        node
    }

    pub fn find(&mut self, node: usize) -> usize {
        let root = self._find_root(node);
        return self.representative[root];
    }

    pub fn find_unmut(&self, node: usize) -> usize {
        let root = self._root(node);
        return self.representative[root];
    }
}

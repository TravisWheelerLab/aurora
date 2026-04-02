#[derive(Debug)]
pub enum RepresentativeType {
    #[allow(dead_code)]
    DEFAULT,
    #[allow(dead_code)]
    SMALLEST,
    LARGEST,
}

#[derive(Debug)]
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
            sizes: vec![1; size],
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
        // Collapse node and all it's parents to point to root.
        while self.parents[node] != root {
            let new_node = self.parents[node];
            self.parents[node] = root;
            node = new_node;
        }
    }

    pub fn union(&mut self, node1: usize, node2: usize) -> usize {
        let mut root1 = self._find_root(node1);
        let mut root2 = self._find_root(node2);

        if root1 == root2 {
            return self.representative[root1];
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

        self.representative[root1]
    }

    fn _find_root(&mut self, node: usize) -> usize {
        let root = self._root(node);
        self._collapse_path(node, root);

        root
    }

    pub fn find(&mut self, node: usize) -> usize {
        let root = self._find_root(node);
        return self.representative[root];
    }

    #[allow(dead_code)]
    pub fn find_unmut(&self, node: usize) -> usize {
        let root = self._root(node);
        return self.representative[root];
    }
}

#[cfg(test)]
mod test {
    use itertools::izip;

    use super::{RepresentativeType, UnionFind};

    #[test]
    fn test_union_find() {
        let uf_links = [[(1, 5), (5, 7), (1, 2)], [(1, 5), (5, 7), (1, 2)]];
        let uf_results = [
            [0, 1, 1, 3, 4, 1, 6, 1, 8, 9],
            [0, 7, 7, 3, 4, 7, 6, 7, 8, 9],
        ];
        let modes = [RepresentativeType::SMALLEST, RepresentativeType::LARGEST];

        for (links, results, mode) in izip!(uf_links, uf_results, modes) {
            let mut uf = UnionFind::new(10, mode);

            for link in links {
                uf.union(link.0, link.1);
            }

            for (src, &dst) in results.iter().enumerate() {
                assert_eq!(uf.find(src), dst);
            }
        }
    }
}

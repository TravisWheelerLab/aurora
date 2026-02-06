use std::{cmp::Ordering, fmt::Debug};

#[derive(Debug)]
struct TreeNode<T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> {
    left: Option<T>,
    right: Option<T>,
    depth: u8,
}

#[derive(Debug)]
pub struct AVLIndexSet<T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> {
    root: T,
    nodes: Vec<TreeNode<T>>,
}

impl<T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> AVLIndexSet<T> {
    pub fn new() -> Self {
        Self {
            root: T::default(),
            nodes: Vec::new(),
        }
    }

    pub fn add(
        &mut self,
        compare: impl Fn(usize) -> Ordering,
    ) -> Result<usize, <T as TryFrom<usize>>::Error> {
        // Create a new node...
        let new_idx: T = self.nodes.len().try_into()?;

        let (new_root, inserted_index) = self._insert_into_tree(self.root, new_idx, compare);
        self.root = new_root;

        if inserted_index.into() == new_idx.into() {
            self.nodes.push(TreeNode {
                left: None,
                right: None,
                depth: 0,
            });
        }

        Result::Ok(new_idx.into())
    }

    fn _depth(&self, index: Option<T>) -> u8 {
        match index {
            Some(idx) => self.nodes[idx.into()].depth,
            None => 0,
        }
    }

    fn _update_depth_from_children(&mut self, root: T) {
        let root_big = root.into();
        self.nodes[root_big].depth = self
            ._depth(self.nodes[root_big].left)
            .max(self._depth(self.nodes[root_big].right))
            + 1
    }

    fn _rotate_left(&mut self, root: T) -> T {
        let root_big = root.into();

        match self.nodes[root_big].right {
            Some(right_idx) => {
                let right_idx_big = right_idx.into();
                self.nodes[root_big].right = self.nodes[right_idx_big].left;
                self.nodes[right_idx_big].left = Some(root);
                self._update_depth_from_children(root);
                self._update_depth_from_children(right_idx);
                right_idx
            }
            _ => {
                // Rotation can't happen as there is no right node, ignore...
                root
            }
        }
    }

    fn _rotate_right(&mut self, root: T) -> T {
        let root_big = root.into();

        match self.nodes[root_big].left {
            Some(left_idx) => {
                let left_idx_big = left_idx.into();
                self.nodes[root_big].right = self.nodes[left_idx_big].right;
                self.nodes[left_idx_big].right = Some(root);
                self._update_depth_from_children(root);
                self._update_depth_from_children(left_idx);
                left_idx
            }
            _ => {
                // Rotation can't happen as there is no right node, ignore...
                root
            }
        }
    }

    fn _insert_into_tree(
        &mut self,
        root: T,
        new_node_idx: T,
        compare: impl Fn(usize) -> Ordering,
    ) -> (T, T) {
        if self.nodes.len() == 0 {
            return (new_node_idx, new_node_idx);
        }

        let root_big: usize = root.into();
        let inserted_index;

        match compare(root.into()) {
            // Case a: if already found, just return the index...
            Ordering::Equal => {
                return (root, root);
            }
            // Smaller, insert if no right node, otherwise
            Ordering::Greater => {
                match self.nodes[root_big].left {
                    Some(idx) => {
                        let root_and_insert = self._insert_into_tree(idx, new_node_idx, compare);
                        self.nodes[root_big].left = Some(root_and_insert.0);
                        inserted_index = root_and_insert.1;
                    }
                    _ => {
                        self.nodes[root_big].left = Some(new_node_idx);
                        inserted_index = new_node_idx;
                        self.nodes[root_big].depth = self.nodes[root_big].depth.max(1);
                    }
                };
            }
            Ordering::Less => {
                match self.nodes[root_big].right {
                    Some(idx) => {
                        let root_and_insert = self._insert_into_tree(idx, new_node_idx, compare);
                        self.nodes[root_big].right = Some(root_and_insert.0);
                        inserted_index = root_and_insert.1;
                    }
                    _ => {
                        self.nodes[root_big].right = Some(new_node_idx);
                        inserted_index = new_node_idx;
                        self.nodes[root_big].depth = self.nodes[root_big].depth.max(1);
                    }
                };
            }
        }

        let left_depth = self._depth(self.nodes[root_big].left);
        let right_depth = self._depth(self.nodes[root_big].right);

        let new_root = if left_depth > right_depth + 1 {
            self._rotate_right(root)
        } else if right_depth > left_depth + 1 {
            self._rotate_left(root)
        } else {
            self._update_depth_from_children(root);
            root
        };

        return (new_root, inserted_index);
    }

    fn search(&self, compare: impl Fn(usize) -> Ordering) -> Option<usize> {
        self._search(self.root.into(), compare)
    }

    fn _search(&self, root: usize, compare: impl Fn(usize) -> Ordering) -> Option<usize> {
        if root >= self.nodes.len() {
            return None;
        }

        let node = &self.nodes[root];

        match compare(root) {
            Ordering::Equal => Some(root),
            Ordering::Greater => match node.left {
                Some(idx) => self._search(idx.into(), compare),
                _ => None,
            },
            Ordering::Less => match node.right {
                Some(idx) => self._search(idx.into(), compare),
                _ => None,
            },
        }
    }
}

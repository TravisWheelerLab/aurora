use std::{cmp::Ordering, collections::VecDeque, fmt::Debug};

use serde_json::map::Iter;

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
            Some(idx) => {
                self.nodes[idx.into()].depth
            }
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

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Get the depth of the tree, including the root node...
    pub fn depth(&self) -> usize {
        if self.root.into() >= self.nodes.len() {
            0
        } else {
            (self._depth(Some(self.root)) as usize) + 1
        }
    }

    pub fn iter(&self) -> AVLInOrderSetIterator<T> {
        let mut stack = Vec::with_capacity(self.depth() + 1);
        stack.push((self.root, 0 as u8));

        AVLInOrderSetIterator { tree: &self, stack }
    }

    pub fn bfs(&self) -> AVLBFSSetIterator<T> {
        AVLBFSSetIterator {
            tree: &self,
            queue: VecDeque::from([Some(self.root)]),
            level: 0,
            offset: 0,
        }
    }
}


pub struct AVLBFSSetIterator<'a, T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> {
    tree: &'a AVLIndexSet<T>,
    queue: VecDeque<Option<T>>,
    level: usize,
    offset: usize,
}

pub struct BFSInfo {
    pub node_index: usize,
    pub level: usize,
    pub offset: usize
}

impl<'a, T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> Iterator
    for AVLBFSSetIterator<'a, T>
{
    type Item = BFSInfo;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(child_elem) = self.queue.pop_front() {
            if let Some(element) = child_elem {
                self.queue.push_front(self.tree.nodes[element.into()].left);
                self.queue.push_front(self.tree.nodes[element.into()].right);
                let result = Some(BFSInfo{
                    node_index: element.into(),
                    level: self.level,
                    offset: self.offset
                });
                self.offset += 1;
                if self.offset >= (1 << self.level) {
                    self.level += 1;
                    self.offset = 0;
                }
                return result;
            }
            self.offset = 0
        } 

        None
    }
}

pub struct AVLInOrderSetIterator<'a, T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> {
    tree: &'a AVLIndexSet<T>,
    stack: Vec<(T, u8)>,
}

impl<'a, T: Into<usize> + TryFrom<usize> + Copy + Default + Debug> Iterator
    for AVLInOrderSetIterator<'a, T>
{
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((idx, phase)) = self.stack.pop() {
            let idx_big = idx.into();
            // Case: Empty tree, just return None...
            if idx_big >= self.tree.len() {
                print!("YEP!");
                break;
            }

            match phase {
                // Attempt to go left...
                0 => {
                    self.stack.push((idx, 1));
                    if let Some(left) = self.tree.nodes[idx_big].left {
                        self.stack.push((left, 0));
                    }
                }
                // Printout this value, add right value to the stack...
                1 => {
                    let depth = self.stack.len();
                    self.stack.push((idx, 2));
                    if let Some(right) = self.tree.nodes[idx_big].right {
                        self.stack.push((right, 0));
                    }
                    return Some((idx_big, depth));
                }
                _ => {}
            }
        }

        None
    }
}

mod tests {
    use std::fmt::Display;

    use anyhow::Ok;
    use itertools::Itertools;

    use crate::balanced_tree::AVLIndexSet;

    struct DummyTree {
        pub tree: AVLIndexSet<u16>,
        pub items: Vec<usize>,
    }

    impl DummyTree {
        fn build_from_list(list: &[usize]) -> Self {
            let mut tree = AVLIndexSet::new();
            let items = list.to_vec();

            for item in items.iter() {
                tree.add(|v| items[v].cmp(item))
                    .expect("Failed to add value!");
            }

            Self { tree, items }
        }
    }

    impl Display for DummyTree {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            writeln!(
                f,
                "{:?}",
                self.tree.iter().map(|v| (self.items[v.0], v.1)).collect_vec()
            )?;

            let max_element_size = self.items.iter().map(|v| v.to_string().len()).max().unwrap_or(1);
            let spacing = 3;
            let depth = self.tree.depth();
            if depth == 0 {
                return Result::Ok(());
            }
            let elements_at_bottom = (1 << (depth - 1));

            let total_width = max_element_size * elements_at_bottom + spacing * elements_at_bottom;

            for level in 0..depth {
                let element_count: usize = (1 << level);
            }            

            Result::Ok(())
        }
    }


    #[test]
    fn test_tree_construction() {
        let tree = DummyTree::build_from_list(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        println!("{:?}", tree.tree);
        println!("{}", tree);
        assert!(1 == 2);
    }
}

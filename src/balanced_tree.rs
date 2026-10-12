use std::{cmp::Ordering, collections::VecDeque, fmt::Debug};

#[derive(Debug)]
struct TreeNode<T: Into<usize> + TryFrom<usize> + Copy + Debug> {
    left: Option<T>,
    right: Option<T>,
    depth: u8,
}

#[derive(Debug)]
pub struct AVLIndexSet<T: Into<usize> + TryFrom<usize> + Copy + Debug> {
    root: Option<T>,
    nodes: Vec<TreeNode<T>>,
}

pub enum SetInsert {
    New(usize),
    Found(usize),
}

#[allow(dead_code)]
impl<T: Into<usize> + TryFrom<usize> + Copy + Debug> AVLIndexSet<T> {
    pub fn new() -> Self {
        Self {
            root: None,
            nodes: Vec::new(),
        }
    }

    pub fn add(
        &mut self,
        compare: impl Fn(usize) -> Ordering,
    ) -> Result<SetInsert, <T as TryFrom<usize>>::Error> {
        // Create a new node...
        let size = self.nodes.len();
        let new_idx: T = size.try_into()?;
        self.nodes.push(TreeNode {
            left: None,
            right: None,
            depth: 1,
        });

        let (new_root, inserted_index) = self._insert_into_tree(self.root, new_idx, compare);
        self.root = Some(new_root);

        // If it found a node that has the same value in the tree, delete the node we just inserted as it already exists...
        if inserted_index.into() != new_idx.into() {
            self.nodes.pop();
            Result::Ok(SetInsert::Found(inserted_index.into()))
        } else {
            Result::Ok(SetInsert::New(inserted_index.into()))
        }
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

    fn _bring_right_node_up(&mut self, root: T) -> T {
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

    fn _bring_left_node_up(&mut self, root: T) -> T {
        let root_big = root.into();

        match self.nodes[root_big].left {
            Some(left_idx) => {
                let left_idx_big = left_idx.into();
                self.nodes[root_big].left = self.nodes[left_idx_big].right;
                self.nodes[left_idx_big].right = Some(root);
                self._update_depth_from_children(root);
                self._update_depth_from_children(left_idx);
                left_idx
            }
            _ => {
                // Rotation can't happen as there is no left node, ignore...
                root
            }
        }
    }

    fn _insert_into_tree(
        &mut self,
        optional_root: Option<T>,
        new_node_idx: T,
        compare: impl Fn(usize) -> Ordering,
    ) -> (T, T) {
        if let Some(root) = optional_root {
            let root_big: usize = root.into();
            let inserted_index;

            match compare(root.into()) {
                // Case a: if already found, just return the index...
                Ordering::Equal => {
                    return (root, root);
                }
                // Smaller, insert if no right node, otherwise continue recursively...
                Ordering::Greater => {
                    let root_and_insert: (T, T) =
                        self._insert_into_tree(self.nodes[root_big].left, new_node_idx, compare);
                    self.nodes[root_big].left = Some(root_and_insert.0);
                    inserted_index = root_and_insert.1;
                }
                Ordering::Less => {
                    let root_and_insert =
                        self._insert_into_tree(self.nodes[root_big].right, new_node_idx, compare);
                    self.nodes[root_big].right = Some(root_and_insert.0);
                    inserted_index = root_and_insert.1;
                }
            }
            self._update_depth_from_children(root);

            let left_depth = self._depth(self.nodes[root_big].left);
            let right_depth = self._depth(self.nodes[root_big].right);

            let new_root = if left_depth > right_depth + 1 {
                self._bring_left_node_up(root)
            } else if right_depth > left_depth + 1 {
                self._bring_right_node_up(root)
            } else {
                root
            };

            (new_root, inserted_index)
        } else {
            (new_node_idx, new_node_idx)
        }
    }

    pub fn search(&self, compare: impl Fn(usize) -> Ordering) -> Option<usize> {
        self._search(self.root, compare)
    }

    fn _search(&self, root: Option<T>, compare: impl Fn(usize) -> Ordering) -> Option<usize> {
        if let Some(some_root) = root {
            let node = &self.nodes[some_root.into()];

            match compare(some_root.into()) {
                Ordering::Equal => Some(some_root.into()),
                Ordering::Greater => self._search(node.left, compare),
                Ordering::Less => self._search(node.right, compare),
            }
        } else {
            None
        }
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Get the depth of the tree, including the root node...
    pub fn depth(&self) -> usize {
        self._depth(self.root) as usize
    }

    /// Iterate over the elements in the tree, in sorted order. Returns the element index and the depth of the element in the tree.
    pub fn iter(&self) -> AVLInOrderSetIterator<'_, T> {
        let mut stack = Vec::with_capacity(self.depth() + 1);
        if let Some(root) = self.root {
            stack.push((root, 0_u8));
        }

        AVLInOrderSetIterator { tree: self, stack }
    }

    pub fn bfs(&self) -> AVLBFSSetIterator<'_, T> {
        let mut queue = VecDeque::new();

        if let Some(root) = self.root {
            queue.push_back(BFSInfo {
                node_index: root.into(),
                level: 0,
                offset: 0,
            });
        }

        AVLBFSSetIterator { tree: self, queue }
    }
}

pub struct AVLBFSSetIterator<'a, T: Into<usize> + TryFrom<usize> + Copy + Debug> {
    tree: &'a AVLIndexSet<T>,
    queue: VecDeque<BFSInfo>,
}

#[derive(Debug)]
pub struct BFSInfo {
    pub node_index: usize,
    pub level: usize,
    pub offset: usize,
}

impl<T: Into<usize> + TryFrom<usize> + Copy + Debug> Iterator for AVLBFSSetIterator<'_, T> {
    type Item = BFSInfo;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(element) = self.queue.pop_front() {
            // Add left and right to queue.
            if let Some(element_left) = self.tree.nodes[element.node_index].left {
                self.queue.push_back(BFSInfo {
                    node_index: element_left.into(),
                    level: element.level + 1,
                    offset: element.offset * 2,
                });
            }
            if let Some(element_right) = self.tree.nodes[element.node_index].right {
                self.queue.push_back(BFSInfo {
                    node_index: element_right.into(),
                    level: element.level + 1,
                    offset: element.offset * 2 + 1,
                });
            }
            // Return current element.
            return Some(element);
        }

        None
    }
}

pub struct AVLInOrderSetIterator<'a, T: Into<usize> + TryFrom<usize> + Copy + Debug> {
    tree: &'a AVLIndexSet<T>,
    stack: Vec<(T, u8)>,
}

impl<T: Into<usize> + TryFrom<usize> + Copy + Debug> Iterator for AVLInOrderSetIterator<'_, T> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((idx, phase)) = self.stack.pop() {
            let idx_big = idx.into();

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
                    // Can be used to pass depth into tree if desired...
                    // let depth = self.stack.len();
                    self.stack.push((idx, 2));
                    if let Some(right) = self.tree.nodes[idx_big].right {
                        self.stack.push((right, 0));
                    }
                    return Some(idx_big);
                }
                _ => {}
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::balanced_tree::AVLIndexSet;
    use std::fmt::Display;

    struct DummyTree {
        pub tree: AVLIndexSet<u16>,
        pub items: Vec<usize>,
    }

    impl DummyTree {
        fn build_from_list(list: &[usize]) -> Self {
            let mut tree = AVLIndexSet::new();
            let mut items: Vec<usize> = Vec::with_capacity(list.len());

            for item in list.iter() {
                if let super::SetInsert::New(_val) = tree
                    .add(|v| items[v].cmp(item))
                    .expect("Failed to add value!")
                {
                    items.push(*item);
                }
            }

            Self { tree, items }
        }

        // Gives a compact, unique representation for the tree we can use for testing...
        fn to_test_string(&self) -> String {
            self.tree
                .bfs()
                .map(|elem| {
                    format!(
                        "{}:{}:{}",
                        self.items[elem.node_index], elem.level, elem.offset
                    )
                })
                .join(",")
        }

        fn search(&self, val: usize) -> Option<usize> {
            self.tree.search(|v| self.items[v].cmp(&val))
        }

        fn iter(&self) -> impl Iterator<Item = usize> + use<'_> {
            self.tree.iter().map(|v| self.items[v])
        }
    }

    impl Display for DummyTree {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let max_element_size = self
                .items
                .iter()
                .map(|v| v.to_string().len())
                .max()
                .unwrap_or(1);
            let spacing = 3;
            let depth = self.tree.depth();
            let width = if depth > 0 {
                max_element_size * (1 << (depth - 1)) + spacing * ((1 << (depth - 1)) - 1)
            } else {
                0
            };

            let mut prior_level = 0;
            let mut prior_offset = 0;

            write!(f, "╔{}╗\n║", "═".repeat(width))?;
            for element in self.tree.bfs() {
                if prior_level != element.level {
                    write!(f, "{}║\n║", " ".repeat(width - prior_offset))?;
                    prior_level = element.level;
                    prior_offset = 0;
                }
                let mul: usize = 1 << ((depth - 1) - element.level);
                let children_start = mul * (max_element_size + spacing) * element.offset;
                let half_offset =
                    (mul * (max_element_size + spacing) - (max_element_size + spacing)) / 2;
                let offset = children_start + half_offset;
                let num_string = self.items[element.node_index].to_string();
                let space_count = offset - prior_offset;

                write!(f, "{}{}", " ".repeat(space_count), num_string)?;
                prior_offset = offset + num_string.len();
            }
            writeln!(
                f,
                "{}║\n╚{}╝",
                " ".repeat(width - prior_offset),
                "═".repeat(width)
            )?;

            Result::Ok(())
        }
    }

    fn get_test_trees() -> Vec<DummyTree> {
        vec![
            DummyTree::build_from_list(&[1, 10, 4, 5, 2, 3, 7, 6, 4, 9, 4, 6, 8, 3, 2, 7, 5, 11]),
            DummyTree::build_from_list(&[8, 6, 7, 5, 3, 0, 9]),
            DummyTree::build_from_list(&[]),
            DummyTree::build_from_list(&[4]),
            DummyTree::build_from_list(&(0..20).collect::<Vec<usize>>()),
        ]
    }

    #[test]
    fn test_tree_structure() {
        let correct_trees = [
            "7:0:0,4:1:0,9:1:1,2:2:0,5:2:1,8:2:2,10:2:3,1:3:0,3:3:1,6:3:3,11:3:7",
            "6:0:0,3:1:0,8:1:1,0:2:0,5:2:1,7:2:2,9:2:3",
            "",
            "4:0:0",
            "7:0:0,3:1:0,15:1:1,1:2:0,5:2:1,11:2:2,17:2:3,0:3:0,2:3:1,4:3:2,6:3:3,9:3:4,13:3:5,16:3:6,18:3:7,8:4:8,10:4:9,12:4:10,14:4:11,19:4:15"
        ];

        for (tree, expected_structure) in get_test_trees().iter().zip(correct_trees) {
            assert_eq!(tree.to_test_string(), expected_structure);
        }
    }

    #[test]
    fn test_tree_properties() {
        let correct_properties = [(4, 11), (3, 7), (0, 0), (1, 1), (5, 20)];

        for (tree, (expected_depth, expected_size)) in
            get_test_trees().iter().zip(correct_properties)
        {
            assert_eq!(tree.tree.depth(), expected_depth);
            assert_eq!(tree.tree.len(), expected_size);
        }
    }

    #[test]
    fn test_tree_inorder_traversal() {
        let correct_ordered_traversal: Vec<Vec<usize>> = vec![
            (1..=11).collect(),
            vec![0, 3, 5, 6, 7, 8, 9],
            vec![],
            vec![4],
            (0..20).collect(),
        ];

        for (tree, expected_result) in get_test_trees().iter().zip(correct_ordered_traversal) {
            assert_eq!(tree.iter().collect::<Vec<usize>>(), expected_result);
        }
    }

    #[test]
    fn test_tree_search() {
        let tree = &get_test_trees()[1];

        assert_eq!(tree.search(8), Some(0));
        assert_eq!(tree.search(7), Some(2));
        assert_eq!(tree.search(3), Some(4));
        assert_eq!(tree.search(0), Some(5));
        assert_eq!(tree.search(10), None);
        assert_eq!(tree.search(2), None);
        assert_eq!(tree.search(4), None);
    }
}

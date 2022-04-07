use super::{
    element::{Element, Materials},
    h_refinement::{HLevels, HRef, HRefError, HRefLoc},
    p_refinement::PolyOrders,
    space::{M2D, V2D},
    EXPECTED_NUM_H_REFINEMENTS,
};
use json::{array, object, JsonValue};
use smallvec::SmallVec;
use std::fmt;
use std::sync::Arc;

/// `Elem`s are the basic geometric unit in the `Mesh` in Parametric Space
///
/// `Elem`s are responsible for keeping track of:
/// * Connections to neighboring `Node`s and `Edge`s
/// * Connections to their parent `Elem` (and their own h-refinement state)
/// * Connections to their child `Elem`s (if h-refined)
/// * Polynomial expansion orders (p-refinement state)
///
/// `Elem`s also maintain a connection to their associated [`Element`] for descriptions of material parameters and mappings to Real Space
///
/// ## Layout
/// The indices of `Node`s and `Edge`s from the perspective of an `Elem` are described as follows:
///
/// ```text
///               N
///         2 --------- 3
///         |     1     |
///         |           |
///      W  |2         3|  E
///         |           |
///         |     0     |
///         0 --------- 1
///               S
/// ```
/// The cardinal directions are also shown. These are used in the h-refinement module to describe different types of h-refinement
///
///
/// ## h-Refinement
///
/// Three variants of h-refinements are supported. The relative indices of the child `Elem`s and their associated `Node`s and `Edge`s are shown below for each type:
///
/// 1. **T-Type**:
/// ```text
///                 1
///     2 --------------------- 3
///     |2    1    3|2    1    3|
///     |           |           |
///     |2    2    3|2    3    3|
///     |           |           |
///     |0    0    1|0    0    1|
///   2 |----------- -----------| 3
///     |2    1    3|2    1    3|
///     |           |           |
///     |2    0    3|2    1    3|
///     |           |           |
///     |0    0    1|0    0    1|
///     0 --------------------- 1
///                 0
/// ```
///
/// 2. **U-Type**:
/// ```text
///                 1
///     2 --------------------- 3
///     |2    1    3|2    1   3 |
///     |           |           |
///     |           |           |
///     |           |           |
///     |           |           |
///   2 |2    0    3|2    1    3| 3
///     |           |           |
///     |           |           |
///     |           |           |
///     |           |           |
///     |0    0    1|0    0    1|
///     0 --------------------- 1
///                 0
/// ```
///
/// 3. **V-Type**:
/// ```text
///                 1
///     2 --------------------- 3
///     |2          1          3|
///     |                       |
///     |2          1          3|
///     |                       |
///     |0          0          1|
///   2 |-----------------------| 3
///     |2          1          3|
///     |                       |
///     |2          0          3|
///     |                       |
///     |0          0          1|
///     0 --------------------- 1
///                 0
/// ```
#[derive(Debug, Clone)]
pub struct Elem {
    pub id: usize,
    pub nodes: [usize; 4],
    pub edges: [usize; 4],
    pub element: Arc<Element>,
    pub h_levels: HLevels,
    pub poly_orders: PolyOrders,
    children: Option<SmallVec<[usize; 4]>>,
    ancestors: SmallVec<[(usize, HRefLoc); EXPECTED_NUM_H_REFINEMENTS]>,
}

impl Elem {
    /// Construct a new Elem from the relevant associative information
    pub fn new(id: usize, nodes: [usize; 4], edges: [usize; 4], element: Arc<Element>) -> Self {
        Self {
            id,
            nodes,
            edges,
            element,
            children: None,
            ancestors: SmallVec::new(),
            h_levels: HLevels::default(),
            poly_orders: PolyOrders::default(),
        }
    }

    /// Construct new 2 or 4 [ElemUninit]'s from an [HRef] of this Elem
    pub(crate) fn h_refine(
        &mut self,
        refinement: HRef,
        id_counter: &mut usize,
    ) -> Result<Vec<ElemUninit>, HRefError> {
        match self.children {
            Some(_) => Err(HRefError::ElemHasChildren(self.id)),
            None => {
                let children = refinement
                    .indices_and_ids(id_counter)
                    .map(|elem_idx_id| {
                        ElemUninit::new(
                            elem_idx_id,
                            refinement,
                            self.element.clone(),
                            self.id,
                            self.loc_stack(),
                            &self.h_levels,
                            self.poly_orders,
                        )
                    })
                    .collect::<Vec<ElemUninit>>();

                self.children = Some(children.iter().map(|ce| ce.id).collect());
                Ok(children)
            }
        }
    }

    /// Id of the Parent Elem if this Elem has a parent
    pub fn parent_id(&self) -> Option<usize> {
        self.ancestors.last().map(|(id, _)| *id)
    }

    /// Get the stack of [HRefLoc]s and Elem-IDs back to this `Elem`s ancestor on the base layer of the mesh
    ///
    /// This is useful for generating mappings between Elem's in parametric space
    pub fn loc_stack(&self) -> &[(usize, HRefLoc)] {
        &self.ancestors
    }

    /// Get the bounds of this `Elem` in parametric space relative to one of its ancestor `Elem`s
    pub fn relative_parametric_range(&self, from_ancestor: usize) -> [[f64; 2]; 2] {
        match self
            .ancestors
            .iter()
            .position(|(ancestor_id, _)| ancestor_id == &from_ancestor)
        {
            Some(starting_index) => self
                .ancestors
                .iter()
                .skip(starting_index)
                .fold([[-1.0, 1.0], [-1.0, 1.0]], |acc, (_, href_loc)| {
                    href_loc.sub_range(acc)
                }),
            None => panic!(
                "Elem {} is not an ancestor of Elem {}; cannot compute relative parametric range!",
                from_ancestor, self.id
            ),
        }
    }

    /// Get the bounds of this `Elem` in parametric space relative to its associated `Element`
    pub fn parametric_range(&self) -> [[f64; 2]; 2] {
        self.ancestors
            .iter()
            .fold([[-1.0, 1.0], [-1.0, 1.0]], |acc, (_, href_loc)| {
                href_loc.sub_range(acc)
            })
    }

    /// Get a reference to the Material properties in this region of the Mesh (Defined by this Elem's associated Element)
    pub fn get_materials(&self) -> &Materials {
        &self.element.materials
    }

    /// Gradients of a parametric point (as a [V2D]) through real space (via this Elem's parent [Element])
    pub fn parametric_mapping(&self, parametric_point: V2D, over_range: [[f64; 2]; 2]) -> M2D {
        self.element
            .parametric_mapping(parametric_point, over_range)
    }

    /// Returns a vector of ids for this Elem's children. Returns `None` if this Elem has no children.
    pub fn child_ids(&self) -> Option<SmallVec<[usize; 4]>> {
        self.children.clone()
    }

    /// Has this `Elem` been h-Refined
    pub fn has_children(&self) -> bool {
        self.children.is_some()
    }

    /// Produce a Json Object that describes this Elem
    #[cfg(feature = "json_export")]
    pub fn to_json(&self) -> JsonValue {
        object! {
            "id": self.id,
            "element_id": self.element.id,
            "parent": self.parent_id(),
            "active": self.children.is_none(),
            "nodes": array![self.nodes[0], self.nodes[1], self.nodes[2], self.nodes[3]],
            "edges": array![self.edges[0], self.edges[1], self.edges[2], self.edges[3]],
            "expansion": self.poly_orders,
            "h_levels": self.h_levels,
            "children": JsonValue::from(
                match &self.children {
                    Some(ids) => ids.to_vec(),
                    None => Vec::new(),
                }
            )
        }
    }
}

// Intermediate data structure used to represent a child [Elem] during the execution of an [HRef]
#[derive(Debug, Clone)]
pub(crate) struct ElemUninit {
    pub id: usize,
    pub nodes: [Option<usize>; 4],
    pub edges: [Option<usize>; 4],
    pub element: Arc<Element>,
    ancestors: SmallVec<[(usize, HRefLoc); EXPECTED_NUM_H_REFINEMENTS]>,
    h_levels: HLevels,
    poly_orders: PolyOrders,
}

impl ElemUninit {
    pub fn new(
        (idx, id): (usize, usize),
        refinement: HRef,
        element: Arc<Element>,
        parent_id: usize,
        parents_ancestors: &[(usize, HRefLoc)],
        parent_h_levels: &HLevels,
        poly_orders: PolyOrders,
    ) -> Self {
        let mut ancestors = SmallVec::from(parents_ancestors);
        ancestors.push((parent_id, refinement.loc(idx)));

        Self {
            id,
            nodes: [None; 4],
            edges: [None; 4],
            element,
            ancestors,
            h_levels: parent_h_levels.refined(refinement),
            poly_orders,
        }
    }

    pub fn set_node(&mut self, node_idx: usize, node_id: usize) {
        assert!(
            node_idx < 4,
            "Node indices must be between 0 and 4; cannot set Node ({}) to {} on ElemUninit {}!",
            node_idx,
            node_id,
            self.id
        );

        if let Some(current_id) = self.nodes[node_idx] {
            assert_eq!(
                current_id,
                node_id,
                "Node ({}) has already been set to {} on ElemUninit {}; Cannot set to {}",
                node_idx,
                self.nodes[node_idx].unwrap(),
                self.id,
                node_id
            );
        } else {
            self.nodes[node_idx] = Some(node_id);
        }
    }

    pub fn set_edge(&mut self, edge_idx: usize, edge_id: usize) {
        assert!(
            edge_idx < 4,
            "Edge indices must be between 0 and 4; cannot set Edge ({}) to {} on ElemUninit {}!",
            edge_idx,
            edge_id,
            self.id
        );

        assert!(
            self.edges[edge_idx].is_none(),
            "Edge ({}) has already been set to {} on ElemUninit {}; Cannot set to {}",
            edge_idx,
            self.edges[edge_idx].unwrap(),
            self.id,
            edge_id
        );

        self.edges[edge_idx] = Some(edge_id);
    }

    pub fn into_elem(self) -> Result<Elem, HRefError> {
        let nodes_init = self.nodes.iter().filter(|n| n.is_some()).count() == 4;
        let edges_init = self.edges.iter().filter(|e| e.is_some()).count() == 4;

        if nodes_init && edges_init {
            Ok(Elem {
                id: self.id,
                nodes: self
                    .nodes
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<usize>>()
                    .try_into()
                    .unwrap(),
                edges: self
                    .edges
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<usize>>()
                    .try_into()
                    .unwrap(),
                element: self.element,
                children: None,
                ancestors: self.ancestors,
                h_levels: self.h_levels,
                poly_orders: self.poly_orders,
            })
        } else {
            Err(HRefError::UninitializedElem(self.id))
        }
    }

    fn fmt_edge(&self, idx: usize) -> String {
        match self.edges[idx] {
            Some(id) => id.to_string(),
            None => String::from("_"),
        }
    }

    fn fmt_node(&self, idx: usize) -> String {
        match self.nodes[idx] {
            Some(id) => id.to_string(),
            None => String::from("_"),
        }
    }
}

impl fmt::Display for ElemUninit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ID: {} \t edges: [{}, {}, {}, {}] \t nodes: [{}, {}, {}, {}]",
            self.id,
            self.fmt_edge(0),
            self.fmt_edge(1),
            self.fmt_edge(2),
            self.fmt_edge(3),
            self.fmt_node(0),
            self.fmt_node(1),
            self.fmt_node(2),
            self.fmt_node(3),
        )
    }
}

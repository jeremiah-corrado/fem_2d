
/// The `Element` and `Materials` structures, which describe the physical properties of the associated `Elem`s
/// 
/// This includes a mapping between Real and Parametric spaces. 
///     (`real <--> parametric` mapping will be more useful once curvilinear Elements have been fully implemented)
/// 
/// JSON mesh files describe the `Element`s in the domain; not the `Elem`s
/// Upon `Mesh` construction, each `Element` has one associated `Elem`, but more can be added through h-Refinements
pub mod element;


/// `Elem`s are the basic finite unit in the `Mesh` 
/// 
/// They defer to their associated `Element` for descriptions of material parameters and mappings to Real Space
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
/// 1. *T-Type*:
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
/// 2. *U-Type*:
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
/// 3. *V-Type*:
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
pub mod elem;

/// A strait line between 2 `Node`s in parametric space
/// 
/// `Edge`s can be connected with 
pub mod edge;

/// `Node`s sit at the junction of at least 2 `Edges` and 2 `Elem`s
/// 
/// Many more `Elem`s can be associated with a `Node` as more refinement layers are added
pub mod node;


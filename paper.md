---
title: 'fem_2d: A Rust Package for 2D Finite Element Method Computations with Extensive Support for *hp*-refinement'
tags:
  - Rust
  - FEM
  - CEM
  - hp-refinement
  - PDE Solvers
authors:
  - name: Jeremiah Corrado^[first author]
    orcid: 0000-0003-2688-0600
    affiliation: 1
affiliations:
 - name: Colorado State University, Department of Electrical and Computer Engineering
   index: 1
date: 29 January 2022
bibliography: paper.bib
---

# Summary

The Finite Element Method (FEM) is a powerful tool used to solve Partial Differential Equations (PDE)s on arbitrary geometries. Although the method has existed for some time, it is constantly being improved upon and applied to new problems in engineering and science. Some common applications include the Navier-Stokes equations which govern fluid dynamics, Schrödinger's equation which governs the evolution of quantum systems, and Maxwells Equations which describe the propagation of Electromagnetic waves.

The `fem_2d` library is a powerful FEM simulation tool implemented in Rust. It was designed within the domain of Computational Electromagnetics (CEM) with a specific focus on solving Maxwells Equations on 2D waveguide cross-sections with very high accuracy. It's functionality is easily extendable to other domains using Rust `Traits` (which are similar to C++20 `Concepts` or Java `Interfaces`). Traits are a highly expressive form of genericism, which allow functions to act on any data structure so long as it implements some shared functionality. In this case, the shared functionality is the ability to compute an integral between a basis function and testing function. As such, the libraries entire functionality can be leveraged against a generic PDE problem simply implementing fem_2d's `Integral` trait to express the variational form of the problem of interest.

`fem_2d` has all the functionality needed to formulate and solve Generalized Eigenvalue Problems. This includes a solver based on the popular pure-Rust linear algebra library Nalgebra[@nalgebra], as well as an external sparse solver implemented using SLEPc [@slepc, @petsc-web-page, @petsc-user-ref, @petsc-efficient]. It also has extensive functionality for generating `.vtk` files of solutions which can be plotted using external tools. 

# Statement of Need
Fem_2D's primary advantage over other FEM libraries, such as the Deal.II library [@dealII93], is its highly dynamic and expressive *hp*-refinement scheme. Unlike many other quadrilateral-element FEM packages, `fem_2d` supports n-irregular anisotropic *h*-refinement as well as anisotropic *p*-refinement. In other words, there is no fundamental limitation on the shape, location, or orientation when adding new elements to the Mesh. The polynomial expansion orders of the Basis Functions associated with each element can also be modified separately in each direction. 

The Isotropic *h*-refinement API alone, is a useful, and even necessary, feature when computing solutions over geometries with sharp edges or material discontinuities, as these tend to introduce very small-scale solution behavior which is otherwise addressed only by excessive amounts of *p*-refinement `[@corrado_harmon_notaros_2021]`. The addition of anisotropic *h*- and *p*-refinement presents an even lager capacity for improved solution efficiency, as small-scale behavior is targeted more directly and redundant Degrees of Freedom can are left out of the system `[@corrado_harmon_notaros_2021]`. 

## Examples of *hp*-Refinement:

The following code sample gives a few examples of how the `Mesh` structures *hp*-refinement methods might be used in practice:
```Rust
use fem_2d::prelude::*;
let mut mesh = Mesh::from_file("./some_mesh.json").unwrap();

// isotopically h-refine all elements
mesh.global_h_refinment(HRef::t()).unwrap();;

// anisotropically p-refine all elements (+2 in the u-direction, +4 in the v-direction)
mesh.global_p_refinement(PRef::from(2, 4)).unwrap();;

// anisotropically h-refine all elements connected to some target_node in the v-direction
let target_node_id = 5;

mesh.h_refine_with_filter(|elem| {
    if elem.nodes.contains(&target_node_id) {
        Some(HRef::v())
    } else {
        None
    }
}).unwrap();

// positively p-refine all elements on the border of the mesh,
// and negatively p-refine all other elements
mesh.p_refine_with_filter(|elem| {
    if elem.edges.iter().any(|edge_id| {
        mesh.edges[*edge_id].is_boundary()
    }) {
        Some(PRef::from(1, 1))
    } else {
        Some(Pref::from(-1, -1))
    }
}).unwrap();

```

## Example of Problem Definition

The Maxwell eigenvalue problem has the following Galerkin Sampled formulation for an arbitrary domain terminated with dirichlet boundary conditions:

>Find a solution: $ \mathrm{U} = \{{\mathbf{u}}, \lambda \} \in B_0 \times \Bbb{R} $ which satisfies:
>
> $$ b(\mathbf{u}, \phi) = \lambda a(\mathbf{u}, \phi) \quad \forall \phi \in B_0 $$ 
>$$\text{where: } \left\{\begin{array}{l}
a(\mathbf{u}, \phi) = \langle \nabla_t \times \mathbf{u}, \nabla_t \times \phi \rangle \cr
b(\mathbf{u}, \phi) = \langle \mathbf{u}, \phi \rangle \cr
B_0 \subset H(\text{curl}; \Omega)
\end{array}\right.$$

The system matrices for this equation are populated in `fem_2d` using a call to `galerkin_sample_gep` with three generic arguments:
```Rust
// Generate a domain (Ω) from a Mesh which has been refined as necessary
let domain = Domain::from(mesh);

// Formulate a generalized Eigenvalue problem
let gep = domain.galerkin_sample_gep::<KOLShapeFn, CurlCurl, L2Inner>(None);
```
 * The first generic argument, which must implement the `ShapeFn` Trait represents the curl-conforming basis $B_0$. In this case `KOLShapeFn` is used.
 * The second generic argument, which implements the `Integral` Trait defines the function $a(\mathbf{u}, \phi)$. The `CurlCurl` Struct implements the L2 inner product of the Curl of two Basis Functions. 
 * The third generic argument, which also implements the `Integral` Trait defines the function $b(\mathbf{u}, \phi)$. The `L2Inner` Struct implements the inner product of two Basis Functions.
 
 The eigenvalue problem can be solved for some target eigenvalue with the following code:
 ```Rust
 // Dense solution (not recommended for large problems)
let solution = nalgebra_solve_gep(gep, target_eigenvalue).unwrap();

// OR: Sparse solution (requires external Slepc solver)
let solution = slepc_solve_gep(gep, target_eigenvalue).unwrap();
 ```

In summary, any eigenvalue problem can be solved by creating a custom implementation of the `Integral` trait which is described in detail in the documentation. A custom Basis set can also be used by implementing the `ShapeFn` trait which is also described in the documentation.

# References

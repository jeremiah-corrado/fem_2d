---
title: "FEM_2D: A Rust Package for 2D Finite Element Method Computations with Extensive Support for *hp*-refinement"
tags:
  - Rust
  - FEM
  - CEM
  - hp-refinement
  - PDE Solvers
authors:
  - name: Jeremiah Corrado
    orcid: 0000-0003-2688-0600
    affiliation: 1
  - name: Jake J. Harmon
    affiliation: 1
  - name: Milan M. Ilic
    affiliation: "1, 2"
  - name: Branislav M. Notaroš
    affiliation: 1
affiliations:
 - name: "Colorado State University; Department of Electrical and Computer Engineering"
   index: 1
 - name: "University of Belgrade; School of Electrical Engineering"
   index: 2
date: 14 February 2022
bibliography: paper.bib
---

# Introduction

The Finite Element Method (FEM) is a powerful computation framework used to solve Partial Differential Equations (PDE)s on arbitrary geometries. In reality, physical systems behave in a continuous manner (both in space and time); however FEM solvers are able to model these dynamics with a high fidelity by decomposing a physical model into a finite set of elements. Each element supports a finite number of degrees of freedom, which are used to describe the behavior of the system. This way, the mathematically continuous dynamics can be expressed in terms of a system of linear equations. Linear algebra tools are then used to solve the problem such that the PDE is satisfied along with some boundary conditions (on the border of the Domain) and some continuity conditions (between neighboring elements). 

Some common PDE's include the Navier-Stokes equations which characterize the behavior of fluids, Schrödinger's equation which governs the evolution of quantum systems, and Maxwell's Equations which are a macroscopic description of essentially all Electromagnetic phenomena. The ability to accurately and efficiently model these differential equations and others is imperative to the success of many engineering projects and scientific endeavors. Most of the technology that engineers are interested in developing has far exceeded the reach of direct mathematical analysis, and thus computational tools such as FEM are used to drive technological development forward.

As such, innovations in FEM have a direct impact on myriad engineering disciplines. The more efficient, accurate, and feature rich, we can make simulation tools, the more beneficial they will be to industrial and scientific applications. This is the motivation force behind academic work in the field of FEM. The `FEM_2D` library is a Rust package that aims to enable further research into a particular FEM innovation called Refinement-by-Superposition (RBS). The related research papers ( [@corrado:2021], [@harmon:2021]) explore benefits of RBS using the 2D Maxwell Eigenvalue Problem as a proving ground.



The `FEM_2D` library is a powerful FEM simulation tool implemented in Rust. It was designed within the domain of Computational Electromagnetics (CEM) with a specific focus on solving Maxwell's Equations over a 2D waveguide cross-sections with very high accuracy. This problem has practical engineering applications for the design of waveguides; however, it is used here as a research tool for exploring and evaluating improvements upon FEM. Particularly, this library employs the Refinement-by-Superposition (RBS) approach to $hp$-refinement. RBS is a novel and relatively unexplored generalization of the more typical Refinement-by-Replacement approach which unlocks a much simpler pathway to enforcing continuity conditions. 

Although `FEM_2D` focuses on the Maxwell Eigenvalue problem, it's functionality is intended to extend easily to other domains using an API based on Rust `Traits` (which are like C++20 `Concepts` or Java `Interfaces`). Traits are a highly expressive form of genericism which allow functions to act on any data structure so long as it implements some shared set of functionality. In this case, the generic functionality is the ability to compute an integral over an expression of two overlapping basis functions. As such, the library's functionality can be leveraged against a given PDE problem by implementing FEM_2D's `Integral` trait to express the variational form of the problem of interest.

`FEM_2D` has all the functionality needed to formulate and solve Generalized Eigenvalue Problems. This includes an easy-to-use Mesh API which handles Mesh instantiation  and refinement. It also includes a Domain API which defines a Generic set of Basis Functions over the Mesh while ensuring adherence to continuity and boundary conditions. Additionally, there are two Eigenvalue-Problem solvers. The first is implemented using Nalgebra [@nalgebra]: a powerful linear algebra library written in Rust. The second is an external sparse solver --for large or poorly conditioned problems-- implemented using SLEPc [@slepc] [@petsc-web-page] [@petsc-user-ref] [@petsc-efficient]. The library also has extensive functionality for generating `.vtk` files of solutions which can be plotted using external tools.

# Statement of Need

Efficiently computing FEM solutions over geometries with sharp edges or stark material discontinuities necessitates *hp*-refinement (whether isotropic or anisotropic). These situations tend to introduce multi-scale solution behavior which is challenging to model with pure *p*- or pure *h*-refinements, motivating combined hp-refinements [@harmon:2021]. Within the class of hp-refinements, the addition of anisotropic *hp*-refinements (over isotropic ones) presents a significantly larger capacity for solution efficiency, as small-scale behavior is targeted more directly and ineffectual Degrees of Freedom are left out of the system [@corrado:2021]. Increased efficiency in terms of the number of degrees of freedom is a key factor in the speed of large scale simulations, as well as the applicability of the method to smaller scale hardware (such as personal computers). Thus, a feature-rich anisotropic *hp*-refinement API is needed to enable efficient solution of challenging FEM problems. This is directly afforded by the underlying RBS methodology.

In addition to the features themselves, a research package is best utilized when the implementation is straightforward and easy to understand. This way, other researchers can use the code to validate the methodology itself, and to serve as a starting point for additional investigation. This is yet another benefit of the RBS approach, as it greatly simplifies the enforcement of continuity conditions, which is typically the most challenging aspect of an *h*-refinement implementation over quadrilateral or hexahedral elements.

# Project Goals

`FEM_2D`: is intended to be the following:

* A proof-of-concept for the RBS method and associated research
* A high quality implementation of the RBS approach to FEM for the purpose of encouraging additional related research
    * This could involve additions to the `FEM_2D` library to support new features or application domains
    * It could also simply serve as a starting point, or example implementation, for new RBS based FEM software
* A powerful solver for the 2D Maxwell Eigenvalue Problem in itself (or any other problems that are implemented by the open-source community)

# Features

## *hp*-Refinement API:

FEM_2D's primary advantage over other FEM libraries, such as the Deal.II library [@dealII93], is its highly dynamic and expressive *hp*-refinement API. Unlike many other quadrilateral-element FEM packages, `FEM_2D` supports n-irregular anisotropic *h*-refinement as well as anisotropic *p*-refinement. In other words, there are far fewer limitations on the shape, location, or orientation of new elements when adding them to the Mesh. The polynomial expansion orders of the Basis Functions associated with each element can also be modified separately in each direction. This level of freedom would not be possible without the underlying RBS methodology. 

The following example shows how some of the $h$-refinement methods may be used to modify a mesh structure. It is important to note that there are three primary $h$-refinement types which are designated by the `HRef` enum:

* T - isotropic: produces 4 child elements
* U - anisotropic in the u-direction: produces 2 child elements
* V - anisotropic in the v-direction: produces 2 child elements

There are also two sub-types associated with the U and V refinements which invoke a subsequent anisotropic refinement on one of the two child elements in the opposite direction. These are constructed with `HRef::U(Some(child_index))` and `HRef::V(Some(child_index))` respectively, where `child_index` must be either 0 or 1. 

It is also important to note that the `global_h_refinement` and `h_refine_with_filter` methods will only apply refinements to Elements that are eligible for $h$-refinement (i.e., they must be leaf elements and the length of each of their edges must be above a minimum threshold). Alternatively, the methods that expose more explicit control (`h_refine_elems` and `execute_h_refinements`) can return an error if one of the specified elements is not eligible for $h$-refinement. A detailed explanation of the error types is provided in the documentation.

```rust
use fem_2d::prelude::*;
use std::error::Error;

fn do_some_h_refinements(mesh_file_path: &str) -> Result<Mesh, Box<dyn Error>> {
    let mut mesh = Mesh::from_file(mesh_file_path)?;

    // isotropically h-refine all elems
    mesh.global_h_refinement(HRef::T);

    // anisotropically h-refine all elems connected to some target node
    let target_node_id = 5;
    mesh.h_refine_with_filter(|elem| {
        if elem.nodes.contains(&target_node_id) {
            Some(HRef::u())
        } else {
            None
        }
    });

    // anisotropically h-refine a list of elems by id
    mesh.h_refine_elems(vec![3, 4, 8, 12], HRef::v())?;

    // directly apply a list of refinements to the mesh
    mesh.execute_h_refinements(vec![
        (1, HRef::T),
        (5, HRef::U(Some(0))),
        (6, HRef::U(Some(1))),
        (10, HRef::V(None)),
    ])?;

    Ok(mesh)
}
```
The following example shows how some of the $p$-refinement methods may be used here, the mesh to be modified is provided as an argument rather than being loaded from a file. 

$P$-refinements are constructed from the `PRef` type using a pair of `i8`'s (8-bit signed integers). As such, any element's u- and v-directed expansion orders can be modified independently in either the positive or negaitve direction.

The behavior of these methods is straightforward with the slight caveat that the `global_p_refinement` and `p_refine_with_filter` methods will guard against any refinement pushing an `Elem` outside of its valid expansion order range. Specifically, refinements are clamped element-wise to ensure that the final expansion order is in the range [1, 20]. The $p$-refinement methods that can return an error (those followed by a `?` in the example) do not exhibit this behavior. This is in keeping with the design of the $h$-refinement API in the sense that methods with less explicit control are safer, while the more explicit methods allow for failure. 

```rust
use fem_2d::prelude::*;

fn do_some_p_refinements(mesh: &mut Mesh) -> Result<(), PRefError> {
    // isotropically p-refine all elems (with a magnitude 2 refinement)
    mesh.global_p_refinement(PRef::from(2, 2));

    // positively p-refine all "leaf" elems (with a magnitude 1 refinement)
    // negatively p-refine all other elems (with a magnitude -1 refinement)
    mesh.p_refine_with_filter(|elem| {
        if elem.has_children() {
            Some(PRef::from(-1, -1))
        } else {
            Some(PRef::from(1, 1))
        }
    });

    // anisotropically p-refine a list of elems by id
    mesh.p_refine_elems(vec![3, 4, 8, 12], PRef::from(4, 2))?;

    // directly apply a list of refinements to the mesh
    mesh.execute_p_refinements(vec![
        (1, PRef::from(3, 2)),
        (5, PRef::from(0, 1)),
        (6, PRef::from(-1, -1)),
        (10, PRef::from(4, -2)),
    ])?;

    Ok(())
}
```


The `Mesh` data structure also has an alternative set of methods to modify expansion orders by setting them directly rather than additively. These methods can be very useful in scenarios where it does not matter what the current expansion orders are, and an element needs to have a specific expansion order which is either known beforehand or computed ad-hoc. The following example juxtaposes some of the functionality with the traditional $p$-refinement API.

Here, both methods can return an error, as it is possible to specify an invalid set of expansion orders. These methods take a length-two array of `u8`'s (8-bit unsigned integers), and thus preemptively remove the possibility of setting negative expansion orders, however, they still Err on expansion orders that are zero or too large. 
```rust
use fem_2d::prelude::*;

fn set_some_expansion_orders(mesh: &mut Mesh) -> Result<(), PRefError> {
    // set the expansion order on all elems to (3, 3)
    mesh.set_global_expansion_orders([3, 3])?;

    // set the expansion orders to (4, 4) on all "leaf" elems
    // set the expansion orders to (2, 2) on all other elems
    mesh.set_expansions_with_filter(|elem| {
        if elem.has_children() {
            Some([2, 2])
        } else {
            Some([4, 4])
        }
    })?;

    Ok(())
}

```

## Problem Formulation and Solution

The following example shows how a simplified formulation of the Maxwell Eigenvalue Problem maps to the corresponding code in the library. This is intended provide a general depiction of how one might translate a mathematical problem into an `FEM_2D` implementation. 

The Maxwell eigenvalue problem has the following Continuous-Galerkin formulation for an arbitrary Domain terminated with Dirichlet boundary conditions, (constraining the solution to TE modes only):

>Find a solution: \begin{equation} \label{eq:solution} \quad \text{U} = \{{\mathbf{u}}, \lambda \} \in B_0 \times \Bbb{R} \quad \end{equation} which satisfies:
> \begin{equation} \label{eq:formulation} b(\mathbf{u}, \phi) = \lambda a(\mathbf{u}, \phi) \quad \forall \phi \in B \end{equation}
>\begin{equation} \label{eq:gen_args} \text{where: } \left\{\begin{array}{l}
B \subset H_0(\text{curl}; \Omega) \cr
a(\mathbf{u}, \phi) = \langle \nabla_t \times \mathbf{u}, \nabla_t \times \phi \rangle \cr
b(\mathbf{u}, \phi) = \langle \mathbf{u}, \phi \rangle
\end{array}\right.\end{equation}

The Generalized Eigenvalue Problem is built from a `Mesh` with the following code.
```rust
use fem_2d::prelude::*;
use rayon::prelude::*;

fn problem_from_mesh(mesh: Mesh) -> Result<GEP, GalerkinSamplingError> {
  // Setup a global thread-pool for parallelizing Galerkin Sampling
  rayon::ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();

  // Generate a Domain (Ω) from a Mesh with H(Curl) Continuity Conditions
  let domain = Domain::from(mesh, ContinuityCondition::HCurl);

  // Compute a Generalized Eigenvalue Problem (composed of the System Matrices A, B)
  let gep = galerkin_sample_gep_hcurl::<
      HierPoly, // Basis Space
      CurlCurl, // Stiffness Integral
      L2Inner,  // Mass Integral
  >(&domain, Some([8, 8]))
}
```

The `Domain` structure represents the entire FEM domain, including the discretization and the basis space which conforms to the provided continuity condition (only H(Curl) is currently implemented; however, a framework is in place for implementing H(Div) and other continuity conditions). 

Galerkin sampling is then executed in parallel over the Domain, generating a Generalized Eigenvalue Problem. The Domain and a Gauss-Legendre-Quadrature grid size are provided as arguments. This function may also return an Error, if the Galerkin Sampling fails due to an ill-posed problem.

The three generic arguments -- designated with the turbofish operator (`::<>`) -- correspond to the three lines of \autoref{eq:gen_args}. The basis space can be swapped for any other space that implements the `HierCurlBasisFnSpace` Trait. `HierPoly` is a relatively simple implementation composed of exponential functions. A more sophisticated basis space: `HierMaxOrtho` can be included using the `max_ortho_basis` Feature Flag. Custom Basis Spaces can also be created by implementing the Trait. 

The `CurlCurl` and `L2Inner` integrals, which correspond to the Stiffness and Mass matrices respectively, can be swapped for any other structure that implements the `HierCurlIntegral` Trait. This generic interface allows users to leverage the galerkin sampling functionality against other curl-conforming problems.^[The provided functionality is obviously somewhat incomplete, as only Curl Conforming problems can be solved; however, the library's module-structure and trait-hierarchy provide a clear template for the analogous H(Div) implementation. There is also room for other galerking sampling and integration functionality associated with alternate continuity conditions. These methods, structures, and traits should require minimal changes to the `Domain` structure, and no changes to the `Mesh` structure.] 

The Generalized Eigenvalue Problem, can then be solved using one of the available solvers:
```rust
// Dense solution (not recommended for large problems)
let eigenpair = nalgebra_solve_gep(gep, target_eigenvalue).unwrap();

// OR: Sparse solution (requires external Slepc solver)
let eigenpair = slepc_solve_gep(gep, target_eigenvalue).unwrap();
```
The sparse [`SLEPc`](https://slepc.upv.es/) solver requires an [external solver](https://github.com/jeremiah-corrado/slepc_gep_solver) to be downloaded and compiled. This is highly recommended for larger problems, or problems with ill-conditioned B-Matrices.

Both solvers look for the eigenvalue closest to the target value. They can return errors if the solution does not converge. Upon success, the returned eigenpair contains the eigenvalue and eigenvector with length equal to the number of degrees of freedom in the domain.

## Field Visualization

# References

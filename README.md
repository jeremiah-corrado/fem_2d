# fem_2d

A Rust library for 2D Finite Element Method computations, featuring:

- Highly felxible *hp*-Refinement
  - Isotropic & Anisotropic *h*-refinements (with support for n-irregularity)
  - Isotropic & Anisotropic *p*-refinements 
- Generic shape function evaluation
  - You can use one of the two built in sets of Shape Functions
  - Or you can define your own by implementing the `ShapeFn` Trait
- Generic eigenvalue Problem construction
  - You can use the built in Integrals (for the Maxwell Eigenvalue Problem)
  - Or you can define your own problem by implementing the `Integral` Trait
- Two Eigensolvers
  - Sparse: Using an external Slepc Solver (installation instructions found [here](https://github.com/jeremiah-corrado/slepc_gep_solver)
  - Dense: Using [Nalgebra](https://nalgebra.org/)'s Eigen Decomposition (not recommended for large problems)
- Expressive Solution Evaluation
  - Field solutions can easily be generated from an eigenvector
  - Functions of solutions can be evaluated
  - Solutions can be printed to `.vtk` files for plotting (via [VISIT](https://visit-dav.github.io/visit-website/index.html) or similar tools)

## Example

Solve the Maxwell Eigenvalue Problem on a standard Waveguide and print the Electric Fields to a VTK file
```Rust
use fem_2d::prelude::*;

// Load a standard air-filled waveguide mesh from a JSON file
let mut mesh = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();

// Set the polynomial expansion order to 4 in both directions on all Elems
        mesh.set_global_expansion_orders([4, 4]).unwrap();

// Isotropically refine all Elems
mesh.global_h_refinement(HRef::t()).unwrap();

// Then anisotropically refine the resultant Elems in the center of the mesh
let center_elem_id = mesh.elems[0].edges[3];
mesh.h_refine_with_filter(|elem| {
    if elem.edges.contains(&center_elem_id) {
        Some(HRef::u())
    } else {
        None
    }
}).unwrap();
        
// Construct a domain with Dirichlet boundary conditions
let domain = Domain::from_mesh(mesh);
println!("Constructed Domain with {} DoFs", domain.dofs.len());

// Construct a generalized eigenvalue problem for the Electric Field
    // (in parallel using the Rayon Global ThreadPool)
let gep = domain.galerkin_sample_gep_parallel::<KOLShapeFn, CurlCurl, L2Inner>(None);

// Solve the generalized eigenvalue problem using Nalgebra's Eigen-Decomposition
    // look for an eigenvalue close to 10.0
let solution = nalgebra_solve_gep(gep, 10.0).unwrap();
println!("Found Eigenvalue: {:.15}", solution.value);

// Construct a solution-field-space over the Domain with 64 samples on each "leaf" Elem
let mut field_space = UniformFieldSpace::new(&domain, [8, 8]);

// Compute the Electric Field in the X- and Y-directions (using the same ShapeFns as above)
let e_field_names = field_space.xy_fields::<KOLShapeFn>("E", solution.normalized_eigenvector()).unwrap();

// Compute the magnitude of the Electric Field
field_space.expression_2arg(e_field_names, "E_mag", |ex, ey| (ex.powi(2) + ey.powi(2)).sqrt()).unwrap();

// Print "E_x", "E_y" and "E_mag" to a VTK file
field_space.print_all_to_vtk("./test_output/electric_field_solution.vtk").unwrap();
```

## Mesh Refinement

***h*-Refinements** are implemented using the Refinement by Superposition (RBS) method
> Technical details can be found in this paper: [A Refinement-by-Superposition Approach to FullyAnisotropichp-Refinement for Improved Efficiencyin CEM](https://www.techrxiv.org/articles/preprint/A_Refinement-by-Superposition_Approach_to_Fully_Anisotropic_hp-Refinement_for_Improved_Efficiency_in_CEM/16695163)

Three types of h-refinement are supported:
* **T**: Elements are superimposed with 4 equally sized child elements
* **U**: Elements are superimposed with 2 child elements, such that the resolution is improved in the x-direction
* **V**: Elements are superimposed with 2 child elements, such that the resolution is improved in the y-direction

These are designated as an Enum: `HRef`, located in the `h_refinements` module

___________________

Mesh coarsening is not currently supported

***p*-Refinements** allow elements to support large expansion orders in the X and Y directions. These can be modulated separately for greater control over resource usage and solution accuracy.

Expansion orders can be increased or reduced by constructing a `PRef`, located in the `p_refinements` module

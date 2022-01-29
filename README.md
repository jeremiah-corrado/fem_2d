# fem_2d

A rust library for 2D Finite Element Method computations, featuring:

- Highly felxible hp-Refinement
  - Isotropic & Anisotropic h-refinements (with support for n-irregularity)
  - Isotropic & Anisotropic p-refinements 
- Generic shape function evaluation
  - You can use one of the two built in sets of Shape Functions
  - Or you can define your own by implementing the `ShapeFn` Trait
- Generic eigenvalue Problem construction
  - You can use the built in Integrals (for the Maxwell Eigenvalue Problem)
  - Or you can define your own problem by implementing the `Integral` Trait
- Two Eigensolvers
  - Sparse: Using an external Slepc Solver (installation instructions found [here](https://github.com/jeremiah-corrado/slepc_gep_solver)]
  - Dense: Using [Nalgebra](https://nalgebra.org/)'s Eigen Decomposition (not recommended for large problems)
- Expressive Solution Evaluation
  - Field solutions can easily be generated from an eigenvector
  - Functions of solutions can be evaluated
  - Solutions can be printed to `.vtk' files for plotting (via [VISIT](https://visit-dav.github.io/visit-website/index.html) or similar tools)


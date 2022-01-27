extern crate cxx;
extern crate rayon;
extern crate bytes;
extern crate nalgebra;

mod eigenproblem;
mod slepc_wrapper;

pub use eigenproblem::{SparseMatrix, GEP, AIJMatrixBinary, retrieve_solution, solve_eigenproblem};
pub use slepc_wrapper::{solve_gep, EigenPair};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_binary_files() {
        let mut sm = SparseMatrix::new(5);

        sm.insert([0, 0], 1.0);
        sm.insert([1, 1], 2.0);
        sm.insert([2, 2], 3.0);
        sm.insert([3, 3], 4.0);
        sm.insert([4, 4], 5.0);
        sm.insert([2, 3], 0.5);
        sm.insert([4, 1], 0.25);

        let sm_aij : AIJMatrixBinary = sm.into();

        sm_aij.to_petsc_binary_format("../../test_output/test_matrix.dat").unwrap();
    }

    #[test]
    fn get_solution_from_files() {
        let epair = retrieve_solution("../../test_input/test").unwrap();
        println!("Eigenvalue: {:15}", epair.value);
        println!("Eigenvector length: {}", epair.vector.len());
    }
}

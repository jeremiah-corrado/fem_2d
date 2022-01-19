extern crate cxx;
extern crate rayon;

mod eigenproblem;
mod slepc_wrapper;

pub use eigenproblem::{GEP, SparseMatrix};
pub use slepc_wrapper::{solve_gep, EigenPair};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

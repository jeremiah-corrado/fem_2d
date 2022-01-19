extern crate cxx;
extern crate rayon;

mod eigenproblem;

pub use eigenproblem::{GEP, AIJMatrix, SparseMatrix};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

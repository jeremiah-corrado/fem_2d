extern crate num_complex;
extern crate smallvec;
extern crate json;

mod domain;

pub use domain::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

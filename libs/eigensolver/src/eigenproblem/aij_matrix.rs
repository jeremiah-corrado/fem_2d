
/// AIJ Sparse Matrix Format native to [PETSC](https://petsc.org/main/docs/manualpages/Mat/MatCreateSeqAIJWithArrays.html#MatCreateSeqAIJWithArrays)
pub struct AIJMatrix {
    /// i[row_idx] = i[row_idx - 1] + number of entries on row i
    pub i: Vec<i32>,
    /// columns of elements in 'a'
    pub j: Vec<i32>,
    /// matrix entries sorted by row
    pub a: Vec<f64>,
    /// symmetric matrix dimension
    pub dim: usize,
}

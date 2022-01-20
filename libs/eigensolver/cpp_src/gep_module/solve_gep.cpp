#include <vector>
#include <array>
#include <mpi.h>

#include <slepceps.h>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "petscmat.h"

void early_exit(int status) {
    SlepcFinalize();
    exit(status);
}

void get_matrix_data(std::string mat_name, std::vector<PetscScalar> &vals, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols, PetscInt dim, PetscInt num_vals);
void build_matrix(Mat &X, int M, std::vector<PetscScalar> &vals, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols);
double solve_eigenproblem(Mat& A, Mat& B, Vec& xr, Vec& xi, double target_eval);
void deposit_solution(std::vector<double> evals);
void deposit_eigenvector(Vec& xr, int size);

int main(int argc, char * argv[]) {
    PetscErrorCode ierr = SlepcInitialize(&argc, &argv, (char *) 0, NULL);
    if (ierr) early_exit(1);

    PetscScalar target_eval = 1.00;
    PetscInt dimension = 0;
    PetscInt num_values = 0;
    PetscOptionsGetScalar(NULL, NULL, "-a", &target_eval, NULL);
    PetscOptionsGetInt(NULL, NULL, "-d", &dimension, NULL);
    PetscOptionsGetInt(NULL, NULL, "-v", &num_values, NULL);

    auto a_vals = std::vector<PetscScalar>();
    auto a_rows = std::vector<PetscInt>();
    auto a_cols = std::vector<PetscInt>();

    auto b_vals = std::vector<PetscScalar>();
    auto b_rows = std::vector<PetscInt>();
    auto b_cols = std::vector<PetscInt>();

    get_matrix_data("A", a_vals, a_rows, a_cols, dimension, num_values);
    get_matrix_data("B", b_vals, b_rows, b_cols, dimension, num_values);

    Mat A, B;
    build_matrix(A, dimension, a_vals, a_rows, a_cols);
    build_matrix(B, dimension, b_vals, b_rows, b_cols);
    
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    MatView(B, PETSC_VIEWER_STDOUT_WORLD);

    Vec xr, xi;
    MatCreateVecs(A, NULL, &xr);
    MatCreateVecs(A, NULL, &xi);
    double best_eval = solve_eigenproblem(A, B, xr, xi, target_eval);

    int thread_id;
    MPI_Comm_rank(PETSC_COMM_WORLD, &thread_id);
    if (thread_id == 0) deposit_solution({best_eval});

    PetscBarrier(PETSC_NULL);
    deposit_eigenvector(xr, dimension);

    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&xr);
    VecDestroy(&xi);

    SlepcFinalize();
    return 0;
}

std::array<std::string, 4> mem_channel_names(std::string name) {
    return {
            name + "_mat_vals",
            name + "_mat_rows",
            name + "_mat_cols",
    };
}

void get_matrix_data(std::string mat_name, std::vector<PetscScalar> &vals, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols, PetscInt dim, PetscInt num_vals) {
    using namespace boost::interprocess;
    auto mem_names = mem_channel_names(mat_name);

    int num_threads;
    int thread_id;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_threads);
    MPI_Comm_rank(PETSC_COMM_WORLD, &thread_id);

    if (num_threads > 1) {
        PetscPrintf(PETSC_COMM_WORLD, "AIJ Solver is only setup to use 1 MPI Thread!");
        early_exit(1);
    }

    shared_memory_object vals_smo (open_only, &mem_names[0][0], read_only);
    shared_memory_object rows_smo (open_only, &mem_names[1][0], read_only);
    shared_memory_object cols_smo (open_only, &mem_names[2][0], read_only);

    mapped_region vals_reg (vals_smo, read_only, 0, num_vals * sizeof(double));
    mapped_region rows_reg (rows_smo, read_only, 0, (dim + 1) * sizeof(int));
    mapped_region cols_reg (cols_smo, read_only, 0, num_vals * sizeof(int));

    auto* vals_array = static_cast<double*>(vals_reg.get_address());
    auto* rows_array = static_cast<int*>(rows_reg.get_address());
    auto* cols_array = static_cast<int*>(cols_reg.get_address());

    vals.assign(vals_array, vals_array + vals_reg.get_size() / sizeof(double));
    rows.assign(rows_array, rows_array + rows_reg.get_size() / sizeof(int));
    cols.assign(cols_array, cols_array + cols_reg.get_size() / sizeof(int));

    // for (int r = 1; r < rows.size(); r++) {
    //     PetscPrintf(PETSC_COMM_WORLD, "%d \t %d\n", rows[r - 1], rows[r]);
    // }

    // PetscPrintf(PETSC_COMM_WORLD, "\n");
    // for (int c = 0; c < cols.size(); c++) {
    //     if (cols[c] >= M || cols[c] < 0) {
    //         PetscPrintf(PETSC_COMM_WORLD, "%d \t %g -------------------------\n", cols[c], vals[c]);
    //     } else {
    //         PetscPrintf(PETSC_COMM_WORLD, "%d \t %g\n", cols[c], vals[c]);
    //     }
    // }
}

void build_matrix(Mat &X, int M, std::vector<PetscScalar> &vals, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols) {
    PetscErrorCode ierr;

    ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, M, M, &rows[0], &cols[0], &vals[0], &X); if (ierr) early_exit(3);
    
    // ierr = MatSetFromOptions(X); if (ierr) early_exit(3);
    // ierr = MatSetUp(X); if (ierr) early_exit(3);

    // ierr =  MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY); if (ierr) early_exit(5);
    // ierr =  MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY); if (ierr) early_exit(5);

    // PetscPrintf(PETSC_COMM_WORLD, "\nMatrix Values: \n");
    // PetscViewer viewer;
    // PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    // PetscViewerSetType(viewer, PETSCVIEWERASCII);
    // MatView(X, viewer);

    // PetscBool sym = PETSC_FALSE;
    // MatIsSymmetric(X, 1e-10, &sym);

    // if (sym) {
    //     PetscPrintf(PETSC_COMM_WORLD, "\nSymmetric");
    // } else {
    //     PetscPrintf(PETSC_COMM_WORLD, "\nAsymmetric");
    // }
}

double solve_eigenproblem(Mat& A, Mat& B, Vec& xr, Vec& xi, double target_eval) {
    PetscErrorCode ierr;
    EPS eps;
    ST st;
    KSP ksp;
    PC pc;
    PetscInt nconv;
    PetscScalar eigr, eigi;
    double eval_solution = 0.0;

    ierr = EPSCreate(PETSC_COMM_WORLD, &eps); if (ierr) early_exit(6);
    ierr = EPSSetOperators(eps, A, B); if (ierr) early_exit(6);
    ierr = EPSSetProblemType(eps, EPS_GHEP); if (ierr) early_exit(6);
    ierr = EPSSetFromOptions(eps); if (ierr) early_exit(6);

    ierr = EPSSetTolerances(eps, 1.0e-15, 100); if (ierr) early_exit(7);
    ierr = EPSSetType(eps, "krylovschur"); if (ierr) early_exit(7);
    
    ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);if (ierr) early_exit(7);
    
    ierr = EPSSetTarget(eps, target_eval);if (ierr) early_exit(8);
    ierr = EPSGetST(eps,&st);if (ierr) early_exit(8);
    ierr = STSetType(st,"sinvert"); if (ierr) early_exit(8);
    ierr = STSetShift(st, target_eval); if (ierr) early_exit(8);

    ierr = STGetKSP(st,&ksp);if (ierr) early_exit(8);
    ierr = KSPSetType(ksp,KSPPREONLY);if (ierr) early_exit(8);
    ierr = KSPGetPC(ksp,&pc);if (ierr) early_exit(8);
    ierr = PCSetType(pc,PCCHOLESKY);if (ierr) early_exit(8);

    // ierr = EPSSetTrueResidual(eps, PETSC_TRUE); if (ierr) early_exit(9);
    // ierr = EPSSetConvergenceTest(eps, EPS_CONV_NORM); if (ierr) early_exit(9);

    ierr = EPSSolve(eps); if (ierr) early_exit(10);

    PetscInt lits;
    KSPGetTotalIterations(ksp,&lits);
    PetscPrintf(PETSC_COMM_WORLD," Number of Eigensolver Iterations: %D\n", lits);

    EPSGetConverged(eps, &nconv);
    if (nconv > 0) {
        ierr = EPSGetEigenpair(eps, 0, &eval_solution, &eigi, xr, xi);
        if (ierr) early_exit(11);
    } else {
        eval_solution = 0.0;
        PetscPrintf(PETSC_COMM_WORLD," Eigensolver Module Failed to Converge!");
    }

    EPSDestroy(&eps);
    return eval_solution;
}

void deposit_solution(std::vector<double> evals) {
    using namespace boost::interprocess;

    shared_memory_object values_smo (create_only, "best_eval_result", read_write);
    values_smo.truncate(evals.size() * sizeof(double));
    mapped_region values_reg(values_smo, read_write);

    std::memcpy(values_reg.get_address(), &evals[0], evals.size() * sizeof(double));
}

void deposit_eigenvector(Vec& xr, int size) {
    using namespace boost::interprocess;

    double* vec_data;
    VecGetArray(xr, &vec_data);

    shared_memory_object vec_smo;

    vec_smo = shared_memory_object(create_only, "best_evec_result", read_write);
    vec_smo.truncate(size * sizeof(double));

    mapped_region vec_reg(vec_smo, read_write, 0, size * sizeof(double));

    std::memcpy(vec_reg.get_address(), vec_data, vec_reg.get_size());

    VecRestoreArray(xr, &vec_data);
}
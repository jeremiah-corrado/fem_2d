#include "eigensolver/cpp_src/slepc_wrapper.h"
#include "eigensolver/src/slepc_wrapper.rs.h"

using namespace boost::interprocess;

// function declarations
void share_mat_data(rust::Vec<double> a, rust::Vec<int32_t> i, rust::Vec<int32_t> j, std::string mat_name);
int call_eigensolver(double target_eigenvalue);
double collect_solution(std::vector<double> &evec);

void clean_memory_channels();
std::array<std::string, 4> mem_channel_names(const std::string &name);

slepc_wrapper::EigenSolutionInternal slepc_wrapper::slepc_eigenproblem(
    double target_eigenvalue,
    AIJMatrix a_mat,
    AIJMatrix b_mat
) {
    clean_memory_channels();

    share_mat_data(a_mat.a, a_mat.i, a_mat.j, "A");
    share_mat_data(b_mat.a, b_mat.i, b_mat.j, "B");

    int status = call_eigensolver(target_eigenvalue);

    std::vector<double> eigenvector;
    double eigenvalue = collect_solution(eigenvector);

    return slepc_wrapper::EigenSolutionInternal {
        status,
        eigenvalue,
        std::make_unique<std::vector<double>>(eigenvector)
    };
}

//https://www.boost.org/doc/libs/1_55_0/doc/html/interprocess/sharedmemorybetweenprocesses.html
void share_mat_data(
    rust::Vec<double> a, 
    rust::Vec<int32_t> i, 
    rust::Vec<int32_t> j, 
    std::string mat_name
) {
    auto mem_names = mem_channel_names(mat_name);

    shared_memory_object values_smo(create_only, &mem_names[0][0], read_write);
    shared_memory_object rows_smo(create_only, &mem_names[1][0], read_write);
    shared_memory_object cols_smo(create_only, &mem_names[2][0], read_write);

    values_smo.truncate(a.size() * sizeof(double));
    rows_smo.truncate(i.size() * sizeof(int32_t));
    cols_smo.truncate(j.size() * sizeof(int32_t));

    mapped_region values_reg(values_smo, read_write);
    mapped_region rows_reg(rows_smo, read_write);
    mapped_region cols_reg(cols_smo, read_write);

    std::memcpy(values_reg.get_address(), &a[0], a.size() * sizeof(double));
    std::memcpy(rows_reg.get_address(), &i[0], i.size() * sizeof(int32_t));
    std::memcpy(cols_reg.get_address(), &j[0], j.size() * sizeof(int32_t));
}

int call_eigensolver(double target_eigenvalue) {
    char* eigensolver_path = std::getenv("EIGSOLVER_PATH");

    if (eigensolver_path) {
        std::stringstream command;
        command << "mpiexec -n 1 " << eigensolver_path << "/solve_gep -a " << target_eigenvalue;
        fflush(stdout);

        int status = system(&command.str()[0]);
        return status / 256;
    } else {
        std::cout << "Environment Variable 'EIGSOLVER_PATH' not set; cannot call Eigensolver!";
        return -1;
    };
}

double collect_solution(std::vector<double> &evec)
{
    // eigenvalue
    shared_memory_object eval_smo(open_only, "best_eval_result", read_only);
    mapped_region eval_reg(eval_smo, read_only);

    auto *best_eval = static_cast<double *>(eval_reg.get_address());
    double eval = best_eval[0];

    // eigenvector
    shared_memory_object evec_smo(open_only, "best_evec_result", read_only);
    mapped_region evec_reg(evec_smo, read_only);

    evec = std::vector<double>(evec_reg.get_size() / sizeof(double));
    std::memcpy(&evec[0], evec_reg.get_address(), evec_reg.get_size());

    return eval;
}

void clean_memory_channels()
{
    auto mem_names_a = mem_channel_names("A");
    auto mem_names_b = mem_channel_names("B");

    for (const auto &n : mem_names_a)
    {
        shared_memory_object::remove(&n[0]);
    }

    for (const auto &n : mem_names_b)
    {
        shared_memory_object::remove(&n[0]);
    }

    shared_memory_object::remove("best_eval_result");
    shared_memory_object::remove("best_evec_result");
}

std::array<std::string, 4> mem_channel_names(const std::string &name)
{
    return {
        name + "_mat_vals",
        name + "_mat_rows",
        name + "_mat_cols",
    };
}
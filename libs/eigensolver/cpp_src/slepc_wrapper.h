#pragma once

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "rust/cxx.h"
#include <memory>
#include <array>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>

namespace slepc_wrapper {
    struct AIJMatrix;
    struct EigenSolutionInternal;

    EigenSolutionInternal slepc_eigenproblem(
        double target_eigenvalue,
        AIJMatrix a_mat,
        AIJMatrix b_mat
    );
}
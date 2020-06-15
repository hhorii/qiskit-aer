/**
 * This code is part of Qiskit.
 *
 * (C) Copyright IBM 2018, 2019.
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */



#ifndef _qv_qvintrin_hpp_
#define _qv_qvintrin_hpp_

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>

namespace QV {

// Type aliases
using uint_t = uint64_t;
using int_t = int64_t;
using reg_t = std::vector<uint_t>;
using indexes_t = std::unique_ptr<uint_t[]>;
template <size_t N> using areg_t = std::array<uint_t, N>;
template <typename T> using cvector_t = std::vector<std::complex<T>>;

bool is_intrinsics();

bool apply_matrix_opt(std::complex<float>* qv_data,
    uint_t data_size,
    const reg_t& qregs,
    const cvector_t<float>& mat,
    uint_t omp_threads);

bool apply_matrix_opt(std::complex<double>* qv_data,
    uint_t data_size,
    const reg_t& qregs,
    const cvector_t<double>& mat,
    uint_t omp_threads);

}
//------------------------------------------------------------------------------
#endif // end module

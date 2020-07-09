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

#if defined(__PPC__) | defined(__PPC64__)
#elif defined(__arm__) || defined(__arm64__)
#else
#ifdef _QV_HEDDER_ONLY
#include <cpuid.h>
// only if x86
#include "qvintrin_avx.hpp"
#endif
#endif

namespace QV {

// Type aliases
using uint_t = uint64_t;
using int_t = int64_t;
using reg_t = std::vector<uint_t>;
using indexes_t = std::unique_ptr<uint_t[]>;
template <size_t N> using areg_t = std::array<uint_t, N>;
template <typename T> using cvector_t = std::vector<std::complex<T>>;

#ifdef _QV_HEDDER_ONLY
bool is_intrinsics () {
#if defined(__PPC__) | defined(__PPC64__)
  return false;
#elif defined(__arm__) || defined(__arm64__)
  return false;
#else
  return is_avx2_supported();
#endif
}

bool apply_matrix_opt(
    float* qv_data,
    const uint64_t data_size,
    const uint64_t* qregs,
    const uint64_t qregs_size,
    const float* fmat,
    const uint_t omp_threads) {
#if defined(__PPC__) | defined(__PPC64__)
  return false;
#elif defined(__arm__) || defined(__arm64__)
  return false;
#else
  return apply_matrix_avx <float> (
      (void *) qv_data, data_size, qregs, qregs_size, fmat, omp_threads);
#endif
}

bool apply_matrix_opt(
    double* qv_data,
    const uint64_t data_size,
    const uint64_t* qregs,
    const uint64_t qregs_size,
    const double* dmat,
    const uint_t omp_threads) {

#if defined(__PPC__) | defined(__PPC64__)
  return false;
#elif defined(__arm__) || defined(__arm64__)
  return false;
#else
  return apply_matrix_avx <double> (
      qv_data, data_size, qregs, qregs_size, dmat, omp_threads);
#endif
}
#else
bool is_intrinsics();

bool apply_matrix_opt(
    float* qv_data,
    const uint64_t data_size,
    const uint64_t* qregs,
    const uint64_t qregs_size,
    float* mat,
    uint_t omp_threads);

bool apply_matrix_opt(
    double* qv_data,
    const uint64_t data_size,
    const uint64_t* qregs,
    const uint64_t qregs_size,
    const double* mat,
    const uint_t omp_threads);
#endif
}
//------------------------------------------------------------------------------
#endif // end module

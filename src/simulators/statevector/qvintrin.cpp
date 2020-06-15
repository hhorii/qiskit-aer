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

#ifdef _MSC_VER
#include <intrin.h>
#elif defined(__GNUC__)
#include <cpuid.h>
#endif

#include "qvintrin.hpp"
#include "qubitvector.hpp"
#include "qvintrin_avx.hpp"

using namespace QV;

inline bool is_intrinsics () {
  return is_avx2_supported();
}

inline bool apply_matrix_opt(std::complex<float>* qv_data,
    uint_t data_size,
    const reg_t& qregs,
    const cvector_t<float>& mat,
    uint_t omp_threads) {

  return apply_matrix_avx <float> (
      qv_data, data_size, qregs, mat, omp_threads) != Avx::NotApplied;
}

inline bool apply_matrix_opt(std::complex<double>* qv_data,
    uint_t data_size,
    const reg_t& qregs,
    const cvector_t<double>& mat,
    uint_t omp_threads) {

  return apply_matrix_avx <double> (
      qv_data, data_size, qregs, mat, omp_threads) != Avx::NotApplied;
}



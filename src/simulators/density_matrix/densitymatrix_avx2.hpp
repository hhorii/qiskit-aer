/**
 * This code is part of Qiskit.
 *
 * (C) Copyright IBM 2020.
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */

#ifndef _qv_density_matrix_avx2_hpp_
#define _qv_density_matrix_avx2_hpp_

#include "densitymatrix.hpp"
#include "simulators/statevector/qv_avx2.hpp"
#include "misc/common_macros.hpp"

namespace AER {
namespace QV {

// Type aliases
using uint_t = uint64_t;
using int_t = int64_t;
using reg_t = std::vector<uint_t>;
using indexes_t = std::unique_ptr<uint_t[]>;
template <size_t N> using areg_t = std::array<uint_t, N>;
template <typename T> using cvector_t = std::vector<std::complex<T>>;

//============================================================================
// DensityMatrixAvx2 class
//============================================================================

template <typename data_t = double>
class DensityMatrixAvx2 : public DensityMatrix<data_t, DensityMatrixAvx2<data_t>> {

  // We need this to access the base class members
  using Base = DensityMatrix<data_t, DensityMatrixAvx2<data_t>>;

public:

  //-----------------------------------------------------------------------
  // Constructors and Destructor
  //-----------------------------------------------------------------------

  DensityMatrixAvx2();
  explicit DensityMatrixAvx2(size_t num_qubits);
  DensityMatrixAvx2(const DensityMatrixAvx2& obj) = delete;
  DensityMatrixAvx2 &operator=(const DensityMatrixAvx2& obj) = delete;

  // Apply a 1-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized 1-qubit matrix.
  void apply_matrix_to_base_vector(const uint_t qubit, const cvector_t<double> &mat);

  // Apply a N-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized N-qubit matrix.
  void apply_matrix_to_base_vector(const reg_t &qubits, const cvector_t<double> &mat);

  protected:
  size_t calculate_num_threads();
};

// We do not define this functions in case we don't use AVX2
// so it can compile, as this class won't be used
#if defined(_MSC_VER) || defined(GNUC_AVX2)
// ostream overload for templated qubitvector
template <typename data_t>
inline std::ostream &operator<<(std::ostream &out, const QV::DensityMatrixAvx2<data_t>&qv) {

  out << "[";
  size_t last = qv.size() - 1;
  for (size_t i = 0; i < qv.size(); ++i) {
    out << qv[i];
    if (i != last)
      out << ", ";
  }
  out << "]";
  return out;
}


//------------------------------------------------------------------------------
// Constructors & Destructor
//------------------------------------------------------------------------------
template <typename data_t>
DensityMatrixAvx2<data_t>::DensityMatrixAvx2(size_t num_qubits) {
  Base::num_qubits_ = 0;
  Base::data_ = nullptr;
  Base::checkpoint_ = 0;
  Base::set_num_qubits(num_qubits);
}

template <typename data_t>
DensityMatrixAvx2<data_t>::DensityMatrixAvx2() : DensityMatrixAvx2(0) {}

template <typename data_t>
void DensityMatrixAvx2<data_t>::apply_matrix_to_base_vector(const uint_t qubit,
                                           const cvector_t<double>& mat) {
  if ((mat[1] == 0.0 && mat[2] == 0.0) || (mat[0] == 0.0 && mat[3] == 0.0)) {
    Base::apply_matrix(qubit, mat);
    return;
  }

  reg_t qubits = {qubit};
  if (apply_matrix_avx<data_t>(reinterpret_cast<data_t*>(Base::data_),
                               Base::data_size_, qubits.data(), qubits.size(),
                               reinterpret_cast<data_t *>(Base::convert(mat).data()),
                               calculate_num_threads()) == Avx::NotApplied) {
    Base::apply_matrix(qubit, mat);
  }
}

template <typename data_t>
void DensityMatrixAvx2<data_t>::apply_matrix_to_base_vector(const reg_t& qubits,
                                           const cvector_t<double>& mat) {
  if (apply_matrix_avx<data_t>(reinterpret_cast<data_t*>(Base::data_),
                               Base::data_size_, qubits.data(), qubits.size(),
                               reinterpret_cast<data_t *>(Base::convert(mat).data()),
                               calculate_num_threads()) == Avx::NotApplied) {
    Base::apply_matrix(qubits, mat);
  }
}

template <typename data_t>
size_t DensityMatrixAvx2<data_t>::calculate_num_threads() {
  if (Base::num_qubits_ > Base::omp_threshold_ && Base::omp_threads_ > 1) {
    return Base::omp_threads_;
  }
  return 1;
}
#endif

template class AER::QV::DensityMatrixAvx2<double>;
template class AER::QV::DensityMatrixAvx2<float>;

}
}
//------------------------------------------------------------------------------
#endif // end module

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


/*
 * This is a simplified version of qubitvector.hpp
 * https://github.com/Qiskit/qiskit-aer/blob/0.5.2/src/simulators/statevector/qubitvector_thrust.hpp
 */


#ifndef _statevector_thrust_hpp
#define _statevector_thrust_hpp

#include <thrust/for_each.h>
#include <thrust/complex.h>
#include <thrust/functional.h>
#include <thrust/tuple.h>

#include <thrust/host_vector.h>

#include <thrust/system/omp/execution_policy.h>


#include <complex>
#include <bitset>
#include <random>

using namespace std;

template <class thrust_vector = thrust::host_vector<thrust::complex<double>>>
class QubitVectorThrust {

public:

  QubitVectorThrust(size_t num_qubits_): num_qubits(num_qubits_) {
    psi.resize(1 << num_qubits);

    zero();
    thrust::complex<double> t;
    t = 1.0;
    psi[0] = t;
  }

  virtual ~QubitVectorThrust() {
    psi.resize(0);
  }

  class _fill {
    public:
    _fill(thrust_vector& _psi, thrust::complex<double> _v) {
      this->psi = thrust::raw_pointer_cast(_psi.data());
      this->v = _v;
    }

    __host__ __device__ void operator()(unsigned i) const {
      psi[i] = v;
    }

    protected:
    thrust::complex<double>* psi;
    thrust::complex<double> v;
  };

  void zero() {
    auto itr = thrust::counting_iterator<unsigned>(0);
    thrust::for_each(itr, itr + (1U << num_qubits), _fill(psi, 0));
  }

  class _unitary_1 {
    public:
    _unitary_1(thrust_vector& _psi, unsigned n,
          thrust::complex<double> _u00, thrust::complex<double> _u01,
          thrust::complex<double> _u10, thrust::complex<double> _u11) {
      this->psi = thrust::raw_pointer_cast(_psi.data());
      this->n = n;
      this->u00 = _u00;
      this->u01 = _u01;
      this->u10 = _u10;
      this->u11 = _u11;
    }

    __host__ __device__ void operator()(unsigned i) const {
      unsigned idx0 = i;
      unsigned low = idx0 & ((1 << n) - 1);
      idx0 >>= n;
      idx0 <<= (n + 1);
      idx0 |= low;

      unsigned idx1 = idx0 | (1 << n);

      thrust::complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1];
      thrust::complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1];

      psi[idx0] = p0;
      psi[idx1] = p1;
    }

    protected:
    thrust::complex<double>* psi;
    unsigned n;
    thrust::complex<double> u00;
    thrust::complex<double> u01;
    thrust::complex<double> u10;
    thrust::complex<double> u11;
  };

  void unitary(unsigned n, //
      complex<double> u00, complex<double> u01, //
      complex<double> u10, complex<double> u11) {

    auto itr = thrust::counting_iterator<unsigned>(0);
    thrust::for_each(itr, itr + (1U << (num_qubits - 1)), _unitary_1(psi, n, u00, u01, u10, u11));
  }

  class _unitary_2 {
    public:
    _unitary_2(thrust_vector& _psi, unsigned n, unsigned m,
          thrust::complex<double> _u00, thrust::complex<double> _u01, thrust::complex<double> _u02, thrust::complex<double> _u03,
          thrust::complex<double> _u10, thrust::complex<double> _u11, thrust::complex<double> _u12, thrust::complex<double> _u13,
          thrust::complex<double> _u20, thrust::complex<double> _u21, thrust::complex<double> _u22, thrust::complex<double> _u23,
          thrust::complex<double> _u30, thrust::complex<double> _u31, thrust::complex<double> _u32, thrust::complex<double> _u33) {
      this->psi = thrust::raw_pointer_cast(_psi.data());
      this->n = n;
      this->m = m;
      this->u00 = _u00; this->u01 = _u01; this->u02 = _u02; this->u03 = _u03;
      this->u10 = _u10; this->u11 = _u11; this->u12 = _u12; this->u13 = _u13;
      this->u20 = _u20; this->u21 = _u21; this->u22 = _u22; this->u23 = _u23;
      this->u30 = _u30; this->u31 = _u31; this->u32 = _u32; this->u33 = _u33;
    }

    __host__ __device__ void operator()(unsigned i) const {
      if (m < n) {
        unsigned idx0 = i;
        unsigned low = idx0 & ((1 << m) - 1);
        idx0 >>= m;
        idx0 <<= (m + 1);
        idx0 |= low;

        unsigned middle = idx0 & ((1 << n) - 1);
        idx0 >>= n;
        idx0 <<= (n + 1);
        idx0 |= middle;

        unsigned idx1 = idx0 | (1 << m);
        unsigned idx2 = idx0 | (1 << n);
        unsigned idx3 = idx0 | (1 << n) | (1 << m);
        thrust::complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1] + u02 * psi[idx2] + u03 * psi[idx3];
        thrust::complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1] + u12 * psi[idx2] + u13 * psi[idx3];
        thrust::complex<double> p2 = u20 * psi[idx0] + u21 * psi[idx1] + u22 * psi[idx2] + u23 * psi[idx3];
        thrust::complex<double> p3 = u30 * psi[idx0] + u31 * psi[idx1] + u32 * psi[idx2] + u33 * psi[idx3];
        psi[idx0] = p0;
        psi[idx1] = p1;
        psi[idx2] = p2;
        psi[idx3] = p3;
      } else {
        unsigned idx0 = i;
        unsigned low = idx0 & ((1 << n) - 1);
        idx0 >>= n;
        idx0 <<= (n + 1);
        idx0 |= low;

        unsigned middle = idx0 & ((1 << m) - 1);
        idx0 >>= m;
        idx0 <<= (m + 1);
        idx0 |= middle;

        unsigned idx1 = idx0 | (1 << n);
        unsigned idx2 = idx0 | (1 << m);
        unsigned idx3 = idx0 | (1 << n) | (1 << m);
        thrust::complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1] + u02 * psi[idx2] + u03 * psi[idx3];
        thrust::complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1] + u12 * psi[idx2] + u13 * psi[idx3];
        thrust::complex<double> p2 = u20 * psi[idx0] + u21 * psi[idx1] + u22 * psi[idx2] + u23 * psi[idx3];
        thrust::complex<double> p3 = u30 * psi[idx0] + u31 * psi[idx1] + u32 * psi[idx2] + u33 * psi[idx3];
        psi[idx0] = p0;
        psi[idx1] = p1;
        psi[idx2] = p2;
        psi[idx3] = p3;
      }
    }

    protected:
    thrust::complex<double>* psi;
    unsigned n;
    unsigned m;
    thrust::complex<double> u00; thrust::complex<double> u01; thrust::complex<double> u02; thrust::complex<double> u03;
    thrust::complex<double> u10; thrust::complex<double> u11; thrust::complex<double> u12; thrust::complex<double> u13;
    thrust::complex<double> u20; thrust::complex<double> u21; thrust::complex<double> u22; thrust::complex<double> u23;
    thrust::complex<double> u30; thrust::complex<double> u31; thrust::complex<double> u32; thrust::complex<double> u33;
  };

  void unitary(unsigned n, unsigned m, //
      complex<double> u00, complex<double> u01, complex<double> u02, complex<double> u03, //
      complex<double> u10, complex<double> u11, complex<double> u12, complex<double> u13, //
      complex<double> u20, complex<double> u21, complex<double> u22, complex<double> u23, //
      complex<double> u30, complex<double> u31, complex<double> u32, complex<double> u33) {

    auto itr = thrust::counting_iterator<unsigned>(0);
    thrust::for_each(itr, itr + (1U << (num_qubits - 2)), _unitary_2(psi, n, m, u00, u01, u02, u03,
                                                                                 u10, u11, u12, u13,
                                                                                 u20, u21, u22, u23,
                                                                                 u30, u31, u32, u33));
  }

  void u3(unsigned n, double theta, double phi, double lambda) {
    complex<double> i(0., 1.);

    complex<double> u00 = cos(theta / 2.);
    complex<double> u01 = -exp(i * lambda) * sin(theta / 2.);
    complex<double> u10 = exp(i * phi) * sin(theta / 2.);
    complex<double> u11 = exp(i * (phi + lambda)) * cos(theta / 2.);

    unitary(n, u00, u01, u10, u11);
  }

  void u2(unsigned n, double phi, double lambda) {
    u3(n, M_PI/2., phi, lambda);
  }

  void u1(unsigned n, double lambda) {
    u3(n, 0., 0., lambda);
  }

  void x(unsigned n) {
    u3(n, M_PI, 0, M_PI);
  }

  void y(unsigned n) {
    u3(n, M_PI, M_PI/2., M_PI/2.);
  }

  void z(unsigned n) {
    u1(n, M_PI);
  }

  void h(unsigned n) {
    u2(n, 0, M_PI);
  }

  void s(unsigned n) {
    u1(n, M_PI/2.);
  }

  void sdg(unsigned n) {
    u1(n, -M_PI/2.);
  }

  void t(unsigned n) {
    u1(n, M_PI/4.);
  }

  void tdg(unsigned n) {
    u1(n, -M_PI/4.);
  }

  void rx(unsigned n, double theta) {
    u3(n, theta, -M_PI/2, M_PI/2);
  }

  void ry(unsigned n, double theta) {
    u3(n, theta, 0., 0.);
  }

  void rz(unsigned n, double phi) {
    u1(n, phi);
  }

  void cnot(unsigned control, unsigned target) {
    unitary(target, control,
        1., 0., 0., 0.,
        0., 1., 0., 0,
        0., 0., 0., 1.,
        0., 0., 1., 0.);
  }

  void dump() {
    for (auto i = 0; i < (1UL << num_qubits); ++i) {
      if (i != 0)
        cout << ", ";
      cout << psi[i];
    }
    cout << endl;
  }

  class _measure {
    public:
    _measure(thrust_vector& _psi, unsigned n) {
      this->psi = thrust::raw_pointer_cast(_psi.data());
      this->n = n;
    }

    __host__ __device__ double operator()(unsigned i) const {
      unsigned idx0 = i;
      unsigned low = idx0 & ((1 << n) - 1);
      idx0 >>= n;
      idx0 <<= (n + 1);
      idx0 |= low;

      return thrust::norm(psi[idx0]);
    }

    protected:
    thrust::complex<double>* psi;
    unsigned n;
  };

  int measure(unsigned n) {
    auto itr = thrust::counting_iterator<unsigned>(0);
    double p = thrust::transform_reduce(itr, itr + (1U << (num_qubits - 1)), _measure(psi, n), 0.0, thrust::plus<double>());

    if (this->random() <= p) {
      unitary(n, std::sqrt(p), 0, 0, 0);
      return 0;
    } else {
      unitary(n, 0, 0, 0, std::sqrt(1. - p));
      return 1;
    }
  }

protected:

  double random() {
    return std::uniform_real_distribution<double>(0., 1.)(rng);
  }

  size_t num_qubits;
  thrust_vector psi;
  mt19937_64 rng;

};

#endif


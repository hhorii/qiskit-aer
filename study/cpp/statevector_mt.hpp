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
 * https://github.com/Qiskit/qiskit-aer/blob/0.5.2/src/simulators/statevector/qubitvector.hpp
 */


#ifndef _statevector_mt_hpp
#define _statevector_mt_hpp

#include <complex>
#include <bitset>
#include <random>

using namespace std;

class QubitVectorMT {

public:

  QubitVectorMT(size_t num_qubits_): num_qubits(num_qubits_) {
    psi = (complex<double>*) malloc(sizeof(complex<double>) * (1 << num_qubits));
#pragma omp parallel for
    for (unsigned i = 0; i < (1 << num_qubits); ++i)
      psi[i] = .0;
    psi[0] = 1.;
  }

  virtual ~QubitVectorMT() {
    free(psi);
  }

  void unitary(unsigned n, //
      complex<double> u00, complex<double> u01, //
      complex<double> u10, complex<double> u11) {

#pragma omp parallel for
    for (unsigned i = 0; i < (1 << (num_qubits - 1)); ++i) {
        unsigned idx0 = i;
        unsigned low = idx0 & ((1 << n) - 1);
        idx0 >>= n;
        idx0 <<= (n + 1);
        idx0 |= low;

        unsigned idx1 = idx0 | (1 << n);

        complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1];
        complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1];

        psi[idx0] = p0;
        psi[idx1] = p1;
    }
  }

  void unitary(unsigned n, unsigned m, //
      complex<double> u00, complex<double> u01, complex<double> u02, complex<double> u03, //
      complex<double> u10, complex<double> u11, complex<double> u12, complex<double> u13, //
      complex<double> u20, complex<double> u21, complex<double> u22, complex<double> u23, //
      complex<double> u30, complex<double> u31, complex<double> u32, complex<double> u33) {

    if (m < n) {
#pragma omp parallel for
      for (unsigned i = 0; i < (1 << (num_qubits - 2)); ++i) {
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
        complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1] + u02 * psi[idx2] + u03 * psi[idx3];
        complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1] + u12 * psi[idx2] + u13 * psi[idx3];
        complex<double> p2 = u20 * psi[idx0] + u21 * psi[idx1] + u22 * psi[idx2] + u23 * psi[idx3];
        complex<double> p3 = u30 * psi[idx0] + u31 * psi[idx1] + u32 * psi[idx2] + u33 * psi[idx3];
        psi[idx0] = p0;
        psi[idx1] = p1;
        psi[idx2] = p2;
        psi[idx3] = p3;
      }
    } else {
#pragma omp parallel for
      for (unsigned i = 0; i < (1 << (num_qubits - 2)); ++i) {
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
        complex<double> p0 = u00 * psi[idx0] + u01 * psi[idx1] + u02 * psi[idx2] + u03 * psi[idx3];
        complex<double> p1 = u10 * psi[idx0] + u11 * psi[idx1] + u12 * psi[idx2] + u13 * psi[idx3];
        complex<double> p2 = u20 * psi[idx0] + u21 * psi[idx1] + u22 * psi[idx2] + u23 * psi[idx3];
        complex<double> p3 = u30 * psi[idx0] + u31 * psi[idx1] + u32 * psi[idx2] + u33 * psi[idx3];
        psi[idx0] = p0;
        psi[idx1] = p1;
        psi[idx2] = p2;
        psi[idx3] = p3;
      }
    }
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

  int measure(unsigned n) {

    double p = .0;
#pragma omp parallel for reduction(+:p)
    for (unsigned i = 0; i < (1 << (num_qubits - 1)); ++i) {
      unsigned idx0 = i;
      unsigned low = idx0 & ((1 << n) - 1);
      idx0 >>= n;
      idx0 <<= (n + 1);
      idx0 |= low;
      p += norm(psi[idx0]);
    }

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

  complex<double>* psi;
  mt19937_64 rng;
  size_t num_qubits;

};

#endif


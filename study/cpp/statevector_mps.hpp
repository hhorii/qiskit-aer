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


#ifndef _statevector_mps_hpp
#define _statevector_mps_hpp

#include <complex>
#include <bitset>
#include <random>

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <sstream>

#include "itensor/all.h"

using namespace std;
using namespace itensor;

class QubitVectorMPS {

public:

  QubitVectorMPS(size_t num_qubits_, bool debug_ = false): num_qubits(num_qubits_), debug(debug_) {
    Index head(1);
    Index prev_rbond = head;
    for (int i = 0; i < num_qubits - 1; ++i) {
      Index qbit(2);
      Index lbond(1);
      itensor::ITensor legMat(qbit, prev_rbond, lbond);
      legMat.set(qbit(1), prev_rbond(1), lbond(1), 1.);
      legMats.push_back(legMat);

      Index rbond(1);
      itensor::ITensor bondMat(lbond, rbond);
      bondMat.set(lbond(1), rbond(1), 1.);
      bondMats.push_back(bondMat);
      prev_rbond = rbond;
    }
    Index qbit(2);
    Index tail(1);
    itensor::ITensor legMat(qbit, prev_rbond, tail);
    legMat.set(qbit(1), prev_rbond(1), tail(1), 1.);
    legMats.push_back(legMat);
  }

  virtual ~QubitVectorMPS() {
  }

  void unitary(int n, //
      complex<double> u00, complex<double> u01, //
      complex<double> u10, complex<double> u11) {

    if (debug) {
      cout << "DEBUG: > unitary(" << n << ", [" << u00 << ", " << u01 << ", " << u10 << ", " << u11 << "]" << endl;
      cout << "DEBUG:   legMats[" << n << "]=" << to_vector_string(legMats[n]) << endl;
    }

    auto ind_in = legMats[n].inds()[0];
    auto ind_out = itensor::Index(2);
    auto gate = itensor::ITensor(ind_in, ind_out);
    if (u00 != 0.)
      gate.set(ind_in(1), ind_out(1), u00);
    if (u01 != 0.)
      gate.set(ind_in(1), ind_out(2), u01);
    if (u10 != 0.)
      gate.set(ind_in(2), ind_out(1), u10);
    if (u11 != 0.)
      gate.set(ind_in(2), ind_out(2), u11);

    legMats[n] = gate * legMats[n];

    if (debug) {
      cout << "DEBUG: < unitary(" << n << ", [" << u00 << ", " << u01 << ", " << u10 << ", " << u11 << "]" << endl;
      cout << "DEBUG:   legMats[" << n << "]=" << to_vector_string(legMats[n]) << endl;
    }
  }

  void unitary(int n, int m, //
      complex<double> u00, complex<double> u01, complex<double> u02, complex<double> u03, //
      complex<double> u10, complex<double> u11, complex<double> u12, complex<double> u13, //
      complex<double> u20, complex<double> u21, complex<double> u22, complex<double> u23, //
      complex<double> u30, complex<double> u31, complex<double> u32, complex<double> u33) {

    if (debug) {
      cout << "DEBUG: > unitary(" << n << ", ["
          << u00 << ", " << u01 << ", " << u02 << ", " << u03 << ", "
          << u10 << ", " << u11 << ", " << u12 << ", " << u13 << ", "
          << u20 << ", " << u21 << ", " << u22 << ", " << u23 << ", "
          << u30 << ", " << u31 << ", " << u32 << ", " << u33 << "]" << endl;
      int from = min(n, m);
      int to = max(n, m);
      for (int i = from; i < to; ++i) {
        cout << "DEBUG:   legMats[" << i << "]=" << to_vector_string(legMats[i]) << endl;
        cout << "DEBUG:   bondMats[" << i << "]=" << to_vector_string(bondMats[i]) << endl;
      }
      cout << "DEBUG:   legMats[" << to << "]=" << to_vector_string(legMats[to]) << endl;
    }

    int iqbit_in0, iqbit_in1;
    if (n < m - 1) {
      permute_to(n, m - 1);
      iqbit_in0 = m - 1;
      iqbit_in1 = m;
    } else if (m < n - 1) {
      permute_to(m, n - 1);
      iqbit_in0 = n;
      iqbit_in1 = n - 1;
    } else {
      iqbit_in0 = n;
      iqbit_in1 = m;
    }

    auto ind_in0 = legMats[iqbit_in0].inds()[0];// control
    auto ind_in1 = legMats[iqbit_in1].inds()[0];
    auto ind_out0 = itensor::Index(2);
    auto ind_out1 = itensor::Index(2);
    Index ind_lower;
    if (iqbit_in0 < iqbit_in1) {
      ind_lower = ind_out0;
    } else {
      ind_lower = ind_out1;
    }

    auto g = itensor::ITensor(ind_in0, ind_in1, ind_out0, ind_out1);
    if (u00 != 0.)
      g.set(ind_out0(1), ind_out1(1), ind_in0(1), ind_in1(1), u00);
    if (u01 != 0.)
      g.set(ind_out0(1), ind_out1(1), ind_in0(1), ind_in1(2), u01);
    if (u02 != 0.)
      g.set(ind_out0(1), ind_out1(1), ind_in0(2), ind_in1(1), u02);
    if (u03 != 0.)
      g.set(ind_out0(1), ind_out1(1), ind_in0(2), ind_in1(2), u03);
    if (u10 != 0.)
      g.set(ind_out0(1), ind_out1(2), ind_in0(1), ind_in1(1), u10);
    if (u11 != 0.)
      g.set(ind_out0(1), ind_out1(2), ind_in0(1), ind_in1(2), u11);
    if (u12 != 0.)
      g.set(ind_out0(1), ind_out1(2), ind_in0(2), ind_in1(1), u12);
    if (u13 != 0.)
      g.set(ind_out0(1), ind_out1(2), ind_in0(2), ind_in1(2), u13);
    if (u20 != 0.)
      g.set(ind_out0(2), ind_out1(1), ind_in0(1), ind_in1(1), u20);
    if (u21 != 0.)
      g.set(ind_out0(2), ind_out1(1), ind_in0(1), ind_in1(2), u21);
    if (u22 != 0.)
      g.set(ind_out0(2), ind_out1(1), ind_in0(2), ind_in1(1), u22);
    if (u23 != 0.)
      g.set(ind_out0(2), ind_out1(1), ind_in0(2), ind_in1(2), u23);
    if (u30 != 0.)
      g.set(ind_out0(2), ind_out1(2), ind_in0(1), ind_in1(1), u30);
    if (u31 != 0.)
      g.set(ind_out0(2), ind_out1(2), ind_in0(1), ind_in1(2), u31);
    if (u32 != 0.)
      g.set(ind_out0(2), ind_out1(2), ind_in0(2), ind_in1(1), u32);
    if (u33 != 0.)
      g.set(ind_out0(2), ind_out1(2), ind_in0(2), ind_in1(2), u33);

    int min_iqbit = min(iqbit_in0, iqbit_in1);
    int max_iqbit = max(iqbit_in0, iqbit_in1);

    auto tobe_svd =
        g * legMats[iqbit_in0] * bondMats[min_iqbit] * legMats[iqbit_in1];

    ITensor legMat(legMats[min_iqbit].inds()[1], ind_lower), bondMat, restTensor;
    itensor::svd(tobe_svd, legMat, bondMat, restTensor, {"Cutoff", svdCutoff});

    legMats[min_iqbit] = legMat;
    bondMats[min_iqbit] = bondMat;

    kickback_ind(restTensor, restTensor.inds()[1]);

    legMats[max_iqbit] = restTensor;
    if (n < m - 1) {
      permute_to(m - 1, n);
    } else if (m < n - 1) {
      permute_to(n - 1, m);
    }

    if (debug) {
      cout << "DEBUG: < unitary(" << n << ", ["
          << u00 << ", " << u01 << ", " << u02 << ", " << u03 << ", "
          << u10 << ", " << u11 << ", " << u12 << ", " << u13 << ", "
          << u20 << ", " << u21 << ", " << u22 << ", " << u23 << ", "
          << u30 << ", " << u31 << ", " << u32 << ", " << u33 << "]" << endl;
      int from = min(n, m);
      int to = max(n, m);
      for (int i = from; i < to; ++i) {
        cout << "DEBUG:   legMats[" << i << "]=" << to_vector_string(legMats[i]) << endl;
        cout << "DEBUG:   bondMats[" << i << "]=" << to_vector_string(bondMats[i]) << endl;
      }
      cout << "DEBUG:   legMats[" << to << "]=" << to_vector_string(legMats[to]) << endl;
    }
  }

  void u3(int n, double theta, double phi, double lambda) {
    complex<double> i(0., 1.);

    complex<double> u00 = cos(theta / 2.);
    complex<double> u01 = -exp(i * lambda) * sin(theta / 2.);
    complex<double> u10 = exp(i * phi) * sin(theta / 2.);
    complex<double> u11 = exp(i * (phi + lambda)) * cos(theta / 2.);

    unitary(n, u00, u01, u10, u11);
  }

  void u2(int n, double phi, double lambda) {
    u3(n, M_PI/2., phi, lambda);
  }

  void u1(int n, double lambda) {
    u3(n, 0., 0., lambda);
  }

  void x(int n) {
    u3(n, M_PI, 0, M_PI);
  }

  void y(int n) {
    u3(n, M_PI, M_PI/2., M_PI/2.);
  }

  void z(int n) {
    u1(n, M_PI);
  }

  void h(int n) {
    u2(n, 0, M_PI);
  }

  void s(int n) {
    u1(n, M_PI/2.);
  }

  void sdg(int n) {
    u1(n, -M_PI/2.);
  }

  void t(int n) {
    u1(n, M_PI/4.);
  }

  void tdg(int n) {
    u1(n, -M_PI/4.);
  }

  void rx(int n, double theta) {
    u3(n, theta, -M_PI/2, M_PI/2);
  }

  void ry(int n, double theta) {
    u3(n, theta, 0., 0.);
  }

  void rz(int n, double phi) {
    u1(n, phi);
  }

  void cnot(int control, int target) {
    unitary(control, target,
        1., 0., 0., 0.,
        0., 1., 0., 0,
        0., 0., 0., 1.,
        0., 0., 1., 0.);
  }

  void swap(int control, int target) {
    unitary(target, control,
        1., 0., 0., 0.,
        0., 0., 1., 0,
        0., 1., 0., 0.,
        0., 0., 0., 1.);
  }

  vector<std::complex<double>> to_vector(itensor::ITensor& tensor) {
    std::vector<std::complex<double>> wf;
    auto store_wf = [&](itensor::Cplx c) {
      auto real = c.real();
      auto imag = c.imag();

      if (std::fabs(real) < 1e-12)
        real = 0.0;
      if (std::fabs(imag) < 1e-12)
        imag = 0.0;

      wf.push_back(std::complex<double>(real, imag));
    };
    tensor.visit(store_wf);
    return wf;
  }

  string to_vector_string(itensor::ITensor& tensor) {
    stringstream ss;
    ss << "[";
    auto print_wf = [&](itensor::Cplx c) {
      auto real = c.real();
      auto imag = c.imag();

      if (std::fabs(real) < 1e-12)
        real = 0.0;
      if (std::fabs(imag) < 1e-12)
        imag = 0.0;

      ss << "[" << real << ", " << imag << "], ";
    };
    tensor.visit(print_wf);
    ss << "]";
    return ss.str();
  }

  std::vector<std::complex<double>> statevector() {
    auto mps = legMats[0];
    for (int i = 1; i < num_qubits; ++i) {
      mps *= bondMats[i - 1];
      mps *= legMats[i];
    }
    auto normed_wf = mps / itensor::norm(mps);
    auto wf = to_vector(normed_wf);

    return wf;
  }

  void dump() {
    auto psi = statevector();
    for (auto i = 0; i < (1UL << num_qubits); ++i) {
      if (i != 0)
        cout << ", ";
      cout << psi[i];
    }
    cout << endl;
  }

  void dump_tensor() {
    auto psi = statevector();
    for (int i = 0; i < num_qubits - 1; ++i) {
      cout << "*****legMats[" << i << "]*****" << endl;
      cout << legMats[i] << to_vector_string(legMats[i]) << endl;
      cout << endl;
      cout << "*****bondMats[" << i << "]*****" << endl;
      cout << bondMats[i] << to_vector_string(bondMats[i]) << endl;
      cout << endl;
    }
    cout << "*****legMats[" << (num_qubits - 1) << "]*****" << endl;
    cout << legMats[(num_qubits - 1)] << to_vector_string(legMats[(num_qubits - 1)]) << endl;
    cout << endl;
  }

  int measure(int n) {

     double p = .0;

     auto ind_measured = legMats[n].inds()[0];
     auto ind_measured_p = ind_measured;
     ind_measured_p.prime();

     auto tMeasure0 = itensor::ITensor(ind_measured, ind_measured_p);
     tMeasure0.set(ind_measured_p(1), ind_measured(1), 1.);
     double p0 = average(n, tMeasure0) / wavefunc_inner();

     double rv = (std::rand() % 1000000) / 1000000.;

     if (rv < p0) {
       legMats[n] =
           tMeasure0 * legMats[n]; // collapse wavefunction
       legMats[n].prime(-1, ind_measured_p);
       return 0;
     } else {
       auto tMeasure1 = itensor::ITensor(ind_measured, ind_measured_p);
       tMeasure1.set(ind_measured_p(2), ind_measured(2), 1.);
       legMats[n] =
           tMeasure1 * legMats[n]; // collapse wavefunction
       legMats[n].prime(-1, ind_measured_p);
       return 1;
     }
  }

//protected:

  double random() {
    return std::uniform_real_distribution<double>(0., 1.)(rng);
  }

  void permute_to(int iqbit, int iqbit_to) {
    int delta = iqbit < iqbit_to ? 1 : -1;
    while (iqbit != iqbit_to) {
      swap(iqbit, iqbit + delta);
      iqbit = iqbit + delta;
    }
  }

  void kickback_ind(ITensor &tensor, const Index &ind) {
    auto ind_p = itensor::prime(ind);
    ITensor identity(ind, ind_p);
    for (int i = 1; i <= ind.dim(); ++i) {
      identity.set(ind(i), ind_p(i), 1.);
    }
    tensor *= identity;
    tensor.prime(-1, ind_p);
  }

  double wavefunc_inner() {
    ITensor inner = itensor::conj(legMats[0] * bondMats[0]) * legMats[0];
    for (int i = 1; i < num_qubits - 1; ++i) {
      inner = inner * itensor::conj(legMats[i] * bondMats[i]) * bondMats[i - 1] *
              legMats[i];
    }
    inner = inner * itensor::conj(legMats[num_qubits - 1]) * bondMats[num_qubits - 2] *
            legMats[num_qubits - 1];

    return inner.cplx().real();
  }

  double average(int iqbit, const ITensor &op_tensor) {
    ITensor inner;
    if (iqbit == 0) {
      auto bra = itensor::conj(legMats[0] * bondMats[0]) * op_tensor;
      bra.noPrime();
      inner = bra * legMats[0];
    } else {
      inner = itensor::conj(legMats[0] * bondMats[0]) * legMats[0];
    }
    for (int i = 1; i < num_qubits - 1; ++i) {
      if (i == iqbit) {
        auto bra = inner * itensor::conj(legMats[i] * bondMats[i]) * op_tensor;
        bra.noPrime();
        inner = bra * bondMats[i - 1] * legMats[i];
      } else {
        inner = inner * itensor::conj(legMats[i] * bondMats[i]) *
                bondMats[i - 1] * legMats[i];
      }
    }
    if (iqbit == num_qubits - 1) {
      auto bra = inner * itensor::conj(legMats[num_qubits - 1]) * op_tensor;
      bra.noPrime();
      inner = bra * bondMats[num_qubits - 2] * legMats[num_qubits - 1];
    } else {
      inner = inner * itensor::conj(legMats[num_qubits - 1]) *
              bondMats[num_qubits - 2] * legMats[num_qubits - 1];
    }
    return inner.cplx().real();
  }

  mt19937_64 rng;
  size_t num_qubits;

  double svdCutoff = 1e-16;
  vector<itensor::ITensor> legMats;
  vector<itensor::ITensor> bondMats;
  bool debug;
};

#endif


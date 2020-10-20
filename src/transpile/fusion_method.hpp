/**
 * This code is part of Qiskit.
 *
 * (C) Copyright IBM 2018, 2019, 2020.
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */

#ifndef _aer_transpile_fusion_method_hpp_
#define _aer_transpile_fusion_method_hpp_

#include <string>

#include "transpile/circuitopt.hpp"
#include "simulators/unitary/unitary_state.hpp"
#include "simulators/superoperator/superoperator_state.hpp"

namespace AER {
namespace Transpile {

using uint_t = uint_t;
using op_t = Operations::Op;
using optype_t = Operations::OpType;
using oplist_t = std::vector<op_t>;
using opset_t = Operations::OpSet;
using reg_t = std::vector<uint_t>;

class FusionMethod {
public:
  FusionMethod(uint_t max_qubit_ = 5, uint_t threshold_ = 15): max_qubit(max_qubit_), threshold(threshold_) {}
  virtual ~FusionMethod() {};

  virtual std::string name() const { return "unitary"; };

  virtual bool can_apply_fusion(const op_t& op) const;

  virtual op_t generate_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                         const reg_t &qubits) const;

  virtual op_t generate_diagonal_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                                  const reg_t &qubits) const;

  uint_t get_default_max_qubit() const { return max_qubit; }

  uint_t get_default_threshold_qubit() const { return threshold; }

  virtual bool support_diagonal() const { return true; }

protected:

  std::vector<op_t> remap_qubits(const std::vector<op_t>& fusioned_ops, const reg_t &qubits) const;

private:

  uint_t max_qubit;
  uint_t threshold;

};

bool FusionMethod::can_apply_fusion(const op_t& op) const {
  if (op.conditional)
    return false;

  switch (op.type) {
    case optype_t::matrix:
      return op.mats.size() == 1;
    case optype_t::gate:
      return QubitUnitary::StateOpSet.contains_gates(op.name);
    default:
      return false;
  }
}

std::vector<op_t> FusionMethod::remap_qubits(const std::vector<op_t>& fusioned_ops, const reg_t &qubits) const {

  std::unordered_map<uint_t, uint_t> qubit_mapping;
  for (size_t i = 0; i < qubits.size(); i++)
    qubit_mapping[qubits[i]] = i;

  std::vector<op_t> ret = fusioned_ops;
  for (auto & op: ret)
    for (size_t i = 0; i < op.qubits.size(); i++)
      op.qubits[i] = qubit_mapping[op.qubits[i]];

  return std::move(ret);
}


op_t FusionMethod::generate_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                                    const reg_t &qubits) const {
  // Run simulation
  RngEngine dummy_rng;
  ExperimentData dummy_data;

  // Unitary simulation
  QubitUnitary::State<> unitary_simulator;
  unitary_simulator.initialize_qreg(qubits.size());
  unitary_simulator.apply_ops(remap_qubits(fusioned_ops, qubits), dummy_data, dummy_rng);
  auto ret = Operations::make_unitary(qubits, unitary_simulator.qreg().move_to_matrix(),
                                  std::string("fusion"));
  return ret;
}

op_t FusionMethod::generate_diagonal_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                                      const reg_t &qubits) const {
  // Run simulation
  RngEngine dummy_rng;
  ExperimentData dummy_data;

  // Unitary simulation
  QubitUnitary::State<> unitary_simulator;
  unitary_simulator.initialize_qreg(qubits.size());
  unitary_simulator.apply_ops(remap_qubits(fusioned_ops, qubits), dummy_data, dummy_rng);

  auto mat = unitary_simulator.qreg().move_to_matrix();
  std::vector<complex_t> vec;
  vec.assign(1UL << qubits.size(), 0);
  for (size_t i = 0; i < vec.size(); ++i)
    vec[i] = mat(i, i);

  auto ret = Operations::make_diagonal(qubits, vec, std::string("fusion"));
  return ret;
}


class UnitaryMatrixFusionMethod : public FusionMethod {
public:
  UnitaryMatrixFusionMethod(uint_t max_qubit = 3, uint_t threshold = 8): FusionMethod(max_qubit, threshold) {}
  virtual ~UnitaryMatrixFusionMethod() {};

  virtual std::string name() const override { return "unitary"; };

};


class SuperopFusionMethod : public FusionMethod {
public:
  SuperopFusionMethod(uint_t max_qubit = 3, uint_t threshold = 8): FusionMethod(max_qubit, threshold) {}
  virtual ~SuperopFusionMethod() {};

  virtual std::string name() const override { return "superop"; };

  virtual bool can_apply_fusion(const op_t& op) const override;

  virtual op_t generate_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                         const reg_t &qubits) const override;

  virtual bool support_diagonal() const override { return false; }

private:
  virtual op_t generate_noise_operation(const reg_t &qubits,
                                        matrix<std::complex<double>>& superop) const;

  const static opset_t noise_opset_;
};

const opset_t SuperopFusionMethod::noise_opset_ = opset_t(
    {Operations::OpType::kraus,
    Operations::OpType::superop,
    Operations::OpType::reset}, {}, {});

bool SuperopFusionMethod::can_apply_fusion(const op_t& op) const {
  if (op.conditional)
    return false;

  switch (op.type) {
    case optype_t::kraus:
    case optype_t::reset:
    case optype_t::superop:
      return true;
    case optype_t::gate:
      return QubitSuperoperator::StateOpSet.contains_gates(op.name);
    default:
      return FusionMethod::can_apply_fusion(op);
  }
}

op_t SuperopFusionMethod::generate_fusion_operation(const std::vector<op_t>& fusioned_ops,
                                              const reg_t &qubits) const {
  bool non_unitary = true;
  for (auto & op: fusioned_ops)
    non_unitary &= (!noise_opset_.contains(op.type));

  if (!non_unitary)
    return FusionMethod::generate_fusion_operation(fusioned_ops, qubits);

  // Run simulation
  RngEngine dummy_rng;
  ExperimentData dummy_data;

  QubitSuperoperator::State<> superop_simulator;
  superop_simulator.initialize_qreg(qubits.size());
  superop_simulator.apply_ops(remap_qubits(fusioned_ops, qubits), dummy_data, dummy_rng);
  auto superop = superop_simulator.qreg().move_to_matrix();
  return Operations::make_superop(qubits, std::move(superop));
}

op_t SuperopFusionMethod::generate_noise_operation(const reg_t &qubits,
                                                   matrix<std::complex<double>>& superop) const {
  return Operations::make_superop(qubits, std::move(superop));
}

class KrausFusionMethod: public SuperopFusionMethod{
public:
  KrausFusionMethod(uint_t max_qubit = 3, uint_t threshold = 8): SuperopFusionMethod(max_qubit, threshold) {}
  virtual ~KrausFusionMethod() {};

  virtual std::string name() const override { return "kraus"; };

private:
  virtual op_t generate_noise_operation(const reg_t &qubits, matrix<std::complex<double>>& superop) const override;
};

op_t KrausFusionMethod::generate_noise_operation(const reg_t &qubits,
                                                 matrix<std::complex<double>>& superop) const {
  size_t dim = 1 << qubits.size();
  return Operations::make_kraus(qubits, Utils::superop2kraus(superop, dim));
}


}
}

#endif

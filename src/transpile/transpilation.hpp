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

#ifndef _aer_transpiler_hpp_
#define _aer_transpiler_hpp_

#include "framework/circuit.hpp"
#include "noise/noise_model.hpp"

namespace AER {
namespace Transpile {

using int_t = int_fast64_t;
using uint_t = uint_fast64_t;
using op_t = Operations::Op;
using optype_t = Operations::OpType;
using oplist_t = std::vector<op_t>;
using opset_t = Operations::OpSet;
using reg_t = std::vector<uint_t>;

#define DEFAULT_THRESHOLD_QUBIT 16

class Transpilation {
public:
  Transpilation();
  ~Transpilation() = default;

  virtual void set_config(const json_t &config) { config_ = config; };
  uint_t threshold() { return threshold_;  }

protected:
  json_t config_;
  uint_t threshold_;
};

Transpilation::Transpilation(): threshold_(DEFAULT_THRESHOLD_QUBIT) {
}

class CircuitTranspilation: public Transpilation {
public:
  CircuitTranspilation() = default;
  virtual ~CircuitTranspilation() = default;

  virtual void transpile_circuit(Circuit& circ,
                                 Noise::NoiseModel& noise,
                                 OutputData& data) const;
};

void CircuitTranspilation::transpile_circuit(Circuit& circ,
                                      Noise::NoiseModel& noise,
                                      OutputData &data) const {
}

class ShotTranspilation: public Transpilation {
public:
  ShotTranspilation() = default;
  virtual ~ShotTranspilation() = default;

  virtual void transpile_shot(std::vector<op_t>& ops,
                              const Operations::OpSet &opset,
                              OutputData& data) const;
};

void ShotTranspilation::transpile_shot(std::vector<op_t>& ops,
                               const Operations::OpSet &opset,
                               OutputData& data) const{
}

//-------------------------------------------------------------------------
} // end namespace transpilation
} // end namespace AER
//-------------------------------------------------------------------------
#endif

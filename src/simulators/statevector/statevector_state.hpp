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


#ifndef _statevector_state_hpp
#define _statevector_state_hpp

#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

#include "framework/json.hpp"
#include "framework/utils.hpp"
#include "simulators/state.hpp"
#include "qubitvector.hpp"
#ifdef AER_THRUST_SUPPORTED
#include "qubitvector_thrust.hpp"
#endif

#include "batched_matrix.hpp"

namespace AER {

//predefinition of StatevectorChunk::State for friend class declaration to access static members
namespace StatevectorChunk {
template <class statevec_t> class State;
}

namespace Statevector {

using OpType = Operations::OpType;

// OpSet of supported instructions
const Operations::OpSet StateOpSet(
    // Op types
    {OpType::gate, OpType::measure,
     OpType::reset, OpType::initialize,
     OpType::snapshot, OpType::barrier,
     OpType::bfunc, OpType::roerror,
     OpType::matrix, OpType::diagonal_matrix,
     OpType::multiplexer, OpType::kraus, OpType::qerror_loc,
     OpType::sim_op, OpType::set_statevec,
     OpType::save_expval, OpType::save_expval_var,
     OpType::save_probs, OpType::save_probs_ket,
     OpType::save_amps, OpType::save_amps_sq,
     OpType::save_state, OpType::save_statevec,
     OpType::save_statevec_dict, OpType::save_densmat,
     },
    // Gates
    {"u1",     "u2",      "u3",  "u",    "U",    "CX",   "cx",   "cz",
     "cy",     "cp",      "cu1", "cu2",  "cu3",  "swap", "id",   "p",
     "x",      "y",       "z",   "h",    "s",    "sdg",  "t",    "tdg",
     "r",      "rx",      "ry",  "rz",   "rxx",  "ryy",  "rzz",  "rzx",
     "ccx",    "cswap",   "mcx", "mcy",  "mcz",  "mcu1", "mcu2", "mcu3",
     "mcswap", "mcphase", "mcr", "mcrx", "mcry", "mcry", "sx",   "sxdg",
     "csx", "mcsx", "csxdg", "mcsxdg",  "delay", "pauli", "mcx_gray", "cu", "mcu", "mcp"},
    // Snapshots
    {"statevector", "memory", "register", "probabilities",
     "probabilities_with_variance", "expectation_value_pauli", "density_matrix",
     "density_matrix_with_variance", "expectation_value_pauli_with_variance",
     "expectation_value_matrix_single_shot", "expectation_value_matrix",
     "expectation_value_matrix_with_variance",
     "expectation_value_pauli_single_shot"});

// Allowed gates enum class
enum class Gates {
  id, h, s, sdg, t, tdg,
  rxx, ryy, rzz, rzx,
  mcx, mcy, mcz, mcr, mcrx, mcry,
  mcrz, mcp, mcu2, mcu3, mcu, mcswap, mcsx, mcsxdg, pauli
};

// Allowed snapshots enum class
enum class Snapshots {
  statevector,
  cmemory,
  cregister,
  probs,
  probs_var,
  densmat,
  densmat_var,
  expval_pauli,
  expval_pauli_var,
  expval_pauli_shot,
  expval_matrix,
  expval_matrix_var,
  expval_matrix_shot
};

// Enum class for different types of expectation values
enum class SnapshotDataType { average, average_var, pershot };

//=========================================================================
// QubitVector State subclass
//=========================================================================

template <class statevec_t = QV::QubitVector<double>>
class State : public Base::State<statevec_t> {
  friend class StatevectorChunk::State<statevec_t>;
public:
  using BaseState = Base::State<statevec_t>;

  State() : BaseState(StateOpSet) {}
  virtual ~State() = default;

  //-----------------------------------------------------------------------
  // Base class overrides
  //-----------------------------------------------------------------------

  // Return the string name of the State class
  virtual std::string name() const override { return statevec_t::name(); }

  // Apply an operation
  // If the op is not in allowed_ops an exeption will be raised.
  virtual void apply_op(const Operations::Op &op,
                        ExperimentResult &result,
                        RngEngine &rng,
                        bool final_op = false) override;

  //applying one operation
  virtual void apply_op_multi_shots(const Operations::Op &op,
                         ExperimentResult &result,
                         std::vector<RngEngine>& rng,
                         bool final_ops = false) override;

  //batched execution
  virtual void apply_batched_ops(const std::vector<Operations::Op> &ops);
  //apply runtime sampled noise in case all inserted ops are Pauli gates
  virtual void apply_batched_pauli(std::vector<std::vector<Operations::Op>> &ops);
  //apply runtime sampled noise in case general ops are inserted 
  virtual void apply_batched_noise_ops(const std::vector<std::vector<Operations::Op>> &op, ExperimentResult &result,
                                               std::vector<RngEngine> &rng);

  virtual bool batchable_op(const Operations::Op& op,bool single_op = true);

  //store asynchronously measured classical bits after batched execution
  virtual void store_measured_cbits(void);

  // Initializes an n-qubit state to the all |0> state
  virtual void initialize_qreg(uint_t num_qubits) override;

  // Initializes to a specific n-qubit state
  virtual void initialize_qreg(uint_t num_qubits,
                               const statevec_t &state) override;

  // Returns the required memory for storing an n-qubit state in megabytes.
  // For this state the memory is independent of the number of ops
  // and is approximately 16 * 1 << num_qubits bytes
  virtual size_t
  required_memory_mb(uint_t num_qubits,
                     const std::vector<Operations::Op> &ops) const override;

  // Load the threshold for applying OpenMP parallelization
  // if the controller/engine allows threads for it
  virtual void set_config(const json_t &config) override;

  // Sample n-measurement outcomes without applying the measure operation
  // to the system state
  virtual std::vector<reg_t> sample_measure(const reg_t &qubits, uint_t shots,
                                            RngEngine &rng) override;
  virtual reg_t local_sample_measure(const reg_t &qubits,
                                            std::vector<double>& rnds);

  virtual reg_t batched_sample_measure(const reg_t &qubits,reg_t& shots,std::vector<RngEngine> &rng);

  virtual double sum(void);

  virtual bool allocate(uint_t num_qubits,uint_t block_bits,uint_t num_parallel_shots = 1) override;
  virtual bool bind_state(State<statevec_t>& state,uint_t ishot,bool batch_enable);

  virtual void end_of_circuit()
  {
    BaseState::qreg_.end_of_circuit();
  }
  
  //cache control for chunks on host
  virtual bool fetch_state(void) const
  {
    return BaseState::qreg_.fetch_chunk();
  }
  virtual void release_state(bool write_back = true) const
  {
    BaseState::qreg_.release_chunk();
  }

  //swap between state
  void apply_state_swap(const reg_t &qubits, State<statevec_t> &chunk, bool write_back = true);
  void apply_state_swap(const reg_t &qubits, uint_t remote_chunk_index);

  //-----------------------------------------------------------------------
  // Additional methods
  //-----------------------------------------------------------------------

  // Initializes to a specific n-qubit state given as a complex std::vector
  virtual void initialize_qreg(uint_t num_qubits, const cvector_t &state);

  // Initialize OpenMP settings for the underlying QubitVector class
  void initialize_omp();

  auto move_to_vector()
  {
    return BaseState::qreg_.move_to_vector();
  }

  //check if this register is on the top of group
  bool top_of_group() override
  {
    return BaseState::qreg_.top_of_group();
  }


  void enable_batch(bool flg)
  {
    BaseState::qreg_.enable_batch(flg);
  }

  //-----------------------------------------------------------------------
  // ClassicalRegister methods
  //-----------------------------------------------------------------------

  // Initialize classical memory and register to default value (all-0)
  virtual void initialize_creg(uint_t num_memory, uint_t num_register);

  // Initialize classical memory and register to specific values
  virtual void initialize_creg(uint_t num_memory,
                       uint_t num_register,
                       const std::string &memory_hex,
                       const std::string &register_hex);

protected:
  //-----------------------------------------------------------------------
  // Apply instructions
  //-----------------------------------------------------------------------

  // Applies a sypported Gate operation to the state class.
  // If the input is not in allowed_gates an exeption will be raised.
  void apply_gate(const Operations::Op &op);

  // Measure qubits and return a list of outcomes [q0, q1, ...]
  // If a state subclass supports this function it then "measure"
  // should be contained in the set returned by the 'allowed_ops'
  // method.
  virtual void apply_measure(const reg_t &qubits, const reg_t &cmemory,
                             const reg_t &cregister, RngEngine &rng);

  // Reset the specified qubits to the |0> state by simulating
  // a measurement, applying a conditional x-gate if the outcome is 1, and
  // then discarding the outcome.
  void apply_reset(const reg_t &qubits, RngEngine &rng);

  // Initialize the specified qubits to a given state |psi>
  // by applying a reset to the these qubits and then
  // computing the tensor product with the new state |psi>
  // /psi> is given in params
  void apply_initialize(const reg_t &qubits, const cvector_t &params,
                        RngEngine &rng);

  // Apply a supported snapshot instruction
  // If the input is not in allowed_snapshots an exeption will be raised.
  virtual void apply_snapshot(const Operations::Op &op, ExperimentResult &result, bool last_op = false);

  // Apply a matrix to given qubits (identity on all other qubits)
  void apply_matrix(const Operations::Op &op);

  // Apply a vectorized matrix to given qubits (identity on all other qubits)
  void apply_matrix(const reg_t &qubits, const cvector_t &vmat);

  // Apply a vector of control matrices to given qubits (identity on all other
  // qubits)
  void apply_multiplexer(const reg_t &control_qubits,
                         const reg_t &target_qubits,
                         const std::vector<cmatrix_t> &mmat);

  // Apply stacked (flat) version of multiplexer matrix to target qubits (using
  // control qubits to select matrix instance)
  void apply_multiplexer(const reg_t &control_qubits,
                         const reg_t &target_qubits, const cmatrix_t &mat);

  // Apply a Kraus error operation
  void apply_kraus(const reg_t &qubits, const std::vector<cmatrix_t> &krausops,
                   RngEngine &rng);

  virtual void apply_bfunc(const Operations::Op &op);

  //-----------------------------------------------------------------------
  // Save data instructions
  //-----------------------------------------------------------------------

  // Save the current state of the statevector simulator
  // If `last_op` is True this will use move semantics to move the simulator
  // state to the results, otherwise it will use copy semantics to leave
  // the current simulator state unchanged.
  void apply_save_statevector(const Operations::Op &op,
                              ExperimentResult &result,
                              bool last_op);

  // Save the current state of the statevector simulator as a ket-form map.
  void apply_save_statevector_dict(const Operations::Op &op,
                                  ExperimentResult &result);

  // Save the current density matrix or reduced density matrix
  void apply_save_density_matrix(const Operations::Op &op,
                                 ExperimentResult &result);

  // Helper function for computing expectation value
  void apply_save_probs(const Operations::Op &op,
                        ExperimentResult &result);

  // Helper function for saving amplitudes and amplitudes squared
  void apply_save_amplitudes(const Operations::Op &op,
                             ExperimentResult &result);

  // Helper function for computing expectation value
  virtual double expval_pauli(const reg_t &qubits,
                              const std::string& pauli) override;
  //-----------------------------------------------------------------------
  // Measurement Helpers
  //-----------------------------------------------------------------------

  // Return vector of measure probabilities for specified qubits
  // If a state subclass supports this function it then "measure"
  // should be contained in the set returned by the 'allowed_ops'
  // method.
  // TODO: move to private (no longer part of base class)
  rvector_t measure_probs(const reg_t &qubits) const;

  // Sample the measurement outcome for qubits
  // return a pair (m, p) of the outcome m, and its corresponding
  // probability p.
  // Outcome is given as an int: Eg for two-qubits {q0, q1} we have
  // 0 -> |q1 = 0, q0 = 0> state
  // 1 -> |q1 = 0, q0 = 1> state
  // 2 -> |q1 = 1, q0 = 0> state
  // 3 -> |q1 = 1, q0 = 1> state
  std::pair<uint_t, double> sample_measure_with_prob(const reg_t &qubits,
                                                     RngEngine &rng);

  void measure_reset_update(const std::vector<uint_t> &qubits,
                            const uint_t final_state, const uint_t meas_state,
                            const double meas_prob);

  //-----------------------------------------------------------------------
  // Special snapshot types
  // Apply a supported snapshot instruction
  //
  // IMPORTANT: These methods are not marked const to allow modifying state
  // during snapshot, but after the snapshot is applied the simulator
  // should be left in the pre-snapshot state.
  //-----------------------------------------------------------------------

  // Snapshot current qubit probabilities for a measurement (average)
  void snapshot_probabilities(const Operations::Op &op, ExperimentResult &result,
                              SnapshotDataType type);

  // Snapshot the expectation value of a Pauli operator
  void snapshot_pauli_expval(const Operations::Op &op, ExperimentResult &result,
                             SnapshotDataType type);

  // Snapshot the expectation value of a matrix operator
  void snapshot_matrix_expval(const Operations::Op &op, ExperimentResult &result,
                              SnapshotDataType type);

  // Snapshot reduced density matrix
  void snapshot_density_matrix(const Operations::Op &op, ExperimentResult &result,
                               SnapshotDataType type);

  // Return the reduced density matrix for the simulator
  cmatrix_t density_matrix(const reg_t &qubits);

  // Helper function to convert a vector to a reduced density matrix
  template <class T> cmatrix_t vec2density(const reg_t &qubits, const T &vec);

  //-----------------------------------------------------------------------
  // Single-qubit gate helpers
  //-----------------------------------------------------------------------

  // Optimize phase gate with diagonal [1, phase]
  void apply_gate_phase(const uint_t qubit, const complex_t phase);

  //-----------------------------------------------------------------------
  // Multi-controlled u3
  //-----------------------------------------------------------------------

  // Apply N-qubit multi-controlled single qubit gate specified by
  // 4 parameters u4(theta, phi, lambda, gamma)
  // NOTE: if N=1 this is just a regular u4 gate.
  void apply_gate_mcu(const reg_t &qubits, const double theta,
                      const double phi, const double lambda,
                      const double gamma);

  //-----------------------------------------------------------------------
  // Config Settings
  //-----------------------------------------------------------------------

  // Apply the global phase
  void apply_global_phase();

  // OpenMP qubit threshold
  int omp_qubit_threshold_ = 14;

  // QubitVector sample measure index size
  int sample_measure_index_size_ = 10;

  // Threshold for chopping small values to zero in JSON
  double json_chop_threshold_ = 1e-10;

  // Table of allowed gate names to gate enum class members
  const static stringmap_t<Gates> gateset_;

  // Table of allowed snapshot types to enum class members
  const static stringmap_t<Snapshots> snapshotset_;
};

//=========================================================================
// Implementation: Allowed ops and gateset
//=========================================================================

template <class statevec_t>
const stringmap_t<Gates> State<statevec_t>::gateset_({
    // 1-qubit gates
    {"delay", Gates::id},// Delay gate
    {"id", Gates::id},   // Pauli-Identity gate
    {"x", Gates::mcx},   // Pauli-X gate
    {"y", Gates::mcy},   // Pauli-Y gate
    {"z", Gates::mcz},   // Pauli-Z gate
    {"s", Gates::s},     // Phase gate (aka sqrt(Z) gate)
    {"sdg", Gates::sdg}, // Conjugate-transpose of Phase gate
    {"h", Gates::h},     // Hadamard gate (X + Z / sqrt(2))
    {"t", Gates::t},     // T-gate (sqrt(S))
    {"tdg", Gates::tdg}, // Conjguate-transpose of T gate
    {"p", Gates::mcp},   // Parameterized phase gate 
    {"sx", Gates::mcsx}, // Sqrt(X) gate
    {"sxdg", Gates::mcsxdg}, // Inverse Sqrt(X) gate
    // 1-qubit rotation Gates
    {"r", Gates::mcr},   // R rotation gate
    {"rx", Gates::mcrx}, // Pauli-X rotation gate
    {"ry", Gates::mcry}, // Pauli-Y rotation gate
    {"rz", Gates::mcrz}, // Pauli-Z rotation gate
    // Waltz Gates
    {"u1", Gates::mcp},  // zero-X90 pulse waltz gate
    {"u2", Gates::mcu2}, // single-X90 pulse waltz gate
    {"u3", Gates::mcu3}, // two X90 pulse waltz gate
    {"u", Gates::mcu3}, // two X90 pulse waltz gate
    {"U", Gates::mcu3}, // two X90 pulse waltz gate
    // 2-qubit gates
    {"CX", Gates::mcx},      // Controlled-X gate (CNOT)
    {"cx", Gates::mcx},      // Controlled-X gate (CNOT)
    {"cy", Gates::mcy},      // Controlled-Y gate
    {"cz", Gates::mcz},      // Controlled-Z gate
    {"cp", Gates::mcp},      // Controlled-Phase gate 
    {"cu1", Gates::mcp},     // Controlled-u1 gate
    {"cu2", Gates::mcu2},    // Controlled-u2 gate
    {"cu3", Gates::mcu3},    // Controlled-u3 gate
    {"cu", Gates::mcu},      // Controlled-u4 gate
    {"cp", Gates::mcp},      // Controlled-Phase gate 
    {"swap", Gates::mcswap}, // SWAP gate
    {"rxx", Gates::rxx},     // Pauli-XX rotation gate
    {"ryy", Gates::ryy},     // Pauli-YY rotation gate
    {"rzz", Gates::rzz},     // Pauli-ZZ rotation gate
    {"rzx", Gates::rzx},     // Pauli-ZX rotation gate
    {"csx", Gates::mcsx},    // Controlled-Sqrt(X) gate
    {"csxdg", Gates::mcsxdg}, // Controlled-Sqrt(X)dg gate
    // 3-qubit gates
    {"ccx", Gates::mcx},      // Controlled-CX gate (Toffoli)
    {"cswap", Gates::mcswap}, // Controlled SWAP gate (Fredkin)
    // Multi-qubit controlled gates
    {"mcx", Gates::mcx},      // Multi-controlled-X gate
    {"mcy", Gates::mcy},      // Multi-controlled-Y gate
    {"mcz", Gates::mcz},      // Multi-controlled-Z gate
    {"mcr", Gates::mcr},      // Multi-controlled R-rotation gate
    {"mcrx", Gates::mcrx},    // Multi-controlled X-rotation gate
    {"mcry", Gates::mcry},    // Multi-controlled Y-rotation gate
    {"mcrz", Gates::mcrz},    // Multi-controlled Z-rotation gate
    {"mcphase", Gates::mcp},  // Multi-controlled-Phase gate 
    {"mcp", Gates::mcp},      // Multi-controlled-Phase gate 
    {"mcu1", Gates::mcp},     // Multi-controlled-u1
    {"mcu2", Gates::mcu2},    // Multi-controlled-u2
    {"mcu3", Gates::mcu3},    // Multi-controlled-u3
    {"mcu", Gates::mcu},      // Multi-controlled-u4
    {"mcswap", Gates::mcswap},// Multi-controlled SWAP gate
    {"mcsx", Gates::mcsx},    // Multi-controlled-Sqrt(X) gate
    {"mcsxdg", Gates::mcsxdg}, // Multi-controlled-Sqrt(X)dg gate
    {"pauli", Gates::pauli},   // Multi-qubit Pauli gate
    {"mcx_gray", Gates::mcx}
});

template <class statevec_t>
const stringmap_t<Snapshots> State<statevec_t>::snapshotset_(
    {{"statevector", Snapshots::statevector},
     {"probabilities", Snapshots::probs},
     {"expectation_value_pauli", Snapshots::expval_pauli},
     {"expectation_value_matrix", Snapshots::expval_matrix},
     {"probabilities_with_variance", Snapshots::probs_var},
     {"density_matrix", Snapshots::densmat},
     {"density_matrix_with_variance", Snapshots::densmat_var},
     {"expectation_value_pauli_with_variance", Snapshots::expval_pauli_var},
     {"expectation_value_matrix_with_variance", Snapshots::expval_matrix_var},
     {"expectation_value_pauli_single_shot", Snapshots::expval_pauli_shot},
     {"expectation_value_matrix_single_shot", Snapshots::expval_matrix_shot},
     {"memory", Snapshots::cmemory},
     {"register", Snapshots::cregister}});

//=========================================================================
// Implementation: Base class method overrides
//=========================================================================

//-------------------------------------------------------------------------
// Initialization
//-------------------------------------------------------------------------
template <class statevec_t>
bool State<statevec_t>::allocate(uint_t num_qubits,uint_t block_bits,uint_t num_parallel_shots)
{
  BaseState::qreg_.set_max_matrix_bits(BaseState::max_matrix_bits_);
  BaseState::num_qubits_ = num_qubits;
  BaseState::num_state_qubits_ = block_bits;
  bool ret = BaseState::qreg_.chunk_setup(block_bits,num_qubits,BaseState::state_index_,num_parallel_shots);

  return ret;
}

template <class statevec_t>
bool State<statevec_t>::bind_state(State<statevec_t>& state,uint_t ishot,bool batch_enable)
{
  //allocate qreg from allocated buffer
  if(BaseState::qreg_.chunk_setup(state.qreg_,ishot)){
    BaseState::qreg_.enable_batch(batch_enable);
    state.qreg_.enable_batch(batch_enable);

    BaseState::state_index_ = ishot;
    BaseState::num_qubits_ = state.num_qubits_;
    BaseState::num_state_qubits_ = state.num_state_qubits_;

    BaseState::max_matrix_bits_ = state.max_matrix_bits_;

    return true;
  }
  return false;
}

template <class statevec_t>
void State<statevec_t>::initialize_qreg(uint_t num_qubits) {
  initialize_omp();
  BaseState::qreg_.set_num_qubits(num_qubits);
  if(BaseState::state_index_ == 0)
    BaseState::qreg_.initialize();
  else
    BaseState::qreg_.zero();
  apply_global_phase();
}

template <class statevec_t>
void State<statevec_t>::initialize_qreg(uint_t num_qubits,
                                        const statevec_t &state) {
  // Check dimension of state
  if (state.num_qubits() != num_qubits) {
    throw std::invalid_argument("QubitVector::State::initialize: initial state "
                                "does not match qubit number");
  }
  initialize_omp();
  BaseState::qreg_.set_num_qubits(num_qubits);
  BaseState::qreg_.initialize_from_data(state.data(), 1ULL << num_qubits);
  apply_global_phase();
}

template <class statevec_t>
void State<statevec_t>::initialize_qreg(uint_t num_qubits,
                                        const cvector_t &state) {
  if (state.size() != 1ULL << num_qubits) {
    throw std::invalid_argument("QubitVector::State::initialize: initial state "
                                "does not match qubit number");
  }
  initialize_omp();
  BaseState::qreg_.set_num_qubits(num_qubits);
  BaseState::qreg_.initialize_from_vector(state);
  apply_global_phase();
}

template <class statevec_t> void State<statevec_t>::initialize_omp() {
  BaseState::qreg_.set_omp_threshold(omp_qubit_threshold_);
  if (BaseState::threads_ > 0)
    BaseState::qreg_.set_omp_threads(
        BaseState::threads_); // set allowed OMP threads in qubitvector
}

template <class state_t>
void State<state_t>::initialize_creg(uint_t num_memory, uint_t num_register) 
{
  BaseState::initialize_creg(num_memory, num_register);

  //initialize creg on GPU memory
  BaseState::qreg_.initialize_creg(num_memory, num_register);
}


template <class state_t>
void State<state_t>::initialize_creg(uint_t num_memory,
                                     uint_t num_register,
                                     const std::string &memory_hex,
                                     const std::string &register_hex) 
{
  BaseState::initialize_creg(num_memory, num_register, memory_hex, register_hex);
  BaseState::qreg_.initialize_creg(num_memory, num_register);
}

//-------------------------------------------------------------------------
// Utility
//-------------------------------------------------------------------------

template <class statevec_t>
void State<statevec_t>::apply_global_phase() {
  if (BaseState::has_global_phase_) {
    BaseState::qreg_.apply_diagonal_matrix({0}, {BaseState::global_phase_, BaseState::global_phase_});
  }
}

template <class statevec_t>
size_t State<statevec_t>::required_memory_mb(uint_t num_qubits,
                                             const std::vector<Operations::Op> &ops)
                                             const {
  (void)ops; // avoid unused variable compiler warning
  return BaseState::qreg_.required_memory_mb(num_qubits);
}

template <class statevec_t>
void State<statevec_t>::set_config(const json_t &config) {
  BaseState::set_config(config);

  // Set threshold for truncating snapshots
  JSON::get_value(json_chop_threshold_, "zero_threshold", config);
  BaseState::qreg_.set_json_chop_threshold(json_chop_threshold_);

  // Set OMP threshold for state update functions
  JSON::get_value(omp_qubit_threshold_, "statevector_parallel_threshold",
                  config);

  // Set the sample measure indexing size
  int index_size;
  if (JSON::get_value(index_size, "statevector_sample_measure_opt", config)) {
    BaseState::qreg_.set_sample_measure_index_size(index_size);
  };
}

template <class statevec_t>
void State<statevec_t>::apply_bfunc(const Operations::Op &op)
{
  if(BaseState::qreg_.batched_optimization_supported()){
    BaseState::qreg_.apply_bfunc(op);
  }
  else
    BaseState::creg_.apply_bfunc(op);
}


//=========================================================================
// Implementation: apply operations
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_op(const Operations::Op &op,
                                  ExperimentResult &result,
                                  RngEngine &rng,
                                  bool final_ops) 
{
  if(BaseState::creg_.check_conditional(op)) {
    switch (op.type) {
      case OpType::barrier:
      case OpType::nop:
      case OpType::qerror_loc:
        break;
      case OpType::reset:
        apply_reset(op.qubits, rng);
        break;
      case OpType::initialize:
        apply_initialize(op.qubits, op.params, rng);
        break;
      case OpType::measure:
        apply_measure(op.qubits, op.memory, op.registers, rng);
        break;
      case OpType::bfunc:
        BaseState::creg_.apply_bfunc(op);
        break;
      case OpType::roerror:
        BaseState::creg_.apply_roerror(op, rng);
        break;
      case OpType::gate:
        apply_gate(op);
        break;
      case OpType::snapshot:
        apply_snapshot(op, result, final_ops);
        break;
      case OpType::matrix:
        apply_matrix(op);
        break;
      case OpType::diagonal_matrix:
        BaseState::qreg_.apply_diagonal_matrix(op.qubits, op.params);
        break;
      case OpType::multiplexer:
        apply_multiplexer(op.regs[0], op.regs[1],
                          op.mats); // control qubits ([0]) & target qubits([1])
        break;
      case OpType::kraus:
        apply_kraus(op.qubits, op.mats, rng);
        break;
      case OpType::sim_op:
        if(op.name == "begin_register_blocking"){
          BaseState::qreg_.enter_register_blocking(op.qubits);
        }
        else if(op.name == "end_register_blocking"){
          BaseState::qreg_.leave_register_blocking();
        }
        break;
      case OpType::set_statevec:
        BaseState::qreg_.initialize_from_vector(op.params);
        break;
      case OpType::save_expval:
      case OpType::save_expval_var:
        BaseState::apply_save_expval(op, result);
        break;
      case OpType::save_densmat:
        apply_save_density_matrix(op, result);
        break;
      case OpType::save_state:
      case OpType::save_statevec:
        apply_save_statevector(op, result, final_ops);
        break;
      case OpType::save_statevec_dict:
        apply_save_statevector_dict(op, result);
        break;
      case OpType::save_probs:
      case OpType::save_probs_ket:
        apply_save_probs(op, result);
        break;
      case OpType::save_amps:
      case OpType::save_amps_sq:
        apply_save_amplitudes(op, result);
        break;
      default:
        throw std::invalid_argument(
            "QubitVector::State::invalid instruction \'" + op.name + "\'.");
    }
  }
}

template <class statevec_t>
void State<statevec_t>::apply_op_multi_shots(const Operations::Op &op,
                                  ExperimentResult &result,
                                  std::vector<RngEngine> &rng,
                                  bool final_ops) 
{
  if(op.conditional){
    BaseState::qreg_.set_conditional(op.conditional_reg);
  }

  switch (op.type) {
    case OpType::barrier:
    case OpType::nop:
      break;
    case OpType::reset:
      BaseState::qreg_.apply_batched_reset(op.qubits,rng);
      break;
    case OpType::initialize:
      BaseState::qreg_.apply_batched_reset(op.qubits,rng);
      BaseState::qreg_.initialize_component(op.qubits, op.params);
      break;
    case OpType::measure:
      BaseState::qreg_.apply_batched_measure(op.qubits,rng,op.memory,op.registers);
      break;
    case OpType::bfunc:
      BaseState::qreg_.apply_bfunc(op);
      break;
    case OpType::roerror:
      BaseState::qreg_.apply_roerror(op, rng);
      break;
    case OpType::gate:
      apply_gate(op);
      break;
    case OpType::snapshot:
      apply_snapshot(op, result, final_ops);
      break;
    case OpType::matrix:
      apply_matrix(op);
      break;
    case OpType::diagonal_matrix:
      BaseState::qreg_.apply_diagonal_matrix(op.qubits, op.params);
      break;
    case OpType::multiplexer:
      apply_multiplexer(op.regs[0], op.regs[1],
                        op.mats); // control qubits ([0]) & target qubits([1])
      break;
    case OpType::kraus:
      BaseState::qreg_.apply_batched_kraus(op.qubits, op.mats,rng);
      break;
    case OpType::sim_op:
      if(op.name == "begin_register_blocking"){
        BaseState::qreg_.enter_register_blocking(op.qubits);
      }
      else if(op.name == "end_register_blocking"){
        BaseState::qreg_.leave_register_blocking();
      }
      break;
    case OpType::set_statevec:
      BaseState::qreg_.initialize_from_vector(op.params);
      break;
    default:
      throw std::invalid_argument(
          "QubitVector::State::invalid instruction \'" + op.name + "\'.");
  }
}

template <class statevec_t>
bool State<statevec_t>::batchable_op(const Operations::Op& op,bool single_op)
{
  if(op.type == OpType::set_statevec || op.type == OpType::save_expval || op.type == OpType::save_expval_var ||
     op.type == OpType::save_densmat || op.type == OpType::save_state || op.type == OpType::save_statevec ||
     op.type == OpType::save_statevec_dict || op.type == OpType::save_probs || op.type == OpType::save_probs_ket ||
     op.type == OpType::save_amps || op.type == OpType::save_amps_sq)
    return false;

  if(single_op)
    return true;

  if(op.type == OpType::bfunc || op.type == OpType::roerror)
    return false;
  if(op.type == OpType::gate && op.name == "pauli")
    return false;   //pauli can be only applied for single_op mode

  return true;
}

template <class statevec_t>
void State<statevec_t>::apply_state_swap(const reg_t &qubits, State<statevec_t> &chunk, bool write_back)
{
  if(qubits[0] < BaseState::num_state_qubits_ && qubits[1] < BaseState::num_state_qubits_){
    //local swap
    BaseState::qreg_.apply_mcswap(qubits);
  }
  else{
    BaseState::qreg_.apply_chunk_swap(qubits,chunk.qreg_,write_back);
  }
}

template <class statevec_t>
void State<statevec_t>::apply_state_swap(const reg_t &qubits, uint_t remote_chunk_index)
{
  BaseState::qreg_.apply_chunk_swap(qubits,remote_chunk_index);
}

//=========================================================================
// Implementation: Save data
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_save_probs(const Operations::Op &op,
                                         ExperimentResult &result) {
  // get probs as hexadecimal
  auto probs = measure_probs(op.qubits);
  if (op.type == OpType::save_probs_ket) {
    // Convert to ket dict
    BaseState::save_data_average(result, op.string_params[0],
                                 Utils::vec2ket(probs, json_chop_threshold_, 16),
                                 op.save_type);
  } else {
    BaseState::save_data_average(result, op.string_params[0],
                                 std::move(probs), op.save_type);
  }
}


template <class statevec_t>
double State<statevec_t>::expval_pauli(const reg_t &qubits,
                                       const std::string& pauli) {
  return BaseState::qreg_.expval_pauli(qubits, pauli);
}

template <class statevec_t>
void State<statevec_t>::apply_save_statevector(const Operations::Op &op,
                                               ExperimentResult &result,
                                               bool last_op) {
  if (op.qubits.size() != BaseState::qreg_.num_qubits()) {
    throw std::invalid_argument(
        op.name + " was not applied to all qubits."
        " Only the full statevector can be saved.");
  }
  std::string key = (op.string_params[0] == "_method_")
                      ? "statevector"
                      : op.string_params[0];
  if (last_op) {
    BaseState::save_data_pershot(result, key,
                                 BaseState::qreg_.move_to_vector(),
                                 op.save_type);
  } else {
    BaseState::save_data_pershot(result, key,
                                 BaseState::qreg_.copy_to_vector(),
                                 op.save_type);
  }
}

template <class statevec_t>
void State<statevec_t>::apply_save_statevector_dict(const Operations::Op &op,
                                                   ExperimentResult &result) {
  if (op.qubits.size() != BaseState::qreg_.num_qubits()) {
    throw std::invalid_argument(
        op.name + " was not applied to all qubits."
        " Only the full statevector can be saved.");
  }
  auto state_ket = BaseState::qreg_.vector_ket(json_chop_threshold_);
  std::map<std::string, complex_t> result_state_ket;
  for (auto const& it : state_ket){
    result_state_ket[it.first] = it.second;
  }
  BaseState::save_data_pershot(result, op.string_params[0],
                               std::move(result_state_ket), op.save_type);
}

template <class statevec_t>
void State<statevec_t>::apply_save_density_matrix(const Operations::Op &op,
                                                  ExperimentResult &result) {
  cmatrix_t reduced_state;

  // Check if tracing over all qubits
  if (op.qubits.empty()) {
    reduced_state = cmatrix_t(1, 1);
    reduced_state[0] = BaseState::qreg_.norm();
  } else {
    reduced_state = density_matrix(op.qubits);
  }

  BaseState::save_data_average(result, op.string_params[0],
                               std::move(reduced_state), op.save_type);
}

template <class statevec_t>
void State<statevec_t>::apply_save_amplitudes(const Operations::Op &op,
                                              ExperimentResult &result) {
  if (op.int_params.empty()) {
    throw std::invalid_argument("Invalid save_amplitudes instructions (empty params).");
  }
  const int_t size = op.int_params.size();
  if (op.type == OpType::save_amps) {
    Vector<complex_t> amps(size, false);
    #pragma omp parallel for if (size > pow(2, omp_qubit_threshold_) &&        \
                                 BaseState::threads_ > 1)                       \
                          num_threads(BaseState::threads_)
      for (int_t i = 0; i < size; ++i) {
        amps[i] = BaseState::qreg_.get_state(op.int_params[i]);
      }
    BaseState::save_data_pershot(result, op.string_params[0],
                                 std::move(amps), op.save_type);
  } else {
    rvector_t amps_sq(size);
    #pragma omp parallel for if (size > pow(2, omp_qubit_threshold_) &&        \
                                 BaseState::threads_ > 1)                       \
                          num_threads(BaseState::threads_)
      for (int_t i = 0; i < size; ++i) {
        amps_sq[i] = BaseState::qreg_.probability(op.int_params[i]);
      }
    BaseState::save_data_average(result, op.string_params[0],
                                 std::move(amps_sq), op.save_type);
  }
}
//=========================================================================
// Implementation: Snapshots
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_snapshot(const Operations::Op &op,
                                       ExperimentResult &result,
                                       bool last_op) {

  // Look for snapshot type in snapshotset
  auto it = snapshotset_.find(op.name);
  if (it == snapshotset_.end())
    throw std::invalid_argument(
        "QubitVectorState::invalid snapshot instruction \'" + op.name + "\'.");
  switch (it->second) {
    case Snapshots::statevector:
      if (last_op) {
        result.legacy_data.add_pershot_snapshot("statevector", op.string_params[0],
                                         BaseState::qreg_.move_to_vector());
      } else {
        result.legacy_data.add_pershot_snapshot("statevector", op.string_params[0],
                                         BaseState::qreg_.copy_to_vector());
      }
      break;
    case Snapshots::cmemory:
      BaseState::snapshot_creg_memory(op, result);
      break;
    case Snapshots::cregister:
      BaseState::snapshot_creg_register(op, result);
      break;
    case Snapshots::probs: {
      // get probs as hexadecimal
      snapshot_probabilities(op, result, SnapshotDataType::average);
    } break;
    case Snapshots::densmat: {
      snapshot_density_matrix(op, result, SnapshotDataType::average);
    } break;
    case Snapshots::expval_pauli: {
      snapshot_pauli_expval(op, result, SnapshotDataType::average);
    } break;
    case Snapshots::expval_matrix: {
      snapshot_matrix_expval(op, result, SnapshotDataType::average);
    } break;
    case Snapshots::probs_var: {
      // get probs as hexadecimal
      snapshot_probabilities(op, result, SnapshotDataType::average_var);
    } break;
    case Snapshots::densmat_var: {
      snapshot_density_matrix(op, result, SnapshotDataType::average_var);
    } break;
    case Snapshots::expval_pauli_var: {
      snapshot_pauli_expval(op, result, SnapshotDataType::average_var);
    } break;
    case Snapshots::expval_matrix_var: {
      snapshot_matrix_expval(op, result, SnapshotDataType::average_var);
    } break;
    case Snapshots::expval_pauli_shot: {
      snapshot_pauli_expval(op, result, SnapshotDataType::pershot);
    } break;
    case Snapshots::expval_matrix_shot: {
      snapshot_matrix_expval(op, result, SnapshotDataType::pershot);
    } break;
    default:
      // We shouldn't get here unless there is a bug in the snapshotset
      throw std::invalid_argument(
          "QubitVector::State::invalid snapshot instruction \'" + op.name +
          "\'.");
  }
}

template <class statevec_t>
void State<statevec_t>::snapshot_probabilities(const Operations::Op &op,
                                               ExperimentResult &result,
                                               SnapshotDataType type) {
  // get probs as hexadecimal
  auto probs =
      Utils::vec2ket(measure_probs(op.qubits), json_chop_threshold_, 16);
  bool variance = type == SnapshotDataType::average_var;
  result.legacy_data.add_average_snapshot("probabilities", op.string_params[0],
                                   BaseState::creg_.memory_hex(),
                                   std::move(probs), variance);
}

template <class statevec_t>
void State<statevec_t>::snapshot_pauli_expval(const Operations::Op &op,
                                              ExperimentResult &result,
                                              SnapshotDataType type) {
  // Check empty edge case
  if (op.params_expval_pauli.empty()) {
    throw std::invalid_argument(
        "Invalid expval snapshot (Pauli components are empty).");
  }

  // Accumulate expval components
  complex_t expval(0., 0.);
  for (const auto &param : op.params_expval_pauli) {
    const auto &coeff = param.first;
    const auto &pauli = param.second;
    expval += coeff * expval_pauli(op.qubits, pauli);
  }

  // Add to snapshot
  Utils::chop_inplace(expval, json_chop_threshold_);
  switch (type) {
  case SnapshotDataType::average:
    result.legacy_data.add_average_snapshot("expectation_value", op.string_params[0],
                              BaseState::creg_.memory_hex(), expval, false);
    break;
  case SnapshotDataType::average_var:
    result.legacy_data.add_average_snapshot("expectation_value", op.string_params[0],
                              BaseState::creg_.memory_hex(), expval, true);
    break;
  case SnapshotDataType::pershot:
    result.legacy_data.add_pershot_snapshot("expectation_values", op.string_params[0],
                              expval);
    break;
  }
}

template <class statevec_t>
void State<statevec_t>::snapshot_matrix_expval(const Operations::Op &op,
                                               ExperimentResult &result,
                                               SnapshotDataType type) {
  // Check empty edge case
  if (op.params_expval_matrix.empty()) {
    throw std::invalid_argument(
        "Invalid matrix snapshot (components are empty).");
  }
  reg_t qubits = op.qubits;
  // Cache the current quantum state
  BaseState::qreg_.checkpoint();
  bool first = true; // flag for first pass so we don't unnecessarily revert
                     // from checkpoint

  // Compute expval components
  complex_t expval(0., 0.);
  for (const auto &param : op.params_expval_matrix) {
    complex_t coeff = param.first;
    // Revert the quantum state to cached checkpoint
    if (first)
      first = false;
    else
      BaseState::qreg_.revert(true);
    // Apply each matrix component
    for (const auto &pair : param.second) {
      reg_t sub_qubits;
      for (const auto &pos : pair.first) {
        sub_qubits.push_back(qubits[pos]);
      }
      const cmatrix_t &mat = pair.second;
      cvector_t vmat =
          (mat.GetColumns() == 1)
              ? Utils::vectorize_matrix(Utils::projector(
                    Utils::vectorize_matrix(mat))) // projector case
              : Utils::vectorize_matrix(mat); // diagonal or square matrix case
      if (vmat.size() == 1ULL << qubits.size()) {
        BaseState::qreg_.apply_diagonal_matrix(sub_qubits, vmat);
      } else {
        BaseState::qreg_.apply_matrix(sub_qubits, vmat);
      }
    }
    expval += coeff * BaseState::qreg_.inner_product();
  }
  // add to snapshot
  Utils::chop_inplace(expval, json_chop_threshold_);
  switch (type) {
  case SnapshotDataType::average:
    result.legacy_data.add_average_snapshot("expectation_value", op.string_params[0],
                              BaseState::creg_.memory_hex(), expval, false);
    break;
  case SnapshotDataType::average_var:
    result.legacy_data.add_average_snapshot("expectation_value", op.string_params[0],
                              BaseState::creg_.memory_hex(), expval, true);
    break;
  case SnapshotDataType::pershot:
    result.legacy_data.add_pershot_snapshot("expectation_values", op.string_params[0],
                              expval);
    break;
  }
  // Revert to original state
  BaseState::qreg_.revert(false);
}

template <class statevec_t>
void State<statevec_t>::snapshot_density_matrix(const Operations::Op &op,
                                                ExperimentResult &result,
                                                SnapshotDataType type) {
  cmatrix_t reduced_state;

  // Check if tracing over all qubits
  if (op.qubits.empty()) {
    reduced_state = cmatrix_t(1, 1);
    reduced_state[0] = BaseState::qreg_.norm();
  } else {
    reduced_state = density_matrix(op.qubits);
  }

  // Add density matrix to result data
  switch (type) {
  case SnapshotDataType::average:
    result.legacy_data.add_average_snapshot("density_matrix", op.string_params[0],
                              BaseState::creg_.memory_hex(),
                              std::move(reduced_state), false);
    break;
  case SnapshotDataType::average_var:
    result.legacy_data.add_average_snapshot("density_matrix", op.string_params[0],
                              BaseState::creg_.memory_hex(),
                              std::move(reduced_state), true);
    break;
  case SnapshotDataType::pershot:
    result.legacy_data.add_pershot_snapshot("density_matrix", op.string_params[0],
                              std::move(reduced_state));
    break;
  }
}

template <class statevec_t>
cmatrix_t State<statevec_t>::density_matrix(const reg_t &qubits) {
  return vec2density(qubits, BaseState::qreg_.data());
}

#ifdef AER_THRUST_SUPPORTED
template <>
cmatrix_t State<QV::QubitVectorThrust<float>>::density_matrix(const reg_t &qubits) {
  return vec2density(qubits, BaseState::qreg_.copy_to_vector());
}

template <>
cmatrix_t State<QV::QubitVectorThrust<double>>::density_matrix(const reg_t &qubits) {
  return vec2density(qubits, BaseState::qreg_.copy_to_vector());
}
#endif

template <class statevec_t>
template <class T>
cmatrix_t State<statevec_t>::vec2density(const reg_t &qubits, const T &vec) {
  const size_t N = qubits.size();
  const size_t DIM = 1ULL << N;
  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());

  // Return full density matrix
  cmatrix_t densmat(DIM, DIM);
  if ((N == BaseState::qreg_.num_qubits()) && (qubits == qubits_sorted)) {
    const int_t mask = QV::MASKS[N];
#pragma omp parallel for if (2 * N > omp_qubit_threshold_ &&                   \
                             BaseState::threads_ > 1)                          \
    num_threads(BaseState::threads_)
    for (int_t rowcol = 0; rowcol < int_t(DIM * DIM); ++rowcol) {
      const int_t row = rowcol >> N;
      const int_t col = rowcol & mask;
      densmat(row, col) = complex_t(vec[row]) * complex_t(std::conj(vec[col]));
    }
  } else {
    const size_t END = 1ULL << (BaseState::qreg_.num_qubits() - N);
    // Initialize matrix values with first block
    {
      const auto inds = QV::indexes(qubits, qubits_sorted, 0);
      for (size_t row = 0; row < DIM; ++row)
        for (size_t col = 0; col < DIM; ++col) {
          densmat(row, col) =
              complex_t(vec[inds[row]]) * complex_t(std::conj(vec[inds[col]]));
        }
    }
    // Accumulate remaining blocks
    for (size_t k = 1; k < END; k++) {
      // store entries touched by U
      const auto inds = QV::indexes(qubits, qubits_sorted, k);
      for (size_t row = 0; row < DIM; ++row)
        for (size_t col = 0; col < DIM; ++col) {
          densmat(row, col) +=
              complex_t(vec[inds[row]]) * complex_t(std::conj(vec[inds[col]]));
        }
    }
  }
  return densmat;
}

//=========================================================================
// Implementation: Matrix multiplication
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_gate(const Operations::Op &op) {
  // Look for gate name in gateset
  auto it = gateset_.find(op.name);
  if (it == gateset_.end())
    throw std::invalid_argument(
        "QubitVectorState::invalid gate instruction \'" + op.name + "\'.");
  switch (it->second) {
    case Gates::mcx:
      // Includes X, CX, CCX, etc
      BaseState::qreg_.apply_mcx(op.qubits);
      break;
    case Gates::mcy:
      // Includes Y, CY, CCY, etc
      BaseState::qreg_.apply_mcy(op.qubits);
      break;
    case Gates::mcz:
      // Includes Z, CZ, CCZ, etc
      BaseState::qreg_.apply_mcphase(op.qubits, -1);
      break;
    case Gates::mcr:
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::r(op.params[0], op.params[1]));
      break;
    case Gates::mcrx:
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::rx(op.params[0]));
      break;
    case Gates::mcry:
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::ry(op.params[0]));
      break;
    case Gates::mcrz:
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::rz(op.params[0]));
      break;
    case Gates::rxx:
      BaseState::qreg_.apply_matrix(op.qubits, Linalg::VMatrix::rxx(op.params[0]));
      break;
    case Gates::ryy:
      BaseState::qreg_.apply_matrix(op.qubits, Linalg::VMatrix::ryy(op.params[0]));
      break;
    case Gates::rzz:
      BaseState::qreg_.apply_diagonal_matrix(op.qubits, Linalg::VMatrix::rzz_diag(op.params[0]));
      break;
    case Gates::rzx:
      BaseState::qreg_.apply_matrix(op.qubits, Linalg::VMatrix::rzx(op.params[0]));
      break;
    case Gates::id:
      break;
    case Gates::h:
      apply_gate_mcu(op.qubits, M_PI / 2., 0., M_PI, 0.);
      break;
    case Gates::s:
      apply_gate_phase(op.qubits[0], complex_t(0., 1.));
      break;
    case Gates::sdg:
      apply_gate_phase(op.qubits[0], complex_t(0., -1.));
      break;
    case Gates::t: {
      const double isqrt2{1. / std::sqrt(2)};
      apply_gate_phase(op.qubits[0], complex_t(isqrt2, isqrt2));
    } break;
    case Gates::tdg: {
      const double isqrt2{1. / std::sqrt(2)};
      apply_gate_phase(op.qubits[0], complex_t(isqrt2, -isqrt2));
    } break;
    case Gates::mcswap:
      // Includes SWAP, CSWAP, etc
      BaseState::qreg_.apply_mcswap(op.qubits);
      break;
    case Gates::mcu3:
      // Includes u3, cu3, etc
      apply_gate_mcu(op.qubits, std::real(op.params[0]), std::real(op.params[1]),
                     std::real(op.params[2]), 0.);
      break;
    case Gates::mcu:
      // Includes u3, cu3, etc
      apply_gate_mcu(op.qubits, std::real(op.params[0]), std::real(op.params[1]),
                      std::real(op.params[2]), std::real(op.params[3]));
      break;
    case Gates::mcu2:
      // Includes u2, cu2, etc
      apply_gate_mcu(op.qubits, M_PI / 2., std::real(op.params[0]),
                     std::real(op.params[1]), 0.);
      break;
    case Gates::mcp:
      // Includes u1, cu1, p, cp, mcp etc
      BaseState::qreg_.apply_mcphase(op.qubits,
                                     std::exp(complex_t(0, 1) * op.params[0]));
      break;
    case Gates::mcsx:
      // Includes sx, csx, mcsx etc
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::SX);
      break;
    case Gates::mcsxdg:
      BaseState::qreg_.apply_mcu(op.qubits, Linalg::VMatrix::SXDG);
      break;
    case Gates::pauli:
        BaseState::qreg_.apply_pauli(op.qubits, op.string_params[0]);
        break;
    default:
      // We shouldn't reach here unless there is a bug in gateset
      throw std::invalid_argument(
          "QubitVector::State::invalid gate instruction \'" + op.name + "\'.");
  }
}

template <class statevec_t>
void State<statevec_t>::apply_batched_ops(const std::vector<Operations::Op> &ops) 
{
  if(BaseState::qreg_.batched_optimization_supported()){
    std::vector<QV::batched_matrix_params> params;
    reg_t qubits;
    cvector_t matrices;
    int_t i,j,nqubit;
    QV::batched_matrix_params param;

    params.reserve(ops.size());
    qubits.reserve(ops.size());
    matrices.reserve(ops.size()*4);

    for(i=0;i<ops.size();i++){
      param.state_index_ = i;
      param.num_qubits_ = 1;
      param.offset_qubits_ = qubits.size();
      param.offset_matrix_ = matrices.size();
      if(ops[i].type == Operations::OpType::gate){
        auto it = State<statevec_t>::gateset_.find(ops[i].name);
        if (it == State<statevec_t>::gateset_.end())
          throw std::invalid_argument("QubitVectorState::invalid gate instruction \'" + 
                                      ops[i].name + "\'.");
        if(it->second != Gates::id && it->second != Gates::pauli){
          if(it->second == Gates::mcswap){
            // Includes SWAP, CSWAP, etc
            param.num_qubits_ = 2;
            param.set_control_mask(ops[i].qubits);
            qubits.push_back(ops[i].qubits[ops[i].qubits.size()-2]);
            qubits.push_back(ops[i].qubits[ops[i].qubits.size()-1]);
            reg_t qubits_sorted(2);
            qubits_sorted[0] = ops[i].qubits[ops[i].qubits.size()-2];
            qubits_sorted[1] = ops[i].qubits[ops[i].qubits.size()-1];
            std::sort(qubits_sorted.begin(), qubits_sorted.end());
            qubits.insert(qubits.end(),qubits_sorted.begin(),qubits_sorted.end());
            params.push_back(param);
            auto mat = Linalg::VMatrix::SWAP;
            matrices.insert(matrices.end(),mat.begin(),mat.end());
          }
          else{
            param.set_control_mask(ops[i].qubits);
            param.qubit_ = ops[i].qubits[ops[i].qubits.size()-1];

            switch (it -> second) {
              case Gates::mcx: {
                // Includes X, CX, CCX, etc
                auto mat = Linalg::VMatrix::X;
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcy: {
                // Includes Y, CY, CCY, etc
                auto mat = Linalg::VMatrix::Y;
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcz: {
                // Includes Z, CZ, CCZ, etc
                auto mat = Linalg::VMatrix::Z;
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcr: {
                auto mat = Linalg::VMatrix::r(ops[i].params[0], ops[i].params[1]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcrx: {
                auto mat = Linalg::VMatrix::rx(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcry: {
                auto mat = Linalg::VMatrix::ry(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcrz: {
                auto mat = Linalg::VMatrix::rz(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::rxx: {
                auto mat = Linalg::VMatrix::rxx(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::ryy: {
                auto mat = Linalg::VMatrix::ryy(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::rzz: {
                auto mat = Linalg::VMatrix::rzz(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::rzx: {
                auto mat = Linalg::VMatrix::rzx(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::h: {
                auto mat = Linalg::VMatrix::H;
                param.set_2x2matrix(mat);
              } break;
              case Gates::s: {
                auto mat = Linalg::VMatrix::S;
                param.set_2x2matrix(mat);
              } break;
              case Gates::sdg: {
                auto mat = Linalg::VMatrix::SDG;
                param.set_2x2matrix(mat);
              } break;
              case Gates::t: {
                auto mat = Linalg::VMatrix::T;
                param.set_2x2matrix(mat);
              } break;
              case Gates::tdg: {
                auto mat = Linalg::VMatrix::TDG;
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcu3: {
                // Includes u3, cu3, etc
                auto mat = Linalg::VMatrix::u3(ops[i].params[0],ops[i].params[1],ops[i].params[2]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcu2: {
                // Includes u2, cu2, etc
                auto mat = Linalg::VMatrix::u2(ops[i].params[0],ops[i].params[1]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcp: {
                // Includes u1, cu1, p, cp, mcp etc
                auto mat = Linalg::VMatrix::phase(ops[i].params[0]);
                param.set_2x2matrix(mat);
              } break;
              case Gates::mcsx: {
                // Includes sx, csx, mcsx etc
                auto mat = Linalg::VMatrix::SX;
                param.set_2x2matrix(mat);
              } break;
              default:
                // We shouldn't reach here unless there is a bug in gateset
                throw std::invalid_argument("QubitVector::State::invalid gate instruction \'" +
                                            ops[i].name + "\'.");
            }
            params.push_back(param);
          }
        }
      }
      else if(ops[i].type == Operations::OpType::matrix){
        param.num_qubits_ = ops[i].qubits.size();
        auto mat = Utils::vectorize_matrix(ops[i].mats[0]);
        if(ops[i].qubits.size() > 1){
          auto qubits_sorted = ops[i].qubits;
          std::sort(qubits_sorted.begin(), qubits_sorted.end());
          qubits.insert(qubits.end(),ops[i].qubits.begin(),ops[i].qubits.end());
          qubits.insert(qubits.end(),qubits_sorted.begin(),qubits_sorted.end());
          matrices.insert(matrices.end(),mat.begin(),mat.end());
        }
        else{
          param.qubit_ = ops[i].qubits[0];
          param.set_2x2matrix(mat);
        }
        param.control_mask_ = 0;
        params.push_back(param);
      }
      else if(ops[i].type == Operations::OpType::diagonal_matrix){
        param.num_qubits_ = ops[i].qubits.size();
        uint_t size = 1ull << ops[i].qubits.size();
        cvector_t mat(size*size,0.0);
        for(j=0;j<size;j++){
          mat[j*(size+1)] = ops[i].params[j];
        }
        if(ops[i].qubits.size() > 1){
          auto qubits_sorted = ops[i].qubits;
          std::sort(qubits_sorted.begin(), qubits_sorted.end());
          qubits.insert(qubits.end(),ops[i].qubits.begin(),ops[i].qubits.end());
          qubits.insert(qubits.end(),qubits_sorted.begin(),qubits_sorted.end());
          matrices.insert(matrices.end(),mat.begin(),mat.end());
        }
        else{
          param.qubit_ = ops[i].qubits[0];
          param.set_2x2matrix(mat);
        }
        param.control_mask_ = 0;
        params.push_back(param);
      }
    }

    BaseState::qreg_.apply_batched_matrix(params,qubits,matrices);
  }
}

template <class statevec_t>
void State<statevec_t>::apply_batched_pauli(std::vector<std::vector<Operations::Op>>& ops)
{
  BaseState::qreg_.apply_batched_pauli(ops);
}

template <class statevec_t>
void State<statevec_t>::apply_batched_noise_ops(const std::vector<std::vector<Operations::Op>> &ops, ExperimentResult &result,
                                                        std::vector<RngEngine> &rng)
{
  int_t i,j,k,count,nop,pos = 0;
  count = ops.size();

  reg_t mask(count);
  std::vector<bool> finished(count,false);
  for(i=0;i<count;i++){
    int_t cond_reg = -1;

    if(finished[i])
      continue;
    if(ops[i].size() == 0 || (ops[i].size() == 1 && ops[i][0].name == "id")){
      finished[i] = true;
      continue;
    }
    mask[i] = 1;

    //find same ops to be exectuted in a batch
    for(j=i+1;j<count;j++){
      if(finished[j]){
        mask[j] = 0;
        continue;
      }
      if(ops[j].size() == 0 || (ops[j].size() == 1 && ops[j][0].name == "id")){
        mask[j] = 0;
        finished[j] = true;
        continue;
      }

      if(ops[i].size() != ops[j].size()){
        mask[j] = 0;
        continue;
      }

      mask[j] = true;
      for(k=0;k<ops[i].size();k++){
        if(ops[i][k].conditional){
          cond_reg = ops[i][k].conditional_reg;
        }
        if(ops[i][k].type != ops[j][k].type || ops[i][k].name != ops[j][k].name){
          mask[j] = false;
          break;
        }
      }
      if(mask[j])
        finished[j] = true;
    }

    //mask conditional register
    int_t sys_reg = BaseState::qreg_.set_batched_system_conditional(cond_reg,mask);
    //batched execution on same ops
    for(k=0;k<ops[i].size();k++){
      Operations::Op cop = ops[i][k];

      //mark op conditional to mask shots
      cop.conditional = true;
      cop.conditional_reg = sys_reg;

      apply_op_multi_shots(cop,result,rng,false);
    }
    mask[i] = 0;
    finished[i] = true;
  }
}

template <class statevec_t>
void State<statevec_t>::apply_multiplexer(const reg_t &control_qubits,
                                          const reg_t &target_qubits,
                                          const cmatrix_t &mat) {
  if (control_qubits.empty() == false && target_qubits.empty() == false &&
      mat.size() > 0) {
    cvector_t vmat = Utils::vectorize_matrix(mat);
    BaseState::qreg_.apply_multiplexer(control_qubits, target_qubits, vmat);
  }
}

template <class statevec_t>
void State<statevec_t>::apply_matrix(const Operations::Op &op) {
  if (op.qubits.empty() == false && op.mats[0].size() > 0) {
    if (Utils::is_diagonal(op.mats[0], .0)) {
      BaseState::qreg_.apply_diagonal_matrix(
          op.qubits, Utils::matrix_diagonal(op.mats[0]));
    } else {
      BaseState::qreg_.apply_matrix(op.qubits,
                                    Utils::vectorize_matrix(op.mats[0]));
    }
  }
}

template <class statevec_t>
void State<statevec_t>::apply_matrix(const reg_t &qubits,
                                     const cvector_t &vmat) {
  // Check if diagonal matrix
  if (vmat.size() == 1ULL << qubits.size()) {
    BaseState::qreg_.apply_diagonal_matrix(qubits, vmat);
  } else {
    BaseState::qreg_.apply_matrix(qubits, vmat);
  }
}

template <class statevec_t>
void State<statevec_t>::apply_gate_mcu(const reg_t &qubits, double theta,
                                       double phi, double lambda, double gamma) {
  BaseState::qreg_.apply_mcu(qubits, Linalg::VMatrix::u4(theta, phi, lambda, gamma));
}

template <class statevec_t>
void State<statevec_t>::apply_gate_phase(uint_t qubit, complex_t phase) {
  cvector_t diag = {{1., phase}};
  apply_matrix(reg_t({qubit}), diag);
}

//=========================================================================
// Implementation: Reset, Initialize and Measurement Sampling
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_measure(const reg_t &qubits, const reg_t &cmemory,
                                      const reg_t &cregister, RngEngine &rng) 
{
  // Actual measurement outcome
  const auto meas = sample_measure_with_prob(qubits, rng);
  // Implement measurement update
  measure_reset_update(qubits, meas.first, meas.first, meas.second);
  const reg_t outcome = Utils::int2reg(meas.first, 2, qubits.size());
  BaseState::creg_.store_measure(outcome, cmemory, cregister);
}

template <class statevec_t>
rvector_t State<statevec_t>::measure_probs(const reg_t &qubits) const {
  return BaseState::qreg_.probabilities(qubits);
}

template <class statevec_t>
void State<statevec_t>::apply_reset(const reg_t &qubits, RngEngine &rng) {
  // Simulate unobserved measurement
  const auto meas = sample_measure_with_prob(qubits, rng);
  // Apply update to reset state
  measure_reset_update(qubits, 0, meas.first, meas.second);
}

template <class statevec_t>
std::pair<uint_t, double>
State<statevec_t>::sample_measure_with_prob(const reg_t &qubits,
                                            RngEngine &rng) {
  rvector_t probs = measure_probs(qubits);
  // Randomly pick outcome and return pair
  uint_t outcome = rng.rand_int(probs);
  return std::make_pair(outcome, probs[outcome]);
}

template <class statevec_t>
void State<statevec_t>::measure_reset_update(const std::vector<uint_t> &qubits,
                                             const uint_t final_state,
                                             const uint_t meas_state,
                                             const double meas_prob) {
  // Update a state vector based on an outcome pair [m, p] from
  // sample_measure_with_prob function, and a desired post-measurement
  // final_state

  // Single-qubit case
  if (qubits.size() == 1) {
    // Diagonal matrix for projecting and renormalizing to measurement outcome
    cvector_t mdiag(2, 0.);
    mdiag[meas_state] = 1. / std::sqrt(meas_prob);
    apply_matrix(qubits, mdiag);

    // If it doesn't agree with the reset state update
    if (final_state != meas_state) {
      BaseState::qreg_.apply_mcx(qubits);
    }
  }
  // Multi qubit case
  else {
    // Diagonal matrix for projecting and renormalizing to measurement outcome
    const size_t dim = 1ULL << qubits.size();
    cvector_t mdiag(dim, 0.);
    mdiag[meas_state] = 1. / std::sqrt(meas_prob);
    apply_matrix(qubits, mdiag);

    // If it doesn't agree with the reset state update
    // This function could be optimized as a permutation update
    if (final_state != meas_state) {
      // build vectorized permutation matrix
      cvector_t perm(dim * dim, 0.);
      perm[final_state * dim + meas_state] = 1.;
      perm[meas_state * dim + final_state] = 1.;
      for (size_t j = 0; j < dim; j++) {
        if (j != final_state && j != meas_state)
          perm[j * dim + j] = 1.;
      }
      // apply permutation to swap state
      apply_matrix(qubits, perm);
    }
  }
}

template <class statevec_t>
std::vector<reg_t> State<statevec_t>::sample_measure(const reg_t &qubits,
                                                     uint_t shots,
                                                     RngEngine &rng) {
  // Generate flat register for storing
  std::vector<double> rnds;
  rnds.reserve(shots);
  for (uint_t i = 0; i < shots; ++i)
    rnds.push_back(rng.rand(0, 1));

  auto allbit_samples = BaseState::qreg_.sample_measure(rnds);

  // Convert to reg_t format
  std::vector<reg_t> all_samples;
  all_samples.reserve(shots);
  for (int_t val : allbit_samples) {
    reg_t allbit_sample = Utils::int2reg(val, 2, BaseState::qreg_.num_qubits());
    reg_t sample;
    sample.reserve(qubits.size());
    for (uint_t qubit : qubits) {
      sample.push_back(allbit_sample[qubit]);
    }
    all_samples.push_back(sample);
  }
  return all_samples;
}

template <class statevec_t>
reg_t State<statevec_t>::local_sample_measure(const reg_t &qubits,
                                            std::vector<double>& rnds)
{
  reg_t samples = BaseState::qreg_.sample_measure(rnds);
  if(BaseState::state_index_ != 0){
    for(int_t i=0;i<samples.size();i++){
      samples[i] += (BaseState::state_index_ << BaseState::num_state_qubits_);
    }
  }

  return samples;
}

template <class statevec_t>
double State<statevec_t>::sum(void)
{
  return BaseState::qreg_.norm();
}

template <class statevec_t>
reg_t State<statevec_t>::batched_sample_measure(const reg_t &qubits,reg_t& shots,std::vector<RngEngine> &rng)
{
  uint_t i,j,pos,n = shots.size();
  uint_t total_shots = 0;

  for(i=0;i<n;i++){
    total_shots += shots[i];
  }
  std::vector<double> rnds(total_shots);

  //offset random values
  pos = 0;
  for(i=0;i<n;i++){
    for(j=0;j<shots[i];j++){
      rnds[pos++] = (double)i + rng[i].rand(0, 1);
    }
  }

  auto samples = BaseState::qreg_.sample_measure(rnds);

  //offset sample values
  pos = 0;
  for(i=0;i<n;i++){
    for(j=0;j<shots[i];j++){
      samples[pos++] -= (i << BaseState::qreg_.num_qubits());
    }
  }

  //return raw sample value, convert in multi-states class
  return samples;
}

template <class statevec_t>
void State<statevec_t>::store_measured_cbits(void)
{
  if(BaseState::qreg_.batched_optimization_supported()){
    uint_t i;
    reg_t pos(1);
    reg_t dummy_pos;

    for(i=0;i<BaseState::creg_.memory_size();i++){
      int bit = BaseState::qreg_.measured_cmemory(i);
      if(bit >= 0){
        const reg_t outcome = Utils::int2reg(bit, 2, 1);
        pos[0] = i;
        BaseState::creg_.store_measure(outcome, pos , dummy_pos);
      }
    }
    for(i=0;i<BaseState::creg_.register_size();i++){
      int bit = BaseState::qreg_.measured_cregister(i);
      if(bit >= 0){
        const reg_t outcome = Utils::int2reg(bit, 2, 1);
        pos[0] = i;
        BaseState::creg_.store_measure(outcome, dummy_pos, pos);
      }
    }
  }
}


template <class statevec_t>
void State<statevec_t>::apply_initialize(const reg_t &qubits,
                                         const cvector_t &params,
                                         RngEngine &rng) {

  if (qubits.size() == BaseState::qreg_.num_qubits()) {
    // If qubits is all ordered qubits in the statevector
    // we can just initialize the whole state directly
    auto sorted_qubits = qubits;
    std::sort(sorted_qubits.begin(), sorted_qubits.end());
    if (qubits == sorted_qubits) {
      initialize_qreg(qubits.size(), params);
      return;
    }
  }
  // Apply reset to qubits
  apply_reset(qubits, rng);
  // Apply initialize_component
  BaseState::qreg_.initialize_component(qubits, params);
}

//=========================================================================
// Implementation: Multiplexer Circuit
//=========================================================================

template <class statevec_t>
void State<statevec_t>::apply_multiplexer(const reg_t &control_qubits,
                                          const reg_t &target_qubits,
                                          const std::vector<cmatrix_t> &mmat) {
  // (1) Pack vector of matrices into single (stacked) matrix ... note: matrix
  // dims: rows = DIM[qubit.size()] columns = DIM[|target bits|]
  cmatrix_t multiplexer_matrix = Utils::stacked_matrix(mmat);

  // (2) Treat as single, large(r), chained/batched matrix operator
  apply_multiplexer(control_qubits, target_qubits, multiplexer_matrix);
}

//=========================================================================
// Implementation: Kraus Noise
//=========================================================================
template <class statevec_t>
void State<statevec_t>::apply_kraus(const reg_t &qubits,
                                    const std::vector<cmatrix_t> &kmats,
                                    RngEngine &rng) 
{
  // Check edge case for empty Kraus set (this shouldn't happen)
  if (kmats.empty())
    return; // end function early

//  if(BaseState::qreg_.batched_optimization_supported()){
//    BaseState::qreg_.apply_kraus(qubits,kmats,rng);
//  }
  // Choose a real in [0, 1) to choose the applied kraus operator once
  // the accumulated probability is greater than r.
  // We know that the Kraus noise must be normalized
  // So we only compute probabilities for the first N-1 kraus operators
  // and infer the probability of the last one from 1 - sum of the previous

  double r = rng.rand(0., 1.);
  double accum = 0.;
  bool complete = false;

  // Loop through N-1 kraus operators
  for (size_t j = 0; j < kmats.size() - 1; j++) {

    // Calculate probability
    cvector_t vmat = Utils::vectorize_matrix(kmats[j]);
    double p = BaseState::qreg_.norm(qubits, vmat);
    accum += p;

    // check if we need to apply this operator
    if (accum > r) {
      // rescale vmat so projection is normalized
      Utils::scalar_multiply_inplace(vmat, 1 / std::sqrt(p));
      // apply Kraus projection operator
      apply_matrix(qubits, vmat);
      complete = true;
      break;
    }
  }

  // check if we haven't applied a kraus operator yet
  if (complete == false) {
    // Compute probability from accumulated
    complex_t renorm = 1 / std::sqrt(1. - accum);
    apply_matrix(qubits, Utils::vectorize_matrix(renorm * kmats.back()));
  }
}

//-------------------------------------------------------------------------
} // namespace Statevector
//-------------------------------------------------------------------------
} // end namespace AER
//-------------------------------------------------------------------------
#endif

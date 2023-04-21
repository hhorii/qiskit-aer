/**
 * This code is part of Qiskit.
 *
 * (C) Copyright IBM 2018, 2019.2023.
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */

#ifndef _shot_branching_hpp
#define _shot_branching_hpp

namespace AER {

namespace CircuitExecutor {

using OpItr = std::vector<Operations::Op>::const_iterator;

class Branch;

// class for shared state for sho-branching
class Branch {
protected:
  uint_t state_index_; // state index

  uint_t shot_index_; // starting shot index

  // creg to be stored to the state
  ClassicalRegister creg_;
  // random generators for shots
  std::vector<RngEngine> shots_;
  // additional operations applied after shot branching
  std::vector<Operations::Op> additional_ops_;

  // mark for control flow
  std::unordered_map<std::string, OpItr> flow_marks_;

  // current iterator of operations
  OpItr iter_;

  // branches from this
  std::vector<std::shared_ptr<Branch>> branches_;

public:
  Branch(void) {}
  ~Branch() {
    shots_.clear();
    additional_ops_.clear();
    branches_.clear();
  }
  Branch(const Branch &src) {
    shots_ = src.shots_;
    creg_ = src.creg_;
    iter_ = src.iter_;
    flow_marks_ = src.flow_marks_;
  }

  uint_t &state_index(void) { return state_index_; }
  uint_t &shot_index(void) { return shot_index_; }
  ClassicalRegister &creg(void) { return creg_; }
  std::vector<RngEngine> &rng_shots(void) { return shots_; }
  OpItr &op_iterator(void) { return iter_; }
  std::unordered_map<std::string, OpItr> &marks(void) { return flow_marks_; }
  uint_t num_branches(void) { return branches_.size(); }
  std::vector<std::shared_ptr<Branch>> &branches(void) { return branches_; }

  uint_t num_shots(void) { return shots_.size(); }
  void clear(void) {
    shots_.clear();
    additional_ops_.clear();
    branches_.clear();
  }
  void clear_branch(void) { branches_.clear(); }

  void set_shots(std::vector<RngEngine> &shots) { shots_ = shots; }
  void initialize_shots(const uint_t nshots, const uint_t seed) {
    shots_.resize(nshots);
    for (int_t i = 0; i < nshots; i++) {
      shots_[i].set_seed(seed + i);
    }
  }

  void add_op_after_branch(Operations::Op &op) {
    additional_ops_.push_back(op);
  }
  void copy_ops_after_branch(std::vector<Operations::Op> &ops) {
    additional_ops_ = ops;
  }
  void clear_additional_ops(void) { additional_ops_.clear(); }

  std::vector<Operations::Op> &additional_ops(void) { return additional_ops_; }

  void branch_shots(reg_t &shots, int_t nbranch);

  bool apply_control_flow(ClassicalRegister &creg, OpItr last) {
    if (iter_->type == Operations::OpType::mark) {
      flow_marks_[iter_->string_params[0]] = iter_;
      iter_++;
      return true;
    } else if (iter_->type == Operations::OpType::jump) {
      if (creg.check_conditional(*iter_)) {
        const auto &mark_name = iter_->string_params[0];
        auto mark_it = flow_marks_.find(mark_name);
        if (mark_it != flow_marks_.end()) {
          iter_ = mark_it->second;
        } else {
          for (++iter_; iter_ != last; ++iter_) {
            if (iter_->type == Operations::OpType::mark) {
              flow_marks_[iter_->string_params[0]] = iter_;
              if (iter_->string_params[0] == mark_name) {
                break;
              }
            }
          }
          if (iter_ == last) {
            std::stringstream msg;
            msg << "Invalid jump destination:\"" << mark_name << "\"."
                << std::endl;
            throw std::runtime_error(msg.str());
          }
        }
      }
      iter_++;
      return true;
    }
    return false;
  }

  void advance_iterator(void);

  bool apply_runtime_noise_sampling(const ClassicalRegister &creg,
                                    const Operations::Op &op,
                                    const Noise::NoiseModel &noise);

  void remove_empty_branches(void);
};

void Branch::branch_shots(reg_t &shots, int_t nbranch) {
  branches_.resize(nbranch);

  for (int_t i = 0; i < nbranch; i++) {
    branches_[i] = std::make_shared<Branch>();
    branches_[i]->creg_ = creg_;
    branches_[i]->iter_ = iter_;
    branches_[i]->flow_marks_ = flow_marks_;
  }
  for (int_t i = 0; i < shots.size(); i++) {
    branches_[shots[i]]->shots_.push_back(shots_[i]);
  }
  // update shot indices
  uint_t index = shot_index_;
  for (int_t i = 0; i < nbranch; i++) {
    branches_[i]->shot_index_ = index;
    index += branches_[i]->shots_.size();
  }
}

void Branch::advance_iterator(void) {
  iter_++;
  for (int_t i = 0; i < branches_.size(); i++) {
    branches_[i]->iter_++;
  }
}

bool Branch::apply_runtime_noise_sampling(const ClassicalRegister &creg,
                                          const Operations::Op &op,
                                          const Noise::NoiseModel &noise) {
  if (op.type != Operations::OpType::sample_noise)
    return false;

  uint_t nshots = num_shots();
  reg_t shot_map(nshots);
  std::vector<std::vector<Operations::Op>> noises;

  for (int_t i = 0; i < nshots; i++) {
    std::vector<Operations::Op> noise_ops =
        noise.sample_noise_loc(op, shots_[i]);

    // search same noise ops
    int_t pos = -1;
    for (int_t j = 0; j < noises.size(); j++) {
      if (noise_ops.size() != noises[j].size())
        continue;
      bool same = true;
      for (int_t k = 0; k < noise_ops.size(); k++) {
        if (noise_ops[k].type != noises[j][k].type ||
            noise_ops[k].name != noises[j][k].name)
          same = false;
        else if (noise_ops[k].qubits.size() != noises[j][k].qubits.size())
          same = false;
        else {
          for (int_t l = 0; l < noise_ops[k].qubits.size(); l++) {
            if (noise_ops[k].qubits[l] != noises[j][k].qubits[l]) {
              same = false;
              break;
            }
          }
        }
        if (!same)
          break;
        if (noise_ops[k].type == Operations::OpType::gate) {
          if (noise_ops[k].name == "pauli") {
            if (noise_ops[k].string_params[0] != noises[j][k].string_params[0])
              same = false;
          } else if (noise_ops[k].params.size() != noises[j][k].params.size())
            same = false;
          else {
            for (int_t l = 0; l < noise_ops[k].params.size(); l++) {
              if (noise_ops[k].params[l] != noises[j][k].params[l]) {
                same = false;
                break;
              }
            }
          }
        } else if (noise_ops[k].type == Operations::OpType::matrix ||
                   noise_ops[k].type == Operations::OpType::diagonal_matrix) {
          if (noise_ops[k].mats.size() != noises[j][k].mats.size())
            same = false;
          else {
            for (int_t l = 0; l < noise_ops[k].mats.size(); l++) {
              if (noise_ops[k].mats[l].size() != noises[j][k].mats[l].size()) {
                same = false;
                break;
              }
              for (int_t m = 0; m < noise_ops[k].mats[l].size(); m++) {
                if (noise_ops[k].mats[l][m] != noises[j][k].mats[l][m]) {
                  same = false;
                  break;
                }
              }
              if (!same)
                break;
            }
          }
        }
        if (!same)
          break;
      }
      if (same) {
        pos = j;
        break;
      }
    }

    if (pos < 0) { // if not found, add noise ops to the list
      shot_map[i] = noises.size();
      noises.push_back(noise_ops);
    } else { // if found, add shot
      shot_map[i] = pos;
    }
  }

  creg_ = creg;
  branch_shots(shot_map, noises.size());
  for (int_t i = 0; i < noises.size(); i++) {
    branches_[i]->copy_ops_after_branch(noises[i]);
  }

  return true;
}

void Branch::remove_empty_branches(void) {
  int_t istart = 0;
  for (int_t j = 0; j < branches_.size(); j++) {
    if (branches_[j]->num_shots() > 0) {
      // copy shots to the root
      shots_ = branches_[j]->rng_shots();
      additional_ops_ = branches_[j]->additional_ops();
      shot_index_ = branches_[j]->shot_index();
      creg_ = branches_[j]->creg();
      branches_[j].reset();
      istart = j + 1;
      break;
    }
    branches_[j].reset();
  }

  std::vector<std::shared_ptr<Branch>> new_branches;

  for (int_t j = istart; j < branches_.size(); j++) {
    if (branches_[j]->num_shots() > 0)
      new_branches.push_back(branches_[j]);
    else
      branches_[j].reset();
  }
  branches_ = new_branches;
}

//-------------------------------------------------------------------------
} // namespace CircuitExecutor
//-------------------------------------------------------------------------
} // end namespace AER
//-------------------------------------------------------------------------
#endif

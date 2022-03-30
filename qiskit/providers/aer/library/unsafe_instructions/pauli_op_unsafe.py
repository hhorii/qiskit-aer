# This code is part of Qiskit.
#
# (C) Copyright IBM 2018, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""
Simulator instruction to save exact operator expectation value.
"""

from numpy import allclose
from qiskit.circuit import Instruction
from qiskit.quantum_info import Pauli, SparsePauliOp, Operator
from qiskit.circuit import QuantumCircuit
from qiskit.extensions.exceptions import ExtensionError

from ..save_instructions.save_expectation_value import _expval_params

class PauliOperatorUnsafe(Instruction):
    """Apply pauli operator."""
    def __init__(self,
                 operator):
        r"""Instruction to apply a Hermitian operator.

        An Hermitian operator :math:`H` is applied to quantum state.
        This instruction does not guarantee that the applied state is valid.

        Args:
            operator (Pauli or SparsePauliOp or Operator): a Hermitian operator.

        Raises:
            ExtensionError: if the input operator is invalid or not Hermitian.

        """
        # Convert O to SparsePauliOp representation
        if isinstance(operator, Pauli):
            operator = SparsePauliOp(operator)
        elif not isinstance(operator, SparsePauliOp):
            operator = SparsePauliOp.from_operator(Operator(operator))
        if not allclose(operator.coeffs.imag, 0):
            raise ExtensionError("Input operator is not Hermitian.")
        params = _expval_params(operator, variance=False)
        super().__init__("pauli_op_unsafe", operator.num_qubits, 0, params)

def pauli_op_unsafe(self,
                    operator,
                    qubits):
    r"""Apply a Hermitian operator.

    Args:
        operator (Pauli or SparsePauliOp or Operator): a Hermitian operator.
        qubits (list): circuit qubits to apply instruction.

    Returns:
        QuantumCircuit: with attached instruction.

    Raises:
        ExtensionError: if the input operator is invalid or not Hermitian.

    .. note::

        This method appends a :class:`SaveExpectationValue` instruction to
        the quantum circuit.
    """
    instr = PauliOperatorUnsafe(operator)
    return self.append(instr, qubits)

QuantumCircuit.pauli_op_unsafe = pauli_op_unsafe

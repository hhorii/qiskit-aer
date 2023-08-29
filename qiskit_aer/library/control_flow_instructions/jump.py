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
Simulator instruction to set a program counter
"""

from qiskit.circuit import Instruction
from qiskit.circuit.classical.expr import Expr


class AerJump(Instruction):
    """
    Jump instruction

    This instruction sets a program counter to specified mark instruction.
    """

    _directive = True

    def __init__(self, jump_to, num_qubits, num_clbits=0):
        super().__init__("jump", num_qubits, num_clbits, [jump_to])
        self.condition_expr = None

    def set_conditional(self, cond):
        if isinstance(cond, Expr):
            self.c_if_expr(cond)
        else:
            self.c_if(*cond)
        return self

    def c_if_expr(self, cond):
        self.condition_expr = cond
        return self

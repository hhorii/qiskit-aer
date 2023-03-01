# This code is part of Qiskit.
#
# (C) Copyright IBM 2018, 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""
AerSimulator Integration Tests
"""
from math import sqrt
from ddt import ddt
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit import CircuitInstruction
from test.terra.reference import ref_algorithms

from test.terra.backends.simulator_test_case import (
    SimulatorTestCase, supported_methods)


@ddt
class TestVariousCircuit(SimulatorTestCase):
    """AerSimulator tests to simulate various types of circuits"""

    @supported_methods(
        ['automatic', 'statevector', 'density_matrix',
         'matrix_product_state', 'extended_stabilizer', 'tensor_network'])
    def test_quantum_register_circuit(self, method, device):
        """Test circuits with quantum registers."""
        
        qubits = QuantumRegister(3)
        clbits = ClassicalRegister(3)

        circuit = QuantumCircuit(qubits, clbits)
        circuit.h(qubits[0])
        circuit.cx(qubits[0], qubits[1])
        circuit.cx(qubits[0], qubits[2])

        for q, c in zip(qubits, clbits):
            circuit.measure(q, c)

        backend = self.backend(
            method=method,
            device=device,
            seed_simulator=1111
        )

        shots=1000
        result = backend.run(circuit, shots=shots).result()
        self.assertSuccess(result)
        self.compare_counts(result, [circuit], [{'0x0':500, '0x7':500}], delta=0.05 * shots)
        
    @supported_methods(
        ['automatic', 'statevector', 'density_matrix',
         'matrix_product_state', 'extended_stabilizer', 'tensor_network'])
    def test_qubits_circuit(self, method, device):
        """Test circuits with quantum registers."""
        
        qubits = QuantumRegister(3)
        clbits = ClassicalRegister(3)

        circuit = QuantumCircuit()
        circuit.add_bits(qubits)
        circuit.add_bits(clbits)
        circuit.h(qubits[0])
        circuit.cx(qubits[0], qubits[1])
        circuit.cx(qubits[0], qubits[2])

        for q, c in zip(qubits, clbits):
            circuit.measure(q, c)

        backend = self.backend(
            method=method,
            device=device,
            seed_simulator=1111
        )

        shots=1000
        result = backend.run(circuit, shots=shots).result()
        self.assertSuccess(result)
        self.compare_counts(result, [circuit], [{'0x0':500, '0x7':500}], delta=0.05 * shots)

    @supported_methods(
        ['automatic', 'statevector', 'density_matrix',
         'matrix_product_state', 'extended_stabilizer', 'tensor_network'])
    def test_qubits_quantum_register_circuit(self, method, device):
        """Test circuits with quantum registers."""
        
        qubits0 = QuantumRegister(2)
        clbits1 = ClassicalRegister(2)
        qubits1 = QuantumRegister(1)
        clbits2 = ClassicalRegister(1)

        circuit = QuantumCircuit(qubits0, clbits1)
        circuit.add_bits(qubits1)
        circuit.add_bits(clbits2)
        circuit.h(qubits0[0])
        circuit.cx(qubits0[0], qubits0[1])
        circuit.cx(qubits0[0], qubits1[0])

        for qubits, clbits in zip([qubits0, qubits1], [clbits1, clbits2]):
            for q, c in zip(qubits, clbits):
                circuit.measure(q, c)

        backend = self.backend(
            method=method,
            device=device,
            seed_simulator=1111
        )

        shots=1000
        result = backend.run(circuit, shots=shots).result()
        self.assertSuccess(result)
        self.compare_counts(result, [circuit], [{'0x0':500, '0x7':500}], delta=0.05 * shots)

        qubits0 = QuantumRegister(1)
        clbits1 = ClassicalRegister(1)
        qubits1 = QuantumRegister(1)
        clbits2 = ClassicalRegister(1)
        qubits2 = QuantumRegister(1)
        clbits3 = ClassicalRegister(1)

        circuit = QuantumCircuit(qubits0, clbits1)
        circuit.add_bits(qubits1)
        circuit.add_bits(clbits2)
        circuit.add_register(qubits2)
        circuit.add_register(clbits3)
        circuit.h(qubits0[0])
        circuit.cx(qubits0[0], qubits1[0])
        circuit.cx(qubits1[0], qubits2[0])

        for qubits, clbits in zip([qubits0, qubits1, qubits2], [clbits1, clbits2, clbits3]):
            for q, c in zip(qubits, clbits):
                circuit.measure(q, c)
 
        backend = self.backend(
            method=method,
            device=device,
            seed_simulator=1111
        )

        shots=1000
        result = backend.run(circuit, shots=shots).result()
        self.assertSuccess(result)
        self.compare_counts(result, [circuit], [{'0x0':500, '0x7':500}], delta=0.05 * shots)

# This code is part of Qiskit.
#
# (C) Copyright IBM 2018, 2019, 2020, 2021, 2022.
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
from qiskit import transpile, QuantumCircuit
from test.terra.reference import ref_algorithms

from test.terra.backends.simulator_test_case import (
    SimulatorTestCase, supported_methods)


@ddt
class TestMetadata(SimulatorTestCase):
    """AerSimulator algorithm tests in the default basis"""

    @supported_methods(
        ['automatic', 'statevector', 'density_matrix',
         'matrix_product_state', 'extended_stabilizer'])
    def test_single_circuit_metadata(self, method, device):
        """Test circuits with object metadata."""
        backend = self.backend(method=method, device=device)
        metadata = {1: object}
        circuit = QuantumCircuit(1, name='circ0', metadata=metadata)
        result = backend.run(circuit).result()
        self.assertSuccess(result)
        self.assertIsNone(result.results[0].header.metadata)
        self.assertEqual(circuit.metadata, {1: object})

    @supported_methods(
        ['automatic', 'statevector', 'density_matrix',
         'matrix_product_state', 'extended_stabilizer'])
    def test_three_circuit_metadata(self, method, device):
        """Test circuits with object metadata."""
        backend = self.backend(method=method, device=device)

        metadata0 = {0: object}
        circuit0 = QuantumCircuit(1, name='circ0', metadata=metadata0)
        
        metadata1 = {1: object}
        circuit1 = QuantumCircuit(1, name='circ1', metadata=metadata1)

        metadata2 = {2: object}
        circuit2 = QuantumCircuit(1, name='circ2', metadata=metadata2)
        
        result = backend.run([circuit0, circuit1, circuit2]).result()
        self.assertSuccess(result)
        self.assertEqual(len(result.results), 3)
        self.assertIsNone(result.results[0].header.metadata)
        self.assertIsNone(result.results[1].header.metadata)
        self.assertIsNone(result.results[2].header.metadata)
        self.assertEqual(circuit0.metadata, {0: object})
        self.assertEqual(circuit1.metadata, {1: object})
        self.assertEqual(circuit2.metadata, {2: object})

    @supported_methods(
        ['automatic', 'statevector', 'density_matrix', 'matrix_product_state'])
    def test_three_parameterized_circuit_metadata(self, method, device):
        """Test circuits with object metadata."""
        backend = self.backend(method=method, device=device)

        metadata0 = {0: object}
        circuit0 = QuantumCircuit(1, name='circ0', metadata=metadata0)
        circuit0.ry(0.1, 0)
        circuit0.measure_all()
        
        metadata1 = {1: object}
        circuit1 = QuantumCircuit(1, name='circ1', metadata=metadata1)
        circuit1.ry(0.1, 0)
        circuit1.measure_all()

        metadata2 = {2: object}
        circuit2 = QuantumCircuit(1, name='circ2', metadata=metadata2)
        circuit2.ry(0.1, 0)
        circuit2.measure_all()
        
        parameterizations=[[[[0, 0], [0, 1]]],
                           [[[0, 0], [0, 1, 2]]],
                           []]
        
        result = backend.run([circuit0, circuit1, circuit2],
                             parameterizations=parameterizations).result()
        self.assertSuccess(result)
        self.assertEqual(len(result.results), 6)
        self.assertIsNone(result.results[0].header.metadata)
        self.assertIsNone(result.results[1].header.metadata)
        self.assertIsNone(result.results[2].header.metadata)
        self.assertIsNone(result.results[3].header.metadata)
        self.assertIsNone(result.results[4].header.metadata)
        self.assertIsNone(result.results[5].header.metadata)
        self.assertEqual(circuit0.metadata, {0: object})
        self.assertEqual(circuit1.metadata, {1: object})
        self.assertEqual(circuit2.metadata, {2: object})
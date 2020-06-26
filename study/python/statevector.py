# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# This is a simplified version of qasm_simulator.py
# https://github.com/Qiskit/qiskit-terra/blob/stable/0.14/qiskit/providers/basicaer/qasm_simulator.py

import numpy as np

class StateVector():
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits;
        self.psi = np.zeros(2 ** self.num_qubits)
        self.psi[0] = 1
        self.psi = np.reshape(self.psi, self.num_qubits * [2])
        self.rng = np.random.RandomState()

    def random(self):
        return self.rng.rand()

    def unitary_2x2 (self, n, unitary):
        
        if unitary.shape != (2, 2):
            raise Exception('invalid shape: {0}'.format(unitary.shape))
        
        index = 'z'
        index += chr(ord('a') + (self.num_qubits - n - 1))
        index += ","
        for i in range(self.num_qubits):
            index += chr(ord('a') + i)
        index += '->'
        for i in range(self.num_qubits):
            if (self.num_qubits - n - 1) == i:
                index += 'z'
            else:
                index += chr(ord('a') + i)
        self.psi = np.einsum(index, unitary, self.psi)
    
    def u3(self, n, theta, phi, lam):
        mat = np.array([[np.cos(theta / 2), -np.exp(1j * lam) * np.sin(theta / 2)],
                        [np.exp(1j * phi) * np.sin(theta / 2), np.exp(1j * phi + 1j * lam) * np.cos(theta / 2)]])
        self.unitary_2x2(n, mat)
    
    def u2(self, n, phi, lam):
        self.u3(n, np.pi/2, phi, lam)
    
    def u1(self, n, lam):
        self.u3(n, 0, 0, lam)
    
    def x(self, n):
        self.u3(n, np.pi, 0, np.pi)
    
    def y(self, n):
        self.u3(n, np.pi, np.pi/2, np.pi/2)
    
    def z(self, n):
        self.u1(n, np.pi)
    
    def z(self, n):
        self.u1(n, np.pi)
        
    def h(self, n):
        self.u2(n, 0, np.pi)
    
    def s(self, n):
        self.u1(n, np.pi / 2)
    
    def sdg(self, n):
        self.u1(n, -np.pi / 2)
    
    def t(self, n):
        self.u1(n, np.pi / 4)
    
    def tdg(self, n):
        self.u1(n, -np.pi / 4)
    
    def rx(self, n, theta):
        self.u3(n, theta, -np.pi/2, np.pi/2)
    
    def ry(self, n, theta):
        self.u3(n, theta, 0, 0)
    
    def rz(self, n, phi):
        self.u1(n, phi)
    
    def unitary_4x4 (self, n, m, unitary):
        
        if unitary.shape != (2, 2, 2, 2):
            raise Exception('invalid shape: {0}'.format(unitary.shape))
        
        index = 'zy'
        index += chr(ord('a') + (self.num_qubits - m - 1))
        index += chr(ord('a') + (self.num_qubits - n - 1))
        index += ","
        for i in range(self.num_qubits):
            index += chr(ord('a') + i)
        index += '->'
        for i in range(self.num_qubits):
            if (self.num_qubits - n - 1) == i:
                index += 'y'
            elif (self.num_qubits - m - 1) == i:
                index += 'z'
            else:
                index += chr(ord('a') + i)
        self.psi = np.einsum(index, unitary, self.psi)
    
    def cnot (self, control, target):
        self.unitary_4x4(target, control, np.array([[[[1, 0], [0, 0]],
                                               [[0, 1], [0, 0]]], 
                                               [[[0, 0], [0, 1]], 
                                               [[0, 0], [1, 0]]]]))
    
    def measure(self, n):
        axis = list(range(self.num_qubits))
        axis.remove(self.num_qubits - 1 - n)
        probabilities = np.sum(np.abs(self.psi) ** 2, axis=tuple(axis))
        p0 = self.random()
        if p0 < probabilities[0]:
            self.unitary_2x2(n, np.array([[1 / np.sqrt(probabilities[0]), 0], [0, 0]]))
            return 0
        else:
            self.unitary_2x2(n, np.array([[0, 1 / np.sqrt(probabilities[1])], [0, 0]]))
            return 1

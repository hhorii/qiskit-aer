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
QasmSimulator Integration Tests
"""

import unittest
from test.terra import common
from test.terra.backends.qasm_simulator.qasm_method import QasmMethodTests
from test.terra.backends.qasm_simulator.qasm_measure import QasmMeasureTests
from test.terra.backends.qasm_simulator.qasm_reset import QasmResetTests
from test.terra.backends.qasm_simulator.qasm_initialize import QasmInitializeTests
from test.terra.backends.qasm_simulator.qasm_conditional import QasmConditionalGateTests
from test.terra.backends.qasm_simulator.qasm_conditional import QasmConditionalUnitaryTests
from test.terra.backends.qasm_simulator.qasm_conditional import QasmConditionalKrausTests
from test.terra.backends.qasm_simulator.qasm_cliffords import QasmCliffordTests
from test.terra.backends.qasm_simulator.qasm_cliffords import QasmCliffordTestsWaltzBasis
from test.terra.backends.qasm_simulator.qasm_cliffords import QasmCliffordTestsMinimalBasis
from test.terra.backends.qasm_simulator.qasm_noncliffords import QasmNonCliffordTests
from test.terra.backends.qasm_simulator.qasm_noncliffords import QasmNonCliffordTestsWaltzBasis
from test.terra.backends.qasm_simulator.qasm_noncliffords import QasmNonCliffordTestsMinimalBasis
from test.terra.backends.qasm_simulator.qasm_algorithms import QasmAlgorithmTests
from test.terra.backends.qasm_simulator.qasm_algorithms import QasmAlgorithmTestsWaltzBasis
from test.terra.backends.qasm_simulator.qasm_algorithms import QasmAlgorithmTestsMinimalBasis
from test.terra.backends.qasm_simulator.qasm_extra import QasmExtraTests
from test.terra.backends.qasm_simulator.qasm_thread_management import QasmThreadManagementTests
from test.terra.backends.qasm_simulator.qasm_fusion import QasmFusionTests
from test.terra.backends.qasm_simulator.qasm_truncate import QasmQubitsTruncateTests
from test.terra.backends.qasm_simulator.qasm_basics import QasmBasicsTests


class TestQasmSimulator(common.QiskitAerTestCase,
                        QasmMethodTests,
                        QasmMeasureTests,
                        QasmResetTests,
                        QasmInitializeTests,
                        QasmConditionalGateTests,
                        QasmConditionalUnitaryTests,
                        QasmConditionalKrausTests,
                        QasmCliffordTests,
                        QasmCliffordTestsWaltzBasis,
                        QasmCliffordTestsMinimalBasis,
                        QasmNonCliffordTests,
                        QasmNonCliffordTestsWaltzBasis,
                        QasmNonCliffordTestsMinimalBasis,
                        QasmAlgorithmTests,
                        QasmAlgorithmTestsWaltzBasis,
                        QasmAlgorithmTestsMinimalBasis,
                        QasmExtraTests,
                        QasmThreadManagementTests,
                        QasmFusionTests,
                        QasmQubitsTruncateTests,
                        QasmBasicsTests):
    """QasmSimulator automatic method tests."""


if __name__ == '__main__':
    unittest.main()

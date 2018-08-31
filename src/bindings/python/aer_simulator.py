"""
Cython quantum circuit simulator.
"""

import sys
import os
import json
import logging
import warnings
import numpy as np

# Import QISKit classes
from qiskit._result import Result
from qiskit.backends import BaseBackend
from qiskit.backends.local.localjob import LocalJob

# Import Simulator tools
from helpers import SimulatorJSONEncoder, qobj2schema
from aer_qv_wrapper import AerSimulatorWrapper

# Logger
logger = logging.getLogger(__name__)


class AerSimulator(BaseBackend):
    """Cython quantum circuit simulator"""

    DEFAULT_CONFIGURATION = {
        'name': 'local_aer_simulator',
        'url': 'NA',
        'simulator': True,
        'local': True,
        'description': 'A C++ statevector simulator for qobj files',
        'coupling_map': 'all-to-all',
        "basis_gates": 'u0,u1,u2,u3,cx,cz,id,x,y,z,h,s,sdg,t,tdg,rzz'
    }

    def __init__(self, configuration=None):
        super().__init__(configuration or self.DEFAULT_CONFIGURATION.copy())
        self.simulator = AerSimulatorWrapper()

    def run(self, qobj):
        """Run a QOBJ on the the backend."""
        return LocalJob(self._run_job, qobj)

    def _run_job(self, qobj):
        self._validate(qobj)
        qobj_str = json.dumps(qobj2schema(qobj), cls=SimulatorJSONEncoder)
        result = json.loads(self.simulator.execute(qobj_str))
        # TODO: get schema in line with result object
        return result  # Result(result)

    def load_noise_model(self, noise_model):
        self.simulator.load_noise_model(json.dumps(noise_model, cls=SimulatorJSONEncoder))

    def clear_noise_model(self):
        self.simulator.clear_noise_model()

    def load_config(self, engine=None, state=None):
        if engine is not None:
            self.simulator.load_engine_config(json.dumps(engine, cls=SimulatorJSONEncoder))
        if state is not None:
            self.simulator.load_state_config(json.dumps(state, cls=SimulatorJSONEncoder))

    def set_max_threads_shot(self, threads):
        """
        Set the maximum threads used for parallel shot execution.

        Args:
            threads (int): the thread limit, set to -1 to use maximum available

        Note that using parallel shot evaluation disables parallel circuit
        evaluation.
        """

        self.simulator.set_max_threads_shot(int(threads))

    def set_max_threads_circuit(self, threads):
        """
        Set the maximum threads used for parallel circuit execution.

        Args:
            threads (int): the thread limit, set to -1 to use maximum available

        Note that using parallel circuit evaluation disables parallel shot
        evaluation.
        """
        self.simulator.set_max_threads_circuit(int(threads))

    def set_max_threads_state(self, threads):
        """
        Set the maximum threads used for state update parallel  routines.

        Args:
            threads (int): the thread limit, set to -1 to use maximum available.

        Note that using parallel circuit or shot execution takes precidence over
        parallel state evaluation.
        """
        self.simulator.set_max_threads_state(int(threads))

    def _validate(self, qobj):
        # TODO
        return

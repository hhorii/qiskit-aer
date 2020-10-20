from benchmark.simulator_benchmark import SimulatorBenchmarkSuite

DEFAULT_APPS = {
#    'fourier_checking',
#    'graph_state',
#    'hidden_linear_function',
#    'iqp',
    'quantum_volume': 1,
#    'phase_estimation'
    }

DEFAULT_QUBITS = [ 25 ]

DEFAULT_RUNTIME = [
    SimulatorBenchmarkSuite.RUNTIME_STATEVECTOR_CPU,
    SimulatorBenchmarkSuite.RUNTIME_STATEVECTOR_GPU,
    SimulatorBenchmarkSuite.RUNTIME_MPS_CPU,
    ]

DEFAULT_MEASUREMENT_COUNTS = [ 1, 10, 100, 1000, 10000 ]

DEFAULT_NOISE_MODELS = [
    SimulatorBenchmarkSuite.NOISE_IDEAL
    ]

class Sampling(SimulatorBenchmarkSuite):

    def __init__(self,
                 apps = DEFAULT_APPS,
                 qubits = DEFAULT_QUBITS,
                 runtime_names = DEFAULT_RUNTIME,
                 measures = [SimulatorBenchmarkSuite.MEASUREMENT_SAMPLING],
                 measure_counts = DEFAULT_MEASUREMENT_COUNTS,
                 noise_model_names = DEFAULT_NOISE_MODELS):
        super().__init__( 'sampling', 
                          apps, 
                          qubits=qubits, 
                          runtime_names=runtime_names, 
                          measures=measures, 
                          measure_counts=measure_counts, 
                          noise_model_names=noise_model_names)
        self.__name__ = 'sampling'

class ExpVal(SimulatorBenchmarkSuite):

    def __init__(self,
                 apps = DEFAULT_APPS,
                 qubits = DEFAULT_QUBITS,
                 runtime_names = DEFAULT_RUNTIME,
                 measures = [SimulatorBenchmarkSuite.MEASUREMENT_EXPVAL],
                 measure_counts = DEFAULT_MEASUREMENT_COUNTS,
                 noise_model_names = DEFAULT_NOISE_MODELS):
        super().__init__( 'expval',
                          apps, 
                          qubits=qubits, 
                          runtime_names=runtime_names, 
                          measures=measures, 
                          measure_counts=measure_counts, 
                          noise_model_names=noise_model_names)
        self.__name__ = 'expval'

    
if __name__ == "__main__":
    Sampling().run_manual()
    ExpVal().run_manual()
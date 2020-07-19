# simplest quantum computing simulators for study of qiskit-aer functionalities


## cpp

move to `cpp/` and compile codes to run variations of qft.
Following shell commands generate: 
- `qft`: single-thread version
- `qft_mt` multi-thread version
- `qft_mt_simd`: multi-thread SIMD version
- `qft_thrust_omp`: thrust-omp version
- `qft_thrust`: thrust-cuda version 


```
g++ -O3 -std=c++11 -o qft qft.cpp -ffast-math

g++ -O3 -std=c++11 -o qft_mt qft_mt.cpp -ffast-math -fopenmp

g++ -O3 -I. -std=c++11 -o qft_simd_mt qft_mt_simd.cpp -mavx2 -mfma -fopenmp

git clone --depth 1 --branch 1.9.5 https://github.com/thrust/thrust.git thrust
g++ -O3 -I. -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -std=c++11 -o qft_thrust_omp qft_thrust.cpp -ffast-math -fopenmp

nvcc -O3 -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA -std=c++11 -o qft_thrust qft_thrust.cu

g++ -c -I . -I /usr/local/opt/openblas/include/ -L /usr/local/opt/openblas/lib/ -lopenblas -ffast-math -fopenmp -O3 -DNDEBUG -fopenmp -std=gnu++17 qft_mps.cpp
g++ -o qft_mps -I . -I /usr/local/opt/openblas/include/ -L /usr/local/opt/openblas/lib/ -lopenblas -ffast-math -fopenmp -O3 -DNDEBUG -fopenmp -std=gnu++17 qft_mps.o libitensor-g.a


```

## python

mote to `python/` and run `qft.py`.


import sys
import numpy as np
import statevector

if len(sys.argv) < 3:
    exit()

N = int (sys.argv[1])
shots = int (sys.argv[2])

results = {}

for shot in range(shots):
    
    sv = statevector.StateVector(N)
    
    for i in range(N):
        sv.h(i);
    
    for i in range(N):
        for j in range(i):
            l = np.pi / 2 ** (i - j)
            sv.u1(i, l/2)
            sv.cnot(i, j)
            sv.u1(j, -l/2)
            sv.cnot(i, j)
            sv.u1(j, l/2)
        sv.h(i);
    
    result = ''
    for i in range(N):
        result = '{0}{1}'.format(result, sv.measure(N - i - 1))
        
    if result in results:
        results[result] += 1
    else:
        results[result] = 1
        
print (results)
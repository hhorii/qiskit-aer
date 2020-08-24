import sys
import numpy as np
import statevector

if len(sys.argv) < 2:
    exit()

N = int (sys.argv[1])

sv = statevector.StateVector(N, unitary_matrix=True)

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

result = sv.get_unitary_matrix()

for i in range(len(result)):
    for j in range(len(result[i])):
        print((round(result[i][j].real, 2) + round(result[i][j].imag, 2)), end=',')
    print()

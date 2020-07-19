#include <map>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "statevector_mps.hpp"

using namespace std;

int main_bell(int argc, char* argv[]) {
//int main(int argc, char* argv[]) {
  QubitVectorMPS qv(2, true);
  qv.h(0);
  qv.cnot(0, 1);

  auto wf = qv.statevector();
  cout << "final state: [";
  for (auto& pa : wf) {
    cout << "[" << pa.real() << ", " << pa.imag() << "], ";
  }
  cout << "]" << endl;

  qv.dump_tensor();

  return 1;
}

int main_ghz(int argc, char* argv[]) {
//int main(int argc, char* argv[]) {
  QubitVectorMPS qv(3);
  qv.h(0);
  qv.cnot(0, 1);
  qv.cnot(1, 2);

  auto wf = qv.statevector();
  cout << "final state: [";
  for (auto& pa : wf) {
    cout << "[" << pa.real() << ", " << pa.imag() << "], ";
  }
  cout << "]" << endl;

  qv.dump_tensor();

  return 1;
}

int main_qv(int argc, char* argv[]) {
//int main(int argc, char* argv[]) {
  int qubits = 10;
  int depth = 10;

  complex<double> u[16];

  u[0] = complex<double>(-0.38422631, 2.83471478e-01);
  u[1] = complex<double>(0.25943494, 7.65221683e-02);
  u[2] = complex<double>(-0.56281274, -3.91037228e-02);
  u[3] = complex<double>(0.54010984, -2.98070487e-01);
  u[4] = complex<double>(-0.34951286, 5.51982645e-01);
  u[5] = complex<double>(-0.1755857 , 1.77050098e-01);
  u[6] = complex<double>(-0.38094516, -5.26882445e-01);
  u[7] = complex<double>(-0.28000149, 9.92657401e-02);
  u[8] = complex<double>(-0.01016025, -4.69180869e-01);
  u[9] = complex<double>(0.05238197-4.75213843e-01);
  u[10] = complex<double>(-0.44501419, -5.84109815e-03);
  u[11] = complex<double>(-0.54594179, 2.34669613e-01);
  u[12] = complex<double>(-0.17806077, -3.05336581e-01);
  u[13] = complex<double>(0.1054005, 7.90556427e-01);
  u[14] = complex<double>(0.00897634, 2.46649688e-01);
  u[15] = complex<double>(-0.42196931, 6.80411917e-04);

  QubitVectorMPS qv(qubits);
  //qv.svdCutoff = 1e-2;

  std::srand(0);

  for (auto i = 0; i < depth; ++i) {
    int q0 = abs(rand() % qubits);
    int q1 = q0;
    while(q0 == q1) q1 = abs(rand() % qubits);

    cout << "SU " << q0 << ", " << q1 << endl;

    qv.unitary(q0, q1,
        u[0], u[1], u[2], u[3],
        u[4], u[5], u[6], u[7],
        u[8], u[9], u[10], u[11],
        u[12], u[13], u[14], u[15]
    );

    qv.dump_tensor();
  }

  return 1;
}

//int main_qft(int argc, char* argv[]) {
int main(int argc, char* argv[]) {

  if (argc < 3) {
    return -1;
  }

  int N = atoi(argv[1]);
  int shots = atoi(argv[2]);
  map<string, unsigned int> results;

  for(auto shot = 0; shot < shots; ++shot) {

    QubitVectorMPS qv(N);

    for (unsigned i =0; i < N; ++i) {
      qv.h(i);
    }
    for (unsigned i =0; i < N; ++i) {
      for (unsigned j =0; j < i; ++j) {
        double l = M_PI / double (1UL << (i - j));
        qv.u1(i, l/2.);
        qv.cnot(i, j);
        qv.u1(j, -l/2.);
        qv.cnot(i, j);
        qv.u1(j, l/2.);
      }
      qv.h(i);
    }
    stringstream ss;
    for (unsigned i = 0; i < N; ++i) {
      ss << qv.measure(N - i - 1);
    }
    string result = ss.str();
    auto itr = results.find(result);
    if (itr == results.end()) {
      results[result] = 1;
    } else {
      ++results[result];
    }
  }
  for(auto itr = results.begin(); itr != results.end(); ++itr) {
    std::cout << itr->first << " : " << itr->second << endl;
  }

  return 1;
}

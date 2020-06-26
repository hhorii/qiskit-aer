#include <map>
#include <sstream>
#include <iostream>
#include "statevector_mt.hpp"

int main(int argc, char* argv[]) {

  if (argc < 3) {
    return -1;
  }

  int N = atoi(argv[1]);
  int shots = atoi(argv[2]);
  map<string, unsigned int> results;

  for(auto shot = 0; shot < shots; ++shot) {

    QubitVectorMT qv(N);

    for (unsigned i =0; i < N; ++i)
      qv.h(i);
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
    for (unsigned i = 0; i < N; ++i)
      ss << qv.measure(N - i - 1);
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

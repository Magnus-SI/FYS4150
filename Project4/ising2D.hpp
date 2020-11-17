#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>

using namespace std;

class Ising {

protected:
  int m_L;
  double m_E, m_M, C_v, m_chi;
  int *m_spin;
  double *m_w;
  int periodic(int i, int add);

public:
  void seed(int s);
  void initialize(int L, double temp, double tol);
  void metropolis();

};
#endif

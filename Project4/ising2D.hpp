#ifndef ISING2D_HPP
#define ISING2D_HPP
#include <fstream>

using namespace std;

class ising2D {

protected:
  int m_L, m_T;
  double m_E, m_M, C_v, m_chi;
  int *m_spin;
  double *m_w;
  int periodic(int i);

public:
  void seed(int s);
  void initialize(int L, double temp, double tol);
  void metropolis();
  void mean_values();
  double *m_mean;

};
#endif

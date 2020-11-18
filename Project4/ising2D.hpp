#ifndef ISING2D_HPP
#define ISING2D_HPP
#include <fstream>

using namespace std;

class ising2D {

protected:
  int m_L, m_T;
  int m_mcs;
  double m_E, m_M, C_v, m_chi;
  int *m_spin;
  double *m_w;
  long idum;


public:
  void seed(int s);
  void initialize(int L, double temp, double tol);
  void metropolis();
  void mean_values();
  void write_to_file(std::ofstream&);
  int periodic(int i, int ixshift, int iyshift);
  //int periodic2(int i);
  double *m_mean;
  int m_accepted;

};
#endif

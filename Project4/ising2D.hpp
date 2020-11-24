#ifndef ISING2D_HPP
#define ISING2D_HPP
#include <fstream>

using namespace std;

class ising2D {

protected:
  int m_L, m_mcs;
  double m_M, m_E, C_v, m_chi, m_T;
  int *m_spin;
  double *m_w;


public:
  void seed(int s);
  void initialize(int L, double temp, double tol);
  void metropolis();
  void mean_values();
  void write_to_file(ofstream&);
  void write_mean(stringstream&);
  int periodic(int i, int ixshift, int iyshift);
  double *m_mean;
  int m_accepted, m_meancount, m_meanstart;
  //int m_0, m_1, m_2, m_3;

};
#endif

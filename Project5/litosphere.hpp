#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP
#include <fstream>

using namespace std;

class Diffusion {

protected:
  double m_alpha, m_dx, m_tc, m_dt;
  int m_nx, m_nL, m_nxb;
  double  *u, *uprev, *uexact, *Q;

public:
  void init(int nx, double dx, double dt, int impmet);
  void init2D(int nxy, double dxy, double dt);
  void EulerForward();
  void EF2D();
  void tridiag_solve();
  void CrankNicolson();
  void GetExact(int);
  void write_andiff(ofstream&, int);
  void write_1D(ofstream&);


};
#endif

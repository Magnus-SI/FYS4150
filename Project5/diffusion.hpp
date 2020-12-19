#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP
#include <fstream>

using namespace std;

class Diffusion {

protected:
  double m_alpha, m_dx, m_dt;
  int m_nx, m_ny, m_n2, m_Qtype;
  double  *u, *uprev, *uexact, *Q;
  double *a, *b, *c, *d; //for the implicit methods

public:
  double m_tc;
  void init(int nx, double dx, double dt, int impmet);
  void init2D(int nx, int ny, double dxy, double dt, int Qtype);
  void update_Q();
  void Q_0();
  void EulerForward();
  void EF2D();
  int xperiodic(int ix, int iy, int ixshift, int iyshift);
  void EF2Dperx();        //with periodic x conditions, used for.
  void tridiag_solve();
  void CrankNicolson();
  void GetExact(int);
  void write_andiff(ofstream&, int);
  double err2D(int);
  double err2Dlocal(int, int);
  void write_analytic(ofstream&, int);
  void write_1D(ofstream&);
  void write_Q(ofstream&);


};
#endif

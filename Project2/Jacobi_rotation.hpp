#ifndef JACOBI_ROTATION_HPP
#define JACOBI_ROTATION_HPP
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

class Jacobi_rotation {

protected:
  mat m_Hamiltonian;
  double m_h, m_rho;
  int m_N;
  vec m_q, m_x, m_v;
  ofstream m_ofile;
  int m_Vchoice;

public:
  mat A, V;
  double omega_r;
  int iters;
  double V_func(double rho);
  void initialize(double a, double b, int N, int Vchoice);
  void print_matrix();
  void rotate(double conv);
  void write_to_file(string filename);
  void rearrange();
  vec return_eig();
  mat test_eig();
  void eigarma(vec &eigval, mat &eigvec);
  float quanteigtest(int method);
};
#endif

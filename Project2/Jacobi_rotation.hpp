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

public:
  void initialize(double a, double b, int N, double V(double rho));
  void print_matrix();
  void rotate(double conv, vec r, mat v);
  void write_to_file(string filename);
};
#endif

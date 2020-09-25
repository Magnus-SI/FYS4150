#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
using namespace arma;

void Jacobi_rotation::initialize(double a, double b, int N, double V(double rho))
{
  /* Initializing the problem, the Hamiltonian matrix, etc. */
  m_N = N;
  m_h = (b - a)/(m_N+1);
  m_Hamiltonian = zeros<mat>(m_N, m_N);

  double  DiagConst, NondiagConst;
  DiagConst = 2.0 / (m_h*m_h);
  NondiagConst =  -1.0 / (m_h*m_h);
  vec x = linspace(a + m_h, b-m_h, m_N);

  // Setting up tridiagonal matrix and diagonalization using Armadillo
  m_Hamiltonian(0,0) = DiagConst;
  m_Hamiltonian(0,1) = NondiagConst;
  for(int i = 1; i < m_N-1; i++) {
    m_Hamiltonian(i,i-1)    = NondiagConst;
    m_Hamiltonian(i,i)    = DiagConst + V(x(i));
    m_Hamiltonian(i,i+1)    = NondiagConst;
  }
  m_Hamiltonian(m_N-1,m_N-2) = NondiagConst;
  m_Hamiltonian(m_N-1,m_N-1) = DiagConst;
}

void Jacobi_rotation::print_matrix()
{
  /* Just printing the matrix. Should only be used for testing on small matrices */
  m_Hamiltonian.print();
}

void Jacobi_rotation::rotate()
{
  /* Method to rotate matrix */
}

void Jacobi_rotation::write_to_file(string filename)
{
  /* Needs to be fixed */
  m_ofile.open(filename);
  for (int i = 0; i < m_N; i++){
    m_ofile << m_x(i) << " " << m_v(i) << endl;
  }
  m_ofile.close();
}

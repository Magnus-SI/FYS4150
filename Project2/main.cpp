#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
using namespace arma;

double V_0(double rho);


int main()
{
  //Defining interval
  double a = 0, b = 1;
  //Number of points
  int N = 10;
  Jacobi_rotation my_solver;
  //Initializing the matrix to diagonalise
  my_solver.initialize(a, b, N, V_0);
  //my_solver.print_matrix();

  double conv = 1e-8;
  vec r = zeros<vec>(N);
  mat V = zeros<mat>(N, N);
  my_solver.rotate(conv, r, V);

}

double V_0(double rho)
{
  /* Zero potential */
  return 0;
}

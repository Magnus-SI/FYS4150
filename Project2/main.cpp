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

  double conv = 1e-8;
  my_solver.rotate(conv);
  //my_solver.test_eig();
  my_solver.rearrange();
}

double V_0(double rho)
{
  /* Zero potential */
  return 0;
}

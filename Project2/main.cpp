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
  double a = 0, b = 1;
  int N = 10;
  Jacobi_rotation my_solver;
  my_solver.initialize(a, b, N, V_0);
  my_solver.print_matrix();
}

double V_0(double rho)
{
  return 0;
}

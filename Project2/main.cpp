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

  double rhomin = 0;
  int N = 10;
  Jacobi_rotation solver;
  for(int i=1; i<6; i++){
    cout<<i<<endl
    double rhomax = i;
    solver.initialize(rhomin, rhomax, N, V_d);
    solver.rotate(conv);
    solver.rearrange();
    string rhomstr = to_string(i);
    rhomstr.append(".csv")
    solver.write_to_file(rhomstr);
  }
}

double V_0(double rho)
{
  /* Zero potential */
  return 0;
}

double V_d(double rho)
{
  return pow(rho, 2);
}

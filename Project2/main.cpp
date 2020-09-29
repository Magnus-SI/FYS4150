#include "Jacobi_rotation.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>


using namespace arma;
using namespace std;

double V_0(double rho);

double V_d(double rho);


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
  my_solver.write_to_file("ref.csv");

  double rhomin = 0;
  int n = 10;
  Jacobi_rotation solver;

  vec rhomax = linspace(1, 10, 10);
  vec Nvals = logspace(1, 2, 3);

  ofstream ofile;
  ofile.open("e2d.csv");
  ofile <<setw(15) << setprecision(8);
  ofile << "rhomax," << "1," << "2," << "3" << endl;

  for (int i=0; i<10; i++){
    ofile<<rhomax(i)<<",";
    for (int j = 0; j<3; j++){
      solver.initialize(rhomin, rhomax(i), Nvals(j), V_d);
      solver.rotate(conv);
      solver.rearrange();
      double maxerr = solver.quanteigtest();
      ofile << maxerr<<",";
    ofile<<endl;
    }
  }
  ofile.close();

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

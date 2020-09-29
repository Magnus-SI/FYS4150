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

int main()
{
  //Defining interval
  double a = 0, b = 1;
  //Number of points
  int N = 10;
  Jacobi_rotation my_solver;
  //Initializing the matrix to diagonalise
  int Vchoice = 0;
  my_solver.initialize(a, b, N, Vchoice);

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

  // ofstream ofile;
  // ofile.open("e2d.csv");
  // ofile <<setw(15) << setprecision(8);
  // ofile << "rhomax," << "1," << "2," << "3," << endl;
  //
  // for (int i=0; i<10; i++){
  //   ofile<<rhomax(i)<<",";
  //   for (int j = 0; j<3; j++){
  //     solver.initialize(rhomin, rhomax(i), Nvals(j), V_d);
  //     solver.rotate(conv);
  //     solver.rearrange();
  //     double maxerr = solver.quanteigtest();
  //     ofile << maxerr<<",";
  //   }
  //   ofile<<endl;
  // }
  // ofile.close();
  Vchoice = 2;
  double rho_max = 4;
  int N_val = 100;
  Jacobi_rotation solver2;
  double omega_rs[] = {0.01, 0.5, 1, 5};
  for(int i =0; i<4; i++){
    solver2.omega_r = omega_rs[i];
    solver2.initialize(rhomin, rho_max, N_val, Vchoice);
    solver2.rotate(conv);
    solver2.rearrange();
    string fname = "el2omega";
    fname.append(to_string(i));
    solver2.write_to_file(fname.append(".csv"));
  }

}

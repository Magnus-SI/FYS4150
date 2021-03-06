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
  //Solve buckling beam problem and save the three lowest eigenstates to file.
  double a = 0, b = 1;
  int N = 100;          //number of points
  Jacobi_rotation my_solver;
  //Initializing the matrix to diagonalise
  int Vchoice = 0;    //buckling beam, no potential
  my_solver.initialize(a, b, N, Vchoice);
  double conv = 1e-8;   //Jacobi tolerance for off-diagonals
  my_solver.rotate(conv);
  my_solver.rearrange();
  my_solver.write_to_file("csv_files/beam.csv");

  //Find optimal rhomax and N for the 1 electron system
  double rhomin = 0;
  vec rhomax = linspace(1, 10, 10);
  vec Nvals = logspace(1, 2, 3);
  Jacobi_rotation solver;

  ofstream ofile;
  ofile.open("csv_files/e2d.csv");
  ofile <<setw(15) << setprecision(8);
  //ofile << "rhomax," << "1," << "2," << "3," <<"4,"<<"5,"<<"6,"<<"7,"<<"8,"<<"9,"<<"10,"<<"11,"<< endl;
  ofile << "rhomax," << "1," << "2," << "3" << endl;
  Vchoice = 1;      //choose potential for the 1 electron quantum situation
  int method = 0;   //method 0 if Jacobi, 1 if armadillo
  for (int i=0; i<10; i++){
    ofile<<rhomax(i)<<",";
    for (int j = 0; j<3; j++){
      solver.initialize(rhomin, rhomax(i), Nvals(j), Vchoice);
      solver.rotate(conv);
      solver.rearrange();
      double maxerr = solver.quanteigtest(method);
      ofile << maxerr<<",";
    }
    ofile<<endl;
  }
  ofile.close();


  // Vchoice = 1;
  // double rho_max = 5;
  // int N_val = 200;
  // Jacobi_rotation solver3;
  // solver3.initialize(rhomin, rho_max, N_val, Vchoice);
  // solver3.rotate(conv);
  // solver3.rearrange();
  // solver3.write_to_file("opt.csv");

  //Save data for the 2 electron system.
  Vchoice = 2;
  double rho_max = 6;
  int N_val = 100;
  Jacobi_rotation solver2;
  double omega_rs[] = {0.01, 0.5, 1, 5};
  for(int i =0; i<4; i++){
    solver2.omega_r = omega_rs[i];
    solver2.initialize(rhomin, rho_max, N_val, Vchoice);
    solver2.rotate(conv);
    solver2.rearrange();
    string fname = "csv_files/el2omega";
    fname.append(to_string(i));
    solver2.write_to_file(fname.append(".csv"));
  }

}

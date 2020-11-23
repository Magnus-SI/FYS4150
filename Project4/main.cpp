#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

void mcice(int L, double temp, double tol, int mcs_max){
  stringstream params;
  params << fixed << L << "_" << setprecision(2) << temp << "_" << tol;
  string filename = "data/mcdep_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << "mcs,E,M,acpt" << endl;
  ofile << setw(15) << setprecision(8);
  ising2D my_ising;
  my_ising.seed(-1);
  my_ising.initialize(L, temp, tol);
  my_ising.write_to_file(ofile);
  for(int mcs=0; mcs<mcs_max; mcs++){
    my_ising.metropolis();
    my_ising.write_to_file(ofile);
  }
  ofile.close();
}

int main(){


  //First run for 2x2 lattice, T = 1.0

  int L = 2; double temp = 1; double tol = 0.0; int mcs_1 = 1000000;
  mcice(L, temp, tol, mcs_1);

  //Then run 20x20 lattice, T = 1.0. Ordered and random initial config
  L = 20; double tol2 = 0.5;
  mcice(L, temp, tol, mcs_1);
  mcice(L, temp, tol2, mcs_1);

  //Then run 20x20 lattice, T = 2.4. Ordered and random initial config
  temp = 2.4;
  mcice(L, temp, tol, mcs_1);
  mcice(L, temp, tol2, mcs_1);

  return 0;
}

#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

int main(){


  //First test for how average values evolve with metropolis cycles
  string filename = "f1.txt";
  ofstream ofile;
  ofile.open(filename);
  ofile << "mc,E,M,Cv,chi" << endl;
  ofile <<setw(15) << setprecision(8);

  ising2D my_ising;
  int L = 2; double temp = 1.0; double tol = 0.0;
  my_ising.initialize(L, temp, tol);
  for(int mcs=1; mcs<1001; mcs++){
    my_ising.metropolis();
    ofile << mcs;
    my_ising.write_to_file(ofile);
    //cout << (double) my_ising.m_mean[0]/mcs << endl;
  }
  ofile.close();

  //Then look at how results depend on temperature

  return 0;
}

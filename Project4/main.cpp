#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

void mcice(int L, double temp, double tol){
  std::stringstream params;
  params << std::fixed << L << "_"<< std::setprecision(2) << temp <<"_" << tol;
  string filename = "data/mcdep_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << "mc,E,M,Cv,chi,acptfrac" << endl;
  ofile <<setw(15) << setprecision(8);
  ising2D my_ising;
  my_ising.initialize(L, temp, tol);
  for(int mcs=1; mcs<10011; mcs++){
    my_ising.metropolis();
    ofile << mcs << ",";
    my_ising.write_to_file(ofile);
    //cout << (double) my_ising.m_mean[0]/mcs << endl;
    //cout << my_ising.m_accepted << mcs << endl;
  }
  ofile.close();
}

int main(){


  //First test for how average values evolve with metropolis cycles
  //this could perhaps be put into a test function

  int L = 2; double temp = 1.0; double tol = 0.5;
  mcice(L, temp, tol);


  ising2D my_ising;
  L = 10;
  my_ising.initialize(L, temp, tol);
  cout << my_ising.periodic(0, 0, 1) << " " << my_ising.periodic(39, 1, 0) << " " << my_ising.periodic(60, -1, 0) << " " << my_ising.periodic(99, 0, 0) << " " << my_ising.periodic(99, 1, 0) << " " << my_ising.periodic(4, 0, -1) << " " << my_ising.periodic(99, 0, 1) <<endl;


  /*
  //Now use study further temp and tol-dependence for L = 20
  L = 20;
  double temps[2] = {1.0, 2.4};
  double tols[2] = {0.0, 0.5};
  for (int i = 0; i<2; i++){
    for (int j = 0; j<2; j++){
      mcice(L, temps[i], tols[j]);
    }
  }

  //Then look at how results depend on temperature
  int Tcount = 7;
  double temps2[Tcount] = {2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3};
  L = 40;

  ofstream ofile;
  string filename = "data/tempvars.csv";
  ofile.open(filename);
  ofile << "T,E,M,Cv,chi,acptfrac" << endl;
  ofile <<setw(15) << setprecision(8);

  ising2D my_ising;
  for (int i = 0; i<Tcount; i++){
    my_ising.initialize(L, temps2[i], tol);
    for(int mcs=1; mcs<4001; mcs++){
      my_ising.metropolis();
    }
    ofile << temps2[i] << ",";
    my_ising.write_to_file(ofile);
  }
  ofile.close();
  */


  return 0;
}

#include "solar_system.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

int main(int argc, char *argv[]){
  int N = atoi(argv[1]);
  int Nt = atoi(argv[2]);
  double beta = atof(argv[3]);
  double T = 1e5*Nt;

  string filename = "solar";
  std::stringstream params;
  params << std::fixed << N << "_"<< std::setprecision(2) << beta <<"_" <<log10(Nt);
  filename.append(params.str()).append(".txt");

  solar_system solar_solver;
  solar_solver.initialize(N, Nt, T, beta);
  solar_solver.remove_drift();
  solar_solver.F_G(0);
  for(int l=0; l<Nt-1; l++){
    solar_solver.velocity_verlet(l);
  }
  solar_solver.write_to_file(filename);
  solar_solver.initialize_mercury_sun(1000000, 1e4*Nt, beta);
  solar_solver.remove_drift();
  solar_solver.F_G_corrected(0);
  for(int k=0; k<Nt-1; k++){
    solar_solver.mercury(k);
  }
  solar_solver.write_to_file("mercury.txt");
}

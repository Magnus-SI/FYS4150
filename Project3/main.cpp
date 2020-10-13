#include "solar_system.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

int main(){
  int N, Nt;
  N = 3;
  Nt = 100000;
  double T = 1e6;
  string filename = "solar.txt";
  solar_system solar_solver;
  solar_solver.initialize(N, Nt, T);
  solar_solver.remove_drift();
  solar_solver.F_G(0);
  for(int l=0; l<Nt-1; l++){
    solar_solver.velocity_verlet(l);
  }
  solar_solver.write_to_file(filename);
}

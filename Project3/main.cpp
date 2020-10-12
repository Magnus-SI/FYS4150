#include "solar_system.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

int main(){
  int N, Nt;
  N = 10;
  Nt = 10000;
  double T = 1e7;
  string filename = "solar.txt";
  solar_system solar_solver;
  solar_solver.initialize(N, Nt, T);
  solar_solver.F_G(0);
  for(int l=0; l<Nt; l++){
    int m = 10*l;
    solar_solver.velocity_verlet(m);
  }
  solar_solver.write_to_file(filename);
}

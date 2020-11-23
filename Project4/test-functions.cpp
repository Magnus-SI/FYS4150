#include "catch.hpp"
#include "ising2D.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "time.h"


TEST_CASE("yes"){
  int a = 1;
}


TEST_CASE("Project 3 example"){
  int Nt = 20000;
  double beta = 2;
  double T = 1e9;
  double N = 10;
  double totE;
  solar_system solar_solver;
  solar_solver.initialize(N, Nt, T, beta, 0);
  solar_solver.remove_drift();
  solar_solver.F_G(0);
  double* E0 = solar_solver.total_energy(0);
  double totE0 = E0[0] - E0[1];

  for (int i = 0; i<Nt-1; i++){
    solar_solver.velocity_verlet(i);
    if(remainder(i,100) == 0){
      double* E = solar_solver.total_energy(i);
      totE = E[0] - E[1];
      REQUIRE(totE/totE0 == Approx(1).epsilon(1e-3));

    }
  }
}

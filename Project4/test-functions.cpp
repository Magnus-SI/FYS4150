#include "catch.hpp"
#include "ising2D.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "time.h"

double Z(double temp, double J), E_func(double temp, double J),
  absM_func(double temp, double J), chi_func(double temp, double J),
  Cv_func(double temp, double J);

TEST_CASE("Project 4 analytic"){
  int mcs_max = 100000;
  double temp = 1.0;
  double tol = 0.0;
  double E, E2, absM, M2;
  ising2D test_ising;
  test_ising.seed(-1);
  test_ising.initialize(2, temp, tol);
  for(int mcs=0; mcs<mcs_max; mcs++){
    test_ising.metropolis();
  }
  int meanc = test_ising.m_meancount;
  E = test_ising.m_mean[0]/meanc; E2 = test_ising.m_mean[1]/meanc;
  absM = test_ising.m_mean[4]/meanc; M2 = test_ising.m_mean[3]/meanc;
  double C_v, chi;
  C_v = E2 - E*E;
  chi = M2 - absM*absM;
  REQUIRE(E/E_func(temp, 1) == Approx(1).epsilon(1e-2));
  REQUIRE(absM/absM_func(temp, 1) == Approx(1).epsilon(1e-2));
  REQUIRE(C_v/Cv_func(temp, 1) == Approx(1).epsilon(1e-2));
  REQUIRE(chi/chi_func(temp, 1) == Approx(1).epsilon(1e-2));
}

double Z(double temp, double J){
  return 4*cosh(8*J/temp) + 12;
}

double E_func(double temp, double J){
  return -32*J/Z(temp, J)*sinh(8*J/temp);
}

double absM_func(double temp, double J){
  return 8/Z(temp, J)*(exp(8*J/temp) + 2);
}

double Cv_func(double temp, double J){
  return 1/pow(temp, 2)*pow(32*J/Z(temp, J), 2)*(3*cosh(8*J/temp) + 1);
}

double chi_func(double temp, double J){
  return 64/(temp*pow(Z(temp, J), 2))*(3*exp(8*J/temp) + exp(-8*J/temp) + 3);
}

#include "catch.hpp"
#include "solar_system.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "time.h"

TEST_CASE("EulerVerlet"){
  int Nt;
  double beta = 2;
  double T = 3e8;   //about 10 years
  int elliptical = 0;

  ofstream ofile;
  ofile.open("timeacc.csv");
  ofile <<setw(15) << setprecision(8);
  ofile << "log10N," << "eultime," << "vertime," << "eulerr," << "vererr" <<endl;
  solar_system verlet_solver;
  solar_system euler_solver;
  double time1, time2, time3, euler_time, verlet_time, euler_err, verlet_err;

  for (int i = 3; i<7; i++){
    Nt = pow(10, i);
    verlet_solver.initialize_earth_sun(Nt, T, beta, elliptical);
    //verlet_solver.remove_drift();
    verlet_solver.F_G(0);
    euler_solver.initialize_earth_sun(Nt, T, beta, elliptical);
    //euler_solver.remove_drift();
    euler_solver.F_G(0);

    time1 = clock();
    for(int m=0; m<Nt-1; m++){
      euler_solver.forward_euler(m);
    }
    time2 = clock();
    for(int m=0; m<Nt-1; m++){
      verlet_solver.velocity_verlet(m);
    }
    time3 = clock();
    euler_time = (time2 - time1)/CLOCKS_PER_SEC;
    verlet_time = (time3 - time2)/CLOCKS_PER_SEC;

    double* cqeuler = euler_solver.conserved_quants(Nt-2);
    double* cqverlet = verlet_solver.conserved_quants(Nt-2);

    euler_err = cqeuler[1];
    verlet_err = cqverlet[1];

    ofile << i << "," << euler_time << verlet_time << euler_err << verlet_err <<endl;

    //verlet_solver.write_to_file("yo.txt");

  }
  ofile.close();


  // for (int i = 0; i<3; i++){
  //   cout << cqeuler[i] << cqverlet[i] << endl;
  //   REQUIRE(cqverlet[i] == Approx(1).epsilon(1e-3));
  //   REQUIRE(cqeuler[i] == Approx(1).epsilon(1e-3));
  //
  // }

}

TEST_CASE("Angular Momentum Conservation"){
  int Nt = 10000;
  double beta = 2;
  double T = 3e8;   //about 10 years
  int elliptical = 1;
  solar_system earth_sun;
  earth_sun.initialize_earth_sun(Nt, T, beta, elliptical);
  earth_sun.F_G(0);
  double dA;
  double dA0 = earth_sun.Kep2Area(0, 1);
  for (int i = 0; i<Nt-1; i++){
    earth_sun.velocity_verlet(i);
    if(remainder(i,100) == 0){
      dA = earth_sun.Kep2Area(i, 1);
      REQUIRE(dA/dA0 == Approx(1).epsilon(1e-2));
    }
  }

}

#include "solar_system.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

int main(int argc, char *argv[]){
  //int N = atoi(argv[1]);
  //int Nt = atoi(argv[2]);
  //double beta = atof(argv[3]);
  //double T = 1e5*Nt;
  int N, Nt;
  double beta, T;

  //First explore different values of beta with a fixed sun
  solar_system earth_sun;
  Nt = 100000;
  T = 12e9;
  double betas[6] = {2,2.01,2.1, 2.2, 2.5, 3};
  for (int i = 0; i<6; i++){
    earth_sun.initialize_earth_sun(Nt, T, betas[i], 0);
    earth_sun.F_G(0);
    for (int j = 0; j<Nt-1; j++){
      earth_sun.velocity_verlet(j);
    }
    earth_sun.write_to_file("data/earth_sun");
  }

  //Then explore different values of beta with elliptical orbit
  solar_system elliptical_earth_sun;
  for (int i = 0; i<6; i++){
    earth_sun.initialize_earth_sun(Nt, T, betas[i], 1.2);
    earth_sun.F_G(0);
    for (int j = 0; j<Nt-1; j++){
      earth_sun.velocity_verlet(j);
    }
    earth_sun.write_to_file("data/elliptical_earth_sun");
  }

  //Then explore different escape velocities with a fixed sun
  solar_system esc_vels;
  double escv[3] = {pow(1.5, 0.5), pow(2, 0.5), pow(2.5, 0.5)};
  string filename0;
  beta = 2;
  for (int i = 0; i<3; i++){
    esc_vels.initialize_earth_sun(Nt, T, beta, escv[i]);
    esc_vels.F_G(0);
    for (int j = 0; j<Nt-1; j++){
      esc_vels.velocity_verlet(j);
    }
    filename0 = "data/";
    std::stringstream params;
    params << std::fixed <<std::setprecision(2) << "v" << escv[i];
    filename0.append(params.str()).append("escvels");
    esc_vels.write_to_file(filename0);
  }
  //Then explore different jupiter masses with a fixed sun
  N = 3;
  Nt = 20000;
  T = 1e9;
  beta = 2;
  double massfactors[3] = {1, 10, 1000};
  solar_system jupiter;
  string filename;
  for (int i = 0; i<3; i++){
    jupiter.initialize(N, Nt, T, beta, 1);
    jupiter.set_jupiter_mass(massfactors[i]);
    jupiter.center_sun();
    jupiter.F_G(0);
    for(int l=0; l<Nt-1; l++){
      jupiter.velocity_verlet(l);
    }
    filename = "data/";
    std::stringstream params;
    params << std::fixed <<std::setprecision(1) << "m" << massfactors[i];
    filename.append(params.str()).append("jupiter");
    jupiter.write_to_file(filename);
  jupiter.save_energies("jupiter");   //only saves the most extreme case
  }

  //Different jupiter masses without a fixed sun
  string filename2;
  for (int i = 0; i<3; i++){
    jupiter.initialize(N, Nt, T, beta, 0);
    jupiter.set_jupiter_mass(massfactors[i]);
    jupiter.remove_drift();
    jupiter.F_G(0);
    for(int l=0; l<Nt-1; l++){
      jupiter.velocity_verlet(l);
    }
    filename2 = "data/nf";
    std::stringstream params;
    params << std::fixed <<std::setprecision(1) << "m" << massfactors[i];
    filename2.append(params.str()).append("jupiter");
    jupiter.write_to_file(filename2);
  jupiter.save_energies("nfjupiter");    //only saves the most extreme case
  }


  //No longer a fixed sun with all planets in the solar system.

  N = 10;
  Nt = 80000;
  T = 1e10;
  solar_system solar_solver;
  solar_solver.initialize(N, Nt, T, beta, 0);
  solar_solver.remove_drift();
  solar_solver.F_G(0);
  for(int l=0; l<Nt-1; l++){
    solar_solver.velocity_verlet(l);
  }
  solar_solver.write_to_file("data/solar");
  solar_solver.save_energies("solar");

  //Investigate mercury perihelion precession (fixed sun again)
  solar_system mercury;
  T = 100*88*24*60*60;
  Nt = 100000;
  //solar_solver.write_to_file(filename);
  mercury.initialize_mercury_sun(Nt, T, beta);
  mercury.remove_drift();
  mercury.F_G_corrected(0);
  for(int k=0; k<Nt-1; k++){
    mercury.mercury(k);
  }
  mercury.write_to_file("data/mercury");
}

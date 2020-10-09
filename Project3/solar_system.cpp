#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

void solar_system::initialize_system(int N, int Nt){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
  m_x = new double[Nparticles];
  m_y = new double[Nparticles];
  m_z = new double[Nparticles];
  m_vx = new double[Nparticles];
  m_vy = new double[Nparticles];
  m_vz = new double[Nparticles];
  mass = new double[Nparticles];
  char* filename_initial = "initial.dat";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_initial, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r") //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i=0; i<m_N; i++){
    fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]); // One %lf (lf=long float or double) for each floating point number on each line of the file.
    fscanf(fp_mass, "%lf", &mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.
}

void solar_system::remove_drift(){
  /*
  removes total momentum of solar system to look at a "stationary" sun
  */
}

/*
void solar_system::F_G(){
  double G = 6.67e-11;
  vec F = -G*ones(3*m_N);
  //Loop over planets
  for(int k=0; k<m_N; k++){
    ind = k/3;
    //Adding force from the sun
    //Maybe there's a way to write this cleaner?
    F(ind) *= m(k);
    F(ind+1) *= m(k);
    F(ind+2) *= m(k);
    //Calculate forces from each planet on this planet
    for(int j=0; j<m_N-1; j++){
      //Do something
    }
  }

}
*/

void solar_system::velocity_verlet(){

}

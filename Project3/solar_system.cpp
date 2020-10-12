#include "solar_system.hpp"
#include <fstream>

using namespace std;

void solar_system::initialize(int N, int Nt){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
  m_x = new double[m_N*m_Nt];
  m_y = new double[m_N*m_Nt];
  m_z = new double[m_N*m_Nt];
  m_vx = new double[m_N*m_Nt];
  m_vy = new double[m_N*m_Nt];
  m_vz = new double[m_N*m_Nt];
  m_mass = new double[m_N*m_Nt];
  char* filename_initial = "initial.txt";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_initial, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i=0; i<m_N; i++){
    fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &m_x[i], &m_y[i], &m_z[i], &m_vx[i], &m_vy[i], &m_vz[i]);
    // One %lf (lf=long float or double) for each floating point number on each line of the file.
    fscanf(fp_mass, "%lf", &m_mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.
}

void solar_system::remove_drift(){
  /*
  removes total momentum of solar system to look at a "stationary" sun
  */
}

void solar_system::F_G(int m){
  double G = 6.67e-11;
  double F[3] = {};
  //Loop over planets
  //F[];
  for(int k=0; k<m_N; k++){
    if(k != m){
      F += ;
    }
  }
  return F;
}


void solar_system::velocity_verlet(){

}

void solar_system::write_to_file(string filename)
{
  /* Stores the positions and velocities */
  m_ofile.open(filename);
  m_ofile << "x,y,z,vx,vy,vz" << endl;
  for (int i=0; i<m_Nt*m_N; i++){
    m_ofile << m_x[i] << "," << m_y[i] << "," << m_z[i] << "," << m_vx[i] << ","
    << m_vy[i] << "," << m_vz[i] << endl;
    }
  m_ofile.close();
}

#include "solar_system.hpp"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void solar_system::initialize(int N, int Nt, double T){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
  //Step-length
  h = T/Nt;
  m_x = new double[m_N*m_Nt];
  m_y = new double[m_N*m_Nt];
  m_z = new double[m_N*m_Nt];
  m_vx = new double[m_N*m_Nt];
  m_vy = new double[m_N*m_Nt];
  m_vz = new double[m_N*m_Nt];
  m_ax = new double[m_N*m_Nt];
  m_ay = new double[m_N*m_Nt];
  m_az = new double[m_N*m_Nt];
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

  //Loop over every acceleration and set it to zero
  for(int k=0; k<m_N*m_Nt; k++){
    m_ax[k] = 0;
    m_ay[k] = 0;
    m_az[k] = 0;
  }
}

void solar_system::remove_drift(){
  /*
  removes total momentum of solar system to look at a "stationary" sun
  */
  double R[3] = {0,0,0};
  double V[3] = {0,0,0};
  for(int i=0; i<m_N; i++){
    R[0] += m_mass[i]*m_x[i];
    R[1] += m_mass[i]*m_y[i];
    R[2] += m_mass[i]*m_z[i];
    V[0] += m_mass[i]*m_vx[i];
    V[1] += m_mass[i]*m_vy[i];
    V[2] += m_mass[i]*m_vz[i];
  }
  double r_den[3] = {0,0,0};
  double v_den[3] = {0,0,0};
  for(int j=0; j<m_N; j++){
    r_den[0] += m_mass[j];
    r_den[1] += m_mass[j];
    r_den[2] += m_mass[j];
    v_den[0] += m_mass[j];
    v_den[1] += m_mass[j];
    v_den[2] += m_mass[j];
  }
  for(int k=0; k<3; k++){
  R[k] /= r_den[k];
  V[k] /= v_den[k];
  }
  for(int l=0; l<m_N; l++){
    m_x[l] -= R[0];
    m_y[l] -= R[1];
    m_z[l] -= R[2];
    m_vx[l] -= V[0];
    m_vy[l] -= V[1];
    m_vz[l] -= V[2];
  }
}

void solar_system::F_G(int m){
  /*
  Method to calculate gravity between objects. Takes the index m to
  keep track of timestep to extract positions
  */
  double G = 6.67e-11;
  double r_norm;
  //Loop over every object
  for(int k=0; k<m_N; k++){
    //Loop over every other object to calculate gravity
    for(int j=0; j<m_N; j++){
      if(j!=k){
        r_norm = pow(((m_x[m+k] - m_x[m+j])*(m_x[m+k] - m_x[m+j]) +
                  (m_y[m+k] - m_y[m+j])*(m_y[m+k] - m_y[m+j]) +
                  (m_z[m+k] - m_z[m+j])*(m_z[m+k] - m_z[m+j])), 1.5);
        m_ax[m+k] += m_mass[j]*(m_x[m+k] - m_x[m+j])/r_norm;
        m_ay[m+k] += m_mass[j]*(m_y[m+k] - m_y[m+j])/r_norm;
        m_az[m+k] += m_mass[j]*(m_z[m+k] - m_z[m+j])/r_norm;
      }
    }
    m_ax[m+k] *= -G;
    m_ay[m+k] *= -G;
    m_az[m+k] *= -G;
  }
}


void solar_system::velocity_verlet(int m){
  /*
  Evolves system one time-step. Takes the index m to
  keep track of timestep to extract positions
  */
  m *= m_N;
  for(int i=0; i<m_N; i++){
    m_x[m+m_N+i] = m_x[m+i] + h*m_vx[m+i] + h*h/2*m_ax[m+i];
    m_y[m+m_N+i] = m_y[m+i] + h*m_vy[m+i] + h*h/2*m_ay[m+i];
    m_z[m+m_N+i] = m_z[m+i] + h*m_vz[m+i] + h*h/2*m_az[m+i];
  }
  F_G(m+m_N);
  for(int j=0; j<m_N; j++){
    m_vx[m+m_N+j] = m_vx[m+j] + h/2*(m_ax[m+m_N+j] + m_ax[m+j]);
    m_vy[m+m_N+j] = m_vy[m+j] + h/2*(m_ay[m+m_N+j] + m_ay[m+j]);
    m_vz[m+m_N+j] = m_vz[m+j] + h/2*(m_az[m+m_N+j] + m_az[m+j]);
  }
}

void solar_system::write_to_file(string filename)
{
  /*
  Stores the positions and velocities in a file
  */
  m_ofile.open(filename);
  m_ofile << "x,y,z,vx,vy,vz" << endl;
  for (int i=0; i<m_Nt*m_N; i++){
    m_ofile << m_x[i] << "," << m_y[i] << "," << m_z[i] << "," << m_vx[i] << ","
    << m_vy[i] << "," << m_vz[i] << endl;
    }
  m_ofile.close();
}

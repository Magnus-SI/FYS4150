#include "solar_system.hpp"
#include <fstream>
#include <cmath>

using namespace std;

void solar_system::initialize(int N, int Nt, double T){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
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
                  (m_z[m+k] - m_z[m+j])*(m_z[m+k] - m_z[m+j])), 3/2);
        m_ax[m+k] += m_mass[j]*(m_x[m+k] - m_x[m+j])/r_norm;
        m_ay[m+k] += m_mass[j]*(m_y[m+k] - m_y[m+j])/r_norm;
        m_az[m+k] += m_mass[j]*(m_z[m+k] - m_z[m+j])/r_norm;
      }
      m_ax[m+k] *= -G;//*m_mass[k/3];
      m_ay[m+k] *= -G;//*m_mass[k/3];
      m_az[m+k] *= -G;//*m_mass[k/3];
      }
    }
}


void solar_system::velocity_verlet(int m){
  for(int i=0; i<m_N; i++){
    m_x[m+m_N+i] = m_x[m+i] + h*m_vx[m+i] + h*h/2*m_ax[m+i];
    m_y[m+m_N+i] = m_y[m+i] + h*m_vy[m+i] + h*h/2*m_ay[m+i];
    m_z[m+m_N+i] = m_z[m+i] + h*m_vz[m+i] + h*h/2*m_az[m+i];
    F_G(m);
    m_vx[m+m_N+i] = m_vx[m+i] + h/2*(m_ax[m+m_N+i] + m_ax[m+i]);
    m_vy[m+m_N+i] = m_vy[m+i] + h/2*(m_ay[m+m_N+i] + m_ay[m+i]);
    m_vz[m+m_N+i] = m_vz[m+i] + h/2*(m_az[m+m_N+i] + m_az[m+i]);
  }
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

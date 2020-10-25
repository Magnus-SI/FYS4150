#include "solar_system.hpp"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void solar_system::initvars(int N, int Nt, double T, double beta){
  /*
  Initializes the class variables
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
  m_beta = beta;
  //Step-length
  h = T/Nt;
  //Positional data
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
}

void solar_system::initialize(int N, int Nt, double T, double beta){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  initvars(N, Nt, T, beta);

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

void solar_system::initialize_earth_sun(int Nt, double T, double beta){
  /*
  Initializes the case of only earth and sun which by the way
  has an analytic solution
  */

  initvars(2, Nt, T, beta);

  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  //Loop over each particle and extract its mass:
  for (int i=0; i<m_N; i++){
    m_x[i] = 0;
    m_y[i] = 0;
    m_z[i] = 0;
    m_vx[i] = 0;
    m_vy[i] = 0;
    m_vz[i] = 0;
    fscanf(fp_mass, "%lf", &m_mass[i]); //Extract mass for particle i.
  }
  fclose(fp_mass); //Close file with masses.
  //Set simple initial conditions for earth, earth is at 1 AU in x direction moving in y direction
  m_x[1] = 1.496e+11; // [meter]
  m_vy[1] = 29789;    // [meter/second]
}

void solar_system::remove_drift(){
  /*
  removes total momentum of solar system to look at a stationary centre of mass system
  */
  double R[3] = {0,0,0};
  double V[3] = {0,0,0};
  double M_tot = 0;
  for(int i=0; i<m_N; i++){
    R[0] += m_mass[i]*m_x[i];
    R[1] += m_mass[i]*m_y[i];
    R[2] += m_mass[i]*m_z[i];
    V[0] += m_mass[i]*m_vx[i];
    V[1] += m_mass[i]*m_vy[i];
    V[2] += m_mass[i]*m_vz[i];
    M_tot += m_mass[i];
  }
  for(int l=0; l<m_N; l++){
    m_x[l] -= R[0]/M_tot;
    m_y[l] -= R[1]/M_tot;
    m_z[l] -= R[2]/M_tot;
    m_vx[l] -= V[0]/M_tot;
    m_vy[l] -= V[1]/M_tot;
    m_vz[l] -= V[2]/M_tot;
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
    //Making sure accerelations are zero before we start adding
    m_ax[m+k] = 0;
    m_ay[m+k] = 0;
    m_az[m+k] = 0;
    //Loop over every other object to calculate gravity
    for(int j=0; j<m_N; j++){
      if(j!=k){
        r_norm = pow(((m_x[m+k] - m_x[m+j])*(m_x[m+k] - m_x[m+j]) +
                  (m_y[m+k] - m_y[m+j])*(m_y[m+k] - m_y[m+j]) +
                  (m_z[m+k] - m_z[m+j])*(m_z[m+k] - m_z[m+j])), (m_beta+1)/(double)2);
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

void solar_system::F_G_corrected(int m){
  /*
  Method to calculate gravity between mercury and sun, corrected for general relativity.
  Takes the index m to keep track of timestep to extract positions
  */
  double G = 6.67e-11;
  double c = 3e8;
  double r_norm, l;
  //Loop over objects
  for(int k=0; k<m_N; k++){
    //Making sure accerelations are zero before we start adding
    m_ax[m+k] = 0;
    m_ay[m+k] = 0;
    m_az[m+k] = 0;
  }
  //Only calculating gravity on mercury, sun is fixed
  r_norm = pow(((m_x[m+6] - m_x[m])*(m_x[m+6] - m_x[m]) +
            (m_y[m+6] - m_y[m])*(m_y[m+6] - m_y[m]) +
            (m_z[m+6] - m_z[m])*(m_z[m+6] - m_z[m])), (m_beta+1)/(double)2);
  //length of cross product between vectors and b
  //l = pow((a2*b3 - a3*b2)**2 + (a3*b1 - a1*b3)**2 + (a1*b2 - a2*b1), 0.5)
  l = pow((m_y[m+6]*m_vz[m+6] - m_z[m+6]*m_vy[m+6])*(m_y[m+6]*m_vz[m+6] - m_z[m+6]*m_vy[m+6]) +
        (m_z[m+6]*m_vx[m+6] - m_x[m+6]*m_vz[m+6])*(m_z[m+6]*m_vx[m+6] - m_x[m+6]*m_vz[m+6]) +
        (m_x[m+6]*m_vy[m+6] - m_y[m+6]*m_vx[m+6])*(m_x[m+6]*m_vy[m+6] - m_y[m+6]*m_vx[m+6]), 0.5);
  m_ax[m+6] += m_mass[0]*(m_x[m+6] - m_x[m])/r_norm*(1 + 3*l*l/pow(r*c,2));
  m_ay[m+6] += m_mass[0]*(m_y[m+6] - m_y[m])/r_norm*(1 + 3*l*l/pow(r*c,2));
  m_az[m+6] += m_mass[0]*(m_z[m+6] - m_z[m])/r_norm*(1 + 3*l*l/pow(r*c,2));
  m_ax[m+6] *= -G;
  m_ay[m+6] *= -G;
  m_az[m+6] *= -G;
}

void solar_system::velocity_verlet(int m){
  /*
  Evolves system one time-step using the velocity-verlet algorithm. Takes the index m to
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

void solar_system::forward_euler(int m){
  /*
  Evolves system one time-step using the forward-euler algorithm. Takes the index m to
  keep track of timestep to extract positions
  */
  m *= m_N;
  for(int i=0; i<m_N; i++){
    //Updating positions
    m_x[m+m_N+i] = m_x[m+i] + h*m_vx[m+i];
    m_y[m+m_N+i] = m_y[m+i] + h*m_vy[m+i];
    m_z[m+m_N+i] = m_z[m+i] + h*m_vz[m+i];
    //Updating velocities
    m_vx[m+m_N+i] = m_vx[m+i] + h*m_ax[m+i];
    m_vy[m+m_N+i] = m_vy[m+i] + h*m_ay[m+i];
    m_vz[m+m_N+i] = m_vz[m+i] + h*m_az[m+i];
  }
  //Updating acceleration
  F_G(m+m_N);
}

void solar_system::mercury(int m){
  /*
  Studying mercury, and adding relativistic correction to gravity
  Sun had index 0, mercury has index 6
  */
  m *= m_N;
  for(int i=0; i<2; i++){
    i*=6;
    m_x[m+m_N+i] = m_x[m+i] + h*m_vx[m+i] + h*h/2*m_ax[m+i];
    m_y[m+m_N+i] = m_y[m+i] + h*m_vy[m+i] + h*h/2*m_ay[m+i];
    m_z[m+m_N+i] = m_z[m+i] + h*m_vz[m+i] + h*h/2*m_az[m+i];
  }
  F_G_corrected(m+m_N);
  for(int j=0; j<2; j++){
    j*=6;
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

#include "solar_system.hpp"
#include <fstream>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
#define G 6.67408e-11
#define c 299792458

using namespace std;

void solar_system::initvars(int N, int Nt, double T, double beta, int fs){
  /*
  Initializes the class variables
  */
  //Number of objects and timesteps.
  m_N = N;
  m_Nt = Nt;
  m_beta = beta;
  //Step-length
  h = T/Nt;
  //Fixed sun (0 if no, 1 if yes)
  fixed_sun = fs;
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
  m_mass = new double[m_N];
}

void solar_system::initialize(int N, int Nt, double T, double beta, int fs){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates.
  */
  initvars(N, Nt, T, beta, fs);

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

void solar_system::initialize_earth_sun(int Nt, double T, double beta, int elliptical){
  /*
  Initializes the case of only earth and sun which by the way
  has an analytic solution
  */

  initvars(2, Nt, T, beta, 1);

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
  m_x[1] = 1.495979e+11; // [meter]
  m_vy[1] = pow(G*m_mass[0]/pow(m_x[1], m_beta-1), 0.5);    // [meter/second]
  if (elliptical == 1){   //elliptical orbit
    //m_vy[1] = 5 * 1.495979e+11/(365*24*3600);
    m_vy[1] *= 1.2;
  }

}

void solar_system::initialize_mercury_sun(int Nt, double T, double beta){
  /*
  Initializes the case of only earth and sun which by the way
  has an analytic solution
  */

  initvars(2, Nt, T, beta, 1);

  char* filename_initial = "initial.txt";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_initial, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &m_x[0], &m_y[0], &m_z[0], &m_vx[0], &m_vy[0], &m_vz[0]);
  fscanf(fp_mass, "%lf", &m_mass[0]);
  //Skip some objects
  char my_string[1000];
  for (int i=0; i<6; i++){
    fgets(my_string, 1000, fp_init);
    fgets(my_string, 1000, fp_mass);
  }
  fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &m_x[1], &m_y[1], &m_z[1], &m_vx[1], &m_vy[1], &m_vz[1]);
  fscanf(fp_mass, "%lf", &m_mass[1]);

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.
}

void solar_system::set_jupiter_mass(int factor){
  m_mass[2] *= factor;
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

void solar_system::center_sun(){
  /*
  Puts the sun at centrum coordinates with no velocity.
  */
  for(int i=0; i<m_N; i++){
    m_x[i] -= m_x[0];
    m_y[i] -= m_y[0];
    m_z[i] -= m_z[0];
    m_vx[i] -= m_vx[0];
    m_vy[i] -= m_vy[0];
    m_vz[i] -= m_vz[0];
  }
}

void solar_system::F_G(int m){
  /*
  Method to calculate gravity between objects. Takes the index m to
  keep track of timestep to extract positions
  */
  double r_norm;
  //Loop over every object
  for(int k=0; k<m_N; k++){
    //Making sure accerelations are zero before we start adding
    m_ax[m+k] = 0;
    m_ay[m+k] = 0;
    m_az[m+k] = 0;
    //Loop over every other object to calculate gravity
    if(k>=fixed_sun){
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
}

void solar_system::F_G_corrected(int m){
  /*
  Method to calculate gravity between mercury and sun, corrected for general relativity.
  Takes the index m to keep track of timestep to extract positions
  */
  double r_norm, v_norm, l_2;
  //Loop over objects
  for(int k=0; k<m_N; k++){
    //Making sure accerelations are zero before we start adding
    m_ax[m+k] = 0;
    m_ay[m+k] = 0;
    m_az[m+k] = 0;
  }
  //Only calculating gravity on mercury, sun is fixed
  r_norm = pow((m_x[m+1] - m_x[m])*(m_x[m+1] - m_x[m]) +
            (m_y[m+1] - m_y[m])*(m_y[m+1] - m_y[m]) +
            (m_z[m+1] - m_z[m])*(m_z[m+1] - m_z[m]), 0.5);
  v_norm = pow((m_vx[m+1] - m_vx[m])*(m_vx[m+1] - m_vx[m]) +
            (m_vy[m+1] - m_vy[m])*(m_vy[m+1] - m_vy[m]) +
            (m_vz[m+1] - m_vz[m])*(m_vz[m+1] - m_vz[m]), 0.5);
  //length of cross product between vectors a and b squared
  l_2 = pow(r_norm, 2) + pow(v_norm,2) - pow(m_x[m+1]*m_vx[m+1] + m_y[m+1]*m_vy[m+1] + m_z[m+1]*m_vz[m+1], 2);
  m_ax[m+1] += m_mass[0]*(m_x[m+1] - m_x[m])/pow(r_norm, 3)*(1 + 3*l_2/pow(r_norm*c,2));
  m_ay[m+1] += m_mass[0]*(m_y[m+1] - m_y[m])/pow(r_norm, 3)*(1 + 3*l_2/pow(r_norm*c,2));
  m_az[m+1] += m_mass[0]*(m_z[m+1] - m_z[m])/pow(r_norm, 3)*(1 + 3*l_2/pow(r_norm*c,2));
  m_ax[m+1] *= -G;
  m_ay[m+1] *= -G;
  m_az[m+1] *= -G;
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
  Sun had index 0, mercury has index 1
  */
  m *= m_N;
  m_x[m+m_N+1] = m_x[m+1] + h*m_vx[m+1] + h*h/2*m_ax[m+1];
  m_y[m+m_N+1] = m_y[m+1] + h*m_vy[m+1] + h*h/2*m_ay[m+1];
  m_z[m+m_N+1] = m_z[m+1] + h*m_vz[m+1] + h*h/2*m_az[m+1];
  m_x[m+m_N] = 0;
  m_y[m+m_N] = 0;
  m_z[m+m_N] = 0;
  F_G_corrected(m+m_N);
  m_vx[m+m_N+1] = m_vx[m+1] + h/2*(m_ax[m+m_N+1] + m_ax[m+1]);
  m_vy[m+m_N+1] = m_vy[m+1] + h/2*(m_ay[m+m_N+1] + m_ay[m+1]);
  m_vz[m+m_N+1] = m_vz[m+1] + h/2*(m_az[m+m_N+1] + m_az[m+1]);
}

double* solar_system::conserved_quants(int m){
  /*
  Returns conserved quantities between the sun and the earth in the 2-planet
  system, given index m, in terms of their original values
  */

  m*= m_N;
  double *cq;
  cq = new double[4];

  double r2init = pow(m_x[1] - m_x[0], 2) + pow(m_y[1] - m_y[0], 2) +
                  pow(m_z[1] - m_z[0], 2);
  double v2init = pow(m_vx[1] - m_vx[0], 2) + pow(m_vy[1] - m_vy[0], 2) +
                  pow(m_vz[1] - m_vz[0], 2);

  double r2 = pow(m_x[m + 1] - m_x[m], 2) + pow(m_y[m + 1] - m_y[m], 2) +
              pow(m_z[m + 1] - m_z[m], 2);
  double v2 = pow(m_vx[m + 1] - m_vx[m], 2) + pow(m_vy[m + 1] - m_vy[m], 2) +
              pow(m_vz[m + 1] - m_vz[m], 2);

  double KEinit = m_mass[1] * v2init;
  double PEinit = G * m_mass[0] * m_mass[1] / r2init;
  double KE = 0.5 * m_mass[1] * v2;
  double PE = G * m_mass[0] * m_mass[1] / r2;

  cq[0] = KE / KEinit;
  cq[1] = PE / PEinit;
  cq[2] = (KE + PE) / (KEinit + PEinit);
  cq[3] = r2/r2init;

  return cq;

}

double solar_system::Kep2Area(int m, int planet_ind){
  /*
  Calculates the Area which should be conserved by Keplers second law.
  Equivalent to conservation of angular momentum.
  At timsetep m for planet_ind. Only works for a fixed sun
  */
  //int ti_min = (m-5) * m_N
  //int ti_max = (m+5) * m_N
  m *= m_N;
  double r2 = pow(m_x[m + planet_ind], 2) + pow(m_y[m + planet_ind], 2) +
                  pow(m_z[m + planet_ind], 2);
  double v2 = pow(m_vx[m + planet_ind], 2) + pow(m_vy[m + planet_ind], 2) +
                  pow(m_vz[m + planet_ind], 2);
  double dotprod =  m_x[m + planet_ind] * m_vx[m + planet_ind] +
                    m_y[m + planet_ind] * m_vy[m + planet_ind] +
                    m_z[m + planet_ind] * m_vz[m + planet_ind];
  double costheta = dotprod/(pow(r2*v2, 0.5));
  double vtheta = pow(1-pow(costheta, 2), 0.5) * pow(v2, 0.5);
  return 0.5 * pow(r2, 0.5) * vtheta;


}

double* solar_system::total_energy(int m){
  /*
  Returns the total energy of the system at the current timestep
  */
  m*= m_N;
  double r2;
  double v2;
  double *tot_E;
  tot_E = new double[4];
  double PE;
  double KE;
  //Loop over every object
  for(int i=0; i<m_N; i++){
    if (i!=0){
      r2 = pow(m_x[m + i] - m_x[m], 2) + pow(m_y[m + i] - m_y[m], 2) +
                  pow(m_z[m + i] - m_z[m], 2);
      PE += G * m_mass[0] * m_mass[i] / pow(r2, 0.5);
    }

    v2 = pow(m_vx[m + i], 2) + pow(m_vy[m + i], 2) + pow(m_vz[m + i], 2);
    KE += 0.5 * m_mass[i] * v2;
  }
  tot_E[0] = PE;
  tot_E[1] = KE;

  return tot_E;
}

void solar_system::write_to_file(string filename)
{
  /*
  Stores the positions and velocities in a file
  */
  std::stringstream params;
  params << std::fixed << m_N << "_"<< std::setprecision(2) << m_beta <<"_" <<log10(m_Nt);
  filename.append(params.str()).append(".txt");

  m_ofile.open(filename);
  m_ofile << "x,y,z,vx,vy,vz" << endl;
  for (int i=0; i<m_Nt*m_N; i++){
    m_ofile << m_x[i] << "," << m_y[i] << "," << m_z[i] << "," << m_vx[i] << ","
    << m_vy[i] << "," << m_vz[i] << endl;
    }
  m_ofile.close();
}

void solar_system::save_energies(string filename){
  /*
  Stores the potential and kinetic energy of the system at certain timesteps
  */
  std::stringstream params;
  params << std::fixed << m_N << "_"<< std::setprecision(2) << m_beta <<"_" <<log10(m_Nt);
  string instring = "data/energy_";
  instring.append(filename).append(params.str()).append(".txt");

  double* E0 = total_energy(0);
  double PE0 = E0[0];
  double KE0 = E0[1];

  ofstream ofile;
  ofile.open(instring);
  ofile <<setw(15) << setprecision(8);
  ofile << "timestep," << "PE," << "KE" <<endl;
  for (int i = 0; i<m_Nt-1; i++){
    if(remainder(i,100) == 0){
      double* E = total_energy(i);
      //totE = E[0] + E[1];
      //REQUIRE(totE/totE0 == Approx(1).epsilon(1e-2));
      ofile << i << "," << E[0] << "," << E[1] <<endl;

    }
  }

}

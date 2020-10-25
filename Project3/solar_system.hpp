#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>

using namespace std;

class solar_system {

protected:
  int m_N, m_Nt; //number of objects and time steps
  int fixed_sun; //whether or not the sun is affected by gravity.
  double *m_x, *m_y, *m_z; //positions
  double *m_vx, *m_vy, *m_vz; //velocities
  double *m_ax, *m_ay, *m_az; //accelerations
  double* m_mass; //masses of sun and planets
  double h; //step size
  double m_beta; //gravity power-law
  ofstream m_ofile;

public:
  void initvars(int N, int Nt, double T, double beta, int fs);
  void initialize(int N, int Nt, double T, double beta, int fs);
  void initialize_earth_sun(int Nt, double T, double beta, int elliptical);
  void initialize_mercury_sun(int Nt, double T, double beta);
  void set_jupiter_mass(int factor);
  void remove_drift();
  void center_sun();
  void velocity_verlet(int m);
  void forward_euler(int m);
  void F_G(int m);
  void F_G_corrected(int m);
  void mercury(int m);
  double* conserved_quants(int m);
  void write_to_file(string filename);
};
#endif

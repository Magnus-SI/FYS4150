#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>

using namespace std;

class solar_system {

protected:
  int m_N, m_Nt; //number of objects and time steps
  double* m_x;
  double* m_y;
  double* m_z; //positions
  double* m_vx;
  double* m_vy;
  double* m_vz; // velocities
  double* m_mass; //masses of sun and planets
  double* F;  //Force between objects (gravity)
  ofstream m_ofile;

public:
  void initialize(int N, int Nt);
  void remove_drift();
  void velocity_verlet();
  void F_G(int m);
  void write_to_file(string filename);
};
#endif

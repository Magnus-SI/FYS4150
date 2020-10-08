#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

class solar_system {

protected:
  int m_N; //number of objects
  mat m_pos, m_vel, m_a;  //positions, velocities, (accelerations ?)
  ofstream m_ofile;

public:
  void initialize(int N);
  void velocity_verlet();
};
#endif

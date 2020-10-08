#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

void solar_system::initialize(str filename){
  /*
  Loads initial positions and velocities for the planets from a file.
  Using	Solar System Barycenter (SSB) coordinates. Ideally finding number
  of planets m_N from the file.
  */
  //m_N = len(filename);
  m_N = 1
  //Store pos and vel in columns because armadillo matrice is column-major
  m_pos = zeros<mat>(3, m_N);
  m_vel = zeros<mat>(3, m_N);
}

void solar_system::velocity_verlet(){

}

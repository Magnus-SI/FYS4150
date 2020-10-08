#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

void solar_system::initialize(int N){
  m_N = N;
  m_pos = zeros<mat>(m_N, m_N);
  m_vel = zeros<mat>(m_N, m_N);
}

void solar_system::velocity_verlet(){

}

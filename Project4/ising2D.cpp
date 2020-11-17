#include "ising2D.hpp"
#include <fstream>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
#define kB 1.38064852e-23     //Boltzmann constant units Joule/Kelvin

using namespace std;

int ising2D::periodic(int i){
  return (m_L*m_L - i)%(m_L*m_L);
}

void ising2D::seed(int s){
  /*
  Gives a random seed
  */
  srand(s);
}

void ising2D::initialize(int L, double temp, double tol){
  /*
  Initializes the 2D lattice in a random state and calculates initial energy,
  magnetisation, and susceptibility.
  tol defines tolerance for spin to point down. If tol = 0, all spins should
  point up.
  */
  m_T = temp; //dimensionless temperature
  m_L = L;
  m_deltaE = 0; m_deltaM = 0;
  m_spin = new int[m_L*m_L];
  m_w = new double[17];
  m_mean = new double[5];
  //
  for(int m=0; m<5; m++) m_mean[m] = 0;
  // setup array for possible energy changes
  for( int de =-8; de <= 8; de++) m_w[de+8] = 0;
  for( int de =-8; de <= 8; de+=4) m_w[de+8] = exp(-de/m_T);
  //Settin up spin matrix and magnetisation
  for(int i=0; i<m_L*m_L; i++){
    m_spin[i] = 1;
    if(rand()/RAND_MAX < tol) m_spin[i] = -1;
    m_deltaM += m_spin[i];
  }
  //Finding energy
  for(int j=0; j<m_L*m_L; j++){
    m_deltaE -= m_spin[j]*(m_spin[periodic(j-1)] +
            m_spin[periodic(m_L - j)]);
  }
  mean_values();
}

void ising2D::metropolis(){
  /*
  Metropolis algorithm
  */
  int ix = rand()/RAND_MAX*m_L;
  int iy = rand()/RAND_MAX*m_L;
  for(int k=0; k<m_L*m_L; k++){
    int ix = rand()/RAND_MAX*m_L;
    int iy = rand()/RAND_MAX*m_L;
    int index = ix + m_L*iy;
    int deltaE = 2*m_spin[index]
      *(m_spin[periodic(index-1)]
      +m_spin[periodic(index+1)]
      +m_spin[periodic(index-m_L)]
      +m_spin[periodic(index+m_L)]);
      // Here we perform the Metropolis test
    if (rand()/RAND_MAX <= m_w[deltaE+8]) {
      m_spin[index] *= -1; // flip one spin and accept new spin config
      // update energy and magnetization
      m_deltaM = (double) 2*m_spin[index];
      m_deltaE = (double) deltaE;
      mean_values();
    }
  }
}

void ising2D::mean_values(){
  /*
  Updates the mean values
  */
  m_mean[0] += m_deltaE; m_mean[1] += m_deltaE*m_deltaE;
  m_mean[2] += m_deltaM; m_mean[3] += m_deltaM*m_deltaM;
  m_mean[4] += fabs(m_deltaM);
}

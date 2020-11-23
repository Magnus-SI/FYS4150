#include "ising2D.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>

using namespace std;

int ising2D::periodic(int i, int ixshift, int iyshift){
  /*
  Ensures periodic boundary conditions. Takes an index from a flattened array,
  converts it to 2D square matrix coordinates, finds correct neighbours and
  converts back to flattened index.
  */
  int ix = i%m_L;
  int iy = i/m_L;
  ix += m_L + ixshift;

  iy += m_L + iyshift;

  return (ix%m_L) + m_L * (iy%m_L);
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
  m_meanstart = 1000000;
  m_meancount = 0;
  /*m_0 = 0;
  m_1 = 0;
  m_2 = 0;
  m_3 = 0;*/
  m_mcs = 0;    //current cycle count
  m_accepted = 0; //accepted spin config count
  m_E = 0; m_M = 0;
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
    if(((float) rand()/RAND_MAX) < tol) m_spin[i] = -1;
    m_M += m_spin[i];
  }
  //Finding energy
  for(int j=0; j<m_L*m_L; j++){
    m_E -= m_spin[j]*(m_spin[periodic(j, -1, 0)] +
            m_spin[periodic(j, 0, -1)]);
  }
}

void ising2D::metropolis(){
  /*
  Metropolis algorithm
  */
  int index;
  for(int k=0; k<m_L*m_L; k++){
    index = ((float) rand()/RAND_MAX) * m_L * m_L;
    //Calculating instances of each index to check uniform distribution
    /*if(index == 0) m_0++;
    if(index == 1) m_1++;
    if(index == 2) m_2++;
    if(index == 3) m_3++;*/
    int deltaE = 2*m_spin[index]
      *(m_spin[periodic(index, -1, 0)]
      +m_spin[periodic(index, 1, 0)]
      +m_spin[periodic(index, 0, -1)]
      +m_spin[periodic(index, 0, 1)]);
      // Here we perform the Metropolis test
    if (((float) rand()/RAND_MAX) <= m_w[deltaE+8]) {
      m_spin[index] *= -1; // flip one spin and accept new spin config
      // update energy and magnetization
      m_M += 2*m_spin[index];
      m_E += deltaE;
      m_accepted++;
    }
  }
  m_mcs++;
  if(m_mcs>=m_meanstart){
    mean_values();
  }
}

void ising2D::mean_values(){
  //Updates the per spin mean values
  m_meancount++;
  m_mean[0] += m_E; m_mean[1] += m_E*m_E;
  m_mean[2] += m_M; m_mean[3] += m_M*m_M;
  m_mean[4] += fabs(m_M);
}

void ising2D::write_to_file(ofstream& ofile){
  /*
  Store quantities in a file
  */
  ofile << m_mcs << "," << m_E << "," << m_M <<  "," << m_accepted << endl;
}

void ising2D::write_mean(stringstream& filedat){

  double Cv = pow(m_T, -2) * (m_mean[1]/m_meancount - pow(m_mean[0]/m_meancount, 2));
  double chi = pow(m_T, -1) * (m_mean[3]/m_meancount - pow(m_mean[4]/m_meancount, 2));
  double ps = (double) 1/(m_L*m_L);  //per spin normalization
  //cout << ps * m_mean[0]/m_meancount << endl;
  filedat << ps * m_mean[0]/m_meancount << "," << ps * m_mean[4]/m_meancount
  << "," << ps * Cv << "," << ps * chi << endl;
}

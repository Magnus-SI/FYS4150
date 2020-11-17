#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

int main(){
  ising2D my_ising;
  int L = 2; double temp = 1.0; double tol = 0.0;
  my_ising.initialize(L, temp, tol);
  for(int mcs=1; mcs<1001; mcs++){
    my_ising.metropolis();
    cout << (double) my_ising.m_mean[0]/mcs << endl;
  }
  return 0;
}

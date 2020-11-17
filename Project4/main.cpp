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
  return 0;
}

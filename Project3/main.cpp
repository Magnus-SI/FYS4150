#include "solar_system.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

int main(){
  string filename = "solar.txt";
  solar_system solar_solver;
  solar_solver.initialize(10, 3);
  solar_solver.write_to_file(filename);

}

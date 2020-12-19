#include "diffusion.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

int main(){

  double dx = 0.025;   //takes almost no time if you want to test code
  //double dx = 0.01;   // used in report, data included

  string filename;
  int nx1 = 1.5/dx - 1;
  int nx2 = 1.5/dx - 1;
  int ny = 1/dx - 1;

  double dt = 0.2 * pow(dx, 2);
  Diffusion diff1, diff2, diff3;
  diff1.init2D(nx1, ny, dx, dt, 1);
  //diff1.update_Q();
  diff2.init2D(nx2, ny, dx, dt, 2);
  //diff2.update_Q();
  diff3.init2D(nx2, ny, dx, dt, 3);
  //diff3.update_Q();
  diff1.Q_0();
  diff2.Q_0();
  diff3.Q_0();

  double t_stable = 0.2;
  int nt_stable = t_stable/dt;

  filename = "data/iso2D_stab.csv";
  ofstream ofile1;
  ofile1.open(filename);
  ofile1 << setw(15) << setprecision(8);

  for (int i = 1; i<nt_stable + 1; i++){
    diff1.EF2Dperx();
    diff2.EF2Dperx();
    diff3.EF2Dperx();
    if (i == nt_stable){
      ofile1 <<i*dt;
      diff1.write_1D(ofile1);
    }
  }
  ofile1.close();

  diff1.m_tc = 0;
  diff2.m_tc = 0;
  diff3.m_tc = 0;
  diff1.update_Q();
  diff2.update_Q();
  diff3.update_Q();


  double t_rad = 0.01877;
  int nt0 = 0.5 * t_rad/dt;
  int nt1 = 1 * t_rad/dt;
  int nt2 = 2 * t_rad/dt;
  int nt3 = 20 * t_rad/dt;

  stringstream params;
  params << fixed << setprecision(2) << log10(dx);
  filename = "data/iso2D_dx_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << setw(15) << setprecision(8);


  for (int i = 1; i<nt3+1; i++){
    diff1.EF2Dperx();
    diff2.EF2Dperx();
    diff3.EF2Dperx();
    diff3.update_Q();
    if (i == nt0 || i == nt1 || i == nt2 || i == nt3){
      ofile << "case1_" <<i*dt;
      diff1.write_1D(ofile);
      ofile << "case2_" <<i*dt;
      diff2.write_1D(ofile);
      ofile << "case3_" <<i*dt;
      diff3.write_1D(ofile);
    }
  }
  ofile.close();
  //dt_err(0.1, 10);
  //dt_err(0.1, 100);

  return 0;
}

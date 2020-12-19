#include "diffusion.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

void comp_methods(double dx, double t1, double t2){
  double dt = 0.4 * pow(dx, 2);
  int nx = 1/dx - 1;

  int nt1 = t1/dt;
  int nt2 = t2/dt;

  Diffusion diff1, diff2, diff3;
  diff1.init(nx, dx, dt, 0);
  diff2.init(nx, dx, dt, 0);
  diff3.init(nx, dx, dt, 1);

  stringstream params;
  params << setprecision(2) << log10(dx);
  string filename = "data/mcdep_dx_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << setw(15) << setprecision(8);

  //ofile << "Efor" << 0;
  //diff1.write_1D(ofile);

  for (int i = 1; i<(nt2+1); i++){
    diff1.EulerForward();
    diff2.tridiag_solve();
    diff3.CrankNicolson();

    if (i == nt1 || i == nt2){
      ofile << "Efor_" << i * dt;
      diff1.write_1D(ofile);
      ofile << "imp1_" << i * dt;
      diff2.write_1D(ofile);
      ofile << "imp2_" << i * dt;
      diff3.write_1D(ofile);
      ofile << "exct_" << i * dt;
      diff3.write_analytic(ofile, 500);
    }
  }

  ofile.close();
}

void dt_err(double dx, double t_end){
  int kmax = 500;
  int nx = 1/dx - 1;
  double dt, nt;
  stringstream params;
  params << fixed << setprecision(3) << log10(dx) <<"_" << t_end;
  string filename = "data/dterr_dxdt_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << setw(15) << setprecision(8);
  ofile << "dtm,Efor,imp1,imp2" << endl;
  Diffusion diff1, diff2, diff3;
  for (double i = -6; i<3; i+= 0.5){
    dt = pow(2, i) * pow(dx, 2);
    nt = t_end/dt;
    //cout << t_end << " " << nt<< endl;
    diff1.init(nx, dx, dt, 0);
    diff2.init(nx, dx, dt, 0);
    diff3.init(nx, dx, dt, 1);

    for (int i = 1; i<nt; i++){
      diff1.EulerForward();
      diff2.tridiag_solve();
      diff3.CrankNicolson();
    }
    ofile << pow(2,i);
    diff1.write_andiff(ofile, kmax);
    diff2.write_andiff(ofile, kmax);
    diff3.write_andiff(ofile, kmax);
    ofile << endl;
  }
  ofile.close();

  //use different dt multipliers, say in 0.05 to 0.5, compare methods.
}

void err2d(double t_end, int gl){

  double dtmults[3] = {0.1, 0.2, 0.4};
  double dxvals[4] = {0.025, 0.05, 0.1, 0.2};

  string filename;
  if (gl!= 0){
    stringstream params;
    params << fixed << setprecision(2) << t_end;
    filename = "data/err2D_tend_";
    filename.append(params.str()).append(".csv");
  }
  else{
    filename = "data/err2D_local.csv";
  }

  ofstream ofile;
  ofile.open(filename);
  ofile << "dt,dx,err" << endl;
  ofile << setprecision(8);


  double dt, dx, err;
  int nxy, nt;
  Diffusion diff2D;

  for (int j = 0; j<4; j++){
    for (int k = 0; k<3; k++){
      dx = dxvals[j];
      dt = dtmults[k] * pow(dx, 2);
      nxy = 1/dx - 1;
      if (gl!= 0){nt = t_end/dt;}
      else {nt = t_end/dt - 1;}    //local error instead
      //cout << nt << ",";


      diff2D.init2D(nxy, nxy, dx, dt, -1);
      diff2D.Q_0();

      for (int i =1; i<(nt+1); i++){
        diff2D.EF2Dperx();

      }
      if (gl!= 0){err = diff2D.err2D(100);}
      else {err = diff2D.err2Dlocal(100, 0);}
      ofile << dt << "," << dx << "," << err << endl;
    }
  }

  ofile.close();
}

void plot2d(){
  double dxy = 0.01;
  double dt = 0.15 * pow(dxy, 2);
  int nxy = 1/dxy - 1;
  Diffusion diff2D;
  diff2D.init2D(nxy, nxy, dxy, dt, 0);
  int nt1 = 40;
  int nt2 = 400;
  int nt3 = 4000;

  stringstream params;
  params << fixed << setprecision(2) << log10(dxy);
  string filename = "data/diff2D_dx_";
  filename.append(params.str()).append(".csv");
  ofstream ofile;
  ofile.open(filename);
  ofile << setw(15) << setprecision(8);


  for (int i = 1; i<nt3+1; i++){
    diff2D.EF2D();
    if (i == nt1 || i == nt2 || i == nt3){
      ofile << i;
      diff2D.write_1D(ofile);
    }
  }
}

int main(){


  //First run for 2x2 lattice, T = 1.0
  double dx1 = 0.1;
  double dx2 = 0.01;
  double t1 = 0.1;
  double t2 = 0.3;


  comp_methods(dx1, t1, t2);
  comp_methods(dx2, t1, t2);
  dt_err(dx1, 0.1);
  dt_err(dx1, 0.3);
  dt_err(dx2, 0.1);
  dt_err(dx2, 0.3);


  err2d(0.2, 1);
  err2d(0.2, 0);
  plot2d();



  return 0;
}

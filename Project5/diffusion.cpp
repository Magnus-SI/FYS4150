#include "diffusion.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>

using namespace std;


void Diffusion::init(int nx, double dx, double dt, int impmet){
  /*
  Initializes the diffusion equation in the 1D case.
  impmet = 0 for euler backward and impmet = 1 for Crank-Nicolson
  */
  u = new double[nx+2];
  uprev = new double[nx+2];
  uexact = new double[nx+2];
  b = new double[nx+2];
  d = new double[nx+2];
  a = new double[nx+1];
  c = new double[nx+1];
  uprev = new double[nx+2];
  m_alpha = dt/pow(dx, 2);

  for (int i = 1; i<nx+1; i++){
    u[i] = 0;
    uprev[i] = 0;
    d[i] = 0;
    a[i] = c[i] = -m_alpha;
    if (impmet == 0){
      b[i] = 1 + m_alpha * 2;
    }
    else if (impmet == 1){
      b[i] = 2 + 2*m_alpha;
    }
  }
  d[0] = d[nx + 1] = b[0] = b[nx+1] = 1;
  //a[0] = c[nx] = -m_alpha;
  a[nx] = c[0] = 0;       //adjust matrix for boundary conditions

  u[0] = uprev[0] = 0.0;
  u[nx+1] = uprev[nx+1] = 1.0;
  m_nx = nx+2;
  m_n2 = nx+2;
  m_dx = dx;
  m_dt = dt;
  m_tc = 0;   //current time

}

void Diffusion::init2D(int nx, int ny, double dxy, double dt, int Qtype){
  /*
  Initializes the diffusion equation in the 2D case. The Qtype parameter
  is related to the initial conditions as well as which type of heat should
  be added in the isothermal case.
  */
  int xdim = nx + 2;
  int ydim = ny + 2;
  int dim = xdim * ydim;
  u = new double[dim];
  uprev = new double[dim];
  uexact = new double[ny+2];    //closed form solution for periodic x boundaries
  Q = new double [dim];
  m_alpha = dt/pow(dxy, 2);

  if (Qtype == 0){
    for (int i = 0; i<dim; i++){
      if (i%xdim == nx+1 || i/xdim == nx + 1) {u[i] = 1;}
      else {u[i] = 0;}
      uprev[i] = u[i];
    }
  }

  else if (Qtype == -1){
    for (int i = 0; i<dim; i++){
      if (i%xdim == 0 || i%xdim == nx + 1) {u[i] = 1;}
      else {u[i] = 0;}
      uprev[i] = u[i];
    }
  }
  else{
    for (int i = 0; i<dim; i++){
      if (i/xdim == 0) {u[i] = 8.0/1300;}
      else if (i/xdim == ny + 1) {u[i] = 1;}
      else {u[i] = 0;}
      uprev[i] = u[i];
    }
  }

  m_nx = nx + 2;
  m_ny = ny + 2;
  m_dx = dxy;
  m_dt = dt;
  m_n2 = dim;
  m_Qtype = Qtype;
  m_tc = 0;
}

void Diffusion::update_Q(){
  /*
  Updates Q according to m_Qtype.
  */
  int ij;
  double Q_base, h, w, Qtot, Qdecay;
  double Qfac = 4.431;
  for (int j = 0; j<m_ny; j++){

    h = j*m_dx;  //dx = dy

    //outside mantle:
    if (h< 4.0/12){
      Q_base = Qfac;
      if (h < 2.0/12){Q_base *= 1.4;}
      else if (h >=2.0/12){Q_base *= 0.35;}
      for (int i = 0; i<m_nx; i++){
        ij = i + j*m_nx;
        Q[ij] = Q_base;
      }
    }
    //inside mantle:
    else{
      Q_base = Qfac;
      Q_base *= 0.05;
      if (m_Qtype == 1){
        for (int i = 0; i<m_nx; i++){
          ij = i + j*m_nx;
          Q[ij] = Q_base;
        }
      }
      else if (m_Qtype == 2 || m_Qtype == 3){
        for (int i = 0; i<m_nx; i++){
          w = i * m_dx;
          ij = i + j*m_nx;
          //Above slab:
          if (w>=0.25 & w<=1.25) {
            //No decay:
            if (m_Qtype == 2){Q[ij] = Q_base + 0.5 * Qfac;}
            //Yes decay:
            else{
              Qdecay = 0.4 * exp(-log(2)*m_tc/0.08391)    //U-decay
                     + 0.4 * exp(-log(2)*m_tc/0.2628)    //Th-decay
                     + 0.2 * exp(-log(2)*m_tc/0.02346);   //K-decay
              Q[ij] = Q_base + 0.5 * Qfac * Qdecay;
            }
          }
          //To the side of the slab:
          else {Q[ij] = Q_base;}
        }
      }
    }
  }
}

void Diffusion::Q_0(){
  /*
  Sets Q to zero.
  */
  for (int i = 0; i<m_n2; i++){Q[i] = 0;}
}

void Diffusion::EulerForward(){
  /*
  Performs the 1D euler forward.
  */
  for (int i = 1; i<m_nx-1; i++){
    u[i] = m_alpha * (uprev[i-1] + uprev[i+1]) + (1 - 2*m_alpha) * uprev[i];
  }
  for (int i = 1; i<m_nx-1; i++){
    uprev[i] = u[i];
  }
  m_tc += m_dt;
}

void Diffusion::EF2D(){
  /*
  Performs the 2D euler forward with fixed boundaries.
  */
  int ij;
  for (int i = 1; i<m_nx-1; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      u[ij] = m_alpha * (uprev[ij-1] + uprev[ij+1] + uprev[ij-m_nx] + uprev[ij+m_nx])
      +(1 - 4*m_alpha) * uprev[ij];
    }
  }
  for (int i = 1; i<m_nx-1; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      uprev[ij] = u[ij];
    }
  }
  m_tc += m_dt;
}

int Diffusion::xperiodic(int ix, int iy, int ixshift, int iyshift){
  /*
  Ensures periodic boundary conditions along x. Takes an index from a flattened array,
  converts it to 2D square matrix coordinates, finds correct neighbours and
  converts back to flattened index.
  */
  int ixs, iys;
  ixs = ix + m_nx + ixshift;

  iys = iy + iyshift;

  return (ixs%m_nx) + m_nx * iys;
}

void Diffusion::EF2Dperx(){
  /*
  Performs the 2D euler forward with periodic boundaries along x
  */
  int ij;
  for (int i = 0; i<m_nx; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      u[ij] = m_alpha * (uprev[xperiodic(i,j,-1,0)] + uprev[xperiodic(i,j,1,0)]
            + uprev[xperiodic(i,j,0,-1)] + uprev[xperiodic(i,j,0,1)])
      +(1 - 4*m_alpha) * uprev[ij] + Q[ij] * m_dt;
    }
  }
  for (int i = 0; i<m_nx; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      uprev[ij] = u[ij];
    }
  }
  m_tc += m_dt;
}

void Diffusion::tridiag_solve(){
  /*
  Solves the tridiagonal matrix equation.
  */
    //forward substitution
  for (int i = 1; i < m_nx; i++){
    d[i] = b[i] - a[i-1]*c[i-1]/d[i-1];
    uprev[i] -= a[i-1]*uprev[i-1]/d[i-1];
    }
    //backward substitution
  u[m_nx-1] = uprev[m_nx-1]/b[m_nx-1];
  for (int i = m_nx-2; i>=0; i--){
    u[i] = (uprev[i]-c[i]*u[i+1])/d[i];
    uprev[i] = u[i];    //reset uprev to u
    }
  uprev[m_nx-1] = u[m_nx-1];
  m_tc += m_dt;
}

void Diffusion::CrankNicolson(){
  /*
  Crank Nicolson is performed by letting B_2*uprev act as the right hand side
  of the tridiagonal solver
  */
  for (int i=1; i<m_nx-1; i++){
    uprev[i] = m_alpha * (u[i-1] + u[i+1]) + (2.0 - 2.0*m_alpha) * u[i];
  }
  tridiag_solve();
}

void Diffusion::GetExact(int kmax){
  /*
  Get the exact solution as described in the report
  */
  int k;
  double uprev;
  for (int i = 1;  i<m_nx - 1; i++){
    uprev = 0;
    uexact[i] = m_dx*i;
    k = 1;
    while (abs(uexact[i] - uprev) > 1e-8 & k<kmax+1){
      uprev = uexact[i];
      uexact[i] += 2 * pow(-1, k)/(k*M_PI) * sin(k*M_PI*m_dx*i)
      * exp(-pow(k*M_PI,2) * m_tc);
      //cout << k << "_" << uexact[i] << uprev << endl;
      k++;
    }
  }
  uexact[0] = 0;
  uexact[m_nx - 1] = 1;
}

void Diffusion::write_andiff(ofstream& ofile, int kmax){
  /*
  Writes difference between exact and current solution. Works for 1d case.
  */
  GetExact(kmax);
  double diff_sum = 0;
  for (int i = 1; i<m_nx - 1; i++){
    diff_sum += abs(u[i] - uexact[i]);
  }
  double diff_mean = diff_sum/(m_nx - 2);
  ofile << "," << diff_mean;
}

void Diffusion::write_1D(ofstream& ofile){
  /*
  Writes a 1d array of u, in the 2d case it is flattened here.
  */
  for (int i = 0; i<m_n2; i++){
    ofile << "," << u[i];
  }
  ofile << endl;
}

double Diffusion::err2D(int kmax){
  /*
  Calculates the mean absolute error in the 2d case.
  */
  int ij;
  GetExact(kmax);
  double diff_sum = 0;
  for (int i = 0; i<m_nx; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      diff_sum += abs(u[ij] - uexact[j]);
    }
  }
  //cout << u[m_nx*m_ny/2 + 5] << "," << uexact[m_ny/2] << endl;
  double diff_mean = diff_sum/((m_nx) * (m_ny-2));
  return diff_mean;
}

double Diffusion::err2Dlocal(int kmax, int test){
  /*
  Calculates the local mean absolute error in the 2d case by setting
  u to its exact value and advancing one step before calculating the error.
  The parameter test is used for unit testing.
  */
  int ij;
  GetExact(kmax);
  for (int i = 0; i<m_nx; i++){
    for (int j = 1; j<m_ny-1; j++){
      ij = i + j*m_nx;
      u[ij] = uexact[j];
      uprev[ij] = uexact[j];
    }
  }
  if (test == 0){EF2Dperx();}

  return err2D(kmax);
}

void Diffusion::write_analytic(ofstream& ofile, int kmax){
  /*
  Stores the analytic solution to a file.
  */
  GetExact(kmax);
  for (int i = 0; i<m_nx; i++){
    ofile << "," << uexact[i];
  }
  ofile << endl;
}

void Diffusion::write_Q(ofstream& ofile){
  /*
  Writes Q instead of u.
  */
  for (int i = 0; i<m_n2; i++){
    ofile << "," << Q[i];
  }
  ofile << endl;
}

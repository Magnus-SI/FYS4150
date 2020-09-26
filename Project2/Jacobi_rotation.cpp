#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
using namespace arma;

void Jacobi_rotation::initialize(double a, double b, int N, double V(double rho))
{
  /* Initializing the problem, the Hamiltonian matrix, etc. */
  m_N = N;
  m_h = (b - a)/(m_N+1);
  m_Hamiltonian = zeros<mat>(m_N, m_N);

  double  DiagConst, NondiagConst;
  DiagConst = 2.0 / (m_h*m_h);
  NondiagConst =  -1.0 / (m_h*m_h);
  vec x = linspace(a + m_h, b-m_h, m_N);

  // Setting up tridiagonal matrix and diagonalization using Armadillo
  m_Hamiltonian(0,0) = DiagConst;
  m_Hamiltonian(0,1) = NondiagConst;
  for(int i = 1; i < m_N-1; i++) {
    m_Hamiltonian(i,i-1)    = NondiagConst;
    m_Hamiltonian(i,i)    = DiagConst + V(x(i));
    m_Hamiltonian(i,i+1)    = NondiagConst;
  }
  m_Hamiltonian(m_N-1,m_N-2) = NondiagConst;
  m_Hamiltonian(m_N-1,m_N-1) = DiagConst;
}

void Jacobi_rotation::print_matrix()
{
  /* Just printing the matrix. Should only be used for testing on small matrices */
  m_Hamiltonian.print();
}

void Jacobi_rotation::rotate(double conv, vec r, mat v)
{
  /* Method to rotate matrix */
  double aip=0, aiq=0, vpi=0, vqi=0;
  double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)
  int count=1;                //count of iterations
  int count_old=count-10;     //keep track of every 10th iteration
  int p=m_N-1, q=m_N-2;           //off diag all same value to start
                              //pick last as first maximum
  clock_t start, end;

  mat a = m_Hamiltonian;

  double app=a(p,p);
  double aqq=a(q,q);
  double apq=a(p,q);

  while(abs(apq)>conv){
      if(count>1){
          apq=0;
          for (int i=0;i<m_N;i++){
               for (int j=0;j<m_N;j++){
                  if(i!=j && abs(a(i,j))>=abs(apq)){
                      apq=a(i,j);
                      p=i;
                      q=j;
                    }
              }
          }
      }

      //calculate sin(theta) and cos(theta)
      aqq=a(q,q);
      app=a(p,p);
      tau=(aqq-app)/(2*apq);
      if(tau>0)
          t=1/(tau+sqrt(1+tau*tau));
      else
          t=-1/(-tau+sqrt(1+tau*tau));
      c=1/sqrt(1+t*t);
      s=c*t;

      //calculate new matrix elements and vectors
      for(int i=0;i<m_N;i++){
          if(i!=p && i!=q){
              aip=a(i,p);
              aiq=a(i,q);
              a(i,p)=aip*c-aiq*s;
              a(p,i)=aip*c-aiq*s;
              a(i,q)=aiq*c+aip*s;
              a(q,i)=aiq*c+aip*s;
          }
          //vpi=v(p,i);
          //vqi=v(q,i);
          vpi=v(i,p);
          vqi=v(i,q);
          //v(p,i)=c*vpi-s*vqi;
         // v(q,i)=c*vqi+s*vpi;
          v(i,p)=c*vpi-s*vqi;
          v(i,q)=c*vqi+s*vpi;
      }
      a(p,p)=app*c*c-2*apq*c*s+aqq*s*s;
      a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
      a(p,q)=0;
      a(q,p)=0;

      count++;
  }
  a.print();
}

void Jacobi_rotation::write_to_file(string filename)
{
  /* Needs to be tailored to our problem */
  m_ofile.open(filename);
  for (int i = 0; i < m_N; i++){
    m_ofile << m_x(i) << " " << m_v(i) << endl;
  }
  m_ofile.close();
}

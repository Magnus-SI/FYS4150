#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
using namespace arma;

double Jacobi_rotation::V_func(double rho)
{
  /* Zero potential */
  if(m_Vchoice == 0){
    return 0;
  }
  if(m_Vchoice ==1){
    return pow(rho, 2); //one electron quantum potential
  }
  if(m_Vchoice == 2){
    return pow(omega_r * rho, 2) + 1/rho; //two electron quantum potential
  }
  else{
    cout<<"whaddup"<<endl;
    return 0;
  }
}

void Jacobi_rotation::initialize(double a, double b, int N, int Vchoice)
{
  /* Initializing the problem, the Hamiltonian matrix, etc. */
  m_N = N;
  m_h = (b - a)/(m_N+1);
  m_Hamiltonian = zeros<mat>(m_N, m_N);
  m_Vchoice = Vchoice;

  double  DiagConst, NondiagConst;
  DiagConst = 2.0 / (m_h*m_h);
  NondiagConst =  -1.0 / (m_h*m_h);
  vec x = linspace(a + m_h, b-m_h, m_N);

  // Setting up tridiagonal matrix and diagonalization using Armadillo
  m_Hamiltonian(0,0) = DiagConst + V_func(x(0));
  m_Hamiltonian(0,1) = NondiagConst;
  for(int i = 1; i < m_N-1; i++) {
    m_Hamiltonian(i,i-1)    = NondiagConst;
    m_Hamiltonian(i,i)    = DiagConst + V_func(x(i));
    m_Hamiltonian(i,i+1)    = NondiagConst;
  }
  m_Hamiltonian(m_N-1,m_N-2) = NondiagConst;
  m_Hamiltonian(m_N-1,m_N-1) = DiagConst + V_func(x(m_N-1));
}

void Jacobi_rotation::print_matrix()
{
  /* Just printing the matrix. Should only be used for testing on small matrices */
  m_Hamiltonian.print();
}

void Jacobi_rotation::rotate(double conv)
{
  /* Method to rotate matrix */
  double aip=0, aiq=0, vpi=0, vqi=0;
  double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)
  int count=1;                //count of iterations
  int count_old=count-10;     //keep track of every 10th iteration
  int p=m_N-1, q=m_N-2;           //off diag all same value to start
  int max_iter = 30000;
                              //pick last as first maximum
  clock_t start, end;

  A = m_Hamiltonian;
  V = mat(m_N, m_N, fill::eye);

  double app=A(p,p);
  double aqq=A(q,q);
  double apq=A(p,q);

  while(abs(apq)>conv && count<max_iter){
      if(count>1){
          apq=0;
          for (int i=0;i<m_N;i++){
               for (int j=0;j<m_N;j++){
                  if(i!=j && abs(A(i,j))>=abs(apq)){
                      apq=A(i,j);
                      p=i;
                      q=j;
                    }
              }
          }
      }

      //calculate sin(theta) and cos(theta)
      aqq=A(q,q);
      app=A(p,p);
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
              aip=A(i,p);
              aiq=A(i,q);
              A(i,p)=aip*c-aiq*s;
              A(p,i)=aip*c-aiq*s;
              A(i,q)=aiq*c+aip*s;
              A(q,i)=aiq*c+aip*s;
          }
          //vpi=v(p,i);
          //vqi=v(q,i);
          vpi=V(i,p);
          vqi=V(i,q);
          //v(p,i)=c*vpi-s*vqi;
         // v(q,i)=c*vqi+s*vpi;
          V(i,p)=c*vpi-s*vqi;
          V(i,q)=c*vqi+s*vpi;
      }
      A(p,p)=app*c*c-2*apq*c*s+aqq*s*s;
      A(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
      A(p,q)=0;
      A(q,p)=0;

      count++;

  iters = count;
  }
  //cout << "Number of iterations at N = " << m_N << " : " << count << endl;
}

void Jacobi_rotation::rearrange()
{
  /* Rearranges order of rotated eigenvectors and eigenvalues with first
  eigenvalues first*/
  vec eigen = zeros(m_N);
  mat V_temp = zeros(m_N, m_N);

  for(int i=0; i<m_N; i++){
      eigen(i) = A(i,i);
  }
  uvec ind = sort_index(eigen);

  for(int i=0; i<m_N; i++){
    A(i,i) = eigen(ind(i));
    for(int j=0; j<m_N; j++){
      V_temp(j,i) = V(j, ind(i));
    }
  }
  V = V_temp;
  cout << endl;
  //V.print();
}

vec Jacobi_rotation::return_eig()
{
  /* Returns a vector of the current eigenvalues */
  vec eigen = zeros(m_N);

  for(int i=0; i<m_N; i++){
      eigen(i) = A(i,i);
  }
  return eigen;
}

mat Jacobi_rotation::test_eig()
{
  /* Returns the eigenvector calculated from the armadillo function eig_sym*/
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, m_Hamiltonian);
  return eigvec;
  //eigvec.print();
}

void Jacobi_rotation::eigarma(vec &eigval, mat &eigvec){
  /* Modifies some defined eigval and eigvec to solve the system using
  the armadillo function eig_sym */
  eig_sym(eigval, eigvec, m_Hamiltonian);
}

void Jacobi_rotation::write_to_file(string filename)
{
  /* Stores the first 3 eigenvectors and the eigenvalue for the current system */
  m_ofile.open(filename);
  m_ofile << "eigenvalues,eigenvector1,eigenvector2,eigenvector3" << endl;
  for (int i = 0; i < m_N; i++){
    m_ofile << A(i,i) << "," << V(i,0) << "," << V(i,1) << "," << V(i,2) << endl;
  }
  m_ofile.close();
}

float Jacobi_rotation::quanteigtest(int method)
{
  /* Tests the quantum eigenvalues by comparing to analytical values. the argument
  method is 0 if the Jacobi rotation should be used, and 1 if the armadillo
  functions should be used. */
  double maxerr = 0;
  vec eigval;
  if(method == 1){
    mat eigvec;
    eig_sym(eigval, eigvec, m_Hamiltonian);
  }
  for (int i = 0; i<5; i++){
    double eigv;
    if(method == 0){
      eigv = A(i,i);
    }
    else{
      eigv = eigval(i);
    }

    double analytic_eigv = 4*i+3;
    if (abs(eigv - analytic_eigv)>maxerr){
      maxerr = abs(eigv - analytic_eigv);
    }
  }
  return maxerr;
}

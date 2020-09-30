#include "catch.hpp"
#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "time.h"

TEST_CASE("Testing eigenvalues/eigenvectors of Toeplitz matrix"){
    //Tests eigenvectors and eigenvalues of a small Toeplitz matrix by
    //comparing to analytical values given in the project description.
    int Dim = 3;
    //    Set up the exact eigenvalues
    vec Exact(Dim);
    mat exactvecs = zeros<mat>(Dim, Dim);
    double pi = acos(-1.0);
    // Integration step length
    //note that example programs didn't use dim+1 here, why?
    double Step    = 1.0/(Dim+1);
    double DiagConst = 2.0 / (Step*Step);
    double NondiagConst =  -1.0 / (Step*Step);
    for(int i = 0; i < Dim; i++) {
      Exact(i) = DiagConst+2*NondiagConst*cos((i+1)*pi/(Dim+1));
      for (int j = 0; j <Dim; j++){
        exactvecs(j,i) = sin((i+1)*(j+1)*pi/(Dim+1)) * pow(2, -0.5); //normalized
        //cout << (i+1)*(j+1)*pi/Dim << endl;
      }
    }

    //get numerical eigenvalues
    Jacobi_rotation my_solver;

    double a = 0;
    double b = 1;
    int Vchoice = 0;
    my_solver.initialize(a, b, Dim, Vchoice);
    double conv = 1e-8;
    my_solver.rotate(conv);
    my_solver.rearrange();
    vec EigvalueNum = my_solver.return_eig();
    mat eigvecNum = my_solver.V;

    for(int i=0; i<Dim; i++){
      REQUIRE(EigvalueNum(i)==Approx(Exact(i)).epsilon(0.00000001));
      for(int j=0; j<Dim; j++){
        REQUIRE(abs(eigvecNum(i,j))==Approx(abs(exactvecs(i,j))).epsilon(0.00000001));
      }
    }

}

TEST_CASE("Eigenvalue deviation for Quantum case"){
  //Tests the maximum eigenvalue deviation function by comparing its value
  //when using armadillo functions with the Jacobi algorithm.
  int Vchoice = 1;
  double rhomin = 0;
  double rhomax = 5;
  vec Nvals = logspace(1,2, 2);
  Jacobi_rotation tester;
  double err1;
  double err2;
  int method1 = 0;    //Jacobi rotation method
  int method2 = 1;    //armadillo eig_sym function
  double conv = 1e-8;
  for (int i=0; i<2; i++){
    tester.initialize(rhomin, rhomax, Nvals(i), Vchoice);
    err2 = tester.quanteigtest(method2);
    tester.rotate(conv);
    tester.rearrange();
    err1 = tester.quanteigtest(method1);
    REQUIRE(err1 == Approx(err2).epsilon(1e-8));
  }

}

TEST_CASE("Timer and tester"){
  //Times and compares the Jacobi algorithm with armadillo functions for different
  //values of N. Also tests that their first eigenvalues are always similar.
  double a = 0, b = 1;
  int Vchoice = 0;
  double conv = 1e-8;
  vec ctimeNs = logspace(1, 2, 17);
  double start1, end1, start2, end2;
  Jacobi_rotation timer;

  ofstream ofile;
  ofile.open("csv_files/timer.csv");
  ofile <<setw(15) << setprecision(8);
  ofile << "log10N," << "Jtime," << "atime," << "iters" <<endl;

  for (int i =0; i<17; i++){
    timer.initialize(a, b, ctimeNs(i), Vchoice);

    vec eigval;
    mat eigvec;
    start2 = clock();
    timer.eigarma(eigval, eigvec);
    end2 = clock();

    start1 = clock();
    timer.rotate(conv);
    end1 = clock();
    timer.rearrange();

    REQUIRE(timer.A(0,0) == Approx(eigval(0)).epsilon(1e-5));

    ofile << log10(ctimeNs(i)) << "," <<(end1-start1)/CLOCKS_PER_SEC<<"," << (end2-start2)/CLOCKS_PER_SEC<< "," << timer.iters << endl;

  }
  ofile.close();
}

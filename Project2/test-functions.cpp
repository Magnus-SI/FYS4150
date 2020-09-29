#include "catch.hpp"
#include "Jacobi_rotation.hpp"
#include <armadillo>
#include <cmath>

double V_0(double rho){
  return 0;
}

TEST_CASE("Testing eigenvalues/eigenvectors of Toeplitz matrix"){

    int Dim = 3;
    //    Set up the exact eigenvalues
    vec Exact(Dim);
    double pi = acos(-1.0);
    // Integration step length
    //note that example programs didn't use dim+1 here, why?
    double Step    = 1.0/(Dim+1);
    double DiagConst = 2.0 / (Step*Step);
    double NondiagConst =  -1.0 / (Step*Step);
    for(int i = 0; i < Dim; i++) {
      Exact(i) = DiagConst+2*NondiagConst*cos((i+1)*pi/(Dim+1));
    }

    //get numerical eigenvalues
    Jacobi_rotation my_solver;

    double a = 0;
    double b = 1;
    my_solver.initialize(a, b, Dim, V_0);
    double conv = 1e-8;
    my_solver.rotate(conv);
    my_solver.rearrange();
    vec EigvalueNum = my_solver.return_eig();
    mat eigvecNum = my_solver.V;
    mat Exactvec = my_solver.test_eig();

    //cout << EigvalueNum(0) << "yo"<< Exact(0) << endl;
    //cout << EigvalueNum(1) << "yo"<< Exact(1) << endl;
    //cout << EigvalueNum(2) << "yo"<< Exact(2) << endl;
    for(int i=0; i<Dim; i++){
      REQUIRE(EigvalueNum(i)==Approx(Exact(i)).epsilon(0.00000001));
      for(int j=0; j<Dim; j++){
        cout << Exactvec(i,j) << " " << eigvecNum(i,j) << endl;
        REQUIRE(abs(eigvecNum(i,j))==Approx(abs(Exactvec(i,j))).epsilon(0.00000001));
      }
    }
    //REQUIRE(EigvalueNum(0)==Approx(Exact(0)).epsilon(0.00000001));
    //REQUIRE(EigvalueNum(1)==Approx(Exact(1)).epsilon(0.00000001));
    //REQUIRE(EigvalueNum(2)==Approx(Exact(2)).epsilon(0.00000001));
}

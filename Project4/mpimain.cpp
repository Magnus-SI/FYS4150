#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "mpi.h"
#include "ising2d.cpp"
using namespace std;


void run_ising(stringstream& filedat, int L, double temp, double tol, int mcs_max, int seed){
  ising2D my_ising;
  my_ising.seed(seed);
  my_ising.initialize(L, temp, tol);
  for(int mcs=0; mcs<mcs_max; mcs++){
    my_ising.metropolis();
  }
  my_ising.write_mean(filedat);

}

stringstream file_stuff(int* Ls, double T, double tol, int mcs_max, int seed);

int main(int argc, char* argv[]){

  //Customizable parameters:
  int Tperproc = 2;
  double tol = 0.5;
  int mcs_max = 4000000;

  //MPI parameters
  int my_rank, numprocs;
  // MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  //Temperaturer and Length values
  double dT, Tmax, Tmin, T;
  //int L;
  Tmin = 2.2; Tmax = 2.35;
  dT = (Tmax - Tmin) / (Tperproc * numprocs - 1);
  int Ls[4] = {40, 60, 80, 100};

  ising2D my_ising;

  ofstream ofile;
  stringstream filedat;

  for (int i = my_rank*Tperproc; i < (my_rank+1)*Tperproc; i++){
    //cout << my_rank << " "<< i << endl;
    T = Tmin + i * dT;
    filedat = file_stuff(Ls, T, tol, mcs_max, i);

    #pragma omp critical
    ofile.open("data/T" + to_string(T) + "multiL.csv"); ofile << filedat.rdbuf(); ofile.close();
    }


  MPI_Finalize();
  return 0;
}

stringstream file_stuff(int* Ls, double T, double tol, int mcs_max, int seed){
  stringstream filedat;
  int L;
  filedat << "L,E,M,Cv,chi" << endl;
  filedat << setw(15) << setprecision(8);
  //cout << i << endl;
  //cout << T << endl;
  for (int j = 0; j<4; j++){
    L = Ls[j];
    cout << L << endl;
    filedat << L << ",";
    run_ising(filedat, L, T, tol, mcs_max, seed);
  }
  return filedat;
}

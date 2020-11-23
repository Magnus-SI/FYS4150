#include "ising2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "mpi.h"
#include "ising2d.cpp"
#include "time.h"

int main(int argc, char* argv[]){
  int Tperproc;
  double tol = 0.5;
  int mcs_max = 10000;
  double time1, time2;

  //MPI parameters
  int my_rank, numprocs;
  // MPI initializations
  time1 = clock();
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  Tperproc = 8/numprocs;

  //Temperaturer and Length values
  double dT, Tmax, Tmin, T;
  int L = atoi(argv[1]);
  //int L;
  Tmin = 2.0; Tmax = 2.3;
  dT = (Tmax - Tmin) / (Tperproc * numprocs - 1);


  ising2D my_ising;

  stringstream filedat;

  for (int i = my_rank*Tperproc; i < (my_rank+1)*Tperproc; i++){
    //cout << my_rank << " "<< i << endl;
    T = Tmin + i * dT;
    my_ising.seed(i);
    my_ising.initialize(L, T, tol);
    for(int mcs=0; mcs<mcs_max; mcs++){
      my_ising.metropolis();
    }
  }

  MPI_Finalize();

  //Write time taken to file:
  time2 = clock();
  ofstream ofile;
  ofile.open("data/timer.csv", ofstream::app);
  ofile << setw(15) << setprecision(8);
  ofile << L <<  ","<< numprocs << "," << (time2-time1)/CLOCKS_PER_SEC << endl;
  ofile.close();
  return 0;
}

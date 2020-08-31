#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

ofstream ofile;

//declare global variables for the constant values
double *a;
double *b;
double *c;
double *x;
double *u;
double *g;
double h;

void init(int N){
  double xmin = 0;
  double xmax = 1;
  h = (xmax - xmin)/(N+1);
  a = new double[N];
  b = new double[N];
  c = new double[N];
  x = new double[N];
  u = new double[N];
  g = new double[N];

  for (int i =0; i<N; i++){
    x[i] = xmin + h*(i+1);
    u[i] = 1 - (1 - exp(-10)) * x[i] - exp(-10 * x[i]);
    g[i] = 100*exp(-10 * x[i]) * pow(h, 2);
    //some if test can be added here for custom abc for testing purposes
    a[i] = -1.; c[i] = -1.;
    b[i] = 2.;

  }
}

void forward_general(double *gt, double *dt, int N){
  dt[0] = b[0];
  gt[0] = g[0];
  for (int i = 1; i<N; i++){
    dt[i] = b[i] - a[i-1]*c[i-1]/dt[i-1];
    gt[i] = g[i] - a[i-1]*gt[i-1]/dt[i-1];
  }
}

void forward_optim(double *gt, double *dt, int N){
  dt[0] = 2.;
  gt[0] = g[0];
  for (int i = 1; i<N; i++){
    dt[i] = 2. - 1/dt[i-1];
    gt[i] = g[i] + gt[i-1]/dt[i-1];
  }
}

void backward_general(double *v, double *gt, double *dt, int N){
  v[N-1] = gt[N-1]/dt[N-1];
  for (int i = N-1; i>=0; i--) {
    v[i] = (gt[i] - c[i]*v[i+1])/dt[i];
  }
}

void backward_optim(double *v, double *gt, double *dt, int N){
  v[N-1] = gt[N-1]/dt[N-1];
  for (int i = N-1; i>=0; i--) {
    v[i] = (gt[i] + v[i+1])/dt[i];
  }
}

int main(int argc, char* argv[]){
  int N = atoi(argv[1]);
  char *ofilename;
  ofilename = argv[2];

  init(N);      //initialize constants
  double *gt = new double[N];
  double *dt = new double[N];
  double *v = new double[N];
  forward_general(gt, dt, N);
  backward_general(v, gt, dt, N);

  ofile.open(ofilename);
  ofile << "x," << "u," << "v" << endl;
  for (int i = 0; i<N; i++){
    ofile << setw(15) << setprecision(8) << x[i] <<"," << u[i] << "," << v[i] << endl;
  }
  ofile.close();

  return 0;
}

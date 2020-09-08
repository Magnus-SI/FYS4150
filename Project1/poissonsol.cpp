#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "armadillo"
#include "time.h"
using namespace arma;
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

double compute_eps(int n, double *v_array){
	double eps_max, eps;
	eps_max = 0;
	for(int i=0; i<n; i++){
		eps = fabs((v_array[i] - u[i])/u[i]);
		if (eps > eps_max){
			eps_max = eps;
			}
		}
  cout << "n = " << n << endl;
	cout << "Max Relative error eps_max: " << eps_max << endl;
  cout << "log10(eps_max) = " << log10(eps_max) << endl;
  return log10(eps_max);
	}

vec LU_decomp(int n){
  mat A = mat(n,n,fill::eye);
  A(0,0) = 2;
  A(0,1) = -1;
  for(int i=1;i<n-1;i++){
    A(i,i) = 2;
    A(i,i-1) = -1;
    A(i,i+1) = -1;
  }
  A(n-1, n-1) = 2;
  A(n-1, n-2) = -1;
  mat L, U, P;
  lu(L,U,P,A);
  vec g_vec(g, n);
  g_vec = P*g_vec;
  vec y = solve(L, g_vec);
  vec x_vec = solve(U, y);
  return x_vec;
}

int main(int argc, char* argv[]){
  if (argc < 3){
    cout << "Call program with exponent and method choice (0 or 1)" <<endl;
    exit(0);
  }

  int max_exponent = atoi(argv[1]);
  int method = atoi(argv[2]);   //0 if optimized, 1 if general
  //ofilename = argv[2];
  double *time = new double[max_exponent];
  double *LUtime = new double[max_exponent];
  double *eps = new double[max_exponent];
  //for loop here
  for (int i = 1; i <=max_exponent; i++){

    int N = (int) pow(10.0, i);
    init(N);      //initialize constantsreturn
    double *gt = new double[N];
    double *dt = new double[N];
    double *v = new double[N];
    vec LU_sol;

    double start, end;
    start = clock();
    if (method){
      cout << "Using general method" << endl;
      forward_general(gt, dt, N);
      backward_general(v, gt, dt, N);
    }
    else{
      cout << "Using optimized method" << endl;
      forward_optim(gt, dt, N);
      backward_optim(v, gt, dt, N);
    }

    end = clock();
    string expstr = to_string(i);

    if (i<=5 && method == 1){
      double startLU = clock();
      LU_sol = LU_decomp(N);
      double endLU = clock();
      LUtime[i-1] = (endLU - startLU)/CLOCKS_PER_SEC;
    }
    else{
      LUtime[i-1] = 0;
    }



    if (i <= 3 && method == 1){
      string vdataname = "values";
      vdataname.append(expstr).append(".csv");
      ofile.open(vdataname);
      ofile << "x," << "u," << "v," << "LUv" << endl;
      for (int i = 0; i<N; i++){
        ofile << setw(15) << setprecision(8) << x[i] <<"," << u[i] << "," << v[i] << "," << LU_sol[i] << endl;
      }
      ofile.close();
    }

    time[i-1] = (end - start)/CLOCKS_PER_SEC;
    eps[i-1] = compute_eps(N, v);
  }

  string extradat = "epsclock";
  extradat.append(to_string(method));
  extradat.append(".csv");
  ofile.open(extradat);
  ofile << "expv," << "time," << "maxeps," << "LUtime" <<endl;
  for (int i = 0; i < max_exponent; i++){
    ofile << setw(15) << setprecision(8) << to_string(i+1) << "," << time[i] << "," << eps[i] << "," << LUtime[i] << endl;
  }

  return 0;
}

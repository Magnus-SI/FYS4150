#include "catch.hpp"
#include "diffusion.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

TEST_CASE("Periodic boundary test"){
  //Tests if the periodic boundary conditions on the x-axis work.

  Diffusion diff;
  int nx = 7;       //nx and ny includes edge points
  int ny = 5;
  double dx = 1;
  double dt = 0.2;
  diff.init2D(nx-2, ny-2, dx, dt, 1);

  REQUIRE(diff.xperiodic(nx-1, 0, 1, 0) == 0);
  REQUIRE(diff.xperiodic(nx-1, 0, -1, 0) == nx-2);
  REQUIRE(diff.xperiodic(nx-1, 0, 0, 2) == nx-1 + nx*2);
  REQUIRE(diff.xperiodic(0, 2, -1, 0) == 2*nx + nx-1);
}

TEST_CASE("Error function"){
  //Checks if the error functions returns no error if u is equal to uexact.
  Diffusion diff;
  int nx = 11;       //nx and ny includes edge points
  int ny = 11;
  double dx = 0.1;
  double dt = 0.2;
  diff.init2D(nx-2, ny-2, dx, dt, 1);

  diff.Q_0();
  for (int i = 0; i<10; i++){
    diff.EF2Dperx();
  }

  REQUIRE(diff.err2Dlocal(20, 1) == 0);
}

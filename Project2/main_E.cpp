#include "classtuff.hpp"
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "armadillo"
#include "time.h"
#include <stdio.h>

using namespace std;
using namespace arma;
double V(double rho){
    double omega = 1./20;
    // Using b since its defined as rho_max.
  return (omega*omega)*(rho*rho) + 1./rho;
}

int main(int argc, char const *argv[]) {

  //Define class object
  classtuff mysolver;

  int c_size = 100;
  //int c_size = 300;
  //int b=5;
  // eig 1: 0.352646 for N = 450 b =14.
  double a = 0;
  double b = 14;
  //Define matrix to solve Ax = lambda x
  //Initialize matrices
  mat A = mysolver.Initialize(a,  b,  V,  c_size);
  vec test_eigvals = mysolver.Jacobi_arm(A);
  mat qen = mysolver.Jacobi(A,1e-16, c_size);
  vec eigs = sort(qen.diag(),"descend");
  int F = eigs.n_elem;
  cout << "eig 1: " << eigs[F-1] << endl;
  return 0;
}

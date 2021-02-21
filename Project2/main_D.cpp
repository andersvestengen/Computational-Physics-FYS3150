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
  return rho*rho;
}

int main(int argc, char const *argv[]) {

  //Define class object
  classtuff mysolver;

  //int size = pow(10,atof(argv[1]));
  int c_size = 300;
  double a = 0;
  double b = 5;
  //Define matrix to solve Ax = lambda x
  //Initialize matrices
  mat A = mysolver.Initialize(a,  b,  V,  c_size);
  vec test_eigvals = mysolver.Jacobi_arm(A);
  mat qen = mysolver.Jacobi(A,1e-16, c_size);
  vec eigs = sort(qen.diag(),"descend");
  int F = eigs.n_elem;
  //cout << F << endl;
  cout << "eig 4: " << eigs[F-4] << endl;
  cout << "eig 3: " << eigs[F-3] << endl;
  cout << "eig 2: " << eigs[F-2] << endl;
  cout << "eig 1: " << eigs[F-1] << endl;
  //cout << eigs << endl;
  return 0;
}

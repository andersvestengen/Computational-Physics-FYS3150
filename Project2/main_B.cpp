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
ofstream ofile;
double V(double x){return 0;}
int main(int argc, char const *argv[]) {
  string outfilename;
  outfilename = "values.txt";
  ofile.open(outfilename);
  //Define class object
  classtuff mysolver;

  //int size = pow(10,atof(argv[1]));
  int c_size = 100;
  //Define matrix to solve Ax = lambda x
  double a = 0;
  double b = 1;
  //Initialize matrices
  mat A = mysolver.Initialize(a,  b,  V,  c_size);
  vec test_eigvals = mysolver.Jacobi_arm(A);
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);
  mat qen = mysolver.Jacobi(A,1e-16, c_size);
  cout << sort(qen.diag()) << endl;
  cout << sort(test_eigvals) << endl;
  for(int i=0; i < c_size; i++){
    if (min(mysolver.S.col(i)) >0){
      for(int k = 0;k< c_size;k++){
        double rho = k*1./(c_size+1);
        ofile << setprecision(15) << rho << " " <<mysolver.S(k,i)<< endl;
      }
    }
  }
  return 0;
}

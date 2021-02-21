#ifndef CLASSTUFF_HPP
#define CLASSTUFF_HPP
#include "armadillo"

using namespace std;
using namespace arma;


class classtuff {
private:
  int c_size;
  double s, c, theta, maxiter;
  int p; int q;

public:
  double maxoff;
  mat S,A,I;
  double eps;
  double rho_max ;
  mat Initialize(double a, double b, double V(double x), double c_size);
  mat Initialize_C(double a, mat ex, double rho_max);
  void Rotate(mat &A, mat &S, int &p, int &q, int n);
  mat Jacobi(mat A, double eps, int c_size);
  vec Jacobi_arm(mat T);
  void offdiag(mat A, int &p, int &q, int n, double &maxoff);
  double rho_func(double h, int i, double a);
  double tau(double theta);
  double tau(double a, double b, double c);
  double zero_test(mat it_A);
  mat returnS();
};

#endif

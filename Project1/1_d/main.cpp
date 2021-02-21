//  Project 1
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "armadillo"

using namespace std;
using namespace arma;

ofstream ofile;

inline double f(double x){return 100*exp(-10*x);}
inline double exactfunc(double x){return 1-(1-exp(-10))*x-exp(-10*x);}
inline double relativeerror(double sol, double exac){return abs((sol-exac)/exac);}

int main(int argc, char const *argv[]) {
  /* code */
  string outfilename;
  outfilename = "values.txt";
  ofile.open(outfilename);
  // Define matrix size

  int ex = atof(argv[1]);
  //Loop over all n to find max error for each n

  vec errlist(ex);
  vec nvals(ex);

  for (int p = 1; p <= ex; p++){

    int n = (int) pow(10,p);

    double h = 1./(n+1);
    cout << "Time step :" << h << endl;


    // Define vectors to solve equation Av = b
    vec v(n);
    vec x(n);
    vec g(n);
    vec gtilde(n);
    vec d(n);
    vec dtilde(n);
    vec e(n);
    vec sol(n);
    vec relerr(n);

    x = linspace(h,1-h,n);
    for (int i = 0; i<n; i++){ e(i) =1.;}
    for (int i =0; i<n; i++){ d(i) = -2.;}

    for (int i =0; i<n; i++){
      gtilde(i) = h*h*f(x(i));
    }

    //Forward Part
    for (int i = 1; i < n; i++)
    {
      gtilde(i) = gtilde(i) + ((double) i/(i+1)*gtilde(i-1));
    }
    //Backward Part
    v(n-1) = (double) gtilde(n-1)*(n)/(n+1);

    //exact solution
    for (int i = 0; i < n; i++){
      sol(i) = exactfunc(x(i));
    }


    cout << "Dimension of vectors:" << v.size() << endl;
    for (int i = n-2; i >= 0; i--)
    {
      v(i) = ((double) (i+1)/(i+2))*(gtilde(i)+v(i+1));
    }

    for (int i = 0; i < n; i++)
    {
      relerr(i) = relativeerror(v(i),sol(i));
    }

    //Print and write out log10 of relative error for a given n
    cout << log10(arma::max(relerr)) << endl;
    ofile << setprecision(15) << n <<" "<< log10(arma::max(relerr)) << endl;
  }

  ofile.close();

  return 0;
}

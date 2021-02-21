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
#include "time.h"
#include <stdio.h>

using namespace std;
using namespace arma;



ofstream ofile;

inline double f(double x){return 100*exp(-10*x);}
inline double exactfunc(double x){return 1-(1-exp(-10))*x-exp(-10*x);}

int main(int argc, char const *argv[]) {
  /* code */

  int ex = atof(argv[1]);

  for (int i=1; i<= ex; i++){
    string outfilename ;

    //Formats name of lists from amount of loops needed
    if (ex < 10){

      char filestuff[9];
      sprintf(filestuff, "valn%d.txt", i );
      outfilename = filestuff;
      ofile.open(outfilename);

    }
    else{

      char filestuff[10];
      sprintf(filestuff, "valn%d.txt", i );
      outfilename = filestuff;
      ofile.open(outfilename);

    }

    int n = (int) pow(10,i);

    double h = 1./(n+2);
    cout << "Time step :" << h << endl;
    cout << "Dimension of vectors:" << n << endl;

    // Define vectors to solve equation Av = b
    mat A = zeros<mat>(n,n);
    vec b(n);
    vec v(n);
    vec x(n);
    vec exac(n);
    vec Y(n);
    vec X(n);
    clock_t start, finish;
    start = clock();
    x = linspace(h,1-h,n);
    A(0,0) = -2; A(0,1) = 1; b(0) = -h*h*f(x(0));
    double hh = h*h;
    b(n-1) = -hh*f(x(n-1));
    for (int i =1; i<n-1; i++){
      b(i) = -hh*f(x(i));
      A(i,i-1) = 1;
      A(i,i) = -2;
      A(i,i+1) = 1;

    }
    A(n-1,n-1) = -2; A(n-1,n-2) = 1; A(n-2,n-1) = 1;

    mat L, U;
    lu(L,U,A);

    Y = solve(L,b);
    X = solve(U,Y);
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << setprecision(10) << "N="<< n<< ":  Time used  for computing=" << setprecision(3) << timeused  << endl;
    for (int i = 0; i<n; i++){
      exac(i) = exactfunc((i+1)*h);
      ofile << setprecision(15) << X(i) << " " << x(i) << " " << exac(i) << endl;
    }
    ofile.close();
  }

  return 0;
}

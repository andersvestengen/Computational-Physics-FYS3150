#include <iostream>
#include "armadillo"
#include "Black_scholes.hpp"
#include <new>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include "time.h"
#include <stdio.h>
#include <tuple>
#include <cmath>

using namespace arma;
using namespace std;
//This main is used to calculate option values for spesific S and \tau values
//and compare them to analytical results.
int main(int argc, char* argv[])
{
   Black_scholes SC;
   clock_t start, finish;
   start = clock();
   double T = 1; double X=0.75; int N=1e3;
   int print_per = 10;
   string filename="u.txt";double r = 0.04; double D=0.12; double sigma=0.4; double E=50;
   SC.Initialize(T,X,N,filename,r,D,sigma,E);
   //SC.D1d_explicit();
   SC.Crank_Nic(print_per);
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
return 0;
}

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
//This main is used to calculate data to visualize the greeks.
int main(int argc, char* argv[])
{
   Black_scholes SC;
   clock_t start, finish;
   start = clock();
   vec sigma_vec =linspace(0.1,1.,1e2); vec r_vec = linspace(0.01,0.5,1e2);
   string rfilename = "greeks_r.txt";string sfilename = "greeks_s.txt";
   SC.Greeks(sigma_vec,r_vec, rfilename, sfilename);
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
return 0;
}

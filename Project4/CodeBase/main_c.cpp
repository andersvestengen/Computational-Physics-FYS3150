#include <iostream>
#include "armadillo"
#include "Ising.hpp"
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
// polluting the namespaces
using namespace arma;
using namespace std;


int main(int argc, char* argv[])
{
   int mcs = 1e8;
   int L = 2;
   double param_1 = 0.;
   Ising Mcint1;
   clock_t start, finish;
   start = clock();
   double T = 1;
   Mcint1.Initialize(L, mcs,T,param_1);
   string filename = "MonteCarloRun.txt";
   Mcint1.init_output(filename);
   Mcint1.MonteCarloV2(filename);
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
   //Calculation and printing of analytical expressions for comparing.
   double Z = 12 + 2*exp(8/((double)T)) + 2*exp(-8/((double)T));
   double E_mean = (-16*exp(8/((double)T)) + 16*exp(-8/((double)T)))/((double) Z);
   double M_mean =  (8*exp(8/((double)T))+16)/((double)Z);
   double EE_mean = (2*(8*8)*exp(8/((double)T)) + 2*(8*8)*exp(-8/((double)T)))/((double) Z);
   double MM_mean = (2*(4*4)*exp(8/((double)T))+8*(2*2))/((double)Z);
   double c_v = (EE_mean-E_mean*E_mean)/((double)T*T);
   double X = (MM_mean-M_mean*M_mean)/((double)T);
   cout << "Mean value of E = " << E_mean<<endl;
   cout << "Mean value of M = " <<M_mean<<endl;
   cout << "C_v = " <<c_v<<endl;
   cout << "X = "<< X<<endl;

return 0;
}

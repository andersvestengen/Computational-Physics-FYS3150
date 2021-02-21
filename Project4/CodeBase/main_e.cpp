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
   Ising Mcint1;
   clock_t start, finish;
   start = clock();
   double L = 20.;
   double T = 1;
   int mcs_max = 1e5;
   int param_1 = 0.;
   Mcint1.Initialize(L, mcs_max,T,param_1);
   string filename = "MonteCarloRun.txt";
   Mcint1.init_output(filename);
   Mcint1.MonteCarloV1();
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
   string filename2 = "E.txt";
   int stabile_indx = 7e4;
   Mcint1.print_E_av(stabile_indx,filename2);
   Mcint1.calc_variance(stabile_indx);
   cout<< setprecision(4)<< "The variance of the energy is "<<Mcint1.m_variance << endl;
   cout<< setprecision(4)<< "The standard deviation of the energy is "<<sqrt(Mcint1.m_variance)<< endl;
return 0;
}

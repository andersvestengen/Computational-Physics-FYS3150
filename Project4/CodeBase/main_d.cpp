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
   int mcs_max = 1e4;
   string filename1 = "T_1_o.txt"; string filename2 = "T_2_o.txt";
   string filename3 = "T_1_n.txt"; string filename4 = "T_2_n.txt";

   Mcint1.Initialize(L, mcs_max,1.,1);Mcint1.init_output(filename1);Mcint1.MonteCarloV1();
   Mcint1.Initialize(L, mcs_max,2.4,1);Mcint1.init_output(filename2);Mcint1.MonteCarloV1();
   Mcint1.Initialize(L, mcs_max,1.,0.);Mcint1.init_output(filename3);Mcint1.MonteCarloV1();
   Mcint1.Initialize(L, mcs_max,2.4,0.);Mcint1.init_output(filename4);Mcint1.MonteCarloV1();
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
return 0;
}

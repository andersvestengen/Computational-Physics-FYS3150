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
#include <omp.h>
#include <tuple>
#include <cmath>
// polluting the namespaces
using namespace arma;
using namespace std;


int main(int argc, char* argv[])
{

   int mcs = 1e3; // Number of Monte Carlo simulations
   double init_temp = 2.1; // start temperature
   double final_temp = 2.4; // last temperature
   double t_step = 0.001; // step size for iterating 
   int param = 0; // ordered lattice structure (all pointing up) = 1, or = 0 for random ordering.

   string name[4] = {"MCL40.txt", "MCL60.txt", "MCL80.txt", "MCL100.txt"};
   int p = 0;

   omp_set_num_threads(8); // this number needs to be optimized for individual pc's !
   int iter = int( (final_temp - init_temp) / t_step );
   double start = omp_get_wtime();


   for(int L = 40; L < 101; L += 20){
      int k = 0;
     #pragma omp parallel for
     for(int i = 0; i <= iter; i++){
       double i_temp = (double) init_temp + i*t_step;
       Ising Mcint1;
       string filename = name[p];

       Mcint1.Initialize(L, mcs, i_temp, param);
       Mcint1.MonteCarloV2(filename);
       cout<< k << " / " << iter << endl;
       k+=1;
     }
     cout << "L = " << L << endl;
     p+=1;

     double finish = omp_get_wtime();
     double timeused = (double) (finish - start);
     cout << setprecision(10) << "Time used  for computing (Multithread) = " << timeused  << " Seconds"<<endl;
   }


return 0;
}

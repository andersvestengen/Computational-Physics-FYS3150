#include <iostream>
#include "armadillo"
#include "Finance.hpp"
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

   int L = 1024;
   double m_0 = 1;
   int mcs = 1e5;
   double tax_or_no = 0;
   double min_tax = 0;
   double alpha05 = 0.5;double alpha1 = 1.;double alpha15 = 1.5;double alpha2 = 2.;
   string filename2 = "V_vis.txt";
   string alpha_05="alpha_05.txt"; string alpha_1="alpha_1.txt";
   string alpha_15="alpha_15.txt";string alpha_2="alpha_2.txt";
   Finance Fc;
   double savings025 = 0;
   clock_t start, finish;
   start = clock();
   Fc.Initialize(mcs, L,m_0,filename2,tax_or_no,min_tax,savings025,alpha05);Fc.MonteCarlo();Fc.print_vec(alpha_05);
   Fc.Initialize(mcs, L,m_0,filename2,tax_or_no,min_tax,savings025,alpha1);Fc.MonteCarlo();Fc.print_vec(alpha_1);
   Fc.Initialize(mcs, L,m_0,filename2,tax_or_no,min_tax,savings025,alpha15);Fc.MonteCarlo();Fc.print_vec(alpha_15);
   Fc.Initialize(mcs, L,m_0,filename2,tax_or_no,min_tax,savings025,alpha2);Fc.MonteCarlo();Fc.print_vec(alpha_2);
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
return 0;
}

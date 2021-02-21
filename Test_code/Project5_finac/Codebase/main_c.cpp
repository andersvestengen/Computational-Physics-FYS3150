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
   int mcs = 1.5e5;
   int L = 500;
   double m_0 = 100;
   double tax_or_no = 0;
   double min_tax = 0;
   double alpha = 0;
   string filename1 = "Savings025.txt";
   string filename2 = "Savings05.txt";
   string filename3 = "Savings09.txt";
   Finance Fc;
   double savings025 = 0.25;double savings05 = 0.5;double savings09 = 0.9;
   clock_t start, finish;
   start = clock();
   Fc.Initialize(mcs, L,m_0,filename1,tax_or_no,min_tax,savings025,alpha,0);Fc.MonteCarlo();Fc.print_vec(filename1);
   Fc.Initialize(mcs, L,m_0,filename2,tax_or_no,min_tax,savings05,alpha,0);Fc.MonteCarlo();Fc.print_vec(filename2);
   Fc.Initialize(mcs, L,m_0,filename3,tax_or_no,min_tax,savings09,alpha,0);Fc.MonteCarlo();Fc.print_vec(filename3);
   finish = clock();
   double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
   cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
return 0;
}

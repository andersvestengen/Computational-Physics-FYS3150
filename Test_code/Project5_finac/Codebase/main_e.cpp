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
#include <omp.h>
//#include <Python/Python.h>

// polluting the namespaces
using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{
   int mcs = 2.5e5;
   int L = 500;
   double m_0 = 100;
   double tax_or_no = 0;
   double min_tax = 0;
   double alpha, gamma;
   string filename = "V_vis.txt";
   string filename1 = "Mon_vis.txt"; string filename6 = "Mon_vis6.txt";
   string filename2 = "Mon_vis2.txt"; string filename7 = "Mon_vis7.txt";
   string filename3 = "Mon_vis3.txt"; string filename8 = "Mon_vis8.txt";
   string filename4 = "Mon_vis4.txt"; string filename9 = "Mon_vis9.txt";
   string filename5 = "Mon_vis5.txt"; string filename10 = "Mon_vis10.txt";

   string files[10] = {filename1,filename2,filename3,filename4,filename5,filename6,
                       filename7, filename8, filename9, filename10};

   //omp_set_num_threads(8); // this number needs to be optimized for individual pc's !


   double savings025 = 0.25;
   double gammaval[5] = {0,1.0,2.0,3.0,4.0};
   double alphaval[2] = {1,2};

   //char filepy[] = "CodeBase/histo.py";
   //double start = omp_get_wtime();
   clock_t start, finish;

   int k = 0;
   //#pragma omp parallel for
   for(int i = 0; i < 2; i++){
     alpha = alphaval[i];


     for(int j = 0; j < 5; j++){
       gamma = gammaval[j];
       Finance Fc;
	     //Py_Initialize();

       start = clock();
       Fc.Initialize(mcs, L,m_0,filename,tax_or_no,min_tax,savings025,alpha,gamma);
       Fc.MonteCarlo();Fc.print_vec(files[k]);

       //double finish = omp_get_wtime();
       //double timeused = (double) (finish - start);
       //cout << setprecision(10) << "Time used  for computing (Multithread) = " << timeused  << " Seconds"<<endl;

       finish = clock();
       double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
       cout << setprecision(10) << "Time used  for computing (single thread) = " << timeused  << " Seconds"<<endl;
       
       // PyObject* fp = PyFile_FromString(filepy, "r");

	     //PyRun_SimpleFileEx(PyFile_AsFile(fp), filepy,1);

       //Py_Finalize();


       cout << "Gamma val = " << gamma << endl;
       k+=1;
     }
   }

return 0;
}

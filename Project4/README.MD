# Project 4

Here are the instructions for operating the numerical machinery.

## preamble:
Due to the nature of the solution there are different mains for different tasks. This is because our solver class is very general, and therefore requires specialized input for different solutions. Therefore, all code is inside the 'CodeBase' folder and you will execute the different aspects by editing and running the makefile. Its intended to run only one main_X at a time. We have not tested with multiple mains, but it should work.

To run main_f.x, you will need to have OpenMP installed. More information on OpenMP can be found at https://openmp.llvm.org/ and https://www.openmp.org/. 

When compiling with OpenMP, iOS differs from other unix based operating systems, as well as Windows. You will need to add or change the following flags in the Makefile to accomodate your operating system:

For Windows/Linux:
```
g++ -fopenmp -o test.cpp
```

For iOS:
```
g++ -Xpreprocessor -fopenmp -lomp -o test.cpp
```

### Makefile setup
The contents of the makefile look like this:
```
all: compile execute

compile:
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_c.x CodeBase/main_c.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_d.x CodeBase/main_d.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_e.x CodeBase/main_e.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -fopenmp -Wall -o CodeBase/main_f.x CodeBase/main_f.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/unit.x CodeBase/Unittest.cpp CodeBase/solver.cpp -larmadillo -llapack -lblas
	# Mac computers needs this instead: -Xpreprocessor -fopenmp -lomp
execute:
	#CodeBase/./main_c.x
	#CodeBase/./main_d.x
	#CodeBase/./main_e.x
	#CodeBase/./main_f.x
	#CodeBase/./unit.x
	#python3 CodeBase/EMC_visual.py
	#python3 CodeBase/E_histogram.py
	#python3 Codebase/T_critical.py
```
#### Instructions for running code
a makefile comment is a '#' character. To configure the Makefile you uncomment the files you would like to run. Do this under both the 'Compile' and 'Execution' sections. For example to run our implementation of c):

  - uncomment the line containing main_c under the 'compile' section
  - uncomment the line containing main_c.x under the 'execute' section
  
to run the makefile open a terminal and go to the 'Project4' folder. Type 'make all' and execute.


### Expected result for c):
  - Terminal statement with time used to run computation
  - Terminal statement with the mean value for energy, maked with E
  - Terminal statement with the mean value for magnetization, marked with M
  - Terminal statement with heat capacity, marked as C_v
  - Terminal statement with susceptibility, marked as X

### Expected result for d):
  - Terminal statement with time used to run computation
  - Plot showing mean magnetization as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 1. 
  - Plot showing mean energy as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 1. 
  - Plot showing accepted configurations as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 1. 
  
  - Plot showing mean magnetization as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 2. 
  - Plot showing mean energy as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 2. 
  - Plot showing accepted configurations as function of Monte Carlo cylces for random and ordered configuration at a given temperature T = 2. 
 
Expected result for e):
  - Terminal statement with time used to run computation
  - Terminal statement with variance of energy
  - Terminal statement with standard deviation of energy
  - Plot showing a histogram of probability as function of energy

### Expected result for f and g):
  - Plot showing the heat capacity as function of temperature for L = 40 with a univariate spline interpolation
  - Plot showing the heat capacity as function of temperature for L = 60 with a univariate spline interpolation
  - Plot showing the heat capacity as function of temperature for L = 80 with a univariate spline interpolation
  - Plot showing the heat capacity as function of temperature for L = 100 with a univariate spline interpolation
  - Plot showing the univariate spline interpolation functions for all four L grids as function of temperature
  
  - Terminal statement with residue for the four interpolation functions
  - Critical temperature for L = &infin;
  
  - Plot showing the susceptibility as function of temperature for L = 40 
  - Plot showing the susceptibility as function of temperature for L = 60 
  - Plot showing the susceptibility as function of temperature for L = 80 
  - Plot showing the susceptibility as function of temperature for L = 100 



#### General example of the code. 
Note that the number of threads possible for a given computer differs with the onboard CPU and you should check different 'omp_set_numthreads()' numbers before running large computations. In brief testing we found there to be as much as 5% difference in computation time. You also have to define the temperature domain, the step size and the number of Monte Carlo cycles. These should be omptimized to your computer, as an increase in Monte Carlo cycles and/or a decrease in temperature step size will significantly increase the computation time. In the for loop, you choose the different spin numbers, and how many different lattices you wish to look at. Also remember to include the OpenMP library, as shown below:
```
#include <omp.h>
int main(int argc, char* argv[])
{
   // Computational variables.
   
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


```


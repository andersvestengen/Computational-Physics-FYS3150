# Project 5 Black-Scholes

Here are the instructions for operating the numerical machinery.

## preamble:
Due to the nature of the solution there are different mains for different tasks. This is because our solver class is very general, and therefore requires specialized input for different solutions. All code is inside the 'CodeBase' folder and you will execute the different aspects by editing and running the makefile. Its intended to run only one main_X at a time. We have not tested with multiple mains, but it should work. This code was tested with Python version 3.7.6, Clang version 4.0.1, and GCC version 9.3.0. Tested on Windows and Mac OS.

### The following C++ and Python packages are required:  
Python:  
  
	numpy 1.19.4  
	matplotlib 3.3.3  
	scipy 1.5.4  
C++:  
  
	armadillo 9.800.4  
	blas 3.9.0  
	lapack 3.9.0  
  
### Makefile setup
The contents of the makefile look like this:
```
all: compile execute

compile:
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_b.x CodeBase/main_b.cpp CodeBase/Black_scholes.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_d.x CodeBase/main_d.cpp CodeBase/Black_scholes.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_c.x CodeBase/main_c.cpp CodeBase/Black_scholes.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/unit.x CodeBase/Unittest.cpp CodeBase/Black_scholes.cpp -larmadillo -llapack -lblas

execute:
	#CodeBase/./main_b.x
	#CodeBase/./main_c.x
	#CodeBase/./main_d.x
	#CodeBase/./unit.x

	#python3 CodeBase/visualize.py 
	#python3 CodeBase/visualize_greeks.py
	#python3 CodeBase/akerbp.py
```

#### Instructions for running code
a makefile comment is a '#' character. To configure the Makefile you uncomment the files you would like to run. Do this under both the 'Compile' and 'Execution' sections. For example to run our implementation of b):

  - Uncomment the line containing main_b under the 'compile' section
  - Uncomment the line containing main_b.x under the 'execute' section
  
to run the makefile open a terminal and go to the 'Project5' folder. Type 'make all' and execute.

you should run all three main cpp files, and then the python files.

### Expected result for b) and c): (also needs to run visualize.py and visualize_greeks.py)
 - Terminal statement with time used to run computation
 - Plot showing numerical solution to Black-Scholes equation for multiple timesteps
 - Plot showing analytical solution to Black-Scholes equation for multiple timesteps
 - Plot showing absolute difference between numerical and analytical solution to Black-Scholes equation for multiple timesteps
 - Plot showing numerical solution to greek Delta for multiple timesteps
 - Plot showing numerical solution to greek gamma for multiple timesteps
 - Plot showing numerical solution to greek Vega for multiple timesteps
 - Plot showing numerical solution to greek Rho for multiple timesteps
 - Plot showing numerical solution to greek Theta for multiple timesteps
 - Plot showing analytical solution to greek Delta for multiple timesteps
 - Plot showing analytical solution to greek gamma for multiple timesteps
 - Plot showing analytical solution to greek Theta for multiple timesteps
 - Plot showing analytical solution to greek Vega for multiple timesteps
 - Plot showing analytical solution to greek Rho for multiple timesteps

### Expected result for d): (also needs to run akerbp.py)
 - Plot showing Aker BP option values vs our numerical models estimated option values.



#### General example of the code. 
Remember to include necessary libraries
```
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

```
Note that the Initializer needs Black-Scholes variables E, sigma, D, and r  
Furthermore, X defines the interval of S values, T is total time, and N is the number of points in the interval.


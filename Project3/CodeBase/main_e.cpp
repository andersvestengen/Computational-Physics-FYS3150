#include "object.hpp"
#include "solving.hpp"
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "armadillo"
#include "time.h"
#include <stdio.h>
#include <tuple>
#include <math.h>
#include <typeinfo>
using namespace std;
using namespace arma;

int main() {

  int Dimension = 3;


  int IntegrationPoints = 1000;
  double FinalTime = 50.0;

  double earth_mass = 3.003e-6;
  double sun_mass = 1.0;

  double beta = 3;
  int alpha =0;
  int fixed =1;
  cout << "Beta " << beta << endl;

  object planet1(earth_mass,1.,0.0,0.0,0,5,0);
  object planet2(sun_mass, 0,0,0,0,0,0);

  solving binary_verlet;
  binary_verlet.add(planet1); binary_verlet.add(planet2);
  binary_verlet.VelocityVerlet(Dimension,IntegrationPoints,FinalTime, beta,fixed,alpha);

  return 0;
}

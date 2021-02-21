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


  int IntegrationPoints = 100000;
  double FinalTime = 18.0;

  double earth_mass = 3.003e-6;
  double sun_mass = 1.0;
  double jup_mass = 0.0009543;
  int fixed = 0;
  double beta = 2;
  int alpha =0;
  cout << "Beta " << beta << endl;

  object planet1(earth_mass,1.,0.0,0.0,0,6.3,0);
  object planet2(jup_mass,5.2,0.0,0.0,0, 3, 0);
  object planet3(sun_mass, 0,0,0,0,0,0);

  solving binary_verlet;
  binary_verlet.add(planet1); binary_verlet.add(planet2);binary_verlet.add(planet3);
  binary_verlet.VelocityVerlet(Dimension,IntegrationPoints,FinalTime, beta,fixed,alpha);


  return 0;
}

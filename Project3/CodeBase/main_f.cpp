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


  int IntegrationPoints = 10000;
  double FinalTime = 200.0;
  double earth_mass = 3.003e-6;
  double sun_mass = 1.0;
  vec pos(3); pos = {1,0,0};
  double r = sqrt( (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2])  );
  int fixed = 1;
  double beta = 2;
  int alpha =0;
  cout << "Beta " << beta << endl;

  object planet1(earth_mass,1.,0.0,0.0,0,8.88577,0);
  object planet2(sun_mass, 0,0,0,0,0,0);

  solving binary_verlet;
  binary_verlet.add(planet1); binary_verlet.add(planet2);

  double vesc = sqrt(2*binary_verlet.G/(r));

  cout << "V_esc = " << vesc << endl;
  binary_verlet.VelocityVerlet(Dimension,IntegrationPoints,FinalTime, beta,fixed,alpha);

  return 0;
}

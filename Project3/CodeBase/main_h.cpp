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


  int IntegrationPoints = 50000;
  double FinalTime = 1;
  double earth_mass = 3.003e-6;
  double sun_mass = 1.0;
  double jup_mass = 0.0009543;
  double mar_mass = 3.213e-7;
  double sat_mass = 0.0002857;
  double nep_mass=0.00005149;
  double ven_mass = 0.000002447;
  double uran_mass = 0.00004365;
  double merc_mass = 1.651e-7;
  double moon_mass = 3.69396868e-8;
  int fixed =1;
  int alpha =0;
  double beta = 2;
  cout << "Beta " << beta << endl;
  //When initializing the planets; make sure to initialize the sun last, and the planet furthest away from the sun second last.
  //Also add earth first.
  //Earth
  object planetearth(earth_mass,9.128884513088843*1e-01,3.928032801600736*1e-01,6.577938183713410*1e-05,-6.957269992142644*1e-03*365, 1.579734315560513*1e-02*365, -2.582593092148153*1e-07*365);
  //Earth moon
  object moon(moon_mass,9.105494529544440E-01,3.923421533389641E-01 ,2.651440460123353E-04,-6.823266197230106E-03*365,1.517678655175216E-02*365,-2.024433736638003E-05*365);
  //Jupiter
  object planetjup(jup_mass,2.556653950007264,-4.428596022378350,-3.882840438937561*1e-02,6.442741439253338*1e-03*365, 4.130146620372741*1e-03*365, -1.612738541610256*1e-04*365);
  //Venus
  object planetvenus(ven_mass,-4.150243727322463e-02,7.249838052804233e-01,1.199344556569153e-02,-2.027816182702826e-02*365,-1.110945668283595e-03*365,1.154794297386271e-03*365);
  //Mercury
  object planetmercury(merc_mass,2.643154303593134e-01, -3.157121726001744e-01,-5.104028200574216e-02 ,1.594913945525348e-02*365,1.943209256008390e-02*365,1.248175570907458e-04*365);
  //Sun
  object sun(sun_mass, -6.107925513172998*1e-03,6.420679726598624*1e-03,8.893727401374147*1e-05,-7.280593132276730*1e-06*365,-5.090234498858063*1e-06*365,2.181619304215098*1e-07*365);
  //Mars
  object planetmar(mar_mass,1.308532937474490E+00,5.379930853282286E-01,-2.102141505220228E-02,-4.717434819666467E-03*365,1.416299295671162E-02*365,4.126758221879923E-04*365);
  //Saturn
  object planetsat(sat_mass,5.144031493282123E+00,-8.563703043404734E+00,-5.588847142609978E-02,4.471296615234172E-03*365,2.858910024157302E-03*365,-2.280276562541721E-04*365);
  //Neptun
  object planetnep(nep_mass,2.941228501033349E+01,-5.465078789506711E+00,-5.652937055144212E-01,5.521932160701419E-04*365,3.104748042343270E-03*365 ,-7.649605995032246E-05*365);
  //Uranus
  object planeturanus(uran_mass,1.555699248858926e+01, 1.221910146632133e+01,-1.561607729538915e-01,-2.458296642987987e-03*365,2.909798437857242e-03*365,4.256035053183806e-05*365);

  solving binary_verlet;
  binary_verlet.add(planetearth);binary_verlet.add(moon);
  binary_verlet.add(planetmercury);binary_verlet.add(planetvenus);
  binary_verlet.add(planetmar); binary_verlet.add(planetjup);
  binary_verlet.add(planetsat); binary_verlet.add(planetnep);
  binary_verlet.add(planeturanus); binary_verlet.add(sun);
  binary_verlet.VelocityVerlet(Dimension,IntegrationPoints,FinalTime, beta, fixed,alpha);

  return 0;
}

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
//#include "classtuff.cpp"
#include "solver.hpp"
#include "object.hpp"
#include "armadillo"
#include <cstdio>
#include <math.h>

using namespace arma ;

TEST_CASE( "Check for errors in code" ) {

  // Make our test system
  //Add objects you wish to look at
  int IntegrationPoints = 50000;
  double FinalTime = 10;
  int Dimension = 3;
  double beta = 2;
  int fixed =0;
  double alpha =0;
  double earth_mass = 3.003e-6;
  double sun_mass = 1.0;

  object planetearth(earth_mass,9.128884513088843*1e-01,3.928032801600736*1e-01,6.577938183713410*1e-05,-6.957269992142644*1e-03*365, 1.579734315560513*1e-02*365, -2.582593092148153*1e-07*365);
  object sun(sun_mass, -6.107925513172998*1e-03,6.420679726598624*1e-03,8.893727401374147*1e-05,-7.280593132276730*1e-06*365,-5.090234498858063*1e-06*365,2.181619304215098*1e-07*365);

  solving unittest_solv;
  unittest_solv.add(planetearth); unittest_solv.add(sun);


  //Run simulation to aquire data
  unittest_solv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime, beta, fixed,alpha);


  int *integrationpoints = nullptr;
  float *number_o_planet = nullptr;
  float *t = nullptr;

  integrationpoints = new int;
  number_o_planet = new float;
  t = new float;

  double *kin, *pot, *l, *x, *y, *z ;
  const char* planet_info = "Planets_pos.txt";
  FILE *fp_init = fopen(planet_info, "r"); //Open file to read, specified by "r".
  fscanf(fp_init, "%d %e", integrationpoints, number_o_planet);

  kin = new double[IntegrationPoints];
  pot = new double[IntegrationPoints];
  l = new double[IntegrationPoints];
  x = new double[IntegrationPoints];
  y = new double[IntegrationPoints];
  z = new double[IntegrationPoints];

  //Waste values for the sun position, just needs to be updated, but are not used in any way
  double *sunx = nullptr; double *suny = nullptr; double *sunz = nullptr;
  sunx = new double[IntegrationPoints]; suny = new double[IntegrationPoints];
  sunz = new double[IntegrationPoints];

  //Read info from file
  int j = 0; int k = 0;
  int o = 0; int u = 0;
  for (int i = 0; i < (3)*IntegrationPoints; i++){

    if (j == 1){
      fscanf(fp_init, "%lf %lf %e %lf", &pot[k], &kin[k], t, &l[k]);
      k += 1;
      j+=1;
    }
    else if (j == 2){
      fscanf(fp_init, "%lf %lf %lf", &sunx[u], &suny[u], &sunz[u] );
      u += 1;
      j = 0;
    }
    else if (j == 0){

      fscanf(fp_init,"%lf %lf %lf ", &x[o], &y[o], &z[o]);
      o += 1;
      j += 1;
    }
  }

  fclose(fp_init);


  SECTION("Check conservation of energy"){
    //Now sum up energy

    double toten[IntegrationPoints];
    for (int i = 0; i < IntegrationPoints; i++){
      toten[i] = kin[i] + pot[i];
    }

    //Find max total energy value for each half of the array.

    double maxval = toten[0];
    for (int i = 1; i < IntegrationPoints/2; i++){
      if (fabs(toten[i]) > fabs(maxval)){
        maxval = toten[i];
      }
    }

    double maxval2 = toten[IntegrationPoints/2];
    for (int k = IntegrationPoints/2; k < IntegrationPoints; k++){
      if (fabs(toten[k]) > fabs(maxval2)){
        maxval2 = toten[k];
      }
    }
    // Checking if the difference between the peak in max-value in the totale-energy
    // is relativley equal for both halves. Therby making the energy periodic,
    // and therfore constant for each period.
    REQUIRE( fabs(maxval - maxval2) < fabs(maxval2)/150 );
  }
  SECTION("Check conservation of angular momentum"){
    REQUIRE(fabs(l[0]-l[int(IntegrationPoints/FinalTime)]) <l[0]/100);
  }

  SECTION("Checking if the earth is still in orbit "){
    //Calculate absolute relative distance for planet in beginning and end of simulation
    double xrel_start = x[0] - sunx[0]; double xrel_end = x[IntegrationPoints-1] - sunx[IntegrationPoints-1];
    double yrel_start = y[0] - suny[0]; double yrel_end = y[IntegrationPoints-1] - suny[IntegrationPoints-1];
    double zrel_start = z[0] - sunz[0]; double zrel_end = z[IntegrationPoints-1] - sunz[IntegrationPoints-1];

    double posrel_start = sqrt(xrel_start*xrel_start + yrel_start*yrel_start + zrel_start*zrel_start);
    double posrel_end = sqrt(xrel_end*xrel_end + yrel_end*yrel_end + zrel_end*zrel_end);
    REQUIRE( posrel_end < 2*posrel_start );
  }
}

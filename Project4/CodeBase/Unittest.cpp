#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
//#include "classtuff.cpp"
#include "armadillo"
#include <cstdio>
#include <math.h>
#include "Ising.hpp"

using namespace arma ;

TEST_CASE( "Check for errors in code" ) {
  Ising Mcint1;
SECTION("Check if the spins ar randomly ordered:"){
  //Checks if the spins are randomly ordered by checking if the avarage spin
  //converges towards 0.
  int mcs =1e2; double T =2.4; double param_1 = 0.;  int L_1 =10;int L_2 = 20;
  Mcint1.Initialize(L_1, mcs,T, param_1);
  int sum_1 = 0;
  for(int i = 0; i <L_1*L_1; i++){sum_1 += Mcint1.m_smatrix[i];}
  double sum_L_10 = sum_1/((double) L_1*L_1);
  Mcint1.Initialize(L_2, mcs,T, param_1);
  int sum_2=0;
  for(int i = 0; i <L_2*L_2; i++){sum_2 += Mcint1.m_smatrix[i];}
  double sum_L_20 = sum_2/((double)L_2*L_2);
  REQUIRE(fabs(sum_L_20) < fabs(sum_L_10));
  REQUIRE(0.0001<fabs(sum_L_20));
}
SECTION("Check if the Variance increases with the temperature:"){
  int mcs =1e5; double T_1 =1; double T_2 = 2.4; double param_1 = 0.;  int L =20;
  int stabile_indx = 7e4;
  Mcint1.Initialize(L, mcs,T_1, param_1);
  Mcint1.MonteCarloV1();
  Mcint1.calc_variance(stabile_indx);
  double variance_1 = Mcint1.m_variance;
  Mcint1.Initialize(L, mcs,T_2, param_1);
  Mcint1.MonteCarloV1();
  Mcint1.calc_variance(stabile_indx);
  double variance_2 = Mcint1.m_variance;
  REQUIRE(fabs(variance_1) < fabs(variance_2));
}
SECTION("Check if numerical values for mean energy and magnetization match analytical up to a certain thershold"){
  int mcs = 2e6;int L = 2;double param_1 = 0.;double T = 1;
  Mcint1.Initialize(L, mcs,T,param_1);
  string filename = "test.txt";
  Mcint1.init_output(filename);
  Mcint1.MonteCarloV2();
  //Calculation and printing of analytical expressions for comparing.
  double Z = 12 + 2*exp(8/((double)T)) + 2*exp(-8/((double)T));
  double E_mean = (-16*exp(8/((double)T)) + 16*exp(-8/((double)T)))/((double) Z);
  double M_mean =  (8*exp(8/((double)T))+16)/((double)Z);
  REQUIRE(fabs(E_mean-Mcint1.m_Etotal_average)<fabs(E_mean)/((double)1000.));
  REQUIRE(fabs(M_mean-Mcint1.m_Mabstotal_average)<fabs(M_mean)/((double)1000.));
}

}

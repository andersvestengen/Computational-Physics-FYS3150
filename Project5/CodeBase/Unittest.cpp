#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
//#include "classtuff.cpp"
#include "armadillo"
#include <cstdio>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <string>
#include "Black_scholes.hpp"
#include <functional>

using namespace arma ;
using namespace std ;

vector<vector<double> > readMatrix(string filename){

  // File pointer
  fstream fin;
  // Open an existing file
  fin.open(filename, ios::in);
  // Read the Data from the file
  // as String Vector
  vector<double> S;
  vector<double> time;
  string line, word, temp;

  int count = 0;
  vector<vector<double> > matrix; //Matrix to hold option value for each time step

  while (fin >> temp){

    getline(fin, line);

    // used for breaking words
    stringstream s(line);

    std::string resultstr = s.str();
    const char* cstr2 = resultstr.c_str();

    vector <double> VperS;

    char* line = (char*) cstr2;
    char* pch = strtok(line, " ");
    if (count == 0){  //If count == 0, line in file is stock price
      while (pch != NULL){
        S.push_back( std::atof(pch) );
        pch = strtok (NULL, " ");
      }
      count += 1;
    }
    else{ //Time and option value
      int hello = 0;
      while (pch != NULL){
        if ( hello == 0 ){
          time.push_back( std::atof(pch) );
          pch = strtok (NULL, " ");
          hello += 1;
        }
        VperS.push_back( std::atof(pch) );
        pch = strtok (NULL, " ");
      }
      matrix.push_back(VperS);
    }
  }

  return matrix;
}

TEST_CASE( "Check for errors in code" ) {
    Black_scholes SC;

  SECTION("Check if value of option is negative at any time"){
    double T = 1; double X=0.75; int N=1e4;
    string filename="u.txt";
    int print_per =10;
    double r = 0.04; double D=0.12;
    double sigma=0.4; double E=50;
    SC.Initialize(T,X,N,filename,r,D,sigma,E);
    SC.Crank_Nic(print_per);


    vector<vector<double> > matrix = readMatrix(filename);

    vector<double> compare_vec;
    for (int i = 0; i < N; i++){
      compare_vec.push_back((double) 0);
    }

    bool ok;
    int number = 0;
    for (int i = 0; i < 11; i++){
        ok = compare_vec <= matrix[i];
        number += ok;
    }

    REQUIRE(number == 11);

  }
  SECTION("Check if Cranck-Nicholsen is indipendent of step-size"){
    double T = 1; double X=0.75;
    double r = 0.04; double D=0.12;
    double sigma=0.4; double E=50;
    int print_per = 10;
    int N1=1e3;
    string filename1="N1e3.txt";
    SC.Initialize(T,X,N1,filename1,r,D,sigma,E);
    SC.Crank_Nic(print_per);
    vector<vector<double> > matrix1 = readMatrix(filename1);

    int N2=1e4;
    string filename2 ="N1e4.txt";
    SC.Initialize(T,X,N2,filename2,r,D,sigma,E);
    SC.Crank_Nic(print_per);

    vector<vector<double> > matrix2 = readMatrix(filename2);
    double sum1=0;double sum2=0;
    int ran_index =10*((int)rand()) / ((int)RAND_MAX);

    for(int i=0;i<N1;i++){
      sum1 += matrix1[ran_index][i];
    }
    for(int i=0;i<N2;i++){
      sum2 += matrix2[ran_index][i];
    }
    sum1/=(double)N1;
    sum2/=(double)N2;
    double tresh = sum2/((double)100);
    REQUIRE(fabs(sum1-sum2)<tresh);

  }
  SECTION("Check that option value increases for increased volatility in early time step"){
    Black_scholes SC;
    //Find random number between 0.1 and 0.5
    double Rnum = ((double)rand()) / ((double)RAND_MAX) / 2.0 + 0.1 ;

    double T = 1; double X=0.75; int N=1e4;
    int print_per = 10;
    string filename="u.txt";
    string filename2 = "u2.txt";
    double r = 0.04; double D=0.12;
    double sigma=Rnum; double sigma2 = Rnum*1.5;
    double E=50;
    //First volatility
    SC.Initialize(T,X,N,filename,r,D,sigma,E);
    SC.Crank_Nic(print_per);

    //Second volatility
    SC.Initialize(T,X,N,filename2,r,D,sigma2,E);
    SC.Crank_Nic(print_per);



    vector<vector<double> > matrix = readMatrix(filename);
    vector<vector<double> > matrix2 = readMatrix(filename2);
    bool ok;
    for (int i = 0; i < N/10; i++){
      ok = matrix[10][i] < matrix2[10][i];
      if(ok == 0){
        ok = 0;
        break;
      }
    }
    REQUIRE(ok > 0);

  }
}

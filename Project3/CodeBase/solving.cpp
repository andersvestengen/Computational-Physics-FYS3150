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
#include <cmath>

solving::solving(){
    total_planets = 0;
    total_mass = 0;
    G = 4*M_PI*M_PI;
}

void solving::add(object newplanet){
  total_planets += 1;
  total_mass += newplanet.mass;
  all_planets.push_back(newplanet);
}

double ** solving::setup_matrix(int height,int width){   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void solving::delete_matrix(double **matrix){   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void solving::delete_matrix3d(double ***matrix,int Integration_points){   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_planets; i++)
      for (int k = 0; k<Integration_points;k++)
        delete [] matrix[i][k];
}

double solving::angularmomentum(object &current){
    double spinvec[3];
    double pos[3]; pos[0] = current.position[0];
    pos[1] = current.position[1]; pos[2] = current.position[2];
    double vel[3]; vel[0] = current.velocity[0];
    vel[1] = current.velocity[1]; vel[2] = current.velocity[2];
    spinvec[0] = pos[1] * vel[2] - pos[2] * vel[1];
    spinvec[1] = pos[0] * vel[2] - pos[2] * vel[0];
    spinvec[2] = pos[0] * vel[1] - pos[1] * vel[0];

    double l = 0;
    for (int i = 0; i<3; i++){
      l += spinvec[i]*spinvec[i];
    }
    l = sqrt(l);
    return l;
    }

void solving::peri(double &thetha,double mon,double tue,double wen,double x_p,double y_p){
  if(mon>tue&&wen>tue){
    thetha = atan2(y_p,x_p);
    std::cout << thetha<<endl;
  }
}
void solving::GravitationalForce(object &current,object &other,double &Fx,double &Fy,double &Fz, double beta,int alpha){   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    double R = pow(r,beta);
    double l = angularmomentum(current);
    double ex;
    if(alpha == 1){
      ex = (double)( 1. + 3.*l*l/(R*cc));
    }
    else{
      ex = 1;
    }
    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/(r*R)*ex;
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/(r*R)*ex;
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/(r*R)*ex;
}

void solving::VelocityVerlet(int dimension, int integration_points, double final_time, double beta, int fixed, int alpha){
  double **acceleration = setup_matrix(total_planets, 3);
  double **acceleration_next = setup_matrix(total_planets, 3);
  double t = 0;
  double Fx, Fy, Fz, Fxnew, Fynew, Fznew, Pot, Kin,mon=0,tue=0,thetha,x_p,y_p;
  double dA1=0;
  double dA2=0;
  double h = final_time/((double) integration_points);
  double hh = h*h;
  double cm[3];
  double dt[4];
  dt[0] = 0.19*final_time;
  dt[1]=0.2*final_time;
  dt[2]=0.69*final_time;
  dt[3] = 0.7*final_time;
  if(fixed==1){center_of_mass(cm,  dimension);}
  std::ofstream ofile;
  std::string outfilename = "Planets_pos.txt";
  ofile.open(outfilename);
  ofile << integration_points << " " << total_planets<< " "<< endl;

  while (t < final_time){
    //Loop over all planets
    for (int nr = 0; nr < total_planets; nr++){

      object &current = all_planets[nr];

      Fx = Fy = Fz = Fxnew = Fynew = Fznew = Pot = Kin = 0;

      //Find forces on other planets
      for (int nr2 = 0; nr2 < total_planets; nr2++){
        if (nr2 != nr && nr < total_planets-1){
          object &other = all_planets[nr2];
          GravitationalForce(current, other, Fx, Fy, Fz, beta,alpha);
          PotentialEnergySystem(current, other, Pot);
          KineticEnergySystem(current, Kin);
          if(nr ==0){
            Delta_A(current,t, dt, dA1, dA2,h);
        }
      }
      else if(fixed==1 && nr==total_planets-1&&nr!=nr2){
        object &other = all_planets[nr2];
        GravitationalForce(current, other, Fx, Fy, Fz, beta,alpha);
      }
      }
      acceleration[nr][0] = Fx/current.mass;
      acceleration[nr][1] = Fy/current.mass;
      acceleration[nr][2] = Fz/current.mass;

      //Calculate position for each planet
      for (int k = 0; k<dimension; k++){
        current.position[k] += current.velocity[k]*h + 0.5*hh*acceleration[nr][k];
      }
      if(fixed==1){for(int j =0;j<dimension;j++){current.position[j]-=cm[j];}}
      //Loop again over all other planets
      for (int nr2 = 0; nr2 < total_planets; nr2++){
        if (nr2 != nr && nr < total_planets-1){
          object &other = all_planets[nr2];
          GravitationalForce(current, other, Fxnew, Fynew, Fznew, beta,alpha);
        }
        else if(fixed==1 && nr==total_planets-1 && nr!=nr2){
          object &other = all_planets[nr2];
          GravitationalForce(current, other, Fxnew, Fynew, Fznew, beta,alpha);
          }
      }
      acceleration_next[nr][0] = Fxnew/current.mass;
      acceleration_next[nr][1] = Fynew/current.mass;
      acceleration_next[nr][2] = Fznew/current.mass;

      //calculate velocity for planets
      for (int y = 0; y < dimension; y++){
        current.velocity[y] += 0.5*h*(acceleration[nr][y] + acceleration_next[nr][y]);
      }
      //Calculate kinetic energy for current object

      print_to_file(all_planets[nr].position, dimension, ofile);

      if(nr==0&&t>0.997*final_time){
      object &other = all_planets[-1];
      double wen = current.distance(other);
      peri(thetha,mon,tue,wen,x_p,y_p);
      x_p = current.position[0];
      y_p = current.position[1];
      mon = tue;
      tue = wen;
      }

      //Calculate angular momentum for print
      double l_print = angularmomentum(current);
      if (nr==0) {
        print_energi(Pot, Kin, t, l_print,ofile);
      }
    }
    t+= h;
    center_of_mass(cm,dimension);


  }
  // Clear memory
  ofile.close();
  delete_matrix(acceleration);
  delete_matrix(acceleration_next);
  std::cout<<setprecision(15)<<"Perihilion angle ="<<thetha*206264.806<<endl;
  std::cout<<"Area of the first time interval is ="<<dA1<<endl;
  std::cout<<"Area of the second time interval is ="<<dA2<<endl;
}

void solving::Forward_Euler(int dimension, int integration_points, double final_time,double beta,int alpha)
{
  double **acceleration = setup_matrix(total_planets, 3);
  double t = 0;
  double Fx, Fy, Fz;
  double h = final_time/((double) integration_points);
  std::ofstream ofile;
  std::string outfilename = "values_Euler.txt";
  ofile.open(outfilename);
  while (t < final_time){
    //Loop over all planets
    for (int nr = 0; nr < total_planets; nr++){

      object &current = all_planets[nr];

      Fx = Fy = Fz;

      //Find forces on other planets
      for (int nr2 = 0; nr2 < total_planets; nr2++){
        if (nr2 != nr && nr < total_planets-1){
          object &other = all_planets[nr2];
          GravitationalForce(current, other, Fx, Fy, Fz, beta,alpha);
      }
      }
      acceleration[nr][0] = Fx/current.mass;
      acceleration[nr][1] = Fy/current.mass;
      acceleration[nr][2] = Fz/current.mass;
      //Calculate position for each planet
      for (int k = 0; k<dimension; k++){
        current.position[k] += current.velocity[k]*h;
      }
      for (int k = 0; k<dimension; k++){
        current.velocity[k] += h*acceleration[nr][k];
      }
      if(nr ==0){print_to_file(all_planets[nr].position, dimension, ofile);}
    }
    t+= h;
  }
  // Clear memory
  ofile.close();
  delete_matrix(acceleration);
}
void solving::print_to_file(double planets[3],int dimension, std::ofstream &ofile){
  ofile << std::setprecision(16)<< planets[0] << " "<< planets[1] << " "<< planets[2] <<endl;
}

void solving::print_energi(double &Pot, double &Kin, double &t, double &l,std::ofstream &ofile){
  ofile <<std::setprecision(20)<<Pot<< " " << Kin <<" " << t<< " " << l << endl;
}

void solving::PotentialEnergySystem(object &current, object &other, double &Pot){
  double r = current.distance(other);
  if ( r != 0){
    Pot = -this->G*current.mass*other.mass/r;
  }
  else {
    Pot = 0;
  }
}

void solving::center_of_mass(double cm[3], int dimension){
  cm[0]=cm[1]=cm[2]=0;
  for(int j =0;j<total_planets;j++){
    object &current = all_planets[j];
    for(int i=0;i<dimension;i++){
      cm[i] += (double)current.mass*current.position[i]/total_mass;
    }
  }
}

void solving::KineticEnergySystem(object &current, double &Kin){
  double velo2 = 0;
  for (int i = 0; i< 3; i++){
    velo2 += (double) current.velocity[i]*current.velocity[i];
  }
  Kin = (double) 1/2*current.mass*velo2;
}

void solving::Delta_A(object &current,double &t, double dt[4],double &dA1, double &dA2,double h){
  if(t>dt[0]&&t<dt[1]){
    double velo2=0;
    double pos1=0;
    for (int i = 0; i< 3; i++){
      velo2 += (double) current.velocity[i]*current.velocity[i];
      pos1 += (double) current.position[i]*current.position[i];
    }
    double pos2 = sqrt(pos1);
    dA1 += (double) 1/2*velo2*pos2*h;
  }
  if(t>dt[2]&&t<dt[3]){
    double velo2=0;
    double pos1=0;
    for (int i = 0; i< 3; i++){
      velo2 += (double) current.velocity[i]*current.velocity[i];
      pos1 += (double) current.position[i]*current.position[i];
    }
    double pos2 = sqrt(pos1);
    dA2 += (double) 1/2*velo2*pos2*h;
  }
}

#include "object.hpp"
#include <cmath>
#include "armadillo"
#include <math.h>

object::object(double M,double x,double y,double z,double vx, double vy,double vz){
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
}

double object::distance(object otherObject){

    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = otherObject.position[0];
    y2 = otherObject.position[1];
    z2 = otherObject.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
}

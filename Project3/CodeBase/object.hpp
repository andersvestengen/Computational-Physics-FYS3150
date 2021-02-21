#ifndef OBJECT_HPP
#define OBJECT_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include "armadillo"

class object
{
    public:
    // Attributes
    double mass;
    double position[3];
    double velocity[3];



    //Initializer
    object(double M,double x,double y,double z,double vx, double vy,double vz);
    //Functions
    double distance(object otherObject);




};

#endif // PLANET_H

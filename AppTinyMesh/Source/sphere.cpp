#include "sphere.h"

Sphere::Sphere(double radius, Vector center) : _radius(radius), _center(center) {};

double Sphere::Radius() const { return _radius; }

Vector Sphere::Center() const { return _center; }

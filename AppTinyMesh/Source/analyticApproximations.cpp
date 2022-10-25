#include "analyticApproximations.h"

const double epsilon = 10e-4;

bool solveQuadratic(double& x1, double& x2, const double a, const double b, const double c)
{
    double delta = b*b - 4 * a *c;

    if(delta < 0)
        return false;
    else
    {
        if(delta > -epsilon && delta < epsilon)
            x1 = x2 = -b / (2 * a);
        else
        {
            double sqrtDelta = std::sqrt(delta);

            x1 = (-b + sqrtDelta) / (2 * a);
            x2 = (-b - sqrtDelta) / (2 * a);
        }

        return true;
    }
}

bool AnalyticSphere::intersect(const Ray& ray, double& t)
{
    double a = ray.Direction() * ray.Direction();
    double b = 2 * ray.Origin() * (ray.Direction() - this->_center);

    Vector oc = (ray.Origin() - this->_center);
    double c = oc*oc - this->_radius * this->_radius;

    double x1, x2;
    if(solveQuadratic(x1, x2, a, b, c))
    {
        if(x1 >= 0 && x2 >= 0)
        {
            t = std::min(x1, x2);

            return true;
        }
        else if(x1 >= 0)
        {
            t = x1;

            return true;
        }
        else if(x2 >= 0)
        {
            t = x2;

            return true;
        }
    }

    return false;
}

#include <cassert>

void AnalyticSphere::intersectionTest()
{
    AnalyticSphere sphere1(Vector(0, 0, 0), 1);

    double t;

    assert(sphere1.intersect(Ray(Vector(0, 0, -1), Vector(0, 0, 1)), t));
    assert(t == 0.0);

    assert(sphere1.intersect(Ray(Vector(0, 0, -1), Vector(0, 0, -2)), t));
    assert(t == 0.0);

    assert(!sphere1.intersect(Ray(Vector(0, 0, -2), Vector(0, 0, -3)), t));

    assert(sphere1.intersect(Ray(Vector(0, 0, -2), Normalized(Vector(1, 0, 0) - Vector(0, 0, -2))), t));
}

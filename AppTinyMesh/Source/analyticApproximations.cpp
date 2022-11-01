#include "analyticApproximations.h"
#include "qmath.h"

const double epsilon = 10e-6;

/*!
 * \brief Computes the closest intersection distance and stores the result in \p t
 * \param x1 First solution of the quadratic equation
 * \param x2 Second solution of the quadratic equation
 * \param t [out] This parameter will receive the closest
 * positive intersection distance if an intersection occured.
 * If no intersection was found, this parameter is unchanged
 * \return True if an intersection occured, false otherwise
 */
bool getTFromX1X2(const double x1, const double x2, double& t)
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

    return false;
}

bool solveQuadratic(double& x1, double& x2, const double a, const double b, const double c)
{
    double delta = b*b - 4 * a *c;

    if(delta < 0)
        return false;
    else
    {
        if(delta == 0)
        {
            double a2 = 2 * a;
            if(a2 == 0)
                return false;

            x1 = x2 = -b / a2;
        }
        else
        {
            double sqrtDelta = std::sqrt(delta);

            x1 = (-b + sqrtDelta) / (2 * a);
            x2 = (-b - sqrtDelta) / (2 * a);
        }

        return true;
    }
}

#include <iostream>

bool intersectWithDisk(const Ray& ray, const double diskRadius, const Vector& diskCenter, const Vector& diskNormal, double& t)
{
    double dotProd = diskNormal * ray.Direction();
    if(dotProd > epsilon || dotProd < -epsilon)//Different from 0
    {
        t = (diskCenter * diskNormal - ray.Origin() * diskNormal) / dotProd;

        if(t < 0)
            return false;

        Vector intersectionPoint = ray(t) - diskCenter;
        if(SquaredNorm(intersectionPoint) <= diskRadius * diskRadius)
            return true;
        else
            return false;
    }

    return false;
}

Vector AnalyticSphere::GetNormalAt(const Vector& vertex)
{
    return Normalized(vertex - _center);
}

bool AnalyticSphere::intersect(const Ray& ray, double& t)
{
    Vector oc = ray.Origin() - this->_center;
    double a = ray.Direction() * ray.Direction();
    double b = 2 * ray.Direction() * oc;
    double c = oc*oc - this->_radius * this->_radius;

    double x1, x2;
    if(!solveQuadratic(x1, x2, a, b, c))
        return false;

    return getTFromX1X2(x1, x2, t);
}

Vector AnalyticCylinder::GetNormalAt(const Vector& vertex)
{
    if(vertex[1] == this->_center[1])//The vertex is on the bottom disk of the cylinder
        return Vector(0, -1, 0);
    else if(vertex[1] == this->_center[1] + _height)//On the top disk of the cylinder
        return Vector(0, 1, 0);
    else//On the body of the cylinder
        return Normalized(vertex - (_center + Vector(0, vertex[1], 0)));
}

bool AnalyticCylinder::intersect(const Ray& ray, double& t)
{
    double centerX = this->_center[0];
    double centerY = this->_center[1];
    double centerZ = this->_center[2];

    double rayOrigX = ray.Origin()[0];
    double rayOrigZ = ray.Origin()[2];

    double rayDirX = ray.Direction()[0];
    double rayDirZ = ray.Direction()[2];

    double a = rayDirX * rayDirX + rayDirZ * rayDirZ;
    double b = 2 * (rayDirX * (rayOrigX - centerX) + rayDirZ * (rayOrigZ - centerZ));
    double c = rayOrigX * rayOrigX + rayOrigZ * rayOrigZ - 2 * (rayOrigX * centerX + rayOrigZ + centerZ) + centerX * centerX + centerZ * centerZ - this->_radius * this->_radius;

    double closestT = std::numeric_limits<double>::max();//This variable will keep the closest intersection we've found (bewteen the body of the cylinder and its disks)

    bool intersectionFound = false;

    //Intersection with the infinite (in height)
    //'body' of the cylinder
    double x1 = -1, x2 = -1;
    double tBody = std::numeric_limits<double>::max();
    if(solveQuadratic(x1, x2, a, b, c))
    {
        if(getTFromX1X2(x1, x2, tBody))
        {
            Vector intersectionPoint = ray(tBody);
            //If we have found an intersection with the
            //infinite cylinder that is below or above
            //the actual cylinder, there won't be an
            //intersection with the actual cylinder,
            //aborting
            if(intersectionPoint[1] < centerY || intersectionPoint[1] > this->_height - centerY)
                return false;

            intersectionFound = true;
            closestT = tBody;
        }
    }

    double tDisk = std::numeric_limits<double>::max();
    if(intersectWithDisk(ray, this->_radius, Vector(0, 0, 0) + this->_center, Vector(0, 1, 0), tDisk))
    {
        closestT = std::min(tDisk, tBody);

        intersectionFound = true;
    }

    if(intersectWithDisk(ray, this->_radius, Vector(0, this->_height, 0) + this->_center, Vector(0, -1, 0), tDisk))
    {
        closestT = std::min(closestT, tDisk);

        intersectionFound = true;
    }

    if(intersectionFound)
        t = closestT;

    return intersectionFound;
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

    assert(sphere1.intersect(Ray(Vector(1, 0, 0), Normalized(Vector(2, 0, 0) - Vector(1, 0, 0))), t));
    assert(t == 0.0);

    Vector direction = Normalized(Vector(2, 0, 0) - Vector(1, 0, 0));
    assert(!sphere1.intersect(Ray(Vector(1, 0, 0) + direction * 10e-4, direction), t));

    AnalyticSphere sphereTranslated(Vector(1, 0, 0), 1);
    assert(!sphereTranslated.intersect(Ray(Vector(3, 0, 0), Vector(1, 0, 0)), t));
    assert(!sphereTranslated.intersect(Ray(Vector(3, 0, 0), Vector(0, 1, 0)), t));
    assert(!sphereTranslated.intersect(Ray(Vector(3, 0, 0), Vector(0, 0, 1)), t));

    assert(sphereTranslated.intersect(Ray(Vector(3, 0, 0), Vector(-1, 0, 0)), t));
    assert(t == 1.0);

    assert(!sphere1.intersect(Ray(Vector(1, 0, 0) + Vector(1, 0, 0) * 1.0e-6, Vector(1, 0, 0)), t));

    //Rays inside the sphere
    assert(sphere1.intersect(Ray(Vector(0, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 1.0);

    assert(sphere1.intersect(Ray(Vector(0, 0, 0), Vector(-1, 0, 0)), t));
    assert(t == 1.0);

    assert(sphere1.intersect(Ray(Vector(0, 0, 0), Vector(0, 1, 0)), t));
    assert(t == 1.0);

    assert(sphere1.intersect(Ray(Vector(0, 0, 0), Vector(0, -1, 0)), t));
    assert(t == 1.0);

    assert(sphere1.intersect(Ray(Vector(0.5, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 0.5);

    assert(sphere1.intersect(Ray(Vector(0.5, 0, 0), Vector(-1, 0, 0)), t));
    assert(t == 1.5);
}

void AnalyticCylinder::intersectionTest()
{
    AnalyticCylinder cylinder1(Vector(0, 0, 0), 1, 2);

    double t;

    //Ray below the cylinder looking straigth up
    assert(cylinder1.intersect(Ray(Vector(0, -1, 0), Vector(0, 1, 0)), t));
    assert(t == 1.0);

    //Ray above the cylinder looking straigth down
    assert(cylinder1.intersect(Ray(Vector(0, 3, 0), Vector(0, -1, 0)), t));
    assert(t == 1.0);

    //Ray next to the cylinder looking straigth up
    assert(!cylinder1.intersect(Ray(Vector(2, -1, 0), Vector(0, 1, 0)), t));
    //Ray nexy the cylinder looking straigth down
    assert(!cylinder1.intersect(Ray(Vector(2, 3, 0), Vector(0, -1, 0)), t));

    //Ray to the right of the bottom disk of the cylinder looking straigth at the cylinder
    assert(cylinder1.intersect(Ray(Vector(2, 0, 0), Vector(-1, 0, 0)), t));
    assert(t == 1.0);

    //Ray to the left of the bottom disk cylinder looking straigth at the cylinder
    assert(cylinder1.intersect(Ray(Vector(-2, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 1.0);

    //Ray to the right of the cylinder, should intersect with the body of the cylinder
    assert(cylinder1.intersect(Ray(Vector(2, 1, 0), Vector(-1, 0, 0)), t));
    assert(t == 1.0);

    //Ray to the left of the cylinder, should intersect with the body of the cylinder
    assert(cylinder1.intersect(Ray(Vector(-2, 1, 0), Vector(1, 0, 0)), t));
    assert(t == 1.0);

    //Ray next to the cylinder but pointing outwards, we can't have a 'negative distance' intersection
    assert(!cylinder1.intersect(Ray(Vector(2, 0, 0), Vector(1, 0, 0)), t));

    AnalyticCylinder cylinderModified(Vector(2, 0, 0), 2, 1);

    //Ray to the left pointing at the cylinder
    assert(cylinderModified.intersect(Ray(Vector(-1, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 1.0);

    //Ray further to the left pointing at the cylinder
    assert(cylinderModified.intersect(Ray(Vector(-8, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 8.0);

    //Ray to the right pointing at the cylinder
    assert(cylinderModified.intersect(Ray(Vector(8, 0, 0), Vector(-1, 0, 0)), t));
    assert(t == 4.0);

    //Ray to the right and above the cylinder, no intersection
    assert(!cylinderModified.intersect(Ray(Vector(8, 2, 0), Vector(-1, 0, 0)), t));
    //Ray to the left nd below the cylinder, no intersection
    assert(!cylinderModified.intersect(Ray(Vector(-1, -2, 0), Vector(1, 0, 0)), t));

    //Ray above pointing down
    assert(cylinderModified.intersect(Ray(Vector(2, 2, 0), Vector(0, -1, 0)), t));
    assert(t == 1.0);

    //Ray below pointing up
    assert(cylinderModified.intersect(Ray(Vector(2, -3, 0), Vector(0, 1, 0)), t));
    assert(t == 3.0);

    //Ray inside the cylinder pointing to the right
    assert(cylinderModified.intersect(Ray(Vector(3, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 1.0);

    //Ray inside the cylinder pointing to the right
    assert(cylinderModified.intersect(Ray(Vector(3.5, 0, 0), Vector(1, 0, 0)), t));
    assert(t == 0.5);

    //Ray inside the cylinder pointing up
    assert(cylinderModified.intersect(Ray(Vector(3, 0.5, 0), Vector(0, 1, 0)), t));
    assert(t == 0.5);

    AnalyticCylinder cylinderTranslatedUp(Vector(0, 2, 0), 1, 2);

    //Ray below pointing up
    assert(cylinderTranslatedUp.intersect(Ray(Vector(0, 0, 0), Vector(0, 1, 0)), t));
    assert(t == 2.0);

    //Ray above pointing down
    assert(cylinderTranslatedUp.intersect(Ray(Vector(0, 5, 0), Vector(0, -1, 0)), t));
    assert(t == 1.0);

    AnalyticCylinder cylinderRad1Height2(Vector(0, 0, 0), 1, 2);
    //Ray's origin at the border of the cylinder, between the top disk and the body. Auto intersection
    assert(cylinderRad1Height2.intersect(Ray(Vector(1, 2, 0), Vector(0.295226, 0.181165, -0.938094)), t));
    assert(t == 0.0);

    //Generating 100 random rays direction and ray's origin around the body of the cylinder. They should all auto-intersect.
    //Also testing for the same ray / ray direction but with the ray origin offset off the cylinder which should prevent
    //auto intersection
    //TODO Not working
//    unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();
//    Vector dump = Normalized(Vector((Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                    (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                    (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1));
//            Normalized(Vector((Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                              (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                              (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1));
//            Normalized(Vector((Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                              (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                              (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1));
//    for(int i = 0; i < 100; i++)
//    {
//        double randomAngle = (Math::xorshift96() / (double)max_unsigned_int) * 2 * M_PI;
//        double randomHeight = (Math::xorshift96() / (double)max_unsigned_int) * 2;

//        double x = std::cos(randomAngle);
//        double y = randomHeight;
//        double z = std::sin(randomAngle);

//        Vector rayOrigin(x, y, z);
//        Vector normalAtRayOrigin = rayOrigin - Vector(0, randomHeight, 0);
//        Vector randomRayDirection = Normalized(Vector((Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                      (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
//                                                      (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1));
//        if(randomRayDirection * normalAtRayOrigin < 0)
//            randomRayDirection *= -1;

//        std::cout << normalAtRayOrigin << std::endl;
//        std::cout << rayOrigin << std::endl;
//        std::cout << randomRayDirection << std::endl;

//        //Auto intersection
//        assert(cylinderRad1Height2.intersect(Ray(rayOrigin, randomRayDirection), t));
//        assert(t == 0.0);

//        //Origin + Offset : shouldn't auto intersect
//        assert(!cylinderRad1Height2.intersect(Ray(rayOrigin + 1.0e-5 * normalAtRayOrigin, randomRayDirection), t));
//    }
}

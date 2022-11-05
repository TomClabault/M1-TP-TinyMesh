#include "analyticApproximations.h"
#include "qmath.h"

#include <vector>

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
    if(dotProd > Math::EPSILON || dotProd < Math::EPSILON)//Different from 0
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
    //First transforming the given vertex to the coordinate system of the unit sphere at the origin
    Vector transformedVertex = vertex - _translation;
    transformedVertex = transformedVertex * _invRotation;
    transformedVertex = transformedVertex * _invScale;

    //For a sphere, the normal is the point itself so there's nothing more
    //to do than transform the normal back to the transformed object space
    Vector normal = transformedVertex * _rotation.GetInverse().GetTranspose() * _scale.GetInverse().GetTranspose();

    return Normalized(normal);
}

bool AnalyticSphere::intersect(const Ray& ray, double& t)
{
    Vector invTransformedOrigin = ray.Origin();
    invTransformedOrigin -= _translation;
    invTransformedOrigin = invTransformedOrigin * _invRotation;
    invTransformedOrigin = invTransformedOrigin * _invScale;

    Vector invTransformedDirection = ray.Direction();
    invTransformedDirection = invTransformedDirection * _invRotation;
    invTransformedDirection = invTransformedDirection * _invScale;

    Ray transformedRay(invTransformedOrigin, Normalized(invTransformedDirection));

    double transformedT;
    if(intersectBasic(transformedRay, transformedT))
    {
        Vector intersectionPoint = transformedRay(transformedT);

        Vector pointWorldCoords = intersectionPoint * _rotation * _scale;
        pointWorldCoords += _translation;

        //Computing the distance to the ray's origin in world coords
        ray.T(pointWorldCoords, t);

        return true;
    }
    else
        return false;
}

bool AnalyticSphere::intersectBasic(const Ray& ray, double& t)
{
    double a = 1;//Ray direction is assumed normalized
    double b = 2 * ray.Direction() * ray.Origin();
    double c = ray.Origin() * ray.Origin() - 1;

    double x1, x2;
    if(!solveQuadratic(x1, x2, a, b, c))
        return false;

    bool tFound = getTFromX1X2(x1, x2, t);

    //Considering there is no intersection below a certain threshold epsilon
    //Intersection distance such as 1.0e-7 are within
    return tFound && (t == 0.0 || t > Math::EPSILON);
}

Vector AnalyticCylinder::GetNormalAt(const Vector& vertex)
{
    //Getting the vertex as if on a unit cylinder
    Vector transformedVertex = vertex - _translation;
    transformedVertex = transformedVertex * _invRotation;
    transformedVertex = transformedVertex * _invScale;

    Vector transformedNormal;

    //Computing the un-transformed normal
    if(transformedVertex[1] < Math::BIGGER_EPSILON)//The transformed vertex is on the bottom disk of the cylinder
        transformedNormal = Vector(0, -1, 0);
    else if(transformedVertex[1] > 1 - Math::BIGGER_EPSILON && transformedVertex[1] < 1 + Math::BIGGER_EPSILON)//On the top disk of the cylinder
        transformedNormal = Vector(0, 1, 0);
    else//On the body of the cylinder
        transformedNormal = Vector(transformedVertex[0], 0, transformedVertex[2]);

    Vector normal = transformedNormal * _rotation.GetInverse().GetTranspose() * _scale.GetInverse().Transpose();

    return Normalized(normal);
}

bool AnalyticCylinder::intersect(const Ray& ray, double& t)
{
    Vector invTransformedOrigin = ray.Origin();
    invTransformedOrigin -= _translation;
    invTransformedOrigin = invTransformedOrigin * _invRotation;
    invTransformedOrigin = invTransformedOrigin * _invScale;

    Vector invTransformedDirection = ray.Direction();
    invTransformedDirection = invTransformedDirection * _invRotation;
    invTransformedDirection = invTransformedDirection * _invScale;

    Ray transformedRay(invTransformedOrigin, Normalized(invTransformedDirection));

    double transformedT;
    if(intersectBasic(transformedRay, transformedT))
    {
        Vector intersectionPoint = transformedRay(transformedT);

        Vector pointWorldCoords = intersectionPoint * _rotation * _scale;
        pointWorldCoords += _translation;

        //Computing the distance to the ray's origin in world coords
        ray.T(pointWorldCoords, t);

        return true;
    }
    else
        return false;
}

bool AnalyticCylinder::intersectBasic(const Ray& ray, double& t)
{
    double rayOrigX = ray.Origin()[0];
    double rayOrigZ = ray.Origin()[2];

    double rayDirX = ray.Direction()[0];
    double rayDirZ = ray.Direction()[2];

    double a = rayDirX * rayDirX + rayDirZ * rayDirZ;
    double b = 2 * (rayDirX * rayOrigX + rayDirZ * rayOrigZ);
    double c = rayOrigX * rayOrigX + rayOrigZ * rayOrigZ - 1;

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
            if(intersectionPoint[1] < 0 || intersectionPoint[1] > 1)
                return false;

            intersectionFound = true;
            closestT = tBody;
        }
    }

    double tDisk = std::numeric_limits<double>::max();
    if(intersectWithDisk(ray, 1, Vector(0, 0, 0), Vector(0, 1, 0), tDisk))
    {
        closestT = std::min(tDisk, tBody);

        intersectionFound = true;
    }

    if(intersectWithDisk(ray, 1, Vector(0, 1, 0), Vector(0, -1, 0), tDisk))
    {
        closestT = std::min(closestT, tDisk);

        intersectionFound = true;
    }

    if(intersectionFound)
        t = closestT;

    return intersectionFound;
}

#include <cassert>
void AnalyticSphere::intersectionTests()
{
    AnalyticSphere sphere1(Vector(0, 0, 0), 1);

    double t;

    assert(sphere1.intersect(Ray(Vector(0, 0, -1), Vector(0, 0, 1)), t));
    assert(t == 0.0);

    assert(sphere1.intersect(Ray(Vector(0, 0, -1), Vector(0, 0, -1)), t));
    assert(t == 0.0);

    assert(!sphere1.intersect(Ray(Vector(0, 0, -2), Vector(0, 0, -1)), t));

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

    AnalyticSphere sphereTranslated2(Vector(2, 0, 0), 1);
    assert(!sphereTranslated2.intersect(Ray(Vector(4, 0, 0), Vector(1, 0, 0)), t));
    assert(!sphereTranslated2.intersect(Ray(Vector(4, 0, 0), Vector(0, 1, 0)), t));
    assert(!sphereTranslated2.intersect(Ray(Vector(4, 0, 0), Vector(0, 0, 1)), t));

    assert(sphereTranslated2.intersect(Ray(Vector(4, 0, 0), Vector(-1, 0, 0)), t));
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

    //Special cases that have caused problems
    AnalyticSphere sphereSpecialCase1(Vector(0, 0, 0), 1);
    assert(!sphereSpecialCase1.intersect(Ray(Vector(0.8506507873535156,-0.525731086730957,0.0001), Vector(0.2680317789104351,0.3653358064721789,0.8914531473966706)), t));

    AnalyticSphere sphereSpecialCase2(Vector(1, 0, std::sqrt(3) * 1), 1);
    t = -INFINITY;
    assert(sphereSpecialCase2.intersect(Ray(Vector(0.8506507873535156,-0.525731086730957,0.0001), Vector(0.2680317789104351,0.3653358064721789,0.8914531473966706)), t));
    assert(t != -INFINITY);

    AnalyticSphere sphereSpecialCase3(Vector(8, 0, 0), 4);
    t = -INFINITY;
    assert(sphereSpecialCase3.intersect(Ray(Vector(3.40261,-2.10293,0), Vector(0.719253,0.690758,0.0743505)), t));
    assert(t != -INFINITY);
//    ray origin: Vector(3.40261,-2.10293,0)
//    ray direction: Vector(0.719253,0.690758,0.0743505)
}

void AnalyticSphere::normalAtTests()
{
    double scale = 2;
    double scaleBig = 100;
    Vector translate(2, 0, 0);

    AnalyticSphere unitSphere(1);
    AnalyticSphere scaledSphere(scale);
    AnalyticSphere scaledBigSphere(scaleBig);
    AnalyticSphere translatedUnitSphere(translate, 1);
    AnalyticSphere translatedScaledSphere(translate, scale);
    AnalyticSphere translatedScaledBigSphere(translate, scaleBig);

    //Testing for a lot of random points.
    //On a unit sphere (translated and/or uniformly scaled or not), the normal should always be the point itself
    for(int i = 0; i < 100000; i++)
    {
        Vector randomPoint = Normalized(Vector(Math::xorshift96(), Math::xorshift96(), Math::xorshift96()));

        //We're using the distance to compare the results as directly comparing floating point numbers is unreliable
        assert(Distance(unitSphere.GetNormalAt(randomPoint), randomPoint) <= Math::EPSILON);
        assert(Distance(scaledSphere.GetNormalAt(randomPoint * scale), randomPoint) <= Math::EPSILON);
        assert(Distance(scaledBigSphere.GetNormalAt(randomPoint * scaleBig), randomPoint) <= Math::EPSILON);
        assert(Distance(translatedUnitSphere.GetNormalAt(randomPoint + translate), randomPoint) <= Math::EPSILON);
        assert(Distance(translatedScaledSphere.GetNormalAt(randomPoint * scale + translate), randomPoint) <= Math::EPSILON);
        assert(Distance(translatedScaledBigSphere.GetNormalAt(randomPoint * scaleBig + translate), randomPoint) <= Math::EPSILON);
    }
}

void AnalyticSphere::tests()
{
    AnalyticSphere::intersectionTests();
    AnalyticSphere::normalAtTests();
}

void AnalyticCylinder::intersectionTests()
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

void AnalyticCylinder::getNormalAtTests()
{
    //For any point on the body of a translated/uniformly scaled cylinder, the normal
    //at this point should be the point itself (normalized)

    double scale = 2;
    double scaleBig = 100;
    Vector translate(2, 0, 0);

    AnalyticCylinder unitCylinder(Vector(0, 0, 0), 1, 1);
    AnalyticCylinder scaledCylinder(Vector(0, 0, 0), scale, scale);
    AnalyticCylinder scaledBigCylinder(Vector(0, 0, 0), scaleBig, scaleBig);
    AnalyticCylinder translatedUnitCylinder(translate, 1, 1);
    AnalyticCylinder translatedScaledCylinder(translate, scale, scale);
    AnalyticCylinder translatedScaledBigCylinder(translate, scaleBig, scaleBig);

    //Testing for a lot of random points.
    for(int i = 0; i < 100000; i++)
    {
        unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

        //Only the y coordinate is generated between 0.1 and the scale -0.1 to be sure the point is on the body of the cylinder
        Vector scaleRandomPoint = Normalized(Vector(Math::xorshift96()/ (double)max_unsigned_int, (Math::xorshift96() / (double)max_unsigned_int) * (scale - 0.2) + 0.1, Math::xorshift96()/ (double)max_unsigned_int));
        Vector scaleBigRandomPoint = Normalized(Vector(Math::xorshift96()/ (double)max_unsigned_int, (Math::xorshift96() / (double)max_unsigned_int) * (scaleBig - 0.2) + 0.1, Math::xorshift96()/ (double)max_unsigned_int));
        Vector randomPoint = Normalized(Vector(Math::xorshift96()/ (double)max_unsigned_int, Math::xorshift96()/ (double)max_unsigned_int, Math::xorshift96()/ (double)max_unsigned_int));

        Vector randomPointTopDisk = Vector(Math::xorshift96()/ (double)max_unsigned_int, 1, Math::xorshift96()/ (double)max_unsigned_int);
        Vector randomPointBottomDisk = Vector(Math::xorshift96()/ (double)max_unsigned_int, 0, Math::xorshift96()/ (double)max_unsigned_int);

        Vector expectedTopDisk = Vector(0, 1, 0);
        Vector expectedBottomDisk = Vector(0, -1, 0);

        //The expected normals will have a y coordinate equal to 0 if the point is on the body
        //of the cylinder and if the cylidner is uniformly scaled (or not scaled at all)
        Vector expectedScaleRandom = scaleRandomPoint;
        Vector expectedScaleBigRandom = scaleBigRandomPoint;
        Vector expectedRandom = randomPoint;
        expectedScaleRandom[1] = 0;
        expectedScaleBigRandom[1] = 0;
        expectedRandom[1] = 0;
        expectedScaleRandom = Normalized(expectedScaleRandom);
        expectedScaleBigRandom = Normalized(expectedScaleBigRandom);
        expectedRandom = Normalized(expectedRandom);

        //We're using the distance to compare the results as directly comparing floating point numbers is unreliable
        //So if the distance from the given normal and the expected one is low enough, it's the same point
        assert(Distance(unitCylinder.GetNormalAt(randomPoint), expectedRandom) <= Math::EPSILON);
        assert(Distance(scaledCylinder.GetNormalAt(scaleRandomPoint * scale), expectedScaleRandom) <= Math::EPSILON);
        assert(Distance(scaledBigCylinder.GetNormalAt(scaleRandomPoint * scaleBig), expectedScaleRandom) <= Math::EPSILON);
        assert(Distance(translatedUnitCylinder.GetNormalAt(randomPoint + translate), expectedRandom) <= Math::EPSILON);
        assert(Distance(translatedScaledCylinder.GetNormalAt(scaleRandomPoint * scale + translate), expectedScaleRandom) <= Math::EPSILON);
        assert(Distance(translatedScaledBigCylinder.GetNormalAt(scaleRandomPoint * scaleBig + translate), expectedScaleRandom) <= Math::EPSILON);

        //Testing for points of the top and bottom disk
        assert(Distance(unitCylinder.GetNormalAt(randomPointTopDisk), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(unitCylinder.GetNormalAt(randomPointBottomDisk), expectedBottomDisk) <= Math::EPSILON);
        assert(Distance(scaledCylinder.GetNormalAt(randomPointTopDisk * scale), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(scaledCylinder.GetNormalAt(randomPointBottomDisk * scale), expectedBottomDisk) <= Math::EPSILON);
        assert(Distance(scaledBigCylinder.GetNormalAt(randomPointTopDisk * scaleBig), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(scaledBigCylinder.GetNormalAt(randomPointBottomDisk * scaleBig), expectedBottomDisk) <= Math::EPSILON);
        assert(Distance(translatedUnitCylinder.GetNormalAt(randomPointTopDisk + translate), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(translatedUnitCylinder.GetNormalAt(randomPointBottomDisk + translate), expectedBottomDisk) <= Math::EPSILON);
        assert(Distance(translatedScaledCylinder.GetNormalAt(randomPointTopDisk * scale + translate), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(translatedScaledCylinder.GetNormalAt(randomPointBottomDisk * scale + translate), expectedBottomDisk) <= Math::EPSILON);
        assert(Distance(translatedScaledBigCylinder.GetNormalAt(randomPointTopDisk * scaleBig + translate), expectedTopDisk) <= Math::EPSILON);
        assert(Distance(translatedScaledBigCylinder.GetNormalAt(randomPointBottomDisk * scaleBig + translate), expectedBottomDisk) <= Math::EPSILON);
    }
}

void AnalyticCylinder::tests()
{
    AnalyticCylinder::intersectionTests();
    AnalyticCylinder::getNormalAtTests();
}

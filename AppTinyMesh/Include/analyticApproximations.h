#pragma once

#include "ray.h"

class AnalyticApproximation
{
public:
    virtual bool intersect(const Ray& ray, double& t) = 0;
};

class AnalyticSphere : public AnalyticApproximation
{
public:
    AnalyticSphere(Vector center, double radius) : _center(center), _radius(radius) {};

    Vector Center() { return _center; }

    virtual bool intersect(const Ray& ray, double& t) override;

    static void intersectionTest();

private:
    Vector _center;

    double _radius;
};

class AnalyticCylinder : public AnalyticApproximation
{
public:
    AnalyticCylinder(Vector center, double radius, double height) : _center(center), _radius(radius), _height(height) {};

    virtual bool intersect(const Ray& ray, double& t) override;

    static void intersectionTest();

private:
    Vector _center;

    double _radius;
    double _height;
};

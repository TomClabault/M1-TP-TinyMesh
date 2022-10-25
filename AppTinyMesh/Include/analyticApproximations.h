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

    virtual bool intersect(const Ray& ray, double& t) override;

    static void intersectionTest();

private:
    Vector _center;

    double _radius;
};

#include "ray.h"

template <typename ShapeType>
class AnalyticShape
{
    bool intersect(const Ray& ray, double& t, double& u, double &v);
};

class AnalyticSphere
{
    bool intersect(const Ray& ray, double& t, double& u, double &v)
    {
        //TODO
        return true;
    }
};

#include "mathematics.h"

class Sphere
{
public:
    /*!
     * \brief Creates a sphere of radius \p radius that has for center \p center
     *
     * \param radius The radius of the sphere
     * \param center The center of the sphere
     */
    Sphere(double radius, Vector center);

    double Radius() const;

    Vector Center() const;

private:
    double _radius;

    Vector _center;
};

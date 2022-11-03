#ifndef __Ray__
#define __Ray__

#include "mathematics.h"

// Ray class
class Ray
{
protected:
  Vector origin; //!< Origin of the ray.
  Vector direction; //!< Direction.
public:
  //! Empty.
  Ray() {}
  explicit Ray(const Vector&, const Vector&);
  //! Empty.
  ~Ray() {}

  Ray Reflect(const Vector&, const Vector&);

  /*!
   * \brief For a ray equation: O + Dt, returns t given a point on the ray
   * \param point The point on the ray
   * \param t [out] The distance from the ray's origin
   * \return True if the point is on the ray, false otherwise
   */
  bool T(const Vector& point, double& t) const;

  Vector operator()(double) const;

  // Functions to access Vector class components
  Vector Origin() const;
  Vector Direction() const;

  friend std::ostream& operator<<(std::ostream&, const Ray&);
};

/*!
\brief Creates a ray.

The direction should be unit:
\code
Ray ray(Vector(0.0,0.0,0.0),Normalized(Vector(2.0,-1.0,3.0)));
\endcode
\param p Origin.
\param d Direction (should be unit vector).
*/
inline Ray::Ray(const Vector& p, const Vector& d)
{
  origin = p;
  direction = d;
}

/*!
\brief Return the origin of the ray.
*/
inline Vector Ray::Origin() const
{
  return origin;
}

/*!
\brief Return the direction of the ray.
*/
inline Vector Ray::Direction() const
{
  return direction;
}

inline bool Ray::T(const Vector& point, double& t) const
{
    Vector cross = direction / (point - origin);
    if(cross[0] > Math::EPSILON || cross[0] < -Math::EPSILON
    || cross[1] > Math::EPSILON || cross[1] < -Math::EPSILON
    || cross[2] > Math::EPSILON || cross[2] < -Math::EPSILON)//Point not on the ray
        return false;
    else
    {
        if(point == origin)
            t = 0.0;
        else
            t = Norm(point - origin);

        return true;
    }
}

/*!
\brief Computes the location of a vertex along the ray.
\param t Parameter.
*/
inline Vector Ray::operator()(double t) const
{
  return origin + t * direction;
}

#endif

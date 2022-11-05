#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "box.h"
#include "mathematics.h"
#include "ray.h"

class Triangle
{
protected:
  Vector p[3] = {Vector(0.0,0.0,0.0),Vector(1.0,0.0,0.0), Vector(0.0,1.0,0.0), }; //!< Array of vertices.
public:
  //! Empty.
  Triangle() {}
  explicit Triangle(const Vector&, const Vector&, const Vector&);

  Triangle(const Triangle& triangle) : p {triangle.p[0], triangle.p[1], triangle.p[2]} {};

  //! Empty.
  ~Triangle() {}

  Vector operator[] (int) const;

  // Point in triangle
  Vector Vertex(double, double) const;

  // Intersection
  bool Intersect(const Ray&, double&, double&, double&) const;

  /*!
  * \brief This function allows the computation of the intersection between a ray and a triangle without having to instantiate a triangle
  */
  static bool IntersectFromPoints(const Vector& a, const Vector& b, const Vector& c, const Ray&, double&, double&, double&);

  void Translate(const Vector&);

  // Geometry
  Vector Normal() const;
  Vector AreaNormal() const;
  Vector Center() const;

  double Area() const;
  double Aspect() const;
  Box GetBox() const;

  Vector Centroid() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Triangle&);

  double InscribedRadius() const;
  double CircumscribedRadius() const;
protected:
  static double epsilon; //!< Internal epsilon constant.
};

/*!
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Triangle::operator[] (int i) const
{
  return p[i];
}

//! Compute the barycenter of the triangle.
inline Vector Triangle::Center() const
{
  return (p[0] + p[1] + p[2]) / 3.0;
}

//! Compute the area of the triangle.
inline double Triangle::Area() const
{
  return 0.5 * Norm((p[0] - p[1]) / (p[2] - p[0]));
}

/*!
\brief Create a triangle.
\param a,b,c Vertices of the triangle.
*/
inline Triangle::Triangle(const Vector& a, const Vector& b, const Vector& c)
{
  p[0] = a;
  p[1] = b;
  p[2] = c;
}

#endif // TRIANGLE_H

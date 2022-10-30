#pragma once

#include "box.h"
#include "color.h"
#include "simpleMeshes.h"
#include "ray.h"
#include "mathematics.h"
#include "matrix.h"
#include "sphere.h"

// Triangle
class Triangle
{
protected:
  Vector p[3] = {Vector(0.0,0.0,0.0),Vector(1.0,0.0,0.0), Vector(0.0,1.0,0.0), }; //!< Array of vertices.
public:
  //! Empty.
  Triangle() {}
  explicit Triangle(const Vector&, const Vector&, const Vector&);

  //TODO remove. On doit faire en sorte que dans le code, on ait plus une seule copie de triangle puisque normalement, une fois que les triangles sont construits, on se balade des références des triangles, on ne les copie pas
  Triangle(Triangle&) {std::cout << "copy constructor triangel";};

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
  static bool IntersectFromPoints(Vector a, Vector b, Vector c, const Ray&, double&, double&, double&);

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


class QString;

class Mesh
{
protected:
  std::vector<Vector> vertices; //!< Vertices.
  std::vector<Vector> normals;  //!< Normals.
  std::vector<int> varray;     //!< Vertex indexes.
  std::vector<int> narray;     //!< Normal indexes.

  std::vector<AnalyticApproximation*> analyticApproximations; //!< The shape that can approximate the current mesh. If no such shape exists,

public:
  explicit Mesh();
  explicit Mesh(const std::vector<Vector>&, const std::vector<int>&);
  explicit Mesh(const std::vector<Vector>&, const std::vector<Vector>&, const std::vector<int>&, const std::vector<int>&);
  ~Mesh();

  void Reserve(int, int, int, int);

  Triangle GetTriangle(int) const;
  std::vector<Triangle*>* GetTriangles() const;
  Vector Vertex(int) const;
  Vector Vertex(int, int) const;

  Vector Normal(int) const;

  int Triangles() const;
  int Vertexes() const;

  std::vector<Vector> Vertices() const;

  std::vector<int> VertexIndexes() const;
  std::vector<int> NormalIndexes() const;

  int VertexIndex(int, int) const;
  int NormalIndex(int, int) const;

  Vector operator[](int) const;

  Box GetBox() const;

  /*!
   * \brief Scales the mesh along the X, Y and Z axis by the given factor
   * \param factor The scaling factor
   */
  void Scale(double factor);

  /*!
   * \brief Scales the mesh using the given homothety matrix
   * \param homothetyMatrix The homothety matrix
   */
  void Scale(const Matrix& homothetyMatrix);

  /*!
   * \brief Rotates the mesh using the given rotation matrix
   * \param rotationMatrix The rotation matrix
   */
  void Rotate(const Matrix& rotationMatrix);

  /*!
   * \brief Translates the mesh using the given translation vector
   * \param translationVector The translation vector
   */
  void Translate(const Vector& translationVector);

  /*!
   * \brief Merges another mesh into the current one
   * \param secondMesh The mesh to merge into the current one
   */
  void Merge(const Mesh& secondMesh);

  /*!
   * \brief Deforms the mesh based on the proximity
   * of the vertices of the mesh to the given sphere.
   * The closer a vertex of the mesh to the center of
   * the sphere, the stronger the deformation.
   *
   * \param sphere The sphere that is going to be used
   * to deform the mesh
   */
  void SphereWarp(Sphere sphere);

  /*!
   * \brief Constructs a vector of color for each
   * vertex of the mesh based on the accessibility
   * of each vertices.
   *
   * \param radius The radius within which to check
   * for intersections with the mesh
   * \param samples How many samples are going to be
   * used to compute the accessibility of each vertex
   * \param accessibilityColors The vector of size
   * this->Vertexes().size() that will receive the color
   * of each vertex (in the same order as the vertices
   * of the this->Vertexes() vector) based on their
   * accessibility. This vector should already be resized
   * to this->Vertexes()
   * \param occlusionStrength A multiplier on the strength
   * of the occlusion. The higher this parameter, the
   * darker the occlusion
   */
  void accessibility(std::vector<Color>& accessibilityColors, double radius, int samples, double occlusionStrength = 1);

  /*!
   * \brief Computes the intersection between a ray
   * and this mesh
   *
   * \param ray The ray
   * \param [out] outT The distance t from the origin of the
   * ray where the intersection has been found asssuming
   * the ray's equation is: orig + dir*t
   *
   * \return True if an intersection was found, false otherwise
   */
  bool intersect(const Ray& ray, double& outT) const;

  /*!
   * \brief Computes the intersection of a ray and all the analytic approximations that compose this mesh
   *
   * \param ray The ray that is going to be tested against the mesh
   * \param t [out] The nearest intersection found with the analytic approximations of this mesh
   *
   * \return True if an intersection was found, false otherwise
   */
  bool intersectAnalytic(const Ray& ray, double& t) const;

  void SmoothNormals();

  // Constructors from core classes
  explicit Mesh(const Box&);
  explicit Mesh(const SimpleMesh&);

  void Load(const QString&);
  void SaveObj(const QString&, const QString&) const;
protected:
  void AddTriangle(int, int, int, int);
  void AddSmoothTriangle(int, int, int, int, int, int);
  void AddSmoothQuadrangle(int, int, int, int, int, int, int, int);
  void AddQuadrangle(int, int, int, int);
};

inline std::vector<Vector> Mesh::Vertices() const
{
    return this->vertices;
}

/*!
\brief Return the set of vertex indexes.
*/
inline std::vector<int> Mesh::VertexIndexes() const
{
  return varray;
}

/*!
\brief Return the set of normal indexes.
*/
inline std::vector<int> Mesh::NormalIndexes() const
{
  return narray;
}

/*!
\brief Get the vertex index of a given triangle.
\param t Triangle index.
\param i Vertex index.
*/
inline int Mesh::VertexIndex(int t, int i) const
{
  return varray.at(t * 3 + i);
}

/*!
\brief Get the normal index of a given triangle.
\param t Triangle index.
\param i Normal index.
*/
inline int Mesh::NormalIndex(int t, int i) const
{
  return narray.at(t * 3 + i);
}

/*!
\brief Get a triangle.
\param i Index.
\return The triangle.
*/
inline Triangle Mesh::GetTriangle(int i) const
{
  return Triangle(vertices.at(varray.at(i * 3 + 0)), vertices.at(varray.at(i * 3 + 1)), vertices.at(varray.at(i * 3 + 2)));
}

inline std::vector<Triangle*>* Mesh::GetTriangles() const
{
    std::vector<Triangle*>* triangles = new std::vector<Triangle*>;
    triangles->reserve(this->VertexIndexes().size() / 3);

    for(int vertexIndex = 0; vertexIndex < VertexIndexes().size(); vertexIndex += 3)
    {
        triangles->push_back(new Triangle(this->vertices.at(this->varray.at(vertexIndex + 0)),
                                      this->vertices.at(this->varray.at(vertexIndex + 1)),
                                      this->vertices.at(this->varray.at(vertexIndex + 2))));
    }

    return triangles;
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector Mesh::Vertex(int i) const
{
  return vertices[i];
}

/*!
\brief Get a vertex from a specific triangle.
\param t The number of the triangle wich contain the wanted vertex.
\param v The triangle vertex: 0, 1, or 2.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector Mesh::Vertex(int t, int v) const
{
  return vertices[varray[t * 3 + v]];
}

/*!
\brief Get the number of vertices in the geometry.
\return The number of vertices in the geometry, in other words the size of vertices.
*/
inline int Mesh::Vertexes() const
{
  return int(vertices.size());
}

/*!
\brief Get a normal.
\param i Index of the wanted normal.
\return The normal.
*/
inline Vector Mesh::Normal(int i) const
{
  return normals[i];
}

/*!
\brief Get the number of triangles.
*/
inline int Mesh::Triangles() const
{
  return int(varray.size()) / 3;
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
\see vertex(int i) const
*/
inline Vector Mesh::operator[](int i) const
{
  return vertices[i];
}


#pragma once

#include "analyticApproximations.h"
#include "box.h"
#include "BVH.h"
#include "color.h"
#include "simpleMeshes.h"
#include "ray.h"
#include "mathematics.h"
#include "matrix.h"
#include "sphere.h"
#include "triangle.h"

#include <map>

class QString;

class Mesh
{
protected:
  std::vector<Vector> vertices; //!< Vertices.
  std::vector<Vector> normals;  //!< Normals.
  std::vector<int> vertexIndices;     //!< Vertex indexes.
  std::vector<int> normalIndices;     //!< Normal indexes.

  std::vector<AnalyticApproximation*> analyticApproximations; //!< The shape that can approximate the current mesh. If no such shape exists,

  //!< The map that is used to store the analytic approximation
  //!< associated with a certain vertex index.
  //!< This map maps the largest vertex number of the mesh to
  //!< its analytic approximation (or nullptr) if the mesh doesn't
  //!< have an analytic approximation
  //!< As an example, if a mesh has 20 vertices, then this map will
  //!< contain an entry (20, approx) that maps 20 (the number of vertices)
  //!< to the approximation of the mesh. If another mesh of 20 vertices
  //!< is merged into this one, another entry (40, approx2) will be added
  //!< to the map. The key of this entry is 40 since the first 20 vertices
  //!< are taken up by the first mesh.
  std::map<int, AnalyticApproximation*> analyticApproxToVertexIndex;

  //!< The Bounding Volume Hierarchy used to optimize mesh intersections
  BVH _bvh;

  unsigned long long int _triangleInterTests = 0;
  unsigned long long int _triangleEffectiveInter = 0;

public:
  explicit Mesh();
  explicit Mesh(const std::vector<Vector>&, const std::vector<int>&);
  explicit Mesh(const std::vector<Vector>&, const std::vector<Vector>&, const std::vector<int>&, const std::vector<int>&);
  ~Mesh();

  void Reserve(int, int, int, int);

  /*!
   * \brief Returns the analytic approximation associated with a certain vertex index
   * \param vertexIndex The vertex index (this should be a number out of the varray vector)
   * \return The analytic approximation to which the vertex "belongs".
   * If there is no analytic approximation for this vertex, returns nullptr.
   */
  AnalyticApproximation* GetAnalyticApproxFromVertexIndex(int vertexIndex) const;

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
   * \param enableAnalyticIntersection Whether or not to allow
   * for the intersection with the analytic approximations
   * of the mesh or not. This parameter has no practical usecases
   * besides benchmarking
   * \param useBVH Whether or not to use the BVH to intersect
   * the mesh. Setting this to false has no practical usecases
   * besides benchmarking
   */
  void accessibility(std::vector<Color>& accessibilityColors, double radius, int samples, double occlusionStrength = 1, bool enableAnalyticIntersection = true, bool useBVH = true);

  /*!
   * \brief Gets the intersection statistics of the mesh.
   * The statistics for this mesh are initialized at 0 at the
   * construction of the mesh and are never reset to 0
   * \param triangleInterTests [out] The number of intersection
   * tests that have been performed between a ray and a triangle
   * \param triangleEffectiveInters The number of triangles that
   * have effectively been intersected by a ray
   */
  void GetInterStats(unsigned long long int& triangleInterTests, unsigned long long int& triangleEffectiveInters) const;

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
  bool intersect(const Ray& ray, double& outT);

  /*!
   * \brief Computes the intersection of a ray and all the
   * analytic approximations that compose this mesh
   *
   * \param ray The ray that is going to be tested against the mesh
   * \param t [out] The nearest intersection found with the
   * analytic approximations of this mesh
   *
   * \return True if an intersection was found, false otherwise
   */
  bool intersectAnalytic(const Ray& ray, double& t) const;

  /*!
   * \brief Computes the intersection of a ray and all the
   * parts of this mesh that don't have an analytic approximation.
   *
   * \param ray The ray
   * \param t [out] The closest intersection found with the
   * analytic approximations of this mesh
   *
   * \return True if an intersection was found, false otherwise
   */
  bool intersectNonAnalytic(const Ray& ray, double& t) const;

  void SmoothNormals();

  // Constructors from core classes
  explicit Mesh(const Box&);
  explicit Mesh(const SimpleMesh& simpleMesh);

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
  return vertexIndices;
}

/*!
\brief Return the set of normal indexes.
*/
inline std::vector<int> Mesh::NormalIndexes() const
{
  return normalIndices;
}

/*!
\brief Get the vertex index of a given triangle.
\param t Triangle index.
\param i Vertex index.
*/
inline int Mesh::VertexIndex(int t, int i) const
{
  return vertexIndices.at(t * 3 + i);
}

/*!
\brief Get the normal index of a given triangle.
\param t Triangle index.
\param i Normal index.
*/
inline int Mesh::NormalIndex(int t, int i) const
{
  return normalIndices.at(t * 3 + i);
}

/*!
\brief Get a triangle.
\param i Index.
\return The triangle.
*/
inline Triangle Mesh::GetTriangle(int i) const
{
  return Triangle(vertices.at(vertexIndices.at(i * 3 + 0)), vertices.at(vertexIndices.at(i * 3 + 1)), vertices.at(vertexIndices.at(i * 3 + 2)));
}

inline std::vector<Triangle*>* Mesh::GetTriangles() const
{
    std::vector<Triangle*>* triangles = new std::vector<Triangle*>;
    triangles->reserve(this->VertexIndexes().size() / 3);

    for(int vertexIndex = 0; vertexIndex < VertexIndexes().size(); vertexIndex += 3)
    {
        triangles->push_back(new Triangle(this->vertices.at(this->vertexIndices.at(vertexIndex + 0)),
                                      this->vertices.at(this->vertexIndices.at(vertexIndex + 1)),
                                      this->vertices.at(this->vertexIndices.at(vertexIndex + 2))));
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
  return vertices[vertexIndices[t * 3 + v]];
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
  return int(vertexIndices.size()) / 3;
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


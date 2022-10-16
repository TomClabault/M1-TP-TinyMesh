#ifndef SIMPLE_MESHES_H
#define SIMPLE_MESHES_H

#include <vector>
#include <iostream>

#include "mathematics.h"

class SimpleMesh
{
public:
    /*!
     * \brief Returns the number of vertices of the mesh
     * \return The number of vertices of the mesh
     */
    unsigned int VerticesCount() const;

    /*!
     * \brief Returns the vertex at index \p
     * index of the mesh
     * \return The vertex at index \p index
     * of the mesh
     */
    Vector Vertex(int index) const;

    /*!
     * \brief Returns the vertex index found at the position
     * \p index of the indices array of the mesh
     * \param index The index
     * \return The vertex index found at the
     * position \p index of the indices
     * array of the mesh
     */
    int VertexIndex(int index) const;

    /*!
     * \brief Returns the number of normals of the mesh
     * \return The number of normals of the mesh
     */
    unsigned int NormalsCount() const;

    /*!
     * \brief Returns the normal at the index \p index of the mesh
     * \param index The index
     * \return The normal at the index \p index of the mesh
     */
    Vector Normal(int index) const;

    /*!
     * \brief Returns the normal index found at the position
     * \p index of the normal indices array of the mesh
     * \param index The index
     * \return The normal index found at the
     * position \p index of the normal indices
     * array of the mesh
     */
    int NormalIndex(int index) const;

    /*!
     * \brief Returns the number of indicies
     * \return The number of indicies
     */
    unsigned int IndicesCount() const;

protected:
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<int> indices;
    std::vector<int> normalIndices;
};

class Icosphere : public SimpleMesh
{
public:
    /*!
     * \brief Icosphere
     * \param radius
     * \param subdivisions
     */
    Icosphere(double radius, int subdivisions);

    double Radius();
    int Subdivisions();

    /*!
     * \brief Subdivides the current icosphere one time (quadruples
     * the number of triangles)
     */
    void subdivide();

public:
    static const Vector baseVertices[12];
    static const int baseIndices[60];

private:
    void initBaseIcosphere(double radius);

private:
    double radius;
    int subdivisions = 0;
};

class Torus : public SimpleMesh
{
public:
    /*!
     * \brief Initializes a torus
     *
     * \param innerRadius The inner radius (radius of the "eatable" part
     * of the doughnut) of the torus
     * \param outerRadius The outer radius (radius of the "empty" part in
     * the middle) of the torus
     * \param ringCount The number of subdivision for the "outer ring"
     * of the torus. The higer the subdivision, the smoother the "empty" ring
     * inside the torus
     * \param ringsSubdivisions The number of subdivision for the inner rings
     * of the torus. The higher the subdivison, the smoother the "eatable"
     * part of the doughnut
     */
    Torus(double innerRadius, double outerRadius, int ringCount = 4, int ringsSubdivisions = 5);
};

class Capsule : public SimpleMesh
{
public:
    Capsule(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions);

private:
    double radius;
    double cylinderHeight;

    int cylinderSubdivisions;
    int sphereHeightSubdivisions;
};

class Cylinder : public SimpleMesh
{
public:
    Cylinder(double radius, double height, int heightSubdivisions, int cylinderSubdivisions);
};

#endif // SIMPLE_MESHES_H

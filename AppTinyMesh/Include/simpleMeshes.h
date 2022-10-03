#ifndef SIMPLE_MESHES_H
#define SIMPLE_MESHES_H

#include <vector>
#include <iostream>

#include "mathematics.h"

class SimpleMesh
{
public:
    unsigned int VerticesCount() const;
    Vector Vertex(int index) const;
    int VertexIndex(int index) const;

    unsigned int NormalsCount() const;
    Vector Normal(int index) const;
    int NormalIndex(int index) const;

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
    Icosphere();
    Icosphere(int subdivisions);

    /*!
     * \brief Subdivides the current icosphere one time (quadruples
     * the number of triangles)
     */
    void subdivide();

public:
    static const Vector baseVertices[12];
    static const int baseIndices[60];
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

#endif // SIMPLE_MESHES_H

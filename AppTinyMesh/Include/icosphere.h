#ifndef ICOSPHERE_H
#define ICOSPHERE_H

#include <vector>
#include <iostream>

#include "mathematics.h"

class Icosphere
{
public:
    Icosphere();
    Icosphere(int subdivisions);

    int VerticesCount() const;
    Vector Vertex(int index) const;
    int VertexIndex(int index) const;

    int NormalsCount() const;
    Vector Normal(int index) const;
    int NormalIndex(int index) const;

    int IndicesCount() const;

private:
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<int> indices;
    std::vector<int> normalIndices;

public:
    static const Vector baseVertices[12];
    static const int baseIndices[60];
};

#endif // ICOSPHERE_H

#include "icosphere.h"

const float X=.525731112119133606f;
const float Z=.850650808352039932f;
const float N=0.f;

//From https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
const Vector Icosphere::baseVertices[12] =
{
    Vector(-X,N,Z), Vector(X,N,Z), Vector(-X,N,-Z), Vector(X,N,-Z),
    Vector(N,Z,X), Vector(N,Z,-X), Vector(N,-Z,X), Vector(N,-Z,-X),
    Vector(Z,X,N), Vector(-Z,X, N), Vector(Z,-X,N), Vector(-Z,-X, N)
};

//From https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
const int Icosphere::baseIndices[60] =
{
  0,4,1,0,9,4,9,5,4,4,5,8,4,8,1,
  8,10,1,8,3,10,5,3,8,5,2,3,2,7,3,
  7,10,3,7,6,10,7,11,6,11,0,6,0,1,6,
  6,1,10,9,0,11,9,11,2,9,2,5,7,2,11
};

Icosphere::Icosphere()
{
    for(int i = 0; i < 12; i++)
        vertices.push_back(Icosphere::baseVertices[i]);

    for(int i = 0; i < 60; i++)
        indices.push_back(Icosphere::baseIndices[i]);

    //Computing normals
    for(int i = 0; i < 60; i += 3)
    {
        int index0 = Icosphere::baseIndices[i];
        int index1 = Icosphere::baseIndices[i + 1];
        int index2 = Icosphere::baseIndices[i + 2];

        Vector normal = (Icosphere::baseVertices[index2] - Icosphere::baseVertices[index0]) / (Icosphere::baseVertices[index1] - Icosphere::baseVertices[index0]);

        normals.push_back(normal);
        normalIndices.push_back(i / 3);
        normalIndices.push_back(i / 3);
        normalIndices.push_back(i / 3);
    }
}

Icosphere::Icosphere(int subdivisions) : Icosphere()
{
    std::cout << "Not implemented";
}

unsigned int Icosphere::VerticesCount() const
{
    return vertices.size();
}

Vector Icosphere::Vertex(int index) const
{
    return vertices.at(index);
}

int Icosphere::VertexIndex(int index) const
{
    return indices.at(index);
}

unsigned int Icosphere::NormalsCount() const
{
    return normals.size();
}

Vector Icosphere::Normal(int index) const
{
    return normals.at(index);
}

unsigned int Icosphere::IndicesCount() const
{
    return indices.size();
}

int Icosphere::NormalIndex(int index) const
{
    return normalIndices.at(index);
}

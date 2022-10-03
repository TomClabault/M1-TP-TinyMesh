#include "icosphere.h"
#include <unordered_map>

//Cache used to avoid the duplication of vertices when subdividing the icosphere
std::unordered_map<int, int> midPointsCache;

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
    midPointsCache = std::unordered_map<int, int>();

    for(int i = 0; i < subdivisions; i++)
        this->subdivide();
}

int getPointKey(const int pointAIndex, const int pointBIndex)
{
    //The Cantor's pairing function allows us to generate a unique key for the hash map
    //based on the two indices of the vertices
    return ((pointAIndex + pointBIndex) * (pointAIndex + pointBIndex + 1) / 2) + std::min(pointAIndex, pointBIndex);
}

/*!
 * \brief Checks wheter the point that is in between the two points
 * represented by the given indices has already been created or not.
 *
 * \param pointAIndex The index of the first point in the 'indices' array
 * \param pointBIndex The index of the second point in the 'indices' array
 * \param midPointIndex The index of the point that is in between the two other
 * points. The index will only be used if the middle point doesn't alredy exist
 * (i.e. if the middle point needs to be added to the cache). If the middle point
 * already exists, this parameter will be ignored
 *
 * \return If the middle point has already been created before, returns the
 * index of the point in the 'indices' array. Returns -1 otherwise.
 */
int midPointExists(const int pointAIndex, const int pointBIndex, const int midPointIndex)
{
    int midPointKey = getPointKey(pointAIndex, pointBIndex);
    auto midPoint = midPointsCache.find(midPointKey);
    if(midPoint != midPointsCache.end())//The middle point already exists
        return midPoint->second;

    //The middle point doesn't exist, we're adding it's index to the cache for
    //further references
    midPointsCache[midPointKey] = midPointIndex;

    //Returning -1 as we haven't found the middle point in the cache
    return -1;
}

void Icosphere::subdivide()
{
    std::vector<int> newIndices;
    std::vector<Vector> newNormals;
    std::vector<int> newNormalsIndices;

    int currentIndicesCount = IndicesCount();
    int currentVerticesCount = VerticesCount();
    int createdVerticesCount = 0;

    for(int i = 0; i < currentIndicesCount; i += 3)
    {
        int vertex1Index = indices[i + 0];
        int vertex2Index = indices[i + 1];
        int vertex3Index = indices[i + 2];

        Vector vertex1 = vertices[vertex1Index];
        Vector vertex2 = vertices[vertex2Index];
        Vector vertex3 = vertices[vertex3Index];

        int vertex12Index = currentVerticesCount + createdVerticesCount;
        int vertex23Index = currentVerticesCount + createdVerticesCount;
        int vertex31Index = currentVerticesCount + createdVerticesCount;

        //vertex12Index = midPointExists(vertex1Index, vertex2Index, vertex12Index);
        //vertex23Index = midPointExists(vertex2Index, vertex3Index, vertex23Index);
        //vertex31Index = midPointExists(vertex3Index, vertex1Index, vertex31Index);
        int vertex12Key = getPointKey(vertex1Index, vertex2Index);
        int vertex23Key = getPointKey(vertex2Index, vertex3Index);
        int vertex31Key = getPointKey(vertex3Index, vertex1Index);

        bool vertex12Exists = false;
        bool vertex23Exists = false;
        bool vertex31Exists = false;

        auto vertex12Find = midPointsCache.find(vertex12Key);
        if(vertex12Find != midPointsCache.end())
        {
            vertex12Index = vertex12Find->second;
            vertex12Exists = true;
        }
        else
        {
            std::cout << vertex12Index << std::endl;
            midPointsCache[vertex12Key] = vertex12Index;

            vertex23Index++;
            vertex31Index++;
        }

        auto vertex23Find = midPointsCache.find(vertex23Key);
        if(vertex23Find != midPointsCache.end())
        {
            vertex23Index = vertex23Find->second;
            vertex23Exists = true;
        }
        else
        {
            std::cout << vertex23Index << std::endl;
            midPointsCache[vertex23Key] = vertex23Index;

            vertex31Index++;
        }

        auto vertex31Find = midPointsCache.find(vertex31Key);
        if(vertex31Find != midPointsCache.end())
        {
            vertex31Index = vertex31Find->second;
            vertex31Exists = true;
        }
        else
        {
            std::cout << vertex31Index << std::endl;
            midPointsCache[vertex31Key] = vertex31Index;
        }

        //Getting the middle points from the already-existing-vertices array or by computing
        //the middle point based on whether we had already computed this point before
        Vector vertex12 = (!vertex12Exists) ? Normalized((vertex1 + vertex2) / 2) : vertices[vertex12Index];
        Vector vertex23 = (!vertex23Exists) ? Normalized((vertex2 + vertex3) / 2) : vertices[vertex23Index];
        Vector vertex31 = (!vertex31Exists) ? Normalized((vertex3 + vertex1) / 2) : vertices[vertex31Index];

        if(!vertex12Exists)//The middle point of vertex1 and vertex2 hadn't been computed before
        {
            createdVerticesCount++;
            vertices.push_back(vertex12);
        }

        if(!vertex23Exists)
        {
            createdVerticesCount++;
            vertices.push_back(vertex23);
        }

        if(!vertex31Exists)
        {
            createdVerticesCount++;
            vertices.push_back(vertex31);
        }

        //Adding the newly created vertices


        //Adding the new connections between the vertices to
        //make the triangles
        newIndices.push_back(vertex1Index);
        newIndices.push_back(vertex12Index);
        newIndices.push_back(vertex31Index);

        newIndices.push_back(vertex12Index);
        newIndices.push_back(vertex2Index);
        newIndices.push_back(vertex23Index);

        newIndices.push_back(vertex23Index);
        newIndices.push_back(vertex3Index);
        newIndices.push_back(vertex31Index);

        newIndices.push_back(vertex12Index);
        newIndices.push_back(vertex23Index);
        newIndices.push_back(vertex31Index);

        //Adding the new normals
        newNormals.push_back((vertex31 - vertex1) / (vertex12 - vertex1));
        newNormals.push_back((vertex12 - vertex2) / (vertex23 - vertex2));
        newNormals.push_back((vertex23 - vertex3) / (vertex31 - vertex3));
        newNormals.push_back((vertex31 - vertex12) / (vertex23 - vertex12));

        //Adding the normals indices
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 3; k++)
                newNormalsIndices.push_back((i/3) * 4 + j);
    }

    indices = newIndices;
    normals = newNormals;
    normalIndices = newNormalsIndices;
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

#include "realtime.h"
#include "simpleMeshes.h"

#define _USE_MATH_DEFINES

#include <cmath>
#include <unordered_map>

unsigned int SimpleMesh::VerticesCount() const
{
    return vertices.size();
}

Vector SimpleMesh::Vertex(int index) const
{
    return vertices.at(index);
}

int SimpleMesh::VertexIndex(int index) const
{
    return indices.at(index);
}

unsigned int SimpleMesh::NormalsCount() const
{
    return normals.size();
}

Vector SimpleMesh::Normal(int index) const
{
    return normals.at(index);
}

unsigned int SimpleMesh::IndicesCount() const
{
    return indices.size();
}

int SimpleMesh::NormalIndex(int index) const
{
    return normalIndices.at(index);
}

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

void Icosphere::initBaseIcosphere(double radius)
{
    for(int i = 0; i < 12; i++)
        vertices.push_back(Icosphere::baseVertices[i] * radius);

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

Icosphere::Icosphere(double radius, int subdivisions)
{
    midPointsCache = std::unordered_map<int, int>();

    this->radius = radius;

    initBaseIcosphere(radius);


    RenderingProfiler profiler;
    profiler.Init();
    for(int i = 0; i < subdivisions; i++)
        this->subdivide();
    profiler.Update();

    std::cout << profiler.msPerFrame << "ms[" << profiler.framePerSecond << "FPS]" << std::endl;
}

double Icosphere::Radius()
{
    return this->radius;
}

int Icosphere::Subdivisions()
{
    return this->subdivisions;
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
    this->subdivisions++;

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
            midPointsCache[vertex31Key] = vertex31Index;

        //Getting the middle points from the already-existing-vertices array or by computing
        //the middle point based on whether we had already computed this point before
        Vector vertex12 = (!vertex12Exists) ? Normalized((vertex1 + vertex2) / 2) : vertices[vertex12Index];
        Vector vertex23 = (!vertex23Exists) ? Normalized((vertex2 + vertex3) / 2) : vertices[vertex23Index];
        Vector vertex31 = (!vertex31Exists) ? Normalized((vertex3 + vertex1) / 2) : vertices[vertex31Index];

        vertex12 *= this->radius;
        vertex23 *= this->radius;
        vertex31 *= this->radius;

        vertex12 = Normalized(vertex12);
        vertex23 = Normalized(vertex23);
        vertex31 = Normalized(vertex31);

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
        newNormals.push_back(Normalized((vertex31 - vertex1) / (vertex12 - vertex1)));
        newNormals.push_back(Normalized((vertex12 - vertex2) / (vertex23 - vertex2)));
        newNormals.push_back(Normalized((vertex23 - vertex3) / (vertex31 - vertex3)));
        newNormals.push_back(Normalized((vertex31 - vertex12) / (vertex23 - vertex12)));

        //Adding the normals indices
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 3; k++)
                newNormalsIndices.push_back((i/3) * 4 + j);
    }

    indices = newIndices;
    normals = newNormals;
    normalIndices = newNormalsIndices;
}

Torus::Torus(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions)
{
    double outerAngleSubdiv = (2 * M_PI) / ringCount;
    double innerAngleSubdiv = (2 * M_PI) / ringsSubdivisions;

    //Computing the vertices
    for(int ring = 0; ring < ringCount; ring++)
    {
        double outerAngle = ring * outerAngleSubdiv;

        for(int ringSubdiv = 0; ringSubdiv < ringsSubdivisions; ringSubdiv++)
        {
            double innerAngle = ringSubdiv * innerAngleSubdiv;

            double x = std::cos(outerAngle) * (outerRadius + std::cos(innerAngle) * innerRadius);
            double y = std::sin(outerAngle) * (outerRadius + std::cos(innerAngle) * innerRadius);
            double z = innerRadius * std::sin(innerAngle);

            vertices.push_back(Vector(x, y, z));
        }
    }

    //Computing the indices and the normals
    int normalsCreated = 0;
    for(int ring = 0; ring < ringCount; ring++)
    {
        for(int ringSubdiv = 0; ringSubdiv < ringsSubdivisions; ringSubdiv++)
        {
            int vertex1Index = ring * ringsSubdivisions + ringSubdiv;
            int vertex2Index = ((ring + 1) % ringCount) * ringsSubdivisions + ringSubdiv;
            int vertex3Index = ring * ringsSubdivisions + (ringSubdiv +  1) % ringsSubdivisions;
            int vertex4Index = ((ring + 1) % ringCount) * ringsSubdivisions + (ringSubdiv + 1) % ringsSubdivisions;

            indices.push_back(vertex1Index);
            indices.push_back(vertex3Index);
            indices.push_back(vertex4Index);

            indices.push_back(vertex1Index);
            indices.push_back(vertex4Index);
            indices.push_back(vertex2Index);

            normals.push_back((vertices[vertex4Index] - vertices[vertex1Index]) / (vertices[vertex3Index] - vertices[vertex1Index]));
            normals.push_back((vertices[vertex2Index] - vertices[vertex1Index]) / (vertices[vertex4Index] - vertices[vertex1Index]));

            normalIndices.push_back(normalsCreated);
            normalIndices.push_back(normalsCreated);
            normalIndices.push_back(normalsCreated);
            normalIndices.push_back(normalsCreated + 1);
            normalIndices.push_back(normalsCreated + 1);
            normalIndices.push_back(normalsCreated + 1);

            normalsCreated += 2;
        }
    }
}

void Capsule::computeSphereRing(double deltaY, int deltaIndex, int ringIndex, double ringIncrement, double ringSubdivIncrement, double localRadius)
{
    //TODO si on est au ring tout en bas, il ne suffit de placer que un seul point vu qu'ils sont tous confondus
    double y = std::sin(M_PI / 2 * ringIndex * ringIncrement) + deltaY;

    for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
    {
        double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;
        double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;

        Vector vertex = Vector(x, y, z);
        this->vertices.push_back(vertex);

        if(ringIndex < sphereHeightSubdivisions - 1)
        {
            int index0 = deltaIndex + ringIndex * cylinderSubdivisions + ringSubdiv;
            int index1 = deltaIndex + ringIndex * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;
            int index2 = deltaIndex + (ringIndex + 1) * cylinderSubdivisions + ringSubdiv;
            int index3 = deltaIndex + (ringIndex + 1) * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;

            this->indices.push_back(index0);
            this->indices.push_back(index1);
            this->indices.push_back(index3);

            this->indices.push_back(index0);
            this->indices.push_back(index3);
            this->indices.push_back(index2);

            this->normalIndices.push_back(0);
            this->normalIndices.push_back(0);
            this->normalIndices.push_back(0);
            this->normalIndices.push_back(0);
            this->normalIndices.push_back(0);
            this->normalIndices.push_back(0);
            this->normals.push_back(Vector(0, 1, 0));
            this->normals.push_back(Vector(0, 1, 0));
            this->normals.push_back(Vector(0, 1, 0));
            this->normals.push_back(Vector(0, 1, 0));
            this->normals.push_back(Vector(0, 1, 0));
            this->normals.push_back(Vector(0, 1, 0));
        }
    }
}

Capsule::Capsule(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions)
{
    this->radius = radius;
    this->cylinderHeight = cylinderHeight;
    this->cylinderSubdivisions = cylinderSubdivisions;
    this->sphereHeightSubdivisions = sphereHeightSubdivisions;

    RenderingProfiler profiler;
    profiler.Init();

    std::cout << "Radius: " << radius << std::endl;
    std::cout << "Cylinder subdiv: " << cylinderSubdivisions << std::endl;
    std::cout << "Sphere height subdiv: " << sphereHeightSubdivisions << std::endl;

    //Generating the first bottom sphere cap of the capsule
    double ringIncrement = 1 / (sphereHeightSubdivisions - 1.0);
    double ringSubdivIncrement = 1 / (double)cylinderSubdivisions;
    double deltaY = -1.0;

    for(int ringIndex = 0; ringIndex < sphereHeightSubdivisions; ringIndex++)
    {
        double localRadius = std::sin(M_PI / 2 * ringIndex * ringIncrement);

        double y = 1 - std::cos(M_PI / 2 * ringIndex * ringIncrement) + deltaY;

        for(int ringSubdiv = 0; ringSubdiv < sphereHeightSubdivisions; ringSubdiv++)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;

            Vector vertex = Vector(x, y, z);
            this->vertices.push_back(vertex);

            if(ringIndex < sphereHeightSubdivisions - 1)
            {
                int index0 = ringIndex * cylinderSubdivisions + ringSubdiv;
                int index1 = ringIndex * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;
                int index2 = (ringIndex + 1) * cylinderSubdivisions + ringSubdiv;
                int index3 = (ringIndex + 1) * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;

                this->indices.push_back(index0);
                this->indices.push_back(index1);
                this->indices.push_back(index3);

                this->indices.push_back(index0);
                this->indices.push_back(index3);
                this->indices.push_back(index2);

                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
            }
        }
    }

    //To account for all the indices we have alredy generated for the first cap sphere
    int deltaIndex = sphereHeightSubdivisions * cylinderSubdivisions;
    double cylinderHeightIncrement = 1 / (cylinderHeightSubdivions - 1.0);
    //Generating the cylinder
    for(int cylinderRingIndex = 0; cylinderRingIndex < cylinderHeightSubdivions; cylinderRingIndex++)
    {
        for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement);
            double y = cylinderRingIndex * cylinderHeightIncrement * cylinderHeight;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement);

            this->vertices.push_back(Vector(x, y, z));

            if(cylinderRingIndex < cylinderHeightSubdivions - 1)
            {
                int index0 = deltaIndex + cylinderRingIndex * cylinderSubdivisions + ringSubdiv;
                int index1 = deltaIndex + cylinderRingIndex * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;
                int index4 = deltaIndex + (cylinderRingIndex + 1) * cylinderSubdivisions + ringSubdiv;
                int index5 = deltaIndex + (cylinderRingIndex + 1) * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;

                this->indices.push_back(index0);
                this->indices.push_back(index1);
                this->indices.push_back(index5);

                this->indices.push_back(index0);
                this->indices.push_back(index5);
                this->indices.push_back(index4);

                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
            }
        }
    }

    deltaY += cylinderHeight + 1.0;
    deltaIndex += cylinderHeightSubdivions * cylinderSubdivisions;
    for(int ringIndex = 0; ringIndex < sphereHeightSubdivisions; ringIndex++)
    {
        double localRadius = std::cos(M_PI / 2 * ringIndex * ringIncrement);

        //TODO si on est au ring tout en bas, il ne suffit de placer que un seul point vu qu'ils sont tous confondus
        double y = std::sin(M_PI / 2 * ringIndex * ringIncrement) + deltaY;

        for(int ringSubdiv = cylinderSubdivisions - 1; ringSubdiv >= 0; ringSubdiv--)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius * radius;

            Vector vertex = Vector(x, y, z);
            this->vertices.push_back(vertex);

            if(ringIndex < sphereHeightSubdivisions - 1)
            {
                int index0 = deltaIndex + ringIndex * cylinderSubdivisions + ringSubdiv;
                int index1 = deltaIndex + ringIndex * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;
                int index2 = deltaIndex + (ringIndex + 1) * cylinderSubdivisions + ringSubdiv;
                int index3 = deltaIndex + (ringIndex + 1) * cylinderSubdivisions + (ringSubdiv + 1) % cylinderSubdivisions;

                this->indices.push_back(index0);
                this->indices.push_back(index1);
                this->indices.push_back(index3);

                this->indices.push_back(index0);
                this->indices.push_back(index3);
                this->indices.push_back(index2);

                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normalIndices.push_back(0);
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
                this->normals.push_back(Vector(0, 1, 0));
            }
        }
    }

    profiler.Update();
    std::cout << profiler.msPerFrame << "ms[" << profiler.framePerSecond << "FPS]" << std::endl;
}

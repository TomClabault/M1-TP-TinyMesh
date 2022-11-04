#include "analyticApproximations.h"
#include "qmath.h"
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

unsigned long long int SimpleMesh::MemorySize() const
{
    unsigned long long int bytes = 0;

    bytes += this->indices.capacity() * sizeof(int);
    bytes += this->normalIndices.capacity() * sizeof(int);
    bytes += this->normals.capacity() * sizeof(double) * 3;//3 double for a Vector
    bytes += this->vertices.capacity() * sizeof(double) * 3;//3 double for a Vector

    return bytes;
}

AnalyticApproximation* SimpleMesh::GetAnalyticApproximation() const
{
    return this->analyticApproximation;
}

void SimpleMesh::transformVerticesAndNormals(const Vector& translation, const Matrix& rotation, const Matrix& scale)
{
    //Translate then rotate then scale.
    for(Vector& vertex : this->vertices)
    {
        vertex = vertex * _scale;
        vertex = vertex * _rotation;
        vertex += _translation;
    }

    for(Vector& normal : this->normals)
    {
        normal = normal * _scale.GetInverse().GetTranspose();
        normal = normal * _rotation.GetInverse().GetTranspose();
    }
}

//Cache used to avoid the duplication of vertices when subdividing the icosphere
std::unordered_map<int, int> midPointsCache;

const double X=.525731112119133606f;
const double Z=.850650808352039932f;
const double N=0.f;

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

Icosphere::Icosphere(const Vector& center, double radius, int subdivisions) : Icosphere(center, Matrix::Homothety(radius, radius, radius), Matrix::RotationX(0), subdivisions) {}

Icosphere::Icosphere(const Vector& translation, const Matrix& rotation, const Matrix& scale, int subdivisions) : SimpleMesh(translation, rotation, scale)
{
    midPointsCache = std::unordered_map<int, int>();

    initBaseIcosphere(1);

    for(int i = 0; i < subdivisions; i++)
        this->subdivide();

    transformVerticesAndNormals(translation, rotation, scale);

    SimpleMesh::analyticApproximation = new AnalyticSphere(translation, rotation, scale);
}

/*!
 * \brief Computes the Cantor's pairing function output of the two given integer.
 * Note that this implementation is meant to work interchangeably for A and B,
 * i.e. cantor(A, B) = cantor(B, A)
 *
 * \param A Integer A
 * \param B Integer B
 *
 * \return
 */
int cantor(const int A, const int B)
{
    //The Cantor's pairing function allows us to generate a unique key for the hash map
    //based on the two indices of the vertices
    return ((A + B) * (A + B + 1) / 2) + std::min(A, B);
}

/*!
 * \brief Checks wheter the vertex in between the two points (middle-vertex)
 * represented by the given indices has already been created or not.
 *
 * \param pointAIndex The index (in the 'indices' array) of the first vertex
 * \param pointBIndex The index (in the 'indices' array) of the second vertex
 * \param[in, out] midPointIndex This parameter is used as an input to the function
 * in case the "middle-vertex" is not found in the cache. In this case
 * the given index will be used to add the middle-vertex index to the cache.
 * If the middle-vertex already exists in the cache when this function is called,
 * this parameter will be used as a container for the index of the already-exsiting
 * middle-vertex
 *
 * \return True if the middle-vertex alredy existed in the cache, false otherwise
 */
bool checkVertexCacheAndUpdate(const int pointAIndex, const int pointBIndex, int& midPointIndex)
{
    int midPointKey = cantor(pointAIndex, pointBIndex);

    auto midPoint = midPointsCache.find(midPointKey);
    if(midPoint != midPointsCache.end())//The middle point already exists
    {
        midPointIndex = midPoint->second;

        return true;
    }

    midPointsCache[midPointKey] = midPointIndex;
    return false;
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

        bool vertex12Exists = checkVertexCacheAndUpdate(vertex1Index, vertex2Index, vertex12Index);
        if(!vertex12Exists)
        {
            //We need to offset the indices of the next vertices that may be created
            //since we just created one. If we don't offset the indices, they're all
            //going to be the same
            vertex23Index++;
            vertex31Index++;
        }
        bool vertex23Exists = checkVertexCacheAndUpdate(vertex2Index, vertex3Index, vertex23Index);
        if(!vertex23Exists)
            vertex31Index++;
        bool vertex31Exists = checkVertexCacheAndUpdate(vertex3Index, vertex1Index, vertex31Index);

        //Getting the middle points from the already-existing-vertices array or by computing
        //the middle point based on whether we had already computed this point before
        Vector vertex12 = Normalized((!vertex12Exists) ? (vertex1 + vertex2) / 2 : vertices[vertex12Index]);
        Vector vertex23 = Normalized((!vertex23Exists) ? (vertex2 + vertex3) / 2 : vertices[vertex23Index]);
        Vector vertex31 = Normalized((!vertex31Exists) ? (vertex3 + vertex1) / 2 : vertices[vertex31Index]);

        //We created one new vertex per each middle-vertex that didn't exist before
        createdVerticesCount += !vertex12Exists + !vertex23Exists + !vertex31Exists;

        if(!vertex12Exists)//The middle point of vertex1 and vertex2 hadn't been computed before
            vertices.push_back(vertex12);

        if(!vertex23Exists)
            vertices.push_back(vertex23);

        if(!vertex31Exists)
            vertices.push_back(vertex31);


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

    this->subdivisions++;
}

Torus::Torus(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions) : Torus(Vector(0, 0, 0), Matrix::RotationX(0), Matrix::Homothety(1, 1, 1), innerRadius, outerRadius, ringCount, ringsSubdivisions) {}

Torus::Torus(const Vector& translation, const Matrix& rotation, const Matrix& scale, double innerRadius, double outerRadius, int ringCount, int ringSubdivisions) : SimpleMesh(translation, rotation, scale)
{
    initBaseTorus(innerRadius, outerRadius, ringCount, ringSubdivisions);

    transformVerticesAndNormals(translation, rotation, scale);
}

void Torus::initBaseTorus(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions)
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

            normals.push_back(Normalized((vertices[vertex4Index] - vertices[vertex1Index]) / (vertices[vertex3Index] - vertices[vertex1Index])));
            normals.push_back(Normalized((vertices[vertex2Index] - vertices[vertex1Index]) / (vertices[vertex4Index] - vertices[vertex1Index])));

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

Capsule::Capsule(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions) : Capsule(Vector(0, 0, 0), Matrix::RotationX(0), Matrix::Homothety(radius, cylinderHeight, radius), cylinderHeightSubdivions, cylinderSubdivisions, sphereHeightSubdivisions) {}

Capsule::Capsule(const Vector& translation, const Matrix& rotation, const Matrix& scale, int cylinderHeightSubdivisions, int cylinderSubdivisions, int sphereHeightSubdivions) : SimpleMesh(translation, rotation, scale)
{
    initBaseCapsule(cylinderHeightSubdivisions, cylinderSubdivisions, sphereHeightSubdivions);

    transformVerticesAndNormals(translation, rotation, scale);
}

void Capsule::initBaseCapsule(int cylinderHeightSubdivisions, int cylinderSubdivisions, int sphereHeightSubdivisions)
{
    this->cylinderSubdivisions = cylinderSubdivisions;
    this->sphereHeightSubdivisions = sphereHeightSubdivisions;

    double ringIncrement = 1 / (sphereHeightSubdivisions - 1.0);
    double ringSubdivIncrement = 1 / (double)cylinderSubdivisions;
    double deltaY = -1.0;

    //Generating the first bottom sphere cap of the capsule
    for(int ringIndex = 0; ringIndex < sphereHeightSubdivisions; ringIndex++)
    {
        double localRadius = std::sin(M_PI / 2 * ringIndex * ringIncrement);

        double y = 1 - std::cos(M_PI / 2 * ringIndex * ringIncrement) + deltaY;

        for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius;

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
            }
        }
    }

    //To account for all the indices we have alredy generated for the first cap sphere
    int deltaIndex = sphereHeightSubdivisions * cylinderSubdivisions;
    double cylinderHeightIncrement = 1 / (cylinderHeightSubdivisions - 1.0);

    //Generating the cylinder
    for(int cylinderRingIndex = 0; cylinderRingIndex < cylinderHeightSubdivisions; cylinderRingIndex++)
    {
        for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement);
            double y = cylinderRingIndex * cylinderHeightIncrement;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement);

            this->vertices.push_back(Vector(x, y, z));

            if(cylinderRingIndex < cylinderHeightSubdivisions - 1)
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
            }
        }
    }

    deltaY += 2.0;
    deltaIndex += cylinderHeightSubdivisions * cylinderSubdivisions;
    for(int ringIndex = 0; ringIndex < sphereHeightSubdivisions; ringIndex++)
    {
        double localRadius = std::cos(M_PI / 2 * ringIndex * ringIncrement);

        //TODO si on est au ring tout en bas, il ne suffit de placer que un seul point vu qu'ils sont tous confondus
        double y = std::sin(M_PI / 2 * ringIndex * ringIncrement) + deltaY;

        for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
        {
            double x = std::cos(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius;
            double z = std::sin(2 * M_PI * ringSubdiv * ringSubdivIncrement) * localRadius;

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
            }
        }
    }

    int normalsCreated = 0;

    //Computing the normals
    for(size_t index = 0; index < this->indices.size(); index += 3)
    {
        int index0 = this->indices[index + 0];
        int index1 = this->indices[index + 1];
        int index2 = this->indices[index + 2];

        Vector normal = Normalized((this->vertices[index2] - this->vertices[index0]) / (this->vertices[index1] - this->vertices[index0]));
        this->normals.push_back(normal);

        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);

        normalsCreated++;
    }
}

Cylinder::Cylinder(const Vector& center, double radius, double height, int heightSubdivisions, int cylinderSubdivisions) : Cylinder(center, Matrix::RotationX(0), Matrix::Homothety(radius, height, radius), heightSubdivisions, cylinderSubdivisions) {}

Cylinder::Cylinder(const Vector& translation, const Matrix& rotation, const Matrix& scale, int heightSubdivisions, int cylinderSubdivisions) : SimpleMesh(translation, rotation, scale)
{
    SimpleMesh::analyticApproximation = new AnalyticCylinder(translation, rotation, scale);

    initBaseCylinder(heightSubdivisions, cylinderSubdivisions);

    transformVerticesAndNormals(translation, rotation, scale);
}

void Cylinder::initBaseCylinder(int heightSubdivisions, int cylinderSubdivisions)
{

    //First point at the middle of the bottom circle of the cylinder
    this->vertices.push_back(Vector(0, 0, 0));

    //Generating the bottom ring of the cylinder
    double ringIncrement = 1.0 / cylinderSubdivisions;
    for(int cylinderSubdiv = 0; cylinderSubdiv < cylinderSubdivisions; cylinderSubdiv++)
    {
        double x = std::cos(2 * M_PI * ringIncrement * cylinderSubdiv);
        double y = 0;
        double z = std::sin(2 * M_PI * ringIncrement * cylinderSubdiv);

        this->vertices.push_back(Vector(x, y, z));
    }


    //Generating the cylinder
    double heightIncrement = 1.0 / (heightSubdivisions + 1.0);
    for(int ringIndex = 0; ringIndex < heightSubdivisions; ringIndex++)
    {
        for(int cylinderRing = 0; cylinderRing < cylinderSubdivisions; cylinderRing++)
        {
            double x = std::cos(2 * M_PI * ringIncrement * cylinderRing);
            double y = heightIncrement * (ringIndex + 1);
            double z = std::sin(2 * M_PI * ringIncrement * cylinderRing);

            this->vertices.push_back(Vector(x, y, z));
        }
    }

    //Copying the bottom circle of the cylinder
    //to have the top circle. Adding in the height
    //of the cylinder to effectively have the
    //"top circle"
    for(int i = 0; i <= cylinderSubdivisions; i++)
        this->vertices.push_back(Vector(this->vertices[i] + Vector(0, 1, 0)));


    //Generating the indices and the normals of both
    //the bottom and top circle of the cylinder
    int verticesCount = this->vertices.size();
    int normalsCreated = 0;
    for(int i = 0; i < cylinderSubdivisions; i++)
    {
        Vector normal;

        int index0 = 0;
        int index1 = i + 1;
        int index2 = (i + 1) % cylinderSubdivisions + 1;

        this->indices.push_back(index0);//Vertex at the center of the bottom ring of the cylinder
        this->indices.push_back(index1);
        this->indices.push_back(index2);
        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);
        normalsCreated++;

        normal = Normalized((this->vertices[index1] - this->vertices[index0]) / (this->vertices[index2] - this->vertices[index0]));
        this->normals.push_back(normal);


        index0 = verticesCount - cylinderSubdivisions - 1;
        index1 = verticesCount - cylinderSubdivisions - 1 + i + 1;
        index2 = index1 + 1;
        if(index2 == verticesCount)
            index2 = index0 + 1;

        this->indices.push_back(index0);//Vertex at the center of the top ring of the cylinder
        this->indices.push_back(index1);
        this->indices.push_back(index2);

        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);
        this->normalIndices.push_back(normalsCreated);
        normalsCreated++;

        normal = Normalized((this->vertices[index2] - this->vertices[index0]) / (this->vertices[index1] - this->vertices[index0]));
        this->normals.push_back(normal);
    }

    //Generating the indices and the normals of
    //the triangles of the 'body' of the cylinder
    //(not the caps of the cylinder)
    for(int i = 0; i < heightSubdivisions + 1; i++)
    {
        for(int ringSubdiv = 0; ringSubdiv < cylinderSubdivisions; ringSubdiv++)
        {
            int index0 = ringSubdiv + 1 + i * cylinderSubdivisions;
            int index1 = index0 + 1;
            if (ringSubdiv == cylinderSubdivisions - 1)
                index1 -= cylinderSubdivisions;

            int index2 = index0 + cylinderSubdivisions + (i == heightSubdivisions ? 1 : 0);
            int index3 = index2 + 1;
            if (ringSubdiv == cylinderSubdivisions - 1)
                index3 -= cylinderSubdivisions;

            this->indices.push_back(index0);
            this->indices.push_back(index1);
            this->indices.push_back(index3);

            Vector normal = Normalized((this->vertices[index3] - this->vertices[index0]) / (this->vertices[index1] - this->vertices[index0]));
            this->normals.push_back(normal);

            this->normalIndices.push_back(normalsCreated);
            this->normalIndices.push_back(normalsCreated);
            this->normalIndices.push_back(normalsCreated);
            normalsCreated++;

            this->indices.push_back(index0);
            this->indices.push_back(index3);
            this->indices.push_back(index2);

            normal = Normalized((this->vertices[index2] - this->vertices[index0]) / (this->vertices[index3] - this->vertices[index0]));
            this->normals.push_back(normal);

            this->normalIndices.push_back(normalsCreated);
            this->normalIndices.push_back(normalsCreated);
            this->normalIndices.push_back(normalsCreated);
            normalsCreated++;

        }
    }
}

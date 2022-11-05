#include "analyticApproximations.h"
#include "BVH.h"
#include "qmath.h" //For M_PI
#include "mesh.h"

#include <chrono>
#include <unordered_set>


/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), vertexIndices(indices)
{
    normals.resize(vertices.size(), Vector::Z);

    //We're adding the analytic approximation that corresponds to the vertices we just added to the mesh.
    //Because we added arbitrary vertices, there is no corresponding approximation so we're adding nullptr
    this->analyticApproxToVertexIndex.insert(std::make_pair(indices.size(), nullptr));
    this->_bvh = BVH(*this->GetTriangles());
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& vIndices, const std::vector<int>& nIndices) :vertices(vertices), normals(normals), vertexIndices(vIndices), normalIndices(nIndices)
{
    this->analyticApproxToVertexIndex.insert(std::make_pair(vIndices.size(), nullptr));
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
    vertices.reserve(nv);
    normals.reserve(nn);
    vertexIndices.reserve(nvi);
    normalIndices.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

AnalyticApproximation* Mesh::GetAnalyticApproxFromVertexIndex(int vertexIndex) const
{
    //Getting the analytic approximation that contains the vertex index number 'vertexIndex' using upper_bound
    auto analyticApproxIte = this->analyticApproxToVertexIndex.upper_bound(vertexIndex);

    if(analyticApproxIte == this->analyticApproxToVertexIndex.end())
        return nullptr;
    else
    {
        AnalyticApproximation* approx = analyticApproxIte->second;

        return approx;
    }
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
    // Initialize
    normals.resize(vertices.size(), Vector::Null);

    normalIndices = vertexIndices;

    // Accumulate normals
    for (size_t i = 0; i < vertexIndices.size(); i += 3)
    {
        Vector tn = Triangle(vertices[vertexIndices.at(i)], vertices[vertexIndices.at(i + 1)], vertices[vertexIndices.at(i + 2)]).AreaNormal();
        normals[normalIndices[i + 0]] += tn;
        normals[normalIndices[i + 1]] += tn;
        normals[normalIndices[i + 2]] += tn;
    }

    // Normalize
    for (size_t i = 0; i < normals.size(); i++)
    {
        Normalize(normals[i]);
    }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
    vertexIndices.push_back(a);
    normalIndices.push_back(na);
    vertexIndices.push_back(b);
    normalIndices.push_back(nb);
    vertexIndices.push_back(c);
    normalIndices.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
    vertexIndices.push_back(a);
    normalIndices.push_back(n);
    vertexIndices.push_back(b);
    normalIndices.push_back(n);
    vertexIndices.push_back(c);
    normalIndices.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
    // First triangle
    AddSmoothTriangle(a, na, b, nb, c, nc);

    // Second triangle
    AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
    AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
    if (vertices.size() == 0)
    {
        return Box::Null;
    }
    return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
    // Vertices
    vertices.resize(8);

    for (int i = 0; i < 8; i++)
    {
        vertices[i] = box.Vertex(i);
    }

    // Normals
    normals.push_back(Vector(-1, 0, 0));
    normals.push_back(Vector(1, 0, 0));
    normals.push_back(Vector(0, -1, 0));
    normals.push_back(Vector(0, 1, 0));
    normals.push_back(Vector(0, 0, -1));
    normals.push_back(Vector(0, 0, 1));

    // Reserve space for the triangle array
    vertexIndices.reserve(12 * 3);
    normalIndices.reserve(12 * 3);

    AddTriangle(0, 2, 1, 4);
    AddTriangle(1, 2, 3, 4);

    AddTriangle(4, 5, 6, 5);
    AddTriangle(5, 7, 6, 5);

    AddTriangle(0, 4, 2, 0);
    AddTriangle(4, 6, 2, 0);

    AddTriangle(1, 3, 5, 1);
    AddTriangle(3, 7, 5, 1);

    AddTriangle(0, 1, 5, 2);
    AddTriangle(0, 5, 4, 2);

    AddTriangle(3, 2, 7, 3);
    AddTriangle(6, 7, 2, 3);

    //Analytic approximation for the boxes not implemented
    this->analyticApproxToVertexIndex.insert(std::make_pair(vertexIndices.size(), nullptr));
}

/*!
\brief Creates a mesh based on a "simple mesh"

\param mesh The simple mesh
*/
Mesh::Mesh(const SimpleMesh& simpleMesh)
{
    // Vertices
    unsigned int verticesCount = simpleMesh.VerticesCount();
    vertices.resize(verticesCount);

    for (unsigned int i = 0; i < verticesCount; i++)
        vertices[i] = simpleMesh.Vertex(i);

    // Normals
    for(unsigned int i = 0; i < simpleMesh.NormalsCount(); i++)
        normals.push_back(simpleMesh.Normal(i));

    // Reserve space for the triangle array
    unsigned int indicesCount = simpleMesh.IndicesCount();
    vertexIndices.reserve(indicesCount);
    normalIndices.reserve(indicesCount);

    for(unsigned int  i = 0; i < indicesCount; i += 3)
        AddTriangle(simpleMesh.VertexIndex(i), simpleMesh.VertexIndex(i + 1), simpleMesh.VertexIndex(i + 2), simpleMesh.NormalIndex(i));

    AnalyticApproximation* meshApprox = simpleMesh.GetAnalyticApproximation();
    if(meshApprox != nullptr)
        this->analyticApproximations.push_back(simpleMesh.GetAnalyticApproximation());
    this->analyticApproxToVertexIndex.insert(std::make_pair(simpleMesh.IndicesCount(), meshApprox));
}

/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (size_t i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (size_t i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}

void Mesh::Scale(const Matrix& homothetyMatrix)
{
    for(Vector& vertex : this->vertices)
        vertex = vertex * homothetyMatrix;
}

void Mesh::Rotate(const Matrix& rotationMatrix)
{
    for(Vector& vertex : this->vertices)
        vertex = vertex * rotationMatrix;

    for(Vector& normal : this->normals)
        normal = Normalized(normal * rotationMatrix.GetInverse().GetTranspose());
}

void Mesh::Translate(const Vector& translationVector)
{
    for(Vector& vertex : this->vertices)
        vertex += translationVector;
}

void Mesh::Merge(const Mesh& secondMesh)
{
    int currentIndicesCount = this->vertexIndices.size();
    int currentVerticesCount = this->Vertexes();
    int currentNormalsCount = this->normals.size();

    for(Vector vertex : secondMesh.vertices)
        this->vertices.push_back(vertex);

    for(Vector normal : secondMesh.normals)
        this->normals.push_back(normal);

    for(int vertexIndex : secondMesh.VertexIndexes())
        this->vertexIndices.push_back(vertexIndex + currentVerticesCount);

    for(int normalIndex : secondMesh.NormalIndexes())
        this->normalIndices.push_back(normalIndex + currentNormalsCount);

    //Also adding all the approximations of the second mesh if they exist
    for(auto secondMeshApproxIte = secondMesh.analyticApproxToVertexIndex.begin();
             secondMeshApproxIte != secondMesh.analyticApproxToVertexIndex.end();
             secondMeshApproxIte++)
        this->analyticApproxToVertexIndex.insert(std::make_pair(secondMeshApproxIte->first + currentIndicesCount, secondMeshApproxIte->second));

    if(secondMesh.analyticApproximations.size() > 0)
        for(AnalyticApproximation* approx : secondMesh.analyticApproximations)
            this->analyticApproximations.push_back(approx);
}

void Mesh::SphereWarp(Sphere sphere)
{
    Vector sphereCenter = sphere.Center();
    double sphereRadius = sphere.Radius();

    for(Vector& vertex : this->vertices)
    {
        Vector directionToCenter = sphereCenter - vertex;
        double distanceToCenter = Distance(sphereCenter, vertex);

        if(distanceToCenter < sphereRadius)
        {
            double distanceToCenterNormalized = distanceToCenter / sphereRadius;

            double& x = distanceToCenterNormalized;
            double intensity = x < 0.5 ? (4 * x * x * x) : (1 - std::pow(-2 * x + 2, 3) / 2);

            vertex += directionToCenter * (1 - intensity);
        }
    }
}

//TODO rajouter un texte dans l'interface qui dit combien de temps ça a pris pour calculer l'AO
//TODO Ajouter une checkbox dans l'interface pour activer/Désactiver la BVH
void Mesh::accessibility(std::vector<Color>& accessibilityColors, double radius, int samples, double occlusionStrength, bool enableAnalyticIntersection, bool useBVH)
{
    std::unordered_set<int> alreadyComputedVertices;//Holds the index of the vertices
    //whose accessibility we already have computed. Useful not to compute
    //several time the accessibility of the same vertex

    double colorIncrement = 1.0 / samples;
    bool analyticIntersection = this->analyticApproximations.size() > 0 && enableAnalyticIntersection;

    size_t vertexIndexCount = this->VertexIndexes().size();
    for(size_t vertexIndex = 0; vertexIndex < vertexIndexCount; vertexIndex++)
    {
        double obstructedValue = 0;

        int vertexNumber = this->vertexIndices.at(vertexIndex);
        if(alreadyComputedVertices.find(vertexNumber) != alreadyComputedVertices.end())//We have already computed this vertex
            continue;//Skipping it
        else
            alreadyComputedVertices.insert(vertexNumber);

        Vector vertex = this->vertices.at(vertexNumber);
        Vector normal;

        if(analyticIntersection)
        {
            AnalyticApproximation* approx = this->GetAnalyticApproxFromVertexIndex(vertexIndex);
            if(approx != nullptr)
                normal = approx->GetNormalAt(vertex);
            else
                normal = this->normals.at(this->normalIndices.at(vertexIndex));
        }
        else
            normal = this->normals.at(this->normalIndices.at(vertexIndex));

        for(int sample = 0; sample < samples; sample++)
        {
            unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

            //Using xorshift96 to generate random numbers is faster than std::rand by 2 orders of magnitude
            Vector randomRayDirection = Normalized(Vector((Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
                                                          (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1,
                                                          (Math::xorshift96() / (double)max_unsigned_int) * 2 - 1));

            if(randomRayDirection * normal < 0)//If the random point we draw is below the surface
                randomRayDirection *= -1;//Flipping the random point for it to be above the surface

            //We're slightly shifting the origin of the ray in the
            //direction of the normal otherwise we will intersect ourselves
            Ray ray(vertex + normal * Math::EPSILON, randomRayDirection);

            double intersectionDistance;
            bool intersectionFound;

            if(analyticIntersection)
            {
                double intersectionDistNonAnalytic = INFINITY;

                //This is only going to intersect the analytic approximations of our mesh
                intersectionFound = this->intersectAnalytic(ray, intersectionDistance);
                //If our mesh is composed of analytic approximations and of other pieces that cannot
                //be approximated, we're going to have to intersect those as well
                intersectionFound  |= this->intersectNonAnalytic(ray, intersectionDistNonAnalytic);

                //We'll keep the smallest intersection distance that we've found between the intersection
                //with the analytic shapes and with the non-analytic shapes
                intersectionDistance = std::min(intersectionDistNonAnalytic, intersectionDistance);
            }
            else
            {
                if(useBVH)
                    intersectionFound = _bvh.intersect(ray, intersectionDistance);
                else
                    intersectionFound = this->intersect(ray, intersectionDistance);
            }

            //In front of the ray and within the given occlusion radius
            if(intersectionFound && intersectionDistance <= radius)
                obstructedValue += colorIncrement;
        }

        accessibilityColors.at(vertexNumber) = Color(1 - obstructedValue * occlusionStrength);
    }
}

void Mesh::GetInterStats(unsigned long long int& triangleInterTests, unsigned long long int& triangleEffectiveInters) const
{
    triangleInterTests = _triangleInterTests;
    triangleEffectiveInters = _triangleEffectiveInter;
}

bool Mesh::intersect(const Ray& ray, double& outT)
{
    bool found = false;
    double closestT = INFINITY;//Used to keep track of the closest

    //Looping through indices 3 by 3 to get the triangles
    for(size_t i = 0; i < this->vertexIndices.size(); i += 3)
    {
        int index1 = this->vertexIndices.at(i + 0);
        int index2 = this->vertexIndices.at(i + 1);
        int index3 = this->vertexIndices.at(i + 2);

        double t, u, v;

        _triangleInterTests++;
        if(Triangle::IntersectFromPoints(this->vertices.at(index1), this->vertices.at(index2), this->vertices.at(index3), ray, t, u, v))
        {
            if(t > 0 && t < closestT)//Intersection in front of the ray
            {
                closestT = t;
                outT = t;

                _triangleEffectiveInter++;
                found = true;
            }
        }
    }

    return found;
}

bool Mesh::intersectAnalytic(const Ray& ray, double& t) const
{
    bool found = false;

    double nearestT = INFINITY;
    for(AnalyticApproximation* approx : this->analyticApproximations)
    {
        if(approx->intersect(ray, t))
        {
            if(t < nearestT && t > 0)
            {
                found = true;

                nearestT = t;
            }
        }
    }

    t = nearestT;

    return found;
}

bool Mesh::intersectNonAnalytic(const Ray& ray, double& t) const
{
    bool found = false;
    double closestT = INFINITY;//Used to keep track of the closest
    //intersection we've found so far

    //Looping through the meshes (possibly merged into this mesh) of this mesh
    int previousEntryLastVertexIndex = 0;//This variable will hold the last vertex
    //index of the previous entry of the map. This will be useful to retrieve
    //the vertices of a mesh because by iterating on the map, we will only
    //have access to the last vertex index of the mesh, not its first one
    //so this is what this variable is here for
    for(auto mapEntry = this->analyticApproxToVertexIndex.begin();
             mapEntry != this->analyticApproxToVertexIndex.end();
             mapEntry++)
    {
        if(mapEntry->second == nullptr) //This means that we found a mesh that
        //doesn't have an analytic approximation
        {
            for(size_t i = previousEntryLastVertexIndex; i < mapEntry->first; i += 3)
            {
                int index1 = this->vertexIndices.at(i + 0);
                int index2 = this->vertexIndices.at(i + 1);
                int index3 = this->vertexIndices.at(i + 2);

                double interT, u, v;
                if(Triangle::IntersectFromPoints(this->vertices.at(index1), this->vertices.at(index2), this->vertices.at(index3), ray, interT, u, v))
                {
                    if(interT > 0 && interT < closestT)//Intersection in front of the ray
                    {
                        closestT = interT;
                        t = interT;

                        found = true;
                    }
                }
            }
        }

        //Keeping the number of the last vertex ofthe mesh we just iterated over
        previousEntryLastVertexIndex = mapEntry->first;
    }

    return found;
}


#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
    vertices.clear();
    normals.clear();
    vertexIndices.clear();
    normalIndices.clear();

    QFile data(filename);

    if (!data.open(QFile::ReadOnly))
        return;
    QTextStream in(&data);

    // Set of regular expressions : Vertex, Normal, Triangle
    QRegularExpression rexV("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
    QRegularExpression rexN("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
    //f 1/1/4 3/7/4 5/8/4 6/6/4
    QRegularExpression rexF4("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
    //f 6/6/5 5/8/5 4/9/5
    QRegularExpression rexF3("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");

    while (!in.atEnd())
    {
        QString line = in.readLine();
        QRegularExpressionMatch matchV = rexV.match(line);
        QRegularExpressionMatch matchN = rexN.match(line);
        QRegularExpressionMatch matchF4 = rexF4.match(line);
        QRegularExpressionMatch matchF3 = rexF3.match(line);
        if (matchV.hasMatch())//rexv.indexIn(line, 0) > -1)
        {
            Vector q = Vector(matchV.captured(1).toDouble(), matchV.captured(2).toDouble(), matchV.captured(3).toDouble());

            vertices.push_back(q);
        }
        else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
        {
            Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());

            normals.push_back(q);
        }
        else if (matchF4.hasMatch())
        {
            //f 3/3/2 2/2/2 4/4/2 5/5/2
            //C B D E
            //--> BEC, BDE --> 2, 4, 1 // 2, 3, 4

            //Quad face, 2 triangles
            vertexIndices.push_back(matchF4.captured(1).toInt() - 1);
            vertexIndices.push_back(matchF4.captured(3).toInt() - 1);
            vertexIndices.push_back(matchF4.captured(5).toInt() - 1);

            vertexIndices.push_back(matchF4.captured(1).toInt() - 1);
            vertexIndices.push_back(matchF4.captured(5).toInt() - 1);
            vertexIndices.push_back(matchF4.captured(7).toInt() - 1);


            normalIndices.push_back(matchF4.captured(2).toInt() - 1);
            normalIndices.push_back(matchF4.captured(4).toInt() - 1);
            normalIndices.push_back(matchF4.captured(6).toInt() - 1);

            normalIndices.push_back(matchF4.captured(2).toInt() - 1);
            normalIndices.push_back(matchF4.captured(6).toInt() - 1);
            normalIndices.push_back(matchF4.captured(8).toInt() - 1);
        }
        else if (matchF3.hasMatch())//rexF3.indexIn(line, 0) > -1)
        {
            vertexIndices.push_back(matchF3.captured(1).toInt() - 1);
            vertexIndices.push_back(matchF3.captured(3).toInt() - 1);
            vertexIndices.push_back(matchF3.captured(5).toInt() - 1);
            normalIndices.push_back(matchF3.captured(2).toInt() - 1);
            normalIndices.push_back(matchF3.captured(4).toInt() - 1);
            normalIndices.push_back(matchF3.captured(6).toInt() - 1);
        }
    }
    data.close();

    this->analyticApproxToVertexIndex.insert(std::make_pair(vertexIndices.size(), nullptr));
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
    QFile data(url);
    if (!data.open(QFile::WriteOnly))
        return;
    QTextStream out(&data);
    out << "g " << meshName << Qt::endl;
    for (size_t i = 0; i < vertices.size(); i++)
        out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
    for (size_t i = 0; i < normals.size(); i++)
        out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
    for (size_t i = 0; i < vertexIndices.size(); i += 3)
    {
        out << "f " << vertexIndices.at(i) + 1 << "//" << normalIndices.at(i) + 1 << " "
            << vertexIndices.at(i + 1) + 1 << "//" << normalIndices.at(i + 1) + 1 << " "
            << vertexIndices.at(i + 2) + 1 << "//" << normalIndices.at(i + 2) + 1 << " "
            << "\n";
    }
    out.flush();
    data.close();
}


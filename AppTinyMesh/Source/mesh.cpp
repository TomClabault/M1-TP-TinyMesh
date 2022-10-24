#include "mesh.h"

//For M_PI
#include "qmath.h"

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
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
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

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
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
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
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
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

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
}

/*!
\brief Creates a mesh based on a "simple mesh"

\param mesh The simple mesh
*/
Mesh::Mesh(const SimpleMesh& mesh)
{
  // Vertices
  unsigned int verticesCount = mesh.VerticesCount();
  vertices.resize(verticesCount);

  for (unsigned int i = 0; i < verticesCount; i++)
      vertices[i] = mesh.Vertex(i);

  // Normals
  for(unsigned int i = 0; i < mesh.NormalsCount(); i++)
      normals.push_back(mesh.Normal(i));

  // Reserve space for the triangle array
  unsigned int indicesCount = mesh.IndicesCount();
  varray.reserve(indicesCount);
  narray.reserve(indicesCount);

  for(unsigned int  i = 0; i < indicesCount; i += 3)
      AddTriangle(mesh.VertexIndex(i), mesh.VertexIndex(i + 1), mesh.VertexIndex(i + 2), mesh.NormalIndex(i));
}

/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
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
    int currentVerticesCount = this->Vertexes();
    int currentNormalsCount = this->normals.size();

    for(Vector vertex : secondMesh.vertices)
        this->vertices.push_back(vertex);

    for(Vector normal : secondMesh.normals)
        this->normals.push_back(normal);

    for(int vertexIndex : secondMesh.VertexIndexes())
        this->varray.push_back(vertexIndex + currentVerticesCount);

    for(int normalIndex : secondMesh.NormalIndexes())
        this->narray.push_back(normalIndex + currentNormalsCount);
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

//TODO remove
#define DEBUG_INDEX_COUNT 4 

void Mesh::accessibility(std::vector<Color>& accessibilityColors, double radius, int samples, double occlusionStrength)
{
    //TODO remove
    /*int debugIndexes[DEBUG_INDEX_COUNT] = { (11 + 20 * 3 - 4) * 3 + 75,
                             30 + 120 * 5,
                             30 + 120 * 10,
                             30 + 120 * 15 };*/

    double colorIncrement = 1.0 / samples;
    const double epsilon = 1.0e-3;

    //TODO remove
    /*Vector* randomSamples = (Vector*)malloc(sizeof(Vector) * samples);
    for (int i = 0; i < samples; i++)
        randomSamples[i] = Normalized(Vector((std::rand() / (double)RAND_MAX) * 2 - 1, (std::rand() / (double)RAND_MAX) * 2 - 1, (std::rand() / (double)RAND_MAX) * 2 - 1));*/

    //TODO ne pas recalculer 50 fois le même vertex. En bouclant sur les indices des vertex comme ça, on va recalculer l'accessibilité même pour des vertex partagés par plusieurs triangles
    for(int vertexIndex = 0; vertexIndex < this->varray.size(); vertexIndex++)
    {
        double obstructedValue = 0;

        Vector vertex = this->vertices.at(this->varray.at(vertexIndex));
        Vector normal = this->normals.at(this->narray.at(vertexIndex));

        //TODO remove
        /*for (int debugIndex = 0; debugIndex < DEBUG_INDEX_COUNT; debugIndex++)
        {
            if (vertexIndex == debugIndexes[debugIndex])
            {
                for (int k = 0; k < 5; k++)
                {
                    this->Merge(Mesh(Box(normal * k/3 + vertex, 0.025)));
                    for (int k2 = 0; k2 < 8; k2++)
                        accessibilityColors.push_back(Color(0, 0, 255));
                }
            }
        }*/

        for(int sample = 0; sample < samples; sample++)
        {
            Vector randomVec = Normalized(Vector((std::rand() / (double)RAND_MAX) * 2 - 1, (std::rand() / (double)RAND_MAX) * 2 - 1, (std::rand() / (double)RAND_MAX) * 2 - 1));
            if(randomVec * normal < 0)//If the random point we draw is below the surface
                randomVec = -1 * randomVec;//Flipping the random point for it to be above the surface

            //TODO remove
            /*for (int debugIndex = 0; debugIndex < DEBUG_INDEX_COUNT; debugIndex++)
            {
                if (vertexIndex == debugIndexes[debugIndex])
                {
                    this->Merge(Mesh(Box((vertex + randomRayDirection) * radius, 0.05)));
                    for (int k = 0; k < 8; k++)
                        accessibilityColors.push_back(Color(0, 255, 0));
                }
            }*/

            //We're slightly shifting the origin of the ray in the
            //direction of the normal otherwise we will intersect ourselves
            Ray ray(vertex + normal * epsilon, randomVec);

            double intersectionDistance;
            bool intersectionFound = this->intersect(ray, intersectionDistance);

            //In front of the ray and within the given occlusion radius
            if(intersectionFound && intersectionDistance > epsilon && intersectionDistance <= radius)
                obstructedValue += colorIncrement;
        }

        //Color occlusion in yellow
        //accessibilityColors.at(this->varray.at(vertexIndex)) = Color((int)(255 * obstructedValue), (int)(255 * obstructedValue), 0);

        accessibilityColors.at(this->varray.at(vertexIndex)) = Color(1 - obstructedValue * occlusionStrength);

        //Red occlusion
        //accessibilityColors.at(this->varray.at(vertexIndex)) = Color(1.0, 1 - obstructedValue * occlusionStrength, 1 - obstructedValue * occlusionStrength);
    }

    //TODO remove
    /*for (int debugIndex = 0; debugIndex < DEBUG_INDEX_COUNT; debugIndex++)
        accessibilityColors.at(varray.at(debugIndexes[debugIndex])) = Color(255, 0, 0);*/
}

bool Mesh::intersect(const Ray& ray, double& outT)
{
    bool found = false;
    double closestT = 10e64;//Used to keep track of the closest

    //Looping through indices 3 by 3 to get the triangles
    for(int i = 0; i < this->varray.size(); i += 3)
    {
        int index1 = this->varray.at(i + 0);
        int index2 = this->varray.at(i + 1);
        int index3 = this->varray.at(i + 2);

        double t, u, v;

        if(Triangle::IntersectFromPoints(this->vertices.at(index1), this->vertices.at(index2), this->vertices.at(index3), ray, t, u, v))
        {
            if(t > 0 && t < closestT)//Intersection in front of the ray
            {
                closestT = t;
                outT = t;

                found = true;
            }
        }
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
  varray.clear();
  narray.clear();

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
        varray.push_back(matchF4.captured(1).toInt() - 1);
        varray.push_back(matchF4.captured(3).toInt() - 1);
        varray.push_back(matchF4.captured(5).toInt() - 1);

        varray.push_back(matchF4.captured(1).toInt() - 1);
        varray.push_back(matchF4.captured(5).toInt() - 1);
        varray.push_back(matchF4.captured(7).toInt() - 1);


        narray.push_back(matchF4.captured(2).toInt() - 1);
        narray.push_back(matchF4.captured(4).toInt() - 1);
        narray.push_back(matchF4.captured(6).toInt() - 1);

        narray.push_back(matchF4.captured(2).toInt() - 1);
        narray.push_back(matchF4.captured(6).toInt() - 1);
        narray.push_back(matchF4.captured(8).toInt() - 1);
    }
    else if (matchF3.hasMatch())//rexF3.indexIn(line, 0) > -1)
    {
      varray.push_back(matchF3.captured(1).toInt() - 1);
      varray.push_back(matchF3.captured(3).toInt() - 1);
      varray.push_back(matchF3.captured(5).toInt() - 1);
      narray.push_back(matchF3.captured(2).toInt() - 1);
      narray.push_back(matchF3.captured(4).toInt() - 1);
      narray.push_back(matchF3.captured(6).toInt() - 1);
    }
  }
  data.close();
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
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}


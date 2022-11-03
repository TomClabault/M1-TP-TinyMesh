#ifndef SIMPLE_MESHES_H
#define SIMPLE_MESHES_H

#include <vector>
#include <iostream>

#include "analyticApproximations.h"
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

    /*!
     * \brief Returns the number of bytes occupied by the instance in memory
     * \return The number of bytes occupied by the instance in memory
     */
    unsigned long long int MemorySize() const;

    /*!
     * \brief Returns the analytic shape that can approximate this simple mesh.
     * \return The analytic mesh that approximates this mesh
     */
    AnalyticApproximation* GetAnalyticApproximation() const;

    /*!
     * \brief Transforms the vertices of the mesh. The vertices are first
     * translated, then rotated then scaled
     * \param translation The translation
     * \param rotation The rotation matrix
     * \param scale The homotethy matrix
     */
    void transformVertices(Vector translation, Matrix rotation, Matrix scale);

protected:
    SimpleMesh(Vector translation, Matrix rotation, Matrix scale) : _translation(translation), _rotation(rotation), _scale(scale) {};

protected:
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<int> indices;
    std::vector<int> normalIndices;

    AnalyticApproximation* analyticApproximation = nullptr;

    Vector _translation;
    Matrix _rotation;
    Matrix _scale;
};

class Icosphere : public SimpleMesh
{
public:
    Icosphere(double radius, int subdivisions) : Icosphere(Vector(0, 0, 0), radius, subdivisions) {};

    /*!
     * \brief Creates an icosphere of radius
     * \param center The center of the icosphere
     * \p radius and that has been subdivided
     * \p subdivisions times
     * \param radius The radius of the icopshere
     * \param subdivisions The number of subdivisions
     * of the icosphere
     */
    Icosphere(const Vector& center, double radius, int subdivisions);

    /*!
     * \brief Initializes an icosphere using different linear transformations
     * \param translation The translation vector
     * \param scale The homothety matrix used to scale the sphere
     * \param rotation The rotation matrix
     * \param subdivisions The number of subdivisions
     * of the icosphere
     */
    Icosphere(const Vector& translation, const Matrix& scale, const Matrix& rotation, int subdivisions);

    /*!
     * \brief Subdivides the current icosphere one time (quadruples
     * the number of triangles)
     */
    void subdivide();

private:
    /*!
     * \brief Creates an icosphere of radius \p
     * radius without any subdivisions steps
     * (only the base vertices are present)
     * \param radius The radius of the icosphere
     * to create
     */
    void initBaseIcosphere(double radius);

public:
    static const Vector baseVertices[12];
    static const int baseIndices[60];

private:
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

    /*!
     * \brief Initializes a torus
     *
     * \param translation The translation that will be applied to the torus
     * \param rotation The rotation that will be applied to the torus
     * \param scale The scaling that will be applied to the torus
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
    Torus(const Vector& translation, const Matrix& rotation, const Matrix& scale, double innerRadius, double outerRadius, int ringCount, int ringSubdivisions);

private:
    /*!
     * \brief Creates a torus at the origin of the coordinates system,
     * not scaled and without any rotation
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
    void initBaseTorus(double innerRadius, double outerRadius, int ringCount, int ringSubdivisions);
};

class Capsule : public SimpleMesh
{
public:
    /*!
     * \brief Create a capsule composed of one cylinder
     * and two caps at its ends.
     *
     * \param radius The radius of the capsule
     * \param cylinderHeight The height of the cylinder
     * of the capsule
     * \param cylinderHeightSubdivions The number of
     * subdivisions in height of the cylinder. This parameter
     * has no visual impact but can be useful for later
     * operations applied to the capsule
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder and the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the capsule.
     * The lower, the rougher. With a value of 4 for example,
     * the capsule will start to look like an elongated cube.
     * \param sphereHeightSubdivisions The number of subdivisions
     * in height of the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the caps
     * (height-wise only).
     */
    Capsule(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions);

    /*!
     * \brief Create a capsule composed of one cylinder
     * and two caps at its ends.
     *
     * \param translation The translation to be applied to the
     * capsule
     * \param rotation The rotation to be applied to the capsule
     * \param scale The scaling to be applied to the capsule
     * \param radius The radius of the capsule
     * \param cylinderHeight The height of the cylinder
     * of the capsule
     * \param cylinderHeightSubdivions The number of
     * subdivisions in height of the cylinder. This parameter
     * has no visual impact but can be useful for later
     * operations applied to the capsule
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder and the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the capsule.
     * The lower, the rougher. With a value of 4 for example,
     * the capsule will start to look like an elongated cube.
     * \param sphereHeightSubdivisions The number of subdivisions
     * in height of the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the caps
     * (height-wise only).
     */
    Capsule(const Vector& translation, const Matrix& rotation, const Matrix& scale, int cylinderHeightSubdivisions, int cylinderSubdivisions, int sphereHeightSubdivions);

    /*!
     * \brief Initializes a capsule of height and radius 1 at the
     * origin of the world
     * \param cylinderHeightSubdivions The number of
     * subdivisions in height of the cylinder. This parameter
     * has no visual impact but can be useful for later
     * operations applied to the capsule
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder and the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the capsule.
     * The lower, the rougher. With a value of 4 for example,
     * the capsule will start to look like an elongated cube.
     * \param sphereHeightSubdivisions The number of subdivisions
     * in height of the caps at the ends of the cylinder.
     * The higher this parameter, the smoother the caps
     * (height-wise only).
     */
    void initBaseCapsule(int cylinderHeightSubdivisions, int cylinderSubdivisions, int sphereHeightSubdivisions);

private:
    int cylinderSubdivisions;
    int sphereHeightSubdivisions;
};

class Cylinder : public SimpleMesh
{
public:
    /*!
     * \brief Creates a cylinder with the given parameters
     *
     * \param radius The radius of the cylinder
     * \param height The height of the cylinder
     * \param heightSubdivisions The number of subdivisions
     * in height of the cylinder. This parameter has no
     * visual impact but can be useful for later operations
     * applied to the cylinder.
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder. The higher this parameter, the
     * smoother the cylinder. With a value of 4 for example,
     * the cylinder will look much like an elongated cube.
     */
    Cylinder(double radius, double height, int heightSubdivisions, int cylinderSubdivisions) : Cylinder(Vector(0, 0, 0), radius, height, heightSubdivisions, cylinderSubdivisions) {};

    /*!
     * \brief Creates a cylinder with the given parameters
     *
     * \param bottomDiskCenter The center of the bottom disk of the cylinder
     * \param radius The radius of the cylinder
     * \param height The height of the cylinder
     * \param heightSubdivisions The number of subdivisions
     * in height of the cylinder. This parameter has no
     * visual impact but can be useful for later operations
     * applied to the cylinder.
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder. The higher this parameter, the
     * smoother the cylinder. With a value of 4 for example,
     * the cylinder will look much like an elongated cube.
     */
    Cylinder(const Vector& bottomDiskCenter, double radius, double height, int heightSubdivisions, int cylinderSubdivisions);

    /*!
     * \brief Creates a cylinder with the given parameters
     *
     * \param translation The translation to be applied to the
     * capsule
     * \param rotation The rotation to be applied to the capsule
     * \param scale The scaling to be applied to the capsule
     * \param heightSubdivisions The number of subdivisions
     * in height of the cylinder. This parameter has no
     * visual impact but can be useful for later operations
     * applied to the cylinder.
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder. The higher this parameter, the
     * smoother the cylinder. With a value of 4 for example,
     * the cylinder will look much like an elongated cube.
     */
    Cylinder(const Vector& translation, const Matrix& rotation, const Matrix& scale, int heightSubdivisions, int cylinderSubdivisions);

    /*!
     * \brief Initializes a base cylinder of radius and height 1 at the
     * origin of the world
     * \param heightSubdivisions The number of subdivisions
     * in height of the cylinder. This parameter has no
     * visual impact but can be useful for later operations
     * applied to the cylinder.
     * \param cylinderSubdivisions The number of subdivisions
     * of the cylinder. The higher this parameter, the
     * smoother the cylinder. With a value of 4 for example,
     * the cylinder will look much like an elongated cube.
     */
    void initBaseCylinder(int heightSubdivisions, int cylinderSubdivisions);
};

#endif // SIMPLE_MESHES_H

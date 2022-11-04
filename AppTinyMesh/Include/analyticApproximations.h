#pragma once

#include "matrix.h"
#include "ray.h"

class AnalyticApproximation
{
public:
    /*!
     * \brief Returns the normal to the surface of the approximation at a certain given vertex
     * \param vertex The vertex
     * \return The normal to the surface of the approximation at the given \p vertex.
     */
    virtual Vector GetNormalAt(const Vector& vertex) = 0;

    /*!
     * \brief Intersects the analytic shape with a given ray
     * \param ray The ray
     * \param t [out] The intersection distance between the
     * ray and the shape
     * \return True if an intersection with the ray occured,
     * false otherwise. This method should return false if an
     * intersection behind the ray occured
     */
    virtual bool intersect(const Ray& ray, double& t) = 0;

protected:
    /*!
     * \brief Protected constructor used by sub classes to initialize the translation, rotation and scale attributes
     * \param translation The translation vector used to translate the analytic approximation
     * \param rotation The rotation matrix
     * \param scale The homothety matrix
     */
    AnalyticApproximation(Vector translation, Matrix rotation, Matrix scale) :
        _translation(translation), _rotation(rotation), _scale(scale),
        _invRotation(rotation.GetInverse()), _invScale(scale.GetInverse()) {};

protected:
    Vector _translation;
    Matrix _rotation;//No rotation by default
    Matrix _scale;//No scale by default

    Matrix _invRotation;//Inverse of the rotate matrix
    Matrix _invScale;//Inverse the of the scale matrix
};

class AnalyticSphere : public AnalyticApproximation
{
public:
    /*!
     * \brief Initializes an analytic sphere
     * \param radius The radius of the sphere
     */
    AnalyticSphere(double radius) : AnalyticApproximation(Vector(0, 0, 0), Matrix::RotationX(0), Matrix::Homothety(Vector(radius, radius, radius))) {};

    /*!
     * \brief Initializes an analytic sphere
     * \param center The center of the sphere
     * \param radius The radius of the sphere
     */
    AnalyticSphere(Vector center, double radius) : AnalyticApproximation(center, Matrix::RotationX(0), Matrix::Homothety(Vector(radius, radius, radius))) {};

    /*!
     * \brief Initializes an analytic sphere
     * \param translation The translation vector used to translate the sphere
     * \param rotation The rotation matrix
     * \param scale The homothety matrix
     */
    AnalyticSphere(Vector translation, Matrix rotation, Matrix scale) : AnalyticApproximation(translation, rotation, scale) {};

    /*!
     * \brief Returns the center of the sphere
     * \return The center of the sphere
     */
    Vector Center() { return _translation; }

    /*!
     * \brief Returns the normal to the sphere at a given vertex.
     * This method assumes that the given vertex actually is on
     * the sphere
     * \param vertex The vertex
     * \return The normal to the sphere at the vertex
     */
    Vector GetNormalAt(const Vector& vertex) override;

    bool intersect(const Ray& ray, double& t) override;

    static void intersectionTests();
    static void normalAtTests();
    static void tests();

private:
    bool intersectBasic(const Ray& ray, double& t);
};

class AnalyticCylinder : public AnalyticApproximation
{
public:
    /*!
     * \brief Initializes an analytic cylinder
     * \param center The center of the bottom disk of the cylinder
     * \param radius The radius of the cylinder
     * \param height The height of the cylinder. This correspond
     * to the height of the top disk of the cylinder.
     */
    AnalyticCylinder(Vector center, double radius, double height) : AnalyticApproximation(center, Matrix::RotationX(0), Matrix::Homothety(radius, height, radius)) {};

    /*!
     * \brief Initializes an analytic cylinder
     * \param translation The translation vector used to translate the cylinder
     * \param rotation The rotation matrix
     * \param scale The homothety matrix
     */
    AnalyticCylinder(Vector translation, Matrix rotation, Matrix scale) : AnalyticApproximation(translation, rotation, scale) {};

    /*!
     * \brief Returns the normal to the cylinder at a given vertex.
     * This method does not check whether the vertex actually is
     * on the cylinder or not
     * \param vertex The vertex
     * \return The normal to the cylinder at the vertex
     */
    Vector GetNormalAt(const Vector& vertex) override;

    bool intersect(const Ray& ray, double& t) override;

    static void intersectionTests();
    static void getNormalAtTests();
    static void tests();

private:
    /*!
     * \brief Intersects a "standard" cylinder that is at the origin,
     * of height 1 and radius 1. This method is mainly used after
     * having inverse transformed the ray with the transformations
     * applied to the actual cylinder.
     * \param ray The inverse transformed ray
     * \param t The intersection distance to the cylinder
     * \return True if an intersection occured, false otherwise
     */
    bool intersectBasic(const Ray& ray, double& t);
};

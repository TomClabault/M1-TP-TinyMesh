#pragma once

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
};

class AnalyticSphere : public AnalyticApproximation
{
public:
    /*!
     * \brief Initializes an analytic sphere
     * \param center The center of the sphere
     * \param radius The radius of the sphere
     */
    AnalyticSphere(Vector center, double radius) : _center(center), _radius(radius) {};

    Vector Center() { return _center; }

    /*!
     * \brief Returns the normal to the sphere at a given vertex.
     * This method does not check whether the vertex actually is
     * on the sphere or not
     * \param vertex The vertex
     * \return The normal to the sphere at the vertex
     */
    Vector GetNormalAt(const Vector& vertex) override;

    bool intersect(const Ray& ray, double& t) override;

    static void intersectionTest();

private:
    Vector _center;

    double _radius;
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
    AnalyticCylinder(Vector center, double radius, double height) : _center(center), _radius(radius), _height(height) {};

    /*!
     * \brief Returns the normal to the cylinder at a given vertex.
     * This method does not check whether the vertex actually is
     * on the cylinder or not
     * \param vertex The vertex
     * \return The normal to the cylinder at the vertex
     */
    Vector GetNormalAt(const Vector& vertex) override;

    bool intersect(const Ray& ray, double& t) override;

    static void intersectionTest();

private:
    Vector _center;// <! The center of the bottom disk of the cylinder

    double _radius;
    double _height;
};

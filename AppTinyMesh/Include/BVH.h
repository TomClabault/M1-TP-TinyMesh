#ifndef BVH_H
#define BVH_H

#include "mathematics.h"
#include "ray.h"
#include "triangle.h"

#include <vector>

class Plane
{
public:
    Plane(const Vector& normal, double d) : _normal(normal), _d(d) {};

    /*!
     * \brief Returns the normal of the plane
     * \return The normal of the plane
     */
    Vector Normal() const;

    /*!
     * \brief Computes the distance to the
     * plane of a point
     * \param point The point
     * \return The distance to the plane
     * of the point
     */
    inline double PointDistance(const Vector& point) const;

    /*!
     * \brief Computes the intersection
     * of a ray and the plane
     * \param ray The ray
     * \param t [out] The distance to the
     * intersection along the ray's
     * direction from its origin
     * \return True if an intersection in
     * front of the ray was found, false
     * otherwise
     */
    inline bool Intersect(const Ray& ray, double& t) const;

private:
    Vector _normal;//Normal to the plane
    double _d;//Distance of the plane from the origin
};

class BoundingVolume
{
public:
    BoundingVolume();
    BoundingVolume(const Triangle& triangle);
    BoundingVolume(const std::vector<Triangle>& triangles);
    BoundingVolume(const std::vector<Vector>& vertexes, const std::vector<int>& indices);

    /*!
     * \brief Extends the current bounds of this volume to fit the
     * given \p volume
     * \param volume The volume that is going to extend this one
     */
    void extend(const BoundingVolume& volume);

    /*!
     * \brief Computes the intersection between a
     * ray and this bounding volume
     * \param ray The ray
     * \param tNear [out] The closest intersection distance in the ray's
     * direction from the ray's origin with the volume
     * \param tFar [out] The furthest intersection distance in the ray's
     * direction from the ray's origin with the volume
     * \return True if an intersection was found
     * in front of the ray, false otherwise
     */
    bool intersect(const Ray& ray, double& tNear, double& tFar) const;

    /*!
     * \brief Returns the centroid of this volume
     * \return The centroid of this volume
     */
    Vector Centroid() const;

    //TODO remove if not used
//    /*!
//     * \brief Computes the minimum and the maximum point of the AABB of
//     * this volume
//     * \param min [out] The minimum point
//     * \param max [out] The maximum point
//     */
//    void GetMinMax(Vector& min, Vector& max) const;

    //TODO remove
//    /*!
//     * \brief Returns the point representing the lower left vertex of the AABB surrounding all the objects of this volume
//     * \return The minimum point of the AABB of this volume
//     */
//    Vector Min() const;

//    /*!
//     * \brief Returns the point representing the upper right vertex of the AABB surrounding all the objects of this volume
//     * \return The maximum point of the AABB of this volume
//     */
//    Vector Max() const;

private:
    /*!
     * \brief Computes all the near and far distances of the 7 planes
     * used by the volume using the given triangles
     *
     * \param triangles The triangles
     */
    void computeDNearDFar(const std::vector<Triangle>& triangles);
    void computeDNearDFar(const Triangle* triangle);

public:
    //The 7 planes that are going to
    //represent the bounding volume.
    //This is basically an extended version of
    //a bounding box. Instead of using only 3
    //planes perpendicular to the X, Y and Z axis,
    //we're using 7 planes to get tighter bounding volumes.
    static const Plane planes[7];

private:
    double _dNears[7];//The nears 'd' (as in the
    //equation of a plane Ax+By+cZ-d=0) of the
    //7 planes (see BoundingVolume::planes) of
    //this bounding volumes
    double _dFars[7];//The fars 'd'

    //TODO remove
    //Vector _min, _max;
};

class OctreeNode
{
public:
    OctreeNode();

    /*!
     * \brief Constructs an octree node using its lower and upper point
     * and its depth
     * \param min The lower point of the AABB of the node
     * \param max The upper point of the AABB of the node
     * \param depth The current depth of the node
     */
    OctreeNode(int centroidNumber, Vector min, Vector max, int depth);//TODO remove centroid number

    /*!
     * \brief Returns the nodes children of this node
     * \return The nodes children of this node
     */
    OctreeNode* ChildrenNodes() const;

    /*!
     * \brief Returns the bounding volume of this node
     * \return The bounding volume of this node
     */
    BoundingVolume GetBoundingVolume() const;

    /*!
     * \brief Is the node a leaf ?
     * \return True if the node is a leaf of the tree, false otherwise
     */
    bool IsLeaf() const;

    /*!
     * \brief Returns the triangles contained in the node
     * \return The triangles contained in the node if it's a leaf. An empty
     * vector otherwise
     */
    std::vector<const Triangle*> Triangles() const;

    /*!
     * \brief Computes the bounding volume of this node and updates the
     * _boundingVolume attribute of this instance
     * \return The bounding volume of the node. This can be useful for
     * further computations
     */
    BoundingVolume& computeVolume();

    /*!
     * \brief Analog to insertTriangles but only inserts one triangle
     */
    //void insertTriangle(const Triangle* triangle, int maxChildren, int maxDepth);
    void insertTriangle(const Triangle* triangle, int maxChildren, int maxDepth, int triangleNumber = -1);//TODO remove triangleNumber

    /*!
     * \brief Inserts the given triangles in the node. If the
     * number of triangles to insert exceeds \p maxChildren,
     * this node will be subdivided into 8 octants and the
     * triangles will be inserted in the newly created nodes
     * \param triangles The triangles to insert into the node.
     * \param maxChildren The maximum number of triangles that
     * can be inserted in this node. If the actual number of
     * triangle exceeds this parameter, the current node will
     * be subdivided into 8 new nodes
     * \param maxDepth The maximum depth
     */
    void insertTriangles(const std::vector<Triangle*>& triangles, int maxChildren, int maxDepth);

    //TODO remove if unused
//    /*!
//     * \brief Computes the intersection between a ray and this octree node
//     * \param ray The ray
//     * \param t The distance to the closest intersection with an object
//     * of the node. Garbage value that shouldn't be considered if no
//     * intersection is found.
//     * \return True if an intersection was found in front of the ray,
//     * false otherwise
//     */
//    bool intersect(const Ray& ray, double& t) const;

private:
    bool _isLeaf;

    int _centroidNumber;//TODO remove

    int _depth;//Current depth of the node in the tree

    Vector _min, _max, _center;

    BoundingVolume _boundingVolume;//Bounding volume of the node itself
    std::vector<const Triangle*> _triangles;//Triangles contained in the node

    //8 children of the node if they exist
    OctreeNode* _childrenNodes = nullptr;
};

class Octree
{
public:
    Octree();
    Octree(Vector min, Vector max);

    /*!
     * \brief Returns the root of this octree
     * \return The root of this octree.
     */
    OctreeNode* Root() const;

    /*!
     * \brief Constructs the octree with the given triangles
     * and the given constraints
     * \param triangles The triangles to insert into the octree
     * \param maxChildren The constraint on the maximum number
     * of children that can be in a leaf of the octree
     * \param maxDepth The constraint on the maximum depth
     * of a leaf in the octree
     */
    void construct(const std::vector<Triangle*>& triangles, int maxChildren, int maxDepth);

    /*!
     * \brief Computes the volumes of all the nodes of the octree
     */
    void computeVolumes();

    /*!
     * \brief Computes the intersection between a ray and this octree
     * \param ray The ray
     * \param t The distance to the closest intersection with an object
     * of the octree. Garbage value that shouldn't be considered if no
     * intersection is found.
     * \return True if an intersection was found in front of the ray,
     * false otherwise
     */
    bool intersect(const Ray& ray, double& t) const;

private:
    Vector _min, _max, _center;

    OctreeNode* _root;

    struct OctreeQueueNode
    {
        OctreeQueueNode(OctreeNode* node, double interDist) : _node(node), _nodeInterDistance(interDist) {};

        //Redefining the < operator for later use in a priority_queue.
        //We're writing it a > b because we want the priority queue to order the smallest distance first in the queue
        friend bool operator < (const OctreeQueueNode& a, const OctreeQueueNode& b) {return a._nodeInterDistance > b._nodeInterDistance;};

        OctreeNode* _node;
        double _nodeInterDistance;
    };
};

class BVH
{
public:

    /*!
     * \brief Empty BVH
     */
    BVH();

    /*!
     * \brief Constructs a BVH from a list of triangles
     * \param mesh The triangles
     * \param maxLeafChildren The maximum number of triangles that is
     * allowed to sit in a volume of the BVH
     * \param The maximum allowed depth of the nodes of the hierarchy
     */
    BVH(const std::vector<Triangle*>& triangles, int maxLeafChildren = 10, int maxDepth = 15);

    /*!
    * \brief Returns the underlying octree
    * \return The octree used by the BVH
    */
   Octree* GetOctree() const;

    /*!
     * \brief Computes the intersection of a ray and the BVH
     * \param ray The ray
     * \param t The distance to the closest intersection with an object
     * of the BVH. Garbage value that shouldn't be considered if no
     * intersection is found.
     * \return True if an intersection was found in front of the ray,
     * false otherwise
     */
    bool intersect(const Ray& ray, double& t) const;

    void static BVHTests();
    void static GetInterStats(unsigned long long int& boundingVolumesTests, unsigned long long int& triangleInterTests, unsigned long long int& triangleEffectiveInters);

private:
    //TODO remove if not used
    std::vector<BoundingVolume> _boundingVolumes;

    Octree* _octree = nullptr;
};

#endif // BVH_H

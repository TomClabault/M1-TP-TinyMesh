#include "BVH.h"

#include <queue>

BoundingVolume::BoundingVolume()
{
    for(int i = 0; i < 7; i++)
        dNears[i] = INFINITY;
    for(int i = 0; i < 7; i++)
        dFars[i] = -INFINITY;
}

BoundingVolume::BoundingVolume(const Triangle& triangle)
{
    _triangles.reserve(1);

    _triangles.push_back(&triangle);

    computeDNearDFar();
}

BoundingVolume::BoundingVolume(const std::vector<Triangle>& triangles)
{
    _triangles.reserve(triangles.size());

    for(const Triangle& triangle : triangles)
        _triangles.push_back(&triangle);

    computeDNearDFar();
}

BoundingVolume::BoundingVolume(const std::vector<Vector>& vertexes, const std::vector<int>& indices)
{
    _triangles.reserve(indices.size() / 3);

    for(int vertexIndex  = 0; vertexIndex < indices.size(); vertexIndex += 3)
        _triangles.push_back(new Triangle(vertexes.at(indices.at(vertexIndex)),
                                          vertexes.at(indices.at(vertexIndex + 1)),
                                          vertexes.at(indices.at(vertexIndex + 2))));

    computeDNearDFar();
}

void BoundingVolume::extend(const BoundingVolume& volume)
{
    for(int i = 0; i < 7; i++)
    {
        dNears[i] = std::min(dNears[i], volume.dNears[i]);
        dFars[i] = std::max(dFars[i], volume.dFars[i]);
    }
}

void BoundingVolume::computeDNearDFar()
{
    for(const Triangle* triangle : _triangles)
    {
        for(int triangleVertex = 0; triangleVertex < 3; triangleVertex++)
        {
            Vector vertex = (*triangle)[triangleVertex];

            //Computing dNears and dFars
            for(int plane = 0; plane < 7; plane++)
            {
                dNears[plane] = BoundingVolume::planes[plane].Normal() * vertex;
                dFars[plane] = BoundingVolume::planes[plane].Normal() * vertex;
            }
        }
    }
}

void BoundingVolume::GetMinMax(Vector& min, Vector& max) const
{
    for(const Triangle* triangle : _triangles)
    {
        for(int triangleVertex = 0; triangleVertex < 3; triangleVertex++)
        {
            Vector vertex = (*triangle)[triangleVertex];

            for(int component = 0; component < 3; component++)
            {
                //Updating the minimum and maximum point of the AABB of this bounding volume
                min[triangleVertex] = std::min(min[triangleVertex], vertex[component]);
                max[triangleVertex] = std::max(max[triangleVertex], vertex[component]);
            }
        }
    }
}

//TODO remove tFar if never used
bool BoundingVolume::intersect(const Ray& ray, double& tNear, double& tFar) const
{
    double outNear = -INFINITY;
    double outFar = INFINITY;

    for(int i = 0; i < 7; i++)
    {
        double ND = BoundingVolume::planes[i].Normal() * ray.Direction();
        double NO = BoundingVolume::planes[i].Normal() * ray.Origin();

        double near = (dNears[i] - NO) / ND;
        double far = (dFars[i] - NO) / ND;

        if(ND < 0)
            std::swap(near, far);

        if(outNear < near)
            outNear = near;
        if(outFar > far)
            outFar = far;

        if(outNear > outFar)
            return false;
    }

    tNear = outNear;
    tFar = outFar;

    if(tFar >= 0)
        return true;
    else
        return false;
}

Vector BoundingVolume::Centroid() const
{
    //The average of the nears and fars of the axis aligned planes
    return Vector((dNears[0] + dFars[0]) / 2, (dNears[1] + dFars[1]) / 2, (dNears[2] + dFars[2]) / 2);
}

//TODO remove
//Vector BoundingVolume::Min() const
//{
//    return _min;
//}

//Vector BoundingVolume::Max() const
//{
//    return _max;
//}

Vector Plane::Normal() const
{
    return _normal;
}

double Plane::PointDistance(const Vector& point) const
{
    return this->_normal * point;
}

bool Plane::Intersect(const Ray& ray, double& t) const
{
    return (t = this->PointDistance(ray.Origin())) < 0;
}

const Plane BoundingVolume::planes[7] = {Plane(Vector(1, 0, 0), 0),
                              Plane(Vector(0, 1, 0), 0),
                              Plane(Vector(0, 0, 1), 0),
                              Plane(Vector(std::sqrt(3) / 3, std::sqrt(3) / 3, std::sqrt(3) / 3), 0),
                              Plane(Vector(-std::sqrt(3) / 3, std::sqrt(3) / 3, std::sqrt(3) / 3), 0),
                              Plane(Vector(-std::sqrt(3) / 3, -std::sqrt(3) / 3, std::sqrt(3) / 3), 0),
                              Plane(Vector(std::sqrt(3) / 3, -std::sqrt(3) / 3, std::sqrt(3) / 3), 0)};

OctreeNode::OctreeNode() {}

OctreeNode::OctreeNode(const Vector min, const Vector max, int depth)
{
    _depth = depth;
    _min = min;
    _max = max;
    _center = (min + max) / 2;

    _isLeaf = true;
};

OctreeNode* OctreeNode::ChildrenNodes() const
{
    return _childrenNodes;
}

BoundingVolume OctreeNode::GetBoundingVolume() const
{
    return _boundingVolume;
}

bool OctreeNode::IsLeaf() const
{
    return _isLeaf;
}

std::vector<const Triangle*> OctreeNode::Triangles() const
{
    return _triangles;
}

BoundingVolume& OctreeNode::computeVolume()
{
    if(_isLeaf)
    {
        for(const Triangle* triangle : _triangles)
            _boundingVolume.extend(BoundingVolume(*triangle));
    }
    else
    {
        //The node is not a leaf, we're going to recursively compute the volume of all the subnodes
        for(int i = 0; i < 8; i++)
            _boundingVolume.extend(_childrenNodes[i].computeVolume());
    }

    return _boundingVolume;
}

void OctreeNode::insertTriangle(const Triangle* triangle, int maxChildren, int maxDepth)
{
    if(_isLeaf)
    {
        //If they are too many objects in this node and that there
        //is still depth left for more nodes in the hierarchy,
        //we're going to subdivide the node into 8
        if(_triangles.size() >= maxChildren && _depth < maxDepth)
        {
            _isLeaf = false;
            _childrenNodes = new OctreeNode[8];

            double halfWidth = (_min[0] + _max[0]) / 2;
            double halfHeight = (_min[1] + _max[1]) / 2;
            double halfDepth = (_min[2] + _max[2]) / 2;

            /**
             *     +-----------+
                  /  4  /  5  /|
                 /-----/-----/ |
                /  0  /  1  /|5|
               +----------+ 1|/|
               |  0  | 1  | /|7|
               |----------|/ | +
               |  2  | 3  |3 /
               |     |    | /
               +----------+/
            */
            _childrenNodes[0] = OctreeNode(Vector(_min[0], halfHeight, _min[2]), _center + Vector(0, halfHeight, 0) , _depth + 1);
            _childrenNodes[1] = OctreeNode(Vector(halfWidth, halfHeight, _min[2]), _center + Vector(0, 0, halfWidth), _depth + 1);
            _childrenNodes[2] = OctreeNode(_min, _center, _depth + 1);
            _childrenNodes[3] = OctreeNode(Vector(halfWidth, _min[1], _min[2]), _center + Vector(0, 0, halfWidth), _depth + 1);
            _childrenNodes[4] = OctreeNode(_center - Vector(halfWidth, 0, 0), _center + Vector(0, halfHeight, halfDepth), _depth + 1);
            _childrenNodes[5] = OctreeNode(_center, _max, _depth + 1);
            _childrenNodes[6] = OctreeNode(_center - Vector(halfWidth, halfHeight, 0), _center + Vector(0, 0, halfDepth), _depth + 1);
            _childrenNodes[7] = OctreeNode(_center - Vector(0, halfHeight, 0), _max - Vector(0, halfHeight, 0), _depth + 1);

            //Now that the node is subdivided and not a leaf anymore, we can reinsert the volumes in the node which will actually insert the volumes in the subnodes
            for(const Triangle* triangle : _triangles)
                this->insertTriangle(triangle, maxChildren, maxDepth);

            //Inserting the volume this function was initially called for
            this->insertTriangle(triangle, maxChildren, maxDepth);

            //This node is not a leaf anymore, we're clearing the volumes
            _triangles.clear();
        }
        else
            _triangles.push_back(triangle);
    }
    else
    {
        //We're going to place the triangle into the octant which contains its centroid
        Vector nodeCentroid = _center;
        Vector triangleCentroid = triangle->Centroid();

        int octant = 0;//The index of the octant in which we're going to place the triangle

        if(triangleCentroid[2] > nodeCentroid[2]) octant += 4;
        if(triangleCentroid[1] > nodeCentroid[1]) octant += 2;
        if(triangleCentroid[0] > nodeCentroid[0]) octant += 1;

        _childrenNodes[octant].insertTriangle(triangle, maxChildren, maxDepth);
    }
}

void OctreeNode::insertTriangles(const std::vector<Triangle*>& triangles, int maxChildren, int maxDepth)
{
    for(const Triangle* triangle : triangles)
        this->insertTriangle(triangle, maxChildren, maxDepth);
}

//TODO remove
//bool OctreeNode::intersect(const Ray& ray, double& t) const
//{
//    double closestT = INFINITY;
//    bool intersectionFound = false;

//    if(_isLeaf)
//    {
//        for(const Triangle* triangle : _triangles)
//        {
//            double u, v;

//            intersectionFound = triangle->Intersect(ray, u, v, t);
//            if(intersectionFound)
//                closestT = std::min(closestT, t);
//        }
//    }
//    else
//    {
//        std::priority_queue<OctreeNode, std::vector<OctreeNode>, std::greater<double>> queue;
//    }

//    return intersectionFound;
//}

Octree::Octree() {};

Octree::Octree(const Vector min, const Vector max) : _min(min), _max(max), _center((min + max) / 2) {};

void Octree::construct(const std::vector<Triangle*>& triangles, int maxChildren, int maxDepth)
{
    _root = new OctreeNode(_min, _max, 0);
    _root->insertTriangles(triangles, maxChildren, maxDepth);
}

void Octree::computeVolumes()
{
    _root->computeVolume();
}

bool Octree::intersect(const Ray& ray, double& t) const
{
    //TODO remove tFar from BoundingVolume.intersect if trash is never used
    double trash;

    double tInterTriangle = INFINITY;//Keeps the closest intersection of the ray with a triangle of the octree
    bool intersectionFound = false;

    //If we're not even intersecting the bounding volumes
    //of the whole scene, we're not going to intersect anything
    if(!_root->GetBoundingVolume().intersect(ray, t, trash))
        return false;
    else
    {
        std::priority_queue<Octree::OctreeQueueNode> queue;
        Octree::OctreeQueueNode _rootQueueNode(_root, t);
        queue.push(_rootQueueNode);

        while(!queue.empty() && tInterTriangle > queue.top()._nodeInterDistance)
        {
            OctreeNode* node = queue.top()._node;
            queue.pop();

            if(node->IsLeaf())
            {
                for(const Triangle* triangle : node->Triangles())
                {
                    if(triangle->Intersect(ray, trash, trash, t))
                    {
                        tInterTriangle = std::min(tInterTriangle, t);
                        intersectionFound = true;
                    }
                }
            }
            else
            {
                OctreeNode* childrenNodes = node->ChildrenNodes();
                for(int i = 0; i < 8; i++)
                {
                    double nodeIntersectionNear;
                    double nodeIntersectionFar;
                    OctreeNode* node = &childrenNodes[i];
                    if(node->GetBoundingVolume().intersect(ray, nodeIntersectionNear, nodeIntersectionFar))
                    {
                        double nodeInterT;

                        //This means that the ray is inside the volume and that the near plane is behind
                        if(nodeIntersectionNear < 0 && nodeIntersectionFar >= 0)
                            nodeInterT = nodeIntersectionFar;
                        else
                            nodeInterT = nodeIntersectionNear;
                        Octree::OctreeQueueNode queueNode(node, nodeInterT);
                        queue.push(queueNode);
                    }
                }
            }

        }
    }

    if(intersectionFound)
    {
        t = tInterTriangle;

        return true;
    }
    else
        return false;
}

BVH::BVH(const std::vector<Triangle*>& triangles, int maxLeafChildren, int maxDepth)
{
    //The minimum and maximum point of the AABB of the
    //root volume, i.e., the points constituting the AABB that can
    //encompass all the triangles
    Vector min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY);

//    //We're going to pack triangles by groups of 8 so that they are not alone in their volume.
//    //If we had 1 bounding volume for each triangle, this would be desastrous in terms of performance
//    for(int i = 0; i < triangles.size() / 8; i += 8)
//    {
//        std::vector<Triangle> trianglePack;
//        trianglePack.reserve(8);

//        for(int j = 0; j < 8; j++)
//        {
//            const Triangle& triangle = triangles.at(i + j);
//            trianglePack.push_back(triangle);

//            //Updating the min and the max of the set of triangles
//            for(int vertex = 0; vertex < 3; vertex++)
//            {
//                for(int component = 0; component < 3; component++)
//                {
//                    min[component] = std::min(min[component], triangle[vertex][component]);
//                    max[component] = std::max(max[component], triangle[vertex][component]);
//                }
//            }
//        }

//        BoundingVolume boundingVolume(trianglePack);
//        _boundingVolumes.push_back(boundingVolume);
//    }

//    //Processing the rest of the triangles that wasn't divisible by the pack size
//    std::vector<Triangle> trianglePack;
//    for(int i = triangles.size() - triangles.size() % 8; i < triangles.size(); i++)
//    {
//        const Triangle& triangle = triangles.at(i);
//        trianglePack.push_back(triangle);

//        //Updating the min and the max of the set of triangles
//        for(int vertex = 0; vertex < 3; vertex++)
//        {
//            for(int component = 0; component < 3; component++)
//            {
//                min[component] = std::min(min[component], triangle[vertex][component]);
//                max[component] = std::max(max[component], triangle[vertex][component]);
//            }
//        }
//    }

    for(const Triangle* triangle : triangles)
    {
        //Updating the min and the max of the set of triangles
        for(int vertex = 0; vertex < 3; vertex++)
        {
            for(int component = 0; component < 3; component++)
            {
                min[component] = std::min(min[component], (*triangle)[vertex][component]);
                max[component] = std::max(max[component], (*triangle)[vertex][component]);
            }
        }

        BoundingVolume boundingVolume(*triangle);
        _boundingVolumes.push_back(boundingVolume);
    }

    _octree = Octree(min, max);
    //Inserts all the triangles in the octree
    _octree.construct(triangles, maxLeafChildren, maxDepth);
    //Computes the volumes of all the nodes of the octree
    _octree.computeVolumes();
}

bool BVH::intersect(const Ray& ray, double& t) const
{
    return _octree.intersect(ray, t);
}

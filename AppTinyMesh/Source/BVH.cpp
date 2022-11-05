#include "BVH.h"

#include <queue>

unsigned long long int _boundingVolumesInterTests;//<! How many bounding volumes were tested against a ray for intersection
unsigned long long int _triangleInterTests;//<! How many triangles were tested against a ray for intersection
unsigned long long int _triangleEffectiveInter;//<! How many triangles were actually intersected

BoundingVolume::BoundingVolume()
{
    for(int i = 0; i < 7; i++)
        _dNears[i] = INFINITY;
    for(int i = 0; i < 7; i++)
        _dFars[i] = -INFINITY;
}

BoundingVolume::BoundingVolume(const Triangle& triangle)
{
    computeDNearDFar(&triangle);
}

BoundingVolume::BoundingVolume(const std::vector<Triangle>& triangles)
{
    computeDNearDFar(triangles);
}

BoundingVolume::BoundingVolume(const std::vector<Vector>& vertexes, const std::vector<int>& indices)
{
    std::vector<Triangle> triangles;
    triangles.reserve(indices.size() / 3);

    for(int vertexIndex  = 0; vertexIndex < indices.size(); vertexIndex += 3)
        triangles.push_back(Triangle(vertexes.at(indices.at(vertexIndex)),
                                          vertexes.at(indices.at(vertexIndex + 1)),
                                          vertexes.at(indices.at(vertexIndex + 2))));

    computeDNearDFar(triangles);
}

void BoundingVolume::extend(const BoundingVolume& volume)
{
    for(int i = 0; i < 7; i++)
    {
        _dNears[i] = std::min(_dNears[i], volume._dNears[i]);
        _dFars[i] = std::max(_dFars[i], volume._dFars[i]);
    }
}

void BoundingVolume::computeDNearDFar(const Triangle* triangle)
{
    for(int triangleVertex = 0; triangleVertex < 3; triangleVertex++)
    {
        Vector vertex = (*triangle)[triangleVertex];

        //Computing dNears and dFars
        for(int plane = 0; plane < 7; plane++)
        {
            double value = BoundingVolume::planes[plane].Normal() * vertex;

            _dNears[plane] = std::min(_dNears[plane], value);
            _dFars[plane] = std::max(_dFars[plane], value);
        }
    }
}

void BoundingVolume::computeDNearDFar(const std::vector<Triangle>& triangles)
{
    for(const Triangle& triangle : triangles)
    {
        for(int triangleVertex = 0; triangleVertex < 3; triangleVertex++)
        {
            Vector vertex = triangle[triangleVertex];

            //Computing dNears and dFars
            for(int plane = 0; plane < 7; plane++)
            {
                _dNears[plane] = BoundingVolume::planes[plane].Normal() * vertex;
                _dFars[plane] = BoundingVolume::planes[plane].Normal() * vertex;
            }
        }
    }
}

//TODO remove
//void BoundingVolume::GetMinMax(Vector& min, Vector& max) const
//{
//    min[0] =
//    return Vector((dNears[0] + dFars[0]) / 2, (dNears[1] + dFars[1]) / 2, (dNears[2] + dFars[2]) / 2);

//    for(const Triangle* triangle : _triangles)
//    {
//        for(int triangleVertex = 0; triangleVertex < 3; triangleVertex++)
//        {
//            Vector vertex = (*triangle)[triangleVertex];

//            for(int component = 0; component < 3; component++)
//            {
//                //Updating the minimum and maximum point of the AABB of this bounding volume
//                min[triangleVertex] = std::min(min[triangleVertex], vertex[component]);
//                max[triangleVertex] = std::max(max[triangleVertex], vertex[component]);
//            }
//        }
//    }
//}

//TODO remove tFar if never used
bool BoundingVolume::intersect(const Ray& ray, double& tNear, double& tFar) const
{
    double outNear = -INFINITY;
    double outFar = INFINITY;

    for(int i = 0; i < 7; i++)
    {
        //TODO precompute ND and NO
        double ND = BoundingVolume::planes[i].Normal() * ray.Direction();
        double NO = BoundingVolume::planes[i].Normal() * ray.Origin();

        double near = (_dNears[i] - NO) / ND;
        double far = (_dFars[i] - NO) / ND;

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

    return true;
}

Vector BoundingVolume::Centroid() const
{
    //The average of the nears and fars of the axis aligned planes
    return Vector((_dNears[0] + _dFars[0]) / 2, (_dNears[1] + _dFars[1]) / 2, (_dNears[2] + _dFars[2]) / 2);
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

            double halfWidth = (_max[0] - _min[0]) / 2;
            double halfHeight = (_max[1] - _min[1]) / 2;
            double halfDepth = (_max[2] - _min[2]) / 2;

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
            _childrenNodes[1] = OctreeNode(Vector(halfWidth, halfHeight, _min[2]), _center + Vector(halfWidth, halfHeight, 0), _depth + 1);
            _childrenNodes[2] = OctreeNode(_min, _center, _depth + 1);
            _childrenNodes[3] = OctreeNode(Vector(halfWidth, _min[1], _min[2]), _center + Vector(halfWidth, 0, 0), _depth + 1);
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

OctreeNode* Octree::Root() const
{
    return _root;
}

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

    double closestIntersection = INFINITY;//Keeps the closest intersection of the ray with a triangle of the octree
    bool intersectionFound = false;

    _boundingVolumesInterTests++;
    if(!_root->GetBoundingVolume().intersect(ray, trash, closestIntersection))
    {
        //If we're not even intersecting the bounding volumes
        //of the whole scene, we're not going to intersect anything
        return false;
    }
    else
    {
        std::priority_queue<Octree::OctreeQueueNode> queue;
        Octree::OctreeQueueNode _rootQueueNode(_root, -1);
        queue.push(_rootQueueNode);

        //As long as there is something to intersect and that
        //the closest intersection we've found so far is not
        //in front of the next volume to test (this would mean
        //that we're going to intersect a volume that is behind
        //the closest intersection we've found so far. We're not
        //going to find a closer intersection than we already
        //have in this case)
        while(!queue.empty() && closestIntersection > queue.top()._nodeInterDistance)
        {
            OctreeNode* node = queue.top()._node;
            queue.pop();

            if(node->IsLeaf())
            {
                for(const Triangle* triangle : node->Triangles())
                {
                    double triangleInterDist = INFINITY;
                    double u, v;

                    _triangleInterTests++;
                    if(triangle->Intersect(ray, triangleInterDist, u, v) && triangleInterDist > 0)
                    {
                        _triangleEffectiveInter++;
                        closestIntersection = std::min(closestIntersection, triangleInterDist);
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

                    _boundingVolumesInterTests++;
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
        t = closestIntersection;

        return true;
    }
    else
        return false;
}

BVH::BVH()
{
    _boundingVolumesInterTests = 0;
    _triangleInterTests = 0;
    _triangleEffectiveInter = 0;
};

BVH::BVH(const std::vector<Triangle*>& triangles, int maxLeafChildren, int maxDepth) : BVH()
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

Octree BVH::GetOctree() const
{
    return _octree;
}

bool BVH::intersect(const Ray& ray, double& t) const
{
    bool returned = _octree.intersect(ray, t);

    return returned;
}

#include <cassert>

void BVH::BVHTests()
{
    //TODO decomment les tests
//    //Test on a simple triangle
//    Triangle triangle(Vector(0, 0, 0), Vector(1, 0, 0), Vector(0.5, 0.5, 0));

//    std::vector<Triangle*> triangles;
//    triangles.push_back(&triangle);

//    BVH bvh(triangles, 8, 2);

//    double t;
//    //Ray pointing towards the center of the triangle
//    assert(bvh.intersect(Ray(Vector(0.5, 0.5/3, 1), Vector(0, 0, -1)), t));
//    assert(t == 1.0);

//    //One ray to each vertex of the triangle
//    //Also one ray to each vertex of the triangle while being inside the volume, directly on the vertex (auto-intersection)
//    //Skipping the last one as the ray/triangle intersection is the problem for this one
//    for(int i = 0; i < 2; i++)
//    {
//        Vector vertex = triangle[i];

//        assert(bvh.intersect(Ray(vertex - Vector(0, 0, 5), Vector(0, 0, 1)), t));
//        assert(t == 5.0);

//        assert(bvh.intersect(Ray(vertex, Vector(0, 0, -1)), t));
//        assert(t == 0.0);
//    }

//    //Ray pointing away from the triangle
//    assert(!bvh.intersect(Ray(Vector(0.5, 0.5/3, 1), Vector(0, 0, 1)), t));
//    //Ray parallel to the triangle
//    assert(!bvh.intersect(Ray(Vector(0.5, 0.5/3, 1), Vector(1, 0, 0)), t));
//    //Ray that intersects the bounding volume of the triangle but not the triangle itself
//    assert(!bvh.intersect(Ray(Vector(0.5, 0.5/3, 0.001), Vector(1, 0, 0)), t));

//    BVH bvhSphere(*Mesh(Icosphere(1, 2)).GetTriangles(), 8, 1000);
//    assert(bvhSphere.intersect(Ray(Vector(0, 0, -1), Vector(0, 0, 1)), t));

    //Testing the repartition of the triangles in the octants of the octree
    std::vector<Triangle*> trianglesRepartition;
    std::vector<Triangle> baseTrianglesRepartition;

    const int NB_TRIANGLES_REPARTITION = 1;

    //Creating triangles close to each other, they will be in the
    //octant 0 and they will be used to create more triangles that
    //will be in the other octants
    for(int i = 0; i < NB_TRIANGLES_REPARTITION; i++)
        baseTrianglesRepartition.push_back(Triangle(Vector(-3, 2, i * 0.05 + 2), Vector(-2, 2, i * 0.05 + 2), Vector(-2.5, 2.5, i * 0.05 + 2)));

    for(int octant = 0; octant < 8; octant++)
    {
        double offsetX = (octant & 1) / 1 * 6;
        double offsetY = (octant & 2) / 2 * -4;
        double offsetZ = (octant & 4) / 4 * -4;
        for(int i = 0; i < NB_TRIANGLES_REPARTITION; i++)
        {
            Triangle* baseTriangleOffset = new Triangle(baseTrianglesRepartition.at(i));
            baseTriangleOffset->Translate(Vector(offsetX, offsetY, offsetZ));

            trianglesRepartition.push_back(baseTriangleOffset);
        }
    }

    //TODO destructeur de la BVH
    BVH bvhRepartition(trianglesRepartition, 1, 10);

    //We're now going to make sure that the triangles have correctly been inserted in the octants
    Octree repartitionOctree = bvhRepartition.GetOctree();
    OctreeNode* octants = repartitionOctree.Root()->ChildrenNodes();
    for(int octant = 0; octant < 8; octant++)
    {
        OctreeNode octreeNode = octants[octant];
        const std::vector<const Triangle*> octantTriangles = octreeNode.Triangles();
        for(const Triangle* triangle : octantTriangles)
        {
            bool triangleFound = false;
            for(int j = octant * NB_TRIANGLES_REPARTITION; j < (octant + 1) * NB_TRIANGLES_REPARTITION; j++)
            {
                bool allVertexEqual = true;

                for(int vertex = 0; vertex < 3; vertex++)
                {
                    if((*triangle)[vertex] == (*trianglesRepartition.at(j))[vertex])
                    {
                        allVertexEqual = false;

                        break;
                    }
                }

                if(allVertexEqual)
                {
                    triangleFound = true;

                    break;
                }
            }

            assert(triangleFound);
        }
    }
}

void BVH::GetInterStats(unsigned long long int& boundingVolumesTests, unsigned long long int& triangleInterTests, unsigned long long int& triangleEffectiveInters)
{
    boundingVolumesTests = _boundingVolumesInterTests;
    triangleInterTests = _triangleInterTests;
    triangleEffectiveInters = _triangleEffectiveInter;
}

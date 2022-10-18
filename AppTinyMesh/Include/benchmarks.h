#ifndef BENCHMARKS_H
#define BENCHMARKS_H

class Benchmarks
{
public:
    void BenchmarkIcosphere(int iterations, int maxSubdivisions);
    void BenchmarkTorus(int iterations, int maxSubdivisions, int ringCount, int ringSubdiv, bool ringCountFixed);
    void BenchmarkCapsule(int iterations, int maxSubdivisions);
    void BenchmarkCylinder(int iterations, int maxSubdivisions);
};

#endif // BENCHMARKS_H

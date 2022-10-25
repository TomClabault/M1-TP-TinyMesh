#include "benchmarks.h"
#include "simpleMeshes.h"

#include <chrono>
#include <iostream>
#include <fstream>

void Benchmarks::BenchmarkIcosphere(int iterations, int maxSubdivisions)
{
    std::ofstream outputDurationFile;
    std::ofstream outputTriangleFile;

    outputDurationFile.open("benchmarkDurationIcosphere.dat");
    outputTriangleFile.open("benchmarkTriangleIcosphere.dat");

    for(int subdiv = 0; subdiv < maxSubdivisions; subdiv++)
    {
        long long int bestTimeForSubdiv = std::numeric_limits<long long int>::max();
        Icosphere icosphereCreated(1, 0);

        for(int i = 0; i < iterations; i++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            icosphereCreated = Icosphere(1, subdiv);
            auto stop = std::chrono::high_resolution_clock::now();

            long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
            if(duration < bestTimeForSubdiv)
                bestTimeForSubdiv = duration;
        }

        std::cout << "Icosphere[" << subdiv << "] --- " << bestTimeForSubdiv << "ms [" << 1.0 / bestTimeForSubdiv << "FPS]" << std::endl;

        outputDurationFile << subdiv << " " << bestTimeForSubdiv << std::endl;
        outputTriangleFile << subdiv << " " << icosphereCreated.IndicesCount() / 3.0 / 1000.0 << std::endl;
    }

    std::cout << std::endl;
}

void Benchmarks::BenchmarkTorus(int iterations, int maxSubdivisions, int ringCount, int ringSubdiv, bool ringCountFixed)
{
    std::ofstream outputDurationFile;
    std::ofstream outputTriangleFile;

    if(ringCountFixed)
    {
        outputDurationFile.open("benchmarkDurationTorus1.dat");
        outputTriangleFile.open("benchmarkTriangleTorus1.dat");
    }
    else
    {
        outputDurationFile.open("benchmarkDurationTorus2.dat");
        outputTriangleFile.open("benchmarkTriangleTorus2.dat");
    }

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        Torus torusCreated(1, 1);

        long long int bestTimeForSubdiv = std::numeric_limits<long long int>::max();

        for(int i = 0; i < iterations; i++)
        {
            if(ringCountFixed)
            {
                auto start = std::chrono::high_resolution_clock::now();
                torusCreated = Torus(1, 2, ringCount, subdiv);
                auto stop = std::chrono::high_resolution_clock::now();

                long long int duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

                if(duration < bestTimeForSubdiv)
                    bestTimeForSubdiv = duration;
            }
            else
            {
                auto start = std::chrono::high_resolution_clock::now();
                torusCreated = Torus(1, 2, subdiv, ringSubdiv);
                auto stop = std::chrono::high_resolution_clock::now();

                long long int duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

                if(duration < bestTimeForSubdiv)
                    bestTimeForSubdiv = duration;
            }
        }

        std::cout << "Torus[" << (ringCountFixed ? ringCount : subdiv) << ", " << (ringCountFixed ? subdiv : ringSubdiv) << "] --- " << bestTimeForSubdiv << "microseconds [" << 1.0 / bestTimeForSubdiv << "FPS]" << std::endl;

        outputDurationFile << subdiv << " " << bestTimeForSubdiv << std::endl;
        outputTriangleFile << subdiv << " " << torusCreated.IndicesCount() / 3 << std::endl;
    }

    outputDurationFile.close();
    outputTriangleFile.close();

    std::cout << std::endl;
    std::cout << std::endl;
}

void Benchmarks::BenchmarkCapsule(int iterations, int maxSubdivisions)
{
    std::ofstream outputDurationFile;
    std::ofstream outputTriangleFile;

    outputDurationFile.open("benchmarkDurationCapsule.dat");
    outputTriangleFile.open("benchmarkTriangleCapsule.dat");

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        Capsule capsuleCreated = Capsule(10, 10, 10, 10, 10);
        long long int bestTimeForSubdiv = std::numeric_limits<long long int>::max();

        for(int i = 0; i < iterations; i++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            capsuleCreated = Capsule(1, 2, 2, subdiv, 30);
            auto stop = std::chrono::high_resolution_clock::now();

            long long int duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

            if(duration < bestTimeForSubdiv)
                bestTimeForSubdiv = duration;
        }

        std::cout << "Capsule[" << subdiv << "] --- " << bestTimeForSubdiv << "microseconds [" << 1.0 / bestTimeForSubdiv << "FPS]" << std::endl;

        outputDurationFile << subdiv << " " << bestTimeForSubdiv << std::endl;
        outputTriangleFile << subdiv << " " << capsuleCreated.IndicesCount() / 3.0 / 1000.0 << std::endl;
    }

    outputDurationFile.close();
    outputTriangleFile.close();

    std::cout << std::endl;
    std::cout << std::endl;
}

void Benchmarks::BenchmarkCylinder(int iterations, int maxSubdivisions)
{
    std::ofstream outputDurationFile;
    std::ofstream outputTriangleFile;

    outputDurationFile.open("benchmarkDurationCylinder.dat");
    outputTriangleFile.open("benchmarkTriangleCylinder.dat");

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        Cylinder cylinderCreated(1, 1, 2, 5);
        long long int bestTimeForSubdiv = std::numeric_limits<long long int>::max();

        for(int i = 0; i < iterations; i++)
        {

            auto start = std::chrono::high_resolution_clock::now();
            cylinderCreated = Cylinder(1, 2, 2, subdiv);
            auto stop = std::chrono::high_resolution_clock::now();

            long long int duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

            if(duration < bestTimeForSubdiv)
                bestTimeForSubdiv = duration;
        }

        std::cout << "Cylinder[" << subdiv << "] --- " << bestTimeForSubdiv << "microseconds [" << 1.0 / bestTimeForSubdiv << "FPS]" << std::endl;

        outputDurationFile << subdiv << " " << bestTimeForSubdiv << std::endl;
        outputTriangleFile << subdiv << " " << cylinderCreated.IndicesCount() / 3 << std::endl;
    }

    outputDurationFile.close();
    outputTriangleFile.close();

    std::cout << std::endl;
    std::cout << std::endl;
}

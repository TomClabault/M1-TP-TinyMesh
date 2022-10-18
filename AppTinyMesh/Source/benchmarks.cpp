#include "realtime.h"
#include "benchmarks.h"

#include <iostream>
#include <fstream>

void Benchmarks::BenchmarkIcosphere(int iterations, int maxSubdivisions)
{
    RenderingProfiler profiler;
    profiler.Update();

    for(int subdiv = 0; subdiv < maxSubdivisions; subdiv++)
    {
        double totalForSubdiv = 0;
        for(int i = 0; i < iterations; i++)
        {
            profiler.InitCPUOnly();
            Icosphere(1, subdiv);
            profiler.Update();

            totalForSubdiv += profiler.msPerFrame;
        }

        std::cout << "Icosphere[" << subdiv << "] --- " << totalForSubdiv / (double)iterations << "ms [" << profiler.framePerSecond << "FPS]" << std::endl;
    }

    std::cout << std::endl;
}

void Benchmarks::BenchmarkTorus(int iterations, int maxSubdivisions, int ringCount, int ringSubdiv, bool ringCountFixed)
{
    RenderingProfiler profiler;
    std::ofstream outputFile;
    if(ringCountFixed)
        outputFile.open("benchmarkTorus1.dat");
    else
        outputFile.open("benchmarkTorus2.dat");

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        //Torus torusCreated = Torus(1, 1, 1, 1);
        double totalForSubdiv = 0;
        for(int i = 0; i < iterations; i++)
        {
            if(ringCountFixed)
            {
                profiler.InitCPUOnly();
                Torus(1, 2, ringCount, subdiv);
            }
            else
            {
                profiler.InitCPUOnly();
                Torus(1, 2, subdiv, ringSubdiv);
            }
            profiler.Update();

            totalForSubdiv += profiler.msPerFrame;
        }

        //std::cout << torusCreated.getMemorySize() / 1000.0 << "KO" << std::endl;
        double time = totalForSubdiv / (double)iterations;

        std::cout << "Torus[" << (ringCountFixed ? ringCount : subdiv) << ", " << (ringCountFixed ? subdiv : ringSubdiv) << "] --- " << time << "ms [" << profiler.framePerSecond << "FPS]" << std::endl;
        outputFile << subdiv << " " << time << std::endl;
    }

    outputFile.close();
    std::cout << std::endl;
    std::cout << std::endl;
}

void Benchmarks::BenchmarkCapsule(int iterations, int maxSubdivisions)
{
    RenderingProfiler profiler;
    std::ofstream outputFile;
    outputFile.open("benchmarkCapsule.dat");

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        double totalForSubdiv = 0;
        Capsule capsuleCreated = Capsule(10, 10, 10, 10, 10);
        for(int i = 0; i < iterations; i++)
        {
            profiler.InitCPUOnly();
            capsuleCreated = Capsule(1, 2, 2, subdiv, 30);
            profiler.Update();

            totalForSubdiv += profiler.msPerFrame;
        }

        double time = totalForSubdiv / (double)iterations;

        //std::cout << capsuleCreated.getMemorySize() /1000.0 << "KO" << std::endl;

        std::cout << "Capsule[" << subdiv << "] --- " << time << "ms [" << profiler.framePerSecond << "FPS]" << std::endl;
        outputFile << subdiv << " " << time << std::endl;
    }

    outputFile.close();
    std::cout << std::endl;
    std::cout << std::endl;
}

void Benchmarks::BenchmarkCylinder(int iterations, int maxSubdivisions)
{
    RenderingProfiler profiler;
    std::ofstream outputFile;
    outputFile.open("benchmarkCylinder.dat");

    for(int subdiv = 3; subdiv < maxSubdivisions; subdiv++)
    {
        double totalForSubdiv = 0;
        for(int i = 0; i < iterations; i++)
        {
            profiler.InitCPUOnly();
            Cylinder(1, 2, 2, subdiv);
            profiler.Update();

            totalForSubdiv += profiler.msPerFrame;
        }

        double time = totalForSubdiv / (double)iterations;

        std::cout << "Cylinder[" << subdiv << "] --- " << time << "ms [" << profiler.framePerSecond << "FPS]" << std::endl;
        outputFile << subdiv << " " << time << std::endl;
    }

    outputFile.close();
    std::cout << std::endl;
    std::cout << std::endl;
}

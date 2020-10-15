#ifndef PARTICLES_CUH
#define PARTICLES_CUH

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../hpp/MT.hpp"
#include "../hpp/conf.hpp"

namespace PhysPeach{
    struct Particles{
        double *diam;
        double packing;
        double *x;
        double *mem;
        double *v;
        double *f;
    };
    double powerParticles(Particles* p);
    void createParticles(Particles*, int);
    void createParticles(Particles*, std::ifstream*);
    void deleteParticles(Particles*);
    void squeezePositions(Particles*, double);
    bool convergedFire(Particles*);
}

#endif
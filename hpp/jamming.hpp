#ifndef JAMMING_HPP
#define JAMMING_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include "../hpp/conf.hpp"
#include "../hpp/particles.hpp"
#include "../hpp/cells.hpp"

namespace PhysPeach{
    struct Jamming {
        double phi;
        double t;
        Particles p;
        Cells cells;
        Lists lists;
    };
    double L(Jamming*);
    void createJamming(Jamming*);
    void deleteJamming(Jamming*);
    int fireJamming(Jamming*);
    void addDphi(Jamming*, double);
    void findJamming(Jamming*);
    void squeezeJamming(Jamming*);
}

#endif
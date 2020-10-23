#include <iostream>

#include "../hpp/conf.hpp"
#include "../hpp/jamming.hpp"

int ID;
int Np;
double Phi_init;
double Dphi;
int main(int argc, char** argv) {
    ID = atoi(argv[1]);
    Np = atoi(argv[2]);
    Phi_init = atof(argv[3]);
    Dphi = atof(argv[4]);

    std::cout << "-- hello jamming --" << std::endl;
    std::cout << "ID       : " << ID << std::endl;
    std::cout << "Np       : " << Np << std::endl;
    std::cout << "Phi_init : " << Phi_init << std::endl;
    std::cout << "Dphi     : " << Dphi << std::endl;
    std::cout << "-------------------" << std::endl << std::endl;

    PhysPeach::Jamming jam;
    PhysPeach::loadSwapMC(&jam);
    PhysPeach::findJamming(&jam);
    PhysPeach::squeezeJamming(&jam);
    PhysPeach::deleteJamming(&jam);

    return 0;
}
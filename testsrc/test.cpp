#include <iostream>
#include "../hpp/MT.hpp"

#include "../testhpp/conf_test.hpp"
#include "../testhpp/particles_test.hpp"

int ID;
int Np;
double Phi_init;

int main(){
    ID = 0;
    Np = 1024;
    Phi_init = 0.3;

    std::cout << "--start test--" << std::endl;
    init_genrand(1);

    //conf_test
    PhysPeach::powIntTest();
    PhysPeach::setZeroTest();

    //particles_test
    PhysPeach::createParticlesTest();

    std::cout << "---finished---" << std::endl;
    return 0;
}
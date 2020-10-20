#include "../testhpp/jamming_test.hpp"

namespace PhysPeach{
    void createJammingTest(){
        Jamming jam;
        createJamming(&jam);
        //std::cout << jam.lists.Nl << std::endl;
        deleteJamming(&jam);
        return;
    }

    void fireJammingTest(){
        Jamming jam;
        createJamming(&jam);
        fireJamming(&jam);

        double fsum = 0.;
        double f2;
        for(int par1 = 0; par1 < Np; par1++){
            f2 = 0.;
            for(int d = 0; d < D; d++){
                f2 += jam.p.f[par1+d*Np] * jam.p.f[par1+d*Np];
            }
            fsum += sqrt(f2);
        }
        assert(fsum < 1.0e-13 * Np);
        deleteJamming(&jam);
        return;
    }
}
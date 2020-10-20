#include "../testhpp/jamming_test.hpp"

namespace PhysPeach{
    void createJammingTest(){
        Jamming jam;
        createJamming(&jam);
        //std::cout << jam.lists.Nl << std::endl;
        deleteJamming(&jam);
        return;
    }
}
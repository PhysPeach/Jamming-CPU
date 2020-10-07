#include <iostream>

#include "../testhpp/conf_test.hpp"

int main(){
    std::cout << "--start test--" << std::endl;

    //conf_test
    PhysPeach::powIntTest();
    PhysPeach::setZeroTest();

    std::cout << "---finished---" << std::endl;
    return 0;
}
#include <iostream>
#include <iomanip>
#include <sstream>

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

    std::cout << "-- p-phi jamming --" << std::endl;
    std::cout << "ID       : " << ID << std::endl;
    std::cout << "Np       : " << Np << std::endl;
    std::cout << "Phi_init : " << Phi_init << std::endl;
    std::cout << "Dphi     : " << Dphi << std::endl;
    std::cout << "---------------------" << std::endl << std::endl;

    PhysPeach::Jamming jam;
    PhysPeach::loadJamming(&jam);

    std::ofstream file;

    std::cout << "Squeeze from phi = " << jam.phi << std::endl;
    double jammingPoint = jam.phi;
    double dphi;
    std::ostringstream PphiName;
    PphiName << "../p-phi/p-phi_N" << Np << "_Phi" << Phi_init << "_Dphi" << Dphi << "_id" << ID <<".data";
    file.open(PphiName.str().c_str());
    file << jammingPoint << std::endl << std::endl;

    double delta = 0.;
    dphi = 1.0e-6;
    for (int i = 0; i < 10; i++){
        if(delta > Dphi){
            break;
        }
        addDphi(&jam, dphi);
        delta = jam.phi - jammingPoint;
        file << delta << " " << P(&jam.p, L(&jam), &jam.lists) << std::endl;
    }
    dphi = 1.0e-5;
    for (int i = 1; i < 10; i++){
        if(delta > Dphi){
            break;
        }
        addDphi(&jam, dphi);
        delta = jam.phi - jammingPoint;
        file << delta << " " << P(&jam.p, L(&jam), &jam.lists) << std::endl;
    }
    dphi = 1.0e-4;
    double OutputAt = jam.phi - jammingPoint;
    while(delta <= Dphi){
        addDphi(&jam, dphi);
        delta = jam.phi - jammingPoint;
        if(delta >= OutputAt){
            file << delta << " " << P(&jam.p, L(&jam), &jam.lists) << std::endl;
            OutputAt = 1.1 * delta;
        }
    }
    file << delta << " " << P(&jam.p, L(&jam), &jam.lists) << std::endl;
    file.close();
    std::cout << "finished!: phi = " << jam.phi << std::endl;

    std::ostringstream squeezeName;
    squeezeName << "../squeeze/sq_N" << Np << "_Phi" << Phi_init << "_Dphi" << Dphi << "_id" << ID <<".data";
    file.open(squeezeName.str().c_str());
    file << std::setprecision(15);
    for(int par1 = 0; par1 < Np; par1++){
        file << jam.p.diam[par1] << " ";
        for(int d = 0; d < D; d++){
            file << jam.p.x[d*Np + par1] << " ";
        }
        file << std::endl;
    }
    file.close();

    PhysPeach::deleteJamming(&jam);

    return 0;
}
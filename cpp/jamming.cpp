#include "../hpp/jamming.hpp"

namespace PhysPeach{
    double L(Jamming* jam){
        return pow(jam->p.packing/jam->phi, 1./(double)D);
    }
    void createJamming(Jamming* jam){
        jam->t = 0.;
        jam->phi = Phi_init;
        createParticles(&jam->p);
        createCells(&jam->cells, L(jam));
        createLists(&jam->lists, &jam->cells);
        updateCellList(&jam->cells, &jam->lists, L(jam), jam->p.x);
        return;
    }
    void deleteJamming(Jamming* jam){
        deleteParticles(&jam->p);
        deleteLists(&jam->lists);
        deleteCells(&jam->cells);
        return;
    }

    int fireJamming(Jamming* jam){
        int loop = 0;
        double dt = dt_init;
        double alpha = alpha_init;
        double power;
        bool mustUpdateCellList;

        bool converged = false;
        int cp = 0;
        setZero(jam->p.v, D*Np);
        while(!converged){
            loop++;
            mustUpdateCellList = updateParticles(&jam->p, L(jam), dt, &jam->lists);
            if(mustUpdateCellList){
                updateCellList(&jam->cells, &jam->lists, L(jam), jam->p.x);
            }
            converged = convergedFire(&jam->p);
            power = powerParticles(&jam->p);
            modifyVelocities(&jam->p, alpha);
            if(power < 0){
                setZero(jam->p.v, D*Np);
                alpha = alpha_init;
                dt *= 0.5;
                cp = 0;
            }else{
                cp++;
                if(cp > 5){
                    dt *= 1.1;
                    if(dt > dt_max){
                        dt = dt_max;
                    }
                    alpha *= 0.99;
                    cp = 0;
                }
            }
            jam->t += dt;
            if(loop == 1000000){
                std::cout << "dt: " << dt << std::endl;
            }
        }
        return loop;
    }

    void addDphi(Jamming* jam, double dphi){
        static double Lstart = L(jam);
        bool mustUpdateCellList;
        int loop;

        jam->phi += dphi;
        double Lend = L(jam);
        squeezePositions(&jam->p, Lend/Lstart);
        Lstart = Lend;

        mustUpdateCellList = updateMem(&jam->p, Lend);
        if((jam->cells.numOfCellsPerSide > 3 && Lend/(double)jam->cells.numOfCellsPerSide < 2 * a_max) || 2 * a_max < Lend/(double)(1+jam->cells.numOfCellsPerSide)){
            deleteCells(&jam->cells);
            createCells(&jam->cells, Lend);
            mustUpdateCellList = true;
        }
        if(mustUpdateCellList){
            updateCellList(&jam->cells, &jam->lists, Lend, jam->p.x);
        }
        loop = fireJamming(jam);
        std::cout << "    " << jam->phi << ", " << U(&jam->p, Lend, &jam->lists) << ", " << P(&jam->p, Lend, &jam->lists) << ", " << loop << std::endl;
        return;
    }

    void findJamming(Jamming* jam){
        std::cout << std::endl << "Find jamming point" << std::endl;

        std::cout << "    Squeeze from phi = " << jam->phi << std::endl;
        std::cout << "    dphi, E, P, loop:" << std::endl;
        fireJamming(jam);
        while (P(&jam->p, L(jam), &jam->lists) < 1.0e-8){
            addDphi(jam, 1.0e-4);
        }
        std::cout << "    Expand from phi = " << jam->phi << std::endl;
        while (P(&jam->p, L(jam), &jam->lists) > 1.0e-8){
            addDphi(jam, -1.0e-5);
        }
        while (P(&jam->p, L(jam), &jam->lists) < 1.0e-8){
            addDphi(jam, 1.0e-6);
        }
        std::cout << "-> Jamming Point: " << jam->phi << std::endl;
        return;
    }
}
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
        bool mustUpdateCells;

        bool converged = false;
        int cp = 0;
        setZero(jam->p.v, D*Np);
        while(!converged){
            loop++;
            mustUpdateCells = updateParticles(&jam->p, L(jam), dt, &jam->lists);
            if(mustUpdateCells){
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
}
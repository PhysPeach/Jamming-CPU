#include "../hpp/particles.hpp"

namespace PhysPeach{
    void createParticles(Particles *p, int numOfCellsPerSide){
        p->packing = 0;

        p->diam = (double*)malloc(Np * sizeof(double));
        p->x = (double*)malloc(D*Np * sizeof(double));
        p->mem = (double*)malloc(D * Np * sizeof(double));
        p->v = (double*)malloc(D * Np * sizeof(double));
        p->f = (double*)malloc(D * Np * sizeof(double));

        for (int par1 = 0; par1 < Np; par1++){
            p->diam[par1] = sqrt(A * a_min * a_min / (A - 2 * a_min * a_min * genrand_real1()));
            p->packing += powInt(p->diam[par1], D);
        }
        p->packing *= 0.5 * pi / D;

        int numOfCells = powInt(numOfCellsPerSide, D);
        int NpC = Np / numOfCells;
        double L = pow(p->packing/Phi_init, 1./(double)D);
        double Lc = L / (double)numOfCellsPerSide;

        for(int i1 = 0; i1 < numOfCellsPerSide; i1++){
            for(int i2 = 0; i2 < numOfCellsPerSide; i2++){
                if(D == 2){
                    for(int m = 0; m < NpC; m++){
                        p->x[(i1*numOfCellsPerSide+i2)*NpC+m] = (i1 + genrand_real1()) * Lc - 0.5 * L;
                        p->x[(i1*numOfCellsPerSide+i2)*NpC+m+Np] = (i2 + genrand_real1()) * Lc - 0.5 * L;
                    }
                }
                else if(D == 3){
                    for(int k = 0; k < numOfCellsPerSide; k++){
                        for(int m = 0; m < NpC; m++){
                            p->x[((i1*numOfCellsPerSide+i2)*numOfCellsPerSide + k)*NpC+m] = (i1 + genrand_real1()) * Lc - 0.5 * L;
                            p->x[((i1*numOfCellsPerSide+i2)*numOfCellsPerSide + k)*NpC+m+Np] = (i2 + genrand_real1()) * Lc - 0.5 * L;
                            p->x[((i1*numOfCellsPerSide+i2)*numOfCellsPerSide + k)*NpC+m+Np] = (k + genrand_real1()) * Lc - 0.5 * L;
                        }
                    }
                }
                else{
                    std::cout << "dimention err" << std::endl;
                    exit(1);
                }
            }
        }
        for(int m = NpC*powInt(numOfCellsPerSide, D); m < Np; m++){
            p->x[m] = (genrand_real1() - 0.5 )* L;
            p->x[m+Np] = (genrand_real1() - 0.5) * L;
            if(D == 3){
                p->x[m+2*Np] = (genrand_real1() - 0.5) * L;
            }
        }
        memcpy(p->mem, p->x, D * Np * sizeof(double));

        for (int par1 = 0; par1 < D*Np; par1++){
            p->v[par1] = 0.;
            p->f[par1] = 0.;
        }
        return;
    }

    void createParticles(Particles *p, std::ifstream *in){
        p->packing = 0;

        p->diam = (double*)malloc(Np * sizeof(double));
        p->x = (double*)malloc(D*Np * sizeof(double));
        p->mem = (double*)malloc(D * Np * sizeof(double));
        p->v = (double*)malloc(D * Np * sizeof(double));
        p->f = (double*)malloc(D * Np * sizeof(double));

        for (int par1 = 0; par1 < Np; par1++){
            *in >> p->diam[par1];
            for(int d = 0; d < D; d++){
                *in >> p->x[par1 + d*Np];
                p->v[par1 + d*Np] = 0;
                p->f[par1 + d*Np] = 0;
            }
            p->packing += powInt(p->diam[par1], D);
        }
        p->packing *= 0.5 * pi / D;

        memcpy(p->mem, p->x, D * Np * sizeof(double));
        return;
    }

    void deleteParticles(Particles *p){
        free(p->diam);
        free(p->x);
        free(p->mem);
        free(p->v);
        free(p->f);
        return;
    }
}
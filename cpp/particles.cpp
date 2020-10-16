#include "../hpp/particles.hpp"

namespace PhysPeach{
    double powerParticles(Particles* p){
        double power = 0.;
        for(int par1 = 0; par1 < D*Np; par1++){
            power += p->v[par1] * p->f[par1];
        }
        return power;
    }
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

    bool updateMem(Particles* p, double L){
        double Lh = 0.5 * L;
        double dx[D];
        double frag = 0.25 * a_max * a_max;
        //fix mem_x by gamma
        /*for(int par1 = 0; par1 < Np; par1++){
            p->mem[par1] += (p->gamma - p->gammaMem) * L;
        }
        p->gammaMem = p->gamma;*/
        for(int par1 = 0 ; par1 < Np; par1++){
            for(int d = 0; d < D; d++){
                dx[d] = p->x[d*Np+par1] - p->mem[d*Np+par1];
            }
            for(int d = D-1; d >= 0; d--){
                if(dx[d] > Lh){dx[d] -= L;}
                if(dx[d] < -Lh){dx[d] += L;}
                if(dx[d] > frag){
                    memcpy(p->mem, p->x, D * Np * sizeof(double));
                    return true;
                }
            }
        }
        return false;
    }

    void squeezePositions(Particles *p, double a){
            for(int par1 = 0; par1 < D * Np; par1++){
            p->x[par1] *= a;
            p->mem[par1] *= a;
        }
        return;
    }

    void modifyVelocities(Particles* p, double a){
        double v, f;
        for(int par1 = 0; par1 < Np; par1++){
            v = 0.;
            f = 0.;
            for(int d = 0; d < D; d++){
                v += p->v[d*Np+par1] * p->v[d*Np+par1];
                f += p->f[d*Np+par1] * p->f[d*Np+par1];
            }
            if(f > 0.){
                v = sqrt(v);
                f = sqrt(f);
                for(int d = 0; d < D; d++){
                    p->v[d*Np+par1] = (1-a) * p->v[d*Np+par1] + a * v * p->f[d*Np+par1] / f;
                }
            }else{
                v = sqrt(v);
                for(int d = 0; d < D; d++){
                    p->v[d*Np+par1] = (1-a) * p->v[d*Np+par1];
                }
            }
        }
        return;
    }

    bool convergedFire(Particles *p){
        double fsum = 0.;
        double fmax = 1.0e-13;
        double f2;
        for(int par1 = 0; par1 < Np; par1++){
            f2 = 0.;
            for(int d = 0; d < D; d++){
                f2 += p->f[par1+d*Np] * p->f[par1+d*Np];
            }
            fsum += sqrt(f2);
        }
        if(fsum > fmax * Np){
            return false;
        }
        return true;
    }
}
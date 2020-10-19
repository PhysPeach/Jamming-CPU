#include "../hpp/cells.hpp"

namespace PhysPeach{
    void createCells(Cells *cells, double L){
        cells->numOfCellsPerSide = (int)(L/(2. * a_max));
        if(cells->numOfCellsPerSide < 3){
            cells->numOfCellsPerSide = 3;
        }
        double buf = 1.5;
        cells->Nc = (int)(buf * (double)Np/ (double)powInt(cells->numOfCellsPerSide, D));

        int NoC = powInt(cells->Nc, D)*cells->Nc;
        cells->cell = (int*)malloc(NoC*sizeof(int));
        setZero(cells->cell, NoC);
        return;
    }
    
    void deleteCells(Cells *cells){
        free(cells->cell);
        return;
    }
    
    void increaseNc(Cells *cells){
        cells->Nc++;
        int NoC = powInt(cells->Nc, D)*cells->Nc;
        free(cells->cell);
        cells->cell = (int*)malloc(NoC*sizeof(int));
        setZero(cells->cell, NoC);
        return;
    }

    bool putParticlesIntoCells(Cells *cells, double L, double* x){
        int i, j, k;
        double Lc = L/(double)cells->numOfCellsPerSide;
        int NoC = powInt(cells->Nc, D)*cells->Nc;
        setZero(cells->cell, NoC);
        for(int n = 0; n < Np; n++){
            i = (x[n] + 0.5 * L)/Lc;
            j = (x[n+Np] + 0.5 * L)/Lc;
            if(i >= cells->numOfCellsPerSide){i -= cells->numOfCellsPerSide;}
            if(i < 0){i += cells->numOfCellsPerSide;}
            if(j >= cells->numOfCellsPerSide){j -= cells->numOfCellsPerSide;}
            if(j < 0){j += cells->numOfCellsPerSide;}
            if (D == 2){
                cells->cell[(i*cells->numOfCellsPerSide+j)*cells->Nc]++;
                if(cells->cell[(i*cells->numOfCellsPerSide+j)*cells->Nc] >= cells->Nc){
                    return false;
                }
                cells->cell[(i*cells->numOfCellsPerSide+j)*cells->Nc + cells->cell[(i*cells->numOfCellsPerSide+j)*cells->Nc]] = n;
            }else{
                std::cout << "dimention err" << std::endl;
                    exit(1);
            }
        }
        return true;
    }

    void updateCells(Cells *cells, double L, double* x){
        bool success = false;
        success = putParticlesIntoCells(cells, L, x);
        while (!success){
            std::cout << cells->Nc << std::endl;
            increaseNc(cells);
            success = putParticlesIntoCells(cells, L, x);
        }

        return;
    }
}
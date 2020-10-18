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
}
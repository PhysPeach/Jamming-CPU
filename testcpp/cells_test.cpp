#include "../testhpp/cells_test.hpp"

namespace PhysPeach{
    void createCellsTest(){

        Cells cells;
        createCells(&cells, 3.);
        assert(cells.numOfCellsPerSide == 3);
        assert(cells.Nc == (int)(1.5 * (double)Np/ (double)powInt(3, D)));
        deleteCells(&cells);

        createCells(&cells, 40.);
        assert(cells.numOfCellsPerSide == 12);
        assert(cells.Nc == (int)(1.5 * (double)Np/ (double)powInt(12, D)));
        deleteCells(&cells);

        return;
    }

    void increaseNcTest(){
        Cells cells;
        createCells(&cells, 40.);

        assert(cells.Nc == (int)(1.5 * (double)Np/ (double)powInt(12, D)));
        increaseNc(&cells);
        assert(cells.Nc == 1 + (int)(1.5 * (double)Np/ (double)powInt(12, D)));
        deleteCells(&cells);

        return;
    }

    void updateCellsTest(){
        //test it in small Np
        double *x;
        x = (double*)malloc(D*Np*sizeof(double));
        for(int par1 = 0; par1 < D*Np; par1++){
            x[par1] = 1.;
        }

        Cells cells;
        createCells(&cells, 10.);
        updateCells(&cells, 10., x);
        assert(cells.Nc == Np+1);
        deleteCells(&cells);

        free(x);
        return;
    }
}
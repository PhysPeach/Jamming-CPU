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
}
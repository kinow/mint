#include "mntGridMover.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>

/**
 * Create uniform VTK grid objects
 * @param nx number of x cells
 * @param ny number of y cells
 * @param lonMin low corner longitude
 * @param lonMax high corner longitude
 * @param latMin low corner latitude
 * @param latMax high corner latitude
 * @param coords vtkDoubleArray object (will be built)
 */
void createUniformGrid(int nx, int ny, double lonMin, double lonMax, double latMin, double latMax,
                       double coords[]) {

    double dx = (lonMax - lonMin) / (double) nx;
    double dy = (latMax - latMin) / (double) ny;
    const int numPointsPerCell = 4;
    const int numDims = 3;
    const int offset = numPointsPerCell * numDims;

    // cell counter
    int k = 0;

    for (int i = 0; i < nx; ++i) {

        double x0 = lonMin + i*dx;
        double x1 = x0 + dx;

        for (int j = 0; j < ny; ++j) {

            double y0 = latMin + j*dy;
            double y1 = y0 + dy;

            // vertices are in counterclockwise order
            coords[k*offset +  0] = x0; coords[k*offset +  1] = y0; coords[k*offset +  2] = 0.0;
            coords[k*offset +  3] = x1; coords[k*offset +  4] = y0; coords[k*offset +  5] = 0.0;
            coords[k*offset +  6] = x1; coords[k*offset +  7] = y1; coords[k*offset +  8] = 0.0;
            coords[k*offset +  9] = x0; coords[k*offset + 10] = y1; coords[k*offset + 11] = 0.0;

            k++;
        }
    }
}


int main() {

    int nx = 10;
    int ny = 5;

    const int fixLonAcrossDateline = 0;
    const int averageLonAtPole = 0;
    const int degrees = 1;
    const int numVertsPerCell = 4;
    const int numDims = 3;

    int ier;

    Grid_t* grid0 = NULL;
    Grid_t* grid1 = NULL;

    ier = mnt_grid_new(&grid0);
    assert(ier == 0);
    ier = mnt_grid_new(&grid1);
    assert(ier == 0);
    vtkIdType numCells = nx * ny;
    std::vector<double> coords0(numCells*numVertsPerCell*numDims);
    createUniformGrid(nx, ny, 0., 360., -90., 90., &coords0[0]);
    // copy
    std::vector<double> coords1 = coords0;
    ier = mnt_grid_setPointsPtr(&grid0, numVertsPerCell, numCells, &coords0[0]);
    assert(ier == 0);
    ier = mnt_grid_setPointsPtr(&grid1, numVertsPerCell, numCells, &coords1[0]);
    assert(ier == 0);

    // flags for lon-lat grid
    // ier = mnt_grid_setFlags(&srcGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);
    // ier = mnt_grid_setFlags(&dstGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);

    // move grid1 along the velocity field
    double dt = 0.1;
    GridMover_t* mover = NULL;
    ier = mnt_gridmover_new(&mover);
    assert(ier == 0);


    // clean up
    ier = mnt_gridmover_del(&mover);
    assert(ier == 0);

}

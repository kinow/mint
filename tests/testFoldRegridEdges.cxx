#include "mntRegridEdges.h"
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

    int srcNx = 10;
    int srcNy = 2;
    int dstNx = 5;
    int dstNy = 3;

    const int fixLonAcrossDateline = 0;
    const int averageLonAtPole = 0;
    const int degrees = 1;
    const int numVertsPerCell = 4;
    const int numDims = 3;

    int ier;

    Grid_t* srcGridObj = NULL;
    Grid_t* dstGridObj = NULL;
    RegridEdges_t* rg = NULL;

    ier = mnt_grid_new(&srcGridObj);
    assert(ier == 0);
    vtkIdType srcNumCells = srcNx * srcNy;
    std::vector<double> srcCoords(srcNumCells*numVertsPerCell*numDims);
    createUniformGrid(srcNx, srcNy, 0., 360., -90., 90., &srcCoords[0]);
    ier = mnt_grid_setPointsPtr(&srcGridObj, numVertsPerCell, srcNumCells, &srcCoords[0]);
    assert(ier == 0);

    // flags for lon-lat grid
    // ier = mnt_grid_setFlags(&srcGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);
    // ier = mnt_grid_setFlags(&dstGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);

    ier = mnt_grid_new(&dstGridObj);
    assert(ier == 0);
    vtkIdType dstNumCells = dstNx * dstNy;
    std::vector<double> dstCoords(dstNumCells*numVertsPerCell*numDims);
    createUniformGrid(dstNx, dstNy, 0., 360., -100., 100., &dstCoords[0]);
    ier = mnt_grid_setPointsPtr(&dstGridObj, numVertsPerCell, dstNumCells, &dstCoords[0]);
    assert(ier == 0);

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    // set the src/dst grids
    rg->srcGridObj = srcGridObj;
    rg->dstGridObj = dstGridObj;

    // compute the regridding weights

    int numCellsPerBucket = 128;
    const double periodX = 360.;
    int debug = 3;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket, periodX, debug);
    assert(ier == 0);

    ier = mnt_grid_dump(&srcGridObj, "foldRegridEdges_srcGrid.vtk");
    assert(ier == 0);

    ier = mnt_grid_dump(&dstGridObj, "foldRegridEdges_dstGrid.vtk");
    assert(ier == 0);

    // clean up, the regridder takes care of freeing the grids
    // ier = mnt_grid_del(&srcGridObj);
    // assert(ier == 0);
    // ier = mnt_grid_del(&dstGridObj);
    // assert(ier == 0);
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

}

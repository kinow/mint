#include "mntGridMover.h"
#include "mntDots.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>

/**
 * Create uniform VTK grid objects with a velocity field attached
 * @param nx number of x cells
 * @param ny number of y cells
 * @param lonMin low corner longitude
 * @param lonMax high corner longitude
 * @param latMin low corner latitude
 * @param latMax high corner latitude
 * @param coords vtkDoubleArray object (will be built)
 * @param velocity vtkDoubleArray object (will be built)
 */
void createUniformGridWithVelocity(int nx, int ny, 
                       double lonMin, double lonMax, double latMin, double latMax,
                       double coords[], double velocity[]) {

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

            // velocity field
            velocity[k*offset +  0] = 0.0; velocity[k*offset +  1] = cos(M_PI*x0/180.); velocity[k*offset +  2] = 0.0;
            velocity[k*offset +  3] = 0.0; velocity[k*offset +  4] = cos(M_PI*x1/180.); velocity[k*offset +  5] = 0.0;
            velocity[k*offset +  6] = 0.0; velocity[k*offset +  7] = cos(M_PI*x1/180.); velocity[k*offset +  8] = 0.0;
            velocity[k*offset +  9] = 0.0; velocity[k*offset + 10] = cos(M_PI*x0/180.); velocity[k*offset + 11] = 0.0;

            k++;
        }
    }
}

/** 
 * Compute the standard deviation distance between two grids
 * @param grid0 first grid
 * @param grid1 second grid
 * @return standard deviation
 */
double computeStandardDeviation(Grid_t* grid0, Grid_t* grid1) {

    int ier;
    double res = 0;
    Vec3 xyz0, xyz1;

    vtkUnstructuredGrid* ugrid0 = NULL;
    vtkUnstructuredGrid* ugrid1 = NULL;

    ier = mnt_grid_get(&grid0, &ugrid0);
    assert(ier == 0);
    ier = mnt_grid_get(&grid1, &ugrid1);
    assert(ier == 0);

    vtkPoints* points0 = ugrid0->GetPoints();
    vtkPoints* points1 = ugrid1->GetPoints();

    vtkIdType numPoints = points0->GetNumberOfPoints();
    for (vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

        points0->GetPoint(pointId, &xyz0[0]);
        points1->GetPoint(pointId, &xyz1[0]);

        Vec3 diff = xyz1 - xyz0;
        res += dot(diff, diff);
    }
    res /= (double) numPoints;
    res = sqrt(res);


    return res;
}

void testInterpVelocity(int nx, int ny) {

    const int fixLonAcrossDateline = 0;
    const int averageLonAtPole = 0;
    const int degrees = 1;
    const int numVertsPerCell = 4;
    const int numDims = 3;

    int ier;

    Grid_t* grid = NULL;
    ier = mnt_grid_new(&grid);

    assert(ier == 0);
    vtkIdType numCells = nx * ny;
    std::vector<double> coords(numCells*numVertsPerCell*numDims);
    std::vector<double> velocity(numCells*numVertsPerCell*numDims);
    createUniformGridWithVelocity(nx, ny, 0., 360., -90., 90., &coords[0], &velocity[0]);

    ier = mnt_grid_setPointsPtr(&grid, numVertsPerCell, numCells, &coords[0]);
    assert(ier == 0);

    // flags for lon-lat grid
    // ier = mnt_grid_setFlags(&srcGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);
    // ier = mnt_grid_setFlags(&dstGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);

    // move grid1 along the velocity field
    GridMover_t* mover = NULL;
    ier = mnt_gridmover_new(&mover);
    assert(ier == 0);

    ier = mnt_gridmover_setGrid(&mover, grid);
    assert(ier == 0);

    int numCellsPerBucket = 128;
    double periodX = 360.;
    ier = mnt_gridmover_build(&mover, numCellsPerBucket, periodX);
    assert(ier == 0);

    ier = mnt_gridmover_setPointVelocityPtr(&mover, numDims, &velocity[0]);
    assert(ier == 0);

    {
        const double targetXyz[] = {180., -90.05, 0.};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }
    {
        const double targetXyz[] = {2.1600000000000e+02, -9.0040450849719e+01, 0.0000000000000e+00};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }
    {
        const double targetXyz[] = {2.5200000000000e+02, -9.0015450849719e+01, 0.0000000000000e+00};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }
    {
        const double targetXyz[] = {2.8800000000000e+02, 9.0015450849719e+01, 0.0000000000000e+00};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }
    {
        const double targetXyz[] = {3.2400000000000e+02, 9.0040450849719e+01, 0.0000000000000e+00};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }
    {
        const double targetXyz[] = {3.6000000000000e+02, 9.0050000000000e+01, 0.0000000000000e+00};
        double vel[3];
        ier = mnt_gridmover_interpVelocity(&mover, targetXyz, vel);
        assert(ier == 0);
        std::cout << "at point " << targetXyz[0] << ',' << targetXyz[1] 
                  << " the velocity is " << vel[0] << ',' << vel[1] << '\n';
    }

    // clean up
    ier = mnt_grid_del(&grid);
    ier = mnt_gridmover_del(&mover);
    assert(ier == 0);

}

void testAdvance(int nx, int ny) {

    const int numVertsPerCell = 4;
    const int numDims = 3;

    const int fixLonAcrossDateline = 0;
    const int averageLonAtPole = 0;
    const int degrees = 1;

    int ier;

    Grid_t* grid = NULL;
    ier = mnt_grid_new(&grid);
    
    assert(ier == 0);
    vtkIdType numCells = nx * ny;
    std::vector<double> coords(numCells*numVertsPerCell*numDims);
    std::vector<double> velocity(numCells*numVertsPerCell*numDims);
    createUniformGridWithVelocity(nx, ny, 0., 360., -90., 90., &coords[0], &velocity[0]);

    ier = mnt_grid_setPointsPtr(&grid, numVertsPerCell, numCells, &coords[0]);
    assert(ier == 0);

    // flags for lon-lat grid
    // ier = mnt_grid_setFlags(&srcGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);
    // ier = mnt_grid_setFlags(&dstGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);

    // move grid1 along the velocity field
    GridMover_t* mover = NULL;
    ier = mnt_gridmover_new(&mover);
    assert(ier == 0);

    ier = mnt_gridmover_setGrid(&mover, grid);
    assert(ier == 0);

    int numCellsPerBucket = 128;
    double periodX = 360.;
    ier = mnt_gridmover_build(&mover, numCellsPerBucket, periodX);
    assert(ier == 0);

    ier = mnt_gridmover_setPointVelocityPtr(&mover, numDims, &velocity[0]);
    assert(ier == 0);

    ier = mnt_gridmover_advance(&mover, 0.1);
    assert(ier == 0);

    // clean up
    ier = mnt_grid_del(&grid);
    ier = mnt_gridmover_del(&mover);
    assert(ier == 0);

}


void testRectilinearGrid2(int nx, int ny) {

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
    std::vector<double> velocity(numCells*numVertsPerCell*numDims);
    createUniformGridWithVelocity(nx, ny, 0., 360., -90., 90., &coords0[0], &velocity[0]);

    // copy
    std::vector<double> coords1 = coords0;

    ier = mnt_grid_setPointsPtr(&grid0, numVertsPerCell, numCells, &coords0[0]);
    assert(ier == 0);

    ier = mnt_grid_setPointsPtr(&grid1, numVertsPerCell, numCells, &coords1[0]);
    assert(ier == 0);

    ier = mnt_grid_dump(&grid0, "testGridMover_grid0.vtk");

    // flags for lon-lat grid
    // ier = mnt_grid_setFlags(&srcGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);
    // ier = mnt_grid_setFlags(&dstGridObj, fixLonAcrossDateline, averageLonAtPole, degrees);

    // move grid1 along the velocity field
    GridMover_t* mover = NULL;
    ier = mnt_gridmover_new(&mover);
    assert(ier == 0);

    ier = mnt_gridmover_setGrid(&mover, grid1);
    assert(ier == 0);

    int numCellsPerBucket = 128;
    double periodX = 360.;
    ier = mnt_gridmover_build(&mover, numCellsPerBucket, periodX);
    assert(ier == 0);

    // DEBUG
    {
        const double point[] = {180., -90.05, 0};
        const double tol = 1.e-12;
        double pcoords[3];
        double weights[12];
        vtkIdType cellId = mover->loc->findCellMultiValued(point, tol, pcoords, weights);
        std::cerr << "can we find the point? point = " << Vec3{point} << " cellId = " << cellId << '\n';
        assert(cellId >= 0);
    }


    ier = mnt_gridmover_setPointVelocityPtr(&mover, numDims, &velocity[0]);
    assert(ier == 0);

    double deltaTime = 50.0;
    ier = mnt_gridmover_advance(&mover, deltaTime);
    assert(ier == 0);

    // // go back
    // ier = mnt_gridmover_advance(&mover, -deltaTime);
    // assert(ier == 0);

    // save grid0 and grid1
    ier = mnt_grid_dump(&grid0, "testGridMover_grid0.vtk");
    assert(ier == 0);
    ier = mnt_grid_dump(&grid1, "testGridMover_grid1.vtk");
    assert(ier == 0);

    double diffGrids = computeStandardDeviation(grid0, grid1);
    std::cout << "diff between grids: " << diffGrids << '\n';

    // clean up
    ier = mnt_grid_del(&grid0);
    assert(ier == 0);
    ier = mnt_grid_del(&grid1);
    assert(ier == 0);
    ier = mnt_gridmover_del(&mover);
    assert(ier == 0);

}


int main() {

    int nx = 10;
    int ny = 5;

    testInterpVelocity(nx, ny);
    testAdvance(nx, ny);
    testRectilinearGrid2(nx, ny);
}

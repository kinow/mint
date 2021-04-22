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
 * @param grid vtkUnstructured grid object (will be built)
 * @param points vtkPoints object (will be built)
 * @param coords vtkDoubleArray object (will be built)
 */
void createUniformGrid(int nx, int ny, double lonMin, double lonMax, double latMin, double latMax,
                       vtkUnstructuredGrid* grid, vtkPoints* points, vtkDoubleArray* coords) {

    double point[3];

    int numCells = nx * ny;
    int numPoints = 4 * numCells;
    double dx = (lonMax - lonMin) / (double) nx;
    double dy = (latMax - latMin) / (double) ny;

    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(numPoints);

    int k = 0;
    for (int i = 0; i < nx; ++i) {
        double x0 = lonMin + i*dx;
        double x1 = x0 + dx;
        for (int j = 0; j < ny; ++j) {
            double y0 = latMin + j*dy;
            double y1 = y0 + dy;
            point[0] = x0; point[1] = y0; point[2] = 0.0;
            coords->SetTuple(k*4 + 0, point);
            point[0] = x1; point[1] = y0; point[2] = 0.0;
            coords->SetTuple(k*4 + 1, point);
            point[0] = x1; point[1] = y1; point[2] = 0.0;
            coords->SetTuple(k*4 + 2, point);
            point[0] = x0; point[1] = y1; point[2] = 0.0;
            coords->SetTuple(k*4 + 3, point);
            k++;
        }
    }

    points->SetData(coords);

    grid->SetPoints(points);
    grid->Allocate(numCells, 1);
    const vtkIdType nptsPerCell = 4;
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(nptsPerCell);
    for (size_t iCell = 0; iCell < numCells; ++iCell) {
        for (size_t i = 0; i < nptsPerCell; ++i) {
            ptIds->SetId(i, nptsPerCell*iCell + i);
        }
        grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    ptIds->Delete();

}


int main() {

    // build the src/dst grids
    vtkUnstructuredGrid* srcGrid = vtkUnstructuredGrid::New();
    vtkPoints* srcPoints = vtkPoints::New();
    vtkDoubleArray* srcArray = vtkDoubleArray::New();

    createUniformGrid(10, 5, 0., 360., -90., 90., srcGrid, srcPoints, srcArray);

    vtkUnstructuredGrid* dstGrid = vtkUnstructuredGrid::New();
    vtkPoints* dstPoints = vtkPoints::New();
    vtkDoubleArray* dstArray = vtkDoubleArray::New();

    createUniformGrid(12, 6, 0., 360., -100., 100., dstGrid, dstPoints, dstArray);

    RegridEdges_t* rg;
    int ier;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    // set the src/dst grids
    rg->srcGrid = srcGrid;
    rg->dstGrid = dstGrid;

    int num_cells_per_bucket = 1;
    double periodX = 360.;
    int debug = 2;
    ier = mnt_regridedges_build(&rg, num_cells_per_bucket, periodX, debug);
    assert(ier == 0);

    ier = mnt_regridedges_print(&rg);
    assert(ier == 0);
    
    std::string outputFilename = "foldRegridEdgesWeights.nc";
    ier = mnt_regridedges_dumpWeights(&rg, outputFilename.c_str(), outputFilename.size());
    assert(ier == 0);

    // clean up
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

    srcGrid->Delete();
    srcPoints->Delete();
    srcArray->Delete();

    dstGrid->Delete();
    dstPoints->Delete();
    dstArray->Delete();

}

#include "mntUgridReader.h"

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>

#ifndef MNT_VTK_UGRID_ADAPTOR
#define MNT_VTK_UGRID_ADAPTOR

class VtkUgridAdaptor {

public:

/**
 * Constructor
 */
VtkUgridAdaptor(const UgridReader& ur);

/**
 * Destructor
 */
~VtkUgridAdaptor() {
    this->vGrid->Delete();
    this->vPoints->Delete();
    this->vPointCoords->Delete();
}

/**
 * Write grid to VTK file
 * @param filename file name
 */
void write(const std::string& filename);


private:

    vtkDoubleArray* vPointCoords;

    vtkPoints* vPoints;

    vtkUnstructuredGrid* vGrid;

};

#endif // MNT_UGRID_READER

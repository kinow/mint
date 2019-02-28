#include <mntVtkUgridAdaptor.h>
#include <vtkCell.h>
#include <vtkUnstructuredGridWriter.h>

VtkUgridAdaptor::VtkUgridAdaptor(const UgridReader& ur) {

    this->vPointCoords = vtkDoubleArray::New();
    this->vPoints = vtkPoints::New();
    this->vGrid = vtkUnstructuredGrid::New();

    this->vPointCoords->SetNumberOfComponents(3);
    // quads
    size_t ncells = ur.getNumberOfFaces();
    const size_t numPointsPerQuad = 4;
    this->vPointCoords->SetNumberOfTuples(4 * ncells);

    // set the coordinates
    size_t pointId = 0;
    for (vtkIdType faceId = 0; faceId < ncells; ++faceId) {
        for (auto point : ur.getFacePointsRegularized(faceId)) {
            this->vPointCoords->SetTuple(pointId, &point[0]);
            pointId++;
        }
    }

    // build the points object
    this->vPoints->SetData(this->vPointCoords);

    // build the grid
    this->vGrid->Allocate(ncells, 1);
    this->vGrid->SetPoints(this->vPoints);
    for (vtkIdType faceId = 0; faceId < ncells; ++faceId) {
        vtkIdType k = 4 * faceId;
        vtkIdType pIdList[] = {k, k + 1, k + 2, k + 3};
        this->vGrid->InsertNextCell(VTK_QUAD, 4, pIdList);
    }

}

void 
VtkUgridAdaptor::write(const std::string& filename) {
    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(this->vGrid);
    writer->Update();
    writer->Delete();    
}
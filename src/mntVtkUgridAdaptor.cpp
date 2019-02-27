#include <mntVtkUgridAdaptor.h>
#include <vtkCell.h>
#include <vtkUnstructuredGridWriter.h>

VtkUgridAdaptor::VtkUgridAdaptor(const UgridReader& ur) {

    this->vPointCoords = vtkDoubleArray::New();
    this->vPoints = vtkPoints::New();
    this->vGrid = vtkUnstructuredGrid::New();

    this->vPointCoords->SetNumberOfComponents(3);
    this->vPointCoords->SetNumberOfTuples(ur.getNumberOfPoints());

    // set the coordinates
    for (vtkIdType i = 0; i < ur.getNumberOfPoints(); ++i) {
        this->vPointCoords->SetTuple(i, ur.getPoint(i));
    }

    // build the points object
    this->vPoints->SetData(this->vPointCoords);

    // build the grid

    int cellType = VTK_QUAD;
    vtkIdList* ptIds = vtkIdList::New();

    const int numPointsPerQuad = 4;
    const int numEdgesPerQuad = 4;
    const int numPointsPerEdge = 2;

    ptIds->SetNumberOfIds(numPointsPerQuad);

    size_t ncells = ur.getNumberOfFaces();
    this->vGrid->Allocate(ncells, 1);

    for (vtkIdType iFace = 0; iFace < ncells; ++iFace) {


        // collect the point Ids of this face
        ptIds->Reset();

        const long long* edgeIds = ur.getFaceEdgeIds(iFace);
        for (size_t iEdge = 0; iEdge < numEdgesPerQuad; ++iEdge) {
            const long long* edgePointIds = ur.getEdgePointIds(edgeIds[iEdge]);
            for (size_t iPoint = 0; iPoint < numPointsPerEdge; ++iPoint) {
                ptIds->InsertUniqueId(edgePointIds[iPoint]);
            }
        }

        // ptIds contains the list of point Ids. Need to order the list
        vtkIdType ptId0 = ptIds->GetId(0); // base
        vtkIdType ptId1 = ptIds->GetId(1);
        vtkIdType ptId2 = ptIds->GetId(2);
        vtkIdType ptId3 = ptIds->GetId(3);
        const double* p0 = ur.getPoint(ptId0);
        const double* p1 = ur.getPoint(ptId1);
        const double* p2 = ur.getPoint(ptId2);
        const double* p3 = ur.getPoint(ptId3);
        Vector<double> point0(p0, p0 + 3);
        Vector<double> point1(p1, p1 + 3);
        Vector<double> point2(p2, p2 + 3);
        Vector<double> point3(p3, p3 + 3);
        point1 -= point0;
        point2 -= point0;
        point3 -= point0;
        vtkIdType pList[] = {ptId0, ptId1, ptId2, ptId3};
        if (point1[0]*point2[1] - point1[1]*point2[2] < 0.) {
            // swap ptId1 and ptId2
            vtkIdType tmpId = ptId1;
            Vector<double> tmpVec = point1;
            pList[1] = pList[2];
            pList[2] = tmpId;
            point1 = point2;
            point2 = tmpVec;
        }
        if (point2[0]*point3[1] - point2[1]*point3[2] < 0.) {
            // swap ptId2 and ptId3
            vtkIdType tmpId = ptId2;
            Vector<double> tmpVec = point2;
            pList[2] = pList[3];
            pList[3] = tmpId;
            point2 = point3;
            point3 = tmpVec;
        }
        this->vGrid->InsertNextCell(VTK_QUAD, 4, pList);

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
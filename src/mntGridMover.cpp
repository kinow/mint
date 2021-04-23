#include <mntGridMover.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <mntVecN.h>

extern "C"
int mnt_gridmover_new(GridMover_t** self) {

    *self = new GridMover_t();
    (*self)->ugrid = NULL;
    (*self)->velArray = vtkDoubleArray::New();

    return 0;
}

extern "C"
int mnt_gridmover_del(GridMover_t** self) {

    (*self)->velArray->Delete();
    delete *self;

    return 0;
}

extern "C"
int mnt_gridmover_setGrid(GridMover_t** self, Grid_t* grid) {

    (*self)->ugrid = grid->grid;

    return 0;
}

extern "C"
int mnt_gridmover_build(GridMover_t** self, int numCellsPerBucket, double periodX, int debug) {

    (*self)->loc->SetDataSet((*self)->ugrid);
    (*self)->loc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    (*self)->loc->setPeriodicityLengthX(periodX);
    (*self)->loc->enableFolding();
    (*self)->loc->BuildLocator();

    return 0;
}

extern "C"
int mnt_gridmover_setPointVelocityPtr(GridMover_t** self, int numComps, double* pointVelocity) {

    if (!(*self)->ugrid) {
        std::cerr << "ERROR: setGrid before calling setPointVelocityPtr\n";
        return 1;
    }

    (*self)->velArray->SetNumberOfComponents(numComps);
    vtkIdType numPoints = (*self)->ugrid->GetNumberOfPoints();
    (*self)->velArray->SetNumberOfTuples(numPoints);
    int save = 1; // caller should take care of disposing the array
    (*self)->velArray->SetVoidArray(pointVelocity, numPoints, save);

    return 0;
}

extern "C"
int mnt_gridmover_interpVelocity(GridMover_t** self, const double xyz[], double velocity[]) {

    const double tol2 = 1.e-12;

    // 3D
    double pcoords[3];

    // hex or quad cell
    double weights[12];

    vtkGenericCell* cell;
    Vec3 uvw;

    int ier = 0;

    // initialize to zero. This will be the value if the target falls outside the domain
    for (int j = 0; j < 3; ++j) {
        velocity[j] = 0;
    }

    // find the cell and interpolation weights
    vtkIdType cellId = (*self)->loc->FindCell(xyz, tol2, cell, pcoords, weights);

    if (cellId >= 0) {

        // found the cell that holds the target point

        // get the vertices of the cell
        vtkIdList* pointIds = cell->GetPointIds();

        for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); ++i) {

            // get the weight for that vertex
            double wght = weights[i];

            // get the point Id for that vertex
            vtkIdType pointId = pointIds->GetId(i);

            // get the velocity attached to that vertex
            (*self)->velArray->GetTuple(pointId, &uvw[0]);

            // add the interpolation contribution from that vertex
            for (int j = 0; j < 3; ++j) {
                velocity[j] += wght * uvw[j];
            }
        }
    }
    else {
        // point is outside of the domain, zero velocity
        ier = 1;
    }

    return ier;
}

extern "C"
int mnt_gridmover_advance(GridMover_t** self, double deltaTime) {

    // error codes
    int ier = 0;
    int status;

    // Runge-Kutta coefficients
    Vec3 k1, k2, k3, k4;

    // positions integrated along the trajectory
    Vec3 xyz0, xyz1, xyz2, xyz3;

    const double dtOver6 = deltaTime / 6.;
    const double dtOver2 = deltaTime / 2.;

    vtkPoints* points = (*self)->ugrid->GetPoints();
    vtkIdType numPoints = points->GetNumberOfPoints();

    for (vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

        // Runge-Kutta 4th order integration

        // compute k1
        points->GetPoint(pointId, &xyz0[0]);
        status = mnt_gridmover_interpVelocity(self, &xyz0[0], &k1[0]);
        ier += status;
        
        // compute k2
        xyz1 = xyz0 + dtOver2*k1;
        status = mnt_gridmover_interpVelocity(self, &xyz1[0], &k2[0]);
        ier += status;

        // compute k3
        xyz2 = xyz0 + dtOver2*k2;
        status = mnt_gridmover_interpVelocity(self, &xyz2[0], &k3[0]);
        ier += status;

        // compute k4
        xyz3 = xyz0 + deltaTime*k3;
        status = mnt_gridmover_interpVelocity(self, &xyz2[0], &k4[0]);
        ier += status;

        // compute and update position
        xyz0 += dtOver6*(k1 + 2.0*k2 + 2.0*k2 + k3);
        points->SetPoint(pointId, &xyz0[0]);
    }

    return ier;
}

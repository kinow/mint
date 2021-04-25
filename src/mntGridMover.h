#include <mntGrid.h>
#include <vmtCellLocator.h>
#include <vtkDoubleArray.h>

#ifndef MNT_GRID_MOVER
#define MNT_GRID_MOVER

/**
 * A class that moves a grid using a velocity field
 */

struct GridMover_t {

    /* unstructured grid */
    vtkUnstructuredGrid* ugrid;

    /* the array that holds the velocity at cell vertices */
    vtkDoubleArray* velArray;

    /** cell locator (octree-based) for fast cell search */
    vmtCellLocator* loc;
};

/**
 * Constructor
 * @param self instance of GridMover_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_new(GridMover_t** self);

/**
 * Destructor
 * @param self instance of GridMover_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_del(GridMover_t** self);

/**
 * Set the grid
 * @param self instance of GridMover_t
 * @param grid instance of Grid_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_setGrid(GridMover_t** self, Grid_t* grid);

/**
 * Build the grid mover
 * @param self instance of GridMover_t
 * @param numCellsPerBucket average number of cells per bucket
 * @param periodX periodicity length (set to 0 if non-periodic)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_build(GridMover_t** self, int numCellsPerBucket, double periodX);

/**
 * Set the pointer to the velocity field
 * @param self instance of GridMover_t
 * @param numComps number of components
 * @param pointVelocity velocity field at the grid vertices [u0, v0, w0, u1, v1, w1, ...]
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_setPointVelocityPtr(GridMover_t** self, int numComps, double* pointVelocity);


/**
 * Interpolate the velocity field at a target point
 * @param self instance of GridMover_t
 * @param xyz target point
 * @param velocity (output)
 * @return error code (0 is OK)
 * @note the velocity is set to zero outside the grid's domain
 */
extern "C"
int mnt_gridmover_interpVelocity(GridMover_t** self, const double xyz[], double velocity[]);

/**
 * Integrate the motion of the grid along the velocity field
 * @param self instance of GridMover_t
 * @param deltaTime time interval
 * @return error code (0 is OK)
 */
extern "C"
int mnt_gridmover_advance(GridMover_t** self, double deltaTime);


#endif // MNT_GRID_MOVER

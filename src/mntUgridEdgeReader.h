#include <vector>
#include <string>
#include <vtkCell.h>
#ifndef MNT_UGRID_EDGE_READER
#define MNT_UGRID_EDGE_READER

class UgridEdgeReader {

/**
 * Constructor
 */
UgridEdgeReader();

/**
 * Destructor
 */
~UgridEdgeReader();

/**
 * Get the number of points
 * @return number
 */
size_t getNumberOfPoints() const;

/**
 * Get the number of edges
 * @return number
 */
size_t getNumberOfEdges() const;

/**
 * Get point/vertex
 * @param point Id 
 * @return point coordinates
 */
const double* getPoint(size_t ptId) const;

/**
 * Get edge to point connectivity
 * @param edge Id 
 * @return point Ids
 */
const vtkIdType* getEdge(size_t edgeId) const;

/**
 * Get edge points, adding/subtracting a longitude period if necessary
 * @param edge Id
 * @param pBeg start edge coordinates
 * @param pEnd end edge coordinatess
 */
void getEdgePoints(size_t edgeId, double pBeg[], double pEnd[]) const;

/**
 * Load from Ugrid file 
 * @param filename file name
 * @return error (0=OK)
 */
int load(const std::string& filename);


private:

	int readPoints(int ncid);

	int readEdgeConnectivity(int ncid);

	int findVariableIdWithCfRole(int ncid, const std::string& cf_role, int* ndims, int dimids[]);

	int findVariableIdWithStandardName(int ncid, const std::string& standard_name, int* ndims, int dimids[]);

    // vertex raw data
    std::vector<double> verts;

    // edge to node connectivity
    std::vector<vtkIdType> edge2Nodes;
};

#endif // MNT_UGRID_EDGE_READER

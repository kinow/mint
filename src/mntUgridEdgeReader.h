#include <vector>
#include <algorithm>
#include "MvVector.h"
#include <map>
#include <string>
#include <limits>
#include <vtkCell.h>
#ifndef MNT_UGRID_EDGE_READER
#define MNT_UGRID_EDGE_READER

#define NUM_PARAM_DIMS 2
#define NUM_SPACE_DIMS 3

class UgridEdgeReader {

public:

/**
 * Constructor
 */
UgridEdgeReader();

/**
 * Destructor
 */
~UgridEdgeReader();


/**
 * Get the number of edges
 * @return number
 */
size_t getNumberOfEdges() const;


/**
 * Get edge
 * @param edge Id 
 * @param pBeg start point of the edge (output)
 * @param pEnd end point of the edge (output)
 */
void getEdge(size_t edgeId, double pBeg[], double pEnd[]) const;

/**
 * Get min/max range of the domain
 * @param xmin low point of the domain (output)
 * @param xmax high point of the domain (output)
 */
void getRange(double xmin[], double xmax[]) const;


/**
 * Load from Ugrid file 
 * @param filename file name
 * @return error (0=OK)
 */
int load(const std::string& filename);

/**
 * Get the edge vertices 
 * @return array [p0x, p0y, p0z, p1x, p1y, p1z, ....]
 */
const std::vector<double>& getEdgePoints() const {
    return this->edge2Points;
}


private:

    // edge to node connectivity
    std::vector<double> edge2Points;

    // domain min/max
    Vector<double> xmin;
    Vector<double> xmax;

    std::vector<double> readPoints(int ncid);

    int readEdgeConnectivity(int ncid, const std::vector<double>& points);

    int findVariableIdWithCfRole(int ncid, const std::string& cf_role, int* ndims, int dimids[]);

    int findVariableIdWithStandardName(int ncid, const std::string& standard_name, int* ndims, int dimids[]);

};

#endif // MNT_UGRID_EDGE_READER

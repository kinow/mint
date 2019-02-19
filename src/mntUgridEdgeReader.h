#include <vector>
#include "MvVector.h"
#include <map>
#include <string>
#include <limits>
#include <vtkCell.h>
#ifndef MNT_UGRID_EDGE_READER
#define MNT_UGRID_EDGE_READER

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
 * Build the edge locator
 * @param numEdgesPerBucket average number of edges per bucket
 */
void buildLocator(int numEdgesPerBucket);

/**
 * Get the edges that likely interesect a line
 * @param pBeg start point of the line
 * @param pEnd end point of the line
 * @return list of edge indices
 */
std::vector<size_t> getEdgesAlongLine(const double pBeg[], const double pEnd[]) const;


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


private:

    // edge to node connectivity
    std::vector<double> edge2Points;

    // domain min/max
    std::vector<double> xmin;
    std::vector<double> deltas;

    std::map<size_t, std::vector<size_t> > buckets;
    size_t nBuckets;

	std::vector<double> readPoints(int ncid);

	int readEdgeConnectivity(int ncid, const std::vector<double>& points);

	int findVariableIdWithCfRole(int ncid, const std::string& cf_role, int* ndims, int dimids[]);

	int findVariableIdWithStandardName(int ncid, const std::string& standard_name, int* ndims, int dimids[]);


	inline Vector<size_t> getBucketLoc(const double p[]) const {

    	Vector<size_t> indexLoc(2);
    	indexLoc[0] = (size_t) std::floor(this->nBuckets * (p[0] - this->xmin[0])/this->deltas[0]);
    	indexLoc[1] = (size_t) std::floor(this->nBuckets * (p[1] - this->xmin[1])/this->deltas[1]);
    	// our locator is in 2d only
    	return indexLoc;
    }

};

#endif // MNT_UGRID_EDGE_READER

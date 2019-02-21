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
 * Build the edge locator
 * @param numEdgesPerBucket average number of edges per bucket
 */
void buildLocator(int numEdgesPerBucket, double tol=1.e-2);

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
    Vector<double> xmin;
    Vector<double> deltas;

    std::map<size_t, std::vector<size_t> > buckets;
    size_t nBuckets;

    std::vector<double> readPoints(int ncid);

    int readEdgeConnectivity(int ncid, const std::vector<double>& points);

    int findVariableIdWithCfRole(int ncid, const std::string& cf_role, int* ndims, int dimids[]);

    int findVariableIdWithStandardName(int ncid, const std::string& standard_name, int* ndims, int dimids[]);

    inline Vector<double> getBucketSpaceLoc(const double p[]) const {
        Vector<double> res(p, p + NUM_SPACE_DIMS);
        res -= this->xmin;
        res /= this->deltas;
        res *= (double) this->nBuckets;
        return res;
    }

    inline Vector<int> getBucketCellLoc(const double p[]) const {
        Vector<double> loc = this->getBucketSpaceLoc(p);
        Vector<int> res(NUM_PARAM_DIMS);
        for (size_t i = 0; i < loc.size(); ++i) {
            res[i] = std::min((int) this->nBuckets - 1, std::max(0, (int) std::floor(loc[i])));
        }
        return res;
    }

};

#endif // MNT_UGRID_EDGE_READER

#include <vector>
#include <set>
#include <algorithm>

#include "MvVector.h"
#include <map>
#include <string>
#include <limits>
#include <vtkCell.h>
#ifndef MNT_EDGES_LOCATOR
#define MNT_EDGES_LOCATOR

#define NUM_PARAM_DIMS 2
#define NUM_SPACE_DIMS 3

class EdgesLocator {

public:

/**
 * Constructor
 */
EdgesLocator();

/**
 * Destructor
 */
~EdgesLocator();


/**
 * Set the domain range
 * @param xmin low end
 * @param xmax high end
 */
void setRange(const double xmin[], const double xmax[]) {
    this->xmin.resize(NUM_SPACE_DIMS);
    this->deltas.resize(NUM_SPACE_DIMS);
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        this->xmin[i] = xmin[i];
        this->deltas[i] = xmax[i] - xmin[i];
    }
}

/**
 * Build the edge locator
 * @param edge2Points array of points
 * @param numEdgesPerBucket average number of edges per bucket
 */
void build(const std::vector<double>& edge2Points, int numEdgesPerBucket);

/**
 * Get the edges that likely interesect a line
 * @param pBeg start point of the line
 * @param pEnd end point of the line
 * @return list of edge indices
 */
std::set<vtkIdType> getEdgesAlongLine(const double pBeg[], const double pEnd[]) const;


private:

    Vector<double> xmin;
    Vector<double> deltas;

    std::map<size_t, std::vector<size_t> > buckets;
    size_t nBuckets;

    std::set<size_t> getBucketsAlongLine(const double pBeg[], const double pEnd[]) const;

    inline Vector<double> getBucketSpaceLoc(const double p[]) const {
        Vector<double> res(NUM_PARAM_DIMS);
        for (size_t i = 0; i < NUM_PARAM_DIMS; ++i) {
            res[i] = ( (p[i] - this->xmin[i]) / this->deltas[i] ) * (double)(this->nBuckets);
        }
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

#endif // MNT_EDGES_LOCATOR

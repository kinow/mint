#include <vector>
#include <set>
#include <algorithm>

#include <vtkQuad.h>
#include <vtkPoints.h>

#include "MvVector.h"
#include "mntUgridReader.h"
#include <map>
#include <string>
#include <limits>

#ifndef MNT_FACES_LOCATOR
#define MNT_FACES_LOCATOR

#define NUM_PARAM_DIMS 2
#define NUM_SPACE_DIMS 3

class FacesLocator {

public:

/**
 * Constructor
 */
FacesLocator();

/**
 * Destructor
 */
~FacesLocator();


/**
 * Build the edge locator
 * @param ur instance of UgridReader
 * @param numFacesPerBucket average number of faces per bucket
 */
void build(const UgridReader& ur, int numFacesPerBucket);

/**
 * Get the edges that likely interesect a line
 * @param pBeg start point of the line
 * @param pEnd end point of the line
 * @return list of edge indices
 */
std::set<vtkIdType> getFacesAlongLine(const double pBeg[], const double pEnd[]) const;

/** 
 * Get the face that contains a point
 * @param point point
 * @param xi parametric coordinates (output)
 * @return face ID
 */
vtkIdType getFace(const double point[], double xi[]) const;


private:

    const UgridReader* ugrid;
    vtkQuad* quad;

    // domain low point and size
    Vector<double> xmin;
    Vector<double> deltas;

    // maps a bucket to a list of faces
    std::map<size_t, std::vector<vtkIdType> > buckets;

    // number of buckets
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

#endif // MNT_FACES_LOCATOR

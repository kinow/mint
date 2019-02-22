#include <mntEdgesLocator.h>
#include <netcdf.h>
#include <iostream>
#include <algorithm>
#include <cmath>

#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

EdgesLocator::EdgesLocator() {
}

/**
 * Destructor
 */
EdgesLocator::~EdgesLocator() {
}

void
EdgesLocator::build(const std::vector<double>& edge2Points, int numEdgesPerBucket, double tol) {

    size_t numEdges = edge2Points.size() / (NUM_SPACE_DIMS * 2);

    // number of cells in x and y directions
    this->nBuckets = std::max(1, (int)(std::sqrt((double)numEdges/(2.0 * (double) numEdgesPerBucket))));

    for (size_t ie = 0; ie < numEdges; ++ie) {

        // get the start/end point of the edge
        const double* pBeg = &edge2Points[0 + (0 + ie*2)*NUM_SPACE_DIMS];
        const double* pEnd = &edge2Points[0 + (1 + ie*2)*NUM_SPACE_DIMS];

        // build vectors
        Vector<double> pa(pBeg, pBeg + NUM_SPACE_DIMS);
        Vector<double> u(pEnd, pEnd + NUM_SPACE_DIMS);
        u -= pa;
        double uNormSquare = dot(u, u);

        // compute the bucket index location
        Vector<int> iBeg = this->getBucketCellLoc(pBeg);
        Vector<int> iEnd = this->getBucketCellLoc(pEnd);
        Vector<int> i0 = std::min(iBeg, iEnd);
        Vector<int> i1 = std::max(iBeg, iEnd);

        // iterate over all the cells in the box
        for (int j = i0[1]; j <= i1[1]; ++j) {
            for (int i = i0[0]; i <= i1[0]; ++i) {

                // average the vertex positions
                Vector<double> quadCentre(NUM_SPACE_DIMS);
                quadCentre[0] = i + 0.5; quadCentre[1] = j + 0.5;

                // check if line intersects with this bucket. Compute the line parameter that
                // minimizes the distance of the cell's vertices to the line. Then check if this
                // point is inside the quad
                double lam = dot(quadCentre - pa, u)/uNormSquare;
                lam = std::min(1.0, lam);
                lam = std::max(0.0, lam);
                Vector<double> quadPt = pa + lam*u;
                // is the point inside the quad?
                if (quadPt[0] >= i - tol && quadPt[0] <= i + 1 + tol &&
                    quadPt[1] >= j - tol && quadPt[1] <= j + 1 + tol) {

                    // flat index
                    size_t k = i + j * this->nBuckets;

                    std::map<size_t, std::vector<size_t> >::iterator it = this->buckets.find(k);

                    if (it != this->buckets.end()) {
                        // add all the edge Ids in this bucket
                        it->second.push_back(ie);
                    }
                    else {
                        // create new entry and add all the edge Ids
                        std::pair<size_t, std::vector<size_t> > kv(k, std::vector<size_t>(1, ie));
                        this->buckets.insert(kv);
                    }
                } // test if closest point is inside bucket
            } // j bucket iteration
        } // i bucket iteration
    } // edge iteration
}

std::vector<vtkIdType> 
EdgesLocator::getEdgesAlongLine(const double pBeg[], const double pEnd[]) const {

    Vector<int> iBeg = this->getBucketCellLoc(pBeg);
    Vector<int> iEnd = this->getBucketCellLoc(pEnd);
    Vector<int> i0 = std::min(iBeg, iEnd);
    Vector<int> i1 = std::max(iBeg, iEnd);

    std::vector<vtkIdType> edgeIds;

    for (int j = i0[1]; j <= i1[1]; ++j) {
        for (int i = i0[0]; i <= i1[0]; ++i) {

            // flat index
            size_t k = i + j * this->nBuckets;

            std::map<size_t, std::vector<size_t> >::const_iterator it = this->buckets.find(k);
            for (auto ie : it->second) {
                // add all the edges that belong to that cell
                edgeIds.push_back(ie);
            }
        }
    }

    return edgeIds;
}


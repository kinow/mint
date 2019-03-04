#include <mntFacesLocator.h>
#include <netcdf.h>
#include <iostream>
#include <algorithm>
#include <cmath>


FacesLocator::FacesLocator() {
    this->ugrid = 0;
    this->quad = vtkQuad::New();
}

FacesLocator::~FacesLocator() {
    this->quad->Delete();
}

void
FacesLocator::build(const UgridReader& ur, int numFacesPerBucket) {

    this->ugrid = &ur;

    this->xmin.resize(NUM_SPACE_DIMS);
    this->deltas.resize(NUM_SPACE_DIMS);
    double xmax[NUM_SPACE_DIMS];
    ur.getRange(&this->xmin[0], &xmax[0]);
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        this->deltas[i] = xmax[i] - this->xmin[i];
    }

    size_t numFaces = ur.getNumberOfFaces();

    // number of buckets in x and y directions
    this->nBuckets = std::max(1, (int)(std::sqrt((double)numFaces/((double)numFacesPerBucket))));

    // create empty entries 
    for (size_t i = 0; i < this->nBuckets; ++i) {
        for (size_t j = 0; j < this->nBuckets; ++j) {
            // flat index
            size_t k = j + i * this->nBuckets;
            std::pair<size_t, std::vector<vtkIdType> > kv(k, std::vector<vtkIdType>());
            this->buckets.insert(kv);
        }
    }

    for (vtkIdType faceId = 0; faceId < numFaces; ++faceId) {

        // get the vertices of the face, adjusting for periodicity
        std::vector< Vector<double> > points = ur.getFacePointsRegularized(faceId);

        // iterate over the edges of the face
        for (size_t ie = 0; ie < 4; ++ie) {
            const double* pBeg = &points[ie][0];
            const double* pEnd = &points[(ie + 1) % 4][0];
            std::set<size_t> ks = this->getBucketsAlongLine(pBeg, pEnd);
            for (size_t k : ks) {
                std::map<size_t, std::vector<vtkIdType> >::iterator it = this->buckets.find(k);
                it->second.push_back(faceId);
            } // bucket iteration
        } // edge iteration
    } // face iteration
}

std::set<vtkIdType> 
FacesLocator::getFacesAlongLine(const double pBeg[], const double pEnd[]) const {

    std::set<vtkIdType> faceIds;

    std::set<size_t> ks = this->getBucketsAlongLine(pBeg, pEnd);
    for (size_t k : ks) {
        std::map<size_t, std::vector<vtkIdType> >::const_iterator it = this->buckets.find(k);
        for (auto faceId : it->second) {
            faceIds.insert(faceId);
        }
    }
    return faceIds;
}

vtkIdType 
FacesLocator::getFace(const double point[], double xi[]) const {

    double* closestPoint = NULL;
    int subId;
    double dist2;
    double weights[4]; // quad

    Vector<int> iPos = this->getBucketCellLoc(point);
    int i = iPos[0];
    int j = iPos[1];
    if (i < 0 || j < 0) {
        // outside domain
        return -1;
    }

    // bucket index 
    size_t k = i * this->nBuckets + j;

    // find the bucket
    std::map<size_t, std::vector<vtkIdType> >::const_iterator it = this->buckets.find(k);
    if (it == this->buckets.end()) {
        // outside 
        return -1;
    }

    vtkIdType res = -1;
    const double tol = 1.e-3;
    // iterate over the faces in this bucket
    for (auto faceId : it->second) {
        if (this->ugrid->containsPoint(faceId, point, tol)) {
            // found!
            res = faceId;
            std::vector< Vector<double> > points = this->ugrid->getFacePointsRegularized(faceId);

            this->quad->GetPoints()->SetPoint(0, &points[0][0]);
            this->quad->GetPoints()->SetPoint(1, &points[1][0]);
            this->quad->GetPoints()->SetPoint(2, &points[2][0]);
            this->quad->GetPoints()->SetPoint(3, &points[3][0]);

            this->quad->EvaluatePosition((double*) point, closestPoint, subId, xi, dist2, weights);

            break;
        }
    }

    return res;
}

std::set<size_t>
FacesLocator::getBucketsAlongLine(const double pBeg[], const double pEnd[]) const {

    // required to test if a line has any chance of intersecting a cell
    // Because the algorithm is based on computing the shortest distance 
    // of the line to the centre of a bucket, we need to account for the cases
    // where the point on the line closest to the centre is outside but the line
    // is still intersecting the bucket. The "halo" parameter takes care of this. 
    const double halo = 0.12;
    std::set<size_t> res;

    // compute the start/end points of line in bucket space
    Vector<double> bBeg = this->getBucketSpaceLoc(pBeg);
    Vector<double> bEnd = this->getBucketSpaceLoc(pEnd);
    Vector<double> bU = bEnd - bBeg;
    double bUNormSquare = dot(bU, bU);

    // compute the bucket index location of the line
    Vector<int> iBeg = this->getBucketCellLoc(pBeg);
    Vector<int> iEnd = this->getBucketCellLoc(pEnd);

    // select the region of interest
    Vector<int> imin = max(0, min(iBeg, iEnd));
    Vector<int> imax = min((int)this->nBuckets - 1, max(iBeg, iEnd));

    // iterate over all the buckets in the region of interest
    for (int i = imin[0]; i <= imax[0]; ++i) {
        for (int j = imin[1]; j <= imax[1]; ++j) {

            // average the vertex positions of the bucket
            Vector<double> bQuadCentre(NUM_PARAM_DIMS, 0.0);
            bQuadCentre[0] = i + 0.5; bQuadCentre[1] = j + 0.5;

            // check if line intersects with this bucket. Compute the line parameter that
            // minimizes the distance of the cell's vertices to the line. Then check if this
            // point is inside the quad to within some tolerance (halo). Now it can happen 
            // that the clostest point is outside but the line still interects the bucket, e.g.
            // by cutting the corner.  
            double lam = dot(bQuadCentre - bBeg, bU)/bUNormSquare;
            lam = std::min(1.0, lam);
            lam = std::max(0.0, lam);

            // point on the line that is closest to the quad
            Vector<double> bClosestPtOnLine = bBeg + lam*bU;

            // is the point inside the bucket to within some tolerance?
            if (bClosestPtOnLine[0] >= i - halo && bClosestPtOnLine[0] <= i + 1 + halo &&
                bClosestPtOnLine[1] >= j - halo && bClosestPtOnLine[1] <= j + 1 + halo) {

                // add the bucket
                size_t k = i * this->nBuckets + j;
                res.insert(k);

            }
        }
    }

    return res;
}


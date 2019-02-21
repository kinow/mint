#include <mntUgridEdgeReader.h>
#include <netcdf.h>
#include <iostream>
#include <algorithm>
#include <cmath>

#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

UgridEdgeReader::UgridEdgeReader() {
}

/**
 * Destructor
 */
UgridEdgeReader::~UgridEdgeReader() {
}

size_t 
UgridEdgeReader::getNumberOfEdges() const {
    return this->edge2Points.size() / (NUM_SPACE_DIMS*2);
}

void
UgridEdgeReader::getEdge(size_t edgeId, double pBeg[], double pEnd[]) const {
    for (size_t j = 0; j < 3; ++j) {
        pBeg[j] = this->edge2Points[j + NUM_SPACE_DIMS*(0 + edgeId*2)];
        pEnd[j] = this->edge2Points[j + NUM_SPACE_DIMS*(1 + edgeId*2)];
    }
}

void
UgridEdgeReader::buildLocator(int numEdgesPerBucket, double tol) {

    // number of cells in x and y directions
    this->nBuckets = std::max(1, (int)(std::sqrt((double) this->edge2Points.size()/(NUM_SPACE_DIMS*2 * 2.0 * (double) numEdgesPerBucket))));

    double pBeg[NUM_SPACE_DIMS], pEnd[NUM_SPACE_DIMS];
    for (size_t ie = 0; ie < this->getNumberOfEdges(); ++ie) {

        // get the start/end point of the edge
        this->getEdge(ie, pBeg, pEnd);

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

std::vector<size_t> 
UgridEdgeReader::getEdgesAlongLine(const double pBeg[], const double pEnd[]) const {

    Vector<int> iBeg = this->getBucketCellLoc(pBeg);
    Vector<int> iEnd = this->getBucketCellLoc(pEnd);
   Vector<int> i0 = std::min(iBeg, iEnd);
    Vector<int> i1 = std::max(iBeg, iEnd);

    std::vector<size_t> edgeIds;

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


void
UgridEdgeReader::getRange(double xMin[], double xMax[]) const {
    for (size_t i = 0; i < this->xmin.size(); ++i) {
        xMin[i] = this->xmin[i];
        xMax[i] = this->xmin[i] + this->deltas[i];
    }
}

int 
UgridEdgeReader::load(const std::string& filename) {

    int ier = 0;
    int ncid;

    // open the file
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open \"" << filename << "\"\n";
        return 1;
    }

    std::vector<double> points = this->readPoints(ncid);
    if (points.size() == 0) {
        std::cerr << "ERROR: cannot read points\n";
        return 1;
    }

    ier = this->readEdgeConnectivity(ncid, points);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot read edge connectivity\n";
        return 1;
    }

    // close the netcdf file
    ier = nc_close(ncid);
}

std::vector<double>
 UgridEdgeReader::readPoints(int ncid) {

    std::vector<double> points;

    int ndims, ier;
    // max number of dimensions is 10
    int dimids[10];
    size_t dims[10];

    // looking for a variable with standard_name "longitude" and "latitude"

    int varidLon = this->findVariableIdWithStandardName(ncid, "longitude", &ndims, dimids);
    if (varidLon < 0) {
        std::cerr << "ERROR: could not find longitude\n";
        return points;
    }
    int varidLat = this->findVariableIdWithStandardName(ncid, "latitude", &ndims, dimids);
    if (varidLon < 0) {
        std::cerr << "ERROR: could not find latitude\n";
        return points;
    }

    // allocate
    size_t nPoints;
    ier = nc_inq_dimlen(ncid, dimids[0], &nPoints);
    std::vector<double> lons(nPoints);
    std::vector<double> lats(nPoints);

    // read
    ier = nc_get_var_double(ncid, varidLon, &lons[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read longitude\n";
        return points;
    }
    ier = nc_get_var_double(ncid, varidLat, &lats[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read latitude\n";
        return points;
    }

    // store
    points.resize(nPoints * 3);
    for (size_t i = 0; i < nPoints; ++i) {
        points[LON_INDEX + 3*i] = lons[i];
        points[LAT_INDEX + 3*i] = lats[i];
        points[ELV_INDEX + 3*i] = 0.0; // no elevation at this point
    }

    return points;
}

int UgridEdgeReader::readEdgeConnectivity(int ncid, const std::vector<double>& points) {

    // looking for a variable with cf_role "edge_node_connectivity"
    int ndims, ier;
    // max number of dimensions is 10
    int dimids[10];
    size_t dims[10];
    int varid = this->findVariableIdWithCfRole(ncid, "edge_node_connectivity", &ndims, dimids);
    if (varid < 0) {
        std::cerr << "ERROR: could not find edge_node_connectivity\n";
        return 1;
    }
    int startIndex;
    int ier2 = nc_get_att_int(ncid, varid, "start_index", &startIndex);

    // allocate
    size_t nEdges;
    ier = nc_inq_dimlen(ncid, dimids[0], &nEdges);

    // read
    std::vector<vtkIdType> edge2Nodes(nEdges * 2);
    ier = nc_get_var_longlong(ncid, varid, &edge2Nodes[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read edge_node_connectivity\n";
        return 1;
    }

    std::vector<double> diffLonMinusZeroPlus(3);
    this->edge2Points.resize(nEdges * NUM_SPACE_DIMS * 2); // 2 points in 3d for each edge

    for (size_t i = 0; i < nEdges; ++i) {

        vtkIdType i0 = edge2Nodes[0 + i*2] - startIndex;
        vtkIdType i1 = edge2Nodes[1 + i*2] - startIndex;
        std::vector<double> pBeg(&points[i0*NUM_SPACE_DIMS], &points[i0*NUM_SPACE_DIMS] + NUM_SPACE_DIMS);
        std::vector<double> pEnd(&points[i1*NUM_SPACE_DIMS], &points[i1*NUM_SPACE_DIMS] + NUM_SPACE_DIMS);

        // add/subtract 360 deg to reduce the length of the edge
        for (size_t j = 0; j < NUM_SPACE_DIMS; ++j) {
            diffLonMinusZeroPlus[j] = pEnd[LON_INDEX] - pBeg[LON_INDEX] + (j - 1) * 360.0;
        }
        std::vector<double>::iterator it = std::min_element(diffLonMinusZeroPlus.begin(), 
                                                            diffLonMinusZeroPlus.end());
        int indexMin = (int) std::distance(diffLonMinusZeroPlus.begin(), it);
        // fix the longitude
        pEnd[LON_INDEX] += (indexMin - 1) * 360.0;

        // if one of the latitudes is +- 90 deg then set the longitude to 
        // be the same as the other longitude
        if (std::abs(pEnd[LAT_INDEX]) == 90.) {
            pEnd[LON_INDEX] = pBeg[LON_INDEX];
        }
        if (std::abs(pBeg[LAT_INDEX]) == 90.) {
            pBeg[LON_INDEX] = pEnd[LON_INDEX];
        }

        this->xmin.resize(NUM_SPACE_DIMS, +std::numeric_limits<double>::max());
        this->deltas.resize(NUM_SPACE_DIMS);
        std::vector<double> xmax(NUM_SPACE_DIMS, -std::numeric_limits<double>::max());
        // store
        for (size_t j = 0; j < NUM_SPACE_DIMS; ++j) {

            this->edge2Points[j + 3*(0 + i*2)] = pBeg[j];
            this->edge2Points[j + 3*(1 + i*2)] = pEnd[j];

            this->xmin[j] = std::min(this->xmin[j], std::min(pBeg[j], pEnd[j]));

            xmax[j] = std::max(xmax[j], std::max(pBeg[j], pEnd[j]));

            this->deltas[j] = xmax[j] - this->xmin[j];
        }
    }

    return 0;
}


int UgridEdgeReader::findVariableIdWithCfRole(int ncid, const std::string& cf_role, int* ndims, int dimids[]) {

    nc_type xtype;
    int natts, ier;
    char varname[NC_MAX_NAME];

    size_t n = cf_role.size();

    int nvars;
    ier = nc_inq_nvars(ncid, &nvars);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: after inquiring the number of variables (ier = " 
                  << ier << ")\n";
        return -1;
    }

    int res = -1;

    std::string cf_role_input;
    cf_role_input.assign(NC_MAX_NAME, ' ');
    for (int ivar = 0; ivar < nvars; ++ivar) {
        ier = nc_inq_var(ncid, ivar, varname, &xtype, ndims, &dimids[0], &natts);
        int ier2 = nc_get_att_text(ncid, ivar, "cf_role", &cf_role_input[0]);
        if (ier2 == NC_NOERR && cf_role_input.substr(0, n) == cf_role) {
            // found the variable
            res = ivar;
        }
    }

    return res;
}

int UgridEdgeReader::findVariableIdWithStandardName(int ncid, const std::string& standard_name, int* ndims, int dimids[]) {

    nc_type xtype;
    int natts, ier;
    char varname[NC_MAX_NAME];

    size_t n = standard_name.size();

    int nvars;
    ier = nc_inq_nvars(ncid, &nvars);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: after inquiring the number of variables (ier = " 
                  << ier << ")\n";
        return -1;
    }

    int res = -1;

    std::string standard_name_input;
    standard_name_input.assign(NC_MAX_NAME, ' ');
    for (int ivar = 0; ivar < nvars; ++ivar) {
        ier = nc_inq_var(ncid, ivar, varname, &xtype, ndims, &dimids[0], &natts);
        int ier2 = nc_get_att_text(ncid, ivar, "standard_name", &standard_name_input[0]);
        if (ier2 == NC_NOERR && standard_name_input.substr(0, n) == standard_name) {
            // found the variable
            res = ivar;
        }
    }

    return res;
}


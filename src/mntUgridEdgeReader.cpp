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
    return this->edge2Points.size() / (3*2);
}

void
UgridEdgeReader::getEdge(size_t edgeId, double pBeg[], double pEnd[]) const {
    for (size_t j = 0; j < 3; ++j) {
        pBeg[j] = this->edge2Points[j + 0*3 + edgeId*3*2];
        pEnd[j] = this->edge2Points[j + 1*3 + edgeId*3*2];
    }
}

void
UgridEdgeReader::buildLocator(int numEdgesPerBucket) {

    // number of cells in x and y directions
    this->nBuckets = std::max(1, (int)(std::sqrt((double) this->edge2Points.size()/(3 * 2 * 2.0 * (double) numEdgesPerBucket))));

    double pBeg[3], pEnd[3];
    for (size_t ie = 0; ie < this->getNumberOfEdges(); ++ie) {

        // get the start/end point of the edge
        this->getEdge(ie, pBeg, pEnd);

        // compute the bucket index location
        Vector<size_t> iBeg = this->getBucketLoc(pBeg);
        Vector<size_t> iEnd = this->getBucketLoc(pEnd);
        Vector<size_t> i0 = std::min(iBeg, iEnd);
        Vector<size_t> i1 = std::max(iBeg, iEnd);

        // add the edge to all the cells spanned by the start/end index locations
        for (size_t j = i0[1]; j <= i1[1]; ++j) {
            for (size_t i = i0[0]; i <= i1[0]; ++i) {

                // flat index
                size_t k = i + j * this->nBuckets;

                std::map<size_t, std::vector<size_t> >::iterator it = this->buckets.find(k);

                if (it != this->buckets.end()) {
                    // append
                    it->second.push_back(ie);
                }
                else {
                    // create new entry
                    std::pair<size_t, std::vector<size_t> > kv(k, std::vector<size_t>(1, ie));
                    this->buckets.insert(kv);
                }
            }
        }
    }
}

std::vector<size_t> 
UgridEdgeReader::getEdgesAlongLine(const double pBeg[], const double pEnd[]) const {

    Vector<size_t> iBeg = this->getBucketLoc(pBeg);
    Vector<size_t> iEnd = this->getBucketLoc(pEnd);
    Vector<size_t> i0 = min(iBeg, iEnd);
    Vector<size_t> i1 = max(iBeg, iEnd);

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

    std::vector<double>&& points = this->readPoints(ncid);
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
    this->edge2Points.resize(nEdges * 3 * 2); // 2 points in 3d for each edge

    for (size_t i = 0; i < nEdges; ++i) {

        vtkIdType i0 = edge2Nodes[0 + i*2] - startIndex;
        vtkIdType i1 = edge2Nodes[1 + i*2] - startIndex;
        std::vector<double> pBeg(&points[i0*3], &points[i0*3] + 3);
        std::vector<double> pEnd(&points[i1*3], &points[i1*3] + 3);

        // add/subtract 360 deg to reduce the length of the edge
        for (size_t j = 0; j < 3; ++j) {
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

        this->xmin.resize(3, +std::numeric_limits<double>::max());
        this->deltas.resize(3);
        std::vector<double> xmax(3, -std::numeric_limits<double>::max());
        // store
        for (size_t j = 0; j < 3; ++j) {

            this->edge2Points[j + 0*3 + i*3*2] = pBeg[j];
            this->edge2Points[j + 1*3 + i*3*2] = pEnd[j];

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


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
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        pBeg[i] = this->edge2Points[i + (0 + edgeId*2)*NUM_SPACE_DIMS];
        pEnd[i] = this->edge2Points[i + (1 + edgeId*2)*NUM_SPACE_DIMS];
    }
}

void
UgridEdgeReader::getRange(double xMin[], double xMax[]) const {
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        xMin[i] = this->xmin[i];
        xMax[i] = this->xmax[i];
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

    const std::vector<double>& points = this->readPoints(ncid);
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
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot close file\n";
        return 1;
    }
}

std::vector<double>
UgridEdgeReader::readPoints(int ncid) {

    this->xmin.resize(NUM_SPACE_DIMS, +std::numeric_limits<double>::max());
    this->xmax.resize(NUM_SPACE_DIMS, -std::numeric_limits<double>::max());

    std::vector<double> points;

    int ndims, ier;
    // max number of dimensions is 10
    int dimids[10];

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
    points.resize(nPoints * NUM_SPACE_DIMS);
    for (size_t i = 0; i < nPoints; ++i) {
        points[LON_INDEX + NUM_SPACE_DIMS*i] = lons[i];
        points[LAT_INDEX + NUM_SPACE_DIMS*i] = lats[i];
        points[ELV_INDEX + NUM_SPACE_DIMS*i] = 0.0; // no elevation at this point
    }

    return points;
}

int UgridEdgeReader::readEdgeConnectivity(int ncid, const std::vector<double>& points) {

    // looking for a variable with cf_role "edge_node_connectivity"
    int ndims, ier;
    // max number of dimensions is 10
    int dimids[10];
    int varid = this->findVariableIdWithCfRole(ncid, "edge_node_connectivity", &ndims, dimids);
    if (varid < 0) {
        std::cerr << "ERROR: could not find edge_node_connectivity\n";
        return 1;
    }
    int startIndex = 0;
    nc_get_att_int(ncid, varid, "start_index", &startIndex);

    // allocate
    size_t nEdges;
    ier = nc_inq_dimlen(ncid, dimids[0], &nEdges);

    // read
    std::vector<vtkIdType> edge2Nodes(nEdges * 2); // two nodes per edge
    ier = nc_get_var_longlong(ncid, varid, &edge2Nodes[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read edge_node_connectivity\n";
        return 1;
    }

    std::vector<double> diffLonMinusZeroPlus(3); // - 360, 0, +360
    this->edge2Points.resize(nEdges * NUM_SPACE_DIMS * 2); // 2 points in 3d for each edge

    for (size_t i = 0; i < nEdges; ++i) {

        vtkIdType i0 = edge2Nodes[0 + i*2] - startIndex;
        vtkIdType i1 = edge2Nodes[1 + i*2] - startIndex;

        const double* p0 = &points[i0*NUM_SPACE_DIMS];
        const double* p1 = &points[i1*NUM_SPACE_DIMS];
        std::vector<double> pBeg(p0, p0 + NUM_SPACE_DIMS);
        std::vector<double> pEnd(p1, p1 + NUM_SPACE_DIMS);

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

        // store
        for (size_t j = 0; j < NUM_SPACE_DIMS; ++j) {

            this->edge2Points[j + NUM_SPACE_DIMS*(0 + i*2)] = pBeg[j];
            this->edge2Points[j + NUM_SPACE_DIMS*(1 + i*2)] = pEnd[j];

            this->xmin[j] = std::min(this->xmin[j], std::min(pBeg[j], pEnd[j]));
            this->xmax[j] = std::max(this->xmax[j], std::max(pBeg[j], pEnd[j]));
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


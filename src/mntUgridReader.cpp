#include <mntUgridReader.h>
#include <netcdf.h>
#include <iostream>
#include <set>
#include <cmath>
#include <algorithm>

#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2


std::set< Vector<double> > 
UgridReader::getFacePoints(long long faceId) const {

    std::set< Vector<double> > res;

    // iterate over the 4 edges
    for (size_t i = 0; i < 4; ++i) {

        // get the edge id
        long long edgeId = this->face2Edges[i + faceId*4];

        // iterate over the points spanning the edge
        for (auto point : this->getEdgePoints(edgeId)) {
            res.insert(point);
        }
    }

    return res;
}


std::vector< Vector<double> > 
UgridReader::getEdgePoints(long long edgeId) const {

    std::vector< Vector<double> > res;

    // itereate over the 2 points spanning the edge
    for (size_t i = 0; i < 2; ++i) {

        // get the point id
        long long pointId = this->edge2Points[i + edgeId*2];

        // get the coordinates of this point
        const double* p = this->getPoint(pointId);

        // add
        res.push_back( Vector<double>(p, p + NUM_SPACE_DIMS) );
    }
    return res;
}


void
UgridReader::getRange(double xMin[], double xMax[]) const {
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        xMin[i] = this->xmin[i];
        xMax[i] = this->xmax[i];
    }
}

int 
UgridReader::load(const std::string& filename) {

    int ier = 0;
    int ncid;

    // open the file
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open \"" << filename << "\"\n";
        return 1;
    }

    ier = this->readPoints(ncid);
    if (ier != 0) {
        std::cerr << "ERROR: cannot read points\n";
        nc_close(ncid);
        return 1;
    }

    size_t n;
    int varid;
    int startIndex;

    varid = this->findVariableIdWithAttribute(ncid, 
                  "cf_role", "edge_node_connectivity", &n);
    this->edge2Points.resize(n);
    ier = nc_get_var_longlong(ncid, varid, &this->edge2Points[0]);
    // we're using zero based indexing
    startIndex = 0;
    nc_get_att_int(ncid, varid, "start_index", &startIndex);
    for (size_t i = 0; i < n; ++i) {
        this->edge2Points[i] -= startIndex;
    }
    this->numEdges = n / 2;

    varid = this->findVariableIdWithAttribute(ncid, 
                  "cf_role", "face_edge_connectivity", &n);
    this->face2Edges.resize(n);
    ier = nc_get_var_longlong(ncid, varid, &this->face2Edges[0]);
    // we're using zero based indexing
    startIndex = 0;
    nc_get_att_int(ncid, varid, "start_index", &startIndex);
    for (size_t i = 0; i < n; ++i) {
        this->face2Edges[i] -= startIndex;
    }
    this->numFaces = n / 4;

    // close the netcdf file
    ier = nc_close(ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot close file\n";
        return 1;
    }

    this->fixPeriodicity();

    return 0;
}

int
UgridReader::readPoints(int ncid) {

    this->xmin.resize(NUM_SPACE_DIMS, +std::numeric_limits<double>::max());
    this->xmax.resize(NUM_SPACE_DIMS, -std::numeric_limits<double>::max());

    int ier;

    // looking for a variable with standard_name "longitude" and "latitude"

    this->numPoints = 0;
    int varidLon = this->findVariableIdWithAttribute(ncid, 
                         "long_name", "longitude of 2D mesh nodes", &this->numPoints);
    if (varidLon < 0) {
        std::cerr << "ERROR: could not find longitude\n";
        return 1;
    }
    std::vector<double> lons(this->numPoints);

    int varidLat = this->findVariableIdWithAttribute(ncid, 
                         "long_name", "latitude of 2D mesh nodes", &this->numPoints);
    if (varidLon < 0) {
        std::cerr << "ERROR: could not find latitude\n";
        return 1;
    }
    std::vector<double> lats(this->numPoints);

    // read
    ier = nc_get_var_double(ncid, varidLon, &lons[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read longitude\n";
        return 1;
    }
    ier = nc_get_var_double(ncid, varidLat, &lats[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read latitude\n";
        return 1;
    }

    // store
    this->points.resize(this->numPoints * NUM_SPACE_DIMS);
    for (size_t i = 0; i < this->numPoints; ++i) {
        this->points[LON_INDEX + NUM_SPACE_DIMS*i] = lons[i];
        this->points[LAT_INDEX + NUM_SPACE_DIMS*i] = lats[i];
        this->points[ELV_INDEX + NUM_SPACE_DIMS*i] = 0.0; // no elevation at this point

        for (size_t j = 0; j < NUM_SPACE_DIMS; ++j) {
            double p = this->points[j + NUM_SPACE_DIMS*i];
            this->xmin[j] = (p < this->xmin[j]? p: this->xmin[j]);
            this->xmax[j] = (p > this->xmax[j]? p: this->xmax[j]);
        }
    }

    return 0;
}



int 
UgridReader::findVariableIdWithAttribute(int ncid, 
                                             const std::string& attrname, 
                                             const std::string& attrvalue, 
                                             size_t* nsize) {

    // initialize
    *nsize = 0;
    int res = -1;
    int ier;

    // find the number of variables
    int nvars;
    ier = nc_inq_nvars(ncid, &nvars);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: after inquiring the number of variables (ier = " 
                  << ier << ")\n";
        return -1;
    }


    size_t attsize = attrvalue.size();
    std::string attv(NC_MAX_NAME, ' '); // attribute value
    for (int ivar = 0; ivar < nvars; ++ivar) {
        int ier2 = nc_get_att_text(ncid, ivar, &attrname[0], &attv[0]);
        if (ier2 == NC_NOERR && attv.substr(0, attsize) == attrvalue) {

            // found the variable, store its id
            res = ivar;

            // get the number of dimensions
            int ndims;
            nc_inq_varndims(ncid, ivar, &ndims);

            // get the dimension ids
            std::vector<int> dimids(ndims);
            nc_inq_vardimid(ncid, ivar, &dimids[0]);

            // number of elements
            *nsize = 1;
            for (size_t i = 0; i < dimids.size(); ++i) {
                size_t ns;
                nc_inq_dimlen(ncid, dimids[i], &ns);
                *nsize *= ns;
            }
        }
    }

    return res;
}

std::vector< Vector<double> > 
UgridReader::getFacePointsRegularized(long long faceId) {

    std::vector< Vector<double> > res(4);

    const long long* edgeIds = this->getFaceEdgeIds(faceId);
    std::set<long long> pointIds;
    for (int i = 0; i < 4; ++i) {
        long long edgeId = edgeIds[i];
        const long long* ptIds = this->getEdgePointIds(edgeId);
        for (int j = 0; j < 2; ++j) {
            pointIds.insert(ptIds[j]);
        }
    }

    int i = 0;
    for (auto pId : pointIds) {
        const double* p = this->getPoint(pId);
        res[i] = Vector<double>(p, p + NUM_SPACE_DIMS);
    }

    // regularize
    for (size_t i = 1; i < 4; ++i) {

        // add/subtract 360 
        double dLon = res[i][LON_INDEX] - res[0][LON_INDEX];
        double dLonsPM360[] = {dLon - 360., dLon, dLon + 360.};
        double* minDLon = std::min_element(&dLonsPM360[0], &dLonsPM360[3]);
        int indexMin = (int) std::distance(dLonsPM360, minDLon);
        res[i][LON_INDEX] += (indexMin - 1)*360.0;

        // if on pole then set the longitude to be that of the base node
        if (std::abs(res[i][LAT_INDEX]) == 90.) {
            res[i][LON_INDEX] = res[0][LON_INDEX];
        }
    }

    // order the points to go anticlockwise
    Vector<double> dp10 = res[1] - res[0];
    Vector<double> dp20 = res[2] - res[0];
    Vector<double> dp30 = res[3] - res[0];
    if (dp10[LON_INDEX]*dp20[LAT_INDEX] - dp10[LAT_INDEX]*dp20[LON_INDEX] < 0.) {
        // swap 1 <-> 2
        Vector<double> tmp = res[1];
        res[1] = res[2];
        res[2] = tmp;
    }
    if (dp20[LON_INDEX]*dp30[LAT_INDEX] - dp20[LAT_INDEX]*dp30[LON_INDEX] < 0.) {
        // swap 2 <-> 3
        Vector<double> tmp = res[2];
        res[2] = res[3];
        res[3] = tmp;
    }

    return res;
}


std::vector< Vector<double> > 
UgridReader::getEdgePointsRegularized(long long edgeId) {

    const long long* ptIds = this->getEdgePointIds(edgeId);
    const double* p0 = this->getPoint(ptIds[0]);
    const double* p1 = this->getPoint(ptIds[1]);

    std::vector< Vector<double> > res(2);
    res[0].assign(p0, p0 + NUM_SPACE_DIMS);
    res[1].assign(p1, p1 + NUM_SPACE_DIMS);

    // fix the longitude to minimize the edge length
    double dLon = p1[LON_INDEX] - p0[LON_INDEX];
    double dLonsPM360[] = {dLon - 360., dLon, dLon + 360.};
    double* minDLon = std::min_element(&dLonsPM360[0], &dLonsPM360[3]);
    int indexMin = (int) std::distance(dLonsPM360, minDLon);
    res[1][LON_INDEX] += (indexMin - 1)*360.0;

    // fix the latitude to minimize the edge length
    if (std::abs(res[0][LAT_INDEX]) == 90.) {
        // lon is not well defined at the pole, choose to be the same as the other lon
        res[0][LON_INDEX] = res[1][LON_INDEX];
    }
    if (std::abs(res[1][LAT_INDEX]) == 90.) {
        // lon is not well defined at the pole, choose to be the same as the other lon
        res[1][LON_INDEX] = res[0][LON_INDEX];
    }

    return res;
}


void
UgridReader::fixPeriodicity() {

    std::vector<double> newPoints;
    newPoints.reserve(this->numPoints/10);

    for (size_t edgeId = 0; edgeId < this->numEdges; ++edgeId) {

        const long long* pointIds = this->getEdgePointIds(edgeId);

        const double* p0 = this->getPoint(pointIds[0]);
        const double* p1 = this->getPoint(pointIds[1]);
        
        double dLon = p1[LON_INDEX] - p0[LON_INDEX];
        const double dLonsPM360[] = {dLon - 360., dLon, dLon + 360};

        const double* minDLon = std::min_element(&dLonsPM360[0], &dLonsPM360[3]);
        int indexMin = (int) std::distance(dLonsPM360, minDLon);

        if (indexMin != 1) {

            // needs a correction

            // point Id of the new vertex
            size_t newPointId = (this->points.size() + newPoints.size()) / NUM_SPACE_DIMS;

            // new vertex
            Vector<double> newPoint(p1, p1 + NUM_SPACE_DIMS);

            // apply periodic correction
            newPoint[LON_INDEX] += (indexMin - 1)*360.0;

            // update min/max of domain
            double lon = newPoint[LON_INDEX];
            this->xmin[LON_INDEX] = (lon < this->xmin[LON_INDEX]? lon: this->xmin[LON_INDEX]);
            this->xmax[LON_INDEX] = (lon > this->xmax[LON_INDEX]? lon: this->xmax[LON_INDEX]);

            for (size_t j = 0; j < NUM_SPACE_DIMS; ++j) {
                newPoints.push_back(newPoint[j]);
            }

            // update the edge to node connectivity
            this->edge2Points[1 + edgeId*2] = newPointId;
        }

    }

    // now insert the new points
    this->points.insert(this->points.end(), newPoints.begin(), newPoints.end());

}


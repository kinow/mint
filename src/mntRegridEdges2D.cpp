#include <mntRegridEdges2D.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <netcdf.h>

/**
 * Compute the inrterpolation weight between a source cell edge and a destination line segment
 * @param srcXi0 start point of src edge
 * @param srcXi1 end point of the src edge
 * @param dstXi0 start point of target line
 * @param dstXi1 end point of target line
 * @return interpolation weight
 */
double computeWeight(const Vector<double>& srcXi0, const Vector<double>& srcXi1,
                     const Vector<double>& dstXi0, const Vector<double>& dstXi1) {

    double weight = 1.0;
    for (size_t d = 0; d < 2; ++d) {

        // mid point of target
        double dstXiMid = 0.5*(dstXi0[d] + dstXi1[d]);
        // length of target
        double dstXiLen = dstXi1[d] - dstXi0[d];

        // mid point of src edge in parameter space
        double x = 0.5*(srcXi0[d] + srcXi1[d]);

        // use Lagrange interpolation to evaluate the basis function integral for
        // any for the 3 possible x values in {0, 0.5, 1}. This formula will make 
        // it easier to extend the code to 3d
        double xm00 = x;
        double xm05 = x - 0.5;
        double xm10 = x - 1.0;
        double lag00 = + 2. * xm05 * xm10;
        double lag05 = - 4. * xm00 * xm10;
        double lag10 = + 2. * xm00 * xm05;

        weight *= (1.0 - dstXiMid)*lag00 + dstXiLen*lag05 + dstXiMid*lag10;

        // fix sign when src edge's direction is negative in src cell param coords
        // expecting srcXi1[d] - srcXi0[d] to be -1, 0 or 1
        weight *= (srcXi1[d] - srcXi0[d] < -0.001? -1.0: 1.0);
    }

    return weight;
}


/**
 * Extract the file name and mesh name from non-zero terminated fortran string 
 * @param fort_filename eg "fileName:meshName"
 * @param n length o filename
 * @param fileName file name (output)
 * @param meshName mesh name (output)
 */
void getFileAndMeshNames(const char* fort_filename, int n, 
                        std::string& fileName, std::string& meshName) {

    std::string fileAndMeshName = std::string(fort_filename, n);

    size_t column = fileAndMeshName.find(':');
    fileName = fileAndMeshName.substr(0, column);

    if (column != std::string::npos) {
        meshName = fileAndMeshName.substr(column + 1);
    }
}

extern "C"
int mnt_regridedges2d_new(RegridEdges2D_t** self) {

    *self = new RegridEdges2D_t();

    (*self)->numSrcEdges = 0;
    (*self)->numSrcFaces = 0;

    (*self)->numDstEdges = 0;
    (*self)->numDstFaces = 0;

    return 0;
}

extern "C"
int mnt_regridedges2d_del(RegridEdges2D_t** self) {

    delete *self;

    return 0;
}

extern "C"
int mnt_regridedges2d_loadEdgeField(RegridEdges2D_t** self,
                                    const char* fort_filename, int nFilenameLength,
                                    const char* field_name, int nFieldNameLength,
                                    size_t ndata, double data[]) {

    std::string fileAndMeshName = std::string(fort_filename, nFilenameLength);

    // filter out the mesh name, if present (not used here)
    size_t columnL = fileAndMeshName.find(':');
    std::string filename = fileAndMeshName.substr(0, columnL);

    std::string fieldname = std::string(field_name, nFieldNameLength);

    // open the file
    int ncid;
    int ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open \"" << filename << "\"\n";
        nc_close(ncid);
        return 1;
    }

    // check if the variable/field exists
    int varId;
    ier = nc_inq_varid(ncid, fieldname.c_str(), &varId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not find variable \"" << fieldname << "\"\n";
        nc_close(ncid);
        return 1;
    }

    // check that the field has the "location" attribute
    size_t nLoc;
    ier = nc_inq_attlen(ncid, varId, "location", &nLoc);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: variable \"" << fieldname << "\" does not appear to have attribute 'location' (ier = " << ier << ")\n";
        nc_close(ncid);
        return 2;
    }
    char location[nLoc + 1];
    ier = nc_get_att_text(ncid, varId, "location", location);
    location[nLoc] = '\0';
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: attribute \"location\" of variable \"" << fieldname << "\" could not be read (ier = " << ier << ")\n";
        nc_close(ncid);
        return 6;
    }
    // check location is set to "edge"
    if (strcmp(location, "edge") != 0) {
        std::cerr << "ERROR: attribute \"location\" of variable " << fieldname << " is not edge  ("
                  << location << ")\n";
        nc_close(ncid);
        return 3;
    }

    // check if the data has the right dimension
    int ndims;
    ier = nc_inq_varndims(ncid, varId, &ndims);
    int dimIds[ndims];
    ier = nc_inq_vardimid(ncid, varId, dimIds);
    size_t n;
    ier = nc_inq_dimlen(ncid, dimIds[0], &n);
    if (n != ndata) {
        std::cerr << "ERROR: size of \"" << fieldname << "\" should be " << n
                  << " but got " << ndata << "\n";
        nc_close(ncid);
        return 5;        
    }

    // TO DO 
    // is there are way to check if a field is an edge integral vs a vector field? 
    // Assume field is a line integral

    // now read
    ier = nc_get_var_double(ncid, varId, data);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: while reading variable '" << fieldname << "'\n";
        nc_close(ncid);
        return 4;
    }

    // close the netcdf file
    ier = nc_close(ncid);

    return 0;
}

extern "C"
int mnt_regridedges2d_dumpEdgeField(RegridEdges2D_t** self,
                                    const char* fort_filename, int nFilenameLength,
                                    const char* field_name, int nFieldNameLength,
                                    size_t ndata, const double data[]) {
    
    std::string filename = std::string(fort_filename, nFilenameLength);
    std::string fieldname = std::string(field_name, nFieldNameLength);

    int ncid, ier;
    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not create file \"" << filename << "\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // create dimensions
    int numEdgesId;

    ier = nc_def_dim(ncid, "num_edges", ndata, &numEdgesId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_edges\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }    

    // create variable
    int dataId;
    int dims[] = {numEdgesId};
    ier = nc_def_var(ncid, fieldname.c_str(), NC_DOUBLE, 1, dims, &dataId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"data\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    // write the data
    ier = nc_put_var_double(ncid, dataId, data);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"data\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    // close the netcdf file
    ier = nc_close(ncid);    

    return 0;
}

extern "C"
int mnt_regridedges2d_checkSrcGrid(RegridEdges2D_t** self, double tol) {

    std::vector<size_t> badFaceIds = (*self)->srcGrid.getNegativeFaces(tol);

    Vector<double> pA(3, 0.0);
    Vector<double> pB(3, 0.0);
    for (const size_t& faceId : badFaceIds) {
        const size_t* pointIds = (*self)->srcGrid.getFacePointIds(faceId);
        const size_t* edgeIds = (*self)->srcGrid.getFaceEdgeIds(faceId);
        std::cout << "Negative area detected for source grid face " << faceId << " with nodes ";
        for (size_t i = 0; i < 4; ++i) { // 4 nodes per face
            size_t pointId = pointIds[i];
            std::cout << pointId << " ";
        }
        std::cout << '\n';
        for (size_t i = 0; i < 4; ++i) { // 4 edges per face
            size_t edgeId = edgeIds[i];
            const size_t* pointIds = (*self)->srcGrid.getEdgePointIds(edgeId);
            (*self)->srcGrid.getEdgeBegEndPointsRegularized(edgeId, faceId, pA, pB);
            std::cout << '\t' << pointIds[0] << '(' << pA << ")----" << edgeId << "---->" << pointIds[1] << '(' << pB << ")\n";
        }
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_checkDstGrid(RegridEdges2D_t** self, double tol) {

    std::vector<size_t> badFaceIds = (*self)->dstGrid.getNegativeFaces(tol);
    
    Vector<double> pA(3, 0.0);
    Vector<double> pB(3, 0.0);
    for (const size_t& faceId : badFaceIds) {
        const size_t* pointIds = (*self)->dstGrid.getFacePointIds(faceId);
        const size_t* edgeIds = (*self)->dstGrid.getFaceEdgeIds(faceId);
        std::cout << "Negative area detected for destination grid face " << faceId << " with nodes ";
        for (size_t i = 0; i < 4; ++i) { // 4 nodes per face
            size_t pointId = pointIds[i];
            std::cout << pointId << " ";
        }
        for (size_t i = 0; i < 4; ++i) { // 4 edges per face
            size_t edgeId = edgeIds[i];
            const size_t* pointIds = (*self)->dstGrid.getEdgePointIds(edgeId);
            (*self)->srcGrid.getEdgeBegEndPointsRegularized(edgeId, faceId, pA, pB);
            std::cout << '\t' << pointIds[0] << '(' << pA << ")----" << edgeId << "---->" << pointIds[1] << '(' << pB << ")\n";
        }
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_loadSrcGrid(RegridEdges2D_t** self, 
		                          const char* fort_filename, int n) {

    int ier;
    std::string fileName, meshName;

    getFileAndMeshNames(fort_filename, n, fileName, meshName);

    if(fileName.size() == 0) return 1; // could not extract the file name
    if(meshName.size() == 0) return 2; // could not extract the mesh name

    ier = (*self)->srcGrid.load(fileName, meshName);

    (*self)->numSrcEdges = (*self)->srcGrid.getNumberOfEdges();
    (*self)->numSrcFaces = (*self)->srcGrid.getNumberOfFaces();

    return ier;
}

extern "C"
int mnt_regridedges2d_loadDstGrid(RegridEdges2D_t** self, 
		                          const char* fort_filename, int n) {

    int ier;
    std::string fileName, meshName;

    getFileAndMeshNames(fort_filename, n, fileName, meshName);

    if(fileName.size() == 0) return 1; // could not extract the file name
    if(meshName.size() == 0) return 2; // could not extract the mesh name

    ier = (*self)->dstGrid.load(fileName, meshName);

    (*self)->numDstEdges = (*self)->dstGrid.getNumberOfEdges();
    (*self)->numDstFaces = (*self)->dstGrid.getNumberOfFaces();

    return ier;
}

extern "C"
int mnt_regridedges2d_build(RegridEdges2D_t** self, int numCellsPerBucket) {

    // checks
    if ((*self)->numSrcEdges == 0) {
        std::cerr << "mnt_regridedges_build: ERROR must load source grid!\n";
        return 1;
    }
    if ((*self)->numDstEdges == 0) {
        std::cerr << "mnt_regridedges_build: ERROR must load destination grid!\n";
        return 2;
    }

    (*self)->weights.clear();

    Vector<double>  dstXi0(3, 0.0);
    Vector<double>  dstXi1(3, 0.0);
    Vector<double>  srcXi0(3, 0.0);
    Vector<double>  srcXi1(3, 0.0);

    // build the source grid locator
    (*self)->srcGrid.buildLocator(numCellsPerBucket);

    // compute the weights

    Vector<double> pA(3, 0.0);
    Vector<double> pB(3, 0.0);
    for (size_t dstCellId = 0; dstCellId < (*self)->numDstFaces; ++dstCellId) {

        // get the edge Ids for this face
        const size_t* dstEdgeIds = (*self)->dstGrid.getFaceEdgeIds(dstCellId);

        // iterate over the four edges
        for (size_t i = 0; i < 4; ++i) { // four edges per face

            size_t dstEdgeId = dstEdgeIds[i];

            (*self)->dstGrid.getEdgeBegEndPointsRegularized(dstEdgeId, dstCellId, pA, pB);
            const Vector<double> u = pB - pA;

            // find all the intersections between the dst edge and the source grid
            std::vector< std::pair<size_t, std::vector<double> > > intersections = 
                (*self)->srcGrid.findIntersectionsWithLine(pA, pB);
            
            for (const std::pair<size_t, std::vector<double> >& srcCellIdLambdas : intersections) {

                size_t srcCellId = srcCellIdLambdas.first;

                // compute and set the regularized nodes of the cell
                (*self)->srcGrid.setCellPoints(srcCellId);
                double lambda0 = srcCellIdLambdas.second[0];
                double lambda1 = srcCellIdLambdas.second[1];

                // start/end points inside src grid cell
                Vector<double> dstPoint0 = pA + lambda0*u;
                Vector<double> dstPoint1 = pB + lambda1*u;

                // compute the src cell parametric coords of the dst edge segment
                bool inside;
                inside = (*self)->srcGrid.getParamCoords(dstPoint0, &dstXi0[0]);
                inside = (*self)->srcGrid.getParamCoords(dstPoint1, &dstXi1[0]);

                // iterate over the edges of the src cell
                const size_t* srcEdgeIds = (*self)->srcGrid.getFaceEdgeIds(srcCellId);

                for (size_t i = 0; i < 4; ++i) { // 2d (4 edges)

                    size_t srcEdgeId = srcEdgeIds[i];

                    // find the the start/end points of srcEdgeId on the srcCellId side
                    Vector<double> srcPoint0(3, 0.0);
                    Vector<double> srcPoint1(3, 0.0);
                    int ier = (*self)->srcGrid.getEdgeBegEndPointsRegularized(srcEdgeId, srcCellId,
                                                                              srcPoint0, srcPoint1);

                    // compute the src cell parametric coords of the src edge
                    inside = (*self)->srcGrid.getParamCoords(srcPoint0, &srcXi0[0]); // need to check
                    inside = (*self)->srcGrid.getParamCoords(srcPoint1, &srcXi1[0]); // need to check

                    // compute the interpolation weight
                    double weight = computeWeight(srcXi0, srcXi1, dstXi0, dstXi1);

                    std::pair<size_t, size_t> ds((long long) dstEdgeId, (long long) srcEdgeId);
                    std::map< std::pair<size_t, size_t>, double >::iterator it = (*self)->weights.find(ds);

                    if (it == (*self)->weights.end()) {
                        // new entry
                        (*self)->weights.insert( std::pair< std::pair<size_t, size_t>, double >(ds, weight) );
                    }
                    else {
                        // update
                        it->second /= 2.0;
                        it->second += 0.5*weight;
                    } 
                } // end of src edge iteration
            } // end of segment iteration
        } // end of dst edge iteration
    } // end of dst cell iteration
    return 0;
}

extern "C"
int mnt_regridedges2d_getNumSrcEdges(RegridEdges2D_t** self, size_t* nPtr) {
    *nPtr = (*self)->numSrcEdges;
    return 0;
}

extern "C"
int mnt_regridedges2d_getNumDstEdges(RegridEdges2D_t** self, size_t* nPtr) {
    *nPtr = (*self)->numDstEdges;
    return 0;
}

extern "C"
int mnt_regridedges2d_getSrcEdgePointsRegularized(RegridEdges2D_t** self, 
                                                  size_t srcEdgeId, size_t srcFaceId,
                                                  double p0[], double p1[]) {

    Vector<double> pA(3, 0.0);
    Vector<double> pB(3, 0.0);
    (*self)->srcGrid.getEdgeBegEndPointsRegularized(srcEdgeId, srcFaceId, pA, pB);
    for (size_t i = 0; i < 3; ++i) {
        p0[i] = pA[i];
        p1[i] = pB[i];
    }
    return 0;
}

extern "C"
int mnt_regridedges2d_getDstEdgePointsRegularized(RegridEdges2D_t** self,
                                                  size_t dstEdgeId, size_t dstFaceId,
                                                  double p0[], double p1[]) {

    Vector<double> pA(3, 0.0);
    Vector<double> pB(3, 0.0);
    (*self)->dstGrid.getEdgeBegEndPointsRegularized(dstEdgeId, dstFaceId, pA, pB);
    for (size_t i = 0; i < 3; ++i) {
        p0[i] = pA[i];
        p1[i] = pB[i];
    }
    return 0;
}

extern "C"
int mnt_regridedges2d_apply(RegridEdges2D_t** self, 
	                        const double src_data[], double dst_data[]) {


    if ((*self)->numSrcEdges == 0 || (*self)->numDstEdges == 0) {
        std::cerr << "ERROR: looks like the src grid connectivity is not set.\n";
        std::cerr << "Typically this would occur if you did not read the grid\n";
        std::cerr << "from the netcdf Ugrid file.\n";
        return 1;
    }

    // initialize the dst data to zero
    for (size_t i = 0; i < (*self)->numDstEdges; ++i) {
        dst_data[i] = 0.0;
    }

    std::cerr << "... size of weights = " << (*self)->weights.size() << '\n';
    // add the contributions from each src edge
    for (const std::pair< std::pair<size_t, size_t>, double >& dsw : (*self)->weights) {

        size_t dstEdgeId = dsw.first.first;
        size_t srcEdgeId = dsw.first.second;
        double weight = dsw.second;

        dst_data[dstEdgeId] += weight * src_data[srcEdgeId];
        if (std::abs(weight) > 1.e-10) std::cerr << "... dstEdgeId = " << dstEdgeId << " srcEdgeId = " << srcEdgeId << " weight = " << weight << '\n';
    }

    return 0;
}


extern "C"
int mnt_regridedges2d_loadWeights(RegridEdges2D_t** self, 
                                const char* fort_filename, int n) {

    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);

    int ncid, ier;
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not open file \"" << filename << "\"!\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // get the sizes
    size_t numWeights;
    int numWeightsId;
    ier = nc_inq_dimid(ncid, "num_weights", &numWeightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not inquire dimension \"num_weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }
    ier = nc_inq_dimlen(ncid, numWeightsId, &numWeights);


    int dstEdgeIdsId, srcEdgeIdsId, weightsId;

    ier = nc_inq_varid(ncid, "dst_edge_ids", &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }
    ier = nc_inq_varid(ncid, "src_edge_ids", &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }
    ier = nc_inq_varid(ncid, "weights", &weightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 7;
    }

    std::vector<long long> weightDstEdgeIds(numWeights);
    std::vector<long long> weightSrcEdgeIds(numWeights);
    std::vector<double> weights(numWeights);

    // read
    ier = nc_get_var_double(ncid, weightsId, &(weights)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read var \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 8;
    }
    ier = nc_get_var_longlong(ncid, dstEdgeIdsId, &weightDstEdgeIds[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_get_var_longlong(ncid, srcEdgeIdsId, &weightSrcEdgeIds[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }

    ier = nc_close(ncid);    
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not close file \"" << filename << "\"!\n";
        std::cerr << nc_strerror (ier);
        return 13;
    }

    (*self)->weights.clear();
    for (size_t i = 0; i < numWeights; ++i) {
        std::pair<size_t, size_t> ds((size_t) weightDstEdgeIds[i], (size_t) weightSrcEdgeIds[i]);
        double w = weights[i];
        (*self)->weights.insert( std::pair< std::pair<size_t, size_t>, double >(ds, w) );
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_dumpWeights(RegridEdges2D_t** self, 
		                         const char* fort_filename, int n) {

    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);

    size_t numWeights = (*self)->weights.size();

    // re-organize the data in three vectors
    std::vector<long long> weightDstEdgeIds(numWeights);
    std::vector<long long> weightSrcEdgeIds(numWeights);
    std::vector<double> weights(numWeights);
    size_t i = 0;
    for (const std::pair< std::pair<size_t, size_t>, double >& dsw : (*self)->weights) {
        weightDstEdgeIds[i] = (long long) dsw.first.first;
        weightSrcEdgeIds[i] = (long long) dsw.first.second;
        weights[i] = dsw.second;
        i++;
    }

    int ncid, ier;
    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not create file \"" << filename << "\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // create dimensions
    int numWeightsId;
    ier = nc_def_dim(ncid, "num_weights", (int) numWeights, &numWeightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_weights\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }

    // create variables
    int numWeightsAxis[] = {numWeightsId};

    int dstEdgeIdsId;
    ier = nc_def_var(ncid, "dst_edge_ids", NC_INT64, 1, numWeightsAxis, &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_edge_ids\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    int srcEdgeIdsId;
    ier = nc_def_var(ncid, "src_edge_ids", NC_INT64, 1, numWeightsAxis, &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    int weightsId;
    ier = nc_def_var(ncid, "weights", NC_DOUBLE, 1, numWeightsAxis, &weightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 7;
    }

    // close define mode
    ier = nc_enddef(ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not end define mode\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 8;
    }

    ier = nc_put_var_longlong(ncid, dstEdgeIdsId, &weightDstEdgeIds[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"dst_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_put_var_longlong(ncid, srcEdgeIdsId, &weightSrcEdgeIds[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"src_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_put_var_double(ncid, weightsId, &weights[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"weights\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 12;
    }

    ier = nc_close(ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not close file \"" << filename << "\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 13;
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_print(RegridEdges2D_t** self) {

    size_t numWeights = (*self)->weights.size();
    std::cout << "Number of weights: " << numWeights << '\n';

    printf("     index  dstEdgeId  srcEdgeId          weight\n");
    size_t i = 0;
    for (const std::pair< std::pair<size_t, size_t>, double >& dsw : (*self)->weights) {
        printf("%10ld %10ld %10ld %15.5lf\n", i, dsw.first.first, dsw.first.second, dsw.second);
        i++;
    }

    return 0;
}


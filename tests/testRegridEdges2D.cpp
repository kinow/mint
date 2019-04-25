#include <mntRegridEdges2D.h>
#include <cmath>
#undef NDEBUG // turn on asserts

double streamFunc(const Vector<double>& p) {
    double lon = p[0];
    double lat = p[1];
    return cos(M_PI*lat/180.0) * sin(M_PI*lon/180.);
}

void regridEdgeFieldTest(const std::string& testName, const std::string& srcFile, const std::string& dstFile) {

    int ier;
    std::string outputFile = testName + "Weights.nc";

    RegridEdges2D_t* rg;

    ier = mnt_regridedges2d_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges2d_loadSrcGrid(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadSrc...OK\n";

    ier = mnt_regridedges2d_loadDstGrid(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadDst...OK\n";

    int numCellsPerBucket = 8;
    
    assert(ier == 0);
    ier = mnt_regridedges2d_build(&rg, numCellsPerBucket);
    std::cerr << testName << ": build...OK\n";

    std::string weightFile = testName + "Weights.nc";
    ier = mnt_regridedges2d_dumpWeights(&rg, weightFile.c_str(), (int) weightFile.size());
    assert(ier == 0);

    // set the source field
    size_t numSrcEdges;
    ier = mnt_regridedges2d_getNumSrcEdges(&rg, &numSrcEdges);
    assert(ier == 0);

    Vector<double> p0(3), p1(3);
    std::vector<double> srcData(numSrcEdges);
    for (size_t srcEdgeId = 0; srcEdgeId < numSrcEdges; ++srcEdgeId) {
        ier = mnt_regridedges2d_getSrcEdgePointsRegularized(&rg, srcEdgeId, &p0[0], &p1[0]);
        assert(ier == 0);
        srcData[srcEdgeId] = streamFunc(p1) - streamFunc(p0);
    }

    size_t numDstEdges;
    ier = mnt_regridedges2d_getNumDstEdges(&rg, &numDstEdges);
    assert(ier == 0);

    std::vector<double> dstData(numDstEdges);

    // regrid
    ier = mnt_regridedges2d_apply(&rg, &srcData[0], &dstData[0]);
    assert(ier == 0);

    // check
    double totError = 0;
    printf("%s\n dstEdgeId   interpVal     exact        error    p0           p1      \n", testName.c_str());
    for (size_t dstEdgeId = 0; dstEdgeId < numDstEdges; ++dstEdgeId) {

        ier = mnt_regridedges2d_getDstEdgePointsRegularized(&rg, dstEdgeId, &p0[0], &p1[0]);
        assert(ier == 0);

        double exact = streamFunc(p1) - streamFunc(p0);

        double interpVal = dstData[dstEdgeId];

        double error = interpVal - exact;

        if (std::abs(error) > 1.e-6) {
            printf("%10ld %10.3lf %10.3lf %12.5le %5.1lf,%5.1lf; %5.1lf,%5.1lf\n", 
                    dstEdgeId, interpVal, exact, error, p0[0], p0[1], p1[0], p1[1]);
        }
        
        totError += std::abs(error);
    }

    std::cout << testName << ": total interpolation |error|: " << totError << '\n';
    assert(totError < 1.e-8);

    // clean up
    ier = mnt_regridedges2d_del(&rg);
    assert(ier == 0);

}


int main() {

    regridEdgeFieldTest("uniqueEdgeIdField_4->4", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics");
    regridEdgeFieldTest("uniqueEdgeIdField_16->4", "@CMAKE_SOURCE_DIR@/data/cs_16.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics"); 
    regridEdgeFieldTest("uniqueEdgeIdField_16->16", "@CMAKE_SOURCE_DIR@/data/cs_16.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_16.nc:physics"); 

    return 0;
}   

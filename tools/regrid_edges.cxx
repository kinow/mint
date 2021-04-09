#include <mntRegridEdges.h>
#include <mntNcAttributes.h>
#include <mntNcDimensions.h>
#include <mntMultiArrayIter.h>
#include <mntNcFieldRead.h>
#include <mntNcFieldWrite.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <cmath>
#include <regex>
#include <netcdf.h>

#define NUM_EDGES_PER_CELL 4

/**
 * Convert a number to a string
 * @param arg number
 * @return string
 */
template <typename T>
std::string toString (T arg) {
    std::stringstream ss;
    ss << arg;
    return ss.str ();
}

/**
 * Split string by separator
 * @param fmname string, eg "filename:meshname"
 * @param separator separator, eg ':'
 * @return filename, meshname pair
 */
std::pair<std::string, std::string> split(const std::string& fmname, char separator) {
    std::pair<std::string, std::string> res;
    size_t pos = fmname.find(separator);
    res.first = fmname.substr(0, pos);
    if (pos < std::string::npos) {
        // file name substring
        // mesh name substring
        res.second = fmname.substr(pos + 1, std::string::npos);
    }
    return res;
}

/**
 * Compute the loop integrals for each cell
 * @param gridObj grid object
 * @param edge data
 * @param avgAbsLoop abs of average loop integral (output)
 * @param minAbsLoop abs of min loop integral (output)
 * @param maxAbsLoop abs of max loop integral (output)
 * @param loop_integrals loop integrals, array must have number of cells elements
 */
void computeLoopIntegrals(Grid_t* grd, const std::vector<double>& edgeData,
                          double* avgAbsLoop, double* minAbsLoop, double* maxAbsLoop,
                          std::vector<double>& loop_integrals) {
    size_t numCells, edgeId;
    int edgeSign, ier;
    mnt_grid_getNumberOfCells(&grd, &numCells);
    *minAbsLoop = + std::numeric_limits<double>::max();
    *maxAbsLoop = - std::numeric_limits<double>::max();
    *avgAbsLoop = 0.0;
    for (size_t cellId = 0; cellId < numCells; ++cellId) {
        double loop = 0.0;
        for (int ie = 0; ie < NUM_EDGES_PER_CELL; ++ie) {

            ier = mnt_grid_getEdgeId(&grd, cellId, ie, &edgeId, &edgeSign);
            assert(ier == 0);

            // +1 for ie = 0, 1; -1 for ie = 2, 3
            int sgn = 1 - 2*(ie/2);
            
            loop += sgn * edgeSign * edgeData[edgeId];
        }

        loop_integrals[cellId] = loop;
        loop = std::abs(loop);
        *minAbsLoop = std::min(loop, *minAbsLoop);
        *maxAbsLoop = std::max(loop, *maxAbsLoop);
        *avgAbsLoop += loop;
    }
    *avgAbsLoop /= double(numCells);
}


class RegridEdgesApp {

public:

    RegridEdgesApp(int argc, char** argv, bool& argParseSucess);

    ~RegridEdgesApp();

    int checkOptions();

    int readGrids();

    int computeWeights();

    int applyWeights();

private:

    int setUpWriter(int srcNdims, const size_t* srcDims);

    CmdLineArgParser args;
    NcFieldRead_t* reader;
    NcFieldWrite_t* writer;
    NcAttributes_t* attrs;
    NcDimensions_t* srcVarDimsObj;
    RegridEdges_t* rg;

    std::string srcFile;
    std::string dstFile;
    std::string weightsFile;
    std::string loadWeightsFile;
    std::string vtkOutputFile;
    std::string dstEdgeDataFile;

    std::string variableName;

    std::vector<double> srcEdgeData;
    std::vector<double> dstEdgeData;


    size_t numSrcEdges;
    size_t numDstEdges;
    size_t numSrcCells;
    size_t numDstCells;

};

RegridEdgesApp::RegridEdgesApp(int argc, char** argv, bool& success) {

    // initialization
    this->reader = NULL;
    this->writer = NULL;
    this->rg = NULL;
    this->attrs = NULL;
    this->srcVarDimsObj = NULL;

    this->numSrcEdges = 0;
    this->numDstEdges = 0;
    this->numSrcCells = 0;
    this->numDstCells = 0;

    // command line options
    this->args.setPurpose("Regrid an edge centred field.");
    this->args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename:meshname\"");
    this->args.set("-v", std::string(""), "Specify edge staggered field variable name in source UGRID file, varname[@filename:meshname]");
    this->args.set("-P", 0.0, "Specify the periodicity length in longitudes (default is non-periodic)");
    this->args.set("-d", std::string(""), "UGRID destination grid file name");
    this->args.set("-w", std::string(""), "Write interpolation weights to file");
    this->args.set("-W", std::string(""), "Load interpolation weights from file");
    this->args.set("-o", std::string(""), "Specify output VTK file where regridded edge data are saved");
    this->args.set("-O", std::string(""), "Specify output 2D UGRID file where regridded edge data are saved");
    this->args.set("-S", 1, "Set to zero to disable source grid regularization, -S 0 is required for uniform lon-lat grid");
    this->args.set("-D", 1, "Set to zero to disable destination grid regularization, -S 0 is required for uniform lon-lat grid");
    this->args.set("-N", 128, "Average number of cells per bucket");
    this->args.set("-debug", 1, "0=no checks, 1=print outside segments, 2=save outside segments");

    success = this->args.parse(argc, argv);
    bool help = this->args.get<bool>("-h");

    if (help) {
        this->args.help();
    }

    if (!success) {
        std::cerr << "ERROR when parsing command line arguments\n";
        this->args.help();
    }

    // get the file names
    this->srcFile = this->args.get<std::string>("-s");
    this->dstFile = this->args.get<std::string>("-d");
    this->weightsFile = this->args.get<std::string>("-w");
    this->loadWeightsFile = this->args.get<std::string>("-W");
    this->vtkOutputFile = this->args.get<std::string>("-o");
    this->dstEdgeDataFile = this->args.get<std::string>("-O");

    // create regridder 
    mnt_regridedges_new(&this->rg);

    std::cerr << "Done with ctor\n";
}

RegridEdgesApp::~RegridEdgesApp() {

    if (this->rg) mnt_regridedges_del(&this->rg);
    std::cerr << "Done with dtor\n";
}

int 
RegridEdgesApp::checkOptions() {

    // run some checks
    if (this->srcFile.size() == 0) {
        std::cerr << "ERROR: must specify a source grid file (-s)\n";
        return 2;
    }
    if (this->dstFile.size() == 0) {
        std::cerr << "ERROR: must specify a destination grid file (-d)\n";
        return 3;
    }

    std::cerr << "Done with checkOptions\n";
    return 0;
}

int
RegridEdgesApp::readGrids() {

    int ier;

    // defaults are suitable for cubed-sphere 
    int fixLonAcrossDateline = 1;
    int averageLonAtPole = 1;
    if (this->args.get<int>("-S") == 0) {
        fixLonAcrossDateline = 0;
        averageLonAtPole = 0;
        std::cout << "info: no regularization applied to source grid\n";
    }
    ier = mnt_regridedges_setSrcGridFlags(&this->rg, fixLonAcrossDateline, averageLonAtPole);

    // ...destination grid
    fixLonAcrossDateline = 1;
    averageLonAtPole = 1;
    if (this->args.get<int>("-D") == 0) {
        fixLonAcrossDateline = 0;
        averageLonAtPole = 0;
        std::cout << "info: no regularization applied to destination grid\n";
    }
    ier = mnt_regridedges_setDstGridFlags(&this->rg, fixLonAcrossDateline, averageLonAtPole);

    // read the source grid
    ier = mnt_regridedges_loadSrcGrid(&this->rg, this->srcFile.c_str(), this->srcFile.size());
    if (ier != 0) {
        std::cerr << "ERROR: could not read file \"" << this->srcFile << "\"\n";
        return 4;
    }

    // read the destination grid
    ier = mnt_regridedges_loadDstGrid(&this->rg, this->dstFile.c_str(), this->dstFile.size());
    if (ier != 0) {
        std::cerr << "ERROR: could not read file \"" << this->dstFile << "\"\n";
        return 5;
    }

    // get the number of edges and allocate src/dst data
    ier = mnt_regridedges_getNumSrcEdges(&this->rg, &this->numSrcEdges);
    ier = mnt_regridedges_getNumDstEdges(&this->rg, &this->numDstEdges);
    std::cout << "info: number of src edges: " << this->numSrcEdges << '\n';
    std::cout << "info: number of dst edges: " << this->numDstEdges << '\n';

    mnt_regridedges_getNumSrcCells(&this->rg, &this->numSrcCells);
    mnt_regridedges_getNumDstCells(&this->rg, &this->numDstCells);
    std::cout << "info: number of src cells: " << this->numSrcCells << '\n';
    std::cout << "info: number of dst cells: " << this->numDstCells << '\n';

    std::cerr << "Done with readGrids\n";
    return 0;
}

int
RegridEdgesApp::computeWeights() {

    int ier;

    if (this->loadWeightsFile.size() == 0) {

        // compute the weights
        std::cout << "info: computing weights\n";
        int numCellsPerBucket = this->args.get<int>("-N");
        int periodicX = this->args.get<double>("-P");
        int debugFlag = this->args.get<int>("-debug");
        ier = mnt_regridedges_build(&rg, numCellsPerBucket, periodicX, debugFlag);
        if (ier != 0) {
            return 6;
        }
    
        // save the weights to file
        if (this->weightsFile.size() != 0) {
            std::cout << "info: saving weights in file " << this->weightsFile << '\n';
            ier = mnt_regridedges_dumpWeights(&this->rg, 
                                              this->weightsFile.c_str(), (int) this->weightsFile.size());
            if (ier != 0) {
                return 7;
            }
        }
    }
    else {
        // weights have been pre-computed, just load them
        std::cout << "info: loading weights from file " << this->loadWeightsFile << '\n';
        ier = mnt_regridedges_loadWeights(&this->rg, 
                                          this->loadWeightsFile.c_str(), (int) this->loadWeightsFile.size());
        if (ier != 0) {
            return 7;
        }

    }

    std::cerr << "Done with computeWeights\n";
    return 0;
}

/**
 * Set up the writer 
 * @param srcNdims number of dimensions
 * @param srcDims  dimensions of the grid
 * @return error (0 is OK)
 */
int 
RegridEdgesApp::setUpWriter(int srcNdims, const size_t* srcDims) {

    int ier;

    std::pair<std::string, std::string> fm = split(this->dstEdgeDataFile, ':');
    // get the dst file name
    std::string dstFileName = fm.first;

    int n1 = dstFileName.size();
    int n2 = this->variableName.size();

    const int append = 0; // new file
    ier = mnt_ncfieldwrite_new(&this->writer, dstFileName.c_str(), n1, 
                               this->variableName.c_str(), n2, append);
    if (ier != 0) {
        std::cerr << "ERROR: create file " << dstFileName << " with field " 
                  << this->variableName << " in append mode " << append << '\n';
        return 14;
    }

    ier = mnt_ncfieldwrite_setNumDims(&this->writer, srcNdims); // matches the number of source field dimensions
    if (ier != 0) {
        std::cerr << "ERROR: cannot set the number of dimensions for field " 
                  << this->variableName << " in file " << dstFileName << '\n';
        ier = mnt_ncfieldwrite_del(&this->writer);
        return 15;
    }

    // add the field's axes. Assume the dst field dimensions are the same as the src field except for the last
    // num edges dimension


    // add num_edges axis. WE SHOULD GET THIS FROM THE DEST FILE?
    std::string axname = "num_edges";
    int n3 = axname.size();
    ier = mnt_ncfieldwrite_setDim(&this->writer, srcNdims - 1, axname.c_str(), n3, numDstEdges);
    if (ier != 0) {
        std::cerr << "ERROR: setting dimension 0 (" << axname << ") to " << numDstEdges
                  << " for field " << this->variableName << " in file " << dstFileName << '\n';
        ier = mnt_ncfieldwrite_del(&this->writer);
        return 16;
    }

    // add the remaining axes, ASSUME THE ADDITIONAL DST AXES TO MATCH THE SRC AXES
    for (int i = 0; i < srcNdims - 1; ++i) {
        axname = "n_" + toString(srcDims[i]);
        ier = mnt_ncfieldwrite_setDim(&this->writer, i, axname.c_str(), axname.size(), srcDims[i]);
    }

    // add the attributes
    ier = mnt_ncattributes_write(&this->attrs, this->writer->ncid, this->writer->varid);
    if (ier != 0) {
        std::cerr << "ERROR: writing attributes for field " << this->variableName << " in file " << dstFileName << '\n';
        ier = mnt_ncfieldwrite_del(&this->writer);
        return 17;
    }

    std::cerr << "Done with setUpWriter\n";
    return 0;
}


int 
RegridEdgesApp::applyWeights() {

    int ier;
    this->srcEdgeData.resize(this->numSrcEdges);
    this->dstEdgeData.resize(this->numDstEdges);

    std::string varAtFileMesh = this->args.get<std::string>("-v");

    if (varAtFileMesh.size() > 0) {

        // get the variable name and the source file/mesh names
        std::pair<std::string, std::string> vfm = split(varAtFileMesh, '@');
        this->variableName = vfm.first;

        // by default the variable is stored in srcFile
        std::string srcFileMeshName = this->srcFile;
        if (vfm.second.size() > 0) {
            srcFileMeshName = vfm.second;
        }
        std::pair<std::string, std::string> fm = split(srcFileMeshName, ':');
        std::string srcFileName = fm.first;

        // get the ncid and varid's so we can read the attributes and dimensions
        int srcNcid;
        ier = nc_open(srcFileName.c_str(), NC_NOWRITE, &srcNcid);
        if (ier != 0) {
            std::cerr << "ERROR: could not open file \"" << srcFileName << "\"\n";
            return 8;
        }

        int srcVarid;
        ier = nc_inq_varid(srcNcid, this->variableName.c_str(), &srcVarid);
        if (ier != 0) {
            std::cerr << "ERROR: could not find variable \"" << this->variableName 
                    << "\" in file \"" << srcFileName << "\"\n";
            return 9;
        }

        // get the attributes of the variable from the netcdf file
        ier = mnt_ncattributes_new(&this->attrs);

        // read the attributes
        ier = mnt_ncattributes_read(&this->attrs, srcNcid, srcVarid);
        if (ier != 0) {
            std::cerr << "ERROR: could not extract attributes for variable \"" 
                      << this->variableName << "\" in file \"" << srcFileName << "\"\n";
            return 11;
        }

        // read the dimensions
        ier = mnt_ncdimensions_new(&this->srcVarDimsObj);
        ier = mnt_ncdimensions_read(&this->srcVarDimsObj, srcNcid, srcVarid);
        int srcNdims;
        ier = mnt_ncdimensions_getNumDims(&this->srcVarDimsObj, &srcNdims);
        size_t srcDims[srcNdims];
        for (int i = 0; i < srcNdims; ++i) {
            ier = mnt_ncdimensions_get(&this->srcVarDimsObj, i, &srcDims[i]); 
        }      
        ier = mnt_ncdimensions_del(&this->srcVarDimsObj);

        // prepare to read the field 
        ier = mnt_ncfieldread_new(&this->reader, srcNcid, srcVarid);


        if (dstEdgeDataFile.size() > 0) {
            // user provided a file name to store the regridded data

            ier = setUpWriter(srcNdims, srcDims);
            if (ier != 0) {
                return ier;
            }

        }

        std::vector<double> loop_integrals(this->numDstCells);
        std::vector<double> dstCellByCellData(this->numDstCells * NUM_EDGES_PER_CELL);

        // attach field to grid so we can save the data to file
        mnt_grid_attach(&this->rg->dstGridObj, this->variableName.c_str(), NUM_EDGES_PER_CELL, 
                        &dstCellByCellData[0]);

        std::string loop_integral_varname = std::string("loop_integrals_of_") + this->variableName;
        mnt_grid_attach(&this->rg->dstGridObj, loop_integral_varname.c_str(), 1, &loop_integrals[0]);

        //
        // iterate over the axes other than edge. ASSUMES num edges is the last dimension!!!
        //

        // leading indices into the src/dst array for each slice
        std::vector<size_t> srcIndices(srcNdims, 0);
        std::vector<size_t> dstIndices(srcNdims, 0);

        // number of data values to read, regrid and write
        std::vector<size_t> srcCounts(srcNdims, 1);
        srcCounts[srcNdims - 1] = this->numSrcEdges;
        std::vector<size_t> dstCounts(srcNdims, 1);
        dstCounts[srcNdims - 1] = this->numDstEdges;
        std::cerr << "**** 4\n";

        MultiArrayIter_t* mai = NULL;
        // iterate over the non-edge indices only, hence srcNdims - 1
        ier = mnt_multiarrayiter_new(&mai, srcNdims - 1, srcDims);

        // total number of elevation * time values
        size_t numIters;
        ier = mnt_multiarrayiter_getNumIters(&mai, &numIters);
        for (size_t iter = 0; iter < numIters; ++iter) {

            // same indices in this version, in later versions srcIndices and dstIndices 
            // might be different
            ier = mnt_multiarrayiter_getIndices(&mai, &srcIndices[0]);
            ier = mnt_multiarrayiter_getIndices(&mai, &dstIndices[0]);

            // read a slice of the data from file
            std::cout << "info: reading slice " << iter << " of field " 
                      << this->variableName << " from file \"" << srcFileName << "\"\n";
            ier = mnt_ncfieldread_dataSlice(&this->reader, &srcIndices[0], &srcCounts[0], 
                                            &this->srcEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: could not read variable \"" << this->variableName 
                          << "\" from file \"" << srcFileName << "\"\n";
                return 12;
            }

            // apply the weights to the src field
            ier = mnt_regridedges_apply(&this->rg, &this->srcEdgeData[0], &this->dstEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: failed to apply weights to dst field \"" 
                          << this->variableName << "\"\n";
                return 13;
            }

            // compute loop integrals for each cell
            double avgAbsLoop, minAbsLoop, maxAbsLoop;
            computeLoopIntegrals(this->rg->dstGridObj, this->dstEdgeData, 
                                 &avgAbsLoop, &minAbsLoop, &maxAbsLoop, loop_integrals);
            std::cout << "Min/avg/max cell loop integrals: " << minAbsLoop 
                      << "/" << avgAbsLoop << "/" << maxAbsLoop << '\n';

            if (this->vtkOutputFile.size() > 0) {

                // new file name for each elevation, time, etc
                std::string vtkFilename = std::regex_replace(vtkOutputFile, 
                                          std::regex(".vtk"), "_" + toString(iter) + ".vtk");

                // compute the cell by cell data
                size_t dstEdgeId;
                int dstEdgeSign;
                for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
                    for (int ie = 0; ie < 4; ++ie) {
                        ier = mnt_grid_getEdgeId(&this->rg->dstGridObj, dstCellId, ie, &dstEdgeId, &dstEdgeSign);
                        size_t k = dstCellId*NUM_EDGES_PER_CELL + ie;
                        dstCellByCellData[k] = this->dstEdgeData[dstEdgeId] * dstEdgeSign;
                    }
                }

                std::cout << "info: writing \"" << this->variableName << "\" to " << vtkFilename << '\n';
                mnt_grid_dump(&this->rg->dstGridObj, vtkFilename.c_str());
            }

            if (dstEdgeDataFile.size() > 0) {

                // write the slice of data to a netcdf file
                ier = mnt_ncfieldwrite_dataSlice(&this->writer, &dstIndices[0], &dstCounts[0], 
                                                 &this->dstEdgeData[0]);
                if (ier != 0) {
                    std::cerr << "ERROR: writing slice " << iter << " of data for field " 
                              << this->variableName << " in file " << dstEdgeDataFile << '\n';
                    ier = mnt_ncfieldwrite_del(&this->writer);
                    return 18;
                }

            }

            // increment the iterator (next time, elevation....)
            ier = mnt_multiarrayiter_next(&mai);

        }

        ier = mnt_multiarrayiter_del(&mai);

        if (dstEdgeDataFile.size() > 0) {
            ier = mnt_ncfieldwrite_del(&this->writer);
        }

        // must destroy before closing the file
        ier = mnt_ncfieldread_del(&this->reader);
        ier = mnt_ncattributes_del(&this->attrs);

        // done with reading the attributes
        ier = nc_close(srcNcid);

    } // has variable 
    else {
        std::cout << "info: no variable name was provided, thus only computing weights\n";
    }

    std::cerr << "**** 18\n";
    std::cerr << "Done with applyWeights\n";
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    int ier;
    bool argParseSucess;

    RegridEdgesApp rea(argc, argv, argParseSucess);
    if(!argParseSucess) {
        std::cerr << "ERROR: invalid command line arguments\n";
        return ier;
    }

    ier = rea.checkOptions();
    if(ier != 0) {
        std::cerr << "ERROR: invalid options\n";
        return ier;
    }

    ier = rea.readGrids();
    if(ier != 0) {
        std::cerr << "ERROR: failed to load the source or destination grid\n";
        return ier;
    }

    ier = rea.computeWeights();
    if(ier != 0) {
        std::cerr << "ERROR: when computing or loading weights\n";
        return ier;
    }


    ier = rea.applyWeights();
    if(ier != 0) {
        std::cerr << "ERROR: when applying the weights\n";
        return ier;
    }

    return 0;
}

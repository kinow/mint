#include <vector>
#include <string>

#ifndef MNT_LAT_LON
#define MNT_LAT_LON

/**
 * A class to generate a uniform lat-lon grid conforming to the output of UM
 */

struct LatLon_t {
    std::vector<double> lats;
    std::vector<double> lons;
    double fillValue;
    double dLat;
    double dLon;
};

/**
 * Constructor
 * @param self instance of LatLon_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_new(LatLon_t** self);

/**
 * Destructor
 * @param self instance of LatLon_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_del(LatLon_t** self);

/**
 * Set the number of latitude cells
 * @param self instance of LatLon_t
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLatCells(LatLon_t** self, std::size_t n);

/**
 * Set the number of longitude cells
 * @param self instance of LatLon_t
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLonCells(LatLon_t** self, std::size_t n);

/**
 * Build the grid
 * @param self instance of LatLon_t
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_build(LatLon_t** self);

/**
 * Load a grid from file
 * @param self instance of LatLon_t
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_load(LatLon_t** self, const std::string& filename);

/**
 * Dump the grid to file
 * @param self instance of LatLon_t
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_dump(LatLon_t** self, const std::string& filename);

#endif // MNT_LAT_LON

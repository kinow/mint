#include <mntUgridEdgeReader.h>
#include <mntEdgesLocator.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test4() {
    UgridEdgeReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");
    std::cout << "Number of edges: " << uer.getNumberOfEdges() << '\n';
    double xmin[3], xmax[3];
    uer.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';
}

void test16() {
    UgridEdgeReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_16.nc");
    std::cout << "Number of edges: " << uer.getNumberOfEdges() << '\n';
    double xmin[3], xmax[3];
    uer.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';
}

int main(int argc, char** argv) {

	test4();
    test16();

    return 0;
}

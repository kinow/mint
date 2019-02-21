#include <mntUgridEdgeReader.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test() {
    UgridEdgeReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");

    std::cout << "Number of edges: " << uer.getNumberOfEdges() << '\n';
    double xmin[3], xmax[3];
    uer.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';

    Vector<double> pBeg(3), pEnd(3);
    for (size_t ie = 0; ie < uer.getNumberOfEdges(); ++ie) {
        uer.getEdge(ie, &pBeg[0], &pEnd[0]);
        std::cerr << "edge " << ie << " connects points " << pBeg << " -> " << pEnd << '\n';
    }

    int numEdgesPerBucket = 10;
    uer.buildLocator(numEdgesPerBucket);

    double pa[] = {0., -90., 0.};
    double pb[] = {360., +90., 0.};
    double p0[3], p1[3];
    std::vector<size_t> edgeIds = uer.getEdgesAlongLine(pa, pb);
    for (size_t ie = 0; ie < edgeIds.size(); ++ie) {
        uer.getEdge(ie, p0, p1);
        std::cout << "edge Id " << ie << " points " << p0[0] << ',' << p0[1] << " -> " << p1[0] << ',' << p1[1] << '\n';
    }
}

int main(int argc, char** argv) {

	test();

    return 0;
}

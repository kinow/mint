#include <mntUgridEdgeReader.h>
#include <mntEdgesLocator.h>
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


    EdgesLocator el;
    el.setRange(xmin, xmax);
    int numEdgesPerBucket = 10;
    el.build(uer.getEdgePoints(), numEdgesPerBucket);

    double pa[] = {0., -90., 0.};
    double pb[] = {0., -90., 0.}; // {360., +90., 0.};
    double p0[3], p1[3];
    std::set<vtkIdType> edgeIds = el.getEdgesAlongLine(pa, pb);

    std::cout << "Found " << edgeIds.size() << " edges that are likely intersecting line "
              << pa[0] << ',' << pa[1] << ',' << pa[2] << " -> "
              << pb[0] << ',' << pb[1] << ',' << pb[2] << '\n';

    for (auto iet = edgeIds.begin(); iet != edgeIds.end(); ++iet) {
        uer.getEdge(*iet, p0, p1);
        std::cout << "edge Id " << *iet << " points " << p0[0] << ',' << p0[1] << " -> " << p1[0] << ',' << p1[1] << '\n';
    }
}

int main(int argc, char** argv) {

	test();

    return 0;
}

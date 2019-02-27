#include <mntVtkUgridAdaptor.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test4() {
    UgridReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");
    std::cout << "Number of points: " << uer.getNumberOfPoints() << '\n';
    std::cout << "Number of  edges: " << uer.getNumberOfEdges() << '\n';
    std::cout << "Number of  faces: " << uer.getNumberOfFaces() << '\n';

    VtkUgridAdaptor vua(uer);
    vua.write("cs_4_ugrid.vtk");
}

int main(int argc, char** argv) {

	test4();

    return 0;
}

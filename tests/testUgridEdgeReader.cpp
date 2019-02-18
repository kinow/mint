#include <mntUgridEdgeReader.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test() {
    UgridEdgeReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");

}

int main(int argc, char** argv) {

	test();

    return 0;
}

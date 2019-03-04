#include <mntUgridReader.h>
#include <mntFacesLocator.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test4() {
    UgridReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");
    std::cout << "Number of points: " << uer.getNumberOfPoints() << '\n';
    std::cout << "Number of  edges: " << uer.getNumberOfEdges() << '\n';
    std::cout << "Number of  faces: " << uer.getNumberOfFaces() << '\n';

    for (size_t pointId = 0; pointId < uer.getNumberOfPoints(); ++pointId) {
        const double* point = uer.getPoint(pointId);
        std::cout << "pointId=" << pointId << " vertex: " << point[0] << ',' << point[1] << ',' << point[2] << '\n';
    }

    for (size_t edgeId = 0; edgeId < uer.getNumberOfEdges(); ++edgeId) {
        std::vector< Vector<double> > points = uer.getEdgePoints(edgeId);
        std::cout << "edgeId=" << edgeId << " vertices: " << points[0] << ',' << points[1] << '\n';
    }

    for (size_t faceId = 0; faceId < uer.getNumberOfFaces(); ++faceId) {
        std::cout << "faceId=" << faceId << " vertices: ";
        for (auto point : uer.getFacePoints(faceId)) {
            std::cout << point << ", ";
        }
        std::cout << '\n';
    }

    double xmin[3], xmax[3];
    uer.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';
}

void test16() {

    UgridReader uer;
    uer.load("${CMAKE_SOURCE_DIR}/data/cs_16.nc");
    std::cout << "Number of points: " << uer.getNumberOfPoints() << '\n';
    std::cout << "Number of  edges: " << uer.getNumberOfEdges() << '\n';
    std::cout << "Number of  faces: " << uer.getNumberOfFaces() << '\n';
    double xmin[3], xmax[3];
    uer.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';

    // locator
    FacesLocator loc;

    // test containsPoint
    vtkIdType faceId = 1;
    double tol = 1.e-6;
    std::vector< Vector<double> > verts = uer.getFacePointsRegularized(faceId);

    Vector<double> pMid = 0.25*(verts[0] + verts[1] + verts[2] + verts[3]);
    assert(uer.containsPoint(faceId, &pMid[0], tol));

    Vector<double> dp = verts[0] - pMid;

    Vector<double> p = pMid + 0.9*dp;
    assert(uer.containsPoint(faceId, &p[0], tol));
    std::cout << "p = " << p << " is in face: " << uer.containsPoint(faceId, &p[0], tol) << '\n';

    p = pMid + 0.99999*dp;
    assert(uer.containsPoint(faceId, &p[0], tol));
    std::cout << "p = " << p << " is in face: " << uer.containsPoint(faceId, &p[0], tol) << '\n';

    p = pMid + 1.00001*dp;
    assert(!uer.containsPoint(faceId, &p[0], tol));
    std::cout << "p = " << p << " is in face: " << uer.containsPoint(faceId, &p[0], tol) << '\n';

    p = pMid + 1.2*dp;
    assert(!uer.containsPoint(faceId, &p[0], tol));
    std::cout << "p = " << p << " is in face: " << uer.containsPoint(faceId, &p[0], tol) << '\n';

}

int main(int argc, char** argv) {

	test4();
    test16();

    return 0;
}

#include <mntUgridReader.h>
#include <mntFacesLocator.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test4() {

    UgridReader ur;
    ur.load("${CMAKE_SOURCE_DIR}/data/cs_4.nc");
    std::cout << "Number of faces: " << ur.getNumberOfFaces() << '\n';
    double xmin[3], xmax[3];
    ur.getRange(xmin, xmax);
    std::cout << "Domain range: " << xmin[0] << ',' << xmin[1] << ',' << xmin[2] << " -> "
                                  << xmax[0] << ',' << xmax[1] << ',' << xmax[2] << '\n';

    FacesLocator fl;
    int numFacesPerBucket = 1;
    fl.build(ur, numFacesPerBucket);

    double p0[3], p1[3];

    size_t ncases = 3;
    // start/end points
    double pab[] = {-50., -90., 0.,  -50., +90., 0.,    // line is outside
                    -90.,   0., 0.,  -90.,   0., 0.,    // line has zero length
                      0., -90., 0.,  360., +90., 0.};   // line traverses domain

    for (size_t icase = 0; icase < ncases; ++icase) {

        double* pa = &pab[icase*6 + 0];
        double* pb = &pab[icase*6 + 3];

        std::set<vtkIdType> faceIds = fl.getFacesAlongLine(pa, pb);

        std::cout << "Found " << faceIds.size() << " faces that likely intersect line "
              << pa[0] << ',' << pa[1] << ',' << pa[2] << " -> "
              << pb[0] << ',' << pb[1] << ',' << pb[2] << '\n';

        for (auto ift = faceIds.begin(); ift != faceIds.end(); ++ift) {
            std::vector< Vector<double> > points = ur.getFacePointsRegularized(*ift);
            std::cout << "face Id " << *ift << " points ";
            for (size_t j = 0; j < points.size(); ++j) {
                std::cout << points[j] << " ; ";
            }
            std::cout << '\n';
        }
    }
}

int main(int argc, char** argv) {

	test4();

    return 0;
}

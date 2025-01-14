import numpy
import vtk
from line_line_intersector import LineLineIntersector


class PolysegmentIter:

    def __init__(self, grid, locator, p0, p1):
        """
        Constructor
        @param grid instance of vtkUnstructuredGrid
        @param locator vtkCellLocator instance attached to the above grid
        @param p0 start point
        @param p1 end point
        """

        # small tolerances 
        self.eps = 1.73654365e-14
        self.eps100 = 100. * self.eps
        self.tol = 1.e-3 # to determine if a point is inside a cell


        self.grid = grid
        self.locator = locator

        self.cellIds = []
        self.xis = []
        self.ts = []
        self.__collectLineGridSegments(p0, p1)

        c2Inds = {}
        for i in range(len(self.cellIds)):
            cId = self.cellIds[i]
            c2Inds[cId] = c2Inds.get(cId, []) + [i]

        self.segCellIds = []
        self.segTas = []
        self.segTbs = []
        self.segXias = []
        self.segXibs = []
        self.segCoeffs = []

        for cId, inds in c2Inds.items():
            # sort by ts values
            indSorted = sorted(inds, key=lambda i: self.ts[i])
            n = len(indSorted)
            for i in range(n - 1):
                i0 = indSorted[i]
                i1 = indSorted[i + 1]
                ta, xia = self.ts[i0], self.xis[i0]
                tb, xib = self.ts[i1], self.xis[i1]
                self.segCellIds.append(cId)
                self.segTas.append(ta)
                self.segTbs.append(tb)
                self.segXias.append(xia)
                self.segXibs.append(xib)
                self.segCoeffs.append(1.0)

        n = len(self.segCellIds)
        inds = [i for i in range(n)]
        # sort by ta values
        indSorted = sorted(inds, key=lambda i: self.segTas[i])
        self.segCellIds = [self.segCellIds[i] for i in indSorted]
        self.segTas = [self.segTas[i] for i in indSorted]
        self.segTbs = [self.segTbs[i] for i in indSorted]
        self.segXias = [self.segXias[i] for i in indSorted]
        self.segXibs = [self.segXibs[i] for i in indSorted]

        # assign coefficients that account for duplicity, ie segments 
        # that are shared between two cells
        self.__assignCoefficientsToSegments()

        self.numSegs = len(self.segCellIds)

        self.totalT = 0
        for i in range(self.numSegs):
            ta, tb = self.segTas[i], self.segTbs[i]
            coeff = self.segCoeffs[i]
            self.totalT += (tb - ta) * coeff

        # reset the iterator
        self.reset()


    def getIntegratedParamCoord(self):
        """
        Get the integrated linear parametric coordinates
        @return value
        """
        return self.totalT


    def reset(self):
        """
        Reset the counter
        """
        self.index = -1


    def __iter__(self):
        return self


    def __next__(self):
        """
        Update iterator
        """
        if self.index < self.numSegs - 1:
            self.index += 1
            return self
        else:
            raise StopIteration()

    # Python2
    def next(self):
        return self.__next__()


    def getCellId(self):
        """
        Get the current cell Id
        @return index
        """
        return self.segCellIds[self.index]


    def getBegCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segXias[self.index]
        

    def getEndCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segXibs[self.index]
 

    def getBegLineParamCoord(self):
        """
        Get the current line parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segTas[self.index]
        

    def getEndLineParamCoord(self):
        """
        Get the current line parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segTbs[self.index]


    def getCoefficient(self):
        """
        Get the coefficient accounting for duplicates
        @return coefficient
        """
        return self.segCoeffs[self.index]
 

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return self.index


    def __assignCoefficientsToSegments(self):

        n = len(self.segCellIds)
        # copy
        sCellIds = self.segCellIds[:]
        sTas = self.segTas[:]
        sTbs = self.segTbs[:]
        sXias = self.segXias[:]
        sXibs = self.segXibs[:]
        sCoeffs = self.segCoeffs[:]
        # remove zero length segments
        self.segCellIds = []
        self.segTas = []
        self.segTbs = []
        self.segXias = []
        self.segXibs = []
        self.segCoeffs = []
        for i in range(n):
            ta, tb = sTas[i], sTbs[i]
            if abs(tb - ta) > self.eps100:
                self.segCellIds.append(sCellIds[i])
                self.segTas.append(sTas[i])
                self.segTbs.append(sTbs[i])
                self.segXias.append(sXias[i])
                self.segXibs.append(sXibs[i])
                self.segCoeffs.append(sCoeffs[i])

        # reduce contribution for overlapping segments. If two 
        # segments overlap then the coefficient of first segment
        # is set to 1.0 - overlap/(tb - ta). Assumes overlap 
        # can only happen for pairs of segment
        n = len(self.segCellIds)
        for i0 in range(n - 1):
            i1 = i0 + 1
            ta0, tb0 = self.segTas[i0], self.segTbs[i0]
            ta1, tb1 = self.segTas[i1], self.segTbs[i1]
            overlap = max(0., min(tb0, tb1) - max(ta1, ta0))
            self.segCoeffs[i0] = 1.0 - overlap/(tb0 - ta0)


    def __collectIntersectionPoints(self, pBeg, pEnd):
        """
        Collect all the intersection points
        @param pBeg starting point
        @param pEnd end point
        @return [(cellId, lambda, point), ...]
        @note lambda is the linear parametric coordinate along the line
        """

        res = []

        intersector = LineLineIntersector()
        cellIds = vtk.vtkIdList()
        ptIds = vtk.vtkIdList()

        dp = pEnd - pBeg

        # find all the cells intersected by the line
        self.locator.FindCellsAlongLine(pBeg, pEnd, self.tol, cellIds)

        # collect the intersection points in between
        for i in range(cellIds.GetNumberOfIds()):

            cId = cellIds.GetId(i)

            self.grid.GetCellPoints(cId, ptIds)

            # iterate over the quads' edges
            for j0 in range(4):

                j1 = (j0 + 1) % 4
                v0 = numpy.array(self.grid.GetPoint(ptIds.GetId(j0)))
                v1 = numpy.array(self.grid.GetPoint(ptIds.GetId(j1)))

                # look for an intersection
                intersector.setPoints(pBeg[:2], pEnd[:2], v0[:2], v1[:2])
                if not intersector.hasSolution(self.eps):
                    continue

                if abs(intersector.getDet()) > self.eps:
                    # normal intersection, 1 solution
                    lambRay, lambEdg = intersector.getSolution()

                    # is it valid? Intersection must be within (p0, p1) and (q0, q1)
                    if lambRay >= 0. - self.eps100 and lambRay <= 1. + self.eps100 and \
                        lambEdg >= 0. - self.eps100 and lambEdg <= 1. + self.eps100:

                        point = pBeg + lambRay*dp
                        res.append( (cId, lambRay, point) )

                else:
                    # det is almost zero
                    # looks like the two lines (p0, p1) and (q0, q1) are overlapping
                    # add the starting/ending points 
                    lama, lamb = intersector.getBegEndParamCoords()
                    pa = pBeg + lama*dp
                    pb = pBeg + lamb*dp
                    res.append( (cId, lama, pa) )
                    res.append( (cId, lamb, pb) )

        return res


    def __collectLineGridSegments(self, p0, p1):
        """
        Collect all the line-grid intersection points
        @param p0 starting point of the line
        @param p1 end point of the line 
        """

        # things we need to define
        ptIds = vtk.vtkIdList()
        cell = vtk.vtkGenericCell()
        cellIds = vtk.vtkIdList()
        subId = vtk.mutable(-1)
        dist = vtk.mutable(0.)
        xi = numpy.zeros((3,), numpy.float64)
        pBeg = numpy.zeros((3,), numpy.float64)
        pEnd = numpy.zeros((3,), numpy.float64)
        point = numpy.zeros((3,), numpy.float64)
        closestPoint = numpy.zeros((3,), numpy.float64)
        weights = numpy.array((4,), numpy.float64)

        # VTK wants 3d positions
        pBeg[:] = p0[0], p0[1], 0.0
        pEnd[:] = p1[0], p1[1], 0.0

        # add starting point
        cId = self.locator.FindCell(pBeg, self.eps, cell, xi, weights)
        if cId >= 0:
            self.cellIds.append(cId)
            self.xis.append(xi[:2].copy())
            self.ts.append(0.)
        else:
            pass
            #print('Warning: starting point {} not found!'.format(p0))

        #
        # find all intersection points in between
        #

        intersections = self.__collectIntersectionPoints(pBeg, pEnd)

        # find the cell id of the neighbouring cells
        for cId, lambRay, point in intersections:

            found = self.grid.GetCell(cId).EvaluatePosition(point, closestPoint, subId, xi, dist, weights)
            if found:
                self.cellIds.append(cId)
                self.xis.append(xi[:2].copy())
                self.ts.append(lambRay)
            else:
                print('Warning: param coord search failed point {} in cell {}'.format(point, cId))

            
        # add last point 
        cId = self.locator.FindCell(pEnd, self.eps, cell, xi, weights)
        if cId >= 0:
            self.cellIds.append(cId)
            self.xis.append(xi[:2].copy())
            self.ts.append(1.)
        else:
            pass
            #print('Warning: end point {} not found!'.format(p1))

#################################################################################################

def main():
    import argparse
    from math import pi
    from ugrid_reader import UgridReader
    from polyline_iter import PolylineIter

    parser = argparse.ArgumentParser(description='Break line into segments')
    parser.add_argument('-i', dest='input', default='mesh_C4.nc', help='Specify input file')
    parser.add_argument('-p', dest='points', default='(0., 0.),(2*pi,0.)', help='Points describing broken line as "(x0, y0),(x1, y1)..."')
    args = parser.parse_args()

    ur = UgridReader(filename=args.input)
    points = eval(args.points)
    pli = PolylineIter(points)
    
    countLine = 0
    tTotal = 0.
    for line in pli:
        t0 = line.getBegParamCoord()
        t1 = line.getEndParamCoord()
        p0 = line.getBegPoint()
        p1 = line.getEndPoint()
        print('line {} t = {} -> {}'.format(countLine, t0, t1))

        psi = PolysegmentIter(ur.getUnstructuredGrid(), 
                              ur.getUnstructuredGridCellLocator(), 
                              p0, p1)

        countSeg = 0
        for seg in psi:
            cellId = seg.getCellId()
            xia = seg.getBegCellParamCoord()
            xib = seg.getEndCellParamCoord()
            ta = seg.getBegLineParamCoord()
            tb = seg.getEndLineParamCoord()
            print('\tseg {} in cell {} t = {} -> {} xi = {} -> {}'.format(countSeg, cellId, ta, tb, xia, xib))

            countSeg += 1

        tTotal += psi.getIntegratedParamCoord()
        countLine += 1


    print('Integrated t = {}'.format(tTotal))

    

if __name__ == '__main__':
    main()

import numpy
import netCDF4
import argparse

parser = argparse.ArgumentParser(description='Generate uniform Ugrid file')
parser.add_argument('-nlon', type=int, default=1, help="Number of longitude cells")
parser.add_argument('-nlat', type=int, default=1, help="Number of latitude cells")
parser.add_argument('-m', type=str, default='physics', help="Mesh name")
parser.add_argument('-o', type=str, default='ugrid.nc', help="Output Ugrid file")
args = parser.parse_args()

numLon = args.nlon
numLat = args.nlat
numLon1 = numLon + 1
numLat1 = numLat + 1


meshName = args.m
numFaces = numLon * numLat
numEdges = numLon * numLat1 + numLon1 * numLat
numPoints = numLon1 * numLat1

# build the connectivity
face2point = numpy.zeros((numFaces, 4), numpy.int)
face2edge = numpy.zeros((numFaces, 4), numpy.int)
edge2point = numpy.zeros((numEdges, 2), numpy.int)

for j0 in range(numLat):
	j1 = j0 + 1
	for i0 in range(numLon):
		i1 = i0 + 1

		iface = i0 + j0*numLon

		p00 = i0 + j0*numLon1
		p10 = i1 + j0*numLon1
		p11 = i1 + j1*numLon1
		p01 = i0 + j1*numLon1

		face2point[iface, :] = p00, p10, p11, p01

		es0 = i0 + j0*numLon
		e1s = i1 + j0*numLon1 + numLon*numLat1
		es1 = i0 + j1*numLon
		e0s = i0 + j0*numLon1 + numLon*numLat1

		face2edge[iface, :] = es0, e1s, es1, e0s

# x edges
for j in range(numLat1):
	for i0 in range(numLon):
		i1 = i0 + 1
		iedge = i0 + j*numLon
		p0 = i0 + j*numLon1
		p1 = i1 + j*numLon1
		edge2point[iedge, :] = p0, p1

# y edges
for j0 in range(numLat):
	j1 = j0 + 1
	for i in range(numLon1):
		iedge = i + j0*numLon1 + numLon*numLat1
		p0 = i + j0*numLon1
		p1 = i + j1*numLon1
		edge2point[iedge, :] = p0, p1

# nodal points
lons = numpy.zeros((numLon1*numLat1), numpy.float64)
lats = numpy.zeros((numLon1*numLat1), numpy.float64)
dLon, dLat = 360./numLon, 180./numLat
for j in range(numLat1):
	for i in range(numLon1):
		k = i + j*numLon1
		lons[k] = 0.0 + i*dLon
		lats[k] = -90.0 + j*dLat

nc = netCDF4.Dataset(args.o, 'w', format="NETCDF4")

# create dimensions
faceDim = nc.createDimension("numFaces", numFaces)
edgeDim = nc.createDimension("numEdges", numEdges)
pointDim = nc.createDimension("numPoints", numPoints)
twoDim = nc.createDimension("two", 2)
fourDim = nc.createDimension("four", 4)

# create variables
mesh = nc.createVariable(args.m, 'int', [])
mesh.mesh_class = "sphere"
mesh.long_name = "Topology data of 2D unstructured mesh"
mesh.topology_dimension = 2 ;
mesh.node_coordinates = "{}_node_x {}_node_y".format(meshName, meshName)
mesh.face_node_connectivity = "{}_face_nodes".format(meshName)
mesh.edge_node_connectivity = "{}_edge_nodes".format(meshName)
mesh.face_edge_connectivity = "{}_face_edges".format(meshName)

face2pointVar = nc.createVariable("{}_face_nodes".format(meshName), "int", ("numFaces", "four"))
face2pointVar.cf_role = "face_node_connectivity"
face2pointVar.start_index = 0
face2pointVar[:] = face2point

face2edgeVar = nc.createVariable("{}_face_edges".format(meshName), "int", ("numFaces", "four"))
face2edgeVar.cf_role = "face_edge_connectivity"
face2edgeVar.start_index = 0
face2edgeVar[:] = face2edge

edge2pointVar = nc.createVariable("{}_edge_nodes".format(meshName), "int", ("numEdges", "two"))
edge2pointVar.cf_role = "edge_node_connectivity"
edge2pointVar.start_index = 0
edge2pointVar[:] = edge2point

xVar = nc.createVariable("{}_node_x".format(meshName), "f8", ("numPoints",))
xVar.standard_name = "longitude"
xVar.units = "degrees_east"
xVar[:] = lons

yVar = nc.createVariable("{}_node_y".format(meshName), "f8", ("numPoints",))
yVar.standard_name = "latitude"
yVar.units = "degrees_north"
yVar[:] = lats

# close
nc.close()




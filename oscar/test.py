import sys,os
sys.path.insert(0, os.getcwd()+"/../src/" )

from oscarH5 import oscarH5
  
#----- Start

h5file  = oscarH5( "yy.h5" )

h5file.setCycle(1)

print("Get coords of nodes 6,7 and 18",h5file.getCoords())
print("FFF" , h5file.getElemIndex(1),h5file.getNodeGroupNames())

h5file.saveAsVTU("zz")

'''
print(h5file.getElemNodes(17))

print(h5file.getElemNodeCount(7))

print(h5file.elemCount())

print(h5file.nodeCount())

print(h5file.nodeDataSets())

print(h5file.getDisplacements([0,4,7]))

print(h5file.getNodeData("S11",[0,4,7]))

print(h5file)

print(h5file.getNodeGroupNames())

print(h5file.getNodeGroup('x0'))

print(h5file.getNodeIndex(6))

print(h5file.getNodeIndex(65))

h5file.saveAsVTU()

h5file.saveAsDat(cycle=2)

h5file  = NDFFile( "silo_test.h5" )

h5file.setCycle(7)

print(h5file.particleCount())

h5file  = NDFFile( "yy.h5" )

h5file.setCycle(1)

print(h5file.nodeDataSets())

print(h5file.elemDataSets())

print(h5file.getElemData("S11",1))
print(h5file.getNodeData("S11",1))

#h5file.saveAsVTU()

'''




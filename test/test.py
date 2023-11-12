from oscar import oscar
  
#----- Start

h5file  = oscar( "test.h5" )

h5file.setCycle(1)

print("Get coords of nodes 6,7 and 18",h5file.getCoords())
print("FFF" , h5file.getElemIndex(1),h5file.getNodeGroupNames())

h5file.saveAsVTU("zz")

print(h5file.getElemNodes(2))

print(h5file.getElemNodeCount(2))

print(h5file.elemCount())

print(h5file.nodeCount())

print(h5file.nodeDataSets())

print(h5file.getDisplacements([0,4,7]))

print(h5file.getNodeData("S11",[0,4,7]))

print(h5file)

print(h5file.getNodeGroupNames())

print(h5file.getNodeGroup('bottom'))

print(h5file.getNodeIndex(6))

h5file.saveAsVTU()

h5file.saveAsDat(cycle=1)





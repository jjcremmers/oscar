"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""

from oscar import oscar
  
#----- Start

h5file  = oscar( "pinched8.h5" )

h5file.setCycle(9)

print(len(h5file.getCoords()))  #578

print(h5file.getCoords(37))      # [9.4705 2.8728 1.5675]

print(len(h5file.getCoords3()))  #578

print(h5file.getCycle())    # 9

print(h5file.getElemNodes(89))    # [ 94  95 112 111 383 384 401 400]

print(h5file.getElemNodeCount(89))    # 8

print(h5file.getElemGroupNames() )  #[]

print(h5file.getElemGroup())        # [0.255]

print(h5file.elemCount())    #256

print(h5file.getElemIndex(78))   #77

print(h5file.getElemIDs())   #1-256

print(h5file.nodeCount())   #578

print(h5file.getElemIndex(131))   #130

print(h5file.getNodeIDs())   #1-578

print(h5file.getNodeGroupNames())  #['loadx', 'loady', 'x0', 'y0', 'z0']

print(h5file.getNodeGroup( 'x0' ) )  #[ 16  33  50  67  84 101 118 135 152 169 186 203 220 237 254 271 288 305 322 339 356 373 390 407 424 441 458 475 492 509 526 543 560 577]

print(h5file.rank())   #3

print(h5file.nodeDataSets() ) # ['S11', 'S12', 'S13', 'S22', 'S23', 'S33', 'displacements']

print(h5file.elemDataSets() ) # ['S11', 'S12', 'S13', 'S22', 'S23', 'S33', 'displacements']

print(h5file.getDisplacements(56)) # [ 0.55154688 -0.14104671 -0.54406478]

print(h5file.getNodeData( "S11" , 78 ) )  # 9681.055732021865

#print(h5file.getMacroFieldData() )


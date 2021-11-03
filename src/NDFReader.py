"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021)
  
"""

import h5py,vtk
import numpy as np
    
def versionCheck( f ):

  if f.attrs['fileFormat'] != 'RNDF':
    print("Error")
    raise RuntimeError
  
  if f.attrs['version'] < 1.0:
    print("Error2")
    raise RuntimeError
    
class NDFFile():

  '''A class to read and access data in the NDF file format.'''  

  def __init__( self , name ):

    '''Inits the class NDFFile.
           
       Args: 
         name:   filename. May be with or without the extension .h5
    '''
    
    if not name.endswith('.h5'):
      name += '.h5'
           
    self.f = h5py.File( name, 'r')
    
    if 'cycleCount' in self.f.attrs.keys():  
      print("Single file with %d datasets" %self.f.attrs['cycleCount'] )
      self.cycle = -1
    else:
      self.data = self.f
      
    if self.f.attrs['version'] < 1.0:
      print("Error2")
      raise RuntimeError  
      
  def __str__( self ):
   
    if 'cycleCount' in self.f.attrs.keys():  
      return "Single file with %d datasets" %self.f.attrs['cycleCount']
  
  
    
  def setCycle( self , iCyc ):
  
    '''
    eee
    '''
    
    self.cycle = iCyc
    self.data = self.f["cycle"+str(self.cycle)]
    print("opened cycle %d." %self.cycle)
    
  def getCoords( self , nodeID = -1 ):
  
    '''
    eee
    '''
    
    grp  = self.data['nodes']
    dset = grp['coordinates']
  
    if nodeID == -1:
      return dset[:]
    else:
      return dset[nodeID]
      
  def getCycle( self ):
    
    '''
    
    '''
    
    return self.cycle  
      
  def getElemNodes( self , elemID ):
  
    '''
    Returns nodes of element ID
    '''
    
    grp   = self.data['elements']
    dset1 = grp['connectivity']
    dset2 = grp['offsets']
  
    if elemID == 0:
      return dset1[0:dset2[elemID]]
    else:
      return dset1[dset2[elemID-1]:dset2[elemID]]  

  def getElemNodeCount( self , elemID ):
  
    '''
    Returns the number of nodes in this element
    '''
    
    return len(self.getElemNodes(elemID))
    
  def getElemGroupNames( self ):
  
    '''
    
    '''
    
    return list(self.data['elementGroups'].keys())
    
  def getElemGroup( self , name ):
  
    '''
    
    '''
    
    grp = self.data['elementGroups']
    return grp[name][:]
            
  def elemCount( self , elemGroup = "all" ):
  
    '''
    Returns the number of elements in this data set
    '''
    
    if elemGroup == 'all':
      grp   = self.data['elements']
      return len(grp['offsets'])
    else:
      grp = self.data['elementGroups']
      return len(grp[elemGroup])


  def nodeCount( self , nodeGroup = 'all' ):
  
    '''
    Returns the number of nodes in this set
    '''
    
    if nodeGroup == 'all':
      grp = self.data['nodes']
      return grp['coordinates'].shape[0]
    else:
      grp = self.data['nodeGroups']
      return len(grp[nodeGroup])
      
  def getNodeGroupNames( self ):
  
    '''
    
    '''
    
    return list(self.data['nodeGroups'].keys())
    
  def getNodeGroup( self , name ):
  
    '''
    
    '''
    
    grp = self.data['nodeGroups']
    return grp[name][:]
                         
  def rank( self ):
  
    grp = self.data['nodes']
    return grp['coordinates'].shape[1]    

  def nodeDataSets( self ):
  
    '''
    Returns the labels of element datasets.
    '''
    grp = self.data['nodeData']
    return list(grp.keys())
    
    
  def elemDataSets( self ):
  
    '''
    Returns the label of element data
    '''
    
    grp = self.data['nodeData']
    return list(grp.keys())
      
  def getDisplacements( self , nodeID ):
  
    '''
    Returns the displacements of nodes nodeID (can be a list or an integer.
    '''       
    
    return self.getNodeData( 'displacements' , nodeID )
    
  def getNodeData( self , label , nodeID ):
  
    '''
    Returns the node data (label) of nodes nodeID (can be a list or an integer).
    '''
    
    grp  = self.data['nodeData']
    dset = grp[label]
    
    if dset.ndim == 1:
      return dset[nodeID]
    else:
      return dset[nodeID,:] 
      
  def getElemData( self , label , elemID ):
  
    '''
    Returns the element data (label) of elements elemID (can be a list or an integer).
    '''
    
    grp  = self.data['elemData']
    dset = grp[label]
    
    if dset.ndim == 1:
      return dset[elemID]
    else:
      return dset[elemID,:]      
      
      
      
      
      
      
      
      
      
      
      
          

def saveAsVTU( h5data , outputName ):
  
  writer = vtk.vtkXMLUnstructuredGridWriter()
  
  writer.SetDataModeToAscii()
  
  writer.SetFileName(outputName)
  
  grid = vtk.vtkUnstructuredGrid()
  
  # -- Read metadata
  
  cycle = h5data.attrs["cycle"]
  
  # -- Read element connectivity
  
  grp  = h5data['elements']
  dset = grp['connectivity']

  connectivity = dset[:]

  dset = grp['offsets']
  pointers = dset[:]
  
  nodeCount = []
     
  i0 = 0  
  for i1 in pointers:
    nodeCount.append(i1-i0)
    i0 = i1
    
  # -- Read nodes

  grp = h5data['nodes']
  dset = grp['coordinates']
  coordinates = dset[:]
  
  nElm = len(pointers)
  nNod = coordinates.shape[0]
  nRan = coordinates.shape[1]
  
  if coordinates.shape[1] == 2:
    newcrd = np.zeros(shape=(nNod,3))
    newcrd[:,:-1] = coordinates
    coordinates   = newcrd
  
  # -- Write nodes
    
  points = vtk.vtkPoints()
   
  for crd in coordinates:
    points.InsertNextPoint(crd)

  grid.SetPoints(points)    
    
  #--Store elements-----------------------------
     
  i0 = 0
  for i1 in pointers:
    nodeIDs = connectivity[i0:i1]
    
    if nRan == 2:   
      if len(nodeIDs) == 3:
        tmp = vtk.vtkTriangle()    
      elif len(nodeIDs) == 4:
        tmp = vtk.vtkQuad()
    elif nRan == 3:
      if len(nodeIDs) == 4:
        tmp = vtk.vtkTetra()
      elif len(nodeIDs) == 6:
        tmp = vtk.vtkWedge()
      elif len(nodeIDs) == 8:
        tmp = vtk.vtkHexahedron()
    else:
      print("Error")
      
    for i,inod in enumerate(nodeIDs):
      tmp.GetPointIds().SetId(i,inod)
    
    grid.InsertNextCell( tmp.GetCellType(),tmp.GetPointIds() );
      
    i0 = i1
          
  # -- Write nodedata

  grp = h5data['nodeData']
  
  for label in grp.keys():
  
    dset = grp[label]
    data = dset[:]
    
    if data.ndim == 2:
      if label == "displacements":
        if data.shape[1] == 2:
          newdata = np.zeros(shape=(nNod,3))
          newdata[:,:-1] = data
          data    = newdata
      
      d = vtk.vtkDoubleArray();
      d.SetName( label );
      d.SetNumberOfComponents(data.shape[1]);
            
      for i,line in enumerate(data):
        for j,l in enumerate(line):         
          d.InsertComponent( i , j , l )
    else:
      d = vtk.vtkDoubleArray();
      d.SetName( label );
      d.SetNumberOfComponents(1);
            
      for i,l in enumerate(data):        
          d.InsertComponent( i , 0 , l )
                   
    grid.GetPointData().AddArray( d )

  writer.SetInputData(grid)
  writer.Write()    
  
  print("  Cycle %6d written as VTU file. " %cycle )
      
      
      
       

  

  
#----- Start

h5file  = NDFFile( "pinched8.h5" )

h5file.setCycle(7)

print("Get coords of nodes 6,7 and 18",h5file.getCoords([6,7,18]))

print(h5file.getElemNodes(7))

print(h5file.getElemNodeCount(7))

print(h5file.elemCount())

print(h5file.nodeCount())

print(h5file.nodeDataSets())

print(h5file.getDisplacements([0,4,7]))

print(h5file.getNodeData("S11",[0,4,7]))

print(h5file)

print(h5file.getNodeGroupNames())

print(h5file.getNodeGroup('x0'))








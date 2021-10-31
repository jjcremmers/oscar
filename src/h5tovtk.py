import h5py,vtk
import numpy as np
    
def versionCheck( f ):

  if f.attrs['fileFormat'] != 'RNDF':
    print("Error")
    raise RuntimeError
  
  if f.attrs['version'] < 1.0:
    print("Error2")
    raise RuntimeError
    
class Hdf5File():

  def __init__( self , name ):
  
    self.data = h5py.File( name, 'r')
    
    versionCheck( self.data )
    
    print("YO")  
    
  def setCycle( self , iCyc ):
  
    self.data = self.f["cycle"+str(iCyc)]
    
  def getCoords( self , nodeID = -1 ):
  
    grp  = self.data['nodes']
    dset = grp['coordinates']
  
    if nodeID == -1:
      return dset[:]
    else:
      return dset[nodeID]

def saveAsVTU( h5data , outputName ):
  
  writer = vtk.vtkXMLUnstructuredGridWriter()
  
  writer.SetDataModeToAscii()
  
  writer.SetFileName(outputName)
  
  grid = vtk.vtkUnstructuredGrid()
  
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
  rank = coordinates.shape[1]
  
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
        
    if len(nodeIDs) == 4:
      tmp = vtk.vtkQuad()
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
      
      
      
       

  

  
#----- Start

h5file = h5py.File( "pinched8.h5" , 'r' )

for i,label in enumerate(h5file.keys()):
  print(i,label)
  saveAsVTU(h5file[label],str("output"+str(i)+".vtu"))



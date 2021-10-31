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

def h5tovtk( name ):

  inputName = name  + '.h5'
  outputName = name + '.vtu'
  
  f = h5py.File(inputName, 'r')
  
  versionCheck( f )
  
  writer = vtk.vtkXMLUnstructuredGridWriter()
  
  #writer.SetDataModeToAscii()
  writer.SetFileName(outputName)
  
  grid = vtk.vtkUnstructuredGrid()
  
  # -- Read element connectivity
  
  grp  = f['elements']
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

  grp = f['nodes']
  dset = grp['coordinates']
  coordinates = dset[:]
  
  nElm = len(pointers)
  nNod = coordinates.shape[0]
  rank = coordinates.shape[1]
  
  if coordinates.shape[1] == 2:
    newcrd = np.zeros(shape=(nNod,3))
    newcrd[:,:-1] = coordinates
    coordinates   = newcrd
  
  # -- Read data
  

        
  # -- Write nodedata
  
  grp  = f['nodeData']
  
  '''
  for label in grp.keys():
  
    dset = grp[label]
    data = dset[:]
    
    if label == "displacements":
      if data.shape[1] == 2:
        newdata = np.zeros(shape=(nNod,3))
        newdata[:,:-1] = data
        data    = newdata
      
    if data.ndim == 1:
      vtkfile.write('<DataArray type="Float64" Name="'+label+'" NumberOfComponents="1" format="ascii" >\n')
    else:
      vtkfile.write('<DataArray type="Float64" Name="'+label+'" NumberOfComponents="'+str(data.shape[1])+'" format="ascii" >\n')
    
    if data.ndim == 1:
      for d in data:
        vtkfile.write(str(d)+'\n')
    else:
      for d in data:    
        for x in d:
          vtkfile.write(str(x)+' ')
        vtkfile.write('\n')
 
    vtkfile.write('</DataArray>\n') 
     
  vtkfile.write('</PointData>\n')
  vtkfile.write('<CellData>\n')
  vtkfile.write('</CellData>\n')'''
  
  points = vtk.vtkPoints()
   
  for crd in coordinates:
    points.InsertNextPoint(crd)

  grid.SetPoints(points)
  
  #--Store elements-----------------------------
     
  i0 = 0
  for i1 in pointers:
    nodeIDs = connectivity[i0:i1]
    
    print(nodeIDs)
    
    if len(nodeIDs) == 4:
      tmp = vtk.vtkQuad()
  
      for i,inod in enumerate(nodeIDs):
        tmp.GetPointIds().SetId(i,inod)
    
      grid.InsertNextCell( tmp.GetCellType(),tmp.GetPointIds() );
      
    i0 = i1

  writer.SetInputData(grid)
  writer.Write()
  
  

        
h5tovtk( "PatchTest4_1" )

h5file = Hdf5File( "PatchTest4_1.h5" )

print(h5file.getCoords([1,3]))

f2 = h5py.File("test.h5", 'r')

print(f2.keys())

grp = f2['step2']
elements = grp['elements']

dset = elements['connectivity']
print(dset[:])

dset = elements['offsets']
print(dset[:])

nodes = grp['nodes']
dset = nodes['coordinates']
coordinates = dset[:]

print(coordinates)

nodeData = grp['nodeData']
dset = nodeData['displacements']

print(dset[:])

print(nodeData.keys())

dset = nodeData['S11']

print(dset[:])

f3 = h5py.File("pinched8.h5", 'r')

print(f3.attrs.keys(),str(f3.attrs["fileFormat"]),f3.attrs["version"])
print(f3.keys())


'''    

print(data)

grp = f['nodes']
dset = grp['coordinates']

data = dset[:]

print(data)

grp = f['nodeData']
dset = grp['displacements']

data = dset[:]

print(data)

dset = grp['S11']

data = dset[:]

print(data)

print(f.attrs.keys())

x = f.attrs['version']

print(x)

print(f['nodeData'].keys())
'''

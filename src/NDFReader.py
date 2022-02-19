"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021)
  
"""

import h5py,vtk
import numpy as np

from VTKutils import insertElement
    
def versionCheck( f ):

  if f.attrs['fileFormat'] != 'RNDF':
    print("Error")
    raise RuntimeError
  
  if f.attrs['version'] < 1.0:
    print("Error2")
    raise RuntimeError
    
class NDFFile():

  '''A class to read and access data in the NDF file format.'''  

  def __init__( self , fileName ):

    '''Inits the class NDFFile.
           
       Args: 
         name:   filename. May be with or without the extension .h5
    '''
    
    if not fileName.endswith('.h5'):
      fileName += '.h5'
      
    self.prefix = fileName.split('.')[0]       
      
    self.f = h5py.File( fileName, 'r')
    
    if 'cycleCount' in self.f.attrs.keys():  
      print("Single file with %d datasets" %self.f.attrs['cycleCount'] )
      self.cycle = -1
    else:
      self.data = self.f
      
    if self.f.attrs['version'] < 1.0:
      print("Error2")
      raise RuntimeError  
      
    self.cycleCount = self.f.attrs['cycleCount']      
      
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
      
  def getCoords3( self ):
  
    '''
    eee
    '''
    
    grp  = self.data['nodes']
    dset = grp['coordinates']
  
    if dset.shape[1] == 2:
      data = np.zeros(shape=(dset.shape[0],3))
      data[:,:2] = dset
      return data
    else:
      return dset[:]
      
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
      
  def getElemIndex( self , elemID ):
   
    '''
    
    '''
    
    
    grp = self.data['elements']
    return list(grp['elementIDs']).index(elemID)      

  def getElemIDs( self ):
   
    '''
    
    '''
        
    grp = self.data['elements']
    return grp['elementIDs'][:]

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
      
  def getNodeIndex( self , nodeID ):
   
    '''
    
    '''
    
    
    grp = self.data['nodes']
    return list(grp['nodeIDs']).index(nodeID)
    
  def getNodeIDs( self ):
   
    '''
    
    '''
        
    grp = self.data['nodes']#nodeIDs']
    return grp['nodeIDs'][:]
      
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
    Returns the labels of element data
    '''
    
    grp = self.data['elementData']
    return list(grp.keys())
      
  def getDisplacements( self , nodeID ):
  
    '''
    Returns the displacements of nodes nodeID (can be a list or an integer.
    '''       
    
    return self.getNodeData( 'displacements' , nodeID )
    
  def getNodeData( self , label , nodeID=-1 ):
  
    '''
    Returns the node data (label) of nodes nodeID (can be a list or an integer).
    '''
    
    grp  = self.data['nodeData']
    dset = grp[label]
    
    if nodeID == -1:
      return dset
    else:
      if dset.ndim == 1:
        return dset[nodeID]
      else:
        return dset[nodeID,:] 
      
  def getElemData( self , label , elemID=-1 ):
  
    '''
    Returns the element data (label) of elements elemID (can be a list or an integer).
    '''
    
    grp  = self.data['elementData']
    dset = grp[label]
    
    if elemID == -1:
      return dset
    else:
      if dset.ndim == 1:
        return dset[elemID]
      else:
        return dset[elemID,:]   
      
  def particleCount( self , particleGroup = 'all' ):
  
    '''
    Returns the number of nodes in this set
    '''
    
    if particleGroup == 'all':
      grp = self.data['particles']
      return grp['position'].shape[0]
    else:
      grp = self.data['particleGroups']
      return len(grp[particleGroup])        

#----------------
#
#---------------
      
  def saveAsVTU( self , prefix = 'None' , cycles = -1 ):
  
    '''
    Saves the data as VTU format
    '''
    
    if prefix == 'None':
      prefix = self.prefix
      
    vtufiles = []
    
    if cycles == -1:
      cycles = np.arange(1,self.cycleCount+1)
      
    for iCyc in cycles:
      self.setCycle( iCyc )

      writer = vtk.vtkXMLUnstructuredGridWriter()
  
      #writer.SetDataModeToAscii()
  
      vtufile = prefix+'-'+str(iCyc)+".vtu"
      
      writer.SetFileName(vtufile)
      
      vtufiles.append(vtufile)
  
      grid = vtk.vtkUnstructuredGrid()
  
      # -- Read metadata
    
      points = vtk.vtkPoints()
   
      coordinates = self.getCoords3()
      
      for crd in coordinates:
        points.InsertNextPoint(crd)

      grid.SetPoints(points)    
    
      #--Store elements-----------------------------
     
      for iElm in range(self.elemCount()):        
        insertElement( grid , self.getElemNodes(iElm) , 3 , 0 )
              
      # -- Write nodedata
  
      labels = self.nodeDataSets()
      
      for label in labels:
        data = self.getNodeData( label )
     
        if data.ndim == 2:
          if label == "displacements":
            if data.shape[1] == 2:
              newdata = np.zeros(shape=(data.shape[0],3))
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
        
      # -- Write elemdata
  
      labels = []#self.elemDataSets()
      
      for label in labels:
        data = self.getElemData( label )
             
        d = vtk.vtkDoubleArray();
        d.SetName( label );
        d.SetNumberOfComponents(1);
            
        for i,l in enumerate(data):        
          d.InsertComponent( i , 0 , l )
                   
        grid.GetCellData().AddArray( d )        

      writer.SetInputData(grid)
      writer.Write()                  
           
    pvdfile = open(prefix+".pvd", 'w')

    pvdfile.write("<VTKFile byte_order='LittleEndian' ")
    pvdfile.write("type='Collection' version='0.1'>\n")
    pvdfile.write("  <Collection>\n")
      
    for iCyc,vtufile in zip(cycles,vtufiles):
      pvdfile.write("    <DataSet file='"+vtufile+"' ")
      pvdfile.write("groups='' part='0' timestep='"+str(iCyc)+"'/>\n")

    pvdfile.write("  </Collection>\n")
    pvdfile.write("</VTKFile>\n")
    
  def saveAsDat( self , fileName = 'None' , cycle = -1 , output = "dawn"):
  
    '''
    Saves the data as dawn data format
    '''
    
    if fileName == 'None':
      fileName = self.prefix + '.dat'
        
    deformed = True    
        
    if cycle == -1:
      deformed = False
      cycle = 1
    
    self.setCycle( cycle )

    datFile = open(fileName, 'w')
       
    coordinates = self.getCoords()
    
    datFile.write("<Nodes>\n")
        
    if deformed:
      coordinates = coordinates + self.getNodeData( "displacements" )
               
    for nodeID,crd in enumerate(coordinates):
      datFile.write("  %d" %nodeID )
      for x in crd:
        datFile.write("  %e" %x )
      datFile.write(" ;\n")
      
    datFile.write("</Nodes>\n\n")      
           
    datFile.write("<Elements>\n\n")
    
    if output == 'dawn':       
      for elemID in np.arange(self.elemCount()):     
        datFile.write("  %d" %elemID )
        elemNodes = self.getElemNodes(elemID)
      
        for iNod in elemNodes:
          datFile.write("  %d" %iNod )
        datFile.write(" ;\n")       
    elif output == 'pyfem':
    
      elemGroups = self.getElemGroupNames()
  
      k = 0
      
      for grp in elemGroups:
        elemIDs = self.getElemGroup( grp )      
        for elemID in elemIDs:
          datFile.write("  %d %s" %(k,grp) )
          elemNodes = self.getElemNodes(elemID)
      
          for iNod in elemNodes:
            datFile.write("  %d" %iNod )
          datFile.write(" ;\n")      
          
          k = k+1 
                  
    datFile.write("</Elements>\n")      
        
    nodeGroups = self.getNodeGroupNames()
  
    for grp in nodeGroups:
      datFile.write("<NodeGroup name = '%s'>\n  {" %grp) 
      nodeIDs = self.getNodeGroup( grp )
      
      for k,nodeID in enumerate(nodeIDs):
        datFile.write(" %d" %nodeID) 
        if (k+1)%10 == 0:
          datFile.write("\n")
                 
      datFile.write(" }\n</NodeGroup>\n\n") 
      
    elemGroups = self.getElemGroupNames()
  
    for grp in elemGroups:
      datFile.write("<ElementGroup name = '%s'>\n  {" %grp) 
      elemIDs = self.getElemGroup( grp )
      
      for k,elemID in enumerate(elemIDs):
        datFile.write(" %d" %elemID) 
        if (k+1)%10 == 0:
          datFile.write("\n")
        
      datFile.write(" }\n</ElementGroup>\n\n")       
  
   

'''
              
      # -- Write nodedata
  
      labels = self.nodeDataSets()
      
      for label in labels:
        data = self.getNodeData( label )
     
        if data.ndim == 2:
          if label == "displacements":
            if data.shape[1] == 2:
              newdata = np.zeros(shape=(data.shape[0],3))
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
        
      # -- Write elemdata
  
      labels = []#self.elemDataSets()
      
      for label in labels:
        data = self.getElemData( label )
             
        d = vtk.vtkDoubleArray();
        d.SetName( label );
        d.SetNumberOfComponents(1);
            
        for i,l in enumerate(data):        
          d.InsertComponent( i , 0 , l )
                   
        grid.GetCellData().AddArray( d )        

      writer.SetInputData(grid)
      writer.Write()                  
           
    pvdfile = open('pinched8.pvd', 'w')

    pvdfile.write("<VTKFile byte_order='LittleEndian' ")
    pvdfile.write("type='Collection' version='0.1'>\n")
    pvdfile.write("  <Collection>\n")
      
    for iCyc,vtufile in zip(cycles,vtufiles):
      pvdfile.write("    <DataSet file='"+vtufile+"' ")
      pvdfile.write("groups='' part='0' timestep='"+str(iCyc)+"'/>\n")

    pvdfile.write("  </Collection>\n")
    pvdfile.write("</VTKFile>\n")    
    
'''      
      
'''  
          

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
      
      
'''      
       

  

  
#----- Start

h5file  = NDFFile( "pinched8.h5" )

h5file.setCycle(7)

print("Get coords of nodes 6,7 and 18",h5file.getCoords([6,7,18]))

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

#h5file.saveAsVTU()

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






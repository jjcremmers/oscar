"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2022)
  
"""

import h5py,vtk
import numpy as np

from VTKutils import insertElement
    
def versionCheck( f ):

  '''
  Checks if the file is an oscar file
  '''
  
  if f.attrs['fileFormat'] != 'RNDF':
    print("Error")
    raise RuntimeError
  
  if f.attrs['version'] < 1.0:
    print("Error2")
    raise RuntimeError

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
def writePVD( prefix , cycles , vtuFiles ):
    
  if len(cycles) != len(vtuFiles):
    print("writePVD: cycles and vtuFiles must have the same length")
    raise RunTimeError
    
  pvdfile = open(prefix+".pvd", 'w')

  pvdfile.write("<VTKFile byte_order='LittleEndian' ")
  pvdfile.write("type='Collection' version='0.1'>\n")
  pvdfile.write("  <Collection>\n")
      
  for iCyc,vtufile in zip(cycles,vtuFiles):
    pvdfile.write("    <DataSet file='"+vtufile+"' ")
    pvdfile.write("groups='' part='0' timestep='"+str(iCyc)+"'/>\n")

  pvdfile.write("  </Collection>\n")
  pvdfile.write("</VTKFile>\n")
            
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class oscarH5():

  '''
    A class to read and access data in the NDF file format.
  '''  

  def __init__( self , fileName ):

    '''Inits the class oscarH5
           
       Args: 
         name:   filename. May be with or without the extension .h5
    '''
    
    if not fileName.endswith('.h5'):
      fileName += '.h5'
      
    self.prefix = fileName.split('.')[0]       
      
    self.f = h5py.File( fileName, 'r')
    
    if 'cycleCount' in self.f.attrs.keys():  
      print("Single file with %d datasets (cycles)" %self.f.attrs['cycleCount'] )
      self.cycle = -1
    else:
      self.data = self.f
      
    if self.f.attrs['version'] < 1.0:
      print("Error2")
      raise RuntimeError  
      
    self.cycleCount = self.f.attrs['cycleCount']      

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def __str__( self ):
  
    '''
    Prints the main contents of the file
    '''
   
    if 'cycleCount' in self.f.attrs.keys():  
      return "Single file with %d datasets (cycles)." %self.f.attrs['cycleCount']
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def setCycle( self , iCyc ):
  
    '''
      Set the cycle for whcih output is read.  
    
      Args:
        iCyc:   Cycle ID number.  
    '''
    
    self.cycle = iCyc
    self.data = self.f["cycle"+str(self.cycle)]
    
    print("opened cycle %d." %self.cycle)
    
#-------------------------------------------------------------------------------
#  getCoords
#-------------------------------------------------------------------------------
    
  def getCoords( self , nodeID = -1 ):
  
    '''
      Gets the coordinates of node(s) nodeID
      
      Args:
        nodeID:   nodeID (list or integrer) according to internal numbering
                  if omitted, all nodes are printed.
    '''
    
    grp  = self.data['nodes']
    dset = grp['coordinates']
  
    if nodeID == -1:
      return dset[:]
    else:
      return dset[nodeID]

#-------------------------------------------------------------------------------
#  getCoords3
#-------------------------------------------------------------------------------
      
  def getCoords3( self ):
  
    '''
      Returns the coordinates always in a 3D format (even if it is 2D).
    '''
    
    grp  = self.data['nodes']
    dset = grp['coordinates']
  
    if dset.shape[1] == 2:
      data = np.zeros(shape=(dset.shape[0],3))
      data[:,:2] = dset
      return data
    else:
      return dset[:]

#-------------------------------------------------------------------------------
#  getCycle
#-------------------------------------------------------------------------------
      
  def getCycle( self ):
    
    '''
      Returns the current cycle ID
    '''
    
    return self.cycle  
      
#-------------------------------------------------------------------------------
#  getElemNodes
#-------------------------------------------------------------------------------
      
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

#-------------------------------------------------------------------------------
#  getElemNodeCount
#-------------------------------------------------------------------------------

  def getElemNodeCount( self , elemID ):
  
    '''
    Returns the number of nodes in this element
    '''
    
    return len(self.getElemNodes(elemID))

#-------------------------------------------------------------------------------
#  getElemGroupNames
#-------------------------------------------------------------------------------
    
  def getElemGroupNames( self ):
  
    '''
      Returns all element group names
    '''
    
    return list(self.data['elementGroups'].keys())

#-------------------------------------------------------------------------------
#  getElemGroup
#-------------------------------------------------------------------------------
    
  def getElemGroup( self , name = "all" ):
  
    '''
      Returns the element ID in a
    '''
    
    if name == "all":
      return range(self.elemCount("all"))
    else:
      grp = self.data['elementGroups']
      return grp[name][:]

#-------------------------------------------------------------------------------
#  elemCount
#-------------------------------------------------------------------------------
            
  def elemCount( self , elemGroup = "all" ):
  
    '''
      Returns the number of elements in this data set (in this cycle).
      
      Args:
        elemGroup:  specify an element group. In that case it returns the 
                    number of elements in that group.      
    '''
    
    if elemGroup == 'all':
      grp   = self.data['elements']
      return len(grp['offsets'])
    else:
      grp = self.data['elementGroups']
      return len(grp[elemGroup])

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def getElemIndex( self , elemID ):
   
    '''
      Return element Index
      
      Args:
        elemID:   integer or list of elemID
    '''
    
    
    grp = self.data['elements']
    return list(grp['elementIDs']).index(elemID)      

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getElemIDs( self ):
   
    '''
      Return element IDs
    '''
        
    grp = self.data['elements']
    return grp['elementIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def getNodeIndex( self , nodeID ):
   
    '''
    
    '''
    
    
    grp = self.data['nodes']
    return list(grp['nodeIDs']).index(nodeID)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getNodeIDs( self ):
   
    '''
    
    '''
        
    grp = self.data['nodes']#nodeIDs']
    return grp['nodeIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def getNodeGroupNames( self ):
  
    '''
      Returns the names of all NodeGroups in this data set.
    '''
    
    return list(self.data['nodeGroups'].keys())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getNodeGroup( self , name ):
  
    '''
    
    '''
    
    grp = self.data['nodeGroups']
    return grp[name][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
                         
  def rank( self ):
  
    grp = self.data['nodes']
    return grp['coordinates'].shape[1]    

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def nodeDataSets( self ):
  
    '''
    Returns the labels of element datasets.
    '''
    grp = self.data['nodeData']
    return list(grp.keys())
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def elemDataSets( self ):
  
    '''
    Returns the labels of element data
    '''
    
    grp = self.data['elementData']
    return list(grp.keys())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def getDisplacements( self , nodeID ):
  
    '''
    Returns the displacements of nodes nodeID (can be a list or an integer.
    '''       
    
    return self.getNodeData( 'displacements' , nodeID )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getNodeData( self , label , nodeID=-1 ):
  
    '''
      Returns the node data (label) of nodes nodeID (can be a list or an integer).
    '''

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    grp  = self.data['nodeData']
    dset = grp[label]
    
    if nodeID == -1:
      return dset
    else:
      if dset.ndim == 1:
        return dset[nodeID]
      else:
        return dset[nodeID,:]
         
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------      

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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def saveAsVTU( self , prefix = 'None' , cycles = -1 ):
  
    '''
      Saves the data as VTU format
      
      Args:
      
        prefix: prefix of the pvd and vtu filenames.
        cycle:  a list (or integer) of cycles that need
                to be written.
    '''
    
    if prefix == 'None':
      prefix = self.prefix
      
    vtufiles = []
    
    if cycles == -1:
      cycles = np.arange(1,self.cycleCount+1)
      
    for iCyc in cycles:
      self.setCycle( iCyc )

      writer = vtk.vtkXMLUnstructuredGridWriter()
  
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
       
    writePVD( prefix , cycles , vtufiles )
        
#-------------------------------------------------------------------------------
#  saveAsDat
#-------------------------------------------------------------------------------
    
  def saveAsDat( self , fileName = 'None' , cycle = -1 , output = "dawn"):
  
    '''
      Saves the data as dawn data format
      
      Args:
      
        fileName:  output filename. If not ending by .dat, it will be added.
        cycle:     cycle number.
        output:    ouput type.
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
    else:
      print("Error") 
                  
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
    
    if output == 'dawn':    
      elemGroups = self.getElemGroupNames()
  
      for grp in elemGroups:
        datFile.write("<ElementGroup name = '%s'>\n  {" %grp) 
        elemIDs = self.getElemGroup( grp )
      
        for k,elemID in enumerate(elemIDs):
          datFile.write(" %d" %elemID) 
          if (k+1)%10 == 0:
            datFile.write("\n")
        
        datFile.write(" }\n</ElementGroup>\n\n")       
        
    datFile.close()
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def saveModes( self , prefix = 'None' ):
  
    '''
      Saves the data as VTU format
    '''
    
    if prefix == 'None':
      prefix = self.prefix
      
    vtufiles = []
    
    if cycles == -1:
      cycles = np.arange(1,self.cycleCount+1)
      
    for iMod in range(5):

      writer = vtk.vtkXMLUnstructuredGridWriter()
  
      vtufile = prefix+'-'+str(iMod+1)+".vtu"
      
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
      '''
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
      '''
           
    pvdfile = open(prefix+".pvd", 'w')

    pvdfile.write("<VTKFile byte_order='LittleEndian' ")
    pvdfile.write("type='Collection' version='0.1'>\n")
    pvdfile.write("  <Collection>\n")
      
    '''cycles = arange(5):
    
    for iCyc,vtufile in zip(cycles,vtufiles):
      pvdfile.write("    <DataSet file='"+vtufile+"' ")
      pvdfile.write("groups='' part='0' timestep='"+str(iCyc+1)+"'/>\n")

    pvdfile.write("  </Collection>\n")
    pvdfile.write("</VTKFile>\n")   
'''
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
      




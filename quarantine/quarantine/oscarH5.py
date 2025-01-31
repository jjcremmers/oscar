"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2022)
  
"""

import h5py,vtk
import numpy as np

import vtk


def setCellNodes( cell , elemNodes ):
        
  for i,inod in enumerate(elemNodes):
    cell.GetPointIds().SetId(i,inod)
          

def insertElement( grid , elemNodes , rank , family ):

  '''
  Inserts an element 
  '''
  
  nNod = len(elemNodes)
  
  if family == 0:   # Continuum
    if rank == 2:
      if nNod == 3:
        cell = vtk.vtkTriangle()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )     
      elif nNod == 4:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )       
      elif nNod == 6:
        cell = vtk.vtkTriangle()      
        setCellNodes( cell , elemNodes[0:6:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
      elif nNod == 8:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes[0:8:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
      elif nNod == 9:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes[0:8:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
      else:
        raise NotImplementedError('Only 3, 4, 6, 8, 9 node continuum elements in 2D.')     
    elif rank == 3:
      if nNod == 4:
        cell = vtk.vtkTetra()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
      elif nNod == 5:
        cell = vtk.vtkPyramid()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
      elif nNod == 6:
        cell = vtk.vtkWedge()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
      elif nNod == 8:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )               
      elif nNod == 16:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , numpy.concatenate(elemNodes[0:8:2],elemNodes[8:16:2] ) ) 
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                         
      else:
        raise NotImplementedError('Only 4, 5, 6, 8 and 16 node continuum elements in 3D.')             
    else:
      raise NotImplementedError('Only 2D and 3D continuum elements.')
  elif family == 1:
    if rank == 2:
      if nNod == 4:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes[0:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[2:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                             
      else:
        raise NotImplementedError('Only 4 node interface elements in 2D.')           
    elif rank == 3:
      if nNod == 6:
        cell = vtk.vtkTria() 
        setCellNodes( cell , elemNodes[0:3] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[3:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
      elif nNod == 8:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes[0:4] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[4:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
      else:
        raise NotImplementedError('Only 6 and 8 node interface elements in 3D.')        
    else:
      raise NotImplementedError('Only 2D and 3D interface elements.')        
  elif family == 2:
    if rank == 2:
      if nNod == 2:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
      else:
        raise NotImplementedError('Only 2 node surface elements in 2D.')              
    elif rank == 3:
      if nNod == 3:
        cell = vtk.vtkTriangle() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
      elif nNod == 4:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
      else:
        raise NotImplementedError('Only 3 and 4 node surface elements in 3D.')              
  elif family == 3:
    if nNod == 2:
      cell = vtk.vtkLine() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )    
    elif nNod == 3:
      cell = vtk.vtkLine() 
      cell.GetPointIds().SetId(0,elemNodes[0])
      cell.GetPointIds().SetId(1,elemNodes[2])      
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    else:
      raise NotImplementedError('Only 2 and 3 node beam elements.')            
  elif family == 4:
    if nNod == 3:
      cell = vtk.vtkTriangle() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
    elif nNod == 4:
      cell = vtk.vtkQuad() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                     
  else:
    raise NotImplementedError('Only 3 and 4 node shell elements.')       

from .VTKutils import insertElement
    
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
    A class to read and access data in the HDF file format.
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
    Set the cycle for which the output is read.  
    
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
    Returns the coordinates always in a 3D format (even if the rank is equal to 2).
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
    
    Args:
    
      elemID:    elementID number.
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
    
    Args:
    
      elemID:    element ID number. This is an integer.
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
    Returns the element IDs in an elementset.
      
    Args:
    
      name:    Name of the elementset. By defauls all element IDs are given.  
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
        elemID:   integer or list of elementIDs
    '''
    
    
    grp = self.data['elements']
    return list(grp['elementIDs']).index(elemID)      

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getElemIDs( self ):
   
    '''
      Return element IDs in the set.
    '''
        
    grp = self.data['elements']
    return grp['elementIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def nodeCount( self , nodeGroup = 'all' ):
  
    '''
    Returns the number of nodes in the nodegroup.
    
    Args:
    
      nodeGroup:    nodeGroup name. By defaults the total number
                    of nodes is given.
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
    Returns the index of nodeID. This is the index of the node in the input file.
    
    Args:
    
      nodeID:   The nodeID. This can be an integer or a list.
    '''
        
    grp = self.data['nodes']
    return list(grp['nodeIDs']).index(nodeID)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getNodeIDs( self ):
   
    '''
    Returns a list of all the nodeIDs.
    '''
        
    grp = self.data['nodes']
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
    Returns the nodeID that are present in a certain NodeGroup.
    
    Args:
      
      name:    Name of the nodegroup.       
    '''
    
    grp = self.data['nodeGroups']
    return grp[name][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
                         
  def rank( self ):
  
    '''
    Returns the rank of the problem. This is the number of spatial
      dimensions.
    '''
  
    grp = self.data['nodes']
    return grp['coordinates'].shape[1]    

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def nodeDataSets( self ):
  
    '''
    Returns all the labels of node datasets.
    '''
    
    grp = self.data['nodeData']
    return list(grp.keys())
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def elemDataSets( self ):
  
    '''
    Returns all the labels of element datasets.
    '''
    
    grp = self.data['elementData']
    return list(grp.keys())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
  def getDisplacements( self , nodeID ):
  
    '''
    Returns the displacements of nodes
    
    Args:
      
      nodeID:   node number. This can be an integer or a list. 
                By default the data of all nodes is given as an array.
    '''       
    
    return self.getNodeData( 'displacements' , nodeID )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getNodeData( self , label , nodeID=-1 ):
  
    '''
    Returns the node data (label) of nodes nodeID (can be a list or an integer).
      
    Args:
      
      label:    Name of the node dataset.
      nodeID:   node number. This can be an integer or a list. 
                By default the data of all nodes is given as an array.
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

  def getMacroFieldData( self , label ):
  
    '''
    Returns the element data (label) of elements elemID (can be a list or an integer).
    '''
    
    grp  = self.data['macrofield']
    return grp[label][:]
         
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
      
        prefix: prefix of the pvd and vtu filenames. If omitted the
                prefix of the original vtu file is used.
        cycle:  a list (or integer) of cycles that need
                to be written. If omitted, all files will be exported.
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
      




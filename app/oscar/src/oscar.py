"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""

import h5py,vtk
import numpy as np

import vtk
    
from .VTKutils import insertElement
 
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
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
#  convertToVTU
#-------------------------------------------------------------------------------

def convertToVTU( fileName , cycles = -1 ):

    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsVTU()
  
#-------------------------------------------------------------------------------
#  convertToDat
#-------------------------------------------------------------------------------

def convertToDat( fileName , cycles = -1 ):

    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsDat()
  
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

class oscar():

    '''
      A class to read and access data in the HDF file format.
    '''  

    def __init__( self , fileName , verbose = False ):

        '''
          Inits the class oscar
           
            Args: 
            name:   filename. May be with or without the extension .h5
        '''
    
        if not fileName.endswith('.h5'):
            fileName += '.h5'
      
        self.prefix = fileName.split('.')[0]       
      
        self.f = h5py.File( fileName, 'r')
    
        self.verbose = verbose
    
        if 'cycleCount' in self.f.attrs.keys():  
            if self.verbose:
                print("Single file with %d datasets (cycles)" %self.f.attrs['cycleCount'] )
            self.cycle = -1
            self.data = self.f["cycle1"]      
        else:
            self.cycle = 1
            self.data = self.f["cycle1"]
      
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
        
        if 'elements' in self.data.keys():
            self.connectivity = self.data['elements']['connectivity']
            self.offsets      = self.data['elements']['offsets']

        if self.verbose:
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
      
      dset = self.data['nodes']['coordinates']
        
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
  
        if elemID == 0:
            return self.connectivity[0:self.offsets[elemID]]
        else:
            return self.connectivity[self.offsets[elemID-1]:self.offsets[elemID]]  

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
            return list(range(self.elemCount("all")))
        else:
            return self.data['elementGroups'][name][:]

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
    
        if 'elementData' in self.data.keys():
            grp = self.data['elementData']
            return list(grp.keys())
        else:
            return []

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
        elif type(cycles) == int:
            cycles = [cycles]  
      
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
   
           
        pvdfile = open(prefix+".pvd", 'w')

        pvdfile.write("<VTKFile byte_order='LittleEndian' ")
        pvdfile.write("type='Collection' version='0.1'>\n")
        pvdfile.write("  <Collection>\n")

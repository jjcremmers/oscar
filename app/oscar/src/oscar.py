"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""

import h5py
import os
import shutil
import vtk
import numpy as np
    
from .VTKutils import insertElement

def apply_rainbow_color_map(mapper):

    rng = mapper.GetScalarRange()

    lookup_table = vtk.vtkLookupTable()
    lookup_table.SetNumberOfTableValues(256)
    lookup_table.Build()
    lookup_table.SetTableRange(rng[0],rng[1])  # Set the range of data values
    lookup_table.SetHueRange(0.0, 0.8)  # Map color hue from blue to red
    mapper.SetLookupTable(lookup_table)
 
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
    
    '''
    Converts the oscar .h5 file to a vtu file. The VTU file is written to the 
    drive as fileName.vtu
    
    Arguments:
        fileName     the name of the h5 file (string)
        cycles       a list of cycles that is converted.
                     If cycles is equal to -1, all cycles
                     are converted.
                     
    Returns:                     
        --
    '''
    
    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsVTU()
  
#-------------------------------------------------------------------------------
#  convertToDat
#-------------------------------------------------------------------------------

def convertToDat( fileName , cycles = -1 ):

    '''
    Converts the oscar .h5 file to a dawn dat file. The dat file is written 
    to the drive as fileName.dat
    
    Arguments:
        fileName     the name of the h5 file (string)
        cycles       a list of cycles that is converted.
                     If cycles is equal to -1, all cycles
                     are converted.
                     
    Returns:                     
        --
    '''
    
    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsDat()
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
def writePVD( prefix , cycles , nProc = 1 ):
      
    '''
    
    '''
    
    if nProc < 1:
        print("writePVD: number of processros should by 1 or more.")
        raise RunTimeError    
    
    pvdfile = open(prefix+".pvd", 'w')

    pvdfile.write("<VTKFile byte_order='LittleEndian' ")
    pvdfile.write("type='Collection' version='0.1'>\n")
    pvdfile.write("  <Collection>\n")
      
    for iCyc in cycles:
        if nProc == 1:
            pvdfile.write("    <DataSet file='"+PVDfileName( prefix , iCyc )+"' ")
            pvdfile.write("groups='' part='0' timestep='"+str(iCyc)+"'/>\n")
        else:
            for iProc in range(nProc):
                pvdfile.write("    <DataSet file='" + 
                                PVDfileName( prefix , iCyc , iProc )+"' ")
                pvdfile.write("groups='' part='0' timestep='"+str(iCyc)+"'/>\n")         
                     
    pvdfile.write("  </Collection>\n")
    pvdfile.write("</VTKFile>\n")
    
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


def PVDfileName( prefix : str , iCyc :int , iProc : int = -1 ) -> str:

    '''
    
    '''
    
    if iProc == -1:
        return prefix + "_t" + str(iCyc) + ".vtu"
    else:
        return prefix + "_p" + str(iProc) + "_t" + str(iCyc) + ".vtu"     
            
            
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


class oscar():

    '''
      A class to read and access data in the HDF file format.
      
          Inits the class oscar
           
            Args: 
            name:   filename. May be with or without the extension .h5
            
    '''  

    def __init__( self , fileName : str , verbose : bool = False ):

        '''
        Constructor
        '''
            
        if not fileName.endswith('.h5'):
            fileName += '.h5'
      
        tempPrefix = fileName.split('.h5')[0]
        self.prefix = tempPrefix.split('/')[-1]

        self.fileDir = fileName.split(self.prefix)[0]
      
        self.f = h5py.File( fileName, 'r')
    
        self.verbose = verbose
    
        if 'cycleCount' in self.f.attrs.keys():  
            if self.verbose:
                print("Single file with %d datasets (cycles)"
                       %self.f.attrs['cycleCount'] )
            self.cycle = -1
            self.data = self.f["cycle1"]      
        else:
            self.cycle = 1
            self.data = self.f["cycle1"]
      
        if self.f.attrs['version'] < 1.0:
            print("Error2")
            raise RuntimeError  
      
        self.cycleCount = self.f.attrs['cycleCount']   
        
        if 'elements' in self.f.keys():
            connectivity = self.f['elements']['connectivity'][:]
            offsets      = self.f['elements']['offsets'][:]
                    
            self.elemNodes = self.unpackElements( connectivity , offsets )
            
        if 'nodes' in self.f.keys():
            self.coordinates  = self.f['nodes']['coordinates']                  

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def __str__( self ) -> str:
  
        '''
         Prints the main contents of the file
        '''
   
        if 'cycleCount' in self.f.attrs.keys():  
            return "Single file with %d datasets (cycles)." %self.f.attrs['cycleCount']
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def setCycle( self , iCyc : int ):
  
        '''
         Set the cycle for which the output is read.  
    
         Args:
           iCyc:   Cycle ID number.  
        '''
    
        self.cycle = iCyc
        self.data = self.f["cycle"+str(self.cycle)]
        
        if 'elements' in self.data.keys():
            connectivity = self.data['elements']['connectivity'][:]
            offsets      = self.data['elements']['offsets'][:]
            
            self.elemNodes = self.unpackElements( connectivity , offsets )            
                        
        if 'nodes' in self.data.keys():
            self.coordinates  = self.data['nodes']['coordinates']                  

        if self.verbose:
            print("opened cycle %d." %self.cycle)
    
#-------------------------------------------------------------------------------
#  getCoords
#-------------------------------------------------------------------------------
    
    def getCoords( self , nodeID : int = -1 ):
  
        '''
          Gets the coordinates of node(s) nodeID
      
          Args:
            nodeID:   nodeID (list or integrer) according to internal numbering
                      if omitted, all nodes are printed.
        '''
      
        if nodeID == -1:
            return self.coordinates[:]
        else:
            return self.coordinates[nodeID]

#-------------------------------------------------------------------------------
#  getCoords3
#-------------------------------------------------------------------------------
      
    def getCoords3( self ):
  
      '''
      Returns the coordinates always in a 3D format (even if the rank is equal to 2).
      '''
              
      if self.coordinates.shape[1] == 2:
          data = np.zeros(shape=(dset.shape[0],3))
          data[:,:2] = self.coordinates
          return data
      else:
          return self.coordinates[:]

#-------------------------------------------------------------------------------
#  getCycle
#-------------------------------------------------------------------------------
      
    def getCycle( self ) -> int:
    
        '''
          Returns the current cycle ID
        '''
    
        return self.cycle  
      
#-------------------------------------------------------------------------------
#  getElemNodes
#-------------------------------------------------------------------------------
      
    def getElemNodes( self , elemID : int ) -> list:
  
        '''
        Returns nodes of element ID
    
        Args:
    
            elemID:    elementID number.
        '''
   
        '''
        if elemID == 0:
            return self.connectivity[0:self.offsets[elemID]]
        else:
            return self.connectivity[self.offsets[elemID-1]:self.offsets[elemID]]  
        '''
        
        return self.elemNodes[elemID]            

#-------------------------------------------------------------------------------
#  getElemNodeCount
#-------------------------------------------------------------------------------

    def getElemNodeCount( self , elemID : int ) -> int:
  
        '''
        Returns the number of nodes in this element
    
        Args:
    
            elemID:    element ID number. This is an integer.
        '''
    
        return len(self.elemNodes[elemID])

#-------------------------------------------------------------------------------
#  getElemGroupNames
#-------------------------------------------------------------------------------
    
    def getElemGroupNames( self ) -> list:
  
        '''
        Returns all element group names
        '''
        
        if 'elementGroups' in self.data.keys():    
            return list(self.data['elementGroups'].keys())
        elif 'elementGroups' in self.f.keys():
            return list(self.f['elementGroups'].keys())
        else:
            return []

#-------------------------------------------------------------------------------
#  getElemGroup
#-------------------------------------------------------------------------------
    
    def getElemGroup( self , name : str = "all" ) -> list[int]:
  
        '''
        Returns the element IDs in an elementset.
      
        Args:
    
            name:    Name of the elementset. By defauls all element IDs are given.  
        '''
    
        if name == "all":
            return list(range(self.elemCount("all")))
        elif 'elementGroups' in self.data.keys(): 
            return self.data['elementGroups'][name][:]
        else:
            return self.f['elementGroups'][name][:]        

#-------------------------------------------------------------------------------
#  elemCount
#-------------------------------------------------------------------------------
            
    def elemCount( self , elemGroup : str = "all" ) -> int:
  
        '''
      Returns the number of elements in this data set (in this cycle).
      
      Args:
        elemGroup:  specify an element group. In that case it returns the 
                    number of elements in that group.      
        '''
    
        if elemGroup == 'all':
            return len(self.elemNodes)
        else:
            return len(self.getElemGroup(elemGroup))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getElemIndex( self , elemID : int | list[int] ) -> list[int]:
   
        '''
          Return element Index
      
          Args:
            elemID:   integer or list of elementIDs
        '''
    
        if 'elements' in self.f.keys():    
            grp = self.f['elements']
        else:
            grp = self.data['elements']
            
        return list(grp['elementIDs']).index(elemID)      

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getElemIDs( self ) -> list[int]:
   
        '''
          Return element IDs in the set.
        '''
 
        if 'elements' in self.f.keys():        
            grp = self.f['elements']
        else:
            grp = self.data['elements']
                        
        return grp['elementIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def nodeCount( self , nodeGroup : str = 'all' ) -> int:
  
        '''
        Returns the number of nodes in the nodegroup.
    
        Args:
    
             nodeGroup:    nodeGroup name. By defaults the total number
                           of nodes is given.
        '''
    
        if nodeGroup == 'all':
            return self.coordinates.shape[0]
        else:
            if 'nodes' in self.f.keys():        
                grp = self.f['nodeGroups']
            else:
                grp = self.data['nodeGroups']
            
            return len(grp[nodeGroup])
            
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def unpackElements( self , a : list , offsets : list ) -> list:
             
        elemNodes = [None] * len(offsets)
        
        elemNodes[0] = a[0:offsets[0]]

        for i, (start, end) in enumerate(zip(offsets[:-1], offsets[1:]), 1):
            elemNodes[i] = a[start:end]
            
        return elemNodes
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getNodeIndex( self , nodeID : int) -> int:
   
        '''
        Returns the index of nodeID. This is the index of the node in the input file.
    
        Args:
    
            nodeID:   The nodeID. This can be an integer or a list.
        '''

        if 'nodes' in self.f.keys():                           
            grp = self.f['nodes']
        else:
            grp = self.data['nodes'] 
                                   
        return list(grp['nodeIDs']).index(nodeID)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getNodeIDs( self ) -> list[int]:
   
        '''
        Returns a list of all the nodeIDs.
        '''

        if 'nodes' in self.f.keys():                           
            grp = self.f['nodes']                
        else:
            grp = self.data['nodes']
            
        return grp['nodeIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getNodeGroupNames( self ) -> list[str]:
  
        '''
        Returns the names of all NodeGroups in this data set.
        '''
    
        if 'nodes' in self.f.keys():                           
            return list(self.f['nodeGroups'].keys())
        else:
            return list(self.data['nodeGroups'].keys())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getNodeGroup( self , name : str ) -> list[int]:
  
        '''
        Returns the nodeID that are present in a certain NodeGroup.
    
        Args:
      
            name:    Name of the nodegroup.       
        '''

        if 'nodes' in self.f.keys():                               
            grp = self.f['nodeGroups']        
        else:
            grp = self.data['nodeGroups']
            
        return grp[name][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
                         
    def rank( self ) -> int:
  
        '''
        Returns the rank of the problem. This is the number of spatial
            dimensions.
        '''
  
        return self.coordinates.shape[1]    

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def nodeDataSets( self ) -> list[str]:
  
        '''
        Returns all the labels of node datasets.
        '''
    
        grp = self.data['nodeData']
        return list(grp.keys())
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def elemDataSets( self ) -> list[str]:
  
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
      
    def getDisplacements( self , nodeID : int ):
  
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
    
    def getNodeData( self , label : str , nodeID : int =-1 ) ->list:
  
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

    def getElemData( self , label : str , elemID : int=-1 ) -> list:
  
        '''
        Returns the element data (label) of elements elemID 
                (can be a list or an integer).
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

    def getMacroFieldData( self , label :str ) -> list:
  
        '''
        Returns the element data (label) of elements 
            elemID (can be a list or an integer).
        '''
    
        grp  = self.data['macrofield']
        return grp[label][:]
         
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def particleCount( self , particleGroup : str  = 'all' ) -> int:
  
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
      
    def saveAsVTU( self , prefix : str = 'None' , cycles : int = -1 ):
  
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
  
            vtufile = prefix+'_t'+str(iCyc)+".vtu"
            
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
     
            for elemNodes in self.elemNodes:  
                insertElement( grid , elemNodes , self.rank() , 0 )
              
            # -- Write nodedata
  
            labels = self.nodeDataSets()
      
            for label in labels:            
                data = self.getNodeData( label )
                if data.ndim == 2:
                    if label == "displacements":
                        if data.shape[1] == 2:
                            newdata = np.zeros(shape=(data.shape[0],3))
                            newdata[:,:-1] = data
                            data = newdata
             												
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
            if iCyc == 11:
            
                # Step 2: Set up a renderer and render window
                stress_array = grid.GetPointData().GetArray("S22")
                grid.GetPointData().SetScalars(stress_array)

                # Create a mapper for the mesh
                mapper = vtk.vtkDataSetMapper()
                mapper.SetInputData(grid)
                mapper.SetScalarRange(stress_array.GetRange())  # Set the color range based on the stress data
                #mapper.SetScalarRange(0,100000000)
	
                # Apply a color map (for example, Jet)
                apply_rainbow_color_map(mapper)

                # Create an actor for the mesh
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)

                # Create a scalar bar actor for the legend
                scalar_bar = vtk.vtkScalarBarActor()
                scalar_bar.SetLookupTable(mapper.GetLookupTable())
                scalar_bar.SetTitle("Stress")
                scalar_bar.SetNumberOfLabels(4)  # Number of labels on the scalar bar

                grid.GetPointData().SetVectors(grid.GetPointData().GetArray("displacements"))
                
                warp_vector = vtk.vtkWarpVector()
                warp_vector.SetInputData(grid)
                warp_vector.SetScaleFactor(1.0)  # Apply magnification factor

                # Update the warp vector filter
                warp_vector.Update()

                mapper = vtk.vtkDataSetMapper()
                #mapper.SetInputData(grid)
                mapper.SetInputConnection(warp_vector.GetOutputPort())

                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                
                
                renderer = vtk.vtkRenderer()
                renderer.AddActor(actor)
                renderer.AddActor(scalar_bar)  # Add the scalar bar to the renderer                
                renderer.SetBackground(1, 1, 1)  # Set background to white

                render_window = vtk.vtkRenderWindow()
                render_window.AddRenderer(renderer)
                render_window.SetSize(800, 800)  # Set the size of the render window

                
                renderer.ResetCamera()
                
                camera = renderer.GetActiveCamera()
                camera.Azimuth(245)  # Rotate the camera around the vertical axis
                camera.Elevation(30)                
                camera.Zoom(0.5)
                
                
                                
                window_to_image_filter = vtk.vtkWindowToImageFilter()
                window_to_image_filter.SetInput(render_window)
                window_to_image_filter.Update()
                
                # Step 6: Save the image as a PNG file
                writer2 = vtk.vtkPNGWriter()
                writer2.SetFileName("output_image.png")
                writer2.SetInputData(window_to_image_filter.GetOutput())
                writer2.Write()
                '''
                                 
        writePVD( prefix , cycles )
        
        return cycles

#-------------------------------------------------------------------------------
#  saveAsDat
#-------------------------------------------------------------------------------
    
    def saveAsDat( self , fileName : str = 'None' , cycle : int = -1 , output : str = "dawn") -> None:
  
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
      
    def saveModes( self , prefix : str = 'None' ) -> None:
  
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

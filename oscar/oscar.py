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
    """Apply a rainbow color map to a VTK mapper.
    
    Creates and applies a rainbow-style lookup table that maps scalar values
    from blue to red (hue range 0.0 to 0.8).
    
    Args:
        mapper: VTK mapper object to which the color map will be applied.
        
    Returns:
        None
        
    Examples:
        >>> mapper = vtk.vtkDataSetMapper()
        >>> mapper.SetInputData(grid)
        >>> apply_rainbow_color_map(mapper)
    """
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
    """Check if the file is a valid oscar HDF5 file.
    
    Verifies that the file format is 'RNDF' and the version is at least 1.0.
    
    Args:
        f: HDF5 file object to check.
        
    Returns:
        None
        
    Raises:
        RuntimeError: If the file format is not 'RNDF' or version is less than 1.0.
        
    Examples:
        >>> import h5py
        >>> f = h5py.File('simulation.h5', 'r')
        >>> versionCheck(f)
    """
  
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
    """Convert an oscar HDF5 file to VTU format.
    
    Converts the specified oscar .h5 file to VTU (VTK Unstructured Grid) format.
    The VTU files are written to disk with the pattern fileName_t<cycle>.vtu and
    a PVD file is created for visualization.
    
    Args:
        fileName (str): The name of the h5 file (with or without .h5 extension).
        cycles (int or list, optional): List or integer of cycle numbers to convert.
            If -1 (default), all cycles are converted.
                     
    Returns:
        None
        
    Examples:
        >>> # Convert all cycles
        >>> convertToVTU('simulation.h5')
        
        >>> # Convert specific cycles
        >>> convertToVTU('simulation.h5', cycles=[1, 5, 10])
        
        >>> # Convert single cycle
        >>> convertToVTU('simulation', cycles=3)
    """
    
    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsVTU()
  
#-------------------------------------------------------------------------------
#  convertToDat
#-------------------------------------------------------------------------------

def convertToDat( fileName , cycles = -1 ):
    """Convert an oscar HDF5 file to Dawn DAT format.
    
    Converts the specified oscar .h5 file to Dawn DAT format.
    The DAT file is written to disk as fileName.dat.
    
    Args:
        fileName (str): The name of the h5 file (with or without .h5 extension).
        cycles (int or list, optional): List or integer of cycle numbers to convert.
            If -1 (default), all cycles are converted.
                     
    Returns:
        None
        
    Examples:
        >>> # Convert all cycles to DAT format
        >>> convertToDat('simulation.h5')
        
        >>> # Convert specific cycles
        >>> convertToDat('simulation', cycles=[1, 10, 20])
    """
    
    oscarFile = oscar( fileName , cycles )
  
    oscarFile.saveAsDat()
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
def writePVD( prefix , cycles , nProc = 1 ):
    """Write a ParaView Data (PVD) file for a collection of VTU files.
    
    Creates a PVD file that references multiple VTU files for time series
    visualization in ParaView. Supports both serial and parallel (multi-processor)
    data files.
    
    Args:
        prefix (str): Prefix for the PVD filename and referenced VTU files.
        cycles (list): List of cycle (timestep) numbers to include.
        nProc (int, optional): Number of processors. If > 1, includes files from
            multiple processors. Default is 1.
            
    Returns:
        None
        
    Raises:
        RuntimeError: If nProc < 1.
        
    Examples:
        >>> # Serial simulation
        >>> writePVD('results', [1, 2, 3, 4, 5])
        
        >>> # Parallel simulation with 4 processors
        >>> writePVD('parallel_results', [1, 2, 3], nProc=4)
    """
    
    if nProc < 1:
        print("writePVD: number of processros should by 1 or more.")
        raise RuntimeError    
    
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
    """Generate a VTU filename for use in PVD files.
    
    Creates a standardized filename for VTU files based on prefix, cycle number,
    and optionally processor number.
    
    Args:
        prefix (str): Filename prefix.
        iCyc (int): Cycle (timestep) number.
        iProc (int, optional): Processor number. If -1 (default), creates a
            serial filename. Otherwise creates a parallel filename.
            
    Returns:
        str: Formatted VTU filename.
        
    Examples:
        >>> # Serial filename
        >>> PVDfileName('simulation', 5)
        'simulation_t5.vtu'
        
        >>> # Parallel filename for processor 2
        >>> PVDfileName('simulation', 5, iProc=2)
        'simulation_p2_t5.vtu'
    """
    
    if iProc == -1:
        return prefix + "_t" + str(iCyc) + ".vtu"
    else:
        return prefix + "_p" + str(iProc) + "_t" + str(iCyc) + ".vtu"     
            
            
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


class oscar():

    '''
    A class to read and access data in the HDF file format. sdd
            
    Args: 
        
        fileName(str): filename May be with or without the extension .h5
        verbose(bool): verbose flag.            
    '''  

    def __init__( self , fileName : str , verbose : bool = False ):
        """Initialize an oscar object for reading HDF5 NDF files.
        
        Opens the specified HDF5 file and initializes data structures for
        accessing mesh, nodes, elements, and simulation results.
        
        Args:
            fileName (str): Path to the HDF5 file (with or without .h5 extension).
            verbose (bool, optional): If True, prints informational messages.
                Default is False.
                
        Raises:
            RuntimeError: If the file version is less than 1.0.
            
        Examples:
            >>> # Open a file silently
            >>> oscarFile = oscar('simulation.h5')
            
            >>> # Open with verbose output
            >>> oscarFile = oscar('simulation', verbose=True)
        """
            
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
        """Return a string representation of the oscar object.
        
        Provides a summary of the file contents, specifically the number
        of datasets (cycles) in the file.
        
        Returns:
            str: Description of the file contents.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> print(oscarFile)
            Single file with 10 datasets (cycles).
        """
   
        if 'cycleCount' in self.f.attrs.keys():  
            return "Single file with %d datasets (cycles)." %self.f.attrs['cycleCount']
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def setCycle( self , iCyc : int ):
        """Set the active cycle for data access.
        
        Changes the active cycle (timestep) from which data will be read.
        Updates internal data structures to point to the specified cycle.
        
        Args:
            iCyc (int): Cycle ID number to activate.
            
        Returns:
            None
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(5)
            >>> # Now all data access is from cycle 5
        """
    
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
        """Get the coordinates of specified node(s).
        
        Retrieves spatial coordinates for one or more nodes.
        
        Args:
            nodeID (int or list, optional): Node ID(s) according to internal numbering.
                If -1 (default), returns coordinates of all nodes.
                
        Returns:
            numpy.ndarray: Node coordinates. Shape is (n_nodes, n_dimensions) for
                multiple nodes or (n_dimensions,) for a single node.
                
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> # Get all coordinates
            >>> all_coords = oscarFile.getCoords()
            >>> # Get single node coordinates
            >>> node_coord = oscarFile.getCoords(5)
            >>> # Get multiple nodes
            >>> nodes_coord = oscarFile.getCoords([0, 5, 10])
        """
      
        if nodeID == -1:
            return self.coordinates[:]
        else:
            return self.coordinates[nodeID]

#-------------------------------------------------------------------------------
#  getCoords3
#-------------------------------------------------------------------------------
      
    def getCoords3( self ):
        """Get node coordinates in 3D format.
        
        Returns the coordinates always in 3D format, padding with zeros if
        the problem is 2D (rank equals 2).
        
        Args:
            None
            
        Returns:
            numpy.ndarray: Node coordinates with shape (n_nodes, 3).
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> # For 2D problems, z-coordinate will be padded with zeros
            >>> coords_3d = oscarFile.getCoords3()
            >>> print(coords_3d.shape)  # (n_nodes, 3)
        """
              
        if self.coordinates.shape[1] == 2:
            data = np.zeros(shape=(self.coordinates.shape[0],3))
            data[:,:2] = self.coordinates
            return data
        else:
            return self.coordinates[:]

#-------------------------------------------------------------------------------
#  getCycle
#-------------------------------------------------------------------------------
      
    def getCycle( self ) -> int:
        """Get the current active cycle number.
        
        Returns:
            int: The currently active cycle ID.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(3)
            >>> current = oscarFile.getCycle()
            >>> print(current)  # 3
        """
    
        return self.cycle  
      
#-------------------------------------------------------------------------------
#  getElemNodes
#-------------------------------------------------------------------------------
      
    def getElemNodes( self , elemID : int ) -> list:
        """Get the node connectivity for a specified element.
        
        Args:
            elemID (int): Element ID number.
            
        Returns:
            list: List of node IDs that form the element.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> nodes = oscarFile.getElemNodes(10)
            >>> print(nodes)  # [23, 45, 67, 89] for a quad element
        """
   
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
        """Get the number of nodes in a specified element.
        
        Args:
            elemID (int): Element ID number.
            
        Returns:
            int: Number of nodes in the element.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> n_nodes = oscarFile.getElemNodeCount(10)
            >>> print(n_nodes)  # 4 for quad, 8 for hex, etc.
        """
    
        return len(self.elemNodes[elemID])

#-------------------------------------------------------------------------------
#  getElemGroupNames
#-------------------------------------------------------------------------------
    
    def getElemGroupNames( self ) -> list:
        """Get the names of all element groups in the dataset.
        
        Returns:
            list: List of element group names as strings.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> groups = oscarFile.getElemGroupNames()
            >>> print(groups)  # ['domain', 'boundary', 'interface']
        """
        
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
        """Get element IDs belonging to a specified element group.
        
        Args:
            name (str, optional): Name of the element group. If "all" (default),
                returns all element IDs in the dataset.
                
        Returns:
            list[int]: List of element IDs in the group.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> # Get all elements
            >>> all_elems = oscarFile.getElemGroup()
            >>> # Get specific group
            >>> boundary_elems = oscarFile.getElemGroup('boundary')
        """
    
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
        """Get the number of elements in the dataset or a specific group.
        
        Args:
            elemGroup (str, optional): Name of an element group. If "all" (default),
                returns the total number of elements in the dataset.
                
        Returns:
            int: Number of elements in the dataset or group.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> total = oscarFile.elemCount()
            >>> print(f"Total elements: {total}")
            >>> boundary_count = oscarFile.elemCount('boundary')
        """
    
        if elemGroup == 'all':
            return len(self.elemNodes)
        else:
            return len(self.getElemGroup(elemGroup))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getElemIndex( self , elemID ) -> list[int]:
        """Get the internal index for a specified element ID.
        
        Converts an element ID to its internal storage index in the HDF5 file.
        
        Args:
            elemID (int or list): Element ID or list of element IDs.
            
        Returns:
            list[int]: Internal index or indices of the element(s).
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> index = oscarFile.getElemIndex(100)
            >>> print(index)  # Internal index for element 100
        """
    
        if 'elements' in self.f.keys():    
            grp = self.f['elements']
        else:
            grp = self.data['elements']
            
        return list(grp['elementIDs']).index(elemID)      

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getElemIDs( self ) -> list[int]:
        """Get all element IDs in the dataset.
        
        Returns:
            list[int]: Array of all element IDs.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> elem_ids = oscarFile.getElemIDs()
            >>> print(f"First 5 elements: {elem_ids[:5]}")
        """
 
        if 'elements' in self.f.keys():        
            grp = self.f['elements']
        else:
            grp = self.data['elements']
                        
        return grp['elementIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def nodeCount( self , nodeGroup : str = 'all' ) -> int:
        """Get the number of nodes in the dataset or a specific group.
        
        Args:
            nodeGroup (str, optional): Node group name. If 'all' (default),
                returns the total number of nodes in the dataset.
                
        Returns:
            int: Number of nodes in the dataset or group.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> total_nodes = oscarFile.nodeCount()
            >>> print(f"Total nodes: {total_nodes}")
            >>> fixed_nodes = oscarFile.nodeCount('fixed')
        """
    
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
        """Unpack element connectivity from flat array using offset indices.
        
        Converts a flat connectivity array and offset array into a list of
        element connectivity lists.
        
        Args:
            a (list): Flat array of node IDs.
            offsets (list): Array of offset indices marking element boundaries.
            
        Returns:
            list: List of element connectivity lists, where each element is
                a list of node IDs.
                
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> connectivity = [0, 1, 2, 3, 3, 4, 5, 6]
            >>> offsets = [4, 8]  # Two elements with 4 nodes each
            >>> elements = oscarFile.unpackElements(connectivity, offsets)
            >>> print(elements)  # [[0, 1, 2, 3], [3, 4, 5, 6]]
        """
        elemNodes = [None] * len(offsets)
        
        elemNodes[0] = a[0:offsets[0]]

        for i, (start, end) in enumerate(zip(offsets[:-1], offsets[1:]), 1):
            elemNodes[i] = a[start:end]
            
        return elemNodes
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getNodeIndex( self , nodeID : int) -> int:
        """Get the internal index for a specified node ID.
        
        Converts a node ID to its internal storage index in the HDF5 file.
        
        Args:
            nodeID (int or list): The node ID. Can be an integer or a list.
            
        Returns:
            int: Internal index of the node.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> index = oscarFile.getNodeIndex(250)
            >>> print(index)  # Internal index for node 250
        """

        if 'nodes' in self.f.keys():                           
            grp = self.f['nodes']
        else:
            grp = self.data['nodes'] 
                                   
        return list(grp['nodeIDs']).index(nodeID)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getNodeIDs( self ) -> list[int]:
        """Get all node IDs in the dataset.
        
        Returns:
            list[int]: Array of all node IDs.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> node_ids = oscarFile.getNodeIDs()
            >>> print(f"First 10 nodes: {node_ids[:10]}")
        """

        if 'nodes' in self.f.keys():                           
            grp = self.f['nodes']                
        else:
            grp = self.data['nodes']
            
        return grp['nodeIDs'][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getNodeGroupNames( self ) -> list[str]:
        """Get the names of all node groups in the dataset.
        
        Returns:
            list[str]: List of node group names.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> groups = oscarFile.getNodeGroupNames()
            >>> print(groups)  # ['fixed', 'loaded', 'free']
        """
    
        if 'nodes' in self.f.keys():                           
            return list(self.f['nodeGroups'].keys())
        else:
            return list(self.data['nodeGroups'].keys())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getNodeGroup( self , name : str ) -> list[int]:
        """Get node IDs belonging to a specified node group.
        
        Args:
            name (str): Name of the node group.
            
        Returns:
            list[int]: Array of node IDs in the group.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> fixed_nodes = oscarFile.getNodeGroup('fixed')
            >>> print(f"Fixed nodes: {fixed_nodes}")
        """

        if 'nodes' in self.f.keys():                               
            grp = self.f['nodeGroups']        
        else:
            grp = self.data['nodeGroups']
            
        return grp[name][:]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
                         
    def rank( self ) -> int:
        """Get the spatial dimension of the problem.
        
        Returns:
            int: Number of spatial dimensions (2 for 2D, 3 for 3D).
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> dim = oscarFile.rank()
            >>> print(f"Problem dimension: {dim}D")
        """
  
        return self.coordinates.shape[1]    

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def nodeDataSets( self ) -> list[str]:
        """Get the names of all available node data fields.
        
        Returns:
            list[str]: List of node data field names (e.g., 'displacements', 'S22').
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(5)
            >>> datasets = oscarFile.nodeDataSets()
            >>> print(datasets)  # ['displacements', 'S11', 'S22', 'S12']
        """
    
        grp = self.data['nodeData']
        return list(grp.keys())
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def elemDataSets( self ) -> list[str]:
        """Get the names of all available element data fields.
        
        Returns:
            list[str]: List of element data field names, or empty list if none exist.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(3)
            >>> datasets = oscarFile.elemDataSets()
            >>> print(datasets)  # ['stress', 'strain', 'damage']
        """
    
        if 'elementData' in self.data.keys():
            grp = self.data['elementData']
            return list(grp.keys())
        else:
            return []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def getDisplacements( self , nodeID : int ):
        """Get displacement values for specified node(s).
        
        Convenience method that calls getNodeData with label='displacements'.
        
        Args:
            nodeID (int or list): Node number(s). Can be an integer, list, or -1
                for all nodes.
                
        Returns:
            numpy.ndarray: Displacement values.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(10)
            >>> # Get all displacements
            >>> all_disp = oscarFile.getDisplacements(-1)
            >>> # Get single node displacement
            >>> node_disp = oscarFile.getDisplacements(5)
            >>> print(node_disp)  # [0.001, -0.002]
        """       
    
        return self.getNodeData( 'displacements' , nodeID )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getNodeData( self , label : str , nodeID : int =-1 ) ->list:
        """Get data values for specified node(s) from a named dataset.
        
        Args:
            label (str): Name of the node dataset (e.g., 'displacements', 'S22').
            nodeID (int or list, optional): Node number(s). Can be an integer,
                list, or -1 (default) for all nodes.
                
        Returns:
            numpy.ndarray: Data values. For scalar data, shape is (n_nodes,) or scalar.
                For vector data, shape is (n_nodes, n_components) or (n_components,).
                
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(5)
            >>> # Get stress S22 for all nodes
            >>> stress = oscarFile.getNodeData('S22')
            >>> # Get displacement for specific nodes
            >>> disp = oscarFile.getNodeData('displacements', [0, 10, 20])
        """
    
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
        """Get data values for specified element(s) from a named dataset.
        
        Args:
            label (str): Name of the element dataset.
            elemID (int or list, optional): Element number(s). Can be an integer,
                list, or -1 (default) for all elements.
                
        Returns:
            numpy.ndarray: Data values. For scalar data, shape is (n_elements,) or scalar.
                For vector data, shape is (n_elements, n_components) or (n_components,).
                
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(3)
            >>> # Get strain for all elements
            >>> strain = oscarFile.getElemData('strain')
            >>> # Get damage for specific element
            >>> damage = oscarFile.getElemData('damage', 25)
        """
    
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
        """Get macroscopic field data for a specified label.
        
        Args:
            label (str): Name of the macrofield dataset.
            
        Returns:
            numpy.ndarray: Macrofield data values.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(5)
            >>> macro_stress = oscarFile.getMacroFieldData('stress')
            >>> print(macro_stress)
        """
    
        grp  = self.data['macrofield']
        return grp[label][:]
         
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
      
    def particleCount( self , particleGroup : str  = 'all' ) -> int:
        """Get the number of particles in the dataset or a specific group.
        
        Args:
            particleGroup (str, optional): Particle group name. If 'all' (default),
                returns the total number of particles.
                
        Returns:
            int: Number of particles in the dataset or group.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> oscarFile.setCycle(10)
            >>> total_particles = oscarFile.particleCount()
            >>> active_particles = oscarFile.particleCount('active')
        """
    
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
        """Export data to VTK Unstructured Grid (VTU) format.
        
        Creates VTU files for visualization in ParaView or other VTK-compatible
        software. Also generates a PVD file for time series visualization.
        
        Args:
            prefix (str, optional): Prefix for output filenames. If 'None' (default),
                uses the prefix of the original HDF5 file.
            cycles (int or list, optional): Cycle number(s) to export. Can be an
                integer, list, or -1 (default) to export all cycles.
                
        Returns:
            list: List of cycle numbers that were exported.
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> # Export all cycles
            >>> oscarFile.saveAsVTU()
            >>> # Export specific cycles with custom prefix
            >>> oscarFile.saveAsVTU(prefix='results', cycles=[1, 5, 10])
            >>> # Export single cycle
            >>> oscarFile.saveAsVTU(cycles=3)
        """
    
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
            
                tt = 0
            
                if len(elemNodes) < 8:
                    tt = 2
                    
                insertElement( grid , elemNodes , self.rank() , tt )
              
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
            writer.SetDataModeToAscii()

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
        """Export data to Dawn DAT format.
        
        Creates a text-based DAT file containing mesh and simulation data.
        
        Args:
            fileName (str, optional): Output filename. If 'None' (default), uses
                the prefix of the HDF5 file. Extension .dat is added if not present.
            cycle (int, optional): Cycle number to export. If -1 (default),
                exports undeformed configuration from cycle 1.
            output (str, optional): Output format type. Either 'dawn' or 'pyfem'.
                Default is 'dawn'.
                
        Returns:
            None
            
        Examples:
            >>> oscarFile = oscar('simulation.h5')
            >>> # Export undeformed mesh
            >>> oscarFile.saveAsDat('mesh.dat')
            >>> # Export deformed configuration at cycle 10
            >>> oscarFile.saveAsDat('deformed.dat', cycle=10)
            >>> # Export in PyFEM format
            >>> oscarFile.saveAsDat('model.dat', cycle=5, output='pyfem')
        """
    
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

        elemGroups = ['all']
             
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
        """Export mode shapes to VTU format.
        
        Creates VTU files for visualization of mode shapes (e.g., from modal analysis).
        
        Args:
            prefix (str, optional): Prefix for output filenames. If 'None' (default),
                uses the prefix of the original HDF5 file.
                
        Returns:
            None
            
        Note:
            This method appears to be incomplete in the implementation.
            
        Examples:
            >>> oscarFile = oscar('modal_analysis.h5')
            >>> oscarFile.saveModes('modes')
            >>> # Creates modes-1.vtu, modes-2.vtu, etc.
        """
    
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

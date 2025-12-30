.. Oscar documentation master file, created by
   sphinx-quickstart on Mon Dec 30 12:05:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Oscar documentation
===================

Oscar is a Python package for reading and converting Finite Element HDF5 files in the NDF format. 
It provides tools for extracting simulation data and converting it to VTK formats for visualization.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   modules
   examples

Command-Line Interface
======================

h5tovtk - HDF5 to VTU Converter
--------------------------------

The ``h5tovtk`` command-line tool converts Oscar HDF5 files to VTK Unstructured Grid (VTU) format
for visualization in ParaView or other VTK-compatible software.

Installation
~~~~~~~~~~~~

After installing the oscar package, the ``h5tovtk`` command becomes available system-wide::

    pip install oscar

Or install from source in development mode::

    cd /path/to/oscar
    pip install -e .

Basic Usage
~~~~~~~~~~~

Convert an entire HDF5 file (all cycles)::

    h5tovtk simulation.h5

The command will create:

* Individual VTU files for each time step: ``simulation_t1.vtu``, ``simulation_t2.vtu``, etc.
* A PVD file for time series animation: ``simulation.pvd``

Command-Line Options
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    usage: h5tovtk [-h] [-c N] [-p PREFIX] [-v] filename

    positional arguments:
      filename              HDF5 file name (with or without .h5 extension)

    optional arguments:
      -h, --help            show this help message and exit
      -c N, --cycle N       Specific cycle number to convert (default: all cycles)
      -p PREFIX, --prefix PREFIX
                            Output prefix for VTU files (default: use input filename)
      -v, --verbose         Increase output verbosity

Examples
~~~~~~~~

**Convert all cycles:**

.. code-block:: bash

    h5tovtk simulation.h5

This converts all time steps in the file.

**Convert a specific cycle:**

.. code-block:: bash

    h5tovtk simulation.h5 -c 5

Converts only cycle number 5 from the simulation.

**Specify custom output prefix:**

.. code-block:: bash

    h5tovtk simulation.h5 -p results

Creates output files with the prefix ``results`` instead of ``simulation``.

**Verbose output:**

.. code-block:: bash

    h5tovtk simulation.h5 -v

Shows detailed progress information during conversion.

**Multi-processor files:**

For simulations run on multiple processors, the files are typically named with a ``_p#`` suffix
(e.g., ``simulation_p0.h5``, ``simulation_p1.h5``, etc.). Simply use the base name::

    h5tovtk simulation

The tool automatically detects all processor files and creates:

* Individual VTU files for each processor and cycle
* A combined PVD file for parallel visualization

File name without extension also works::

    h5tovtk simulation.h5

Output Files
~~~~~~~~~~~~

The converter creates the following output files:

**VTU Files**
  Individual unstructured grid files for each time step. Format: ``<prefix>_t<cycle>.vtu``
  
  Example: ``simulation_t1.vtu``, ``simulation_t2.vtu``, etc.

**PVD Files**
  ParaView Data files that reference all VTU files and provide time step information.
  Format: ``<prefix>.pvd``
  
  Example: ``simulation.pvd``

**Multi-processor outputs**
  For parallel simulations, files are organized by processor:
  
  * ``<prefix>_p<proc>_t<cycle>.vtu`` - Individual processor/cycle files
  * ``<prefix>.pvd`` - Combined time series for all processors

Opening in ParaView
~~~~~~~~~~~~~~~~~~~

After conversion, open the PVD file in ParaView:

1. Launch ParaView
2. File → Open → Select the ``.pvd`` file
3. Click "Apply" in the Properties panel
4. Use the time controls to animate through cycles

The PVD file automatically loads all time steps and provides proper time information
for animation and analysis.

Error Handling
~~~~~~~~~~~~~~

**File not found:**

.. code-block:: text

    Error: File 'simulation.h5' not found

Ensure the file path is correct and the file exists.

**Invalid HDF5 format:**

.. code-block:: text

    Error during conversion: ...

The file may not be a valid Oscar HDF5 file. Check that the file format is 'RNDF'
and version is at least 1.0.

**Cycle out of range:**

If you specify a cycle number that doesn't exist in the file, the conversion will fail.
Use ``-v`` flag for more detailed error information.

Troubleshooting
~~~~~~~~~~~~~~~

**Command not found after installation:**

If ``h5tovtk`` is not found after installation, ensure:

1. The package is properly installed: ``pip list | grep oscar``
2. Your Python scripts directory is in PATH
3. Try reinstalling: ``pip install --force-reinstall oscar``

**Permission errors:**

If you get permission errors when writing output files:

1. Check write permissions in the current directory
2. Try running from a different directory
3. Specify an output prefix in a writable location: ``-p /tmp/results``

**Memory issues with large files:**

For very large HDF5 files:

* Convert specific cycles instead of all at once
* Process multi-processor files separately
* Ensure adequate RAM is available

Python API
~~~~~~~~~~

The CLI internally uses the ``oscar`` class. For programmatic access, use the Python API directly:

.. code-block:: python

    from oscar import oscar
    
    # Open HDF5 file
    h5file = oscar('simulation.h5')
    
    # Convert all cycles
    h5file.saveAsVTU()
    
    # Convert specific cycles
    h5file.saveAsVTU(cycles=[1, 5, 10])
    
    # Custom output prefix
    h5file.saveAsVTU(prefix='results')

See the API documentation for more details on the ``oscar`` class and its methods.

Python API Reference
====================

The oscar package provides a Python API for programmatic access to HDF5 files in the NDF format.
This section documents the main ``oscar`` class and its methods.

The oscar Class
---------------

The ``oscar`` class provides methods to read and manipulate Finite Element HDF5 files.

Initialization
~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    # Open an HDF5 file
    oscarFile = oscar('simulation.h5')
    
    # Open with verbose output
    oscarFile = oscar('simulation.h5', verbose=True)

**Constructor:**

.. py:class:: oscar(fileName, verbose=False)

   Initialize an oscar object for reading HDF5 NDF files.
   
   :param fileName: Path to the HDF5 file (with or without .h5 extension)
   :type fileName: str
   :param verbose: If True, prints informational messages
   :type verbose: bool
   :raises RuntimeError: If the file version is less than 1.0

Cycle Management
~~~~~~~~~~~~~~~~

Methods for managing time steps (cycles) in the simulation data.

.. py:method:: setCycle(iCyc)

   Set the active cycle for data access.
   
   :param iCyc: Cycle ID number to activate
   :type iCyc: int
   
   Example::
   
       oscarFile.setCycle(5)  # Access data from cycle 5

.. py:method:: getCycle()

   Get the current active cycle number.
   
   :return: Currently active cycle ID
   :rtype: int
   
   Example::
   
       current_cycle = oscarFile.getCycle()
       print(f"Current cycle: {current_cycle}")

Node Operations
~~~~~~~~~~~~~~~

Methods for accessing node coordinates and information.

.. py:method:: getCoords(nodeID=-1)

   Get the coordinates of specified node(s).
   
   :param nodeID: Node ID(s). If -1 (default), returns all node coordinates.
   :type nodeID: int, list, or -1
   :return: Node coordinates array
   :rtype: numpy.ndarray
   
   Example::
   
       # Get all coordinates
       all_coords = oscarFile.getCoords()
       
       # Get single node
       node_coord = oscarFile.getCoords(5)
       
       # Get multiple nodes
       coords = oscarFile.getCoords([0, 5, 10])

.. py:method:: getCoords3()

   Get node coordinates in 3D format (padding with zeros for 2D problems).
   
   :return: Node coordinates with shape (n_nodes, 3)
   :rtype: numpy.ndarray
   
   Example::
   
       coords_3d = oscarFile.getCoords3()

.. py:method:: nodeCount(nodeGroup='all')

   Get the number of nodes in the dataset or a specific group.
   
   :param nodeGroup: Node group name or 'all' for total count
   :type nodeGroup: str
   :return: Number of nodes
   :rtype: int
   
   Example::
   
       total_nodes = oscarFile.nodeCount()
       fixed_nodes = oscarFile.nodeCount('fixed')

.. py:method:: getNodeIDs()

   Get all node IDs in the dataset.
   
   :return: Array of all node IDs
   :rtype: list[int]
   
   Example::
   
       node_ids = oscarFile.getNodeIDs()

.. py:method:: getNodeGroupNames()

   Get the names of all node groups.
   
   :return: List of node group names
   :rtype: list[str]
   
   Example::
   
       groups = oscarFile.getNodeGroupNames()
       print(f"Available node groups: {groups}")

.. py:method:: getNodeGroup(name)

   Get node IDs belonging to a specified node group.
   
   :param name: Name of the node group
   :type name: str
   :return: Array of node IDs in the group
   :rtype: list[int]
   
   Example::
   
       fixed_nodes = oscarFile.getNodeGroup('fixed')

Element Operations
~~~~~~~~~~~~~~~~~~

Methods for accessing element connectivity and information.

.. py:method:: getElemNodes(elemID)

   Get the node connectivity for a specified element.
   
   :param elemID: Element ID number
   :type elemID: int
   :return: List of node IDs forming the element
   :rtype: list
   
   Example::
   
       nodes = oscarFile.getElemNodes(10)
       print(f"Element 10 nodes: {nodes}")

.. py:method:: getElemNodeCount(elemID)

   Get the number of nodes in a specified element.
   
   :param elemID: Element ID number
   :type elemID: int
   :return: Number of nodes in the element
   :rtype: int
   
   Example::
   
       n_nodes = oscarFile.getElemNodeCount(10)

.. py:method:: elemCount(elemGroup='all')

   Get the number of elements in the dataset or a specific group.
   
   :param elemGroup: Element group name or 'all' for total count
   :type elemGroup: str
   :return: Number of elements
   :rtype: int
   
   Example::
   
       total_elements = oscarFile.elemCount()
       boundary_elements = oscarFile.elemCount('boundary')

.. py:method:: getElemIDs()

   Get all element IDs in the dataset.
   
   :return: Array of all element IDs
   :rtype: list[int]
   
   Example::
   
       elem_ids = oscarFile.getElemIDs()

.. py:method:: getElemGroupNames()

   Get the names of all element groups.
   
   :return: List of element group names
   :rtype: list[str]
   
   Example::
   
       groups = oscarFile.getElemGroupNames()

.. py:method:: getElemGroup(name='all')

   Get element IDs belonging to a specified element group.
   
   :param name: Element group name or 'all' for all elements
   :type name: str
   :return: List of element IDs in the group
   :rtype: list[int]
   
   Example::
   
       all_elems = oscarFile.getElemGroup()
       boundary_elems = oscarFile.getElemGroup('boundary')

Data Access
~~~~~~~~~~~

Methods for accessing node and element field data.

.. py:method:: nodeDataSets()

   Get the names of all available node data fields.
   
   :return: List of node data field names
   :rtype: list[str]
   
   Example::
   
       fields = oscarFile.nodeDataSets()
       print(f"Available node fields: {fields}")

.. py:method:: elemDataSets()

   Get the names of all available element data fields.
   
   :return: List of element data field names
   :rtype: list[str]
   
   Example::
   
       fields = oscarFile.elemDataSets()
       print(f"Available element fields: {fields}")

.. py:method:: getNodeData(label, nodeID=-1)

   Get node field data for specified node(s).
   
   :param label: Name of the data field
   :type label: str
   :param nodeID: Node ID or -1 for all nodes
   :type nodeID: int
   :return: Node field data
   :rtype: numpy.ndarray
   
   Example::
   
       # Get displacement for all nodes
       displacements = oscarFile.getNodeData('displacement')
       
       # Get data for specific node
       temp = oscarFile.getNodeData('temperature', nodeID=5)

.. py:method:: getElemData(label, elemID=-1)

   Get element field data for specified element(s).
   
   :param label: Name of the data field
   :type label: str
   :param elemID: Element ID or -1 for all elements
   :type elemID: int
   :return: Element field data
   :rtype: numpy.ndarray
   
   Example::
   
       # Get stress for all elements
       stress = oscarFile.getElemData('stress')
       
       # Get data for specific element
       strain = oscarFile.getElemData('strain', elemID=10)

.. py:method:: getDisplacements(nodeID)

   Get displacement data for specified node(s).
   
   :param nodeID: Node ID or list of node IDs
   :type nodeID: int or list
   :return: Displacement values
   :rtype: numpy.ndarray
   
   Example::
   
       disp = oscarFile.getDisplacements(5)

.. py:method:: getMacroFieldData(label)

   Get homogenized macro field data from RVE simulations.
   
   :param label: Name of the macro field
   :type label: str
   :return: Macro field data
   :rtype: numpy.ndarray
   
   Example::
   
       thermal_cond = oscarFile.getMacroFieldData('thermalCond')
       density = oscarFile.getMacroFieldData('density')

Export Methods
~~~~~~~~~~~~~~

Methods for exporting data to various file formats.

.. py:method:: saveAsVTU(prefix='None', cycles=-1)

   Export data to VTK Unstructured Grid (VTU) format.
   
   :param prefix: Output filename prefix (default: use input filename)
   :type prefix: str
   :param cycles: Cycle number(s) to export. Can be int, list, or -1 for all
   :type cycles: int or list
   :return: List of exported cycle numbers
   :rtype: list
   
   Example::
   
       # Export all cycles
       oscarFile.saveAsVTU()
       
       # Export specific cycles
       oscarFile.saveAsVTU(cycles=[1, 5, 10])
       
       # Custom prefix
       oscarFile.saveAsVTU(prefix='results')

.. py:method:: saveAsDat(fileName='None', cycle=-1, output='dawn')

   Export data to DAT format.
   
   :param fileName: Output filename (default: use input filename)
   :type fileName: str
   :param cycle: Cycle number to export or -1 for current cycle
   :type cycle: int
   :param output: Output format type
   :type output: str
   
   Example::
   
       oscarFile.saveAsDat(fileName='output.dat', cycle=5)

.. py:method:: saveModes(prefix='None')

   Export modal analysis data.
   
   :param prefix: Output filename prefix
   :type prefix: str
   
   Example::
   
       oscarFile.saveModes(prefix='modes')

Utility Methods
~~~~~~~~~~~~~~~

.. py:method:: rank()

   Get the spatial dimension of the problem.
   
   :return: Number of spatial dimensions (2 or 3)
   :rtype: int
   
   Example::
   
       dim = oscarFile.rank()
       print(f"Problem is {dim}D")

.. py:method:: particleCount(particleGroup='all')

   Get the number of particles in the dataset or a specific group.
   
   :param particleGroup: Particle group name or 'all'
   :type particleGroup: str
   :return: Number of particles
   :rtype: int
   
   Example::
   
       total_particles = oscarFile.particleCount()

Complete API Example
~~~~~~~~~~~~~~~~~~~~

Here's a complete example demonstrating common API usage patterns:

.. code-block:: python

    from oscar import oscar
    import numpy as np
    
    # Open the HDF5 file
    oscarFile = oscar('simulation.h5', verbose=True)
    
    # Get basic information
    print(f"Total cycles: {oscarFile.cycleCount}")
    print(f"Total nodes: {oscarFile.nodeCount()}")
    print(f"Total elements: {oscarFile.elemCount()}")
    print(f"Problem dimension: {oscarFile.rank()}D")
    
    # Iterate through cycles
    for cycle in range(1, oscarFile.cycleCount + 1):
        oscarFile.setCycle(cycle)
        
        # Get node data
        coords = oscarFile.getCoords()
        displacements = oscarFile.getNodeData('displacement')
        
        # Get element data
        stress = oscarFile.getElemData('stress')
        
        # Process data
        max_disp = np.max(np.linalg.norm(displacements, axis=1))
        print(f"Cycle {cycle}: Max displacement = {max_disp:.6f}")
    
    # Access specific node groups
    fixed_nodes = oscarFile.getNodeGroup('fixed')
    print(f"Number of fixed nodes: {len(fixed_nodes)}")
    
    # Access element groups
    boundary_elems = oscarFile.getElemGroup('boundary')
    print(f"Number of boundary elements: {len(boundary_elems)}")
    
    # List available data fields
    print("Available node fields:", oscarFile.nodeDataSets())
    print("Available element fields:", oscarFile.elemDataSets())
    
    # Export to VTU for visualization
    oscarFile.saveAsVTU(prefix='results')
    print("Exported to VTU format")

RVE Utilities API
-----------------

The ``RVEutils`` module provides utilities for Representative Volume Element (RVE) 
homogenization analysis.

.. code-block:: python

    from oscar.RVEutils import getConductivity, getMechanical, getCTE

Thermal Properties
~~~~~~~~~~~~~~~~~~

.. py:function:: getConductivity(baseName, props=None)

   Get thermal conductivity tensor from RVE simulations.
   
   :param baseName: Base name of HDF5 files
   :type baseName: str
   :param props: Additional properties (optional)
   :type props: str
   :return: Thermal conductivity components [k_x, k_y, k_z]
   :rtype: list
   
   Example::
   
       cond = getConductivity('rve_thermal')
       print(f"Conductivity: {cond}")

.. py:function:: getCTE(baseName)

   Get coefficients of thermal expansion from RVE analysis.
   
   :param baseName: Base name of HDF5 file
   :type baseName: str
   :return: Thermal expansion coefficients
   :rtype: list
   
   Example::
   
       cte = getCTE('rve_thermal')

.. py:function:: getCapacity(baseName)

   Get specific heat capacity from RVE homogenization.
   
   :param baseName: Base name of HDF5 file
   :type baseName: str
   :return: Specific heat capacity
   :rtype: float
   
   Example::
   
       capacity = getCapacity('rve_thermal')

.. py:function:: getDensity(baseName)

   Get effective density from RVE homogenization.
   
   :param baseName: Base name of HDF5 file
   :type baseName: str
   :return: Density value
   :rtype: float
   
   Example::
   
       rho = getDensity('rve_thermal')

Mechanical Properties
~~~~~~~~~~~~~~~~~~~~~

.. py:function:: getMechanical(baseName, eps=0.01)

   Get effective mechanical properties from RVE homogenization.
   
   :param baseName: Base name of HDF5 files
   :type baseName: str
   :param eps: Strain magnitude used in simulations
   :type eps: float
   :return: Tuple of (E, G, nu, cmat) where:
   
       - E: Young's moduli [E_x, E_y, E_z]
       - G: Shear moduli [G_xy, G_xz, G_yz]
       - nu: Poisson's ratios [nu_xy, nu_yz, nu_xz]
       - cmat: 6x6 stiffness matrix
   
   :rtype: tuple
   
   Example::
   
       E, G, nu, C = getMechanical('rve_mechanical')
       print(f"Young's moduli: {E}")
       print(f"Shear moduli: {G}")
       print(f"Poisson's ratios: {nu}")

Complete RVE Example
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar.RVEutils import *
    
    # Base name for RVE simulations
    base = 'composite_rve'
    
    # Get thermal properties
    conductivity = getConductivity(base)
    cte = getCTE(base)
    capacity = getCapacity(base)
    density = getDensity(base)
    
    print("Thermal Properties:")
    print(f"  Conductivity: {conductivity}")
    print(f"  CTE: {cte}")
    print(f"  Capacity: {capacity}")
    print(f"  Density: {density}")
    
    # Get mechanical properties
    E, G, nu, C = getMechanical(base, eps=0.01)
    
    print("\nMechanical Properties:")
    print(f"  Young's moduli: {E}")
    print(f"  Shear moduli: {G}")
    print(f"  Poisson's ratios: {nu}")
    print(f"  Stiffness matrix shape: {C.shape}")

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

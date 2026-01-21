Python API Reference
====================

The oscar package provides a Python API for programmatic access to HDF5 files in the NDF format.
This documentation is automatically generated from the source code docstrings.

Quick Start
-----------

.. code-block:: python

    from oscar import oscar
    
    # Open an HDF5 file
    oscarFile = oscar('simulation.h5')
    
    # Get basic information
    print(f"Total cycles: {oscarFile.cycleCount}")
    print(f"Total nodes: {oscarFile.nodeCount()}")
    print(f"Total elements: {oscarFile.elemCount()}")
    
    # Set active cycle and get data
    oscarFile.setCycle(5)
    coords = oscarFile.getCoords()
    displacements = oscarFile.getNodeData('displacement')
    
    # Export to VTU format
    oscarFile.saveAsVTU()

Core Module
-----------

Main Functions
~~~~~~~~~~~~~~

.. autofunction:: oscar.oscar.convertToVTU

.. autofunction:: oscar.oscar.writePVD

.. autofunction:: oscar.oscar.PVDfileName

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: oscar.oscar.apply_rainbow_color_map

.. autofunction:: oscar.oscar.versionCheck

The oscar Class
---------------

.. autoclass:: oscar.oscar.oscar
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource

Complete API Example
--------------------

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

RVE Utilities Module
--------------------

The ``RVEutils`` module provides utilities for Representative Volume Element (RVE) 
homogenization analysis to extract effective material properties.

.. code-block:: python

    from oscar.RVEutils import getConductivity, getMechanical, getCTE

Functions
~~~~~~~~~

.. autofunction:: oscar.RVEutils.getFileName

.. autofunction:: oscar.RVEutils.getConductivity

.. autofunction:: oscar.RVEutils.getCTE

.. autofunction:: oscar.RVEutils.getCapacity

.. autofunction:: oscar.RVEutils.getDensity

.. autofunction:: oscar.RVEutils.getMechanical

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

VTK Utilities Module
--------------------

The ``VTKutils`` module provides utilities for working with VTK data structures.

Functions
~~~~~~~~~

.. autofunction:: oscar.VTKutils.insertElement

Examples
========

This section provides practical examples of using the Oscar package for various tasks.

Command-Line Examples
---------------------

Using h5tovtk
~~~~~~~~~~~~~

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

Python API Examples
-------------------

Basic File Access
~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    # Open HDF5 file
    h5file = oscar('simulation.h5')
    
    # Get basic information
    print(f"Total cycles: {h5file.cycleCount}")
    print(f"Total nodes: {h5file.nodeCount()}")
    print(f"Total elements: {h5file.elemCount()}")
    print(f"Problem dimension: {h5file.rank()}D")

Converting to VTU Format
~~~~~~~~~~~~~~~~~~~~~~~~

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

Accessing Node Data
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    # Open file and set cycle
    h5file = oscar('simulation.h5')
    h5file.setCycle(5)
    
    # Get all node coordinates
    coords = h5file.getCoords()
    
    # Get coordinates for specific node
    node_coord = h5file.getCoords(nodeID=10)
    
    # Get displacement data for all nodes
    displacements = h5file.getNodeData('displacement')
    
    # Get node groups
    groups = h5file.getNodeGroupNames()
    fixed_nodes = h5file.getNodeGroup('fixed')

Accessing Element Data
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    # Open file and set cycle
    h5file = oscar('simulation.h5')
    h5file.setCycle(5)
    
    # Get element connectivity
    elem_nodes = h5file.getElemNodes(elemID=0)
    
    # Get stress data for all elements
    stress = h5file.getElemData('stress')
    
    # Get element groups
    groups = h5file.getElemGroupNames()
    boundary_elems = h5file.getElemGroup('boundary')

Processing Multiple Cycles
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    import numpy as np
    
    # Open the HDF5 file
    h5file = oscar('simulation.h5', verbose=True)
    
    # Initialize storage for results
    max_displacements = []
    
    # Iterate through cycles
    for cycle in range(1, h5file.cycleCount + 1):
        h5file.setCycle(cycle)
        
        # Get node data
        displacements = h5file.getNodeData('displacement')
        
        # Process data
        max_disp = np.max(np.linalg.norm(displacements, axis=1))
        max_displacements.append(max_disp)
        
        print(f"Cycle {cycle}: Max displacement = {max_disp:.6f}")
    
    # Plot or save results
    print(f"Maximum displacement overall: {max(max_displacements):.6f}")

Working with Node and Element Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    h5file = oscar('simulation.h5')
    h5file.setCycle(1)
    
    # List available groups
    node_groups = h5file.getNodeGroupNames()
    elem_groups = h5file.getElemGroupNames()
    
    print(f"Node groups: {node_groups}")
    print(f"Element groups: {elem_groups}")
    
    # Access specific node groups
    fixed_nodes = h5file.getNodeGroup('fixed')
    loaded_nodes = h5file.getNodeGroup('loaded')
    
    print(f"Number of fixed nodes: {len(fixed_nodes)}")
    print(f"Number of loaded nodes: {len(loaded_nodes)}")
    
    # Access specific element groups
    boundary_elems = h5file.getElemGroup('boundary')
    domain_elems = h5file.getElemGroup('domain')
    
    print(f"Number of boundary elements: {len(boundary_elems)}")
    print(f"Number of domain elements: {len(domain_elems)}")

Extracting Field Data
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    
    h5file = oscar('simulation.h5')
    h5file.setCycle(10)
    
    # List available data fields
    node_fields = h5file.nodeDataSets()
    elem_fields = h5file.elemDataSets()
    
    print(f"Available node fields: {node_fields}")
    print(f"Available element fields: {elem_fields}")
    
    # Extract specific field data
    temperature = h5file.getNodeData('temperature')
    stress = h5file.getElemData('stress')
    strain = h5file.getElemData('strain')

RVE Homogenization Examples
----------------------------

Extracting Thermal Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar.RVEutils import getConductivity, getCTE, getCapacity, getDensity
    
    # Base name for RVE simulation files
    base = 'composite_rve'
    
    # Get thermal properties
    conductivity = getConductivity(base)
    cte = getCTE(base)
    capacity = getCapacity(base)
    density = getDensity(base)
    
    print("Thermal Properties:")
    print(f"  Conductivity [W/mK]: {conductivity}")
    print(f"  CTE [1/K]: {cte}")
    print(f"  Specific heat capacity [J/kgK]: {capacity}")
    print(f"  Density [kg/m³]: {density}")

Extracting Mechanical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar.RVEutils import getMechanical
    
    # Base name for RVE simulation files
    base = 'composite_rve'
    
    # Get mechanical properties (assuming 1% strain was applied)
    E, G, nu, C = getMechanical(base, eps=0.01)
    
    print("Mechanical Properties:")
    print(f"  Young's moduli [GPa]:")
    print(f"    E_x = {E[0]/1e9:.2f}")
    print(f"    E_y = {E[1]/1e9:.2f}")
    print(f"    E_z = {E[2]/1e9:.2f}")
    print(f"  Shear moduli [GPa]:")
    print(f"    G_xy = {G[0]/1e9:.2f}")
    print(f"    G_xz = {G[1]/1e9:.2f}")
    print(f"    G_yz = {G[2]/1e9:.2f}")
    print(f"  Poisson's ratios:")
    print(f"    nu_xy = {nu[0]:.3f}")
    print(f"    nu_yz = {nu[1]:.3f}")
    print(f"    nu_xz = {nu[2]:.3f}")
    
    # Stiffness matrix
    print(f"\nStiffness matrix shape: {C.shape}")
    print("Stiffness matrix [GPa]:")
    print(C / 1e9)

Complete Homogenization Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar.RVEutils import *
    import json
    
    # Base name for RVE simulations
    base = 'composite_rve'
    
    # Extract all properties
    properties = {
        'thermal': {
            'conductivity': getConductivity(base),
            'expansion_coefficient': getCTE(base),
            'specific_heat': getCapacity(base),
            'density': getDensity(base)
        }
    }
    
    # Get mechanical properties
    E, G, nu, C = getMechanical(base, eps=0.01)
    properties['mechanical'] = {
        'youngs_moduli': E,
        'shear_moduli': G,
        'poissons_ratios': nu,
        'stiffness_matrix': C.tolist()
    }
    
    # Save to file
    with open('material_properties.json', 'w') as f:
        json.dump(properties, f, indent=2)
    
    print("Material properties saved to material_properties.json")

Advanced Examples
-----------------

Custom Data Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    import numpy as np
    import matplotlib.pyplot as plt
    
    h5file = oscar('simulation.h5')
    
    # Track maximum stress over time
    times = []
    max_stresses = []
    
    for cycle in range(1, h5file.cycleCount + 1):
        h5file.setCycle(cycle)
        
        # Get stress data
        stress = h5file.getElemData('stress')
        
        # Calculate von Mises stress
        von_mises = np.sqrt(0.5 * (
            (stress[:, 0] - stress[:, 1])**2 +
            (stress[:, 1] - stress[:, 2])**2 +
            (stress[:, 2] - stress[:, 0])**2 +
            6 * (stress[:, 3]**2 + stress[:, 4]**2 + stress[:, 5]**2)
        ))
        
        max_stress = np.max(von_mises)
        max_stresses.append(max_stress)
        times.append(cycle)
    
    # Plot results
    plt.figure()
    plt.plot(times, max_stresses)
    plt.xlabel('Cycle')
    plt.ylabel('Maximum von Mises Stress [Pa]')
    plt.title('Stress Evolution')
    plt.grid(True)
    plt.savefig('stress_evolution.png')
    print("Plot saved to stress_evolution.png")

Batch Processing
~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar
    import os
    import glob
    
    # Find all HDF5 files in directory
    h5_files = glob.glob('*.h5')
    
    for filename in h5_files:
        print(f"Processing {filename}...")
        
        try:
            h5file = oscar(filename)
            
            # Convert to VTU
            h5file.saveAsVTU()
            
            # Extract summary information
            summary = {
                'file': filename,
                'cycles': h5file.cycleCount,
                'nodes': h5file.nodeCount(),
                'elements': h5file.elemCount(),
                'dimension': h5file.rank()
            }
            
            print(f"  Converted {summary['cycles']} cycles")
            print(f"  Mesh: {summary['nodes']} nodes, {summary['elements']} elements")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("Batch processing complete!")

Parallel File Processing
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from oscar import oscar, writePVD
    import os
    
    # Base name for parallel simulation
    prefix = 'parallel_sim'
    
    # Detect number of processors
    nProc = 0
    while os.path.exists(f'{prefix}_p{nProc}.h5'):
        nProc += 1
    
    print(f"Found {nProc} processor files")
    
    # Process each processor file
    all_cycles = None
    for iProc in range(nProc):
        proc_file = f'{prefix}_p{iProc}.h5'
        print(f"Processing processor {iProc}...")
        
        h5file = oscar(proc_file)
        cycles = h5file.saveAsVTU()
        
        if all_cycles is None:
            all_cycles = cycles
    
    # Create combined PVD file
    writePVD(prefix, all_cycles, nProc)
    print(f"Created combined PVD file: {prefix}.pvd")

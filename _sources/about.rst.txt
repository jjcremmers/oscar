Data Structure
==============

Oscar HDF5 File Format (NDF/RNDF)
----------------------------------

The oscar package reads HDF5 files in the NDF (Numeric Data Format) or RNDF format.
This section describes the hierarchical data structure within these files.

File Structure Overview
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    simulation.h5 (HDF5 File)
    │
    ├── Attributes (File-level metadata)
    │   ├── fileFormat: "RNDF"
    │   ├── version: 1.0+
    │   └── cycleCount: N
    │
    ├── nodes (Optional, static mesh)
    │   ├── coordinates: [n_nodes × n_dim] array
    │   └── nodeIDs: [n_nodes] array
    │
    ├── elements (Optional, static mesh)
    │   ├── connectivity: [n_connections] flat array
    │   ├── offsets: [n_elements] offset indices
    │   └── elementIDs: [n_elements] array
    │
    ├── nodeGroups (Optional)
    │   ├── <group_name_1>: [node_ids] array
    │   ├── <group_name_2>: [node_ids] array
    │   └── ...
    │
    ├── elementGroups (Optional)
    │   ├── <group_name_1>: [element_ids] array
    │   ├── <group_name_2>: [element_ids] array
    │   └── ...
    │
    ├── cycle1 (First time step)
    │   ├── nodes (Optional, dynamic mesh)
    │   │   ├── coordinates: [n_nodes × n_dim] array
    │   │   └── nodeIDs: [n_nodes] array
    │   │
    │   ├── elements (Optional, dynamic mesh)
    │   │   ├── connectivity: [n_connections] flat array
    │   │   ├── offsets: [n_elements] offset indices
    │   │   └── elementIDs: [n_elements] array
    │   │
    │   ├── nodeGroups (Optional)
    │   │   └── <group_name>: [node_ids] array
    │   │
    │   ├── elementGroups (Optional)
    │   │   └── <group_name>: [element_ids] array
    │   │
    │   ├── NodeData
    │   │   ├── <field_name_1>: [n_nodes × n_components] array
    │   │   ├── <field_name_2>: [n_nodes × n_components] array
    │   │   └── ... (e.g., displacement, velocity, temperature)
    │   │
    │   ├── ElemData
    │   │   ├── <field_name_1>: [n_elements × n_components] array
    │   │   ├── <field_name_2>: [n_elements × n_components] array
    │   │   └── ... (e.g., stress, strain, damage)
    │   │
    │   └── MacroFields (Optional, for RVE homogenization)
    │       ├── <field_name_1>: scalar or array
    │       └── ... (e.g., thermalCond, stresses, density)
    │
    ├── cycle2 (Second time step)
    │   └── ... (same structure as cycle1)
    │
    ├── cycle3
    │   └── ...
    │
    └── cycleN
        └── ...

Detailed Data Structure
~~~~~~~~~~~~~~~~~~~~~~~~

Root Level Attributes
^^^^^^^^^^^^^^^^^^^^^

The HDF5 file contains global attributes that describe the file format and contents:

* **fileFormat** (string): Must be "RNDF" for oscar compatibility
* **version** (float): File format version, must be ≥ 1.0
* **cycleCount** (int): Total number of time steps (cycles) in the file

Mesh Data Structure
^^^^^^^^^^^^^^^^^^^

Mesh data can be stored at two levels:

1. **Static Mesh** (root level): Used when the mesh doesn't change between cycles
2. **Dynamic Mesh** (cycle level): Used when the mesh changes with each time step

**Nodes Group:**

* **coordinates**: NumPy array of shape ``[n_nodes, n_dim]`` where ``n_dim`` is 2 or 3
* **nodeIDs**: Array of unique node identifiers

**Elements Group:**

* **connectivity**: Flat array containing all node IDs for all elements
* **offsets**: Array marking where each element's connectivity ends in the flat array
* **elementIDs**: Array of unique element identifiers

**Example:**

.. code-block:: python

    # Element connectivity stored as:
    connectivity = [0, 1, 2, 3,  # Element 0 (quad)
                    3, 4, 5, 6]  # Element 1 (quad)
    offsets = [4, 8]             # Element 0 ends at index 4, Element 1 at 8

Node and Element Groups
^^^^^^^^^^^^^^^^^^^^^^^^

Groups provide a way to organize and select subsets of nodes or elements:

* **nodeGroups**: Named collections of node IDs (e.g., "fixed", "loaded", "boundary")
* **elementGroups**: Named collections of element IDs (e.g., "domain", "interface")

Field Data Structure
^^^^^^^^^^^^^^^^^^^^

**NodeData Group:**

Contains field data defined at nodes (e.g., displacements, temperatures):

* Each field is stored as an array of shape ``[n_nodes, n_components]``
* Common fields: "displacement", "velocity", "temperature", "pressure"

**ElemData Group:**

Contains field data defined at element integration points:

* Each field is stored as an array of shape ``[n_elements, n_components]``
* Common fields: "stress", "strain", "damage", "plastic_strain"

**MacroFields Group (RVE Analysis):**

For Representative Volume Element homogenization:

* Stores effective/homogenized properties
* Examples: "thermalCond", "stresses", "density", "specificCapac"
* Can be scalar values or arrays

Multi-Processor Files
~~~~~~~~~~~~~~~~~~~~~~

For parallel simulations, data is distributed across multiple files:

.. code-block:: text

    simulation_p0.h5  ← Processor 0 partition
    simulation_p1.h5  ← Processor 1 partition
    simulation_p2.h5  ← Processor 2 partition
    ...
    simulation_pN.h5  ← Processor N partition

Each file contains the same structure but with a subset of the mesh and data.

Data Access Pattern
~~~~~~~~~~~~~~~~~~~

The oscar class provides a unified interface to access this hierarchical structure:

.. code-block:: text

    oscar Object
    │
    ├── File Management
    │   ├── prefix: Base filename
    │   ├── cycleCount: Total cycles
    │   └── cycle: Current active cycle
    │
    ├── Mesh Access
    │   ├── coordinates: Node positions
    │   ├── elemNodes: Element connectivity
    │   ├── getCoords(nodeID)
    │   └── getElemNodes(elemID)
    │
    ├── Group Access
    │   ├── getNodeGroupNames()
    │   ├── getNodeGroup(name)
    │   ├── getElemGroupNames()
    │   └── getElemGroup(name)
    │
    ├── Data Access
    │   ├── nodeDataSets()
    │   ├── elemDataSets()
    │   ├── getNodeData(label, nodeID)
    │   ├── getElemData(label, elemID)
    │   └── getMacroFieldData(label)
    │
    └── Export
        ├── saveAsVTU()
        ├── saveAsDat()
        └── saveModes()

Workflow Diagram
~~~~~~~~~~~~~~~~

.. code-block:: text

    ┌─────────────────┐
    │   HDF5 File     │
    │  (simulation.h5)│
    └────────┬────────┘
             │
             │ open
             ▼
    ┌─────────────────┐
    │  oscar Object   │
    │                 │
    │  setCycle(n)    │◄─── Select time step
    └────────┬────────┘
             │
             ├─────► getCoords()          ─► Node coordinates
             │
             ├─────► getNodeData(label)   ─► Field data at nodes
             │
             ├─────► getElemData(label)   ─► Field data at elements
             │
             ├─────► getNodeGroup(name)   ─► Node groups
             │
             ├─────► getElemGroup(name)   ─► Element groups
             │
             └─────► saveAsVTU()          ─► Export to VTK format
                     └─► .vtu files
                     └─► .pvd file

Memory Layout
~~~~~~~~~~~~~

**Efficient Data Storage:**

.. code-block:: text

    Node Coordinates (2D Example):
    ┌──────────┬──────────┐
    │    X     │    Y     │
    ├──────────┼──────────┤
    │  0.0     │   0.0    │  ← Node 0
    │  1.0     │   0.0    │  ← Node 1
    │  1.0     │   1.0    │  ← Node 2
    │  0.0     │   1.0    │  ← Node 3
    └──────────┴──────────┘
    
    Element Connectivity (Flat Array):
    ┌───┬───┬───┬───┬───┬───┬───┬───┐
    │ 0 │ 1 │ 2 │ 3 │ 1 │ 4 │ 5 │ 2 │
    └───┴───┴───┴───┴───┴───┴───┴───┘
      └─────────┘   └─────────────┘
        Element 0      Element 1
    
    Offsets Array:
    ┌───┬───┐
    │ 4 │ 8 │  ← Element boundaries
    └───┴───┘
    
    Field Data (Displacement, 2D):
    ┌──────────┬──────────┐
    │   U_X    │   U_Y    │
    ├──────────┼──────────┤
    │  0.001   │  0.002   │  ← Node 0
    │  0.003   │  0.001   │  ← Node 1
    │  0.004   │  0.005   │  ← Node 2
    └──────────┴──────────┘

Type Definitions
~~~~~~~~~~~~~~~~

The oscar package uses the following data types:

* **Coordinates**: ``numpy.ndarray`` of shape ``(n, 2)`` or ``(n, 3)``
* **Connectivity**: ``list[list[int]]`` or flat array with offsets
* **Field Data**: ``numpy.ndarray`` of shape ``(n_entities, n_components)``
* **Node/Element IDs**: ``numpy.ndarray`` or ``list[int]``
* **Groups**: ``list[int]`` containing entity IDs

API Reference
=============

.. toctree::
   :maxdepth: 4

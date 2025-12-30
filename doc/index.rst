.. Oscar documentation master file, created by
   sphinx-quickstart on Mon Dec 30 12:05:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Oscar documentation
===================

Oscar is a Python package for reading and converting Finite Element HDF5 files in the NDF format. 
It provides tools for extracting simulation data and converting it to VTK formats for visualization.

For installation instructions, see :doc:`install`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   api
   examples

Command-Line Interface
======================

h5tovtk - HDF5 to VTU Converter
--------------------------------

The ``h5tovtk`` command-line tool converts Oscar HDF5 files to VTK Unstructured Grid (VTU) format
for visualization in ParaView or other VTK-compatible software.

For detailed usage examples, see :doc:`examples`.

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

For complete API documentation, see :doc:`api`.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

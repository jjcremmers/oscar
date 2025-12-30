Installation
============

Oscar is a Python package for reading and converting Finite Element HDF5 files in the NDF format.

Requirements
------------

* Python >= 3.10
* h5py >= 3.6
* vtk >= 9.2
* numpy

Installation from PyPI
----------------------

The easiest way to install oscar is using pip::

    pip install oscar

This will install the package and all required dependencies, making the ``h5tovtk`` 
command-line tool available system-wide.

Installation from Source
------------------------

To install from source or for development purposes:

1. Clone the repository::

    git clone https://github.com/jjcremmers/oscar.git
    cd oscar

2. Install in development mode::

    pip install -e .

This installs the package in editable mode, allowing you to make changes to the 
source code without reinstalling.

Verifying Installation
----------------------

After installation, verify that the package is correctly installed:

**Check the package version:**

.. code-block:: bash

    pip list | grep oscar

**Test the CLI tool:**

.. code-block:: bash

    h5tovtk --help

You should see the help message for the h5tovtk command.

**Test the Python API:**

.. code-block:: python

    from oscar import oscar
    print("Oscar package imported successfully!")

Troubleshooting
---------------

**Command not found:**

If ``h5tovtk`` is not found after installation:

1. Ensure your Python scripts directory is in your PATH
2. Try reinstalling: ``pip install --force-reinstall oscar``
3. Check if the package is installed: ``pip show oscar``

**ImportError:**

If you get import errors:

1. Verify all dependencies are installed: ``pip install h5py vtk numpy``
2. Check your Python version: ``python --version`` (must be >= 3.10)
3. Try reinstalling with ``pip install --upgrade oscar``

**Permission errors:**

If you get permission errors during installation:

* Use a virtual environment (recommended)
* Install for your user only: ``pip install --user oscar``
* Use sudo (not recommended): ``sudo pip install oscar``

Using Virtual Environments
---------------------------

It's recommended to use a virtual environment to avoid conflicts with other packages.

**Using venv:**

.. code-block:: bash

    # Create virtual environment
    python -m venv oscar-env
    
    # Activate (Linux/Mac)
    source oscar-env/bin/activate
    
    # Activate (Windows)
    oscar-env\Scripts\activate
    
    # Install oscar
    pip install oscar

**Using conda:**

.. code-block:: bash

    # Create conda environment
    conda create -n oscar-env python=3.10
    
    # Activate environment
    conda activate oscar-env
    
    # Install oscar
    pip install oscar

Upgrading
---------

To upgrade to the latest version::

    pip install --upgrade oscar

To check your current version::

    pip show oscar

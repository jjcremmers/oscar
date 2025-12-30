import numpy as np

from .oscar import oscar
   
#-------------------------------------------------------------------------------
#  getFileName
#-------------------------------------------------------------------------------

def getFileName( name : str , prefix : str = None ) -> str:

    '''
    Returns the filename in HDF5 format with proper prefix and extension.
    
    Constructs a filename by adding an optional prefix and ensuring the .h5
    extension is present. If a prefix is provided and the name doesn't start
    with it, the prefix is prepended with an underscore separator.
  
    Args:
        name (str): Base name of the HDF5 file (oscar format)
        prefix (str, optional): Prefix to prepend to the filename. Defaults to None.
        
    Returns:
        str: Formatted filename with prefix and .h5 extension
        
    Example:
        >>> getFileName("mydata", prefix="sim")
        "sim_mydata.h5"
    '''
  
    if prefix != None:
        if not name.startswith(prefix):
            output = prefix + "_" + name
        else:
            output = name
   
    if output.endswith("h5"):
        output = output + ".h5"
  
    return output
  
#-------------------------------------------------------------------------------
#  getConductivity
#-------------------------------------------------------------------------------
      
def getConductivity( baseName : str , props : str = None ) -> float:

    '''
    Returns the thermal conductivity tensor components from RVE simulations.
    
    Reads thermal conductivity data from three separate HDF5 files (one for each
    direction) and returns the conductivity components as a list.
  
    Args:
        baseName (str): Base name of the HDF5 files (oscar format)
        props (str, optional): Additional properties parameter. Defaults to None.
        
    Returns:
        list: Thermal conductivity components [k_x, k_y, k_z]
        
    Note:
        Expects files named "cond_{baseName}_0.h5", "cond_{baseName}_1.h5",
        and "cond_{baseName}_2.h5" to exist.
    '''
  
    cond = []
    
    for i in range(3):
  
        fileName = "cond_" + baseName + "_"+str(i)+".h5"    
    
        h5file  = oscar( fileName )  
                      
        cond.append(h5file.getMacroFieldData("thermalCond")[0])
    
    return cond
  
#-------------------------------------------------------------------------------
#  getCTE
#-------------------------------------------------------------------------------

def getCTE( baseName : str ) -> float:

    '''
    Returns the coefficients of thermal expansion (CTE) as an array.
    
    Reads thermal expansion coefficients from an HDF5 file containing
    RVE homogenization results.
  
    Args:
        baseName (str): Base name of the HDF5 file (oscar format)
        
    Returns:
        list: Coefficients of thermal expansion as a list
        
    Note:
        Expects a file named "cte_{baseName}.h5" to exist with macro field
        data for "thermalExpansion".
    '''

    fileName = "cte_" + baseName +".h5"    
   
    h5file  = oscar( fileName )  
                  
    return h5file.getMacroFieldData("thermalExpansion").tolist()
  
#-------------------------------------------------------------------------------
#  getMechanical
#-------------------------------------------------------------------------------

def getMechanical( baseName : str , eps : float  = 0.01 ) -> list:

    '''
    Returns the mechanical properties from RVE homogenization analysis.
    
    Computes effective mechanical properties (Young's moduli, shear moduli,
    Poisson's ratios, and stiffness matrix) by collecting stress responses
    from six load cases corresponding to unit strains in each direction.
  
    Args:
        baseName (str): Base name of the HDF5 files (oscar format)
        eps (float, optional): Strain magnitude used in the simulations. Defaults to 0.01.
        
    Returns:
        tuple: A tuple containing:
            - E (list): Young's moduli [E_x, E_y, E_z]
            - G (list): Shear moduli [G_xy, G_xz, G_yz]
            - nu (list): Poisson's ratios [nu_xy, nu_yz, nu_xz]
            - cmat (numpy.ndarray): 6x6 stiffness matrix
            
    Note:
        Expects files named "mech_{baseName}_0.h5" through "mech_{baseName}_5.h5"
        corresponding to six independent strain states.
    '''
  
    stress = np.zeros(shape=(6,6))  
    E  = []
    G  = []
    nu = []  
  
    for i in range(6):
  
        fileName = "mech_" + baseName + "_" +str(i)+".h5"
  
        h5file  = oscar( fileName )  
                  
        stress[i,:] = h5file.getMacroFieldData("stresses")
       
    cmat = stress / eps
    smat = np.linalg.inv(cmat)
  
    E.append(1.0/smat[0,0])
    E.append(1.0/smat[1,1])
    E.append(1.0/smat[2,2])
  
    G.append(1.0/smat[3,3])
    G.append(1.0/smat[4,4])
    G.append(1.0/smat[5,5])
  
    nu.append(-smat[0,1]/smat[0,0])
    nu.append(-smat[1,2]/smat[1,1])
    nu.append(-smat[0,2]/smat[2,2])
  
    return E,G,nu,cmat

#-------------------------------------------------------------------------------
#  getCapacity
#-------------------------------------------------------------------------------
  
def getCapacity( baseName : str ) -> float:
  
    '''
    Returns the specific heat capacity from RVE homogenization.
    
    Reads the volumetric average specific heat capacity from an HDF5 file
    containing homogenized thermal properties.
  
    Args:
        baseName (str): Base name of the HDF5 file (oscar format)
        
    Returns:
        float: Specific heat capacity value
        
    Note:
        Expects a file named "capac_{baseName}.h5" to exist with macro field
        data for "specificCapac".
    '''
  
    fileName = "capac_" + baseName + ".h5"
  
    h5file  = oscar( fileName )  
                        
    return h5file.getMacroFieldData("specificCapac")[0]  

#-------------------------------------------------------------------------------
#  getDensity
#-------------------------------------------------------------------------------

def getDensity( baseName : str ) -> float:
  
    '''
    Returns the effective density from RVE homogenization.
    
    Reads the volumetric average density from an HDF5 file containing
    homogenized material properties.
  
    Args:
        baseName (str): Base name of the HDF5 file (oscar format)
        
    Returns:
        float: Density value
        
    Note:
        Expects a file named "capac_{baseName}.h5" to exist with macro field
        data for "density".
    '''
  
    fileName = "capac_" + baseName + ".h5"
  
    h5file  = oscar( fileName )  
                        
    return h5file.getMacroFieldData("density")[0]

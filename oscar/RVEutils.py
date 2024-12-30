import numpy as np

from .oscar import oscar
   
#-------------------------------------------------------------------------------
#  getFileName
#-------------------------------------------------------------------------------

def getFileName( name : str , prefix : str = None ) -> str:

    '''
    Returns the filename in h5 format
  
    Args: name of h5 file (oscar format)
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
    Returns the themral conductivity
  
    Args: name of h5 file (oscar format)
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
    Returns the coefficients of thermal expansion as an array
  
    Args: name of h5 file (oscar format)
    '''

    fileName = "cte_" + baseName +".h5"    
   
    h5file  = oscar( fileName )  
                  
    return h5file.getMacroFieldData("thermalExpansion").tolist()
  
#-------------------------------------------------------------------------------
#  getMechanical
#-------------------------------------------------------------------------------

def getMechanical( baseName : str , eps : float  = 0.01 ) -> list:

    '''
    Returns the mechanical properties by collecting all stresses
    due to a given strain with magnitude eps.
  
    Args: baseName name of h5 file (oscar format)
  
          eps      strain magnitude
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
    Returns the heat capacity
  
    Args: name of h5 file (oscar format)
    '''
  
    fileName = "capac_" + baseName + ".h5"
  
    h5file  = oscar( fileName )  
                        
    return h5file.getMacroFieldData("specificCapac")[0]  

#-------------------------------------------------------------------------------
#  getDensity
#-------------------------------------------------------------------------------

def getDensity( baseName : str ) -> float:
  
    '''
    Returns the coefficients of thermal expansion as an array
  
    Args: name of h5 file (oscar format)
    '''
  
    fileName = "capac_" + baseName + ".h5"
  
    h5file  = oscar( fileName )  
                        
    return h5file.getMacroFieldData("density")[0]

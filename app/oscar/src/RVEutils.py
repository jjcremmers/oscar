import numpy as np

from .oscar import oscar
   
def getFileName( name , prefix = None ):

  if prefix != None:
    if not name.startswith(prefix):
      output = prefix + "_" + name
    else:
      output = name
   
  if output.endswith("h5"):
    output = output + ".h5"
  
  return output
  
#
#
#    
      
def getConductivity( baseName , props = None ):

  cond = []
    
  for i in range(3):
  
    fileName = "cond_" + baseName + "_"+str(i)+".h5"    
    
    h5file  = oscar( fileName )  
    h5file.setCycle(1)
                      
    cond.append(h5file.getMacroFieldData("thermalCond")[0])
    
  return cond
  
#-----------------------
#
#---------------------------------

def getCTE( baseName ):

  fileName = "cte_" + baseName +".h5"    
   
  h5file  = oscar( fileName )  
  h5file.setCycle(1)
                  
  return h5file.getMacroFieldData("thermalExpansion")
  
#--------------------
#
#-------------------------------------------

def getMechanical( baseName , eps = 0.01 ):

  stress = np.zeros(shape=(6,6))  
  E  = []
  G  = []
  nu = []  
  
  for i in range(6):
  
    fileName = "mech_" + baseName + "_" +str(i)+".h5"
  
    h5file  = oscar( fileName )  
    h5file.setCycle(1)
                  
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
  
  return E,G,nu

#---------------------------------------------------
#
#---------------------------------------------------
  
def getCapacity( baseName ):
  
  fileName = "capac_" + baseName + ".h5"
  
  h5file  = oscar( fileName )  
  h5file.setCycle(1)
                        
  return h5file.getMacroFieldData("specificCapac")[0]  

#------------------------------------------
#
#---------------------------------------------------

def getDensity( baseName ):
  
  fileName = "capac_" + baseName + ".h5"
  
  h5file  = oscar( fileName )  
  h5file.setCycle(1)
                        
  return h5file.getMacroFieldData("density")[0]

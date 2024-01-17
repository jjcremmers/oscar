"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2022)
  
"""

import numpy as np

import vtk


def setCellNodes( cell , elemNodes ):

  '''
  
  '''
          
  for i,inod in enumerate(elemNodes):
    cell.GetPointIds().SetId(i,inod)
          
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
          
def insertElement( grid , elemNodes , rank , family ):

  '''
  Inserts an element 
  '''
  
  nNod = len(elemNodes)
  
  if family == 0:   # Continuum
    if rank == 2:
      if nNod == 2:
        cell = vtk.vtkLine()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
      elif nNod == 1:
        print('pass node element')
      elif nNod == 3:
        cell = vtk.vtkTriangle()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )     
      elif nNod == 4:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )       
      elif nNod == 6:
        cell = vtk.vtkTriangle()      
        setCellNodes( cell , elemNodes[0:6:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
      elif nNod == 8:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes[0:8:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
      elif nNod == 9:
        cell = vtk.vtkQuad()      
        setCellNodes( cell , elemNodes[0:8:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
      elif nNod > 1:
        print(nNod)
        raise NotImplementedError('Only 2, 3, 4, 6, 8, 9 node continuum elements in 2D.')     
    elif rank == 3:
      if nNod == 2:
        cell = vtk.vtkLine()    
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() ) 
      elif nNod == 3:
        cell = vtk.vtkLine()    
        setCellNodes( cell , elemNodes[0:3:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
      elif nNod == 4:
        cell = vtk.vtkTetra()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
      elif nNod == 5:
        cell = vtk.vtkPyramid()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
      elif nNod == 6:
        cell = vtk.vtkWedge()      
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )             
      elif nNod == 8:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )               
      elif nNod == 16:
        cell = vtk.vtkHexahedron()
        setCellNodes( cell , numpy.concatenate(elemNodes[0:8:2],elemNodes[8:16:2] ) ) 
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                         
      else:
        raise NotImplementedError('Only 2, 3, 4, 5, 6, 8 and 16 node continuum elements in 3D.')             
    else:
      raise NotImplementedError('Only 2D and 3D continuum elements.')
  elif family == 1:
    if rank == 2:
      if nNod == 4:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes[0:2] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[2:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                             
      else:
        raise NotImplementedError('Only 4 node interface elements in 2D.')           
    elif rank == 3:
      if nNod == 6:
        cell = vtk.vtkTria() 
        setCellNodes( cell , elemNodes[0:3] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[3:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
      elif nNod == 8:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes[0:4] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
        setCellNodes( cell , elemNodes[4:] )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
      else:
        raise NotImplementedError('Only 6 and 8 node interface elements in 3D.')        
    else:
      raise NotImplementedError('Only 2D and 3D interface elements.')        
  elif family == 2:
    if rank == 2:
      if nNod == 2:
        cell = vtk.vtkLine() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
      else:
        raise NotImplementedError('Only 2 node surface elements in 2D.')              
    elif rank == 3:
      if nNod == 3:
        cell = vtk.vtkTriangle() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )              
      elif nNod == 4:
        cell = vtk.vtkQuad() 
        setCellNodes( cell , elemNodes )  
        grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )            
      else:
        raise NotImplementedError('Only 3 and 4 node surface elements in 3D.')              
  elif family == 3:
    if nNod == 2:
      cell = vtk.vtkLine() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )    
    elif nNod == 3:
      cell = vtk.vtkLine() 
      cell.GetPointIds().SetId(0,elemNodes[0])
      cell.GetPointIds().SetId(1,elemNodes[2])      
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )
    else:
      raise NotImplementedError('Only 2 and 3 node beam elements.')            
  elif family == 4:
    if nNod == 3:
      cell = vtk.vtkTriangle() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )        
    elif nNod == 4:
      cell = vtk.vtkQuad() 
      setCellNodes( cell , elemNodes )  
      grid.InsertNextCell( cell.GetCellType(),cell.GetPointIds() )                     
  else:
    raise NotImplementedError('Only 3 and 4 node shell elements.')

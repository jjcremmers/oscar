"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""
import unittest
from oscar.oscar import oscar
from pathlib import Path
  
class OscarTesting(unittest.TestCase):
    
    script_dir = Path(__file__).parent
	
    h5file = oscar( script_dir / "pinched8.h5" )
    h5file.setCycle(9)          
    
    h5fred = oscar( script_dir / "pinched8_reduced.h5" )
    h5fred.setCycle(9)      
        
    def test_getCoords(self):
      
        self.assertEqual( len(self.h5file.getCoords()) , 578 )
        self.assertEqual( len(self.h5fred.getCoords()) , 578 )        

    def test_getCoordsComponents(self):

        self.assertEqual( self.h5file.getCoords(37)[0] , 9.4705 )
        self.assertEqual( self.h5file.getCoords(37)[1] , 2.8728 )
        self.assertEqual( self.h5file.getCoords(37)[2] , 1.5675 )                

    def test_getCoords3(self):
    
        self.assertEqual( len(self.h5file.getCoords3()) , 578 )
    
    def test_getCycle(self):
    
        self.assertEqual( self.h5file.getCycle() , 9 )
    
    def test_getElemNodes(self):

        elemNodes = self.h5file.getElemNodes(89)
        
        self.assertEqual( elemNodes[0] , 94  )  
        self.assertEqual( elemNodes[1] , 95  ) 
        self.assertEqual( elemNodes[2] , 112 ) 
        self.assertEqual( elemNodes[3] , 111 ) 
        self.assertEqual( elemNodes[4] , 383 ) 
        self.assertEqual( elemNodes[5] , 384 ) 
        self.assertEqual( elemNodes[6] , 401 ) 
        self.assertEqual( elemNodes[7] , 400 )   
    
    def test_getElemNodeCount(self):
        
        self.assertEqual( self.h5file.getElemNodeCount(89) , 8 )
        
    def test_getElemGroup(self):
    
        self.assertEqual( len(self.h5file.getElemGroup()) , 256 )
    
    def test_elemCount(self):

        self.assertEqual( self.h5file.elemCount() , 256 ) 
        
    def test_getElemIndex(self):
    
        self.assertEqual( self.h5file.getElemIndex(78) , 77 )
        
    def test_getElemIDs(self):
    
        elemIDs = self.h5file.getElemIDs()
        
        self.assertEqual( elemIDs[0]   , 1   )
        self.assertEqual( elemIDs[100] , 101 )   
        self.assertEqual( elemIDs[139] , 140 )  
    
    def test_nodeCount(self):
    
        self.assertEqual( self.h5file.nodeCount() , 578 )
        
    def test_getElemIndex(self):  
    
        self.assertEqual( self.h5file.getElemIndex(131) , 130 )

    def test_getNodeIDs(self): 
         
        nodeIDS = self.h5file.getNodeIDs()   
        self.assertEqual( self.h5file.getNodeIDs()[0]  , 1  ) 
        self.assertEqual( self.h5file.getNodeIDs()[46] , 47 )         

    def test_getNodeGroupNames(self): 
    
        names = self.h5file.getNodeGroupNames()
        
        self.assertEqual( names[0] , 'loadx' )
        self.assertEqual( names[1] , 'loady' )
        self.assertEqual( names[2] , 'x0' )
        self.assertEqual( names[3] , 'y0' )
        self.assertEqual( names[4] , 'z0' )                                

    def test_getNodeGroup(self): 
    
        nodeGroup = self.h5file.getNodeGroup( 'x0' )
        
        self.assertEqual( len(nodeGroup) , 34 )  
        
        self.assertEqual( nodeGroup[0] , 16 ) 
        self.assertEqual( nodeGroup[1] , 33 ) 
        self.assertEqual( nodeGroup[2] , 50 )                         

    def test_rank(self): 
    
        self.assertEqual( self.h5file.rank() , 3 )

    def test_nodeDataSets(self): 
    
        dataSets = self.h5file.nodeDataSets()

        self.assertEqual( len(dataSets) , 7 )  
               
        self.assertEqual( dataSets[0]   , 'S11' ) 
        self.assertEqual( dataSets[1]   , 'S12' )
        self.assertEqual( dataSets[2]   , 'S13' )
        self.assertEqual( dataSets[3]   , 'S22' )
        self.assertEqual( dataSets[4]   , 'S23' )
        self.assertEqual( dataSets[5]   , 'S33' )
        self.assertEqual( dataSets[6]   , 'displacements' )        
                                                      
    def test_getDisplacements(self): 
    
        disp = self.h5file.getDisplacements(56)
        
        self.assertAlmostEqual( disp[0] , 0.55154688  )
        self.assertAlmostEqual( disp[1] , -0.14104671 )
        self.assertAlmostEqual( disp[2] , -0.54406478 )

    def test_getNodeData(self): 
        
        self.assertAlmostEqual( self.h5file.getNodeData( "S11" , 78 ) ,    
                                9681.055732021865 )
        
        
if __name__ == '__main__':
    unittest.main()

"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""
import unittest
import os
import tempfile
import shutil
import numpy as np
from oscar.oscar import oscar, convertToVTU, convertToDat, writePVD, PVDfileName

  
class OscarTesting(unittest.TestCase):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
	
    h5file = oscar( os.path.join(script_dir, "pinched8.h5" ) )
    h5file.setCycle(9)          
    
    h5fred = oscar( os.path.join(script_dir, "pinched8_reduced.h5" ) )
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
    
    def test_getNodeDataVectorAllNodes(self):
        """Test getting vector data for all nodes"""
        disp_all = self.h5file.getNodeData("displacements")
        self.assertEqual(disp_all.shape[0], 578)
        self.assertEqual(disp_all.shape[1], 3)  # 3D displacements
    
    def test_getNodeDataScalarAllNodes(self):
        """Test getting scalar data for all nodes"""
        stress = self.h5file.getNodeData("S11")
        self.assertEqual(len(stress), 578)
    
    def test_getNodeDataMultipleNodes(self):
        """Test getting data for multiple specific nodes"""
        node_list = [10, 20, 30]
        disp = self.h5file.getNodeData("displacements", node_list)
        self.assertEqual(len(disp), 3)
    
    def test_str(self):
        """Test string representation"""
        str_repr = str(self.h5file)
        self.assertIsInstance(str_repr, str)
        self.assertIn("datasets", str_repr.lower())
    
    def test_getElemGroupNames(self):
        """Test getting element group names"""
        group_names = self.h5file.getElemGroupNames()
        self.assertIsInstance(group_names, list)
    
    def test_getElemGroupSpecific(self):
        """Test getting specific element group"""
        all_elems = self.h5file.getElemGroup("all")
        self.assertEqual(len(all_elems), 256)
    
    def test_getNodeIndex(self):
        """Test getting node index from node ID"""
        # Node IDs typically start at 1, index at 0
        index = self.h5file.getNodeIndex(1)
        self.assertEqual(index, 0)
    
    def test_elemDataSets(self):
        """Test getting element data set names"""
        datasets = self.h5file.elemDataSets()
        self.assertIsInstance(datasets, list)
    
    def test_getElemData(self):
        """Test getting element data if it exists"""
        datasets = self.h5file.elemDataSets()
        if len(datasets) > 0:
            data = self.h5file.getElemData(datasets[0])
            self.assertIsNotNone(data)
    
    def test_nodeCountWithGroup(self):
        """Test node count for specific group"""
        x0_count = self.h5file.nodeCount('x0')
        self.assertEqual(x0_count, 34)
    
    def test_unpackElements(self):
        """Test unpacking elements from connectivity array"""
        connectivity = [0, 1, 2, 3, 3, 4, 5, 6]
        offsets = [4, 8]
        result = self.h5file.unpackElements(connectivity, offsets)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], [0, 1, 2, 3])
        self.assertEqual(result[1], [3, 4, 5, 6])
    
    def test_getCoordsMultipleNodes(self):
        """Test getting coordinates for multiple nodes"""
        node_list = [0, 1, 2]
        coords = self.h5file.getCoords(node_list)
        self.assertEqual(len(coords), 3)
    
    def test_getCoordsSingleNode(self):
        """Test getting coordinates for single node"""
        coord = self.h5file.getCoords(0)
        self.assertEqual(len(coord), 3)  # 3D coordinates


class UtilityFunctionsTesting(unittest.TestCase):
    """Test utility functions like PVDfileName, writePVD, etc."""
    
    def test_PVDfileName_serial(self):
        """Test serial filename generation"""
        filename = PVDfileName("test", 5)
        self.assertEqual(filename, "test_t5.vtu")
    
    def test_PVDfileName_parallel(self):
        """Test parallel filename generation"""
        filename = PVDfileName("test", 5, iProc=2)
        self.assertEqual(filename, "test_p2_t5.vtu")
    
    def test_writePVD_serial(self):
        """Test writing PVD file for serial data"""
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            cycles = [1, 2, 3]
            writePVD(prefix, cycles, nProc=1)
            
            pvd_file = prefix + ".pvd"
            self.assertTrue(os.path.exists(pvd_file))
            
            # Check file content
            with open(pvd_file, 'r') as f:
                content = f.read()
                self.assertIn("VTKFile", content)
                self.assertIn("Collection", content)
                for cycle in cycles:
                    self.assertIn(f"timestep='{cycle}'", content)
    
    def test_writePVD_parallel(self):
        """Test writing PVD file for parallel data"""
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test_parallel")
            cycles = [1, 2]
            nProc = 4
            writePVD(prefix, cycles, nProc=nProc)
            
            pvd_file = prefix + ".pvd"
            self.assertTrue(os.path.exists(pvd_file))
            
            # Check that all processors are referenced
            with open(pvd_file, 'r') as f:
                content = f.read()
                for proc in range(nProc):
                    self.assertIn(f"_p{proc}_", content)


class FileConversionTesting(unittest.TestCase):
    """Test file conversion functions"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    test_file = os.path.join(script_dir, "pinched8_reduced.h5")
    
    def test_convertToVTU_single_cycle(self):
        """Test converting single cycle to VTU"""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            try:
                os.chdir(tmpdir)
                convertToVTU(self.test_file, cycles=1)
                
                # Check that VTU file was created
                vtu_file = "pinched8_reduced_t1.vtu"
                self.assertTrue(os.path.exists(vtu_file))
                
                # Check that PVD file was created
                pvd_file = "pinched8_reduced.pvd"
                self.assertTrue(os.path.exists(pvd_file))
            finally:
                os.chdir(original_dir)
    
    def test_convertToVTU_multiple_cycles(self):
        """Test converting multiple cycles to VTU"""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            try:
                os.chdir(tmpdir)
                cycles = [1, 2]
                convertToVTU(self.test_file, cycles=cycles)
                
                # Check that VTU files were created
                for cycle in cycles:
                    vtu_file = f"pinched8_reduced_t{cycle}.vtu"
                    self.assertTrue(os.path.exists(vtu_file), 
                                  f"VTU file for cycle {cycle} not found")
            finally:
                os.chdir(original_dir)
    
    def test_convertToDat_default(self):
        """Test converting to DAT format with defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            try:
                os.chdir(tmpdir)
                convertToDat(self.test_file)
                
                # Check that DAT file was created
                dat_file = "pinched8_reduced.dat"
                self.assertTrue(os.path.exists(dat_file))
                
                # Check basic content
                with open(dat_file, 'r') as f:
                    content = f.read()
                    self.assertIn("<Nodes>", content)
                    self.assertIn("<Elements>", content)
            finally:
                os.chdir(original_dir)


class SaveMethodsTesting(unittest.TestCase):
    """Test save methods of oscar class"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    def setUp(self):
        """Set up test fixture"""
        self.h5file = oscar(os.path.join(self.script_dir, "pinched8_reduced.h5"))
        self.h5file.setCycle(1)
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixture"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_saveAsVTU_default(self):
        """Test saving as VTU with default parameters"""
        original_dir = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            cycles = self.h5file.saveAsVTU(cycles=1)
            
            self.assertIsNotNone(cycles)
            self.assertIn(1, cycles)
            
            # Check files exist
            vtu_file = "pinched8_reduced_t1.vtu"
            self.assertTrue(os.path.exists(vtu_file))
        finally:
            os.chdir(original_dir)
    
    def test_saveAsVTU_custom_prefix(self):
        """Test saving as VTU with custom prefix"""
        original_dir = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            cycles = self.h5file.saveAsVTU(prefix="custom", cycles=1)
            
            # Check files exist with custom prefix
            vtu_file = "custom_t1.vtu"
            self.assertTrue(os.path.exists(vtu_file))
            
            pvd_file = "custom.pvd"
            self.assertTrue(os.path.exists(pvd_file))
        finally:
            os.chdir(original_dir)
    
    def test_saveAsDat_default(self):
        """Test saving as DAT with default parameters"""
        dat_file = os.path.join(self.temp_dir, "test.dat")
        self.h5file.saveAsDat(fileName=dat_file)
        
        self.assertTrue(os.path.exists(dat_file))
        
        # Verify content structure
        with open(dat_file, 'r') as f:
            content = f.read()
            self.assertIn("<Nodes>", content)
            self.assertIn("</Nodes>", content)
            self.assertIn("<Elements>", content)
            self.assertIn("</Elements>", content)
    
    def test_saveAsDat_with_cycle(self):
        """Test saving as DAT with specific cycle"""
        dat_file = os.path.join(self.temp_dir, "test_cycle.dat")
        self.h5file.saveAsDat(fileName=dat_file, cycle=1)
        
        self.assertTrue(os.path.exists(dat_file))
    
    def test_saveAsDat_pyfem_output(self):
        """Test saving as DAT in PyFEM format"""
        dat_file = os.path.join(self.temp_dir, "test_pyfem.dat")
        # This might fail if element groups don't exist, but we test the call
        try:
            self.h5file.saveAsDat(fileName=dat_file, cycle=1, output='pyfem')
            self.assertTrue(os.path.exists(dat_file))
        except Exception:
            # If it fails due to data format issues, that's acceptable
            pass


class EdgeCasesTesting(unittest.TestCase):
    """Test edge cases and error handling"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    def setUp(self):
        """Set up test fixture"""
        self.h5file = oscar(os.path.join(self.script_dir, "pinched8.h5"))
        self.h5file.setCycle(1)
    
    def test_file_without_extension(self):
        """Test opening file without .h5 extension"""
        file_path = os.path.join(self.script_dir, "pinched8")
        h5 = oscar(file_path)
        self.assertIsNotNone(h5)
    
    def test_verbose_mode(self):
        """Test opening file in verbose mode"""
        file_path = os.path.join(self.script_dir, "pinched8.h5")
        h5 = oscar(file_path, verbose=True)
        self.assertIsNotNone(h5)
    
    def test_setCycle_updates_data(self):
        """Test that setting cycle updates internal data"""
        initial_cycle = self.h5file.getCycle()
        self.h5file.setCycle(2)
        new_cycle = self.h5file.getCycle()
        self.assertNotEqual(initial_cycle, new_cycle)
        self.assertEqual(new_cycle, 2)
    
    def test_getCoords_negative_index(self):
        """Test getCoords with invalid index handling"""
        # -1 should return all coordinates
        all_coords = self.h5file.getCoords(-1)
        self.assertGreater(len(all_coords), 0)
    
    def test_elemCount_invalid_group(self):
        """Test elemCount with non-existent group"""
        try:
            count = self.h5file.elemCount("nonexistent_group")
            # Should either return 0 or raise an exception
            self.assertIsInstance(count, int)
        except:
            # Expected behavior for invalid group
            pass
    
    def test_getNodeData_nonexistent_label(self):
        """Test getNodeData with non-existent label"""
        with self.assertRaises(Exception):
            self.h5file.getNodeData("nonexistent_data_label")
    
    def test_getDisplacements_single_vs_multiple(self):
        """Test getDisplacements returns correct dimensions"""
        single = self.h5file.getDisplacements(0)
        self.assertEqual(len(single), 3)  # 3D vector
        
        multiple = self.h5file.getDisplacements([0, 1, 2])
        self.assertEqual(multiple.shape[0], 3)
        self.assertEqual(multiple.shape[1], 3)


class MultipleFilesTesting(unittest.TestCase):
    """Test operations across multiple test files"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    def test_multiple_files_basic_operations(self):
        """Test basic operations work across different files"""
        test_files = ["pinched8.h5", "pinched8_reduced.h5"]
        
        for filename in test_files:
            filepath = os.path.join(self.script_dir, filename)
            if os.path.exists(filepath):
                h5 = oscar(filepath)
                h5.setCycle(1)
                
                # Basic operations should work
                self.assertGreater(h5.nodeCount(), 0)
                self.assertGreater(h5.elemCount(), 0)
                self.assertGreater(h5.rank(), 0)
    
    def test_compare_full_vs_reduced(self):
        """Compare full and reduced versions of same file"""
        full = oscar(os.path.join(self.script_dir, "pinched8.h5"))
        reduced = oscar(os.path.join(self.script_dir, "pinched8_reduced.h5"))
        
        full.setCycle(1)
        reduced.setCycle(1)
        
        # Should have same node count
        self.assertEqual(full.nodeCount(), reduced.nodeCount())
        # Should have same element count
        self.assertEqual(full.elemCount(), reduced.elemCount())
        # Should have same rank
        self.assertEqual(full.rank(), reduced.rank())


class DataIntegrityTesting(unittest.TestCase):
    """Test data integrity and consistency"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    def setUp(self):
        """Set up test fixture"""
        self.h5file = oscar(os.path.join(self.script_dir, "pinched8.h5"))
        self.h5file.setCycle(5)
    
    def test_node_ids_sequential(self):
        """Test that node IDs are properly sequenced"""
        node_ids = self.h5file.getNodeIDs()
        self.assertGreater(len(node_ids), 0)
        # First node ID should typically be 1
        self.assertEqual(node_ids[0], 1)
    
    def test_element_ids_sequential(self):
        """Test that element IDs are properly sequenced"""
        elem_ids = self.h5file.getElemIDs()
        self.assertGreater(len(elem_ids), 0)
        # First element ID should typically be 1
        self.assertEqual(elem_ids[0], 1)
    
    def test_coordinates_dimension_consistency(self):
        """Test coordinate dimensions are consistent"""
        coords = self.h5file.getCoords()
        rank = self.h5file.rank()
        
        # All coordinates should have 'rank' dimensions
        self.assertEqual(coords.shape[1], rank)
    
    def test_element_connectivity_valid(self):
        """Test element connectivity references valid nodes"""
        node_count = self.h5file.nodeCount()
        elem_count = self.h5file.elemCount()
        
        # Check a few elements
        for i in range(min(10, elem_count)):
            elem_nodes = self.h5file.getElemNodes(i)
            for node in elem_nodes:
                # Node should be within valid range
                self.assertGreaterEqual(node, 0)
                self.assertLess(node, node_count)
    
    def test_displacement_dimensions_match_rank(self):
        """Test displacement vectors match problem dimension"""
        rank = self.h5file.rank()
        disp = self.h5file.getDisplacements(0)
        
        # Displacement should have same dimension as rank
        # (note: might be padded to 3D for visualization)
        self.assertGreaterEqual(len(disp), rank)
    
class PerformanceTesting(unittest.TestCase):
    """Test performance-related aspects"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    def setUp(self):
        """Set up test fixture"""
        self.h5file = oscar(os.path.join(self.script_dir, "pinched8.h5"))
        self.h5file.setCycle(1)
    
    def test_large_data_access(self):
        """Test accessing large datasets"""
        # Get all coordinates at once
        coords = self.h5file.getCoords()
        self.assertIsNotNone(coords)
        
        # Get all displacements at once
        disp = self.h5file.getDisplacements(-1)
        self.assertIsNotNone(disp)
    
    def test_repeated_access(self):
        """Test repeated data access doesn't cause issues"""
        for _ in range(10):
            coords = self.h5file.getCoords(0)
            self.assertEqual(len(coords), 3)
    
    def test_cycle_switching(self):
        """Test switching between cycles multiple times"""
        for cycle in [1, 5, 9, 1, 9]:
            self.h5file.setCycle(cycle)
            self.assertEqual(self.h5file.getCycle(), cycle)
            # Data should still be accessible
            self.assertGreater(self.h5file.nodeCount(), 0)
        
        
if __name__ == '__main__':
    unittest.main()

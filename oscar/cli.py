#!/usr/bin/env python
"""
Command-line interface for Oscar HDF5 to VTU conversion.

This module provides CLI tools to convert Oscar HDF5 files to VTK Unstructured
Grid (VTU) format for visualization in ParaView or other VTK-compatible software.
"""

from .oscar import oscar, writePVD
import argparse
import os
import sys


def main():
    """Main CLI entry point for h5tovtk converter.
    
    Parses command-line arguments and converts Oscar HDF5 files to VTU format.
    Supports both single files and multi-processor outputs with automatic
    prefix detection.
    
    Usage:
        h5tovtk <filename> [options]
        
    Examples:
        # Convert entire file
        h5tovtk simulation.h5
        
        # Convert specific cycle
        h5tovtk simulation.h5 -c 5
        
        # Convert multi-processor files
        h5tovtk simulation
        
        # With verbose output
        h5tovtk simulation.h5 -v
    """
    
    parser = argparse.ArgumentParser(
        description='h5tovtk: A tool to convert Oscar H5 files to VTK files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s simulation.h5              Convert all cycles
  %(prog)s simulation.h5 -c 5         Convert cycle 5
  %(prog)s simulation_p0.h5           Convert multi-processor files
  %(prog)s simulation -v              Verbose output
        '''
    )
    
    parser.add_argument('filename', 
                        help='HDF5 file name (with or without .h5 extension)')
    
    parser.add_argument('-c', '--cycle', 
                        type=int, 
                        default=-1, 
                        metavar='N',
                        help='Specific cycle number to convert (default: all cycles)')
    
    parser.add_argument('-p', '--prefix',
                        type=str,
                        default='None',
                        help='Output prefix for VTU files (default: use input filename)')
    
    parser.add_argument('-v', '--verbose', 
                        help='Increase output verbosity',
                        action='store_true')
    
    args = parser.parse_args()
    fileName = args.filename
    nProc = 0
    
    # Handle file name and check for multi-processor files
    if not fileName.endswith('.h5'):
        if os.path.exists(fileName + '_p0.h5'):
            prefix = fileName
            while os.path.exists(fileName + '_p' + str(nProc) + '.h5'):
                nProc = nProc + 1
            
            if args.verbose:
                print(f"Processing multi-processor files with prefix '{prefix}' "
                      f"for {nProc} processors")
        else:
            fileName = fileName + '.h5'
            if args.verbose:
                print(f"Processing file: {fileName}")
    else:
        if args.verbose:
            print(f"Processing file: {fileName}")
    
    # Check if file exists
    if nProc == 0 and not os.path.exists(fileName):
        print(f"Error: File '{fileName}' not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Single file processing
        if nProc == 0:
            if args.verbose:
                print(f"Opening file: {fileName}")
            
            h5file = oscar(fileName)
            
            # Convert specified cycle(s)
            if args.cycle == -1:
                if args.verbose:
                    print(f"Converting all cycles ({h5file.cycleCount} total)")
                cycles = h5file.saveAsVTU(prefix=args.prefix)
            else:
                if args.verbose:
                    print(f"Converting cycle {args.cycle}")
                cycles = h5file.saveAsVTU(prefix=args.prefix, cycles=[args.cycle])
            
            if args.verbose:
                print(f"Successfully converted {len(cycles)} cycle(s)")
        
        # Multi-processor file processing
        else:
            all_cycles = None
            for iProc in range(nProc):
                proc_file = prefix + '_p' + str(iProc) + '.h5'
                
                if args.verbose:
                    print(f"Processing processor {iProc}/{nProc-1}: {proc_file}")
                
                h5file = oscar(proc_file)
                
                if args.cycle == -1:
                    cycles = h5file.saveAsVTU()
                else:
                    cycles = h5file.saveAsVTU(cycles=[args.cycle])
                
                if all_cycles is None:
                    all_cycles = cycles
            
            # Create combined PVD file for multi-processor output
            if all_cycles is not None:
                if args.verbose:
                    print(f"Creating combined PVD file for {nProc} processors")
                writePVD(prefix, all_cycles, nProc)
            
            if args.verbose:
                print(f"Successfully converted {len(all_cycles)} cycle(s) "
                      f"for {nProc} processors")
        
        if args.verbose:
            print("Conversion complete!")
            
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

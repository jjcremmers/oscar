from oscar import oscar
import os
  
#----- Start

script_dir = os.path.dirname(os.path.abspath(__file__))

h5file  = oscar( os.path.join(script_dir, "pinched8.h5" ) )
h5file.saveAsVTU("pinched8")

h5file  = oscar( os.path.join(script_dir, "pinched8_reduced.h5" ) )
h5file.saveAsVTU("pinched8_reduced")

h5file  = oscar( os.path.join(script_dir, "pinched8_stateonly.h5" ) )
h5file.saveAsVTU()



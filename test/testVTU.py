from oscar import oscar
from pathlib import Path
  
#----- Start

script_dir = Path(__file__).parent 

h5file  = oscar( script_dir / "pinched8.h5" )
h5file.saveAsVTU("pinched8")

h5file  = oscar( script_dir / "pinched8_reduced.h5" )
h5file.saveAsVTU("pinched8_reduced")

h5file  = oscar( script_dir / "pinched8_stateonly.h5" )
h5file.saveAsVTU()



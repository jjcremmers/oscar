from oscar import oscar
  
#----- Start

h5file  = oscar( "pinched8.h5" )
h5file.saveAsVTU("pinched8")

h5file  = oscar( "pinched8_reduced.h5" )
h5file.saveAsVTU("pinched8_reduced")

h5file  = oscar( "pinched8_stateonly.h5" )
h5file.saveAsVTU()



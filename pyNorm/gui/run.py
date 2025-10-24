from astropy.io import fits
from pynorm.gui import Absorber_pn as A
from pynorm.gui import Metal_Plot_pn as M   
from linetools.spectra.xspectrum1d import XSpectrum1D

# Read in the 1D spectrum to be analyzed
filename='SDSS2338+1504_f.fits'

a=XSpectrum1D.from_file(filename)
wave=a.wavelength.value
flux=a.flux.value
error=a.sig.value
#------------------------------

instrument = 'HIRES'
#Specify redshift at which to perform analysis
z=2.23275

#Give approximate absorption line rest frame wavelengths to be analyzed
lines = [1036.3367,1334,1191,1193.29,1206,1394,1402,1549,1550,1304,1808]


#Preprocessing:
# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines)   
Abs=absys.ions
M.Transitions(Abs, instrument=instrument)

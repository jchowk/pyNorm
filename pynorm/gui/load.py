import pickle
pfile='Spectrum_Analysis_z_2.23275.p'
with open(pfile,'rb') as pklfile:
    absys=pickle.load(pklfile)

#Run the Master GUI
from pynorm.gui import Metal_Plot_pn as M   
M.Transitions(absys)

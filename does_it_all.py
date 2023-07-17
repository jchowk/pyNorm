from pyNorm.io import read_rbcodes
from pyNorm.aod import pyn_batch
from pyNorm.aod import xlimit
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import glob
fi = input("ENTER SIGHTLINE:")                  #user input for a given sightline directory
ra = float(input("ra of QSO:"))
dec = float(input("dec of QSO:"))
redshift = input("ENTER REDSHIFT OF THE ABSORBER:")
fil = fi+"/Spectrum_Analysis_z_"+redshift+"*.p"              #string that contains the address
file_list = glob.glob(fil)                      #obtains all frequency files and creates a list of the files
tarname = fi
#ra = 209.6324424188400 
#dec = 05.0896574614700
spec = []
spec_in = []
ion_list=[]
print(fil)
pdf_pages_flux = PdfPages(fi+'/'+fi+'_'+redshift+'_Nav.pdf')
pdf_pages_blem = PdfPages(fi+'/'+fi+'_'+redshift+'_blemish.pdf')
for file in file_list:
    print(file)
    with (open(file, "rb")) as openfile:
        while True:
            try:
                spec_in=pickle.load(openfile)
                spec.append(spec_in)
            except EOFError:
                break
# reading ions in the list
    print("Redshift of the absorber:",spec_in['Target']['z'])
    keys_dict = list(spec_in.keys())
    ion_in_file = keys_dict[0:len(keys_dict)-1]
    print('Ions in file '+file+' are: ',ion_in_file) 
    for j in ion_in_file:
        spec_obj = read_rbcodes(file,fi,ra,dec,j)
        if(spec_obj['flag_blemish']==True):
            fig = plt.figure(figsize=(10,8))
            plt.title(j)
 #           plt.plot(spec_obj['vel'],spec_obj['flux']/spec_obj['contin'],linewidth=2,color='blue')
            plt.plot(spec_obj['vel'],spec_obj['fnorm'],alpha=0.5,color='red',label='Uncorrected')
            plt.plot(spec_obj['vel'],spec_obj['eflux'],color='black')
            plt.axvline(spec_obj['v1'],linestyle='--',color='blue')
            plt.axvline(spec_obj['v2'],linestyle='--',color='blue')
            plt.scatter(spec_obj['vel'][spec_obj['eflux']==-0.9],spec_obj['flux'][spec_obj['eflux']==-0.9]/spec_obj['contin'][spec_obj['eflux']==-0.9],color='red')
            plt.axhline(0.0,color='black')
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('Normalized flux')
            plt.legend()
            pdf_pages_blem.savefig(fig)
            fig = plt.figure(figsize=(10,8))
            plt.title(j)
            plt.plot(spec_obj['vel'],spec_obj['flux']/spec_obj['contin'],color='blue',label='Corrected')
#            plt.plot(spec_obj['vel'],spec_obj['fnorm'],alpha=0.5,color='red')
            plt.plot(spec_obj['vel'],spec_obj['eflux'],color='black')
            plt.axvline(spec_obj['v1'],linestyle='--',color='blue')
            plt.axvline(spec_obj['v2'],linestyle='--',color='blue')
            plt.scatter(spec_obj['vel'][spec_obj['eflux']==-0.9],spec_obj['flux'][spec_obj['eflux']==-0.9]/spec_obj['contin'][spec_obj['eflux']==-0.9],color='red')
            plt.axhline(0.0,color='black')
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('Normalized flux')
            plt.legend()
            pdf_pages_blem.savefig(fig)
        fig = plt.figure(figsize=(10,8))
        plt.title(j)
        plt.plot(spec_obj['vnorm'],spec_obj['Nav'],drawstyle='steps-mid')
        plt.axhline(0,linestyle=':',color='k',linewidth=1,zorder=0)
        plt.axhline(1,linestyle=':',color='k',linewidth=1,zorder=0)
        plt.axvline(spec_obj['v1'],linestyle='--',color='blue')
        plt.axvline(spec_obj['v2'],linestyle='--',color='blue')
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('$N_a(v)$ (cm$^{-2}$ (km/s)$^{-1}$)');
        plt.xlim(-250,250)
        plt.legend()
        pdf_pages_flux.savefig(fig)
    ion_list.extend(keys_dict[0:len(keys_dict)-1])
pdf_pages_blem.close()
pdf_pages_flux.close()
print(ion_list)
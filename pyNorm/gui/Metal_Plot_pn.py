import numpy as np
from numpy.polynomial.legendre import legfit, legval
import numpy.polynomial.legendre as L

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.constants import c

from matplotlib.figure import Figure
from matplotlib import rcParams
import matplotlib.pyplot as plt
import pickle
import sys
from collections import OrderedDict
import pandas as pd
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QStyleFactory, QPushButton,QLineEdit,QMainWindow,QInputDialog,QLabel,QMessageBox,QScrollBar,QVBoxLayout,QHBoxLayout
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QPalette, QColor
import traceback
from PyQt5.QtWidgets import QMessageBox
from . import Absorber_pn
from matplotlib.axes import Axes
from pyNorm.aod import pyn_batch
from pyNorm.continuum import continuum_fit

rcParams['lines.linewidth'] = .9

USE_COMPUTE_EW = False  # Turn to True to re-enable compute_EW (Rongmon's function)

def handle_exception(exc_type, exc_value, exc_traceback):
    if not issubclass(exc_type, KeyboardInterrupt):
        print("Uncaught exception:", exc_type, exc_value)
        
        tb_str = ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))
        QMessageBox.critical(None, "Python Error", tb_str)


def compute_EW(lam,flx,wrest,lmts,flx_err,plot=False,**kwargs):
    """
    ------------------------------------------------------------------------------------------
       Function to compute the equivalent width within a given velocity limits lmts=[vmin,vmax]
               [Only good for high resolution spectra]
      Caveats:- Not automated, must not include other absorption troughs within the velocity range.
     
       Parameters:- 
               lam         :- Observed Wavelength vector (units of Angstrom)
               flx         :- flux vector ( same length as wavelgnth vector, preferably continuum normalized)
               wrest       :- rest frame wavelength of the line [used to make velcity cuts]
               lmts        :- [vmin,vmax], the velocity window within which equivalent width is computed.
               flx_err     :- error spectrum [same length as the flux vector]
    
       OPTIONAL :-
               f0=f0       :- fvalue of the transition 
               zabs=zabs   :- absorber redshift
               plot        :- plot keyword, default = no plots plot=0
                               plot=1 or anything else will plot the corresponding spectrum 
                                and the apparent optical depth of absorption. 
    
    
    
        Returns:-  In a Python dictionary format
               output['ew_tot']      :- rest frame equivalent width of the absorpiton system [Angstrom]
               output['err_ew_tot']  :- error on rest fram equivalent width 
               output['col']         :- AOD column denisty 
               output['colerr']      :- 1 sigma error on AOD column density 
               output['n']           :- AOD column density as a function of velocity
               output['Tau_a']       :- AOD as a function of velocity
               output['med_vel']     :- velocity centroid (Median Equivalent Width weighted velocity within lmts)
               output['vel_disp']    : 1 sigma velocity dispersion
               output['vel50_err']   : error on velocity centroid
    
    
       Written :- Rongmon Bordoloi                             2nd November 2016
    -  I translated this from my matlab code compute_EW.m, which in turn is from Chris Thom's eqwrange.pro. 
       This was tested with COS-Halos/Dwarfs data. 
       Edit:  RB July 5 2017. Output is a dictionary. Edited minor dictionary arrangement
              RB July 25 2019. Added med_vel
              RB April 28, 2021, changed med_vel to weight be EW & vel_disp
    ------------------------------------------------------------------------------------------
    """
    defnorm=1.0;
    spl=2.9979e5;  #speed of light
    if 'zabs' in kwargs:
        zabs=kwargs['zabs']
    else:
        zabs=0.

    if 'sat_limit' in kwargs:
        sat_limit=kwargs['sat_limit']
    else:
        sat_limit=0.10 #  Limit for saturation (COS specific). Set to same as fluxcut for now. WHAT SHOULD THIS BE???
    vel = (lam-wrest*(1.0 + zabs))*spl/(wrest*(1.0 + zabs));
    lambda_r=lam/(1.+zabs);

    

    norm=defnorm

    norm_flx=flx/norm;
    flx_err=flx_err/norm;
    sq=np.isnan(norm_flx);
    tmp_flx=flx_err[sq]
    norm_flx[sq]=tmp_flx
    #clip the spectrum. If the flux is less than 0+N*sigma, then we're saturated. Clip the flux array(to avoid inifinite optical depth) and set the saturated flag
    q=np.where(norm_flx<=sat_limit);
    tmp_flx=flx_err[q]
    norm_flx[q]=tmp_flx
    q=np.where(norm_flx<=0.);
    tmp_flx=flx_err[q]+0.01
    norm_flx[q]=tmp_flx;


    del_lam_j=np.diff(lambda_r);
    del_lam_j=np.append([del_lam_j[0]],del_lam_j);


    pix = np.where( (vel >= lmts[0]) & (vel <= lmts[1]));
    Dj=1.-norm_flx

    # Equivalent Width Per Pixel
    ew=del_lam_j[pix]*Dj[pix];


    sig_dj_sq=(flx_err)**2.;
    err_ew=del_lam_j[pix]*np.sqrt(sig_dj_sq[pix]);
    err_ew_tot=np.sqrt(np.sum(err_ew**2.));
    ew_tot=np.sum(ew);

    #compute the velocity centroid of ew weighted velcity.
    ew50=np.cumsum(ew)/np.max(np.cumsum(ew))
    vel50=np.interp(0.5,ew50,vel[pix])
    vel16=np.interp(0.16,ew50,vel[pix])
    vel_disp=np.abs(vel50-vel16)
    vel50_err = vel_disp/np.sqrt(len(ew))



    print('W_lambda = ' + str('%.3f' % ew_tot) + ' +/- ' + str('%.3f' % err_ew_tot)  +'  \AA   over [' + str('%.1f' % np.round(lmts[0]))+' to ' +str('%.1f' % np.round(lmts[1])) + ']  km/s')
    output={}
    output["ew_tot"]=ew_tot
    output["err_ew_tot"]=err_ew_tot
    output["vel_disp"]=vel_disp
    output['vel50_err']=vel50_err


    if 'f0' in kwargs:
        f0=kwargs['f0']
        #compute apparent optical depth
        Tau_a =np.log(1./norm_flx);
        

        # REMEMBER WE ARE SWITCHING TO VELOCITY HERE
        del_vel_j=np.diff(vel);
        del_vel_j=np.append([del_vel_j[0]],del_vel_j)
        
        # Column density per pixel as a function of velocity
        nv = Tau_a/((2.654e-15)*f0*lambda_r);# in units cm^-2 / (km s^-1), SS91 
        n = nv* del_vel_j# column density per bin obtained by multiplying differential Nv by bin width 
        tauerr = flx_err/norm_flx;
        nerr = (tauerr/((2.654e-15)*f0*lambda_r))*del_vel_j; 
        col = np.sum(n[pix]);
        colerr = np.sum((nerr[pix])**2.)**0.5; 
        print('Direct N = ' + str('%.3f' % np.log10(col))  +' +/- ' + str('%.3f' % (np.log10(col+colerr) - np.log10(col))) + ' cm^-2')
        output["col"]=col
        output["colerr"]=colerr
        output["Tau_a"]=Tau_a
        output["med_vel"]=vel50
        
    # If plot keyword is  set start plotting
    if plot is not False:
        fig = plt.figure()
        ax1=fig.add_subplot(211)
        ax1.step(vel,norm_flx)
        ax1.step(vel,flx_err,color='r')
        #plt.xlim([lmts[0]-2500,lmts[1]+2500])
        plt.xlim([-600,600])
        plt.ylim([-0.02,1.8])
        ax1.plot([-2500,2500],[0,0],'k:')
        ax1.plot([-2500,2500],[1,1],'k:')       
        plt.plot([lmts[0],lmts[0]],[1.5,1.5],'r+',markersize=15)        
        plt.plot([lmts[1],lmts[1]],[1.5,1.5],'r+',markersize=15)    
        plt.title(r' $W_{rest}$= ' + str('%.3f' % ew_tot) + ' $\pm$ ' + str('%.3f' % err_ew_tot) + ' $\AA$')
        ax1.set_xlabel('vel [km/s]')
    
        ax2=fig.add_subplot(212)
        ax2.step(vel,n)
        ax2.set_xlabel('vel [km/s]')
        ax2.plot([-2500,2500],[0,0],'k:')
        #plt.xlim([lmts[0]-2500,lmts[1]+2500])
        plt.xlim([-600,600])
        plt.show()

    
    return output

def rb_set_color():
    """
    Defines a set of colors
       
       Parameters 
       ----------
       None

       Returns
       -------
       clr: a dictionary with different colors
    """
    clr={}
    clr['black']=[  0.,   0.,   0.]
    clr['white'] =[ 255./ 255., 255./ 255., 255./ 255.]
    clr['red'] =[ 255./ 255.,   0./ 255.,   0./ 255.]
    clr['green'] =[   0., 255./ 255.,   0.]
    clr['blue'] =[   0.,   0., 255./ 255.]
    clr['dark_orange']=[1., 0.5, 0.2]
    clr['orange'] =[ 230./ 255., 159./ 255.,   0./ 255.]
    clr['sky_blue'] =[  86./ 255., 180./ 255., 233./ 255.]
    clr['bluish_green'] =[   0./ 255., 158./ 255., 115./ 255.]
    clr['yellow'] =[ 240./ 255., 228./ 255.,  66./ 255.]
    clr['blue2'] =[   0./ 255., 114./ 255., 178./ 255.]
    clr['vermillion'] =[ 213./ 255.,  94./ 255.,   0.]
    clr['reddish_purple'] =[ 204./ 255., 121./ 255., 167./ 255.] 
    clr['cream'] =[ 248./ 255., 248./ 255., 248./ 255.]
    clr['cyan'] =[   0./ 255., 255./ 255., 255./ 255.]
    clr['light_lime_green'] =[ 153./ 255., 255./ 255., 153./ 255.]
    clr['pale_lime_green'] =[ 200./ 255., 255./ 255., 200./ 255.]
    clr['purple_wordle'] =[ 204./ 255.,   0., 255./ 255.]
    clr['light_purple'] =[ 231./ 255., 113./ 255., 255./ 255.]
    clr['orange2'] =[ 217./ 255., 110./ 255.,   0./ 255.]
    clr['light_orange'] =[ 233./ 255., 172./ 255.,  55./ 255.]
    clr['teal'] =[  62./ 255., 133./ 255., 181./ 255.]
    clr['pale_red'] =[ 255./ 255., 118./ 255., 110./ 255.]
    clr['pale_cyan'] =[ 133./ 255., 249./ 255., 255./ 255.]
    clr['dark_red'] =[ 100./ 255.,   0.,   0.]
    clr['dark_green'] =[   0., 100./ 255.,   0.] 
    clr['dark_blue'] =[   0.,   0., 100./ 255.]
    clr['gray']=[.5,.5,.5]
    clr['light_gray']=[.8,.8,.8]

    return clr

clr=rb_set_color()

HELP =  '''
        ---------------------------------------------------------------------------
        pyNorm interactive 1D absorption line measurement toolbox.
        This allows for interactive continuum fitting and equivalent width measurement of CGM/IGM/ISM absorption lines. [Based on rb_codes written by Rongmon Bordoloi.]

        ---------------------
        Screen Layout:
        ---------------------
            LHS/RHS = Left Hand Side/Right Hand Side
            
            LHS shows spectrum with overlaid legendre polynomial continuum estimate; grayed regions indicate masked regions.
            
            RHS shows normalized spectrum with velocity limits

        ---------------------    
        Mouse button Layout:
        ---------------------
            LMB/RMB = Left Mouse Button/Right Mouse Button

            
        ---------------------
        Useful Mouse Clicks:
        ---------------------
        *LHS:
            LMB     : Remove wavelengths from continuum fit (click left/right).
            RMB     : Add wavelengths to continuum fit (click left/right).
        *RHS
            LMB     : Set lower velocity limit
            RMB     : Set upper velocity limit
            
        ---------------------
        Useful Keystrokes:            
        ---------------------

            v           : place mouse on desired subplot
                                   LHS: manually enter regions to mask continuum 
                                   RHS: manually enter EW intergration limits
            
            V (RHS only): Use active subplot velocity limits for all lines
            
            Up arrow    : Increase Polynomial Order [default 4]
            Down arrow  : Decrease Polynomial Order [default 4]

            m           : Measure EW/N for active subplot
            M           : Measure EW/N for ALL subplots
            x/X         : Zoom in/out along x-axis
            y/Y         : Zoom in/out along y-axis
            [,]         : Move left and right in velocity on LHS
            w,s         : Move up and down in flux on LHS
            W,S         : Move up and down in flux on LHS by larger steps

            1/2/0 (RHS only): flag absorber as
                              (0) positive detection
                              (1) upper limit 
                              (2) lower limit

            t (RHS only): Cycle text printed on absorbers. Displays logN, or EW
            
            Q           : Exit the GUI (close all windows)
            ------------------------------------------------------------------------------
            Each tab displays up to 6 transitions. There are maximum 5 tabs allowed. This limits the total number of transitions that can be simultaneously analyzed to 30.
            '''

# Autocontinuum

def velocity(wave, w0):
    return (wave - w0) / w0 * c / 1e3

def fit_legendre(x, y, order, mask=None):
    if mask is not None:
        x_fit, y_fit = x[mask], y[mask]
    else:
        x_fit, y_fit = x, y
    coeffs = legfit(x_fit, y_fit, order)
    model = legval(x, coeffs)
    return coeffs, model

def estimate_continuum(xdata, ydata, yerr, wline, zabs, vwin, order):
    wobs = wline * (1 + zabs)
    v = velocity(xdata, wobs)
    in_window = (v > -vwin) & (v < vwin)

    xcut = xdata[in_window]
    ycut = ydata[in_window]
    yerrcut = np.abs(yerr[in_window])
    xnorm = 2 * (xcut - xcut.min()) / (xcut.max() - xcut.min()) - 1

    good = (ycut / yerrcut > 3)
    coeffs1, cont1 = fit_legendre(xnorm, ycut, order, mask=good)
    residual = np.abs(ycut - cont1)

    bad_resid = (residual >= 2.5 * yerrcut)
    padded = np.pad(bad_resid.astype(int), (1, 1), constant_values=(0, 0))
    diff = np.diff(padded)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]

    absorption_mask = np.ones_like(residual, dtype=bool)
    for s, e in zip(starts, ends):
        if (e - s) > 10:
            s_clip = max(s - 20, 0)
            e_clip = min(e + 20, len(residual))
            absorption_mask[s_clip:e_clip] = False

    good2 = (ycut / yerrcut > 3) & absorption_mask
    coeffs2, final_cont = fit_legendre(xnorm, ycut, order, mask=good2)
    return xcut, ycut, yerrcut, final_cont, coeffs2, good2, absorption_mask

# --- Function to compute Na(v) from in-memory data --- #
def compute_nav_profile(ion_data):
    spec = OrderedDict()

    vel = ion_data['vel']
    wave = ion_data['wave']
    flux = ion_data['flux']
    error = ion_data['error']
    cont = ion_data['cont']
    zabs = ion_data['z']
    fval = ion_data['f']
    lam0 = ion_data['lam_0']
    EWlims = ion_data['EWlims']

    norm_flux = flux / cont
    norm_error = error / cont
    wavelength = lam0 * (vel / 2.998e5) + lam0

    spec['vel'] = vel
    spec['flux'] = flux
    spec['eflux'] = error
    spec['wave'] = wavelength
    spec['contin'] = cont
    spec['contin_err'] = np.zeros_like(cont)
    spec['contin_order'] = ion_data.get('order', 0)
    spec['fval'] = fval
    spec['wavc'] = lam0
    spec['z'] = zabs
    spec['ion'] = ion_data['name']
    spec['wni'] = str(lam0)
    spec['v1'] = EWlims[0]
    spec['v2'] = EWlims[1]
    spec['object'] = 'QSO'

    spec['RA'], spec['Dec'] = 0.0, 0.0
    coords = SkyCoord(0.0, 0.0, unit="deg", frame="icrs")
    spec['gl'] = coords.galactic.l.value
    spec['gb'] = coords.galactic.b.value

    spec['contin_coeff'] = np.array([0.0, 0.0])
    spec['mask_cont'] = np.ones_like(vel)
    spec['contin_mask_bits'] = np.ones_like(vel)

    spec = pyn_batch(spec, verbose=False, partial_pixels=True, blemish_correction=True)
    #spec = continuum_fit(spec, minord=spec['contin_order'], maxord=spec['contin_order'])
    if spec['flag_sat']:
        flag = -2
    elif not spec['detection_3sig']:
        flag = -1
    else:
        flag = 0

    spec['pyn_flag'] = flag  # Add flag to spec for storage
    return spec

# --- Na(v) Panel Plotting Function --- #
def plot_nav_panel(parent, key_idx, ii):
    ion_data = parent.ions[parent.keys[key_idx]]

    nav = ion_data['Nav']
    nav_err = ion_data['Nav_err']
    nav_vel = ion_data['Nav_vel']
    vmin = ion_data['EWlims'][0]
    vmax = ion_data['EWlims'][1]

    ax = parent.axesN[parent.page][ii]
    ax.clear()

    if ii == 0:
        ax.set_title('Na(v) Profile')
        ax.set_xlabel('Velocity (km/s)')

    ax.step(nav_vel, nav, where='mid', color='purple')
    ax.fill_between(nav_vel, nav - nav_err, nav + nav_err, color='violet', alpha=0.3)

    ax.set_xlim(ion_data['window_lim_p'])
    ymax = np.nanmax(nav)
    ymin = -0.05 * ymax  # Pad slightly below zero
    ax.set_ylim(ymin, ymax * 1.1)

    ax.axvline(vmin, ls='--', color='b')
    ax.axvline(vmax, ls='--', color='r')
    ax.axvline(0, ls='--', color='k')
    ax.axhline(0, color='k')

    if ii == parent.nions[parent.page] - 1:
        ax.set_xlabel('Velocity (km/s)')
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])

    # Flag and Column Density Text
    flag = ion_data.get('nav_flag', 0)
    flag_text = {0: "OK", -1: "ND", -2: "SAT"}.get(flag, "NA")

    logN = ion_data.get('ncol_pyn', np.nan)
    logN_lo = ion_data.get('ncol_err_lo', np.nan)
    logN_hi = ion_data.get('ncol_err_hi', np.nan)

    if flag == -1 and 'ncol_linear2sig' in ion_data:
        text = r"$\log N < {0:.2f}$".format(ion_data['ncol_linear2sig'])
    elif np.isfinite(logN) and np.isfinite(logN_lo) and np.isfinite(logN_hi):
        text = r"$\log N = {0:.2f}^{{+{1:.2f}}}_{{{2:.2f}}}$".format(logN, logN_hi, logN_lo)
    else:
        text = r"$\log N$: NA"

    ax.text(0.05, 0.85, text, transform=ax.transAxes, fontsize=9, color='indigo')
    ax.text(0.95, 0.95, f"Flag: {flag_text}", transform=ax.transAxes,
            fontsize=8, color='crimson', ha='right', va='top')

    parent.figs[parent.page].canvas.draw()

def shift2vel(z1,z2):
    # z1 at rest
    # z2 for which relative velocity is computed
    spl=2.9979e5;  #speed of light
    vel=((1215.67*(1.+z2))-(1215.67*(1.0 + z1)))*spl/(1215.67*(1.0 + z1))#(self.wrest-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))
    return vel

def grab_intervening_linelist(filename,z_gal,wrest_galaxy,wavelength):
    
    #s=ascii.read(filename)
    s=pd.read_csv(filename,sep=',',encoding='latin1')
    ion=s['Name'].values#s['col1']
    wobs=s['Wave_obs'].values #s['col3']
    zobs=s['Zabs'].values#s['col4']
    wrest=wobs/(1.+zobs)#s['col2']

    spl=2.9979e5;  #speed of light

    #compute relative velocity difference with the host galaxy
    delv=shift2vel(z_gal,zobs)
    # Now select only lines within the given wavelength range
    window_max=np.max(wavelength)#*(1.+z_gal)
    window_min=np.min(wavelength)#*(1.+z_gal)
    q=np.where( (wobs >= window_min) & (wobs <= window_max))#& ( np.abs(delv) >= 200.) )
    outlist={}
    if np.sum(q)>0:
        #If there are lines create a new list
        outlist['ion']=ion[q]
        outlist['wrest']=wrest[q]
        outlist['wobs']=wobs[q]
        outlist['zobs']=zobs[q]
        wobs_small=wobs[q]
        outlist['number']=len(wobs_small)
        vel=np.zeros((len(wobs_small),))
        outlist['delv']=delv[q]
        #Now computing velocity for each
        for i in range(0,len(wobs_small)):
            vel[i] = (outlist['wobs'][i]-wrest_galaxy*(1.0 + z_gal))*spl/(wrest_galaxy*(1.0 + z_gal))
        outlist['vel']=vel
    else:
        outlist['number']=0

    return outlist



def plot_intervening_lines(ax,outlist,delv):
    #Plot all intervening lines. 
    # Other lines associated with current absorber = blue
    # All other intervening lines = Red
    if outlist['number']>0:
        #print('There are Intervening Lines!')
        vellist=outlist['vel']
        relative_vel_z=outlist['delv']
        for index in range(0,outlist['number']):
            if (np.abs(vellist[index]) <delv):
                #print(vellist[index])

                if np.abs(relative_vel_z[index]) >200.:
                    color =clr['pale_red']
                else:
                    color='b'
                ax.text(vellist[index],1.05, str(outlist['ion'][index])+' '+ str('%.0f' % outlist['wrest'][index]),
                    fontsize=8,rotation=90, rotation_mode='anchor',color=color)
                ax.text(vellist[index]+50.,1.05, 'z = '+str('%.3f' % outlist['zobs'][index]),
                    fontsize=8,rotation=90, rotation_mode='anchor',color=color)


class mainWindow(QtWidgets.QTabWidget):
    
    def __init__(self,ions, parent=None,intervening=False,instrument = None):
        #-----full spectra properties---------#
        self.instrument = instrument or 'GENERIC'
        self.z = ions['Target']['z']; self.flux = ions['Target']['flux']
        self.wave = ions['Target']['wave']; self.error = ions['Target']['error']
        #------------------------------------#
        self.old_axes = None
        self.tab_names = [f'Ions {i+1}' for i in range(20)]  # supports up to 20 tabs

        self.ions = ions
        self.keys = list(self.ions.keys())[:-1] # last item is the full target spectrum
        # --- Ensure contamination_mask exists for all ions ---
        for key in self.keys:
            if 'contamination_mask' not in self.ions[key]:
                wave = self.ions[key]['wave']
                self.ions[key]['contamination_mask'] = np.zeros_like(wave, dtype=bool)

        # -- Auto-continuum fitting for all ions --
        for key in self.keys:
            ion_data = self.ions[key]
            wave = ion_data['wave']
            flux = ion_data['flux']
            error = ion_data['error']
            zabs = ion_data['z']
            fval = ion_data['f']
            wline = ion_data['lam_0']
            order = ion_data.get('order', 4)
            if 'cont' not in ion_data or not np.any(np.isfinite(ion_data['cont'])):
                # Do auto-continuum only if no valid continuum exists
                try:
                    xcut, ycut, yerrcut, cont_model, coeffs, final_mask, _ = estimate_continuum(
                        wave, flux, error, wline, zabs, vwin=1000, order=order
                    )
                    cont = np.full_like(flux, np.nan)
                    mask_full = np.zeros_like(flux, dtype=bool)
                    in_window = (wave >= xcut.min()) & (wave <= xcut.max())
                    cont[in_window] = cont_model
                    mask_full[in_window] = final_mask

                    self.ions[key]['cont'] = cont
                    self.ions[key]['order'] = order
                    self.ions[key]['pco'] = coeffs
                    self.ions[key]['wc'] = mask_full
                    # --- Normalize flux and set initial y-limits for plotting ---
                    norm_flux = flux / cont
                    ylims = [0.0, np.nanmax(norm_flux) * 1.2]
                    self.ions[key]['y_lim'] = ylims
                    self.ions[key]['y_lim_init'] = ylims.copy()

# --- Compute velocity array for this ion (needed for axis limits, masks) ---
                    vel = (wave - wline * (1 + zabs)) / wline / (1 + zabs) * 2.9979e5
                    self.ions[key]['vel'] = vel
                    self.ions[key]['contamination_mask'] = np.zeros_like(wave, dtype=bool)

# --- Set default window limits if not already present ---
                    if 'window_lim' not in self.ions[key]:
                        window_lim = [-1000, 1000]
                        self.ions[key]['window_lim'] = window_lim
                        self.ions[key]['window_lim_p'] = window_lim.copy()

                except Exception as e:
                    print(f"[AutoFit] Skipped {key} due to: {e}")

        self.wc = None #initialize overall parameter
        self.ylims = []
        self.vclim = None
        self.vclim_contam = None

        self.vclim_all = []
        self.name = None
        self.EWlim = [None,None] #left,right
        self.event_button = None
        self.Manual_Mask = None
        self.pFlag = 1
        self.intervening=intervening
        self.save = False
        super(mainWindow,self).__init__(parent)
        
#---------------Initial page setup------------------# 
        self.page = 0
        self.tabs = [QtWidgets.QWidget()]
        self.addTab(self.tabs[self.page],self.tab_names[self.page])
        self.figs = [Figure()]
        self.canvas = [FigureCanvasQTAgg(self.figs[self.page])]

        self.nions = [len(self.keys)]
        if self.nions[0] > 6: self.nions[0]=6
        self.page=0
        #Need to set click focus for keyboard functionality
        self.setParent(parent)
        self.figs[0].canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.figs[0].canvas.setFocus()

        #LAYOUT
        #Top layer is a horizontal layout with help|save|add_ion; middle is canvas; bottom is Page number
        self.add_ion_button = QPushButton("Add Ion",self)
        self.add_ion_button.setGeometry(630,30,200,30)
        self.add_ion_button.clicked.connect(lambda: self.NewTransition(self))
        
        self.remove_ion_button = QPushButton("Remove Ion",self)
        self.remove_ion_button.setGeometry(1030,30,200,30)
        self.remove_ion_button.clicked.connect(lambda: self.removeIon(self))

        self.openButton = QPushButton("Help",  self)
        self.openButton.setGeometry(830,30,200,30)
        self.openButton.clicked.connect(lambda: self.opensub(self))
        
        self.saveButton = QPushButton("Save",  self)
        self.saveButton.setGeometry(430,30,200,30)
        self.saveButton.clicked.connect(lambda: self.onsave(self))
        
        self.PageLabel = QtWidgets.QLabel("Page: " + str(self.currentIndex()+1)+"/"+str(len(self.figs)),self)
        self.PageLabel.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
        
        self.zoom_combo = QtWidgets.QComboBox(self)
        self.zoom_combo.setGeometry(340, 30, 180, 30)  # Position it to the left of Save
        self.zoom_combo.addItems(["Zoom ±1000 km/s", "Zoom ±500 km/s", "Zoom ±250 km/s"])
        self.zoom_combo.setCurrentIndex(0)
        self.zoom_combo.currentIndexChanged.connect(self.apply_zoom_from_combo)

        self.main_layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        self.bot_layout = QHBoxLayout()
        
        self.spacerItem = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        
        self.top_layout.addItem(self.spacerItem)
        self.top_layout.addWidget(self.add_ion_button)
        self.top_layout.addWidget(self.saveButton)
        self.top_layout.addWidget(self.openButton)
        self.top_layout.addWidget(self.remove_ion_button)
        self.top_layout.addItem(self.spacerItem)
        
        self.bot_layout.addItem(self.spacerItem)
        self.bot_layout.addWidget(self.PageLabel)
        self.status_label = QLabel(" ", self)
        self.status_label.setStyleSheet("color: white; font: 20pt;")
        self.bot_layout.addWidget(self.status_label)

        self.bot_layout.addItem(self.spacerItem)
        
        self.main_layout.addLayout(self.top_layout,stretch=1)
        self.main_layout.addWidget(self.canvas[0],stretch=14)
        self.main_layout.addLayout(self.bot_layout,stretch=1) #0.5
        self.tabs[self.page].setLayout(self.main_layout)
        
        #initializing left and right axes
        self.axesL = [list(range(6))]; self.axesR = [list(range(6))]; self.axesN = [list(range(6))]
        for ii in range(self.nions[0]):
            self.axesL[0][ii] = self.figs[0].add_subplot(6, 3, 3 * ii + 1)
            self.axesR[0][ii] = self.figs[0].add_subplot(6, 3, 3 * ii + 2)
            self.axesN[0][ii] = self.figs[0].add_subplot(6, 3, 3 * ii + 3)

            self.figs[self.page].subplots_adjust(hspace=0.01, left=0.05, right=0.98, wspace=0.1)
                
            Plotting(self,ii,modify=True)
            
         # Set up connectivity
        self.cid1 = self.figs[0].canvas.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.figs[0].canvas.mpl_connect("key_press_event", self.onpress)
        self.cid3 = self.figs[0].canvas.mpl_connect("motion_notify_event",self.onmotion)
        
#----------------Setup for Additional pages-----------#  
        AddPage(self)

                       
#---------------------Save Button/function----------------# 
        

#--------------------------------------------------------------#

#         self.PageLabel = QtWidgets.QLabel("Page: " + str(self.currentIndex()+1)+"/"+str(len(self.figs)),self)
#         self.PageLabel.setStyleSheet("font: 16pt;color: black;background-color:white")
#         self.PageLabel.setGeometry(630,850,200,30)
        
    def apply_zoom_from_combo(self):
        selected = self.zoom_combo.currentText()
        if "1000" in selected:
            vmax = 1000
        elif "500" in selected:
            vmax = 500
        elif "250" in selected:
            vmax = 250
        else:
            vmax = 1000  # default fallback

        for page_idx in range(len(self.figs)):
            for ii in range(self.nions[page_idx]):
                key_idx = page_idx * 6 + ii
                if key_idx >= len(self.keys):
                    continue
                key = self.keys[key_idx]
                self.ions[key]['window_lim_p'] = [-vmax, vmax]

                try:
                    self.axesL[page_idx][ii].set_xlim(-vmax, vmax)
                    self.axesR[page_idx][ii].set_xlim(-vmax, vmax)
                    self.axesN[page_idx][ii].set_xlim(-vmax, vmax)
                except Exception:
                    pass

            self.figs[page_idx].canvas.draw_idle()


        
    def getPage(self):
        self.page = self.currentIndex()
        self.PageLabel.setText("Page: " + str(self.page+1)+"/"+str(len(self.figs)))
        self.currentChanged.connect(lambda: getPage(self))
            
#-------------------Add Ion Button------------------------------# 
    def NewTransition(self,parent):
    #will run back through absorber class to identify line, obtain slice window, 
        new_line,ok3 = QInputDialog.getDouble(self,'Add Line','Enter new transition:')#,decimals=4)
        if ok3:
            #add new line
            new_abs = Absorber_pn.Absorber(self.z,self.wave,self.flux,self.error,[new_line])
            #update initial dictionary
            self.ions.update(new_abs.ions); self.ions['Target'] = self.ions.pop('Target')# moves Target to last index
            for key in list(new_abs.ions.keys())[:-1]: self.keys.append(key)

#                 self.page = len(self.nions) - 1 #finds current page max
            if self.nions[self.page] < 6: #proceed with filling page
                ii = self.nions[self.page]
                self.axesL[self.page][ii] = self.figs[self.page].add_subplot(6,2,2*(ii)+1)
                self.axesR[self.page][ii] = self.figs[self.page].add_subplot(6,2,2*(ii+1))
                self.figs[self.page].subplots_adjust(hspace=0.01, left=0.05, right=0.98, wspace=0.1)
                self.nions[self.page] = self.nions[self.page]+1
                Plotting(self,ii,modify=True)

            else: # need to add a new page
                AddPage(self)
                    
#         self.add_ion_button = QPushButton("Add Ion",self)
#         self.add_ion_button.setGeometry(630,30,200,30)
#         self.add_ion_button.clicked.connect(lambda: NewTransition(self))
        
#--------------------------------------------------------------------# 

    def handle_nav_keypress(self):
        if self.Ridx is not None:
            key_idx = self.page * 6 + self.Ridx
            ion_key = self.keys[key_idx]
            ion_data = self.ions[ion_key]

        # --- Compute Na(v) ---
            nav_spec = compute_nav_profile(ion_data)
            self.ions[ion_key]['flux'] = nav_spec['flux']
            self.ions[ion_key]['error'] = nav_spec['eflux']
            self.ions[ion_key]['cont'] = nav_spec['contin']  # optional: update continuum too
            self.ions[ion_key]['tau_array'] = nav_spec.get('tau_array')
            self.ions[ion_key]['tau_array_err'] = nav_spec.get('tau_array_err')
            self.ions[ion_key]['fnorm_err'] = nav_spec.get('fnorm_err')
            self.ions[ion_key]['contin_err'] = nav_spec.get('contin_err')

            self.ions[ion_key]['integration_weights'] = nav_spec['integration_weights']
            self.ions[ion_key]['Nav'] = nav_spec['Nav']
            self.ions[ion_key]['Nav_err'] = nav_spec['Nav_err']
            self.ions[ion_key]['Nav_vel'] = nav_spec['vel']
            self.ions[ion_key]['EW'] = nav_spec['EW']
            self.ions[ion_key]['EW_err'] = nav_spec['EW_err']
            self.ions[ion_key]['N'] = nav_spec['ncol']
            self.ions[ion_key]['Nsig'] = nav_spec['ncol_err_hi']  # or average of hi/lo

            self.ions[ion_key]['ncol_pyn'] = nav_spec['ncol']
            self.ions[ion_key]['ncol_err_lo'] = nav_spec['ncol_err_lo']
            self.ions[ion_key]['ncol_err_hi'] = nav_spec['ncol_err_hi']
            self.ions[ion_key]['va'] = nav_spec['va']
            self.ions[ion_key]['va_err'] = nav_spec['va_err']
            self.ions[ion_key]['nav_flag'] = nav_spec['pyn_flag']
            self.ions[ion_key]['ba'] = nav_spec.get('ba', np.nan)
            self.ions[ion_key]['ba_err'] = nav_spec.get('ba_err', np.nan)
            self.ions[ion_key]['dv90'] = nav_spec.get('dv90', np.nan)
            self.ions[ion_key]['dv90_err'] = nav_spec.get('dv90_err', np.nan)
            self.ions[ion_key]['ncol_linear2sig'] = nav_spec.get('ncol_linear2sig', np.nan)
        # --- Plot Na(v) ---
            plot_nav_panel(self, key_idx, self.Ridx)

        # --- Also compute EW/AOD and populate those ---
            if USE_COMPUTE_EW:
                try:
                    EW(self, self.page, self.Ridx, ion_data['EWlims'])
                except Exception as e:
                    print(f"[Na(v)+EW] Failed EW computation for {ion_key}: {e}")
            else:
                Plotting(self, self.Ridx, modify=False, Print=True)


#-------HELP------#
    def onsave(self,parent):
            self.savepg = SavePage(self)
            if self.savepg.closewin == None:
                return None
            else:
                self.savepg.show()
                
    def opensub(self,parent):
        self.sub = HelpWindow()
        self.sub.show()

# ----------------- Updated by saloni to remove ions from the GUI ------------------------ #
    def removeIon(self, parent):
        ion_to_remove, ok4 = QInputDialog.getText(self, 'Remove Ion', 'Enter ion to be removed:')
        if ok4:
            if ion_to_remove in self.ions:
                index_removed = self.keys.index(ion_to_remove)

                # Remove from data structures
                self.ions.pop(ion_to_remove)
                self.keys.remove(ion_to_remove)

                total_ions = len(self.keys)
                total_pages = (total_ions + 5) // 6  # ceil division

                # If no ions remain
                if total_ions == 0:
                    self.page = 0
                    self.Ridx = None
                    self.Lidx = None
                else:
                    # Adjust current page to the last valid one if needed
                    self.page = min(self.page, total_pages - 1)
                    new_index = min(index_removed, total_ions - 1)
                    self.Ridx = new_index % 6
                    self.Lidx = self.Ridx

                # Trim unused pages if needed
                if len(self.figs) > total_pages:
                    for pg in range(total_pages, len(self.figs)):
                        self.removeTab(pg)
                    self.figs = self.figs[:total_pages]
                    self.canvas = self.canvas[:total_pages]
                    self.tabs = self.tabs[:total_pages]
                    self.axesL = self.axesL[:total_pages]
                    self.axesR = self.axesR[:total_pages]
                    self.axesN = self.axesN[:total_pages]
                    self.nions = self.nions[:total_pages]

                self.updatePlotAfterIonRemoval()
            else:
                QMessageBox.warning(self, 'Warning', 'Ion not found.')

    def updatePlotAfterIonRemoval(self):
        """Redraw plots and update state after ion removal."""
        if self.page >= len(self.figs):  # avoid out-of-bounds if page count changed
            self.page = len(self.figs) - 1

        # Clear figure and axes on current page
        self.figs[self.page].clf()
        self.axesL[self.page] = []
        self.axesR[self.page] = []
        self.axesN[self.page] = []

        if len(self.keys) == 0:
            self.nions[self.page] = 0
            self.canvas[self.page].draw()
            self.PageLabel.setText(f"Page: {self.page + 1}/{len(self.figs)}")
            return

        start_idx = self.page * 6
        end_idx = min(start_idx + 6, len(self.keys))
        self.nions[self.page] = end_idx - start_idx

        for ii in range(self.nions[self.page]):
            actual_ion_idx = start_idx + ii
            self.axesL[self.page].append(self.figs[self.page].add_subplot(6, 3, 3 * ii + 1))
            self.axesR[self.page].append(self.figs[self.page].add_subplot(6, 3, 3 * ii + 2))
            self.axesN[self.page].append(self.figs[self.page].add_subplot(6, 3, 3 * ii + 3))

            try:
                Plotting(self, ii, modify=True)
            except Exception as e:
                print(f"Plotting failed for index {ii}: {e}")

        self.figs[self.page].subplots_adjust(hspace=0.01, left=0.05, right=0.98, wspace=0.1)
        self.canvas[self.page].draw()
        self.PageLabel.setText(f"Page: {self.page + 1}/{len(self.figs)}")


    def onmotion(self, event):
        if event.inaxes is None:
            self.status_label.setText(" ")
            return
    
    # Find which page this event came from by checking the figure canvas
        event_page = None
        for page_idx, fig in enumerate(self.figs):
            if event.canvas == fig.canvas:
                event_page = page_idx
                break
    
        if event_page is None:
            self.status_label.setText(" ")
            return
    
    # Update the active page to match where the mouse is
        old_page = self.page
        self.page = event_page
    
    # Find which subplot contains the current axes (only check the event page)
        found = False
        if event.inaxes in self.axesL[event_page]:
            self.Lidx = self.axesL[event_page].index(event.inaxes)
            self.Ridx = None
            found = True
        elif event.inaxes in self.axesR[event_page]:
            self.Ridx = self.axesR[event_page].index(event.inaxes)
            self.Lidx = None
            found = True
        elif event.inaxes in self.axesN[event_page]:
            self.Lidx = None
            self.Ridx = None
            found = True
    
        if not found:
            self.status_label.setText(" ")
            return
    
    # Highlight active axes - clear old highlighting first
        if hasattr(self, 'old_axes') and self.old_axes and self.old_axes != event.inaxes:
            for pos in ['top', 'bottom', 'left', 'right']:
                self.old_axes.spines[pos].set_edgecolor('black')
                self.old_axes.spines[pos].set_linewidth(0.5)
        # Redraw the figure that contained the old axes
            if hasattr(self, 'old_page') and self.old_page is not None:
                self.figs[self.old_page].canvas.draw()
    
    # Highlight new axes
        if not hasattr(self, 'old_axes') or self.old_axes != event.inaxes:
            for pos in ['top', 'bottom', 'left', 'right']:
                event.inaxes.spines[pos].set_edgecolor('#01DF01')
                event.inaxes.spines[pos].set_linewidth(2)
            self.figs[self.page].canvas.draw()
            self.old_axes = event.inaxes
            self.old_page = self.page  # Store which page had the old axes
    
            self.status_label.setText(" ")

#----------------------key button events-----------------------------#            
    
    '''on press is used to reduce the order for the polyfit, toggle measurement displays, measure properties, and assist selecting EW vel bounds'''
    def onpress(self,event):
        if event.key == 'v':
            if self.old_axes in self.axesL[self.page]:
                mask_reg,ok = QInputDialog.getText(self,'Manual Mask Limits','Input Region to Mask (e.g. 200,250)')

                if ok:
                    key_idx = (self.page*6)+self.Lidx
                    vel = self.ions[self.keys[key_idx]]['vel']
                    wc = self.ions[self.keys[key_idx]]['wc']

                    mask = mask_reg.split(',')
                    mask = np.array(mask).astype('float32')

                    wc=((vel<mask[0]) | (vel>mask[1])) & wc
                    self.ions[self.keys[key_idx]]['wc'] = wc
                    Plotting(self,self.Lidx,modify=True)
                    
            elif self.old_axes in self.axesR[self.page]:
                integ_lims,ok = QInputDialog.getText(self,'Manual EW Limits','Input integration region (eg -100,100)')
                if ok:
                    key_idx = (self.page*6)+self.Ridx
                    integ_lims = integ_lims.split(',')
                    integ_lims = np.array(integ_lims).astype('float32')

                    self.ions[self.keys[key_idx]]['EWlims'][0] = integ_lims[0]
                    self.ions[self.keys[key_idx]]['EWlims'][1] = integ_lims[1]
                    Plotting(self,self.Ridx,modify=False,Print=False)

        
        #right arrow key (directional) increases poly order 
        if event.key == 'up':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx
                
                self.ions[self.keys[key_idx]]['order'] = self.ions[self.keys[key_idx]]['order']+1
                Plotting(self,self.Lidx,modify=True)
            else:
                print('click on a left transition window first')
                
        # Reduce polynomial, but don't allow order < 0
        if event.key == 'down':
            if self.Lidx is not None:
                key_idx = self.page * 6 + self.Lidx
                current_order = self.ions[self.keys[key_idx]]['order']
                self.ions[self.keys[key_idx]]['order'] = max(0, current_order - 1)
                Plotting(self, self.Lidx, modify=True)
            else:
                print('click on a transition window first')

        # Added by saloni
        if event.key == ']':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_xlim = self.ions[self.keys[key_idx]]['window_lim_p']
                new_xlim = [current_xlim[0] + 25, current_xlim[1] + 25]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['window_lim_p']=new_xlim
                Plotting(self,self.Lidx,modify=True)   

        if event.key == '[':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_xlim = self.ions[self.keys[key_idx]]['window_lim_p']
                new_xlim = [current_xlim[0] - 25, current_xlim[1] - 25]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['window_lim_p']=new_xlim
                Plotting(self,self.Lidx,modify=True)   

        if event.key == 'x':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_xlim = self.ions[self.keys[key_idx]]['window_lim_p']
                new_xlim = [current_xlim[0] * 0.9, current_xlim[1] * 0.9]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['window_lim_p']=new_xlim
                Plotting(self,self.Lidx,modify=True)  

        if event.key == 'X':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_xlim = self.ions[self.keys[key_idx]]['window_lim_p']
                new_xlim = [current_xlim[0] * 1.1, current_xlim[1] * 1.1]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['window_lim_p']=new_xlim
                Plotting(self,self.Lidx,modify=True)  

        if event.key == 'y':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] * 0.9, current_ylim[1] * 0.9]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True)         

        if event.key == 'Y':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] * 1.1, current_ylim[1] * 1.1]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True) 

        if event.key == 'w':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] + 0.2, current_ylim[1] + 0.2]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True) 

        if event.key == 's':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] -0.2, current_ylim[1] -0.2]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True) 

        if event.key == 'W':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] + 10, current_ylim[1] + 10]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True) 

        if event.key == 'S':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx   

                current_ylim = self.ions[self.keys[key_idx]]['y_lim']
                new_ylim = [current_ylim[0] -10, current_ylim[1] -10]  # Adjust zoom factor as needed
                self.ions[self.keys[key_idx]]['y_lim']=new_ylim
                Plotting(self,self.Lidx,modify=True) 

        if event.key == 'm':
            key_idx = self.page*6+self.Ridx
            EW(self,self.page,self.Ridx,self.ions[self.keys[key_idx]]['EWlims'])

                
        #evaluate all ions with limits on page
        if event.key == 'M':
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):
                    EW(self,self.page,ii,self.ions[self.keys[ii+self.page*6]]['EWlims'])

            
        
        #use same range for all ion Velocity limits and directly measured by following with clicks bounds on a single subplot
        if event.key == 'V':
            EWlims = self.ions[self.keys[self.Ridx + 6 * self.page]]['EWlims']
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):
                    self.ions[self.keys[ii + self.page * 6]]['EWlims'] = EWlims.copy()
                    Plotting(self, ii, modify=False, Print=False)


        if event.key in ['0','1','2']: #detection, upperlimit, lower limit
            key_idx = self.page*6 +self.Ridx
            self.ions[self.keys[key_idx]]['flag'] = int(event.key)
            Plotting(self,self.Ridx,modify=False,Print=True)
            
        if event.key == 't':#toggle the EW/N display
            self.pFlag = (self.pFlag+1)%3
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):
                    Plotting(self,ii,modify=False,Print=True)
            
        if event.key == 'N':
            self.handle_nav_keypress()

        if event.key == 'r':
            if self.Lidx is not None:
                key_idx = self.page * 6 + self.Lidx
                ion = self.ions[self.keys[key_idx]]

                ion['window_lim_p'] = ion['window_lim'].copy()

                if 'y_lim_init' in ion:
                    ion['y_lim'] = ion['y_lim_init'].copy()
                else:
                    ion['y_lim'] = [0, 1.5]  # or some fallback default

                Plotting(self, self.Lidx, modify=True)

            elif self.Ridx is not None:
                key_idx = self.page * 6 + self.Ridx
                ion = self.ions[self.keys[key_idx]]

                ion['window_lim_p'] = ion['window_lim'].copy()

                if 'y_lim_init' in ion:
                    ion['y_lim'] = ion['y_lim_init'].copy()

                Plotting(self, self.Ridx, modify=False)

                # Na(v) reset if available
                if 'nav_ylim_init' in ion:
                    self.axesN[self.page][self.Ridx].set_ylim(ion['nav_ylim_init'])
                    self.figs[self.page].canvas.draw()

        if event.key == 'Q':  # Quit/exit the GUI
            import matplotlib.pyplot as plt
            plt.close('all')


#------------------------------click button events----------------------------#        
    def onclick(self, event):
        if event.button in [1, 3]:
            # --- Left panel: continuum masking --- #
            if self.Lidx is not None:
                key_idx = self.page * 6 + self.Lidx
                vel = self.ions[self.keys[key_idx]]['vel']
                name = self.ions[self.keys[key_idx]]['name']

                if self.vclim is None:
                    self.vclim = [event.xdata]
                    self.axesL[self.page][self.Lidx].plot(event.xdata, event.ydata, 'ro', ms=5)
                    self.figs[self.page].canvas.draw()
                else:
                    if event.xdata is not None:
                        vclim = np.sort(np.append(self.vclim, event.xdata))
                        self.vclim = None
                        wc = self.ions[self.keys[key_idx]]['wc']
                        if event.button == 1:
                            wc = ((vel < vclim[0]) | (vel > vclim[1])) & wc
                        else:
                            wc = ((vel > vclim[0]) & (vel < vclim[1])) | wc
                        self.ions[self.keys[key_idx]]['wc'] = wc
                        Plotting(self, self.Lidx, modify=True)
                    else:
                        self.vclim = None

        # --- Right panel: set EW limits --- #
            if self.Ridx is not None:
                key_idx = self.page * 6 + self.Ridx
                if event.button == 1:
                    self.EWlim[0] = event.xdata
                    self.ions[self.keys[key_idx]]['EWlims'][0] = event.xdata
                    Plotting(self, self.Ridx, modify=False, Print=False)
                elif event.button == 3:
                    self.EWlim[1] = event.xdata
                    self.ions[self.keys[key_idx]]['EWlims'][1] = event.xdata
                    Plotting(self, self.Ridx, modify=False, Print=False)
                       
        # --- Bottom panel: contamination masking --- #
            if event.inaxes in self.axesN[self.page]:
                self.Nidx = self.axesN[self.page].index(event.inaxes)
            
            # Initialize contamination mask selection state
                if not hasattr(self, 'contam_vclim'):
                    self.contam_vclim = None

                if not hasattr(self, 'contam_vclim_marker'):
                    self.contam_vclim_marker = None

                key_idx = self.page * 6 + self.Nidx
                vel = self.ions[self.keys[key_idx]]['vel']

            # Start new selection or complete existing selection
                if self.contam_vclim is None:
                # First click - start selection
                    self.contam_vclim = [event.xdata]
                    self.contam_vclim_marker = self.axesN[self.page][self.Nidx].plot(event.xdata, event.ydata, 'ro', ms=5)[0]
                    self.figs[self.page].canvas.draw()
                else:
                # Second click - complete selection and apply mask
                    if event.xdata is not None:
                        vclim = np.sort(np.append(self.contam_vclim, event.xdata))
                        self.contam_vclim = None
                    
                    # Remove the marker
                        if self.contam_vclim_marker:
                            self.contam_vclim_marker.remove()
                            self.contam_vclim_marker = None
                    
                    # Get current contamination mask
                        contam_mask = self.ions[self.keys[key_idx]]['contamination_mask']
                    
                    # Apply mask logic (same as left panel but for contamination_mask)
                        if event.button == 1:
                        # Left click: mark selected region as contaminated (True)
                            new_mask = (vel >= vclim[0]) & (vel <= vclim[1])
                            contam_mask[new_mask] = True
                        else:
                        # Right click: mark selected region as clean (False)
                            new_mask = (vel >= vclim[0]) & (vel <= vclim[1])
                            contam_mask[new_mask] = False
                    
                    # Update the mask
                        self.ions[self.keys[key_idx]]['contamination_mask'] = contam_mask
                    
                        # Redraw
                        Plotting(self, self.Nidx, modify=True)
                        # Force canvas refresh to ensure red shading appears immediately
                        self.figs[self.page].canvas.draw()
                        self.figs[self.page].canvas.flush_events()
                    else:
                    # Cancel selection if click is outside axes
                        self.contam_vclim = None
                        if self.contam_vclim_marker:
                            self.contam_vclim_marker.remove()
                            self.contam_vclim_marker = None


def get_mask_regions(vel, mask):
    from itertools import groupby
    from operator import itemgetter
    idxs = np.where(mask)[0]
    if len(idxs) == 0:
        return []
    grouped = []
    for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
        group = list(map(itemgetter(1), g))
        grouped.append((vel[group[0]], vel[group[-1]]))
    return grouped

class Plotting:
    def __init__(self, parent, ii, modify=False, Print=False, **kwargs):
        key_idx = ii + 6 * parent.page

        # ----- Remove deleted ions if needed -----
        for key, value in kwargs.items():
            if value in parent.keys:
                parent.keys.remove(value)

        # ----- Define variables for clarity -----
        vel = parent.ions[parent.keys[key_idx]]['vel']
        wave = parent.ions[parent.keys[key_idx]]['wave']
        error = parent.ions[parent.keys[key_idx]]['error']
        flux = parent.ions[parent.keys[key_idx]]['flux']
        weight = parent.ions[parent.keys[key_idx]]['weight']
        name = parent.ions[parent.keys[key_idx]]['name']
        wc = np.array(parent.ions[parent.keys[key_idx]]['wc'] & (error != 0))
        cont = parent.ions[parent.keys[key_idx]]['cont']
        window_lim = parent.ions[parent.keys[key_idx]]['window_lim']
        window_lim_p = parent.ions[parent.keys[key_idx]]['window_lim_p']
        order = parent.ions[parent.keys[key_idx]]['order']
        EWlims = parent.ions[parent.keys[key_idx]]['EWlims']
        lam_0 = parent.ions[parent.keys[key_idx]]['lam_0']
        fvals = parent.ions[parent.keys[key_idx]]['f']
        ylims = parent.ions[parent.keys[key_idx]]['y_lim']

        if not Print:
            if modify:
                try:
                    if len(wave[wc]) > 0 and order >= 0:
                        parent.ions[parent.keys[key_idx]]['pco'] = L.Legendre.fit(wave[wc], flux[wc], order, w=weight[wc])
                        parent.ions[parent.keys[key_idx]]['cont'] = parent.ions[parent.keys[key_idx]]['pco'](wave)
                    else:
                        raise ValueError("Invalid data for fitting")
                except Exception as e:
                    print(f"[PolyFit] Fallback due to error: {e}")
                    parent.ions[parent.keys[key_idx]]['pco'] = L.Legendre([1.0])  # flat continuum
                    parent.ions[parent.keys[key_idx]]['cont'] = np.ones_like(wave)

                cont = parent.ions[parent.keys[key_idx]]['cont']

                if wc[0]:
                    gray_idx = np.where(np.diff(wc, prepend=np.nan, append=np.nan))[0][1:]
                else:
                    gray_idx = np.where(np.diff(wc, prepend=np.nan, append=np.nan))[0]

                parent.axesL[parent.page][ii].clear()
                parent.axesL[parent.page][ii].step(vel[(flux > 0) & (error > 0)], flux[(flux > 0) & (error > 0)], color='k', where='mid')
                parent.axesL[parent.page][ii].step(vel[(flux > 0) & (error > 0)], error[(flux > 0) & (error > 0)], color='r', where='mid')
                parent.axesL[parent.page][ii].step(vel[(flux > 0) & (error > 0)], cont[(flux > 0) & (error > 0)], color='b', where='mid')
                parent.axesL[parent.page][ii].text(
                        0.9, 0.12, f"Order: {order}", 
                        transform=parent.axesL[parent.page][ii].transAxes,
                        fontsize=8, color='blue', ha='left', va='top',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
                                                )

                for zz in range(int(len(gray_idx) / 2)):
                    vel_gray = vel[gray_idx[zz * 2]:gray_idx[2 * zz + 1]]
                    flux_gray = flux[gray_idx[zz * 2]:gray_idx[2 * zz + 1]]
                    error_gray = error[gray_idx[zz * 2]:gray_idx[2 * zz + 1]]
                    parent.axesL[parent.page][ii].step(
                        vel_gray[(flux_gray > 0) & (error_gray > 0)],
                        flux_gray[(flux_gray > 0) & (error_gray > 0)],
                        where='mid',
                        color='lightgray',
                        linewidth=1.3,
                        alpha=1
                    )

            parent.axesR[parent.page][ii].clear()
            parent.axesR[parent.page][ii].step(
                vel[(flux > 0) & (error > 0)],
                flux[(flux > 0) & (error > 0)] / cont[(flux > 0) & (error > 0)],
                color='k',
                where='mid'
            )
            parent.axesR[parent.page][ii].step(
                vel[(flux > 0) & (error > 0)],
                error[(flux > 0) & (error > 0)] / cont[(flux > 0) & (error > 0)],
                color='r',
                where='mid'
            )
            parent.axesR[parent.page][ii].axhline(y=1, xmin=window_lim[0], xmax=window_lim[1], ls='--', c='b')

            parent.axesL[parent.page][ii].set_yticks([])
            parent.axesL[parent.page][ii].set_ylabel(name)
            parent.axesR[parent.page][ii].set_ylabel(name)

            parent.axesL[parent.page][ii].set_xlim(window_lim_p)
            parent.axesR[parent.page][ii].set_xlim(window_lim_p)
            parent.axesL[parent.page][ii].set_ylim(ylims)
            parent.axesR[parent.page][ii].set_ylim([0, 2.2])

            # Show tick marks on all axes, but label only the bottom one
            for ax in [parent.axesL[parent.page][ii], parent.axesR[parent.page][ii]]:
                ax.tick_params(axis='x', which='both', bottom=True, top=False)
                if ii == parent.nions[parent.page] - 1:
                    ax.set_xlabel('Velocity (km/s)')
                    ax.tick_params(labelbottom=True)
                else:
                    ax.tick_params(labelbottom=False)


            if ii == 0:
                parent.axesL[parent.page][0].set_title('Continuum Fitter')
                parent.axesR[parent.page][0].set_title('Normalized Spectra')
                parent.axesN[parent.page][0].set_title('Na(v) Profile')

            parent.axesL[parent.page][ii].axvline(0, ls='--', c='k')
            parent.axesR[parent.page][ii].axvline(0, ls='--', c='k')

            if EWlims[0] is not None:
                parent.axesR[parent.page][ii].axvline(EWlims[0], ymin=0, ymax=2.5, ls='--', c='b')
            if EWlims[1] is not None:
                parent.axesR[parent.page][ii].axvline(EWlims[1], ymin=0, ymax=2.5, ls='--', c='r')

            if parent.intervening is not False:
                outlist = grab_intervening_linelist(parent.intervening, np.double(parent.z), lam_0, wave)
                plot_intervening_lines(parent.axesR[parent.page][ii], outlist, np.max(vel))

            parent.axesR[parent.page][ii].text(
                0.8,
                0.85,
                r'$\log(f\lambda):\,$' + str('%.2f' % (np.log10(fvals * lam_0))),
                transform=parent.axesR[parent.page][ii].transAxes,
                color=clr['teal']
            )

            parent.figs[parent.page].canvas.draw()

            # ---- Auto-plot Na(v) if already available ----
            if all(k in parent.ions[parent.keys[key_idx]] for k in ['Nav', 'Nav_err', 'Nav_vel']):
                plot_nav_panel(parent, key_idx, ii)
            # --- Overlay contamination mask regions ---
            key = parent.keys[parent.page * 6 + ii]
            contam_mask = parent.ions[key].get('contamination_mask', None)

            if contam_mask is not None and np.any(contam_mask):
                vel = parent.ions[key]['vel']
                regions = get_mask_regions(vel, contam_mask)
                for vmin, vmax in regions:
                    parent.axesN[parent.page][ii].axvspan(vmin, vmax, color='red', alpha=0.3)

        if Print:
            if parent.ions[parent.keys[key_idx]]['EW_text'] is not None:
                try:
                    parent.ions[parent.keys[key_idx]]['EW_text'].remove()
                except (NotImplementedError, ValueError, AttributeError):
                    # Text artist might not be removable or already removed
                    pass

            plotText(parent, parent.ions[parent.keys[key_idx]])
            text = parent.ions[parent.keys[key_idx]]['text']
            EWtoggle = parent.axesR[parent.page][ii].text(.05, 0.85, text, transform=parent.axesR[parent.page][ii].transAxes)
            parent.ions[parent.keys[key_idx]]['EW_text'] = EWtoggle
            parent.figs[parent.page].canvas.draw()

#Calculates N,EW,V,          
class EW:
    def __init__(self,parent,page,ii,lims):

        #determine which page is being accessed
        key_idx = ii+6*parent.page
            
        #define variables for readability
        vel = parent.ions[parent.keys[key_idx]]['vel']
        wave = parent.ions[parent.keys[key_idx]]['wave']
        error = parent.ions[parent.keys[key_idx]]['error']
        flux = parent.ions[parent.keys[key_idx]]['flux']
        name = parent.ions[parent.keys[key_idx]]['name']
        zabs = parent.ions[parent.keys[key_idx]]['z']
        f0 = parent.ions[parent.keys[key_idx]]['f']
        lam_0 = parent.ions[parent.keys[key_idx]]['lam_0']
        cont = parent.ions[parent.keys[key_idx]]['cont']

        #compute EW/N/med_vel
        output = compute_EW(wave,flux/cont,lam_0,lims,error/cont,plot=False,zabs=zabs,f0=f0)
        #save variables in ion's respective dictionary
        parent.ions[parent.keys[key_idx]]['N'] = output['col']
        parent.ions[parent.keys[key_idx]]['Nsig'] = output['colerr']
        parent.ions[parent.keys[key_idx]]['EW'] = output['ew_tot']*1000
        parent.ions[parent.keys[key_idx]]['EWsig'] = output['err_ew_tot']*1000
        parent.ions[parent.keys[key_idx]]['med_vel'] = output['med_vel']
        parent.ions[parent.keys[key_idx]]['EWlims'] = lims
        
        #plot EW on page
        Plotting(parent,ii,modify=False,Print=True)

        
# plots measurements, with toggle will display EW/N 
class plotText:
    def __init__(self,parent,line):
        ew = line.get('EW')
        ew_err = line.get('EW_err')  #  pyn_batch

        # Handle None or missing values
        if ew is None or ew_err is None:
            EW_det_text = "N/A"
            EW_limit_text = "N/A"
            logN_det_text = "N/A"
        elif ew is None:
            EW_det_text = "N/A"
            EW_limit_text = "<{:.0f} m$\AA$".format(2. * ew_err)  # upper limit
            logN_det_text = "N/A"            
        else:
            EW_det_text = str('%.0f' % ew) + ' $\pm$ ' + str('%.0f' % ew_err) + ' m$\AA$'
            EW_limit_text = "<{:.0f} m$\AA$".format(2. * ew_err)  # upper limit
            logN_det_text= str('%.2f' % np.log10(line['N'])) +' $\pm$ ' + str('%.3f' % (np.log10(line['N']+line['Nsig']) - np.log10(line['N']))) + ' /cm$^2$'

        #line.flag is the line specific upper/lower/detections
        #pflag is the toggle button for which to show
        # Measurement: 
        if line['flag']==0:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_det_text
            elif parent.pFlag==2: text=logN_det_text

        # Upper Limit
        elif line['flag']==1:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_limit_text
            elif parent.pFlag==2: text="<{:.2f}".format(np.log10(line['Nsig']))

        elif line['flag']==2:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_det_text
            elif parent.pFlag==2: text=">{:.2f}".format(np.log10(line['N']))
        line['text'] = text
        
class AddPage:
    def __init__(self, parent):
        ions_per_page = 6
        total_ions = len(parent.keys)
        existing_pages = len(parent.figs)
        total_pages_needed = (total_ions + ions_per_page - 1) // ions_per_page

        # Only add more pages if needed
        while len(parent.figs) < total_pages_needed:
            parent.page = len(parent.figs)

            initialize(parent)

            canvas = parent.figs[parent.page].canvas
            canvas.mpl_connect("button_press_event", parent.onclick)
            canvas.mpl_connect("key_press_event", parent.onpress)
            canvas.mpl_connect("motion_notify_event", parent.onmotion)

            
class initialize:
    def __init__(self,parent):
        parent.tabs.append(QtWidgets.QWidget())
        parent.addTab(parent.tabs[parent.page], parent.tab_names[parent.page])
        parent.figs.append(Figure())
        parent.canvas.append(FigureCanvasQTAgg(parent.figs[parent.page]))
        layout = QtWidgets.QVBoxLayout()
#         layout.addWidget(parent.canvas[parent.page])
#         parent.tabs[parent.page].setLayout(layout)
        
        self.add_ion_button = QPushButton("Add Ion",parent)
        self.add_ion_button.setGeometry(630,30,200,30)
        self.add_ion_button.clicked.connect(lambda: parent.NewTransition(parent))
        
        self.openButton = QPushButton("Help",  parent)
        self.openButton.setGeometry(830,30,200,30)
        self.openButton.clicked.connect(lambda: parent.opensub(parent))
        
        self.saveButton = QPushButton("Save",  parent)
        self.saveButton.setGeometry(430,30,200,30)
        self.saveButton.clicked.connect(lambda: parent.onsave(parent))
        
        self.remove_ion_button = QPushButton("Remove Ion",  parent)
        self.remove_ion_button.setGeometry(1030,30,200,30)
        self.remove_ion_button.clicked.connect(lambda: parent.removeIon(parent))


        self.PageLabel = QtWidgets.QLabel("Page: " + str(parent.page+1)+"/"+str(len(parent.figs)),parent)
        self.PageLabel.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
        
        self.main_layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        self.bot_layout = QHBoxLayout()
        
        self.spacerItem = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.top_layout.addItem(self.spacerItem)
        self.top_layout.addWidget(self.add_ion_button)
        self.top_layout.addWidget(self.saveButton)
        self.top_layout.addWidget(self.openButton)
        self.top_layout.addWidget(self.remove_ion_button)
        self.top_layout.addItem(self.spacerItem)
        
        self.bot_layout.addItem(self.spacerItem)
        self.bot_layout.addWidget(self.PageLabel)
        self.bot_layout.addItem(self.spacerItem)
                                   
        self.main_layout.addLayout(self.top_layout,stretch=1)
        self.main_layout.addWidget(parent.canvas[parent.page],stretch=14)
        self.main_layout.addLayout(self.bot_layout,stretch=1)
        parent.tabs[parent.page].setLayout(self.main_layout)

        ions_per_page = 6
        start_idx = parent.page * ions_per_page
        end_idx = min(start_idx + ions_per_page, len(parent.keys))
        n_ions_on_page = end_idx - start_idx
        parent.nions.append(n_ions_on_page)

        #initializing left and right axes
        parent.axesL.append(list(range(6))); parent.axesR.append(list(range(6))); parent.axesN.append(list(range(6)))                   
        for ii in range(n_ions_on_page):
            parent.axesL[parent.page][ii] = parent.figs[parent.page].add_subplot(6,3,3*ii+1)
            parent.axesR[parent.page][ii] = parent.figs[parent.page].add_subplot(6,3,3*ii+2)
            parent.axesN[parent.page][ii] = parent.figs[parent.page].add_subplot(6,3,3*ii+3)

            parent.figs[parent.page].subplots_adjust(hspace=0.01, left=0.05, right=0.98, wspace=0.1)

            Plotting(parent,ii,modify=True)


        # Set up connectivity and Need to set click focus for keyboard functionality
        parent.figs[parent.page].canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        parent.figs[parent.page].canvas.setFocus()            
        parent.figs[parent.page].canvas.mpl_connect('motion_notify_event', parent.onmotion)



class HelpWindow(QtWidgets.QWidget):

    def __init__(self,parent=None):
        super(HelpWindow, self).__init__(parent)
        self.setWindowTitle("PyNorm GUI Help")
        self.resize(600, 700)
        
        # Create layout
        layout = QtWidgets.QVBoxLayout(self)
        
        # Use QTextBrowser for better text rendering with wrapping and scrolling
        text_browser = QtWidgets.QTextBrowser(self)
        text_browser.setMarkdown(False)  # Use plain text mode
        text_browser.setText(HELP)
        text_browser.setReadOnly(True)
        text_browser.setStyleSheet("QTextBrowser { font-family: monospace; font-size: 10pt; }")
        
        layout.addWidget(text_browser)
        self.setLayout(layout)
        
class SavePage(QtWidgets.QWidget):
    def __init__(self,parentvals,parent=None):
        super(SavePage, self).__init__(parent)
        self.resize(700,400)
        self.closewin = False

        def onpdf(self,parentvals):
            # need to eliminate highlighted axes
            for pos in ['top','bottom','left','right']:
                parentvals.old_axes.spines[pos].set_edgecolor('black')
                parentvals.old_axes.spines[pos].set_linewidth(0.5)
            figurefile = self.pdfline.text()
            i = 1
            for figures in parentvals.figs:
                figures.savefig(figurefile[:-4]+str(i)+'.pdf',bbox_inches='tight')
                i = i+1
            self.pdfsave.setStyleSheet('background-color : green')
            parentvals.save = True
        def ontable(self,parentvals,Table_e):
            savefile = self.tableline.text()
            ascii.write(Table_e,savefile,overwrite=True)
            #how to display to user it has been saved?
            self.tablesave.setStyleSheet("background-color : green")
            parentvals.save = True
            
        def onpickle(self,parentvals):
            pfile = self.pickleline.text()
            with open(pfile,'wb') as pklfile:
                pickle.dump(parentvals.ions,pklfile,protocol=pickle.HIGHEST_PROTOCOL)
            self.picklesave.setStyleSheet('background-color : green')
            parentvals.save = True
        # want to print out table to the user regardless of whether it is saved
        # Build measurement table
        Table_e = Table()
        Table_e['Transitions'] = parentvals.keys

        EW = []; EW_err = []; N = []; Nsig = []; Vel = []
        EWlims_low = []; EWlims_high = []

        N_pyn = []; N_lo_pyn = []; N_hi_pyn = []; pyn_flags = []

        va_list = []; va_err_list = []
        ba_list = []; ba_err_list = []
        dv90_list = []; dv90_err_list = []

        unevaluated_keys = ['EW', 'EW_err', 'N', 'Nsig', 'med_vel']
        unevaluated_ions = []

        # First pass: check if anything is missing
        for ion in parentvals.keys:
            this_ion = parentvals.ions[ion]
            if any(this_ion.get(k) is None for k in unevaluated_keys):
                unevaluated_ions.append(ion)

        # Warn if unevaluated
        if unevaluated_ions:
            reply = QMessageBox.question(
                self,
                'Warning',
                f"{len(unevaluated_ions)} ions are unevaluated. Proceed and fill with NaNs?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                self.closewin = None
                return self.closewin

        # Second pass: fill everything properly
        for ion in parentvals.keys:
            this_ion = parentvals.ions[ion]

        # Fill missing keys with NaN if needed
            for k in unevaluated_keys:
                if this_ion.get(k) is None:
                    this_ion[k] = np.nan

            EW.append(np.round(this_ion.get('EW', np.nan), 2))
            EW_err.append(np.round(this_ion.get('EW_err', np.nan), 2))
    
            N_val = this_ion.get('N', np.nan)
            Nsig_val = this_ion.get('Nsig', np.nan)
            N.append(np.round(np.log10(N_val), 2) if N_val > 0 else np.nan)
            Nsig.append(np.round(np.log10(Nsig_val), 2) if Nsig_val > 0 else np.nan)
            Vel.append(np.round(this_ion.get('med_vel', np.nan), 2))

            lims = this_ion.get('EWlims', [np.nan, np.nan])
            EWlims_low.append(np.round(lims[0], 2))
            EWlims_high.append(np.round(lims[1], 2))

            # pynorm values
            flag = this_ion.get('nav_flag', np.nan)
            pyn_flags.append(flag)

            if flag == -1:  # non-detection
                N_pyn.append(np.round(this_ion.get('ncol_linear2sig', np.nan), 2))
                N_lo_pyn.append(np.nan)
                N_hi_pyn.append(np.nan)
            elif 'ncol_pyn' in this_ion:
                N_pyn.append(np.round(this_ion['ncol_pyn'], 2))
                N_lo_pyn.append(np.round(this_ion['ncol_err_lo'], 2))
                N_hi_pyn.append(np.round(this_ion['ncol_err_hi'], 2))
            else:
                N_pyn.append(np.nan)
                N_lo_pyn.append(np.nan)
                N_hi_pyn.append(np.nan)

            # velocity structure parameters
            va_list.append(np.round(this_ion.get('va', np.nan), 2))
            va_err_list.append(np.round(this_ion.get('va_err', np.nan), 2))
            ba_list.append(np.round(this_ion.get('ba', np.nan), 2))
            ba_err_list.append(np.round(this_ion.get('ba_err', np.nan), 2))
            dv90_list.append(np.round(this_ion.get('dv90', np.nan), 2))
            dv90_err_list.append(
                np.round(this_ion.get('va_err', np.nan) * np.sqrt(2), 2) if this_ion.get('va_err') else np.nan
            )

        # Assign columns to the table
        Table_e['EW'] = EW
        Table_e['EW_err'] = EW_err
        Table_e['Vmin'] = EWlims_low
        Table_e['Vmax'] = EWlims_high
        Table_e['logN'] = N
        Table_e['logNsig'] = Nsig
        Table_e['logN_pyn'] = N_pyn
        Table_e['logN_lo_pyn'] = N_lo_pyn
        Table_e['logN_hi_pyn'] = N_hi_pyn
        Table_e['pyn_flag'] = pyn_flags
        Table_e['Vel'] = Vel
        Table_e['va'] = va_list
        Table_e['va_err'] = va_err_list
        Table_e['ba'] = ba_list
        Table_e['ba_err'] = ba_err_list
        Table_e['dv90'] = dv90_list
        Table_e['dv90_err'] = dv90_err_list


        print(Table_e)
        
        #pdf save
        pdflabel = QLabel("Enter path and filename: (e.g. pathname\Ions.pdf)",self)
        pdflabel.setGeometry(100,100,400,30)
        
        self.pdfline = QLineEdit(self)
        self.pdfline.setText(parentvals.instrument+'_z_'+str(parentvals.z)+'_Ions.pdf')
        self.pdfline.setGeometry(100,125,300,30)
        
        self.pdfsave = QPushButton("save PDF",self)
        self.pdfsave.setGeometry(410,125,200,30)
        self.pdfsave.clicked.connect(lambda: onpdf(self,parentvals))
        #table save
        tablelabel = QLabel("Enter path and filename: (e.g. pathname\Table.dat)",self)
        tablelabel.setGeometry(100,175,400,30)
        
        self.tableline = QLineEdit(self)
        self.tableline.setText(parentvals.instrument+"_z_"+str(parentvals.z)+"_Measurement_Table.dat")
        self.tableline.setGeometry(100,200,300,30)
        
        self.tablesave = QPushButton("Save Table",self)
        self.tablesave.setGeometry(410,200,200,30)
        self.tablesave.clicked.connect(lambda: ontable(self,parentvals,Table_e))
        
        #pickle save
        picklelabel = QLabel("Enter path and filename: (e.g. pathname\Table.p)",self)
        picklelabel.setGeometry(100,250,400,30)
        
        self.pickleline = QLineEdit(self)
        self.pickleline.setText(parentvals.instrument+"_z_"+str(parentvals.z)+".p")
        self.pickleline.setGeometry(100,275,300,30)
        
        self.picklesave = QPushButton("Save Progress",self)
        self.picklesave.setGeometry(410,275,200,30)
        self.picklesave.clicked.connect(lambda: onpickle(self,parentvals))

        
#Initial inputs and callable class to run proram        
class Transitions:
    def __init__(self,Abs,intervening=False,instrument=None):
        if not QtWidgets.QApplication.instance():
            # Set the exception hook BEFORE launching the app
            sys.excepthook = handle_exception
            app = QtWidgets.QApplication(sys.argv)
            app.setStyle("Fusion")

            # Now use a palette to switch to dark colors:
            palette = QPalette()
            palette.setColor(QPalette.Window, QColor(53, 53, 53))
            palette.setColor(QPalette.WindowText, QtCore.Qt.white)        
            palette.setColor(QPalette.Base, QColor(25, 25, 25))
            palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
            palette.setColor(QPalette.Button, QColor(53, 53, 53))
            palette.setColor(QPalette.ButtonText, QtCore.Qt.white)
            palette.setColor(QPalette.BrightText, QtCore.Qt.red)
            palette.setColor(QPalette.Link, QColor(42, 130, 218))
            palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
            palette.setColor(QPalette.Text, QtCore.Qt.white)
    
            app.setPalette(palette)

        else:
            app = QtWidgets.QApplication.instance() 




        #app = QtWidgets.QApplication(sys.argv)
        # Force the style to be the same on all OSs:
        #app.setStyle("Fusion")

        # Now use a palette to switch to dark colors:
        #palette = QPalette()
        #palette.setColor(QPalette.Window, QColor(53, 53, 53))
        #palette.setColor(QPalette.WindowText, QtCore.Qt.white)        
        #palette.setColor(QPalette.Base, QColor(25, 25, 25))
        #palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        #palette.setColor(QPalette.Button, QColor(53, 53, 53))
        #palette.setColor(QPalette.ButtonText, QtCore.Qt.white)
        #palette.setColor(QPalette.BrightText, QtCore.Qt.red)
        #palette.setColor(QPalette.Link, QColor(42, 130, 218))
        #palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        #palette.setColor(QPalette.Text, QtCore.Qt.white)

        #app.setPalette(palette)
        main = mainWindow(Abs,intervening=intervening,instrument=instrument)
        main.resize(1800, 900)
        main.show()
        QtWidgets.QApplication.setQuitOnLastWindowClosed(True)

        try:
            exit_code = app.exec_()
            main.deleteLater()         # Schedule proper deletion
            app.processEvents()        # Handle pending events
            sys.exit(exit_code)        # Clean exit
        except Exception as e:
            print(f"Error during shutdown: {e}")
            sys.exit(1)
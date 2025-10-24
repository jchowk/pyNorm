


'''
Absorber 

Inputs:
flux; wave; error; linelist; redshift; bin

1st. Asborber will bin the flux,wave and error to clean the data
2nd. will pull the actual lamd_rest from the atom.dat file for all lines
-->this (2nd) will return a dictionary of lamd_rest,ion_name,fval,gamma

3rd. Initialized entries for the Vstack plotting
'''
import numpy as np
import numpy.polynomial.legendre as L
import numpy as np
from astropy.io import ascii
from pkg_resources import resource_filename

def rb_setline(lambda_rest,method,linelist='atom',target_name=None):
    """
    Function to read in atomic line information for a given rest frame  wavelength.
                           Or 
    For the line matching the closest wavelength. 

    Parameters
    ----------
    lambda_rest :-  Rest Frame wavelength (in \AA) of the line to match
    method     :-   'closest' ->  If set will match the closest line.
                    'Exact'  ->  If set will match the exact wavelength.
                    'Name'   -> Match by name, USE WITH CARE. MUST INPUT OPTIONAL NAMELIST
 
    Returns
    ----------
    
    dic :- Dictionary with fval,lambda and species name.

    Example
    -------

       str=rb_setline(2796.3,'closest')


    Written By: Rongmon Bordoloi                Jan 2018, Python 2.7
    Edit:       Rongmon Bordoloi                            Sep 2018, Depreciated kwargs to be compatible with python 3
   """
    
    #if kwargs.has_key('linelist'):
    #   linelist=kwargs['linelist']
    #else:
    #   linelist='LLS'
    
    line_str=read_line_list(linelist)
    wavelist=np.zeros((len(line_str),))
    name = np.empty(len(line_str), dtype='object')
    fval=np.zeros((len(line_str),))
    if linelist=='atom':
        gamma=np.zeros((len(line_str),))


    for i in range(0,len(wavelist)):
        wavelist[i]=np.double(line_str[i]['wrest'])
        fval[i]=float(line_str[i]['fval'])
        name[i]=np.str_(line_str[i]['ion'])
        if linelist=='atom':
            gamma[i]=np.str_(line_str[i]['gamma'])

    if method=='Exact':
        q= np.where( (np.abs(lambda_rest-wavelist) < 1e-3))
        if linelist=='atom':
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q],'gamma':gamma[q]}
        else:
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q]}

    if method=='Name':
        #USE INPUT NAME LIST TO MATCH

        q= np.where(name == target_name)
        if linelist=='atom':
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q],'gamma':gamma[q]}
        else:
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q]}


    elif method=='closest':
        idx=(np.abs(lambda_rest-wavelist)).argmin()
        if linelist=='atom':
            outstr={'wave':wavelist[idx],'fval':fval[idx],'name':name[idx],'gamma':gamma[idx]}  

        else:

            outstr={'wave':wavelist[idx],'fval':fval[idx],'name':name[idx]} 
    else:
        raise NameError('Specify the matching method, closest or Exact')

    return outstr



def read_line_list(label):
    """Module to read a linelist defined by the label

    Parameters
    ----------
    lable : Label string [e.g. atom, LLS, LLS Small, LBG, Gal, Eiger_Strong]
      Must include redshift

    Returns
    ----------
    a dictionary with wrest, ion name and fvalues

    """
    

    if label=='atom':
        filename=resource_filename('pynorm.gui','lines/atom_full.dat')


    data = []

    if label=='atom':

        s=ascii.read(filename)

        for line in range(0,len(s['col1'])):
            source = {}
            source['wrest'] = float(s['col2'][line])
            ion_label = s['col1'][line]
            wrest_val = float(s['col2'][line])

            if (ion_label.startswith("NI")) and (ion_label!='NII'):  # Apply only to nitrogen NI species
                source['ion'] = f"{ion_label} {wrest_val:.2f}"
            else:
                source['ion'] = f"{ion_label} {int(wrest_val)}"

            source['fval']=float(s['col3'][line])
            source['gamma']=float(s['col4'][line])

            data.append(source)

    elif ((label =='LBG') | (label =='Gal')):

        s=ascii.read(filename)

        for line in range(0,len(s['wrest'])):
            source = {}
            source['wrest'] = float(s['wrest'][line])
            source['ion'] = s['name'][line]+' '+s['transition'][line]
            source['fval']=float(s['ID'][line])
            source['gamma']=float(s['ID'][line])

            data.append(source)

    elif (label =='Eiger_Strong') |(label =='Gal_Em') | (label =='Gal_Abs') |(label =='Gal_long') | (label =='AGN'):

        s=ascii.read(filename)

        for line in range(0,len(s['wrest'])):
            source = {}
            source['wrest'] = float(s['wrest'][line])
            source['ion'] = s['name'][line]#+' '+s['transition'][line]
            source['fval']=float(0)#s['ID'][line])
            source['gamma']=float(0)#s['ID'][line])

            data.append(source)

    elif (label =='HI_recomb') |((label =='HI_recomb_light')):
        s=ascii.read(filename)

        for line in range(0,len(s['wrest'])):
            source = {}
            source['wrest'] = float(s['wrest'][line]*10**4)
            source['ion'] = s['name'][line]#+' '+s['transition'][line]
            source['fval']=float(0)#s['ID'][line])
            source['gamma']=float(0)#s['ID'][line])

            data.append(source)

    else:       
        f=open(filename,'r')
        header1 = f.readline()
        for line in f:
            line = line.strip()
            columns = line.split()
            source = {}
            source['wrest'] = float(columns[0])
            source['ion'] = columns[1]+' '+columns[2]
            source['fval']=float(columns[3])
            data.append(source)


    return data

c =  299792.458

class Absorber:
    
    #defining variables to be used in the transition plotting GUI
    def Transition(self,ion_dict,line_dat,wave,flux,error,z,mask,window_lim,nofrills=False):
            # Edit RB May21, 2020: added nofrills keyword to toggle between initializing the continuum fields. 
            #        Added to avoid issues of calling Absorber class very near to the edge of the detector.
            #        Does not apply when calling abstools.

            # VARIABLE INITIALIZATION
            ion_dict['f'] = line_dat['fval']
            ion_dict['lam_0'] = line_dat['wave']
            ion_dict['name'] = line_dat['name']
            ion_dict['gamma'] = line_dat['gamma']
            ion_dict['z'] = z
            ion_dict['window_lim'] = window_lim
            ion_dict['window_lim_p'] = [-1000,1000]

            '''Shifting to vel-space centered on lam_o'''
            ion_dict['lam_0_z'] = ion_dict['lam_0']*(1+z)
            ion_dict['vel'] = (wave-ion_dict['lam_0_z'])/ion_dict['lam_0_z']*c

            '''create window for flux,wave,error based on max and min velocity'''
            window = (ion_dict['vel']>window_lim[0]) & (ion_dict['vel']<window_lim[1])
            ion_dict['flux'] = flux[window]; ion_dict['wave']=wave[window]
            ion_dict['error'] = error[window]; ion_dict['vel'] = ion_dict['vel'][window]
            min_flux = max(min(ion_dict['flux']),0.25*np.mean(ion_dict['flux'])); max_flux= min(max(ion_dict['flux']),2*np.mean(ion_dict['flux']))
            #ion_dict['y_lim'] = [min(ion_dict['flux']),max(ion_dict['flux'])]
            ion_dict['y_lim'] = [min_flux,max_flux]
            '''Initial Polyfit assuming a masked region of -200:200 and polyorder=4
            cont= continuum, pco= polynomial coefficients for cont fitting; weight= parameter to fit polys
            order = order of polynomial

            lets also give each ion object the Legendre function for ease of use during plotting'''

            if nofrills==False:
                ion_dict['wc'] = ((ion_dict['vel']<mask[0])|(ion_dict['vel']>mask[1]))&(ion_dict['vel']>-500)&(ion_dict['vel']<500)&(ion_dict['error'] != 0) #error != 0 is a bad pixel mask
                ion_dict['weight'] = 1/(ion_dict['error']**2)
                ion_dict['order'] = 4 #order of poly fit
                ion_dict['pco']=L.Legendre.fit(ion_dict['wave'][ion_dict['wc']],ion_dict['flux'][ion_dict['wc']],ion_dict['order'],w=ion_dict['weight'][ion_dict['wc']])
                ion_dict['cont'] = np.full_like(ion_dict['wave'][ion_dict['wc']], np.nan)


            '''Property initializations:'''
            ion_dict['N']=None; ion_dict['Nsig']=None; ion_dict['EW']=None; ion_dict['EWsig']=None
            ion_dict['med_vel'] = None; ion_dict['EWlims'] = [mask[0],mask[1]]; ion_dict['flag'] = 0
            #for text
            ion_dict['EW_text'] = None
    
    
    def __init__(self,z,wave,flux,error,lines=None,mask_init=[-200,200],window_lim=[-1000,1000],load_file = False,nofrills=False):
        mask = mask_init
        self.z =z
        self.ions = {}

        if lines:
            for line in lines:
                line_dat = rb_setline(line,method='closest')
                
                #if using Transition class uncomment below line. Also comment transition def, while uncommenting transition class, comment out lines 80-82
                #self.ions[line_dat['name']] =Transition(line_dat,wave,flux,error,self.z,mask_init,window_lim)
        
                self.ions[line_dat['name']] = {}
                ion_dict = self.ions[line_dat['name']]
                self.Transition(ion_dict,line_dat,wave,flux,error,z,mask,window_lim,nofrills=nofrills)
                                
            #last dictionary item for full spectra data                    
            self.ions['Target'] = {}
            self.ions['Target']['flux'] = flux
            self.ions['Target']['wave'] = wave
            self.ions['Target']['error'] = error
            self.ions['Target']['z'] = z
        else:
            print('Input Linelist and rerun')
            
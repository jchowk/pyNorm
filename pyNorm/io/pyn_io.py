def read_rbcodes(input_filename, targname, ra, dec, ion, partial_pixels=True):
    import numpy as np
    from collections import OrderedDict
    from scipy.io import readsav
    from pyNorm.aod import pyn_batch
    from pyNorm.continuum import continuum_fit
    import pickle
    from astropy.coordinates import SkyCoord

    # Read the saved pickle file:
    spec_in = []
    with (open(input_filename, "rb")) as openfile:
        while True:
            try:
                spec_in.append(pickle.load(openfile))
            except EOFError:
                break
    #adjust to the level of the specified ion in the pickle file
    spec_in = spec_in[0][ion]

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec_in['vel'].flags.writeable:
        spec_in = __fix_unwriteable_spec(spec_in)


    # Create the dictionary
    spec = OrderedDict()

    # Observation information
    try:
        spec['ion'] = spec_in['name'].decode('utf-8')
    except:
        spec['ion'] = spec_in['name']


    spec['wni'] = spec_in['name'].split( )[-1]  #should be iNorm wavelength label

    spec['wavc'] = spec_in['lam_0']
    spec['fval'] = spec_in['f']
    spec['gamma'] = spec_in['gamma']
    spec['redshift'] = spec_in['z']

    #this whole section through the vlsr is not adjusted for rbcodes yet
    #don't have the info in the pickle file right now
    # Sometimes targname key is missing:
    try:
        spec['targname'] = targname#spec_in['targname']
    except:
        spec['targname'] = 'NoTargName'#spec_in['object']

    # Fix binary format
    try:
        spec['object'] = spec_in['object'].decode('utf-8')
    except:
        spec['object'] = 'NoTargName'#spec_in['object']

    try:
        spec['targname'] = spec['targname'].decode('utf-8')
    except:
        spec['targname'] = spec['targname']

    # Coordinates
    spec['RA'] =ra
    spec['Dec'] =dec
    coords  = SkyCoord(ra, dec, unit="deg",frame = 'icrs')
    #fill in the lat and long
    spec['gl'] = coords.galactic.l.value
    spec['gb'] = coords.galactic.b.value

    # Sometimes LSR shift is missing
    try:
        spec['vlsr'] = lsrvel(spec['gl'],spec['gb'])
    except:
        spec['vlsr'] = float("nan")


    #
    # Raw data
    spec['vel'] = spec_in['vel']
    spec['flux'] = spec_in['flux']
    spec['eflux'] = spec_in['error']
    spec['wave'] = spec['wavc']*(spec['vel']/2.998e5)+spec['wavc']
    #
    # Continuum definition
    try:
        spec['contin'] = spec_in['cont']
        spec['contin_err'] = np.zeros(len(spec_in['cont']))
    except:
        spec['contin'] = spec_in['ycon']
        spec['contin_err'] = spec_in['ycon_sig']

    try:
        spec['contin_order'] = spec_in['order']
    except:
        spec['contin_order'] = 0

    spec['contin_coeff'] = np.array([0,0.])
    #
    # Construct continuum masks / velocity range
    #rb_codes saves the mask as a boolean array, need to convert to 0s and 1s before using the
    #__convert_inorm_mask function
    masked_area = (spec_in['wc'] == False) #these are the areas masked from rb_codes
    spec['contin_mask_bits'] = np.zeros(len(spec_in['wc'])) #True = included
    spec['contin_mask_bits'][masked_area] = 0
    spec['contin_mask_bits'][~masked_area] = 1
    spec['mask_cont'] = np.zeros(len(spec_in['wc']))
    spec['mask_cont'][masked_area] = 0
    spec['mask_cont'][~masked_area] = 1
    spec = __convert_inorm_mask(spec)
    #
    # Normalized spectrum
    spec['vnorm'] = spec_in['vel']
    spec['fnorm'] = spec_in['flux']/spec_in['cont']
    spec['fnorm_err'] = spec_in['error']/spec_in['cont']
    spec['fnorm_err_contin'] = 0.
    spec['fnorm_err_stat'] = 0.
    #
    # Na(v) data
    spec['Nav'] = float('nan')
    spec['Nav_err'] = float('nan')
    spec['Nav_sat'] = float('nan')
    #
    # Na(v) quantities
    try:
        spec['SNR'] = 0.0
    except:
        spec['SNR'] = 0.0

    spec['v1'] = spec_in['EWlims'][0]
    spec['v2'] = spec_in['EWlims'][1]
    spec['EW'] = spec_in['EW']
    spec['EW_err'] = spec_in['EWsig']
    spec['EW_err_stat'] = 0.0
    spec['EW_err_cont'] = 0.0
    spec['EW_err_zero'] = 0.0
    spec['ncol_linearCoG'] = 1.13e17*spec['EW']/(spec['fval']*spec['wavc']**2)
    spec['ncol_linear2sig'] = 0.0
    spec['ncol_linear3sig'] = 0.0
    spec['EW_cumulative'] = np.zeros_like(spec['vel'])

    # Detection flags
    if spec['EW'] >= 2*spec['EW_err']:
        print(spec['EW'],spec['EW_err'],type(spec['EW_err']))
        spec['detection_2sig'] = True
    else:
        spec['detection_2sig'] = False
    if spec['EW'] >= 3*spec['EW_err']:
        spec['detection_3sig'] = True
    else:
        spec['detection_3sig'] = False

    spec['ncol'] = spec_in['N']
    spec['ncol_err_lo'] = spec_in['Nsig']
    spec['ncol_err_hi'] = spec_in['Nsig']
    spec['flag_sat'] = False
    #
    spec['va'] = 0.0#spec_in['va'] #average velocity
    spec['va_err'] = 0.0#spec_in['vaerr']
    spec['ba'] = 0.0#spec_in['ba'] #velocity width (b value)
    spec['ba_err'] = 0.0#spec_in['baerr']
    spec['m3'] = 0.0#spec_in['m3'] #skewness
    spec['m3_err'] = 0.0#spec_in['m3err']
    spec['dv90'] = 0.0#spec_in['dv90']
    spec['v90a'] = 0.0#spec_in['v90a']
    spec['v90b'] = 0.0#spec_in['v90b']

    # Create an integration limit if not already available
    if spec['v1'] == spec['v2']:
        spec['v1'] = -100.
        spec['v2'] = +100.
    
    spec = continuum_fit(spec,minord=spec['contin_order']-1,maxord=spec['contin_order']+1)
    spec = pyn_batch(spec, verbose=False, partial_pixels=partial_pixels)

    return spec

def read_inorm(input_filename, partial_pixels=True):
    import numpy as np
    from collections import OrderedDict
    from scipy.io import readsav
    from pyNorm.aod import pyn_batch

    # Read the save file:
    spec_in = readsav(input_filename)

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec_in['vel'].flags.writeable:
        spec_in = __fix_unwriteable_spec(spec_in)


    # Create the dictionary
    spec = OrderedDict()

    # Observation information
    try:
        spec['ion'] = spec_in['ion'].decode('utf-8')
    except:
        spec['ion'] = spec_in['ion']

    try:
        spec['wni'] = spec_in['wni'].decode('utf-8')
    except:
        spec['wni'] = '{0:0.1f}'.format(spec_in['wavc'])

    spec['wavc'] = spec_in['wavc']
    spec['fval'] = spec_in['fval']
    spec['gamma'] = spec_in['gam']
    spec['redshift'] = spec_in['redshift']

    # Sometimes targname key is missing:
    try:
        spec['targname'] = spec_in['targname']
    except:
        spec['targname'] = spec_in['object']

    # Fix binary format
    try:
        spec['object'] = spec_in['object'].decode('utf-8')
    except:
        spec['object'] = spec_in['object']

    try:
        spec['targname'] = spec['targname'].decode('utf-8')
    except:
        spec['targname'] = spec['targname']

    # Coordinates
    spec['RA'] = spec_in['ra']
    spec['Dec'] = spec_in['dec']
    spec['gl'] = spec_in['gl']
    spec['gb'] = spec_in['gb']

    # Sometimes LSR shift is missing
    try:
        spec['vlsr'] = spec_in['vlsr']
    except:
        spec['vlsr'] = lsrvel(spec['gl'],spec['gb'])

    #
    # Raw data
    spec['vel'] = spec_in['vel']
    spec['flux'] = spec_in['flux']
    spec['eflux'] = spec_in['eflux']
    spec['wave'] = spec['wavc']*(spec['vel']/2.998e5)+spec['wavc']
    #
    # Continuum definition
    try:
        spec['contin'] = spec_in['ycon']
        spec['contin_err'] = spec_in['ycon_sig']
    except:
        spec['contin'] = spec_in['cont']
        spec['contin_err'] = spec_in['econt']

    try:
        spec['contin_order'] = spec_in['maxord']
    except:
        spec['contin_order'] = 0

    spec['contin_coeff'] = np.array([0,0.])
    #
    # Construct continuum masks / velocity range
    spec['contin_mask_bits'] = spec_in['mask_cont']
    spec['mask_cont'] = spec_in['mask_cont']
    spec = __convert_inorm_mask(spec)
    #
    # Normalized spectrum
    spec['vnorm'] = spec_in['vnorm']
    spec['fnorm'] = spec_in['fnorm']
    spec['fnorm_err'] = spec_in['efnorm']
    spec['fnorm_err_contin'] = spec_in['efnorm1']*0.
    spec['fnorm_err_stat'] = spec_in['efnorm2']*0.
    #
    # Na(v) data
    spec['Nav'] = np.zeros_like(spec_in['fnorm'])
    spec['Nav_err'] = np.zeros_like(spec_in['fnorm'])
    spec['Nav_sat'] = np.repeat(False,np.size(spec['Nav']))
    #
    # Na(v) quantities
    try:
        spec['SNR'] = spec_in['sn']
    except:
        spec['SNR'] = spec_in['snr']

    spec['v1'] = spec_in['v1']
    spec['v2'] = spec_in['v2']
    spec['EW'] = spec_in['w']
    spec['EW_err'] = spec_in['w_et']
    spec['EW_err_stat'] = spec_in['w_et']
    spec['EW_err_cont'] = spec_in['w_ec']
    spec['EW_err_zero'] = spec_in['w_ez']
    spec['ncol_linearCoG'] = 1.13e17*spec['EW']/(spec['fval']*spec['wavc']**2)
    spec['ncol_linear2sig'] = spec_in['col2sig']
    spec['ncol_linear3sig'] = spec_in['col3sig']
    spec['EW_cumulative'] = np.zeros_like(spec['vel'])

    # Detection flags
    if spec['EW'] >= 2*spec['EW_err']:
        spec['detection_2sig'] = True
    else:
        spec['detection_2sig'] = False
    if spec['EW'] >= 3*spec['EW_err']:
        spec['detection_3sig'] = True
    else:
        spec['detection_3sig'] = False

    spec['ncol'] = spec_in['ncol']
    spec['ncol_err_lo'] = spec_in['ncole1']
    spec['ncol_err_hi'] = spec_in['ncole2']
    spec['flag_sat'] = False
    #
    spec['va'] = spec_in['va']
    spec['va_err'] = spec_in['vaerr']
    spec['ba'] = spec_in['ba']
    spec['ba_err'] = spec_in['baerr']
    spec['m3'] = spec_in['m3']
    spec['m3_err'] = spec_in['m3err']
    spec['dv90'] = spec_in['dv90']
    spec['v90a'] = spec_in['v90a']
    spec['v90b'] = spec_in['v90b']

    # Create an integration limit if not already available
    if spec['v1'] == spec['v2']:
        spec['v1'] = -100.
        spec['v2'] = +100.

    spec = pyn_batch(spec, verbose=False, partial_pixels=partial_pixels)

    return spec

def __fix_unwriteable_spec(spec):
    import numpy as np

    # FIX NON-WRITEABLE ARRAYS due to discontiguous memory
    for kkk in spec.keys():
        if isinstance(spec[kkk],(np.ndarray)):
            spec[kkk] = spec[kkk].copy()

    return spec


def __convert_inorm_mask(spec_in):
    import numpy as np

    spec = spec_in.copy()

    # Extract the mask
    try:
        mask = spec['mask_cont']
    except:
        return spec

    # Create a boolean mask from original
    gdcont = (mask == 1)
    # Find the places where the mask changes.
    delta = mask-np.roll(mask,1)

    # Identify where the mask transitions are:
    vstarts = np.where(delta == -1)[0]
    vstops = np.where(delta == 1)[0]
    # If the first data point is fitted, adjust the continuum boundaries.
    if mask[0] == 1:
        strt = vstops
        stps = vstarts
        vstarts = np.concatenate([np.array([0]), strt])
        vstops = np.concatenate([stps, np.array([np.size(mask)-1])])
    # Transform the results to velocity ranges:
    vstarts = spec['vel'][vstarts]
    vstops = spec['vel'][vstops]

    # # Create a boolean mask from velocity ranges
    # gdcont = np.repeat(False,np.size(spec['vel']))
    # for j in np.arange(np.size(vstarts)):
    #     gdnew = (spec['vel'] >= vstarts[j]) & (spec['vel'] <= vstops[j])
    #     gdcont = gdcont | gdnew

    spec['contin_mask_bits'] = mask
    spec['contin_mask_bool'] = gdcont
    spec['contin_v1'] = vstarts
    spec['contin_v2'] = vstops

    try:
        del spec['mask_cont']
    except:
        pass

    return spec


def lsrvel(long, lat, radec=False, mihalas=False, silent=True):
    """delta_v = lsrvel(long, lat, mihalas=False, SILENT=False):

       This program calculates the projection of the velocity vector of the local standard of
       rest on the sightline specified by (l,b)=(long,lat) used to calculate the shift from
       heliocentric to LSR velocities: v(LSR) = v(helio) + LSR

          Assumes v(LSR) = 20   km/sec to (l,b)=(56, 22) or
                  v(LSR) = 16.5 km/sec to (l,b)=(53, 25)
                                  from Mihalas & Binney

          Created by JCH 9/27/99

    :param long: Longitude [Galactic unless radec=True]
    :param lat: Latitude [Galactic unless radec=True]
    :param radec: input coordinates are RA/Dec? [default: False]
    :param mihalas: use the Mihalas & Binney definition of LSR [default: False]
    :param SILENT: suppress printing (default: False)
    :return: LSR shift, where v(LSR) = v(helio) + LSR
    """

    import numpy as np
    from astropy.coordinates import SkyCoord

    # Radio definition coordinates
    llsr = 56.
    blsr = 22.
    lsr_coords = SkyCoord(llsr, blsr, frame='galactic', unit='deg')
    # Radio definition velocity
    vlsr = 20.0

    # M&B coordinates
    lmb = 53.
    bmb = 25.
    mb_coords = SkyCoord(lmb, bmb, frame='galactic', unit='deg')
    # M&B velocity
    vmb = 16.5

    if radec == False:
        input_coords = SkyCoord(long, lat, frame='galactic', unit='deg')
    else:
        input_coords = SkyCoord(long, lat, frame='icrs', unit='deg')

    # Calculate the separations on the sky [given in degrees]:
    dlsr = input_coords.separation(lsr_coords)
    dmb = input_coords.separation(mb_coords)

    # Calculate the projected velocities:
    vlsr_out = vlsr * np.cos(dlsr)
    vmb_out = vmb * np.cos(dmb)

    if silent == False:
        ## Print output...
        print("\n LSR Correction: ")
        print("     LSR     = {0:0.2f} km/s".format(vlsr_out))
        print("     LSR(MB) = {0:0.2f} km/s.".format(vmb_out))
        print("\n v(LSR) = v(helio) + LSR \n")

    def _ret():
        if mihalas == True:
            deltav_out = vmb_out
        else:
            deltav_out = vlsr_out
        return deltav_out

    return _ret()

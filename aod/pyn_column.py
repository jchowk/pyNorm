import numpy as np

# Assumptions:
#
#  -- Arrays are passed as velocity, flux
#  -- Velocity spacing is constant (enough)

def fix_unwriteable_spec(spec):
    # FIX NON-WRITEABLE ARRAYS due to discontiguous memory
    for kkk in spec.keys():
        if isinstance(spec[kkk],(np.ndarray)):
            spec[kkk] = spec[kkk].copy()

    return spec

def xlimit(x, limits):

   def _ret():  return idx1, idx2

   # x1 = where(ravel(x >= xpos1))[0]  ;  x1 = x1[0]
   # x2 = where(ravel(x > xpos2))[0]  ;  x2 = x2[0] - 1

   idx1 = (np.abs(x - limits[0])).argmin()
   idx2 = (np.abs(x - limits[1])).argmin()

   return _ret()

def integrate_column(velocity, flux, flux_err,
                continuum, continuum_err, wavc, fval,
                integration_limits = [-100,100]):
    # TODO: Asymmetrical error bars?

    # Some constants and flags
    column_factor = 2.654e-15
    flag_sat = False

    # TODO: Adjust integration limits to reflect whole pixels used.
    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        int_idx[xlim1:xlim2] = True
    else:
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        int_idx[xlim1:xlim2] = True

    #Define the velocity spacing
    delv = np.median(velocity[1:]-velocity[:-1])

    # Test for clearly saturated pixels:
    idx_saturation = (flux <= 0.)
    if idx_saturation.sum() > 0: flag_sat = True

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    tau_int = np.sum(tau_array[int_idx]*delv)
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err[int_idx]*delv)**2))

    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err = tau_array_err/(wavc*fval*column_factor)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    column_err_cont = \
        np.sum(((continuum_err[int_idx]/continuum[int_idx])*delv)) /\
         (wavc*fval*column_factor)

    # Background uncertainty
    z_eps = 0.01  # Fractional bg error
    yc1 = continuum[int_idx]*(1.-z_eps)
    y1  = flux[int_idx]-continuum[int_idx]*z_eps
    tau1 = np.sum(np.log(yc1/y1)*delv)
    col1 = tau1 / (wavc*fval*column_factor)
    column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    return column, column_err_total


def pyn_column(spec_in,integration_limits = None):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Some constants and flags
    column_factor = 2.654e-15
    flag_sat = False

    velocity = spec['vel'].copy()
    flux = spec['flux'].copy()
    flux_err = spec['eflux'].copy()
    wavc=spec['wavc'].copy()
    fval=spec['fval'].copy()

    # Deal with the continuum:
    try:
        continuum=spec['ycon'].copy()
        continuum_err = spec['ycon_sig'].copy()
    except:
        continuum=spec['cont']
        continuum_err = spec['econt'].copy()


    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        int_idx[xlim1:xlim2] = True
    else:
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        int_idx[xlim1:xlim2] = True

    #Define the velocity spacing
    delv = np.median(velocity[1:]-velocity[:-1])

    # Test for clearly saturated pixels:
    idx_saturation = (flux <= 0.)
    if idx_saturation.sum() > 0: flag_sat = True

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    tau_int = np.sum(tau_array[int_idx]*delv)
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err[int_idx]*delv)**2))

    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err = tau_array_err/(wavc*fval*column_factor)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    column_err_cont = \
        np.sum(((continuum_err[int_idx]/continuum[int_idx])*delv)) /\
         (wavc*fval*column_factor)

    # Background uncertainty
    z_eps = 0.01  # Fractional bg error
    yc1 = continuum[int_idx]*(1.-z_eps)
    y1  = flux[int_idx]-continuum[int_idx]*z_eps
    tau1 = np.sum(np.log(yc1/y1)*delv)
    col1 = tau1 / (wavc*fval*column_factor)
    column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    spec['v1'] = integration_limits[0]
    spec['v2'] = integration_limits[1]
    spec['vaod1'] = integration_limits[0]
    spec['vaod2'] = integration_limits[1]

    spec['ncol'] = np.log10(column)
    spec['necol1'] = column_err_total/column*np.log10(np.e)
    spec['necol2'] = column_err_total/column*np.log10(np.e)

    return spec


def pyn_eqwidth(spec_in,integration_limits = None):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Some constants and flags
    column_factor = 2.654e-15
    lightspeed = 2.998e5 # km/s
    flag_sat = False

    velocity = spec['vel'].copy()
    flux = spec['flux'].copy()
    flux_err = spec['eflux'].copy()
    wavc=spec['wavc'].copy()
    fval=spec['fval'].copy()

    # Deal with the continuum:
    try:
        continuum=spec['ycon'].copy()
        continuum_err = spec['ycon_sig'].copy()
    except:
        continuum=spec['cont']
        continuum_err = spec['econt'].copy()


    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        int_idx[xlim1:xlim2] = True
    else:
        int_idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        int_idx[xlim1:xlim2] = True

    #Define the velocity spacing
    delv = np.median(velocity[1:]-velocity[:-1])

    # Create the wavelength array
    try:
        wave=spec['wave'].copy()
    except:
        wave = spec['wavc']*(velocity/lightspeed)+wavc
        spec['wave'] = wave

    # Wavelength spacing
    delw = np.median(wave[1:]-wave[:-1])

    # Calculate the equivalent width
    eqw_int = np.sum((1.-flux[int_idx]/continuum[int_idx])*delw)

    # Random flux errors
    eqw_stat_err = \
        np.sqrt(np.sum((flux_err[int_idx]/continuum[int_idx]*delw)**2))
    # Continuum errors
    eqw_cont_err = \
        np.sum(continuum_err[int_idx]*(flux[int_idx]/continuum[int_idx]**2)*delw)
    # Zero point error
    # TODO: Check this calculation
    z_eps = 0.01
    eqw_zero_err = z_eps*eqw_int

    # Combine errors
    eqw_err = np.sqrt(eqw_stat_err**2 \
        +eqw_cont_err**2 + eqw_zero_err**2)

    spec['v1'] = integration_limits[0]
    spec['v2'] = integration_limits[1]

    # Store the EW in milliAngstrom
    spec['EW'] = eqw_int*1000.
    spec['EW_err'] = eqw_err*1000.
    spec['EW_stat_err'] = eqw_stat_err*1000.
    spec['EW_cont_err'] = eqw_cont_err*1000.
    spec['EW_zero_err'] = eqw_zero_err*1000.

    spec['EW_cumulative'] = \
      np.cumsum((1.-flux[int_idx]/continuum[int_idx])*delw)*1000.

    # Delete the old versions of the EW quantities
    try:
        del spec['w']
        del spec['w_es']
        del spec['w_ec']
        del spec['w_et']
        del spec['w_ez']
    except:
        pass

    return spec

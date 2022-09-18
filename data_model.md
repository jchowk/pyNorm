* `ion` - Ion
* `wni` - iNorm wavelength label
* `wavc` - Central wavelength
* `fval` - f-value
* `gamma` - Damping constant $\gamma$
* `redshift` - Assumed absorber redshift [redundant with `z`]
* `targname` - Object name used for observations [in header]
* `object` - Object name [preferred, e.g., from Simbad; may be different than targname]
* `vlsr` - LSR velocity correction
* `RA` - Right ascension
* `Dec` - Declination
* `gl` - Galactic latitude
* `gb` - Galactic longitude
* `vel` - Velocity [LSR?]
* `flux` - Flux
* `eflux` - Flux error
* `wave` - Wavelength array [LSR frame?]
* `contin` - Continuum
* `contin_err` - Error in the continuum points
* `contin_order` - Order of fitted Legendre polynomial
* `contin_coeff` - Coefficients of fitted Legendre polynomial
* `contin_mask_bits` - Continuum mask (bits: 1 = use)
* `contin_mask_bool` - Continuum mask (Boolean)
* `contin_v1` - Starting velocities of continuum regions
* `contin_v2` - Stopping velocities of continuum regions
* `vnorm` - Velocity for normalized spectrum
* `fnorm` - Flux for normalized spectrum
* `fnorm_err` - Error of normalized spectrum
* `fnorm_err_contin` - Continuum error of normalized spectrum
* `fnorm_err_stat` - Statistical error of normalized spectrum
* `Nav` - Apparent column density profile, Na(v) [cm^-2]
* `Nav_err` - Apparent column density error
* `Nav_sat` - Apparent column density saturation flag (True = saturated)
* `SNR` - Median signal-to-noise ratio (SNR) in continuum regions.
* `v1` - Lower [left] integration limit
* `v2` - Upper [right] integration limit
* `EW` - Equivalent width [mA]
* `EW_err` - Equivalent width total error [mA]     
* `EW_err_stat` - Equivalent width statistical error [mA]
* `EW_err_cont` - Equivalent width continuum error [mA]
* `EW_err_zero` - Equivalent width zero point error [mA]
* `EW_cumulative` - Cumulative EW over the line profile [mA]
* `ncol_linearCoG` - log N from linear CoG
* `ncol_linear2sig` - Detection limit at 3$\sigma$ from linear CoG
* `ncol_linear3sig` - Detection limit at 3$\sigma$ from linear CoG
* `detection_2sig` - Flag indicating EW is detected at 2sigma [True/False]
* `detection_3sig` - Flag indicating EW is detected at 3sigma [True/False]
* `ncol` - Log10 Na (apparent column density)
* `ncol_err_lo` - Negative error in log10 Na (apparent column density)
* `ncol_err_hi` - Positive error in log10 Na (apparent column density)
* `flag_sat` - Flag denoting presence of obvious saturation [True/False]
* `va` - Average velocity [first moment] [km/s]
* `va_err` - Average velocity [first moment] error [km/s]
* `ba` - Velocity width [sqrt(2)*second moment] b-value [km/s]
* `ba_err` - Velocity width [sqrt(2)*second moment] b-value error [km/s]
* `m3` - Skewness [third moment]
* `m3_err` - Skewness [third moment] error
* `dv90` - $\Delta v_{90}$: 90% of EW [km/s]
* `v90a` - $v_{90,a}$: 5% of EW velocity limit [km/s]
* `v90b` - $v_{90,b}$: 95% of EW velocity limit [km/s]

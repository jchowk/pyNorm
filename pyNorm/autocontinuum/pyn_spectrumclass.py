from pyNorm.io import read_inorm

# Read data from IDL save file
input_filename = 'Data/CIV1548.2i_o.save'
spec=read_inorm(input_filename)

# Convert spec dictionary to a class containing all the same information. 
# This is necessary for the GP continuum fitting.
class Spectrum:
    def __init__(self, spec):
        self.meta = {key: value for key, value in spec.items() if key not in 
                     ['vel',
                    'flux',
                    'eflux',
                    'wave',
                    'contin',
                    'contin_err',
                    'contin_mask_bits',
                    'contin_mask_bool',
                    'vnorm',
                    'fnorm',
                    'fnorm_err',
                    'fnorm_err_contin',
                    'fnorm_err_stat',
                    'Nav',
                    'Nav_err',
                    'Nav_sat',
                    'EW_cumulative',
                    'integration_weights']}
        self.vel = spec['vel']
        self.flux = spec['flux']
        self.eflux = spec['eflux']
        self.wave = spec['wave']
        self.contin = spec['contin']
        self.contin_err = spec['contin_err']
        self.contin_mask_bits = spec['contin_mask_bits']
        self.contin_mask_bool = spec['contin_mask_bool']
        self.vnorm = spec['vnorm']
        self.fnorm = spec['fnorm']
        self.fnorm_err = spec['fnorm_err']
        self.fnorm_err_contin = spec['fnorm_err_contin']
        self.fnorm_err_stat = spec['fnorm_err_stat']
        self.Nav = spec['Nav']
        self.Nav_err = spec['Nav_err']
        self.Nav_sat = spec['Nav_sat']
        self.EW_cumulative = spec['EW_cumulative']
        self.integration_weights = spec['integration_weights']

# Convert spec dictionary to a Spectrum object
spectrum = Spectrum(spec)




"""
                               iFIND.PRO
                               Version 1.0;
Program Description:
 This procedure finds the ion and atomic parameters given a wavelength

Restrictions:
      Need to define the path to the atomic file.

Screen Output:
       Error text


On Input:
   wave  :== wavelength to search for
   tol :== wavelength tolerance (optional)
#On Output:
   wave_out  :== wavelength of adopted f-value
   ion_out :==  ion for given wavelength
   fval_out :== f-value for given wavelength
   gam_out :== gamma for given wavelength

Common Blocks / Structures:
       None

Latest Update Comments:
       04/12/13  NL: Version 1.0

External Routines called:
       None
"""

from __future__ import print_function

from numpy import *

def iFIND(wave, tol, wave_out, ion_out, fval_out, gam_out):
    """
    Error control.
    """

    n_params = 6
    def _ret():  return (wave, tol, wave_out, ion_out, fval_out, gam_out)
    
    # ON_IOERROR, ESCAPE
    #
    #Find the fvalue of line within a delta wavelength defined by variable tol.
    #
    ## LOOP:
    if bitwise_not((tol is not None)):    
        tol = 0.001
    
    path = '~/Dropbox/IDL/iNorm/lists/'
    RESTORE(path + 'ilines.save')
    
    gg = where(ravel(absolute(wavc - wave) <= tol))[0]
    if ct != 0:    
        wave_out = wavc[gg]
        fval_out = fval[gg]
        ion_out = ion[gg]
        gam_out = gam[gg]
        if ion_out.size > 1:    
            print( 'iFIND::  WARNING: Found more than one transition!')
            for i in arange(0, (ion_out.size - 1)+(1)):
                print( ion_out(i), ' ', array(wave_out(i), copy=0).astype(float64))
            READ('iFIND::  Enter EXACT wavelength....   ', wave)
            tol = 1.e-4
            ## GOTO, ;; LOOP
        # we do not want an array, but a float or string.
    else:    
        print( 'iFIND::  Unable to locate line ', str(wave, '(f8.3)'))
        print( 'iFIND::  WARNING: Cannot find wavelength!')
        print( 'iFIND::  Enter another (w)avelength or (t)olerance, or (q)uit')
        choice = GET_KBRD(1)
        if choice == 'w':    
            READ('iNORM:: Enter wavelength: ', wave)
            ## GOTO, ;; LOOP
        if choice == 't':    
            READ('iNORM:: Enter tolerance: ', tol)
            ## GOTO, ;; LOOP    #else begin
        ## GOTO, ESCAPE
        #endelse
    wave_out = array(wave_out(0), copy=0).astype(float64)
    fval_out = array(fval_out(0), copy=0).astype(float64)
    ion_out = ion_out(0)
    gam_out = array(gam_out(0), copy=0).astype(float64)
    print( 'iFIND::  Located line ', ion_out[0], ' ', str(wave_out[0], '(f8.3)'))
    return _ret()
    #------------------------------------------------------------------------------
    # ESCAPE:
    #   PRINT,'iFIND:: '+!err_string
    #   PRINT,'iFIND:: No wavelength defined! Press (Q)uit!'
    #   RETURN  & 


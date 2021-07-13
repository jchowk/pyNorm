"""
Kenneth Sembach
                               LEGERR.PRO

Created: Unknown
Last Revised:

Program Description:
	This procedure calculates errors for use with LEGFIT.

Restrictions:
	Errors treated with uniform weighting ala LEGFIT.

Screen Output:
       None

Use:
       LEGERR,x,y,a,eps,chi2,error

On Input:
		x       :== abcissa array
		y       :== ordinate array
		a	:== Legendre polynomial coefficient array (from LEGFIT)
		eps	:== Error matrix (from LEGFIT)
		chi2	:== Variance returned by LEGFIT (uniform weighting)

On Ouptut:
               error  :== error array for y

Common Blocks / Structures:
       None

Latest Update Comments:

External Routines Called:
       None
"""

from __future__ import print_function

from numpy import *

def LEGERR(x, y, a, eps, chi2, error):

##  ;; ;; ;; IF N_PARAMS() EQ 0 THEN BEGIN MAN,'legerr' & RETURN & ENDIF

    n_params = 6
    def _ret():  return (x, y, a, eps, chi2, error)
    
    ncoeff = a.size
    nx = x.size
    ix = arange(nx, dtype=float32)
    nord = a.size - 1
    
    p = zeros([ncoeff, nx], "float32")
    p(ix) = 1.
    p(ix + nx) = x
    for j in arange(2, (nord)+(1)):
        p(ix + j * nx) = ((2. * j - 1.) * x * p[j - 1,:] - (j - 1) * p[j - 2,:]) / j
    
    eps1 = chi2 * eps
    tot = zeros([nx], "float32")
    for i in arange(0, (nx - 1)+(1)):
        tot(i) = 0
        for l in arange(0, (nord)+(1)):
            for k in arange(0, (nord)+(1)):
                tot(i) = tot(i) + eps1[k,l] * p[l,i] * p[k,i]
    error = SQRT(tot)
    
    return _ret()


"""
Kenneth Sembach
				LEGPOLY.PRO

Created: 04/23/91
Last Revised: 10/08/92

Program Description:
	This function calcuates the Legendre polynomial corresponding
	to an abcissa vector and a coefficient vector.

Restrictions:
	None

Screen Output:
	None

Use:
	result = LEGPOLY(x,a)

On Input:
		x	:== abcissa array
		a	:== coefficient array

On Ouptut:
		result	:== array containing Legendre polynomial construction

Common Blocks / Structures:
	None

Latest Update Comments:
	10/08/92  KRS	- Updated to run under Version 2 IDL.

External Routines Called:
	None
"""

from __future__ import print_function

from numpy import *

def LEGPOLY(x, coeff):
    """
    
    Array length variables.
    """

    n_params = 2
    def _ret():  return None
    
    nx = x.size
    ix = arange(nx, dtype=float32)
    ncoeff = coeff.size
    nord = ncoeff - 1
    nc = choose(ncoeff < 2, (ncoeff, 2))
    #
    #Form legendre polynomial pieces.
    #
    p = zeros([ncoeff, nx], "float32")
    p(ix) = 1.
    p(ix + nx) = x
    for j in arange(2, (nord)+(1)):
        p[ix + j * nx] = ((2. * j - 1.) * x * p[j - 1,:] - (j - 1) * p[j - 2,:]) / j
    #
    #Add the pieces to form entire polynomial.
    #
    y = zeros([nx], "float32")
    for k in arange(0, (nord)+(1)):
        y = y + coeff[k] * p[k,:]
    #
    #Return to caller.
    #
    return y


"""
Kenneth Sembach
				LEGFIT.PRO

Created: 04/23/91
Last Revised: 05/02/99

Program Description:
	This procedure calculates a Legendre polynomial fit to data.  Derived
       from procedure similar to that given by Bevington (1969).

Restrictions:
	None

Screen Output:
	None

Use:
	LEGFIT,x,y,minord,maxord,yfit,a,eps,chi2

On Input:
		x	:== abscissa array
		y	:== ordinate array
		minord	:== minimum order to start fit
		maxord	:== maximum order to terminate fit

On Output:
		yfit	:== fitted array (over x)
		a	:== coefficients of fit
		eps	:== error matrix
		chi2	:== chi-squared for fit (by definition = sigma^2 here)
                           We use uniform weighting = 1 here, and let sigma
                           be goodness of fit.

Common Blocks / Structures:
	None

Latest Update Comments:
	10/08/92  KRS	- Now runs under Version 2 IDL.
	05/02/99  KRS	- Documentation updates

External Routines Called:
	FTEST		- to check for need for another term in polynomial
"""

from __future__ import print_function

from numpy import *

def LEGFIT(x, y, minord, maxord, yfit, a, eps, chi2):
    """

    Flags and counters.
    """

    n_params = 8
    def _ret():  return (x, y, minord, maxord, yfit, a, eps, chi2)

    nflag = 0
    nord = maxord
    #
    #Array subscript length and vector.
    #
    nx = x.size
    ix = arange(nx, dtype=int32)
    #
    #Form legendre polynomial.
    #
    p = zeros([maxord + 1, nx], "float32")
    p[ix] = 1.
    p[ix + nx] = x
    for j in arange(2., (maxord)+(1)):
        p(ix + j * nx) = ((2. * j - 1.) * x * p[j - 1,:] - (j - 1) * p[j - 2,:]) / j
    #
    #Begin loop to do fit.
    #
    for nord in arange(minord, (maxord)+(1)):

    ## LOOP:
        ncoeff = nord + 1
        #
        #Form alpha and beta matrices.
        #
        beta = zeros([nord + 1], "float32")
        alpha = zeros([nord + 1, nord + 1], "float32")
        for k in arange(0, (nord)+(1)):
            beta(k) = TOTAL(y * p[k,:])
        for k in arange(0, (nord)+(1)):
            for j in arange(0, (nord)+(1)):
                alpha(j, k) = TOTAL(p[j,:] * p[k,:])
        #
        #Invert alpha matrix ==> error matrix eps.
        #
        eps = INVERT(alpha)
        #
        #Calculate coefficients and fit.
        #
        a = transpose(matrixmultiply(transpose(beta), transpose(eps)))
        yfit = zeros([nx], "float32")
        for j in arange(0, (nord)+(1)):
            yfit = yfit + a(j) * p[j,:]
        #
        #Calculate chi squared of fit - uniform weighting=1.
        #
        #		sigma = SQRT(TOTAL((y-yfit)^2)/(nx-ncoeff-1))
        chisq = TOTAL((y - yfit) ** 2)
        chi2 = chisq / (nx - ncoeff - 1)
        #
        #Check chi2 against previous chi2 to see if additional term should
        #be kept.  Check probability for "95% confidence".
        #
        # IF nflag EQ 1 THEN  GOTO,OUT
        if nord > minord:
            f = (chisq1 - chisq) / chi2
            fcutoff = FTEST(nx - ncoeff - 1, 0.05)
            if f < fcutoff:
                nord = nord - 1
                nflag = 1
                ## GOTO,;; LOOP  ;back up to top to do it over again
        chisq1 = chisq
    # OUT:
    # 	maxord = maxord < nord



    #
    #Calculate errors on coefficients - not used here since uniform weighting of
    #data, but could be used someday.
    #
    #	siga = FLTARR(maxord+1)
    #	FOR j=0,maxord DO siga(j) = SQRT(chi2*eps(j,j))

    # RETURN

    return _ret()

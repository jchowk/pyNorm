"""
				iYFIT.PRO
				Version 1.0
Program Description:
	This procedure fits a Legendre polynomial continuum to a spectrum
	and calculates associated errors as outline in Sembach & Savage
	1992.

Restrictions:
       Requires cursor (3 button mouse) input

Screen Output:
	Text  &  Graphics

Use:
	IYFIT,x,y,xarray,yarray,store,ycon,a,sigma,ycon_sig,ftflag

On Input:
		x	:== x coordinate array
		y	:== y coordinate array
	 ey :== ey coordinate array
		xarray	:== x defined continuum array
		yarray  :== y defined continuum array
		store   :== stored continuum region array (2xn)
On Output:
		ycon	:== y calculated continuum array
		a	:== coefficient array
		sigma   :== sigma of fit
		ycon_sig:== error on continuum points

Common Blocks / Structures:
	None

Latest Update Comments:
	04/10/13  NL	- Version 1.0

External Routines Called:
	LEGERR		- to calculate Legenedre polynomial errors
	LEGFIT		- to calculate Legendre polynomial fit
	LEGPOLY		- to construct Legendre polynomial fit
"""

from __future__ import print_function

from numpy import *

def iYFIT(x, y, ey, xarray, yarray, store, ycon, a, sigma, ycon_sig, ftflag):
    """
    
    Error control.
    """

    n_params = 11
    def _ret():  return (x, y, ey, xarray, yarray, store, ycon, a, sigma, ycon_sig, ftflag)
    
    # ON_IOERROR, ESCAPE
    #
    #Overplot defined continuum regions.
    #
    print( 'iYFIT::  Current fit order: ', str(a.size - 1, '(I1)'))
    ## LOOP:
    loadct(39, silent=True)
    liney1 = array([[_sys_cymax / 30.0, -_sys_cymax / 30.0]]) + 0.9 * _sys_cymax
    for k in arange(0, (store.size / 2 - 1)+(1)):
        OPLOT(array([[store[k,0], store[k,0]]]), liney1, color=200)
        OPLOT(array([[store[k,1], store[k,1]]]), liney1, color=200)
        OPLOT(array([[store[k,0], store[k,1]]]), array([[1.0, 1.0]]) * 0.9 * _sys_cymax, color=200)
    #
    #Print heading and instructions.
    #
    READ('iYFIT::  Minimum order of fit (-1=clear,0=quit): ', minord)
    minord = array(minord, copy=0).astype(int32)
    # IF minord EQ 0 THEN RETURN
    if minord < 0:    
        PLOT(x, y)
        ## GOTO,;; LOOP
    READ('iYFIT::  Maximum order of fit: ', maxord)
    maxord = choose(array(maxord, copy=0).astype(int32) < minord, (array(maxord, copy=0).astype(int32), minord))
    print( 'iYFIT::  Working...')
    #
    #Check to be sure nord is not too large for array.
    #
    # IF maxord GE N_ELEMENTS(xarray) THEN GOTO,ESCAPE
    maxord = choose(maxord > 9, (maxord, 9))
    minord = choose(minord > 9, (minord, 9))
    #
    #Scale xarray and x into appropriate range so that Legendre polynomial can be
    #fit over the range from -1 to + 1.
    #
    maxx = array(absolute(xarray), copy=0).max()  ;  _sys_c = 0
    xarray = xarray / maxx
    x = x / maxx
    #
    #Get Legendre polynomial fit and form continuum over all x.
    #
    LEGFIT(xarray, yarray, minord, maxord, yfit, a, eps, chi2)
    ycon = LEGPOLY(x, a)
    #
    #Get sigma from Legendre polynomial fit.  It is just square root of chi2 here.
    #
    sigma = SQRT(chi2)
    #
    #Calculate mean signal to noise (added 09/01/89).
    #
    sn = 0.0
    if sigma != 0.0:    
        sn = AVG(yfit) / sigma
    #
    #Calculate signal to noise at 0 km/s.
    #
    loc0 = SORT(absolute(x)) ; loc0 = loc0(0)
    sn0 = ycon(loc0) / sigma
    #
    #Overplot the derived continuum as a dashed curve.
    #
    _sys_p.linestyle = 2  ;  OPLOT(x * maxx, ycon, color=190)	;  _sys_p.linestyle = 0
    #
    #Calculate the error bars on the continuum due to errors in the coefficients
    #of the Legendre polynomial fit.
    #
    LEGERR(x, y, a, eps, chi2, ycon_sig)	#chi2 = variance (=sigma^2) here
    #
    #Convert x and xarray back into the correct units.
    #
    x = x * maxx
    xarray = xarray * maxx
    #
    #Bound ycon so that plotting errors and arithmetic errors aren't important.
    #
    maxy = array(y, copy=0).max() ; miny = array(y, copy=0).min()
    ycon = choose(choose(ycon < (-absolute(miny) * 2), (ycon, (-absolute(miny) * 2))) > (maxy * 2), (choose(ycon < (-absolute(miny) * 2), (ycon, (-absolute(miny) * 2))), (maxy * 2)))
    #
    #Print the fitting statistics to upper corner of screen.
    #
    print( 'iYFIT::  Information Dump Follows')
    print( '----------------------------------------')
    print( 'Order used = ', maxord)
    print( '# points   = ', str(yarray.size, '(I8)'))
    print( 'RMS sigma  = ', str(sigma, '(E9.3)'))
    print( 'Mean S/N   = ', str(sn, '(F8.2)'))
    print( 'Mean S/N(0)= ', str(sn0, '(F8.2)'))
    print( '----------------------------------------')
    print( 'iYFIT::  End of Information Dump')
    print( 'iYFIT::  Hit any key to continue')
    choice = GET_KBRD(1)
    ftflag = 1
    return _ret()
    #------------------------------------------------------------------------------
    # ESCAPE:
    # 	PRINT,'iYFIT::  '+!err_string
    # 	RETURN  & 
    


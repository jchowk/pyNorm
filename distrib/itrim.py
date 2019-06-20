"""
				ITRIM.PRO
				Version 1.0

Program Description:
	This program trims a spectrum and associated continuum and error bars.

Screen Output:
	Text & Graphics

Use:
	IMTRIM,x,y,ycon,lbar,ubar,coflag,ebflag,updates

On Input:
		x	:== x coordinate array
   y :== y coordinate array
   ey :== ey coordinate array
   ycon  :== continuum coordinate array
   ycon_sig  :== continuum coordinate array
		coflag	:== continuum fitting flag (0=no,1=yes)
On Output:
		x	:== trimmed x coordinate array
		y	:== trimmed y coordinate array
   ey :== ey coordinate array
		ycon 	:== trimmed continuum coordinate array
   ycon_sig  :== continuum coordinate array

Common Blocks / Structures:
	None

Latest Update Comments:
	04/12/13  NL	- Version 1.0

External Routines called:
	XLIMIT		- to determine trimming elements
"""

from __future__ import print_function

from numpy import *

def iTRIM(x, y, ey, ycon, ycon_sig, coflag):
    """
    
    Error control.
    """

    n_params = 6
    def _ret():  return (x, y, ey, ycon, ycon_sig, coflag)
    
    # ON_IOERROR, ESCAPE
    #
    #Print heading and get left trim limit.
    #
    print( 'iTrim::  (c)ursor trim   (k)eyboard trim   (q)uit')
    ## LOOP:
    choice = GET_KBRD(1)
    #
    #Quit if user wants to do so.
    #
    if choice == 'q':    
        return _ret()
    #
    #Do the trim in the x-coordinate.
    #
    loadct(39, silent=True)
    ## LOOP1:
    if choice == 'c':    
        print( 'iTrim::  Mark (C1,C2,C3)')
        CURSOR(xpos1, ypos1, DOWN=True)
        print( 'iTrim::  Marked left limit:  ', xpos1)
        _sys_p.psym = 1 ; OPLOT(array([[xpos1, xpos1]]), array([[ypos1, ypos1]]), color=230, thick=3, symsize=2) ; _sys_p.psym = 10
        CURSOR(xpos2, ypos2, DOWN=True)
        print( 'iTrim::  Marked right limit: ', xpos2)
    else:    
        if choice == 'k':    
            READ('iTrim::  Enter left limit:  ', xpos1)
            READ('iTrim::  Enter right limit: ', xpos2)
        else:    
            print( 'iTrim::  (c)ursor trim  (k)eyboard trim  (q)uit')
            ## GOTO,;; LOOP
    
    
    xpos1 = choose(choose(xpos1 < array(x, copy=0).min(), (xpos1, array(x, copy=0).min())) > array(x, copy=0).max(), (choose(xpos1 < array(x, copy=0).min(), (xpos1, array(x, copy=0).min())), array(x, copy=0).max()))
    xpos2 = choose(choose(xpos2 > array(x, copy=0).max(), (xpos2, array(x, copy=0).max())) < array(x, copy=0).min(), (choose(xpos2 > array(x, copy=0).max(), (xpos2, array(x, copy=0).max())), array(x, copy=0).min()))
    
    if xpos1 > xpos2:    
        print( 'iTrim::  ' + 'Limits will be reversed for trimming')
        xtmp = xpos1 ; xpos1 = xpos2 ; xpos2 = xtmp
    #
    #Find the trim limits.
    #
    XLIMIT(x, xpos1, xpos2, x1, x2)
    
    if absolute(x1 - x2) <= 1:    
        print( 'iTrim::  ' + 'Insufficient spectral range remaining.')
        print( 'iTrim::  Please re-enter limits.')
        ## GOTO,;; LOOP1
    #
    #Ask if limits are acceptable.  If not, return.
    #
    PLOT(x[x1:(x2)+1], y[x1:(x2)+1])
    READ('iTrim::  Is the trimmed spectrum acceptable? ', choice)
    if STRMID(choice, 0, 1) != 'y':    
        print( 'iTrim::  Spectrum untrimmed')
        return _ret()
    #
    #Trim the spectrum.
    #
    x = x[x1:(x2)+1]
    y = y[x1:(x2)+1]
    ey = ey[x1:(x2)+1]
    #
    #Trim the continuum and error bar arrays if they exist.  Return if only two
    #parameters are passed (ie., called outside IMNORM).
    #
    if n_params == 4:    
        return _ret()
    if coflag == 1:    
        ycon = ycon[x1:(x2)+1]
        ycon_sig = ycon_sig[x1:(x2)+1]
    #
    #Update message to be put into file header and return.
    #
    npts = x.size
    comment = 'iTrim::  Spectrum trimmed to ' + str(npts, '(I5)') + ' points  ' + _sys_stime
    return _ret()
    #------------------------------------------------------------------------------
    # ESCAPE:
    # 	PRINT,'iTrim::  '+!err_string
    # 	RETURN  & 


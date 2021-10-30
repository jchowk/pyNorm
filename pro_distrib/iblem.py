"""
Kenneth Sembach
                               iBLEM.PRO
                               Version 5.2
Created: 09/01/89
Last Revision: 05/02/99

Program Description:
       This program linearly interpolates across spectral regions identified
	by the user with the mouse.

Restrictions:
       3-button mouse assumed

Screen Output:
       Text & Graphics

Use:
       IMBLEM,x,y,yorig

On Input:
               x       :== x coordinate array
               y       :== y coordinate array
               ey       :== ey coordinate array
On Output:
               y       :== blemish corrected y coordinate array
               yorig   :== original y coordinate array

Common Blocks / Structures:
       None

Latest Update Comments:
       04/11/13  NL   - Version 1.0

External Routines called:
	XLIMIT
"""

from __future__ import print_function

from numpy import *

def iBLEM(x, y, ey, yorig):
    """
    
    Error control.
    
    """

    n_params = 4
    def _ret():  return (x, y, ey, yorig)
    
    yorig = y
    ## LOOP:
    #
    #Get input from user.  If right mouse button is pushed, the return.
    #If middle button is pushed, then clear.
    #
    print( 'IMBLEM(v5.2)::  Mark (C1)  Reset (C2)  Quit (C3)')
    CURSOR(xpos1, ypos1, DOWN=True)
    xpos1 = choose(xpos1 < array(x, copy=0).min(), (xpos1, array(x, copy=0).min()))  ;  _sys_c = 0
    if _sys_err == 4:    
        if n_params == 4:    
            loc = where(ravel(y != yorig))[0]
            if loc(0) != -1:    
                npts = loc.size
                comment = 'IMBLEM(v5.2)::  ' + STRTRIM(npts, 2) + ' blemished points removed  ' + _sys_stime
        # RETURN
    
    if _sys_err == 2:    
        y = yorig
        PLOT(x, y)
        ## GOTO,;; LOOP
    print( 'IMBLEM(v5.2)::  Left limit:  ', xpos1)
    OPLOT(array([[1, 1]]) * xpos1, array([[1, 1]]) * ypos1, psym=1)
    #
    #Get and mark right of region.
    #
    CURSOR(xpos2, ypos2, DOWN=True)
    xpos2 = choose(xpos2 > array(x, copy=0).max(), (xpos2, array(x, copy=0).max()))  ;  _sys_c = 0
    print( 'IMBLEM(v5.2)::  Right limit: ', xpos2)
    OPLOT(array([[1, 1]]) * xpos2, array([[1, 1]]) * ypos2, psym=1)
    #
    #Check to be sure range is not reversed.  If reversed, get again.
    #
    if xpos2 <= xpos1:    
        print( 'IMBLEM(v5.2)::  Intended action unclear - try again.')
    #
    #Determine which elements to interpolate and do interpolation.
    #
    XLIMIT(x, xpos1, xpos2, x1, x2)
    xint = x[x1:(x2)+1]
    yint = INTERPOL(array([[y(x1), y(x2)]]), array([[x(x1), x(x2)]]), xint)
    eyint = INTERPOL(array([[ey(x1), ey(x2)]]), array([[x(x1), x(x2)]]), xint)
    y[x1:(x2)+1] = yint
    ey[x1:(x2)+1] = eyint
    OPLOT(xint, yint)
    ## GOTO,;; LOOP
    #----------------------------------------------------------------------------
    # ESCAPE:
    # PRINT,'IMBLEM(v5.2)::  '+!err_string
    
    return _ret()


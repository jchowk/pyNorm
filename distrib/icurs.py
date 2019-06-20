"""
				iCURS.PRO

Program Description:
	This procedure reads the cursor position from the current X window.
	Coordinates are read in data coordinates.

Restrictions:
	Only works with 3 button mouse.  Use right mouse button to exit.

Screen Output:
	Graphics  &  Text

Use:
	iCURS

On Input:
	None

On Output:
	None

Common Blocks / Structures:
	None

Latest Update Comments:
	04/11/13  NL	- Version 1.0

External Routines called:
	None
"""

from __future__ import print_function

from numpy import *

def iCURS():
    """
    Error control
    
     ON_IOERROR,ESCAPE
    
    Print heading.
    """

    n_params = 0
    def _ret():  return None
    
    print( 'iCURS::  Mark (C1,C2)   Quit (C3)')
    print( 'iCURS::  Information Dump Follows')
    print( '----------------------------------------')
    print( '          X             Y             ')
    
    ## LOOP:
    #
    #Get cursor position in data coordinates.  Use the /DOWN qualifier.
    #
    ## ;; ;; IF N_PARAMS() EQ 0 THEN CURSOR,xpos,ypos,/DOWN
    if n_params != 0:    
        CURSOR(xpos, ypos, DOWN=True, DEVICE=True)
    #
    #If the right mouse button is pressed, then exit.  Otherwise print position
    #and plot a crosshair on the plot.
    #
    if _sys_err == 4:    
        print( '----------------------------------------')
        print( 'iCURS::  End of Information Dump')
        return _ret()
    if absolute(xpos) <= 1.e4:    
        format = '$(4x,f10.4,e)'
        print( format, xpos, ypos)
    else:    
        print( xpos, ypos)
    
    ##IF N_PARAMS() EQ 0 THEN BEGIN
    _sys_p.linestyle = 2
    nsum = _sys_nsum ; _sys_nsum = 1
    OPLOT(_sys_x.crange, array([[ypos, ypos]]))
    OPLOT(array([[xpos, xpos]]), _sys_y.crange)
    _sys_p.linestyle = 0
    _sys_nsum = nsum
    # ENDIF
    
    ## GOTO,;; LOOP
    #------------------------------------------------------------------------------
    # ESCAPE:
    # 	PRINT,'iCURS:: '+!err_string
    
    return _ret()


"""
Kenneth Sembach
				XLIMIT.PRO
				Version 6.0
Created: 09/01/89
Last Revision:	02/27/95

Program Description:
	This procedure finds the end pixels (x1 & x2) of array x between the
	values xpos1 & xpos2.

Screen Output: None

Use:
	XLIMIT,x,xpos1,xpos2,x1,x2

On Input:
		x	:== xarray
		xpos1	:== left data limit
		xpos2	:== right data limit

On Output:
		x1	:== left x pixel
		x2	:== right x pixel

Common Blocks / Structures:
	None

Latest Update Comments:
	10/15/92  KRS	- Version 5.0, IMLIMIT renamed to XLIMIT
	02/27/95  KRS	- Stupid conditional removed.  No more changing of
			   xpos1 and xpos2.

External Routines called:
	None
"""

from __future__ import print_function

from numpy import *

def XLIMIT(x, xpos1, xpos2, x1, x2):

##  ;; ;; ;; IF N_PARAMS() EQ 0 THEN BEGIN MAN,'xlimit' & RETURN & ENDIF

    n_params = 5
    def _ret():  return (x, xpos1, xpos2, x1, x2)
    
    x1 = where(ravel(x >= xpos1))[0]  ;  x1 = x1(0)
    x2 = where(ravel(x > xpos2))[0]  ;  x2 = x2(0) - 1
    
    if x1 <= -1:    
        x1 = 0
    if x2 <= -1:    
        x2 = x.size - 1
    
    return _ret()  ;

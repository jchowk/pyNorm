from numarray import *

#+
#Kenneth Sembach
#				XLIMIT.PRO
#				Version 6.0
#Created: 09/01/89
#Last Revision:	02/27/95
#
#Program Description:
#	This procedure finds the end pixels (idx1 & x2) of array x between the
#	values xpos1 & xpos2.
#
#Screen Output: None
#
#Use:
#	XLIMIT,x,xpos1,xpos2,idx1,x2
#
#On Input:
#		x	:== xarray
#		xpos1	:== left data limit
#		xpos2	:== right data limit
#
#On Output:
#		idx1	:== left x pixel
#		x2	:== right x pixel
#
#Common Blocks / Structures:
#	None
#
#Latest Update Comments:
#	10/15/92  KRS	- Version 5.0, IMLIMIT renamed to XLIMIT
#	02/27/95  KRS	- Stupid conditional removed.  No more changing of
#			   xpos1 and xpos2.
#
#External Routines called:
#	None
#----------------------------------------------------------------------------
def xlimit(x, xpos1, xpos2):

   def _ret():  return idx1, idx2

   # x1 = where(ravel(x >= xpos1))[0]  ;  x1 = x1[0]
   # x2 = where(ravel(x > xpos2))[0]  ;  x2 = x2[0] - 1

   idx1 = (np.abs(x - xpos1)).argmin()
   idx2 = (np.abs(x - xpos2)).argmin()

   return _ret()

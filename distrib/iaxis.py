"""
Kenneth Sembach
				iAXIS.PRO
				Version 5.0

Program Description:
	This function converts between wavelength and velocity space.
	If an error is encountered, no message is printed but in_flag is
	set to zero.

Restrictions:
	None

Screen Output:
	None

Use:
	result = iAXIS(x,wavc,in_flag)

On Input:
		x	:== x coordinate array
		wavc	:== laboratory wavelength of line
		in_flag	:== (-1=lambda->velocity, +1=velocity->lambda)
On Output:
		result	:== converted x coordinate array
		in_flag :== 0 if wavc = 0

Common Blocks / Structures:
	None

External Routines called:
	None
"""

from __future__ import print_function

from numpy import *

def iAXIS(x, wavc, in_flag):
    """
    Speed of light in km/sec.
    """

    n_params = 3
    def _ret():  return None
    
    c = 2.997925e5
    #
    #Return if no wavelength is defined .
    #
    if (bitwise_or((wavc == 0), (absolute(in_flag) != 1))):    
        in_flag = 0  ;  return x
    #
    #Covert between lambda and velocity.
    #
    if in_flag == -1:    
        return c * (x - wavc) / wavc	#Lambda to velocity.
    if in_flag == +1:    
        return wavc * (x / c) + wavc 	#Velocity to lambda.
    
    return x  ;

"""
				IMATH.PRO
				Version1.0

Description:
	This routine performs various mathematical operations on a
	spectrum. Limited use.

Screen Output: Graphics text

Use:
	IMMATH,x,y

On Input:
		x	:== x coordinate array
		y	:== y coordinate array

On Output:
		y	:== y coordinate array after mathematical operation

Common Blocks:
	None

Latest Update Comments:
	04/13/13  NL	- Version 1.0

External Routines called:

"""

from __future__ import print_function

from numpy import *

def IMATH(x, y):
    """
    
    Error control.
    
     ON_IOERROR,ESCAPE
    
    Initialize
    """

    n_params = 2
    def _ret():  return (x, y)
    
    com = STRARR(80, 1) ; j = 0
    #
    #Print heading.
    #
    ## LOOP:
    print( 'iMATH::  Information Dump Follows')
    print( 'iMATH::  ----------------------------------------')
    print( 'iMATH::  x = (' + STRTRIM(array(x, copy=0).min(), 2) + ' , ' + STRTRIM(array(x, copy=0).max(), 2) + ')')
    print( 'iMATH::  y = (' + STRTRIM(array(y, copy=0).min(), 2) + ' , ' + STRTRIM(array(y, copy=0).max(), 2) + ')')
    print( 'iMATH::  AVG = ' + STRTRIM(AVG(y), 2))
    print( 'iMATH;;  MED = ' + STRTRIM(MEDIAN(y), 2))
    print( 'iMATH::  ----------------------------------------')
    print( 'iMATH::  End of Information Dump')
    print( 'iMATH::  (a)bs  (e)xp  (i)nverse  (l)og  (n)atural log')
    print( 'iMATH::  (s)quare root (r)eplace value  (R)eset')
    print( 'iMATH::  (+)Add constant  (*)Multiply constant  (q)uit')
    #
    #Get operator.
    #
    ## LOOP1:
    operator = GET_KBRD(1)
    ysave = y
    loc = where(ravel(y != 0))[0]
    _expr = operator
    #
    #Absolute value of y.
    #
    if _expr == a:    
        com(j) = 'iMATH::  Absolute value of spectrum taken  ' + _sys_stime
        y = absolute(y)
        #
        #Base e exponential of y.
        #
    elif _expr == e:    
        com(j) = 'iMATH::  Exponential of spectrum taken  ' + _sys_stime
        y = EXP(y)
        #
        #Reciprocal of y.
        #
    elif _expr == i:    
        com(j) = 'iMATH::  Reciprical of spectrum taken  ' + _sys_stime
        y(loc) = 1.0 / y(loc)
        #
        #Base 10 log of y > 0.
        #
    elif _expr == l:    
        com(j) = 'iMATH::  Base 10 log of spectrum taken ' + _sys_stime
        y = choose(y < 0, (y, 0))
        y(loc) = log10(y(loc))
        #
        #Base e log of y > 0
        #
    elif _expr == n:    
        com(j) = 'iMATH::  Base e log of spectrum taken ' + _sys_stime
        y = choose(y < 0, (y, 0))
        y(loc) = log(y(loc))
        #
        #Replace values.
        #
    elif _expr == r:    
        cutoff = 0.0
        READ('iMATH::  Enter lower cutoff: ', cutoff)
        y = choose(y < cutoff, (y, cutoff))
        com(j) = 'iMATH::  Y-Cutoff of ' + STRTRIM(cutoff, 2) + ' imposed  ' + _sys_stime
        #
        #Square root of y.
        #
    elif _expr == s:    
        com(j) = 'iMATH::  Square root of spectrum taken ' + _sys_stime
        y = choose(y < 0, (y, 0))
        y = SQRT(y)
        
    elif _expr == *:    
        constant = 0.0
        READ('iMATH::  Enter multiplication constant: ', constant)
        y = y * constant		#Multiply y by constant.
        com(j) = 'iMATH::  Constant ' + STRTRIM(constant, 2) + ' mulitplied ' + _sys_stime
        
    elif _expr == +:    
        constant = 0.0
        READ('iMATH::  Enter additive constant: ', constant)
        y = y + constant		#Add constant to y.
        com(j) = 'iMATH::  Constant ' + STRTRIM(constant, 2) + ' added ' + _sys_stime
        
    elif _expr == q:    
        return _ret()
        
    elif _expr == R:    
        y = ysave
        print( 'iMATH::  Spectrum restored')
        com = STRARR(80, 1)
        j = -1
        
    else:    
        print( 'iMATH::  Invalid command: ' + operator)
        ## GOTO,;; LOOP1
    
    print( ' ')
    PLOT(x, y)
    j = j + 1
    ## GOTO,;; LOOP
    #----------------------------------------------------------------------------
    # ESCAPE:
    #         PRINT,'iMATH::  '+!err_string
    #         RETURN
    
    
    return _ret()


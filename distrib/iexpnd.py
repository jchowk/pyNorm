"""
				IEXPND.PRO
				Version 1.0

Program Description:
	This procedure allows interactive expansion of a plot in the current
	X window. Limited use for analyzing the spectrum.

Screen Output:
	Text  &  Graphics

Use:
	IMEXPND,x,y

On Input:
		x	:== x coordinate array
		y 	:== y coordinate array
On Output:
	None

Common Blocks / Structures:
	None

Latest Update Comments:
	04/13/13  - NL: version 1.0

External Routines Called:
	None
"""

from __future__ import print_function

from numpy import *

def IEXPND(x, y):
    """
    
    Error control.
    
     ON_IOERROR,ESCAPE
    
    Print heading and commands.
    
    """

    n_params = 2
    def _ret():  return (x, y)
    
    print( 'IExpnd::  (c)ursor  (k)eyboard   (u)nexpand   (q)uit')
    print( 'IExpnd::  [1-0] fast-x  [!-)] fast-y')
    choice = GET_KBRD(1)
    #
    #Fast scaling of x-axis.
    #
    if choice == '1':    
        xpos1 = -100  ;  xpos2 = 100  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '2':    
        xpos1 = -200  ;  xpos2 = 200  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '3':    
        xpos1 = -300  ;  xpos2 = 300  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '4':    
        xpos1 = -400  ;  xpos2 = 400  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '5':    
        xpos1 = -500  ;  xpos2 = 500  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '6':    
        xpos1 = -600  ;  xpos2 = 600  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '7':    
        xpos1 = -700  ;  xpos2 = 700  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '8':    
        xpos1 = -800  ;  xpos2 = 800  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '9':    
        xpos1 = -900  ;  xpos2 = 900  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    if choice == '0':    
        xpos1 = -1000  ;  xpos2 = 1000  ;  ypos1 = 0  ;  ypos2 = 0
        ## GOTO,FASTOUT
    #
    #Fast scaling of y-axis.
    #
    if choice == '!':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -100  ;  ypos2 = 100
        ## GOTO,FASTOUT
    if choice == '@':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -200  ;  ypos2 = 200
        ## GOTO,FASTOUT
    if choice == '#':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -300  ;  ypos2 = 300
        ## GOTO,FASTOUT
    if choice == '$':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -400  ;  ypos2 = 400
        ## GOTO,FASTOUT
    if choice == '%':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -500  ;  ypos2 = 500
        ## GOTO,FASTOUT
    if choice == '^':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -600  ;  ypos2 = 600
        ## GOTO,FASTOUT
    if choice == '&':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -700  ;  ypos2 = 700
        ## GOTO,FASTOUT
    if choice == '*':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -800  ;  ypos2 = 800
        ## GOTO,FASTOUT
    if choice == '(':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -900  ;  ypos2 = 900
        ## GOTO,FASTOUT
    if choice == ')':    
        xpos1 = 0   ;  xpos2 = 0  ;  ypos1 = -1000  ;  ypos2 = 1000
        ## GOTO,FASTOUT
    #
    #Quit if asked to do so.
    #
    if choice == 'q':    
        return _ret()
    #
    #Unexpand spectrum if asked to do so.
    #
    if choice == 'u':    
        print( 'IExpnd::  Axes will be unexpanded')
        _sys_x.style = 0  ;  _sys_y.style = 0
        _sys_x.range = 0  ;  _sys_y.range = 0
        return _ret()
    #
    #Use cursor or keyboard to mark limits.
    #
    if choice == 'c':    
        print( 'IExpnd::  Axes will be expanded by cursor')
        print( 'IExpnd::  Mark lower left limit')
        CURSOR(xpos1, ypos1, DOWN=True)
        print( 'IExpnd::  Marked:    (', xpos1, ypos1, '    )')
        _sys_p.psym = 1 ; OPLOT(array([[xpos1, xpos1]]), array([[ypos1, ypos1]])) ; _sys_p.psym = 0
        print( 'IExpnd::  Mark upper right limit')
        CURSOR(xpos2, ypos2, DOWN=True)
        print( 'IExpnd::  Marked:    (', xpos2, ypos2, '    )')
    else:    
        if choice == 'k':    
            print( 'IExpnd::  Axes will be expanded by keyboard')
            READ('IExpnd::  Enter xmin,ymin: ', xpos1, ypos1)
            _sys_p.psym = 1 ; OPLOT(array([[xpos1, xpos1]]), array([[ypos1, ypos1]])) ; _sys_p.psym = 0
            READ('IExpnd::  Enter xmax,ymax: ', xpos2, ypos2) #ELSE ;; GOTO,;; LOOP
    #
    #Check to be sure expansion is ok.
    #
    loc = where(ravel(bitwise_and((x >= xpos1), (x <= xpos2))))[0]
    _sys_psym = 10
    
    if (bitwise_and((loc(0) == -1), (xpos1 != xpos2))):    
        print( 'IExpnd::  Unable to x-expand plot')
        _sys_x.range = array([[0, 0]])
        _sys_psym = 10
        ## GOTO,;; LOOP
    loc = where(ravel(bitwise_and((y >= ypos1), (y <= ypos2))))[0]
    if (bitwise_and((loc(0) == -1), (ypos1 != ypos2))):    
        print( 'IExpnd::  Unable to y-expand plot')
        _sys_y.range = array([[0, 0]])
        _sys_psym = 10
        ## GOTO,;; LOOP
    # FASTOUT:
    #
    #Do the expansion according to limits defined above and return to caller.
    #
    minx = array(x, copy=0).min()
    maxx = array(x, copy=0).max()
    
    if (bitwise_or((bitwise_and((xpos1 < minx), (xpos2 < minx))), (bitwise_and((xpos1 > maxx), (xpos2 > maxx))))):    
        print( 'IExpnd:: X-axis expansion impossible')
        return _ret()
    else:    
        _sys_x.style = 1  ;  _sys_y.style = 1
        _sys_x.range = array([[xpos1, xpos2]])
        _sys_y.range = array([[ypos1, ypos2]])
        _sys_psym = 10
    # RETURN
    #------------------------------------------------------------------------------
    # ESCAPE:
    # 	PRINT,'IExpnd::  '+!err_string
    # 	RETURN  &  END
    
    return _ret()


"""
                               iCOL_PREP.PRO
                                 Version 1.0

Program Description:
       This procedure calculates the apparent column density, equivalent width,
       and line statistics.

Restrictions:
       None

Screen Output:
       Text  &  Graphics

Use:
       IMCOL_PREP,x,y,ex,ycon,sigma,ycon_sig,wavc,fval

On Input:
               x       :== velocity array
               y       :== intensity array
               ey      :== intensity error array
               ycon    :== continuum array
               sigma   :== RMS sigma of continuum fit
               ycon_sig:== continuum fit error array
               wavc    :== laboratory wavelength of line
               fval    :== laboratory f-value of line

On Ouptut:
               storing inputs and outputs with possibility to save them at
               the end the session

Common Blocks / Structures:
       SHARE: inputs and outputs.
       NOISE
       FLAGSAT

Latest Update Comments:
 04/10/13  NL    - Version 1.0
 11/16/16  NL    - included errors on lower limits.

External Routines Called:
       iCOL            - to calculate column density and errors
       iEQW            - to calculate equivalent width and errors
       iSTAT2          - to calculate line statistics and errors
       XLIMIT          - to determine elements of integration
       PLOTCOLORFILL   - to plot nice fill colored histograms
"""

from __future__ import print_function

from numpy import *

def iCOL_PREP(x, y, ey, ycon, sigma, ycon_sig, wavc, fval):
    n_params = 8
    def _ret():  return (x, y, ey, ycon, sigma, ycon_sig, wavc, fval)
    
    global ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err, v1, v2, col3sig, dv90, v90a, v90b
    global set_noise, eyflag
    global flag_sat
    
    # loadct,39,/silent
    ## IF N_PARAMS() EQ 0 THEN BEGIN MAN,'iCcol_prep' & RETURN & ENDIF
    #
    #Error control
    #
    # ON_IOERROR,ESCAPE
    #
    #Find the appropriate fvalue from a standard line list.
    #
    ## LOOP:
    #
    #Plot the column density versus velocity.
    #
    _sys_ytitle = 'Apparent Column Density (cm!u-2!n/(km/s))'
    ycol = zeros([y.size], "float32")
    for i in arange(0, (y.size - 1)+(1)):
        if bitwise_and((y[i] > 0), (ycon[i] > 0)):    
            ycol[i] = log(ycon[i] / y[i]) / (wavc * fval * 2.654e-15)
    PLOT(x, ycol)
    #
    #Overplot zero line.
    #
    OPLOT(_sys_x.crange, array([[0, 0]]), LINESTYLE=2, color=190, thick=2)
    #
    #Print wavelength, f-value, and error type to screen for a user check.
    #
    print( 'iCOL_PREP::  Wavelength = ', str(wavc, '(f8.3)'))
    print( 'iCOL_PREP::  f-value = ', str(fval, '(f7.5)'))
    #
    #Type of input.
    #
    print( 'iCOL_PREP::  (c)ursor   (k)eyboard')
    choice1 = GET_KBRD(1)
    #
    #Get limits of integration regardless of which mode is used.
    #
    if choice1 == 'c':    
        print( 'iCOL_PREP::  Mark (C1)   Clear (C2)   Quit(c3)')
        
        CURSOR(xpos1, ypos1, DOWN=True)  ;  if _sys_err == 4:    
            return _ret()
        if _sys_err == 2:    
            PLOT(x, y)  ;  ;; GOTO,;; LOOP
        xpos1 = choose(xpos1 < array(x, copy=0).min(), (xpos1, array(x, copy=0).min()))  ;  _sys_c = 0
        print( 'iCOL_PREP::  Left limit:  ', xpos1, IAXIS(xpos1, wavc, +1))
        CURSOR(xpos2, ypos2, DOWN=True)  ;  if _sys_err == 4:    
            return _ret()
        if _sys_err == 2:    
            PLOT(x, y)  ;  ;; GOTO,;; LOOP
        xpos2 = choose(xpos2 > array(x, copy=0).max(), (xpos2, array(x, copy=0).max()))  ;  _sys_c = 0
        print( 'iCOL_PREP::  Right limit: ', xpos2, IAXIS(xpos2, wavc, +1))
    else:    
        READ('iCOL_PREP::  Enter left limit (v):  ', xpos1)
        READ('iCOL_PREP::  Enter right limit (v): ', xpos2)
    
    if xpos1 > xpos2:    
        print( 'iCOL_PREP::  ' + 'Limits will be reversed for integration')
        xtmp = xpos1 ; xpos1 = xpos2 ; xpos2 = xtmp
    #
    #Compute which range of elements should be included in integration and fill
    #in line.  XLIMIT will find points in between limits, but be sure to set
    #endpoints of integration equal to the limits.
    #
    v1 = xpos1
    v2 = xpos2
    
    XLIMIT(x, xpos1, xpos2, x1, x2)
    if absolute(x1 - x2) <= 1:    
        print( 'iCOL_PREP::  ' + 'Insufficient spectral range specified.')
        print( 'iCOL_PREP::  Please re-enter limits.')
        ## GOTO,;; LOOP
    xwork = array([[xpos1, x[x1:(x2)+1], xpos2]])
    ywork = array([[INTERPOL(y, x, xpos1), y[x1:(x2)+1], INTERPOL(y, x, xpos2)]])
    eywork = array([[INTERPOL(ey, x, xpos1), ey[x1:(x2)+1], INTERPOL(ey, x, xpos2)]])
    yconwork = array([[INTERPOL(ycon, x, xpos1), ycon[x1:(x2)+1], INTERPOL(ycon, x, xpos2)]])
    ycolwork = array([[INTERPOL(ycol, x, xpos1), ycol[x1:(x2)+1], INTERPOL(ycol, x, xpos2)]])
    ycon_sig_work = array([[INTERPOL(ycon_sig, x, xpos1), ycon_sig[x1:(x2)+1]]])
    ycon_sig_work = array([[ycon_sig_work, INTERPOL(ycon_sig, x, xpos2)]])
    
    #
    #Ask about what kind of statistics should be used (poisson or fixed pattern)
    #to describe the noise characteristics of the data.
    #
    #
    if eyflag != 0:    
        print( 'iCOL_PREP::  Error Vector provided')
        y_sig_work = eywork
        kind = 'u'
    else:    
        if sigma != 0:    
            print( 'iCOL_PREP::  Noise:  (p)oisson   (f)ixed pattern')
            kind = GET_KBRD(1)
            if kind == 'f':    
                print( 'iCOL_PREP::  Fixed pattern noise assumed')
                y_sig_work = sigma + zeros([ywork.size], "float32")
            else:    
                print( 'iCOL_PREP::  Poisson noise assumed')
                y_sig_work = sigma * SQRT(absolute(ywork) / yconwork)
        else:    
            print( 'iCOL_PREP:: The errors are not defined')
            print( 'iCOL_PREP:: No error will be estimated')
            y_sig_work = sigma + zeros([ywork.size], "float32")
            kind = 'f'
    # save noise variable
    set_noise = kind
    
    #
    #Calculate the column density and error by calling iCOL.
    #
    iCOL(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, wavc, fval, col, y_err, ycon_err, zero_err)
    tot_err = SQRT(y_err ** 2 + ycon_err ** 2)
    if col <= 0.0:    
        col = 1.0
    #
    #FIll in the integrated area for the viewer to see.
    #
    plotcolorfill(xwork, ycolwork, midpoint=True, noerase=True, col=234, bottom=0)
    
    #
    #Print information dump.
    #
    if flag_sat == 0:    
        print( '--------------------------------------------')
        print( 'log N (best) =   ', str(log10(col), '(F9.3)'))
        print( '(+1 sig)     =   ', str(log10(col + tot_err) - log10(col), '(F9.3)'))
        print( '(-1 sig)     =   ', str(-log10(col - tot_err) + log10(col), '(F9.3)'))
        #	PRINT,'2% zero err  =   ',STRING(ALOG10(col+zero_err)-ALOG10(col),'(F9.3)')
        print( '--------------------------------------------')
        ncol = log10(col)
        ncole1 = log10(col + tot_err) - log10(col)
        ncole2 = -log10(col - tot_err) + log10(col)
        ncolez = log10(col + zero_err) - log10(col)
    if flag_sat == 1:    
        print( '--------------------------------------------')
        print( 'log N   > ', str(log10(col), '(F9.3)'))
        print( '(+1 sig)     =   ', str(log10(col + tot_err) - log10(col), '(F9.3)'))
        print( '(-1 sig)     =   ', str(-log10(col - tot_err) + log10(col), '(F9.3)'))
        print( '--------------------------------------------')
        ncol = log10(col)
        ncole1 = log10(col + tot_err) - log10(col)
        ncole2 = -log10(col - tot_err) + log10(col)
        ncolez = 1
    #------------------------------------------------------------------------------
    #
    xpos1 = v1
    xpos2 = v2
    XLIMIT(x, xpos1, xpos2, x1, x2)
    xwork = x[x1:(x2)+1]
    ywork = y[x1:(x2)+1]
    eywork = ey[x1:(x2)+1]
    yconwork = ycon[x1:(x2)+1]
    ycon_sig_work = ycon_sig[x1:(x2)+1]
    
    #
    #Ask about what kind of statistics should be used (poisson or fixed pattern)
    #to describe the noise characteristics of the data.
    #
    
    if kind == 'f':    
        y_sig_work = sigma + zeros([ywork.size], "float32")
    else:    
        if kind == 'u':    
            y_sig_work = eywork
        else:    
            if kind == 'p':    
                y_sig_work = sigma * SQRT(absolute(ywork) / yconwork)
    
    
    #Calculate line statistics.
    #
    iSTAT2(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, m1, m1err, m2, m2err, m3, m3err, m4, m4err, dv90, v90a, v90b)
    
    # show the v90 interval
    
    OPLOT(array([[v90a, v90a]]), array([[0, 0]]), color=40, psym=7, symsize=1.5, thick=2)
    OPLOT(array([[v90b, v90b]]), array([[0, 0]]), color=40, psym=7, symsize=1.5, thick=2)
    #
    #Print information dump.
    #
    
    va = m1
    vaerr = m1err
    ba = m2 * sqrt(2)
    baerr = m2err * sqrt(2)
    
    print( '----------------------------------------------')
    print( "$('<v>       = ',f8.3,'  +/- ',f7.3)", m1, m1err)
    print( "$('<b>       = ',f8.3,'  +/- ',f7.3)", m2 * sqrt(2), m2err * sqrt(2))
    print( "$('dv90       = ',f8.3,'  +/- ',f7.3)", dv90, m1err * sqrt(2))
    print( "$('Skew      = ',f8.3,'  +/- ',f7.3)", m3, m3err)
    print( '----------------------------------------------')
    #------------------------------------------------------------------------------
    
    #------------------------------------------------------------------------------
    #
    #Calculate equivalent width and associated errors (cont ,stat, tot).
    # change first velocity to wavelength
    #
    xpos1 = v1
    xpos2 = v2
    xpos1 = iAXIS(xpos1, wavc, 1)
    xpos2 = iAXIS(xpos2, wavc, 1)
    x = iAXIS(x, wavc, 1)
    XLIMIT(x, xpos1, xpos2, x1, x2)
    xwork = array([[xpos1, x[x1:(x2)+1], xpos2]])
    ywork = array([[INTERPOL(y, x, xpos1), y[x1:(x2)+1], INTERPOL(y, x, xpos2)]])
    eywork = array([[INTERPOL(ey, x, xpos1), ey[x1:(x2)+1], INTERPOL(ey, x, xpos2)]])
    yconwork = array([[INTERPOL(ycon, x, xpos1), ycon[x1:(x2)+1], INTERPOL(ycon, x, xpos2)]])
    ycon_sig_work = array([[INTERPOL(ycon_sig, x, xpos1), ycon_sig[x1:(x2)+1]]])
    ycon_sig_work = array([[ycon_sig_work, INTERPOL(ycon_sig, x, xpos2)]])
    
    if kind == 'f':    
        y_sig_work = sigma + zeros([ywork.size], "float32")
    else:    
        if kind == 'u':    
            y_sig_work = eywork
        else:    
            if kind == 'p':    
                y_sig_work = sigma * SQRT(absolute(ywork) / yconwork)
    
    iEQW(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, ew, y_err, ycon_err, zero_err)
    toterr = SQRT(y_err ** 2 + ycon_err ** 2)
    #
    #Calculate linear column density and error.
    #
    col = (ew / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval
    colerr = col - ((ew - toterr) / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval
    col3sig = log10(3 * (toterr / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval)
    col2sig = log10(2 * (toterr / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval)
    
    #
    #Print information dump.  Equivalent widths in mA.
    #
    w = ew * 1000.0
    w_es = y_err * 1000.0
    w_ec = ycon_err * 1000.0
    w_et = toterr * 1000.0
    w_ez = zero_err * 1000.0
    print( '--------------------------------------------')
    print( 'EW           = ', str(ew * 1000.0, '(F9.2)'))
    print( 'Stat Error   = ', str(y_err * 1000.0, '(F9.2)'))
    print( 'Cont Error   = ', str(ycon_err * 1000.0, '(F9.2)'))
    print( 'Tot Error    = ', str(toterr * 1000.0, '(F9.2)'))
    #     PRINT,'2% Zero Err  = ',STRING(zero_err * 1000.0,'(F9.2)')
    print( 'Linear COG N = ', str(log10(col), '(F9.4)'))
    print( '3sigma EW    < ', str(toterr * 3000.0, '(F9.2)'))
    print( '3sigma N     < ', str(col3sig, '(F9.2)'))
    print( '2sigma EW    < ', str(toterr * 2000.0, '(F9.2)'))
    print( '2sigma N     < ', str(col2sig, '(F9.2)'))
    print( '--------------------------------------------')
    #return to velocity
    x = iAXIS(x, wavc, -1)
    print( 'iCOL_PREP::  Press ENTER to continue....')
    pause()
    _sys_ytitle = 'Flux'
    # RETURN
    # ESCAPE:
    #         PRINT,'iCOL_PREP:: '+!err_string
    #         RETURN  &  END
    
    return _ret()


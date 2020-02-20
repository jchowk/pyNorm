from numarray import *
from numpy import *
import numpy as np

#+
#                               iCOL_BATCH.PRO
#                                 Version 1.0
#
#Program Description:
#       This procedure calculates the apparent column density, equivalent width,
#       and line statistics.
#
#Restrictions:
#       None
#
#Screen Output:
#       Text  &  Graphics
#
#Use:
#       IMCOL_PREP,x,y,ex,ycon,sigma,ycon_sig,wavc,fval
#
#On Input:
#               x       :== velocity array
#               y       :== intensity array
#               ey      :== intensity error array
#               ycon    :== continuum array
#               sigma   :== RMS sigma of continuum fit
#               ycon_sig:== continuum fit error array
#               wavc    :== laboratory wavelength of line
#               fval    :== laboratory f-value of line
#
#On Ouptut:
#               storing inputs and outputs with possibility to save them at
#               the end the session
#
#Common Blocks / Structures:
#       SHARE: inputs and outputs.
#       FLAGSAT
#
#Latest Update Comments:
# 04/24/15  NL    - Version 1.0
#
#External Routines Called:
#       iCOL            - to calculate column density and errors
#       iEQW            - to calculate equivalent width and errors
#       iSTAT2          - to calculate line statistics and errors
#       XLIMIT          - to determine elements of integration
#------------------------------------------------------------------------------
def icol_batch(x, y, ey, ycon, sigma, ycon_sig, wavc, fval, vbatch1, vbatch2,
        ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez,
        va, vaerr, ba, baerr, m3, m3err, col3sig, col2sig, dv90, v90a, v90b, eyflag, kind):

   n_params = 32

   def _ret():  return (x, y, ey, ycon, sigma, ycon_sig, wavc, fval,
                vbatch1, vbatch2, ncol, ncole1, ncole2, ncolez,
                w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err,
                col3sig, col2sig, dv90, v90a, v90b, eyflag, kind)

   global flag_sat
   v90a = 0
   v90b = 0

   ycol = zeros([array(y, copy=0).nelements()], Float32)
   print 'iCOL_BATCH::  Wavelength = ', string(wavc, '(f8.3)')
   print 'iCOL_BATCH::  f-value = ', string(fval, '(f7.5)')
   #
   #
   #Get limits of integration, needs to be velocities
   #
   xpos1 = maximum(vbatch1, array(x, copy=0).min())  ;  _sys_c = 0
   xpos2 = minimum(vbatch2, max(x))  ;  _sys_c = 0


   if xpos1 > xpos2:
      print 'iCOL_PREP_BATCH::  ' + 'Limits will be reversed for integration'
      xtmp = xpos1 ; xpos1 = xpos2 ; xpos2 = xtmp
   #
   #Compute which range of elements should be included in integration and fill
   #in line.  XLIMIT will find points in between limits, but be sure to set
   #endpoints of integration equal to the limits.
   #
   v1 = xpos1
   v2 = xpos2

   x1, x2 = xlimit(x, xpos1, xpos2)
   xwork = concatenate([xpos1, x[x1:(x2)+1], xpos2])
   ywork = concatenate([interpol(y, x, xpos1), y[x1:(x2)+1], interpol(y, x, xpos2)])
   eywork = concatenate([interpol(ey, x, xpos1), ey[x1:(x2)+1], interpol(ey, x, xpos2)])
   yconwork = concatenate([interpol(ycon, x, xpos1), ycon[x1:(x2)+1], interpol(ycon, x, xpos2)])
   ycolwork = concatenate([interpol(ycol, x, xpos1), ycol[x1:(x2)+1], interpol(ycol, x, xpos2)])
   ycon_sig_work = concatenate([interpol(ycon_sig, x, xpos1), ycon_sig[x1:(x2)+1]])
   ycon_sig_work = concatenate([ycon_sig_work, interpol(ycon_sig, x, xpos2)])
   # what error to use....
   if eyflag != 0:
      print 'iCOL_PREP_BATCH::  Error Vector provided'
      y_sig_work = eywork
      kind = 'u'
   if eyflag == 0:
      if kind != _sys_null:
         if kind == 'f':
            print 'iCOL_PREP_BATCH::  Fixed pattern noise assumed'
            y_sig_work = sigma + zeros([array(ywork, copy=0).nelements()], Float32)
         if kind == 'p':
            print 'iCOL_PREP_BATCH::  Poisson noise assumed'
            y_sig_work = sigma * sqrt(abs(ywork) / yconwork)
      if kind == _sys_null:
         print 'iCOL_PREP_BATCH:: The errors are not defined'
         print 'iCOL_PREP_BATCH:: FP will be assumed'
         y_sig_work = sigma + zeros([array(ywork, copy=0).nelements()], Float32)
         kind = 'f'
   # save noise variable
   set_noise = kind

   #
   #Calculate the column density and error by calling iCOL
   icol(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, wavc, fval, col, y_err, ycon_err, zero_err)
   tot_err = sqrt(y_err ** 2 + ycon_err ** 2)
   if col <= 0.0:
      col = 1.0

   #
   #Print information dump.
   #
   if flag_sat == 0:
      print '--------------------------------------------'
      print 'log N (best) =   ', string(alog10(col), '(F9.3)')
      print '(+1 sig)     =   ', string(alog10(col + tot_err) - alog10(col), '(F9.3)')
      print '(-1 sig)     =   ', string(-alog10(col - tot_err) + alog10(col), '(F9.3)')
      #	PRINT,'2% zero err  =   ',STRING(ALOG10(col+zero_err)-ALOG10(col),'(F9.3)')
      print '--------------------------------------------'
      ncol = alog10(col)
      ncole1 = alog10(col + tot_err) - alog10(col)
      ncole2 = -alog10(col - tot_err) + alog10(col)
      ncolez = alog10(col + zero_err) - alog10(col)
   if flag_sat == 1:
      print '--------------------------------------------'
      print 'log N   > ', string(alog10(col), '(F9.3)')
      print '--------------------------------------------'
      ncol = alog10(col)
      ncole1 = 1
      ncole2 = 1
      ncolez = 1
   #------------------------------------------------------------------------------
   #
   xpos1 = v1
   xpos2 = v2
   x1, x2 = xlimit(x, xpos1, xpos2)
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
      y_sig_work = sigma + zeros([array(ywork, copy=0).nelements()], Float32)
   else:
      if kind == 'u':
         y_sig_work = eywork
      else:
         if kind == 'p':
            y_sig_work = sigma * sqrt(abs(ywork) / yconwork)



   #Calculate line statistics.
   #
   istat2(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, m1, m1err, m2, m2err, m3, m3err, m4, m4err, dv90, v90a, v90b)
   #
   #Print information dump.
   #

   va = m1
   vaerr = m1err
   ba = m2 * sqrt(2)
   baerr = m2err * sqrt(2)

   print '----------------------------------------------'
   print "$('<v>       = ',f8.3,'  +/- ',f7.3)", m1, m1err
   print "$('<b>       = ',f8.3,'  +/- ',f7.3)", m2 * sqrt(2), m2err * sqrt(2)
   print "$('dv90       = ',f8.3,'  +/- ',f7.3)", dv90, m1err * sqrt(2)
   print "$('Skew      = ',f8.3,'  +/- ',f7.3)", m3, m3err
   print '----------------------------------------------'
   #------------------------------------------------------------------------------

   #------------------------------------------------------------------------------
   #
   #Calculate equivalent width and associated errors (cont ,stat, tot).
   # change first velocity to wavelength
   #
   xpos1 = v1
   xpos2 = v2
   xpos1 = iaxis(xpos1, wavc, 1)
   xpos2 = iaxis(xpos2, wavc, 1)
   x = iaxis(x, wavc, 1)
   x1, x2 = xlimit(x, xpos1, xpos2)
   xwork = concatenate([xpos1, x[x1:(x2)+1], xpos2])
   ywork = concatenate([interpol(y, x, xpos1), y[x1:(x2)+1], interpol(y, x, xpos2)])
   eywork = concatenate([interpol(ey, x, xpos1), ey[x1:(x2)+1], interpol(ey, x, xpos2)])
   yconwork = concatenate([interpol(ycon, x, xpos1), ycon[x1:(x2)+1], interpol(ycon, x, xpos2)])
   ycon_sig_work = concatenate([interpol(ycon_sig, x, xpos1), ycon_sig[x1:(x2)+1]])
   ycon_sig_work = concatenate([ycon_sig_work, interpol(ycon_sig, x, xpos2)])

   #
   #Ask about what kind of statistics should be used (poisson or fixed pattern)
   #to describe the noise characteristics of the data.
   #

   if kind == 'f':
      y_sig_work = sigma + zeros([array(ywork, copy=0).nelements()], Float32)
   else:
      if kind == 'u':
         y_sig_work = eywork
      else:
         if kind == 'p':
            y_sig_work = sigma * sqrt(abs(ywork) / yconwork)

   #
   ieqw(xwork, ywork, yconwork, y_sig_work, ycon_sig_work, ew, y_err, ycon_err, zero_err)
   toterr = sqrt(y_err ** 2 + ycon_err ** 2)
   #
   #Calculate linear column density and error.
   #
   col = (ew / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval
   colerr = col - ((ew - toterr) / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval
   col3sig = alog10(3 * (toterr / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval)
   col2sig = alog10(2 * (toterr / wavc) / 8.85e-13 / (wavc * 1.e-8) / fval)

   #
   #Print information dump.  Equivalent widths in mA.
   #
   w = ew * 1000.0
   w_es = y_err * 1000.0
   w_ec = ycon_err * 1000.0
   w_et = toterr * 1000.0
   w_ez = zero_err * 1000.0
   print '--------------------------------------------'
   print 'EW           = ', string(ew * 1000.0, '(F9.2)')
   print 'Stat Error   = ', string(y_err * 1000.0, '(F9.2)')
   print 'Cont Error   = ', string(ycon_err * 1000.0, '(F9.2)')
   print 'Tot Error    = ', string(toterr * 1000.0, '(F9.2)')
   #     PRINT,'2% Zero Err  = ',STRING(zero_err * 1000.0,'(F9.2)')
   print 'Linear COG N = ', string(alog10(col), '(F9.4)')
   print '3sigma EW    < ', string(toterr * 3000.0, '(F9.2)')
   print '3sigma N     < ', string(col3sig, '(F9.2)')
   print '2sigma N     < ', string(col2sig, '(F9.2)')
   print '--------------------------------------------'
   #return to velocity
   x = iaxis(x, wavc, -1)


   return _ret()
   # escape:
   print 'iCOL_BATCH:: ' + _sys_err_string
   return _ret()  ;

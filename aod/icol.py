from numarray import *

#+
#				ICOL.PRO
#				Version 1.0
#
#Created: Unknown
#Last Revised: 05/02/99
#
#Program Description:
#       This procedure calculates column densities and errors for
#       spectral lines given a fitted continuum with errors.  Column
#	densities and errors are calculated using the integrated apparent
#	column density technique.
#
#Restrictions:
#       Error points in the ordinate array are assumed to be uncorrelated.
#       Error points in the continuum array are assumed to be correlated.
#
#;
#Use:
# 	IMCOL,x,y,ycon,y_sig,ycon_sig,wavc,fval,col,y_err,ycon_err,zero_err
#
#On Input:
#               x       :== abcissa array
#               y       :== ordinate array
#               ycon    :== fitted continuum array
#               y_sig   :== error array for ordinate
#               ycon_sig:== error array for fitted continuum
#		wavc	:== wavelength of line
#		fval	:== fvalue of line
#
#On Ouptut:
#               col     :== integrated column density
#               y_err   :== column density error due to y_sig
#               ycon_err:== column density error due to ycon_sig
#               zero_err:== column density error due to 2% background err
#
#Common Blocks / Structures:
#       None
#
#Latest Update Comments:
#       04/12/13  KRS   - Version 1.0
#
#External Routines Called:
#       None
#------------------------------------------------------------------------------
def icol(x, y, ycon, y_sig, ycon_sig, wavc, fval, col, y_err, ycon_err, zero_err, root):
   n_params = 12
   def _ret():  return (x, y, ycon, y_sig, ycon_sig, wavc, fval, col, y_err, ycon_err, zero_err, root)
   
   global flag_sat
   
   # IF N_PARAMS() EQ 0 THEN BEGIN MAN,'imcol' & RETURN & ENDIF
   #
   #Calculate dx.
   #
   nx = array(x, copy=0).nelements()
   dx = zeros([nx], Float32)
   for j in arange(1, (nx - 2)+(1)):
      dx(j) = (x(j + 1) - x(j - 1)) / 2.0
   dx(0) = x(1) - x(0)               #Updated over dx(0)=dx(1)
   dx(nx - 1) = x(nx - 1) - x(nx - 2)      #Updated over dx(nx-1) = dx(nx-2)
   #
   #Calculate integrated apparent column density.  Units are atoms/cm2/(km/s).
   #
   tau = total(alog(ycon / y) * dx)
   col = tau / (wavc * fval * 2.654e-15)
   flag_sat = 0
   test = where(ravel(y <= 0))[0]
   yy1 = y
   if ct != 0:   
      # artificially remove 0 and negative flux
      yy1[test] = abs(yy1[test])
      yy2 = yy1[test]
      test1 = where(ravel(yy2 == 0))[0]
      if ct1 != 0:   
         yy1[test1] = 2 * abs(y_sig[test1])
      flag_sat = 1
      tau = array(total(alog(ycon / yy1) * dx), copy=0).astype(Float64)
      col = array(tau / (wavc * fval * 2.654e-15), copy=0).astype(Float64)
   #
   #Calculate error on integrated apparent optical depth due to intensity errors.
   #
   y_err = sqrt(total(y_sig ** 2 * (-1. / y) ** 2 * dx ** 2)) / (wavc * fval * 2.654e-15)
   #
   #Calculate error on integrated apparent optical depth due to continuum errors.
   #Add the errors by straight addition since errors are correlated => cannot add
   #in quadrature.
   ycon_err = total(ycon_sig * (1. / ycon) * dx) / (wavc * fval * 2.654e-15)
   #
   #Calculate error on integrated apparent optical depth due to 2% zero level
   #shift
   #
   z_eps = 0.02
   zero_err = total(alog(1 + z_eps * alog(ycon / y)) * dx) / (wavc * fval * 2.654e-15)
   #
   #Empirical estimate of error due to 2% zero level shift.
   #
   z_eps = 0.02
   yc1 = ycon * (1 - z_eps)
   y1 = y - ycon * z_eps
   tau1 = total(alog(yc1 / y1) * dx)
   col1 = tau1 / (wavc * fval * 2.654e-15)
   zero_err = abs(col1 - col)
   
   return _ret()  ;

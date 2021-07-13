"""
Kenneth Sembach
				iSTAT2.PRO
				Version 1.0

Program Description:
	This procedure calculates line statistics weighted by the data
	values (mean, width, skewness, etc.).  It is called by IMSTAT_PREP.

Screen Output:
	Text

Use:
	iSTAT2,x,y,ycon,y_sig,ycon_sig,m1,m1err,m2,m2err,m3,m3err,m4,m4err

On Input:
		x       :== wavelength array
		y       :== intensity array
		ycon    :== continuum array
		y_sig   :== statistical errors in y
		ycon_sig:== continuum fit error array

On Output:
		m1	:== average velocity
		m1err	:== error on m1
		m2	:== width
		m2err	:== error on m2
		m3	:== skewness
		m3err	:== error on m3
		m4	:== extent (same as m2 with m1=0)
		m3err	:== error on m4

Common Blocks / Structures:
	None

Latest Update Comments:
	04/13/13  NL	- Version 1.0

External Routines called:
	None
"""

from __future__ import print_function

from numpy import *

def iSTAT2(x, y, ycon, y_sig, ycon_sig, m1, m1err, m2, m2err, m3, m3err, m4, m4err, dv90, v90a, v90b):
    """
    
    Calculate dx.
    """

    n_params = 16
    def _ret():  return (x, y, ycon, y_sig, ycon_sig, m1, m1err, m2, m2err, m3, m3err, m4, m4err, dv90, v90a, v90b)
    
    nx = x.size
    dx = zeros([nx], "float32")
    for j in arange(1, (nx - 2)+(1)):
        dx(j) = (x(j + 1) - x(j - 1)) / 2.0
    dx(0) = dx(1)  ;  dx(nx - 1) = dx(nx - 2)
    #
    #Calculate zeroth moment (equivalent width).
    #
    tau = array(log(absolute(ycon / y)), copy=0).astype(float64)
    tautot = array(TOTAL(tau * dx), copy=0).astype(float64)
    #
    #Calculate first moment (average velocity).
    #
    a = array(TOTAL(tau * x * dx), copy=0).astype(float64)
    m1 = array(a / tautot, copy=0).astype(float64)
    #
    #Calculate second moment (width).
    #
    b = array(TOTAL(tau * (x - m1) ** 2 * dx), copy=0).astype(float64)
    m2 = array(SQRT(b / tautot), copy=0).astype(float64)
    #
    #Calculate third moment (skewness).
    #
    f = array(TOTAL(tau * ((x - m1) / m2) ** 3 * dx), copy=0).astype(float64)
    m3 = array(f / tautot, copy=0).astype(float64)
    #
    #Calculate extent (same as m2, except that m1 is assumed to be 0).
    #
    b4 = array(TOTAL(tau * (x - 0.0) ** 2 * dx), copy=0).astype(float64)
    m4 = array(SQRT(b4 / tautot), copy=0).astype(float64)
    
    # Estimate the velocities at 5% and 95% to estimate dv90 (90% optical-depth width)
    
    tautot5 = 0.05 * tautot
    tautot95 = 0.95 * tautot
    
    
    k1 = 0
    k2 = 0
    tautot0 = tau(0) * (x(2) - x(0)) / 2.
    for i in arange(1, (x.size - 2)+(1)):
        tautot0 = tautot0 + tau(i) * (x(i + 1) - x(i - 1)) / 2.
        if bitwise_and(tautot0 >= tautot5, k1 == 0):    
            v90a = x(i)
            k1 = 1
        if bitwise_and(tautot0 >= tautot95, k2 == 0):    
            v90b = x(i)
            k2 = 1
    dv90 = v90b - v90a
    
    
    #
    #Calculate error on m1.
    #
    dadi = -1. / y * x * dx
    dadc = 1 / ycon * x * dx
    dwdi = -1. / y * dx
    dwdc = 1 / ycon * dx
    
    dm1di = (tautot * dadi - a * dwdi) / tautot ** 2
    dm1dc = (tautot * dadc - a * dwdc) / tautot ** 2
    
    q1 = SQRT(TOTAL(y_sig ** 2 * dm1di ** 2))
    q2 = TOTAL(SQRT(ycon_sig ** 2 * dm1dc ** 2))
    
    m1err = SQRT(q1 ** 2 + q2 ** 2)
    #
    #Calculate error on m2.
    #
    dbdi = -1. / y * (x - m1) ** 2 * dx
    dbdc = 1. / ycon * (x - m1) ** 2 * dx
    dbdm1 = -2. * tau * (x - m1) * dx
    
    dm2di = (tautot * dbdi - b * dwdi) / tautot ** 2
    dm2dc = (tautot * dbdc - b * dwdc) / tautot ** 2
    dm2dm1 = dbdm1 / tautot
    
    q1 = SQRT(TOTAL(y_sig ** 2 * dm2di ** 2))
    q2 = TOTAL(SQRT(ycon_sig ** 2 * dm2dc ** 2))
    q3 = SQRT(TOTAL(m1err ** 2 * dm2dm1 ** 2))
    
    m2err = SQRT(q1 ** 2 + q2 ** 2 + q3 ** 2)
    m2err = m2err / (2. * m2)
    #
    #Calculate error on m3.
    #
    dfdi = -1. / y * ((x - m1) / m2) ** 3 * dx
    dfdc = 1. / ycon * ((x - m1) / m2) ** 3 * dx
    dfdm1 = tau * 3. * ((x - m1) / m2) ** 2 * (-1. / m2) * dx
    dfdm2 = tau * 3. * ((x - m1) / m2) ** 2 * (m1 - x) / m2 ** 2 * dx
    
    dm3di = (tautot * dfdi - f * dwdi) / tautot ** 2
    dm3dc = (tautot * dfdc - f * dwdc) / tautot ** 2
    dm3dm1 = dfdm1 / tautot
    dm3dm2 = dfdm2 / tautot
    
    q1 = SQRT(TOTAL(y_sig ** 2 * dm3di ** 2))
    q2 = TOTAL(SQRT(ycon_sig ** 2 * dm3dc ** 2))
    q3 = SQRT(TOTAL(m1err ** 2 * dm3dm1 ** 2))
    q4 = SQRT(TOTAL(m2err ** 2 * dm3dm2 ** 2))
    
    m3err = SQRT(q1 ** 2 + q2 ** 2 + q3 ** 2 + q4 ** 2)
    #
    #Calculate error on m4.
    #
    dbdi = -1. / y * (x - 0.0) ** 2 * dx
    dbdc = 1. / ycon * (x - 0.0) ** 2 * dx
    
    dm4di = (tautot * dbdi - b4 * dwdi) / tautot ** 2
    dm4dc = (tautot * dbdc - b4 * dwdc) / tautot ** 2
    
    q1 = SQRT(TOTAL(y_sig ** 2 * dm4di ** 2))
    q2 = TOTAL(SQRT(ycon_sig ** 2 * dm4dc ** 2))
    
    m4err = SQRT(q1 ** 2 + q2 ** 2)
    m4err = m4err / (2. * m4)
    #
    #Return to caller.
    #
    return _ret()



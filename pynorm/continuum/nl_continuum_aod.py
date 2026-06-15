"""
   Version 3.0

 PURPOSE:
 This program allows you to automatically model a continuum near a given transition using a Legendre polynome and estimate the AOD parameters.

 CALLING SEQUENCE:

   nl_continuum_leg, '*i.save'

 OPTIONAL KEYWORDS:

 EXAMPLE:
   IDL> cbw_continuum_aod, '*i.save'

 OUTPUTS

 Continuum, AOD results

 External program:
   iCOL_BATCH
 OPTIONAL OUTPUTS:
 ;

 REVISION HISTORY:
   v1.0 - 25-April-2015 Written by NL;
   v1.1 - 2016-10-21 Updated by CBW to fix a typo. When saving, v1.0 was
                     saving ..."efnorm1,efnorm1"... but we want
                     to save ..."efnorm1,efnorm2"...
   v2.0 - 2016-10-21 Updated by CBW to double-fit the continuum.
                     This is a huge update. The new approach:
                     (1) Separates out the continuum fitting portion (into a
                     separate function at the top). I did NOT change anything
                     here, just separated it out because we use it twice.
                     (2) We now double-fit the continuum, meaning that, once
                     the continuum is fit the first time, we bin the data, and
                     then, if the mean if each bin lies a certain "sigma" BELOW
                     (we're only looking to reject ABSORPTION lines the second time)
                     the first continuum fit, we reject it. This method
                     much more effectively rejects absorption lines.

                     The variables that can be changed in this section of the
                     code are the binning size (pix2) and the "sigma" (nsig2)
                     below which the data within the bin are rejected.
   v3.0 - 2016-11-21 Updated by NL to improve the continuum placement by
   	    	      clipping the flux in different sections to estimate the continuum
   	    	      region (new keyword vclip1)
   	    	      Added a pdf file flag to create PDF files from the EPS file; however
   	    	      this does slow down things.
-
"""

from __future__ import print_function

from numpy import *

def cbw_contfit(x, y, continuum, continuum_err, mask, minord, maxord, sn, ftflag):

    n_params = 9
    def _ret():  return (x, y, continuum, continuum_err, mask, minord, maxord, sn, ftflag)
    
    g = where(ravel(mask != 0))[0]
    xarray = x[g]
    yarray = y[g]
    
    
    minord = array(minord, copy=0).astype(int32)
    maxord = choose(array(maxord, copy=0).astype(int32) < minord, (array(maxord, copy=0).astype(int32), minord))
    maxord = choose(maxord > 11, (maxord, 11))
    minord = choose(minord > 11, (minord, 11))
    #
    
    #Scale xarray and x into appropriate range so that Legendre polynomial can be
    #fit over the range from -1 to + 1.
    #
    maxx = array(absolute(xarray), copy=0).max()  ;  _sys_c = 0
    xarray = xarray / maxx
    x = x / maxx
    #
    #Get Legendre polynomial fit and form continuum over all x.
    #
    LEGFIT(xarray, yarray, minord, maxord, yfit, a, eps, chi2)
    continuum = LEGPOLY(x, a)
    #
    #Get sigma from Legendre polynomial fit.  It is just square root of chi2 here.
    #
    sigma = SQRT(chi2)
    #
    #Calculate mean signal to noise
    #
    sn = 0.0
    if sigma != 0.0:    
        sn = mean(yfit) / sigma
    #
    
    #Calculate the error bars on the continuum due to errors in the coefficients
    #of the Legendre polynomial fit.
    #
    LEGERR(x, y, a, eps, chi2, continuum_err)    #chi2 = variance (=sigma^2) here
    #
    #Convert x and xarray back into the correct units.
    #
    x = x * maxx
    xarray = xarray * maxx
    #
    #Bound continuum so that plotting errors and arithmetic errors aren't important.
    #
    maxy = array(y, copy=0).max() ; miny = array(y, copy=0).min()
    continuum = choose(choose(continuum < (-absolute(miny) * 2), (continuum, (-absolute(miny) * 2))) > (maxy * 2), (choose(continuum < (-absolute(miny) * 2), (continuum, (-absolute(miny) * 2))), (maxy * 2)))
    #
    #Print the fitting statistics to upper corner of screen.
    #
    print( '----------------------------------------')
    print( 'Order used = ', maxord)
    print( '# points   = ', str(yarray.size, '(I8)'))
    print( 'RMS sigma  = ', str(sigma, '(E9.3)'))
    print( 'Mean S/N   = ', str(sn, '(F8.2)'))
    print( '----------------------------------------')
    ftflag = 1
    
    
    
    
    return _ret()



################################################################################################
################################################################################################
################################################################################################


def nl_continuum_aod(filein, maxord=None, minord=None, nsig1=None, nsig2=None, vl1=None, vl2=None, vaod1=None, vaod2=None, contplot=None, aodplot=None, vclip1=None, blemish=None, vshift=None, flgaod=None, pix1=None, pix2=None, fmin=None, nocont=None, noerror=None, kind=None, pdfflg=None, changename=None):

    n_params = 1
    _opt = (maxord, minord, nsig1, nsig2, vl1, vl2, vaod1, vaod2, contplot, aodplot, vclip1, blemish, vshift, flgaod, pix1, pix2, fmin, nocont, noerror, kind, pdfflg, changename)
    def _ret():
        _optrv = zip(_opt, [maxord, minord, nsig1, nsig2, vl1, vl2, vaod1, vaod2, contplot, aodplot, vclip1, blemish, vshift, flgaod, pix1, pix2, fmin, nocont, noerror, kind, pdfflg, changename])
        _rv = [filein]
        _rv += [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    if bitwise_not((filein is not None)):    
        print( 'No input file...')
        return _ret()
    if bitwise_not((minord is not None)):    
        minord = 0
    if bitwise_not((maxord is not None)):    
        maxord = 10
    if bitwise_not((NSIG1 is not None)):    
        nsig1 = 3
    if bitwise_not((NSIG2 is not None)):    
        nsig2 = 0.75
    if bitwise_not((vl1 is not None)):    
        vl1 = -1000
    if bitwise_not((vl2 is not None)):    
        vl2 = +1000
    if bitwise_not((vaod1 is not None)):    
        vaod1 = -100
    if bitwise_not((vaod2 is not None)):    
        vaod2 = +100
    if bitwise_not((pix1 is not None)):    
        pix1 = 4
    if bitwise_not((pix2 is not None)):    
        pix2 = 8
    if bitwise_not((fmin is not None)):    
        fmin = 0.5e-15
    if bitwise_not((vclip1 is not None)):    
        vclip1 = 5.
    if bitwise_not((vshift is not None)):    
        vshift = 0
    if (changename is not None):    
        str = strjoin(strsplit(filein, 'i.save', regex=True, extract=True))
    
    # kind will define error if /noerror is selected. By default FP. choices are kind = 'P' (poisson noise) or kind = 'F' (fixed pattern noise)
    
    restore(filein)
    
    ns = f.size
    
    v = v + vshift
    
    # define region where the continuum will be fitted
    g = where(ravel(bitwise_and(v >= vl1, v <= vl2)))[0]
    fi1 = f[g]
    vi1 = v[g]
    ei1 = ef[g]
    #===========================================================================
    # remove blemish
    if (blemish is not None):    
        for j in arange(1, (fi1.size - 2)+(1)):
            if (fi1(j) <= 0):    
                if (bitwise_and(absolute(fi1(j + 1)) > 0, absolute(fi1(j - 1)) > 0)):    
                    fi1(j) = (fi1(j - 1) + fi1(j + 1)) / 2.
    #===========================================================================
    
    #g= where(fi1 gt fmin,test)
    #if test ne 0 then fm = median(fi1[where(fi1 gt fmin)])
    #if test eq 0 then fm = median(fi1)
    
    ns1 = fi1.size
    mask_cont = ((1)*ones(ns1, dtype="int32"))
    
    # cut the spectrum in sub-sections to check regions to be masked. This should help when
    # flux changes by a large amount over the considered velocity interval.
    
    diffv = (vl2 - vl1) / vclip1
    diffv1 = vl1
    diffv2 = vl1 + diffv
    j = 0
    for k in arange(0, (array(vclip1, copy=0).astype(int32) - 1)+(1)):
        tt = where(ravel(bitwise_and(vi1 >= diffv1, vi1 < diffv2)))[0]
        ftemp = fi1[tt]
        etemp = ei1[tt]
        g = where(ravel(ftemp > fmin))[0]
        if test != 0:    
            fm = median(ftemp[where(ravel(ftemp > fmin))[0]])
        if test == 0:    
            fm = median(ftemp)
        sigp = fm + nsig1 * etemp
        sigm = fm - nsig1 * etemp
        #    print, diffv1, diffv2, fm,pix1
        for i in arange(pix1, (ftemp.size - 1)+(pix1), pix1):
            j = j + pix1
            if (mean(ftemp[i - pix1:(i)+1]) < mean(sigm[i - pix1:(i)+1])):    
                mask_cont[j - pix1:(j)+1] = 0
            if (mean(ftemp[i - pix1:(i)+1]) > mean(sigp[i - pix1:(i)+1])):    
                mask_cont[j - pix1:(j)+1] = 0
            #     print, i,j,ns1, mean(ftemp[i-pix1:i]), mean(sigp[i-pix1:i]), mean(sigm[i-pix1:i])
        diffv1 = diffv2
        diffv2 = diffv1 + diffv
    
    
    #for i = 0,ns1 -1 do print, vi1[i],fi1[i],mask_cont[i]
    
    
    
    
    if bitwise_not((nocont is not None)):    
        
        
        ########################
        ## CBW 2016-10-21 begin
        ########################
        
        ##continuum fit the first time
        vi1, fi1, ycon, ycon_sig, mask_cont, minord, maxord, sn, ftflag = cbw_contfit(vi1, fi1, ycon, ycon_sig, mask_cont, minord, maxord, sn, ftflag)
        
        ##copy this, just in case. We don't use it elsewhere at the moment.
        mask_cont_firstfit = mask_cont
        
        
        ############
        ##Now reassess the errors as "sigma away from the continuum fit"
        ##That is, recreate the mask after a first continuum fit
        ##But, only look for absorption lines the second time,
        ##  not positive spikes!!
        ############
        
        ##comment out; don't overwrite the mask!
        # mask_cont = make_array(ns1,/integer,value=1)
        
        
        ######
        ## first method: pixel-by-pixel
        ######
        
        ##comment out; let's do this using a binned method like Nicolas did above
        # bcontrefit=where(fi1/ycon lt (1-1.5*stdev(fi1/ycon)))
        # mask_cont[bcontrefit] = 0
        
        
        ######
        ## second method: binned
        ######
        
        ##I have two ways to do this:
        ##  using mean(efnorm)
        ##  using stdev(fnorm)
        ##The first uses the error vector, which may be proper
        ##The second uses the stdev of the flux itself, which may
        ##  produce better results, depending??
        
        
        ##OK, *these are equivalent* as long as the coefficients (nsig2) are scaled based on "mean(efnorm)/stdev(fnorm)"
        ##  i.e., 1-0.75*mean(efnorm) rejects the same number of points as
        ##        1-0.593243*stdev(fnorm)
        ##  where 0.593243 = 0.75*mean(efnorm)/stdev(fnorm)
        
        
        ##first way -- this is safer and rejects more consistently, whereas
        ##  with the other way, nsig2 has to be re-calculated each time
        ##using mean(efnorm)
        # pix2 = 2*pix1
        #    for i = pix2,ns1-1,pix2 do begin
        ##if the mean of the *binned* pixels is less than
        ##  half of stdev(fnorm) below 1.0, then mask it
        #        if (mean(fi1[i-pix2:i]/ycon[i-pix2:i]) lt (1-nsig2*mean(ei1/ycon))) then mask_cont[i-pix2:i]=0
        #    endfor
        
        
        # ;;second way
        # ;;using stdev(fnorm)
        # ; pix2 = 2*pix1
        # pix2 = 8.
        # nsig2=0.50
        # stop;;cbw
        for i in arange(pix2, (ns1 - 1)+(pix2), pix2):
        #     ;;if the mean of the *binned* pixels is less than
        #     ;;  half of stdev(fnorm) below 1.0, then mask it
            if (mean(fi1[i - pix2:(i)+1] / ycon[i - pix2:(i)+1]) < (1 - nsig2 * stdev(fi1 / ycon))):    
                mask_cont[i - pix2:(i)+1] = 0
        
        
        
        
        ##continuum fit the second time, hopefully with fewer absorption lines
        vi1, fi1, ycon, ycon_sig, mask_cont, minord, maxord, sn, ftflag = cbw_contfit(vi1, fi1, ycon, ycon_sig, mask_cont, minord, maxord, sn, ftflag)
        
        
        
        ########################
        ## CBW 2016-10-21 end
        ########################
        
        
        vnorm = vi1
        fnorm = fi1 / ycon
        eform = 0 * fnorm
        eform1 = 0 * fnorm
        eform2 = 0 * fnorm
        sigma0 = 0 * fnorm
        if bitwise_not((noerror is not None)):    
            efnorm = sqrt((ycon_sig / ycon) ** 2.0 + (ei1 / ycon) ** 2.)
            efnorm1 = sqrt((ycon_sig / ycon) ** 2.0 + (ei1 / ycon) ** 2.)
            efnorm2 = (ycon_sig / ycon) + (ei1 / ycon)
            efnorm = (efnorm1 + efnorm2) / 2
            sigma0 = ei1 / ycon
        if bitwise_not((name is not None)):    
            if (object is not None):    
                targname = object
            else:    
                targname = 'file_out'
        else:    
            targname = name
        
        if (changename is not None):    
            outname = 'plot_cont_' + str + '.eps'
        if bitwise_not((changename is not None)):    
            outname = 'plot_cont_' + ion + wni + '.eps'
        
        if (contplot is not None):    
            
            SET_PLOT('PS')
            DEVICE(FILENAME=outname, FONT_SIZE=11, INCHES=True, XSIZE=8, YSIZE=5, ENCAPSULATED=True, COLOR=True, portrait=True)
            
            _sys_P.MULTI = array([[0, 1, 1]])
            #!y.margin = [4,2]
            _sys_P.CHARSIZE = 1.6
            _sys_p.psym = 10
            loadct(40, silent=True)
            
            
            yax = array(fi1 * 1e15, copy=0).max() + array(fi1 * 1e15, copy=0).max() * 0.05
            
            loadct(40, silent=True)
            plot(vi1, fi1 * 1e15, xr=array([[vl1, vl2]]), yr=array([[-0.2, yax]]), psym=10, thick=2, charthick=2, ys=True, xs=True, xtitle='!6Velocity (km s!e-1!n)', ytitle='Flux (10!e-15!n erg cm!e-2!n s!e-1!n !n!sA!r!u!9 %!6!N!e-1 !N)', position=array([[0.11, 0.14, 0.94, 0.98]]))
            plotzero()
            
            
            g = where(ravel(mask_cont == 0))[0]
            if test != 0:    
                oplot(vi1[g], fi1[g] * 1e15, psym=7, thick=4, color=250)
            
            oplot(vi1, ycon * 1e15, line=0, thick=4, color=200)
            ver(vaod1, lines=1, color=40, thick=4)
            ver(vaod2, lines=1, color=40, thick=4)
            
            xyouts(0.12, 0.16, ion + ' !7k!6' + STRTRIM(str(wavc, FORMAT='(f7.1)'), 1), al=0, charsize=2., charthick=2, normal=True, color=79)
            
            DEVICE(CLOSE=True)
            
            set_plot('x')
            cleanplot(silent=True)
            if (pdfflg is not None):    
                spawn('epstopdf ' + outname)
            
        
    else:    
        vnorm = vi1
        fnorm = fi1
        if bitwise_not((noerror is not None)):    
            efnorm = ei1
            efnorm1 = ei1
            sigma0 = ei1
        if bitwise_not((name is not None)):    
            if (object is not None):    
                targname = object
            else:    
                targname = 'file_out'
        else:    
            targname = name
        ycon = ((1.0)*ones(vnorm.size, dtype="float32"))
        ycon_sig = ((0.0)*ones(vnorm.size, dtype="float32"))
        sn = 1. / mean(efnorm)
    
    
    if (aodplot is not None):    
        
        if (changename is not None):    
            outname = 'plot_cont_aod_' + str + '.eps'
        if bitwise_not((changename is not None)):    
            outname = 'plot_cont_aod_' + ion + wni + '.eps'
        
        
        SET_PLOT('PS')
        
        DEVICE(FILENAME=outname, FONT_SIZE=11, INCHES=True, XSIZE=8, YSIZE=9, ENCAPSULATED=True, COLOR=True, portrait=True)
        
        _sys_P.MULTI = array([[0, 1, 2]])
        #!y.margin = [4,2]
        _sys_P.CHARSIZE = 1.6
        _sys_p.psym = 10
        loadct(40, silent=True)
        
        if bitwise_not((nocont is not None)):    
            yax = array(fi1 * 1e15, copy=0).max() + array(fi1 * 1e15, copy=0).max() * 0.1
            
            fp = fi1 * 1e15
            yp = ycon * 1e15
            
            plot(vi1, fp, xr=array([[vl1, vl2]]), yr=array([[-0.2, yax]]), psym=10, thick=2, charthick=2, ys=True, xs=True, xtitle='', ytitle='Flux (10!e-15!n erg cm!e-2!n s!e-1!n !n!sA!r!u!9 %!6!N!e-1 !N)', position=array([[0.13, 0.58, 0.94, 0.99]]))
            plotzero()
            
        else:    
            yax = array(fi1, copy=0).max() + array(fi1, copy=0).max() * 0.1
            fp = fi1
            yp = ycon
            #!y.tickinterval =0.5
            
            plot(vi1, fp, xr=array([[vl1, vl2]]), yr=array([[-0.2, yax]]), psym=10, thick=2, charthick=2, ys=True, xs=True, xtitle='', ytitle='Normalized Flux', position=array([[0.13, 0.58, 0.94, 0.99]]))
            plotzero()
        
        
        oplot(vi1, yp, line=0, thick=4, color=200)
        ver(vaod1, lines=1, color=40, thick=4)
        ver(vaod2, lines=1, color=40, thick=4)
        g = where(ravel(mask_cont == 0))[0]
        if test != 0:    
            oplot(vi1[g], fp[g], psym=7, thick=4, color=250)
        
        xyouts(0.15, 0.60, ion + ' !7k!6' + STRTRIM(str(wavc, FORMAT='(f7.1)'), 1), al=0, charsize=2., charthick=2, normal=True, color=79)
        
        
        g = where(ravel(bitwise_and(vi1 >= vaod1, vi1 <= vaod2)))[0]
        ycol = -log(fnorm) / (wavc * fval * 2.654e-15)
        yax = array(ycol[g] * 1e-13, copy=0).max() + array(ycol[g] * 1e-13, copy=0).max() * 0.1
        
        plot(vi1, ycol * 1e-13, xr=array([[vaod1 - 150, vaod2 + 150]]), yr=array([[-0.05, yax]]), psym=10, thick=2, charthick=2, ys=True, xs=True, xtitle='!6Velocity (km s!e-1!n)', ytitle='N!da!n(v) (10!e13!n cm!e-2!n (km s!e-1!n)!e-1!n)', position=array([[0.13, 0.13, 0.94, 0.54]]))
        plotzero()
        plotcolorfill(vi1[g], ycol[g] * 1e-13, midpoint=True, noerase=True, col=234, bottom=0)
        
        ver(vaod1, lines=1, color=40, thick=4)
        ver(vaod2, lines=1, color=40, thick=4)
        
        if (pdfflg is not None):    
            spawn('epstopdf ' + outname)
        
        
        DEVICE(CLOSE=True)
        
        set_plot('x')
        cleanplot(silent=True)
        
    
    flux = fi1
    eflux = 0 * flux
    if bitwise_not((noerror is not None)):    
        eflux = ei1
    if (noerror is not None):    
        efnorm = eflux * 0
        efnorm1 = eflux * 0
        efnorm2 = eflux * 0
    
    vel = vi1
    
    if (flgaod is not None):    
        print( 'iCOL_BATCH::  Ion = ', ion)
        
        if (changename is not None):    
            outname = str + 'i_o.save'
        if bitwise_not((changename is not None)):    
            outname = ion + wni + 'i_o.save'
        
        
        if bitwise_not((noerror is not None)):    
            eyflag = 1
            iCOL_BATCH(vi1, fi1, ei1, ycon, sigma, ycon_sig, wavc, fval, vaod1, vaod2, ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err, col3sig, col2sig, dv90, v90a, v90b, eyflag, kind)
            
            v1 = vaod1
            v2 = vaod2
            save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, name, ra, dec, gl, gb, z, ycon, ycon_sig, vnorm, fnorm, efnorm, efnorm1, efnorm2, sigma0, targname, sn, maxord, v1, v2, vaod1, vaod2, ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err, col3sig, col2sig, dv90, v90a, v90b, mask_cont, vel, flux, eflux, fil=outname)
            
        
        if (noerror is not None):    
            eyflag = 0
            iCOL_BATCH(vi1, fi1, ei1, ycon, sigma, ycon_sig, wavc, fval, vaod1, vaod2, ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err, col3sig, col2sig, dv90, v90a, v90b, eyflag, kind)
            v1 = vaod1
            v2 = vaod2
            
            save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, name, ra, dec, gl, gb, z, ycon, ycon_sig, vnorm, fnorm, efnorm, efnorm1, efnorm2, sigma0, targname, sn, maxord, v1, v2, vaod1, vaod2, ncol, ncole1, ncole2, ncolez, w, w_es, w_ec, w_et, w_ez, va, vaerr, ba, baerr, m3, m3err, col3sig, col2sig, dv90, v90a, v90b, mask_cont, vel, flux, eflux, fil=outname)
        
    else:    
        save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, name, ra, dec, gl, gb, z, ycon, ycon_sig, vnorm, fnorm, efnorm, targname, sn, maxord, mask_cont, vel, flux, eflux, fil=outname)
    
    
    
    
    return _ret()
    


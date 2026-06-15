; +
;   Version 3.0

; PURPOSE:
; This program allows you to automatically model a continuum near a given transition using a Legendre polynome and estimate the AOD parameters.
;
; CALLING SEQUENCE:
;
;   nl_continuum_leg, '*i.save'
;
; OPTIONAL KEYWORDS:
;
; EXAMPLE:
;   IDL> cbw_continuum_aod, '*i.save'
;
; OUTPUTS
;
; Continuum, AOD results
;
; External program:
;   iCOL_BATCH
; OPTIONAL OUTPUTS:
; ;
;
; REVISION HISTORY:
;   v1.0 - 25-April-2015 Written by NL;
;   v1.1 - 2016-10-21 Updated by CBW to fix a typo. When saving, v1.0 was
;                     saving ..."efnorm1,efnorm1"... but we want
;                     to save ..."efnorm1,efnorm2"...
;   v2.0 - 2016-10-21 Updated by CBW to double-fit the continuum.
;                     This is a huge update. The new approach:
;                     (1) Separates out the continuum fitting portion (into a
;                     separate function at the top). I did NOT change anything
;                     here, just separated it out because we use it twice.
;                     (2) We now double-fit the continuum, meaning that, once
;                     the continuum is fit the first time, we bin the data, and
;                     then, if the mean if each bin lies a certain "sigma" BELOW
;                     (we're only looking to reject ABSORPTION lines the second time)
;                     the first continuum fit, we reject it. This method
;                     much more effectively rejects absorption lines.
;                     
;                     The variables that can be changed in this section of the
;                     code are the binning size (pix2) and the "sigma" (nsig2)
;                     below which the data within the bin are rejected.
;   v3.0 - 2016-11-21 Updated by NL to improve the continuum placement by 
;   	    	      clipping the flux in different sections to estimate the continuum
;   	    	      region (new keyword vclip1)
;   	    	      Added a pdf file flag to create PDF files from the EPS file; however
;   	    	      this does slow down things.    	    	
;-
;------------------------------------------------------------------------------



PRO cbw_contfit,x,y,continuum,continuum_err,mask,minord,maxord,sn,ftflag

g =where(mask ne 0)
xarray = x[g]
yarray = y[g]


minord = FIX(minord)
maxord = FIX(maxord) > minord
maxord = maxord < 11
minord = minord < 11
;

;Scale xarray and x into appropriate range so that Legendre polynomial can be
;fit over the range from -1 to + 1.
;
maxx = MAX(ABS(xarray))  &  !c = 0
xarray = xarray/maxx
x = x / maxx
;
;Get Legendre polynomial fit and form continuum over all x.
;
LEGFIT,xarray,yarray,minord,maxord,yfit,a,eps,chi2
continuum = LEGPOLY(x,a)
;
;Get sigma from Legendre polynomial fit.  It is just square root of chi2 here.
;
sigma = SQRT(chi2)
;
;Calculate mean signal to noise
;
sn = 0.0
IF sigma NE 0.0 THEN sn = mean(yfit)/sigma
;

;Calculate the error bars on the continuum due to errors in the coefficients
;of the Legendre polynomial fit.
;
LEGERR,x,y,a,eps,chi2,continuum_err    ;chi2 = variance (=sigma^2) here
;
;Convert x and xarray back into the correct units.
;
x = x * maxx
xarray = xarray * maxx
;
;Bound continuum so that plotting errors and arithmetic errors aren't important.
;
maxy = MAX(y) & miny = MIN(y)
continuum = continuum > (-ABS(miny)*2) < (maxy*2)
;
;Print the fitting statistics to upper corner of screen.
;
PRINT,'----------------------------------------'
PRINT,'Order used = ',maxord
PRINT,'# points   = ',STRING(N_ELEMENTS(yarray),'(I8)')
PRINT,'RMS sigma  = ',STRING(sigma,'(E9.3)')
PRINT,'Mean S/N   = ',STRING(sn,'(F8.2)')
PRINT,'----------------------------------------'
ftflag=1



END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO nl_continuum_aod,filein,maxord=maxord,minord=minord,nsig1=nsig1,nsig2=nsig2,vl1=vl1,$
    vl2=vl2,vaod1=vaod1,vaod2=vaod2,contplot=contplot,aodplot=aodplot,vclip1=vclip1,blemish=blemish,vshift=vshift,$
    flgaod=flgaod,pix1=pix1,pix2=pix2,fmin=fmin,nocont=nocont,noerror=noerror,kind=kind,pdfflg=pdfflg,changename=changename

IF not keyword_set(filein)THEN begin
    print,'No input file...'
    RETURN
endif
if not keyword_set(minord) then minord=0
if not keyword_set(maxord) then maxord=10
if not keyword_set(NSIG1) then nsig1 = 3
if not keyword_set(NSIG2) then nsig2 = 0.75
if not keyword_set(vl1) then vl1 = -1000
if not keyword_set(vl2) then vl2 = +1000
if not keyword_set(vaod1) then vaod1 = -100
if not keyword_set(vaod2) then vaod2 = +100
if not keyword_set(pix1) then pix1 = 4
if not keyword_set(pix2) then pix2 = 8
if not keyword_set(fmin) then fmin = 0.5e-15
if not keyword_set(vclip1) then vclip1 = 5.
if not keyword_set(vshift) then vshift=0
if keyword_set(changename) then str = strjoin(strsplit(filein,'i.save',/regex,/extract))

    ; kind will define error if /noerror is selected. By default FP. choices are kind = 'P' (poisson noise) or kind = 'F' (fixed pattern noise)

restore,filein

ns = n_elements(f)

v = v + vshift 

; define region where the continuum will be fitted
g = where(v ge vl1 and v le vl2)
fi1 = f[g]
vi1 = v[g]
ei1= ef[g]
;===========================================================================
; remove blemish
if keyword_set(blemish) then begin 
for j = 1, n_elements(fi1) - 2 do begin
    if (fi1(j) le 0 ) then begin
	if (abs(fi1(j+1)) gt 0 and abs(fi1(j-1)) gt 0) then fi1(j) = (fi1(j-1) + fi1(j+1))/2.
    endif
endfor
endif 
;===========================================================================

;g= where(fi1 gt fmin,test)
;if test ne 0 then fm = median(fi1[where(fi1 gt fmin)])
;if test eq 0 then fm = median(fi1)

ns1 = n_elements(fi1)
mask_cont = make_array(ns1,/integer,value=1)

; cut the spectrum in sub-sections to check regions to be masked. This should help when
; flux changes by a large amount over the considered velocity interval.

diffv = (vl2-vl1)/vclip1
diffv1  = vl1
diffv2   = vl1 + diffv
j=0
for k =0,fix(vclip1)-1 do begin 
tt = where(vi1 ge diffv1 and vi1 lt diffv2)
ftemp = fi1[tt]
etemp = ei1[tt]
g= where(ftemp gt fmin,test)
if test ne 0 then fm = median(ftemp[where(ftemp gt fmin)])
if test eq 0 then fm = median(ftemp)
sigp = fm + nsig1*etemp
sigm = fm - nsig1*etemp
;    print, diffv1, diffv2, fm,pix1
    for i = pix1,n_elements(ftemp)-1,pix1 do begin
    j = j +pix1
     if (mean(ftemp[i-pix1:i]) lt mean(sigm[i-pix1:i])) then mask_cont[j-pix1:j]=0
     if (mean(ftemp[i-pix1:i]) gt mean(sigp[i-pix1:i])) then mask_cont[j-pix1:j]=0
;     print, i,j,ns1, mean(ftemp[i-pix1:i]), mean(sigp[i-pix1:i]), mean(sigm[i-pix1:i])
    endfor
diffv1  = diffv2
diffv2   = diffv1 + diffv
endfor


;for i = 0,ns1 -1 do print, vi1[i],fi1[i],mask_cont[i]




if not keyword_set(nocont) then begin


    ;;;;;;;;;;;;;;;;;;;;;;;;
    ;; CBW 2016-10-21 begin
    ;;;;;;;;;;;;;;;;;;;;;;;;

    ;;continuum fit the first time
    cbw_contfit,vi1,fi1,ycon,ycon_sig,mask_cont,minord,maxord,sn,ftflag

    ;;copy this, just in case. We don't use it elsewhere at the moment.
    mask_cont_firstfit=mask_cont


    ;;;;;;;;;;;;
    ;;Now reassess the errors as "sigma away from the continuum fit"
    ;;That is, recreate the mask after a first continuum fit
    ;;But, only look for absorption lines the second time,
    ;;  not positive spikes!!
    ;;;;;;;;;;;;

    ;;comment out; don't overwrite the mask!
    ; mask_cont = make_array(ns1,/integer,value=1)


    ;;;;;;
    ;; first method: pixel-by-pixel
    ;;;;;;

    ;;comment out; let's do this using a binned method like Nicolas did above
    ; bcontrefit=where(fi1/ycon lt (1-1.5*stdev(fi1/ycon)))
    ; mask_cont[bcontrefit] = 0


    ;;;;;;
    ;; second method: binned
    ;;;;;;

    ;;I have two ways to do this:
    ;;  using mean(efnorm)
    ;;  using stdev(fnorm)
    ;;The first uses the error vector, which may be proper
    ;;The second uses the stdev of the flux itself, which may
    ;;  produce better results, depending??


    ;;OK, *these are equivalent* as long as the coefficients (nsig2) are scaled based on "mean(efnorm)/stdev(fnorm)"
    ;;  i.e., 1-0.75*mean(efnorm) rejects the same number of points as
    ;;        1-0.593243*stdev(fnorm)
    ;;  where 0.593243 = 0.75*mean(efnorm)/stdev(fnorm)


    ;;first way -- this is safer and rejects more consistently, whereas
    ;;  with the other way, nsig2 has to be re-calculated each time
    ;;using mean(efnorm)
    ; pix2 = 2*pix1
;    for i = pix2,ns1-1,pix2 do begin
        ;;if the mean of the *binned* pixels is less than
        ;;  half of stdev(fnorm) below 1.0, then mask it
;        if (mean(fi1[i-pix2:i]/ycon[i-pix2:i]) lt (1-nsig2*mean(ei1/ycon))) then mask_cont[i-pix2:i]=0
;    endfor


    ; ;;second way
    ; ;;using stdev(fnorm)
    ; ; pix2 = 2*pix1
    ; pix2 = 8.
    ; nsig2=0.50
    ; stop;;cbw
     for i = pix2,ns1-1,pix2 do begin
    ;     ;;if the mean of the *binned* pixels is less than
    ;     ;;  half of stdev(fnorm) below 1.0, then mask it
         if (mean(fi1[i-pix2:i]/ycon[i-pix2:i]) lt (1-nsig2*stdev(fi1/ycon))) then mask_cont[i-pix2:i]=0
     endfor




    ;;continuum fit the second time, hopefully with fewer absorption lines
    cbw_contfit,vi1,fi1,ycon,ycon_sig,mask_cont,minord,maxord,sn,ftflag



    ;;;;;;;;;;;;;;;;;;;;;;;;
    ;; CBW 2016-10-21 end
    ;;;;;;;;;;;;;;;;;;;;;;;;


    vnorm = vi1
    fnorm = fi1/ycon
    eform = 0*fnorm
    eform1 = 0*fnorm
    eform2 = 0*fnorm
    sigma0 = 0*fnorm
    if not keyword_set(noerror) then begin
        efnorm =  sqrt((ycon_sig/ycon)^2.0 + (ei1/ycon)^2.)
        efnorm1 =  sqrt((ycon_sig/ycon)^2.0 + (ei1/ycon)^2.)
        efnorm2 =  (ycon_sig/ycon) + (ei1/ycon)
        efnorm =  (efnorm1+efnorm2)/2
        sigma0 = ei1/ycon
    endif
    if not keyword_set(name) then begin 
    	if keyword_set(object) then begin 
	    targname = object
	 endif else begin 
	    targname = 'file_out'
	 endelse 
    endif else begin     
    targname = name
    endelse 

if keyword_set(changename) then outname = 'plot_cont_'+str+'.eps'
if not keyword_set(changename) then outname = 'plot_cont_'+ion+wni+'.eps'

    if keyword_set(contplot) then begin

        SET_PLOT, 'PS'
        DEVICE, FILENAME=outname, FONT_SIZE=11, /INCHES, XSIZE=8, YSIZE=5, $
            /ENCAPSULATED, /COLOR, /portrait

        !P.MULTI = [0, 1, 1]
        ;!y.margin = [4,2]
        !P.CHARSIZE = 1.6
        !p.psym=10
        loadct,40,/silent


        yax = max(fi1*1e15) + max(fi1*1e15) * 0.05

        loadct, 40, /silent
        plot,vi1,fi1*1e15, xr= [vl1,vl2],yr=[-0.2,yax],psym=10,thick=2, charthick=2,/ys,/xs,$
        xtitle  = '!6Velocity (km s!e-1!n)', ytitle='Flux (10!e-15!n erg cm!e-2!n s!e-1!n !n!sA!r!u!9 %!6!N!e-1 !N)',$
             position=[0.11,0.14,0.94,0.98]
             plotzero


        g =where(mask_cont eq 0,test)
        if test ne 0 then oplot, vi1[g],fi1[g]*1e15,psym=7,thick=4,color=250

        oplot, vi1,ycon*1e15,line=0,thick=4,color=200
        ver,vaod1,lines=1,color=40, thick=4
        ver,vaod2,lines=1,color=40, thick=4

        xyouts,0.12,0.16,ion+' !7k!6'+STRTRIM(string(wavc,format='(f7.1)'),1), al=0, charsize = 2., charthick = 2,/normal,color=79

        DEVICE, /CLOSE

        set_plot,'x'
        cleanplot, /silent
    if keyword_set(pdfflg) then    spawn,'epstopdf '+outname

    endif

endif else begin
    vnorm = vi1
    fnorm = fi1
    if not keyword_set(noerror) then begin
        efnorm =  ei1
        efnorm1 = ei1
        sigma0 = ei1
    endif
    if not keyword_set(name) then begin 
    	if keyword_set(object) then begin 
	    targname = object
	 endif else begin 
	    targname = 'file_out'
	 endelse 
    endif else begin     
    targname = name
    endelse 
    ycon = make_array(n_elements(vnorm),/float,value=1.0)
    ycon_sig = make_array(n_elements(vnorm),/float,value=0.0)
    sn = 1./mean(efnorm)
endelse


if keyword_set(aodplot) then begin

if keyword_set(changename) then outname = 'plot_cont_aod_'+str+'.eps'
if not keyword_set(changename) then outname = 'plot_cont_aod_'+ion+wni+'.eps'


    SET_PLOT, 'PS'

    DEVICE, FILENAME=outname, FONT_SIZE=11, /INCHES, XSIZE=8, YSIZE=9, $
        /ENCAPSULATED, /COLOR, /portrait

    !P.MULTI = [0, 1, 2]
    ;!y.margin = [4,2]
    !P.CHARSIZE = 1.6
    !p.psym=10
    loadct,40,/silent

    if not keyword_set(nocont) then begin
        yax = max(fi1*1e15) + max(fi1*1e15) * 0.1

        fp = fi1*1e15
        yp = ycon*1e15

        plot,vi1,fp, xr= [vl1,vl2],yr=[-0.2,yax],psym=10,thick=2, charthick=2,/ys,/xs,$
        xtitle  = '', ytitle='Flux (10!e-15!n erg cm!e-2!n s!e-1!n !n!sA!r!u!9 %!6!N!e-1 !N)',$
             position=[0.13,0.58,0.94,0.99]
             plotzero

    endif else begin
        yax = max(fi1) + max(fi1) * 0.1
        fp = fi1
        yp = ycon
        ;!y.tickinterval =0.5

        plot,vi1,fp, xr= [vl1,vl2],yr=[-0.2,yax],psym=10,thick=2, charthick=2,/ys,/xs,$
        xtitle  = '', ytitle='Normalized Flux',$
             position=[0.13,0.58,0.94,0.99]
             plotzero
    endelse


    oplot, vi1,yp,line=0,thick=4,color=200
    ver,vaod1,lines=1,color=40, thick=4
    ver,vaod2,lines=1,color=40, thick=4
    g =where(mask_cont eq 0, test)
    if test ne 0 then oplot, vi1[g],fp[g],psym=7,thick=4,color=250

    xyouts,0.15,0.60,ion+' !7k!6'+STRTRIM(string(wavc,format='(f7.1)'),1), al=0, charsize = 2., charthick = 2,/normal,color=79


    g =where(vi1 ge vaod1 and vi1 le vaod2)
    ycol = -ALOG(fnorm) / (wavc*fval*2.654e-15)
    yax = max(ycol[g]*1e-13) + max(ycol[g]*1e-13) * 0.1

    plot,vi1,ycol*1e-13, xr= [vaod1-150,vaod2+150],yr=[-0.05,yax],psym=10,thick=2, charthick=2,/ys,/xs,$
    xtitle  = '!6Velocity (km s!e-1!n)', ytitle='N!da!n(v) (10!e13!n cm!e-2!n (km s!e-1!n)!e-1!n)',$
         position=[0.13,0.13,0.94,0.54]
         plotzero
              plotcolorfill, vi1[g], ycol[g]*1e-13,/midpoint, /noerase, col=234, bottom=0

    ver,vaod1,lines=1,color=40, thick=4
    ver,vaod2,lines=1,color=40, thick=4

    if keyword_set(pdfflg) then    spawn,'epstopdf '+outname


    DEVICE, /CLOSE

    set_plot,'x'
    cleanplot, /silent

endif

flux = fi1
eflux = 0 * flux
if not keyword_set(noerror) then eflux = ei1
if keyword_set(noerror) then begin
    efnorm= eflux * 0
    efnorm1= eflux * 0
    efnorm2= eflux * 0
endif

vel = vi1

if keyword_set(flgaod) then begin
    PRINT,'iCOL_BATCH::  Ion = ',ion

if keyword_set(changename) then outname = str+'i_o.save'
if not keyword_set(changename) then outname = ion+wni+'i_o.save'


    if not keyword_set(noerror) then begin
        eyflag=1
        iCOL_BATCH,vi1,fi1,ei1,ycon,sigma,ycon_sig,wavc,fval,vaod1,vaod2, ncol,ncole1,$
            ncole2,ncolez, w,w_es,w_ec,w_et,w_ez,va,vaerr,ba,baerr,m3,m3err,col3sig,col2sig,dv90,v90a,v90b,eyflag,kind

v1 = vaod1
v2 = vaod2
        save,fil=outname,v,f,ef,ion,wni,wavc,fval,gam,redshift,object,vlsr,name,ra,dec,gl,gb,z,$
            ycon,ycon_sig,vnorm,fnorm,efnorm,efnorm1,efnorm2,sigma0,targname,sn,maxord,v1,v2,vaod1,vaod2, ncol,ncole1,$
            ncole2,ncolez, w,w_es,w_ec,w_et,w_ez,va,vaerr,ba,baerr,m3,m3err,col3sig,col2sig,dv90,v90a,v90b,mask_cont,vel,flux,eflux

    endif

    if keyword_set(noerror) then begin
        eyflag =  0
        iCOL_BATCH,vi1,fi1,ei1,ycon,sigma,ycon_sig,wavc,fval,vaod1,vaod2, ncol,ncole1,$
            ncole2,ncolez, w,w_es,w_ec,w_et,w_ez,va,vaerr,ba,baerr,m3,m3err,col3sig,col2sig,dv90,v90a,v90b,eyflag,kind
v1 = vaod1
v2 = vaod2

           save,fil=outname,v,f,ef,ion,wni,wavc,fval,gam,redshift,object,vlsr,name,ra,dec,gl,gb,z,$
            ycon,ycon_sig,vnorm,fnorm,efnorm,efnorm1,efnorm2,sigma0,targname,sn,maxord,v1,v2,vaod1,vaod2, ncol,ncole1,$
            ncole2,ncolez, w,w_es,w_ec,w_et,w_ez,va,vaerr,ba,baerr,m3,m3err,col3sig,col2sig,dv90,v90a,v90b,mask_cont,vel,flux,eflux
    endif

endif else begin
    save,fil=outname,v,f,ef,ion,wni,wavc,fval,gam,redshift,object,vlsr,name,ra,dec,gl,gb,z,$
        ycon,ycon_sig,vnorm,fnorm,efnorm,targname,sn,maxord,mask_cont,vel,flux,eflux
endelse




return

end

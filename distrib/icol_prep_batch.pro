;+
;                               iCOL_PREP_BATCH.PRO
;                                 Version 1.0
;
;Program Description:
;       This procedure calculates the apparent column density, equivalent width,
;       and line statistics. 
;
;Restrictions:
;       None
;
;Screen Output:
;       Text  &  Graphics
;
;Use:
;       IMCOL_PREP,x,y,ex,ycon,sigma,ycon_sig,wavc,fval
;
;On Input:
;               x       :== velocity array
;               y       :== intensity array
;               ey      :== intensity error array
;               ycon    :== continuum array
;               sigma   :== RMS sigma of continuum fit
;               ycon_sig:== continuum fit error array
;               wavc    :== laboratory wavelength of line
;               fval    :== laboratory f-value of line
;
;On Ouptut:
;               storing inputs and outputs with possibility to save them at 
;               the end the session
;
;Common Blocks / Structures:
;       SHARE: inputs and outputs. 
;       NOISE
;       FLAGSAT
;
;Latest Update Comments:
; 04/10/13  NL    - Version 1.0
;
;External Routines Called:
;       iCOL            - to calculate column density and errors
;       iEQW            - to calculate equivalent width and errors
;       iSTAT2          - to calculate line statistics and errors
;       XLIMIT          - to determine elements of integration
;       PLOTCOLORFILL   - to plot nice fill colored histograms
;------------------------------------------------------------------------------
PRO iCOL_PREP_BATCH,x,y,ey,ycon,sigma,ycon_sig,wavc,fval,vbatch1,vbatch2
	common SHARE, ncol,ncole1,ncole2,ncolez, w,w_es,w_ec,w_et,w_ez,va,vaerr,ba,baerr,m3,m3err,v1,v2,col3sig,dv90,v90a,v90b
  common NOISE, set_noise,eyflag
  common FLAGSAT, flag_sat
  
loadct,39,/silent
        IF N_PARAMS() EQ 0 THEN BEGIN MAN,'iCcol_prep' & RETURN & ENDIF
;
;Error control
;
	ON_IOERROR,ESCAPE
;
;Find the appropriate fvalue from a standard line list.
;
;Print wavelength, f-value, and error type to screen for a user check.
;
	ycol = FLTARR(N_ELEMENTS(y))
	PRINT,'iCOL_PREP_BATCH::  Wavelength = ',STRING(wavc,'(f8.3)')
	PRINT,'iCOL_PREP_BATCH::  f-value = ',STRING(fval,'(f7.5)')
;
;
;Get limits of integration, needs to be velocities
;
    xpos1 = vbatch1 > MIN(x)  &  !c = 0
    xpos2 = vbatch2< MAX(x)  &  !c = 0
    

        IF xpos1 GT xpos2 THEN BEGIN
                PRINT,'iCOL_PREP_BATCH::  '$
                        +'Limits will be reversed for integration'
                xtmp = xpos1 & xpos1 = xpos2 & xpos2 = xtmp
        ENDIF
;
;Compute which range of elements should be included in integration and fill
;in line.  XLIMIT will find points in between limits, but be sure to set
;endpoints of integration equal to the limits.
;
    v1 = xpos1 
    v2 = xpos2  

        XLIMIT,x,xpos1,xpos2,x1,x2	
        xwork = [xpos1,x(x1:x2),xpos2]
        ywork = [INTERPOL(y,x,xpos1),y(x1:x2),INTERPOL(y,x,xpos2)]
	eywork = [INTERPOL(ey,x,xpos1),ey(x1:x2),INTERPOL(ey,x,xpos2)]
        yconwork = [INTERPOL(ycon,x,xpos1),ycon(x1:x2),INTERPOL(ycon,x,xpos2)]
        ycolwork = [INTERPOL(ycol,x,xpos1),ycol(x1:x2),INTERPOL(ycol,x,xpos2)]
        ycon_sig_work = [INTERPOL(ycon_sig,x,xpos1),ycon_sig(x1:x2)]
        ycon_sig_work = [ycon_sig_work,INTERPOL(ycon_sig,x,xpos2)]
;	
     if eyflag ne 0 then begin
        PRINT,'iCOL_PREP_BATCH::  Error Vector provided'
        y_sig_work = eywork
        kind = 'u'
        endif else begin 
        kind = 'f'
    	PRINT,'iCOL_PREP_BATCH::  Fixed pattern noise assumed'
    	y_sig_work = sigma + FLTARR(N_ELEMENTS(ywork))
         endelse 
; save noise variable 
         set_noise = kind          

;
;Calculate the column density and error by calling iCOL.
;
	iCOL,xwork,ywork,yconwork,y_sig_work,ycon_sig_work,wavc,fval,col,$
		y_err,ycon_err,zero_err
	tot_err = SQRT(y_err^2 + ycon_err^2)
	IF col LE 0.0 THEN col=1.0

;
;Print information dump.
;
if flag_sat eq 0 then begin
        PRINT,'--------------------------------------------'
	PRINT,'log N (best) =   ',STRING(ALOG10(col),'(F9.3)')
  PRINT,'(+1 sig)     =   ',STRING(ALOG10(col+tot_err)-ALOG10(col),'(F9.3)')
  PRINT,'(-1 sig)     =   ',STRING(-ALOG10(col-tot_err)+ALOG10(col),'(F9.3)')
;	PRINT,'2% zero err  =   ',STRING(ALOG10(col+zero_err)-ALOG10(col),'(F9.3)')
        PRINT,'--------------------------------------------'
	ncol   = alog10(col)
	ncole1 = ALOG10(col+tot_err) - ALOG10(col)
	ncole2 = -ALOG10(col-tot_err) + ALOG10(col)
	ncolez = ALOG10(col+zero_err)-ALOG10(col)
endif 
if flag_sat eq 1 then begin
        PRINT,'--------------------------------------------'
    PRINT,'log N   > ',STRING(ALOG10(col),'(F9.3)')
   PRINT,'(+1 sig)     =   ',STRING(ALOG10(col+tot_err)-ALOG10(col),'(F9.3)')
  PRINT,'(-1 sig)     =   ',STRING(-ALOG10(col-tot_err)+ALOG10(col),'(F9.3)')   
        PRINT,'--------------------------------------------'
  ncol   = alog10(col)
  ncole1 = ALOG10(col+tot_err) - ALOG10(col)
  ncole2 = -ALOG10(col-tot_err) + ALOG10(col)
  ncolez = 1
endif	
;------------------------------------------------------------------------------
;
xpos1 = v1
xpos2 = v2 
        XLIMIT,x,xpos1,xpos2,x1,x2
        xwork    = x(x1:x2)
        ywork    = y(x1:x2)
        eywork   = ey(x1:x2)
        yconwork = ycon(x1:x2)
        ycon_sig_work = ycon_sig(x1:x2)
 
;
;Ask about what kind of statistics should be used (poisson or fixed pattern)
;to describe the noise characteristics of the data.
;
 
        IF kind EQ 'f' THEN BEGIN
                y_sig_work = sigma + FLTARR(N_ELEMENTS(ywork))
        endif else if kind eq 'u' then begin 
        y_sig_work = eywork
        endif else if kind eq 'p' then begin 
             y_sig_work = sigma * SQRT(abs(ywork)/yconwork)
         endif  


;Calculate line statistics.
;
  iSTAT2,xwork,ywork,yconwork,y_sig_work,ycon_sig_work,m1,m1err, m2,m2err,m3,m3err,m4,m4err,dv90,v90a,v90b
;
;Print information dump.
;

  va = m1
  vaerr = m1err
  ba = m2 * sqrt(2)
  baerr = m2err * sqrt(2)

        PRINT,'----------------------------------------------'
  PRINT,"$('<v>       = ',f8.3,'  +/- ',f7.3)",m1,m1err
  PRINT,"$('<b>       = ',f8.3,'  +/- ',f7.3)",m2* sqrt(2),m2err* sqrt(2)
  PRINT,"$('dv90       = ',f8.3,'  +/- ',f7.3)",dv90,m1err* sqrt(2)
  PRINT,"$('Skew      = ',f8.3,'  +/- ',f7.3)",m3,m3err
  PRINT,'----------------------------------------------'
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;
;Calculate equivalent width and associated errors (cont ,stat, tot).
; change first velocity to wavelength 
;
xpos1 = v1
xpos2 = v2 
               xpos1 = iAXIS(xpos1,wavc,1) 
               xpos2 = iAXIS(xpos2,wavc,1)
               x  = iAXIS(x,wavc,1) 
        XLIMIT,x,xpos1,xpos2,x1,x2
        xwork = [xpos1,x(x1:x2),xpos2]
        ywork = [INTERPOL(y,x,xpos1),y(x1:x2),INTERPOL(y,x,xpos2)]
        eywork = [INTERPOL(ey,x,xpos1),ey(x1:x2),INTERPOL(ey,x,xpos2)]
        yconwork = [INTERPOL(ycon,x,xpos1),ycon(x1:x2),INTERPOL(ycon,x,xpos2)]
        ycon_sig_work = [INTERPOL(ycon_sig,x,xpos1),ycon_sig(x1:x2)]
        ycon_sig_work = [ycon_sig_work,INTERPOL(ycon_sig,x,xpos2)]
 
         IF kind EQ 'f' THEN BEGIN
                y_sig_work = sigma + FLTARR(N_ELEMENTS(ywork))
        endif else if kind eq 'u' then begin 
        y_sig_work = eywork
        endif else if kind eq 'p' then begin 
             y_sig_work = sigma * SQRT(abs(ywork)/yconwork)
         endif  
        
        iEQW,xwork,ywork,yconwork,y_sig_work,ycon_sig_work,ew,y_err,ycon_err,zero_err
        toterr = SQRT(y_err^2 + ycon_err^2)
;
;Calculate linear column density and error.
;
        col = (ew/wavc)/8.85e-13/(wavc*1.e-8)/fval
        colerr = col - ((ew-toterr)/wavc)/8.85e-13/(wavc*1.e-8)/fval
        col3sig = alog10(3* (toterr/wavc)/8.85e-13/(wavc*1.e-8)/fval)

;
;Print information dump.  Equivalent widths in mA.
; 
  w    = ew * 1000.0
  w_es =  y_err * 1000.0
  w_ec =  ycon_err * 1000.0
  w_et =  toterr * 1000.0
  w_ez = zero_err* 1000.0
        PRINT,'--------------------------------------------'
        PRINT,'EW           = ', STRING(ew * 1000.0,'(F9.2)')
        PRINT,'Stat Error   = ',STRING(y_err * 1000.0 ,'(F9.2)')
        PRINT,'Cont Error   = ',STRING(ycon_err * 1000.0,'(F9.2)')
        PRINT,'Tot Error    = ',STRING(toterr * 1000.0 ,'(F9.2)')
   ;     PRINT,'2% Zero Err  = ',STRING(zero_err * 1000.0,'(F9.2)')
        PRINT,'Linear COG N = ',STRING(ALOG10(col),'(F9.4)')
        PRINT,'3sigma EW    < ',STRING(toterr * 3000.0,'(F9.2)')
        PRINT,'3sigma N     < ',STRING(col3sig,'(F9.2)')
        PRINT,'--------------------------------------------'
 ;return to velocity                    
  x  = iAXIS(x,wavc,-1) 
 !ytitle = 'Flux'
         RETURN 
ESCAPE:
        PRINT,'iCOL_PREP_BATCH:: '+!err_string
        RETURN  &  END


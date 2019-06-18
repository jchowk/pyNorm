; +
; nl_do_cont_aod
;   Version 1.0

; PURPOSE:
; This program automatically continuum fits a bunch of transtion and estimate their AOD parameters
;
; CALLING SEQUENCE:
;  
;   nl_continuum_leg, '*i.save',order=[low,high]
;
; OPTIONAL KEYWORDS:
;
; EXAMPLE:
;   IDL> nl_continuum_leg, '*i.save',order=[low,high]
;
; OUTPUTS 
;
;;External Routines Called:
;   nl_continuum_aod.pro
;   
; REVISION HISTORY:
;   v1.0 - 25-April-2015 Written by NL;
;-
;------------------------------------------------------------------------------

pro nl_do_cont_aod,filein=filein,maxord=maxord,minord=minord,nsig1=nsig1,nsig2=nsig2,vl1=vl1,vclip1=vclip1,$
    	    	    vl2=vl2,vaod1=vaod1,vaod2=vaod2,contplot=contplot,aodplot=aodplot,flgaod=flgaod,blemish=blemish,vshift=vshift,$
		    pix1=pix1,pix2=pix2,fmin=fmin,nocont=nocont,noerror=noerror,kind=kind,pdfflg=pdfflg,changename=changename
    
    IF not keyword_set(filein)THEN begin 
    print,'No input file...'
    RETURN
    endif


; choice: filein can be a single *i.save file or a file that contains a list of *i.save file.     
  g = where(strmatch(filein,'*i.save') eq 1 ,test)
  
  if test eq 0 then begin 
    readcol,filein,files,format='a'
    	for i =0,n_elements(files) -1 do begin
    	nl_continuum_aod,files[i],maxord=maxord,minord=minord,nsig1=nsig1,nsig2=nsig2,vl1=vl1,vclip1=vclip1,$
    	    	vl2=vl2,vaod1=vaod1,vaod2=vaod2,contplot=contplot,aodplot=aodplot,blemish=blemish,vshift=vshift,$
		flgaod=flgaod,pix1=pix1,pix2=pix2,fmin=fmin,nocont=nocont,noerror=noerror,kind=kind,pdfflg=pdfflg,changename=changename

    endfor 

    endif
    
  if test eq 1 then begin 
    nl_continuum_aod,filein,maxord=maxord,minord=minord,nsig1=nsig1,nsig2=nsig2,vl1=vl1,vclip1=vclip1,$
    	    	vl2=vl2,vaod1=vaod1,vaod2=vaod2,contplot=contplot,aodplot=aodplot,blemish=blemish,vshift=vshift,$
		flgaod=flgaod,pix1=pix1,pix2=pix2,fmin=fmin,nocont=nocont,noerror=noerror,kind=kind,pdfflg=pdfflg,changename=changename
    endif
    
    


return

end

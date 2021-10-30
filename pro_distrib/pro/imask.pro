;+
;				iMask.PRO
;				Version 1.0
;
;Program Description:
;	This procedure serves to mask part of the spectrum
;
;Restrictions:
;       Requires cursor (3 button mouse) input
;
;Screen Output:
;	Graphics  &  Text
;
;Use:
;	imask,x,y,ey,xarray,yarray,mask,maskflg
;
;On Input:
;		x	:== x coordinate array
;   	    	y   	:== y coordinate array
;   	    	ey  	:== ey error coordinate array
;		mask 	:==  mask region
;		maskflg	:== mask flag
;
;On Output:
;		mask     :==  mask region
;		maskflg  :== mask flag
;
;Common Blocks / Structures:
;	None
;
;Latest Update Comments:
; 04/17/15  NL  - Version 1.0
;
;External Routines Called:
;----------------------------------------------------------------------------
PRO imask,x,y,ey,mask,maskflg

  	;; ;; ;; IF N_PARAMS() EQ 0 THEN BEGIN MAN,'imask' & RETURN & ENDIF
;
;Error control.
;
	; ON_IOERROR,ESCAPE
	; IF N_PARAMS() NE 7 THEN maskflg = 0
;
;Print instructions.
;
	PRINT,'iMask::  Mark (C1)   Clear (C2)   Quit (C3)'


;; LOOP:
;
;Determine how many regions have already been masked and show them.
;
loadct,39,/silent
	liney1 = [1,-1]*!y.crange(1)/30. + 0.9*!y.crange(1)
	IF maskflg EQ 0 THEN mask = x/x * 1.
	IF maskflg EQ 1 THEN BEGIN
		bad = where(mask eq 0,ct)
	    if ct ne 0 then oplot,x[bad], y[bad], color=150,psym=7, thick=3
	ENDIF
;
;Get input from user.  If right mouse button is pushed, the return.
;If middle button is pushed, then clear.
;
	CURSOR,xpos1,ypos1,/DOWN
	xpos1 = xpos1 > MIN(x)  &  !c = 0
	; IF !err EQ 4 THEN RETURN
	IF !err EQ 2 THEN BEGIN
		maskflg = 0
		PLOT,x,y
		;; GOTO,;; LOOP
	ENDIF
;
;Mark left of region.
;
	PRINT,'iMask::  Left limit:  ',xpos1
	OPLOT,[xpos1,xpos1],liney1,color=200
;
;Get and mark right of region.
;
	CURSOR,xpos2,ypos2,/DOWN
	xpos2 = xpos2 < MAX(x)  &  !c = 0
  	PRINT,'iMask::  Right limit: ',xpos2
	OPLOT,[xpos2,xpos2],liney1,color=200
;
;Check to be sure range is not reversed.  If reversed, get again.
;
	IF xpos2 GT xpos1 THEN BEGIN
		OPLOT,[xpos1,xpos2],[1,1]*0.9*!y.crange(1)	,color=200
	ENDIF ;ELSE ;; GOTO,;; LOOP
;
;Detrmine which regions is masked
;
	XLIMIT,x,xpos1,xpos2,x1,x2
  bad = where(x ge xpos1 and x le xpos2)
  mask[bad] = 0
  maskflg  = 1

	;; GOTO,;; LOOP
;----------------------------------------------------------------------------
; ESCAPE:
; 	PRINT,'iMask::  '+!err_string
; 	RETURN
END

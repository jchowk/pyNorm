"""
 nl_do_cont_aod
   Version 1.0

 PURPOSE:
 This program automatically continuum fits a bunch of transtion and estimate their AOD parameters

 CALLING SEQUENCE:

   nl_continuum_leg, '*i.save',order=[low,high]

 OPTIONAL KEYWORDS:

 EXAMPLE:
   IDL> nl_continuum_leg, '*i.save',order=[low,high]

 OUTPUTS

#External Routines Called:
   nl_continuum_aod.pro

 REVISION HISTORY:
   v1.0 - 25-April-2015 Written by NL;
-
"""

from numpy import *

def nl_do_cont_aod(filein=None, maxord=None, minord=None, nsig1=None, nsig2=None, vl1=None, vclip1=None, vl2=None, vaod1=None, vaod2=None, contplot=None, aodplot=None, flgaod=None, blemish=None, vshift=None, pix1=None, pix2=None, fmin=None, nocont=None, noerror=None, kind=None, pdfflg=None, changename=None):

    n_params = 0
    _opt = (filein, maxord, minord, nsig1, nsig2, vl1, vclip1, vl2, vaod1, vaod2, contplot, aodplot, flgaod, blemish, vshift, pix1, pix2, fmin, nocont, noerror, kind, pdfflg, changename)
    def _ret():
        _optrv = zip(_opt, [filein, maxord, minord, nsig1, nsig2, vl1, vclip1, vl2, vaod1, vaod2, contplot, aodplot, flgaod, blemish, vshift, pix1, pix2, fmin, nocont, noerror, kind, pdfflg, changename])
        _rv = [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    if bitwise_not((filein is not None)):    
        print 'No input file...'
        return _ret()
    
    
    # choice: filein can be a single *i.save file or a file that contains a list of *i.save file.
    g = where(ravel(strmatch(filein, '*i.save') == 1))[0]
    
    if test == 0:    
        readcol(filein, files, format='a')
        for i in arange(0, (files.size - 1)+(1)):
            nl_continuum_aod(files[i], maxord=maxord, minord=minord, nsig1=nsig1, nsig2=nsig2, vl1=vl1, vclip1=vclip1, vl2=vl2, vaod1=vaod1, vaod2=vaod2, contplot=contplot, aodplot=aodplot, blemish=blemish, vshift=vshift, flgaod=flgaod, pix1=pix1, pix2=pix2, fmin=fmin, nocont=nocont, noerror=noerror, kind=kind, pdfflg=pdfflg, changename=changename)
            
        
    
    if test == 1:    
        nl_continuum_aod(filein, maxord=maxord, minord=minord, nsig1=nsig1, nsig2=nsig2, vl1=vl1, vclip1=vclip1, vl2=vl2, vaod1=vaod1, vaod2=vaod2, contplot=contplot, aodplot=aodplot, blemish=blemish, vshift=vshift, flgaod=flgaod, pix1=pix1, pix2=pix2, fmin=fmin, nocont=nocont, noerror=noerror, kind=kind, pdfflg=pdfflg, changename=changename)
    
    
    
    
    return _ret()
    


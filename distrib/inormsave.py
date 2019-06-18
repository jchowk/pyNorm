"""
                               iNORMSAVE.PRO
                               Version 1.0;
Program Description:
 This procedure will create a list of inorm save files from
 a list of ions.

Latest Update Comments:
       04/15/13  NL: Version 1.0

External Routines called:
       None
"""

from __future__ import print_function

from numpy import *

def iNORMSAVE(root=None, redshift=None, vlsr=None, object=None):

    n_params = 0
    _opt = (root, redshift, vlsr, object)
    def _ret():
        _optrv = zip(_opt, [root, redshift, vlsr, object])
        _rv = [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    if bitwise_not((object is not None)):    
        object = ''
    if bitwise_not((redshift is not None)):    
        redshift = 0
    if bitwise_not((vlsr is not None)):    
        vlsr = 0
    if bitwise_not((root is not None)):    
        root = ''
    
    ion = ''
    vlsr = 0.
    wavc = 0.
    fval = 0.
    gam = 0.
    wni = ''
    rext = ''
    flgf = -1
    if (root is not None):    
        root = STRTRIM(root, 2)
        # detemine file type
        rext1 = array([['.save', '.sav', '.fits', '.txt', '.dat', '.asc']])
        flgfile = array([[0, 0, 1, 2, 2, 2]])
        for j in arange(0, (rext1.size - 1)+(1)):
            test = file_test(root + rext1[j])
            if (test == 1):    
                rext = rext1[j]
                flgf = flgfile[j]
    
    if flgf == -1:    
        print( 'iSAVE::  Enter  type of files')
        print( 'iSAVE::  (s)ave, (x).fits, (a)scii or (q)uit')
        kind = GET_KBRD(1)
        rext = ''
        if kind == 'q':    
            return _ret()
        if kind == 's':    
            flgf = 0
        if kind == 'x':    
            flgf = 1
        if kind == 'a':    
            flgf = 2
    #
    if flgf == 0:    
        print( 'iNORM::  Reading: ', root + rext)
        restore(root + rext)
    if flgf == 1:    
        # assume xidl fits file
        fluxi = file_search('*f.fits')
        fluxi = fluxi[0]
        erri = file_search('*e.fits')
        erri = erri[0]
        err = xmrdfits(erri, 0, hdr, silent=True)
        flux = xmrdfits(fluxi, 0, hdr, silent=True)
        object = sxpar(hdr, 'TARGNAME')
        wave = 10 ** (sxpar(hdr, 'CRVAL1') + arange(flux.size, dtype=float32) * sxpar(hdr, 'cdelt1'))
    if flgf == 2:    
        print( 'iNORM::  Reading: ', root + rext)
        readcol(root + rext, wave, flux, err, silent=True)
    
    readcol('~/Dropbox/IDL/iNorm/lists/lineslls.dat', elem, w, gv, fv, format='A,D,D,D')
    
    for i in arange(0, (fv.size - 1)+(1)):
        wni = STRTRIM(str(w(i), '(f8.1)'), 2)
        ion = STRTRIM(elem(i), 2)
        wavc = w(i)
        fval = fv(i)
        gam = gv(i)
        z = redshift
        
        vlsr1 = vlsr
        if z > 0:    
            vlsr1 = 0 # helio frame for z>0 absorbers
        
        if (bitwise_and(wavc * (1 + z) >= array(wave, copy=0).min(), wavc * (1 + z) <= array(wave, copy=0).max())):    
            
            vel = ((wave - wavc) / wavc * 2.9979e5 - 2.9979e5 * z) / (1. + z) + vlsr1
            
            index = where(ravel(bitwise_and(vel >= -2000, vel <= 2000)))[0]
            if ct != 0:    
                v = vel[index]
                f = flux[index]
                ef = err[index]
                save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, fil=ion + wni + 'i.save')
            
        else:    
            v = 0. + zeros([1], "float32")
            f = 0. + zeros([1], "float32")
            ef = 0. + zeros([1], "float32")
            save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, fil=ion + wni + 'i.save')
        
    
    ##plot
    
    print( 'Done......................')
    
    
    return _ret()


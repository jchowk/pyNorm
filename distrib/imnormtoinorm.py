"""
                               iMNORMTOINORM.PRO
                               Version 1.0;
Program Description:
 This procedure will update files *o.save files to *i.save files
 with updated variable names.

Latest Update Comments:
       04/15/13  NL: Version 1.0

External Routines called:
       None
"""

from __future__ import print_function

from numpy import *

def iMNORMTOINORM(name):
    n_params = 1
    def _ret():  return name
    
    if bitwise_not((name is not None)):    
        name = ''
    
    file = file_search('*o.save')
    
    
    for i in arange(0, (file.size - 1)+(1)):
        restore(file(i))
        if vlsr.size == 0:    
            vlsr = 0.
        ion = el
        gam = gv
        fval = fv
        object = name
        redshift = z
        wavc = wc
        wni = STRTRIM(str(wavc, '(f8.1)'), 2)
        
        save(v, f, ef, ion, wni, wavc, fval, gam, redshift, object, vlsr, fil=ion + wni + 'i.save')
        
    
    print( 'Done......................')
    
    
    return _ret()


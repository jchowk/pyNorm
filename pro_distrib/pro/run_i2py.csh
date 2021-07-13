#!/bin/csh

set files = ('*pro')

foreach ff ($files)
        echo "---  "+$ff+"  ---"
        idl2python $ff

end

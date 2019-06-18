#!/bin/csh

set files = ('*pro')
set mancall_in = "IF N_PARAMS() EQ 0 THEN"
set mancall_out = ";; IF N_PARAMS() EQ 0 THEN"

foreach ff ($files)
        gsed -i "s/$mancall_in/$mancall_out/g" $ff


end

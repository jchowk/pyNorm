#!/bin/csh

set files = ('*pro')
set mancall_in = "IF N_PARAMS() EQ 0 THEN"
set mancall_out = ";; IF N_PARAMS() EQ 0 THEN"

set goto_in = "GOTO"
set goto_out = ";; GOTO"

set loop_in = "LOOP"
set loop_out = ";; LOOP"

foreach ff ($files)
        gsed -i "s/$mancall_in/$mancall_out/g" $ff
        gsed -i "s/$goto_in/$goto_out/g" $ff
        gsed -i "s/$loop_in/$loop_out/g" $ff


end

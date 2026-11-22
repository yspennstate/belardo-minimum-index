#!/bin/bash
# Data of the last four rows of Table 2: all connected graphs with 13 vertices and 15 edges (geng -c), and for the other
# three pairs only those with all degrees in {2,3} (geng -c -d2 -D3). Set GENG to the geng binary if it is not on the PATH.
GENG=${GENG:-geng}
PY=${PY:-python3}
e=15
"$GENG" -cq 13 "$e:$e" | "$PY" check_min_index.py exhaustive 13 2
for ns in "14 4" "15 2" "17 2"; do
    set -- $ns
    e=$(( $1 + $2 ))
    echo "degrees 2..3 only:"
    "$GENG" -cq -d2 -D3 "$1" "$e:$e" | "$PY" check_min_index.py exhaustive "$1" "$2"
done

#!/bin/bash
# Data of the table: all connected graphs with n vertices and n+s edges (geng -c), and for the last four rows only
# those with all degrees in {2,3} (geng -c -d2 -D3). Set GENG to the geng binary if it is not on the PATH.
GENG=${GENG:-geng}
PY=${PY:-python3}
for ns in "4 1" "4 2" "5 1" "5 2" "6 1" "6 2" "6 3" "7 1" "7 2" "7 3" "8 1" "8 2" "8 3" "8 4" "9 1" "9 2" "9 3" "9 4" \
          "10 1" "10 2" "10 3" "10 4" "10 5" "11 1" "11 2" "11 3" "12 2" "12 3" "12 4"; do
    set -- $ns
    e=$(( $1 + $2 ))
    "$GENG" -cq "$1" "$e:$e" | "$PY" check_min_index.py exhaustive "$1" "$2"
done
for ns in "14 2" "15 3" "16 4" "14 7"; do
    set -- $ns
    e=$(( $1 + $2 ))
    echo "degrees 2..3 only:"
    "$GENG" -cq -d2 -D3 "$1" "$e:$e" | "$PY" check_min_index.py exhaustive "$1" "$2"
done

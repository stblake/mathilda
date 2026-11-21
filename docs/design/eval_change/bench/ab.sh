#!/bin/bash
# usage: CASE=name ab.sh reps [binary-for-000]
D=$(cd $(dirname $0); pwd); B=$D/../inst/Mathilda
for r in $(seq $1); do
  for v in "0 0 0" "1 1 0" "1 1 1"; do set -- $v
    out=$(MTH_P_SYMSUM=$1 MTH_P_RBSUM=$2 MTH_P_NOBUMP=$3 /usr/bin/time -l $B -file $D/ab_$CASE.m 2>&1)
    echo "[$1$2$3] $(echo "$out" | grep cpu | cut -c1-110) | Ginstr $(echo "$out" | awk '/instructions retired/{printf "%.3f", $1/1e9}') Gcyc $(echo "$out" | awk '/cycles elapsed/{printf "%.3f", $1/1e9}')"
  done
done

#!/bin/bash
# corpus_ab.sh VARIANT(4 digits: symsum rbsum nobump binder) -> one line per case
D=$(cd $(dirname $0); pwd); B=$D/../inst/Mathilda; V=$1; set -- $(echo $1 | sed 's/./& /g')
for f in $D/corpus/*.m; do
  out=$(MTH_P_SYMSUM=$1 MTH_P_RBSUM=$2 MTH_P_NOBUMP=$3 MTH_P_BINDER=$4 /usr/bin/time -l $B -file $f 2>&1)
  echo "[$V] $(echo "$out" | grep ' cpu ' | cut -c1-120) | Ginstr $(echo "$out" | awk '/instructions retired/{printf "%.3f", $1/1e9}')"
done

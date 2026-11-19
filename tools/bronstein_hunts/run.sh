#!/bin/bash
# Run every Bronstein / Risch hunt, instrumented, in parallel, against a binary
# (default ./Mathilda), and write per-hunt verdict summaries and the ranking of
# rde_cert_mark reasons over the UNDECIDED cases (blocking.py) and their exact
# mark sets (solo.py; "NO MARK" = a silent decline the completeness ratchet
# cannot see).  Usage: tools/bronstein_hunts/run.sh [binary] [outdir]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="${1:-$HERE/../../Mathilda}"
OUT="${2:-$HERE/out}"
mkdir -p "$OUT"
HUNTS="hunt hunt2 hunt_exp f8hunt f9hunt2 f12hunt dephunt rdehunt primhunt primhunt2 tanhhunt tanihunt multitanhunt"
WRAP='SetAttributes[eiq, HoldAll]; eiq[f_, x_] := Module[{v}, WriteString["stderr", "CASE\n"]; v = Risch`ElementaryIntegralQ[f, x]; WriteString["stderr", "VERDICT ", If[Head[v] === Risch`ElementaryIntegralQ, "undec", ToString[v]], "\n"]; v];'
for h in $HUNTS; do
  { echo "$WRAP"; sed 's/Risch`ElementaryIntegralQ\[/eiq[/g' "$HERE/$h.m"; } > "$OUT/blk_$h.m"
  ( ulimit -v 16000000 2>/dev/null; MATHILDA_RDE_CERT_DEBUG=1 timeout 3000 "$BIN" -file "$OUT/blk_$h.m" \
      > "$OUT/$h.out" 2> "$OUT/$h.err" ) &
done
wait
for h in $HUNTS; do
  echo "== $h"
  grep -vE "::undec|^$" "$OUT/$h.out"
  python3 "$HERE/solo.py" < "$OUT/$h.err"
done

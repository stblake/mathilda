#!/usr/bin/env bash
# Benchmark NDArray reductions/sort in Mathilda against NumPy on a 1e7 float64
# array. Run from the repo root: comparisons/ndarray_reductions_bench.sh
# Requires: ./Mathilda built, and a python with numpy (edit PY below if needed).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3.11}
N=${N:-10000000}

echo "== Mathilda ($N float64, best of 5; Median/Sort best of 3) =="
printf '%s\n' \
  "{\"id\":0,\"expr\":\"nd = NDArray[RandomReal[{0,1}, $N]]; NDArrayQ[nd]\"}" \
  '{"id":1,"expr":"{\"Total\",  Min[Table[First[Timing[Total[nd]]], {5}]]}"}' \
  '{"id":2,"expr":"{\"Mean\",   Min[Table[First[Timing[Mean[nd]]], {5}]]}"}' \
  '{"id":3,"expr":"{\"StdDev\", Min[Table[First[Timing[StandardDeviation[nd]]], {5}]]}"}' \
  '{"id":4,"expr":"{\"Max\",    Min[Table[First[Timing[Max[nd]]], {5}]]}"}' \
  '{"id":5,"expr":"{\"Median\", Min[Table[First[Timing[Median[nd]]], {3}]]}"}' \
  '{"id":6,"expr":"{\"Sort\",   Min[Table[First[Timing[Sort[nd];]], {3}]]}"}' \
  | ./Mathilda 2>/dev/null | grep -oE '"payload":"[^"]*"'

echo
echo "== NumPy =="
"$PY" - "$N" <<'PY'
import sys, time, numpy as np
n = int(sys.argv[1]); x = np.random.random(n)
def b(f, r=5):
    return min((lambda t0: (f(), time.perf_counter()-t0)[1])(time.perf_counter()) for _ in range(r))
for name, f, r in [("Total", x.sum, 5), ("Mean", x.mean, 5),
                   ("StdDev", lambda: x.std(ddof=1), 5), ("Max", x.max, 5),
                   ("Median", lambda: np.median(x), 3), ("Sort", lambda: np.sort(x), 3)]:
    print(f"{name:8s} {b(f, r)*1e3:8.2f} ms")
PY

# Mathematica (optional): set WOLFRAMSCRIPT, or auto-detect the app bundle.
WS=${WOLFRAMSCRIPT:-/Applications/Mathematica.app/Contents/MacOS/wolframscript}
if [ -x "$WS" ]; then
  echo
  echo "== Mathematica =="
  "$WS" -code "
    n = $N; x = RandomReal[1, n];
    best[f_, r_] := Min[Table[First[AbsoluteTiming[f[]; ]], {r}]];
    Print[\"Total    \", best[Total[x]&, 5]*1000., \" ms\"];
    Print[\"Mean     \", best[Mean[x]&, 5]*1000., \" ms\"];
    Print[\"StdDev   \", best[StandardDeviation[x]&, 5]*1000., \" ms\"];
    Print[\"Max      \", best[Max[x]&, 5]*1000., \" ms\"];
    Print[\"Median   \", best[Median[x]&, 3]*1000., \" ms\"];
    Print[\"Sort     \", best[Sort[x]&, 3]*1000., \" ms\"];
  " 2>/dev/null
else
  echo; echo "(wolframscript not found at \$WS=$WS; skipping Mathematica)"
fi

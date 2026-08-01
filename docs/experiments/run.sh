#!/bin/bash
# Run one experiment's kernels in all three systems and print the three
# outputs side by side.
#
#     docs/experiments/run.sh 12-graph-and-sparse            # all three
#     docs/experiments/run.sh 12-graph-and-sparse mathilda   # just one
#
# Every experiment folder holds the same three files:
#
#     README.md      the write-up: what was measured, what it found, the
#                    roadmap for any row where Mathilda is not the fastest
#     <slug>.m       the kernels, runnable unmodified in Mathilda AND in
#                    Mathematica -- one source text, so the two CAS columns
#                    cannot silently be timing different programs
#     <slug>.py      the same algorithm in NumPy, in the same order
#
# All three systems take a script file directly: Mathilda and Mathematica both
# spell it `-file <path>`, python takes the path on its own.
#
# Environment:
#     MATHILDA_BIN   path to the Mathilda binary (default: ./Mathilda)
#     WOLFRAMSCRIPT  path to wolframscript
#     HPC_PYTHON     a python3 that has numpy and scipy
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="$ROOT/docs/experiments"

MATHILDA="${MATHILDA_BIN:-$ROOT/Mathilda}"
WOLFRAM="${WOLFRAMSCRIPT:-/Applications/Mathematica.app/Contents/MacOS/wolframscript}"
PYTHON="${HPC_PYTHON:-/usr/local/bin/python3.11}"

if [ $# -lt 1 ]; then
  echo "usage: run.sh <experiment-folder> [mathilda|mathematica|python]" >&2
  echo "" >&2
  echo "available:" >&2
  ls -1 "$EXP" | grep -E '^[0-9]{2}-' | sed 's/^/  /' >&2
  exit 2
fi

DIR="$EXP/$1"
WHICH="${2:-all}"
[ -d "$DIR" ] || { echo "no such experiment: $1" >&2; exit 2; }

MFILE=$(ls "$DIR"/*.m 2>/dev/null | head -1)
PYFILE=$(ls "$DIR"/*.py 2>/dev/null | head -1)

run_mathilda() {
  [ -n "$MFILE" ] || { echo "(no .m file)"; return; }
  [ -x "$MATHILDA" ] || { echo "(no Mathilda binary at $MATHILDA -- run make)"; return; }
  echo "===== Mathilda ====="
  "$MATHILDA" -file "$MFILE" 2>&1
}

run_mathematica() {
  [ -n "$MFILE" ] || { echo "(no .m file)"; return; }
  [ -x "$WOLFRAM" ] || { echo "(no wolframscript at $WOLFRAM)"; return; }
  echo "===== Mathematica ====="
  "$WOLFRAM" -file "$MFILE" 2>&1
}

run_python() {
  [ -n "$PYFILE" ] || { echo "(no .py file)"; return; }
  command -v "$PYTHON" >/dev/null 2>&1 || { echo "(no python at $PYTHON)"; return; }
  echo "===== Python / NumPy ====="
  "$PYTHON" "$PYFILE" 2>&1
}

case "$WHICH" in
  mathilda)     run_mathilda ;;
  mathematica)  run_mathematica ;;
  python|numpy) run_python ;;
  all)          run_mathilda; echo; run_mathematica; echo; run_python ;;
  *) echo "unknown system: $WHICH" >&2; exit 2 ;;
esac

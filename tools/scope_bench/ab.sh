#!/bin/bash
# Instruction-count A/B over the scoping corpus (docs/design/scoping_substitution.md).
#
#   tools/scope_bench/ab.sh [binary] [ENV=VAL ...]
#
# Runs every tools/scope_bench/corpus/*.m case with the given environment and
# prints one line per case: CPU, wall, answer hash, and instructions retired
# (billions). Compare two runs, e.g. the default against MATHILDA_NO_SYMSUM=1
# (summaries off) or MATHILDA_MODULE_CLOCK_BUMP=1 (the old Module clock bump);
# an unchanged hash means an unchanged evaluation path, so the instruction
# delta is the change's own. Run it on a quiet machine: the Integrate cascade is
# wall-clock budgeted, and under load a case can take a different path.
# macOS: /usr/bin/time -l; Linux: perf stat -e instructions.
D=$(cd "$(dirname "$0")" && pwd)
B=${1:-$D/../../Mathilda}; shift
for f in "$D"/corpus/*.m; do
  if [ "$(uname)" = Darwin ]; then
    out=$(env "$@" /usr/bin/time -l "$B" -file "$f" 2>&1)
    gi=$(echo "$out" | awk '/instructions retired/{printf "%.3f", $1/1e9}')
  else
    out=$(env "$@" perf stat -x, -e instructions "$B" -file "$f" 2>&1)
    gi=$(echo "$out" | awk -F, '/instructions/{printf "%.3f", $1/1e9}')
  fi
  echo "$(echo "$out" | grep ' cpu ' | cut -c1-110) | check $(echo "$out" | awk '/^check /{print $2}') | Ginstr $gi"
done

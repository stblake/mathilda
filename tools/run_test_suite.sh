#!/bin/sh
# run_test_suite.sh -- build and run EVERY test binary, and report against a baseline.
#
# WHY THIS EXISTS. The suite is 400+ separate binaries and nothing ran all of them: CI compiled the
# tree and ran two source-level gates, and a human deciding to "run the tests" ran whichever subset
# they had built. A PR then claimed the tests passed on the strength of 40 of 426 binaries. The
# subset is the problem, so this script removes the choice: it builds every target and runs every
# binary, once, with the same configuration each time.
#
# THE BASELINE IS WHAT MAKES IT USABLE. Some tests fail for reasons unrelated to whatever change is
# being made, and a gate that is red on arrival gets ignored within a week -- so `tests/known_
# failures.txt` lists them with a one-line reason each. A failure NOT in that list fails this
# script; a listed test that has started passing is reported so the line can be deleted. That is the
# same ratchet the repo already uses for OFF_BUFFER and BASELINE in the packing gates.
#
# CONFIGURATION IS PINNED, NOT INHERITED. The failures that wasted the most time were configuration
# artifacts: a CMake tree left at USE_FLINT=OFF reported 36 failures, of which 32 were the missing
# library rather than a defect, and one target could not even LINK. So the configuration is spelled
# out here and printed at the top of the run -- a result is only meaningful alongside the flags that
# produced it.
#
# Usage:  tools/run_test_suite.sh [build-dir]
# Exit:   0 = no unexpected failures; 1 = at least one; 2 = the build itself failed.

set -e

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD_DIR=${1:-"$ROOT/tests/build-ci"}
BASELINE="$ROOT/tests/known_failures.txt"
# Per-test ceiling. A hung test would otherwise stall the whole run, and a timeout is a failure:
# "it did not finish" is not a pass.
TIMEOUT=${TEST_TIMEOUT:-300}

# The pinned configuration. Anything optional that is OFF here is off for everyone reading the
# result, which is the point.
CMAKE_FLAGS="-DUSE_FLINT=ON -DUSE_MPFR=ON -DUSE_LAPACK=ON -DUSE_REGEX=ON -DUSE_FFTW=ON"

echo "=== configuration ==="
echo "build dir : $BUILD_DIR"
echo "cmake     : $CMAKE_FLAGS"
echo "timeout   : ${TIMEOUT}s per binary"
echo "commit    : $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
# shellcheck disable=SC2086
cmake $CMAKE_FLAGS "$ROOT/tests" > cmake.log 2>&1 || { echo "cmake failed:"; tail -20 cmake.log; exit 2; }

JOBS=$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 4 )
echo "=== building all test targets (-j$JOBS) ==="
if ! make -j"$JOBS" > build.log 2>&1; then
    echo "BUILD FAILED — the suite cannot be run:"
    grep -E "error:|Error [0-9]" build.log | head -20
    exit 2
fi

# A PORTABLE TIMEOUT. `timeout` is GNU coreutils and does not exist on macOS, where an unguarded
# `timeout 300 ./x` is `command not found` -- exit 127, which a pass/fail loop reads as a FAILED
# TEST. That mistake produced a "0 passed, 444 failed" run during this script's own development,
# and the number looked like a catastrophe rather than a missing binary. Homebrew coreutils
# provides `gtimeout`; with neither, tests run unbounded and the ceiling is documented as absent
# rather than silently turning every result red.
if command -v timeout > /dev/null 2>&1; then
    RUN="timeout $TIMEOUT"
elif command -v gtimeout > /dev/null 2>&1; then
    RUN="gtimeout $TIMEOUT"
else
    RUN=""
    echo "note: no timeout(1) or gtimeout(1) — tests run without a time limit"
fi

echo "=== running every test binary ==="
pass=0
fail=0
failed=""
for t in *_tests; do
    [ -x "$t" ] || continue
    # shellcheck disable=SC2086
    if $RUN "./$t" > "run-$t.log" 2>&1; then
        pass=$((pass + 1))
    else
        fail=$((fail + 1))
        failed="$failed $t"
    fi
done

# The baseline: first field of each non-comment line.
known=""
if [ -f "$BASELINE" ]; then
    known=$(grep -vE '^\s*(#|$)' "$BASELINE" | awk '{print $1}')
fi

unexpected=""
for t in $failed; do
    echo "$known" | grep -qx "$t" || unexpected="$unexpected $t"
done

fixed=""
for k in $known; do
    if [ -x "./$k" ]; then
        echo "$failed" | grep -qw "$k" || fixed="$fixed $k"
    fi
done

echo
echo "=== result ==="
echo "$pass passed, $fail failed, $(echo "$known" | grep -c . ) baselined"

if [ -n "$fixed" ]; then
    echo
    echo "BASELINED TESTS THAT NOW PASS — delete these lines from tests/known_failures.txt:"
    for f in $fixed; do echo "  $f"; done
fi

if [ -n "$unexpected" ]; then
    echo
    echo "UNEXPECTED FAILURES:"
    for t in $unexpected; do
        echo "--- $t"
        tail -5 "run-$t.log" | sed 's/^/    /'
    done
    echo
    echo "Each of these is either a real regression or a test that needs a line in"
    echo "tests/known_failures.txt with the reason. Do not add a line without one:"
    echo "a baseline whose entries are unexplained is a list of tests nobody will ever fix."
    exit 1
fi

echo "no unexpected failures"
exit 0

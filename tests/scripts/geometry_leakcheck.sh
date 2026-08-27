#!/bin/sh
# geometry_leakcheck.sh -- leak gate for the GEO-1 geometry builtins.
#
# The C test binary cannot detect a leak on its own, so this is the real gate.
# Drives every geometry head (exact and machine paths) under a leak checker and
# FAILS on any leaked byte.
#
#   macOS : leaks --atExit
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 -- an unavailable checker is not
#   a pass, so it says so rather than printing nothing (docs/adr/0004).
#
# Usage: bash tests/scripts/geometry_leakcheck.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "geometry_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t geoleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
Do[Area[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]], {200}];
Do[Area[Polygon[{{0,0},{1/3,0},{1/7,1/5}}]], {200}];
Do[Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]], {200}];
Do[Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]], {200}];
Do[Perimeter[Polygon[{{0,0},{3.,0},{3.,4.}}]], {200}];
Do[RegionCentroid[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]], {200}];
Do[RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {1/2,1/2}], {200}];
Do[ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}], {200}];
Do[ConvexHullRegion[{{0,0},{1,1},{2,2}}], {200}];
Do[ConvexHullRegion[{{1,2}}], {200}];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    OUT="$(leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null | grep 'total leaked bytes' | tail -1)"
    rm -f "$WORK"
    echo "geometry_leakcheck: $OUT"
    case "$OUT" in
        *"0 leaks for 0 total leaked bytes"*) echo "geometry_leakcheck: PASS"; exit 0 ;;
        "") echo "geometry_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "geometry_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>/tmp/geoleak.valgrind
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { echo "geometry_leakcheck: PASS (valgrind)"; exit 0; }
    echo "geometry_leakcheck: FAIL (valgrind) — see /tmp/geoleak.valgrind"; exit 1
else
    rm -f "$WORK"
    echo "geometry_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "geometry_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi

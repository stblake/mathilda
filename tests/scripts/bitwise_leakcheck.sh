#!/bin/sh
# bitwise_leakcheck.sh -- leak gate for the bit-level integer builtins.
#
# The C test binary cannot detect a leak on its own, so this is the real gate.
# Drives BitLength (machine-int, bignum, negative, INT64_MIN, listable, packed,
# and error paths) under a leak checker and FAILS on any leaked byte.
#
#   macOS : leaks --atExit
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 -- an unavailable checker is not
#   a pass, so it says so rather than printing nothing (docs/adr/0004).
#
# Usage: bash tests/scripts/bitwise_leakcheck.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "bitwise_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t bitleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
Do[BitLength[255], {200}];
Do[BitLength[-8], {200}];
Do[BitLength[0], {200}];
Do[BitLength[-9223372036854775808], {200}];
Do[BitLength[2^100], {200}];
Do[BitLength[2^100 - 1], {200}];
Do[BitLength[100!], {200}];
Do[BitLength[{0, 1, 2, 7, 8, 255, 256}], {200}];
Do[BitLength[{-1, -2, -8, -256}], {200}];
Do[BitLength[Range[64]], {200}];
Do[BitLength[NDArray[{-9223372036854775808, -1, 0, 1, 255}, DataType -> "int64"]], {200}];
Do[BitLength[3.5], {200}];
Do[BitLength[], {200}];
Do[BitLength[1, 2], {200}];
Do[BitLength[x], {200}];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    OUT="$(leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null | grep 'total leaked bytes' | tail -1)"
    rm -f "$WORK"
    echo "bitwise_leakcheck: $OUT"
    case "$OUT" in
        *"0 leaks for 0 total leaked bytes"*) echo "bitwise_leakcheck: PASS"; exit 0 ;;
        "") echo "bitwise_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "bitwise_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>/tmp/bitleak.valgrind
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { echo "bitwise_leakcheck: PASS (valgrind)"; exit 0; }
    echo "bitwise_leakcheck: FAIL (valgrind) — see /tmp/bitleak.valgrind"; exit 1
else
    rm -f "$WORK"
    echo "bitwise_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "bitwise_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi

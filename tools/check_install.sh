#!/bin/sh
# check_install.sh — prove the installed embedding API works (issue #88).
#
# Run via `make check-install` (which builds the executable and both libraries
# first). Runs `make install` into a scratch prefix, then, using ONLY the
# installed files (never the source tree):
#
#   0. CLI     — the installed executable starts and loads bootstrap modules;
#   1. shared  — a C host linked with -lmathilda evaluates through the shared
#                library with no $MATHILDA_HOME and from a CWD outside the repo,
#                so the library must locate share/mathilda/internal itself
#                (src/loadmodule.c: library-relative lookup);
#   2. static  — the same host linked from libmathilda.a with the flags the
#                installed mathilda.pc reports (`pkg-config --static`), pointed at
#                the module tree via `pkg-config --variable=internaldir`;
#   3. PIC     — libmathilda.a linked into a SHARED object, which is what a
#                Cython/Python extension does (fails without -fPIC on Linux);
#   4. exports — the shared library dynamically exports only mathilda_ffi_*.
#
# The probe is HeavisideTheta[0]: its value 0 comes from distributions.m, which
# init.m loads, so a library that silently failed to find its module tree
# answers HeavisideTheta[0] unevaluated and the check fails.
#
# Environment: CC (default gcc), MAKE (default make). POSIX sh; Linux + macOS.
set -eu

CC=${CC:-gcc}
MAKE=${MAKE:-make}
PKG_CONFIG=${PKG_CONFIG:-pkg-config}

case "$(uname -s)" in
    Darwin) OS=Darwin ;;
    Linux)  OS=Linux ;;
    *) echo "check-install: skipped (Linux/macOS only)"; exit 0 ;;
esac

ROOT=$(pwd)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/mathilda-check-install.XXXXXX")
trap 'rm -rf "$WORK"' EXIT INT TERM
# A real prefix rather than DESTDIR staging: pkg-config's SYSROOT would also
# re-root the dependency -L paths in Libs.private, which are host paths.
PREFIX="$WORK/inst"
LIBDIR="$PREFIX/lib"
fail=0

say()  { printf '  %-8s %s\n' "$1" "$2"; }
bad()  { say "$1" "FAIL: $2"; fail=1; }

echo "check-install: make install prefix=$PREFIX"
"$MAKE" -s install prefix="$PREFIX" >"$WORK/install.log" 2>&1 || {
    cat "$WORK/install.log"; echo "check-install: make install failed"; exit 1; }

for f in include/mathilda_ffi.h lib/libmathilda.a lib/pkgconfig/mathilda.pc \
         share/mathilda/internal/init.m; do
    [ -f "$PREFIX/$f" ] || bad files "missing $PREFIX/$f"
done

cat >"$WORK/startup.m" <<'EOF'
Print[2 + 2];
Print[HeavisideTheta[0]];
Print[$InstallationDirectory];
EOF
if out=$(cd "$WORK" && env -u MATHILDA_HOME \
         "$PREFIX/bin/Mathilda" -file "$WORK/startup.m" 2>&1); then
    expected=$(printf '4\n0\n%s' "$(cd "$PREFIX" && pwd -P)")
    if [ "$out" = "$expected" ]; then
        say CLI "ok (startup, bootstrap, installation directory)"
    else
        bad CLI "unexpected output: $out"
    fi
else
    bad CLI "startup failed: $out"
fi

# Host program: optional argv[1] = module-tree directory (static case).
cat >"$WORK/host.c" <<'EOF'
#include <stdio.h>
#include "mathilda_ffi.h"
int main(int argc, char** argv) {
    if (argc > 1) mathilda_ffi_set_home(argv[1]);
    char* r = mathilda_ffi_eval("HeavisideTheta[0]");
    printf("%s\n", r);
    mathilda_ffi_free(r);
    return 0;
}
EOF

# pkg-config against the installed .pc only.
pc() {
    PKG_CONFIG_PATH="$LIBDIR/pkgconfig" PKG_CONFIG_LIBDIR="$LIBDIR/pkgconfig" \
    "$PKG_CONFIG" "$@" mathilda
}

expect_closed_form() {   # $1 = label, $2 = output
    case "$2" in
        *HeavisideTheta*|"") bad "$1" "module tree not loaded (got: $2)" ;;
        0)            say "$1" "ok ($2)" ;;
        *)            bad "$1" "unexpected output: $2" ;;
    esac
}

# 1. shared, no MATHILDA_HOME, CWD outside the repo.
if "$CC" -I"$PREFIX/include" "$WORK/host.c" -L"$LIBDIR" -lmathilda \
        -o "$WORK/host_shared" 2>"$WORK/shared.log"; then
    out=$(cd "$WORK" && env -u MATHILDA_HOME \
          LD_LIBRARY_PATH="$LIBDIR" DYLD_LIBRARY_PATH="$LIBDIR" \
          ./host_shared 2>&1) || true
    expect_closed_form shared "$out"
else
    cat "$WORK/shared.log"; bad shared "link against -lmathilda failed"
fi

# 2. static, flags from the installed mathilda.pc.
if command -v "$PKG_CONFIG" >/dev/null 2>&1; then
    cflags=$(pc --cflags)
    # Force the ARCHIVE (the shared library sits beside it and would win).
    libs=$(pc --static --libs | sed "s#-lmathilda#$LIBDIR/libmathilda.a#")
    home=$(pc --variable=internaldir)
    # shellcheck disable=SC2086
    if "$CC" $cflags "$WORK/host.c" $libs -o "$WORK/host_static" \
            2>"$WORK/static.log"; then
        out=$(cd "$WORK" && env -u MATHILDA_HOME ./host_static "$home" 2>&1) || true
        expect_closed_form static "$out"
    else
        cat "$WORK/static.log"; bad static "link from pkg-config --static failed"
    fi

    # 3. libmathilda.a inside a shared object (the Cython extension case).
    # shellcheck disable=SC2086
    if "$CC" -fPIC -c $cflags "$WORK/host.c" -o "$WORK/ext.o" \
       && "$CC" -shared "$WORK/ext.o" $libs -o "$WORK/ext.so" \
            2>"$WORK/pic.log"; then
        say PIC "ok (libmathilda.a links into a shared object)"
    else
        cat "$WORK/pic.log"; bad PIC "libmathilda.a cannot be linked into a shared object"
    fi
else
    bad static "$PKG_CONFIG not found"
fi

# 4. dynamic exports are exactly the mathilda_ffi_* API.
if [ "$OS" = Darwin ]; then
    # GCC on Darwin links its runtime shims (-lemutls_w, -lheapt_w) into every
    # image as WEAK exports, deliberately, so the process coalesces them to one
    # copy; they are the toolchain's, not ours, and an export list cannot drop
    # them.
    syms=$(nm -gU "$LIBDIR/libmathilda.0.dylib" | awk '{print $3}' | sed 's/^_//' |
           grep -v -E '^_(_emutls_|_gcc_nested_func_ptr_)')
else
    syms=$(nm -D --defined-only "$LIBDIR/libmathilda.so.0" | awk '{print $3}' |
           grep -v -E '^(_init|_fini|_edata|_end|__bss_start)$')
fi
extra=$(printf '%s\n' "$syms" | grep -v '^mathilda_ffi_' | grep -v '^$' || true)
napi=$(printf '%s\n' "$syms" | grep -c '^mathilda_ffi_' || true)
if [ -n "$extra" ]; then
    bad exports "non-API symbols exported: $(printf '%s\n' "$extra" | head -5 | tr '\n' ' ')"
elif [ "$napi" -eq 0 ]; then
    bad exports "no mathilda_ffi_* symbols exported"
else
    say exports "ok ($napi mathilda_ffi_* symbols, nothing else)"
fi

if [ "$fail" -ne 0 ]; then
    echo "check-install: FAILED"; exit 1
fi
echo "check-install: all checks passed"

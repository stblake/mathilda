#!/usr/bin/env python3
"""Guard against constructs that compile on macOS and break the Linux build.

`src/` is compiled with `-std=c99`. Under that flag glibc defines
`__STRICT_ANSI__` and hides everything outside ISO C99 -- `M_PI`, `jn`, `yn`,
`strdup`, `fileno`, `clock_gettime`, ... Darwin's headers expose those symbols
regardless of the standard, so an unguarded use compiles cleanly on macOS and
reaches users as a Linux build failure. That has now happened twice: issue #36
(`M_PI`) and issue #37 (`jn` / `yn`, an *error* rather than a warning under
GCC 14, which promoted implicit-function-declaration to an error).

There are two sanctioned ways to use such a symbol, and this script enforces
that one of them is present.

**Constants** get a per-file fallback right after `#include <math.h>` (see
`src/trig.c`, `src/numeric.c`)::

    #include <math.h>
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif

**Functions** cannot be re-implemented by a macro, so they instead request the
namespace they live in with a feature-test macro, placed *before any include*
(see `src/core.c`, `src/repl.c`, `src/loadmodule.c`, `src/ndkernels.c`)::

    #ifndef _XOPEN_SOURCE
    #define _XOPEN_SOURCE 600
    #endif
    #include <math.h>

Ordering matters: a feature-test macro defined after the first `#include` has
no effect, because the header it was meant to unlock has already been parsed
with the wrong namespace. That case is reported separately.

The third check covers a different platform difference with the same shape.
`int64_t` is `long long` on Darwin and `long` under glibc -- the same width,
but not the same type. `src/checked_int.h` therefore ships two families of
overflow-checked helpers, and mixing them is invisible on macOS and a hard
error under GCC 14 on Linux, which is issue #40: `ci_powi` was handed an
`int64_t*` out-pointer. Only `src/compile/` may use the `long long` family --
its Slot register is a `long long`; everything else holds `int64_t` buffers and
must use the `_i64` spelling.

Note that a compile flag cannot catch this one on the development machine: on
Darwin the two types coincide, so there is nothing for the compiler to object
to. A source-level rule is the only thing that fails locally.

Run via `make check-c99`; it exits nonzero and names every offending file.

Scope
-----
Constants are checked across `src/` and `tests/`. Functions and the
checked-integer families are checked in `src/` only: the test suite is built by
CMake at the compiler's default standard (gnu17), where the whole POSIX surface
is visible, so flagging it would be noise rather than signal.

Extending
---------
`FUNCTIONS` is a curated table, not the entire POSIX standard -- it covers the
symbols this project plausibly reaches for. Adding an entry is one line: map
the name to the feature-test macros (and minimum values) that expose it in
glibc. `_GNU_SOURCE` is accepted for every symbol and need not be listed.

`WIDE_INT_FAMILY` is the `long long` half of `src/checked_int.h`; a new helper
added there in both spellings gets one more name in that tuple.
"""

import os
import re
import subprocess
import sys

# ---------------------------------------------------------------------------
# Check 1: POSIX <math.h> constants
# ---------------------------------------------------------------------------

# The full set of POSIX <math.h> constants. All carry the same portability
# hazard, so all are checked -- not just the ones that have bitten us so far.
CONSTANTS = (
    "M_E", "M_LOG2E", "M_LOG10E", "M_LN2", "M_LN10",
    "M_PI", "M_PI_2", "M_PI_4", "M_1_PI", "M_2_PI",
    "M_2_SQRTPI", "M_SQRT2", "M_SQRT1_2",
)

# Reference values, so the error message can be pasted straight into the fix.
VALUES = {
    "M_E":        "2.71828182845904523536",
    "M_LOG2E":    "1.44269504088896340736",
    "M_LOG10E":   "0.434294481903251827651",
    "M_LN2":      "0.693147180559945309417",
    "M_LN10":     "2.30258509299404568402",
    "M_PI":       "3.14159265358979323846",
    "M_PI_2":     "1.57079632679489661923",
    "M_PI_4":     "0.785398163397448309616",
    "M_1_PI":     "0.318309886183790671538",
    "M_2_PI":     "0.636619772367581343076",
    "M_2_SQRTPI": "1.12837916709551257390",
    "M_SQRT2":    "1.41421356237309504880",
    "M_SQRT1_2":  "0.707106781186547524401",
}

USE_RE = re.compile(r"\b(" + "|".join(CONSTANTS) + r")\b")
GUARD_RE = re.compile(r"^\s*#\s*(?:ifndef|define)\s+(M_\w+)")

# ---------------------------------------------------------------------------
# Check 2: POSIX-only libc functions
# ---------------------------------------------------------------------------

# Each symbol maps to the feature-test macros that expose it in glibc, as
# (macro, minimum value) alternatives -- satisfying ANY one of them is enough.
# A minimum of 0 means "defined at all, value irrelevant". `_GNU_SOURCE` turns
# on everything and is accepted implicitly for every symbol.
_XSI = (("_XOPEN_SOURCE", 500), ("_DEFAULT_SOURCE", 0), ("_BSD_SOURCE", 0))
_P1 = (("_POSIX_C_SOURCE", 1), ("_XOPEN_SOURCE", 0), ("_DEFAULT_SOURCE", 0))
_P199309 = (("_POSIX_C_SOURCE", 199309), ("_XOPEN_SOURCE", 500))
_P199506 = (("_POSIX_C_SOURCE", 199506), ("_XOPEN_SOURCE", 500))
_P200112 = (("_POSIX_C_SOURCE", 200112), ("_XOPEN_SOURCE", 600))
_P200809 = (("_POSIX_C_SOURCE", 200809), ("_XOPEN_SOURCE", 700),
            ("_DEFAULT_SOURCE", 0))
_GNU_ONLY = ()

FUNCTIONS = {
    # <math.h> X/Open extensions -- the issue #37 family. `gamma` and `finite`
    # are legacy XSI names, easy to reach for by accident.
    "j0": _XSI, "j1": _XSI, "jn": _XSI,
    "y0": _XSI, "y1": _XSI, "yn": _XSI,
    "scalb": _XSI, "drem": _XSI, "significand": _XSI,
    "finite": _XSI, "gamma": _XSI, "lgamma_r": _XSI,
    # <math.h> GNU extensions.
    "sincos": _GNU_ONLY, "exp10": _GNU_ONLY, "pow10": _GNU_ONLY,
    # <string.h> / <strings.h>
    "strdup": _P200809, "strndup": _P200809,
    "strtok_r": _P199506,
    "strcasecmp": _P200809, "strncasecmp": _P200809,
    "strsep": _GNU_ONLY, "strcasestr": _GNU_ONLY, "memmem": _GNU_ONLY,
    # <stdio.h>
    "fileno": _P1, "fdopen": _P1, "popen": _P200112, "pclose": _P200112,
    "getline": _P200809, "getdelim": _P200809,
    "asprintf": _GNU_ONLY, "vasprintf": _GNU_ONLY,
    # <stdlib.h>
    "setenv": _P200112, "unsetenv": _P200112, "putenv": _XSI,
    "realpath": _P200112, "mkstemp": _P200112, "mkdtemp": _P200809,
    "random": _XSI, "srandom": _XSI, "drand48": _XSI, "lrand48": _XSI,
    # <unistd.h> -- glibc hides essentially all of it without _POSIX_C_SOURCE.
    "access": _P1, "isatty": _P1, "unlink": _P1, "rmdir": _P1,
    "chdir": _P1, "getcwd": _P1, "dup": _P1, "dup2": _P1, "pipe": _P1,
    "fork": _P1, "execvp": _P1, "getpid": _P1, "ftruncate": _P1,
    "sleep": _P1, "usleep": _XSI,
    # <time.h> / <sys/time.h> / <sys/stat.h> / <dirent.h>
    "clock_gettime": _P199309, "nanosleep": _P199309,
    "localtime_r": _P199506, "gmtime_r": _P199506,
    "gettimeofday": _P1, "mkdir": _P1,
    "opendir": _P1, "readdir": _P1, "closedir": _P1,
}

# A call, but not `p->jn(...)` / `s.jn(...)` -- a struct member of that name is
# the project's own, not libc's.
CALL_RE = {
    sym: re.compile(r"(?<![\w.])(?<!->)\b" + sym + r"\s*\(")
    for sym in FUNCTIONS
}

# ---------------------------------------------------------------------------
# Check 3: the `long long` checked-integer family outside src/compile/
# ---------------------------------------------------------------------------

# The `long long` half of src/checked_int.h. Each has an `_i64` twin with
# identical semantics on int64_t; see this module's docstring for why the two
# are not interchangeable off Darwin.
WIDE_INT_FAMILY = ("ci_add", "ci_sub", "ci_mul", "ci_neg", "ci_abs", "ci_powi")

# The compiler's Slot register genuinely is a `long long`, so src/compile/ is
# the one legitimate caller. checked_int.h itself defines them.
WIDE_INT_HOME = "src/compile/"
WIDE_INT_DEFN = "src/checked_int.h"

# `ci_add(` but not `ci_add_i64(` -- the `(` must follow the name directly.
WIDE_INT_RE = {
    sym: re.compile(r"(?<![\w.])(?<!->)\b" + sym + r"\s*\(")
    for sym in WIDE_INT_FAMILY
}

# `#define _POSIX_C_SOURCE 200809L` -- the trailing L is optional, and a bare
# `#define _GNU_SOURCE` carries no value at all.
FTM_RE = re.compile(r"^\s*#\s*define\s+(_[A-Z0-9_]+)(?:\s+(\d+)L?)?\s*$")
INCLUDE_RE = re.compile(r"^\s*#\s*include\b")
FTM_NAMES = ("_GNU_SOURCE", "_DEFAULT_SOURCE", "_BSD_SOURCE",
             "_POSIX_C_SOURCE", "_XOPEN_SOURCE")


def strip_comments_and_strings(text):
    """Blank out comments and string/char literals.

    A symbol named only in prose ("M_PI is POSIX, not C99", "the candidate set
    is finite (bounded by deg h)") is not a usage, and a checker that cannot
    tell the difference gets ignored. Replacing with spaces rather than
    deleting keeps line numbers intact for reporting.
    """
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        two = text[i:i + 2]
        if two == "/*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        elif two == "//":
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'":
            quote, j = c, i + 1
            while j < n and text[j] != quote:
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def source_files(repo_root, dirs=("src", "tests")):
    """Git-tracked C sources under `dirs`, excluding vendored code."""
    try:
        tracked = subprocess.check_output(
            ["git", "-C", repo_root, "ls-files"] + list(dirs),
            text=True,
        ).split("\n")
    except (subprocess.CalledProcessError, OSError) as exc:
        sys.stderr.write("check_c99_portability: cannot list git files: %s\n"
                         % exc)
        sys.exit(2)

    ext = (".c", ".h", ".inc")
    return [
        f for f in tracked
        if f.endswith(ext) and not f.startswith("src/external/")
    ]


def read(repo_root, path):
    try:
        with open(os.path.join(repo_root, path), encoding="utf-8",
                  errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


def scan_constants(text):
    """Return the sorted constants used in `text` that it does not guard."""
    used, guarded = set(), set()
    for line in text.splitlines():
        guard = GUARD_RE.match(line)
        if guard:
            # `#ifndef M_PI` / `#define M_PI ...` -- the guard itself, not a use.
            guarded.add(guard.group(1))
            continue
        used.update(USE_RE.findall(line))
    return sorted(used - guarded)


def feature_test_macros(text):
    """Return (effective, late) feature-test macros defined by this file.

    `effective` maps macro -> int value for those defined before the first
    `#include`; `late` is the set defined after it, which the preprocessor has
    already ignored by the time the system headers were read.
    """
    effective, late = {}, set()
    seen_include = False
    for line in text.splitlines():
        if INCLUDE_RE.match(line):
            seen_include = True
            continue
        m = FTM_RE.match(line)
        if not m or m.group(1) not in FTM_NAMES:
            continue
        name, value = m.group(1), int(m.group(2) or 0)
        if seen_include:
            late.add(name)
        else:
            effective[name] = max(effective.get(name, 0), value)
    return effective, late


def declares_locally(text, sym):
    """True if this file declares or defines its own `sym`.

    Guards against a project function that happens to share a POSIX name. The
    `=` test keeps a genuine call in an initialiser (`double v = jn(n, x);`)
    from reading as a declaration just because a type precedes it.
    """
    pat = re.compile(
        r"^[ \t]*(?:static[ \t]+|extern[ \t]+|inline[ \t]+)*"
        r"[A-Za-z_][A-Za-z0-9_]*(?:[ \t]+[A-Za-z_][A-Za-z0-9_]*)*"
        r"[ \t\*]+" + sym + r"[ \t]*\(", re.M)
    for m in pat.finditer(text):
        end = text.find("\n", m.start())
        line = text[m.start():len(text) if end < 0 else end]
        if "=" not in line.split(sym, 1)[0]:
            return True
    return False


def scan_functions(text):
    """Return (unguarded, late) POSIX functions used by `text`.

    `unguarded` is a list of (symbol, alternatives); `late` names symbols whose
    feature-test macro exists but sits after the first `#include`.
    """
    effective, late_macros = feature_test_macros(text)
    if "_GNU_SOURCE" in effective:
        return [], []          # turns on everything glibc has

    unguarded, late = [], []
    for sym in sorted(FUNCTIONS):
        alternatives = FUNCTIONS[sym]
        if not CALL_RE[sym].search(text) or declares_locally(text, sym):
            continue
        if any(effective.get(macro, -1) >= minimum
               for macro, minimum in alternatives):
            continue
        accepted = set(m for m, _ in alternatives) | {"_GNU_SOURCE"}
        if late_macros & accepted:
            late.append(sym)
        else:
            unguarded.append((sym, alternatives))
    return unguarded, late


def scan_wide_int_family(text):
    """Return the sorted `long long` checked-int helpers called in `text`."""
    return sorted(sym for sym in WIDE_INT_FAMILY if WIDE_INT_RE[sym].search(text))


def suggest(alternatives):
    """The feature-test macro we recommend for a symbol, as source lines."""
    if not alternatives:
        return "#ifndef _GNU_SOURCE\n      #define _GNU_SOURCE\n      #endif"
    macro, minimum = alternatives[0]
    if macro == "_POSIX_C_SOURCE":
        value = "200809L" if minimum <= 200809 else "%dL" % minimum
    elif macro == "_XOPEN_SOURCE":
        value = "600" if minimum <= 600 else str(minimum)
    else:
        value = ""
    return ("#ifndef %s\n      #define %s%s\n      #endif"
            % (macro, macro, (" " + value) if value else ""))


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    const_offenders, func_offenders, late_offenders = [], [], []
    wide_offenders = []

    for path in source_files(repo_root):
        raw = read(repo_root, path)
        if raw is None:
            continue
        text = strip_comments_and_strings(raw)

        missing = scan_constants(text)
        if missing:
            const_offenders.append((path, missing))

        # Functions: `src/` only. The CMake test build uses the compiler's
        # default standard, where the POSIX surface is visible anyway.
        if path.startswith("src/"):
            unguarded, late = scan_functions(text)
            if unguarded:
                func_offenders.append((path, unguarded))
            if late:
                late_offenders.append((path, late))

            if not path.startswith(WIDE_INT_HOME) and path != WIDE_INT_DEFN:
                wide = scan_wide_int_family(text)
                if wide:
                    wide_offenders.append((path, wide))

    if not (const_offenders or func_offenders or late_offenders
            or wide_offenders):
        return 0

    if const_offenders:
        sys.stderr.write(
            "error: POSIX <math.h> constants used without a C99 fallback.\n"
            "       glibc hides these under -std=c99; the build fails on Linux\n"
            "       even though macOS compiles it. Add, after #include <math.h>:\n\n"
        )
        for path, missing in const_offenders:
            sys.stderr.write("  %s\n" % path)
            for name in missing:
                sys.stderr.write(
                    "      #ifndef %s\n      #define %s %s\n      #endif\n"
                    % (name, name, VALUES[name]))
        sys.stderr.write("\n")

    if func_offenders:
        sys.stderr.write(
            "error: POSIX-only functions used without a feature-test macro.\n"
            "       glibc hides these under -std=c99, so GCC reports an implicit\n"
            "       declaration -- an ERROR since GCC 14 -- on Linux while macOS\n"
            "       compiles clean. Add, BEFORE the first #include:\n\n"
        )
        for path, syms in func_offenders:
            sys.stderr.write("  %s  (%s)\n"
                             % (path, ", ".join(s for s, _ in syms)))
            for text in sorted(set(suggest(a) for _, a in syms)):
                sys.stderr.write("      %s\n" % text)
        sys.stderr.write("\n")

    if late_offenders:
        sys.stderr.write(
            "error: feature-test macro defined AFTER the first #include, where\n"
            "       the preprocessor has already read the header it was meant\n"
            "       to unlock. Move the #define to the top of the file:\n\n"
        )
        for path, syms in late_offenders:
            sys.stderr.write("  %s  (%s)\n" % (path, ", ".join(syms)))
        sys.stderr.write("\n")

    if wide_offenders:
        sys.stderr.write(
            "error: the `long long` half of src/checked_int.h used outside\n"
            "       src/compile/. Only the compiler's Slot register is a\n"
            "       `long long`; every other caller holds int64_t buffers and\n"
            "       must use the _i64 spelling. Darwin typedefs int64_t to\n"
            "       `long long`, so a mix-up compiles clean here and is an\n"
            "       ERROR under GCC 14 on glibc, where it is `long`\n"
            "       (issue #40). Rename the call:\n\n"
        )
        for path, syms in wide_offenders:
            sys.stderr.write("  %s\n" % path)
            for sym in syms:
                sys.stderr.write("      %s(...)   ->   %s_i64(...)\n"
                                 % (sym, sym))
        sys.stderr.write("\n")

    total = len(set(p for p, _ in
                    const_offenders + func_offenders + late_offenders
                    + wide_offenders))
    sys.stderr.write("%d file(s) affected.\n" % total)
    return 1


if __name__ == "__main__":
    sys.exit(main())

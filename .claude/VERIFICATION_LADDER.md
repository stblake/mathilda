# Verification ladder configuration

One file, read by `skills/verification-ladder/scripts/ladder.py`. Change a value here instead
of editing that script.

Written by `/setup-kit` on 2026-08-22. DETECT (`kit-setup/scripts/detect_ladder.py`) found
**zero signal for all four phases** against this repo: `CMakeLists.txt` lives at
`tests/CMakeLists.txt`, not repo root (the detector only checks the root), and there is no
`c-cmake`/`make` toolchain entry that recognizes a hand-written top-level `makefile` at all —
so a `make ...`-only CI (`.github/workflows/build.yml`) produced no CI-marker hits either.
Every value below comes from `SPEC.md` §9 and `.github/workflows/build.yml`, confirmed with
the human directly, not from DETECT.

## Static analysis

```
static = make check-c99 && make check-packed-aware
```

`make check-c99` (`tools/check_c99_portability.py`) flags POSIX-only symbols glibc hides
under `-std=c99` (missing `#ifndef` guards on `<math.h>` constants, missing feature-test
macros on POSIX functions) — this repo's actual lint tier, and the one CI runs first
(`.github/workflows/build.yml:73`). `make check-packed-aware` is a static structural audit:
diffs the NDArray dispatch sites in source against `src/pack.c`'s `AWARE` list.

## Type checking

```
typecheck = make -j$(nproc)
```

C99 has no separate type-check phase from compilation — `SPEC.md`'s own words: "for
C/C++/Go/Java, this is often just 'does it build.'" `gcc -std=c99 -Wall -Wextra`
(`-Werror` on several diagnostics GCC 14 promotes) is the closest equivalent.

## Unit tests

```
unit = cd tests && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc) && for t in *_tests; do ./$t || exit 1; done
```

The CMake-built `test_*.c` suite under `tests/` (SPEC.md §9).

## Integration tests

```
integration = make check-nd-surfaces
```

Runs the built `./Mathilda` binary over probe expressions across three representations
(plain `List`, packed `List`, visible `NDArray`) and requires them to agree — components
(evaluator, pack transparency gate, ND kernels) exercised together, with real execution, not
isolated unit logic. `valgrind --leak-check=full ./Mathilda` (SPEC.md §4) is real practice
here but was left out of this phase: it needs a script argument with no single canonical
choice, so scripting it as the ladder's integration command would be a guess, not a
confirmed answer.

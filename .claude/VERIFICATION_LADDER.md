# Verification ladder configuration — mathilda (C99, GCC-only)

Filled in manually for this repo's real toolchain (per VERIFICATION_LADDER.example.md;
/setup-kit's DETECT pass would find the same commands — makefile:469 `check-c99`,
README.md:230-237 CMake test build). GCC is mandatory: both the makefile (clang → hard
error) and tests/CMakeLists.txt:20-25 (FATAL_ERROR on Apple clang) enforce it.

## Static analysis

```
static = python3 tools/check_c99_portability.py
```

The repo's own C99/POSIX portability gate (`make check-c99`) — catches the
implicit-declaration class that is a warning on macOS and a wrong answer on Linux.

## Type checking

```
typecheck = make -j8
```

For this C toolchain, type checking IS the build: gcc-16 with
-Werror=implicit-function-declaration/-incompatible-pointer-types/-int-conversion/
-implicit-int promoted (makefile:44-67).

## Unit tests

```
unit = cmake -S tests -B tests/build -DCMAKE_C_COMPILER=gcc-16 && cmake --build tests/build -j8 && ctest --test-dir tests/build --output-on-failure
```

CMake-only test suite; ctest runs every add_test-registered binary.

## Integration tests

```
integration = not-configured
```

Honest gap: the repo has no separate integration tier (benchmarks/ exist but compare
against golden outputs out-of-band). Multi-statement `./Mathilda -file` scripts
mis-parse (observed 2026-08-27), so scripted end-to-end runs are not a trustworthy
rung yet.

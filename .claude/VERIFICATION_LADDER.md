# Verification ladder configuration — mathilda (C99, GCC-only)

Read by the ais kit's `skills/verification-ladder/scripts/ladder.py`. Format: fenced
`phase = command` lines. Written for GEO-1 (2026-08-27); commands are the repo's own
real toolchain, mirroring CI order (.github/workflows/build.yml: check-c99 ->
check-packed-aware -> build -> stdin smoke test).

## Static analysis

C99/POSIX portability sweep (tools/check_c99_portability.py) — the class of bug that
compiles clean on macOS and breaks on glibc (repo issues #36/#37/#40):

```
static = make check-c99
```

## Type checking

C has no separate typecheck rung; compilation with the repo's promoted
-Werror set (implicit-function-declaration, incompatible-pointer-types,
int-conversion, implicit-int) IS the type gate, per the kit example's own
note for C toolchains. SDKROOT must be exported for Homebrew GCC on macOS:

```
typecheck = SDKROOT=$(xcrun --show-sdk-path) make -j8
```

## Unit tests

The C test binaries under tests/, built and run via ctest. NOTE the baseline,
measured on unmodified main (commit 15b088da, 2026-08-27): `bench_pack` and
`interp_tests` fail there and are inherited, NOT introduced. A third test,
`primenu_tests`, fails only under full parallel load (its ECM factor search
exhausts its budget: "FactorInteger::nofac ... no factor was found within the
search bounds", 7 prime factors instead of 8) and PASSES standalone in ~8s —
environment-dependent, not deterministic. So a full-suite run on this repo
shows 2 or 3 reds depending on machine load, and none of them are geometry.
Use ladder.py --baseline to attribute rungs honestly:

```
unit = cd tests/build && ctest --output-on-failure
```

## Integration tests

Full REPL pipeline (parse -> evaluate -> print) over the GEO-1 acceptance rows.
The grep is load-bearing: Quit[1] does not propagate an exit code through -file
mode (fault-injection receipt 2026-08-27), so the marker line is the verdict:

```
integration = ./Mathilda -file tests/scripts/geometry_e2e.m | grep -q "geometry_e2e: ALL PASS"
```

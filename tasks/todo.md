# Task: Fix Simplify hang + add Simplify TimeConstraint option

Reported: `Simplify[(x+Sqrt[u]) u^(-1/2)]` hangs. Root cause = unguarded
`while(true)` in `pseudo_rem` (poly GCD) reached via Simplify's Factor-polish.

## Implementation
- [x] Fix 1: guard `pseudo_rem` loop (degree-monotonicity) — `src/poly/poly.c`
- [x] Fix 2: `has_compound_radicand` + narrow the Factor-polish gate — `src/simp/simp_builtins.c`
- [x] Feature: `simp_mono_seconds()` helper (+ POSIX guard) — `src/simp/simp_util.c` / `simp_internal.h`
- [x] Feature: dynamically-scoped `time_budget` (save/restore, not param threading) + deadline checks in `simp_search`
- [x] Feature: parse `TimeConstraint` option in `builtin_simplify` (before positional fallback)
- [x] Feature: register `TimeConstraint -> Infinity` default — `src/options_builtin.c`
- [x] Docstring: mention `TimeConstraint` on Simplify — `src/info.c`
- [x] Tests: new `tests/test_simplify_hang.c` (24 termination battery + option tests) + CMake register (ctest)
- [x] Docs: `docs/spec/builtins/simplification.md` + weekly changelog `docs/spec/changelog/2026-09-07.md`

## Verify
- [x] `make -j` clean (0 warnings), `make check-c99` green
- [x] reported case returns `1 + x/Sqrt[u]` promptly; `Factor[1 + x/Sqrt[a]]` terminates
- [x] `simplify_hang_tests` PASS; existing simplify/fullsimplify/radical/trigrat/trig_canon/logexp/invtrig/factorial + facpoly/factorlist/factor_baseline/groebner/intrat all PASS
- [x] valgrind: no Mathilda frames in leak stacks (only macOS libobjc/dyld baseline noise) — bail path leak-free

## Review

**Root cause was NOT what the report suggested.** The user hypothesised a
*cycle in Simplify*; tracing (`$SimplifyDebug`, lldb, `sample`) showed the
search finishes instantly and the hang is an unguarded `while(true)` in
`pseudo_rem` (multivariate GCD) reached via Simplify's post-search `Factor`
polish. A Simplify-layer cycle/deadline guard could not have caught it (control
never returns to the Simplify loop). Confirmed scope with the user before
coding.

**Delivered:** (1) `pseudo_rem` degree-monotonicity guard — fixes the hang for
every Factor/GCD caller; (2) narrowed the Simplify Factor-polish gate
(`has_compound_radicand`) so it skips bare-symbol radicals; (3) new synchronous,
per-sub-expression, leak-free `TimeConstraint` option on Simplify (the user's
explicit feature request), implemented via a dynamically-scoped budget +
local deadline in `simp_search` rather than parameter threading (far less churn,
re-entrancy-safe). All output-neutral except the previously-hanging class.

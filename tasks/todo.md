# Task: overflow-robust Schwartz–Zippel sampler in PossibleZeroQ / zero_test

Fixes POSSIBLE_ZEROQ_IMPROVEMENTS.md #2 (wide dynamic range): overflow sample points
made the sampler abort to UNKNOWN→True for nowhere-zero functions. Re-draw/skip overflow
points instead. Scope: overflow-robust sampler only (item #1 residual out of scope).

## Implementation (all src/zero_test.c unless noted)
- [ ] Add `ZT_OVERFLOW_SHELL_BITS[]` ladder constant {6,4,2,1,0}
- [ ] `evaluate_rung`: add `bool* out_overflow`; ok&&!isfinite(mag) → flag + UNKNOWN; !ok → plain UNKNOWN
- [ ] `decide_numeric`: add `bool* out_overflow` (+ fwd decl :159); propagate from rung-0 call, NULL to higher rungs
- [ ] `screen_point`: add `bool* out_overflow`; propagate
- [ ] `sample_random_value(int num_bits)`; `sample_channel(...,int num_bits)`; `sample_random_value_spec(...,int num_bits)` (two-sided range branch ignores num_bits)
- [ ] `sz_trial`: add `int num_bits`, `bool* out_overflow`; thread through draw + eval calls
- [ ] `decide_schwartz_zippel_core`: shell-ladder retry on overflow-UNKNOWN in screen + confirm loops; skip if all shells overflow; screen returns UNKNOWN if zero informative points; confirm retains TRUE
- [ ] `zt_decide_core` Stage 2 call (:1716): pass NULL

## Tests
- [ ] Group 18 in tests/test_zero_test.c: E^(-10t)+E^(-100t), Gamma[x^2]+1, Gamma[x+1]-x Gamma[x]+E^(x^2) → False (timed + stable); overflow-prone identities stay True; matched non-identities stay False
- [ ] Run zero_test_tests, trigexp_zero_tests, possiblezeroq_* (grep FAIL), + downstream: refine, comparisons/Equal, dsolve, solve, integrate, limit, interp, radical_simplify, nullspace
- [ ] Inspect every True→False flip

## Docs / release
- [ ] POSSIBLE_ZEROQ_IMPROVEMENTS.md #2 status
- [ ] docs/spec/builtins/expression-information.md PossibleZeroQ note
- [ ] docs/spec/changelog/2026-09-21.md entry
- [ ] src/version.h → 0.206; tag v0.206
- [ ] valgrind spot-check; rebuild code-review graph

## Review (completed 2026-09-26, v0.208)

**Outcome.** Fixed the wrong/flaky `True` in `PossibleZeroQ` for wide-dynamic-range
exponential sums (POSSIBLE_ZEROQ_IMPROVEMENTS.md #2). Root cause: the Schwartz–Zippel
sampler aborted the whole test on the first IEEE-overflow sample point (`E^(large) → ±Inf`),
collapsing `UNKNOWN → True` for nowhere-zero functions, with the verdict decided by hash-seeded
draw order. Fix (`src/zero_test.c`): `evaluate_rung` flags an overflow (`ok && !isfinite(mag)`)
distinctly from a symbolic residue; `sz_trial_shelled` re-draws an overflow point through a
shrinking magnitude-shell ladder (`2^6→2^4→2^2→2^1→2^0`, `|value|>=1` floor kept), skipping only
if every shell overflows; a whole screen of overflow-only points returns an honest `UNKNOWN`.

**Scope narrowed during implementation.** An initial broader version (re-drawing symbolic
residues and sampled poles too) regressed `test_battery_weierstrass_cosh_product_roundtrip` by
exposing a *separate* pre-existing deep-cancellation false-negative in the screen phase. Narrowed
to the IEEE-overflow class only, which is byte-for-byte behaviour-preserving for non-overflowing
inputs — so no regressions. `Gamma[x^2]+1` (special-function magnitude-decline residue) and the
deep-cancellation screen weakness are logged as follow-ups #4/#5 in POSSIBLE_ZEROQ_IMPROVEMENTS.md.

**Verification.** All green: `zero_test_tests` (incl. new Group 18) + `trigexp_zero_tests`,
`possiblezeroq_{expcombine,assumptions,stress}`, `refine_tests`, `interp`/`nullspace`/
`solve_radicals_reals`/`limit`, `integrate_goursat`/`integrate_jeffrey`, DSolve corpus §2.2.7
(wide-spectrum systems, within baseline). Determinism confirmed (stable ×8). Clean GCC build
(`-std=c99 -Wall -Wextra`, no warnings). No new heap allocations (leak risk nil). `$VersionNumber
= 0.208`.

**Not committed yet** — awaiting go-ahead to commit to main + tag `v0.208`.

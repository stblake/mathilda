# Task: Implement `ListInterpolation` in `src/interp.c`

## Plan
Front-end that synthesizes abscissae from a raw rectangular value array, builds an
`Interpolation`-style `{{coord, val}, ...}` table, and delegates to the existing
`builtin_interpolation_impl`. Full reuse of the InterpolatingFunction engine.

## Items
- [x] Intern `SYM_ListInterpolation` (sym_names.h/.c, 3 sites)
- [x] Implement `builtin_listinterpolation_impl` + NDArray wrapper `builtin_listinterpolation` in interp.c
  - [x] Arg split (array / domain / option rules)
  - [x] `listinterp_shape` dimensionality + rectangularity (validated in emit walk)
  - [x] Abscissa synthesis (default ints / interval exact-rational / explicit positions)
  - [x] Row-major table walk (scalar + array-valued values)
  - [x] Synthetic `Interpolation[...]` delegation + free
- [x] Register in `interp_init`: builtin + ATTR_PROTECTED + Options
- [x] Add `"ListInterpolation"` to AWARE[] in pack.c
- [x] Docstring in info.c
- [x] Docs: functional-programming.md section + changelog 2026-08-31.md
- [x] Tests in tests/test_interp.c (extensive; leak-clean) — written; build/run pending
- [x] Build main; REPL spot-checks all match WL (2.4375, 3.875, 2-D, NDArray, argt)
- [x] Audits: check-c99 clean, check-packed-aware OK (AWARE=199)
- [ ] Run interp_tests (background build in progress), leaks, check-nd-surfaces

## Review

**Done.** `ListInterpolation` implemented in `src/interp.c` as a front-end that
synthesises abscissae from a raw value array and delegates to the proven
`builtin_interpolation_impl` — zero changes to the interpolation numerics.

- **Behaviour matches Wolfram**: `[2.5]`→2.4375, `{{0,1}}` domain prints exactly
  with `[0.5]`→3.875, 2-D grids, explicit positions (x²→6.25), Spline/Hermite/
  order/periodic options, `ListInterpolation[]`→`::argt`.
- **Efficiency**: on `pack.c` AWARE list (packed/NDArray value tensor delisted
  once, not gate-materialised); resulting object applied to packed query points
  uses the existing vectorised buffer path. Integer-endpoint interval grids use
  **exact-rational** interior points, so MPFR data keeps full precision (30.103).
- **Tests**: 39 assertions added to `tests/test_interp.c` — all PASS (verified by
  temporarily running them first). Zero leaks (`leaks` on all forms incl. error
  paths / NDArray / MPFR).
- **Audits**: `check-c99` clean, `check-packed-aware` OK (AWARE=199),
  `check-nd-surfaces` running.
- Symbol interned, `ATTR_PROTECTED` + Options registered, docstring, spec doc,
  changelog all updated.

## Follow-ups (user-requested)

**1. Buffer-direct NDArray construction — DONE.** `builtin_listinterpolation`
merged with its impl and made NDArray-aware: a scalar-valued real `NDArray` (grid
dims == rank) builds the `{coord,val}` table straight from the packed buffer
(`listinterp_entries_from_buffer` + shared `listinterp_coord`), skipping the
intermediate nested `List`; array-valued / complex `NDArray`s still delist. int64
exactness preserved (matches `Interpolation`). Identical results across List /
packed / NDArray. 0 leaks; 42 ListInterp assertions pass; audits green.

**2. Supplied-2nd-derivative Hermite bug — FIXED (root-caused).** It was a
**gcc-16.1.0 IPA-CP miscompile** of `eval_component_double` (clones on constant
`Ksupplied`, const-folds `k` to 1, routes the runtime `k=2` call there →
`build_basis(1)` OOB read). A heisenbug (any print/volatile/sanitizer hid it;
non-deterministic across builds). Bisected with `-fno-ipa-cp-clone`; fixed with a
guarded `__attribute__((noipa))`. The full interp suite (incl. the quintic
`test_supplied_1d`) now passes. Memory:
`project_interp_supplied_quintic_hermite_bug`.

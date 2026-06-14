# Task: Implement UnitStep (piecewise.c)

## Spec
- `UnitStep[x]` = 0 for x<0, 1 for x>=0 (value at 0 is **1**).
- `UnitStep[x1,...,xn]` = 1 iff none of the xi are negative, else 0 (acts as product of UnitStep[xi]).
- `UnitStep[]` = 1.
- Always returns an exact result (0 or 1) for real numeric input.
- Exact symbolic reals (Pi, Sqrt[2], ...) resolved by numerical certification (MPFR precision-doubling).
- Non-real / symbolic args -> leave unevaluated (drop resolved-nonneg args, keep unknowns).
- Attributes: Listable, NumericFunction, Orderless, Protected.

## Plan
- [ ] piecewise.c: `builtin_unitstep` + helpers (real_sign, ustep_class, MPFR certify).
- [ ] piecewise.h: prototype.
- [ ] piecewise_init: register builtin + attributes.
- [ ] info.c: terse docstring.
- [x] sym_names.c: SYM_UnitStep — ALREADY PRESENT (verified line 400/811).
- [ ] deriv.c: D[UnitStep[a1..an],x] via product rule -> Piecewise[{{Indeterminate, ai==0}},0]*UnitStep[rest].
- [ ] tests/test_piecewise.c: extensive UnitStep cases.
- [ ] docs/spec/builtins/elementary-functions.md + changelog 2026-06-08.md.
- [ ] Build clean, run piecewise_tests, valgrind clean.

## Review — DONE 2026-06-14

All explicit deliverables complete and verified.

- **piecewise.c**: `builtin_unitstep` + helpers `ustep_real_sign`,
  `ustep_is_pos_infinity`, `ustep_certify` (MPFR precision-doubling), `ustep_class`.
  Multidim = product semantics; drops proven-non-negative args; non-real/unresolved
  left unevaluated. Registered with attributes Listable|NumericFunction|Orderless|Protected.
- **info.c**: terse docstring (no examples, per house style).
- **deriv.c**: product-rule derivative; `D[UnitStep[x],x]` = `Piecewise[{{Indeterminate,x==0}},0]`,
  multidim matches MMA (`UnitStep[x,y] Piecewise[...]`).
- **tests**: `test_unitstep` (~55 cases) in test_piecewise.c — all pass (exit 0).
  deriv_tests / deriv_symbolic_order_tests / deriv_array_tests all still pass.
- **build**: clean under `-std=c99 -Wall -Wextra`.
- **valgrind**: only the documented macOS baseline (12,800 B / 400 blocks from
  dyld/libobjc/Accelerate); zero leak stacks reference Mathilda source.
- **docs**: elementary-functions.md UnitStep section + changelog 2026-06-08.md entry.

All 23 spec example outputs reproduced exactly (incl. boundary case Sqrt[2]-99/70 → 0).

### Deferred (noted in changelog)
- `Integrate[UnitStep[x],x]` = `x UnitStep[x]` — needs an integral-table hook.
- `UnitStep'[x] === UnitStep''[x]` — needs cleaner `D` of `Piecewise`.

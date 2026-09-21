# Task: Implement `AlgebraicNumberNorm`

Field norm of an algebraic number (product of conjugates). Option `Extension -> None`
(default) gives the relative norm over `Q(theta)`. Attributes `Listable, Protected`.
FLINT-backed, exact, leak-free. Plan: `~/.claude/plans/continuing-on-from-our-crystalline-yeti.md`.

## Steps
- [x] Engine: `qqbar_abs_norm` + `flint_qqbar_algebraic_number_norm` in `flint_qqbar.c` (+ `.h` proto + `#else` stub)
- [x] Leverage: teach `to_qqbar` + `flint_qqbar_is_constant_algebraic` about `GoldenRatio -> (1+Sqrt[5])/2`
- [x] Wrapper: `src/poly/algebraicnumbernorm.c` / `.h` (Extension option, tri-state mapping)
- [x] Symbols: `sym_names.h` / `sym_names.c` (`SYM_AlgebraicNumberNorm`)
- [x] Register: `core.c` init call; attributes `ATTR_LISTABLE | ATTR_PROTECTED`
- [x] Docstring: `info.c` (terse, no examples)
- [x] Options default: `options_builtin.c` (`Extension -> None`)
- [x] Tests: `tests/test_algebraicnumbernorm.c` + `tests/CMakeLists.txt` — all passed
- [x] Docs: `docs/spec/builtins/algebra.md` + `docs/spec/changelog/2026-09-21.md`
- [x] Version: `src/version.h` 0.158 -> 0.159 (tag `v0.159` at commit time)
- [x] Build clean; check-c99 PASS; check-packed-aware OK
- [x] No-drift suite (rootreduce/algebraicnumber/numberfield/zero_test/core...): all PASS
- [x] GoldenRatio leverage spot-check (RootReduce, comparisons, AlgebraicIntegerQ): improved, no drift
- [x] check-array-exactness: 346 probes, 0 MIXED (OK); valgrind: no Mathilda/FLINT leak frames (macOS objc/dyld baseline only)
- [x] Rebuilt code-review-graph
- [ ] Commit + tag `v0.159` (awaiting user — not auto-committed)

## Review

Implemented `AlgebraicNumberNorm[a]` (+ `Extension -> theta` relative norm), Listable/Protected,
FLINT-backed, exact, leak-free. Key design: both cases reduce to one minimal-polynomial read —
absolute norm `(-1)^n c_0/c_n` off the qqbar minpoly, relative norm by norm-in-towers
`(absolute norm)^{[Q(theta):Q(a)]}` with membership via `qqbar_express_in_field`. No new
resultant/`nf_elem_norm`/`mp_core` code.

Bonus (leverage, user-approved): `to_qqbar` + `flint_qqbar_is_constant_algebraic` now recognise
`GoldenRatio = (1+Sqrt[5])/2`, so `RootReduce`, `AlgebraicNumber`, `ToNumberField`, and algebraic
comparisons handle it for free (`RootReduce[GoldenRatio]` → `1/2 (1+Sqrt[5])`, `GoldenRatio < 2`
→ `True`, `AlgebraicIntegerQ[GoldenRatio]` → `True`).

Verification: all spec examples match; new suite `test_algebraicnumbernorm.c` passes; 10 no-drift
suites pass; build clean; check-c99 / check-packed-aware / check-array-exactness green; valgrind
clean (baseline noise only). Two spec examples using `NumberFieldFundamentalUnits` were omitted —
that head is not yet implemented in Mathilda (independent of `AlgebraicNumberNorm`).

Note: `AlgebraicNumberNorm[-2/3]` prints as `-2/3` (Mathilda's rational convention), where the
Mathematica reference shows `-(2/3)`; the value is identical.

# Reduce — Phase 7: Quantifier Elimination (Exists / ForAll / Resolve)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-temporal-cook.md`
Status: **v1 landed, v0.094 (2026-08-26).**

## Tasks

- [x] Add `SYM_Exists`, `SYM_ForAll`, `SYM_Resolve` (sym_names.{h,c}, 3 sites each)
- [x] Add `reduce_cad_qe` public seam in reduce_cad.c + declare in reduce_cad.h
- [x] Create `src/solve/reduce_qe.{c,h}` — front-end, builtins, Case A/B/C dispatch
- [x] Wire `reduce_qe_init()` into `reduce_init`; add Exists/ForAll peel in `builtin_reduce`
- [x] Add `reduce_qe.c` to tests/CMakeLists.txt COMMON_SRC
- [x] Tests: `test_quantifiers_{decision,parametric,decline}` — all pass
- [x] Build clean (make, 0 errors/warnings); reduce_tests green (exit 0)
- [x] `make check-c99` clean; `make check-packed-aware` clean (heads exempt)
- [x] valgrind: no QE-attributable leak (residual = macOS libobjc/dyld init noise)
- [x] Docs (solutions-of-equations.md: Exists/ForAll/Resolve sections + intro)
- [x] Changelog (2026-08-24.md) + version bump 0.093→0.094 + REDUCE_PLAN.md status
- [x] Rebuild code-review graph

## Review

Implemented quantifier elimination for `Reduce` in three regimes:

- **Case A (fully quantified)** — a real-closed-field decision procedure that reuses
  the whole engine: `Exists[{v},φ] = Reduce[φ,v,Reals] =!= False`,
  `ForAll[{v},φ] = ... === True`. No new CAD code.
- **Case B (one free variable)** — new public seam `reduce_cad_qe` in reduce_cad.c:
  build CAD with the free var outermost, read each cell's `Exists`/`ForAll` verdict
  (`!empty` / `all_true` — already computed by `cad_build`), emit via
  `rru_emit_sign_diagram`.
- **Case C (≥2 free vars / alternating / non-Reals / algebraic boundary)** — declines
  (NULL), sound.

Key design wins vs. the original sketch: the fold was already computed by `cad_build`
(no new projection machinery); `ForAll` done directly via `all_true` (no
`Not[Exists[Not]]`, no `rform_not_*`); the fully-quantified case reuses `builtin_reduce`
wholesale. Soundness invariant preserved everywhere: undecidable/out-of-scope → NULL.

## Deferred (future work)

- Phase 6b (real-algebraic-coefficient fibres) unblocks: ≥2-free-var parametric QE and
  the algebraic-boundary `Resolve[Exists[x, x^2+bx+c==0]]→b^2-4c>=0` example.
- Alternating quantifier prefixes (`ForAll[x, Exists[y, ...]]`) in the parametric case.
- R6 completeness: `cad_leaf` sometimes over-declines a verdict-decidable fibre
  (e.g. `Exists[y, 0<y<x]`) — sound decline, could add a verdict-only leaf later.
- Phase 8 companions: `LogicalExpand`, `FindInstance`, `CylindricalDecomposition`.

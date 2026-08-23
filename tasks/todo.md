# Reduce Phase 6d — n-variable CAD over the Reals

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-fizzy-giraffe.md`

## Stage 0 — Baseline (golden-master must be green first)
- [x] Build main binary + `reduce_tests` + `reduce_corpus_tests`; confirm all green (2-var pins + corpus) — green (94/94)

## Stage A — recursive engine (hybrid: nu>=3 new path, nu==2 untouched)  ✅ DONE
- [x] `PolySet` + `cad_project_out`; iterated McCallum stack `pstack[0..d-1]`
- [x] `_n` helpers: `subst_n`, `is_poly_n`, `atom_truth_n`, `form_truth_n`, `cell_dead_n`
- [x] `cad_leaf` = parameter-generalized `lift_fiber`; `symbolic_branch_lvl`/`bound_expr_lvl`
- [x] `cad_recurse` + `cad_leaf_emit` (per-level rational gate; nullification bail; partial-CAD prune; flat `parts` DNF)
- [x] `reduce_cad_nvar` driver; gate routes `nu>=3` there (nu==2 old path UNTOUCHED = byte-identical)
- [x] Bugfix: only fibre-var-bearing factors trigger nullification (a lower-var factor vanishing at its section is skipped)
- [x] **Parity gate**: 2-var pins byte-identical; corpus 103/103 (form-agnostic oracle verifies all closed/verbose cases); leak-clean; check-c99 clean
- [ ] (stretch) unify nu==2 onto the recursive path only if it reproduces the pins byte-identically — deferred

## Stage B — n-D boundary merge  ✅ DONE (v0.088)
- [x] Restructure emission to a cell TREE (`CADRegion`/`CADCell`, `cad_build`) + merge (`cad_region_expr`)
- [x] `cad_templates_equal` (structural) + `cad_absorbable` (sampling-based equality via `cad_sample_cell`/`formula_truth_at`)
- [x] Closed regions close outer ranges: `x^2+y^2+z^2<=1 -> -1<=x<=1 && ...`; sphere surface, half-ball, 4-var closed ball; strict stay open
- [x] Sound fallback: undecidable comparison leaves verbose form; corpus oracle certifies (105/105)
- [x] Leak-clean (valgrind), check-c99 clean, solve_tests green, version 0.088

## Tests / docs / hygiene
- [ ] `tests/test_reduce.c`: flip `x^2+y^2+z^2<1` decline pin (:616-617) → solved; add <=1, >=0→True, <0→False, x y z>0, half-ball, nu==4 ball, `<=2`→decline
- [ ] `tests/reduce_corpus.m`: flip `dec-nl-multivar3` → solved; add cad3-*/cad4-* rows
- [ ] `reduce_cad.h` header prose (recursive engine, d>=3 rational-fibre scope)
- [ ] docs/spec/builtins/solutions-of-equations.md (Reduce bullets + deferral paragraph)
- [ ] changelog `docs/spec/changelog/2026-08-24.md` (verify Monday-of-ISO-week); `src/version.h` 0.086 → 0.087
- [ ] `make check-c99`; `valgrind --leak-check=full` on corpus

## Review

**Delivered (Stage A):** `Reduce[..., {x1..xn}, Reals]` now solves nonlinear real
inequalities in 3+ effective variables via a recursive McCallum-projection CAD
(`reduce_cad_nvar` + `cad_recurse`/`cad_leaf` in `src/solve/reduce_cad.c`). The
2-variable path is byte-identical (untouched); the new engine is a hybrid that
reuses the shared primitives and the generalized leaf.

- Strict inequalities → clean nested form (`x^2+y^2+z^2<1`, box, octant, planes).
- Closed regions → correct but verbose union of cells (Stage-B boundary merge
  deferred — cosmetic only; sound + complete either way).
- v1 rational-fibre regime: irrational non-innermost breakpoints decline
  (`x^2+y^2+z^2<=2`); interval nullification declines (6e deferred).
- Key bugfix: only fibre-variable-bearing factors trigger the nullification bail
  (a lower-variable factor vanishing at its own section is skipped) — this is why
  `x y z > 0` initially declined.

**Verification:** reduce_tests green (new `test_cad_nvar`); corpus 103/103
(form-agnostic sample-point oracle certifies the verbose closed cases);
solve_tests green (no collateral); `make check-c99` clean; valgrind leak-clean
(only macOS ObjC-runtime baseline noise); version 0.087; docs + changelog updated.

**Deferred / follow-ups:** Stage-B n-D boundary merge (close outer ranges for
closed regions); Phase 6b (real-algebraic-coefficient fibre isolation to widen
past the rational-fibre regime); Phase 6e (McCallum well-orientedness
augmentation); unifying nu==2 onto the recursive path.

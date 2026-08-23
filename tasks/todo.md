# Reduce — Phase 6: Cylindrical Algebraic Decomposition (Reals)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-fluffy-kite.md`

Goal: 2-variable McCallum CAD for `Reduce[..., Reals]` on multivariate nonlinear
inequalities. Soundness over completeness — undecidable → `NULL` (unevaluated).

## Step 1 — Extract shared real-algebraic primitives (refactor)  ✅
- [x] `src/solve/reduce_real_util.h` — public API (prefixed `rru_`)
- [x] `src/solve/reduce_real_util.c` — moved sign/sample/root helpers +
      provenance-aware `rru_collect_roots`
- [x] `reduce_univar.c` — includes the header, moved statics dropped, rewired
- [x] both added to `tests/CMakeLists.txt` COMMON_SRC (+ reduce_cad.c)
- [x] main + reduce_tests (all pass) + reduce_corpus_tests (0/67) → no regression

## Step 2 — CAD engine  ✅
- [x] `src/solve/reduce_cad.{c,h}` — `reduce_cad(F, vars, nv)`
- [x] collect + gate (PolynomialQ[{x,y}]) + factor atoms → squarefree basis B
- [x] projection (disc_y, ldcf_y, pairwise res_y; factor; keep x-only factors)
- [x] base x-cells (roots, order/dedup, samples; irrational-section decline)
- [x] lift per x-cell (substitute, isolate y-roots w/ provenance, y-cells)
- [x] truth at full sample via qqbar; DNF fold; `-2` bail; all-cells-true → True
- [x] partial-CAD pruning (xcell_dead)
- [x] symbolic sector emission + branch matching + decline on non-clean Solve
- [x] soundness bails (nullification, irrational section, non-clean fiber, nv≥3)

## Step 3 — Dispatch  ✅
- [x] `reduce.c` — `reals && nv>=2`: FM then CAD fallback; effective-var pruning in reduce_cad

## Step 4 — Tests  ✅
- [x] flipped 3 CAD placeholders; new `test_cad_real` block (pinned FullForm)
- [x] Phase-6 corpus rows in `reduce_corpus.m` (11 new; total 77)
- [x] reduce_tests + reduce_corpus_tests (77/77) green; solve_tests unregressed

## Step 5 — Verify + docs  ✅
- [x] `make check-c99` clean; both new TUs compile FLINT-off and MPFR-off (bail stubs)
- [x] valgrind CAD probe — zero reduce_cad/reduce_real_util frames in lost blocks
- [x] docs: solutions-of-equations.md + changelog 2026-08-17.md + REDUCE_PLAN status
- [x] rebuild code-review graph; memory: project_reduce_cad_2var_engine

## Review
- **Engine**: 2-var McCallum CAD. Factor atoms → squarefree basis (confines
  nullification to sections). Projection = disc/ldcf/pairwise-res, factored.
  Base sign diagram lifted fibre-by-fibre; truth via exact qqbar oracle; symbolic
  sector emission with per-root provenance + branch matching. Partial-CAD prune.
- **Reuse**: extracted `reduce_real_util.{c,h}` shared by the univariate engine
  and CAD (one exact sign/root/sample implementation) — no duplication.
- **Soundness**: every undecidable path (qqbar -2, irrational section, interval
  nullification, non-clean fibre Solve, nv≥3) declines to NULL. Verified sound on
  `x^2==2 && y<x` (irrational section → decline) and nv=3 (→ decline).
- **Correctness**: 77/77 corpus (sample-point + witness), exact unit pins, plus
  hand-checked conics/hyperbola/parabola/products/cubic fibres. Leak-clean.
- **Deferred (documented)**: n-var (6d), well-orientedness augmentation (6e),
  Complexes nonlinear equations.

## Step 6 — Boundary-merge pass  ✅
- [x] structured per-cell `YRegion` (symbolic-bound segments) from `lift_fiber`
- [x] `yregion_at_point` (template limit), `regions_equal` (sampled, sound),
      `templates_equal`, `breakpoint_absorbable`
- [x] merge runs of same-template interval cells; absorb interior/boundary
      breakpoints; standalone sections for non-absorbed non-empty fibres
- [x] closed disk → `-1<=x<=1 && ...`; half-disk → `0<=x<=1`; circle → arcs;
      point → `x==0&&y==0`; strict stays open; punctured hole preserved
- [x] updated unit pins + 4 new corpus rows (0/81); solve unregressed;
      valgrind leak-clean; build + check-c99 clean

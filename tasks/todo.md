# Reduce Phase 6e: McCallum well-orientedness augmentation

## Goal
n-var CAD over Reals currently bails to NULL on positive-dim fibre nullification.
Fix soundly via McCallum augmentation: add the nullified poly's coefficients to
the projection level below, rebuild, retry (bounded). Reproducing case:
`Reduce[x^2 - y^2 z == 0, {x,y,z}, Reals]` (currently declines; MMA answers it).

Hard invariant: undecidable / still-nullifying-after-augment => NULL.

## Tasks
- [x] reduce_cad.c: `NullReport` struct (owned poly, level, hit)
- [x] reduce_cad.c: thread `NullReport* nz` through cad_build + cad_leaf; set at
      the two positive-dim bails
- [x] reduce_cad.c: augment+restart loop in reduce_cad_nvar (max 4 rounds);
      re-project lower stack, re-alloc caches per round; reduce_cad_qe passes NULL
- [x] tests/test_reduce.c: test_cad_well_oriented (3 solve; 4-var still declines)
- [x] tests/reduce_corpus.m: three 6e-* rows (sample-oracle equivalence, 170/170)
- [x] Build + reduce_tests (pass) + corpus (170/170) + solve (pass/99)
- [x] Efficiency: 4 common 3-/4-var cases in 0.14s (augment loop = 0 extra rounds)
- [x] valgrind: macOS baseline 13,440/420, none attributable
- [x] make check-c99 (clean)
- [x] docs + changelog + REDUCE_PLAN.md; version 0.125

## Review

**Shipped (v0.125) — Phase 6e, closing the CAD roadmap.** McCallum
well-orientedness augmentation in `reduce_cad.c`:
- Detect a positive-dimensional fibre nullification, report the offending factor +
  fibre level via `NullReport` (threaded through `cad_build`/`cad_leaf`).
- In `reduce_cad_nvar`, add all that factor's coefficients (w.r.t. its fibre var) to
  the projection level below via `add_proj`, re-derive the lower projection stack,
  rebuild. Bounded to 4 rounds; persistent nullification → sound decline.
- `Reduce[x^2 - y^2 z == 0, {x,y,z}, Reals]` now solves (was declining).

**Key design decision:** resampling the interior point to avoid the nullification
would be UNSOUND for equation atoms (drops the on-locus solution slice, e.g. `x=0`).
Only augmentation (making the locus a cell boundary) is sound.

**Zero-cost for the common case:** the augmentation path is entered only on an
actual nullification; every non-nullifying CAD runs the identical code. Verified:
common 3-/4-var cases unchanged and fast (0.14s).

**Verification:** test_cad_well_oriented; three 6e-* corpus rows (oracle certifies
input≡output); independent grid equivalence (0 mismatches); reduce+solve
suites/corpora green; check-c99 clean; valgrind at baseline (restart frees each
superseded pstack/cache generation).

**Only open Reduce item:** 7-extended's algebraic free-variable-boundary sub-case
(free leading coefficient) — sound-declines.

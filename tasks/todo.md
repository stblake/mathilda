# Reduce Phase 6b — Real-Algebraic-Coefficient Fibre Isolation (core + tower)

## Goal
Lift the CAD "rational-fibre" restriction so non-innermost breakpoints may be
irrational algebraic. Unlocks 3-var sphere octant, closed ball, `x^2==2 && y<x`.

## Tasks
- [x] Read the key reduce_cad.c regions + reduce_real_util + flint_qqbar signatures
- [x] New file `src/solve/reduce_algfiber.{c,h}` — `rru_algebraic_fiber_roots`
      (iterated-resultant tower projection + integer-isolate + qqbar exact filter
      + order + resource guard). USE_FLINT-guarded.
- [x] Thread `asgdef[]` (defining pstack factor per irrational level) through
      cad_build / cad_leaf / lift_fiber / reduce_cad_qe; populate from
      rru_collect_roots provenance.
- [x] Relax rational gates (reduce_cad.c:960-962 and :1448); route fibre isolation
      through the new primitive when asg has ≥1 irrational coord; keep all-rational
      fast path unchanged.
- [x] Verify emission at algebraic samples: octant, closed 3-/4-ball, cubic,
      hyperbola, x^2==2&&y<x all correct; Pi declines soundly.
- [x] Tests: test_cad_algebraic_fibre in tests/test_reduce.c (added, passing);
      stale decline assertions updated (cad_nvar, cylindrical_decomposition).
- [x] cad-alg-* corpus rows (7 added; corpus 166/166 pass).
- [x] check-c99 clean; reduce suite + corpus no-regression; valgrind macOS
      baseline 13,440/420 (identical at 1x and 25x — zero per-call leak).
- [x] Update reduce_cad.h scope doc.
- [x] Docs: solutions-of-equations.md, changelog 2026-08-24.md, version 0.113,
      REDUCE_PLAN.md status, memory notes.
- [x] Rebuild code-review-graph.

## Review — Phase 6b landed (v0.113), all verification green

- New `src/solve/reduce_algfiber.{c,h}` — `rru_algebraic_fiber_roots`:
  iterated-resultant tower projection to ℚ + integer-coeff isolation + exact qqbar
  spurious-root filter + order + resource ceiling. USE_FLINT-guarded (off→decline).
- `reduce_cad.c` — removed both rational-fibre gates; shared `isolate_fiber_at`
  dispatch (all-rational fast path unchanged, ≥1 algebraic → new primitive);
  borrowed `asgdef[]` (section defining factor via provenance) threaded through
  cad_build/cad_leaf/lift_fiber/reduce_cad_qe; 2-var driver gained base-root
  provenance (`bxfac`). `reduce_cad.h` scope doc updated.

Verification: reduce_tests PASS (new test_cad_algebraic_fibre; stale decline
assertions updated); reduce corpus 166/166 (7 new alg-* rows); solve_tests PASS +
solve corpus 99/99 (no regression); check-c99 clean; valgrind macOS baseline
13,440/420 — IDENTICAL at 1× and 25× octant, zero per-call leak, no leak stack
references the new code.

Newly solved (declined before): x^2==2 && y<x; x^2+y^2+z^2<=2; sphere positive
octant; x^3+y^3==1 && x>0 && y>0; hyperbola branches; closed 3-/4-balls; CD too.
Transcendental (Sqrt[Pi]) still declines soundly.

Deferred (documented): 6e (well-orientedness); 7-extended (≥2-free-var / algebraic
discriminant-variety QE — needs a nested-QE emitter). Closed-region bound cosmetic
asymmetry is pre-existing (pinned in the shipped x^2+y^2+z^2<=1 test), not new.

Not committed — awaiting user review.

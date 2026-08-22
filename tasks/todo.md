# Close 4 Diophantine gaps in Solve[…, Integers] vs sympy

Plan: /Users/user/.claude/plans/cozy-imagining-pizza.md

## Gap 2 — General Pythagorean `x1²+…+xk² == y²` (k≥3), unbounded  [DONE]
- [x] New `src/solve/solveint_pythag.c`: `si_solve_general_pythagorean(SICtx*)`
- [x] Prototype in solveint_internal.h; wire into solveint.c chain
- [x] Build; verify `x²+y²+z²==w²`, `x²+y²+z²+w²==v²` return a family (residual→0)
- [ ] Unit tests (substitution → identity 0; known solution reproduced)

## Gap 1 — General homogeneous ternary quadratic (cross-term / general coeff)  [DONE]
- [x] New `src/solve/solveint_ternary_general.c`: `si_solve_ternary_general(SICtx*)`
- [x] Diagonalize (mpq) + Legendre + Holzer witness + chord family in original coords
- [x] Wire after si_solve_ternary_quadratic
- [x] Verify `4x²−5y²+z²==0`, `x²+3xy+2y²−z²==0`; PROJECTIVE completeness (matches sympy)
- [x] Family is projectively complete (up-to-scaling), same representation as sympy/MMA
- [ ] Unit tests

## Gaps 3 & 4 — General binary quadratic parametric (parabolic + hyperbolic)  [DONE]
- [x] New `src/solve/solveint_bqf_parametric.c`: `si_solve_bqf_parametric(SICtx*)`
- [x] δ==0 parabolic families (congruence loop, x/y quadratic in t) — 0 missing
- [x] δ>0 non-square hyperbolic: base search + automorphism M closed-form family — complete
- [x] (did NOT need to de-static genpell — hyperbolic uses M_aut + si_pell_cf directly)
- [x] Wire into chain; verify `x²−4xy+4y²−3x==0` (2 fam), `x²−3xy+y²==1` (6 fam)
- [x] Unit tests

## Wrap-up
- [x] tests/CMakeLists.txt: add 3 files to COMMON_SRC
- [x] make check-c99 clean; solve_integers_tests + solve_tests pass; valgrind = baseline only
- [x] Fixed real Holzer-bound bug (x²+y²==65z² now solves projectively)
- [ ] Docs: solutions-of-equations.md + changelog 2026-08-17.md
- [ ] Review section

## Review

All four unbounded Diophantine gaps vs sympy are closed. Three new files, wired
into the `si_solve_*` unbounded-family chain in `solveint.c`; no new builtin
symbol (all internal `Solve` paths), no packed/ND/Compile surface.

- **Gap 2 (general Pythagorean)** — `solveint_pythag.c`. Stereographic family,
  O(k), pure symbolic. `x²+y²+z²==w²` etc. residual→0.
- **Gap 1 (general ternary quadratic)** — `solveint_ternary_general.c`. mpq
  congruent diagonalisation → Legendre decision + Holzer witness → chord family
  in original coords (only the witness is mapped back). Projectively complete
  (matches sympy/MMA up-to-scaling); anisotropic→trivial `{{0,0,0}}`. Also picks
  up multi-rep `k` (65) and non-symmetric diagonals the symmetric solver declines.
- **Gaps 3&4 (binary quadratic)** — `solveint_bqf_parametric.c`. Parabolic
  (δ=0) congruence families; hyperbolic (δ>0 non-square) base-search +
  automorphism `M=[[t-Bu,-2Cu],[2Au,t+Bu]]` closed-form families. Both verified
  0-missing against brute force.

Bug found & fixed along the way: the ternary witness search used the wrong
Holzer-bound coefficient pair (each searched var's bound must involve `C[si]`,
the solved-for coefficient), which spuriously declined solvable forms such as
`x²+y²==65z²`.

Behaviour changes to existing tests (both are improvements — the engine now
solves cases it used to defer):
- `x²+y²==65z²`: symmetric solver declines → general solver returns a projective
  family. `test_ternary_quadratic` updated (was "Solve", now "List" + residual 0).
- `y==x²`: now solved as the parabola `x=t, y=t²`. Moved out of
  `test_deferred_unevaluated` into `test_bqf_parabolic`.

Verification: `make check-c99` clean; `solve_integers_tests` (incl. 4 new tests)
and `solve_tests` pass; valgrind shows only the documented macOS libobjc baseline
(identical to a trivial script — zero leaks from new code). Completeness checked
by brute force for every gap.

Scope declines (sound — unevaluated, never wrong): inhomogeneous ternary
quadratics (linear/const term); rank-deficient ternary forms; hyperbolic BQF with
linear D/E terms; weighted general-Pythagorean coefficients; witnesses beyond the
Holzer box.

## Follow-up (benchmark-driven, done)

Benchmarked Mathilda vs sympy vs PARI on all four gaps (simple→hard). Mathilda is
sub-ms across the board and 50–1000× faster than sympy where both engage; PARI's
qfsolve (gaps 1–2 only) is µs but returns one witness, not the family. Two fixes
fell out and were applied:
- **Ternary small-witness pre-scan** — `99991x²−99989y²−2z²==0` (witness (1,1,1))
  used to decline; now instant.
- **Gap-3 large-δ decline closed** — the hyperbolic base search was bounded by
  `s_min·unit`, so δ=61 (unit ≈3.5e9) and δ=10⁶ declined. Rewrote to reduce to
  `X²−δY²=N` and scan `Y` over the Nagell bound (unit-size-independent) via new
  `bqf_pos_base`. Both now sub-ms; small cases verified unchanged + complete.
  Tests added; docs/changelog/memory updated; valgrind = baseline only.

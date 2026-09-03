# Task: Generalized Trig & Hyperbolic Equations in Solve[] / Reduce[]

Plan: `/Users/user/.claude/plans/we-also-need-solve-linear-phoenix.md`

## Stages

- [x] Stage 0 — New file `src/solve/solvetrigpair.{c,h}`, build wiring, init hook
- [x] Stage 1 — Engine 1 normalize: to_sincos + Together + numer/denom split
- [x] Stage 2 — Engine 1 factor: sum-to-product / reverse-angle-addition rule tables
- [x] Stage 3 — Engine 1 solve atoms via solveinv (+ poly_in_kernel fallback)
- [x] Stage 4 — Shared pole/spurious gate (zero_test_decide, C[k]->0 representative)
- [x] Stage 5 — Degenerate cases: N≡0 -> True ({{}}); N var-free const -> False ({})
- [x] Stage 6 — Engine 2: self-contained rational exp-generator (try_exp_generator)
- [x] Stage 7 — Dispatch insertion in solve.c (Engine 1 -> solvetrig fallback)
- [x] Stage 8 — Reduce routing (reduce_has_trig_pair) + fix 2 Reals bugs (family collapse)
- [x] Stage 9 — Init/registration/docs/spec/changelog/memory
- [x] Verify — build, behavioral corpus, unit tests, valgrind, check-c99

## Review

Delivered generalized two-argument trig/hyperbolic solving in Solve[] and
Reduce[] via a new `src/solve/solvetrigpair.c`.

**Coverage (both user corpora, ~41 cases): 39 solved/True/False correctly**;
the 2 declines (`Sinh[a x+b y]+Cosh[c x+d y]==0` and the Sech/Csch analogue)
are symbolic-coefficient mixed Sinh/Cosh — no Engine-1 identity and no
Engine-2 generator, documented non-goals.

**Engine 1 (primary, clean output):** reciprocal-normalize (Tan/Sec/... ->
Sin/Cos) under `trig_canon` suppression -> Together -> numerator factored by an
explicit sum-to-product / reverse-angle-addition rule table -> each single-arg
atom solved by the existing `solveinv` (falling back to
`solvetrig_solve_poly_in_kernel` for `p(Cosh[y])`-shaped atoms) -> pole gate
(drop families where a denominator factor is identically zero). Handles
same-head ±, mixed cofunctions, reciprocals, SYMBOLIC coefficients,
higher-multiplicity args, parity->True, spurious traps->False, and
constant-RHS cases that collapse via product-to-sum (Cosh[x+y]+Cosh[x-y]==2).

**Engine 2 (fallback, Root/Log output):** a self-contained single rational
exponential generator `u=E^(sigma y)` (sigma = I/q circular, 1/q hyperbolic;
q from the rational gcd of the var-coefficients) for non-collapsing
constant-RHS (Tan-Tan==1) and numeric mixed Sinh/Cosh. Declines symbolic
coefficients (no generator).

**Reduce:** new `reduce_has_trig_pair`/`reduce_is_trig_resid` gate routes
Complexes AND Reals; two pre-existing Reals bugs fixed
(`Reduce[Exp[y]==3,y,Reals]` False->`y==Log[3]`, `Reduce[Sinh[y]==0,y,Reals]`
True->`y==0`) via a family-collapse in `reduce_eq_transcendental` that keeps
only the real members of a periodic family (real period -> keep all, imaginary
period -> isolated member at C[k]=0).

**Verification:** all Engine-1 answers back-substitute to 0; Engine-2 verified
numerically (x subst -> N -> Chop); solve/reduce/radicals/dsolve/trig_canon/
simplify suites all pass; `make check-c99` clean; valgrind at baseline for
Engine 1 and one-time +40 bytes for Engine 2 (a lazy-init cache, not per-call).

**Key gotchas hit:** (1) evaluator auto-folds Sin/Cos->Tan and 1/Cos->Sec, so
the reciprocal pipeline must run under `trig_canon_suppress_inc/dec`;
(2) `Together` leaves exponents unexpanded (net-zero var terms), so Engine 2
expands each E-exponent before scanning; (3) `Simplify` has a pre-existing
per-call leak on trig-of-fractional-argument inputs, so the degenerate check
uses the leak-free `zero_test`+`var_in` instead of a Simplify probe.

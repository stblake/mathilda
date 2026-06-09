# Task: Jeffrey–Rich continuous Weierstrass substitution integrator

Implement the continuous extension of the Weierstrass (tan(x/2)) substitution for
rational functions of trigonometric functions per Jeffrey & Rich, *The Evaluation
of Trigonometric Integrals Avoiding Spurious Discontinuities* (ACM TOMS 1994), and
generalise it to hyperbolic functions.

File: `src/calculus/integrate_jeffrey.c` (+ `.h`).
Exposed as `Method -> "Weierstrass"` and `Integrate`Weierstrass`. Wired into the
main `Integrate[]` cascade.

## Algorithm (paper §5)
1. Detect: integrand is a rational function of {Sin,Cos,Tan,Cot,Sec,Csc}[x]
   (trig mode) OR {Sinh,Cosh,Tanh,Coth,Sech,Csch}[x] (hyperbolic mode), with the
   kernel argument exactly the integration variable `x`. Every occurrence of `x`
   must be inside such a kernel. Reject mixed trig+hyperbolic.
2. Substitute (Table I choice (a) for trig, u = tan(x/2); tanh(x/2) for hyp):
   replace each kernel by its u-rational form, multiply by the Jacobian dx, then
   Cancel[Together[.]] -> rational function g(u), free of x.
3. Recurse: Integrate[g, u] (closes via BronsteinRational). ghat_u(u).
   - Trig: if (a) does not close, retry with choice (c), u = cot(x/2), b = 0.
4. Continuity correction (TRIG ONLY -- tanh(x/2) is a monotone bijection
   R->(-1,1), so the hyperbolic substitution introduces no spurious
   discontinuity):
   K = lim_{u->+inf} ghat_u - lim_{u->-inf} ghat_u  (choice (a); negated for (c)).
   Computed by a small linearity-aware evaluator that splits Plus terms, strips
   u-free constant factors, and resolves ArcTan/Log limits via the limit of their
   inner argument (the core Limit engine cannot distribute over Times/Plus at inf,
   but resolves rational arg-limits and bare ArcTan[n u]). If any term diverges /
   the limit is non-finite (genuine singularity), K is dropped (no correction) --
   matches paper section 4.
   Result: g(x) = ghat(x) + K*Floor[(x - b)/p],  b = Pi (a) / 0 (c),  p = 2 Pi.
5. ghat(x) = ghat_u with u -> Tan[x/2] (a) / Cot[x/2] (c) / Tanh[x/2] (hyp).

Correct by construction (exact substitution identity + verified rational
sub-integral + Floor' = 0 a.e.); no differentiate-back verification (matches
integrate_linrad.c rationale; the Floor term defeats symbolic D anyway).

## Steps
- [ ] src/calculus/integrate_jeffrey.h -- public API (try / builtin / init).
- [ ] src/calculus/integrate_jeffrey.c -- detection, substitution table,
      jump evaluator, core driver, builtin, init + docstring + attributes.
- [ ] Wire into integrate.c: METHOD_WEIERSTRASS enum, parse_method_option
      string "Weierstrass", try_weierstrass wrapper, cascade insertion
      (after derivdivides, before risch), dispatch case, integrate_jeffrey_init().
- [ ] Update Integrate docstring with the new method name.
- [ ] tests/test_integrate_jeffrey.c -- trig + hyperbolic cases, continuity
      (Floor) cases, strict no-match, method plumbing. Verify via
      Simplify[D[Int /. Floor[_]->0, x] - f] === 0.
- [ ] CMakeLists.txt: add file to COMMON_SRC; add integrate_jeffrey_tests target.
- [ ] Build (make + cmake), run new tests + integration regression suite.
- [ ] valgrind the new test binary -- no leaks.
- [ ] docs/spec/builtins/calculus.md + weekly changelog (Mon 2026-06-08).

## Review

Done. `src/calculus/integrate_jeffrey.{c,h}` implement the continuous
Weierstrass substitution (trig + hyperbolic), exposed as
`Method -> "Weierstrass"` and `Integrate`Weierstrass`, wired into the Automatic
cascade **ahead of** `DerivativeDivides`/`RischNorman` (per user feedback:
deterministic, correct-by-construction methods go first).

Verified:
- Headline trig cases match the paper: `Integrate[3/(5-4Cos[x]),x]` = eq. (10);
  `Integrate[1/(2+Cos[x]),x]` continuous with `Floor` correction.
- Hyperbolic: `Integrate[1/(2+Cosh[x]),x]` (was unevaluated) now closes, NO
  `Floor` term (tanh(x/2) is a pole-free monotone bijection — no spurious jump).
- `TrigExpand` pre-pass handles multiple-angle args (`Cosh[x] Cosh[2 x]`).
- Automatic gate leaves polynomial trig clean (`Integrate[Sin[x],x]` = -Cos[x]).
- New suite `integrate_jeffrey_tests` passes; regression suite
  (`integrals_tests`, `integrate_dispatch/derivdivides/linrad/quadrad/
  linratiorad/unknown_tests`) all green.
- valgrind: no leak stacks reference `integrate_jeffrey` (only macOS baseline
  noise ~12.8 KB/402 blocks).

Design notes:
- Jump `K` computed in u-space (limits at ±Infinity), summed term-by-term with a
  linearity/ArcTan/Log-aware evaluator because the core `Limit` engine does not
  distribute over `Plus`/`Times` at infinity. Divergent limit → no correction
  (genuine singularity), graceful degradation to the standard antiderivative.
- Correct by construction (exact substitution + verified rational sub-integral +
  `Floor' = 0` a.e.); no differentiate-back gate, matching `integrate_linrad.c`.
- Docs: `docs/spec/builtins/calculus.md` (cascade list, method list, dedicated
  section) + `docs/spec/changelog/2026-06-08.md`.

# Task: Fix post-integration simplification (radical antiderivatives)

## Problem
`Integrate[1/(x^3 (a + b x)^(1/3)), x]` returned a correct but verbose result:
unsimplified rational pieces over denominators in the compound `(a+b x)` and
stray `(1/a)^(k/3)` constants. Target = Mathematica's compact form with
`6 a^2 x^2` denominator and `a^(7/3)` coefficients.

## Plan (done)
- [x] New `intsimp_finalize(r, x)` in `src/calculus/intsimp.c`, wired at the
      single `builtin_integrate` chokepoint (`integrate.c`). Takes ownership.
- [x] Scope to **radical-bearing** results (`it_has_x_radical`); non-radical /
      nested-`Integrate` results pass through untouched.
- [x] x-free, positive-base `PowerExpand` (`it_normalize_xfree_powers` +
      `it_pos_base`): `(1/a)^(2/3) -> a^(-2/3)`.
- [x] Algebraic recombine: `Cancel[Together]` then `Numerator`/`Denominator`
      `Expand` + `Cancel` so `(a+b x)`-poly denominators collapse in x
      (`it_recombine_algebraic`); only kept when it strictly shrinks.
- [x] Inverse-trig arg distribute + sign-pull (`it_tidy_invtrig_args`).
- [x] Docs: changelog `2026-06-15.md` + `docs/spec/builtins/calculus.md`.

## Review
- **Result now matches the target** (algebraic part collapses to `6 a^2 x^2`,
  `ArcTan[1/Sqrt[3] + ...]` with sign pulled, `a^(7/3)` denominators, clean
  `Log` args) up to an immaterial `Log` sign (constant of integration).
  `PossibleZeroQ[D[r,x] - integrand] = True`.
- **No new builtins/attributes.** One new function + helpers in an existing
  module; prototype in `intsimp.h`.
- **Tests:** `integrate_linrad/quadrad/linratiorad`, `integrate_unknown`
  (previously segfaulted — fixed by the nested-`Integrate` guard + radical-only
  scoping), `intrat`, `integrate_dispatch`, `integrals`, `intrischnorman`,
  `integrate_derivdivides` — all pass. `IntegrateRationalTests.m` corpus passes
  (1 diff ≤ baseline 2, 0 crashes; was 9 diffs + 2 crashes mid-development).
- **valgrind:** no new leaks — only the documented 12,800 B / 400 blocks macOS
  dyld/Accelerate baseline; no `intsimp`/`it_*` frames in any leak stack.

## Key gotchas hit (see lessons.md)
1. Broad x-free `PowerExpand` on *all* results changed negative/complex constant
   values across branch cuts → corpus DIFF regressions. Gated to positive bases
   AND radical-bearing results only.
2. Re-evaluating results that contain nested `Integrate` re-enters the
   integrator → recursion-limit segfault. Skip those.
3. Pulling an ArcTan sign by `Times[-1, NegativePlus]` is undone by ArcTan's own
   oddness; must `Expand` the negation into a genuinely positive sum first.

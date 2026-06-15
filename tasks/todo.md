# Task: Simplify must recover the integrand from D[Integrate[…radical…]]

## Problem
`Simplify[D[Integrate[1/(x^3 (a + b x)^(1/3)), x], x]]` returned a 4-term sum of
fractions over `a^(k/3)` and `(a+b x)^(k/3)` instead of the integrand
`1/(x^3 (a + b x)^(1/3))`. The integral is correct
(`PossibleZeroQ[D[r,x]-integrand] = True`); the gap is in `Simplify`: `Together`/
`Cancel` handle one radical generator but not two, and the collapse needs
reduction modulo the generator relations `s^3 = a`, `t^3 = a + b x`.

## Plan (done)
- [x] New `simp_radical_rational` (`src/simp/radrat.c` + `radrat.h`), wired into
      `simp_dispatch` (`src/simp/simp_search.c`) like `simp_trig_rational`.
- [x] Detect ≥2 distinct, provably-positive radical bases; substitute each base
      (outer-first) to a fresh generator (`poly_subst_radical_to_gen`).
- [x] Build relations `g_k^{q_k} - V_k` (compound bases only; ring-membership
      guard drops degenerate / non-nesting cases).
- [x] Together → reduce num/den mod the ideal (`PolynomialRemainder`) →
      rationalise denominator (`PolynomialExtendedGCD`) → substitute back →
      `Cancel[N / Factor[D]]`.
- [x] Strict `SimplifyCount` gate: return only when it beats the input.
- [x] Branch safety: non-positive (negative/complex numeric) bases left untouched.
- [x] Build wiring: makefile auto-globs `src/simp/*.c`; added to
      `tests/CMakeLists.txt` COMMON_SRC.
- [x] Docs: `docs/spec/builtins/simplification.md` + changelog `2026-06-15.md`.

## Review
- **Target fixed:** `Simplify[D[r,x]] = 1/(x^3 (a+b x)^(1/3))`,
  `PossibleZeroQ[Simplify[D[r,x]] - integrand] = True`. Sibling
  `1/(x^2 (a+b x)^(1/3))` also recovers (PossibleZeroQ True).
- **No new builtins/attributes** — one internal pass, gated and inert by default
  (returns NULL unless ≥2 positive radical generators AND a strict improvement).
- **Algorithm validated in-REPL before coding:** numerator of
  `(D[r,x] − integrand)` reduces to **0** mod `{t^3 − s^3 − b x}`; the full
  pipeline (Together → reduce → ExtendedGCD rationalise → back-substitute →
  Cancel[N/Factor[D]]) reproduces the integrand. The naive variant
  (reduce-then-Cancel without rationalisation, or rationalise the *un-reduced*
  numerator) was rejected — the first leaves a relation-dependent factor, the
  second blows up the multivariate GCD (runaway processes; matches the
  "Simplify multi-generator explosion" lesson). Reduce-first then rationalise
  stays small.
- **Tests:** all pass — `simplify`, `simp`, `trigrat`, `radical_simplify`,
  `simp_algebraic_cuberoot`, `intrat`, `rat`, `integrals`, `integrate_*`
  (linrad/quadrad/linratiorad/dispatch/derivdivides/jeffrey/unknown),
  `intrischnorman`, `factorial_simplify`, `logexp_simplify`, `invtrig_simplify`,
  `simp_log`, `rationalize`; `fullsimplify_corpus` 9/9.
- **Pre-existing FAILs (proven, not regressions):**
  - `simplify_tests` 2 soft-assert FAILs (`Sqrt[x^2+6]`/`Sqrt[6]` multi-surd) —
    one symbolic radical base only (`6` is constant), so the pass is inert
    (n=1 → NULL). Confirmed by neutralising the hook and rebuilding:
    identical 2 FAILs.
  - `integrate_jeffrey_tests` 1 soft-assert FAIL (`Cosh[x] Cosh[2 x]`
    Weierstrass) — no radicals (n=0 → inert); known Limit/Weierstrass limit.
  - `intrat_corpus`/`fullsimplify_corpus` initial runs only failed to *load*
    their `.m` data via a `../`-relative path (cwd issue), not a result DIFF.
- **Out of scope (intentional):** a direct two-generator sum like
  `a^(1/3)/(a+b x)^(2/3) + (a+b x)^(1/3)/a^(2/3)` is left unchanged — there the
  free variables `b, x` are absorbed into the generators so no usable relation
  survives; combining it would require changing `Together`, which the user
  scoped out.

## Key gotchas hit (see lessons.md)
1. Reduce numerator/denominator mod the ideal **before** rationalising; rationalising
   the un-reduced numerator detonates the multivariate GCD (runaway CPU).
2. The relation `t^3 = a+b x` is essential — plain `Cancel` after reduction cannot
   use it, leaving a `(s^2+s t+t^2)` factor; rationalise the denominator via
   `PolynomialExtendedGCD` instead.
3. FRAMING-2 (eliminate base symbols) only yields a usable relation when bases
   **nest** (inner base's symbols survive); guard relations by ring membership.
4. Prove pre-existing test FAILs by neutralising the new hook and rebuilding —
   don't `git stash` on `main` (untracked files + stray-stash hazard).

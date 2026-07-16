# Bug: `Together`/`SolveAlways` stack overflow over `Q(Sqrt d)` (d ≥ 5)

**Discovered:** 2026-07-16, while landing the Cherry `ExpIntegralEi` algebraic-constant
layer (`src/calculus/cherry_ei.c`). **Status:** OPEN (worked around in cherry_ei; see below).
**Severity:** crash (stack overflow → SIGSEGV), data-dependent.

## Symptom

`Together` (and hence `SolveAlways`, `Cancel`, `Simplify`) **stack-overflows** when a rational
function mixes a **reducible-over-`Q(Sqrt d)` denominator** `x^2 - d` with the corresponding
**algebraic linear factors** `x +- Sqrt d`, for `d >= 5`. Small radicands `d = 2, 3` do NOT
crash.

## Minimal repro

```
(* CRASH (SIGSEGV, stack overflow) for d = 5, 7, 11, ... *)
SolveAlways[
  Numerator[Together[
    1/(x^2 - 5)
    - (D[(a0 + a1 x + a2 x^2 + a3 x^3 + a4 x^4)/(x^2 - 5), x]
       + (a0 + a1 x + a2 x^2 + a3 x^3 + a4 x^4)/(x^2 - 5)
       + c1/(x + Sqrt[5]) + c2/(x - Sqrt[5]))]] == 0, x]

(* OK for d = 2, 3 (returns the solution) *)
(* OK for ALL d when the polynomial-over-(x^2-d) term is absent: *)
SolveAlways[Numerator[Together[
   1/(x^2 - 5) - (c1/(x + Sqrt[5]) + c2/(x - Sqrt[5]))]] == 0, x]   (* fine *)
```

The crash needs BOTH a `.../(x^2 - d)` rational term AND the `Sqrt d` partial-fraction terms
in the same `Together`. Either alone is fine. A `bt` shows an unwindable (overflowed) stack
crashing inside a `libgmp` division — i.e. unbounded recursion in the algebraic
GCD/normalization path, not a genuine GMP fault.

## Likely area

`Together`/`Cancel` over an algebraic extension (`src/simp/`, `flint_rational_normalize_core`
/ `flint_algebraic_field_normalize`, and the `Q(Sqrt d)` GCD path). The `d`-dependence points
at a discriminant-driven recursion that fails to bottom out for larger radicands.

## Workaround (in cherry_ei.c)

The Cherry ei engine keeps algebraic-root factors OUT of the elementary part's denominator:
`y = Y/s` with `s` = product of the **degree-1 (rational-root)** factors of `den(g) q` only
(`rational_pole_denominator`, via a `Factor` walk). An algebraic ei pole `x +- Sqrt d` is then
never reduced against a `y`-denominator carrying `x^2 - d`, so the crashing `Together`
configuration never arises. (Real algebraic ei answers are still emitted and diff-back
verified.) This is a workaround, not a fix — the underlying `Together` recursion should be
bounded.

---

## Related algebraic-normaliser gaps (found 2026-07-16, Cherry C2 stress)

These decline cleanly (no crash) but block otherwise-valid answers:

### G1. `Together` does not reduce conjugate linear factors to a rational quadratic

```
Together[1/(x^2 - 2 x - 1) - (c1/(x - 1 - Sqrt[2]) + c2/(x - 1 + Sqrt[2]))]
```
keeps a degree-4 common denominator (numerator degree 3) instead of recognising
`(x - 1 - Sqrt[2])(x - 1 + Sqrt[2]) = x^2 - 2 x - 1`, so a downstream `Solve` for
`c1, c2` is overdetermined and finds no solution. (`x^2 - d` and `x^2 - x - 1`
happen to reduce.) Blocks Cherry ei over `Q(sqrt d)` for quadratics with a
non-zero linear term / mixed algebraic+rational denominators.

### G2. `PolynomialSqrt` fails on numeric radical coefficients

```
PolynomialSqrt[x^4 + 2 Sqrt[2] x^2 + 2]   (* -> $Failed, should be x^2 + Sqrt[2] *)
PolynomialSqrt[x^4 + 2 Sqrt[a] x^2 + a]   (* -> Sqrt[a] + x^2  (symbolic works)   *)
```
`ps_is_numeric` / the content-handling in `builtin_polynomialsqrt`
(`src/poly/facpoly.c`) does not take the square root of a coefficient that is an
exact numeric radical. Blocks Cherry erf when the completing-square constant
`beta` is a numeric irrational (e.g. `E^((x^4+2)/x^2)`, `beta = +-2 Sqrt[2]`).

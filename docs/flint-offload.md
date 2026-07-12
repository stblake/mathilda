# Offloading hot loops from the evaluator to FLINT

A reusable optimisation tactic in Mathilda: when an algorithm performs many
polynomial/rational operations in a loop, **do not** implement each step by
building an `Expr` tree and calling a builtin (`Expand`, `Together`, `Cancel`,
`PolynomialGCD`, `PolynomialQuotient`, `CoefficientList`, `D`, …). Every such
call re-enters the generic evaluator, and the evaluator's per-call overhead —
attribute lookup, `Orderless` canonical sort, own-value application, symbol
interning, deep `expr_free` — dwarfs the actual arithmetic. Instead, convert to
a packed FLINT representation **once**, run the whole loop in native FLINT
kernels, and convert back **once**.

## Why the round-trip is the cost (measured)

Profiling `Integrate[((x-99)/(x+2)^102) Exp[x], x, Method->"RischTranscendental"]`
(In17) at ~17 s, the heaviest self-time functions were:

```
evaluate_step   expr_free   builtin_plus   symtab_get_def   expr_eq
intern_symbol   evaluate    apply_own_values   is_zero_poly_depth
builtin_coefficientlist   builtin_times   builtin_power
```

**None is polynomial math.** The base-field Risch differential equation was
solved by an SPDE ladder (~100 recursion levels for a degree-102 integrand),
each level doing a handful of `rde_mul`/`rde_gcd`/`rde_quot`/`rde_rem`/`rde_dx`
helpers — and every one of those built an `Expr` and evaluated it. The
constant-factor evaluator overhead, multiplied across ~100 levels × several ops
× degree-100 operands, was the entire runtime. The polynomial arithmetic itself
was noise.

## The tactic

1. **Detect the fast-path domain.** The packed types are exact and narrow:
   `fmpq_poly` (univariate over Q), `fmpq_mpoly` (multivariate over Q),
   `fmpz_mpoly_q` (rational functions over Q), number-field / tower variants.
   Gate on the operands being in that domain (e.g. "univariate in `x` with
   rational coefficients") and return a sentinel to fall back otherwise.
2. **Convert once, at the boundary.** Walk the `Expr` straight into the FLINT
   type without going through `Together`/`Numerator`/`Denominator`. See
   `expr_accum_fmpq_poly` (polynomial) and `expr_to_fmpq_ratfunc` (rational
   function) in `src/poly/flint_bridge.c` — the latter combines a `Plus` of
   fractions and `Power[·,-k]` denominators directly in `fmpq_poly`, so no
   evaluator `Together` is ever called.
3. **Run the whole loop in FLINT.** All the kernels are native and O(M(n)):
   `fmpq_poly_gcd`, `_xgcd`, `_divrem`, `_rem`, `_mul`, `_sub`, `_derivative`,
   `_degree`, `_get_coeff_fmpq`. No `Expr` is allocated inside the loop.
4. **Convert back once.** `fmpq_poly_to_expr_x` emits an unsimplified
   `Plus[Times[coeff, x^k], …]`; the evaluator canonicalises it after the
   calling builtin returns — one canonicalisation, not one per operation.
5. **Preserve semantics exactly.** A decision procedure's "no solution" verdict
   must survive: the FLINT path returns the *same* NULL/decline as the Expr
   path (it is a faithful reimplementation, not an approximation), and any
   downstream verification gate (e.g. the Risch diff-back) still runs.

## Worked example — the base-field RDE

`flint_rde_base_solve_fg` (`src/poly/flint_bridge.c`) solves the exponential
tower's base-field Risch equation `Dq + f q = g` entirely in `fmpq_poly`:
RdeNormalDenominator + RdeBoundDegreeBase + Bronstein's SPDE ladder +
PolyRischDENoCancel, with `f`, `g` converted directly from `Expr`. It is
dispatched at the top of `rde_base` in
`src/calculus/integrate_risch_transcendental.c`; the Expr implementation
(`rde_spde`, `rde_polyrischde_nocancel1`, …) remains as the fallback for
tower/parametric base fields where the coefficients carry other symbols.

Result: the base-field SPDE solve for In17 went from being the bottleneck to
**~0.002 s**, and In16 / the `poly·e^x` family from seconds to milliseconds,
with the diff-back gate confirming correctness and the non-elementary siblings
(`∫e^x/(x+2)^2` → `Ei`) still declining.

## Candidate sites (same tactic, not yet applied)

The Risch engine still has Expr-level hot spots that fit this pattern exactly:

- **`rt_exp_poly_case` Laurent splitting** — `Together`/`Numerator`/
  `Denominator`/`Coefficient`/`PolynomialQ` on the integrand's `t = e^u`
  Laurent form, whose `x`-coefficients are dense rationals. This is In17's
  *remaining* cost (~7 s): the `CoefficientList` storm visible in the profile.
  A bivariate `fmpq_mpoly` (or per-coefficient `fmpq_poly`) representation of
  the Laurent polynomial would remove it.
- Any `rt_eval1("Together"/"Cancel"/"Expand", …)` inside a loop over tower
  levels or Laurent powers.

When adding one, follow the five steps above and keep the Expr path as the
out-of-domain fallback.

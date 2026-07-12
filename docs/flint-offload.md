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

## Second worked example — `Together` on a plain rational

After the RDE core was offloaded, In17's residual (~4-7 s) turned out to be a
*different* consequence of the same problem: `rt_exp_poly_case`'s Laurent split
and `rde_base`'s zero-test called `Together` (and `rt_is_zero`, which is built
on `Together`) several times on the degree-102 rational, and `Together` itself
cost ~0.5 s each. Two fixes, both instances of the tactic:

- **`flint_rational_together`** (`src/poly/flint_bridge.c`) — `Together` for a
  plain rational function over Q: convert to `fmpz_mpoly_q` (which stores the
  fraction reduced automatically) and read back num/den. Dispatched in
  `builtin_together_compute` after the algebraic/parametric normalizers decline.
  Gated to fire only when the input actually has a denominator (a
  denominator-free product stays factored) and only for genuine rationals (a
  symbolic/fractional power such as `a^x` or `Sqrt[2]` makes the conversion
  decline). Output matches the classical `Together` exactly (expanded, reduced).
  `Together[(x-99)/(x+2)^102]`: **0.47 s → 0.0016 s** (~290×).

- **structural `rt_is_zero`** (`integrate_risch_transcendental.c`) — a value is
  zero iff its numerator is; a recursion over `Times`/`Power` never combines
  over the big denominator (a reciprocal factor `1/g` is never zero), so
  `Together` is confined to the small numerator polynomials, not the degree-2n
  denominator.

Net: In17 **~17 s → 0.66 s** (~26×), and the whole `Together`-heavy path across
Mathilda (Simplify, Apart, every integrator) benefits from the same fast path.

## More candidate sites (same tactic)

- Any remaining `rt_eval1("Cancel"/"Expand", …)` inside a loop over tower levels
  or Laurent powers.
- `Cancel` could take the same plain-rational `fmpz_mpoly_q` fast path
  (`flint_rational_together` is essentially a reduced-fraction normaliser that
  Cancel could reuse).

When adding one, follow the five steps above and keep the Expr path as the
out-of-domain fallback.

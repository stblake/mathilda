# FLINT-native fast paths — session writeup

This document records a body of work that moved several hot polynomial /
rational routines off Mathilda's generic evaluator and onto native FLINT
kernels. The motivating trail was the transcendental Risch integrator, but the
fixes landed in core builtins (`Cancel`, `Together`, `PolynomialGCD`,
`PolynomialQuotient`, `Expand`) and speed up the whole system.

The reusable *tactic* is documented separately in
[`docs/flint-offload.md`](docs/flint-offload.md); this file is the concrete
inventory of what was built, why, and the measured effect.

---

## 1. The core insight

Profiling `Integrate[((x-99)/(x+2)^102) Exp[x], x,
Method->"RischTranscendental"]` (call it **In17**) at ~17 s, the heaviest
self-time functions were, every time:

```
evaluate_step   expr_free   builtin_plus   symtab_get_def   expr_eq
intern_symbol   evaluate    apply_own_values   is_zero_poly_depth
builtin_coefficientlist   builtin_times   builtin_power
```

**None is polynomial arithmetic.** The cost was the *constant factor* of routing
each intermediate polynomial back through the evaluator — attribute lookup,
`Orderless` canonical sort, own-value application, symbol interning, deep
`expr_free` — multiplied across an `O(n)`-deep algorithm on degree-100 operands.

The fix pattern (see `docs/flint-offload.md` for the full recipe): detect the
packed domain, convert `Expr → FLINT` **once**, run the whole loop in native
FLINT kernels, convert back **once**, and keep the Expr path as the
out-of-domain fallback.

---

## 2. What was built

| Area | Function(s) | File | Commit |
|------|-------------|------|--------|
| Base-field RDE | `flint_rde_base_solve_fg`, `expr_to_fmpq_ratfunc`, `expr_accum_fmpq_poly`, `fq_spde`, `fq_polyrischde_nocancel1`/`_integrate` | `src/poly/flint_bridge.c` | `9316eb5` |
| `Together` (plain rational) | `flint_rational_together` | `src/poly/flint_bridge.c` | `3204f87` |
| `Cancel` (plain rational) | `flint_rational_cancel` (+ shared `flint_rational_normalize_core`) | `src/poly/flint_bridge.c` | `a4dc856` |
| `Cancel` exact division | FLINT path in `cancel_exact_div_wrapper` / `_strict` | `src/rat.c` | `b8eb456` |
| `PolynomialGCD`/`Quotient`/`Remainder`, `Expand` | FLINT gating on dense high-degree input | `src/poly/poly.c`, `src/expand.c` | `e248cf3` |
| `rt_is_zero` | structural (numerator-only) zero test | `src/calculus/integrate_risch_transcendental.c` | `3204f87` |

All FLINT entry points are `#ifdef USE_FLINT`; each has a no-FLINT stub that
returns `NULL`/`-1` so the build degrades to the classical path.

---

## 3. Component detail

### 3.1 Native base-field RDE (`9316eb5`)

The base-field Risch differential equation `Dq + f q = g` (the core of the
`poly·e^x` case) was solved by an SPDE ladder whose every `rde_mul`/`rde_gcd`/
`rde_quot`/`rde_rem`/`rde_dx` step round-tripped through the evaluator.
`flint_rde_base_solve_fg` runs the *whole* base-field solve — RdeNormalDenominator
+ RdeBoundDegreeBase + Bronstein's SPDE ladder + PolyRischDENoCancel — in
`fmpq_poly`, converting `f` and `g` straight from `Expr` (no evaluator
`Together`). Dispatched at the top of `rde_base`; the Expr implementation stays
as the fallback for tower/parametric base fields.

Key correctness point: for the exponential tower `f = i·u'` is always a
*polynomial* (`u` a polynomial in `x`), so `den(f) = 1` and WeakNormalizer is
**provably** a no-op — the native path asserts that (`FD == 1`, else it declines
to the Expr path, which is the primitive/log case).

- In16 `((x-100)/x^101) e^x`: **4.8 s → 0.006 s**
- `(x^101+1) e^x`: **→ 0.002 s**
- In17 base-field SPDE solve: bottleneck → **0.002 s**

### 3.2 `Together` / `Cancel` for plain rationals (`3204f87`, `a4dc856`)

After the RDE core was native, In17's residual was `Together` (and
`rt_is_zero`, built on it) called several times on the degree-102 rational at
~0.5 s each — `Together` itself was slow. `flint_rational_together` combines a
plain rational over Q into a single reduced fraction via `fmpz_mpoly_q` (stored
in lowest terms automatically). `flint_rational_cancel` reuses the same core.

The subtlety is the **per-builtin structural gate**, not the arithmetic:

- Both require an actual denominator (`expr_has_denominator`) — a
  denominator-free product/polynomial is left *factored*, matching the classical
  builtins (`Together[(x+1)(x+2)]` must **not** expand).
- Both decline symbolic/fractional powers (`a^x`, `Sqrt[2]`), which
  `expr_to_mpolyq` rejects, so those keep their existing paths.
- `Cancel` additionally requires **no denominator inside a `Plus`**
  (`denom_inside_plus`), because `Cancel` — unlike `Together` — leaves a sum of
  fractions uncombined (`Cancel[1/(x+1)+1/(x+2)]` stays `1/(1+x)+1/(2+x)`).

Output matches the classical builtins **exactly** (both expand num/den when a
denominator is present; verified across ~20 fraction / product / sum cases).

- `Together[(x-99)/(x+2)^102]`: **0.47 s → 0.0016 s** (~290×)
- `Cancel[(x-99)/(x+2)^102]`: **0.24 s → 0.0013 s** (~185×)

These are system-wide — `Simplify`, `Apart`, and every integrator that combines
fractions ride the same fast path.

### 3.3 `Cancel` exact division (`b8eb456`)

`Cancel` already computed the num/den GCD via FLINT, but the subsequent exact
divisions still used the classical `exact_poly_div` (dense `Expr` arithmetic).
Both exact-division sites in `src/rat.c` now try FLINT's `fmpq_mpoly_divides`
first, gated on a max-integer-exponent test (`≥ 32`) so low-degree/factored
outputs are undisturbed.

- `Cancel[Expand[(x+3)^120 (2x+5)^60]/Expand[(2x+5)^60]]`: **70 s → 1.6 s** (~44×)
- multivariate analogue: **2.2 s → 0.2 s** (~11×)

### 3.4 Core-builtin fixes (`e248cf3`)

Two latent defects surfaced by the dense high-degree cases: `PolynomialGCD`/
`Quotient`/`Remainder` returned wrong results / ran for tens of seconds on dense
large-coefficient polynomials (int64 overflow) → routed through FLINT
(`fmpq_mpoly`), gated on exponent ≥ 32; and `Expand` had an arbitrary
`exponent < 100` ceiling → replaced with a result-size estimate.

### 3.5 Structural `rt_is_zero` (`3204f87`)

`rt_is_zero` combined the whole expression with `Together` to test for zero. A
value is zero iff its *numerator* is; a structural recursion over `Times`/`Power`
(a product is zero iff a factor is zero; a reciprocal `1/g` is never zero)
confines `Together` to the small numerator polynomials, never the degree-2n
denominator. Only a `Plus` (whose terms might cancel) falls back to the exact
test.

---

## 4. Net effect

| | before | after |
|---|--------|-------|
| **In17** `∫((x-99)/(x+2)^102) e^x` | ~17 s | **~0.5 s** (~34×) |
| In16 `∫((x-100)/x^101) e^x` | 4.8 s | 0.006 s |
| `Together[(x-99)/(x+2)^102]` | 0.47 s | 0.0016 s (~290×) |
| `Cancel[(x-99)/(x+2)^102]` | 0.24 s | 0.0013 s (~185×) |
| `Cancel[Expand[(x+3)^120(2x+5)^60]/…]` | 70 s | 1.6 s (~44×) |
| `integrate_risch_transcendental_tests` | ~80 s | ~15 s |

---

## 5. Correctness methodology

- **Decision-procedure semantics preserved.** The FLINT paths are faithful
  reimplementations, not approximations. The Risch "no solution" verdict still
  holds — `∫e^x/(x+2)^2`, `∫e^x/(x-3)^2` (non-elementary `Ei`) still decline —
  and every returned antiderivative is still diff-back verified
  (`Simplify[D[∫f]-f]==0`).
- **Output form preserved.** For `Cancel`/`Together` the FLINT output was checked
  against the classical output on a spread of fraction / product / sum / edge
  cases before shipping; the per-builtin gates exist precisely to preserve the
  cases where the classical builtin leaves things factored or uncombined.
- **No regressions, verified by A/B.** Every "new" test failure encountered was
  A/B-checked against the committed baseline and found pre-existing (simplify 4
  radical-Sqrt cases; intrat 2; intrat_corpus 5 BronsteinRational/radical;
  intrischnorman 2; parfrac 2). Suites confirmed green: rat, together,
  fullsimplify, expandfrac, both risch, goursat, derivdivides, trigrat,
  sum_rational, radical_simplify.
- **Memory-clean.** Valgrind on representative cases sits at the macOS `Sin[1.0]`
  dyld/Accelerate baseline (13,440 B / 420 blocks); no allocations traced to the
  new code.

---

## 6. Candidate future work

Same tactic, not yet applied:

- **`Apart`** and other `PolynomialGCD` / `Resultant`-heavy passes that still
  route through the classical multivariate Euclid on high-degree inputs.
- Any remaining `rt_eval1("Expand", …)` inside a loop over tower levels or
  Laurent powers in the Risch engine.

When adding one, follow the five-step recipe in `docs/flint-offload.md` and keep
the Expr path as the out-of-domain fallback; the correctness work is the
per-builtin domain/structural gate, not the FLINT arithmetic.

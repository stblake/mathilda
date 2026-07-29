# RDE SPDE optimisation — scope

**Goal:** cut the transcendental Risch RDE solver's cost on high-degree
polynomial×exponential integrands (In16/In17 and the `poly·e^x` family) from
seconds to milliseconds, without touching the decision-procedure correctness
guarantees.

**Measured target:** In17 `∫((x−99)/(x+2)^102) e^x dx` — currently ~17.5 s,
target < 0.5 s. Ripple effect: the `integrate_risch_transcendental_tests`
binary drops from ~80 s back under the 60 s harness alarm.

---

## 1. Root cause (measured, not assumed)

`Integrate[..., Method->"RischTranscendental"]` on `p(x)·e^x` solves the Risch
differential equation `Dq + f·q = g` in the base field via Bronstein's SPDE
ladder (`rde_spde`) + `PolyRischDENoCancel` (`rde_polyrischde_nocancel1` /
`_integrate`) in `src/calculus/integrate_risch_transcendental.c`.

For In17 the base field is **Q(x)** (univariate over the rationals): `f = 1`,
`g = (x−99)/(x+2)^102`. `rde_base` forms `a = dn·h = (x+2)^101`,
`b = a − D((x+2)^101)`, `c = dn·h²·g`, degree bound `n ≈ 100`, then SPDE reduces
the bound by `deg(a1)` per level. Each level does one `PolynomialExtendedGCD`,
`PolynomialGCD`, a couple of `PolynomialQuotient/Remainder`, `Expand`-based
`mul/add/sub`, and `D`, on degree-~100 polynomials — **~100 levels deep**.

`sample` of the live 17.5 s run (`/tmp/samp2.txt`), heaviest self time:

```
 366 evaluate_step        204 builtin_plus       113 is_polynomial
 296 expr_free            193 symtab_get_def      111 is_zero_poly_depth
 189 intern_symbol        190 expr_eq             99 builtin_times
 176 evaluate             123 apply_own_values    82 expr_expand_patt
                                                   69 builtin_coefficientlist
```

**None of the top functions is polynomial math.** The entire cost is the
*constant-factor overhead of routing every intermediate polynomial back through
the generic evaluator*: each `rde_mul/rde_add/rde_sub/rde_quot/rde_rem/rde_gcd/
rde_dx/rt_degree/rt_is_poly/rt_is_zero` builds an `Expr` tree and calls
`Expand`/`Plus`/`Times`/`PolynomialGCD`/`PolynomialQuotient`/`CoefficientList`/`D`
through `evaluate()` — attribute lookup, Orderless canonical sort, own-value
application, `intern_symbol`, deep `expr_free`. Multiply by ~100 levels ×
several ops × degree-100 operands and that is the 17.5 s.

The FLINT GCD/exact-division fast paths already added help the *individual*
`PolynomialGCD`/quotient calls, but cannot remove the per-op `evaluate()`
round-trip, the `Expand`-based `mul/add/sub`, or the `CoefficientList`-based
degree/zero tests. The win requires taking the *whole ladder* off the evaluator.

## 2. Fix — run the SPDE ladder natively over `fmpq_poly`

FLINT's univariate rationals `fmpq_poly` provide every operation the ladder
needs as a native O(M(n)) kernel: `fmpq_poly_gcd`, `_xgcd`, `_divrem`, `_rem`,
`_mul`, `_sub`, `_add`, `_derivative`, `_degree`, `_get_coeff_fmpq`,
`_is_zero`. (Confirmed present in flint 3.6.0.)

**Plan:** port the hot core — `rde_spde`, `rde_polyrischde_nocancel1`,
`rde_polyrischde_integrate` — to `fmpq_poly`, converting Expr→`fmpq_poly` once
on the way in and `fmpq_poly`→Expr once on the way out. The ~100-level ladder
then runs entirely in packed FLINT arithmetic: microseconds, no evaluator.

### Dispatch / scope boundary
Engage the native path in `rde_base` **iff `a`, `b`, `c` are univariate
polynomials in `x` over Q** — i.e. the Expr→`fmpq_poly` conversion succeeds
(returns NULL on any non-x symbol, tower variable, or non-rational/algebraic
coefficient). This is exactly the slow family:

| Case | base field | native path? |
|------|-----------|--------------|
| In16 `((x−100)/x^101) e^x` | Q(x) | **yes** (fast) |
| In17 `((x−99)/(x+2)^102) e^x` | Q(x) | **yes** (fast) |
| `(x^101+1) e^x`, `((x^101+1)/(x+1)) e^x` | Q(x) | **yes** |
| multi-pole product · e^x, `x^100/100` case | Q(x) | **yes** |
| In2/In5/In8/In9 (nested/log/exp towers) | C(x, t₀…) multivariate | no — keep Expr path (already fast, small degree) |

The existing Expr `rde_spde`/`rde_polyrischde_*` stay as the fallback for
multivariate/tower base fields; the native path is a drop-in for the univariate
common case only.

### What stays on the Expr path (deliberately)
`rde_weak_normalizer` and the RdeNormalDenominator prep in `rde_base` are
**one-shot O(1) passes**, not in the ~100-level loop — the profile shows the
cost is the ladder, not the prep. Port the ladder first; leave prep on Expr.
`a`, `b`, `c` are already materialised as Expr polynomials at the SPDE call
site, so conversion happens there.

## 3. Correctness (unchanged guarantees)

- **Decision procedure preserved.** `fmpq_poly` SPDE/PolyRischDE are exact; they
  return the same "no solution" verdict when no bounded-degree `q` exists, so
  `∫e^x/(x+2)^2` and `∫e^x/(x−3)^2` (genuinely non-elementary, `Ei`) still
  decline. The native path returns NULL in exactly the cases the Expr path does.
- **Diff-back gate still runs.** `rt_verify_antideriv` (Simplify-to-0) is
  unchanged and independently re-checks every returned antiderivative.
- **No arbitrary caps** (project rule): the degree bound is still the derived
  `RdeBoundDegreeBase` value; the native ladder just executes the same
  recursion faster.
- **Round-trip fidelity:** Expr↔`fmpq_poly` must be exact; unit-tested directly.

## 4. Work items

1. **Conversion helpers** (new, in `src/poly/flint_bridge.{c,h}` under `USE_FLINT`):
   `Expr* flint_expr_to_fmpq_poly(const Expr* e, const char* xvar, fmpq_poly_t out)`
   returning success/fail, and `flint_fmpq_poly_to_expr(const fmpq_poly_t p, const char* xvar)`.
   (Coefficient extraction can reuse the existing `to_mpoly` var machinery or a
   direct `CoefficientList`-free walk.)
2. **Native ladder** (new, either a new file `src/calculus/integrate_risch_rde_flint.c`
   — needs a `tests/CMakeLists.txt` COMMON_SRC line — or inline in the existing
   engine under `#ifdef USE_FLINT` to avoid the CMake edit):
   `rde_spde_fmpq`, `rde_polyrischde_nocancel1_fmpq`, `rde_polyrischde_integrate_fmpq`,
   mirroring the Expr versions one-for-one with `fmpq_poly` ops.
3. **Dispatch** in `rde_base` (§2): try converting `aa,bb,cc`; on success run the
   native ladder and convert `q` back; else fall through to the current code.
4. **Tests:** `test_bronstein_rde_examples` already asserts correctness for all
   these integrals (numeric single-integration diff-back). Add a machine-
   independent perf gate (In17 normalised against a calibration op, mirroring
   `bench_assoc.c`) so the ladder blow-up can't silently return. Restore the
   test alarm to 60 s once In17 is sub-second.

Estimated size: ~300–400 LoC (conversion + three ported functions + dispatch).

## 5. Alternative considered — and why rejected

*Fast-path each `rde_*` helper to bypass the evaluator* (do the poly arithmetic
directly without the `Expand`/`Together` round-trip, keep the Expr rep). This
removes the `builtin_*`/`Orderless`-sort overhead but **still allocates and
frees an Expr tree per operation per level** (`expr_free` is already the #2 self
cost) and still walks Expr for degree/zero tests. `fmpq_poly` removes all of it
in one representation change and is the natural substrate for univariate-over-Q.

## 6. Expected outcome

- In16/In17/`poly·e^x` family: seconds → milliseconds (100–1000×).
- `integrate_risch_transcendental_tests`: ~80 s → ~30 s; alarm back to 60 s.
- Nested-tower examples (In2/In5/In8/In9): unchanged (Expr fallback).
- Correctness and "no solution" verdicts: unchanged; diff-back gate intact.

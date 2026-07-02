# Plan: closing the Goursat named integrals fast

**Goal.** The two classic Goursat pseudo-elliptic named integrals must close
end-to-end, fast (Mathematica does integral 1 in **0.37 s**):

```
Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[x (1-x)(1-k x)]), x]        -> -2 ArcTan[...]
Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[(1-x^2)(1-k^2 x^2)]), x]    -> ArcTanh[...]
```

The integrand ring is `Q(a,b,k)(√k)[x]` (parameters a,b,k; a single symbolic
radical √k). The V4 (Goursat) descent reduces each to an inner
`∫ G(u)/√(quadratic(u)) du` over that field.

## Status (2026-07-01)

**Done — every polynomial-algebra layer is FLINT-fast over `Q(a,b,k)(√k)`.**
All routed through `src/poly/flint_bridge.c` via the `√k → S, k → S²` collapse:

| Layer | Route | Notes |
|-------|-------|-------|
| `PolynomialGCD` | `flint_extension_gcd` (`fmpq_mpoly_gcd`) | |
| `Cancel` / `Together` (1- & 2-arg `Extension`) | `flint_cancel_fraction` + `flint_parametric_field_normalize` | the latter (fmpz_mpoly_q) handles a **Plus of fractions**, which the num/den path dropped to the slow QA path (3.5 s → 0.4 ms) |
| `Resultant` | `flint_parametric_sqrt_resultant` (`fmpq_mpoly_resultant`) | 90 s → 0.5 ms |
| `Factor` / `FactorSquareFree` | `flint_parametric_sqrt_factor` (`fmpq_mpoly_factor`) | parametric field IS a rational function field, so factoring works (unlike constant-radicand number fields) |
| `Apart` | Bézout modular-inverse (`parfrac.c`) | replaced the exponential `RowReduce`; handles the ArcTan quadratics |
| `PolynomialExtendedGCD`, `PolynomialQuotient`/`Remainder` | `flint_parametric_field_{xgcd,divrem}` (`gr_poly` over `fmpz_mpoly_q` = `Q(a,b,S)[x]`) | field coeffs incl. denominators (`allow_neg_pow`) |

**Result:** the **inner integral now closes** (`ArcTan`), where it previously
hung. Typical parametric integrals close cleanly and sub-second
(`∫1/((x-√k)(x+√k)) = -ArcTanh[x/√k]/√k`, etc.).

**Remaining gap:** the inner integral is **~17.5 s** (was ~200 s; the full named
integral runs it ~3× plus descent overhead). Target is ~0.37 s, so **~50× more**
is needed. It is no longer a hang or a single slow primitive — it is
*cumulative* cost across many operations on the large `QuadraticRadicals`-
rationalised intermediate. Profiling (`intrat_trace` phase timing) attributes the
17.5 s roughly: traced LRT phases (`IntRationalLogPart`, `LogToReal`,
`LogToAtan`) ~4 s; ~13 s untraced in BronsteinRational's non-LRT steps
(`intrat_integrate_summands`, the six `intsimp_*` post-passes on the huge result).

## Next modular improvements (in priority order)

Each is independently testable and independently landable.

1. **Profile & attack the ~13 s untraced BronsteinRational post-processing.**
   Add step timing around `intrat_integrate_summands` (intrat.c:3325) and the
   `intsimp_distribute_plus` / `intsimp_strip_log_constants` /
   `intsimp_log_to_arctanh` / `intsimp_normalize_inverse_trig_signs` chain
   (intrat.c:3339-3353). These walk/rewrite a multi-thousand-term result;
   likely quadratic in expression size. Fix: operate on the already-canonical
   FLINT form, or skip redundant re-distribute/re-strip passes.

2. **Route the remaining `intsimp_*`/`intrat_canonic` radical simplifications
   through `flint_parametric_field_normalize`.** Any place that Cancels/Togethers
   a scalar field expression (log arguments, `A²+B²`, residues) should use the
   fmpz_mpoly_q normalizer, not classical Simplify/Cancel. Grep intrat.c/intsimp.c
   for `intrat_canonic`, `internal_cancel`, `internal_together`, `Simplify` on
   radical-bearing scalars.

3. **Avoid the degree-doubling `QuadraticRadicals` rationalisation for the
   `∫rational/√quadratic` case.** The Euler substitution turns a degree-2 radical
   problem into a **degree-4 rational** with `k^(1/2)…k^(9/2)` coefficients — the
   root of the expression swell. A direct `∫rational/√quadratic → ArcTan/ArcTanh`
   reduction (Hermite for the √quadratic + one arctan/log per simple pole) stays
   low-degree. This is the biggest structural win but the largest change; do it
   after 1–2 exhaust the cheap wins.

4. **Constant-algebraic parametric regime (`fq_nmod` outer loop).** For
   `Q(t)(α)` with `α` a genuine algebraic number over `Q(t)` (constant-coefficient
   minimal polynomial defining `GF(p^d)`) — not needed for √k, but the last
   general regime. Prime selection + `fq_nmod_mpoly_gcd` + CRT/reconstruction per
   ALGEBRAIC_EXTENSION_ARITHMETIC_PLAN.md §3.4.

5. **Number-field factoring over towers** (`Q(√2,√3)` etc.): FLINT's `gr` layer
   is `GR_UNABLE`; stays in `qafactor.c` (Trager). Only blocks Apart-tower /
   Integrate-tower for *constant* radicands, not the parametric target.

## Acceptance

- Both named integrals close under the normal `TimeConstrained` budget with
  `D[result,x] − integrand ≡ 0` (verify via `PossibleZeroQ` at sample points —
  the answers land in nested-radical `ArcTan` forms `Cancel` cannot reduce to 0).
- Inner integral wall-clock within a small factor of Pari/Mathematica.
- All regressions green; valgrind-clean; `USE_FLINT=0` fallback stays correct.

## Test coverage added this pass

`tests/test_flint_bridge.c`: `test_parametric_field_ops` (xgcd Bézout identity
incl. the degree-2 `k^(3/2)` case, divrem round-trip, field-normalize of a Plus
of fractions), `test_parametric_resultant`, `test_parametric_factor`,
`test_goursat_descent_operands` (14 real descent Cancel/Together operands).
`tests/test_extension_auto_builtins.c`: `test_apart_parametric_sqrt_k`.
`tests/test_poly.c`, `tests/test_qafactor.c`: radical-coefficient div/factor.

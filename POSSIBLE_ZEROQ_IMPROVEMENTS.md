# PossibleZeroQ / `zero_test` — improvement backlog

Follow-up items for the hybrid symbolic-numeric zero recogniser
(`src/zero_test.c`, `src/zero_test.h`; design in `ZERO_RECOGNISE_PLAN.md`).
Each entry: a minimal repro, the diagnosis, the observed impact, and a
concrete direction. These are **deferred** — logged here so a later session can
pick them up without re-deriving the analysis.

---

## 1. Spins on a Gaussian × Erf residual (the exact-ODE integrating-factor form)

**Status:** RESOLVED (2026-09-08). Fixed by the core recogniser change described in
the "Suggested direction" below: an exponential-combining normalisation (Stage 0.5)
runs `ExpandAll` on any input carrying a constant-base exponential with a *non-linear*
symbol-dependent exponent, collapsing same-base exponentials before the numeric ladder.
`PossibleZeroQ[res]` now returns `True` at machine speed (the flagship residual reduces
to literal `0`). The gate is narrowed to *non-linear* exponents so an affine `E^x`
(representable, ladder-friendly) is left untouched — a blanket gate regressed a
constant-coefficient-ODE verify whose residual has an incidental affine `E^x` and a trig
core that `ExpandAll` distributed into a sampler-hostile form. The local DSolve-side
workaround in `dsolve_verify_body` was removed as redundant. See
`src/zero_test.c` (`exp_exponent_is_nonlinear`, `expr_has_symbolic_exp_kernel`,
`zt_normalize_exp_kernels`), `tests/test_zero_test.c` (Group 16),
`tests/test_possiblezeroq_expcombine_stress.c`, and the 2026-09-07 changelog.

Note the pre-existing (unrelated) sampler limitation this exposed: a genuine
non-identity that differs from a zero residual only *inside* the `E·Erf[imaginary]`
tiny·huge terms (e.g. `D[b,x]-b`) is still mis-classified `True` — the imaginary-argument
`Erf` overflows to `Inf` at the sampler's moderate range, degrading to `UNKNOWN → True`.
That is independent of this fix (which only ever produces a value-equal form).

_Original report (retained for context):_

### Repro (minimal)

```
b = E^(-1/2 x^2) (C[1] + I (C[2] Sqrt[Pi] Erf[-(I x)/Sqrt[2]]) / Sqrt[2]);
res = D[b, {x, 2}] + x D[b, x] + b;        (* identically 0 — b solves y''+x y'+y==0 *)

PossibleZeroQ[res]        (* SPINS: > 20 s, floods "$IterationLimit exceeded", never usable *)
Simplify[res]            (* 0, fast *)
ExpandAll[res] // PossibleZeroQ   (* True, fast *)
```

Contrast — these all decide **fast** (so it is not "Erf" per se):
```
PossibleZeroQ[Erf[I x] - Erf[I x]]                                    (* True *)
PossibleZeroQ[D[Erf[x], x] - 2/Sqrt[Pi] Exp[-x^2]]                    (* True *)
PossibleZeroQ[D[Erf[-I x/Sqrt[2]], x] - 2/Sqrt[Pi] Exp[x^2/2]*(-I/Sqrt[2])]  (* True *)
```

### Diagnosis

The residual is a sum of products in which `E^(-x^2/2)` (from `b`) multiplies
`E^(x^2/2)` (from differentiating `Erf[-I x/Sqrt[2]]`, whose derivative is
`∝ E^(-(-I x/Sqrt2)^2) = E^(x^2/2)`). These two powers of `E` are **not adjacent
factors** — they sit in different summands — so `Times` never combines them to
`E^0 = 1`. `x` is free, so `zero_test_decide` reaches Stage 3 (Schwartz–Zippel):
for each of k=4 random rationals it runs the Stage 2 precision ladder
(machine → 200 → 500 → 1000 bits). At a random `x`, `E^(-x^2/2)·E^(x^2/2)`
numericalises as `tiny · huge`, i.e. **catastrophic cancellation**, so the ladder
keeps climbing to 1000 bits, and numericalising `Erf` at 1000 bits (× several
points × 4 samples) is what floods `$IterationLimit` and blows the time budget.
The expression *is* zero; the recogniser just cannot afford to confirm it in this
form.

`ExpandAll` fixes it because distributing the sums makes the `E`-powers adjacent
factors, which `Times` then collapses to `E^0` — removing the cancellation before
the numeric ladder ever runs. A bare `E^p_ E^q_ :> E^(p+q)` rule does **not**
help (the powers are not adjacent until the sums are distributed).

### Impact

- **DSolve `ExactODE`** (`dsolve_verify_body` in `src/calculus/dsolve_common.c`):
  §2.2.5 case 428 (`y''+x y'+y==0`) produces the correct Erf closed form but the
  verify's `zero_test_decide` on the 2nd-order residual spins → whole solve times
  out. The existing Erf guard in `dsolve_verify_body` (the `Power[E, ztexp_] :>
  Power[E, Expand[ztexp]]` re-check) only fires on the `ZERO_TEST_FALSE` path — it
  cannot help here, because the *first* `zero_test_decide` call never returns in
  time.
- More generally, any solution whose residual is a Gaussian × (Erf/Erfi/Fresnel)
  product hits this (integrating-factor solutions of `y' + x y == g`, parabolic-
  cylinder normal forms, etc.).

### Suggested direction (core fix, deferred)

Add an **exponential-combining normalisation to Stage 1** of `zero_test`, before
the numeric ladder: distribute products/powers enough to combine same-base
exponentials (`E^a · E^b → E^(a+b)`), e.g. an `ExpandAll`-class pass gated to
inputs containing `Power[E, _]` (or `Exp`). This is cheap for the elementary
inputs Stage 1 already Expands, and it removes exactly the tiny·huge cancellation
that defeats Stage 2/3. Alternatives worth weighing:
- Cap the Stage 2 precision ladder (or Stage 3 sample count) when the integrand
  numericalisation cost per point is high, returning `UNKNOWN` rather than
  climbing to 1000 bits — callers that only need "not-decidably-nonzero" (like
  `dsolve_verify_body`) then proceed cheaply.
- A `zero_test` entry point / flag for "quick-reject only": return `FALSE` only on
  a confident machine-precision nonzero, else `UNKNOWN` fast — for verify callers
  that reject on `FALSE` and keep on `TRUE`/`UNKNOWN`.

### Cross-references

- `src/zero_test.c` (Stages 1–3), `src/zero_test.h`, `ZERO_RECOGNISE_PLAN.md`.
- `src/calculus/dsolve_common.c` — `dsolve_verify_body` (the affected caller and
  the site of the M25 local workaround).
- DSolve §2.2.5 corpus cases 428, 482 (`DSolve_test_status/DE_examples_225.m`).

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

---

## 2. `Simplify` hangs on a sum of exponentials with widely-separated real rates

**Status:** WORKED AROUND (2026-09-09, M27); the **`PossibleZeroQ` overflow-abort
subset is now RESOLVED** (2026-09-26, v0.208); a residual `Simplify`-search
deficiency remains OPEN.

**Resolved (v0.208): the wrong/flaky `True` in `PossibleZeroQ`.** Investigation
showed the "widely-separated rates" family was not only slow but produced a
**wrong, draw-order-dependent `True`** for nowhere-zero sums like
`E^(-10 t) + E^(-100 t)`: they reach the Schwartz–Zippel sampler (the Stage-0b
certificate does not fire on a `Plus`), a sample at negative `t` overflows `E` to
IEEE `±Inf`, and the sampler **aborted the whole test on that first `UNKNOWN`**
→ `UNKNOWN` → `True`. Fix (`src/zero_test.c`): distinguish an IEEE overflow
(`ok && !isfinite(mag)` in `evaluate_rung`) from a symbolic residue and **re-draw**
the overflow point from a shrinking magnitude shell (`2^6 → … → 2^0`, `|value| >= 1`
floor kept), skipping only if every shell overflows. `PossibleZeroQ[E^(-10 t) +
E^(-100 t)]` and `PossibleZeroQ[Gamma[x+1] - x Gamma[x] + E^(x^2)]` now decide
`False` at machine speed. Tests: `test_zero_test.c` Group 18. See the 2026-09-21
changelog. **The M27 `dsolve_linsys_tidy` `Expand`-not-`Simplify` workaround stays**
(harmless, and the underlying `Simplify` search — separate from `PossibleZeroQ`,
which must never call `Simplify` — is unchanged).

_Original report (retained for context):_

**Minimal repro:** `Simplify[Exp[-10 t] + Exp[-100 t]]` does not return (an 8 s
`TimeConstrained` aborts it); `Simplify[Exp[-t] + Exp[-2 t]]` returns instantly.

**Diagnosis.** `Simplify` invokes a zero-test (`PossibleZeroQ` / `simp_search`
equivalence checks) that numericalises the two terms at sample points. `E^(-10 t)`
against `E^(-100 t)` spans a huge dynamic range at a generic `t` (e.g. at `t=1.1`,
`e^{-11}` vs `e^{-110}` — a ratio of ~10^43), so any difference/equivalence probe
sees `small ± tiny` and the precision ladder climbs to its ceiling without deciding.
The separation, not the magnitude, is the trigger: rates within ~1 order (e.g.
-1,-2) never provoke it.

**Impact.** Constant-coefficient linear ODE *systems* whose spectrum is real but
spread — e.g. `x'=-50x+20y, y'=100x-60y` (eigenvalues -10, -100) — produce a
fundamental-matrix body that is exactly such a sum, and the tidy `Simplify` in
`dsolve_linsys_tidy` hung the whole solve (uninterruptibly). §2.2.7 corpus systems
636, 650 (and more broadly any wide-spectrum system) are affected.

**Worked around (M27):** `dsolve_linsys_tidy` (`src/calculus/dsolve_linsys.c`) now
routes *any* exponential-carrying body through `Expand` rather than `Simplify`
(previously only complex-spectrum / large bodies took `Expand`); the result is
back-substitution-verified regardless, so the cosmetic loss of combination is
harmless. This does not fix the underlying `Simplify`/`zero_test` deficiency — a
direct `Simplify` of a wide-spectrum exponential sum still hangs.

**Suggested direction (core fix, deferred).** Same family as #1: before the numeric
ladder, factor out the dominant exponential (`E^(-10 t) + E^(-100 t) = E^(-100 t)
(E^(90 t) + 1)`) or compare terms in log-magnitude so a term negligible at the
sample is not differenced against a large one; or cap the ladder / return `UNKNOWN`
when per-point dynamic range exceeds the working precision.

### Cross-references (#2)

- `src/calculus/dsolve_linsys.c` — `dsolve_linsys_tidy` (the `Expand`-not-`Simplify`
  workaround and its rationale comment).
- DSolve §2.2.7 corpus systems 636, 650 (`DSolve_test_status/DE_examples_227.m`).

---

## 3. `PossibleZeroQ[-E^(-a x)]` false positive (`True` for a nowhere-zero function)

**Status:** RESOLVED (2026-09-26, v0.204).

**Minimal repro:** `PossibleZeroQ[-E^(-a x)]` returned `True`; `-e^{-ax}` is never zero.

**Diagnosis.** The expression has free symbols, so it reached the Schwartz–Zippel sampler.
The sampler draws each symbol with `|value| ∈ [1, 64]`, so the exponent `-a x` ranges over
`±[1, 4096]`. At every sample `E^(-a x)` is therefore either an IEEE overflow (`Inf`) or an
underflow (`0`) — the sampler never obtains a finite, moderate, obviously-non-zero value, and
the near-zero points drove the "residual never shrank → lenient machine zero" branch to a
`True`. Same family as #1/#2 (a symbol-dependent exponent defeats the numeric ladder), but
here the outcome was a *wrong verdict*, not a hang.

**Fix (core, sound).** A structural non-zero certificate `provably_nonzero` (Stage 0b in
`src/zero_test.c`), run right after `decide_structural` and before the guards that route to
sampling in both `zt_decide_core` and `zt_decide_assuming_core`. The exponential is entire and
has no zeros (`E^z ≠ 0` for all complex `z`), so a product of exponentials and finite non-zero
constants is provably non-zero and settles `FALSE` with no sampling. Certified: a finite
non-zero literal; a known non-zero constant; `Exp[_]`/`Power[E, _]`; `Power[b, _]` with `b`
certified; `Times[…]` of certified factors. Sound and unconditional (holds under assumptions);
it never fires on `Plus`, a possibly-zero/infinite base, or a bare symbol, so no identity is
affected. Tests: `tests/test_zero_test.c` Group 17. See the 2026-09-21 changelog.

---

## 4. Special-function *magnitude-decline* residue defeats the sampler (`Gamma[x^2]+1` → `True`)

**Status:** OPEN (deferred from the v0.208 overflow fix, which is scoped to the IEEE-overflow
class only).

**Minimal repro:** `PossibleZeroQ[Gamma[x^2] + 1]` returns `True`; `Gamma[x^2]+1` is nowhere
zero. Same family as #2/#3 (a symbol-dependent magnitude defeats the numeric ladder), but the
failure mode is different from the IEEE-overflow one the v0.208 fix cures.

**Diagnosis.** At a sampled `x` with `x^2 > ~171`, `N[Gamma[x^2]]` does not overflow to `±Inf`
— the machine/MPFR `Gamma` **declines and returns the unevaluated head** (`N[Gamma[256.0]]` →
`Gamma[256.0]`). That surfaces as a **symbolic residue** (`is_pure_numeric == false`), which the
sampler treats as "genuinely undecidable at this point" and aborts to `UNKNOWN` → `True`. The
v0.208 overflow re-draw fires only on `ok && !isfinite(mag)` (a real numeric out-of-range value),
so it deliberately does **not** re-draw this residue class.

**Why not fixed in v0.208.** Making the shell ladder also re-draw magnitude-decline residues (so
`Gamma[x^2]+1` shrinks to a resolvable small-argument point) is sound in isolation, but it makes
the sampler *reach* points it used to abort before — which re-exposes deficiency #5 below (a true
identity then returns a wrong `False`). The two must be fixed together.

**Suggested direction (deferred).** Either (a) make `N[Gamma[big]]`, `N[Erf[I·big]]`, … overflow
to `Infinity` rather than declining (turns this into the already-handled overflow class, but is a
`src/numeric.c` change with wider blast radius), or (b) treat a *numeric-head-of-large-numeric-arg*
residue as re-drawable in the sampler — but only alongside a fix for #5.

### Cross-references (#4)

- `src/zero_test.c` — `evaluate_rung` (the `is_pure_numeric` residue path), `sz_trial_shelled`.
- `src/numeric.c` — the machine/MPFR `Gamma`/`Erf` decline path.

---

## 5. Machine-precision deep-cancellation false-negative in the Stage-3 *screen* phase

**Status:** OPEN (pre-existing; surfaced while scoping the v0.208 fix).

**Minimal repro:** the Weierstrass antiderivative round-trip
`D[Integrate[Cosh[x] Cosh[2 x], x, Method -> "Weierstrass"] /. Floor[_] -> 0, x] -
Cosh[x] Cosh[2 x]` is identically `0` (`Simplify` → `0`), but a thorough sampler can return
`False`. (Today it still returns `True` — the sampler aborts on an earlier `Indeterminate` point
before reaching the bad ones — but that is luck of draw order, not robustness.)

**Diagnosis.** The Weierstrass form is a rational function of `Tanh[x/2]` whose
`(1 - Tanh[x/2]^2)^k` denominators produce terms of magnitude `~10^15`–`10^22` at moderate
samples (`x ≈ 11`–`18`). These huge terms must cancel to `0`, but machine-precision `Tanh` is not
accurate enough for the cancellation to occur, so the residual is a *large* fraction of the
operand scale (ratio `~0.4`, far above the `2^-12` "obvious non-zero" gate). The **screen** phase
therefore reports a decisive `False` at machine precision, even though the confirm ladder (MPFR)
would resolve it to `0`. The screen trusts a machine-precision `False`; deep cancellation
(`> ~52` bits) breaks that trust.

**Suggested direction (deferred).** Verify a screen `False` against catastrophic cancellation
before trusting it — e.g. when the operand scale is huge (terms `≫ 2^40`), climb one MPFR rung and
require the residual to *persist* (not shrink) before concluding `False`; a shrinking residual is a
cancellation zero. This is a change to the core `decide_numeric` gate and must be weighed against
the perf cost on the verify-reject path (a large non-zero would then climb before rejecting), so it
was left out of the focused v0.208 overflow fix.

### Cross-references (#5)

- `src/zero_test.c` — `screen_point`, `decide_numeric` (the `ZT_OBVIOUS_NONZERO_BITS` gate).
- `tests/test_zero_test.c` — `test_battery_weierstrass_cosh_product_roundtrip` (the fragile case).

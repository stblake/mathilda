# Q(i√d) Strengthening Scope — making the declined complex-quadratic integrals evaluate efficiently

**Date:** 2026-07-21
**Context:** While fixing the Cherry Gaussian integral `∫ x⁴ e^{-c x²}(a x⁴ − b) dx`, we found
that the complex-conjugate `ExpIntegralEi` cases whose constants live in **Q(i√d)** (d not a
perfect square) were **hanging**. They are now *contained* — the Cherry engine gates the
radical-bearing complex roots out (`cherry_ei.c` `expr_has_radical`) so `Integrate` returns
unevaluated instead of hanging. This document scopes what Mathilda must strengthen to make them
**evaluate** (fast, correct) instead of decline.

Companions: `CHERRY_PLAN.md`, `CHERRY_BLOCKERS.md` (§A1 "general complex case"), `CHERRY_GAUSSIAN_PLAN.md`.

---

## 1. The declined examples (the review set)

All are single-kernel `g·E^x` integrands whose only pole is an irreducible real quadratic with
**complex-irrational** roots. Each currently returns unevaluated (clean decline, **no hang**):

| Integrand | ei-argument field | root form |
|---|---|---|
| `E^x/(x²+x+1)` | Q(i√3) | (−1 ± i√3)/2 |
| `E^x/(x²+3)` | Q(i√3) | ± i√3 |
| `x E^x/(x²+x+1)` | Q(i√3) | (−1 ± i√3)/2 |
| `(x²+1) E^x/(x²+x+1)` (Cherry **d12**) | Q(i√3) | (−1 ± i√3)/2 |
| `E^x/(x²+2x+3)` | Q(i√2) | −1 ± i√2 |
| `E^x/(x²−2x+3)` | Q(i√2) | 1 ± i√2 |
| `(3x+1) E^x/(x²+2x+3)` | Q(i√2) | −1 ± i√2 |
| `E^x/(x²+2x+7)` | Q(i√6) | −1 ± i√6 |
| `(2x+1) E^x/(x²+x+3)` | Q(i√11) | (−1 ± i√11)/2 |

**Contrast — these SOLVE** (kept working): the **Q(i)** pairs whose roots are `p + q·i` with `q`
rational (no radical): `E^x/(x²+1)`, `E^x/(x²+9)`, `E^x/(x²+2x+2)`, `E^x/(x²+2x+5)`,
`E^(2x)/(x²+1)`. And all **real-algebraic** ei constants (`E^x/(x²−2)` over Q(√2), golden ratio,
√5, …) still solve. The dividing line is exactly: **a complex root that also carries a radical**
(⇒ Q(i√d), d not a perfect square).

**The class generalises beyond Cherry.** Any integrator whose answer needs Q(i√d) constants (a
rational-function log-part with an irreducible quadratic of negative non-square discriminant, a
Risch residue over Q(i√d), etc.) hits the same wall. The Cherry set is the concrete, reproducible
sample; the fix is shared infrastructure, not a Cherry patch.

---

## 2. Root cause — sharply localised

Every one of these routes into the shared normaliser
`builtin_together`/`builtin_cancel` → `cancel_recursive` (`rat.c:711`) → generic
`builtin_polynomialgcd` (`poly.c:2602`) → `poly_gcd_internal` → **`exact_poly_div`
(`poly.c:~1300`)**, which suffers **super-exponential expression swell** (effectively
non-terminating) when the polynomial's coefficients live in Q(i√d) *and* a second generator (the
`E^x` kernel, or a substituted radical symbol) is present. Confirmed by repeated backtraces:
`exact_poly_div` spins in `expr_expand`/`estimate_terms`/`is_zero_poly`, with `var_count = 2`
(`{x, kernel-or-generator}`).

### 2.1 What we already fixed (the constant sub-case)

`poly_find_radical_gen` (`poly.c:756`) now **declines on a pure numeric constant** (no polynomial
variable): substituting a radical for a fresh generator on `(1 − i√3)/2` fabricated a bogus
polynomial that blew up the GCD. With that guard **plus** a numeric-constant early-out in
`cancel_recursive` (`rat.c`), the constant layer is healthy:

```
Simplify[(1 − I Sqrt[3])/2]   → 1/2 (1 − I Sqrt[3])      (was: HANG)
Cancel[(1 − I Sqrt[3])/2]     → 1/2 (1 − I Sqrt[3])      (was: HANG)
Apart[1/((x−a)(x−ā))]         → 1/(1 − x + x²)            (a = (1+I√3)/2)
```

This is why the Cherry **nf-fallback setup** (`center`/`hd`/`disc` = `Simplify`/`Together` of
constants) no longer hangs.

### 2.2 The remaining gap — x-dependent + transcendental kernel

The sharp, load-bearing finding:

| operation over Q(i√d), with `x` present | status |
|---|---|
| `Together[E^x/(x + (1−i√3)/2)]` | ✅ `(2 E^x)/(1 − i√3 + 2x)` |
| `Together[1/(x−a) + 1/(x−ā)]` | ✅ `(−4 + 8x)/(4 − 4x + 4x²)` |
| `Cancel[(x²+x+1)/(x−a)]` | ✅ |
| **`Simplify[E^x/(x + (1−i√3)/2)]`** | ❌ **HANG** |
| **`Simplify[1/(x−a) − 1/(x−ā)]`** | ❌ **HANG** |

So **`Together` and `Cancel` are already fine over Q(i√d) with a variable — only `Simplify`
hangs.** A `Simplify` sub-pass substitutes the Q(i√d) radical into a fresh polynomial generator
(so the working set becomes `{x, gen}` or `{x, gen, E^x}`) and *then* the multivariate GCD blows
up — whereas plain `Together`/`Cancel` never fire that substitution. (The constant guard in §2.1
does not apply here because these expressions *do* contain a genuine variable `x`.)

This is exactly the wall the Cherry engine hits: the coefficient solve and the diff-back verifier
(`rt_verify_antideriv` → `Simplify`) both call `Simplify` on x-dependent Q(i√d) expressions.

---

## 3. Where Mathilda must be strengthened (two independent tracks)

> **Update (2026-07-21): Track B LANDED — the Q(i√d) ei cases now EVALUATE.** All nine §1
> examples close to a complex-conjugate `ExpIntegralEi` pair (diff-back 0). The fix was *not*
> a full number-field GCD but a targeted, exact **diff-back zero-test**: the Cherry
> number-field fallback already solved the coefficients over `Q` (symbolic generator `chs`);
> the only remaining blocker was `rt_verify_antideriv` failing to certify `E^x·R(x) = 0` over
> `Q(i√d)`. Adding `Together[ComplexExpand[TrigToExp[diff]]]` *before* the hang-prone
> `Simplify` certifies it — `ComplexExpand` splits `I·Sqrt[d]` into real/imaginary parts over
> ℝ, where `Together` decides exactly (it is additive — real-radical diffs fall through to
> `Simplify` — and sound — `ComplexExpand` is an exact identity). The full
> `integrate_risch_transcendental` battery stays green. **Track B (a general FLINT
> number-field GCD) is therefore no longer needed for these integrals**; it remains the
> durable hardening only if the *generic* `exact_poly_div` must be made safe for arbitrary
> algebraic-coefficient inputs elsewhere.
>
> **Update (2026-07-21): Track A LANDED.** The general `Simplify`/`FullSimplify`
> hang over Q(i√d) is **fixed** — the denominator-rationalisation pass
> (`simp_rationalize.c` `denom_compute_inverse`) now bails when the denominator carries an
> explicit complex literal (it is meant for *real* radical denominators; a complex one left
> `I` in the coefficients and blew up the extended-GCD over Q(i) in `{gen, x}`). Verified:
> `Simplify[E^x/(x + (1−I√3)/2)]`, `Simplify[1/(x−a) − 1/(x−ā)]`,
> `FullSimplify[1/(x−a) − 1/(x−ā)] = (I√3)/(1−x+x²)` all return fast; real-radical
> rationalisation (`1/(1+√2)→−1+√2`) is unchanged.
>
> **The Cherry complex-ei cases still decline**, though — un-gating them re-hangs, now in the
> **diff-back verifier** `rt_verify_antideriv` (`risch_util.c`). Traced precisely: the
> number-field fallback *solves* (its residual `Together` over the symbolic generator `chs`
> and the coefficient `Solve` both complete); the answer is correct
> (`E^{−a}E^{x+a}` collapses to `E^x`); but the verifier's `Together[TrigToExp[diff]]` does
> **not** reduce the `E^x·R(x)` difference to 0 (leaves it uncombined, `zero=0`) and then the
> fallback `Simplify[diff]` hits a *different* Simplify sub-pass over Q(i√d) and hangs. So the
> remaining blocker is a **number-field-aware zero-test** for the diff-back — squarely Track B.

### Track A — the targeted fix (small, high-leverage): don't substitute radicals into a blow-up

**Status: the `denom_compute_inverse` complex-denominator guard is landed (above).** It fixed
the standalone `Simplify` hang but *not* the Cherry solve, because the diff-back verifier reaches
the blow-up through a different route. The remaining Track-A-shaped options for a full Cherry
solve without Track B:

The `Simplify` pass that substitutes a Q(i√d) radical for a polynomial generator should **decline
the substitution when it would create a multivariate GCD over a transcendental kernel** — i.e.
when the substitution introduces a generator alongside an existing variable *and* a transcendental
kernel (E^x, Log, …) is present. This mirrors the guard already inside `poly_find_radical_gen`
(`poly.c:781`, `has_other_var` → "leave B^(1/m) as an opaque coefficient so the GCD stays in
K[x]"), which today fires for polynomial second variables but not for the transcendental-kernel
case.

- **Task A1.** Identify the exact `Simplify` sub-pass. Backtraces show it is *not*
  `simp_algebraic` (bails on `contains_explicit_complex`, `simp_algebraic.c:486`) nor
  `simp_radical_rational` (needs ≥2 radicals, `radrat.c:302`). Instrument `builtin_together`/the
  simp search to capture the first substituted (`$…gen`) input; likely a radical/denominator
  canonicalisation in `simp_rationalize.c` or a `simp_search` node transform.
- **Task A2.** Guard that pass with the same predicate `poly_find_radical_gen` uses: if the
  expression carries a transcendental kernel (or the substitution yields ≥2 live generators),
  keep the radical as an opaque algebraic coefficient (K[x] over K = Q(i√d)); let `Together`/
  `Cancel` — which already work — do the cancellation.
- **Payoff.** `Simplify` over Q(i√d) with a kernel stops hanging. The Cherry diff-back verifier
  then certifies, and **all §1 examples can be un-gated and SOLVE** via the existing
  `rt_cherry_ei_conjpair_nf` number-field fallback (which already solves over Q with a symbolic
  generator — it was only ever blocked by the `Simplify` hang, not by its own logic).
- **Effort:** medium. **Risk:** medium (touches the Simplify search) — gate behind the
  kernel/second-generator predicate so ordinary radical Simplify is unchanged; verify against
  `simplify_tests`, `radical_simplify_tests`, `simp_algebraic_cuberoot_tests`.

### Track B — the robust fix (larger): exact number-field-aware GCD / cancellation

Even with Track A, the *generic* `exact_poly_div` remains a latent hazard for any input that
reaches it with Q(i√d) coefficients over ≥2 generators (other integrators, user `Simplify`). The
principled cure is to **never fall to the subresultant PRS over an algebraic-number coefficient
ring** — route through exact number-field arithmetic instead.

- **Task B1.** In `cancel_recursive` (`rat.c`, after the FLINT-over-Q decline at `:810`), detect
  algebraic-number coefficients (Complex + radical / `Root[]`) and route the GCD through the
  FLINT `qqbar` / number-field layer (`flint_algebraic_field_normalize`, cyclotomic `Together`;
  see memories *FLINT extension engine*, *RootReduce qqbar*, *Cyclotomic extension support*)
  rather than the generic `poly_gcd_internal`.
- **Task B2.** Give the Cherry coefficient solve and `rt_verify_antideriv` a number-field zero-test
  that does not depend on generic `Simplify` — the exact `flint_algebraic_field_normalize` over
  Q(i√d) (this is the `CHERRY_BLOCKERS.md` §A1 "FLINT number-field linear solve" already named as
  the resolution for the general complex case).
- **Task B3.** Safety net: a structural termination guard in `exact_poly_div` — exact division
  strictly reduces `deg(R)` each step, so a non-decreasing step means the leading coefficient did
  not cancel in this ring ⇒ decline (`return NULL`, leave uncancelled) instead of spinning. This
  is a *structural* bound, not an arbitrary cap; it converts any *future* mis-routed input from a
  hang into a clean decline. (Prototyped here; did not fix the §2.2 case because that hang is
  within a single expansion step, but it is a cheap, correctness-preserving backstop worth
  landing with B1.)
- **Effort:** large. **Risk:** medium — correct-by-construction if the field arithmetic is exact;
  guard against timeouts; prove non-regression against the whole `rat`/`simp`/Risch battery.

**Recommended order:** Track A first (un-gates all §1 examples at medium cost by reusing the
number-field fallback that already exists), then Track B as the durable hardening (removes the
generic-GCD hazard system-wide and lets the coefficient solve itself run in the number field).

---

## 4. Verification pins (when the work lands)

- **Un-gate:** remove `expr_has_radical` gate in `cherry_ei.c` `gen_alpha_candidates`; each §1
  example must return a closed `E^{−α} ei(x+α)` conjugate pair with **exact diff-back = 0** and a
  wall-clock ceiling (no hang). Flip the `assert_declines` back to `assert_ei` in
  `tests/test_cherry_ei.c::test_complex_ei`.
- **Primitive pins (Track A/B):** `Simplify[E^x/(x + (1−I Sqrt[3])/2)]`,
  `Simplify[1/(x−a) − 1/(x−ā)]`, and `RootReduce`/`FullSimplify` siblings must return fast.
- **Regression:** full `simplify_tests`, `radical_simplify_tests`, `simp_algebraic_cuberoot_tests`,
  `rat_tests`, `cherry_ei_tests`, `integrate_risch_transcendental_tests` unchanged; the Q(i) and
  real-algebraic ei families stay green; `valgrind` clean on the new paths.

---

## 5. One-line summary

The Q(i√d) declines are one shared infrastructure gap: a `Simplify` pass substitutes the complex
radical into a polynomial generator and the generic multivariate `exact_poly_div` then blows up —
even though `Together`/`Cancel` over Q(i√d) already work. Guard that substitution (Track A) to
un-gate every §1 example via the existing number-field fallback; then route algebraic-coefficient
cancellation through exact FLINT number-field arithmetic (Track B) to remove the hazard for good.

# Knowles erf/li Integration — Future Extensions

**Status:** the `knowles_erf.c` engine (K2) closes transcendental **Liouvillian**
integrands in terms of elementary functions plus `Erf`/`Erfi`, over the K0
`RT_PRIM` tower, via the perfect-square gate + undetermined-coefficient elementary
part + `SolveAlways` + diff-back gate. The **radical / quasiquadratic** case
(Part I §6, Ex 8.1) landed 2026-07-21 (`collapse_exp_of_log` + `s_k = Sqrt[g_k]`
solve). See [`KNOWLES_DESIGN.md`](KNOWLES_DESIGN.md) for the full spec and
[`docs/spec/changelog/2026-07-20.md`](docs/spec/changelog/2026-07-20.md) for the
landed increments.

This document records the **remaining** extensions, each with a concrete integral
it would unlock. Every "declines today" claim below was verified against the
current build; every target antiderivative diff-backs to its integrand (checked
symbolically or numerically). The engine is **sound by construction** (diff-back
gates every emission), so none of these are correctness bugs — they are
*completeness* gaps: integrals with an erf/li-elementary antiderivative that the
current bounded search does not find, and so declines.

> Convention reminder: Mathilda's classical `Erf` has `d/dx Erf[u] = (2/√π)e^{−u²}`;
> Knowles' `erf(u) = ∫e^{−u²}du = (√π/2) Erf[u]`. Targets below are in Mathilda's
> `Erf`.

---

## 1. `x`-rational (non-constant) elementary-part coefficients

**What's missing.** The elementary part `v` in `∫g = v + Σ kᵢ erf(uᵢ)` is currently
a bounded-degree polynomial in the tower monomials with **constant** coefficients
(`kerf$c` symbols in `knowles_erf.c`). Knowles/SSC allow `v ∈ F` — any element of
the differential field, i.e. coefficients **rational in `x`** (and in the tower
variables). Any integrand whose antiderivative has an `x`-dependent elementary part
declines.

**Concrete unlock** (declines today):

```
∫ [ E^(−Erf[x]²) − (4x/√π) Erf[x] E^(−x²−Erf[x]²) ] dx  =  x · E^(−Erf[x]²)
```

The elementary part `x·E^(−Erf[x]²)` carries the coefficient `x`, which the
constant-coefficient ansatz cannot represent, so the linear solve has no solution
and the engine declines. (Verified: the integrand is `d/dx[x·E^(−Erf[x]²)]`, and
`Integrate[…]` returns unevaluated.)

**Implementation sketch.** Replace each constant `kerf$c_b` with an *undetermined
rational function* `c_b(x) = p_b(x)/q_b(x)` of bounded degree — i.e. run the
elementary-part step through the field-`F` limited-integration oracle
(`rt_limited_integrate` / `Integrate\`BronsteinRational`) that the recursive Risch
already uses for the base field, rather than a pure `SolveAlways` over constants.
The bound on `deg v` comes from the same rank analysis as the erf-term count.

---

## 2. Additive decomposition across special-function families (mixed erf + li)

**What's missing.** Dispatch (`RT_SPECIAL_FORMS`) routes an integrand to **one**
special-form engine by its top monomial kind (`RT_SF_TOP_EXP` → erf/ei,
`RT_SF_TOP_LOG` → li/dilog). An integrand whose antiderivative needs an `Erf` term
**and** a `LogIntegral` term simultaneously is never split, so it declines even
though each summand integrates on its own. This is the outward face of the deferred
**deep-tower recursion** (Part II Lem 2.1/2.2, Prop 3.4): fire the argument
generators on the inner peeled coefficients, not just the outermost integrand
(`CHERRY_BLOCKERS.md` B2).

**Concrete unlock** (declines today; each summand works alone):

```
∫ [ 2 E^(−x²) + 1/(Log[x] · Log[LogIntegral[x]]) ] dx  =  √π Erf[x] + LogIntegral[LogIntegral[x]]
```

Verified: `Integrate[2 E^(−x²), x] = √π Erf[x]` and
`Integrate[1/(Log[x] Log[LogIntegral[x]]), x] = LogIntegral[LogIntegral[x]]` both
evaluate, but their **sum** returns unevaluated.

**Implementation sketch.** At the `rt_field_integrate` / `extended_liouville_solve`
decline, peel the top monomial by partial fractions (`Risch\`CanonicalRepresentation`)
and recurse `extended_liouville_solve(remainder, x, top−1)` on each coefficient
(the K0 "deep-tower recursion hook", already stubbed in `tasks/todo.md` C.1). Sum
the erf-part and li-part results. Soundness stays free (final diff-back gate).

---

## 3. Completing-the-square erf arguments (affine / linear-shift arguments)

**What's missing.** The perfect-square gate accepts a top exponential `e^w` only
when `−w` is *literally* a perfect square (root `u`). Knowles Part I (Lemma 6.2,
the `Π(s+β)=R²` gate) allows an additive **completing-the-square** step: `−w = u² + c`
with `c` a constant folded into the elementary part, so the erf argument is an
**affine shift** `u = (kernel + const)`. The current `rt_expand_exp_sums` split
even *separates* the base-field `−x²` term from the kernel terms, hiding the cross
term needed to complete the square.

**Concrete unlock** (declines today):

```
∫ E^(−x² − 2 Erf[x] − Erf[x]²) dx  =  (e·π/4) · Erf[Erf[x] + 1]
```

Here `−x² − 2Erf[x] − Erf[x]² = −x² − (Erf[x]+1)² + 1`, so the erf argument is the
shifted `Erf[x] + 1` and the constant `+1` becomes the `e` factor. (Verified
numerically: `D[(eπ/4)Erf[Erf[x]+1]] − integrand ≈ 2.8×10⁻¹⁷`; `Integrate[…]`
returns unevaluated.)

**Implementation sketch.** Before the perfect-square gate, run the completing-square
β-finder (the one already in `cherry_ei.c` for the Gaussian `∫P(x)e^{−cx²}` case)
over the *combined* exponent (do not pre-split the base-field term away from the
kernel terms for this test): find `β` constant with `−w − β = u²`, emit `erf(u)`,
and absorb `e^{β}` into `v`.

---

## 4. Certified non-existence (the "NO" verdict)

**What's missing.** Declines are **sound but not certified**. When the engine
declines, it means "the bounded candidate search (≤ 2 erf terms, degree-bounded `v`,
the perfect-square/Σ-decomposition gate) found no antiderivative" — not a *proof*
that none exists. Knowles Part I Cases 1–2.3.2 + Prop 3.4, with the structural
bounds implemented faithfully, upgrade this to a genuine decision procedure that
returns a certified **NO**.

**Concrete example** (correctly declines today, and a real NO):

```
∫ x · E^(−x² − Erf[x]²) dx        (Part II Ex 4.2 — provably no erf-elementary antiderivative)
```

The engine declines, which is the right answer — but only Part I's rank analysis
turns that decline into the theorem "no `v + Σ kᵢ erf(uᵢ)` exists." Until then, a
decline could in principle also be a *completeness* miss (items 1–3 above), so the
two are indistinguishable to a caller. Certifying NO requires proving the candidate
set is exhaustive: the erf-term count bound, Σ-decomposition termination, and the
perfect-square gate finiteness (Lemma 6.2).

---

## 5. Algebraic constants in the coefficient solve (`C = C̄`)

**What's missing.** Algebraic constants in the erf **argument** already work — e.g.
`∫E^(−2x²) = (√π/(2√2)) Erf[√2 x]` closes with `√2` in the argument. The gap is when
the constant **solve** itself lands on an algebraic number: a coupled multi-`erf`
system whose coefficients `kᵢ` are roots of an irreducible quadratic (or higher)
over `Q`. `SolveAlways` runs over `Q`; the algebraic case needs the FLINT
number-field solve shared with **Cherry A1** (`RootReduce`/`qqbar`, off the critical
path for the rational-constant pins).

**Illustrative shape.** Any target of the form `k₁ Erf[u₁] + k₂ Erf[u₂]` where
`{k₁, k₂}` are the two roots of `k² − k − 1 = 0` (so `k = (1±√5)/2`) — the golden-
ratio coefficients only appear after solving the coupling, and a rational
`SolveAlways` returns no solution. (This item shares its infrastructure dependency
with Cherry A1; see `CHERRY_BLOCKERS.md`.)

---

## Summary

| # | Extension | Concrete unlock (target antiderivative) | Status today |
|---|-----------|------------------------------------------|--------------|
| 1 | `x`-rational elementary part | `∫[E^(−Erf²x) − (4x/√π)Erf x·E^(−x²−Erf²x)] = x·E^(−Erf²x)` | declines |
| 2 | Mixed erf + li (deep-tower recursion) | `∫[2E^(−x²) + 1/(Log x·Log li x)] = √π Erf[x] + li(li x)` | declines |
| 3 | Completing-square (affine erf arg) | `∫E^(−x²−2Erf x−Erf²x) = (eπ/4)Erf[Erf x + 1]` | declines |
| 4 | Certified non-existence | `∫x·E^(−x²−Erf²x)` — certify the sound decline (Part II Ex 4.2) | sound, not certified |
| 5 | Algebraic constant solve (`C=C̄`) | coupled `k₁ Erf[u₁]+k₂ Erf[u₂]`, `kᵢ` algebraic (Cherry A1 infra) | declines |

Items 1–3 are self-contained completeness increments (each reuses existing
machinery — the field-`F` oracle, the deep-tower hook, the `cherry_ei.c`
completing-square β-finder — behind the same `rt_tower_solve` / diff-back seam).
Item 4 is the decision-procedure hardening (Part I Cases 1–2.3.2). Item 5 is
infra-gated on the shared Cherry A1 FLINT number-field solve.

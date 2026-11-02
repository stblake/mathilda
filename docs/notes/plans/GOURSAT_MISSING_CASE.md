# Why Mathilda's Goursat integrator cannot (yet) close these cube-root integrals

A mathematical summary of the obstruction, for

```
(a)  ∫ ((x-1)² (x+1))^(1/3) / x²  dx
(b)  ∫ (1 - k x) / ((1 + (k-2) x) · (x(1-x)(1-k x))^(2/3))  dx      (k a parameter)
```

Both are **elementary** (this is not in dispute), yet the Goursat cube-root
integrator returns them unevaluated. The reasons are different for (a) and (b),
but they share one root cause: the algorithm's *elementarity certificate* is too
narrow for integrands whose rational cofactor has a pole **off the branch locus**.

---

## 1. The geometry: an elliptic curve

Write the integrand as `F(t) · R(t)^(-p)`, `F` rational, `R` a cubic, `p ∈ {1/3, 2/3}`.
The natural home of the differential `ω = F R^(-p) dt` is the cyclic trigonal curve

```
        C :  y³ = R(t),        y = R(t)^(1/3).
```

For a **squarefree** cubic `R`, `C` is ramified (index 3) over each of the three
roots of `R` and unramified over `∞` (since `deg R = 3 ≡ 0 mod 3`). Riemann–Hurwitz:

```
   2g - 2 = 3·(2·0 - 2) + Σ (e_P - 1) = -6 + 3·(3-1) = 0   ⟹   g = 1.
```

**`C` is an elliptic curve.** This single fact drives everything below.

Case (b) has this generic genus-1 geometry. Case (a) is special: `R = (x-1)²(x+1)`
is **not** squarefree (a double root at `x=1`), which collapses the genus to `0` and
also breaks a structural assumption of the algorithm — treated separately in §5.

---

## 2. What the algorithm actually certifies (Goursat eigendescent)

The implemented cube-root method (`goursat_cubic`, mirroring `CubicTest`/`CubicTest23`
of `GoursatAppendix.wl`) is:

1. Build the order-3 Möbius map `S` cyclically permuting the three roots of `R`.
   Its fixed points `α, β` give the coordinate `z = (t-α)/(t-β)`, in which `S` acts
   as `z ↦ ω z` (`ω = e^{2πi/3}`) and `R(t(z))·(1-z)³ = c (z³ - K)`.
2. Transport the cofactor: `H̃(z) = (α-β) F(t(z)) / c^{2/3}` (the `p=2/3` case), and
   split it into the three `ω`-eigencomponents under `z ↦ ω z`:
   `H̃ = H̃₀ + H̃₁ + H̃₂`, with `H̃_k(z) = z^k φ̃_k(z³)`.
3. **Elementarity criterion.** The integral is declared elementary iff the
   *obstruction eigencomponent* vanishes:
   ```
        p = 1/3 :  H₁ ≡ 0            p = 2/3 :  H̃₀ ≡ 0.
   ```
   When it holds, the substitution `x = z³` reduces each surviving component to a
   **rational** integral (Corollary 4.8(i) of `..._30042026B.tex`):
   ```
     ∫ F R^(-2/3) dt = −∫ φ̃₁(K s³/(s³-1)) · s ds/(s³-1)  +  ∫ φ̃₂(u³+K) du.
   ```

The content of the criterion, in classical language, is **Liouville's theorem for
abelian integrals**:

> `∫ ω` is elementary  ⟺  `ω = dg + Σ_i c_i · du_i/u_i`
> for algebraic functions `g, u_i` on `C`, with constant `c_i`.

Decompose `ω` into its three Liouville pieces:

* **first kind** — holomorphic (no poles). On a genus-1 curve the space of these is
  1-dimensional and **spanned by a non-elementary period** (the elliptic integral).
* **second kind** — poles but *no residues*. Always elementary (it is `dg`).
* **third kind** — simple poles with nonzero residues.

The eigendescent criterion `H̃₀ ≡ 0` is exactly the statement **"the first-kind part
of `ω` vanishes."** That is *necessary and sufficient* for elementarity **only when
`ω` has no third-kind part** — i.e. when `F` is a polynomial, or more generally when
`F`'s poles all sit at branch points of `C` (roots of `R`). The reference proves this
sufficiency precisely under the hypothesis that `φ̃₁` (equivalently the relevant
projection) is a **polynomial** (Theorem 4.3(iii)).

---

## 3. Why case (b) escapes the criterion

In (b),

```
   F(t) = (1 - k t)/(1 + (k-2) t),      R(t) = t(1-t)(1-k t).
```

`F` has a **simple pole at** `t₀ = -1/(k-2)`, and

```
   R(t₀) ≠ 0   ⟹   t₀ is NOT a branch point of C.
```

So `ω = F R^(-2/3) dt` has genuine **third-kind** behaviour: over `t₀` sit three
points `P₀, P₁, P₂` of `C` (the three cube-root sheets), and `ω` has a simple pole
with nonzero residue at each. Two consequences:

1. **The criterion no longer decides.** Measuring the obstruction directly gives
   `H̃₀ ≠ 0` (numerically `H̃₀(z=0.37) = 0.193 − 0.805 i`, verified to 30 digits — not
   a zero-test artefact). By the reference's own residue analysis, `H̃₀ ≠ 0` proves
   non-elementarity **only** when the obstruction is second-kind (`φ̃₁` polynomial).
   Here `F` is a *rational* cofactor, so `H̃₀ ≠ 0` carries a third-kind contribution
   and **certifies nothing**. The algorithm, seeing `H̃₀ ≠ 0`, declines — a false
   negative.

2. **Elementarity is now an Abel–Jacobi torsion condition.** A third-kind
   differential with residues `c_i` at points `P_i` integrates to elementary form iff
   the divisor `D = Σ_i c_i P_i` is (after clearing denominators of the `c_i`) the
   divisor of an algebraic function `u` on `C` — i.e. iff `D` is **principal**, which
   for the genus-1 curve `C` means its image under the Abel–Jacobi map is a **torsion
   point of the elliptic curve's Jacobian**:
   ```
       Σ_i c_i · [P_i]  =  0   in   Pic⁰(C)   (a torsion relation).
   ```
   Then `ω = dg + c · du/u` and `∫ω = g + c Log u`.

The integrand (b) is **elementary precisely because its pole was placed to make this
divisor torsion** — the factor `1 + (k-2)x` in the denominator is tuned (uniformly in
`k`) so that `P₀+P₁+P₂` over `t₀` lands on a torsion point of `C`. That is the whole
art of Goursat's examples. But **nothing in the current algorithm computes residues,
forms the divisor, or tests torsion in `Jac(C)`.** The `H̃₀ ≡ 0` gate sees only the
holomorphic part; it is blind to the torsion condition that actually governs the
third-kind case. This is exactly the step the preprint (`..._30042026B.tex`,
lines 838–840) **explicitly defers to a Risch–Trager–Bronstein analysis**.

### Why a naive log-ansatz patch does not work

One might try to skip the torsion machinery and directly fit the log part

```
   ∫ F R^(-2/3) dt  =?  C · Σ_j ω^{j} Log[ R^{1/3} - ω^j (κ t + m) ],
   κ = (lead R)^{1/3},   m = R(t₀)^{1/3} - κ t₀   (centre on the pole t₀).
```

This *does* pass the local test (its differential is proportional to `y = R^{1/3}`,
so the two non-matched eigen-slots vanish — confirmed numerically). But the required
scale `C = F·N/(R·p_match)` comes out **not constant in `t`** (`Craw(1/7)=0.686` vs
`Craw(2/9)=0.701`; `|D[G]-f| = 7.24`). The reason is structural: the log arguments
of the *true* answer are the back-substitution variables

```
   s = c^{1/3}(t-α) / (R^{1/3}(α-β)),        u = R^{1/3}(α-β) / (c^{1/3}(t-β)),
```

i.e. `R^{1/3}` sits in the **denominator** — the arguments are genuine algebraic
*curve functions*, not `radical − (linear in t)`. Moreover `F` splits as

```
   F = (1-kt)/(1+(k-2)t) =  −k/(k-2)   +   (2k-2)/((k-2)((k-2)t+1)),
                            └ constant ┘   └──── pole part at t₀ ────┘
```

whose two pieces land in **different** eigencomponents: the constant part integrates
to an **algebraic** term (a rational function of `R^{1/3}`, the `x=u³+K` branch), the
pole part to a **Log** (the `x=Ks³/(s³-1)` branch). The antiderivative is genuinely
`algebraic + logarithmic`, which no single linear-argument log family can span. The
non-constant `C` is the concrete symptom.

---

## 4. The performance layer (a solved sub-problem)

Independently of the certificate gap, even the cases the algorithm *should* close
were failing for a mundane reason: the eigendescent lives over an **algebraic tower**
`Q(k)(k^{1/3}, R(t₀)^{1/3}, ω)`, and Mathilda's native `Cancel`/`Together` over that
tower **does not terminate** (still running at 300 s per pass in these examples).

This part is fixable and was verified: routing the tower zero-tests and
normalisations through the FLINT genuine-algebraic field engine
(`flint_algebraic_field_normalize` / `flint_algebraic_field_canonical`) reduces the
slot decisions from "non-terminating" to `< 1 s` and returns numerically correct
verdicts. So the tower blow-up is **not** the fundamental obstruction — the missing
elementarity certificate (§3) is.

---

## 5. Case (a): a different, structural failure

```
   ((x-1)² (x+1))^(1/3) / x²   →   F = (x-1)²(x+1)/x²,   R = (x-1)²(x+1),   p = 2/3.
```

Two independent problems:

* **Non-squarefree radicand.** `R` has a *double root* at `x=1`, so it has only two
  distinct roots. The order-3 cyclic Möbius construction assumes **three distinct
  branch points**; the root-count gate (`dR=3, nr=2`) rejects the integrand before any
  descent. (Geometrically the curve is now genus 0, hence everything is elementary —
  but the algorithm never gets that far.)
* **A double pole off the branch locus.** `F = x - 1 - 1/x + 1/x²` has a *second-order*
  pole at `x=0`, which is not a root of `R`. That is a **second-kind** part on top of
  the third-kind pole — outside the simple-pole scope of both the eigendescent and any
  third-kind log construction.

So (a) needs (i) squarefree/double-root handling in the Möbius step and (ii) a
second-kind reduction for the higher-order off-branch pole — neither present.

---

## 6. Summary

| | curve | `F` pole | why it fails today |
|---|---|---|---|
| (a) | genus 0 (double root) | double pole at `x=0` (off branch) | Möbius step rejects non-squarefree `R`; higher-order off-branch pole is second-kind, unhandled |
| (b) | genus 1 (elliptic) | simple pole at `t₀=-1/(k-2)` (off branch) | criterion `H̃₀≡0` tests only the holomorphic part; the third-kind pole's elementarity is an **Abel–Jacobi torsion** condition the algorithm never computes |

The essential gap is **§3**: for a rational cofactor with an off-branch pole, the
implemented criterion `H̃₀ ≡ 0` is the wrong test. Elementarity is then decided by
whether the residue divisor of the off-branch pole is torsion in `Jac(C)`, and the
answer is `algebraic part + Σ c_i Log(u_i)` with the `u_i` the algebraic curve
functions realising that torsion relation. Constructing them is precisely the
cube-root **third-kind Risch–Trager–Bronstein** step that the reference defers; the
two-substitution *splitting* alternative (Welz) is asserted in the literature but
never given a criterion or a formula. Until one of those is implemented, (b) is out
of reach on principle (not merely on performance), and (a) additionally needs
non-squarefree + second-kind handling.

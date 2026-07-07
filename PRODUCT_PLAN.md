# `Product` — Survey of Symbolic Product Algorithms

> A research survey of the known algorithms for computing **indefinite** and
> **definite** symbolic products — Mathematica's `Product[]`, Maple's
> `product()`, Maxima's `product()`, and SymPy's `Product` — together with the
> decision-procedure literature behind them, their reference papers,
> implementation complexity for a small C/GMP CAS, and effectiveness on
> well-known problems.
>
> This document is the *symbolic* counterpart to the already-shipped numerical
> `NProduct` (see `NPRODUCT_PLAN.md`). `NProduct` evaluates a product to a number
> via `Exp[NSum[Log[f], …]]`; `Product` must return an **exact closed form**
> (`n!`, `Pochhammer[a, n]`, `Gamma`-ratios, `QPochhammer`, …) and is therefore a
> different problem requiring different machinery.
>
> Scope note: this is a survey + roadmap, not a build sheet. A per-stage build
> sheet (in the style of `SUM_DEVEL_PLAN.md`) is the natural follow-on once the
> staging in §8 is chosen.

---

## 1. The problem and its duality with summation

A **product** is the multiplicative analogue of a sum. Everything in the
summation literature has a product image under the `log`/`exp` correspondence,
so the cleanest way to understand product algorithms is as the multiplicative
dual of the `Sum` machinery Mathilda already has in `src/sum/`.

| Summation                              | Product (multiplicative dual)                  |
|----------------------------------------|------------------------------------------------|
| indefinite sum `S` with `S(k+1)−S(k)=f`| indefinite product `P` with **`P(k+1)/P(k)=f`** |
| telescoping `Σ f = S(n+1) − S(a)`      | telescoping **`∏ f = P(n+1)/P(a)`**            |
| antidifference                         | **anti-quotient** (multiplicative antidifference)|
| closed forms in polynomials, harmonic # | closed forms in **factorials / Pochhammer / Gamma** |
| Gosper's algorithm (have it: `sum_gosper.c`) | Gosper-Petkovšek normal form, read off multiplicatively |
| `Sum[…]` cascade (`sum.c`)             | proposed `Product[…]` cascade (`product.c`)    |

Two structural facts drive the whole design:

1. **A term is *hypergeometric* iff its term-ratio `f(k+1)/f(k)` is rational in
   `k`.** This rationality is exactly what makes closed forms decidable. The
   natural domain of `Product` is hypergeometric terms, whose products live in
   Gamma/Pochhammer space.

2. **Naïve `∏ f = exp(Σ log f)` does *not* work symbolically** — `log f(k)` is
   generally *not* hypergeometric even when `f` is, so the sum has no closed
   form. (This is why it works *numerically* for `NProduct` but not here.)
   Products must be handled in their own multiplicative terms.

The two problems are:

- **Indefinite product** `∏ f(k)` with a *symbolic* upper index `n` → find an
  anti-quotient `P(n)` (e.g. `Product[k+a, {k,1,n}] = Pochhammer[a+1, n]`).
- **Definite product** `∏_{k=a}^{b} f(k)` → evaluate the anti-quotient at the
  endpoints, `P(b+1)/P(a)`; the **infinite** case is the `b→∞` limit of that
  closed form, when it converges (§6).

---

## 2. What the reference systems actually do

A short orientation before the algorithm-by-algorithm survey. The systems agree
on the *answers* (the Gamma/Pochhammer closed forms) and differ mostly in how
hard they try.

### Mathematica `Product[]`

The strongest of the four. Returns minimal factorial/Gamma closed forms for
hypergeometric terms, recognises `QPochhammer` for `q`-factors, and evaluates
convergent infinite products through Gamma reflection, the Barnes `G` function,
and special constants (Glaisher–Kinkelin). Worked outputs (verified against the
Wolfram reference page):

| Input | Output |
|---|---|
| `Product[k, {k,1,n}]` | `n!` |
| `Product[k^2, {k,1,n}]` | `(n!)^2` |
| `Product[x+i, {i,0,n-1}]` | `Pochhammer[x, n]` |
| `Product[x-i, {i,0,n-1}]` | `FactorialPower[x, n]` |
| `Product[1 - 1/k^2, {k,2,n}]` | `(n+1)/(2 n)` (telescopes) |
| `Product[1 - 1/k^2, {k,2,∞}]` | `1/2` |
| `Product[1 + 1/k, {k,1,n}]` | `n+1` |
| `Product[1 + 1/i, {i,1,∞}]` | divergent (doc's own example) |
| `Product[k/(k+1), {k,1,n}]` | `1/(n+1)` |
| `Product[1 - a q^i, {i,0,n-1}]` | `QPochhammer[a, q, n]` |
| `Product[Gamma[i], {i, n-1}]` | `BarnesG[n+1]` |

Limitations: divergent products, multiplicand limit ≠ 1, non-integer-spaced
irrational/transcendental roots, and non-hypergeometric terms (`∏(k!+1)`) are
returned unevaluated; Barnes `G`/Glaisher results require those special
functions to exist in the system.

### Maple `product()` / `Product()`

`product` is the active form, `Product` the inert (pretty-printed Π) form. For a
finite integer range it unrolls; for a symbolic bound it produces **`GAMMA`-ratio
/ Pochhammer** closed forms by the same root-matching construction. Verified
documented outputs: `product(k^2, k=1..n)` → `GAMMA(n+1)^2` (i.e. `(n!)^2`),
`Product(n+k, k=0..m)` → `GAMMA(n+m+1)/GAMMA(n)`, `product(k, k=x..5x)` →
`GAMMA(5x+1)/GAMMA(x)`, `product(k*x, k=1..5)` → `120 x^5`. Irreducible
denominators are expressed via `RootOf` (product over polynomial roots). The
product internals are not separately documented; the companion *summation*
machinery is `SumTools[IndefiniteSum]` (Abramov's algorithm) and
`SumTools[Hypergeometric]` (`Gosper`, `Zeilberger`, `RationalCanonicalForm`, …).

### Maxima `product()`

The weak baseline. From the manual: if the bounds differ by an integer the
product is unrolled, "otherwise the range of the index is indefinite" and a noun
(Π) form is returned. The option flag `simpproduct:true` enables limited
simplification (enough for `product(k,k,1,n) → n!`) but there is **no general
rational→Gamma engine**. Note `nusum`/`unsum` are *summation* (Gosper) tools, not
products — Maxima has no Gosper-style product engine.

### SymPy `Product` / `concrete/products.py`

The cleanest **port target**, because `_eval_product` is an explicit cascade we
can mirror almost directly (dispatch order verified from source):

1. term free of `k` → `term**(n-a+1)`; `term ≡ 1` → `1`;
2. trivial `a == n` → substitute;
3. `KroneckerDelta` → `deltaproduct`;
4. small definite integer range → direct unroll;
5. **polynomial term → roots → `RisingFactorial(a-r, n-a+1)^m` per root `r` of
   multiplicity `m`, times `LC^(n-a+1)`** (the core rational→Pochhammer path);
6. `Add` → factor and recurse;
7. `Mul` → split factors with/without `k`, recurse (rational = numerator-product
   / denominator-product);
8. `Power base^exp`: base free of `k` → `base ** Sum(exp)` (geometric, handled by
   the *sum* engine); exp free of `k` → product of base, powered;
9. nested `Product` → inner first;
10. fallback direct computation.

It follows **Karr's convention** for empty/reversed ranges. SymPy's `Product`
produces the Wallis product in closed form via path 5/7, e.g.
`Product[(2i/(2i-1))(2i/(2i+1)), {i,1,n}]` →
`2^{-2n} 4^n (n!)^2 / ((1/2)_n (3/2)_n)`.

### mpmath `nprod` (for completeness)

Purely numerical: extrapolated partial products, or `exp(nsum(log f))` with
`nsum=True`. This is exactly Mathilda's existing `NProduct`. It is the
**validation reference** for the symbolic engine, not new work.

---

## 3. Indefinite products — the practical algorithms

These are the algorithms a `Product` builtin should actually implement first.
They cover the overwhelming majority of textbook inputs and rest only on
machinery Mathilda already has (`src/poly/` factorization, `Pochhammer`,
`Gamma`).

### 3.1 Constant / trivial / finite-unroll dispatch

**Overview.** The cheap front of the cascade, mirroring SymPy paths 1–4 and
`sum.c`'s surface handling: term free of the index → power; `a == b` →
substitute; a finite integer range below a size threshold → multiply out
directly (reusing the evaluator). Empty/reversed ranges follow Karr's
convention (`∏_{k=a}^{a-1} = 1`).

**Reference.** M. Karr, "Summation in finite terms," *J. ACM* 28(2):305–350
(1981) — the source of the empty-range convention SymPy adopts.

**Complexity.** Trivial. ~100 LoC; pure expression plumbing, reuses
`sum_stage_args`-style iterator parsing.

**Effectiveness.** Handles every concrete finite product and the degenerate
cases; produces nothing symbolic on its own but gates the harder stages.

### 3.2 Multiplicative telescoping (rational anti-quotient)

**Overview.** The product analogue of summation telescoping. Detect
`f(k) = g(k+1)/g(k)` for a rational `g`; then `∏_{k=a}^{b} f = g(b+1)/g(a)`. This
catches the "obvious" collapses without any factorial machinery:
`∏ (k+1)/k = (n+1)`, `∏ k/(k+1) = 1/(n+1)`, `∏ (1-1/k^2) = ∏ (k-1)(k+1)/k^2`. It
is the multiplicative image of the `c(k+1)/c(k)` shell in the Gosper-Petkovšek
form (§5.1), so in a mature implementation it falls out of that normal form
rather than being a separate detector — but a standalone telescoping pass is a
high-value, low-risk first increment.

**Reference.** The product specialisation of Gosper (1978; §5.1); textbook
treatment in Graham–Knuth–Patashnik, *Concrete Mathematics*, 2nd ed. (1994),
ch. 2 (finite calculus / products).

**Complexity.** Low (~200–300 LoC) given polynomial GCD and shift. Detect
shift-equivalence of numerator and denominator factors (the dispersion idea of
§5.1 in miniature).

**Effectiveness.** Decides every product whose anti-quotient is itself rational
— a large, common class. Cannot express products that genuinely need factorials
(`∏ k = n!`): those need §3.3.

### 3.3 Rational-function products → Gamma / Pochhammer / Factorial

**The single most valuable stage.** Shared by Mathematica, Maple, and SymPy; it
is the workhorse.

**Overview.** Write the term as
`f(k) = c · ∏_i (k - α_i) / ∏_j (k - β_j)` by factoring numerator and
denominator. Each linear factor telescopes into a shifted factorial via the
atomic identity

```
∏_{k=1}^{n} (k + a) = Pochhammer[a+1, n] = Gamma[n+a+1] / Gamma[a+1],
```

so the whole finite product becomes `c^n` times a ratio of
Pochhammer/Gamma/Factorial values. The single-shift template that the
implementation keys on:

```
∏_{k=a}^{b} (k+α)/(k+β) = Gamma[b+1+α] Gamma[a+β] / ( Gamma[a+α] Gamma[b+1+β] ).
```

The algorithmic refinement that produces *minimal* output: **match numerator
and denominator roots that differ by an integer.** A pair `α, β` with
`α − β ∈ ℤ` telescopes into a rational (polynomial) part rather than residual
Gammas, cutting the number of Gamma factors. Wolfram states products of rational
functions are returned "with a minimal number of factorial functions"; this
root-difference matching is how. Half-integer shifts stay inside the framework:
`∏_{k=1}^{n} (2k-1)/(2k) = (1/2)_n / (1)_n = Binomial[2n,n]/4^n`.

**References.**
- Foundational *summation* analogue (often mis-cited for products): R. Moenck,
  "On computing closed forms for summations," *Proc. 1977 MACSYMA Users' Conf.*,
  pp. 225–236. **⚠ This is a summation paper** (Hermite-style rational
  summation, polygamma remainder); it does **not** contain the product→Gamma
  algorithm. Cite it for the rational-*summation* lineage only.
- The actual product→Gamma normal form: W. Koepf, "Algorithms for the indefinite
  and definite summation," arXiv:math/9412227 (1995) — the Γ-ratio normalisation
  (build `a_k/a_{k-1}`, convert factorial/binomial/Pochhammer to Γ, rewrite
  `Γ(a+k) = (a)_k Γ(a)` for positive-integer argument differences, cancel).
- Pirastu & Strehl, "Rational summation and Gosper–Petkovšek representation"
  (the integer-root-difference matching).
- Textbook: Graham–Knuth–Patashnik, *Concrete Mathematics*, 2nd ed. (1994),
  ch. 2; Andrews–Askey–Roy, *Special Functions*, CUP (1999), ch. 1 (the
  `(a)_n = Γ(a+n)/Γ(a)` Pochhammer machinery).

**Complexity.** Moderate (~400–700 LoC) on top of existing infrastructure:
factor numerator & denominator over ℚ (have it: `src/poly/`), collect roots as
GMP rationals, bucket by fractional part and pairwise-match integer-spaced roots
(`O(p·q)` naïve, `O(n log n)` with sorting), telescope matched families to a
rational part, emit `Pochhammer`/`Gamma`/`Factorial` for the residue. Risk areas:
factor multiplicities, the `c^n` leading-coefficient term, and lower-limit
normalisation. Irreducible quadratic/irrational roots need the algebraic-number
layer (or stay as `RootOf`-style Gamma products); restricting to ℚ-linear roots
is the pragmatic first cut.

**Effectiveness.** Complete for rational `f` factoring into integer-spaced
ℚ-linear factors — i.e. essentially all textbook `Product` inputs:
`Product[k,{k,1,n}]=n!`, `Product[k^2,…]=(n!)^2`, `Product[(k+a),…]=Pochhammer[a+1,n]`,
`Product[k/(k+1),…]=1/(n+1)`. Degrades to an unreduced Gamma product when roots
are irrational and don't pair into integer-difference families.

### 3.4 Geometric and polynomial-exponential factors

**Overview.** A factor `r^k` with `r` free of `k` contributes
`r^{Σ k} = r^{(geometric/arithmetic sum of the exponent)}`, i.e. it routes the
*exponent* through the existing `Sum` engine (SymPy path 8). `Product[2^k,{k,1,n}]
= 2^{n(n+1)/2}`. A `Mul` of a geometric factor and a rational factor splits and
each part is handled independently (SymPy path 7).

**Reference.** SymPy `concrete/products.py` path 8; the exponent sum is exactly
Mathilda's `Sum` (`sum_polynomial.c` / `sum_geometric.c`).

**Complexity.** Low (~150 LoC). The work is recognising the factor split and
delegating the exponent to `Sum`; reuse is the whole point.

**Effectiveness.** Closes all `p(k) r^k`-shaped products whose rational part is
already handled by §3.3. Direct reuse of shipped `Sum` stages.

---

## 4. `q`-products

**Overview.** The `q`-analogue of §3.3. A term is `q`-hypergeometric when
`t(k+1)/t(k)` is rational in `q^k`. Such a factor splits into linear pieces
`(1 - a_i q^k)`, and

```
∏_{k=0}^{n-1} (1 - a q^k) = QPochhammer[a, q, n],
```

so definite `q`-products collapse to **ratios of `q`-Pochhammer symbols** — the
exact mirror of §3.3 with `Gamma → QGamma`, `Pochhammer → QPochhammer`. The
machinery for *deciding* the indefinite case is the **`q`-Gosper algorithm**
(the multiplicative side is the `q`-telescoping dual); `q`-numbers, `q`-factorial,
`q`-Gamma, and Gaussian binomials are thin wrappers over `QPochhammer`:
`[n]_q! = QPochhammer[q,q,n]/(1-q)^n`.

**References.**
- H. Böing & W. Koepf, "Algorithms for q-hypergeometric summation in computer
  algebra," *J. Symbolic Computation* 28(4–5):777–799 (1999), DOI
  10.1006/jsco.1998.0339.
- W. Koepf, *Hypergeometric Summation*, 2nd ed., Springer Universitext (2014) —
  Gosper/Zeilberger/Petkovšek and their `q`-analogues.
- P. Paule & A. Riese, "A Mathematica q-analogue of Zeilberger's algorithm…,"
  *Fields Institute Communications* 14 (1997), 179–210 (the `qZeil` package).

**Complexity.** Moderate, **gated on a `QPochhammer` builtin** (prerequisite for
*any* symbolic `q`-product result). Once `QPochhammer` exists, the
factor-and-map path mirrors §3.3. A full `q`-Gosper/`q`-Zeilberger decision
procedure is a separate, larger project; the factor-and-map path covers explicit
`q`-products without it.

**Effectiveness.** `∏_{k=1}^{n}(1-q^k) = QPochhammer[q,q,n]`,
`∏(1+q^k) = QPochhammer[-q,q,n]`, Gaussian binomials. Requires `|q|<1` for the
infinite/convergence side; symbolic `q` stays a free variable.

---

## 5. Decision-procedure theory (the rigorous foundation)

The §3 algorithms are the *constructive read-off*; the papers below are the
*theory* that says when a closed form exists and makes the read-off canonical and
minimal. Mathilda already implements the additive side of the first one
(`sum_gosper.c`), so this is well-trodden ground.

### 5.1 Gosper–Petkovšek normal form (the engine of everything)

**Overview.** For a hypergeometric certificate `r(k) = t(k+1)/t(k)`, the
**Gosper-Petkovšek (GP) normal form** is the unique (up to constants) decomposition

```
r(k) = z · (a(k)/b(k)) · (c(k+1)/c(k)),   gcd(a(k), b(k+h)) = 1  ∀ h ∈ ℕ.
```

In the *additive* (Gosper summation) problem you then solve a degree-bounded
linear equation for a polynomial `x`; `t` is Gosper-summable iff it has a
solution. In the **multiplicative** problem the GP form is read off *directly*:
`c(k+1)/c(k)` contributes the rational shell `c(k)` to the anti-quotient
(telescoping, §3.2); the shift-coprime kernel `a(k)/b(k)`, factored into linear
pieces, contributes the **Pochhammer/Gamma factors** (§3.3); and `z` gives `z^n`.
So §3.2 + §3.3 *are* the GP read-off; building them on a proper GP-form pass is
the elegant unification. Computing the GP form needs the **dispersion set** (the
non-negative integer shifts `h` where `a(k)` and `b(k+h)` share a root — the
integer roots of `Res_k(a(k), b(k+h))` in `h`, or, more practically,
factor-difference matching).

**Reference.** R. W. Gosper, "Decision procedure for indefinite hypergeometric
summation," *Proc. Natl. Acad. Sci. USA* 75(1):40–42 (1978). Textbook: Petkovšek,
Wilf & Zeilberger, *A=B*, A K Peters (1996), ch. 5 (Gosper, where the GP form is
built) — note *A=B* is a *summation* text, the product reading is the dual.

**Complexity.** Mathilda has the additive engine (`sum_gosper.c`). The extra
work for products is the factorisation read-off (needs `src/poly/` factoring,
which the additive path avoids) and dispersion via resultant or factor-matching.
Moderate; the factor-matching route reuses §3.3's factoring.

**Effectiveness.** Decides whether a hypergeometric term has a hypergeometric
anti-quotient and constructs it. Closed-form class fixed to hypergeometric;
algebraic-extension roots are the main completeness risk; `q`/multibasic need the
§4 analogue.

### 5.2 Abramov–Petkovšek rational normal forms & minimal decomposition

**Overview.** Refines the GP form to a *canonical, minimal* multiplicative shape.
The **rational normal form (RNF)** writes the certificate's kernel as
**shift-reduced** (no numerator factor shift-equivalent to a denominator factor,
for *all* `k ∈ ℤ`, not just `k ≥ 0`), giving the **minimal multiplicative
representation** `t(n) = V(n) · ∏_k F(k)` where `∏ F` telescopes to the minimal
set of Γ/Pochhammer factors and `V` is the residual rational shell. The companion
minimal additive decomposition `t = ΔT₁ + T₂` makes "Gosper-summable iff `T₂=0`"
a clean criterion. The 2010 `(w,σ)`-canonical form reduces *minimality* to a
bipartite **assignment problem** (min-cost matching of shift-equivalent factors).

**References.**
- S. A. Abramov & M. Petkovšek, "Minimal decomposition of indefinite
  hypergeometric sums," *ISSAC 2001*, pp. 7–14.
- S. A. Abramov & M. Petkovšek, "Rational normal forms and minimal
  decompositions of hypergeometric terms," *J. Symbolic Computation*
  33(5):521–543 (2002); **erratum** *JSC* 38(3):1165 (2004).
- S. A. Abramov & M. Petkovšek, "Polynomial ring automorphisms, rational
  (w,σ)-canonical forms, and the assignment problem," *JSC* 45(6):684–708 (2010).

**Complexity.** Moderate-to-high, dominated by dispersion + ℚ-factorisation.
Primitives: univariate GCD (have), dispersion via subresultant resultants,
ℚ-factorisation (the hard dependency — have it in `src/poly/`), and (for the 2010
minimality) an `O(n³)` Hungarian assignment solver. Exact arithmetic only — no
floating dispersion. Mainly relevant if `Product` output must be *provably
minimal*; §3.3's integer-root matching already produces near-minimal output for
common inputs.

**Effectiveness.** Decides Gosper-summability and produces the minimal
multiplicative representation; underpins modern creative-telescoping. Strictly
univariate single-shift hypergeometric.

### 5.3 Petkovšek's `Hyper` algorithm

**Overview.** Finds all **hypergeometric solutions** of a linear recurrence with
polynomial coefficients `Σ p_i(n) y(n+i) = 0` — a finite basis (possibly empty).
Relevant because **every hypergeometric term *is* a product** (`y(n) =
y(0)∏ t(k)` with `t` rational), so `Hyper`'s GP-form factoring is the canonical
multiplicative decomposition, and the algorithm is what lets a CAS *recognise* a
product as a known hypergeometric closed form (or prove none exists). Enumerates
monic `a | p₀` and `b | p_r(n-r+1)`, solves a leading-degree auxiliary equation
for `z`, then a degree-bounded polynomial solve for `c`.

**Reference.** M. Petkovšek, "Hypergeometric solutions of linear recurrences
with polynomial coefficients," *J. Symbolic Computation* 14(2–3):243–264 (1992);
*A=B* ch. 8.

**Complexity.** High — dominated by **full factorisation over ℚ** (Yun → mod-`p`
→ Hensel lift → recombination), the single largest subsystem. The "Hyper shell"
(enumeration + auxiliary equation + GP coprimality) on top of a working
factoriser is comparatively small. Not needed for the §3 constructive path; only
for closed-form *recognition* of arbitrary recurrence-defined products.

**Effectiveness.** Complete decision procedure for hypergeometric (and,
self-layered, d'Alembertian) solutions; correctly reports "none of this class"
for non-hypergeometric sequences.

### 5.4 Karr / Schneider difference-field & difference-ring theory

**Overview.** The most general framework — the discrete analogue of Risch's
integration algorithm. A **difference field** `(𝔽, σ)` with shift `σ` builds a
tower of **Σ-extensions** (indefinite sums, `σ(t)=t+a`) and **Π-extensions**
(indefinite products, `σ(t)=α·t`), a **ΠΣ-field**, in which telescoping and
creative telescoping are *decidable* (Karr 1981). **Products are exactly the
Π-extensions**: the Gosper-for-products test is whether a candidate Π-extension is
genuinely new, i.e. whether *no* `g ∈ 𝔽` satisfies `σ(g)/g = α` (if one exists,
the product already has a closed form one level down). Schneider's **ΠΣ\***- and
**RΠΣ\***-extensions add depth-optimality and algebraic root-of-unity generators
(fixing Karr's inability to model `(-1)^k` in a field), and the Ocansey–Schneider
work canonically represents nested and `q`-/multibasic products and solves
product **zero-recognition**.

**References.**
- M. Karr, "Summation in finite terms," *J. ACM* 28(2):305–350 (1981); "Theory
  of summation in finite terms," *JSC* 1(3):303–315 (1985).
- C. Schneider, "Symbolic summation assists combinatorics," *Sém. Lothar.
  Combin.* 56 (2007), Art. B56b; "A refined difference field theory for symbolic
  summation," *JSC* 43(9):611–644 (2008); E. D. Ocansey & C. Schneider,
  "Representing (q-)hypergeometric products … in difference rings," *Springer
  Proc. Math. & Stat.* 226 (2018), 175–213.

**Complexity.** Research-grade — **out of scope for a pico-CAS.** A faithful
ΠΣ-tower needs multivariate difference-field arithmetic, recursive
constant-field/transcendence decisions at each adjunction, a parameterized
linear difference-equation solver with per-level bounds, and (for products)
cyclotomic + algebraic-independence machinery. Reference implementations
(Schneider's `Sigma`) are tens of thousands of lines on top of a full algebra
stack. The honest in-scope subset is "Gosper + its product analogue over ℚ(n)"
(§5.1) — the single-Π specialisation.

**Effectiveness.** Complete telescoping/creative-telescoping decision within the
constructed tower; canonical, algebraically-independent representations and
zero-recognition for nested products. "No closed form" means "none in this
tower," not absolute non-existence.

---

## 6. Definite & infinite products

### 6.1 Definite = endpoint evaluation of the indefinite anti-quotient

Once §3 yields an anti-quotient `P(n)`, a definite product is just
`∏_{k=a}^{b} f = P(b+1)/P(a)` with the bound `b` (symbolic or numeric) substituted
into the Gamma/Pochhammer form. The multiplicative fundamental theorem of
difference calculus; no extra algorithm. Care at Gamma poles (non-positive
integers) and exact-zero factors (a single zero factor forces the whole product
to 0 regardless of the closed form).

### 6.2 Infinite products via special functions and constants

The convergent `b→∞` limit of the finite closed form, recognised against a table
of known infinite-product identities (best implemented as `.m` rewrite rules,
mirroring `deriv.m` / the integral tables):

- **Wallis:** `∏ 4k²/(4k²-1) = π/2`.
- **Euler sine / cosine:** `sin(πz) = πz ∏(1 - z²/k²)`,
  `cos(πz) = ∏(1 - 4z²/(2k-1)²)`.
- **The `(k²+a²)/(k²+b²)` family:**
  `∏_{k≥1}(k²+a²)/(k²+b²) = b sinh(πa) / (a sinh(πb))`, with special cases
  `∏(1 - z²/k²) = sin(πz)/(πz)`, `∏(1 + 1/k²) = sinh(π)/π`,
  and (telescoping, index from `k=2`) `∏(1 - 1/k²) = 1/2`.
- **Gamma as an infinite product** (Euler/Weierstrass, DLMF §5.8) — mainly for
  numericalisation.
- **Hyperfactorial / Barnes `G` / Glaisher–Kinkelin:** `∏_{k=1}^{n} k^k =
  Hyperfactorial[n]`, `∏_{k=1}^{n-1} k! = BarnesG[n+1]`, with the
  Glaisher–Kinkelin constant `A = exp(1/12 - ζ'(-1)) ≈ 1.2824271291` as the
  asymptotic normaliser. **These require Barnes `G`, `Hyperfactorial`, and the
  Glaisher constant as prerequisites** (none exist in Mathilda yet — cf. the
  special-functions subsystem). Until then `∏ k^k` / `∏ k!` stay unevaluated.

**References.** DLMF §4.22 (sine product), §5.5 (reflection), §5.8 (Gamma
products), §17.2 (`q`-Pochhammer); Whittaker & Watson, *A Course of Modern
Analysis*, 4th ed. (1927); MathWorld "Hyperfactorial", "Barnes G-Function",
"Glaisher–Kinkelin Constant".

### 6.3 Convergence gate (must run before asserting an infinite value)

`∏(1+a_k)` converges to a *nonzero* limit iff `Σ a_k` converges (for eventually
one-signed `a_k`), absolutely iff `Σ|a_k| < ∞`. The CAS checklist:

1. **No zero factor** in range (else product is exactly 0).
2. **Term → 1** (`a_k = f(k)-1 → 0`); else divergent — reject. (`∏(1+1/k)`
   fails here.)
3. **Rational closed-form gate** for `f = P/Q`: require `deg P = deg Q`, **equal
   leading coefficients** (so `f→1`), **and equal next-to-leading coefficients**
   (root-sum of `P` = root-sum of `Q`). The third is decisive: it forces
   `a_k ~ c/k²` (convergent) rather than `a_k ~ c/k` (divergent). This is exactly
   why `∏ k/(k+1) → 0` and `∏(k+1)/k → ∞` diverge but `∏(k²+a)/(k²+b)` converges.

**References.** DLMF §1.9(vi); Whittaker & Watson §2.7; Rudin, *Real and Complex
Analysis*, 3rd ed., Thm 15.4–15.6; W. F. Trench, "Conditional convergence of
infinite products," *Amer. Math. Monthly* 106 (1999), 646–651 (the `Σa_k` +
`Σ|a_k|²` refinement).

---

## 7. Effectiveness summary (the benchmark inputs)

| Problem | Closed form | Stage |
|---|---|---|
| `Product[k, {k,1,n}]` | `n!` | §3.3 |
| `Product[k^2, {k,1,n}]` | `(n!)^2` | §3.3 |
| `Product[k+a, {k,1,n}]` | `Pochhammer[a+1, n]` | §3.3 |
| `Product[(2k-1)/(2k), {k,1,n}]` | `Binomial[2n,n]/4^n` | §3.3 |
| `Product[1+1/k, {k,1,n}]` | `n+1` | §3.2 |
| `Product[k/(k+1), {k,1,n}]` | `1/(n+1)` | §3.2 |
| `Product[1-1/k^2, {k,2,n}]` | `(n+1)/(2n)` | §3.2 |
| `Product[2^k, {k,1,n}]` | `2^(n(n+1)/2)` | §3.4 |
| `Product[1-a q^k, {k,0,n-1}]` | `QPochhammer[a,q,n]` | §4 |
| `Product[4k^2/(4k^2-1), {k,1,∞}]` | `π/2` (Wallis) | §6.2 |
| `Product[1-1/k^2, {k,2,∞}]` | `1/2` | §6.2 |
| `Product[1+1/k^2, {k,1,∞}]` | `Sinh[π]/π` | §6.2 |
| `Product[(k^2+a^2)/(k^2+b^2), {k,1,∞}]` | `b Sinh[πa]/(a Sinh[πb])` | §6.2 |
| `Product[k^k, {k,1,n}]` | `Hyperfactorial[n]` | §6.2 † |
| `Product[1+1/k, {k,1,∞}]` | divergent (rejected) | §6.3 |
| `Product[k!+1, {k,1,n}]` | unevaluated (non-hypergeometric) | — |

† requires the `Hyperfactorial`/`BarnesG`/Glaisher prerequisites of §6.2.

---

## 8. Recommended staging for Mathilda

`Product` should mirror `Sum`'s architecture: a dispatcher (`product.c`) plus one
file per algorithm in a new `src/product/`, each algorithm also exposed as a
context-qualified builtin (``Product`Rational``, ``Product`Gosper``, …), exactly
as `Sum` does (`Sum`Gosper`` etc.) and `Integrate` does
(`Integrate`BronsteinRational``). `Product` is `HoldAll | Protected`; the index
is `Block`-localised. Adding a stage stays additive: new `src/product/product_*.c`,
one `try_*` line in the dispatch cascade, one `*_init()` call in `product_init()`.

| Stage | Algorithm | Effort | Depends on |
|---|---|---|---|
| 0 | iterator surface, finite unroll, constant/trivial, Karr empty-range, multi-index nesting, Method cascade (§3.1) | Low | evaluator, `Block` |
| 1 | multiplicative telescoping (§3.2) | Low | poly GCD/shift |
| 2 | **rational → Pochhammer/Gamma** (§3.3) — *the workhorse* | Moderate | `src/poly/` factoring, `Pochhammer`, `Gamma` |
| 3 | geometric / poly-exponential factors (§3.4) | Low | reuse `Sum` |
| 4 | infinite-product identity table + convergence gate (§6.2–6.3) | Moderate | `.m` rules; `Sin`/`Sinh`/`π` |
| 5 | `q`-products (§4) | Moderate | **new `QPochhammer` builtin** |
| 6 | Barnes `G` / Hyperfactorial / Glaisher results (§6.2 †) | Moderate | **new `BarnesG`, `Hyperfactorial`, Glaisher** |
| 7 | GP-form unification / Abramov–Petkovšek minimality (§5.1–5.2) | High | dispersion, factoring |
| — | full Karr/Schneider ΠΣ towers (§5.4) | Research-grade | out of scope |

**Stages 0–3 deliver the bulk of practical `Product`** on machinery Mathilda
already has, and stand alone. Stages 4–6 each unlock a recognised family but gate
on a new prerequisite (an infinite-product rule table, `QPochhammer`, Barnes
`G`/Glaisher respectively). Stage 7 is the "provably minimal output" refinement.
Stage 5.4 is explicitly out of scope for a pico-CAS.

### The recurring hard primitive

Across every Gamma-producing stage the gating dependency is **univariate
polynomial factorisation over ℚ** (and, for completeness, over algebraic
extensions). Mathilda's `src/poly/` already provides this, including an
algebraic-number layer — the single most leverageable asset for the whole
subsystem. The additive Gosper path uniquely *avoids* factoring (GCD + linear
solve only); everything that reads off Γ/Pochhammer factors depends on it.
Secondary engineering risks: exact dispersion-set computation (never floating
point) and intermediate coefficient swell in resultants (keep contents factored,
work with primitive parts).

---

## 9. Verification notes & caveats from the source survey

- **Moenck (1977) is a *summation* paper**, not the product→Gamma algorithm — it
  is widely mis-cited for products. Cite Koepf (arXiv:math/9412227) and the
  Gosper–Petkovšek lineage for the product construction.
- *A=B* (Petkovšek–Wilf–Zeilberger, 1996) is a *summation* text; the product
  results are the multiplicative dual, not a dedicated chapter.
- There is **no single canonical "multiplicative Gosper" paper** — it is the
  product specialisation of the GP normal form; the rigorous modern treatments
  are in the Schneider/Bauer–Petkovšek (multibasic) literature.
- `∏_{k≥2}(1-1/k²) = 1/2` is **telescoping**, not the sine product — it must
  index from `k=2` (the sine product at `z=1` has a vanishing `k=1` factor).
- The Jacobi triple product and `q`-series have multiple incompatible conventions
  (`q` vs `q²` in the exponent); a CAS must fix and document one.
- All Gamma/Pochhammer identities use the **rising-factorial** convention
  `(a)_n = Γ(a+n)/Γ(a)`; combinatorics texts use `(a)_n` for the *falling*
  factorial — fix the convention internally to avoid shift/sign bugs.

---

## 10. References (consolidated)

**Decision procedures**
- R. W. Gosper, "Decision procedure for indefinite hypergeometric summation,"
  *PNAS* 75(1):40–42 (1978).
- M. Petkovšek, "Hypergeometric solutions of linear recurrences with polynomial
  coefficients," *J. Symbolic Computation* 14(2–3):243–264 (1992).
- M. Petkovšek, H. S. Wilf, D. Zeilberger, *A=B*, A K Peters (1996).
- S. A. Abramov & M. Petkovšek, "Rational normal forms and minimal decompositions
  of hypergeometric terms," *JSC* 33(5):521–543 (2002); erratum *JSC* 38(3):1165
  (2004); "Minimal decomposition of indefinite hypergeometric sums," *ISSAC 2001*,
  7–14; "Polynomial ring automorphisms, rational (w,σ)-canonical forms, and the
  assignment problem," *JSC* 45(6):684–708 (2010).
- M. Karr, "Summation in finite terms," *J. ACM* 28(2):305–350 (1981); "Theory of
  summation in finite terms," *JSC* 1(3):303–315 (1985).
- C. Schneider, "Symbolic summation assists combinatorics," *Sém. Lothar.
  Combin.* 56 (2007), Art. B56b; "A refined difference field theory…," *JSC*
  43(9):611–644 (2008); E. D. Ocansey & C. Schneider, "Representing
  (q-)hypergeometric products … in difference rings," *Springer Proc. Math. &
  Stat.* 226 (2018), 175–213.

**Rational/Gamma construction & `q`-analogues**
- R. Moenck, "On computing closed forms for summations," *Proc. 1977 MACSYMA
  Users' Conf.*, 225–236 (⚠ summation, not products).
- W. Koepf, "Algorithms for the indefinite and definite summation,"
  arXiv:math/9412227 (1995); *Hypergeometric Summation*, 2nd ed., Springer (2014).
- H. Böing & W. Koepf, "Algorithms for q-hypergeometric summation in computer
  algebra," *JSC* 28(4–5):777–799 (1999).
- P. Paule & A. Riese, "A Mathematica q-analogue of Zeilberger's algorithm…,"
  *Fields Inst. Commun.* 14 (1997), 179–210.

**Textbooks & tables**
- R. Graham, D. Knuth, O. Patashnik, *Concrete Mathematics*, 2nd ed.,
  Addison-Wesley (1994).
- G. Andrews, R. Askey, R. Roy, *Special Functions*, CUP (1999).
- E. Whittaker & G. Watson, *A Course of Modern Analysis*, 4th ed., CUP (1927).
- W. Rudin, *Real and Complex Analysis*, 3rd ed., McGraw-Hill (1987).
- W. F. Trench, "Conditional convergence of infinite products," *Amer. Math.
  Monthly* 106 (1999), 646–651.
- NIST DLMF, §§1.9(vi), 4.22, 5.5, 5.8, 17.2.

**Reference systems**
- Wolfram Language `Product` reference.
- Maple `product` / `SumTools` documentation.
- Maxima manual, "Sums, Products and Series."
- SymPy `sympy/concrete/products.py` (`_eval_product`); mpmath `nprod`.

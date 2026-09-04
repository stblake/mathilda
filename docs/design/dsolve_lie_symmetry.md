# DSolve`LieSymmetry — heuristic Lie point-symmetry method (M10)

`DSolve`LieSymmetry` integrates a **first-order** ODE `y' = ω(x, y)` by finding a
one-parameter Lie group of point symmetries and reducing the equation to a
quadrature. It is the general first-order backstop in the `DSolve` cascade: it
runs after the deterministic specialists (which give tidier closed forms) and
before the implicit / series fallbacks.

It is the **one deliberately heuristic method** in an otherwise
decision-procedure cascade. For a first-order ODE the symmetry-determining
equation is a *single* PDE in *two* unknown functions `ξ(x,y)`, `η(x,y)` — it is
underdetermined, has infinitely many solutions, and finding a *useful* one is as
hard as solving the ODE. There is no algorithm; SymPy and Maple both use a fixed
table of ansätze (functional forms) for `ξ, η`, each of which collapses the
determining PDE to something solvable. This note follows that table.

## 1. The symmetry (determining) condition

The generator `V = ξ ∂ₓ + η ∂_y` is a point symmetry of `y' = ω(x,y)` iff its
first prolongation annihilates `y' − ω`, which reduces to the linearized
invariance condition

```
S(ξ, η) ≡ η_x + (η_y − ξ_x) ω − ξ_y ω²  − ξ ω_x − η ω_y  =  0.        (†)
```

(Subscripts are partial derivatives; `ω_x = ∂ω/∂x`, `ω_y = ∂ω/∂y`.) A candidate
`(ξ, η)` is accepted only when `S` is provably zero (`zero_test_decide`) — the
`checkinfsol` gate. In the implementation `y[x]` is replaced by a plain symbol
`Y` and `x` is the independent-variable symbol, so `(†)` is an ordinary algebraic
expression differentiated with `D`.

## 2. Integrating the ODE from a known symmetry

Two classical routes; the implementation uses the first.

**Integrating factor (used).** Lie's theorem: if `(ξ, η)` is a symmetry of
`y' = ω`, then

```
μ(x, y) = 1 / (η − ω ξ)                                                (‡)
```

is an integrating factor of the 1-form `ω dx − dy = 0`. The condition `(†)` is
exactly the exactness condition `∂(μω)/∂y = ∂(−μ)/∂x`, so a first integral
`F(x, y) = C₁` exists with `F_x = μ ω`, `F_y = −μ`. Build it by the standard
exact-equation potential:

```
Φ  = ∫ μ ω  dx           (y held constant)
F  = Φ + ∫ ( −μ − ∂Φ/∂y ) dy
```

Return `F(x, y[x]) == C[1]` through `dsolve_run_implicit` (which verifies by
implicit differentiation `y' == −F_x/F_y` and fits an IVP constant), then attempt
`Solve[F == C[1], y]` for an explicit `y[x] -> …` branch; keep the implicit form
when it does not invert. Decline the branch if either integral is non-elementary
(`ds_has_head("Integrate")`) or `η − ω ξ ≡ 0` (the symmetry is trivial /
tangent to the direction field).

**Canonical coordinates (fallback, not yet used).** `r` = group invariant
(`ξ r_x + η r_y = 0`, i.e. `dx/ξ = dy/η`) and `s = ∫ dx/ξ`; in `(r, s)` the ODE
becomes `ds/dr = f(r)`, a quadrature. Heavier (two PDEs) — only needed when `(‡)`
yields a non-integrable `μ`.

## 3. The heuristics (ansatz table)

Each heuristic substitutes its ansatz into `(†)`, then splits the result into a
*determining system* by collecting coefficients in powers of `y` (or of `ω`) with
`Coefficient` / `Collect`, giving ODEs / algebraic equations for the unknown
functions, solved by `Solve` / `Integrate`. `ω` is first put over a common
denominator (`Together`) so `(†)` clears to a polynomial identity.

| Heuristic | Ansatz `(ξ, η)` | Source |
|---|---|---|
| `abaco1_simple` | one nonzero, one variable: `(0, f(x))`, `(0, f(y))`, `(f(x), 0)`, `(f(y), 0)` | CT–Duarte–da Mota, CPC 101 (1997), p.8 |
| `abaco1_product` | `(f(x)·g(y), 0)`, `(0, f(x)·g(y))` | CT–Roche, CPC 113 (1998), pp.7–8 |
| `function_sum` | `(f(x)+g(y), 0)`, `(0, f(x)+g(y))` | CT–Roche, pp.7–8 |
| `abaco2_similar` | `(f(x), g(x))`, `(g(y), f(y))` | CT–Roche, pp.10–12 |
| `linear` | `(a x + b y + c, e x + f y + g)` → linear algebraic system in the constants | CT–Kolokolnikov, math-ph/0007023 |
| `bivariate` | `ξ, η` bivariate polynomials in `x, y`, degree-bounded → linear system | Olver, *Applications of Lie Groups to Differential Equations*, pp.327–329 |
| `chi` | `η = ξ ω + χ(x, y)`, `χ` solving an auxiliary linear PDE | CT–Duarte–da Mota, p.8 |
| `abaco2_unique_unknown` | forms for `ω` carrying an unknown function / non-integer power | CT–Roche, pp.10–12 |
| `abaco2_unique_general` | `(g(y), f(x))` with no assumption on `ω` | CT–Roche, pp.10–12 |

Heuristics are tried in the order above; the first that yields a `(ξ, η)` passing
the `checkinfsol` gate wins (matching SymPy's speed-first `default`). Note that
`abaco1_simple` largely overlaps the earlier `Separable` / `LinearFirstOrder`
methods (a one-variable symmetry ⇔ separable or linear), so in the automatic
cascade Lie earns its keep mainly through `linear`, `abaco2_*`, and `bivariate`,
which catch equations the specialists miss (e.g. the linear-coefficients family
`y' = (a₁x+b₁y+c₁)/(a₂x+b₂y+c₂)`).

## 4. Implementation staging

- **L1** — substrate (`lie_symmetry_residual`, `lie_check_infinitesimals`,
  `lie_first_integral`), the integrating-factor → first-integral →
  `dsolve_run_implicit` (+ explicit `Solve`) pipeline, and `abaco1_simple`. Proves
  the pipeline end-to-end via the pinned `DSolve`LieSymmetry[eqn, y, x]` builtin.
- **L2** — `linear` ✅ (the coefficient-splitting determining-system machinery);
  `abaco1_product` ✅; `abaco2_similar` ✅ (§4.3, below); `function_sum` ✅ (§4.2, below).
- **L3** — `bivariate` ✅ (degree-bounded: `lie_poly_symmetry` tries degree 2 then
  3, the exact generalization of `linear`'s degree-1 NullSpace determining system —
  builds `ξ, η` as general bivariate polynomials over interned `DSolve\`lieB*`
  coefficients, then `Numerator[Together[S]]` → `CoefficientList[·,{x,y}]` →
  `Outer[Coefficient, forms, coeffs]` → `NullSpace`; catches quadratic/projective
  symmetries such as `ξ=x², η=x y` that the affine ansatz misses);
  `abaco2_unique_unknown` ✅ (§4.4.1, below, incl. the order-zero extension Eqs 73–81);
  `chi` ✅ (CPC 101, 5th algorithm; §4.5, below). The formal §4.4.2 Case I/II
  (`abaco2_unique_general`) is a **documented exemption** — the authors call the
  closed-form route "very inefficient, if not just impractical", and every
  `[F(x),G(y)]`-symmetric ODE it targets is already caught by Riccati/Bernoulli/kernel
  methods, so it cannot be tested non-vacuously.

### 4.1 `abaco1_product` (CT–Roche §4.1) — the first quadrature ansatz

Symmetry `[ξ = F(x) G(y), η = 0]` and its inverse `[0, F(x) G(y)]`. Necessary
condition (their Eq 19): `L = (ω_xy ω − ω_x ω_y)/ω⁴` is product-separable in `x, y`;
then `F(x)` is its x-factor and `g(y) = F · ∂_x(1/(F ω))` (Eq 20) must be free of
`x`, giving `G = Exp[∫ g dy]`. The inverse pattern is found by running the same
extractor on the inverse ODE `y' = 1/ω(y, x)` and mapping the symmetry back
(`[swap(η), swap(ξ)]`), so one routine covers both. Reaches rational-but-non-
polynomial infinitesimals (`ξ = y/x`, …) that `linear`/`bivariate` cannot.

Shared substrate introduced here and reused by the later quadrature heuristics:
- `lie_sep_xfactor(L)` — the product-separable x-factor `Exp[∫ (L_x/L) dx]`, gated by
  a **fast rational free-of test** `lie_free_of_var` (for `P/Q`, free of a variable
  ⇔ the polynomial `P_v Q − P Q_v` expands to 0 — a polynomial zero-test, ~2 ms,
  vs ~12 ms for the general `zero_test` on the rational derivative).
- `lie_ratsimp` — `Cancel[Together[·]]`, used in place of the far costlier `Simplify`
  on the hot path (all these forms are rational).
- `lie_swap_xy` / `lie_inverse_omega` — the `x↔y` swap and the inverse ODE.

The Lie sub-cascade is ordered **cheapest-first**: `abaco1_simple` → `linear` →
`abaco1_product` → `bivariate`, so the degree-2/3 NullSpace is a last resort. Every
returned first integral back-substitutes to zero (`dsolve_run_implicit`), so a
mis-derived candidate can only decline, never mislead. **Validation:** SymPy 1.14's
`lie_group` (`infinitesimals`) *times out* (>25 s) on the pure product family
`2xy/(x²+2y⁴+c)` that this solves in ~30 ms; it agrees with SymPy on the ODEs SymPy
can handle. Per-solve valgrind is leak-flat (identical `definitely lost` for `1+1`,
1×, and 8×).

### 4.2 `function_sum` (CT–Roche §4.2) — the additive quadrature ansatz

Symmetry `[ξ = F(x) + G(y), η = 0]` and its inverse `[0, F(x) + G(y)]` — the additive
counterpart of `abaco1_product`. The classifying quantity is the **rational** factor
`ω · ∂²ₓ(1/ω) = F''(x)/(F(x)+G(y))` (their Eq 27; the leading `ω` is essential — it
cancels the transcendental part of `1/ω`, which for a genuine member is Log/ArcTan, so
`∂²ₓ(1/ω)` alone still carries it). `∂_y` of its reciprocal is product-separable
(Eq 28) with x-factor `1/F''(x)` (via `lie_sep_xfactor`), from which
`ξ = F+G = 1/(xfactor · factor)`; `factor ≡ 0` marks the invert-linear case (Eq 27
footnote) and declines. The inverse pattern reuses the inverse-ODE driver
(`lie_run_with_inverse`). Note the paper (§3) remarks these patterns rarely help
Kamke-book ODEs — the value is completeness of the L2 ansatz table, not new coverage of
named equations. **Validation:** `t_lie_function_sum` and the `t_stress_lie_function_sum`
Log-form family, each a member of the §4.2 invariant family built with `F=1/x`, `G=y`;
per-solve valgrind leak-flat.

### 4.3 `abaco2_similar` (CT–Roche §4.3) — both components single-variable

Symmetry `[ξ = F(x), η = H(x)]` and its inverse `[F(y), H(y)]` — the first ansatz
where *both* infinitesimals are one-variable functions. From `Q = ω_y/ω_yy` (Eq 39),
`T = Q_x/Q_y` must be free of `y` (Eq 44a) and the integrand
`(T·ω_y − T_x − ω_x)/(ω + T)` must be free of `y` (Eq 44b); then
`F = Exp[∫(integrand) dx]` (Eq 43) and `H = −T·F` (Eq 40). The rarer `Q_y = 0`
sub-case (Eqs 47/51) is deferred (declined). The inverse pattern reuses the
inverse-ODE driver.

The point of this heuristic is **reach**, not a new mechanism: it is the first to
solve ODEs with an **irrational** `ω`, where every rational ansatz
(`abaco1_simple`/`linear`/`abaco1_product`) structurally declines. It handles the
substitution family `y' = (a x + b y + c)^p` for non-integer `p` — e.g.
`y' = Sqrt[x + y]` returns the verified first integral
`x − 2 Sqrt[x+y] + 2 Log[1 + Sqrt[x+y]] == C[1]` (the antiderivative uses a CRC
integral-table rule, so this path needs `init.m` loaded — which the DSolve test
suites now do, matching production). Sits after `abaco1_product` and before
`bivariate` in the cheapest-first sub-cascade. **Validation:** `t_lie_abaco2_similar`
+ the `t_stress_lie_similar` radical family; per-solve valgrind leak-flat (the
`ds_simplify` finalize step exposed — and this work fixed — a pre-existing
`simp_power.c` power-distribution leak, unrelated to DSolve).

### 4.4.1 `abaco2_unique_unknown` (CT–Roche §4.4.1) — from functions/non-integer powers

Symmetries `[ξ=F(x), η=G(y)]` and `[ξ=G(y), η=F(x)]` (not inverses of each other; one
scheme finds both). Prop 7: if `ω` has such a symmetry and contains a function `M` of
both `x` and `y`, that `M` must be a function of `f(x)+g(y)`, so `R = M_y/M_x = g_y/f_x`
(Eq 63; the outer function's derivative cancels) **separates by product**. A recursive
tree-walk `lie_collect_kernels` gathers each candidate mapping `M` — a non-arithmetic
function application, or a power with a non-integer / variable exponent, that holds both
variables. For each, with `X` the x-factor of `R` (via `lie_sep_xfactor`), scheme (iib)
tries the two candidates `[X, −X/R]` (type `[F(x),G(y)]`) and `[−R/X, 1/X]` (type
`[G(y),F(x)]`), each gated by `lie_check`. A purely rational `ω` has no kernels → instant
decline. **Reach vs. limit:** it solves ODEs with a non-integer power of both variables —
`y' = (x/y)(x²+y²)^p` → symmetry `[1/x, −1/y]`, first integral elementary. For a genuinely
*arbitrary* function `F` (Kamke-85 form `y'=(x/y)F[(x²+y²)/2]`) the symmetry is still
found but the Lie integrating-factor quadrature `∫F/(1+F)` is non-elementary, so the
branch declines at integration (correctly — no inert head). **Validation:**
`t_lie_abaco2_unique_unknown` and the `t_stress_lie_unique_unknown` non-integer-power
family; per-solve valgrind adds only the pre-existing FLINT `rat_canon` leak
(`tasks/flint_ratcanon_leak.md`), its own allocations balanced. `chi` is now implemented
(§4.5); the general case §4.4.2 (Cases I/II) is a documented exemption (see §4).

### 4.4.1-general `abaco2_unique_unknown` — differential invariant of order zero (Eqs 73–81)

Beyond the two separable candidates, for each mapping `M` with `R = M_y/M_x` the same
scheme admits — *without* a separability test — the "order-zero" candidates
`[-R, 1]` (Eq 75, pattern `[f(x)g(y), 1]`) and `[1, -R]` / `[1, -1/R]` (the family-77
patterns). These reach ODEs whose `R` does **not** separate by product, e.g. Kamke's
first-order ODE 433, `(x y' + y + 2x)² = 4(x y + x² + a)`: the mapping
`M = Sqrt[x y + x² + a]` gives `R = x/(2x+y)` (not separable), and the order-zero
candidate `[1, -R]` yields the first integral `x − Sqrt[x² + x y + a] == C[1]`. Each
candidate is gated by `lie_check` and integrated by `lie_first_integral` exactly as the
separable ones, so a mis-derived candidate can only decline. This is the practically
useful part of the `[F(x),G(y)]`-general territory; the *formal* §4.4.2 Case I/II route
(closed-form `η`/`ξ` from derivatives of `φ = Log ω`) is, in the authors' own words,
"very inefficient, if not just impractical" and rarely helps a named ODE, so it is not
implemented (documented exemption).

### 4.5 `chi` (CPC 101 1997, 5th algorithm) — the `η = ξω + χ` reformulation

The last and richest heuristic. Because the determining residual is linear and its
tangent solution `η = ξω` is trivial, `S(ξ, ξω + χ) = S(0, χ)`, so a `χ(x,y)` solving
the **linear PDE**

```
χ_x + ω χ_y − ω_y χ = 0                                                    (Eq 10)
```

gives the symmetry `[0, χ]` for *any* `ξ`. `χ` is sought by a **rich-basis** ansatz:

```
χ = ( Σ_k c_k m_k(x,y) ) / Dtrans
```

where the `m_k` are monomials (up to total degree 2 then 3, each times a `{1, x, y}`
factor) in the **transcendental atoms** of `ω` — the `Sin/Cos` and `Sinh/Cosh` base
pairs, `Log`, the inverse-trig functions — and `Dtrans` is the product of those atoms
(supplying the reciprocal / denominator structure a plain polynomial cannot reach). The
`c_k` are undetermined constants. Substituting into Eq 10 and:

1. **clearing** — multiplying by `x^(d+3) ∏ atom^(d+3)`, which both clears every
   denominator *and* folds the reciprocal trig `D` produces (`Csc·Sin → 1`,
   `Cot·Sin → Cos`, …) back to `Sin/Cos`, so the subsequent substitution never needs a
   reciprocal (which would send `Together` into a blow-up — the naive route measured a
   timeout);
2. **substituting** each atom by an independent generator, then `Numerator[Together[·]]`
   and a `PolynomialQ` gate (a leftover `Exp` / special-function makes it non-polynomial
   → decline);
3. **reducing** modulo the Pythagorean relation of each pair (`gc² → 1 − gs²`,
   `gch² → 1 + gsh²`) so identity-based cancellations are found;
4. `CoefficientList` over `{x, y, generators}` → `Outer[Coefficient, ·, {c_k}]` →
   `NullSpace`

gives a homogeneous linear determining system whose null space is a basis of `χ`'s.
Each nonzero `χ` is integrated by `lie_first_integral`, **simplest-first** (a smaller
`χ` has a simpler, faster, more-likely-elementary quadrature; a messy sibling can send
`Integrate` into a long trig search). `lie_check` is skipped here — the null vectors are
exact Eq-10 solutions by construction, and the returned first integral is still verified
by `dsolve_run_implicit`'s back-substitution.

This is the one heuristic whose `χ` may be a **genuine transcendental** beyond the
polynomial reach of `bivariate` — e.g. Kamke's first-order ODE 357,
`x y' ln(x) sin(y) + cos(y)(1 − x cos(y)) = 0`, has the symmetry
`[0, cos²(y)/(ln(x) sin(y))]`, from which `chi` returns the verified first integral
`−x + Log[x] Sec[y[x]] == C[1]`. First cut: elementary-transcendental atoms only (an
undefined function → decline, the `abaco2_unique_*` domain; `Exp` and special functions
are not folded, so the `PolynomialQ` gate declines them); atom degree ≤ 3, `m ≤ 5`
atoms, ≤ 90 monomials (cost caps). Since `chi` is built for the trig case, the
cascade **skips the rational/algebraic classifiers** (`abaco1_product`, `function_sum`,
`abaco2_similar`) when `ω` contains trig — they cannot solve it and their construction
is slow on trig — and `abaco1_simple` uses the bounded `lie_ratsimp` in place of the
full `Simplify`. Runs **last** (richest, biggest system). Source: Cheb-Terrab, Duarte &
da Mota, CPC 101 (1997), Eqs 9–10 and the 5th algorithm.

### 4.6 Robustness on transcendental / undefined-function `ω`

The quadrature heuristics build their classifiers by repeated symbolic differentiation
of `ω`; on an `ω` carrying an **undefined function** of both variables — e.g.
`y' = Tan[ArcTan[y] + F[x²+y²]]` — those intermediates balloon (a measured pre-`Expand`
numerator hit 217 k nodes, blowing past a million after `Expand`) and the general
`zero_test` / `Integrate` engine hangs on the arbitrary-function atoms. This used to
hang `function_sum` before `abaco2_unique_*` could run. The method is now bounded:

- **Node budget.** `lie_ratsimp`, `lie_free_of_var`, `lie_is_zero`, `lie_sep_xfactor`
  bail the instant an input crosses `LIE_EXPR_BUDGET` (6000 nodes) — a genuine
  rational/algebraic target keeps every intermediate to a few hundred nodes, so the
  budget only ever fires on a blow-up. A deterministic node budget is the
  machine-independent analogue of the wall-clock timeout SymPy/Maple use here.
- **Polynomial zero-test.** The fast helpers decide "is this zero" by `Expand` +
  a structural literal-`0` test (`lie_lit_zero`), never the general `zero_test`
  (which hangs on undefined-function atoms). `lie_check` keeps the full `zero_test`
  as a fallback for identity-based zeros, but only when the residual is free of
  undefined functions.
- **Undefined-function gate.** `abaco1_product` / `function_sum` / `abaco2_similar`
  (rational/algebraic ansätze that structurally cannot solve an arbitrary-function
  ODE) are skipped when `ω` has an undefined function; `lie_first_integral` declines
  such `ω` up front (its Lie quadrature is non-elementary anyway).

Effect: `y' = Tan[ArcTan[y] + F[x²+y²]]` and the CT–Roche Eq-70 example
`y' = -Tan[ArcTan[x/y] + H[x²+y²]]` decline in ~1 s instead of hanging, with no inert
head and no wrong answer (verified by `t_lie_undefined_function_declines`).

## 5. References

- E. S. Cheb-Terrab, A. D. Roche, *Symmetries and First Order ODE Patterns*,
  Comput. Phys. Commun. 113 (1998) 239. Preprint: UW CS-98-11.
- E. S. Cheb-Terrab, L. G. S. Duarte, L. A. C. P. da Mota, *Computer Algebra
  Solving of First Order ODEs Using Symmetry Methods*, Comput. Phys. Commun. 101
  (1997) 254.
- E. S. Cheb-Terrab, T. Kolokolnikov, *First-order ordinary differential
  equations, symmetries and linear transformations*, Eur. J. Appl. Math. (2003);
  arXiv:math-ph/0007023.
- E. S. Cheb-Terrab, L. G. S. Duarte, L. A. C. P. da Mota, *Computer Algebra
  Solving of Second Order ODEs Using Symmetry Methods*, Comput. Phys. Commun. 108
  (1998) 90; arXiv:gr-qc/9703082 (background; the algorithmic 2nd-order route, not
  used by this first-order method).
- P. J. Olver, *Applications of Lie Groups to Differential Equations*, GTM 107.
- J. Starrett, *Solving Differential Equations by Symmetry Groups*, Amer. Math.
  Monthly 114 (2007) 778 (pedagogical; cited by SymPy's docs).
- Cross-reference: SymPy `sympy/solvers/ode/lie_group.py`
  (`lie_heuristic_*`, `infinitesimals`, `checkinfsol`).

Working PDF copies of the freely-available preprints are in `lie_refs/`
(git-ignored — not redistributed with the source).

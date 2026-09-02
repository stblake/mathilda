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
  `abaco1_product` ✅; `function_sum`, `abaco2_similar` remain.
- **L3** — `bivariate` ✅ (degree-bounded: `lie_poly_symmetry` tries degree 2 then
  3, the exact generalization of `linear`'s degree-1 NullSpace determining system —
  builds `ξ, η` as general bivariate polynomials over interned `DSolve\`lieB*`
  coefficients, then `Numerator[Together[S]]` → `CoefficientList[·,{x,y}]` →
  `Outer[Coefficient, forms, coeffs]` → `NullSpace`; catches quadratic/projective
  symmetries such as `ξ=x², η=x y` that the affine ansatz misses); `chi`,
  `abaco2_unique_unknown`, `abaco2_unique_general` remain.

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

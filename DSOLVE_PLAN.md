# DSolve Implementation Plan

`DSolve` is Mathilda's symbolic differential-equation solver — the symbolic
analogue of `Integrate`. It is a **cascade polyalgorithm in C**: a dispatcher
(`builtin_dsolve`, `src/calculus/dsolve.c`) tries a sequence of methods until one
succeeds; each method is *also* a REPL-callable builtin `DSolve`<Method>[...]`
with its own docstring and attributes, exactly like the `Integrate`<Method>`
family.

- **Language:** C throughout, one file per method (`src/calculus/dsolve_<m>.c`).
  No `.m` tier; the special-function recognizer table (Phase 1c) is a static C
  table.
- **Missing special functions:** recognizers emit the head inertly (e.g.
  `MathieuC[...]`, `WeierstrassP[...]`) even when the function has no evaluator.
  Inert-head solutions cannot be verified by back-substitution, so their
  correctness rests on recognizer structure + the stress corpus.
- **Scope:** Phase 1 = ODEs, Phase 2 = PDEs. DAEs, delay DEs, integral /
  integro-differential, and hybrid (`WhenEvent`) equations are future work.
- **Determinism:** every method is a deterministic decision procedure *except*
  `LieSymmetry` (M10), which is heuristic by necessity — for a first-order ODE the
  linearized symmetry condition is one PDE in two unknowns (`ξ`, `η`), i.e.
  underdetermined, so no decision procedure exists and a fixed table of ansätze
  (Cheb-Terrab et al., as in SymPy/Maple) is the state of the art.

## Architecture

| Piece | File | Role |
|---|---|---|
| Dispatcher + cascade + fail-memo + depth + init | `dsolve.c` | mirrors `integrate.c` |
| Problem substrate (parse / verify / fit / assemble + helpers) | `dsolve_common.{c,h}` | `DSolveProblem`, `dsolve_run` |
| One method each | `dsolve_<method>.c` | `try` + `builtin` + `init` |

Each method file exposes the three-function contract:
`dsolve_<m>_try(P, &nbranch)` (cascade routine, returns branch bodies or NULL),
`builtin_dsolve_<m>(res)` (REPL entry, via `dsolve_method_builtin`), and
`dsolve_<m>_init()` (registers `DSolve`<Name>` with
`ATTR_PROTECTED` + docstring; chained from `dsolve_init`).

`DSolve` is **not** `HoldAll` (attributes `Protected`, matching
Mathematica): with the dependent function undefined, a symbolic equation
(`y'[x] == a Sin[x]`) and its point conditions (`y[0] == 5`) evaluate to
unevaluated `Equal[...]` on their own, so they reach the solver formal without
holding — and an equation stored in a variable (`eq = …; DSolve[eq, y, x]`) is
solved rather than declined. The parser
folds any `D[...]` into `Derivative[...]`, splits equations from point
conditions, detects the output form (`u` → `Function`, `u[x]` → expression), and
parses the options `GeneratedParameters` (default `C`), `Assumptions`, `Method`,
`IncludeSingularSolutions`. Every returned branch is verified by substituting
`u -> Function[{x}, body]` into each residual and requiring the result not to be
a *decidable* non-zero (undecidable is kept, matching Solve); constants are then
fitted to conditions via `Solve`, and `C[k]` renamed to the requested head last.

Reuse: `solvepoly_solve_polynomial_equality` (characteristic roots), `Eigenvalues/
Eigenvectors/JordanDecomposition` + `MatrixPower`/`DiagonalMatrix` (systems — the
fundamental matrix `e^{Ax}` is assembled from the Jordan form, as symbolic
`MatrixExp` is currently inert), `Series`/`SeriesData` (Frobenius),
`Eliminate` (PDE characteristics), `Piecewise`/`UnitStep` (BVP/wave/heat),
`zero_test_decide` (verify), `Integrate`/`D`/`Solve` (throughout).

## Milestones

- **M0 — substrate + first-order core.** ✅ DONE. Dispatcher, `dsolve_common`,
  `DSolve`Quadrature`, `DSolve`LinearFirstOrder`, `DSolve`Separable`, IVP fitting,
  both output forms, verification, fail-memo, tests (`tests/test_dsolve.c`).
- **M1 — rest of first-order (1a).** ✅ DONE. Bernoulli, Homogeneous, Exact
  (with μ(x)/μ(y) search), Clairaut (with singular solutions). Added the shared
  helpers `dsolve_linear_factor_solve`, `dsolve_algebraic_residual`,
  `dsolve_extract_solutions`, and robust `ds_free_of` (derivative-based
  free-of-variable test) / `ds_simplify`.
- **M2 — linear constant-coefficient (1b) + BVP.** ✅ DONE. `DSolve`LinearConstantCoefficients`
  (any order; real/complex/repeated roots; variation of parameters for forcing);
  BVP constant-fitting works through the substrate (`y''+y==0, y[0]==0, y[π/2]==1`
  → `Sin[x]`).
- **M3 — variable-coefficient (1c).** IN PROGRESS. `[✓] EulerCauchy` done
  (matches the reference `x²y''+4xy'+7y==0`). Also refactored the linear-ODE
  substrate into shared helpers `dsolve_linear_coeffs`, `dsolve_analyze_roots`,
  `dsolve_variation_of_parameters` (used by both const-coeff and Euler). Still to
  do: reduction of order, normal form, special-function recognizer table.
- **M4 — systems (1e) + nonlinear higher-order (1d).** ✅ DONE. `LinearFirstOrderSystem`
  (constant-coefficient, eigen-based, real + complex eigenvalues, constant
  forcing — *diagonalizable only*; the defective / singular-`A` / general-forcing
  cases are lifted in **M8**) and `DecoupleSystem`; multi-function
  verify/fit/assemble added to the substrate; the `nfun>1` dispatch route added.
  Nonlinear higher-order (1d):
  `ReductionOfOrder` (2nd-order missing y) and `AutonomousReduction` (2nd-order
  missing x, two-stage recursion) both done; `EnergyIntegral` subsumed by the
  latter for elementary cases.
- **M5 — 2nd-order linear variable-coefficient: Kovacic + Frobenius.** ✅ DONE.
  The target is the rational-coefficient equation `y'' + P(x) y' + Q(x) y == 0`,
  which previously only `SpecialFunctionForm` skimmed (Airy / Bessel patterns).
  *Landed:* `DSolve`NormalForm` (substrate `dsolve_second_order_PQ` /
  `dsolve_normal_form`, reused by `SpecialFunctionForm`), `DSolve`Kovacic`
  (Cases 1 & 2; Case 3 declines), and `DSolve`FrobeniusSeries`/`DSolve`PowerSeries`
  (auto last-resort fallback). Kovacic uses a **Riccati/undetermined-coefficient**
  construction (transparent, back-substitution-verifiable) rather than the
  classical exponent tables; Case 1b adds the polynomial-`P` step for apparent
  singularities (`z''=(x²+3)z → x Exp[x²/2]`); Case 2 numerically verifies its
  algebraic candidates. The `dsolve_verify_body` substrate now substitutes
  derivative terms as `D[body,{x,k}]` directly so `SeriesData` bodies verify
  (the pure-function derivative of a SeriesData evaluates to 0). The later-M5
  items (`ExactODE`, `OperatorFactor`/`DFactor`, Sturm–Liouville
  `EigenvalueProblem`) remain future work. The two flagship methods slot into the
  scalar cascade **after** the recognizers and **before** a generic decline
  (closed form preferred over series):

  - **`NormalForm`** — shared prerequisite, land first. Substitute
    `y = z·Exp[-∫P/2]` to kill the `y'` term, giving the reduced form
    `z'' == r z` with `r = P²/4 + P'/2 − Q` (rational). Kovacic *and* the
    Airy/Bessel recognizers want this form; expose it as
    `dsolve_normal_form(P, Q, x) → r` in the substrate and route the recovered
    `z`-solution back through the `Exp[-∫P/2]` factor.

  - **`Kovacic`** — Liouvillian solutions of `z'' == r z`, `r ∈ C(x)`, by the
    classical three-case algorithm driven by the poles of `r` (located + ordered
    via `Apart` / `FactorList`) and its order at infinity. Staged by rising cost:
      1. **Case 1** — a solution `Exp[∫ω]` with `ω ∈ C(x)`: assemble the
         candidate `ω` from the per-pole and at-∞ local data, fix its polynomial
         part by a degree bound `d`, and solve the resulting linear system
         (`Solve`); the second solution follows by `ReductionOfOrder`. Covers the
         common `Exp[rational]` / `Exp[polynomial]` / rational-power Liouvillian
         solutions.
      2. **Case 2** — a solution algebraic of degree 2 over `C(x)` (`ω` a root of
         a quadratic over `C(x)`): the same pole-exponent + `d`-bound search over
         the sign choices, same solver.
      3. **Case 3** — algebraic of degree 4/6/12 (tetrahedral / octahedral /
         icosahedral): recognized but **gated/optional** — return the structural
         `Root`-form or decline with a message rather than run the heavy search by
         default.
    Local exponents are algebraic numbers → reuse `RootReduce` / qqbar.
    Elementary answers verify by back-substitution; the algebraic cases rest on
    recognizer structure + the stress corpus (as with the inert special-function
    heads).

  - **`FrobeniusSeries` / `PowerSeries`** — series solutions about `x0`, the
    always-available fallback when no closed form is found. Classify `x0`:
    **ordinary** → two `PowerSeries`; **regular singular** → Frobenius
    `y = (x−x0)^s Σ aₙ (x−x0)ⁿ` with indicial quadratic `s(s−1) + p₀ s + q₀ == 0`
    (`p₀ = lim (x−x0)P`, `q₀ = lim (x−x0)²Q`; solved by `Solve`, reusing
    `dsolve_analyze_roots` as `EulerCauchy` does); **irregular** → decline. The
    root-difference `s₁−s₂` picks the sub-case: non-integer → two independent
    series; equal → second solution carries a `Log`; positive integer → a
    coefficient-obstruction test decides whether a `Log` appears. Output is a
    truncated `SeriesData` to the requested order, verified by requiring the
    truncated residual to be `O[(x−x0)^N]`. Reuse `Series`/`SeriesData` arithmetic
    throughout.

  Later within M5: `ExactODE` (higher-order exact equations) and
  `OperatorFactor`/`DFactor` (factor a higher-order linear operator into
  lower-order factors and compose their solution sets) are both **done** — see §1c.
  Still future: the Sturm–Liouville `EigenvalueProblem` (1f). Reuse hooks across M5: `dsolve_linear_coeffs`
  (extract P, Q), `Apart`/`Together`/`FactorList` (pole structure),
  `Solve`/`solvepoly` (indicial + coefficient systems), `Series`/`SeriesData`
  (Frobenius), `RootReduce`/qqbar (algebraic exponents), `ReductionOfOrder`
  (second solution).
- **M6 — Phase 2 PDEs.** STARTED. First-order linear constant-coefficient PDE
  (method of characteristics) done — transport, `3u_x+5u_y==x`, `u_x+3u_y+u==1`.
  The `is_pde` dispatch route + PDE verify/assemble (2-variable `Function`) added.
  Quasilinear/nonlinear first-order and 2nd-order (wave/heat) still to do.
- **M7 — first-order substitution + attribute cleanup.** ✅ DONE.
  `DSolve`FirstOrderSubstitution` (`y'==F(a x + b y + c)`, completing the 1a
  first-order family bar Riccati/Lagrange/Abel/Chini) and `AutonomousReduction`
  (1d, above). `DSolve` and every `DSolve`<Method>` are now `Protected` only:
  `HoldAll` was dropped (equations survive evaluation as formal `Equal[...]`, so an
  equation held in a variable now solves), and `ReadProtected` was removed
  system-wide (Mathilda is fully open source).
- **M8 — general linear systems (defective, singular, arbitrary forcing) +
  triangular systems.** ✅ DONE. Two principled generalizations that together
  retire the M4 "diagonalizable only" restriction and the `DecoupleSystem`
  "one function per equation" restriction, closing the `{y'==0, x'+y==0}` class —
  a coupled constant-coefficient system whose matrix `A = {{0,0},{-1,0}}` is
  *defective* (one Jordan block, eigenvalue 0 doubled) **and** singular — in full
  generality rather than via a triangular special-case hack.
  - **Fundamental-matrix rework of `LinearFirstOrderSystem`.** Solve
    `Y' == A Y + b(x)` for **any** constant `A` (diagonalizable *or* defective)
    via `Φ(x) = e^{Ax}` built from the Jordan form: `{S,J} = JordanDecomposition[A]`;
    split `J = D + N` (diagonal eigenvalues `D` + strictly-upper nilpotent `N`);
    `e^{Jx} = e^{Dx} · e^{Nx}` with `e^{Dx} = DiagonalMatrix[Exp[λ_i x]]` and
    `e^{Nx} = Σ_{m<n} N^m x^m / m!` (a **finite** sum — `N` is nilpotent, so
    `N^n == 0`, verified via `MatrixPower`); then `Φ = S · e^{Jx} · S^{-1}`.
    Homogeneous solution `Y = Φ · {C[1],…,C[n]}`; forcing `b(x)` by variation of
    parameters `Y = Φ·(C + ∫ Φ^{-1} b dx)`, which **subsumes** the old `-A^{-1}b`
    particular and — crucially — stays valid when `A` is singular (exactly the
    failing example). `D` and `N` commute (each Jordan block is scalar on its
    diagonal), so the split is exact. Complex eigenvalues make `J`/`S` complex;
    the resulting complex-exponential body is reduced to real `e^{αx}Cos/Sin[βx]`
    form via `ComplexExpand`/`Simplify` (kept in complex form if it does not
    reduce, matching the "structurally exact" policy). The zero-eigenvector guard
    and the defective decline in `dsolve_linsys.c` are deleted; the diagonalizable
    case is now just `N == 0`. Symbolic `MatrixExp` is inert, so `Φ` is built from
    Jordan directly — factor the `Φ`-builder (`dsolve_fundamental_matrix`) so a
    future symbolic `MatrixExp[m]` can reuse it.
  - **`TriangularSystem`** (new method; generalizes `DecoupleSystem`). When the
    inter-function dependency graph is a DAG, topologically sort it, solve the
    functions in order, substitute each solved function forward into the
    still-unsolved equations, recurse into the scalar engine per function, and
    renumber constants (reusing `renumber`/`extract_body` from
    `dsolve_decouple.c`). This covers coupled-but-triangular systems at **any**
    coefficient — constant *or* variable (`{y'==x^2 y, x'==y}`) — the class that
    neither `DecoupleSystem` (needs zero edges) nor the constant-`A` matrix
    exponential (needs constant coefficients) can reach. Cascade order for
    `nfun>1`: `DecoupleSystem` → `TriangularSystem` → `LinearFirstOrderSystem`
    (cheapest / cleanest constants first; the matrix exponential is the general
    backstop for irreducibly-coupled constant systems).
  - **Higher-order linear systems** reduce to first order by state augmentation
    (introduce `y_i == u^(i)` auxiliary functions) in the substrate, then feed
    the reworked core — lands with, or immediately after, the above.

  Union of the two methods leaves exactly one honest gap: genuinely coupled,
  *non-triangular*, *variable*-coefficient systems (`LinearSystemVarCoeff`,
  Floquet/Magnus — still future).

  *Landed:* `dsolve_linsys.c` reworked around `mat_exp` (Jordan + finite
  nilpotent series) with the `Simplify[ComplexExpand[·]] //. Cosh[a]+Sinh[a]->E^a`
  realifier; new `dsolve_triangular.c`; shared `dsolve_renumber_constants` /
  `dsolve_extract_system_body` factored out of `dsolve_decouple.c`; cascade is
  `DecoupleSystem → TriangularSystem → LinearFirstOrderSystem`. **Constant-
  namespace hazard** (fixed): forward-substituting a solved function's `C[k]`
  into a later equation collides with the fresh `C[k]` the scalar engine emits
  for that equation (`GeneratedParameters` does not help — DSolve renames *all*
  `C[k]`, merging them); `TriangularSystem` parks solved constants in the private
  head `DSolve\`sysK` during the peel and remaps `sysK[k] -> C[k]` only at the
  end. Verified: `{y'==0, x'+y==0}` → `{y->C[1], x->C[2]-C[1] t}`; defective
  non-triangular, complex, forced-singular, and variable-coefficient triangular
  IVPs all back-substitute to zero (`tests/test_dsolve.c` t_sys_*).
- **M9 — SymPy parity gaps (deterministic).** ✅ DONE. `Factorable`, `NthAlgebraic`,
  `AlmostLinear`, `LinearCoefficients`, `SeparableReduced`, `Liouville`,
  `UndeterminedCoefficients`, `FirstOrderPowerSeries` — all eight implemented as
  self-contained `dsolve_<m>.c` (three-function contract), reducing to existing
  methods with one shared helper added (`dsolve_homog_basis`, moved from
  `dsolve_constcoeff.c` into the substrate). Cascade: `Factorable` + `NthAlgebraic`
  run EARLY (front, matching SymPy — a product/power-in-`y'` form is split before the
  specialists match it); `UndeterminedCoefficients` before `LinearConstantCoefficients`;
  the substitution reductions (`LinearCoefficients`/`AlmostLinear`/`SeparableReduced`)
  after the named first-order specialists; `Liouville` in the 2nd-order group;
  `FirstOrderPowerSeries` pinned-only (not auto — opt-in, matching SymPy/MMA). Unit +
  forward-generator stress families; all three DSolve ctest suites + `make check-c99`
  green; seven of eight per-call valgrind-flat (`AlmostLinear` inherits the
  pre-existing `Integrate`-engine per-call leak that `LinearFirstOrder` also has).
  Robustness notes: `Factorable` factors over plain-symbol substitutes with a
  `PolynomialQ` gate (raw funcapp `FactorList` hangs/misfactors) and keeps only
  differential factors; `UndeterminedCoefficients` uses `Expand` not `Simplify` before
  `Coefficient[·,Cos[b x]]`.
- **M10 — heuristic Lie point-symmetry (`lie_group`).** The one deliberately
  heuristic method in the cascade (first-order symmetry finding is an
  underdetermined PDE — no decision procedure exists). Nine ansatz heuristics from
  Cheb-Terrab et al.; a general first-order backstop inserted after `Abel` and
  before the implicit / series fallbacks. Substrate: the linearized symmetry
  residual `S(ξ,η) = η_x + (η_y − ξ_x)ω − ξ_y ω² − ξ ω_x − η ω_y` and a
  `checkinfsol` gate; integration by the Lie integrating factor `μ = 1/(η − ω ξ)`
  → first integral → `dsolve_run_implicit` (+ explicit `Solve` inversion), so every
  branch verifies with no inert heads. Staged **L1** ✅ (substrate +
  `abaco1_simple` end-to-end), **L2** (`linear` ✅ — affine ansatz →
  linear-coefficients class via determining-system `NullSpace`; `abaco1_product` ✅;
  `abaco2_similar` ✅ — §4.3 similarity ansatz `[F(x), H(x)]`, the first to reach
  **irrational** ODEs (`y' == Sqrt[a x + b y + c]`, `(a x + b y + c)^p`);
  `function_sum` ✅ — §4.2 additive ansatz `[F(x)+G(y), 0]`, classified by the
  rational factor `ω·∂²ₓ(1/ω) = F''/(F+G)`), **L3** (`bivariate` ✅ — general
  degree-2/3 bivariate-polynomial ansatz, the exact generalization of `linear`'s
  NullSpace determining system, catching genuinely quadratic/projective symmetries
  the affine ansatz misses; `abaco2_unique_unknown` ✅ — §4.4.1, the `[F(x),G(y)]`/
  `[G(y),F(x)]` ansätze from a function/non-integer-power of both variables in `ω`
  (`y' == (x/y)(x²+y²)^(1/3)`); `chi`, `abaco2_unique_general`
  remain). The quadrature ansätze are transcribed directly
  from the Cheb-Terrab & Roche (1998) invariant-family necessary conditions (the
  paper is in the repo root). `abaco1_product` (§4.1) uses the Eq-19 separability of
  `L = (ω_xy ω − ω_x ω_y)/ω⁴`; its shared substrate — `lie_sep_xfactor` (the
  product-separable x-factor via a fast rational free-of test), `lie_ratsimp`
  (`Cancel[Together]` in place of the far costlier `Simplify` on the hot path),
  `lie_inverse_omega`/`lie_swap_xy` (the inverse-ODE, so one extractor covers a
  pattern and its inverse) — is reused by the remaining quadrature heuristics. `abaco2_similar` (§4.3) finds `[ξ = F(x), η = H(x)]` from `Q = ω_y/ω_yy`,
  `T = Q_x/Q_y` (free of `y`), `F = Exp[∫((T ω_y − T_x − ω_x)/(ω+T)) dx]`, `H = −T F`
  — an irrational-`ω` reach the rational ansätze structurally lack. Lie
  sub-cascade is ordered **cheapest-first** (`abaco1_simple` → `linear` →
  `abaco1_product` → `abaco2_similar` → `bivariate`), so the degree-2/3 NullSpace runs
  only as a last resort. Reuses `D`, `Coefficient`/`Collect`, `Solve`, `Integrate`, `Together`,
  `dsolve_run_implicit`, `zero_test_decide`. `dsolve_lie.c`. References:
  Cheb-Terrab & Roche (CPC 113, 1998); Cheb-Terrab, Duarte & da Mota (CPC 101,
  1997); Cheb-Terrab & Kolokolnikov (math-ph/0007023); see
  `docs/design/dsolve_lie_symmetry.md`.

## Phase 1 — ODE method catalog

Cascade order: cheap deterministic recognizers first. `[✓]` implemented,
`[ ]` planned.

### 1a. First order
- `[✓] Quadrature` — `y^(n)==f(x)`, `f` free of `y`: integrate `n` times + constant polynomial.
- `[✓] LinearFirstOrder` — `y'+p(x)y==q(x)`: integrating factor `Exp[∫p]`.
- `[✓] Separable` — `y'==g(x)h(y)`: `∫dy/h==∫g dx + C[1]`, solved for `y`.
- `[✓] Bernoulli` — `y'==A y + B y^n` (n≠0,1): substitution `v=y^(1-n)` (exponent recovered robustly, incl. fractional/negative n).
- `[✓] Homogeneous` — `y'==F(y/x)`: substitution `y=v x` → separable. Direct
  log-form inversion, with an exponentiate-and-clear-radicals fallback
  (`homog_exp_log_invert`: `Prod g_i^{c_i} == C[1] x` raised to power `d` → `Solve`
  Root branches) for the pure-log (real-root) rational family, e.g.
  `(x+2y)/(2x+y)`. The transcendental (ArcTan log-spiral) subset has no explicit
  inverse and is returned as the **implicit first integral** `G(x,y[x]) == C[1]`
  via the `dsolve_run_implicit` path (`dsolve_homogeneous_implicit_try`), verified
  by implicit differentiation (`y' == -G_x/G_y` satisfies the ODE).
- `[✓] Exact` — `M+N y'==0`, `M_y==N_x`; + integrating-factor search `μ(x)`, `μ(y)`.
- `[✓] Clairaut` — `y==x y'+f(y')`: general line `y=C[1]x+f(C[1])` + singular envelope (`IncludeSingularSolutions`).
- `[✓] Riccati` — `y'==q0+q1 y+q2 y^2` (`q2≠0`): linearise `y=-u'/(q2 u)` →
  2nd-order linear `u'' - (q1+q2'/q2) u' + q0 q2 u == 0`, solved by recursing into
  the scalar cascade (const-coeff / Euler / Airy-Bessel / Kovacic / Frobenius);
  the redundant constant is collapsed (`C[2]->1`) to the single Riccati parameter.
  Runs after `FirstOrderSubstitution` so `fos` keeps its cleaner `-x-Tan[…]` for
  `y'==(a x+b y+c)^2`. Solves the Airy-linearised `y'==y^2+x` that previously
  declined; declines (no wrong answer) when the linearisation has no closed form.
- `[✓] Lagrange` (d'Alembert) — `y==x φ(y')+ψ(y')` (`φ(y')≠y'`): the general
  solution is **parametric** `{x=X(t,C), y=X φ+ψ}`, `t=y'`, where `X(t)` solves the
  linear ODE `dx/dt − [φ'/(t−φ)]x = ψ'/(t−φ)` (via `dsolve_linear_factor_solve`).
  New parametric substrate path (`dsolve_run_parametric` / `_verify_parametric` /
  `_assemble_parametric` / `_method_builtin_parametric`, mirroring the implicit
  path); output `{{x->Function[{t},X], y->Function[{t},Y]}}`, verified by
  `y'=Y'(t)/X'(t)`. Declines Clairaut (`φ≡p`) and genuinely-linear equations. Runs
  after `Clairaut` in the cascade. *Deferred (future):* singular-line solutions
  (roots of `φ(p)=p`) and parametric IVP constant-fitting (an IVP declines).
- `[✓] Chini` — `y'==f(x) y^n+g(x) y+h(x)` (n≠0,1,2): the reducible-to-autonomous
  sub-class, via `y=f^(-1/(n-1)) u` → `u'==u^n+B u+C` (B,C constant); implicit first
  integral `∫du/(u^n+Bu+C)−x==C[1]` (rational integrand, always elementary) returned
  through `dsolve_run_implicit`. Non-reducible cases decline. `dsolve_chini.c`
  (shared `dsolve_chini_first_integral`).
- `[✓] Abel` — `y'==f3 y^3+f2 y^2+f1 y+f0` (f3,f2≠0): remove the y² term
  (`z=y+f2/(3 f3)`) → Chini n=3 → same implicit first integral. `dsolve_abel.c`
  (thin front-end over the shared Chini helper). The fuller constant-invariant
  class (with an x-rescaling) is future.
- `[✓] FirstOrderSubstitution` — `y'==F(a x + b y + c)`: detect the constant ratio
  `r = F_x/F_y`, substitute `v = y + r x` → autonomous separable `v'==r+H(v)`,
  solved inline; declines (stays symbolic) when the antiderivative does not invert.
- `[✓] Factorable` — factor the equation as a polynomial in the highest derivative
  (`F1·F2·…==0`); recurse the cascade on each factor, union the branch bodies.
  (SymPy `factorable`.) Factors over plain-symbol substitutes under a `PolynomialQ`
  gate (raw funcapp `FactorList` hangs/misfactors); keeps only differential factors.
  Runs at the front. `dsolve_factorable.c`.
- `[✓] NthAlgebraic` — algebraic (degree ≥ 2) in the top derivative `y^(n)`: `Solve`
  for `y^(n)`, recurse each root branch (a branch free of `y` hits `Quadrature`);
  also the degenerate no-derivative case. Runs at the front. (SymPy `nth_algebraic`.)
  `dsolve_nth_algebraic.c`.
- `[✓] AlmostLinear` — `f(x)g(y)y' + k(x)l(y) + m(x)==0`: substitution `u=∫g dy` →
  `u'+P u==Q` (integrating factor), then `l(y)==U(x)` solved for `y`. (SymPy
  `almost_linear`.) `dsolve_almostlinear.c`.
- `[✓] LinearCoefficients` — `y'==(a1 x+b1 y+c1)/(a2 x+b2 y+c2)`: det≠0 → shift to the
  lines' intersection → `Homogeneous` (explicit / implicit); det=0 (parallel) →
  substitute `v=a1 x+b1 y` → `Separable`. Explicit + implicit two entries. (SymPy
  `linear_coefficients`.) `dsolve_lincoeff.c`.
- `[✓] SeparableReduced` — `x y'/y == G(x^n y)`: `n = x r_x/(y r_y)`, substitution
  `w=x^n y` → `Separable`; implicit first integral. (SymPy `separable_reduced`.)
  `dsolve_sepreduced.c`.
- `[~] LieSymmetry` (`DSolve`LieGroup`/`LieSymmetry`) — heuristic infinitesimal
  point-symmetry method; the general first-order backstop, run after the specialists
  and before the series fallback. *Implemented:* `abaco1_simple` (L1) + `linear`
  (L2, affine ansatz → linear-coefficients class via determining-system `NullSpace`)
  + `bivariate` (L3, general degree-2/3 polynomial ansatz — same NullSpace machinery
  at higher degree, `lie_poly_symmetry`; catches quadratic/projective symmetries the
  affine case misses, e.g. ξ=x², η=xy for ω = y/x + A(y/x)/x)
  + `abaco1_product` (§4.1, the symmetry [F(x)G(y), 0] and its inverse [0, F(x)G(y)]
  — a rational-but-non-polynomial infinitesimal, e.g. ξ=y/x, that the polynomial
  ansätze miss; found by the Cheb-Terrab & Roche Eq-19 product-separability of
  L = (ω_xy ω − ω_x ω_y)/ω⁴, then Eq-20 for G(y); the inverse pattern via the
  inverse ODE 1/ω(y,x)). SymPy's `lie_group` times out (>25 s) on the pure product
  family 2xy/(x²+2y⁴+c), which this solves in ~30 ms.
  + `abaco2_similar` (§4.3, the symmetry [F(x), H(x)] and its inverse [F(y), H(y)] —
  both components single-variable functions; from Q = ω_y/ω_yy, T = Q_x/Q_y free of y,
  F = Exp[∫((T ω_y − T_x − ω_x)/(ω+T)) dx], H = −T F. First heuristic to reach
  irrational ω: solves y' = Sqrt[a x + b y + c] and (a x + b y + c)^p, p non-integer).
  + `function_sum` (§4.2, the additive symmetry [F(x)+G(y), 0] and its inverse
  [0, F(x)+G(y)]; classified by the rational factor ω·∂²ₓ(1/ω) = F''/(F+G) whose
  reciprocal's ∂_y separates by product, x-factor 1/F'').
  + `abaco2_unique_unknown` (§4.4.1, the symmetries [F(x),G(y)] and [G(y),F(x)] from a
  function/non-integer-power M of both variables in ω: R = M_y/M_x separates by product
  with x-factor X, candidates [X,−X/R] and [−R/X,1/X]; solves y' = (x/y)(x²+y²)^(1/3)).
  Remaining: `chi`, `abaco2_unique_general`.
  Nine ansatz heuristics (`abaco1_simple`,
  `abaco1_product`, `function_sum`, `abaco2_similar`, `linear`, `bivariate`, `chi`,
  `abaco2_unique_unknown`, `abaco2_unique_general`); each candidate `(ξ,η)` is gated
  through the symmetry condition, then integrated by the Lie integrating factor
  `μ = 1/(η − ω ξ)` → first integral (reuses the exact-equation quadrature) →
  `dsolve_run_implicit` + explicit `Solve` inversion. (SymPy `lie_group`.) See M10.
  `dsolve_lie.c`.

### 1b. Linear constant-coefficient
- `[✓] LinearConstantCoefficients` — one method for homogeneous + inhomogeneous.
  Characteristic polynomial `Σ a_k λ^k`; roots via `Solve` with derivative-based
  multiplicity + dedup; complex-conjugate pairs → `e^(ax)(Cos,Sin)`; repeated
  roots → `x^k e^(rx)`. Inhomogeneous by variation of parameters (Wronskian /
  Cramer + `Integrate`, particular `Simplify`d).
  *Cosmetic gap:* for simple forcing the var-params particular can carry a
  homogeneous component (`7/2 Cos^2 x` for `7/4`) — correct and verified, less
  tidy than undetermined coefficients (a future refinement).
- `[✓] UndeterminedCoefficients` — tidy particular for polynomial / exponential /
  sinusoid forcing (incl. resonance) for **constant-coefficient** ODEs; runs before
  `LinearConstantCoefficients` (which stays the var-params fallback for other
  forcing). Superposition over `Expand[g]`; resonance shift `s` found by
  incrementing until the coefficient system solves. New file `dsolve_undetcoeff.c`
  (reuses the shared `dsolve_homog_basis`). Euler variant is future. (SymPy
  `nth_linear_constant_coeff_undetermined_coefficients`.)

### 1c. Linear variable-coefficient
- `[✓] EulerCauchy` — `a_n x^n y^(n)+…+a_0 y == g`: indicial polynomial
  `Σ a_k (r)_k` (falling factorial) via trial `x^r`; real roots → `x^r (Log x)^j`,
  complex pairs → `x^a Cos/Sin[b Log x]`, repeated → `Log x` powers; forcing via
  variation of parameters (real-root forcing works; complex-root forcing with a
  hard Wronskian integral declines gracefully).
- `[→M5] ReductionOfOrder` — second solution from one known (also the second-
  solution engine reused by Kovacic Case 1).
- `[✓] ExactODE` — higher-order exact linear equations: `L[y] == d/dx(M[y])`
  (exactness `Σ(-1)^k a_k^(k) == 0`, tested as `a_0 == b_0'` via the first-integral
  recurrence `b_{n-1}=a_n`, `b_{k-1}=a_k-b_k'`). Integrate once to the first
  integral `M[y] == ∫g + C[n]` and recurse into the scalar cascade on the
  order-(n-1) equation (constant `C[n]` contiguous with the sub-solve's
  `C[1..n-1]`, no renumbering; iterated exactness free via the recursion). Runs
  after `EulerCauchy`, before `SpecialFunctionForm`. First cut linear/order≥2/
  genuinely exact; integrating-factor (adjoint) exactness and nonlinear
  total-derivative detection are future. `dsolve_exactode.c`.
- `[✓] NormalForm` — reduce `y''+P y'+Q y` to `z''==r z` via `y=z Exp[-∫P/2]`,
  `r=P²/4+P'/2−Q`; prerequisite for Kovacic and the special-function recognizers
  (`dsolve_normal_form` substrate helper).
- `[~] SpecialFunctionForm` — matches the normalised 2nd-order form
  `y''+P y'+Q y==0`: **Airy** (`P=0`, `Q` linear → AiryAi/AiryBi, verifies),
  **Bessel / modified Bessel** (`P=1/x`, `Q=±1−ν²/x²` → BesselJ/Y, BesselI/K —
  correct heads; residual reduces only via Bessel recurrences zero_test can't
  decide, so kept as structurally exact), **Kummer** confluent hypergeometric
  (`x y''+(b−x)y'−a y==0` → `Hypergeometric1F1[a,b,x]` + `x^(1−b) 1F1[a−b+1,2−b,x]`)
  and **Gauss** (`x(1−x)y''+(c−(a+b+1)x)y'−ab y==0` → `Hypergeometric2F1[a,b,c,x]` +
  `x^(1−c) 2F1[a−c+1,b−c+1,2−c,x]`; `a,b` recovered from `a+b`, `ab` via a
  quadratic whose linear factors give radical-free roots). The hypergeometric
  heads auto-rewrite to `HypergeometricPFQ`, which has a `deriv.c` z-derivative
  rule, so the branches verify. The second solution carries `x^(1−b)`/`x^(1−c)`
  and is emitted only when that exponent parameter is a **non-integer number**
  (an integer makes the pair dependent / the pFq lower parameter singular; a
  symbolic exponent makes the verify residual a symbolic-power+pFq sum on which
  zero_test currently hangs). Both degenerate cases decline to the Frobenius
  series fallback; the other parameters (`a`; `a,b`) may stay symbolic.
  TODO: canonical-singular-point restriction lifted by an affine/Möbius change of
  variable. NOTE: the integer-degree **Legendre / Chebyshev / Gegenbauer / Jacobi**
  family (both solutions elementary) is now solved *algorithmically* by Kovacic
  Case 1 (below), not by a recognizer — no LegendreP/Q head needed. A recognizer is
  only wanted for the **non-integer** degree (genuinely hypergeometric) cases and
  for inert heads: LegendreQ, HermiteH, Laguerre, Whittaker, Mathieu, Spheroidal,
  Kelvin, ParabolicCylinder, Struve, Weierstrass.
- `[✓] Kovacic` — Liouvillian solutions of the reduced form `z''==r z`
  (`r∈C(x)`); three-case algorithm on the poles of `r` (`Apart`/`FactorList`) +
  order at ∞. Staged: Case 1 (`z==P Exp[∫θ]`, `θ,P` rational) → Case 2 (degree-2
  algebraic) → Case 3 (degree 4/6/12, gated). **Case 1 apparent-singularity
  completion (`kovacic_case1_general`):** the pole-only Riccati ansatz misses a
  solution whose `z1` has zeros off the poles of `r`; the classical fix builds
  `θ = Σ_c α_c^{s}/(x−c)` from the local pole exponents `α_c=(1±√(1+4 b_c))/2`
  (`b_c = lim(x−c)²r`; poles located incl. complex via `dsolve_analyze_roots`) and
  a monic `P` of degree `d = α_∞ − Σα_c` (classical degree bound, tested
  *numerically* so complex-α combinations reject instantly), then `z1 = P Exp[∫θ]`.
  The second solution is taken at the **y-level** (`y2 = y1 ∫ w²/y1²`, reduction of
  order on the original ODE) so the recovery radical cancels up front — the z-level
  `z1∫1/z1²` instead sends Simplify into a minutes-long blow-up. Solves integer-degree
  Legendre/Chebyshev/Gegenbauer/Jacobi in elementary form. Case 2 is **skipped when a
  denominator factor has degree ≥ 2** (complex poles): its σ-`Solve` blows up there and
  Case 1c already covers those poles. Algebraic exponents via `RootReduce`/qqbar;
  second solution via `ReductionOfOrder`.
- `[✓] FrobeniusSeries` / `[✓] PowerSeries` — series about `x0`: ordinary →
  power series, regular-singular → Frobenius (indicial quadratic + `Log`-term
  sub-cases by root difference), irregular → decline; truncated `SeriesData`
  verified as `O[(x−x0)^N]`. Reuse `Series`/`SeriesData`.
- `[✓] FirstOrderPowerSeries` — order-1 power series about the ordinary point x0=0
  (`dsolve_frobenius.c`, `dsolve_first_order_series_try`; a₀=C[1], one Taylor read
  per order). Pinned-only (not auto — opt-in, matching SymPy/MMA). Solves nonlinear
  `y'==x+y²` etc. (SymPy `1st_power_series`.)
- `[✓] OperatorFactor` (`DSolve`DFactor`) — factor a homogeneous linear operator
  (order ≥ 3) by finding a first-order right factor `(D − r)`, `r ∈ C(x)` (a
  hyperexponential solution `Exp[∫r]`, via a rational Riccati `Σ a_k P_k(r) == 0`
  undetermined-coefficient search); peel via operator right-division, recurse `DSolve`
  on the order-(n−1) quotient, close with the trailing first-order solve. Reaches
  reducible variable-coefficient operators the earlier methods miss (shifted-Euler at
  a pole ≠ 0). `DSolve`DFactor[eqn,y,x]` returns `{Dx − r1, Dx − r2, …}`. Runs after
  Kovacic (order 2), before the reduction/series methods. First cut: first-order
  **right** factors, homogeneous, rational coefficients; irregular-singular / 2nd-order
  right factors (Beke) are future. `dsolve_operator_factor.c` (self-contained; no
  changes to the Kovacic engine).

### 1d. Nonlinear higher-order
- `[✓] ReductionOfOrder` — `y''==F(x,y')` missing y: reduce to first order in
  p=y' (recurse into the scalar engine), then `y=∫p dx + C[2]`. Guards against a
  wrong `Integrate` antiderivative (requires `D[∫p]==p`, decided by `zero_test`
  then `PossibleZeroQ` sampling so a correct-but-unsimplified antiderivative — a
  multi-`Log` `∫Tan`, an `ArcTan[x/Sqrt[C]]` — is accepted while a degenerate
  `y=const` is still rejected) so it declines instead of shipping a degenerate
  solution. Solves `y''==(y')^2`, the autonomous `a+b(y')^2` (→ Tan/Tanh), and the
  Riccati-in-p `c x (y')^2` families.
- `[✓] AutonomousReduction` — `y''==f(y,y')` missing `x`: `p=y'(y)`, `p p'(y)==f`
  (recurse), then `y'==p(y)` separable (recurse), constants renumbered across the
  two stages; final body required to depend on `x` (rejects the degenerate
  `y=const` that trivially back-substitutes). Solves `y y''==(y')^2 → C[2] E^(C[1] x)`.
- `[~] EnergyIntegral` — `y''==f(y)`: subsumed by AutonomousReduction for the
  elementary cases (`f` free of `y'` is a special case); genuinely elliptic ones
  (`y''==2y^3`, `y''==-Sin[y]`) still decline (`WeierstrassP` inert head: future).
- `[✓] Liouville` — `y'' + g(y)(y')^2 + h(x)y' == 0`: Liouville's transformation →
  two quadratures `∫Exp[∫g dy] dy == C[1] ∫Exp[-∫h dx] dx + C[2]`, solved for `y`.
  Distinct from `AutonomousReduction` (missing-`x`) / `ReductionOfOrder` (missing-`y`).
  After `AutonomousReduction` in the cascade. (SymPy `Liouville`.) `dsolve_liouville.c`.

### 1e. Systems

Cascade order (`nfun>1`): `DecoupleSystem` → `TriangularSystem` →
`LinearFirstOrderSystem`.

- `[✓] DecoupleSystem` — each equation involves one function: recurse into the scalar
  engine per function, renumber the generated constants. Handles variable-coefficient
  components (`y' == x^2 y`).
- `[✓] TriangularSystem` — inter-function dependency graph is a DAG: topologically
  sort, solve in order substituting each solved function forward, recurse into the
  scalar engine + renumber constants (reuses `dsolve_renumber_constants` /
  `dsolve_extract_system_body`, with solved constants parked in `DSolve\`sysK` to
  avoid colliding with the scalar engine's fresh `C[k]`). Coupled-but-triangular at
  **any** coefficient — constant or variable
  (`{y'==0, x'+y==0}` → `y=C[1], x=C[2]-C[1]x`;  `{y'==y/x, z'==y}`).
- `[✓] LinearFirstOrderSystem` — `Y' == A Y + b(x)`, constant `A`, **any** spectrum.
  Fundamental matrix `Φ = e^{Ax} = S · e^{Jx} · S^{-1}` from `JordanDecomposition`
  (diagonalizable → `C e^{λx} v`; **defective → `x^k e^{λx}`** generalized-eigenvector
  terms; complex pairs → real `e^{αx}Cos/Sin[βx]` via `ComplexExpand`). Forcing `b(x)`
  by variation of parameters `Φ·(C + ∫ Φ^{-1} b dx)` — subsumes `-A^{-1}b` and stays
  valid for singular `A`. Multi-function verify + IVP/BVP fit in the substrate
  (`dsolve_verify_system`/`dsolve_fit_system`/`dsolve_assemble_system`).
  *Was (M4): eigen-only, diagonalizable, `-A^{-1}b` forcing, defective decline.*
- `[ ] LinearSystemVarCoeff` — genuinely coupled, *non-triangular*, *variable* `A`
  (Floquet/Magnus); future. This is SymPy's non-constant-coefficient
  `linear_neq_order1` (types with a commutative antiderivative of the coefficient
  matrix) and is the systems parity target.

### 1f. Conditions
- `[✓]` IVP (fit constants at one point) — in the substrate.
- `[ ]` Linear BVP (multiple points) — in the substrate + `Solve`.
- `[ ] EigenvalueProblem` — Sturm–Liouville: `Piecewise` with parameter conditions.

## Phase 2 — PDE method catalog

### 2a. First order
- `[✓] PDELinearFirstOrder` — constant-coefficient `a u_{v1}+b u_{v2}+c u==f`
  (a,b constant; c,f functions) by the method of characteristics: invariant
  `ξ = a v2 - b v1`, then the linear ODE `u_{v1}+(c/a)u=f/a` along the
  characteristic gives `u = Exp[-∫c/a]( C[1][ξ] + ∫Exp[·]f/a )`. Solves the
  transport equation, `3u_x+5u_y==x`, `u_x+3u_y+u==1`. PDE verify substitutes a
  concrete test function (`C[1][z_]:>Sin[z]`) — zero_test cannot sample an
  arbitrary function, and D-of-a-2-var-Function-with-arbitrary-function crashes
  the evaluator (both pre-existing).
- `[ ] PDEQuasilinear` (Lagrange), `[ ] PDECharpit` (nonlinear complete
  integral), `[ ] PDEClairaut`.

### 2b. Second order
- `[ ] PDEHyperbolicGeneral` (operator factoring → arbitrary functions),
  `[ ] WaveEquation` (d'Alembert; inhomogeneous; half-line; `Piecewise`),
  `[ ] HeatEquation` (heat kernel / `Erf`), `[ ] SeparationOfVariables`,
  `[ ] PDEClassify` (discriminant).

## Testing

- **In-engine self-verification:** every returned branch back-substitutes to a
  residual that is not decidably non-zero (`dsolve_run` → `zero_test_decide`).
- **Unit tests** (`tests/test_dsolve.c`): the Wolfram reference examples, checked
  by `PossibleZeroQ` of the residual / of the difference from the known form.
  Every scalar method with a backtick builtin has a pinned-method test
  (`DSolve`<Method>[...]`); systems/PDE (no backtick builtin) are covered through
  the automatic dispatch.
- **Stress tests:** parametrized forward-generator families per method group,
  back-substitution verified, each guarding `Head === List` first so a declined
  solve cannot pass vacuously. `tests/test_dsolve_m5_stress.c` covers Kovacic and
  Frobenius/PowerSeries; `tests/test_dsolve_stress.c` covers the rest of the
  cascade (LinearFirstOrder, Separable, Bernoulli, Homogeneous, Exact,
  LinearConstantCoefficients, EulerCauchy, ReductionOfOrder, 2×2 + triangular
  systems, first-order PDE). A generator builds the equation from parameters whose
  closed form is guaranteed — a chosen spectrum (ConstCoeff/Euler), a potential
  (Exact), or a target solution `yt` with `q = yt' + p·yt` (LinearFirstOrder, so
  the integrating-factor integral is elementary by construction).
- **Gates:** `dsolve_tests` (ctest), `valgrind --leak-check=full`, `make
  check-c99`, and REPL spot-checks (`DSolve[...]`, `?DSolve`, `?DSolve`Separable`).

DSolve is a symbolic/structural head (no element-wise numeric mapping); it is a
documented exemption from the packed-aware / Compile surfaces.

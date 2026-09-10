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
- **M6 — Phase 2 PDEs.** IN PROGRESS. First-order linear constant-coefficient PDE
  (method of characteristics) done — transport, `3u_x+5u_y==x`, `u_x+3u_y+u==1`.
  The `is_pde` dispatch route + PDE verify/assemble (2-variable `Function`) added.
  **`PDELinearSecondOrder`** (`dsolve_pde2.c`) done — the homogeneous, principal-
  part-only, constant-coefficient 2nd-order linear PDE `A u_{v1 v1}+B u_{v1 v2}+
  C u_{v2 v2}==0` via operator factoring (trial `f(v2+λ v1)`, characteristic
  quadratic `A λ²+B λ+C==0`). One method covers all three discriminant signs —
  hyperbolic/distinct roots (the **wave equation** `u_tt==c² u_xx →
  C[1][x-c t]+C[2][x+c t]`, d'Alembert), elliptic/complex roots (Laplace
  `u_xx+u_yy==0 → C[1][y-I x]+C[2][y+I x]`, matching Mathematica's complex-
  characteristic form), parabolic/repeated root (`C[1][w]+v1 C[2][w]`) — reusing
  `dsolve_analyze_roots` (distinct/complex/repeated λ uniformly). This realizes
  the §2b `PDEHyperbolicGeneral` item in full generality, so no separate elliptic/
  parabolic method is needed. The shared `dsolve_verify_pde` was generalized to
  arbitrary derivative order (scanned from the residual — `max_order` is 0 for
  PDEs) and multiple arbitrary functions (`C[1..4]` → distinct concrete test
  functions). **`SeparationOfVariables`** (`dsolve_pdesep.c`, pinned-only) done —
  the separated product mode `u == X(v1) Y(v2)` for a homogeneous, constant-
  coefficient, no-mixed-term linear PDE (heat, Helmholtz), splitting into two
  constant-coefficient ODEs in a separation constant λ and recursing into the
  scalar cascade; back-substitution verified. **`PDEClassify`**
  (`dsolve_pdeclassify.c`) done — the discriminant classifier (Hyperbolic /
  Parabolic / Elliptic) over the principal part. **`WaveDAlembert`**
  (`dsolve_wave.c`) done — the d'Alembert initial-value problem for the 1-D wave
  equation on the whole line (`u_tt==c² u_xx` with `u(x,t0)==f(x)`,
  `u_t(x,t0)==g(x)`), auto-dispatched and pinned; the method sorts the PDE from its
  two initial conditions (which arrive as ordinary equations, `neq==3`) and carries
  its own multi-equation verify. **`HeatKernel`** (`dsolve_heat.c`) done — the
  Cauchy problem for the 1-D heat equation by the heat-kernel convolution (auto +
  pinned; `neq==2`, verified at the kernel level). **`PDEQuasilinear`**
  (`dsolve_pdequasi.c`) done — Lagrange's method of characteristics for the
  quasilinear `P u_{v1}+Q u_{v2}==R` (semilinear → explicit `C[1][ξ]`; conservation
  law R==0 → implicit `φ1==C[1][u]`, e.g. inviscid Burgers), via a new
  `dsolve_run_pde_implicit` substrate (implicit relation verified by the
  implicit-function rule). **`PDEClairaut`** (`dsolve_pdeclairaut.c`) done — the
  Clairaut form `u==v1 u_{v1}+v2 u_{v2}+f(u_{v1},u_{v2})` → complete integral +
  singular envelope. **`PDECharpit`** (`dsolve_pdecharpit.c`) done — first-order
  fully nonlinear `F(v1,v2,u,p,q)==0` by Charpit's method (the three standard forms
  `F(p,q)` / `F(u,p,q)` / separable `f(v1,p)==g(v2,q)`), via the new `PDERelation`
  implicit-constant verify path. **`PDELinearSecondOrder` lower-order terms** done —
  constant-coefficient 2nd-order with `D u_{v1}+E u_{v2}+F u` when the symbol factors
  into two first-order operators → exponential-damped `Σ e^{−m_i v1} C[i][ξ_i]`
  (distortionless telegraph, damped/convection). Still to do: inhomogeneous
  2nd-order forcing; the non-factorable general telegraph (Bessel); Charpit's general
  integrable-combination + singular solutions; genuinely-quasilinear non-conservation
  (P/Q depends on u with R≠0); inhomogeneous / half-line / `Piecewise` wave data; the
  Erf-producing heat data; finite-interval Fourier-series problems.
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
  (`y' == (x/y)(x²+y²)^(1/3)`); `abaco2_unique_unknown` also carries the §4.4.1
  "differential invariant of order zero" extension (Eqs 73–81: the non-separable
  candidates `[-R,1]`/`[1,-R]`/`[1,-1/R]`, catching Kamke 433 `(x y'+y+2x)²=4(xy+x²+a)`
  → `x−Sqrt[x²+xy+a]==C[1]`). `chi` ✅ (CPC 101 1997, 5th algorithm: the `η=ξω+χ`
  reformulation with a rich transcendental-atom basis for `χ`, solving Kamke 357
  `x y' ln(x)sin(y)+cos(y)(1−x cos(y))=0` → `−x+Log[x]Sec[y[x]]==C[1]`, the first
  genuinely-transcendental `χ` beyond `bivariate`). The formal §4.4.2 Case I/II
  (`abaco2_unique_general`) is a **documented exemption** — the authors call it "very
  inefficient, if not just impractical" and every `[F(x),G(y)]`-symmetric target is
  already caught by Riccati/Bernoulli/kernel methods, so it cannot be tested
  non-vacuously). **Robustness:** an `ω` carrying an undefined function of both
  variables (`Tan[ArcTan[y]+F[x²+y²]]`) used to hang the quadrature heuristics; a node
  budget on the hot-path helpers, a structural polynomial zero-test in place of the
  hanging general `zero_test`, and an undefined-function gate now make it a ~1 s clean
  decline. The quadrature ansätze are transcribed directly
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

- **M11 — finish Phase-1 ODE gaps.** ✅ DONE. Three principled additions closing the
  remaining §1e/§1f gaps, each back-substitution verified with unit + forward-generator
  stress tests.
  - **Refactor:** the constant-coefficient system core is factored out of
    `dsolve_linsys.c` into `dsolve_linsys.h` — `dsolve_linsys_matexp` (`e^{Mt}` Jordan
    builder), `dsolve_linsys_tidy` (realifier, now applying `ComplexExpand` only when
    the body carries the imaginary unit, avoiding a gratuitous `Log[x]→Log[Abs]+I Arg`
    split), `dsolve_linsys_extract_Ab` (constant-`A` guard moved to callers), and
    `dsolve_linsys_assemble(M,t,xvar,b,b_zero,n)`. `LinearFirstOrderSystem` output
    unchanged.
  - **`LinearSystemVarCoeff`** (§1e, `dsolve_linsys_varcoeff.c`) — coupled variable-
    coefficient systems of the scalar-factor class `A(x)=f(x)B`, via `t=∫f dx` and the
    reused assembler. Cascade after `LinearFirstOrderSystem`.
  - **Formal Linear BVP soundness** (§1f, `dsolve_common.c`) — an inconsistent /
    over-determined BVP returns `{}` (no solution) instead of the silent unfitted
    general solution: the constant-fitters gained a `no_solution` out-param (set on an
    empty-`Solve` result) threaded through `dsolve_run` / `dsolve_run_system` as the
    concrete empty list.
  - **`EigenvalueProblem`** (§1f, `dsolve_eigenvalue.c`) — pinned-only Sturm-Liouville
    `y''+λy==0` on `[a,b]` with two homogeneous BCs, eigenvalue family + eigenfunctions
    verified under the integer index.
  All three DSolve ctest suites + `make check-c99` green; valgrind leak-clean; SymPy
  cross-validated the coupled variable-coefficient systems.

- **M12 — nonlinear second-order Lie point symmetry (`SecondOrderSymmetry`).** ✅ DONE.
  The general **nonlinear**-second-order backstop, targeting the ~110 nonlinear
  2nd-order ODEs in the 12000.org "SymPy-failed" corpus (Kamke/Murphy; Maple's
  `_with_linear_symmetries` on nonlinear equations + `_reducible,_mu_*`). Extends the
  M10 first-order Lie machinery to the second prolongation. `dsolve_lie2.c`.
  - **Symmetry search.** A point symmetry `X = ξ(x,y)∂ₓ + η(x,y)∂_y` satisfies the
    second-prolongation determining equation (Cheb-Terrab, Duarte & da Mota,
    physics/9703082 — the PDF is in the repo root; Eq. 2, here in the standard
    prolongation form with `p=y'`):
    `η_xx + (2η_xy−ξ_xx)p + (η_yy−2ξ_xy)p² − ξ_yy p³ + (η_y−2ξ_x−3ξ_y p)Φ − ξΦ_x −
    ηΦ_y − (η_x+(η_y−ξ_x)p−ξ_y p²)Φ_p == 0`, `Φ = y''`. With `ξ,η` general bivariate
    polynomials (degree ≤ 2) and `Φ` rational in `(x,y,p)`, this clears to a
    polynomial identity in `(x,y,p)` whose coefficients are linear/homogeneous in the
    ansatz coefficients; the determining system's `NullSpace` is a basis of admissible
    `(ξ,η)` — the same machinery `dsolve_lie.c::lie_poly_symmetry` uses at first order,
    lifted to `{x,y,p}`. (Maple/SymPy `symgen way=3`.)
  - **Order reduction.** From one symmetry, build canonical coordinates `(r,s)` with
    `Xr=0`, `Xs=1` (zero-component symmetries → `r,s` explicit; both-nonzero → `r`
    from the characteristic `y'=η/ξ` via the cascade, `s=∫dx/ξ`). The ODE becomes a
    first-order ODE `dq/dr==F(r,q)` in `q=ds/dr`, solved by recursing into the scalar
    cascade; then `s=∫q dr + C[2]` and `s(x,y)==S(r(x,y))+C[2]` is solved for `y`.
  - **Two verification gates (both essential).** (1) Every inversion branch is tried
    and **numerically** back-substituted (`l2_num_ok`, 5 sample points) — the dsolve_run
    symbolic verify KEEPS an undecidable residual (Solve policy), so a wrong reduction
    branch with a symbolically-intractable residual (the seeds carry logs/radicals)
    would otherwise slip through as a WRONG answer; N1/N2 were exactly such false
    positives before this gate. (2) Only nonlinear ODEs: a **linearity gate** declines
    `Φ` affine in `(y,y')` (linear ODEs are the domain of Euler/Kovacic/
    SpecialFunction/Frobenius) — this also avoids the 8-dimensional-algebra symmetry
    search on a Kovacic-declined high-degree rational coefficient.
  - **Robustness.** lie2 chains recursive `DSolve`/`Solve`/`Integrate`/`NullSpace`
    whose cost is data-dependent and can hit a **pre-existing** slow path (e.g. the
    explicit inversion inside `DSolve\`Separable` hangs on a plain separable reduced
    ODE). Every such sub-call is `TimeConstrained`-bounded (2–3 s), the whole attempt
    carries a ~6 s wall-clock deadline, an undefined-function `Φ` is gated, and a
    per-top-level **decline memo** collapses the evaluator's fixed-point re-invocation
    (a decline is otherwise recomputed ~3×). A decline is a clean bounded fall-through
    to the series fallback; a solve is fast.
  - **Cascade slot:** after the nonlinear-2nd-order specialists (`ReductionOfOrder`,
    `AutonomousReduction`, `Liouville`) and before the Frobenius series fallback,
    mirroring how first-order Lie sits before series.
  - *Solves:* N3 `x³y''=(y−xy')²`, N4 `x²yy''+x²y'²−5xyy'=4y²`, N5
    `2x²yy''+y²=x²y'²`, the projective family `x²(x+y)y''=(xy'−y)²` and Eq.3
    `y''=(xy'−y)²/x³`, and the arbitrary-coefficient generalizations (`dsolve_m12_stress`).
    *Anti-overfit:* three forward-generator families all 100% back-substitution
    verified — projective `x³y''=a(y−xy')²` (a∈[−4,6]), projective+linear
    `x³y''=(y−xy')²+d x(y−xy')`, scaling `k x²yy''+y²=x²y'²`. *Declines* (bounded,
    no wrong answer): genuinely harder nonlinear equations whose reduced ODE the
    first-order cascade cannot close. Per-call leak is the inherited Integrate/Solve-
    engine leak (as `AlmostLinear`/`LinearFirstOrder`), amplified by the sub-call
    count; ownership within `dsolve_lie2.c` is clean. All DSolve ctest suites +
    `dsolve_m12_stress` + `make check-c99` green.

- **M13 — Abel invariant (AIR).** ⏸ DEFERRED (research-grade). Investigation
  established that the 12000.org `[_rational, _Abel]` corpus (≈9/10 of the
  first-kind examples are `y' = f3(x)y³ + f2(x)y²`) is NOT in the tractable
  constant-invariant class (the canonical `G0/G3` invariant is non-constant; the
  reciprocal `y=1/v → Chini(n=−1)` forms also decline), and needs the full Abel
  Invariant Rational method (Cheb-Terrab's multi-year algorithm; Nasser's own
  solver *and* SymPy both fail on precisely these). A faithful implementation is
  out of a single milestone's reach and a partial one risks wrong answers /
  overfit, so it is deferred rather than hacked. Revisit with the dedicated AIR
  invariant hierarchy + a checked-in solvable-class table.

- **M14 — linear change-of-variable to canonical form (`ChangeOfVariable`).** ✅ DONE.
  Second-order **linear** ODEs with TRANSCENDENTAL coefficients that the direct
  rational methods (Euler/Kovacic/SpecialFunctionForm/Frobenius) miss, because they
  are rational only after a change of the INDEPENDENT variable `t = phi(x)`.
  `dsolve_changevar.c`.
  - For `y'' + P(x) y' + Q(x) y == 0`, the substitution gives `Y'' + A(t) Y' +
    B(t) Y == 0` with `A = (phi''+P phi')/phi'^2`, `B = Q/phi'^2` re-expressed in
    `t`. Candidate `phi ∈ {Cos, Sin, Tan}[x]`, each with the trig-identity rewrite
    (`Sin[x]→√(1−t²)`, … for `t=Cos[x]`) that rationalizes the coefficients; when
    `A,B` come out rational in `t` the transformed equation is handed back to the
    scalar cascade and the solution composed with `t=phi(x)`.
  - Flagship: `y'' + Cot[x] y' + k(k+1) y == 0 --(t=Cos[x])--> (1−t²)Y''−2tY'+
    k(k+1)Y==0` (Legendre), solved by Kovacic (integer degree).
  - Verified NUMERICALLY on the ORIGINAL equation (the composed solution carries
    Legendre-Q `Log` terms whose residual `zero_test` cannot decide, so dsolve_run's
    symbolic verify alone would keep a bad transform). Bounded exactly as M12
    (TimeConstrained sub-solves + wall-clock deadline + decline memo + a re-entry
    guard, since it recurses `DSolve`); transcendental-coefficient gate so it never
    touches an already-rational ODE.
  - Cascade slot: after `Kovacic`, before the reduction/series methods.
  - *Anti-overfit:* two forward-generator families, 9/9 verified — Legendre
    `y''+Cot[x]y'+k(k+1)y==0` (k=1..5) and the Sturm-Liouville spelling
    `y''Sin[x]+y'Cos[x]+m(m+1)ySin[x]==0` (m=1..4) — `tests/test_dsolve_m14_stress.c`
    (+ pinned unit in `test_dsolve.c`). Valgrind leak-flat vs baseline (unlike M12,
    changevar makes fewer sub-calls). *Future:* the associated-Legendre / Pöschl-
    Teller potential `y''−(a+n(n−1)/sin²x)y==0` (transforms rationally but the
    transformed hypergeometric solve needs the non-integer recognizer), Bessel/
    Lommel via power substitution `t=x^k`, and the `E^x`/`Log` (Euler) substitutions.
  All DSolve ctest suites + `make check-c99` green; no regression.

- **M15 — external corpus harness + §2.1.2 baseline + crash hardening.** ✅ DONE.
  The first time DSolve is measured against a large *external* reference set:
  Nasser Abbasi's 12000.org "Solving ODEs" §2.1.2 — 1204 ODEs Maple **and**
  Mathematica both solve (SymPy only 157). Everything lives in
  **`DSolve_test_status/`** (self-contained cross-session dashboard): the converter
  `tools/latex_ode_to_mathilda.py` (tex4ht LaTeX → Mathilda syntax; function/indvar
  detection, subscripts, `\operatorname`/`\textit` Maple heads; all 1391 equations
  round-trip through the parser), the corpus `DE_examples_2.m` (1000 scalar + 204
  systems), the fork-per-case self-verifying harness `test_dsolve_corpus.c` +
  `dsolve_corpus_prelude.m` (runs `DSolve` under `TimeConstrained`, numerically
  back-substitutes each explicit branch; systems skipped — scalar-first), the
  ranked-gap reporter `tools/dsolve_corpus_report.py`, and `STATUS.md`. Registered
  as the per-section ctest `dsolve_corpus_2_1_2_tests` (a progress dashboard gated at
  the checked-in non-PASS baseline, argv[3]; each wave lowers it).
  - **Baseline: 385/1000 scalar (38.5%), 0 FAIL, 6 crashes.** Measurement corrected
    the roadmap: the orthogonal-polynomial (54/54), elliptic (7/7) and Bessel (2/2)
    buckets are **already fully solved** by Kovacic — no recognizer wave needed.
    Ranked scalar gaps: 2nd-order linear 222, 3rd/high-order linear 110,
    2nd-order reducible-μ 70, **Abel 64**, 1st-order symmetry 38, solvable-for-y/x 25.
  - **Two root-cause crash fixes** the corpus surfaced (both protect every caller):
    (A) `dsolve_linear_normalize` now gates its polynomial/rational normalisation to
    **rational-in-x coefficients** (new `ds_is_rational_in`) — a transcendental
    coefficient (`a(λe^{λx}−a e^{2λx})`, `(1+e^{t²/2})²`) drove `PolynomialGCD`/`Cancel`
    to a non-finite content and `GCD(−∞,1)` stack-overflow (SIGSEGV on 298/876/879,
    which now **solve** via Frobenius); (B) `expr_to_mpolyq` guards `base^(−k)` with a
    zero base, which called `fmpz_mpoly_q_inv(0)` and hard-ABORTed FLINT (SIGABRT on
    the nonlinear 3rd-order 269). All 6 crashes now run clean; 0 crashes remain.
  - Method waves M16+ target the ranked gaps (measured-biggest-first). All existing
    DSolve ctest suites + `make check-c99` green; no regression.

- **M16 — 2nd-order linear: change-of-variable + special-function recognizer.** ✅ DONE.
  First measured method wave on the §2.1.2 corpus's largest bucket (2nd-order linear,
  222 gap). Two verified increments, **388 → 396 scalar solved** (gap 612 → 604), 0
  wrong answers, all DSolve suites + `make check-c99` green.
  - **`cv_num_ok` symbolic-parameter verification** (`dsolve_changevar.c`). M14's numeric
    back-substitution guard instantiated only `C[1]`/`C[2]`, so a transform whose
    composed solution carried a **symbolic parameter** (e.g. `LegendreP[k, Cos x]` for
    symbolic `k`) never numericized — every sample was NaN and the correct transform was
    rejected. It now collects the residual's free parameters (argument-position symbols,
    never function heads) and instantiates each at a distinct generic real, so
    symbolic-degree Legendre-Cot (`y''+Cot x y'+k(k+1)y==0`) and its kin now solve
    through M14's existing `t=Cos x` transform. Guard: `t_m16_legendre_symbolic`.
  - **Pöschl-Teller / trigonometric-potential recognizer** (`dsolve_specialform.c`).
    New row for `y'' == (a + p(p−1)Csc²x + q(q−1)Sec²x) y` (and the Csc-only and
    `(a Cos²+b Sin²+c)/Sin²` spellings): extract `Q Sin²x Cos²x` as an even quadratic in
    `Cos x` (rewriting `Sin^{2k}→(1−Cos²)^k`, then **Expand** — `Simplify` reverts it),
    giving `a=−c0`, `p(p−1)=−c1`, `q(q−1)=−c2`; emit the verifiable
    `Sin^p Cos^q ₂F₁((p+q±√(−a))/2, p+½, Sin²x)` and its `p→1−p` partner. Because the
    ₂F₁ residual is undecidable by `zero_test`, emission is gated by an **in-method
    numeric self-verify** (`sf_num_ok`, samples in `(0,π/2)`) — it can never ship a wrong
    answer; a degenerate parameter set (singular ₂F₁ lower parameter) declines to
    Frobenius. Guard: `t_m16_poschl_teller` (3 forms × 2 instantiations). BesselJ/Y,
    LegendreP, Hypergeometric2F1 all numericize, so these solutions are corpus-verifiable;
    GegenbauerC/HermiteH are inert (not used).
  - **Measurement corrected the roadmap:** naive change-of-*independent*-variable
    substitutions (`t=x^m`, `e^x`, `ln x`, `1/x`) had ~0 yield on the homogeneous-rational
    residue (those need special-function recognisers, often at symbolic exponents), and the
    canonical Gegenbauer/Jacobi bucket is already Kovacic-solved. Residue for M17:
    Gegenbauer/Jacobi via an affine→Gauss ₂F₁ change of variable (math validated),
    parabolic-cylinder, power→Bessel, and the Abel/operator-factoring buckets.

- **M17 — 2nd-order-linear: affine→Gauss ₂F₁ recognizer + normal-form pre-pass.** ✅ DONE.
  The largest §2.1.2 bucket (2nd-order linear). Two increments, both in
  `dsolve_specialform.c` (no new method file, no cascade edit), both emitting
  heads that NUMERICIZE (`Hypergeometric2F1`, `BesselJ/Y/I/K`, `AiryAi/Bi`) so the
  in-method numeric self-verify AND the corpus back-substitution are both genuine —
  the 0-FAIL invariant is preserved by construction, not by trust.
  - **`Fuchsian → Gauss ₂F₁` (affine map + local-exponent shift).** A rational-
    coefficient equation with exactly two finite regular singular points `{x1,x2}`
    (the poles of `P`,`Q`: `L = ` squarefree `PolynomialLCM[denom P, denom Q]`,
    `deg L == 2`) is mapped affinely `x = x1 + h s` (`h = x2 − x1`) onto the canonical
    interval `x(1-x)` — `P̃ = h P(x1+h s)`, `Q̃ = h² Q(x1+h s)` — then the local
    exponents at `s=0,1` (indicial roots) are pulled out by the **F-homotopy**
    `Y = s^{r0}(1-s)^{r1} F`, so `F` satisfies the canonical Gauss equation and solves
    as `Hypergeometric2F1`. Reaches **Gegenbauer / Jacobi / associated Legendre at
    symbolic (non-integer) degree/order** — exactly the residue Kovacic declines
    (no Liouvillian solution). The canonical Gauss row was factored into
    `specialform_gauss_basis` and reused by the F-homotopy; the exponent-extraction and
    the mapped coefficients use `Cancel[Together[·]]`, since `Simplify` can mis-reduce a
    constant-over-quadratic pole factor (`3/(4s(s-1)) → −3/(4s)`, dropping a singular
    point → the mapped equation is no longer canonical). Emission gated by `hgc_num_ok`,
    a numeric back-substitution on the ORIGINAL equation sampled INSIDE the mapped
    interval (`s∈(0,1)`), the sole guard against a mis-mapped 2F1 whose contiguous-
    relation residual `zero_test` cannot disprove.
  - **Ordinary vs associated Legendre.** The Legendre row now emits `LegendreP`/`LegendreQ`
    only for the ordinary case (`μ == 0`, which numericizes for symbolic-then-instantiated
    degree); the associated case (`μ ≠ 0`) declines HERE so it falls through to the
    F-homotopy row and is emitted as verifiable `Hypergeometric2F1` — the 3-arg
    `LegendreP[ν,μ,x]` does not numericize for symbolic ν,μ (only integer ν, integer
    μ≥0), so its residual verify and the corpus check both fail. This is what closes the
    corpus trig-potential family (`(p(1+p) − k²Csc²x)y + Cot x y' + y'' == 0`): M14's
    `t=Cos x` reduces it to a rational associated-Legendre equation, which re-enters the
    cascade → this recognizer → 2F1, and M14's `cv_num_ok` verifies.
  - **Liouville normal-form pre-pass.** An equation with a `y'` term is also tried through
    its normal form `z'' = r z` (`y = z Exp[−∫P/2]`, when `∫P/2` is elementary) against the
    y'-free Airy/Bessel recognisers (factored into `specialform_reduced_basis`); a match
    multiplies the recovered basis by the recovery factor. Gated by `sf_num_ok`. Solves the
    Bessel/Airy-reducible equations that carry a first-derivative term (e.g. spherical
    Bessel `y'' + (2/x)y' + y == 0 → Sin[x]/x, Cos[x]/x`).
  - **Cascade slot:** the whole wave lives inside `dsolve_specialform_try`, which already
    sits after the recognizers/Euler and before Kovacic (`dsolve.c:341`). The canonical
    Gauss row runs first, then the affine/F-homotopy row (before Pöschl-Teller), then the
    normal-form pre-pass.
  - *Anti-overfit:* `tests/test_dsolve_m17_stress.c` — forward-generator grids over
    Gegenbauer's λ, Jacobi's (α,β), associated Legendre's order, shifted-interval Gauss
    endpoints, and the normal-form Bessel power, each back-substitution verified with the
    symbolic degree instantiated. Plus pinned units in `test_dsolve.c` (`t_m17_*`).
    *Declines (no wrong answer):* integer 2F1 lower-parameter cases (dependent solutions),
    two-power / Heun potentials, the quartic parabolic-cylinder potential (no native
    `ParabolicCylinderD`), and three-finite-RSP Heun/Möbius equations. All DSolve ctest
    suites + `make check-c99` green; refactor byte-for-byte behavior-preserving on the
    existing Airy/Bessel/Gauss rows.

- **M18 — 2nd-order reducible-μ integrating factor (`ReducibleIntegratingFactor`).** IN
  PROGRESS. The §2.1.2 corpus's largest *well-defined algorithmic* gap: the
  `[_2nd_order, _reducible, _mu_*]` class (80 UNEVAL, all nonlinear). Method:
  Cheb-Terrab & Roche, *Integrating Factors for Second-order ODEs*, J. Symb. Comput. 27
  (1999) 501–519 (PDF in repo). A function μ(x,y,y') is an integrating factor of
  `y''==Φ(x,y,y')` when `μ(y''−Φ)` is a total x-derivative `dR/dx`; from (2.8)–(2.10)
  `μ=R_{y'}`, `R=∫μ dy'+G(x,y)`, `G` fixed by `R_x+y'R_y+ΦR_{y'}==0`. Then `R==C[1]` is a
  first-order ODE → cascade. Two verify gates (symbolic `A(R)==0` = the paper's exactness
  test (2.3), + numeric back-substitution) ⇒ a wrong μ is always a clean decline.
  `dsolve_ifactor.c`, before `SecondOrderSymmetry` in the cascade; reuses the M12 lie2
  robustness kit (deadline / TimeConstrained sub-solves / decline memo / node budget /
  undefined-fn gate) and a **linearity gate** (declines linear ODEs — not the method's
  domain, and it terminates the Case-B ν-ODE recursion).
  - **Stage 1 — μ(x,y) (Section 2.1). ✅ DONE.** Φ a degree-≤2 poly in y'
    (`y''=a y'²+b y'+c`), branch on `2a_x−b_y`: Case A closed-form μ (2.16–2.17); Case B
    `μ=ν(x)e^{−∫a dy}` with ν from one linear ODE (2.18–2.21). Solves `y y'+y''==1`,
    `y''−y y'==6` (→ Airy) and the Coth-coefficient Case-B family
    `−(y'²/y²)+y''/y+2Coth[2x]y'/y==2`. *Side-fix:* `TrigToExp[Coth]` had a sign-flipped
    denominator (=−Coth), surfaced while canonicalizing mixed hyperbolic/exp coefficients;
    fixed in `src/simp/trigsimp.c`. Anti-overfit `tests/test_dsolve_m18_stress.c` (Case-A
    `y''+k y y'==c`, Case-B `y y''−y'²+h(x)y²==0` grids); units `t_m18_*` in
    `test_dsolve.c`. All DSolve ctest suites + trig/hyperbolic suites + `make check-c99`
    green; no regression. **Measured: +7 §2.1.2 reducible-μ solves** (184, 207, 693, 906,
    1013, 1014, 1094), 0 FAIL. Side-fix + symbolic-parameter verify gate (M16 lesson:
    instantiate the residual's free *argument-position* params at generic reals — not
    heads — so symbolic-coefficient ODEs numericize; this unlocked 184).
  - **Stage 2 — μ(x,y') (Section 2.2, Lemma 3).** NEXT INCREMENT. `μ = 𝓕(x,y')·μ̃(x)`
    (2.24-2.25). **𝓕 by Lemma 3** from `Υ = Φ_y` (2.35), six branches:
    **A** `∂_{y'}(Υ_y/Υ)≠0` (2.36) → `𝓕 = 1/(y'-only-not-y factors of Υ)` (2.40);
    **B** those factors free of y' → try A, test μ̃; **C** `G_xy/G_yy≠0` indep of y →
    `w`=y-not-y' factors of Υ, `𝓗=∂_y ln w=𝒢''/𝒢'` (2.47), `p'=𝓗_x/𝓗_y=(w_xy w−w_x
    w_y)/(w_yy w−w_y²)` (2.48), `𝓕=(p'+y')w/Υ` (2.50); **D** `𝓗=0`: `Λ=1/Υ`, `Ψ=Φ/Υ−y`,
    diff (2.62-64) → linear-alg `p'`; **E** `𝓗'=0,𝓗≠0`: `Λ,Ψ` (2.80-83) → linear-alg
    `p'`; **F** (2.83)/Λ_{y'} indep of y' → `β,γ,(2.94)` → linear-alg `p'`. **μ̃ by
    Lemma 2** from `φ1=Φ_y𝓕−y'∂_{y'}(Φ_y𝓕)`, `φ2=∂_{y'}(Φ_y𝓕)`, `φ3=−∂_{y'}(Φ𝓕)`,
    `φ4=∂_{y'}𝓕` (2.30-31): if φ2≠0 `μ̃=Exp∫(φ1_y−φ2_x)/φ2 dx` (2.33), else φ4≠0
    `μ̃=Exp∫(φ3_{y'}−φ4_x)/φ4 dx` (2.34) — the integrand being x-only is the existence
    condition. **Case discrimination:** the Lemma-2 μ̃-existence check is necessary but
    NOT sufficient (a wrong case passes it); the `A(R)=0` gate per candidate is what
    picks the right case. **BLOCKED finding (measured this session):** Cases A/C/D + Lemma
    2 were implemented and VERIFIED to find valid μ (Kamke 226 → μ=y'; Kamke 136 →
    (y'−1)/h(y'); Kamke 66 → (y'+b)/(a(1+y'²)^{3/2}) — all with `A(R)=0` holding), but the
    wave yields **0 new corpus solves**: the reduced first integrals `R==C[1]` are
    NON-ELEMENTARY first-order ODEs (`y'=√(x²y²+2C)`, `y'=Tan[C+Log[x−y]]`) that neither
    our cascade nor — verified directly — Maple/Mathematica close in elementary explicit
    form (those CAS return them implicitly). So the Cases-A/C/D μ-search was reverted, and
    Stage 2 is **blocked on** either (a) a non-elementary/implicit first-order ODE solver,
    or (b) a policy decision to emit the reduced first integral `R(x,y[x],y'[x])==C[1]` as
    an implicit answer (as the existing chini/abel/homogeneous-implicit methods do for
    first-order ODEs). Cases E/F (`𝓗'=0` exponential; the general p'(x)-elimination) were
    not needed. The verified Cases-A/C/D code + μ̃ recovery are recoverable from this
    session's history.
  - **Stage 3 — μ(y,y') (Section 2.3).** Point-swap `y↔x`, reuse Stage 2, back-transform
    `μ=μ_swapped(x,1/y')/y'²` (2.95). Pending Stage 2.
  - *Future:* the 3rd-order `_mu_y2`/`_mu_poly_yn` integrating factors (6 corpus cases);
    a first-order-cascade path for the radical reduced ODEs Stage-2 Case A produces.

- **M19 — confluent Whittaker / ₁F₁ recognizer + §2.1.2 re-baseline.** ✅ DONE. The
  confluent twin of M17's affine→Gauss ₂F₁ row, closing the 2nd-order-linear bucket's
  confluent residue as verified `₁F₁` closed forms. A 2nd-order linear ODE whose
  Liouville normal form `z''==r z` has ONE finite regular singular point `x0` (a double
  pole of `r`) and a rank-1 irregular point at ∞ (`r → b2 ≠ 0`) is Whittaker's equation.
  New `specialform_whittaker_basis()` in `dsolve_specialform.c` (no new method file, no
  cascade edit).
  - **Re-baseline (the plan's pending re-run).** The `reports/2.1.2.tsv` scoreboard was
    stale (pre-M17/M18). A full re-run put the true post-M17/M18 baseline at **424/1000
    scalar (42.4%), 0 FAIL** — not the projected ~403; M17 gained more than estimated.
    This is the honest denominator for the M19 delta.
  - **Method.** `Qc = b2 + b1/(x−x0) + b0/(x−x0)²` matched to `W'' + (−1/4 + κ/z +
    (1/4−μ²)/z²)W==0` under `z = c(x−x0)`: `c = 2√(−b2)`, `μ = √(1/4−b0)`, `κ = b1/c`.
    Emit `Exp[−z/2] z^(1/2±μ) Hypergeometric1F1[1/2±μ−κ, 1±2μ, z]` (numericizes, auto-
    rewrites to `HypergeometricPFQ`, so the residual back-substitutes; the inert
    `WhittakerM/W` heads are never emitted). Single finite double pole isolated by the
    squarefree-part trick (`den/gcd(den,den')` degree 1) — the two-finite-pole Gauss
    case (M17) and the pole-free Airy/polynomial cases decline here. Declines `2μ ∈ ℤ`
    (dependent partners / singular ₁F₁ lower parameter → Frobenius log solution); keeps
    symbolic `μ`. `base` self-verifies against the reduced equation `w'' + Qc w == 0`
    before returning (0-FAIL by construction).
  - **Scope: the y'-free (`P == 0`) surface only.** Deliberately NOT run in the Liouville
    normal-form pre-pass (`P ≠ 0`): there the recovery factor `Exp[−∫P/2]` shares the
    finite-pole base `(x−x0)` with the Whittaker `z^(1/2±μ)` factors, so the composed
    candidate stacks two same-base symbolic-radical powers whose verify/`zero_test`/
    `HypergeometricPFQ`-numeric can hit `$IterationLimit` and STARVE the Frobenius series
    fallback — measured as a **regression** on ≈5 corpus cases (94/470/472/806/811).
    Restricting to `P == 0` (no recovery, no stacking) keeps **0 regressions**.
  - *Solves* the y'-free confluent cases 2.1.2-102 (`x²y''+(cx²+bx+a)y`), -568 and the
    Whittaker/Coulomb normal forms (previously series-fallback). *Future work (biggest
    residue, ~13 corpus cases):* the `P ≠ 0` confluent family (2.1.2-97/-101/-104 and kin)
    — the pre-pass Whittaker with the recovery factor — pending an evaluator-robustness fix
    for the same-base symbolic-radical-exponent verify; also parabolic-cylinder/Hermite
    (`POLY_r ndeg=2`) and the trig/hyperbolic-potential-with-`y'` residues.
  - Anti-overfit `tests/test_dsolve_m19_stress.c` (P==0 Whittaker `(κ,μ)` grid, shifted
    pole, nonzero energy); units `t_m19_*` in `test_dsolve.c`. All DSolve ctest suites +
    `make check-c99` green; Whittaker path valgrind leak-flat.

- **M20 — first-order symmetry gap, Stage 1: `PolynomialShiftSubstitution`.** ✅ DONE.
  First wave on the §2.1.2 `1st_with_symmetry` bucket (65 tot, 38 UNEVAL — Mathilda's
  clearest first-order upside, where SymPy is closest). The 38 subgroup as
  `[F(x),G(x)y+H(x)]` (17), `[F(x),G(y)]` (12), `[F(x),G(x)]` (7), `[F(x)G(y),0]` (2) —
  all *existing* M10 Lie-ansatz classes that decline/abort. Stage 1 takes the largest
  deterministic ABORTing sub-cluster: the **radical substitution** family — `y' = R(x) +
  g(x)(φ(x)+c y)^p` (`p` non-integer, `φ` poly-in-x, `c` const, `R = −φ'/c`), where
  `u = φ+c y` reduces to the separable `u' = c g(x) u^p`, first integral `u^(1−p)/(1−p)
  − ∫c g dx == C[1]` returned **implicitly** (branch-safe; verified by the implicit-
  function rule). New `dsolve_polyshift.c` (`DSolve\`PolynomialShiftSubstitution`), the
  x-dependent-shift generalisation of `FirstOrderSubstitution`.
  - These are Maple's `[F(x),G(x)]`-symmetry cases; the heuristic `abaco2_similar`
    (`dsolve_lie.c`) targets them but its `Q=ω_y/ω_yy`,`T=Q_x/Q_y` differentiates the
    radical `ω`, blows past the node budget, and ABORTs. The deterministic substitution
    sidesteps the symmetry machinery.
  - **Cascade slot:** after `Separable`, before the heavier `Linearizable`/`Exact`/
    `Lagrange` searches (which otherwise spin on the parameter-laden radical before the
    late implicit slot is reached). Its detection needs a fractional-power-of-(linear-in-y)
    atom that the standard methods never produce, so early placement is safe and steals
    nothing (a `y' = a y + b√y` Bernoulli has an additive `a y` → reduced form is not the
    pure `k(x)u^p` → declines).
  - **Robustness:** the reduced form is built by plain evaluation, NOT `Simplify` —
    `Simplify[(x+Sqrt[u]) u^(-1/2)]` loops (radical rationalisation), which the
    `nth_algebraic` branches of a Lagrange equation (`y=2xy'+y'²` → `y'=−x±√(x²+y)`) hit;
    the evaluator's same-base-power/additive cancellation exposes u-freeness without it.
  - *Measured:* **427→432 scalar (+6), 0 FAIL, 0 real regressions.** *Solves*
    2.1.2-402/-371/-372/-376/-403/-424 (previously `abaco2_similar` aborts).
    *Declines (correctly):* 2.1.2-365 (general separable `u'=k(x)(1+I√u)`, not pure power),
    2.1.2-378 (x-dependent `c` → substitution reintroduces y). *Future stages of the
    symmetry gap:* quadratic-in-y' (48/342/981/992 → Factorable/NthAlgebraic robustness);
    arbitrary-function families (~11, need undefined-function symmetry support); the
    transcendental `u=y^{3/2}`/`u=e^{y/x}` analogues; general-separable reduced forms.
  - Anti-overfit `tests/test_dsolve_m20_stress.c` (`(φ,c,g,p)` grids, implicit-function-
    rule verified); units `t_m20_*` in `test_dsolve.c`. All DSolve ctest suites +
    `make check-c99` green.

- **M21 — §2.2.1 corpus (Problems 1–100) + initial-condition coverage.** ✅ DONE. The
  first corpus section carrying **initial value problems**, and the first time the harness
  verifies initial conditions rather than only the general solution's ODE residual. Nasser
  Abbasi §2.2.1 "Table 2.19, Problems 1 to 100" — 100 elementary ODEs (quadrature / linear
  / separable / homogeneous / Riccati), 63 IVPs (4 symbolic `y(a)=b`, 3 swapped-variable
  `x=x(y)`), **zero overlap** with §2.1.2. **86 → 96 / 100 scalar, 0 FAIL, 0 regression.**
  - **Corpus infrastructure (IC support).** The converter `tools/latex_ode_to_mathilda.py`
    now (a) emits initial conditions instead of discarding them — an IVP record's equation
    slot is the DSolve-native list `{ode, ic1, …}` (each ic a point equation `y[x0]==v` /
    `y'[x0]==v`); (b) detects a swapped independent variable (`x=x(y)`, `y` un-primed) and
    an autonomous IC-only parameter (`y(a)=b` → fresh `x`); (c) is section-agnostic
    (auto-selects the problems table, reads columns from the header — §2.2.1 is `TBL-12`
    with ODE at col 2, vs §2.1.2's `TBL-4`/col 3). The prelude `dsolve_corpus_prelude.m`
    passes the list to `DSolve` and verifies **every member** — the ODE residual swept over
    `iv`, and each IC (free of `iv`); an IVP whose solved branch still carries a generated
    constant `C[k]` is scored UNEVAL (general solution, IC unfitted — not a wrong answer).
  - **Verifier bug fixed (all sections).** `dsFreeParams` used `Cases[…, Heads->True]`,
    collecting operator heads (`Plus`, `Times`, `Tan`, `Sec`) as "parameters" and
    substituting numbers for them, so every residual became non-numericizable → a vacuous
    `UNK` (trusted). Removing `Heads->True` (matching the internal `l2_num_ok` "never heads"
    policy) makes the numeric back-substitution real for the first time; §2.1.2 re-verified
    with **0 new FAIL**.
  - **Three solver fixes (chase full coverage).**
    1. **Swapped-variable `A/y'==B`** (2.2.1-98/99/100): `dsolve_nth_algebraic.c` clears a
       top-derivative-bearing denominator (`Numerator[Together[·]]`, top derivative
       restored) and recurses on the cleared ODE, so the linear normalizer owns a
       `y`-dependent leading coefficient. Gated by `!PolynomialQ[Rsub,Dn]` so a normal
       polynomial ODE is never hijacked; presence of the derivative in the cleared form is
       tested with `FreeQ`, **not** `Exponent` (which mis-returns 0 whenever a funcapp is
       present — a separate `src/poly/exponent.c` bug, worked around here).
    2. **Transcendental-inverse IC fit** (2.2.1-29/30/33/34/40/61): `dsolve_fit_constants`
       fits a single-condition/single-constant IVP through Solve's **scalar** form
       `Solve[eq, C[1]]`, since only the scalar form applies inverse-function inversion —
       the list form `Solve[{eq},{C}]` bubbles back on `Sqrt[C]==1`, `Log[3+C]==0`, an Airy
       Möbius ratio, etc., leaking the general solution. A scalar-form empty result is
       treated as "couldn't fit, keep general" (not no-solution: Solve returns `{}` at a
       singular fit point, e.g. a Riccati→Bessel IVP fitted at `x0=0`).
    3. **`ConditionalExpression` principal branch** (2.2.1-60): a multivalued `Tan` inverse
       fits as `ConditionalExpression[…, Element[C[1],Integers]]`; strip it and collapse the
       family index `C[_]→0`, guarded to the ConditionalExpression case so a clean fit is
       untouched → `y'=1+y², y(0)=0` gives `Tan[x]`.
  - **Residue (4, bounded UNEVAL, no wrong answers):** 2.2.1-35 (`y'=Log[1+y²]`,
    non-elementary + missed equilibrium `y≡0`), -47 (homogeneous class-G degree-12 `Root`),
    -48 (slow `ArcSin` separable, >8 s), -67 (`Solve[E^y==Q,y]` reuses `C[1]` as the Log
    branch-index colliding with the integration constant — a `Solve` generated-constant bug,
    separate follow-up). *Follow-ups filed:* the `Exponent`-with-funcapp bug and the
    `Solve` generated-constant collision.
  - New corpus `DSolve_test_status/DE_examples_221.m`; ctest `dsolve_corpus_2_2_1_tests`
    (gate baseline 4); `reports/2.2.1.{tsv,md}`; STATUS.md §2.2.1 block. All DSolve ctest +
    stress suites and `make check-c99` green; §2.1.2 gate held.

- **M22 — §2.2.2 corpus (Problems 101–200) + homogeneous-correctness / exact-transcendental
  / FOS-implicit waves.** ✅ DONE. The continuation of §2.2.1's Table 2.19 (100 elementary
  ODEs, 9 IVPs, 0 systems). **82 → 92 / 100 scalar, 0 FAIL, 0 regression** (§2.1.2 and §2.2.1
  ctests held). Every solver fix reuses the verified implicit first-integral substrate
  (`dsolve_run_implicit` + `dsolve_verify_implicit`), so no wave can ship a wrong answer.
  - **Converter fix (`tools/latex_ode_to_mathilda.py`, all sections).** `is_condition_row`
    matched `<main>·(…)` multiplication (`y²(y'x+y)`, `x(5−x)`) as a `y(P)` initial condition
    and silently dropped 8 ODE rows (109/119/129/166/175–178). Anchored to the LHS: a
    condition row's LHS is *solely* the application `y(…)`/`y'(…)`. Regenerating §2.2.1 is
    byte-for-byte identical (no regression); +6 of the 8 immediately PASS.
  - **Homogeneous correctness (`dsolve_homogeneous.c`).** Retired latent WRONG answers (117
    `y'x=y+√(x²+y²)` returned a spurious extra `x`; 112 `x²y'=xy+x²E^(y/x)` a spurious
    `2 I C[1]π` inverse-branch) and a radical hang (107) — all previously masked as UNEVAL by
    a timeout, i.e. FAILs waiting to surface. Root cause: the reduced RHS was built as
    `F(x,v·x)`, whose `x` does not cancel for a radical/exp `F` (`√(x²+v²x²)` needs `x>0`) and
    leaks into the `v`-integral. Now built as `F(1,v)` (degree-0 homogeneity makes them equal,
    but `x→1` collapses radicals textually). Plus: `homog_exp_log_invert` gated to a log-sum
    antiderivative (the rational-`F` case it was written for); a per-branch **numeric**
    verification drops a spurious inverse-function branch the symbolic verify keeps as
    undecidable (→ implicit fallback); `$rad` internal-placeholder leaks rejected (118). A
    first attempt rejecting `ConditionalExpression[…,C∈ℤ]` in the *shared*
    `dsolve_extract_solutions` was reverted — it broke legitimate periodic `Tan`-inverse
    general solutions (121) and a §2.2.1 case — in favour of the method-local numeric check.
  - **Exact transcendental (`dsolve_exact.c`, `dsolve.c`).** The potential `F(x,Y)` was built
    correctly but `Solve[F==C,Y]` returns unevaluated for a transcendental `F`. Factored the
    build into `exact_potential()` and added `dsolve_exact_implicit_try` returning
    `F(x,y[x])==C[1]` via `dsolve_run_implicit`, wired **immediately after** explicit Exact in
    the cascade so it claims 140/141/142/182/195 *before* a downstream method hangs. A
    `Linearizable` Bernoulli-shape recursion gate (skip a reduced eqn carrying a transcendental
    of `u`, e.g. the `E^(u+E^u)` the Log candidate builds on an `E^y` coefficient) removes the
    pre-Exact hang that otherwise swallowed 141.
  - **FirstOrderSubstitution implicit (`dsolve_fos.c`, `dsolve.c`).** `dsolve_fos_implicit_try`
    returns the inert-integral relation `∫dv/(r+H(v))−x==C[1]` for `y'==f[a x+b y+c]` with an
    arbitrary `f` (159) — the form Mathematica also returns; verified by the implicit-function
    rule, kept on an undecidable inert-integral residual (the keep-the-undecidable policy).
  - **Bernoulli mixed-radical gate (`dsolve_bernoulli.c`).** `(x y)^p` is not the pure
    `B(x) y^n` form (a genuine Bernoulli term's y-power base is `Y` alone), and the exponent
    detector's radical zero-test spun on it; decline early so Homogeneous owns 107.
  - **Residue (8, bounded UNEVAL, no wrong answers):** 133 (`y=G(x,y')` trig), 160 (Bernoulli
    with symbolic exponent `n`), 165 (correct but slow `Sin[x−y]` explicit — implicit is
    future), 170 (elastica `r y''=(1+y'²)^{3/2}`), 175/176 (logistic IVP, flaky under the
    forked 8 s limit), 177/178 (a **pre-existing** general-solution hang on the autonomous
    quadratic `x'=a x(b−x)`; separate follow-up).
  - New corpus `DSolve_test_status/DE_examples_222.m`; ctest `dsolve_corpus_2_2_2_tests`
    (gate baseline 8); `reports/2.2.2.{tsv,md}`; STATUS.md §2.2.2 block. All DSolve ctest +
    stress suites and `make check-c99` green; §2.1.2 / §2.2.1 gates held.

- **M23 — §2.2.3 corpus (Problems 201–300) + exact-radical-potential → implicit wave.** ✅
  DONE. The continuation of §2.2.1/§2.2.2's Table 2.19 (100 elementary ODEs, 35 IVPs, 0
  systems), skewed toward higher-order **constant-coefficient linear** (2nd/3rd/4th order,
  homogeneous + forced, real/repeated/complex roots), **Euler–Cauchy / Emden–Fowler**, and a
  handful of elementary first-order. **98 → 99 / 100 scalar, 0 FAIL, 0 regression** (§2.1.2,
  §2.2.1, §2.2.2 ctests held). The section is near-fully covered out of the box by the
  existing const-coeff (any order), Euler, and first-order specialists; the one solver fix
  reuses the verified implicit first-integral substrate, so it cannot ship a wrong answer.
  - **Corpus infrastructure (converter UNCHANGED).** `DE_examples_223.m` generated by the
    section-agnostic `tools/latex_ode_to_mathilda.py` from `indexsubsection12.htm`
    ("Problems 201 to 300"); §2.2.1/§2.2.2 regenerate **byte-for-byte identical** (no
    converter edit this wave → no regression). New ctest `dsolve_corpus_2_2_3_tests`
    (gate baseline 1); `reports/2.2.3.{tsv,md}`; STATUS.md §2.2.3 block; README row.
  - **Exact radical potential → implicit (`dsolve_exact.c`).** 2.2.3-204
    (`9√x y^(4/3) − 12 x^(1/5) y^(3/2) + (8 x^(3/2) y^(1/3) − 15 x^(6/5)√y) y' == 0`) is EXACT
    (`M_y == N_x`); its potential `F = 6 x^(3/2) y^(4/3) − 10 x^(6/5) y^(3/2)` is built fast,
    but the explicit `Solve[F == C[1], y]` on the mixed FRACTIONAL powers (4/3, 3/2) does not
    terminate — hanging the whole solve. The explicit Exact entry `dsolve_exact_try` now gates
    the inversion on `ds_is_rational_in(F, y)`: a rational-in-`y` potential (Solve reduces to a
    terminating polynomial solve) keeps the explicit path; a radical/transcendental potential
    declines to the existing implicit entry `dsolve_exact_implicit_try`, which returns
    `F(x, y[x]) == C[1]` verbatim (as Maple/Mathematica do), verified by the implicit-function
    rule. The gate never demotes an invertible case — the Laurent-in-`y` `x^a y^b`-factor exact
    (`t_exact_xayb`) and `2xy+1+x²y'==0` still solve explicitly. This is the radical twin of
    M22's transcendental-exact-implicit wave (M22 handled `Solve` returning UNEVALUATED; M23
    handles `Solve` not terminating), and needs no cascade edit — the implicit fallback was
    already wired after explicit Exact (`dsolve.c:351`). Anti-overfit unit `t_m23_exact_radical`
    (implicit-function-rule numeric verify + explicit-preservation guard).
  - **Residue (1, bounded UNEVAL, 0 wrong answers):** 2.2.3-232 `y y'' == 6 x^4` — an
    Emden–Fowler `_with_linear_symmetries` whose (correct) scaling-symmetry reduction
    `r = y/x³, s = ln x` lands on the autonomous `r r'' + 5 r r' + 6 r² == 6`, whose first-order
    reduction `r p p' == 6 − 5 r p − 6 r²` is an **Abel equation of the 2nd kind** —
    non-elementary for the cascade and squarely in the deferred-M13 (Abel Invariant Rational)
    territory. Declines cleanly (verified: the reduced autonomous and intermediate first-order
    both decline in isolation).
  - All DSolve ctest + stress suites and `make check-c99` green; §2.1.2/§2.2.1/§2.2.2 gates
    held.

- **M24 — §2.2.4 corpus (Problems 301–400) + trig-power forcing / numeric-root IVP /
  harness-collision waves.** ✅ DONE. The continuation of §2.2.1/§2.2.2/§2.2.3's Table 2.19
  (100 elementary ODEs, 24 IVPs, 0 systems), skewed toward higher-order **constant-coefficient
  linear** (2nd/3rd/**5th** order, homogeneous + forced), **missing-x/missing-y** reductions,
  nonlinear **`_with_linear_symmetries`**, plus a few Euler/Emden–Fowler/exact/quadrature.
  **93 → 99 / 100 scalar, 0 FAIL (down from 1), 0 regression** (§2.1.2, §2.2.1, §2.2.2, §2.2.3
  ctests held). Two converter fixes made the corpus faithful before any solver work; the
  solver fixes reuse existing verified machinery, so no wave can ship a wrong answer.
  - **Converter fixes (`tools/latex_ode_to_mathilda.py`; §2.2.1/2/3 regenerate byte-for-byte
    identical).** (1) **Imaginary unit `i`** (309/310/311, `y''+2 i y'+3 y=0`,
    `y''=(-2+2 i√3)y`): a constant-coefficient **complex** ODE with no explicit independent
    variable made `detect_symbols` pick the imaginary unit `i` as the indep var and keep it a
    plain symbol — `i`/`I` are now excluded from indep-var candidates (as `e` already was) and
    a standalone `i` maps to Mathilda's `I`. (2) **`y^{(n)}` derivative notation** (336/340/343,
    5th-order): `y^{(5)}` was converted to `y[x]^((5))` (a power of y) instead of the 5th
    derivative; the mains substitution now recognises the parenthesized-order superscript and
    emits `y'''''[x]` (`Derivative[5]`).
  - **Harness verify variable-capture (`dsolve_corpus_prelude.m`, benefits every section).**
    The sole FAIL (387, `m x''+k x==F0 Cos[om t]`, a *correct* fitted solution) was a false
    "BAD": `dsResidVerdict` swept the residual over a loop variable `k` that **collided with the
    ODE parameter `k`** (the spring constant), forcing it to 0…5 instead of its generic sample
    value → a bogus nonzero residual. The sweep index / value holder are now `$`-prefixed
    (`$dsSweep`/`$dsVal`) — names the converter can never emit — so it is strictly more correct
    (fixes false FAILs only; can never create one).
  - **Trig-power/product forcing (`dsolve_undetcoeff.c`, 326/362/363/365).**
    `UndeterminedCoefficients` now `TrigReduce`-linearises the forcing (`Sin[x]^2→(1−Cos2x)/2`,
    `Cos[x]^3→(3Cosx+Cos3x)/4`, `Sin[3x]Sin[x]→(Cos2x−Cos4x)/2`, `x Cos[x]^3→(3x Cosx+x Cos3x)/4`)
    into first-harmonic sinusoids — each a UC function — so a trig power/product forcing solves
    tidily instead of declining to the (hanging) variation-of-parameters fallback. TrigReduce
    preserves value and never turns a UC function into a non-UC one, so it can only help.
  - **Numeric complex roots concretized (`dsolve_common.c` `dsolve_homog_basis`, 312).**
    `y'''==y`'s complex cube roots were emitted as `Re[-(-1)^(1/3)]`/`Im[-(-1)^(1/3)]` (Re/Im do
    not auto-evaluate on a radical power), blocking the IVP constant-fit. The basis builder now
    `ComplexExpand`s the real/imag parts of a **numeric** complex root (gated by `NumericQ`, so a
    symbolic-parameter root — where ComplexExpand could introduce Abs/Sign — is untouched), so
    `y'''==y, y(0)=1, y'(0)=0, y''(0)=0` fits to a concrete solution.
  - **Residue (1, bounded UNEVAL, 0 wrong answers):** 2.2.4-381 `(x²−1)y''−2x y'+2y==x²−1`, a
    variable-coefficient Legendre-type (`_with_linear_symmetries`, SymPy-failed) whose homogeneous
    solves via Kovacic (`y1=x`) but whose nonhomogeneous particular needs variation-of-parameters
    on a Kovacic basis with a rational forcing — a genuine new capability, not a reuse tweak.
    Declines cleanly (Maple solves it; SymPy does not).
  - New corpus `DSolve_test_status/DE_examples_224.m`; ctest `dsolve_corpus_2_2_4_tests` (gate
    baseline 1); `reports/2.2.4.{tsv,md}`; STATUS.md §2.2.4 block; README row. Anti-overfit units
    `t_m24_trig_power_forcing`, `t_m24_complex_cuberoot_ivp` in `test_dsolve.c`. All DSolve ctest +
    stress suites and `make check-c99` green; §2.1.2/§2.2.1/§2.2.2/§2.2.3 gates held.

- **M25 — §2.2.5 corpus (Problems 401–500) + Erf-verify / Kovacic fundamental-set /
  transcendental-Frobenius waves.** ✅ DONE. The continuation of the Table 2.19 sequence, but a
  **series-solution-heavy** chunk (Edwards & Penney, Ch. 8): 2nd-order linear (constant- and
  variable-coefficient), Airy/Emden–Fowler, **Gegenbauer** (already Kovacic), **Bessel**,
  **Jacobi/₂F₁**, one **Liénard**, one 3rd-order (100 records, 15 IVPs, 0 systems).
  **95 → 99 / 100 scalar, 0 FAIL, 0 regression** (§2.1.2, §2.2.1–§2.2.4 ctests held). The
  converter needed no change (§2.2.1–4 regenerate byte-for-byte identical). All three solver fixes
  reuse verified machinery or return a series the Frobenius recurrence gates, so no wave can ship a
  wrong answer.
  - **Erf integrating-factor verify (`dsolve_common.c`, 428).** `y''+x y'+y==0` is exact, reducing
    to the first-order linear `y'+x y==C[2]` whose integrating-factor solution carries
    `Erf[-I x/Sqrt[2]]`. The closed form was produced correctly, but `dsolve_verify_body`'s
    `zero_test` spun on the Gaussian×Erf residual: `E^(-x²/2)·E^(x²/2)` products sit in separate
    summands, never combine to `E^0`, and numericalise as tiny·huge (a catastrophic cancellation),
    so the numeric precision ladder climbs to 1000 bits on every Schwartz–Zippel sample and
    effectively hangs. The verify now `ExpandAll`-normalises a residual **that contains Erf/Erfi**
    before the zero-test (distributing the sums collapses the exponentials); gated to Erf/Erfi so
    every other verify path is byte-for-byte unchanged. The underlying `zero_test` deficiency is
    logged in `POSSIBLE_ZEROQ_IMPROVEMENTS.md` #1 for a later core fix.
  - **Kovacic fundamental-set guard (`dsolve_kovacic.c`) + ExactODE non-elementary decline
    (`dsolve_exactode.c`, 482).** `2x y''+(1-2x²)y'-4x y==0` has one Liouvillian solution
    `√x E^(x²/2)`; its reduction-of-order second solution is non-elementary. Kovacic's
    coincident-exponent path (`Sqrt[D]==0`) collapsed the two basis solutions to the same function
    and returned the rank-deficient `(C[1]+C[2])√x E^(x²/2)` — it verifies (it does solve the ODE)
    but is not a general solution. A final independence guard (`kovacic_body_independent`) extracts
    `y1 = body/.{C[1]->1,C[2]->0}`, `y2 = body/.{C[1]->0,C[2]->1}` and rejects a dependent pair, so
    the cascade falls through to Frobenius, whose two-parameter series is the correct general
    solution. ExactODE also now declines when its reduced sub-solve leaves an unevaluated
    `Integrate` (482's `Integrate[x^(-3/2) E^(-x²/2), x]` is non-elementary here), so 482 reaches
    Kovacic/Frobenius rather than hanging on the junk body.
  - **Transcendental-coefficient Frobenius (`dsolve_frobenius.c`, 463/490).** At a regular singular
    point, `frobenius_regsing` forms `xP = x·P` and `x²Q`, then reads Taylor coefficients by direct
    `x->0` substitution. When P,Q carry an analytic transcendental coefficient this leaves a
    **removable singularity** (`6 Sin[x]/x` is 6 at 0 but substitutes to `6 Sin[0]/0 =
    Indeterminate`), so the indicial roots came out garbage and Frobenius declined. A new
    `normal_series` helper replaces `xP`, `x²Q` by their Taylor polynomials
    (`Normal[Series[·,{x,0,N}]]`, a no-op for genuine polynomials) before the coefficient read, so
    analytic transcendental coefficients yield a verified Frobenius series (463 → indicial roots
    −2,−3; 490 → −1/2, 1).
  - **Residue (1, bounded UNEVAL, 0 wrong answers):** 2.2.5-459 `x² y''+Cos[x] y'+x y==0`, an
    **irregular** singular point at x=0 (`x·P = Cos[x]/x` is not analytic) whose only analytic
    solution is a one-parameter formal power series (the second has an essential singularity). A
    transcendental-coefficient equation Mathilda leaves unevaluated, **matching Mathematica** — the
    shifted-Frobenius path deliberately declines transcendental coefficients rather than expand
    about an arbitrary ordinary point. Declines cleanly.
  - New corpus `DSolve_test_status/DE_examples_225.m`; ctest `dsolve_corpus_2_2_5_tests` (gate
    baseline 1); `reports/2.2.5.{tsv,md}`; STATUS.md §2.2.5 block; README row. Anti-overfit units
    `t_m25_exact_erf`, `t_m25_kovacic_fundamental_set`, `t_m25_transcendental_frobenius` in
    `test_dsolve.c`. `POSSIBLE_ZEROQ_IMPROVEMENTS.md` created (zero_test Gaussian×Erf follow-up).
    All DSolve ctest + stress suites and `make check-c99` green; §2.1.2/§2.2.1–§2.2.4 gates held.

- **M26 — §2.2.6 corpus (Problems 501–600) + general-forcing / DiracDelta wave.** ✅ DONE. Table
  2.29 (Edwards & Penney 6th ed.): 100 records — **74 scalar (47 IVP) + 26 systems** (systems
  skipped by the scalar harness). A **forced-linear** chunk: constant-coefficient 2nd/high-order
  IVPs with **general forcing f(t)** and **DiracDelta impulses**, plus variable-coefficient
  series/Bessel/Emden–Fowler/Liénard and the special Riccati `y'=x²+y²`. **57 → 72 / 74 scalar,
  0 FAIL, 0 regression** (§2.1.2, §2.2.1–§2.2.5, series/reduce/integrate ctests + DSolve stress all
  held). The converter needed a one-line, paren-guarded fix (`\delta(arg)→DiracDelta[arg]`, leaving
  a bare Greek `\delta` parameter — e.g. §2.1.2-544 Heun — as the symbol `delta`). The whole
  forcing family (561–575) now solves via **Green's-function variation of parameters** (no Laplace
  transform); each branch is probe-verified so no wave ships a wrong answer.
  - **Definite-integral variation of parameters (`dsolve_common.c`).** `dsolve_variation_of_parameters`
    kept its indefinite (closed-form) attempt but, when the integral does not close, now builds the
    causal Duhamel convolution `x_p = Integrate[Σ_i basis_i(t)·cof_i(s)·g(s)/(a_n W(s)), {s,0,x}]`
    over a fresh dummy `DSolve`vpS`, where `cof_i = Det(vp_matrix(dv,n,i,e_n))`. `K(t,t)=0` makes
    `x_p` and its first n−1 derivatives vanish at the base point, so a zero-IC IVP fits its constants
    to 0. `TrigReduce`+`Expand` normalise the kernel so a resonant cos/sin forcing closes to the
    clean `t Sin[w t]` (569 → `t Sin[3t]/6 − Sin[3t] HeavisideTheta[t−3π]/3`) and an exponential
    kernel integrates termwise. Covers 561–563, 572–575 (arbitrary f → convolution integral).
  - **DiracDelta sifting under a definite integral (`integrate_dirac.c`/`.h`, new).**
    `integrate_dirac_try`, called first in `integrate_definite`'s real-axis branch:
    `∫ DiracDelta[αx+β] h(x) dx over [lo,hi] → (h/.x→x0)/|α|·B`, `x0=−β/α`. Boundary `B`:
    `HeavisideTheta[hi−x0]` for a symbolic upper limit with `x0≥lo` (causal convention — the FULL
    step even at `x0=lo`, matching MMA), and `1 / 0 / ½` (interior / outside / endpoint) for numeric
    limits (`Positive`/`Negative` so a `Pi` shift is decided, not the inert `Sign[Pi]`). A mixed
    integrand (1+δ, t+δ, δ+cos) is `ExpandAll`-split by linearity; the delta is located anywhere in
    a (possibly nested) product. Covers 564–571 and the standalone `Integrate[δ·f]` gap.
  - **HeavisideTheta / DiracDelta rules (`distributions.m` new, loaded from `init.m`; `deriv.c`).**
    `H(0)=0`, `H(x>0)=1`, `H(x<0)=0`, `δ(x≠0)=0` on definite-sign numeric arguments (symbolic left
    inert), and `d/dg HeavisideTheta[g] = DiracDelta[g]` in `elementary_fprime`. Left-continuous
    `H(0)=0` makes a causal impulse response satisfy its pre-impulse ICs, and the numeric rules let
    an IC fit at the base point resolve a shifted `H[t−a]` / `δ[t−a]`. `dsolve_constcoeff.c` gates
    its "homogeneous" branch off a DiracDelta forcing (which samples numerically to 0, so the
    numeric zero-test would wrongly drop it).
  - **`D[]` Leibniz rule for a variable-limit integral (`deriv.c`).** `D[Integrate[e,{u,a,b}],x] =
    (e/.u→b)·D[b,x] − (e/.u→a)·D[a,x] + Integrate[D[e,x],{u,a,b}]` (guarded to a 3-element List spec,
    bound var ≠ x); with the equal-limits rule `Integrate[_,{s,a,a}]→0` (`integrate.c`) it lets the
    convolution IVP fit cleanly. Replaces the prior garbage output that leaked the bound variable.
  - **Anti-hang guards.** VoP skips its indefinite attempt for an arbitrary/undefined forcing (it
    never closes and can hang — `ds_has_undefined_function`), and `integrate_definite` skips the
    improper/parametric methods (residue/Ramanujan/differentiation-under-the-integral) on an
    undefined-function integrand: they cannot close a `K(t,s) f(s)` convolution and churned for
    seconds (an exp kernel dropped 6 s → 0.3 s per eval, so the DSolve fixed-point no longer times
    out). `dsolve_verify_body` keeps (never rejects) a distributional residual (definite Integrate /
    DiracDelta / HeavisideTheta) instead of driving `zero_test` into a spin.
  - **Residue (2, bounded UNEVAL, 0 wrong answers):** 2.2.6-524 `y''+x⁴ y==0` (Emden–Fowler; the
    correct `√x BesselJ/Y[1/6, x³/3]` form is returned but exceeds the 8 s per-case DSolve budget —
    a performance residue) and 2.2.6-555 `t x''+(t-2)x'+x==0` (exact → the regular-singular
    reduction `t x'+(t-3)x==C[2]` whose integrating-factor quadrature `∫E^t/t⁴` is non-elementary
    (`ExpIntegralEi`); a series residue, declines via timeout). Both pre-existing (UNEVAL in the
    baseline), no wrong answers.
  - New corpus `DSolve_test_status/DE_examples_226.m`; ctest `dsolve_corpus_2_2_6_tests` (gate
    baseline 2); `reports/2.2.6.{tsv,md}`; STATUS.md §2.2.6 block; README row. Anti-overfit units
    `t_m26_distributions`, `t_m26_impulse_forcing`, `t_m26_general_forcing` in `test_dsolve.c`.
    Converter `\delta` fix in `tools/latex_ode_to_mathilda.py`. All DSolve ctest + stress suites,
    series/reduce/integrate suites, and `make check-c99` green; all prior corpus gates held.

- **M27 — §2.2.7 corpus (Problems 601–700) + corpus-harness SYSTEM VERIFICATION.** ✅ DONE.
  Table 2.31: 100 records — **50 scalar (23 IVP) + 50 systems** (25 2-D, 18 3-D, 7 4-D; 47
  constant-coefficient + 3 variable-coefficient). The first section that is **half systems**,
  and the wave that taught the corpus harness to **verify systems** rather than skip them
  ("scalar-first" ran M15–M26). The scalar half is elementary first-order (separable /
  quadrature / linear / homogeneous class-G / Riccati). **90 → 93 / 100** after two engine
  fixes, **0 FAIL, 0 regression** (§2.2.1–§2.2.5 gates held; §2.1.2 and §2.2.6 re-baselined).
  The converter needed no change (the `--label 2.2.7` run of the section-agnostic
  `latex_ode_to_mathilda.py` produced a clean 100-record file; §2.2.1–§2.2.6 regenerate
  byte-for-byte identical in their records — the only header edit is a wording change from
  "the scalar harness skips them" to a system-shape description).
  - **System verification in the harness (`dsolve_corpus_prelude.m`).** The numeric
    back-substitution machinery (`dsResidVerdict` / `dsBranchVerdict`) already substituted a
    whole rule-list `/. br` into each residual, so it was already multi-function-capable; only
    `dsExplicitQ` hard-coded a single function symbol and `dsolveCheckCode` returned SKIP for a
    system. `dsExplicitQ` now accepts a List function slot (every rule resolves to a
    `Function`, and every dependent function is present — a partially-solved system is not
    "explicit"), and the system-skip is removed. A system branch
    `{x->Function[…], y->Function[…], …}` is back-substituted per equation exactly like a
    scalar ODE; a general system solution's `C[1..n]` behave as free constants under the sweep
    (residual ~0 for all sampled values). The verifier's `leaked→UNFIT` rule (a demonstrably-
    nonzero residual with `C[k]` still present scores UNEVAL, not FAIL) means a wrong *general*
    system solution cannot become a FAIL — the FAIL surface is only fitted IVP systems, of
    which §2.2.7 has none. **Prior sections carrying systems were re-baselined:** §2.1.2 (204
    systems) and §2.2.6 (26 systems, 72/74 scalar → 95/100 total = 72 scalar + 23 systems,
    baseline 2→5); §2.2.1–§2.2.5 are pure scalar and unchanged. `tools/dsolve_corpus_report.py`
    reports scalar and system populations separately.
  - **Separated-exponent Simplify hang (`dsolve_linsys.c`).** A constant-coefficient system as
    ordinary as `x'=-50x+20y, y'=100x-60y` (real eigenvalues -10, -100) HUNG uninterruptibly:
    `dsolve_linsys_tidy` gave a small real-exponential body the "pretty" `Simplify`, and
    `Simplify`'s zero-test spins on a **sum of exponentials with widely-separated decay rates**
    (`Simplify[Exp[-10 t]+Exp[-100 t]]` does not return — it numericises `E^(-10 t)` against
    `E^(-100 t)`, whose dynamic range drives the precision ladder to its ceiling; close rates
    like -1,-2 are fine). The `heavy → Expand` gate (already used for complex-spectrum and
    large bodies) now also covers **any** exponential body, so a real-exponential body never
    reaches `Simplify`; the result is back-substitution-verified regardless. +2 systems
    (636, 650). Documented in `POSSIBLE_ZEROQ_IMPROVEMENTS.md`.
  - **Solve periodicity-index collision (`solveinv.c` + `dsolve_common.c`).** `y'=2x Sec[y]`
    shipped a WRONG general solution (masked as UNEVAL by leaked→UNFIT): `DSolve\`Separable`
    feeds Solve an equation already carrying the integration constant `C[1]`, and the
    inverse-trig peeler minted `C[1]` again as the `2πk` periodicity index (its counter started
    at 0, ignoring existing constants); once `dsolve_extract_solutions` stripped the
    `Element[C[1],Integers]` constraint, the shared `C[1]` shifted `y` by a non-multiple of 2π
    for non-integer values. Two root-cause fixes: (a) `solveinv` **seeds its mint counter with
    the largest `C[k]` index already in the equation** (`max_param_index`), so a family
    parameter is always fresh; (b) `dsolve_extract_solutions` **collapses the integer family to
    its principal branch** (`Element[C[k],Integers]`-constrained `C[k] → 0`), scoped strictly to
    the `Element[…,Integers]` atom — a *range* condition on the integration constant (e.g. the
    `-π/2 < x³+C[1] ≤ π/2` from a `Tan`/`ArcTan` inversion) legitimately mentions `C[1]` and
    must NOT be collapsed (that regressed `y'=3x²(1+y²), y(0)=1` before the scoping was
    tightened). +1 (684). This is the pre-existing Solve generated-constant collision M21 filed
    as a follow-up (also 2.2.1-67). `dsolve_tests`/`solve_tests`/`reduce_tests`/`solve_corpus`
    all green (no regression from the minting change).
  - **Residue (7, bounded UNEVAL, 0 wrong answers):** 603/606/607 (forced constant-coefficient
    systems whose closed form with irrational/complex eigenvalues exceeds the 8 s per-case
    budget — a performance residue); 604/608 (variable-coefficient non-triangular systems —
    the honest `LinearSystemVarCoeff` gap, matching the plan's "one honest gap"); 675
    (`y'=Log[1+y²]`, non-elementary, matches Mathematica); 683 (`y'=4(x y)^(1/3)`, homogeneous
    class-G whose implicit inversion exceeds the 8 s budget).
  - New corpus `DSolve_test_status/DE_examples_227.m`; ctest `dsolve_corpus_2_2_7_tests` (gate
    baseline 7); `reports/2.2.7.{tsv,md}`; STATUS.md §2.2.7 block + M26/M27 wave-history bullets;
    README row + systems-now-verified note. Anti-overfit units `t_m27_system_verify`,
    `t_m27_separable_inverse_constant`, `t_m27_ivp_family_intact` in `test_dsolve.c`. All DSolve
    ctest + stress suites (`dsolve`, m5/m12/m14/m17/m18/m19/m20), `solve`/`reduce`/`solve_corpus`,
    and `make check-c99` green.

- **M28 — §2.2.8 corpus (Problems 701–800) + Bernoulli cascade-hang fix.** ✅ DONE. Table 2.33
  (Edwards & Penney): 100 records — **100 scalar (20 IVP), 0 systems**, a return to elementary
  first-order after §2.2.7's half-systems chunk. The same families as §2.2.1–§2.2.4 (23 linear,
  16 separable, ~25 homogeneous class A/G/C, 11 exact incl. two fractional-power potentials, 19
  Bernoulli several with fractional exponent, 3 quadrature, plus a few `y'=F(ax+by+c)` /
  Riccati). **98 → 99 / 100 scalar, 0 FAIL, 0 crashes, 0 regression** (§2.1.2, §2.2.1–§2.2.7 gates
  held). The converter needed no change (`--label 2.2.8` on the section-agnostic
  `latex_ode_to_mathilda.py` produced a clean 100-record file; §2.2.1–§2.2.7 regenerate
  byte-for-byte identical). One root-cause engine fix reusing verified machinery, so no wave ships
  a wrong answer.
  - **Bernoulli fast-decline on a transcendental-in-`y` RHS (`dsolve_bernoulli.c`, 757).**
    `2 x Sin[y]Cos[y] y' == 4 x² + Sin[y]²` reduces (`u = Sin[y]²`) to the linear `u' − u/x == 4 x`,
    which `Linearizable` solves in one step — but solved for `y'` the RHS is a rational function of
    `Sin[y]`/`Cos[y]`, transcendental in `y`, NOT the Bernoulli form `A(x) y + B(x) y^n`. The
    exponent detector (`Q = Y − Y F_Y`, then `n = Y Q_Y/Q` via `Cancel`, then `ds_free_of`) spun for
    seconds on the trig-rational expression, timing out the whole cascade (`$Aborted`) before
    `Linearizable` was reached. A new early-decline gate `bern_Y_nonalgebraic` returns NULL
    immediately when `y` appears inside a non-`Power` function head (`Sin[y]`, `Exp[y]`, …) or in a
    power exponent — mirroring the existing `bern_mixed_radical` guard; a genuine Bernoulli
    (algebraic in `y`, e.g. 752 `F = y − E^(-2x)/(2x) y³`) is untouched (the reconstruction check
    already rejected any transcendental F, so this only makes the decline FAST, never changes an
    answer). +1 (757).
  - **Residue (1, bounded UNEVAL, 0 wrong answers):** 2.2.8-783 `y'=1+x²+y²+x²y⁴`, a quartic-in-`y`
    equation (beyond Abel) that Maple/Mma/SymPy solve only via the general "solve-for-`y` then
    differentiate" method (Maple's `y=_G(x,y')` class) — a `SolvableForY`/generalised-d'Alembert
    method Mathilda does not have. Declines cleanly (no hang, no wrong answer); future work.
  - New corpus `DSolve_test_status/DE_examples_228.m`; ctest `dsolve_corpus_2_2_8_tests` (gate
    baseline 1); `reports/2.2.8.{tsv,md}`; STATUS.md §2.2.8 block + M28 wave-history bullet; README
    row. Anti-overfit unit `t_m28_bernoulli_hang_trig_substitution` in `test_dsolve.c` (757 solves +
    back-substitutes, the general `a`-coefficient family solves, and a genuine Bernoulli 752 still
    solves via Bernoulli). All DSolve ctest + stress suites, `solve`/`reduce`, and `make check-c99`
    green; §2.1.2/§2.2.1–§2.2.7 gates held.
- **M29 — §2.2.9 corpus (Problems 801–900) + piecewise rounding-function derivatives.** ✅ DONE.
  Edwards & Penney, Problems 801–900: 100 records — **100 scalar (35 IVP), 0 systems**, dominated
  by second-order linear — 38 constant-coefficient homogeneous, 35 constant-coefficient
  nonhomogeneous (undetermined coefficients / variation of parameters), 11 Euler/Emden–Fowler —
  plus 7 with `x(t)` as the dependent variable (`x''` notation, 862–868), 3 complex-coefficient
  (`i` in the equation: 857/858/859), and 6 first-order (separable / homogeneous / Bernoulli /
  Abel). **100 / 100 scalar, 0 FAIL, 0 crashes, 0 regression, baseline 0** — Mathilda's strongest
  DSolve territory, fully covered out of the box by the existing const-coeff (any order), Euler,
  variation-of-parameters/Green's-function and first-order specialists. The converter needed no
  change (`--label 2.2.9 --url …indexsubsection18.htm` on the section-agnostic
  `latex_ode_to_mathilda.py`; complex `i`→`I` and `x''[t]` independent-variable inference already
  handled; §2.2.1–§2.2.8 regenerate byte-for-byte identical in their records). **No ODE-solver fix
  was required.** Instead the wave landed one general engine feature the corpus made visible:
  - **Piecewise differentiation of the integer-rounding functions (`src/calculus/deriv.c`).**
    `D` of `Floor`/`Ceiling`/`Round`/`IntegerPart`/`FractionalPart` now returns the Mathematica
    `Piecewise[{{v, cond}}, Indeterminate]` form (value 0 off the jump set, 1 for `FractionalPart`;
    `Round`'s jump set is the half-integers, `IntegerPart`/`FractionalPart` guard `Re`/`Im`
    non-integrality), multiplied by `D[g,x]` so the chain rule composes. Modeled on the existing
    `UnitStep` derivative handler. Problem `2.2.9-898` (`y''+9y == 2 Sec[3x]`) is solved by
    variation of parameters and its verified solution carries a `Floor` branch-tracking term;
    before this, `D[Floor[u],x]` was the inert `Derivative[1][Floor][u]`, so the ODE residual never
    numericized and the corpus harness passed 898 only under the "non-numericizable ⇒ trust DSolve"
    policy. The residual now reduces to a genuine numeric ~0 — 898 is verified, not merely trusted.
  - **Residue: none** (0 UNEVAL, 0 FAIL).
  - New corpus `DSolve_test_status/DE_examples_229.m`; ctest `dsolve_corpus_2_2_9_tests` (gate
    baseline 0); `reports/2.2.9.{tsv,md}`; STATUS.md §2.2.9 block + M29 wave-history bullet; README
    row. Anti-overfit units `t_m29_sec_floor_verifies` (`test_dsolve.c`: the Floor derivative
    numericizes, 898 + a sibling `Sec`-forced equation solve and back-substitute to numeric zero)
    and `test_rounding_deriv` (`test_deriv.c`: the five exact `Piecewise` forms + chain rule). All
    DSolve ctest + stress suites, `deriv`, and `make check-c99` green; §2.1.2/§2.2.1–§2.2.8 gates
    held.
- **M30 — §2.2.10 corpus (Problems 901–1000) + forced-system / Kovacic-inhomogeneous fixes.** ✅
  DONE. Edwards & Penney 901–1000: 100 records — **56 scalar (15 IVP) + 44 first-order linear
  systems** (the most systems-heavy section; 2×2/3×3/4×4 constant matrices). **100/100, 0 FAIL,
  0 regression, baseline 0.** Two root-cause fixes, both about *forced* linear ODEs: (1) a
  real-irrational-spectrum forced system (924, `(1±√89)/2`) ran >90 s because `Integrate`
  rationalised the `1/λᵏ` variation-of-parameters coefficient into a hundreds-of-digit integer —
  `dsolve_linsys.c` now abstracts a real-irrational eigenvalue to a symbol before the integral
  (complex/rational spectra stay concrete); (2) Kovacic gained an inhomogeneous closure
  (`dsolve_kovacic.c` + `dsolve_second_order_PQ_forced`: accept a forcing, de-obfuscate the
  fundamental set, add a variation-of-parameters particular), solving the Legendre-type 907.
  Anti-overfit units `t_m30_linsys_irrational_forcing`, `t_m30_kovacic_inhomogeneous`. See the
  §2.2.10 block in `DSolve_test_status/STATUS.md`.
- **M31 — §2.2.11 corpus (Problems 1001–1100) + linear-integrand combine / cancellation-robust
  verifier.** ✅ DONE. Edwards & Penney 1001–1100: 100 records — **59 scalar (13 IVP) + 41
  first-order constant-coefficient linear systems** (2×2 … 6×6, defective/repeated/complex
  spectra). The scalar half (separable / quadrature / first-order-linear / 2nd-order const-coeff
  + Euler + Gegenbauer + Emden–Fowler + Liénard + Airy) solves entirely out of the box; both
  gaps were systems. **§2.2.11 100/100, 0 FAIL, 0 crashes, 0 regression, baseline 0.** The
  converter needed no change (the §2.2.10 `indexsubsection→section` scheme extends to
  `indexsubsection20.htm`; §2.2.1–§2.2.10 regenerate byte-for-byte identical). Two root-cause
  fixes:
  - **`Integrate` linearity over a distributed product (`src/calculus/integrate.c`,
    `try_linearity`).** `2.2.11-1014` (`{x1'=2x1, x2'=−7x1+9x2+7x3, x3'=2x3}`, a DAG solved by
    `TriangularSystem`) asked `Integrate` for `Integrate[e^{−9x}(7 C[k]e^{2x}−7 C[j]e^{2x}), x]` —
    a product of an exponential with a **sum** of exponentials (`Times[c, Plus[…]]`, exponents
    uncombined) — which took the exponential-substitution path: **55 s** and a **branch-wrong**
    `(−1)^{1/9}` antiderivative (accepted because its residual is zero-test-*undecidable*). The
    **root fix** is in `Integrate`: `try_linearity` (the Plus-splitting stage, ahead of the
    substitution stages) now distributes a product over a sum factor — `Integrate` is linear, so
    `c(g+h) → cg+ch`, committed only if every term closes elementary (else the whole-integrand
    cascade still runs). The exponentials then collapse (`e^{−9x}e^{2x}→e^{−7x}`) — fast and
    correct. This *also* repairs the user-reported **direct** bug
    `Integrate[e^{−9x}(a e^{2x}−b e^{2x}), x]` (was `−(a−b)/7·(e^{2x})^{−7/2}` in ~9 s, and
    genuinely branch-wrong for symbolic/funcapp coefficients; now `−(a−b)/7·e^{−7x}` in ~4 ms).
    Guarded by `test_linearity_distributes_product` (`test_integrate_dispatch.c`); all 22 integrate
    unit suites stay green (a fully-elementary split is the only thing committed).
  - **Corpus verifier made cancellation-robust (`dsolve_corpus_prelude.m`, `dsResidVerdict`).**
    `2.2.11-1001` (4×4, eigenvalues {16,32,48,64}) solves **correctly** (`Simplify[residual]≡0`)
    but back-substitutes to a difference of `e^{64x}`-scale terms that, at the prelude's 20-digit
    sweep over `x≈1.1…3`, looks large (catastrophic cancellation) → a false "BAD" → UNEVAL. The
    shared verifier now re-evaluates a not-small residual sample at 200-digit precision; the
    change is **monotone** (can only turn a spurious "not small" into "small" — never introduces
    a FAIL, never raises a section's non-PASS count). It also lifted four earlier sections whose
    correct-but-cancellation-heavy answers now verify, re-baselined to match.
  Anti-overfit units `t_m31_triangular_exp_forcing`, `t_m31_linsys_large_eigenvalue`
  (`tests/test_dsolve.c`). All DSolve ctest + stress suites and `make check-c99` green;
  §2.1.2/§2.2.1–§2.2.10 gates held (four improved). See the §2.2.11 block in
  `DSolve_test_status/STATUS.md`.

- **M32 — §2.2.12 corpus (Problems 1101–1200) + IVP-fitter / Separable-implicit /
  Integrate-Erf-variable fixes.** ✅ DONE. Edwards & Penney 1101–1200: 100 records —
  **100 scalar (38 IVP), 0 systems**, elementary first-order (separable / linear / quadrature /
  homogeneous / Bernoulli / exact, plus homogeneous-class-A "Abel" rationals and two
  solvable-for-y/x forms). **The first §2.2.x section NOT green out of the box:** baseline
  **80/100 with a wrong answer (1 FAIL) + 19 UNEVAL → 97/100, 0 FAIL, 0 crashes, 0 regression.**
  Three shared-substrate root-cause fixes (each lifts earlier sections too, none regress):
  - **IVP constant-fitting (`src/calculus/dsolve_common.c`).** `dsolve_fit_constants` now reports a
    per-branch `fit_state` (OK / EMPTY / UNDEF) and `dsolve_run` drops, WITH sibling context, a
    branch an initial condition cannot be met on: one whose fit substituted `Undefined`/`$Failed`
    (Solve had no consistent constant — `2.2.12-1147`, the section's sole **FAIL/wrong answer**),
    or a scalar unsatisfiable inverse branch (wrong `±`/`Root` index) **when a sibling actually
    fits** — so a lone basis singularity is still kept and a genuinely UNDER-determined BVP
    (`y''+y==0, y[0]==0, y[π]==0 → C[2] Sin[x]`) keeps its free constant. `dsolve_bernoulli.c` emits
    BOTH real signs for an even `1−n` root (`y'==(1−2x)/y, y[1]==−2` needs `−√`). Fixed the FAIL +
    11 UNEVAL separable IVPs, including the Root-form cubics 1149/1150 (constant fitted on the
    implicit first integral `G(x0,y0)` with no inversion, via the new twin below).
  - **Separable recognizer + implicit twin (`src/calculus/dsolve_separable.c`).** The 36-sample
    split search is factored into `sep_find_split`; the gate now rejects a sample only when the
    denominator is PROVABLY zero (accepting generic-parameter splits like `(a y+b)/(c y+d)`); and a
    new `dsolve_separable_implicit_try` (dispatched via `dsolve_run_implicit`, mirroring
    Exact/ExactImplicit) returns the first integral `∫dy/g == ∫f dx + C[1]` — keeping a
    non-elementary integral **unevaluated** — when the relation does not invert for y. Solves
    `Cot[t]y/(1+y)`, `Cos²x Cos²2y`, and the autonomous non-elementary `−2 ArcTan[y]/(1+y²)`.
  - **Integrate Gaussian→Erf/Ei/PolyLog recognizer variable (`src/calculus/risch_special.c`).** The
    completing-the-square templates (`rt_try_erf`/`rt_try_ei`/dilog) emitted the antiderivative in a
    **literal `x` for every integration variable** — `Integrate[E^(a^2/2), a]` came back
    `… Erf[… x …]` — so any Bernoulli/linear DSolve over a non-`x` variable produced a stray-variable,
    non-verifying answer. The real variable is now threaded through the template (ReplaceAll does not
    re-traverse substitutions and the coefficients are free of `x`, so a parameter named `x` is safe).
    Fixed 1182/1190.
  - **Residue (3, research-grade, bounded declines — no wrong answers):** `1135` (`y=_G(x,y')`,
    transcendental generalised-d'Alembert — induced `p`-ODE not cascade-solvable), `1200`
    (`x=_G(y,y')` — cannot be solved for `x`: `e^x` + linear `x`), `1157` (`(a y+b)/(c y+d)` with `a`
    the independent variable — Abel 2nd kind, the M13-deferred class). A general
    generalised-d'Alembert method was prototyped and confirmed to close NONE of these, so it is
    deferred to its own milestone rather than half-built here. Anti-overfit units `t_m32_*`
    (`tests/test_dsolve.c`); `make check-c99` green; §2.1.2/§2.2.1–§2.2.11 corpus gates held. See the
    §2.2.12 block in `DSolve_test_status/STATUS.md`. (`t_rischnorman_enum_cap_no_crash` exceeds the
    120 s whole-binary alarm in `tests/test_utils.h` on slower hardware — a pre-existing local-only
    condition identical on pristine `main`; CI does not run `dsolve_tests`.)

- **M33 — §2.2.13 corpus (Problems 1201–1300) + exact-method robustness / separable
  fast-decline / exact-vs-homogeneous IVP fall-through.** ✅ DONE. Edwards & Penney 1201–1300:
  100 records — **100 scalar (32 IVP), 0 systems**, a MIX of first-order (exact / linear /
  separable / homogeneous / Abel / symmetry) and 2nd-order linear (const-coeff,
  Euler–Cauchy / "Emden–Fowler", reducible). Baseline **91/100 (0 FAIL, 9 UNEVAL) → 99/100,
  0 FAIL, 0 regression.** Residue: `2.2.13-1203`, an Abel-2nd-kind (class B) — the
  research-grade `[_Abel]` class deferred at M13. Five shared-substrate root-cause fixes
  (each lifts earlier sections too; none regress):
  - **Exact potential Path-1/Path-2 + syntactic-denominator clearing + robust `mu(y)`
    (`dsolve_exact.c`).** (a) Build the potential from whichever coefficient integrates
    cleanly — `∫M dx` (Path 1) or `∫N dy` (Path 2), NO `Simplify` on the hot path — so the
    E^(x y) exact family (1201), whose `∫M dx` lands in a Tan half-angle form that leaves
    `g'(y)` only *cancelling* to a function of y, is solved via the clean Path 2. (b) Clear an
    inexact rational form (`y'==P/Q`, `1/x`/`1/y` poles) by its common denominator D (an
    integrating factor), trying two exactness-gated candidates — the SYNTACTIC denominators
    (handles a NEGATIVE exponential `E^(-x)`, which `Together` mis-factors as a spurious `E^x`
    — 1233), then the `Together` denominator (summed `1/x`,`1/y` — 1216/1238); restricted to
    the condition-free general solve (an IVP's cleared Root form does not inverse-fit). (c)
    `mu(y)` is tried whenever `mu(x)` yields no factor (its free-of test can mis-decide on a
    trig-rational derivative), fixing `mu = Sin y` equations (1214).
  - **Implicit verify `Together`s the residual (`dsolve_common.c`).**
    `dsolve_verify_implicit` combines `F_x + N(-F_x/F_y)` over a common denominator before the
    zero-test: a `mu = Sin y` exact equation's telescoping residual has uncancelled `Csc`/`Cot`
    poles that the numeric zero-test FALSE-NEGATIVES (proved "nonzero"), wrongly REJECTING a
    correct branch (1214). Value-preserving, so it cannot mask a real nonzero.
  - **Separable fast numeric pre-filter (`dsolve_separable.c`).** `sep_find_split` numerically
    samples the separability check `F - g·h` at a generic real point and skips the expensive
    symbolic zero-test when it is clearly nonzero — a non-separable transcendental RHS
    (1201/1233) fast-declines in ~0 s not ~15 s (the pre-Exact cascade cost that pushed those
    solvable equations past the 8 s per-case timeout). Never rejects a genuine split.
  - **First-order IVP undecided-fit fall-through (`dsolve_common.c`).** New `FIT_UNDECIDED`
    state: when a scalar FIRST-ORDER IVP's fit bubbles back unevaluated (Solve could not
    resolve the single constant) and no branch achieved a real fit, `dsolve_run` DECLINES so
    the cascade continues — closing the exact/homogeneous overlap 1205/1231, where Homogeneous
    runs first and returns a transcendental log-form whose constant does not inverse-fit,
    shadowing Exact's fitting polynomial first integral. Gated to `nfun==1 && max_order==1`:
    a legitimately UNDER-determined BVP keeps its free constant (FIT_OK), and a higher-order
    series IVP whose SeriesData fit bubbles yet verifies is untouched (else it would decline
    correct §2.2.5/§2.2.6/§2.2.11 2nd-order IVPs — caught and gated during the wave).
  Anti-overfit units `t_m33_*` (`tests/test_dsolve.c`); `make check-c99` green; §2.1.2 /
  §2.2.1–§2.2.12 corpus gates all held at baseline. See the §2.2.13 block in
  `DSolve_test_status/STATUS.md`. Version 0.130 → 0.131.

- **M34 — §2.2.14 corpus (Problems 1301–1400, Boyce & DiPrima) + VoP verify short-circuit /
  robust variation of parameters / bounded-Kovacic complex-pole gate / SeriesData IVP fit /
  IC-point Frobenius series.** ✅ DONE. The first Boyce & DiPrima section: 2nd-order-linear
  dominated (38 with-symmetry, 19 const-coeff, 13 nonhomogeneous, 12 Emden–Fowler, 11 exact),
  99 scalar (25 IVP) + 1 system. **90/100 → 99/100, 0 FAIL, 0 regression** via four shared-
  substrate root-cause fixes (they lift earlier sections, none regress):
  1. **Numeric-zero verify short-circuit + robust VoP** (`dsolve_common.c`). `dsolve_verify_body`
     keeps a branch whose residual is NUMERICALLY zero at a spread of clean real points before the
     symbolic `zero_test_decide`, whose precision ladder climbs for >8 s on a residual that IS zero
     but carries Log/ArcTan branch cuts (`y''+y==Tan[x]` / `2 Sec[x/2]`, 1337/1341); the reject
     path is unchanged. `dsolve_variation_of_parameters` reserves the symbolic-limit definite
     convolution for DiracDelta forcing and keeps a per-term INDEFINITE Wronskian integral (inert
     when non-elementary) otherwise — matching Mathematica's integral form and never entering the
     parametric DiffUnderInt escalation that blows up on `Tan`/`Sec`/arbitrary `g` (1350/1354).
  2. **Bounded-Kovacic complex-pole gate** (`dsolve_kovacic.c`). Case-1c declines a NON-REAL pole
     (a `time()` wall-clock backstop, no alarm) — the complex-conjugate pole pair of
     `(x^3+1)y''+4x y'+y==0`, whose `ds_simplify` on the complex radicals spins; a genuinely Heun
     equation with no Liouvillian solution that then falls to the Frobenius series (1392/1393). The
     forcing closure accepts an arbitrary-`g` inert-Integrate VoP particular (skips the
     un-numericizable `numeric_verify`), closing the forced Bessel operator 1350.
  3. **Exact-ODE plain-symbol first integral + SeriesData IVP fit** (`dsolve_exactode.c`,
     `dsolve_common.c`). The exact reduction uses a PLAIN symbol (not `C[n]`) for the first-integral
     constant — `C[n]` reads as a parametric function and drives the reduced integrating-factor
     quadrature `Integrate[C[2] E^(-Cos[x]),x]` into a >8 s DiffUnderInt hang; it declines BEFORE
     renaming (the rename re-evaluates and would re-trigger the hang inside the inert integral),
     falling to the Frobenius series (1384). `dsolve_fit_constants` takes `Normal[body]` of a
     SeriesData body so a series IVP fits `a[0]=C[1], a[1]=C[2]`, and the FIT_UNDECIDED fall-through
     extends to a 2nd-order IVP (a special-function solution singular at the IC point, Bessel at
     x=0, 1381, declines to the origin series).
  4. **IC-point Frobenius series** (`dsolve_frobenius.c`). `dsolve_frobenius_shifted_try` prefers
     the IVP's IC point as the expansion center when ordinary and lifts the rational-only gate
     there, so a transcendental-coefficient equation Taylor-expands about x0
     (`x^2 y''+(x+1)y'+3 Log[x] y==0, y[1]==2, y'[1]==0`, 1385).
  Residue 1 (NOT a wrong answer): `1360` `u''+u'+u^3/5==Cos[t]` — a forced Duffing oscillator with
  no closed form in Maple/Mathematica/SymPy. Anti-overfit units `t_m34_*` (`tests/test_dsolve.c`);
  `make check-c99` green; all prior corpus gates held. See the §2.2.14 block in
  `DSolve_test_status/STATUS.md`. Version 0.131 → 0.132.

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
- `[✓] PolynomialShiftSubstitution` — `y'==R(x)+g(x)(φ(x)+c y)^p` (`p` non-integer, `φ`
  polynomial in `x`, `c` constant, `R=−φ'/c`): the x-dependent-shift generalisation of
  `FirstOrderSubstitution` (whose argument is only the constant-coefficient linear form).
  The substitution `u=φ+c y` reduces it to the separable `u'==c g(x) u^p`; the first
  integral `u^(1−p)/(1−p)−∫c g dx==C[1]` is returned implicitly (branch-safe, verified by
  the implicit-function rule). The deterministic replacement for the radical
  `[F(x),G(x)]`-symmetry cases `abaco2_similar` aborts on. Detection/reduction is
  Simplify-free (plain evaluation) to avoid a radical-rationalisation loop. Runs after
  `Separable`, before the heavier substitution/parametric searches. See M20.
  `dsolve_polyshift.c`.
- `[✓] Factorable` — factor the equation as a polynomial in the highest derivative
  (`F1·F2·…==0`); recurse the cascade on each factor, union the branch bodies.
  (SymPy `factorable`.) Factors over plain-symbol substitutes under a `PolynomialQ`
  gate (raw funcapp `FactorList` hangs/misfactors); keeps only differential factors.
  Runs at the front. `dsolve_factorable.c`.
- `[✓] NthAlgebraic` — algebraic (degree ≥ 2) in the top derivative `y^(n)`: `Solve`
  for `y^(n)`, recurse each root branch (a branch free of `y` hits `Quadrature`);
  also the degenerate no-derivative case. Also clears a top-derivative-bearing
  **denominator** (`A/y'==B`, the 12000.org `x=x(y)` spelling): `Numerator[Together[·]]`
  restores a polynomial ODE that the linear normalizer then owns (M21). Runs at the front.
  (SymPy `nth_algebraic`.) `dsolve_nth_algebraic.c`.
- `[✓] AlmostLinear` — `f(x)g(y)y' + k(x)l(y) + m(x)==0`: substitution `u=∫g dy` →
  `u'+P u==Q` (integrating factor), then `l(y)==U(x)` solved for `y`. (SymPy
  `almost_linear`.) `dsolve_almostlinear.c`.
- `[✓] LinearCoefficients` — `y'==(a1 x+b1 y+c1)/(a2 x+b2 y+c2)`: det≠0 → shift to the
  lines' intersection → `Homogeneous` (explicit / implicit); det=0 (parallel) →
  substitute `v=a1 x+b1 y` → `Separable`. Explicit + implicit two entries. (SymPy
  `linear_coefficients`.) `dsolve_lincoeff.c`.
- `[✓] SeparableReduced` — `x y'/y == G(x^n y)`: `n = x r_x/(y r_y)`, substitution
  `w=x^n y` → `Separable`; implicit first integral. (SymPy `separable_reduced`.)
  `dsolve_sepreduced.c`. Gated to `F` rational-in-`y` (a trig/irrational
  dependence on `y` made the `n` Simplify diverge).
- `[✓] Linearizable` — first-order ODEs that become **linear or Bernoulli in a
  new variable `u = phi(y)`** for elementary `phi` (Log/Exp/Sin/Cos/Tan): with
  `u' = phi'(y) F(x,y)` re-expressed through `y = phi^{-1}(u)`, the RHS collapses
  to `G(x,u)` linear/Bernoulli in `u`. Linear `G` is integrated directly
  (`dsolve_linear_factor_solve`); a Bernoulli `G` recurses into the cascade under
  a re-entry guard; then `y = phi^{-1}(H)`. The Sin/Cos/Tan cases introduce a
  `Sqrt[1-u^2]`/`Sqrt[1+u^2]` from `Cos[ArcSin[u]]` etc. — a candidate is kept
  only when that radical cancels (cheap eval-then-check, full Simplify only when
  an inverse-trig remains that it can collapse, e.g. `Cos[2 ArcCos u] -> 2u^2-1`).
  Gate: `F` transcendental in `y`; a per-`phi` kernel pre-filter avoids the
  costly Simplify on inapplicable substitutions. Runs after the y-linear/
  y-Bernoulli/separable specialists and **before Exact** (whose integrating-factor
  search otherwise spins on a transcendental-in-`y` equation). Solves the Maple
  `[F(x),G(x)y+H(x)]`-symmetry log family and the trig/exp-linearizable family
  (`y'Cos[y] = Cos[x]Sin[y]^2 + Sin[y]`, `y' = e^{x-y}(e^x-e^y)`, …) with a clean
  explicit `y = phi^{-1}(H)`. `dsolve_linearizable.c`.
- `[✓] LieSymmetry` (`DSolve`LieGroup`/`LieSymmetry`) — heuristic infinitesimal
  point-symmetry method; the general first-order backstop, run after the specialists
  and before the series fallback. All ansatz heuristics implemented except the
  documented-exempt formal §4.4.2 Case I/II (`abaco2_unique_general`).
  *Implemented:* `abaco1_simple` (L1) + `linear`
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
  with x-factor X, candidates [X,−X/R] and [−R/X,1/X]; solves y' = (x/y)(x²+y²)^(1/3));
  plus the §4.4.1 order-zero extension (Eqs 73–81: non-separable [−R,1]/[1,−R]/[1,−1/R],
  catching Kamke 433). `chi` ✅ (CPC 101 1997, 5th algorithm: `η=ξω+χ`, rich
  transcendental-atom basis for `χ`; solves Kamke 357). The formal §4.4.2 Case I/II
  (`abaco2_unique_general`) is a documented exemption (impractical per the paper;
  untestable non-vacuously).  Undefined-function `ω` no longer hangs the heuristics
  (node budget + polynomial zero-test + undefined-function gate); a trig `ω` skips the
  rational/algebraic classifiers and goes straight to `chi`.
  Nine ansatz heuristics (`abaco1_simple`,
  `abaco1_product`, `function_sum`, `abaco2_similar`, `linear`, `bivariate`, `chi`,
  `abaco2_unique_unknown`, `abaco2_unique_general`); each candidate `(ξ,η)` is gated
  through the symmetry condition, then integrated by the Lie integrating factor
  `μ = 1/(η − ω ξ)` → first integral (reuses the exact-equation quadrature) →
  `dsolve_run_implicit` + explicit `Solve` inversion. (SymPy `lie_group`.) See M10.
  `dsolve_lie.c`.

### 1b. Linear constant-coefficient

*Coefficient normalization (shared by const-coeff / Euler / undetermined-coeff):*
`dsolve_linear_normalize` clears denominators (so a rational-RHS form
`y''' == (24x+24y)/x^3` becomes the polynomial-coefficient Euler `x^3 y''' - 24y
== 24x`) and divides by the polynomial GCD of the coefficients (so
`x(y'''+2y''-y'-2y) == 1` becomes the CONSTANT-coefficient `y'''+2y''-y'-2y ==
1/x`). Depth-gated to the outermost call so it does not slow `OperatorFactor`'s
recursive sub-solves.

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
  quadratic whose linear factors give radical-free roots), and **Pöschl-Teller**
  (M16, `P=0`, `Q=c0+c1 Csc²x+c2 Sec²x` — extracted as an even quadratic in `Cos x`
  after `Sin^{2k}→(1−Cos²)^k`; → `Sin^p Cos^q 2F1((p+q±√(−a))/2, p+½, Sin²x)` and the
  `p→1−p` partner, `p(p−1)=−c1`, `q(q−1)=−c2`, `a=−c0`; gated by an in-method NUMERIC
  self-verify since the 2F1 residual is undecidable — covers the Kamke/Murphy
  Csc²/Sec²/(a Cos²+b Sin²+c)/Sin² trig-potential family, and the `y''+Cot x y'+…`
  spellings once M14's t=Cos transform reaches this form). The hypergeometric
  heads auto-rewrite to `HypergeometricPFQ`, which has a `deriv.c` z-derivative
  rule, so the branches verify. The second solution carries `x^(1−b)`/`x^(1−c)`
  and is emitted only when that exponent parameter is a **non-integer number**
  (an integer makes the pair dependent / the pFq lower parameter singular; a
  symbolic exponent makes the verify residual a symbolic-power+pFq sum on which
  zero_test currently hangs). Both degenerate cases decline to the Frobenius
  series fallback; the other parameters (`a`; `a,b`) may stay symbolic.
  The canonical-singular-point restriction is **lifted (M17)** by an affine map plus
  a local-exponent shift `Y=s^{r0}(1-s)^{r1}F`: a rational-coefficient equation with
  two finite regular singular points `{x1,x2}` off `{0,1}` is mapped onto `x(1-x)` and
  solved as `Hypergeometric2F1` (Gegenbauer/Jacobi/associated Legendre at symbolic
  degree; numeric-self-verified). Ordinary Legendre keeps `LegendreP`; the associated
  case (μ≠0) declines to the affine path (3-arg `LegendreP` does not numericize). M17
  also adds a Liouville normal-form pre-pass so the P==0 Airy/Bessel recognizers fire
  on an equation with a y' term. **Confluent Whittaker / ₁F₁ (M19):** the confluent twin
  of the affine→Gauss row — a normal form `z''==r z` with ONE finite regular singular
  point (double pole of `r`) and a rank-1 irregular point at ∞ (`r → b2 ≠ 0`) is emitted
  as verifiable `Exp[−z/2] z^(1/2±μ) Hypergeometric1F1[1/2±μ−κ, 1±2μ, z]` (`z=c(x−x0)`,
  `c=2√(−b2)`, `μ=√(1/4−b0)`, `κ=b1/c`); declines `2μ ∈ ℤ`. Run on the **y'-free (P==0)
  surface only** (the P≠0 recovery-factor path stacks same-base radical powers and can
  `$IterationLimit` → future work); covers 2.1.2-102/-568 and the Whittaker/Coulomb normal
  forms. NOTE: the integer-degree **Legendre / Chebyshev / Gegenbauer / Jacobi**
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
- `[✓] ChangeOfVariable` (`DSolve\`ChangeOfVariable`) — 2nd-order linear with
  TRANSCENDENTAL coefficients `y''+P y'+Q y==0`: try a change of the independent
  variable `t=phi(x)` (`Cos`/`Sin`/`Tan`) that rationalizes the coefficients
  (`A=(phi''+P phi')/phi'²`, `B=Q/phi'²` become rational in `t`), recurse the
  cascade on the transformed equation, and back-substitute `t=phi(x)`. Flagship:
  `y''+Cot[x]y'+k(k+1)y==0 → (1−t²)Y''−2tY'+k(k+1)Y==0` (Legendre) via `t=Cos[x]`.
  NUMERICALLY verified on the original (Legendre-Q `Log` terms make the residual
  undecidable); bounded sub-solves + deadline + decline memo + re-entry guard;
  transcendental-coefficient gate. Runs after `Kovacic`, before series. See M14.
  `dsolve_changevar.c`.
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
- `[✓] SecondOrderSymmetry` (`DSolve\`SecondOrderSymmetry`) — the general **nonlinear**
  2nd-order backstop: find a Lie **point** symmetry `X=ξ∂ₓ+η∂_y` of `y''==Φ(x,y,y')`
  by a polynomial-ansatz determining system (2nd-prolongation `NullSpace`, the M10
  first-order machinery lifted to `{x,y,p}`), then reduce the order via canonical
  coordinates `(r,s)` → first-order `dq/dr==F(r,q)` → cascade → invert for `y`. Every
  inversion branch is NUMERICALLY back-substitution verified (the symbolic verify keeps
  undecidable residuals, so a wrong branch must be caught here); linear ODEs are gated
  out (their domain is Euler/Kovacic/SpecialFunction/Frobenius). All recursive
  sub-solves are `TimeConstrained`-bounded with a wall-clock deadline + decline memo so
  a decline is a clean bounded fall-through. Solves the scaling/projective/`_mu_*`
  reducible families (Kamke/Murphy nonlinear 2nd-order). Runs after `Liouville`, before
  the series fallback. See M12. `dsolve_lie2.c`.
- `[~] ReducibleIntegratingFactor` (`DSolve\`ReducibleIntegratingFactor`) — nonlinear
  2nd-order `y''==Φ(x,y,y')` admitting an integrating factor μ of a restricted form
  (Cheb-Terrab & Roche 1999; Maple `_reducible,_mu_*`): `μ=R_{y'}`, `R=∫μ dy'+G` with `G`
  from `R_x+y'R_y+ΦR_{y'}==0` (2.9–2.10), then `R==C[1]` → first-order cascade. Symbolic
  `A(R)==0` + numeric verify ⇒ wrong μ always declines. **Stage 1 done: μ(x,y)** (Φ deg-≤2
  poly in y'; Case A closed-form / Case B linear-ν-ODE). Stages 2/3 (μ(x,y'), μ(y,y')) +
  the 3rd-order forms are pending. Runs before `SecondOrderSymmetry` (cheaper algebraic
  search, reaches ODEs with no point symmetry). Linearity gate. See M18. `dsolve_ifactor.c`.

### 1e. Systems

Cascade order (`nfun>1`): `DecoupleSystem` → `TriangularSystem` →
`SystemReduce` → `LinearFirstOrderSystem` → `LinearSystemVarCoeff` →
`LinearSystemCommutative` → `AutonomousSystem`.

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
- `[✓] SystemReduce` — **higher-order coupled linear systems** by STATE
  AUGMENTATION. For functions of orders `m_j`, the state
  `Y = (u_1,u_1',…,u_1^(m_1-1), u_2,…)` (`|Y| = Σ m_j`) advances as
  `s_{j,k}' = s_{j,k+1}`, with the top rows `s_{j,m_j-1}' = u_j^(m_j)` obtained
  by solving the `n` original equations for the `n` top derivatives (`LinearSolve`
  on the leading matrix `L = ∂R/∂(u_j^(m_j))`). This yields `Y' == A Y + b(x)`;
  when `A` is constant it is fed straight to `dsolve_linsys_assemble` (the same
  Jordan → `e^{At}` → variation-of-parameters engine as `LinearFirstOrderSystem`,
  so defective/complex spectra and forcing come for free), then the `u_j = s_{j,0}`
  components are read back. Solves the constant-coefficient second-/higher-order
  systems (`{x''+x'+y'-2y==0, x'+x-y'==0}`, `{x''==4y, y''==4x}`, mixed orders,
  forced). Declines a **variable** `A` (variable-coefficient higher-order system)
  or a **singular** leading matrix `L` (a differential-algebraic leading form) —
  both routed to the future operator-determinant elimination fallback (P1b).
  `dsolve_sys_reduce.c`. *(Realification note: the shared `dsolve_linsys_tidy` is
  size/content-adaptive — a large or mixed exp+trig augmented body gets `Expand`
  instead of a full `Simplify`, which otherwise spends tens of seconds in
  `Together`; the result is back-substitution-verified either way.)* This is the
  work `M8` promised at `:247` but never landed.
- `[✓] LinearFirstOrderSystem` — `Y' == A Y + b(x)`, constant `A`, **any** spectrum.
  Fundamental matrix `Φ = e^{Ax} = S · e^{Jx} · S^{-1}` from `JordanDecomposition`
  (diagonalizable → `C e^{λx} v`; **defective → `x^k e^{λx}`** generalized-eigenvector
  terms; complex pairs → real `e^{αx}Cos/Sin[βx]` via `ComplexExpand`). Forcing `b(x)`
  by variation of parameters `Φ·(C + ∫ Φ^{-1} b dx)` — subsumes `-A^{-1}b` and stays
  valid for singular `A`. Multi-function verify + IVP/BVP fit in the substrate
  (`dsolve_verify_system`/`dsolve_fit_system`/`dsolve_assemble_system`).
  *Was (M4): eigen-only, diagonalizable, `-A^{-1}b` forcing, defective decline.*
- `[✓] LinearSystemVarCoeff` — genuinely coupled, *non-triangular*, *variable* `A`,
  the **scalar-factor class** `A(x) == f(x) B` (`B` constant): `t = ∫f dx` reduces
  `Y' == f(x) B Y + b(x)` to `dY/dt == B Y`, so `Φ == e^{B t}` reuses the constant
  fundamental-matrix builder (`dsolve_linsys_assemble`); forcing by variation of
  parameters. Cascade slot after `LinearFirstOrderSystem` (declines a constant `A`).
  `dsolve_linsys_varcoeff.c`. *Future (still `[ ]`):* the wider commutative-
  antiderivative class (`[A,∫A]==0` but not scalar×constant) and the genuinely
  non-commuting Floquet/Magnus case — both need a symbolic `MatrixExp` of a variable
  matrix. (SymPy's non-constant `linear_neq_order1`.)
- `[✓] AutonomousSystem` — 2-D **autonomous** systems `x'==f(x,y)`, `y'==g(x,y)`
  (`f,g` free of `t`), including **nonlinear**, by the phase-plane reduction:
  eliminate `t` via the orbit ODE `dy/dx == g/f` (a scalar first-order ODE the
  cascade solves), then reconstruct `x(t)` from `x'==f(x,Y(x))` along the orbit;
  `y(t)=Y(x(t))`. The orbit constant is renumbered to `C[2]` so the
  reconstruction's fresh `C[1]` cannot collide. Solves `{x'=y, y'=y^2/x}` (orbit
  `y=Cx`), `{x'=-1/y, y'=1/x}` (`xy=C`), `{x'=x/y, y'=y/x}` (`1/x-1/y=C`),
  `{x'=1/y, y'=1/x}`. Restricted to a **rational** orbit + field: a radical orbit
  (e.g. `Sqrt[x^2+C]` from `{x'=y/(x-y), y'=x/(x-y)}`) gives a nonelementary
  reconstruction and is declined (that path also currently exposes a pre-existing
  NULL deref in the radical Risch–Norman `risch_squarefree_t`). Last in the system
  cascade. `dsolve_autosys.c`. *Future:* implicit-orbit / radical-orbit
  reconstruction once the elliptic quadratures are supported.
  Also: `dsolve_linsys_extract_Ab` now admits a **t-dependent leading-derivative
  coefficient** (`t x' + y == 0`), dividing through to a variable `A` — so the
  scalar-factor class is reached from the natural `t x' = …` spelling, not only
  the pre-divided `x' = …/t` one.
- `[✓] LinearSystemCommutative` — the **2x2 commutative class**
  `A(t) == a(t) I + b(t) K0` (`K0` constant), i.e. `[A, ∫A] == 0`, where the
  fundamental matrix is the ordinary exponential
  `Phi = Exp[∫a] · Exp[(∫b) K0]` with, for constant traceless `K0` and
  `mu^2 = -det K0`, `Exp[s K0] = Cosh[mu s] I + (Sinh[mu s]/mu) K0` (nilpotent
  `mu=0`: `I + s K0`; complex `mu` realifies to the rotation Cos/Sin). Decomposed
  by trace (`a`) and one scalar factor for the traceless remainder (`b K0`, `K0`
  constant), which is robust where a symbolic `Eigenvectors[A]` returns a
  non-constant `Sign[t]`-scaled basis. Covers the rotation/decay families the
  scalar-factor reduction misses (`{x'=-x+t y, y'=t x-y}`, `{x'=x Cos t-y Sin t,
  y'=x Sin t+y Cos t}`, `{x'=x/t+y, y'=-x+y/t}`, `(t^2+1)`-scaled). Forcing by
  variation of parameters; after `LinearSystemVarCoeff` in the cascade; declines
  constant or non-commutative `A`. `dsolve_linsys_commutative.c`. *Future:* the
  n×n commutative case and the genuinely non-commuting Floquet/Magnus case.

### 1f. Conditions
- `[✓]` IVP (fit constants at one point) — in the substrate.
- `[✓]` Linear BVP (multiple points) — in the substrate + `Solve`; **sound**: an
  inconsistent / over-determined BVP now returns `{}` (no solution) rather than the
  silent unfitted general solution (the constant-fitters signal an empty-`Solve`
  result up through `dsolve_run` / `dsolve_run_system` as the concrete empty list,
  distinct from a decline). Under-determined BVPs keep the free constant.
- `[✓] EigenvalueProblem` — Sturm–Liouville `y'' + λ y == 0` on `[a,b]` with two
  homogeneous BCs (Dirichlet/Neumann/mixed), first cut: eigenvalue family + eigen-
  functions, verified under `C[1] ∈ Integers`; pinned-only `DSolve`EigenvalueProblem`
  (`dsolve_eigenvalue.c`). *Future:* non-constant weight, Robin/periodic BCs, the
  `λ==0` Neumann mode, and the `DEigensystem`/`DEigenvalues` surfaces.

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
- `[✓] PDEQuasilinear` (Lagrange) — `P u_{v1}+Q u_{v2}==R`, linear in the first
  derivatives (P,Q,R may depend on u), by the method of characteristics
  `dv1/P=dv2/Q=du/R`. First cut, two bounded classes: **semilinear** (P,Q free of
  u; base characteristic `dv2/dv1=Q/P` decouples → first integral ξ via recursion
  into the scalar ODE cascade; linear u-ODE along the characteristic by the
  integrating factor, with `Solve[ξ==const, v2]` for the y-along-characteristic
  when the coefficients need it) → **explicit** `u==body` with C[1][ξ]
  (`x u_x+y u_y==u → x C[1][y/x]`), generalizing `PDELinearFirstOrder` to variable
  coefficients; **conservation law** (R==0, coefficients touching u; u invariant,
  φ2=u, φ1 from `Q dv1-P dv2=0` with u a parameter) → **implicit**
  `φ1(v1,v2,u)==C[1][u]` (inviscid Burgers `u u_x+u_y==0`). Declines coefficients
  depending on u with R≠0 (full characteristic system: future) and non-elementary
  characteristic integrals. New substrate `dsolve_run_pde_implicit`
  (`PDEImplicit`/`PDEExplicit`/`PDEBranches` wrappers; implicit relation verified by
  the implicit-function rule with concrete test functions). `dsolve_pdequasi.c`.
- `[✓] PDECharpit` (nonlinear complete integral) — `F(v1,v2,u,p,q)==0` (p=u_{v1},
  q=u_{v2}) by Charpit's method, the three classic **standard forms**: `F(p,q)`
  (explicit `u=C[1] v1 + q v2 + C[2]`, q from `F(C[1],q)==0`); `F(u,p,q)` (implicit
  `∫du/P == v1 + C[1] v2 + C[2]`, P from `q==C[1] p`); separable `f(v1,p)==g(v2,q)`
  (explicit `u=∫P dv1 + ∫Q dv2 + C[2]`). Explicit forms verify through the ordinary
  explicit PDE path (bare constants survive); the implicit `F(u,p,q)` form via the
  new `PDERelation` wrapper (implicit-function-rule verify). Runs after
  quasilinear/Clairaut (nonlinear-in-derivatives gate). Declines equations in all
  of v1,v2,u (the general integrable-combination search + singular solutions are
  future) and non-elementary integrals. `dsolve_pdecharpit.c`.
- `[✓] PDEClairaut` — `u==v1 u_{v1}+v2 u_{v2}+f(u_{v1},u_{v2})`; complete integral
  `u==C[1] v1+C[2] v2+f(C[1],C[2])`, plus the singular envelope (eliminate the
  constants from `{v1+f_{C1}==0, v2+f_{C2}==0}`) under `IncludeSingularSolutions`.
  Nonlinear in the derivatives (linear/quasilinear methods decline; a degenerate
  linear-f Clairaut is semilinear, so `PDEQuasilinear` runs first for the general
  arbitrary-function solution). `dsolve_pdeclairaut.c`.

### 2b. Second order
- `[✓] PDELinearSecondOrder` (the `PDEHyperbolicGeneral` item, generalized) —
  homogeneous, principal-part-only, constant-coefficient 2nd-order linear PDE
  `A u_{v1 v1}+B u_{v1 v2}+C u_{v2 v2}==0` by operator factoring: the trial
  `u==f(v2+λ v1)` gives the characteristic quadratic `A λ²+B λ+C==0`, so the
  principal operator factors over ℂ and one method covers all three discriminant
  signs — distinct real roots (hyperbolic) → `C[1][v2+λ1 v1]+C[2][v2+λ2 v1]` (the
  wave equation `u_tt==c² u_xx → C[1][x-c t]+C[2][x+c t]`, d'Alembert); complex
  roots (elliptic) → the complex-characteristic form (Laplace `u_xx+u_yy==0 →
  C[1][y-I x]+C[2][y+I x]`, matching Mathematica — no realification needed); a
  repeated root (parabolic) → `C[1][w]+v1 C[2][w]`, `w=v2+λ v1`. Reuses
  `dsolve_analyze_roots` (distinct/complex/repeated λ uniformly); back-substitution
  verified. `A==0` handled by swapping `v1↔v2`; pure mixed `B u_{v1 v2}==0 →
  C[1][v1]+C[2][v2]`. The shared `dsolve_verify_pde` was generalized to arbitrary
  order (scanned from the residual — `max_order` is 0 for PDEs) and multiple
  arbitrary functions (`C[1..4]` → distinct test functions). **Lower-order terms**
  `D u_{v1}+E u_{v2}+F u` are now handled when the FULL symbol
  `A ξ²+B ξη+C η²+D ξ+E η+F` factors into two first-order operators
  `(∂_v1−λ_i ∂_v2+m_i)` — the shifts `m_i` solve `m1+m2=D/A`, `λ2 m1+λ1 m2=−E/A`
  with the factorability check `m1 m2==F/A` (repeated λ: `E+λD==0`, then the
  `m²−(D/A)m+F/A` quadratic) — giving the exponential-damped form
  `Σ e^{−m_i v1} C[i][v2+λ_i v1]`. Covers the **distortionless telegraph**
  `u_tt−c² u_xx+a u_t+(a²/4)u==0 → e^{−a t/2}(C[1][x+c t]+C[2][x−c t])` and
  damped/convection factorable cases; `D=E=F=0` reproduces the principal-part form.
  `dsolve_pde2.c`. *Future:* inhomogeneous forcing, the non-factorable general
  telegraph (`b≠a²/4` → Bessel functions).
- `[✓] SeparationOfVariables` — separated product solution `u == X(v1) Y(v2)` of a
  homogeneous, constant-coefficient linear PDE with **no mixed derivative term**.
  Dividing by `X Y` separates into two constant-coefficient ODEs in a separation
  constant λ — `Σ a_i X^(i) − λ X == 0`, `Σ b_j Y^(j) + (e+λ) Y == 0` — each solved
  by recursing into the scalar cascade; λ becomes a generated constant. The single
  redundant overall scale is absorbed (fixing a first-order side's lone constant to
  1 loses no generality). Example: the heat equation `u_t == u_xx →
  E^(λ t)(C[1] E^(−√λ x) + C[2] E^(√λ x))`; also Helmholtz `u_xx+u_yy+u==0`. Returns
  the *representative product mode* (the general solution is a superposition over λ),
  so **pinned-only** (not in the automatic cascade — matching `FirstOrderPowerSeries`
  / `EigenvalueProblem`). Back-substitution verified. `dsolve_pdesep.c`. *Future:*
  variable (product-separable) coefficients, BC-driven eigenfunction expansions.
- `[✓] PDEClassify` — `PDEClassify[eqn, u, {v1, v2}]` classifies a 2nd-order
  linear PDE by the discriminant `Δ = B² − 4 A C` of its principal part
  (`A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2}`; only the highest-order terms count):
  `"Hyperbolic"` (Δ>0, wave), `"Parabolic"` (Δ=0, heat), `"Elliptic"` (Δ<0,
  Laplace). A standalone builtin (not a solver). A discriminant whose sign is not
  a decidable constant — a mixed-type / parameter-dependent equation such as
  Tricomi's `y u_xx + u_yy == 0` (Δ = −4y) — leaves the call unevaluated (an
  honest decline, not a region-blind label). `dsolve_pdeclassify.c`.
- `[✓] WaveDAlembert` — the initial-value problem for the 1-D wave equation on
  the whole line: `u_tt == c² u_xx`, `u(x,t0) == f(x)`, `u_t(x,t0) == g(x)` →
  d'Alembert `u = ½(f(x−cτ)+f(x+cτ)) + 1/(2c)∫_{x−cτ}^{x+cτ} g(s)ds`, `τ = t−t0`
  (the velocity integral kept unevaluated for undefined `g`, dummy `K`). The two
  initial conditions are two-argument point conditions the shared parser does not
  recognise, so they arrive as ordinary equations (`neq==3`); the method sorts
  PDE-vs-ICs itself (the fixed variable is time, the free one space), and has its
  own multi-equation runner + verify (rebuild with `f==Cos, g==Sin` and require the
  PDE residual and both ICs to vanish — the unevaluated integral of an undefined
  `g` cannot be differentiated). Solved automatically by `DSolve` and via the
  pinned `DSolve`WaveDAlembert`. `dsolve_wave.c`. *Future:* inhomogeneous forcing,
  half-line / boundary problems, `Piecewise` data.
- `[✓] HeatKernel` — the Cauchy (initial-value) problem for the 1-D heat equation
  on the whole line: `u_t == k u_xx`, `u(x,t0) == f(x)` → the heat-kernel
  convolution `u = 1/Sqrt[4πkτ] ∫_{-∞}^{∞} f(K) Exp[-(x−K)²/(4kτ)] dK`, `τ = t−t0`.
  The single initial condition arrives as an ordinary equation (`neq==2`); the
  method sorts PDE-vs-IC itself (shared with the wave IVP's approach). The Gaussian
  convolution is nonelementary, so the integral is returned unevaluated (matching
  Mathematica for general `f`) and built WITHOUT attempting integration (the
  `Function` body holds it, avoiding a hang on the improper integral). Verified at
  the KERNEL level — the heat kernel `G` is a decidable Gaussian solving
  `G_t == k G_xx` (so the convolution does, by differentiation under the integral),
  which also checks `k` was read correctly; the initial condition holds by the
  kernel's convergence to `δ(x−K)` as `τ→0+` (a distributional limit, not a
  back-substitution). Solved automatically and via pinned `DSolve`HeatKernel`.
  Requires `k>0` (forward diffusion); declines backward-heat, second-order-in-time,
  advection (`u_x`), and reaction (`u`) terms. `dsolve_heat.c`. *Future:* the
  Erf-producing step/box data (where the integral evaluates), advection–diffusion,
  reaction, and finite-interval Fourier-series problems.

## Testing

- **In-engine self-verification:** every returned branch back-substitutes to a
  residual that is not decidably non-zero (`dsolve_run` → `zero_test_decide`).
- **Unit tests** (`tests/test_dsolve.c`): the Wolfram reference examples, checked
  by `PossibleZeroQ` of the residual / of the difference from the known form.
  Every method — scalar, system, and PDE — is now a pinned backtick builtin
  (`DSolve`<Method>[...]`) with a docstring and a pinned-method test (systems/PDE via
  `dsolve_method_builtin_system` / `_pde`); the automatic dispatch is exercised too.
- **Stress tests:** parametrized forward-generator families per method group,
  back-substitution verified, each guarding `Head === List` first so a declined
  solve cannot pass vacuously. `tests/test_dsolve_m5_stress.c` covers Kovacic and
  Frobenius/PowerSeries; `tests/test_dsolve_m12_stress.c` covers the nonlinear
  2nd-order `SecondOrderSymmetry` families (verified NUMERICALLY — the solutions
  carry logs/radicals PossibleZeroQ cannot decide, and lie2 itself relies on a
  numeric back-substitution guard); `tests/test_dsolve_m14_stress.c` covers the
  `ChangeOfVariable` Legendre families (verified via vanishing C[1]/C[2] residual
  coefficients); `tests/test_dsolve_stress.c` covers the rest of
  the cascade (LinearFirstOrder, Separable, Bernoulli, Homogeneous, Exact,
  LinearConstantCoefficients, EulerCauchy, ReductionOfOrder, 2×2 + triangular
  systems, first-order PDE). A generator builds the equation from parameters whose
  closed form is guaranteed — a chosen spectrum (ConstCoeff/Euler), a potential
  (Exact), or a target solution `yt` with `q = yt' + p·yt` (LinearFirstOrder, so
  the integrating-factor integral is elementary by construction).
- **Gates:** `dsolve_tests` (ctest), `valgrind --leak-check=full`, `make
  check-c99`, and REPL spot-checks (`DSolve[...]`, `?DSolve`, `?DSolve`Separable`).

DSolve is a symbolic/structural head (no element-wise numeric mapping); it is a
documented exemption from the packed-aware / Compile surfaces.

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
  - **Stage 2 — μ(x,y') (Section 2.2, Lemma 3).** ✅ DONE in **M56** (below). `μ = 𝓕(x,y')·μ̃(x)`
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
    picks the right case. **RESOLVED in M56 via path (b):** Cases A/C/D + Lemma 2 find valid
    μ (Kamke 226 → μ=y'; Kamke 136 → (y'−1)/h(y'); Kamke 66 → (y'+b)/(a(1+y'²)^{3/2}) — all
    with `A(R)=0` holding), but the reduced first integrals `R==C[1]` are NON-ELEMENTARY
    first-order ODEs (`y'=√(x²y²+2C)`, `y'=Tan[C+Log[x−y]]`) that no CAS closes in elementary
    explicit form (they return them implicitly). M18 had therefore reverted the search for lack
    of *explicit* yield; **M56 keeps the search and emits the first integral
    `R(x,y[x],y'[x])==C[1]` as a reduction-of-order answer** (the new `dsolve_run_first_integral`
    runner + `DSolve\`ReducibleFirstIntegral`), as chini/abel/homogeneous-implicit already do for
    first-order ODEs. Cases E/F (`𝓗'=0` exponential; the general p'(x)-elimination) and the
    general (non-degenerate) Case D remain future.
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

- **M35 — §2.2.15 corpus (Problems 1401–1500, Boyce & DiPrima) + `DSolve\`PiecewiseForcing`
  (step / piecewise / Heaviside forcing).** ✅ DONE. The Boyce & DiPrima Laplace-transform
  chapter: 100 records — **39 scalar (18 IVP), 61 systems** (53 2×2 + 8 3×3 constant-coefficient
  linear). Scalars are high-order constant-coefficient linear (1462–1489) plus **step-forced
  2nd-order IVPs** (1492–1500). **98/100, 0 FAIL, 0 regression.** The 1492–1500 family was a set
  of FALSE passes at baseline (DSolve returned an inert `Integrate[UnitStep[…]·Cos,t]` the prelude
  could not numericize → UNK → trusted); they now GENUINELY solve.
  - **`DSolve\`PiecewiseForcing`** (`dsolve_piecewise.c`, new method + `ATTR_PROTECTED` builtin).
    Mathilda has `Piecewise`/`UnitStep` but no `LaplaceTransform`, so step-forced linear IVPs are
    solved by **interval continuation**: normalize the forcing's finite breakpoints, solve on each
    interval (recursing the scalar cascade), match the `C^(n-1)` continuity data `y,…,y^(n-1)`
    across each breakpoint, and assemble a verified `Piecewise` closed form (Mathematica's own
    output form). Routes through the standard `dsolve_run` (the constant-free `Piecewise` body
    verifies; fit is a no-op), with its OWN per-interval + per-IC numeric guard `pw_num_ok` (since
    `dsolve_run`'s probe samples only the first interval). Cascade slot: BEFORE
    `UndeterminedCoefficients`/`LinearConstantCoefficients` (whose VoP leaves the inert integral),
    gated to step/piecewise-forced linear IVPs so smooth forcing falls through; bounded exactly as
    M12/M14 (re-entry guard + `TimeConstrained` sub-solves + wall-clock deadline + decline memo).
    Solves e.g. `y''+4y==Piecewise[{{1,0<=t<Pi}},0], y(0)=1,y'(0)=0` →
    `Piecewise[{{1/4+3/4Cos[2t],t<Pi},{Cos[2t],t>=Pi}},0]`.
  - **`dsolve_verify_body` distributional-keep** (`dsolve_common.c`). A residual carrying
    `UnitStep`/`Piecewise` joins `DiracDelta`/`HeavisideTheta`/definite-`Integrate` as
    accepted-on-construction: the numeric probe cannot run on it (differentiating a `Piecewise`
    yields a boundary term the probe reads as an undefined function → it bails) and `zero_test`
    evaluates a piecewise condition on the wrong branch — it spuriously rejected the correct
    step-forced answer of 1497 (`y''+4y==Sin[t]−UnitStep[t−2π]Sin[t]`). The producing method
    self-verifies every interval and IC.
  - **Converter** (`tools/latex_ode_to_mathilda.py`): the `\left\{…cases…\right.` environment now
    converts to `Piecewise[{{v,c},…}, default]` (`otherwise`→default, `±∞` bounds dropped — a
    comparison against `Infinity` does not reduce to `True`, which would leave the tail clause
    unevaluated inside a residual); `\le`/`\ge`/`\infty` map to `<=`/`>=`/`Infinity`; Maple
    `Heaviside` → Mathilda `UnitStep`. Symbol detection reads the UNprotected tex so the forcing's
    independent variable is not hidden.
  - **Residue (2, NOT wrong answers — bounded declines, both also ✗ in SymPy):** `1463`
    `t(t−1)y''''+E^t y''+4t²y==0` (4th-order transcendental) and `1469` `t y'''+2y''−y'+t y==0`
    (3rd-order variable-coefficient — the `3rd_high_linear` operator-factoring bucket, future work).
  - Anti-overfit unit test `t_m35_piecewise_forcing` (`tests/test_dsolve.c`); `make check-c99`
    green; §2.1.2/§2.2.1–§2.2.14 corpus gates all held. See the §2.2.15 block in
    `DSolve_test_status/STATUS.md`. Version 0.132 → 0.133.

- **M36 — §2.2.16 corpus (Problems 1501–1600, Boyce & DiPrima) + cubic-log separable → implicit
  first integral.** ✅ DONE. 100 records, **all scalar (56 IVP), 0 systems**: 77 elementary
  first-order (39 separable, 25 linear, 13 quadrature), 18 step/impulse-forced linear IVPs
  (1501–1518: `UnitStep`/Heaviside → `PiecewiseForcing`, `DiracDelta` → variation of parameters,
  one mixed `DiracDelta`+`UnitStep`), 5 specials (Clairaut 1536, Riccati 1577, first-order symmetry
  1575/1576, class-A 1561). **Baseline 96/100 → 97/100, 0 FAIL, 0 regression.** Two increments:
  - **Converter `Abs[…]`** (`tools/latex_ode_to_mathilda.py`). `\left|…\right|` / `\lvert…\rvert` /
    `\vert…\vert` / bare `|…|` now convert to `Abs[…]` right after the `\left`/`\right` strip and
    before the mains/`{}`→`()` passes (so the wrapped body still picks up its `[x]`). Without it the
    `|` bars in 1535 `y'==|y|+1` survived into `(|y[x]|)` and the WHOLE corpus file failed to parse
    as a `List` literal (blocking the section entirely).
  - **Cubic-log separable → implicit first integral** (`src/calculus/dsolve_separable.c`). A
    separable ODE whose y-side antiderivative `Integrate[1/h, y]` carries **≥3 distinct `Log[…]`
    arguments** (the partial-fraction integral of a cubic-or-higher `h(y)`, e.g.
    `−½Log[y−1]+⅙Log[y+1]+⅓Log[y−2]`) has no elementary explicit inverse for `y` that `ds_solve`
    closes — it churns, and the implicit twin (next in the cascade) never runs. New
    `sep_noninvertible_logsum` (a recursive distinct-`Log`-argument counter, robust to whether the
    coefficient sits outside the `Log` or is absorbed as a `Log[(y−1)^(−1/2)]` power) makes the
    explicit path decline, so the implicit first-integral twin returns `G(x,y)==C[1]`, verified by
    the implicit-function rule and fitted on an IVP as `C=G(x0,y0)`. Closes **1590**
    `y'+(1+y)(y−1)(y−2)/(x+1)==0, y(1)=0` — which **SymPy does not solve**. A one- or two-log
    relation is untouched (single `Log` inverts by `Exp`; two logs merge to a Möbius
    `Log[(y−a)/(y−b)]` → `Tanh`/logistic), so no separable regresses (§2.2.8/§2.2.12 held).
  - *Residue (3, NOT wrong answers — bounded declines):* `1508`
    `y''+3y'+2y==DiracDelta[t−5]+UnitStep[t−10]` and `1509` `y''+2y'+3y==Sin[t]+DiracDelta[t−3π]` —
    the impulse Green's-function convolution churns on the irrational-frequency `−1±I√2` kernel; its
    `TrigReduce` linearisation is REQUIRED for the `DiracDelta` sift to be correct (skipping it
    produces a spurious ramp — verified on 1506), so it cannot simply be gated. `1534`
    `y'==a y^((a−1)/a)` — the converter reads the parameter `a` as the independent variable (no
    explicit indvar present); with that reading the equation is genuinely non-elementary, and
    forcing `x` yields only a messy implicit form the prelude cannot verify. (All three would need a
    substantial substrate rewrite — a direct Green's-kernel impulse sift, a dispatcher-level forcing
    superposition, or an indep-var heuristic change risking all 15 prior corpora — so deferred.)
  - Anti-overfit unit test `t_m36_separable_cubic_log` (`tests/test_dsolve.c`); `make check-c99`
    green; §2.1.2/§2.2.1–§2.2.15 corpus gates all held. See the §2.2.16 block in
    `DSolve_test_status/STATUS.md`. Version 0.133 → 0.134.

- **M37 — §2.2.17 corpus (Problems 1601–1700, Nasser Abbasi) + fractional-power branch
  correctness fix + `DSolve\`AbelAIR` (Abel-2nd-kind scaling reduction).** ✅ DONE. 100 records,
  **all scalar (28 IVP), 0 systems**: first-order-dominated (27 homogeneous, 23 Abel-tagged —
  most overlap `_homogeneous`/`_exact` and already solve — 14 separable, 8 Bernoulli, 8 exact,
  6 quadrature, 4 linear, 2 Riccati, 8 `y=_G(x,y')` no-method). **Baseline 79/100 with 5 FAIL →
  87/100, 0 FAIL, 0 regression.** Two increments:
  - **Fractional-power branch correctness (the 5 FAILs).** Five `√y`/`y^(1/3)`/`y^(3/2)` ODEs
    shipped a demonstrably wrong explicit branch (e.g. `y'−2y==2√y, y(0)=1` → spurious `y=1`
    instead of `(2Eˣ−1)²`; `y'=3x(y−1)^(1/3)` → a spurious complex twin beside the correct
    `1+x³`). `dsolve_run`'s symbolic verify keeps an undecidable branch-cut residual (Solve
    policy) and the corpus prelude's numeric check then flags it. Fixed with a shared numeric
    back-substitution filter matching the prelude's verdict, at three points
    (`src/calculus/dsolve_common.c`, `dsolve_bernoulli.c`, `dsolve_separable.c`): the IVP
    constant-fitter (`dsolve_fit_constants`) tries all `Solve` roots and keeps the one that
    numerically satisfies the ODE; Separable/Bernoulli drop a robustly-nonzero `±`/principal-root
    branch (`ds_branch_num_ok`, gated by a cheap `ds_has_radical_power` pre-check); and a
    prelude-matching post-fit gate (`ds_branch_corpus_verifiable`) in `dsolve_run` drops any
    first-order explicit branch the harness would score BAD, letting the cascade fall through to a
    verifiable (often implicit) form. The gate reuses the prelude's grid+majority, so it removes
    only would-be-FAIL branches, never a would-PASS one.
  - **`DSolve\`AbelAIR`** (`src/calculus/dsolve_abel_air.c`) — a down payment on the deferred AIR
    method (M13). Solves the Abel-2nd-kind / rational-in-y forms `y' == N/D` with `D` a single
    linear y-factor `c(x)(P1 y+P0)^k` (`k∈{1,2}`) — the shapes Riccati/Chini/Abel-1st reject — by
    the scaling `u = y/s` (`s = −P0/P1`) that makes the equation separable; returns the implicit
    first integral `∫1/B du|_{u→y/s} − ∫A dx`, verified by the implicit-function rule with the
    separation checked exactly before integrating. Closes 1604/1606/1607/1676. The full
    rational-invariant solvable-class table stays deferred. Cascade after Chini/Abel-1st, before
    the substitution reductions and the Lie backstop.
  - *Residue (13, bounded declines, no wrong answers):* the eight `y=_G(x,y')` nonelementary
    equations (SymPy also fails), `1624` (cube-root real-branch), `1601` (garbled `Solve`
    inversion), `1673` (Riccati), `1681`/`1689` (nonelementary integrating factor), `1691`.
  - Anti-overfit unit tests `t_m37_fractional_power_branch` / `t_m37_abel_air`
    (`tests/test_dsolve.c`); `make check-c99` green; §2.1.2/§2.2.1–§2.2.16 corpus gates all held.
    See the §2.2.17 block in `DSolve_test_status/STATUS.md`. Version 0.134 → 0.135.

- **M38 — §2.2.18 corpus (Problems 1701–1800, Nasser Abbasi) + `zero_test` decay false-positive →
  dropped-forcing correctness fix.** ✅ DONE. 100 records, **all scalar (14 IVP), 0 systems**; 45
  first-order / 55 second-order (no order ≥ 3), mixed classes (25 2nd-order `_with_linear_symmetries`,
  17 2nd-order linear exact/(non)homog, 10 separable, 8 2nd-order `_missing_x`, 11 Abel-2nd-kind, 8
  quadrature, 5 Emden–Fowler linear/Euler subtype, plus Bernoulli/Riccati/linear/homogeneous). The
  second-order stack (Kovacic/NormalForm/SpecialFunctionForm/change-of-variable/Frobenius) already
  solved the bulk: **baseline 95/100, 0 FAIL → 96/100, 0 FAIL, 0 regression.** One increment:
  - **`zero_test` decay false-positive → dropped forcing** (`src/calculus/dsolve_common.{c,h}`,
    `dsolve_kovacic.c`). `PossibleZeroQ`/`zero_test` return **True** for a genuinely-nonzero *decaying*
    expression (documented sampler limitation; `PossibleZeroQ[E^(-2x²)] === True`). Two second-order
    gates were corrupted by it: (1) the Wronskian nonzero-check in `dsolve_variation_of_parameters`
    (a decaying `W = E^(-∫P)` read as a degenerate basis → VoP declined), and (2) `DSolve\`Kovacic`'s
    forcing detection (a decaying forcing read as homogeneous → the particular dropped and a
    **homogeneous-only WRONG answer** shipped, masked as UNEVAL by the corpus prelude's leaked-`C[k]`→
    UNFIT leniency). Fix: new shared `ds_is_structural_zero` (`Expand[·]===0`) gates both — strictly
    safer (a truly-zero forcing → `yp==0`; a truly-dependent basis is still structurally 0) — with the
    sensitive `zero_test` sampler left untouched (per the documented "don't gate on `PossibleZeroQ` for
    decaying exprs"). `1763` `y''+4x y'+(4x²+2)y==8E^(-x(x+2))` now solves as
    `(C[1]+C[2]x)E^(-x²)+2E^(-x²-2x)`. **Bonus (0 regression):** the same fix raised §2.1.2 (−5),
    §2.2.1 (−1), §2.2.2 (−2), §2.2.6 (−2), §2.2.7 (−1) — forced equations with decaying Wronskians that
    previously wrongly declined now solve.
  - *Residue (4, bounded declines, all `sympySolved=False`, no wrong answers):* `1708`/`1709`
    (Abel-2nd-kind class B, M13-deferred AIR class), `1729` (`∫Sin[x]/(b Cos x−x Sin x)` non-elementary
    → no elementary integrating factor), `1769` (VoP needs `∫E^x/√x→√π Erfi[√x]`; the Erf/Erfi Risch
    tower handles only `∫x^(1/2)E^x`, not negative half-integer powers — a deep-Risch extension, future).
  - Anti-overfit unit test `t_m38_forced_decay_wronskian` (`tests/test_dsolve.c`); `make check-c99`
    green; §2.1.2/§2.2.1–§2.2.17 corpus gates all held (several improved). See the §2.2.18 block in
    `DSolve_test_status/STATUS.md`. Version 0.135 → 0.136.

- **M39 — §2.2.19 corpus (Problems 1801–1900, Nasser Abbasi) + `DSolve\`VariationOfParameters`
  (transcendental nonhomogeneous backstop).** ✅ DONE. 100 records, **all scalar (27 IVP), 0 systems**;
  second-order-linear dominated (3 Riccati, ~30 nonhomogeneous 2nd-order linear with transcendental
  forcing, a large `_with_linear_symmetries` variable-coefficient homogeneous family, Emden–Fowler,
  Gegenbauer). The 2nd-order stack (Kovacic/SpecialFunctionForm/change-of-variable/Frobenius + M38)
  already solved the bulk: **baseline 94/100, 0 FAIL → 95/100, 0 FAIL, 0 regression.** One increment:
  - **`DSolve\`VariationOfParameters`** (`dsolve_nonhomog_vop.c`) — the general nonhomogeneous
    backstop for `L[y] == g(x)` (`g ≢ 0`) whose homogeneous part is solved by a method that does not
    itself carry forcing (chiefly `ChangeOfVariable`). Recurses `DSolve` on `L[y]==0`, extracts and
    normalises the fundamental set (`PowerExpand[Simplify[D[hom, C[k]]]]` — Simplify canonicalises
    `1+Tan²→Sec²`, PowerExpand then reduces the residual radicals `1/Sqrt[Sec²x]→Cos x` /
    `Sqrt[1/(Pi x)]·x→Sqrt[x]/Sqrt[Pi]` so the VoP integrals close), builds the particular via
    `dsolve_variation_of_parameters`, returns `hom + yp` (numerically self-verified by
    `ds_branch_num_ok`; rejected on an inert `Integrate`). **Gated to transcendental (trigonometric)
    coefficients** — a rational-coefficient nonhomogeneous equation is Kovacic's/Euler's domain (each
    with its own forcing closure), so the gate keeps the class the backstop uniquely reaches while
    sparing every rational case a redundant recursive re-solve on the large corpus. Cascade slot:
    after `SpecialFunctionForm`/`Kovacic`/`ChangeOfVariable`, so it is a pure backstop that never
    preempts their results (hence 0-regression by construction). Deliberately uses **no nested
    `TimeConstrained`** (the documented no-nest hazard aborts the whole subtree); bounded by a
    `time()`-deadline + decline memo + re-entry guard. `SpecialFunctionForm`'s depth-1-gated
    normal-form pre-pass is allowed to fire during this method's recursion (exposed via
    `dsolve_nh_vop_active`), since the backstop extracts the basis rather than composing the `μ`
    recovery factor back into a reduction. Solves `1822`
    `Sin[x] y'' + (2Sin−Cos) y' + (Sin−Cos) y == E^-x → −E^-x Sin[x] + {E^-x, E^-x Cos[x]}`.
    Anti-overfit `t_m39_nonhomog_vop` (`tests/test_dsolve.c`): 1822, a different forcing on the same
    transcendental operator, and the rational-coefficient decline.
  - *Residue (5, bounded declines, no wrong answers):* `1817` (non-elementary Struve/Lommel
    particular over a Bessel homogeneous set), `1823` (elementary `−√x/2` but Kovacic's constant-`r`
    path churns before the transcendental-gated backstop is reached — a Kovacic-robustness fix,
    future), `1836` (Solve leaves the `ExpIntegralEi`-constant IVP fit system unevaluated), `1876`
    (Kovacic `E^(nLog)` IC-fit branch artifacts — a form-normalisation gap), `1900` (exact →
    non-elementary integrating factor; complex-singular ₂F₁, future).
  - `make check-c99` green; §2.1.2/§2.2.1–§2.2.18 corpus gates all held (0 regression — the backstop
    only fires on prior-declined transcendental cases). See the §2.2.19 block in
    `DSolve_test_status/STATUS.md`. Version 0.136 → 0.137.

- **M40 — §2.2.20 (Problems 1901–2000) + §2.2.21 (Problems 2001–2100) corpus baselines + converter
  LaTeXML migration.** ✅ DONE. Two 100-ODE series-solution sections, each **96/100, 0 FAIL, 0 crash**,
  all homogeneous 2nd-order linear with a regular singular point ("series expansion around x0"; §2.2.20
  has 34 IVPs at x0 ∈ [-4,3], §2.2.21 all general at x0=0). Same territory as the series-heavy §2.2.5
  (99/100); the 2nd-order stack (Kovacic / SpecialFunctionForm / Frobenius) solves the bulk out of the
  box, and a verified truncated Frobenius `SeriesData` scores PASS by back-substitution.
  - **Converter LaTeXML migration** (`tools/latex_ode_to_mathilda.py`). The 12000.org site moved its
    HTML generator **tex4ht → LaTeXML ("oxide")**; the old `indexsubsectionN.htm` URLs the earlier
    corpora were built from are now 1.4 KB redirect stubs, and the current source is the
    `Ch2.S2.SSN.htm` "sorted sequentially" pages (LaTeXML). The converter now auto-detects the format
    (LaTeXML = no tex4ht `id='TBL-'` cells) and adds `_parse_table_latexml`: select the ODE `<table>`,
    read the column layout from its header `<tr>` legend (`# | ODE | classification | Solved? | Maple |
    Mma | Sympy | time`) over `ltx_td` cells, and pull each equation from the ODE cell's first
    `<math alttext>` (a multi-row `\begin{array}` — row 1 the ODE, rows 2+ the ICs — read with `re.S`).
    Two core fixes: `strip_array` now consumes the `\begin{array}[]` optional arg, and a new
    `_implicit_mult` pass inserts spaces at implicit-multiplication boundaries (LaTeXML juxtaposes
    products with no delimiter — `8y`, `yx`, `3y'x` — where tex4ht was space-delimited; a strict no-op
    on the old format, macro-protected so `\prime`/`\operatorname{…}`/`\mathrm{…}` survive). The
    `convert_row` LaTeX→Mathilda core and its IVP IC-splitting are unchanged; both sections convert
    100/100 and round-trip through the parser.
  - **Kovacic churn — investigated, no safe fix (bounded declines kept).** The 8 non-PASS across the
    two sections (§2.2.20: `1941`/`1964`/`1976` + symbolic-Heun `1916`; §2.2.21: `2003`/`2005`/`2006`/
    `2080`) are non-Liouvillian regular-singular ODEs whose `DSolve\`Kovacic` Case-1 Riccati solve
    (`ds_solve` on `ω'+ω²==r`, cleared against the `(leading)²` denominator into a heavy high-degree
    system) runs 9–18 s and then declines — overrunning the corpus harness's 8 s `TimeConstrained`
    before the Frobenius fallback (correct series, <0.1 s) runs. This is the `1823` "Kovacic
    constant-`r` churn" class flagged in M39. A bounded fix was attempted three ways and none is safe:
    (a) a nested `TimeConstrained` around the solve — the documented no-nest hazard, which standalone
    blew up to >90 s; (b)/(c) two structural pre-gates dropping the churning complex-pole ansatz term
    (`df≥2`, then `df≥2 && order≥2`, then coexisting-high-real-pole). **No cheap structural discriminant
    separates a fruitless churn from a real solve** — the churn *is* the Liouvillian-existence
    decision, and every gate that speeds up `2003`/`2005`/`2006`/`1941` also drops the order-2
    quadratic-pole term that `2.2.19-1822`'s `t=Tan[x]` transform genuinely needs, regressing that
    documented solve and its pinned `t_m39_nonhomog_vop`. Reverted to 0-regression clean main. A real
    fix needs a fail-fast Riccati coefficient solver (or an up-front Liouvillian/regular-singular
    classifier), not a heuristic gate — future work, tracked with `1823`.
  - No C/behavior change (Kovacic reverted); the deliverable is the converter migration + two corpora +
    two `dsolve_corpus_2_2_2{0,1}_tests` gates (baseline 4) + dashboards. Version 0.137 → 0.138. See the
    §2.2.20 / §2.2.21 blocks in `DSolve_test_status/STATUS.md`.

- **M41 — §2.2.22 corpus (Problems 2101–2200, Nasser Abbasi) + Frobenius integer-root-difference /
  logarithmic second solution.** ✅ DONE. First §2.2.x section dominated by **higher-order (3rd–6th)
  constant-coefficient linear ODEs** (85 of 100: 36 homogeneous `_missing_x`, 49 nonhomogeneous with
  `E^x·poly` / `E^x·(poly Cos + poly Sin)` / mixed forcing), plus a 3rd-order Euler–Cauchy IVP family
  (2106–2111, incl. an x=−1 point and a fully symbolic-IC problem 2111) and 5 2nd-order
  variable-coefficient regular-singular series ODEs (2101–2105). The constant-coefficient bulk — incl.
  resonant nonhomogeneous, `Root`-object spectra, and mixed forcing — solves out of the box; the wave
  drives the section to **100/100, 0 FAIL, 0 crash** with three general fixes and **0 regression**
  (§2.1.2, §2.2.20, §2.2.21 and all `series`/`nseries` ctests held). New gate
  `dsolve_corpus_2_2_22_tests` (baseline 0). Version 0.138 → 0.139.
  - **`DSolve\`FrobeniusSeries` positive-integer indicial-root-difference / logarithmic second
    solution** (`dsolve_frobenius.c`). The header contract's "positive-integer difference → decline
    when a genuine Log is required" gap is closed. For `d = r1 − r2` a positive integer with the
    smaller-root recurrence obstructed, `y2` is built by the **(s − r2)-modified Frobenius derivative
    method**: the symbolic-exponent coefficients `a_n(s)` carry a simple pole at `s = r2` for `n ≥ d`,
    so `b̄_n = lim (s−r2)a_n(s)` (the `Log` coefficient, ∝ y1) and `b̄_n' = lim d/ds[(s−r2)a_n(s)]`
    (the algebraic correction) are removable-singularity limits — a clean generalisation of the
    existing equal-root d/ds path, declining (NULL) only if a limit is non-finite so no wrong answer
    ships. The truncation window is sized to `ceil(r1−r2) + FROB_ORDER` so both series' leading terms
    survive the merged `SeriesData` (a fixed 6-term window had dropped the larger-root `C[1]`). Solves
    obstructed `2104` (roots {0,4}) and wide-window `2101` (roots {5/2,−7/2}).
  - **General `series.c` fix — `scalar · Laurent-SeriesData` order truncation.** Multiplying a
    `SeriesData` with `nmin < 0` by any scalar truncated the product order to `order_num + nmin`
    (`so_from_constant` allocated only `order_num` coefficients; `so_mul`'s
    `min(const.order + series.nmin, …)` then lost `|nmin|` orders). `series_combine` now passes the
    most-negative operand `nmin`, and the constant spans `order_num − min_nmin` coefficients (extends
    only in the Laurent case). This is a general `Series[]` correctness fix and is what recovered
    Frobenius `2103` (roots {3,−7}, difference 10, whose merged series had collapsed to a single
    constant).
  - **Converter `<var>\cos(` juxtaposition** (`tools/latex_ode_to_mathilda.py`). A trig function
    directly after a variable (`15x\cos(2x)`) glued into the bogus symbol `xCos`; `_implicit_mult`
    now inserts the multiplication space before a FUNCS head (`2185`). §2.2.20/§2.2.21 regenerate
    byte-identically.
  - `§2.2.21-2080` (the integer-difference member of that section's residue) stays UNEVAL: its Kovacic
    Case-1 solve overruns the 8 s harness window before Frobenius runs (the documented churn class,
    orthogonal to this Frobenius capability fix). See the §2.2.22 block in `DSolve_test_status/STATUS.md`.

- **M42 — §2.2.23 corpus (Problems 2201–2300, Nasser Abbasi) + undetermined-coefficients hyperbolic
  forcing + converter `\frac`/`\sqrt` juxtaposition.** ✅ DONE. The first **systems-heavy** §2.2.x
  section since §2.2.15: 45 scalar (17 IVP) + **55 constant-coefficient homogeneous 2×2/3×3 linear
  first-order systems** (2238–2292; real / complex / repeated-defective / irrational-`Root` spectra).
  Scalar half is higher-order (3rd–6th) constant-coefficient linear nonhomogeneous (2201–2220),
  3rd/4th-order Euler–Cauchy IVPs (2221–2233), general undefined-`F(x)` forcing (2234–2237), and
  elementary first-order (2293–2300). The constant-coefficient bulk — scalars via
  `LinearConstantCoefficients` / `UndeterminedCoefficients` / `EulerCauchy` / `ExactODE`, systems via
  the `LinearFirstOrderSystem` Jordan matrix exponential — solves out of the box; the wave drives the
  section to **98/100, 0 FAIL, 0 crash** with two general fixes and **0 regression** (§2.1.2, §2.2.4,
  §2.2.9, §2.2.20–22 and all `series`/`nseries` ctests held). New gate `dsolve_corpus_2_2_23_tests`
  (baseline 2). Version 0.139 → 0.140.
  - **`DSolve\`UndeterminedCoefficients` — hyperbolic (`Sinh`/`Cosh`) forcing**
    (`dsolve_undetcoeff.c`). A forcing carrying a hyperbolic×trig/poly/exp product
    (`Sinh[x] Cos[x] − Cosh[x] Sin[x]`, `2203`) is not a UC function as written, so UC declined and
    the equation fell through to the (here **hanging**) variation-of-parameters fallback (timeout).
    The forcing normalisation now expands `Sinh`/`Cosh` to **real** exponentials
    (`Sinh[u]→(E^u−E^{−u})/2`, keeping `Cos`/`Sin`; `TrigToExp` is deliberately avoided so the trig
    atoms do not become complex exponentials the `Coefficient[·,Cos]` matcher loses), so the forcing
    becomes a sum of genuine `E^{a x}{1|Cos|Sin}` atoms — including the resonant `x·E^{±x}(…)` shift.
    A structural no-op when the forcing has no hyperbolic head (non-hyperbolic forcing byte-identical);
    the strict in-method residual gate is unchanged, so no wrong particular ships. Closes `2203` and
    intercepts the `y''−4y == Sinh/Cosh[..]` const-coeff family (§2.2.4-328/329/369, §2.2.9-875/876/
    895) with a cleaner particular.
  - **Converter — `\frac`/`\sqrt` math brace-arg juxtaposition** (`tools/latex_ode_to_mathilda.py`).
    `_implicit_mult` protected a backslash macro *together with its first brace argument*, shielding a
    `\frac` numerator / `\sqrt` radicand from the digit/letter split — so a juxtaposed
    coefficient×function inside a fraction (`\frac{2ty}{t^2+1}`, `\frac{4y_1}{3}`) glued into a bogus
    symbol whose dependent function never got its `[x]`/`[t]` (systems `2239`/`2240`/`2261`, `2300`).
    It now protects only the *name* of the math macros `\frac`/`\sqrt` so their braces stay exposed,
    while name-macros (`\operatorname{Heaviside}`, `\mathrm{e}`) keep their brace arg protected.
    §2.2.20/§2.2.21/§2.2.22 regenerate **byte-identical**.
  - **Residue 2 (both `sympySolved=False`, bounded declines — no wrong answers):** `2289` (3×3 system,
    irreducible-cubic `t³−14t+40` complex-`Root` spectrum — the `Simplify[ComplexExpand[·]]` realifier
    churns on the `Root`-heavy matrix exponential, no solve in 60 s; the M40/M41 irrational-spectrum-
    churn class, where every attempted structural gate was reverted as unsafe) and `2220` (4th-order
    IVP whose homogeneous basis carries the same cubic-`Root` spectrum, so the `Solve` IC-fit cannot
    close over the `Root`-object exponential basis). See the §2.2.23 block in
    `DSolve_test_status/STATUS.md`.

- **M43 — §2.2.24 corpus (Problems 2301–2400, Nasser Abbasi) + IVP inverse-branch correctness +
  homogeneous-Root→implicit + Lagrange linear-coefficients deferral.** ✅ DONE. A **mixed**
  first-order + second-order section: 100 scalar (51 IVP), 0 systems — first-order separable /
  linear / homogeneous (class A/C) / exact / Bernoulli / Riccati / Abel / Lagrange–d'Alembert,
  and second-order `_missing_x` / linear / Emden–Fowler / `_with_linear_symmetries` /
  Gegenbauer–Legendre. Most solves out of the box (the special-function Riccati `y'=t+y²`→Airy
  and `y'=t²+y²`→Bessel, every Abel-tagged case, and Legendre included). Baseline **90/100 with
  1 FAIL** → **92/100, 0 FAIL, 0 crash** with two general root-cause fixes and **0 regression**
  (all §2.2.x + §2.1.2 gates held; §2.2.1/2/6/7 improved 1–2 as a bonus). New gate
  `dsolve_corpus_2_2_24_tests` (baseline 8). Version 0.140 → 0.141.
  - **IVP inverse-branch correctness — a wrong-answer fix** (`dsolve_common.c`). The IVP
    constant-fitter collapsed a multivalued inverse fit's `ConditionalExpression` integer-family
    to its principal member (`C[_]→0`) only at the FINAL step, not while *choosing* among Solve's
    inverse branches. For `t y'==y+√(t²+y²)`, `y[1]==0` (2329), `Solve[Sinh[C]==0,C]` returns an
    odd family `Iπ+2Ikπ` and an even family `2Ikπ`; the candidate loop substitutes a generic
    non-integer for the free family index, so BOTH look complex-nonzero and it fell back to
    `args[0]` (odd/wrong), shipping `y=(1−t²)/2` — meets the IC, fails the ODE (correct
    `(t²−1)/2`), and the transcendental residual is undecidable by `zero_test`. Fix: a new
    `ds_collapse_principal` collapses each candidate to its principal member BEFORE the numeric
    `ds_branch_num_ok` check (correct even family wins), plus a scoped final numeric-verify gate
    (first-order scalar, radical bodies excluded) that drops any confidently-wrong fit. Bonus:
    closes first-order IVPs in §2.2.1/2/6/7.
  - **`DSolve\`Lagrange` defers the linear-coefficients class** (`dsolve_lagrange.c` + new public
    predicate `dsolve_is_linear_coefficients_form` in `dsolve_lincoeff.c`). Every affine-ratio
    `y'==(a1 x+b1 y+c1)/(a2 x+b2 y+c2)` also matches Lagrange as `y==x F(y')+G(y')` with `F`
    rational in `y'`, but Lagrange's integrating-factor linear ODE spins uninterruptibly
    (`TimeConstrained` cannot bound an inner `Integrate` — see
    [[project_timeconstrained_no_nest]]); it timed out on 2335, which is NOT a genuine d'Alembert
    equation. Lagrange now defers the class to `LinearCoefficients`, which solves it in ~1 s.
    Genuine d'Alembert (polynomial `F`, e.g. `y==2x y'+y'²`) is unaffected — it still solves via
    Lagrange's parametric path.
  - **Deferred-class attempt** (per request): the special-function Riccati (Airy/Bessel) and every
    Abel-tagged case in this section ALREADY solve — no gap, no new method needed. **Residue 8**
    (all `sympySolved=False`, bounded declines — no wrong answers): `2304` (non-elementary
    integrating-factor integral `∫(t²+1)^{1/2}(t⁴+1)^{1/4}dt`), `2327` (symbolic-coefficient
    Erf-Riccati — `u''+k(t+b)u'+k²tb u==0` factors `(D+kt)(D+kb)` with an Erf second solution;
    solves for concrete k,b, symbolic needs a symbolic 2nd-order operator-factoring method),
    `2349`/`2350`/`2351` (`y'=e^{−t²}+y²` → `u''=e^{−t²}u`, non-elementary), `2352`/`2355`/`2356`
    (`y=_G(x,y')` non-elementary). See the §2.2.24 block in `DSolve_test_status/STATUS.md`.

- **M44 — §2.2.25 corpus (Problems 2401–2500, Nasser Abbasi) + VariationOfParameters
  fractional-power Simplify-hang fix + three converter root-cause fixes.** ✅ DONE. A
  **second-order-linear-heavy** section: 100 scalar (32 IVP), 0 systems — constant-coefficient
  nonhomogeneous (incl. `sec(t)` and resonant fractional-power `t^{5/2}e^{-2t}` forcing), the
  orthogonal-polynomial / special-function families (Airy `y''=ty`, Hermite, Legendre, Chebyshev,
  Laguerre, Gegenbauer, Bessel at symbolic order `v`), Euler–Cauchy, regular-singular (Frobenius),
  Emden–Fowler nonlinear, and first-order separable / linear / exact / homogeneous. The 2nd-order
  stack (Kovacic / SpecialFunctionForm hypergeometric+Bessel+Airy recognizers / Frobenius fallback
  / UndeterminedCoefficients / VariationOfParameters / EulerCauchy) solves the bulk out of the box.
  Baseline **95/100, 0 FAIL, 0 crash** → **96/100** with one general root-cause fix and **0
  regression**. New gate `dsolve_corpus_2_2_25_tests` (baseline 4). Version 0.141 → 0.142.
  - **`dsolve_variation_of_parameters` fractional-power Simplify hang** (`dsolve_common.c`). A
    resonant fractional-power forcing (`y''+4y'+4y == t^{5/2}e^{-2t}`, 2406) closes to the
    elementary `4/63 t^{9/2}e^{-2t}`, but the routine ended with an unconditional
    `ds_simplify(yp)`, and `Simplify[t^{5/2}e^{-2t}]` itself HANGS (the documented
    radical/pseudo-remainder Simplify pathology; the harness's own `TimeConstrained` forbids a
    nested one — see [[project_timeconstrained_no_nest]]). The VoP body is already auto-evaluated to
    a clean form, so the final Simplify is now SKIPPED when the answer carries a fractional power
    (new static `ds_has_fractional_power`). Correct answer, no hang, no nesting; intercepts the
    whole resonant-fractional-forcing class. Non-fractional answers are byte-identical.
  - **Three converter root-cause fixes** (`tools/latex_ode_to_mathilda.py`), each a latent bug the
    prior LaTeXML sections had masked; **§2.2.20–23 regenerate byte-identically**, and §2.2.24 was
    regenerated (7 records corrected — see its STATUS block; 2327 now solves, honest **91/100**,
    gate 8→9): (a) **independent-variable juxtaposition** — a preferred indvar letter glued to a
    digit or the dependent letter (`ty`, `2t`, `3t²`) was invisible to the word-boundary scan, so it
    defaulted to `x`, silently generating an `x`/`t`-mixed WRONG equation (8 recs here, 6 in
    §2.2.24); fixed by a fallback re-scan on the implicit-multiplication-separated body. (b)
    **autonomous parameter as indvar** — `y'=k(a−y)(b−y)` picked the parameter `a` (a spurious
    `dy/da` Riccati that ABORTs); a unique-candidate rule now falls through to a fresh `x` (2498;
    §2.2.24-2327). (c) **piecewise/cases forcing** — the `PWFORCE` sentinel was shredded by the
    implicit-multiplication pass and the `\begin{array}[]{cc}` optional arg leaked; a
    whitespace-tolerant expand + optional-arg strip now yields a clean `Piecewise[…]` (2487).
  - **Residue 4** (all `sympySolved=False` — SymPy fails them too, bounded declines, no wrong
    answers): `2409` (`y''+t²y/4 == f cos t`, parabolic-cylinder homogeneous part + forcing —
    needs a ParabolicCylinderD recognizer), `2410` (`(t²+1)y''−2ty'+2y = t²+1`; the homogeneous
    solves via Kovacic but with a complex-radical basis `√(t−i)√(t+i)` that does not reduce to
    `√(t²+1)`, so VoP for the forcing cannot close — a clean-basis / polynomial-solution VoP would
    give the elementary answer), `2444` (transcendental-coefficient power series, slower than the
    8 s harness bound), `2477` (non-elementary integrating factor, same class as §2.2.24-2304). All
    §2.2.x + §2.1.2 gates held; see the §2.2.25 block in `DSolve_test_status/STATUS.md`.
- **M45 — §2.2.26 corpus (Problems 2501–2600, Nasser Abbasi) + Bernoulli integrating-factor
  symbol-leak fix + new polynomial-solution 2nd-order method.** ✅ DONE. A Braun-textbook
  MIXED section: 100 scalar (50 IVP), 0 systems — a first-order-nonlinear-heavy front half
  (homogeneous class A/C, dAlembert, Bernoulli, Abel, exact, separable, Riccati) and a
  second-order block (constant-coefficient `_missing_x`, Euler–Cauchy, Emden–Fowler,
  Gegenbauer, `_with_linear_symmetries`). Baseline **87/100, 0 FAIL, 0 crash** → **89/100**
  with two general root-cause fixes and **0 regression**. New gate `dsolve_corpus_2_2_26_tests`
  (baseline 11). Version 0.142 → 0.143.
  - **Bernoulli integrating-factor `DSolve\`Y` leak** (`dsolve_bernoulli.c`). The linearised
    Bernoulli coefficients A, B are mathematically free of the reduction variable Y but stored
    TEXTUALLY carrying a frozen `Q = FY − Y F_Y` (`Times` does not distribute over `Plus`).
    When the reduced integrating-factor integral is ELEMENTARY the evaluator collapses it, but
    when it is NON-elementary (`y' == (1+cos 4t)/4·y − (1−cos 4t)/800·y²`, 2532) the frozen
    `DSolve\`Y` leaked into the unevaluated `Integrate`, so back-substitution scored UNEVAL.
    `Cancel` each coefficient to lowest terms in Y before the linear solve — the same cheap
    rational-GCD already used for the exponent n (never `Simplify`, which hangs on radical
    coefficients). Intercepts the whole non-elementary-integrating-factor Bernoulli class;
    Y-free coefficients are unchanged.
  - **Polynomial-solution 2nd-order method** (new `src/calculus/dsolve_ratsol2.c`,
    `DSolve\`PolynomialSolution`). A homogeneous 2nd-order linear ODE with rational coefficients
    whose fundamental set is polynomial (`(t²+1)y''−2ty'+2y==0` → `{t, t²−1}`) is solved by a
    degree-bounded undetermined-coefficient search (ansatz `Σ a_k x^k`, clear the denominator,
    `Solve` the coefficient equations for the a_k null space); when forced, variation of
    parameters over that clean basis. Runs BEFORE Kovacic, which returns the same set wrapped in
    a complex-radical basis `√(t−i)√(t+i)` that does not reduce and blocks the VoP integrals.
    Every returned solution is numerically self-verified (`ds_branch_num_ok`). Fixes **2592 here
    AND §2.2.25-2410** (that section 96→**97/100**, gate 4→3). Const-coefficient / Euler /
    special-function families are still claimed by their specialists earlier in the cascade
    (verified: `y''−y==0`, Euler, `y''==0` unchanged).
  - **Residue 11** (all `sympySolved=False` — SymPy/Mathematica fail them too, bounded declines,
    no wrong answers): non-integrable Riccati `y'=e^{−t²}+y²` (`2524`/`2525`/`2526`), Abel
    `y'=y³+e^{−5t}` (`2528`), implicit `y=G(x,y')` / `x=G(y,y')` forms not solvable for `y'`
    (`2514`/`2527`/`2530`/`2531`/`2537`), rational `y'=(t²+y²)/(1+t+y²)` (`2539`), and the
    abstract-coefficient `y''+p(t)y'+q(t)y==1+t` (`2591`, unsolvable for arbitrary p, q). All
    §2.2.x + §2.1.2 gates held; see the §2.2.26 block in `DSolve_test_status/STATUS.md`.

- **M46 — §2.2.27 corpus (Problems 2601–2700, Nasser Abbasi) + corpus-harness
  chained-inequality crash fix.** ✅ DONE. A 2nd-order-LINEAR-heavy section: 96 scalar
  (39 IVP) + 4 systems — constant-coefficient nonhomogeneous
  (UndeterminedCoefficients/VariationOfParameters), `_missing_y` reduction, Emden–Fowler
  (Airy/Bessel-reducible `y''±tⁿy`, 13/13), 2nd-order exact-linear, orthogonal-polynomial
  (Gegenbauer/Jacobi/Laguerre) + Bessel + Lienard special forms, `_with_linear_symmetries`,
  and 4 constant-coefficient 2×2 linear systems (two forced, 4/4). Baseline **97/100, 0 FAIL,
  1 CRASH** → **98/100, 0 FAIL, 0 crash** with one general root-cause fix and **0 regression**.
  New gate `dsolve_corpus_2_2_27_tests` (baseline 2). Version 0.143 → 0.144.
  - **Corpus-harness chained-inequality crash** (`DSolve_test_status/dsolve_corpus_prelude.m`,
    `dsFreeParams`). A constant-coefficient IVP with Piecewise/step forcing
    (`y''+y'+7y == Piecewise[{{t, 0<=t<2},{0, 2<=t}}]`, 2690) back-substitutes to a residual
    carrying the **chained** inequality `0<=t<2`, represented as
    `Inequality[0, LessEqual, t, Less, 2]` with the comparison operators in **argument**
    position. The verifier's free-parameter collector (`Cases[resid, s_Symbol, {0,∞}]`, whose
    "instantiate params, never Heads" reasoning assumes operators occupy head position) then
    treated `LessEqual`/`Less` as free parameters and substituted numbers for them, corrupting
    the `Piecewise` into garbage that **SIGSEGV**ed when the residual was differentiated
    (`Integrate` derivative-divides → `D`/`higher_order_partial` → NULL-arg deref in
    `evaluate_step`). Fix: exclude relational/logical/piecewise operator symbols (`Less`,
    `LessEqual`, `Greater`, `GreaterEqual`, `Equal`, `Unequal`, `Inequality`, `And`, `Or`,
    `Not`, `Xor`, `Piecewise`, `UnitStep`, `HeavisideTheta`, `DiracDelta`) from the collector —
    they are never ODE parameters (the converter emits value symbols only), so the exclusion is
    **monotone-safe**. Single-sided conditions (`Less[t,2]`, operator as head) never triggered
    it, which is why the piecewise-heavy §2.2.15 stayed green — only a chained inequality
    reaches it. The CAS itself is robust on legitimate chained-inequality `Piecewise`
    (`D`/`Integrate` verified clean); the fault was entirely in the test harness's substitution.
    With the fix, 2690's variation-of-parameters answer verifies (the residual's inert integrals
    cancel by the fundamental theorem) → PASS. **Regression:** the 3 sections whose corpus carries
    a chained-inequality Piecewise/step condition were re-run — §2.2.15 (2/2), §2.2.16 (3/3),
    §2.2.25 (3/3) all hold at baseline; the remaining sections contain no such operators, so
    `dsFreeParams` is provably unchanged for them. No C code changed.
  - **Residue 2** (both `sympySolved=False` — SymPy fails them, bounded UNEVAL declines, no
    wrong answers): the 2nd-order **exact** nonhomogeneous `y''+t³y'+3t²y==eᵗ` (`2621` —
    integrates once to the linear `y'+t³y==eᵗ+C`, whose integrating-factor integral
    `∫e^{t⁴/4}eᵗ dt` is non-elementary), and the **transcendental-coefficient**
    `(1−t²)y''+y'/Sin[1+t]+y==0` (`2641`, no closed form). All §2.2.x + §2.1.2 gates held; see the
    §2.2.27 block in `DSolve_test_status/STATUS.md`.

- **M47 — §2.2.28 corpus (Problems 2701–2800, Nasser Abbasi) + radical-root VoP integrand
  simplify + homogeneous zero-IVP shortcut.** ✅ DONE. A SYSTEMS-heavy section: 18 scalar
  (3 IVP) + 82 systems — constant-coefficient linear systems (2×2/3×3/4×4, homogeneous +
  `E^t`/trig/impulse/step forcing), 8 nonlinear systems (Lotka–Volterra / epidemic /
  competition), and 18 higher-order (3rd–6th) constant-coefficient scalar linear ODEs
  (homogeneous `_missing_x` + arbitrary-forcing `_missing_y`/VoP). Baseline **84/100, 0 FAIL,
  0 crash** → **86/100, 0 FAIL, 0 crash** (scalars **18/18**) with two general root-cause fixes
  and **0 regression**. New gate `dsolve_corpus_2_2_28_tests` (baseline 14). Version 0.144 → 0.145.
  - **Converter — subscripted arbitrary forcing functions** (`tools/latex_ode_to_mathilda.py`).
    The source's `f_1(t)`, `f_2(t)` (and the author's inconsistent text form `\textit{f\_1}(t)`)
    normalise to the symbols `f1`,`f2` but were then read as multiplication (`f1 t`) — a fidelity
    bug like §2.2.24's `t`-as-constant. Generalised the arbitrary-function detection/application
    from single letters `f`/`g`/`h` to subscripted `[fgh][0-9]+` (and fold `\_`→`_`), so
    `f1[t]`,`f2[t]` are emitted as genuine function applications (2708/2762/2780). Byte-identical
    on prior LaTeXML sections (§2.2.27 regenerates identically).
  - **Constant-base-radical VoP integrand simplify** (`dsolve_common.c`,
    `dsolve_variation_of_parameters`). An irrational-root fundamental set — `y''''+y==g[t]`
    (`2719`), roots `(±1±i)/√2`, basis `E^(±t/√2) Cos/Sin[t/√2]` — makes each per-term Cramer
    integrand hundreds of leaves, so the integrator churns ~2 s/term and four terms overrun the
    8 s solve budget (→ UNEVAL), while the sibling `y''''−y==g[t]` (`2718`, clean basis
    `{E^t,E^-t,Cos,Sin}`) closes fast. Fix: `Simplify` the per-term integrand before integrating
    **when it carries a constant-base radical** (`2^(-1/2)` in the exponent) and **no** `t^(p/q)`
    of the variable (collapses ~485 → ~44 leaves; the inert convolution integral then closes
    instantly). New `ds_has_var_fractional_power` distinguishes the safe constant-base radical
    from the documented `t^(5/2) E^(-2t)` Simplify hang (§2.2.25-`2406`), which is left untouched;
    a clean-root VoP has no fractional power and is byte-identical.
  - **Homogeneous-linear-IVP zero-ladder → `y≡0`** (`dsolve_common.c`, `dsolve_fit_constants`).
    A 4th-order homogeneous IVP with an *irreducible* characteristic polynomial (`2713`:
    `r⁴+4r³+14r²−20r+25`, irreducible over ℚ) has a `Root[]`-object fundamental set; `Solve`
    bubbles unevaluated on the Root-coefficient constant-fit system, so the IVP kept its
    constants → UNEVAL. But the conditions are the **complete** zero derivative-ladder
    `y^(k)(x0)==0`, k=0..n−1, at a single point — a well-posed IVP whose unique solution is
    `y≡0` (nonsingular Wronskian forces `C=0`). Fix: substitute every generated constant with 0
    under a tight guard (one condition per constant, all values the literal 0, one base point,
    order set exactly `{0..n−1}`), which excludes a BVP (`y[0]==0, y[Pi]==0 → C[2] Sin[x]`) and
    any under/over-determined set, so it can never turn a genuine non-trivial answer to 0.
  - **Residue 14** (all systems, no scalar gap; `sympySolved=False` where noted): 8 nonlinear
    systems `2788`–`2795` (Lotka–Volterra / epidemic / competition, no closed form), the 3×3
    `E^t` system `2785` (`sympy=False`), system arbitrary forcing `f1[t]`/`f2[t]`
    (`2708`/`2762`/`2780`, needs system-level VoP over arbitrary functions), and system
    DiracDelta/UnitStep forcing (`2781`/`2782`, needs the systems Laplace/Green's path) — deep
    systems-solver work, deferred. All §2.2.x + §2.1.2 gates held; see the §2.2.28 block in
    `DSolve_test_status/STATUS.md`.

- **M48 — §2.2.29 corpus (Problems 2801–2900, Nasser Abbasi) + three general converter
  transcription fixes.** ✅ DONE. A first-order-nonlinear + linear-systems section: 75 scalar
  (21 IVP) + 25 systems — 2×2/3×3/4×4 constant-coefficient linear systems (incl. subscripted
  `x1`/`x2` and 4-variable `x,y,z,h`), autonomous 2nd-order `_missing_x` reducibles
  (`z''+g(z)==0`), Sturm–Liouville eigenvalue BVPs (`y''+λy==0` with symbolic boundary `L`), and
  a large first-order block (separable / linear / homogeneous class A/C / dAlembert / Abel /
  Bernoulli / exact). Baseline **88/100, 0 FAIL, 0 crash** (scalars **70/75**, 93.3%) — **no
  solver change**: the first-order + linear-system stack solves the section out of the box. The
  wave was three general **converter** fixes (`tools/latex_ode_to_mathilda.py`), each a
  transcription-fidelity bug found by spot-audit and each verified **byte-identical** on the
  prior LaTeXML sections (§2.2.27/§2.2.28 full old-vs-new diff) and on §2.1.2's arbitrary-
  function records (direct `detect_symbols` comparison). New gate `dsolve_corpus_2_2_29_tests`
  (baseline 12). Version 0.145 → 0.146.
  - **`\sqrt` letter-juxtaposition glue** (`_implicit_mult`). A variable letter juxtaposed with
    `\sqrt` (`x-k\sqrt{x²+y²}`, `2890`; `t\sqrt{1-y²}`, latent prior `2.2.24-2360` / `2.2.26-2536`)
    glued to the `Sqrt` head after `replace_sqrt` as the bogus single symbol `kSqrt[...]` /
    `tSqrt[...]`, because `\sqrt` is protected (to expose its radicand) before the letter-split
    and is not a FUNCS head. Fix: split a **letter** before `\sqrt`; a **digit** (`2\sqrt{x}` →
    `2Sqrt[x]`) already parses as multiplication and is left byte-identical. (`\frac` needs no
    such rule — it lowers to `((n)/(d))`, not a named head.) This also corrected the two latent
    prior records; both were already PASS (formal implicit separable solutions even on the
    garbled head) and remain PASS with the faithful equation, so those gates are unchanged (9/11).
  - **`\textit{x\_}N` italic-glued subscript** (`normalize_subscripts`). One system (`2824`)
    rendered its subscripted variables as `\textit{x\_}1` (italic core, index *outside* the
    brace) instead of the clean `x_{1}` used by siblings `2811`/`2825`; the trapped underscore
    defeated subscript+prime handling, so `^{\prime}` degraded to a literal `^(prime)` power and
    the dependent functions were never detected (function list defaulted to `{y}`). Fix: fold
    `\textit{X\_}N → X_{N}` early. Distinct from §2.2.28's `\textit{f\_1}` (index *inside* the
    brace), which is untouched.
  - **System dependent variable `h` misread as an arbitrary function** (`detect_symbols`).
    `2806`/`2807` are 4-variable linear systems `x,y,z,h`, but `h` (a conventional arbitrary-
    function letter, `ARBFUN={f,g,h}`) was dropped from the function list. Fix: a symbol that
    **heads its own derivative row** (start-anchored `h^{\prime}&=−2z`) is a dependent variable
    even when its letter is in `ARBFUN`; an arbitrary function differentiated only *inside*
    another equation's body (`f'(x)`/`g'(x)` in a scalar Abel/Riccati ODE, §2.1.2) is not a row
    head and stays arbitrary — the promotion is gated on a start-anchored match, never a
    substring hit (verified §2.1.2 arbitrary-function records unchanged).
  - **Residue 12** (all no elementary closed form): 7 nonlinear systems `2811`/`2813`–`2818`
    (Riccati-type / Lotka–Volterra / coupled-nonlinear, `sympy=False`); 4 autonomous 2nd-order
    `2820`/`2821`/`2822`/`2823` (`z''+g(z)==0`, Duffing/hyperelliptic, `sympy=False`); and `2819`
    (`z''+z³==0`, `sympy=True` — energy first integral gives an *elliptic-integral* implicit
    solution; DSolve times out at the 8 s bound and the form is not back-substitution-verifiable
    — deep elliptic-ODE work, deferred). All §2.2.x + §2.1.2 gates held; see the §2.2.29 block in
    `DSolve_test_status/STATUS.md`.

- **M49 — §2.2.30 corpus (Problems 2901–3000, Nasser Abbasi) + Greek-variable converter fix +
  inverse-hyperbolic integrating-factor solver fix.** ✅ DONE. An all-first-order section: 100
  scalar (24 IVP), 0 systems — a large homogeneous class A/C/G + Abel(2nd type) + rational +
  Bernoulli + linear + exact + separable block, plus symmetry / dAlembert cases. Baseline
  **93/100, 0 FAIL, 0 crash**. Two general root-cause fixes, each verified 0-regression. New gate
  `dsolve_corpus_2_2_30_tests` (baseline 7). Version 0.146 → 0.147.
  - **Greek-letter macro as a first-class variable** (`latex_ode_to_mathilda.py`). The
    single-letter (`[A-Za-z]`) symbol-detection layer could not see `θ` (`\theta`) as a variable.
    2972/2973/2992 are `dr/dθ` ODEs — θ the **independent** variable — so the converter defaulted
    the indvar to `x` and read θ as a constant; 2984 is `sin(θ)·θ'(t)+…=0` — θ the **dependent**
    function of `t` — so the derivative `\theta^{\prime}` mangled to a literal `^(prime)` power and
    the function list defaulted to a fallback letter. Fix, four coordinated parts: (a) primed-symbol
    detection recognises a Greek macro as a dependent variable (`NAME_RE`); (b) a lone present Greek
    letter is adopted as the independent variable **only for a first-order ODE** — the essential
    gate, so an autonomous 2nd-order eigenvalue problem `y''+λy=0` (§2.2.29-2834ff, `_missing_x`)
    keeps `λ` a parameter and a fresh `x` as indvar; (c) `convert_side` canonicalises Greek
    *variable* macros to ASCII before the mains/derivative pass; (d) placeholder-expansion accepts
    multi-character identifiers. **Byte-identical** verification: OLD-vs-NEW converter on the same
    fetched HTML for all of §2.2.20–§2.2.29 (IDENTICAL) and §2.1.2's arbitrary-function
    `detect_symbols` (direct comparison) — the only Greek variable names anywhere in the corpus are
    the four new records.
  - **Inverse-hyperbolic integrating factor** (`dsolve_linear_factor_solve`,
    `src/calculus/dsolve_common.c`). 2980 is the *linear* ODE `(x²−1)y'+4y=−(x²−1)²`, whose
    integrating factor is `μ=Exp[∫4/(x²−1)dx]=Exp[−4 ArcTanh[x]]`. The existing IF cleanup
    (`PowerExpand[Simplify[μ]]`, tuned for the trig `Tan→Sec` case) leaves `Exp[c ArcTanh]` intact,
    so `Integrate[μ q]` spun past the 8 s harness bound → UNEVAL. Fix: when the IF still carries an
    `ArcTanh`/`ArcCoth` factor **and** no `ArcTan`/`ArcSin`/`ArcCos` (those rationalise to *complex*
    powers, not simpler), `TrigToExp` rewrites the inverse hyperbolics to Logs and `Simplify`
    collapses `Exp[…]` to the real algebraic `(x−1)^a(x+1)^b`, after which the `μ q` integral is
    elementary. General solution + `y[0]=−6` IVP both verify by back-substitution (residual 0);
    every trig / polynomial-exponent IF the pipeline already handled is byte-identical.
  - **Residue 7**: 4 `sympy=False` with no elementary closed form (`2923`/`2944`/`2948`/`2955`);
    3 deferred `sympy=True` gaps — `2933` (dAlembert with a non-elementary intermediate integral),
    `2970` (`_with_symmetry_[F(x)*G(y),0]` Lie-symmetry case, DSolve spins to the bound), `2979`
    (correct 4-branch nested-radical general solution, but IVP fit over the radical branches fails).
    All §2.2.x + §2.1.2 gates held; see the §2.2.30 block in `DSolve_test_status/STATUS.md`.

- **M50 — §2.2.31 corpus (Problems 3001–3100, Nasser Abbasi) + corpus-verifier exact-radical
  hang fix.** ✅ DONE. 100 scalar (20 IVP), 0 systems, in two blocks: a first-order block
  (3001–3056 — homogeneous class A/C/G + Abel(2nd type) + rational + Bernoulli + linear + exact +
  separable) and a constant-coefficient **linear** block (3057–3100 — 2nd/3rd/high-order
  `_missing_x` + quadrature). Baseline **99/100, 0 FAIL, 0 crash** — **no solver change** (the
  first-order stack + `LinearConstantCoefficients` solve the section out of the box). New gate
  `dsolve_corpus_2_2_31_tests` (baseline 1). Version 0.147 → 0.148. One general **harness**
  root-cause fix, verified 0-regression:
  - **Real-valued sample point in the numeric back-substitution verifier**
    (`DSolve_test_status/dsolve_corpus_prelude.m`, `dsResidVerdict`). 3043
    (`y'x = 2y + 2x⁴y³`, `y[1]=1`) is solved instantly and correctly by DSolve
    (`y = 1/√((3/2 − x⁸/2)/x⁴)`; residual `Simplify`s to 0, IC fits), yet scored a spurious
    **UNEVAL**: the verifier substituted an **exact rational** sample point into the radical
    residual, forming `(hugeRational)^(3/2)`, whose square-factor extraction factors a ~20-digit
    integer and hung `Power` for minutes (at x = 193/130, 243/130 the exact substitution never
    returned; the per-case fork alarm then scored it non-PASS). The fix substitutes a
    high-**precision real** for the independent variable — `resid /. iv -> N[pt, 24]` on the fast
    tier, `N[pt, 210]` on the high-precision escalation tier — so every power is a fast floating
    power. The verifier's job is *numeric* back-substitution, so it must numericize the point, not
    do exact radical algebra; the residual stays symbolic, so the two-tier catastrophic-cancellation
    robustness (eigenvalue-64 systems, `2.2.11-1001`) is preserved. General across every section —
    it also incidentally fixed the same-class hang in prior sections (e.g. `2.2.1-35`).
  - **Residue 1**: `3049` (`3xy + (3x²+y²)y' = 0`, `y[0]=1`) — DSolve returns the correct 4-branch
    homogeneous general solution, but every explicit radical branch factors `x` (so `y(0)=0`
    structurally) and the IVP fitter cannot fit `y[0]=1`; the clean IVP answer is the *implicit*
    first integral `y⁴+6x²y² = 1`. Same deferred class as `2.2.30-2979`. Honest UNEVAL, not a wrong
    answer. All §2.2.x + §2.1.2 gates re-run (standalone) — all hold. See the §2.2.31 block in
    `DSolve_test_status/STATUS.md`.

- **M51 — §2.2.32 corpus (Problems 3101–3200, Nasser Abbasi); symbolic complex-root basis fix
  investigated + reverted.** ✅ DONE. 100 scalar (8 IVP), 0 systems — entirely constant-coefficient
  **linear** ODEs (2nd/3rd/high-order nonhomogeneous) with polynomial / exponential / sinusoid
  forcing + resonance, plus variation-of-parameters-only forcing (`Sec`/`Tan`/`Csc`/`Log`). Baseline
  **96/100, 0 FAIL, 0 crash** — **no solver change** (solved out of the box by
  `UndeterminedCoefficients` + `LinearConstantCoefficients`). New gate `dsolve_corpus_2_2_32_tests`
  (baseline 4). Version 0.148 → 0.149.
  - **Symbolic complex-root wrong answer (root-caused, fix DEFERRED).** `3155` (`y''+a²y==Sec[a x]`)
    is a **masked wrong answer** for a symbolic coefficient: `dsolve_homog_basis`
    (`src/calculus/dsolve_common.c`) realifies the complex roots `±√(-a²)` to `Exp[Re x](Cos,Sin)[Im x]`
    but only concretizes `Re`/`Im` when the root is **numeric**; for a symbolic coefficient `Re`/`Im`
    do not resolve to explicit reals (`ComplexExpand[Re[-½√(-4a²)]]` leaves `Arg[-4a²]`/`(a⁴)^(1/4)`),
    so the basis does not back-substitute and the VoP particular is wrong — DSolve returns
    `Sec[a x]/a² + Re/Im` mush (a genuine REPL wrong answer for the whole `y''+(symbol)²y` class, only
    *masked* as UNEVAL in the corpus by the verifier's leaked-`C[k]` rule). The complex-**exponential**
    basis `Exp[r x]` for the symbolic case is mathematically correct (`(√(-a²))² = -a²` exactly, so it
    back-substitutes) and was implemented, but the `√(-a²)`-laden exponential form **slows the
    2nd-order cascade past the 8 s cold budget on 13 symbolic-coefficient §2.1.2 cases** (58/439/877/…
    Pöschl-Teller, shifted-Euler, exponential-potential — 557→544 §2.1.2 PASS), a net regression, so it
    was **reverted** (M51 note in `dsolve_common.c`). A narrower fix that does not add cascade latency
    is future work.
  - **Residue 4:** `3155` (masked wrong, deferred as above); `3165` (`y''+y==Tan[x/3]²`, correct VoP
    but cold `DSolve` > 8 s; solves warm); `3161`/`3164` (`y''+4y==Sec[x]Tan[x]`, `y''+9y==Csc[2x]`,
    correct VoP whose cold solve sits at the 8 s boundary — timing-sensitive). Because the fix was
    reverted, **no solver behavior change lands** — every prior gate is behaviorally unchanged from its
    confirmed baseline (§2.2.1–31 additionally re-run green; §2.1.2's code is byte-identical to the
    M50-confirmed 647 ≤ 655, a clean re-run being impossible under heavy external machine load that by
    itself swings §2.1.2 647→667). See the §2.2.32 block in `DSolve_test_status/STATUS.md`.

- **M52 — symbolic complex-root wrong answer fixed (narrow, latency-safe).** ✅ DONE. The
  M51-deferred wrong answer for the `y'' + (symbol)² y` variation-of-parameters class
  (`DSolve[y''+a²y==Sec[a x]] → Sec[a x]/a² + Re/Im` mush) is resolved. Root-caused to **two
  independent** defects, each fixed narrowly so the M51 latency regression cannot recur:
  - **Homogeneous basis** (`dsolve_homog_basis`, `src/calculus/dsolve_common.c`). A symbolic
    complex conjugate pair realified to `Exp[Re x](Cos,Sin)[Im x]` with `Re`/`Im` left
    **unevaluated** (they do not concretize on a radical root), so the basis did not
    back-substitute. For the **pure-imaginary** pair (`r_c == -r`) the real frequency is now
    `β = Simplify[PowerExpand[√(-r²)]]` (field arithmetic; `β² = -r²` exactly, and
    `PowerExpand` collapses `√(a²)→a` so `β` matches the forcing and VoP integrates), giving
    the real `Cos[a x]`/`Sin[a x]` basis. Restricted to the pure-imaginary pair, so general
    complex-root cases (α≠0) stay byte-identical — the 13 symbolic-coefficient §2.1.2 cases
    (58/439/877/… Pöschl-Teller, shifted-Euler, exponential-potential) that the
    M51-reverted whole-branch `Exp[r x]` pushed past the 8 s cold budget are untouched
    (re-timed cold ≈5–6 s, all still solving). Numeric path (`NumericQ`→`ComplexExpand`)
    unchanged, so the M24 complex-cube-root IVP still fits.
  - **Undetermined coefficients** (`dsolve_undetcoeff.c`). Independently, UC **accepted** the
    non-UC forcing `Sec[a x]` (`y_p = Sec[a x]/a²`), winning the cascade ahead of
    `LinearConstantCoefficients`' correct VoP. UC's residual zero-test gate (its intended
    safety net) was defeated by a `zero_test` **false positive** — `PossibleZeroQ` of the raw
    residual `D[Sec[a x],{x,2}]/a²` is `True` for a symbolic parameter though genuinely
    nonzero (the numeric case already declined). The gate now `Simplify`s the residual first,
    so the non-UC term declines to VoP; genuine UC residuals are provably zero either way, so
    real accepts (poly/exp/sinusoid, incl. resonance) are unchanged.

  Verified: new unit `t_m52_symbolic_complex_root` (`tests/test_dsolve.c`); all DSolve ctest
  suites (unit + M12/M14/M18 stress) + `make check-c99` green; the three named §2.1.2 latency
  cases re-timed within budget. Corpus §2.2.32-3155 flips from masked-UNEVAL to PASS. Version
  0.149 → 0.150.

- **M53 — §2.2.33 corpus (Problems 3201–3300) + IVP condition-verification.** ✅ DONE. 12000.org
  §2.2.33: 100 records (93 scalar [14 IVP] + 7 systems), converted via
  `tools/latex_ode_to_mathilda.py`. Baseline **77/100 PASS, 0 FAIL, 23 UNEVAL, 0 crash** (scalar
  72/93, systems 5/7). New gate `dsolve_corpus_2_2_33_tests` (baseline 23); report
  `DSolve_test_status/reports/2.2.33.md`. Dominant gap the `2nd_reducible_mu` class (13 UNEVAL —
  the **M18 Stage-2** target), plus quadrature / dAlembert / rational and one converter miss
  (3296, a multi-line `array` ODE row the converter cannot extract).
  - **IVP condition-verification (0-FAIL restored).** 3279 (`y''==(y')² Sin[x]`, `y[0]==0`,
    `y'[0]==1/2`) initially FAILed — DSolve shipped `y==0`, satisfying the ODE and `y[0]==0` but
    **violating `y'[0]==1/2`** (correct is `Tan[x/2]`). Pre-existing bug: `Solve` returns the
    degenerate `C[1]→1` at a *removable singularity* of the fit system and the constant-fitter
    accepted it; the spurious-fit backstop (`dsolve_fit_constants`, `src/calculus/dsolve_common.c`)
    checked only the ODE residual (which `y==0` passes) and only at first order. New helper
    `ds_fit_meets_conditions` numerically verifies a fitted body against every point condition at
    **any** order; the backstop now rejects a condition-violating fit (→ honest decline) rather
    than shipping a wrong answer. Conservative (rejects only a robustly-nonzero condition residual),
    so genuine fits are untouched. Regression-checked: control IVPs (incl. nonhomogeneous zero-IC
    `y''+y==x`, BVP, M24) unchanged; IVP-heavy §2.2.12 (3≤3) and §2.2.14 (1≤1) within baseline;
    `make check-c99` green. Version 0.150 → 0.151.

- **M54 — Kovacic `Q==0` early-decline (reducible-ODE hang fix).** ✅ DONE. Targeting the
  §2.2.33 reducible-μ gap, investigation showed the dominant concrete failure is not the (M18
  note's own 0-yield) integrating-factor method but a **Kovacic churn**: `(1-x²)y''+xy'==1`
  (§2.2.33-3256) **hung** because it has no y-term (`Q≡0`) — first order in `y'`, solved by
  `ReductionOfOrder` in ~0.05 s — yet `DSolve\`Kovacic` (earlier in the cascade) has a nonzero
  normal form and its Case-1 search churns to its 5 s budget before declining, starving the
  cascade past budget. Fix (`src/calculus/dsolve_kovacic.c`): decline at once when `Q≡0`
  (Kovacic is never uniquely needed for a y-free equation). 3256 now solves and back-substitutes;
  genuine Kovacic cases (nonzero y-term) are untouched (the gate keys on `Q`, not the search).
  New unit `t_m54_kovacic_missing_y_no_churn`; §2.2.33 78/100 PASS (3256 UNEVAL → PASS, 0 FAIL,
  no other change), gate baseline 23 → 22; §2.1.2 (Kovacic's home turf) re-run within baseline;
  `make check-c99` green. Version 0.151 → 0.152.

- **M55 — §2.1.2 method wave: generalised power-potential recogniser.** ✅ DONE. First
  gap-driven method wave on the master §2.1.2 corpus's largest bucket (2nd-order linear). Extends
  `DSolve\`SpecialFunctionForm` (`dsolve_specialform.c`, new `specialform_power_potential` +
  `whittaker_M_1F1`) to two power-potential families the numeric-exponent Bessel row missed, both
  at **symbolic exponent** (the corpus shape) and both emitting heads that NUMERICIZE, so the
  in-method `sf_num_ok` gate and the corpus back-substitution are genuine (0-FAIL by construction):
  - **Single power** `y'' + A x^m y == 0` → `Sqrt[x] Z_{1/(m+2)}(κ x^((m+2)/2))`, `Z = BesselJ/Y`,
    `κ = 2 Sqrt[A]/(m+2)`. Generalises the pre-existing `NumberQ[m]`-gated pure-power row to a
    symbolic `m` (which previously fell to the series fallback and emitted `0^n`/`ComplexInfinity`
    garbage).
  - **Two-term** `y'' + (α x^(2c) + β x^(c-1)) y == 0`. This is NOT Bessel: under `ξ = x^(c+1)` it
    becomes the Coulomb equation `w_ξξ + (c/d)(1/ξ) w_ξ + (α/d² + (β/d²)/ξ) w == 0` (`d = c+1`),
    solved by `y = x^((1-d)/2) WhittakerM[κ, ±μ, z]` with `μ = 1/(2d)`, `κ = β/(2 d Sqrt[-α])`,
    `z = (2 Sqrt[-α]/d) x^d`, emitted as the verifiable confluent form
    `Exp[-z/2] z^(1/2±μ) Hypergeometric1F1[1/2±μ-κ, 1±2μ, z]`. Extraction groups the potential's
    additive terms by their (symbolic) exponent `e = x T'/T` after `Expand` (so a factored
    normal-form potential and un-combined like powers both reduce cleanly) and checks the
    structural constraint `P_big − 2 P_sm == 2`. The `y'`-carrying members (`y'' + a x^n y' +
    b x^(n-1) y == 0`, whose normal form is the two-term family) are reached through the existing
    Liouville normal-form pre-pass (`y = z Exp[-∫P/2]`), whose monomial `z = κ x^d` carries no
    finite-pole radical, so the pre-pass's Whittaker-exclusion NOTE does not apply.
  - **Root-cause verify fix** (`dsolve_common.c`, `ds_residual_numeric_zero`, protects every
    method): the generic verify's numeric KEEP short-circuit declined *any* residual flagged by
    `ds_has_undefined_function`, which is true for a `Derivative` head — including the derivative
    of a DEFINED special function (Bessel/pFq/Airy) produced by differentiating a special-function
    solution. That false decline routed a symbolic-exponent special-function residual to
    `zero_test`, whose precision ladder spins for many seconds (name-sensitively: it hung on the
    exponent symbol `n`/`nn`, was fine on `c`/`k`/`m`). It now declines only inert `Integrate`; a
    genuinely arbitrary forcing `f[x]` still fails to numericize (NaN samples → skipped → the same
    `zero_test` fallback, which is fast for arbitrary forcing). This removed the hang and made 792
    and the `y'`-carrying members solvable.
  - *Anti-overfit:* four forward-generator families (`tests/test_dsolve_m55_stress.c`) — single
    power (a>0/a<0, integer + fractional exponent), two-term Whittaker, `y'`-carrying (pre-pass),
    and an `n`/`nn`-named exponent guard for the verify fix — all `Head === List` + numeric-
    back-substitution verified; pinned unit `t_m55_generalized_power_potential`.
  - §2.1.2: **557 → 564 PASS (+7), 640 non-PASS, 0 FAIL**; gate baseline 655 → 648. The
    2nd-order-linear bucket rises to 239/414 (57.7%). The nine flagship power-potential cases
    (2.1.2-72/81/103/434/791/792/803/804/805) all verify. Remaining apparent flips vs the older
    checked-in report are pre-existing ~6-7 s timing-boundary cases (58/375/439/877) that flip
    P↔U on load and are untouched by M55 (the verify fix keys on a symbolic x-exponent, which
    those trig/radical residuals lack). All DSolve ctest suites (unit + M12/M14/M18/M55 stress)
    green; `make check-c99` green; no regression. Version 0.152 → 0.153.
  - **Kovacic purely-imaginary-pole fix** (follow-on, `dsolve_kovacic.c`): `kovacic_case1_general`
    blanket-declined every non-real pole (a guard against a `ds_simplify(theta)` hang on Heun
    complex-pole equations, §2.2.14-1392/1393), which also broke the ±i case
    `y''−((3+2x²)/(1+x²)²)y==0 → x√(1+x²)` (regression masked for days by the pre-existing
    Risch-Norman unit-suite hang, so `t_kovacic_complex_poles` never ran). Now declines only a
    complex pole with a NONZERO real part (roots of `x²+bx+c`, `b≠0`, whose √-discriminant radical
    is what spins the simplify); a purely-imaginary pair (`x²+c`) is solved. Heun cases still
    decline fast (no hang) → Frobenius. §2.1.2 **564 → 568 PASS, 0 FAIL** (gate 648 → 644);
    §2.2.14 unchanged 99/100. `t_kovacic_complex_poles` passes. Version 0.153 → 0.154.

- **M56 — second-order integrating factors, Stage 2: μ(x, y′) → reduction-of-order first
  integral.** ✅ DONE. Completes the "NEXT INCREMENT" the plan marked open inside M18: the
  μ(x, y′) integrating-factor search (Cheb-Terrab & Roche 1999, Section 2.2, Lemma 3) — which
  M18 had implemented-and-reverted because the reduced first integrals `R==C[1]` are
  non-elementary first-order ODEs — now ships, **unblocked via path (b)**: emit
  `R(x, y[x], y'[x]) == C[1]` as a Maple-style reduction-of-order answer.
  - **μ(x, y′) search (`dsolve_ifactor.c`).** From `Υ = Φ_y`: **Case A** (`∂_{y′}(Υ_y/Υ)≠0`)
    → `𝓕 = 1/(y′-only factors of Υ)`; **Case C** (`Υ_y≠0`, `∂_{y′}(Υ_y/Υ)=0`) → `w` the
    y-only factors, `𝓗=w_y/w`, `p'=𝓗_x/𝓗_y`, `𝓕=(p'+y')w/Υ`; **Case D** (`Υ_y=0`, narrow
    degenerate family) → `p'=Ψ_x`, `Ψ=Φ/Υ−y`, `𝓕=(p'+y')/Υ`. `μ̃(x)` by Lemma 2 (φ₁…φ₄,
    x-only integrand existence). `if_factor_select` (via `FactorList`) does the split; a
    mis-extraction only ever declines. (Case B → A/C; Cases E/F never occur in Kamke — future.)
  - **Correctness gate unchanged.** Each candidate μ reconstructs `R=∫μ dy'+G(x,y)` (the
    existing general `ifactor_build_R`) and must pass `A(R)=R_x+y'R_y+Φ R_{y'}==0` **symbolic
    (`ifactor_R_ok`) + numeric (`ifactor_R_num_ok`, new)** before use → a wrong μ is a clean
    decline, never a wrong answer.
  - **Emit (`dsolve_run_first_integral` + `dsolve_method_builtin_first_integral`,
    `dsolve_common.c`).** New runner parallel to `dsolve_run_implicit` but 2nd-order-aware:
    verifies `d/dx(R)` vanishes modulo `y''==Φ` and assembles `{{ Rf == C[1] }}`; **declines an
    IVP** (one constant cannot fit two conditions). The search is shared by two try-fns — the
    unchanged explicit `dsolve_ifactor_try` (which still wins with a full two-constant closed
    form when the reduced ODE is solvable) and the new `dsolve_ifactor_first_integral_try`
    (reduction fallback), each with its own decline memo. Cascade: new
    `DS_FIRSTINTEGRAL` after `SecondOrderSymmetry`, so full solutions always win.
  - **Solves** Kamke 226 `y''=(x²yy'+xy²)/y'` → `y'[x]²/2 − x²y[x]²/2 == C[1]` (Case A),
    Kamke 136 `y''=(1+y'²)/(x−y)` (Case C), Kamke 66 `y''=a(c+bx+y)(1+y'²)^{3/2}` (Case D).
    New pinned builtin `DSolve\`ReducibleFirstIntegral` (docstring + `ATTR_PROTECTED`); units
    `t_m56_*` (`tests/test_dsolve.c`); anti-overfit `tests/test_dsolve_m56_stress.c` (A/C/D
    forward-generator grids, intrinsic first-integral verify). §2.1.2: **568 → 585 (+17), 0
    FAIL** — the entire net gain is the `2nd_reducible_mu` bucket (38 → 55, gap 64 → 47),
    verified first integrals for 1156 (`_mu_x_y1`), 14/198/897/1157 (`_mu_xy`, all
    `sympy=False`); measurement deterministic across two runs. The 7 P→U flips are the
    pre-existing 8 s timing-boundary cluster (linear/symmetry cases M56 declines; 112 needs
    18.5 s), not regressions. Gate baseline 644 → 631. All DSolve stress suites
    (m5/m12/m14/m18/m55/m56) + `make check-c99` green.
    *Stage 3 (μ(y,y′), the point-swap) and Cases E/F remain future.* Version 0.171 → 0.172.

- **M57 — `SolvableForY` / `SolvableForX`: the `y=G(x,y′)` differentiation method.** ✅ DONE.
  The highest-yield tractable, deterministic §2.1.2 bucket (`1st_solvable_for_yx`) — a capability
  Maple/Mathematica have and Mathilda lacked (Maple's `dp`/`dp2`). `dsolve_solvefor.c`.
  - **Method.** For `F(x,y,y')==0` polynomial in `y` (`SolvableForY`) isolate `y = G(x,p)` (`p=y'`),
    differentiate w.r.t. `x` (`p = G_x + G_p p'`), and recurse the cascade on the induced first-order
    ODE `dx/dp = G_p/(p − G_x)` for `x = X(p,C)`; the general solution is **parametric**
    `{{x->Function[{p},X], y->Function[{p},G(X,p)]}}`. `SolvableForX` is the `x`-mirror
    (`x = H(y,p)`, `dy/dp = H_p/(1/p − H_y)`). This GENERALIZES `DSolve`Lagrange`, which is exactly
    the special case where `y = x φ(p)+ψ(p)` makes the induced ODE linear.
  - **Substrate reuse (zero substrate edits).** The whole parametric path
    (`dsolve_run_parametric` / `dsolve_verify_parametric` / `dsolve_assemble_parametric`) is reused
    verbatim, so every branch is back-substitution verified — **0 FAIL by construction**. The new
    logic is: solve `R==0` for `y`/`x` generally (Lagrange did only R-linear-in-`y`), build the
    induced ODE, and recurse.
  - **Latency.** A cheap `PolynomialQ[R, solve_var]` pre-gate (mirroring `NthAlgebraic`) declines the
    transcendental time-burners (`Tan[x y]`, `Log[Log[y]]`) that never close — `y'==Sin[x y]` went
    20 s → 0.24 s. Bounded exactly as M14 (`TimeConstrained` sub-solves + wall-clock deadline +
    re-entry guard + decline memo + bounded Simplify).
  - **`SolvableForX` is pinned-only** (opt-in, like `FirstOrderPowerSeries`/`EigenvalueProblem`):
    automatic corpus yield ~0, and a degree-1-in-`x` form with a cubic-denominator induced ODE
    produces a transcendental branch whose back-substitution verify is slow (the risk to the Lie
    backstop's time budget). `SolvableForY` is automatic, slotted after every named first-order
    specialist and before the Lie backstop (Maple's late `dp` ordering).
  - **Solves** the flagship `x^(n-1)(y')^n − n x y' + y == 0` (2.1.2-347 family; `y = n x y' −
    x^(n-1)(y')^n` is not affine in `x`, so Lagrange declines), verified n=3..6 to ~1e-16.
    *Anti-overfit:* `tests/test_dsolve_m57_stress.c` (347 + sign-variant 347+ forward-generator grids
    + pinned SolvableForX, parametric numeric verify). Units `t_m57_*` (`test_dsolve.c`, verified in
    isolation — pre-existing `t_m19` in-suite abort). §2.1.2 (clean re-run): **585 → 590 (net +5),
    0 FAIL**; the **deterministic** M57 gain is **+2** (347 `y=_G(x,y')`, 352 `[F(x),G(y)]`-symmetry,
    both via `SolvableForY`), the other +3 being the `_with_linear_symmetries` timing cluster
    (2nd/high-order, declined instantly by `SolvableForY`) oscillating; 350/351 solve but exceed the
    8 s forked budget. Gate baseline 631 → 629.
    All DSolve stress suites (m5/m12/m14/m18/m55/m56/m57) + `make check-c99` green. Version 0.172 →
    0.173. *Future:* explicit (non-parametric) `p(x)` route; parametric IVP fitting; parameter
    elimination to an implicit `Φ(x,y)=0`; the `y⁽ⁿ⁾`-solvable higher-order generalization; the
    radical / degree-≥4 bucket cases needing `Root`-object handling.

- **M58 — higher-order autonomous reduction (missing-x, order n ≥ 3).** ✅ DONE.
  Lifted `DSolve`AutonomousReduction` (`dsolve_autonomous.c`) from order-2-only to any
  order n ≥ 2, completing the missing-x / missing-y reduction pair (the missing-y half at
  order 3+ was already `LowerDerivativeReduction`). For `y⁽ⁿ⁾ == f(y,…,y⁽ⁿ⁻¹⁾)` (no explicit
  x), the derivative chain `D₁=p, D_{k+1}=p·d/dy(D_k)` (`y''=p·p_y`, `y'''=p²·p_yy+p·p_y²`, …)
  turns the ODE into an **(n−1)-order ODE in p(y)**, solved by recursion, then the separable
  `y'==p(y)` closes it. Constants `C[1..n−1]` from stage 1 are frozen to `C[2..n]` via
  `dsolve_renumber_constants` before stage 2 mints its fresh `C[1]`. For n=2 the chain
  reproduces the classical `p·p_y == f(y,p)` exactly (order-2 behavior unchanged; the
  Tan/Tanh and exponential families still solve, the elliptic case still declines).
  - **Zero cascade/wiring changes** (the method was already registered); the change is
    confined to `dsolve_autonomous.c` internals + reuse of `dsolve_renumber_constants`.
  - **Stage-2 is the yield gate.** Stage 1 (the reduction) always closes, but the separable
    stage-2 quadrature `∫dy/p(y)` is **non-elementary for most order-3 corpus cases** (a
    `Sqrt` of a `Log`, a high-degree/hyperelliptic radicand, or a rational-under-radical),
    and Integrate SPINS uninterruptibly on those (`TimeConstrained` cannot preempt it). A
    **stage-2 spin-guard** declines before that spin — when `p` carries a `Log` or a
    y-dependent denominator — keeping the elementary-quadrature cases (`y y'''==y'y''` →
    `p=Sqrt[C₁y²+C₂]`, elementary) and the order-2 forms. A per-top-level decline memo +
    a 5 s wall-clock deadline bound the recursion.
  - **Yield:** §2.1.2 clean re-run **590 → 591 (deterministic +1), 0 FAIL** (`1143`
    `y y'''==y'y''`); the other Group-A autonomous 3rd-order cases (263/264/267/268/1167/1168)
    reduce correctly but their stage-2 quadrature is non-elementary, so they decline (matching
    what an elementary-quadrature engine can do). The net P↔U movements are the
    `_with_linear_symmetries` timing cluster oscillating. Gate baseline 629 → 628. The capability
    is nonetheless real — Mathilda solved NO 3rd+ order autonomous ODE before, and future Integrate
    improvements extend the set automatically. Version 0.173 → 0.174. *Future (the big lever):* return the **implicit first integral** `∫dy/p(y)==x+C`
    (inert quadrature, verified by implicit differentiation) for non-elementary stage-2, which
    would unlock the full Group-A (~+7) — matching Mathematica — and is a substrate-level task.

- **M59 — `Inactive` primitive + implicit-first-integral autonomous companion.** ✅ DONE.
  Realises M58's deferred "big lever". Two parts.
  - **`Inactive` / `Activate` primitives** (`src/core.c`, `src/sym_names.{c,h}`, `src/calculus/deriv.c`).
    `Inactive[f]` is an inert head-wrapper: `Inactive[f][args]` evaluates its args but does not fire
    `f`'s rules, so `Inactive[Integrate][g,x]` holds the integral WITHOUT running the (uninterruptibly
    spinning) integration cascade — the key finding was that this compound-head form is **already a
    fixed point** (like `Derivative[n][f][x]`), so registration is just `ATTR_PROTECTED` and inertness
    is automatic. The only real additions: a **`D`-FTC rule** `D[Inactive[Integrate][f,u],u]==f`
    (`deriv.c`, at the compound-head dispatch — the different-variable case is handled by the generic
    free-of-x short-circuit), and **`Activate`** (`Inactive[h]→h` + re-evaluate).
  - **`dsolve_autonomous_implicit_try`** (`dsolve_autonomous.c`). Where the M58 explicit method
    declines a non-elementary stage-2 quadrature, this returns the inert first integral
    `Inactive[Integrate][1/p, y[x]] − x == C[1]` via `dsolve_run_implicit` — whose verify computes
    `y'=−G_x/G_y=p` through the new FTC rule with no integration. Correctness is a **numeric
    self-verify** (reconstruct `y'..y⁽ⁿ⁾` from the reduction chain, check the original ODE ≈0), since
    `dsolve_run_implicit`'s verify passes vacuously at order n≥2. The M58 stage-1 reduction is factored
    into a shared `ar_reduce`; the explicit method (order-2, 1143) is unchanged. Wired as a second
    cascade slot + explicit-then-implicit pinned method (mirroring `Homogeneous`).
  - *Solves:* the order-3 autonomous residue 263/264/267/268/1167/1168 (`p=Sqrt[…]`, + dups
    710/711/712) now returns a verified inert first integral instead of declining;
    `y'==Sqrt[y Log y+…]`'s 45 s spin → 0.0 s. §2.1.2 clean re-run **591 → 595 (deterministic +9),
    0 FAIL** (`3rd_high_reducible` bucket 8 → 17 PASS; the net P↔U is the x-dependent timing cluster);
    gate baseline 628 → 619. New `tests/test_inactive.c` + M59 units/stress;
    core/deriv/integrate regression suites unaffected. Version 0.174 → 0.175. *Future:* upgrade the
    `separable`/`fos`/`chini`/`exact` implicit paths to emit inert first integrals for their own
    non-elementary quadratures (they currently decline), now that the mechanism exists.

- **M60 — higher-order linear: reducibility completion + generalised-Airy recogniser.** ✅ DONE.
  Gap-driven wave on the largest measured §2.1.2 bucket, `3rd_high_linear` (142 total, 41 PASS,
  101 UNEVAL). Probing the binary against the real corpus equations split that bucket into three
  independently addressable mechanisms plus one outright defect, and M60 lands all of them. Every
  new branch is gated by an in-method numeric back-substitution **and** the substrate's own verify,
  so the 0-FAIL invariant holds by construction.
  - **`DSolve\`GeneralizedAiry`** (new `dsolve_genairy.c`) — the n-th order **pure-power
    potential** `u^(n) == A x^m u`, the higher-order analogue of Airy's equation, whose fundamental
    set is `x^j ₀F_{n−1}(; {1 + (j−i)/p : i ≠ j}; A x^p/pⁿ)` for `j = 0..n−1`, `p = m + n`. Derived
    from the balance `x^{j+pk−n}` against `A x^{m+j+p(k−1)}`: the `i == j` factor of the recurrence
    is exactly the `k!` of a `₀F_{n−1}` and the other `n−1` are its lower parameters. The equation
    reaches that shape through a **depression (gauge) pre-pass** generalising `SymmetricSquare`'s
    order-3 gauge to arbitrary `n`: with `v = −c_{n−1}/n` and `W_0 = 1, W_{i+1} = W_i' + v W_i`
    (Bell polynomials in `v`, so no `Exp` ever has to cancel), `L[w u]/w = Σ_j d_j u^(j)` with
    `d_j = Σ_{k≥j} c_k C(k,j) W_{k−j}`; `d_n == 1` and `d_{n−1} == 0` by construction, and when
    every remaining intermediate `d_j` vanishes the potential is read off `−d_0` by the M55
    `m = x T'/T` extraction (so a **symbolic** `A` and `m` work — 2.1.2-595 is `y''' == a x^b y`).
    Forcing is added by variation of parameters. Declines `p == 0` (Euler), a lower parameter that
    is a non-positive integer (exponents differing by a multiple of `p` — the logarithmic Frobenius
    case, which is why `x w''' + w == 0` and `x² y'''' == A y` correctly fall through), and any
    non-power potential. Cascade slot: after `SpecialFunctionForm`, before `SymmetricSquare`, gated
    to order ≥ 3. *Solves* 2.1.2-230/292/311/1201 (`y''' == x y` → `₀F₂({},{1/2,3/4},x⁴/64)` and its
    `x`, `x²` partners), -1073, -595 (symbolic), and -241/-604 through the gauge (`w = 1/x`).
  - **`OperatorFactor` at order 2 + forcing.** The method was gated to order ≥ 3 and homogeneous.
    Both gates are lifted. Order 2 matters because `Kovacic` runs first and *owns* the tidy answers
    there, so what reaches `OperatorFactor` is the rational-Riccati Case-1 residue Kovacic declined —
    and a closed form from it beats the truncated Frobenius series that was winning. That is exactly
    what unblocks the order-3 peel: 2.1.2-253's quotient `x²(1+x)u'' + 2x(2+x)u' + 2u == 0` has the
    hyperexponential solution `u = 1/x`, which `DSolve\`DFactor` found all along while `DSolve`
    returned `O[x]^6`. Forcing is carried by handing the monic-normalised right-hand side to the
    recursive quotient solve, so the `_linear,_nonhomogeneous` members ride along.
  - **Adjoint (left-factor) peel = the "Beke / 2nd-order right factor" target.** When no first-order
    RIGHT factor exists, the same search runs on the adjoint `L* = Σ (−1)^k D^k ∘ a_k`. Since
    `(A∘B)* = B*∘A*` and `(D − s)* = −(D + s)`, a first-order right factor of `L*` is a first-order
    LEFT factor `(D + s)` of `L`, i.e. an order-(n−1) RIGHT factor `Q` — the classical case for
    `n = 3`, reached **without exterior powers**. Left division `q_{n−1} = a_n`,
    `q_{m−1} = a_m − q_m' − s q_m`, remainder `a_0 − q_0' − s q_0`, is the exact acceptance test, so
    a mis-found `s` can only decline. `Q[z] == 0` is solved by recursion and `L[y] == g` closed by
    variation of parameters over `Q`'s basis for `W = Exp[−∫s](∫ g Exp[∫s] dx + C[n])`.
    `DSolve\`DFactor` reports the left factor too, emitted **last** (the list is innermost-first) —
    e.g. `(D + 1/x) ∘ (D² − x)` → `{Dx² − x, Dx + 1/x}`, an order-2 right factor of an operator with
    no first-order right factor at all. **Measurement corrected the roadmap here:** a `Q` reachable
    only through the adjoint is by definition free of hyperexponential solutions, so its basis is
    Kovacic-case-2 / special-function, and the VoP particular is then essentially never elementary —
    2.1.2-294's quotient is `Sqrt[x] BesselJ[I Sqrt[3], 2 Sqrt[x]]` and the integrals do not close.
    So the left peel's *factorisation* is the real capability and its *solve* yield is small; the
    named next step is to emit the particular with M59's `Inactive[Integrate]` instead of declining.
  - **Generalised Bessel row** (`dsolve_specialform.c`) — `Q = A x^m + B x^(−2)` →
    `Sqrt[x] Z_ν(κ x^((m+2)/2))`, `ν = Sqrt[1−4B]/(m+2)`, `κ = 2 Sqrt[A]/(m+2)`. A strict
    generalisation of the single-power row (`B == 0` gives the same `ν = 1/(m+2)`) that also reaches
    the inverse-power pair `{−1,−2}` the Whittaker condition `P_big − 2 P_sm == 2` misses. Complex
    `ν` is admissible: `BesselJ`/`BesselY` numericize there, so `sf_num_ok` stays genuine.
  - **Three latency root-causes, all pre-existing but newly exposed by the wider reach.** (1) The
    trailing integrand was built as a raw `Exp[−(Log[x] + Log[1+x²]/2)]`, which sends `Integrate`
    down the transcendental-tower path: 3.98 s **and it fails**, where the same integrand simplified
    to `1/(x Sqrt[1+x²])` is algebraic, closes, and costs 0.015 s. The gauge factor is now simplified
    — and then *checked* (`w' == r w`) rather than trusted, falling back to the raw form if the
    simplifier picks a branch that breaks it. (2) The algebraic integrator is two orders of magnitude
    slower on an integrand carrying symbolic `C[k]` coefficients (2.06 s vs 0.011 s), so the trailing
    integral and the left peel's VoP are now split along the constants by linearity. (3)
    `of_find_factor` bounded the pole *count* but not the ansatz *width*: two double poles plus a
    degree-2 polynomial part is 11 unknowns in a cubic determining system, which does not return.
    New `OF_MAX_UNKNOWNS` bound (every factor the search has ever found needed ≤ 4), and the ansatz
    pole order is capped at 1: a DOUBLE pole in `r` is an irregular singularity, which the file's
    documented scope already defers, so searching for one only doubled the width of every
    combination — it earned nothing on the corpus (all 30 gains verified individually at pole order
    1) while costing the most (a failing search on a two-pole order-3 operator: **60.2 s → 0.72 s**).
    Together these took 2.1.2-250 from a >120 s non-answer to a **2.6 s solve**. A structural gate also declines the
    left peel when the quotient basis is a series or a special function, where the VoP integrals
    churn for seconds without closing (2.1.2-294: 11.9 s → 0.37 s).
  - *Tests:* `tests/test_dsolve_m60_stress.c` — eight forward-generator families (pure power over an
    (n, m, a) grid; the gauge family `D[x^k y, {x,n}] == a x^(m+k) y`; symbolic exponent with a
    post-instantiation numeric check; order-2 and order-3 factorable operators over a rational `r`
    grid, requiring a closed form and not a `SeriesData`; the adjoint family, whose returned
    factorisation is **reconstructed and compared against the original operator**; forcing; and
    bounded declines). Units `t_m60_*` in `tests/test_dsolve.c`.
  - §2.1.2: **605 → 634 PASS (+30 gained, 1 lost, net +29), 570 non-PASS, 0 FAIL**, measured against
    a fresh same-day baseline on the same tree (the checked-in M59 row of 595 had drifted +10 on
    intervening non-DSolve work, so it is not a valid comparator) and **deterministic across two full
    runs**. Bucket attribution: the targeted **`3rd_high_linear` 42 → 64 (+22)** (29.6% → 45.1%),
    `2nd_linear` 243 → 249 (+6), `1st_solvable_for_yx` +1, `Emden_Fowler` +1; the one loss (342) is a
    first-order *nonlinear* ODE `OperatorFactor` declines in 0.4 ms and which solves in 5.25 s
    standalone — the documented 8 s timing cluster. Gate baseline **619 → 580**. All 14 DSolve ctest
    stress suites pass (including the new `dsolve_m60_stress_tests`, 19 s) and `make check-c99` is
    green. *Note:* `dsolve_tests` (the unit suite) is **pre-existing red** on this machine — it hits
    the `alarm(120)` self-kill in `test_utils.h` at test 48 of 266, so the `t_m60_*` group registered
    at the end never executes there; all 16 of its assertions were verified by direct evaluation, and
    the coverage that actually runs is the stress binary. *Future:* emit the left peel's
    variation-of-parameters particular with M59's `Inactive[Integrate]` when it is non-elementary
    (2.1.2-294/296 factor correctly but their `Sqrt[x] BesselJ[I Sqrt[3], 2 Sqrt[x]]` quotient basis
    gives no elementary particular), and the genuine `m`-th exterior power for right factors of order
    `2 ≤ m ≤ n−2`.

- **M61 — §2.2.34 corpus (Problems 3301–3400) + inert VoP particular, autonomous IVP constant
  fit, separable singular solution, two latency root causes, nonlinear Clairaut.** ✅ DONE.
  12000.org's next hundred, measured end to end and then driven up by fixes that are general
  rather than case-shaped: **80/100 → 91/100 (+11), 0 FAIL, 0 crash**, deterministic across two
  per-case-identical runs, every gain attributable to a named defect. New gate
  `dsolve_corpus_2_2_34_tests` (baseline 9); report `DSolve_test_status/reports/2.2.34.md`.
  - **Upstream renumbered.** The site regenerated 2026-09-28 and SWAPPED its chapter-2 section
    numbers: the sequential pages moved §2.2.N → **§2.1.N** (`Ch2.S1.SSN.htm`), and the master
    corpus this plan calls §2.1.2 is now upstream §2.2.2, paginated `Ch2.S2.SS2.SSS1…13.htm`.
    Internal names keep the `2.2.34` spelling (continuity with 33 sections and every milestone
    here); the mapping and the consequences are in `DSolve_test_status/README.md`. The converter
    needed no format work, *verified* not assumed: re-converting §2.2.33 from the new page gives
    the same record counts and **100/100 semantically identical** records (the 90 textual diffs
    are upstream LaTeX reformatting). One lasting consequence: the converter's byte-identity
    contract can no longer be exercised against the live site, so **semantic** equivalence
    replaces it.
  - **`DSolve\`VariationOfParameters` — the inert `Inactive[Integrate]` particular.** This is the
    step M60's closing note named twice, and measurement changed its shape entirely: the
    mechanism already existed (the shared VoP helper keeps a non-closing Wronskian integral) and
    the real blocker was a missing call site — `dsolve_specialform.c` is **homogeneous-only by
    construction** (it gates on `dsolve_second_order_PQ`), so for a 2nd-order equation at a
    regular singular point with polynomial forcing it found the Bessel/hypergeometric
    fundamental set and threw it away because the equation had a right-hand side, while
    `Kovacic` declines (Bessel is not Liouvillian). A second mode (`nh_try_core(…, inert)`)
    occupies a new **LAST** cascade slot, after both Frobenius fallbacks. The slot position is a
    correctness property, not a convenience: from there it can only turn UNEVAL into an answer,
    so no existing PASS can be lost to it and its latency lands only on equations that were
    returning nothing. The homogeneous part is solved by the **pinned** `SpecialFunctionForm`,
    which structurally guarantees a closed-form basis — so a truncated-series basis, over which
    an inert integral is meaningless, can never reach the quadrature (that is the 3393/3394 gate,
    held by construction rather than by a post-hoc `SeriesData` check).
    Two mechanisms carry it. `vp_integral_hopeless` **never hands `Integrate` an integrand whose
    DENOMINATOR carries a special function** — the failing attempt on a Bessel Wronskian quotient
    costs up to **47.8 s**, which blows every solve budget, while the screen costs microseconds
    (and it is a denominator test, not a blanket one: `Integrate` genuinely closes a special
    function in the numerator). And `ds_inert_vop_verified` is the **only** correctness barrier,
    because nothing downstream can catch a wrong inert answer: `ds_residual_numeric_zero` and
    `ds_branch_num_ok` both bail to KEEP on an `Integrate` head (and `ds_has_head` is a *name*
    test, so it cannot even tell the inert head from the active one — hence the new
    `ds_has_active_integrate` / `ds_has_inactive_integrate`), and `PossibleZeroQ` answers `True`
    for any residual containing an inert integral. The gate exploits the particular being LINEAR
    in its inert integrals: replace each by a symbol, split with `Coefficient`, and require every
    coefficient of the linear form to vanish — the FTC part must cancel the forcing and each
    `L[basis_i]` must be zero. Sampling is at **exact rationals** (`ds_subst_generics_exact`):
    the same Bessel residual reads 1e-7 at machine reals, where twelve digits go to cancellation
    between terms of magnitude 1e4, and 1e-28 exactly. An **IVP declines** — an inert particular
    has no value at a point, so its constants cannot be fitted, and the prelude would score the
    unfitted general solution as solved. *Solves* 3387/3388/3389/3392/3395. `deriv.c` also gains
    the **Leibniz rule for the inert DEFINITE integral**, the companion of M59's indefinite rule.
  - **`AutonomousReduction` — fit the stage-1 constants from the point conditions.** For an IVP
    the conditions determine them exactly and independently of the quadrature
    (`D_k(y0) == y⁽ᵏ⁾(x0)`, the same chain the reduction is built from), and fitting them FIRST
    is what closes the class: `y''+2yy'==0, y(0)=0, y'(0)=1` reduces to `p == C[2]−y²`, whose
    quadrature is a symbolic-parameter `ArcTanh` the stage-2 spin guard turns away — with
    `C[2]==1` fitted it is `ArcTanh[y]` and the answer is `Tanh[x]`. The value condition is then
    handed to the stage-2 sub-solve too, so the body comes back constant-free (otherwise the
    substrate must invert `{Tanh[C[1]]==0, Sech[C[1]]²==1}`, which it cannot decide, and declines
    the correct branch). Where the quadrature stays non-elementary the relation is fitted all the
    same, in **definite** form `Inactive[Integrate][1/p, {t, y0, y[x]}]`: the indefinite form
    could not be fitted at all, because `dsolve_implicit_rhs`'s `y[x] -> y0` is a blind
    `ReplaceAll` that also rewrote the integration-variable slot, producing a meaningless
    `Inactive[Integrate][1, 0]` which — once no free constant remained — would have scored as a
    solved answer. A radical over a transcendental function of `y` now declines the explicit path
    fast (that stage-2 sub-solve costs 90 s and returns an inert relation anyway), where the
    existing denominator test could not see it. *Solves* 3345/3347; 3346 stays elliptic.
  - **`Separable` — the singular (equilibrium) solution of an IVP.** Dividing by `h(y)` drops the
    constant solutions, so `y'==x²y², y(1)=0` shipped `1/(C[1]−x³/3)` with its constant unfitted.
    Placed in the substrate (`dsolve_run`), keyed on **every surviving branch being FIT_EMPTY** —
    i.e. `Solve` PROVED the family cannot reach the initial point — and not on `h(y0)==0`, which
    also holds for `y'==y, y(0)=0`, `y'==y(1−y), y(0)=1` and `y'==Sqrt[y], y(0)=0`, all three
    already correct (testing it would emit a duplicate branch on two and a Mathematica-divergent
    extra branch on the third). Covers Bernoulli/Homogeneous/Chini/Abel too. *Solves* 3336.
  - **Two latency root causes, both with reach far beyond this corpus.** (1)
    `dsolve_verify_parametric` substituted an **uncancelled** `dY/dX` into a residual that raises
    it to a power (`x y'³ == y y' + 1` cubes it), and `zero_test`'s canonicalisation of the
    degree-exploded rational did not return — a **>120 s hang** reachable from every parametric
    answer (`Lagrange`, `SolvableForY`, `SolvableForX`), i.e. any user ODE of d'Alembert shape.
    `Together` on the quotient before substitution collapses it to `t`. **Gated to a quotient
    rational in the parameter**, which the wave's own measurement forced: ungated, `Together` is
    itself the expensive step on a radical-carrying candidate (2.1.2-980, `(x²+y²)^(3/2)`: 6 s →
    51 s plus a recursion-limit blowup). Cancelling the residual instead does not work — it has
    to happen before the power. (2) `NthAlgebraic` is the **second** method tried on every scalar
    ODE and its per-root recursion was an unbounded `DSolve` whose *implicit* results it then
    discards: on `(y − x y')² == 1 + y'²` that was two 5.9 s solves thrown away, repeated by the
    evaluator's fixed-point re-invocation for ~25 s, on an equation `Clairaut` answers in 0.01 s.
    It now carries the standard kit (wall-clock deadline, per-branch `TimeConstrained`, decline
    memo). *Solves* 3312 and, with the next item, 3331.
  - **`Clairaut` — a Clairaut equation written NONLINEARLY in `y`.** `(y − x y')² == 1 + y'²` is
    textbook Clairaut (its roots `y = x p ± Sqrt[1+p²]` are two Clairaut equations) and was
    reachable by no method: this one required the algebraic residual to be *linear* in `y`, and
    `SolvableForY`, which does isolate the roots, discards them deliberately — for a Clairaut
    equation its denominator `p − G_x` vanishes identically, and its comment says the family is
    "owned earlier", which it was not. The residual now need only be **polynomial** in `y`
    (`PolynomialQ` gate, so a transcendental residual never reaches `Solve`); each root runs
    through the same `d/dx Yexpr == p` test via the factored `clairaut_emit`. The linear path is
    byte-identical, envelope included.
  - *Tests:* six `t_m61_*` units in `tests/test_dsolve.c` and `tests/test_dsolve_m61_stress.c` —
    seven families: an **18-member forward generator** over the generalised-Bessel grid
    `x²y''+xy'+(a x^m+b)y == p(x)` (whose homogeneous basis is Bessel by construction, so nothing
    is hand-picked), an elementary-still-wins family, the **gate accept/reject margin measured on
    planted wrong bases** (correct 5.4e-51 against 0.18, 0.18 and 0.49 for an extra `Sqrt[x]`, a
    mismatched Bessel order and a wrong exponent — the margin is *tested*, not merely observed),
    series-basis and IVP declines, a **latency bound** that notices if the denominator screen
    stops firing, and bounded declines. Each verified independently of the in-method gate, at
    different sample points and precision.
  - *Regression:* measured as a **same-machine A/B against a HEAD binary built in an isolated git
    worktree** — the only honest comparator, since the checked-in per-section reports are stale
    and several section gates are **already red on main today**. Over 11 exposed sections
    (separable/IVP-heavy, 2nd-order-linear-heavy, dAlembert-heavy): M61 **better on 6, equal on
    5, 0 FAIL in all 22 runs** (2212 5→4, 2214 3→1, 2219 11→10, 2225 4→2, 2227 5→2, 2233 18→12).
    All 15 DSolve unit/stress suites report no failure; `make check-c99` and `check-messages`
    green. `dsolve_tests` remains **pre-existing red** (it hits `alarm(120)` in `test_utils.h`
    before any test runs, so the `t_m61_*` group never executes there) — its 30 assertions were
    verified by direct evaluation instead, all passing; two stress suites truncate at the same
    alarm, and their stalling families were A/B-confirmed identical on the HEAD binary.
  - *Future, named by this wave:* the **series particular** — `FrobeniusSeries`/`PowerSeries`
    both call the homogeneous-only extractor and so decline ANY forced equation, which is exactly
    what loses 3393/3394 and is the mechanism this whole "series expansion" block is about. And a
    **nonlinear Taylor-series IVP** method would close the six no-closed-form cases
    (3319/3338/3341/3342/3348/3349) but only as PASS-on-trust, since the verifier samples far
    from `x0`; it belongs in its own milestone behind an `O[(x−x0)^N]` residual gate, not bolted
    on here. Also unclaimed: `SolvableForX` is pinned-only and solves 3319 in 0.088 s — the
    documented objection to admitting it to the cascade was the slow parametric verify this wave
    just fixed, so that decision is now worth re-measuring. Separately, and larger than this
    wave: **`TimeConstrained` refunds nested time** (`builtin_time_constrained` saves and restores
    the outer `ITIMER_PROF` instead of clamping to it), so every nested sub-solve extends the
    user's deadline — which is why a `TimeConstrained[…, 12]` here returned at 20 s. That is a
    global fix (Integrate, Simplify, NIntegrate, the corpus harness itself) and deserves its own
    change and measurement.

- **M62 — §2.2.35 corpus (Problems 3401–3500) + `TimeConstrained` nesting, symbolic-exponent
  incomplete Gamma, sequential constant fit, nonlinear exact ODEs.** ✅ DONE. The next hundred,
  measured end to end and then driven up by four fixes that are general rather than case-shaped:
  **95/100 → 98/100 (+3), 0 FAIL, 0 crash, 0 timeout**, deterministic across two per-case-identical
  runs. New gate `dsolve_corpus_2_2_35_tests` (baseline 2); report
  `DSolve_test_status/reports/2.2.35.md`. Two of M61's closing notes are discharged here: the
  global `TimeConstrained` fix it named, and the latency half of its series-particular item.
  - **`TimeConstrained` CLAMPS a nested budget instead of refunding it** (`src/core.c`) — the
    global fix M61's closing paragraph named, and the largest single lever here. Both enforcement
    layers lifted the caller's deadline merely by entering an inner scope: the `ITIMER_PROF` layer
    reinstalled the outer timer at the value it held when the inner call *started* (a full refund
    of everything the inner call spent, and the inner was armed for its full request even when
    that exceeded the outer's remaining time), and the cooperative wall-clock layer set an inner
    ABSOLUTE deadline that could be later than the outer's — plus it dropped the outer deadline
    entirely when `clock_gettime` failed, i.e. on exactly the hosts that layer exists for. Both
    had to be fixed; either alone breaks the property. The repro is two lines and needs no
    subsystem: `TimeConstrained[TimeConstrained[<loop>, 30], 3]` ran the loop **to completion in
    24.5 s**. Because Integrate (three call sites), Simplify, NIntegrate and DSolve's own
    per-method kit all bound sub-steps this way, *no* budget in the system was an upper bound —
    which is why the measured cost of a failing attempt could be 46 s under a 3 s constraint. Now
    `min(inner, outer remaining)` is armed and the outer is reinstalled charged for what the inner
    consumed (1 µs when exhausted, so the caller aborts at its next step). Consequences beyond
    this corpus: `DSolve\`Kovacic`'s forced closure on 3401 went **375 s → 9.85 s**, and the two
    §2.2.32 cases documented for waves as "correct VoP at the 8 s boundary" (3164/3165) are now
    inside it. A companion, smaller lever landed with it: `dsolve_vp_set_integral_budget` lets a
    variation-of-parameters caller that needs an ELEMENTARY closure bound each Wronskian integral,
    which `DSolve\`Kovacic` now does (its own Liouvillian basis is invisible to the structural
    `vp_integral_hopeless` screen — `Exp[c ArcTanh[radical]]` carries no special-function head yet
    still sends `Integrate` on a search that does not return). A timed-out term must NOT come back
    as a raw unevaluated `Integrate`: that re-enters the integration cascade on every later
    re-evaluation of the body, which turned a 3 s budget into 375 s over the re-evaluations — so
    the elementary mode bails out of the particular instead, and Kovacic declines a concrete
    forcing whose integral did not close rather than handing `numeric_verify` a residual it cannot
    decide.
  - **Symbolic-exponent `x^p E^(a x^m)` → incomplete Gamma** (new
    `src/calculus/integrate_gammapower.c`, one cascade line in `integrate.c` after the Fresnel
    recogniser). With `p` symbolic the integrand is outside every elementary stage, and those
    stages do not merely decline — they SEARCH: `Integrate[x^n E^(-x), x]` cost a measured
    **12.9 s** to come back unevaluated, and 12.4 s for the Gaussian sibling `x^n E^(-x²)`. Two of
    those are the entire 26 s cost of 3495 `y'' - y == x^n`, whose variation-of-parameters
    particular is exactly that pair of integrals. The answer is
    `-(1/m)(-a)^(-s) Gamma[s, -a x^m]` with `s = (p+1)/m`, emitted only behind an exact
    differentiate-back certificate — two passes, because the closed form is stated in the
    principal-branch convention `(-a x^m)^k == (-a)^k x^(m k)` which `Simplify` does not apply on
    its own, so `a > 0` or `m > 1` needs `PowerExpand`; the second pass quotients out exactly that
    convention and nothing else. Gated to a symbolic exponent: for a non-negative integer `p` the
    answer is elementary and the existing stages give it, and for a numeric non-integer they give
    the cleaner Erf form, so firing there would be a quality regression dressed as a speedup. The
    gate needed a `NumberQ` test, not a structural symbol scan — `Rational[1,2]` is an
    `EXPR_FUNCTION` with a symbol head, so `Sqrt[x]`'s exponent reads as "symbolic" to a bare
    walk, and the recogniser silently took over `Integrate[Sqrt[x] E^(-x), x]`. *Solves* 3495
    (8 s abort → 0.09 s, closed form).
  - **Sequential scalar constant fit** (`ds_fit_sequential`, `src/calculus/dsolve_common.c`). Only
    Solve's SCALAR form applies inverse-function inversion — the substrate already relied on that
    for the single-condition case — and every multi-constant fit uses the LIST form, so a constant
    nested inside a transcendental bubbles back unevaluated and takes the whole method down with
    it. 3482 `y'' + y'² + y' == 0, y(0) = 0` is the small example: the general solution
    `C[2] + Log[C[1] - E^(-x)]` is found in 30 ms, each condition is individually invertible, and
    the IVP was declined. The fallback fits one condition at a time in the scalar form, accepting
    a substitution only when the constant COUNT drops (so a bubbled scalar Solve cannot smuggle an
    unfitted body through), and leaves surplus constants free — correct for an under-determined
    IVP, and tested as such: one condition on a second-order equation must leave exactly one free
    constant, not zero and not two. It runs only after the list form produced no fit, so no fit
    that already worked can change. *Solves* 3482.
  - **`DSolve\`ExactODE`: a nonlinear total derivative, and a per-level first-integral constant**
    (`src/calculus/dsolve_exactode.c`). Two independent repairs in one file.
    *(a) An incomplete general solution that scored PASS.* The first-integral constant is carried
    as a plain private symbol rather than a `C[k]` (a documented `DiffUnderInt` avoidance), and it
    was ONE fixed name — so on a **doubly**-exact equation the outer and inner constants were the
    same symbol and merged inside the inner sub-solve, before any rename could tell them apart.
    `x y''' + 2 y'' == A x` came back as `C[1] + C[2] x + A x³/18 + C[2] Log[x]`: two constants
    for a third-order ODE, with `x` and `Log[x]` sharing one. Nothing downstream can catch that —
    the residual of an incomplete family is still exactly zero, so the corpus harness scored 3498
    as solved. Fixed with one name per nesting level plus numbering the constant one past the
    largest `C[k]` the sub-solve actually used (identical to the old `C[n]` whenever the recursion
    does not nest). The regression test is a COUNT: an order-*n* general solution must carry
    exactly *n* distinct constants.
    *(b) Nonlinear exactness.* The linear path matches coefficients; a nonlinear left side needs
    the general construction over the jet variables `j_k = y^(k)`. Since
    `dF/dx = F_x + Σ F_{j_k} j_{k+1}`, we have `dL/dj_n = F_{j_{n-1}}`, so one integration recovers
    `F`'s dependence on the top jet and subtracting that piece's total derivative leaves a shorter
    expression to which the same step applies; peeling `k = n … 1` either exhausts `L` or leaves a
    remainder still carrying a derivative, which PROVES non-exactness. One exact `dF/dx == L`
    comparison is the acceptance certificate. 3497 `2 y y''' + 2(y+3y')y'' + 2 y'² == Sin[x]`
    peels to `2 y y'' + 2 y y' + 2 y'²` — which is `(y²)'' + (y²)'` — and the recursion closes the
    chain. Note the certificate is self-referential if the operator part is wrong: an early cut
    formed `R + g0` instead of `R - g0` and verified happily against its own wrong `L`, producing
    a first integral off by `2 Cos[x]`; that is why the sign is now commented.
    **Placement cost two measured corrections.** The peel is a general backstop — it claims ANY
    equation whose left side happens to be a total derivative — so sharing ExactODE's early linear
    slot made it preempt methods that answer the same equations directly. It now has its own LATE
    slot (after `ifactor_first_integral`, before the series fallbacks) and runs for a general
    solution only (`ncond == 0`). Both restrictions came from the corpus, not from reasoning:
    ungated for IVPs it claimed and lost §2.2.34-3345/3347, which `AutonomousReduction` answers
    `Tanh[x]` by fitting its stage-1 constant from the point conditions; and from the early slot it
    took **§2.1.2-1143** (`y y''' == y' y''`, whose first integral is `y y'' - y'^2`) from 0.08 s to
    9.96 s for the *same answer* — a 125x latency regression only the master-corpus run found.
    Pinning `DSolve\`ExactODE` tries both paths, so the pinned method stays complete.
    *Solves* 3497, and §2.2.34-3346 (the elliptic stage-2 case M61 left as residue) as a
    side-effect.
  - *Withdrawn, deliberately:* a **hygiene gate** rejecting any answer that carries one of the
    solver's own private `DSolve\`` symbols. The leak it targets is real and is a silent wrong
    answer — `DSolve\`Y` escaping the Bernoulli linearisation reads as a free parameter, so
    neither the symbolic nor the numeric verify rejects it — but a blanket test is wrong: M61's
    inert DEFINITE integral legitimately carries its bound integration variable `DSolve\`impT` in
    the answer, and the gate cost **nine** §2.2.33 passes. The correct version needs
    bound-variable analysis (a private symbol is a leak only where it occurs FREE) and is its own
    change. Built, measured, removed — recorded here because the measurement is the useful part.
  - *Tests:* `tests/test_dsolve_m62_stress.c`, five families, each shaped to the failure mode it
    guards. The `TimeConstrained` family asserts an INEQUALITY on a pure CPU loop at one, two and
    three levels of nesting (no subsystem heuristics in the way). The constant-count family runs a
    forward generator of doubly-exact equations and counts distinct `C[k]` — the only test that
    can see an incomplete general solution. The Gamma family verifies by the recogniser's own
    certificate over a `(p, a, m)` grid, and its NEGATIVE controls (integer and rational exponents
    must keep their elementary / Erf answers) matter more than its positives; it also carries a
    latency bound, since a recogniser that stops firing shows up as 13 s rather than as a wrong
    answer. The fit family generates nonlinear second-order IVPs over `y'' + y'² + k y' == 0`,
    checks the ODE *and* every condition, and pins the under-determined member's constant count.
    The nonlinear-exact family generates `(y²)'' == g` and `(y²)''' + (y²)'' == g` over six
    forcings — exact by construction, so nothing is hand-picked — and asserts both gates.
  - *Regression:* same-machine A/B against a HEAD binary built in an isolated git worktree, over
    ten exposed sections: **better on 4 (2223 +1, 2227 +2, 2232 +2, 2234 +1), equal on 6, zero
    cases lost, 0 FAIL in all 22 runs** — and then the §2.1.2 master corpus, 1204 records through
    both binaries, which is what a latency change has to answer to: **635 PASS against HEAD's 614
    (+21), 0 FAIL on both, ZERO cases lost, timeouts 8 → 4, crashes 2 → 1.** §2.1.2-225
    (`A y + (a+2bx+cx²+y²)² y'' == 0`) is the shape of the gain: both binaries decline it, but
    HEAD's `TimeConstrained[…, 8]` overran to **48 s** where M62 returns at 7.48 s. The one
    remaining crash is the documented intermittent macOS libmalloc/GMP-lock SIGILL the
    `siglongjmp` abort carries (HEAD crashed too, on a different case; neither reproduces in
    isolation over three trials) — and the clamp makes SIGPROF-driven aborts land where the
    child `alarm(20)` used to, so that pre-existing hazard is now exercised more often, which is
    worth a look on its own. Valgrind on the five changed paths is **better** than baseline, not
    merely flat: 13.6 KB definitely lost against HEAD's 74.8 KB, and 6.7 KB indirect against
    2.66 MB — an abandoned unbounded search allocates. **THREE regressions of this
    wave's own were found by measurement and none by reasoning:** the hygiene gate's nine §2.2.33
    losses, the nonlinear peel's two §2.2.34 IVP losses, and — only in the master run — the peel's
    §2.1.2-1143 slowdown, which is why its cascade slot moved twice before it was right.
    `make check-c99` and `check-messages` green; `dsolve_m62_stress_tests` and every other DSolve
    stress suite pass. `dsolve_stress_tests` and `dsolve_m34_stress_tests` fail one assertion each
    — **verified pre-existing**, both failing identically on the HEAD worktree binary
    (`DSolve\`UndeterminedCoefficients[y''-2y'+y == Cos[2x]]` declines on both, and
    `(x^3+2)y''+4xy'+y == 0` times out on both even at 40 s); `dsolve_tests` remains the documented
    `alarm(120)` casualty, so the `t_m62_*` units in it were verified by direct evaluation.

- **M63a — the parameter-as-indvar converter repair: eighteen corpus records were the WRONG
  equation.** ✅ DONE. Converting §2.2.36 turned up one mis-transcribed record (3570,
  `y'' - 2a y' + a²y == 0` read as an ODE *in* `a`), and auditing the class turned up seventeen
  more already sitting in the gated corpora. This is the sixth appearance of the "parses but is
  the WRONG equation" family (M42, M44 ×3, M48 ×3, M49, M61) and the worst-behaved: a
  mis-transcribed record is back-substituted into the *garbled* equation, so it scores PASS or
  UNEVAL for a question the book never asked, and neither the residual, the verify gate, nor the
  ctest baseline can see it. One of the eighteen (`§2.2.16-1534`) had been standing in
  `README.md` as an unexplained residue since M36.
  - **Two gates in `detect_symbols` (`tools/latex_ode_to_mathilda.py`).** *(a)* The lone-letter
    step adopted ANY present non-main letter; its Latin candidate set is now `{y}` alone. The
    argument is complete rather than heuristic: a genuine independent variable with any other
    spelling is in `INDVAR_PREF` and was already claimed by the preferred-letter step, so the one
    thing this step can legitimately contribute is the swapped-variable reading `x = x(y)` — and
    the complete set of legitimate adoptions across all 36 corpora is nine records, every Latin
    one of them `y` (§2.2.1-98/99/100, §2.2.30-2961/2962/2966; the other three are the M49 Greek
    `dr/dθ` path, untouched). *(b)* The record's CAS classification is now plumbed through
    `convert_row`, and a `_missing_x` tag — upstream *stating* that the equation carries no
    explicit dependence on its independent variable — restricts the indvar to a function-argument
    position (`y''(t)`), falling through to a fresh letter otherwise. That is what catches
    `§2.2.2-170` (`r y'' == (1+y'²)^(3/2)`, where `r` is the radius of curvature but is also in
    `INDVAR_PREF`, so the *preferred-letter* step took it — the one victim gate (a) cannot see).
  - **Measured, not argued, in both directions.** Over the 495 `_missing_x` records in the
    committed corpora exactly 8 had the independent variable free in the equation body and all 8
    were mis-transcriptions — zero false positives. The old-vs-new converter run on nineteen
    re-fetched pages moves **13 records, every one a known victim**, and leaves the other ten
    sections byte-identical (§2.2.1/8/17/20/21/24/25/29/30/35 as negative controls); of the 168
    `_missing_x` records on those pages, 159 are byte-identical and the 9 that move each move to
    the right letter (`x` for a `y`-ODE, `t` where `x` is the dependent function). The
    argument-position branch has never fired upstream — not one of the 168 writes `y''(t)` — but
    it keeps the rule honest rather than lucky. A blanket `_quadrature` veto was *rejected* by
    measurement: 14 of those 157 records have the dependent variable undifferentiated with the
    indvar free, 8 are victims and **6 are legitimate** factorable `x`-equations
    (§2.2.8-746, §2.2.17-1682/1684, §2.2.33-3285/3292/3293).
  - **The audit, `make check-corpus-indvar`** (`tools/check_corpus_indvar.py`). The class had been
    found by hand five times and each new section re-opened it, so it is mechanical now: two rules
    read off the record (a `_missing_x` classification forbids the indvar in the body; the indvar
    must be a letter that names a variable by convention), assert-empty, `EXEMPT` for a deliberate
    case. **It earned its place on first run** — it found the eighteenth victim, `§2.2.28-2789`,
    which every hand scan had missed because it is a *system* record: the SIR model
    `x' = -bxy + m, y' = bxy - gy`, whose `m` is the immigration RATE, read as an ODE in `m` while
    every sibling in its block uses `t`.
  - **Four of the eighteen needed no code change at all** (§2.2.12-1157/1182, §2.2.16-1537,
    §2.2.17-1603): their parameter is glued to the dependent letter (`ay`), which the lone-letter
    scan cannot see, so converter work that landed *after* those sections were generated had
    already repaired them — the committed records were stale artifacts. That four wrong equations
    survived many waves *because nobody re-ran the converter* is the argument for the audit, and
    the reason `README.md` now says to run it before trusting a conversion.
  - **Also corrected in `README.md`**, both measured wrong: the fetch URL (`Ch2.S2.SSN.htm` has
    404'd since the 2026-09-28 upstream regeneration; the live path is `Ch2.S1.SSN.htm`) and the
    claim that §2.2.1–19 can no longer be regenerated — every section is back, which is what made
    repairing the pre-LaTeXML records possible at all. Upstream now paginates to §2.1.145.
  - *Corpus effect, measured before/after on the same machine, back to back, per section:*
    **+3 (§2.2.12-1157, §2.2.16-1534, §2.2.33-3247), zero cases lost, 0 FAIL in all 18 runs**;
    fourteen records keep their verdict but now hold the right equation, which is the point — the
    garbled reading of a logistic IVP or a constant-coefficient equation is *also* solvable, and
    that is exactly why the class stayed invisible. §2.2.16-1534 went `$Aborted` (7.9 s) →
    PASS in 0.02 s. Non-PASS over the nine sections 67 → 64. Gate baselines: §2.2.2 7 → 6,
    §2.2.16 3 → 2, §2.2.18 4 → 3, §2.2.26 11 → 9, §2.2.33 22 → 10, §2.2.17/§2.2.28 unchanged —
    and §2.2.12 3 → **4** / §2.2.13 1 → **3**, which are *raised*: both were already red on main
    (the pristine runs measure 5 and 3), so they are now honest instead of aspirational. The
    report refresh surfaced that drift explicitly — eight pre-existing gains and six pre-existing
    losses, every one present in the *pristine* run, so none of it belongs to this change.
    §2.2.12-1133, §2.2.13-1201 and §2.2.13-1219 are the losses worth chasing next.

- **M63b — §2.2.36 corpus (Problems 3501–3600) + a shared integral budget for `Separable`.**
  ✅ DONE. The next hundred, measured end to end: **96/100 → 99/100, 0 FAIL, 0 crash,
  0 timeout**, two per-case-identical runs on each side. New gate
  `dsolve_corpus_2_2_36_tests`; report `DSolve_test_status/reports/2.2.36.md`. 100 scalar
  records, 9 IVPs, from Goode & Annin (44), Goode 2nd ed. (43) and Riley–Hobson–Bence (13).
  Three blocks: **3501–3513** second-order linear with `z` as the independent variable
  (series / special-function — Jacobi, Laguerre, Gegenbauer at a symbolic eigenvalue,
  Liénard, Emden–Fowler), **3514–3556** elementary first order, **3557–3600** second-order
  `_missing_x` + Euler–Cauchy + quadrature, with 3592–3600 repeating 3514–3522 so a fix
  there scores twice. Every `z`-block answer verifies, including the five returned as
  truncated `SeriesData` — the question §2.2.20/21 settled, re-confirmed here.
  - **One root-cause fix, and the diagnosis was not the one the symptom suggested.**
    `3521`/`3599` (`y' == Cos[x−y]/(Sin x Sin y) − 1`) and `3527` (the `Sin[x+y]` sibling,
    as an IVP) are upstream `_separable` with `sympy=True` and each spun to the 8 s wall.
    The obvious reading — Separable needs the angle-addition identity to see the split — is
    wrong **twice**, and bisecting the cascade with the pinned builtins is what showed it.
    (1) `DSolve\`Linearizable` already answers all three from the **original** equation in
    0.06–0.11 s (`ArcCos[C[1] Csc[x]]`, `ArcCos[C[1] Sec[x]]`, `ArcCos[Sec[x]/2]` for the
    IVP), three cascade slots after Separable; the records were never missing a capability.
    (2) Separable's cost is not in its split search — the separability zero test on this
    shape is 2.4 ms — but in `Integrate` on its own **sampled** integrands: sampling at
    integer points turns a mixed-angle kernel into a shifted one, `Csc[2] Csc[x] Cos[x−2]`
    for `g` and the quotient `Cot[2]²/(Csc[2] Csc[Y] Cos[2−Y] − 1)` for `1/h`, each ~20 s
    to come back unevaluated. So one method at cascade slot 7 of 52 was spending the entire
    solve on behalf of every method behind it. An answer-only test would have passed before
    the fix as well as after; this is a latency property and its test asserts a time.
  - **The fix is a budget, and it is SHARED across both Separable entries** (`sep_remaining`
    / `sep_integrate_bounded`, `src/calculus/dsolve_separable.c`): one 6 s wall-clock
    deadline per top-level solve, armed on first use and keyed on `eval_toplevel_id()`, so
    the explicit path and the implicit twin cannot each charge the cascade a full budget.
    Every part of that was forced by a measurement, and two earlier designs were built and
    discarded:
      1. *A `TrigExpand`-before-sampling rewrite.* Faster (0.08 s, because Separable then
         receives `Cot[x] Cot[y]` and both integrals close at once) and **discarded**: it
         makes Separable *claim* these records, and its implicit twin's relation carries the
         sampling artefact `Cot[2]` where the cascade left alone returns the explicit
         `ArcCos` that Mathematica gives. It also does not generalise — `TrigExpand` expands
         multiple angles as well as sums and cannot reach into a denominator, so on
         `Cos[x−2y]/(Sin x Sin 2y) − 1` it half-expands into a form strictly worse to work
         with (Separable's search 12 s → 8 s on it, the full cascade 12.2 s → 19.5 s).
         A node-count "must shrink" guard rescued that, which is when the simpler answer
         became obvious: bound the attempt instead of teaching the method a new identity.
      2. *A per-integral bound.* It cannot serve two measured requirements at once.
         §2.1.2-1134's answer legitimately carries an unevaluated `Integrate` and needs ~5 s
         to get it *decided*: a 1 s bound **loses the record** (an answer at 5.1 s became a
         16 s abort) and a 2 s bound leaves it on the boundary, passing or failing with the
         load. Raising the bound walks the repaired records toward the wall instead, since
         each path pays its own — 3521 closes at 4.1 s under 2 s and 6.2 s under 3 s. The
         shared deadline separates the cases on the quantity that distinguishes them (1134's
         integrals cost seconds, 3521's cost twenty) and charges a hopeless integrand once.
      3. *Six seconds*, where both hold with margin: 1134 unchanged at 5.3 s (twice), 3521
         at 6.15 s against the 8 s bound, and §2.2.36 per-case identical across two runs.
    A timeout is a **decline on both paths**: a decided non-elementary integrand is an
    answer the implicit twin keeps, a timed-out one has decided nothing, and a raw
    unevaluated `Integrate` re-enters the integration cascade on every re-evaluation — the
    trap that turned a 3 s budget into 375 s in M62.
  - *Reach beyond the three records,* because the budget bounds ANY integrand rather than
    one trig shape. Measured A/B on the separable controls: the autonomous non-elementary
    `y' == −2 ArcTan[y]/(1+y²)` **7.76 s → 1.39 s**, `y' == Cos[x]² Cos[2y]²`
    **14.2 s → 8.0 s**, `y' == Cot[x] y/(1+y)` 0.035 s → 0.014 s, and §2.1.2-1134,
    §2.2.2-165, the symbolic-parameter split and the trivial separables all unchanged.
  - *Residue 1, honest:* `3579`, a `y=_G(x,y')` equation with `sympy=False` — the documented
    non-elementary class (§2.2.17 carries eight). Its sibling `3576` solves.
  - **Named by this wave, not done — three unbounded steps the generator exposed.**
    `DSolve\`Exact` and `DSolve\`LieSymmetry` carry the same defect on the same shape
    (~10 s each, measured); no corpus record needs them for this family because
    Linearizable claims it first. `DSolve\`Homogeneous` spins **16 s** on a *scaled*
    mixed-angle RHS — `k(Cos[x−y]/(Sin x Sin y) − 1)` for `k ≥ 2`, identical on both
    binaries — which is why the `k`-varying members of the stress generator never reach
    Separable at all, and is the reason that family is tested through the pinned entry.
    And the deeper cause of the original spin is in `Integrate`, not DSolve:
    `wj_has_kernel_in_denominator` (`src/calculus/integrate_jeffrey.c`) admits Weierstrass
    to the cascade only on a literal `Power[base, negative]` with `x` in the base, but the
    canonicaliser writes `Cos[x−a]/Sin[x]` as `Times[Csc[x], Cos[…]]`, so **no
    `Csc`/`Sec`/`Cot` integrand ever reaches that stage** — `Integrate[Cos[x−a]/Sin[x], x]`
    is elementary (`Cos[a] Log[Sin x] + Sin[a] x`) and the system either spins or falsely
    reports non-elementarity, while `Integrate\`Weierstrass` pinned returns it in 0.007 s.
    Hoisting the admission newly admits multiple-angle integrands whose answers change
    spelling (`Integrate[Csc[x], x]` would move from `½(Log[2−2Cos x] − Log[2+2Cos x])` to
    `Log[Tan[x/2]]`), so it is a milestone with its own measurement.
  - *Also corrected:* `dsolve_separable.c`'s header advertised a `DSolve\`SeparableImplicit`
    builtin that was never registered — both paths are reached through `DSolve\`Separable`,
    and the phantom name is a trap when bisecting with pinned methods.

- **M64 — §2.2.37 corpus (Problems 3601–3700) + `DSolve\`Bernoulli` handles an irrational /
  transcendental constant exponent.** ✅ DONE. The next hundred, measured end to end:
  **96/100 → 98/100, 0 FAIL, 0 crash, 0 timeout**. New corpus `DE_examples_2237.m` (upstream
  §2.1.37, 100 scalar records, 21 IVPs; `make check-corpus-indvar` green); new gate
  `dsolve_corpus_2_2_37_tests` at baseline **2**; report `DSolve_test_status/reports/2.2.37.md`.
  Heavily first-order (20 linear, 10 separable, 8 Bernoulli, and the homogeneous-classA /
  dAlembert / Abel-2nd-type / Riccati families) with a tail of missing-x constant-coefficient
  2nd/3rd order — Mathilda's strong suit, hence the high baseline.
  - **One general fix, diagnosed from inside the cascade with the pinned builtins.** `3666`
    (`y' − y/((π−1)x) == 3x y^π/(1−π)`, `n = π`) and `3668` (`(1−√3)y' + y sec x == y^√3 sec x`,
    `n = √3`) are textbook Bernoulli equations, and `v = y^(1−n)` linearises them for ANY
    constant `n ≠ 1`. Pinned `DSolve\`Bernoulli` **declined** both: its exponent detector read
    `n` off the whole-`Q` logarithmic derivative `Y Q_Y/Q` (`Q = F − Y F_Y`), and `Cancel`
    cannot reduce that ratio when `Y^π` is a non-polynomial kernel, so `n` came back carrying
    `Y` and the free-of guard rejected it. (`3668` additionally spun in a *later* cascade method
    to the wall — UNEVAL by timeout, not by a missing method; owning it in Bernoulli at slot 5
    makes it return fast, the recurring "a later method burns the budget" shape.)
  - **The obstruction is structural, and the fix steps around it in two places.** The evaluator
    collapses `Y^a·Y^b → Y^(a+b)` only for *direct* `Times` factors, and neither `Plus` nor
    `Cancel` will merge `c1 Y^n + c2 Y^n` once the coefficients carry the symbolic `n` — so the
    `Y·(… + B Y^(n−1) + …)` left by `Y F_Y`, and the like-term `Y^n` sums in `A`/`B`, never
    reduce for an irrational `n`. (1) `n` is now read off a **single term** of `Expand[Q]` as
    the per-monomial log-derivative `Y t_Y/t`, where a lone monomial's `Y^(n−1)·Y·Y^(−n)`
    collapses to `Y^0` for every exponent. (2) `A`,`B` are extracted by **abstracting `Y^n` to a
    fresh symbol `W`** (so `F = A Y + B W` is linear) and taking `Coefficient`, replacing the
    `(F − B Y^n)/Y` division that left a spurious `Y^(n−1)` residue for irrational `n`.
    `src/calculus/dsolve_bernoulli.c`. The `recon` reconstruction check and the `bern_mixed_radical`
    / `bern_Y_nonalgebraic` / `bern_Y_in_sum_power` early-decline guards are **unchanged**, so no
    non-Bernoulli form is newly claimed; integer / rational / negative Bernoulli is unaffected
    (verified against the in-tree Bernoulli forms and the full corpus — 0 FAIL, no section regressed).
  - *Residue 2, honest:* `3650` (`y' == (−2x+4y)/(x+y)`, `y(0)=2`) is a Root-object
    homogeneous/Abel IVP — DSolve returns `y = x·Root[cubic]` but the IC cannot fit it (that form
    forces `y(0)=0`), so `C[k]` is left unfitted → UNEVAL by rule; Root-object IC fitting is a
    separate, delicate piece of work. `3662` (`(x−a)(x−b)(y′−√y) == 2(b−a)y`) is solved correctly,
    but its `Sqrt`-branch general solution is not confirmable by the prelude's numeric sampler —
    a verification limitation, not a DSolve gap. v0.266→0.267.

- **M65 — §2.2.38 corpus (Problems 3701–3800, Goode & Annin 4th ed.) + `Integrate`
  linear-argument substitution.** ✅ DONE. The next hundred, measured end to end:
  **98/100 → 99/100, 0 FAIL, 0 crash, 0 timeout**, per-case identical on two runs. New corpus
  `DE_examples_2238.m` (upstream §2.1.38, 100 scalar records, 9 IVPs; `make check-corpus-indvar`
  green); new gate `dsolve_corpus_2_2_38_tests` at baseline **1**; report
  `DSolve_test_status/reports/2.2.38.md`. Second-order-linear dominated (67 `_linear`, 18
  `_with_linear_symmetries`, 7 `_missing_x`, 6 `_missing_y`, 5 `_Emden`, 4 `_exact`, 1 `_Gegenbauer`):
  higher-order constant-coefficient (UndeterminedCoefficients / VariationOfParameters), Euler–Cauchy,
  abstract `F(x)` forcing, and a few Bessel/Legendre/Gegenbauer specials — Mathilda's strong suit,
  hence the high baseline.
  - **One general fix, in `Integrate`, diagnosed from inside the cascade.** `3746`
    (`y'' + 9y == 18 Sec[3x]^3`) was a >25 s timeout → UNEVAL. Its variation-of-parameters
    particular needs `∫ Sin[3x] Sec[3x]^3 dx` and `∫ Cos[3x] Sec[3x]^3 dx` — both elementary
    (`Sec[3x]^2/6`, `Tan[3x]/3`) — but the Jeffrey–Rich Weierstrass stage always substitutes
    `t = Tan[x/2]`, so a **scaled** argument like `Sec[3x]^2` is first multiple-angle expanded into
    a degree-12 rational in `Tan[x/2]`; DSolve then spun trying to simplify it. The cascade order is
    the cause: Weierstrass runs *before* the CRC table and derivative-divides, so it grabs the scaled
    integrand and makes a mess the clean rules never get to answer.
  - **The fix is a new cascade stage, not a change to Weierstrass** (`src/calculus/integrate_linarg.c`,
    stage `IP_LINARG`, run immediately before Weierstrass). A rational trig/hyperbolic integrand
    whose kernels all share one **non-trivial linear argument** `w = a·x + b` (`a` a non-zero number,
    `a ≠ 1` or `b ≠ 0`, with a kernel in a denominator) is reduced via `u = a·x + b` to the
    bare-argument integral `(1/a)·(Integrate[f(u), u] /. u → a·x+b)`, which the recursive cascade
    closes cleanly (reusing the identical bare-argument machinery, so clean cases keep their spelling
    by construction). It **declines** on a bare argument (trivial `w`), on polynomial trig
    (`Sin[2x]^3`, no denominator kernel), and if the bare sub-integral does not close — so Weierstrass
    still gets its turn and nothing else is touched. Leaving Weierstrass's internals alone was the
    low-risk choice: the only behaviour change across the entire integrate / risch / trig / simp / CRC
    suite was one *improvement* (below).
  - *Measured wins, all differentiate-back verified:* `3746` now solves in **0.4 s**
    (`y = C[1] Cos[3x] − Cos[6x] Sec[3x] + C[2] Sin[3x]`); `Int[Sec[3x]^2] = Tan[3x]/3`,
    `Int[Csc[3x]^2] = −Cot[3x]/3`, `Int[Tan[2x]Sec[2x]^2] = Sec[2x]^2/4`, `Int[Sec[x+1]^2] = Tan[1+x]`
    (previously degree-12 `Tan[x/2]` rationals). *Reach beyond trig:* the definite Laplace–Bessel
    transform `∫₀^∞ e^{−c x} J₀(a x) dx` now closes to `ConditionalExpression[1/√(a²+c²),
    Re[c] > 0 ∧ a > 0]` (matching Mathematica) where it had stayed unevaluated — the stage closes a
    trig sub-integral inside the integral-representation path. `test_integrate_intrep` was updated
    from a decline-assertion to that value (the sole test change; the previous decline encoded a
    limitation, not a desired behaviour). No DSolve corpus section regressed.
  - *Residue 1, honest:* `3764` (`y''' + 3y'' + 3y' + y == 2 e^{−x}/(x²+1)`, `sympy=False`) is solved
    **correctly** — `y = (C[1] + C[2] x + C[3] x²) e^{−x} + e^{−x}(x + ArcTan[x](x²−1) − x Log[1+x²])`,
    its integrals already clean — but cold DSolve is ~12.6 s, over the prelude's 8 s wall, so the
    harness scores it UNEVAL by timeout. A latency residue, the same class as §2.2.32's 3161/3164/3165,
    not a missing capability. v0.340→0.341.

- **M66 — §2.2.39 corpus (Problems 3801–3900, Nasser Abbasi) + `Integrate` `Log·trig`
  integration-by-parts stage.** ✅ DONE. The next hundred, measured end to end:
  **95/100 → 96/100, 0 FAIL, 0 crash**. New corpus `DE_examples_2239.m` (upstream §2.1.39, 100
  records — **8 scalar, 0 IVP, 92 systems**; `make check-corpus-indvar` green); new gate
  `dsolve_corpus_2_2_39_tests` at baseline **4**; report `DSolve_test_status/reports/2.2.39.md`.
  **Systems-heavy**: 92 constant-coefficient linear systems (53 2×2, 32 3×3, 7 4×4), including a
  forced block 3822/3870–3876 (`E^{kt}` / trig / `t E^{3t}` forcing) and a few variable-coefficient
  systems; the 8 scalars are second-order linear nonhomogeneous (VoP forcings `Tan`, `Log`, and a
  variable-coefficient `y''+x y`) — the constant-coefficient system + VoP stack solves the bulk out
  of the box, hence the high baseline.
  - **One general fix, in `Integrate`, diagnosed from inside the cascade with the profiler.** `3805`
    (`y'' + 4y == Log[x]`) was a latency UNEVAL: cold DSolve ~9.4 s, over the prelude's 8 s wall. Its
    variation-of-parameters particular needs `∫ Log[x] Sin[2x] dx` and `∫ Log[x] Cos[2x] dx` — both
    elementary (`CosIntegral`/`SinIntegral` forms) — but a `Log·trig` integrand is outside every cheap
    cascade stage (`MATHILDA_INTEGRATE_PROFILE=1` shows DerivativeDivides / RischTranscendental /
    CRCTable all decline) and falls through to the special-function stage `ParallelMixedSpecial`,
    which closed each correctly but at ~4.7 s (the sibling `x³ Sin[x] Log[x]²` costs ~12 s there, per
    an existing comment in `integrate.c`). The two slow searches were the whole ~9.4 s cost of the ODE.
  - **The fix is a new cascade stage, `LogByParts`** (`src/calculus/integrate_logbyparts.c`, stage
    `IP_LOGBYPARTS`), run after the cheap elementary stages and before the `ParallelMixedTower` /
    `ParallelMixedSpecial` tail. It recognises `c Log[g(x)] K(x)` with exactly one `Log` factor, `K`
    free of `Log` and carrying a trig/hyperbolic kernel of `x`, and does one integration by parts:
    `∫ Log[g] K dx = Log[g] V − ∫ V (g'/g) dx`, `V = ∫ K dx`. On this family `V` is an elementary trig
    antiderivative and `V (g'/g) = (trig)/x` closes to `Si`/`Ci`. The cascade's `if (!result)`
    short-circuit means the stage only ever sees integrands the cheap stages declined — so no
    already-fast integral is perturbed (it was the clean low-risk placement) — and its recursive
    sub-integrals run with a new `g_integrate_no_special` counter raised, so a `Log·trig` case whose
    IBP residual is itself non-elementary declines promptly rather than paying `ParallelMixedSpecial`
    twice. Acceptance is an exact `Simplify` diff-back (the sole test, as in `integrate_gammapower.c`;
    `Sinc[z]` from `D[SinIntegral[z]]` is rewritten to `Sin[z]/z` first), so a mis-recognition can
    only decline, never emit a wrong closed form.
  - *Measured wins, all differentiate-back verified:* `∫ Log[x] Sin[x]` **4.9 s → 0.10 s**
    (`CosIntegral[x] − Cos[x] Log[x]`), `∫ Log[x] Sin[2x]` / `∫ Log[x] Cos[2x]` ~4.7 s → ~0.09 s
    (identical forms to the previous `ParallelMixedSpecial` output), `DSolve[y''+4y==Log[x]]`
    **9.4 s → 0.35 s** (`y = C[1] Cos[2x] − C[2] Sin[2x] + ¼(Log[x] − Cos[2x] CosIntegral[2x] −
    Sin[2x] SinIntegral[2x])`). The full integrate/risch suite (35 tests) and every prior DSolve
    corpus section are unchanged (0 regressions).
  - *Residue 4, honest:* `3804` (`y'' + x y == Sin[x]`) is a variable-coefficient Airy-inhomogeneous
    equation whose particular is non-elementary — DSolve declines instantly; `3832`/`3891` are
    variable-coefficient 2×2 / 3×3 systems, outside the constant-coefficient system solver; `3886`
    is a constant-coefficient 3×3 system with an irreducible cubic-`Root` spectrum (`l³−5l²+8l−8`),
    the known spectrum-churn class (>45 s, cf. §2.2.23's 2289/2220). v0.342→0.343.

## Phase 1 — ODE method catalog

Cascade order: cheap deterministic recognizers first. `[✓]` implemented,
`[ ]` planned.

### 1a. First order
- `[✓] Quadrature` — `y^(n)==f(x)`, `f` free of `y`: integrate `n` times + constant polynomial.
- `[✓] LinearFirstOrder` — `y'+p(x)y==q(x)`: integrating factor `Exp[∫p]`.
- `[✓] Separable` — `y'==g(x)h(y)`: `∫dy/h==∫g dx + C[1]`, solved for `y`.
- `[✓] Bernoulli` — `y'==A y + B y^n` (n≠0,1): substitution `v=y^(1-n)` (exponent recovered robustly for any constant n — integer, fractional/negative, AND irrational/transcendental such as Pi or Sqrt[3]; see M64).
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
- `[✓] SolvableForY` — the **"dp" differentiation method** (Maple's `dp`), of which
  `Lagrange` is the linear special case. For `F(x,y,y')==0` **polynomial in `y`**
  (`PolynomialQ` gate), isolate `y=G(x,p)` (`p=y'`), differentiate w.r.t. `x`, and
  recurse the cascade on the induced first-order ODE `dx/dp=G_p/(p−G_x)` for
  `x=X(p,C)` → parametric `{x=X(p,C), y=G(X,p)}` through `dsolve_run_parametric`
  (back-substitution verified). Solves the `y=_G(x,y')` class, e.g.
  `x^(n-1)(y')^n−n x y'+y==0`. Declines the Clairaut/singular case (`p−G_x==0`), a
  `G` free of `p`, a `Root`/implicit branch. Runs late (after every named first-order
  specialist, before the Lie backstop); bounded (M14 kit). See M57.
  `dsolve_solvefor.c`.
- `[✓] SolvableForX` — the `x`-mirror of `SolvableForY` (`x=H(y,p)`, `dx/dy=1/p` →
  `dy/dp=H_p/(1/p−H_y)`). **Pinned-only** (opt-in; automatic yield ~0 and a
  cubic-denominator induced ODE gives a slow-to-verify transcendental branch),
  matching `FirstOrderPowerSeries`/`EigenvalueProblem`. Solves e.g.
  `x−y y'−(y')^2==0`. See M57. `dsolve_solvefor.c`.
- `[✓] Chini` — `y'==f(x) y^n+g(x) y+h(x)` (n≠0,1,2): the reducible-to-autonomous
  sub-class, via `y=f^(-1/(n-1)) u` → `u'==u^n+B u+C` (B,C constant); implicit first
  integral `∫du/(u^n+Bu+C)−x==C[1]` (rational integrand, always elementary) returned
  through `dsolve_run_implicit`. Non-reducible cases decline. `dsolve_chini.c`
  (shared `dsolve_chini_first_integral`).
- `[✓] Abel` — `y'==f3 y^3+f2 y^2+f1 y+f0` (f3,f2≠0): remove the y² term
  (`z=y+f2/(3 f3)`) → Chini n=3 → same implicit first integral. `dsolve_abel.c`
  (thin front-end over the shared Chini helper). The fuller constant-invariant
  class (with an x-rescaling) is future.
- `[✓] AbelAIR` — `y'==N(x,y)/D(x,y)` with `D=c(x)(P1 y+P0)^k` a SINGLE linear y-factor
  (`k∈{1,2}`): the Abel-2nd-kind / rational-in-y forms Riccati/Chini/Abel-1st reject
  (they need `D` free of y). Scaling `u=y/s` (`s=−P0/P1`, the denominator's y-root) makes
  it separable → implicit first integral `∫1/B du|_{u→y/s}−∫A dx`, implicit-function-rule
  verified with the separation checked exactly before integrating. Runs after Chini/Abel-1st,
  before the substitution reductions + Lie backstop. A down payment on the deferred full AIR
  method (M13; the rational-invariant solvable-class table stays future). `dsolve_abel_air.c`.
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
  *Fixed (M52):* the M51 symbolic-coefficient wrong answer (`y''+a²y==Sec[a x] →
  Sec[a x]/a²` mush) is resolved. Two narrow, latency-safe fixes: `dsolve_homog_basis`
  now realifies a **pure-imaginary** symbolic pair by field arithmetic
  `β = Simplify[PowerExpand[√(-r²)]]` (real `Cos[a x]`/`Sin[a x]`, back-substitutable;
  general complex pairs α≠0 untouched, so the 13 §2.1.2 cold-budget cases are unchanged),
  and `UndeterminedCoefficients` now `Simplify`s its residual before the zero-test gate so
  the non-UC `Sec` forcing declines to VoP instead of shipping `Sec[a x]/a²` (a `zero_test`
  false positive on the raw symbolic residual had defeated the gate).
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
- `[✓] OperatorFactor` (`DSolve`DFactor`) — factor a linear operator of order ≥ 2,
  homogeneous **or forced**, two ways (M60). **Right factor:** find a first-order
  `(D − r)`, `r ∈ C(x)` (a hyperexponential solution `Exp[∫r]`, via a rational Riccati
  `Σ a_k P_k(r) == 0` undetermined-coefficient search); peel via operator
  right-division, recurse `DSolve` on the order-(n−1) quotient (which keeps the
  forcing), close with the trailing first-order solve. Order 2 is admitted because
  Kovacic runs first and owns the tidy answers there — the residue it declines is
  reached here, and a closed form beats the Frobenius series that would otherwise win
  (this is also what lets the order-3 peel's own order-2 quotient close).
  **Left factor (Beke):** failing that, run the same search on the adjoint
  `L* = Σ (−1)^k D^k ∘ a_k`; a first-order right factor of `L*` is a first-order LEFT
  factor `(D + s)` of `L`, i.e. an order-(n−1) RIGHT factor `Q`, solved by recursion
  and closed by variation of parameters. `DSolve`DFactor[eqn,y,x]` returns
  `{Dx − r1, Dx − r2, …}` **innermost first**, with any left factor emitted last.
  Runs after Kovacic, before the reduction/series methods. Still future: the
  irregular-singular (double-pole) `r`, symbolic parameters in the coefficients, and
  right factors of order `m` with `2 ≤ m ≤ n−2` (which genuinely need the `m`-th
  exterior power). `dsolve_operator_factor.c` (self-contained; no changes to the
  Kovacic engine).
- `[✓] GeneralizedAiry` — linear of order n ≥ 3 that becomes the pure-power potential
  `u^(n) == A x^m u` after the depression `y = w u`, `w = Exp[−∫c_{n−1}/n]`: the
  fundamental set is `x^j ₀F_{n−1}(; {1 + (j−i)/p : i ≠ j}; A x^p/pⁿ)`, `p = m + n`.
  Symbolic `A` and `m` supported; forcing by variation of parameters; numeric
  self-verify gate. Declines `p == 0`, a non-positive-integer lower parameter (the
  logarithmic Frobenius case), and any non-power potential. See M60.
  `dsolve_genairy.c`.

### 1d. Nonlinear higher-order
- `[✓] ReductionOfOrder` — `y''==F(x,y')` missing y: reduce to first order in
  p=y' (recurse into the scalar engine), then `y=∫p dx + C[2]`. Guards against a
  wrong `Integrate` antiderivative (requires `D[∫p]==p`, decided by `zero_test`
  then `PossibleZeroQ` sampling so a correct-but-unsimplified antiderivative — a
  multi-`Log` `∫Tan`, an `ArcTan[x/Sqrt[C]]` — is accepted while a degenerate
  `y=const` is still rejected) so it declines instead of shipping a degenerate
  solution. Solves `y''==(y')^2`, the autonomous `a+b(y')^2` (→ Tan/Tanh), and the
  Riccati-in-p `c x (y')^2` families.
- `[✓] AutonomousReduction` — `y⁽ⁿ⁾==f(y,y',…,y⁽ⁿ⁻¹⁾)` missing `x`, **any order n ≥ 2**
  (M58): `p=y'(y)`, the derivative chain `D_{k+1}=p·d/dy(D_k)` reduces to an
  order-(n−1) ODE in `p(y)` (recurse), then `y'==p(y)` separable (recurse); stage-1
  constants `C[1..n−1]` frozen to `C[2..n]` (via `dsolve_renumber_constants`) before
  stage 2; final body required to depend on `x` (rejects the degenerate `y=const`).
  For n=2 this is the classical `p p_y==f` (`y y''==(y')² → C E^(C x)`). A stage-2
  spin-guard declines when `∫dy/p` would be non-elementary (a `Log`- or y-denominator
  radicand — Integrate spins uninterruptibly there). Solves the 3rd-order
  `y y'''==y'y''` (`p=Sqrt[C₁y²+C₂]`); the elliptic/hyperelliptic 3rd-order cases
  reduce but decline (non-elementary quadrature). See M58.
- `[~] EnergyIntegral` — `y''==f(y)`: subsumed by AutonomousReduction for the
  elementary cases (`f` free of `y'` is a special case); genuinely elliptic ones
  (`y''==2y^3`, `y''==-Sin[y]`) still decline (`WeierstrassP` inert head: future).
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

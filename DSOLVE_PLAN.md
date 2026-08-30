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

## Architecture

| Piece | File | Role |
|---|---|---|
| Dispatcher + cascade + fail-memo + depth + init | `dsolve.c` | mirrors `integrate.c` |
| Problem substrate (parse / verify / fit / assemble + helpers) | `dsolve_common.{c,h}` | `DSolveProblem`, `dsolve_run` |
| One method each | `dsolve_<method>.c` | `try` + `builtin` + `init` |

Each method file exposes the three-function contract:
`dsolve_<m>_try(P, &nbranch)` (cascade routine, returns branch bodies or NULL),
`builtin_dsolve_<m>(res)` (REPL entry, via `dsolve_method_builtin`), and
`dsolve_<m>_init()` (registers `DSolve`<Name>` with `ATTR_HOLDALL |
ATTR_PROTECTED | ATTR_READPROTECTED` + docstring; chained from `dsolve_init`).

`DSolve` is **HoldAll** (equations and conditions must stay formal); the parser
folds any `D[...]` into `Derivative[...]`, splits equations from point
conditions, detects the output form (`u` → `Function`, `u[x]` → expression), and
parses the options `GeneratedParameters` (default `C`), `Assumptions`, `Method`,
`IncludeSingularSolutions`. Every returned branch is verified by substituting
`u -> Function[{x}, body]` into each residual and requiring the result not to be
a *decidable* non-zero (undecidable is kept, matching Solve); constants are then
fitted to conditions via `Solve`, and `C[k]` renamed to the requested head last.

Reuse: `solvepoly_solve_polynomial_equality` (characteristic roots), `Eigenvalues/
Eigenvectors/JordanDecomposition` (systems), `Series`/`SeriesData` (Frobenius),
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
- **M4 — systems (1e).** ✅ DONE (systems half). `LinearFirstOrderSystem`
  (constant-coefficient, eigen-based, real + complex eigenvalues, constant
  forcing) and `DecoupleSystem`; multi-function verify/fit/assemble added to the
  substrate; the `nfun>1` dispatch route added. Nonlinear higher-order (1d) still
  open.
- **M5** — Kovacic / Frobenius / operator factoring / eigenvalue problems.
- **M6** — Phase 2 PDEs.

## Phase 1 — ODE method catalog

Cascade order: cheap deterministic recognizers first. `[✓]` implemented,
`[ ]` planned.

### 1a. First order
- `[✓] Quadrature` — `y^(n)==f(x)`, `f` free of `y`: integrate `n` times + constant polynomial.
- `[✓] LinearFirstOrder` — `y'+p(x)y==q(x)`: integrating factor `Exp[∫p]`.
- `[✓] Separable` — `y'==g(x)h(y)`: `∫dy/h==∫g dx + C[1]`, solved for `y`.
- `[✓] Bernoulli` — `y'==A y + B y^n` (n≠0,1): substitution `v=y^(1-n)` (exponent recovered robustly, incl. fractional/negative n).
- `[✓] Homogeneous` — `y'==F(y/x)`: substitution `y=v x` → separable.
- `[✓] Exact` — `M+N y'==0`, `M_y==N_x`; + integrating-factor search `μ(x)`, `μ(y)`.
- `[✓] Clairaut` — `y==x y'+f(y')`: general line `y=C[1]x+f(C[1])` + singular envelope (`IncludeSingularSolutions`).
- `[ ] Riccati` — `y'==q0+q1 y+q2 y^2`: reduce to 2nd-order linear (needs the M2/M3 linear engine; may branch).
- `[ ] Lagrange` (d'Alembert) — `y==x f(y')+g(y')`.
- `[ ] Abel` / `[ ] Chini` — implicit solutions via `Solve`.
- `[ ] FirstOrderSubstitution` — `y'==F(ax+by+c)` and related heuristic substitutions.

### 1b. Linear constant-coefficient
- `[✓] LinearConstantCoefficients` — one method for homogeneous + inhomogeneous.
  Characteristic polynomial `Σ a_k λ^k`; roots via `Solve` with derivative-based
  multiplicity + dedup; complex-conjugate pairs → `e^(ax)(Cos,Sin)`; repeated
  roots → `x^k e^(rx)`. Inhomogeneous by variation of parameters (Wronskian /
  Cramer + `Integrate`, particular `Simplify`d).
  *Cosmetic gap:* for simple forcing the var-params particular can carry a
  homogeneous component (`7/2 Cos^2 x` for `7/4`) — correct and verified, less
  tidy than undetermined coefficients (a future refinement).

### 1c. Linear variable-coefficient
- `[✓] EulerCauchy` — `a_n x^n y^(n)+…+a_0 y == g`: indicial polynomial
  `Σ a_k (r)_k` (falling factorial) via trial `x^r`; real roots → `x^r (Log x)^j`,
  complex pairs → `x^a Cos/Sin[b Log x]`, repeated → `Log x` powers; forcing via
  variation of parameters (real-root forcing works; complex-root forcing with a
  hard Wronskian integral declines gracefully).
- `[ ] ReductionOfOrder` — second solution from one known.
- `[ ] ExactODE` — higher-order exact equations.
- `[ ] NormalForm` — remove the `y'` term (2nd order).
- `[~] SpecialFunctionForm` — matches the normalised 2nd-order form
  `y''+P y'+Q y==0`: **Airy** (`P=0`, `Q` linear → AiryAi/AiryBi, verifies) and
  **Bessel / modified Bessel** (`P=1/x`, `Q=±1−ν²/x²` → BesselJ/Y, BesselI/K —
  correct heads; residual reduces only via Bessel recurrences zero_test can't
  decide, so kept as structurally exact). TODO: LegendreP, Hypergeometric/Kummer
  (real heads); inert heads for LegendreQ, HermiteH, Chebyshev, Laguerre,
  Gegenbauer, Jacobi, Whittaker, Mathieu, Spheroidal, Kelvin, ParabolicCylinder,
  Struve, Weierstrass.
- `[ ] Kovacic` — Liouvillian solutions of 2nd-order rational-coefficient ODEs.
- `[ ] FrobeniusSeries` / `[ ] PowerSeries` — series solutions (reuse `Series`).
- `[ ] OperatorFactor` (DFactor) — factor higher-order linear operators.

### 1d. Nonlinear higher-order
- `[✓] ReductionOfOrder` — `y''==F(x,y')` missing y: reduce to first order in
  p=y' (recurse into the scalar engine), then `y=∫p dx + C[2]`. Guards against a
  wrong `Integrate` antiderivative (requires `D[∫p]==p`) so it declines instead of
  shipping a degenerate `y=const`. Solves `y''==(y')^2 → C[2]-Log[C[1]-x]`.
- `[ ] AutonomousReduction` — `y''==f(y,y')` missing `x`: `p=y'(y)`.
- `[ ] EnergyIntegral` — `y''==f(y)`: `(y')^2=2∫f dy`; elliptic / `WeierstrassP` (inert).

### 1e. Systems
- `[✓] LinearFirstOrderSystem` — `Y' == A Y + b`, constant A: eigen-decomposition
  (real λ → `C e^{λx} v`; complex pairs → real `e^{αx}Cos/Sin[βx]` combinations of
  Re/Im of the eigenvector); constant forcing → particular `-A^{-1}b`. Diagonalizable
  case; defective/non-constant-forcing decline. Multi-function verify + IVP fit in the
  substrate (`dsolve_verify_system`/`dsolve_fit_system`/`dsolve_assemble_system`).
- `[✓] DecoupleSystem` — each equation involves one function: recurse into the scalar
  engine per function, renumber the generated constants. Handles variable-coefficient
  components (`y' == x^2 y`).
- `[ ] LinearSystemVarCoeff` (coupled, variable A) — limited/future.

### 1f. Conditions
- `[✓]` IVP (fit constants at one point) — in the substrate.
- `[ ]` Linear BVP (multiple points) — in the substrate + `Solve`.
- `[ ] EigenvalueProblem` — Sturm–Liouville: `Piecewise` with parameter conditions.

## Phase 2 — PDE method catalog

### 2a. First order
- `[ ] PDELinearFirstOrder` (transport / characteristics), `[ ] PDEQuasilinear`
  (Lagrange), `[ ] PDECharpit` (nonlinear complete integral), `[ ] PDEClairaut`.

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
- **Stress tests:** randomized families per method group, back-substitution
  verified (added alongside each method group as it lands).
- **Gates:** `dsolve_tests` (ctest), `valgrind --leak-check=full`, `make
  check-c99`, and REPL spot-checks (`DSolve[...]`, `?DSolve`, `?DSolve`Separable`).

DSolve is a symbolic/structural head (no element-wise numeric mapping); it is a
documented exemption from the packed-aware / Compile surfaces.

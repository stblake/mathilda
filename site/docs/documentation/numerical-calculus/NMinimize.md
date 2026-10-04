# NMinimize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NMinimize[f, x]`**

searches for a global minimum of f with respect to x.

**`NMinimize[f, {x, y, ...}]`**

global minimum with respect to several variables.

**`NMinimize[{f, cons}, vars]`**

global minimum of f subject to the constraints cons (a single constraint, an And of constraints, or additional list elements {f, c1, c2, ...} that are implicitly And-ed).

<details>
<summary>Notes</summary>

Variables may be given as bare symbols, {x, lo, hi} search-interval specs, or indexed variables x\[i\]; a held generator such as Table\[x\[i\], {i, 1, n}\] or Array\[x, n\] for the variables (and Table\[...\] for the constraints) is expanded automatically. Constraints may be equalities (==), inequalities (\<, \<=, \>, \>=), chained inequalities, and their And combinations.  Disjunctive (Or) constraints c1 || c2 are also supported: a point is feasible if it satisfies at least one branch.  Scalar integer variables are declared with Element\[x, Integers\] in the variable list or the constraints.  An empty feasible set returns {Infinity, {x -\> Indeterminate, ...}}. Methods (Method -\> ...): Automatic            uses DifferentialEvolution. "DifferentialEvolution" DE/rand/1/bin with Deb feasibility rules; the default global engine. "NelderMead"          downhill-simplex search with random restarts. "RandomSearch"        multiple random starts refined by a local solver. "SimulatedAnnealing"  Metropolis search with geometric cooling. "SHGO"                Simplicial Homology Global Optimization (Endres, Sandrock & Focke 2018; cf. scipy.optimize.shgo): samples the bounded box, builds a graph, and starts one local search from each "minimizer pool" vertex (a point better than all its graph neighbours), so every basin is reached with few local searches. "DualAnnealing"       Generalized Simulated Annealing (Tsallis & Stariolo 1996; cf. scipy.optimize.dual\_annealing): a heavy-tailed visiting distribution and a generalized Metropolis acceptance rule with reannealing, plus a local search after each Markov chain. "DIRECT"              DIviding RECTangles (Jones, Perttunen & Stuckman 1993; locally-biased DIRECT-L of Gablonsky & Kelley 2001; cf. scipy.optimize.direct): a deterministic Lipschitzian search that normalizes the box to the unit hypercube and repeatedly subdivides the "potentially optimal" cells (those on the lower-right hull of the (size, value) trade-off), needing no derivatives or random seed. "BasinHopping"        Monte-Carlo minimization (Wales & Doye 1997; cf. scipy.optimize.basinhopping): each step randomly perturbs the current point, LOCALLY MINIMIZES it (the "quench"), and accepts the move by a Metropolis rule on the two locally-minimized energies, with an adaptive step size targeting a fixed acceptance rate; strong on funnel-shaped landscapes. The global best is polished with the exact local optimizer.  A method may be given with sub-options as {"Name", "SearchPoints" -\> n, "ScalingFactor" -\> F, "CrossProbability" -\> cr, "RandomSeed" -\> s}.  "NelderMead" also takes the simplex coefficients "ReflectRatio" (default 1), "ExpandRatio" (default 2), "ContractRatio" (default 0.5), "ShrinkRatio" (default 0.5), and "Tolerance" (convergence threshold), and "InitialPoints" -\> {{x1,...}, ...} to seed the initial simplex.  "SimulatedAnnealing" takes "SearchPoints" -\> K (number of annealing chains, default 1), "PerturbationScale" -\> s (trial-step scale, default 1), and "BoltzmannExponent" -\> f (uphill acceptance probability Exp\[f\[i, df, f0\]\]; Automatic keeps the built-in -df/T).  "SHGO" takes "SamplingMethod" -\> "Simplicial" | "Sobol" | "Halton" (default "Simplicial": the exact Kuhn-triangulation graph, best for low dimension and falling back to "Sobol" above 7 variables; "Sobol" and "Halton" are low-discrepancy point sets whose connectivity is a k-nearest-neighbour graph, a documented approximation of scipy's Delaunay pool), "SearchPoints" -\> n (number of sampling points, default 100), "Iterations" -\> k (sampling/refinement rounds, default 1; stops early once no new local minimum appears), and "RandomSeed" -\> s (the QMC digital/fractional shift).  "DualAnnealing" takes "VisitingParameter" -\> qv (visiting-distribution shape in (1, 3\], default 2.62), "AcceptanceParameter" -\> qa (acceptance shape in \[-1e4, -5\], default -5), "InitialTemperature" -\> T0 (default 5230), "RestartTemperatureRatio" -\> r (reanneal once the temperature falls below T0 r, default 2\*^-5), "LocalSearch" -\> True | False (run the per-chain local search, default True), "SearchPoints" -\> K (independent chains, default 1), and "RandomSeed" -\> s; MaxIterations is the per-chain temperature-step budget (default 1000).  "DIRECT" takes "LocallyBiased" -\> True | False (True, the default, is DIRECT-L, biased toward the incumbent and best for few minima; False is the original unbiased DIRECT, better for many minima), "Epsilon" -\> e (potentially-optimal slack, default 1\*^-4), "MaxFunctionEvaluations" -\> m (objective-evaluation budget, default 1000 n), "MaxIterations" -\> k (division-round budget, default 1000), "VolumeTolerance" -\> v and "LengthTolerance" -\> l (stop once the incumbent cell is smaller than these fractions of the box, defaults 1\*^-16 and 1\*^-6), and "MinValue" -\> f\* with "MinValueTolerance" -\> rt (stop once within relative rt of a known optimum f\*).  "BasinHopping" takes "Temperature" -\> T (Metropolis temperature, default 1), "StepSize" -\> s (initial random-displacement half-width, default 0.5), "StepInterval" -\> k (hops between step-size adaptations, default 50), "TargetAcceptanceRate" -\> r (the acceptance rate the adaptation aims for, default 0.5), "StepFactor" -\> a (step-size adjustment factor in (0,1), default 0.9), "SuccessIterations" -\> m (stop a run once the best stalls for m hops; Automatic disables it), "SearchPoints" -\> K (independent multi-start runs, default 1), and "RandomSeed" -\> seed; MaxIterations is the hop count (default 100).  "PostProcess" (any method) controls the final exact local polish: True | Automatic | a named local method ("InteriorPoint", "FindMinimum", "KKT", ...) turn it on; False | None return the raw global-search point.  "PenaltyFunction" (any method) is the function applied to each constraint's violation when scoring infeasible points during the search: Automatic | None keep the built-in squared penalty; a pure function or function symbol (#^2 &, (10 #) &, Sqrt, ...) replaces it (Automatic is #^2 &). Options: Method               global-search selector (see above). WorkingPrecision     MachinePrecision, or a positive digit count (MPFR refinement of unconstrained/box continuous problems). MaxIterations        Automatic or a positive integer cap on generations; default 100. AccuracyGoal         Automatic | Infinity | digits. PrecisionGoal        Automatic | Infinity | digits. EvaluationMonitor    :\> body run on every objective evaluation. StepMonitor          :\> body (accepted). NMinimize is Protected but NOT HoldAll (matching Mathematica); its variables should be unbound symbols, which evaluate to themselves.  During the search their values are set and restored Block-style, so an unbound variable is not left modified.  The search is deterministic for a fixed RandomSeed.  At MachinePrecision the objective and constraints are auto-compiled to bytecode for the trial-point loop, falling back to the interpreter where a construct cannot be compiled.  Expected numeric-domain messages (e.g. Power::infy from a 1/0 in a gradient on a non-differentiable ridge) are quieted during the search.  Returns {fmin, {x -\> xmin, ...}}.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= NMinimize[x^4 - 3 x^2 - x, x]
Out[1]= {-3.51391, {x -> 1.30084}}

In[2]:= NMinimize[{x + y, x^2 + y^2 <= 9}, {x, y}]
Out[2]= {-4.24264, {x -> -2.12132, y -> -2.12132}}

In[3]:= NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, x >= 1}, {x, y}]
Out[3]= {2.33333, {x -> 1.66667, y -> 0.333333}}

In[4]:= NMinimize[{x + y, x + 2 y >= 3, x >= -2}, {Element[x, Integers], Element[y, Integers]}]
Out[4]= {1.0, {x -> -2, y -> 3}}

In[5]:= NMinimize[{x, x > 2 && x < 1}, x]
Out[5]= {Infinity, {x -> Indeterminate}}

In[6]:= NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[6]= {1.41421, {x -> 0.707107, y -> 0.707107}}
```

### Applications (4)

A convex parabola, minimum at the vertex

```mathematica
In[7]:= NMinimize[x^2 - 4 x + 7, x]
Out[7]= {3.0, {x -> 2.0}}
```

An equality constraint: smallest x + y on the unit circle

```mathematica
In[8]:= NMinimize[{x + y, x^2 + y^2 == 1}, {x, y}]
Out[8]= {-1.41421, {x -> -0.707107, y -> -0.707107}}
```

A small linear program with inequality constraints

```mathematica
In[9]:= NMinimize[{2 x + 3 y, x >= 1, y >= 1, x + y >= 4}, {x, y}]
Out[9]= {9.0, {x -> 3.0, y -> 1.0}}
```

Restricting x to the integers

```mathematica
In[10]:= NMinimize[{(x - 3)^2 + 1, Element[x, Integers]}, x]
Out[10]= {1.0, {x -> 3}}
```

## Implementation notes

**Algorithm.** `nm_minimize_driver` (shared with `NMaximize`, registered in
`findmin.c`) parses the objective and variables: bare symbols, `{x, y, …}`
lists, `{x, lo, hi}` search boxes, indexed `x[i]` (rewritten to fresh scalar
symbols so the symbol-keyed solver runs unchanged), and `Element[x, Integers]`
domain declarations. A `{f, c1, c2, …}` first argument is the objective with its
constraints implicitly `And`-ed; `fm_collect_constraints` recurses through
`And`, expands chained inequalities, and routes each leaf to a simple-variable
**box**, a general constraint `g(x) ≤ 0` / `h(x) == 0`, or an `Or` disjunction.
Constraints are enforced two ways at once: a quadratic exterior **penalty**
`Σ max(0, gᵢ)² + Σ hⱼ²` (a disjunction contributes its minimum branch, zero iff
one branch holds), and **Deb's feasibility rules** for *selection* (feasible
beats infeasible; then smaller objective; then smaller violation), with a loose
threshold for search ranking and a tight one for what the caller receives.

The default `Method -> Automatic` is **differential evolution** (`nm_de.c`),
DE/rand/1/bin with Deb-rule selection. Automatic mode strengthens it to
current-to-best/1/bin with a Latin-hypercube initial population (`15n`, clamped),
per-generation dithered `F ∈ [0.5, 1.0]`, `CR = 0.9`, and bounce-back bound
handling; an explicit `"DifferentialEvolution"` keeps the classic `F = 0.6`,
`NP = 10n` path. It stops when the feasible sub-population's objective spread
collapses, then polishes up to `Min[2n, 50]` distinct basins. The other engines
(`nm_neldermead.c`, `nm_randomsearch.c`, `nm_sa.c`, `nm_shgo.c`,
`nm_dual_annealing.c`, `nm_direct.c`, `nm_basin_hopping.c`) implement
restart Nelder–Mead, multi-start local search, geometric-cooling simulated
annealing, SHGO, Tsallis dual annealing, deterministic DIRECT/DIRECT-L, and
basin hopping respectively. The global best is then refined by the **reused
FindMinimum local solver**: plain BFGS (box bounds only) or BFGS inside a PHR
augmented-Lagrangian outer loop when there are general constraints — not SLSQP
or interior-point, which are selectable only for `FindMinimum`. Mixed-integer
problems add integer coordinate descent plus a continuous-relaxation/round/pin
recovery. `NMaximize` builds a synthetic `NMinimize` of `−f` (negating only the
objective element, leaving constraints intact) and negates the reported optimum.

**Data structures.** Population and working buffers are flat row-major `double*`
(`NP×n` in DE, the `(n+1)×n` simplex in Nelder–Mead, chain buffers in SA); the
context `NmDriver` carries the objective program, variables, region bounds,
boxes, general constraints, disjunctions, and integer flags. At
`MachinePrecision` the objective and every general constraint are
**auto-compiled to bytecode** (`compile_expr_ex`, `COMPILE_FOLD_GLOBALS`) and
evaluated by the register machine per trial point, falling back to the
interpreter where a body cannot be lowered. Variables are bound and restored
`Block`-style (`FmVarBind` saves OwnValues + attributes, clears them so the
symbol is a free argument during the solve, and restores on every exit). A
`WorkingPrecision > machine` request refines a continuous, general-constraint-free
result with an MPFR BFGS at the requested bits.

**Complexity / limits.** The RNG is SplitMix64, so the search is **deterministic
for a fixed `RandomSeed`** (default seed fixed); per-attempt seeds derive from
it. Options: `Method` (with per-method sub-options), `MaxIterations` (default
100; setting it disables the Automatic-only best-of-K runs and region growth),
`WorkingPrecision`, `AccuracyGoal`/`PrecisionGoal` (default `WorkingPrecision/2`,
feeding the DE convergence tolerance), `EvaluationMonitor`/`StepMonitor` (fire
inside the local polish, not the global loop). An empty or unmet feasible set
returns `{Infinity, {x -> Indeterminate, …}}`; a success returns
`{fmin, {x -> xmin, …}}` with integer coordinates as exact `Integer`s. Per-trial
evaluation runs under `arith_warnings_mute` and the local solvers under
`g_fm_quiet`, so numeric-domain chatter from probing bad points is suppressed
while driver diagnostics still honour `Quiet[]`/`Check[]`. `NMinimize` is
`Protected` but **not** `HoldAll`.

**Attributes:** `Protected`.

## References

**See also:** [NMaximize](../../numerical-calculus/NMaximize/), [FindMinimum](../../numerical-calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [HoldAll](../../expression-information/HoldAll/), [Rule](../../assignment-and-rules/Rule/), [Sqrt](../../arithmetic/Sqrt/), [Round](../../arithmetic/Round/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- R. Storn and K. Price, *Differential Evolution — a simple and efficient heuristic for global optimization over continuous spaces*, J. Global Optim. **11** (1997) 341–359.
- K. Deb, *An efficient constraint handling method for genetic algorithms*, Comput. Methods Appl. Mech. Engrg. **186** (2000) 311–338 — the feasibility-rule selection.
- M. J. D. Powell, *A fast algorithm for nonlinearly constrained optimization calculations*, Lecture Notes in Math. **630** (1978) — the PHR augmented-Lagrangian polish.
- Source: [`src/numerical_calculus/nm_driver.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/nm_driver.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_direct.c)
- Tests: [`tests/test_dual_annealing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dual_annealing.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)

## Notes & additional examples

### Notes

`NMinimize[f, vars]` searches for a global minimum and returns
`{fmin, {x -> xmin, ...}}`. Variables may be bare symbols, `{x, lo, hi}` search
intervals, or indexed variables, and the constraint set may mix equalities,
inequalities, chained inequalities, `And`, and `Or`. A scalar integer variable
is declared with `Element[x, Integers]`.

The default `Method -> Automatic` is `"DifferentialEvolution"` (DE/rand/1/bin
with Deb's feasibility rules); `"NelderMead"`, `"RandomSearch"`,
`"SimulatedAnnealing"`, `"SHGO"`, `"DualAnnealing"`, `"DIRECT"`, and
`"BasinHopping"` are also available. The global-search point is polished by an
exact local optimizer, so smooth convex and polynomial problems return clean
values. The search is **deterministic for a fixed `RandomSeed`** (the default
seed is fixed), so these results are reproducible. At `MachinePrecision` the
objective and constraints are auto-compiled to bytecode for the trial-point
loop; `WorkingPrecision -> d` refines a continuous box/unconstrained problem to
`d` digits with MPFR. An empty feasible set returns
`{Infinity, {x -> Indeterminate, ...}}`. `NMinimize` is `Protected` but not
`HoldAll`; its variables should be unbound symbols, and their values are set and
restored `Block`-style during the search.

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

## Examples (6)

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
Out[4]= {1.0, {x -> -1, y -> 2}}

In[5]:= NMinimize[{x, x > 2 && x < 1}, x]
Out[5]= {Infinity, {x -> Indeterminate}}

In[6]:= NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[6]= {1.41421, {x -> 0.707107, y -> 0.707107}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NMaximize](../../numerical-calculus/NMaximize/), [FindMinimum](../../numerical-calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [HoldAll](../../expression-information/HoldAll/), [Rule](../../assignment-and-rules/Rule/), [Sqrt](../../arithmetic/Sqrt/), [Round](../../arithmetic/Round/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_direct.c)
- Tests: [`tests/test_dual_annealing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dual_annealing.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)

# Experiment 86 — Basin Hopping across the NMinimize test corpus

This experiment runs the **prior `NMinimize` stress / unit-test problems** (from
`tests/test_nminimize.c`) through the Basin Hopping engine
(`Method -> {"BasinHopping", …}`) and compares each to
`scipy.optimize.basinhopping` on the identical problem, with an **appropriate,
matched parameterisation** on both sides.

Basin Hopping is a **Monte-Carlo minimization** (Wales & Doye 1997): perturb →
local-minimize ("quench") → Metropolis-accept on the minimized energies, with an
adaptive step size. Its reach depends on two knobs that a single run at the scipy
defaults (`stepsize 0.5`, `niter 100`) does not always set well:

- **Multi-start.** A single run crosses widely-separated basins only through the
  displacement, so both implementations benefit from restarts. The matched knob is
  `K = 8` (Mathilda's `"SearchPoints" -> 8`; on the scipy side, 8 seeded
  `basinhopping` runs kept at the min).
- **Step size vs box.** The uniform `0.5` displacement cannot traverse a
  `[-500, 500]` box in 100 hops, so the wide-box deceptive functions get a
  box-scaled `"StepSize"` (Schwefel `150`, Schaffer `40`) — the same on both sides.

Both knobs are applied **identically** to Mathilda and scipy, so the race stays
like-for-like on the algorithm. The harness files (`basin_hopping_testbed.m` /
`.py`) contain only the cleanly-comparable regime; this README documents all five.

Run the harness comparison:

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 86 --system mathilda,python
```

---

## 1. Cleanly comparable — the harness cases (both reach the same global)

Continuous, box-bounded stress functions, 1-D to 10-D, with the matched `K = 8`
parameterisation. The check is the objective at the optimum, rounded per problem;
**0 CHECK-FAIL**. Times are min-of-3 (K = 8 multi-start on both sides).

| # | Problem (corpus test) | Dim | Global | Mathilda | scipy BH | Ratio |
|---|-----------------------|-----|-------:|---------:|---------:|------:|
| 01 | Chained trig (`test_chained_inequality`) | 1 | −1.7602 | **1.9 ms** | 689 ms | **363× faster** |
| 02 | Rugged sine (`test_sa_suboptions`) | 2 | 0 | **4.8 ms** | 1.36 s | **283× faster** |
| 03 | Schwefel-2D `step 150` (`test_sa_deceptive_landscapes`) | 2 | ≈0 | **28.7 ms** | 1.28 s | **45× faster** |
| 04 | Schaffer N.2 `step 40` (`test_schaffer2_simulatedannealing`) | 2 | 0 | **31.3 ms** | 22.68 s | **725× faster** |
| 05 | Ackley-3D `[-32,32]` (`test_indexed_real_coefficient`) | 3 | 0 | **21.2 ms** | 7.67 s | **362× faster** |
| 06 | Rosenbrock-10D (`test_indexed_rosenbrock`) | 10 | 0 | **204.8 ms** | 8.98 s | **44× faster** |

Both implementations reach the same global on every case (**0 CHECK-FAIL**);
Mathilda is **44×–725× faster**. A basin-hopping run spends its time inside the
quenches (K × 100 local minimizations, the same count on both sides), and
Mathilda's objective is auto-compiled to bytecode where scipy calls back into
Python for every evaluation each quench needs.

---

## 2. Mathilda reaches the global where scipy's multi-start stalls

The one deliberate difference between the two engines is the **quench**: scipy uses
L-BFGS-B, Mathilda its own BFGS / augmented-Lagrangian polish. On these **funnel**
landscapes (many local minima organised around a central global) Mathilda's walk
reaches the origin global where scipy's does not — at the *same* `K` on both sides.
Kept out of the harness (the check values diverge, in Mathilda's favour) and
recorded here. `RandomSeed -> 1` / matched seeds.

| Problem (test) | Dim | Global | Mathilda BH | scipy BH (K=8) | scipy BH (single) |
|----------------|-----|-------:|------------:|---------------:|------------------:|
| Rastrigin-5D (`test_de_options_effective`) | 5 | 0 | **0.0** | 0.995 | 1.99 |
| Rastrigin-8D (`test_search_points_honored`) | 8 | 0 | **0.0** | 2.985 | 3.98 |

The gap is even wider at a **single** run (scipy's documented default): on Ackley-3D
scipy's lone run stalls in the flat outer region at **19.76** and on Chained-trig at
**−0.369**, where Mathilda's single run already reaches `0` and `−1.7602`. Those two
recover once scipy is given the same `K = 8` restarts (so they sit in the harness
above); Rastrigin-5D/8D do **not** recover — scipy's L-BFGS-B quench keeps landing on
a lattice ring the `0.5` step will not walk off, even across eight restarts.

---

## 3. Both miss — deceptive high-D without a funnel

Basin hopping walks locally; a function whose global sits far from any descent
funnel defeats it in either implementation.

| Problem (test) | Dim | Global | Mathilda BH (K=8, step 150) | scipy BH (K=8, step 150) |
|----------------|-----|-------:|----------------------------:|-------------------------:|
| Schwefel-10D (`test_de_boundary_no_stagnation`) | 10 | ≈0 | 1126.7 | 236.9 — both far from 0 |

Schwefel's global sits in one corner of a `[-500, 500]¹⁰` box with no funnel toward
it; neither basin-hopping run resolves 10 dimensions of it (scipy edges closer this
time). `"DifferentialEvolution"` / `"DualAnnealing"` — whose heavy-tailed or
population moves are box-scale-independent — are the engines for this shape
(benchmarks `63-global-optimization`, `82-dual-annealing-testbed`).

---

## 4. Constrained problems — `scipy.optimize.basinhopping` is box-only

`scipy.optimize.basinhopping(func, x0, …)` optimizes over the box its inner
minimizer enforces; its Metropolis accept has **no feasibility notion**, so
constrained problems are outside its standard envelope. Mathilda's Basin Hopping
routes each quench through the shared **augmented-Lagrangian polish**
(`nm_local_polish`), so it solves box + linear/nonlinear inequality + **equality** +
**disjunctive (Or)** constrained problems with no engine-specific code. There is no
`basinhopping` baseline for these; Mathilda's results (`"SearchPoints"` as noted):

| Problem (test) | Dim | Constraint kind | Global | Mathilda BH |
|----------------|-----|-----------------|-------:|------------:|
| Disk linear (`test_disk_linear`) | 2 | nonlinear ineq | −4.2426 | **−4.2426** ✓ |
| Quadratic + linear (`test_quadratic_linear`) | 2 | linear ineq | 0.2 | **0.2** ✓ |
| Linear program (`test_linear_program`) | 2 | LP | 3.25 | **3.25** ✓ |
| Equality (`test_equality_constraint`) | 2 | **equality** | 2.3333 | **2.3333** ✓ |
| Equation system (`test_equation_system`) | 3 | **2 equalities** | 0.4286 | **0.4286** ✓ |
| Mishra's Bird disk (`test_mishra_bird_randomsearch`) | 2 | nonlinear ineq | −106.76 | **−106.76** ✓ |
| NMaximize on disk (`test_nmaximize_constrained`) | 2 | nonlinear ineq | 1.4142 | **1.4142** ✓ |
| Disjunctive Himmelblau (`test_disjunctive_constraints`) | 2 | **Or** of disks | 0 | **≈0** at (3,2) ✓ |
| Risk-parity (`test_risk_parity`) | 8 | **equality** `Σw=1` | ≈0 | **≈1e-18** ✓ |

---

## 5. Not applicable — combinatorial / mixed-integer

`scipy.optimize.basinhopping` optimizes over a **continuous** box and has no integer
support, so there is no basinhopping-vs-basinhopping comparison for the combinatorial
testbed:

> mixed-integer LP / quadratic (`test_integer_domain_*`, `test_mixed_integer`), QAP,
> single-machine scheduling, multiple-knapsack, maximum independent set, 4×4
> Latin-square / sudoku, 3-SAT, fixed-charge flow, and the cardinality /
> transaction-cost / adjacency portfolios (`test_qap_assignment`,
> `test_job_scheduling`, `test_multiknapsack`, `test_max_independent_set`,
> `test_sudoku_latin`, `test_3sat_feasibility`, `test_fixed_charge_flow`,
> `test_cardinality_portfolio`, `test_txncost_portfolio`,
> `test_adjacency_assignment`).

Mathilda's Basin Hopping *can* run these (integer coordinates flow through
`nm_local_polish`'s integer descent — e.g. `min (x−3)² + (y+2)²` over integers → `0`
at `(3, −2)`), but the engines built for combinatorial problems are
`"DifferentialEvolution"` / `"RandomSearch"`, whose scipy comparison lives in
`benchmarks/63-global-optimization`.

---

## Summary

Across the `NMinimize` testbed, Basin Hopping's envelope and the comparison:

- **Inside the envelope** (bounded continuous, with the matched `K = 8` multi-start
  and box-scaled step): Mathilda reaches the same global as scipy on every harness
  case and is **44×–725× faster** per solve (compiled objective vs Python callback).
- **On funnel multimodal problems** (Rastrigin-5D/8D) Mathilda reaches the global
  where scipy's basinhopping stalls in a ring minimum *even at the same `K = 8`*, and
  its single run is far more reliable than scipy's single run (Ackley-3D, Chained
  trig, Rastrigin) — the quench matters.
- **Deceptive, funnel-free high-D** (Schwefel-10D): both miss — basin hopping is the
  wrong tool there; DE / DualAnnealing are the engines.
- **Constrained problems**: Mathilda is **strictly more capable** — its
  augmented-Lagrangian quench solves inequality, equality, and disjunctive
  constraints that `scipy.optimize.basinhopping` cannot express at all.
- **Combinatorial / mixed-integer**: no basinhopping baseline exists (scipy is
  continuous-only); DE / RandomSearch are the engines (benchmark 63).

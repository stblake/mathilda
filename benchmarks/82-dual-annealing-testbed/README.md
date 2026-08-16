# Experiment 82 — Dual Annealing across the NMinimize test corpus

This experiment runs the **prior `NMinimize` stress / unit-test problems** (from
`tests/test_nminimize.c`) through the Dual Annealing engine
(`Method -> {"DualAnnealing", ...}`) and compares each to
`scipy.optimize.dual_annealing` on the identical problem — same bounded box,
scipy's default parameters (`visit 2.62`, `accept -5`, `initial_temp 5230`,
`maxiter 1000`), and a matched seed on both sides.

Dual Annealing is a **general-purpose bounded global optimizer** (Generalized
Simulated Annealing + a local search per Markov chain). Unlike SHGO
(`benchmarks/80-shgo-testbed`), its envelope extends to higher dimension, but
`scipy.optimize.dual_annealing` is **box-only** — it has no `constraints`
argument at all — so the constrained, equality-constrained, disjunctive, and
combinatorial parts of the corpus are not dual-annealing-vs-dual_annealing races.
This README documents all five regimes; the harness files
(`dual_annealing_testbed.m` / `.py`) contain only the first, where a like-for-like
timing race is meaningful.

Run the harness comparison:

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 82 --system mathilda,python
```

> **Budget note.** The corpus run surfaced a real bug: `NMinimize`'s generic
> `MaxIterations` default is 100, but scipy's `dual_annealing` default is 1000, and
> the engine was inheriting 100 — a 10×-too-small temperature-step budget that
> starved the anneal on the larger boxes (e.g. Ackley over `[-32, 32]³` stalled at
> 2.12 instead of 0). Fixed so an unset `MaxIterations` uses Dual Annealing's own
> default of 1000 (scipy parity); an explicit `MaxIterations` still wins. The
> results below are post-fix.

---

## 1. Cleanly comparable — the harness cases (both reach the same global)

Continuous, box-bounded multimodal stress functions, spanning 1-D to 10-D. The
check is the objective at the optimum, rounded per problem; **0 CHECK-FAIL**.
Times are min-of-3 (ms).

| # | Problem (corpus test) | Dim | Global | Mathilda | scipy DA | Ratio |
|---|-----------------------|-----|-------:|---------:|---------:|------:|
| 01 | Chained trig (`test_chained_inequality`) | 1 | −1.7602 | **0.61 ms** | 66.3 ms | **109× faster** |
| 02 | Schwefel-2D (`test_sa_deceptive_landscapes`) | 2 | 0 | **1.9 ms** | 124.4 ms | **65× faster** |
| 03 | Schaffer N.2 (`test_schaffer2_simulatedannealing`) | 2 | 0 | **2.3 ms** | 175.2 ms | **76× faster** |
| 04 | Rugged sine (`test_sa_suboptions`) | 2 | 0 | **1.7 ms** | 127.1 ms | **75× faster** |
| 05 | Ackley-3D `[-32,32]` (`test_indexed_real_coefficient`) | 3 | 0 | **4.2 ms** | 288.8 ms | **69× faster** |
| 06 | Rastrigin-5D (`test_de_options_effective`) | 5 | 0 | **7.6 ms** | 394.8 ms | **52× faster** |
| 07 | Rastrigin-8D (`test_search_points_honored`) | 8 | 0 | **16.8 ms** | 635.3 ms | **38× faster** |
| 08 | Rosenbrock-10D (`test_indexed_rosenbrock`) | 10 | 0 | **25.8 ms** | 871.7 ms | **34× faster** |
| 09 | Schwefel-10D (`test_de_boundary_no_stagnation`) | 10 | ≈0 | **29.1 ms** | 790.7 ms | **27× faster** |

Both implementations reach the same global on every case (**0 CHECK-FAIL**);
Mathilda is **27×–109× faster**. Its objective is auto-compiled to bytecode, so
its per-evaluation cost is a fraction of scipy's Python callback over the ~10⁴
evaluations a run needs.

---

## 2. Mathilda reaches a deeper basin than scipy — same algorithm, luckier walk

On these deceptive functions **both** default runs are approximate, but Mathilda's
DualAnnealing reaches a strictly deeper basin than `scipy.optimize.dual_annealing`
at the matched seed. They are kept out of the harness (the check values diverge —
in Mathilda's favour) and recorded here. Deterministic, `RandomSeed -> 1` /
`seed=1`.

| Problem (test) | Dim | Global | Mathilda DA | scipy DA |
|----------------|-----|-------:|------------:|---------:|
| Eggholder (`test_sa_deceptive_landscapes`) | 2 | −959.64 | **−956.9** | −894.6 |
| Drop-wave | 2 | −1 | **−1.000** | −0.936 |
| Griewank-5D (`test_randomsearch_searchpoints_verbatim`) | 5 | 0 | **0.0148** | 0.0320 |
| Modified Ackley-10D (`test_modified_ackley`) | 10 | ≈−0.98 | **−0.975** | −0.957 |

`test_modified_ackley` is the case whose source comment already earmarks it for a
head-to-head against scipy `dual_annealing`; Mathilda reaches the deeper value.

---

## 3. Sharp needles — both implementations are approximate

Functions with a near-measure-zero global basin. Neither default dual annealing
resolves them (both find a good but not exact optimum, or tie); `NMinimize`'s
`"DifferentialEvolution"` is the engine for the ones that are solvable at all.

| Problem (test) | Dim | Global | Mathilda DA | scipy DA |
|----------------|-----|-------:|------------:|---------:|
| Bukin N.6 ridge (`test_bukin6_no_warning`) | 2 | 0 | 0.0162 | 0.0107 — both approximate the ridge |
| Easom spike (`test_initial_points`) | 2 | −1 | **−1.000** | −1.000 — both hit the (π,π) spike |
| Griewank-10D (`test_griewank_*`) | 10 | 0 | 0.0246 | ≈0 — scipy edges it |
| Gaussian well-10D (`test_gaussian_well`) | 10 | 0 | ≈1.0 | ≈1.0 — both miss the narrow well (DE solves it) |

---

## 4. Constrained problems — `scipy.optimize.dual_annealing` is box-only

`scipy.optimize.dual_annealing(func, bounds, …)` takes **no constraints
argument**: it optimizes over a box and nothing else. Mathilda's DualAnnealing
routes its local phase through the shared **augmented-Lagrangian polish**
(`fm_run_penalty`), so it solves box + linear/nonlinear inequality + **equality**
+ **disjunctive (Or)** constrained problems that scipy's dual_annealing cannot even
express. There is no dual_annealing baseline for these; Mathilda's results:

| Problem (test) | Dim | Constraint kind | Global | Mathilda DA |
|----------------|-----|-----------------|-------:|------------:|
| Disk linear (`test_disk_linear`) | 2 | nonlinear ineq | −4.2426 | **−4.2426** ✓ |
| Quadratic + linear (`test_quadratic_linear`) | 2 | linear ineq | 0.2 | **0.1999** ✓ |
| Linear program (`test_linear_program`) | 2 | LP | 3.25 | **3.25** ✓ |
| Equality (`test_equality_constraint`) | 2 | **equality** | 2.3333 | **2.3333** ✓ |
| Equation system (`test_equation_system`) | 3 | **2 equalities** | 0.4286 | **0.4286** ✓ |
| Mishra's Bird disk (`test_mishra_bird_randomsearch`) | 2 | nonlinear ineq | −106.76 | **−106.76** ✓ |
| NMaximize on disk (`test_nmaximize_constrained`) | 2 | nonlinear ineq | 1.4142 | **1.4142** ✓ |
| Disjunctive (`test_disjunctive_constraints`) | 1 | **Or** | 4.0 | **4.0** ✓ |
| Risk-parity (`test_risk_parity`) | 3 | **equality** `Σw=1` | ≈0 | **4.8e-21** ✓ |

**Limit.** On the hardest constrained case — refinery pooling
(`test_refinery_pooling`, **bilinear equalities**, global −3900) — single-chain
DualAnnealing reaches only a feasible −2820; the bilinear-equality feasibility
restoration wants a population, so `"DifferentialEvolution"` is the engine there
(it reaches −3900, benchmarked in `benchmarks/63-global-optimization`).

---

## 5. Not applicable — combinatorial / mixed-integer

`scipy.optimize.dual_annealing` optimizes over a **continuous** box and has no
integer support, so there is no dual-annealing-vs-dual_annealing comparison for the
combinatorial testbed:

> mixed-integer LP / quadratic (`test_integer_domain_*`, `test_mixed_integer`),
> QAP, single-machine scheduling, multiple-knapsack, maximum independent set,
> 4×4 Latin-square / sudoku, 3-SAT, fixed-charge flow, and the cardinality /
> transaction-cost / adjacency portfolios (`test_qap_assignment`,
> `test_job_scheduling`, `test_multiknapsack`, `test_max_independent_set`,
> `test_sudoku_latin`, `test_3sat_feasibility`, `test_fixed_charge_flow`,
> `test_cardinality_portfolio`, `test_txncost_portfolio`,
> `test_adjacency_assignment`).

Mathilda's DualAnnealing *can* run these (integer coordinates flow through
`nm_local_polish`'s integer descent), but the engines built for combinatorial
problems are `"DifferentialEvolution"` / `"RandomSearch"`, whose scipy comparison
lives in `benchmarks/63-global-optimization`.

---

## Summary

Across the `NMinimize` testbed, Dual Annealing's envelope and the comparison:

- **Inside the envelope** (bounded continuous multimodal, 1-D to 10-D): Mathilda's
  DualAnnealing reaches the same global as scipy's on every harness case and is
  faster per solve (compiled objective vs Python callback).
- **On the deceptive functions** where both are approximate, Mathilda reaches a
  **strictly deeper basin** than scipy at the matched seed (Eggholder, drop-wave,
  Griewank-5, modified Ackley-10).
- **Constrained problems**: Mathilda is **strictly more capable** — its
  augmented-Lagrangian polish solves inequality, equality, and disjunctive
  constraints that `scipy.optimize.dual_annealing` cannot express at all. The one
  gap is bilinear-equality feasibility restoration (a population method's job).
- **Combinatorial / mixed-integer**: no dual_annealing baseline exists (scipy is
  continuous-only); DE / RandomSearch are the engines (benchmark 63).

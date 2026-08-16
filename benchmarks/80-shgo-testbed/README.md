# Experiment 80 — SHGO across the NMinimize test corpus

This experiment runs the **prior `NMinimize` stress / unit-test problems** (from
`tests/test_nminimize.c`) through the new SHGO engine
(`Method -> {"SHGO", ...}`) and compares each to `scipy.optimize.shgo` on the
identical problem — same bounded box, same `sampling_method`, same sample count
`n`.

SHGO (Simplicial Homology Global Optimization) is a **low-to-moderate-dimensional
bounded global optimizer**: it samples the box, builds a graph, and starts one
local search from each minimizer-pool vertex. The corpus spans three regimes
relative to that envelope, and this README documents all three; the harness
files (`shgo_testbed.m` / `.py`) contain only the first, where a SHGO-vs-shgo
timing race is meaningful.

Run the harness comparison:

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 80 --system mathilda,python
```

---

## 1. Cleanly comparable — the harness cases (both reach the same global)

Multimodal stress functions and inequality-constrained problems. The check is
the objective at the optimum, rounded per problem; **0 CHECK-FAIL, 9 AHEAD,
2 SLOWER**. Times are min-of-3 (ms).

| # | Problem | Dim | Sampling | Global | Mathilda | scipy shgo | Ratio |
|---|---------|-----|----------|--------|---------:|-----------:|------:|
| 01 | Eggholder | 2 | simplicial | −959.64 | **2.7** | 117.3 | **43× faster** |
| 02 | Schwefel | 2 | simplicial | 0 | **1.6** | 64.0 | **40× faster** |
| 03 | Schaffer N.2 | 2 | simplicial | 0 | **3.0** | 376.7 | **126× faster** |
| 04 | Mishra's Bird (disk) | 2 | simplicial | −106.76 | **4.3** | 27.4 | **6.4× faster** |
| 05 | Rugged sine | 2 | simplicial | 0 | **0.45** | 43.0 | **96× faster** |
| 06 | Ackley | 3 | simplicial | 0 | 22.1 | **11.3** | 2.0× slower |
| 07 | Rastrigin | 5 | simplicial | 0 | 291.2 | **26.1** | 11.1× slower |
| 08 | Griewank | 5 | simplicial | 0 | **5.6** | 31.7 | **5.7× faster** |
| 09 | Disk (`x+y`, `x²+y²≤9`) | 2 | simplicial | −4.243 | **1.0** | 2.2 | **2.2× faster** |
| 10 | Chained trig | 1 | simplicial | −1.760 | **0.067** | 1.9 | **28× faster** |
| 11 | Quadratic half-plane | 2 | simplicial | 0.5 | **0.68** | 0.92 | **1.35× faster** |

Both implementations reach the same global on every case. Mathilda wins 9 of 11,
several by 40–126× — the compiled-objective trial-point loop dominates on the
smooth 2-D functions.

**The two slower cases** are Ackley-3D and Rastrigin-5D, both `"Simplicial"`.
Mathilda's simplicial refinement uses longest-edge bisection with a local-star
reconstruction whose cost grows with the edge count; scipy's `Complex` class is
more optimized in 3–5 dimensions. Both find the global — the gap is construction
overhead, not search quality, and is the clear optimisation target for the
simplicial sampler. (`"Sobol"`/`"Halton"` avoid it but are weaker on these
rugged separable functions.)

---

## 2. Constrained problems `scipy.optimize.shgo` mishandles — Mathilda solves

These are simple / convex constrained test problems. Mathilda's SHGO routes the
local phase through the shared **augmented-Lagrangian polish** (`fm_run_penalty`),
so it reaches the constrained optimum; `scipy.optimize.shgo` returns a
non-optimal vertex, reports no feasible point, or does not finish. They are kept
out of the harness (a divergent or non-terminating baseline yields no valid
timing) and recorded here.

| Problem (test) | Dim | Global | Mathilda SHGO | scipy shgo |
|----------------|-----|-------:|--------------:|-----------|
| Linear program (`test_linear_program`) | 2 | 3.25 | **3.25** ✓ | 20.0 ✗ (non-optimal feasible vertex; a linear objective has no interior minimizer pool) |
| Quadratic + equality (`test_equality_constraint`) | 2 | 2.333 | **2.333** ✓ | no feasible point ✗ (missed the thin `x+y==2` sliver) |
| Minimax Chebyshev (`test_minimax_chebyshev`) | 15 | 0.12512 | **0.12512** ✓ | did not return within the time budget |
| Almgren–Chriss liquidation (`test_optimal_liquidation`) | 16 | 0.03521 | **0.03521** ✓ | did not return within the time budget |
| Risk-parity portfolio (`test_risk_parity`) | 8 | ≈0 | **3.1e-16** ✓ | — |
| VWAP tracking (`test_vwap_tracking`) | 20 | 2.25e-4 | **2.25e-4** ✓ | did not return within the time budget |

**Finding:** on equality-constrained and degenerate (linear-objective) problems,
Mathilda's SHGO is markedly more robust than scipy's — the augmented-Lagrangian
inner solver restores feasibility that sampling alone does not reach.

---

## 3. Outside SHGO's envelope — both implementations miss the global

High-dimensional multimodal functions and sharp-needle functions. SHGO is not the
right tool here in *either* implementation; `NMinimize`'s
`"DifferentialEvolution"` / `"SimulatedAnnealing"` are (they solve these in
`tests/test_nminimize.c`, and are benchmarked in `benchmarks/63-global-optimization`).

| Problem (test) | Dim | Global | Mathilda SHGO | scipy shgo |
|----------------|-----|-------:|--------------:|-----------|
| Bukin N.6 (`test_bukin6_no_warning`) | 2 | 0 | 0.026 | 0.0034 — both approximate (the ridge defeats sampling) |
| Easom (`test_initial_points`) | 2 | −1 | ≈0 | ≈0 — both miss the narrow spike at (π,π) |
| Griewank-10 (`test_griewank_*`) | 10 | 0 | 30.3 | misses; also >120 s (times out) |
| Schwefel-10 (`test_de_boundary_no_stagnation`) | 10 | ≈0 | 1086 | misses |
| Rastrigin-8 (`test_search_points_honored`) | 8 | −80 | −66 | misses |
| Katsuura-8, Lennard-Jones-18, rotated Rastrigin-20 | 8–20 | — | miss (high-dim) | miss |

The global-basin volume shrinks exponentially with dimension, so a bounded sample
budget cannot resolve the optimum past ~6–8 variables — the documented reason
`"Simplicial"` caps at 7 variables and falls back to `"Sobol"`. On the needles
(Bukin ridge, Easom spike) both implementations are approximate.

---

## 4. Not applicable — combinatorial / integer problems

`scipy.optimize.shgo` optimizes over a **continuous** box and has no integer
support, so there is no SHGO-vs-shgo comparison for the combinatorial testbed:

> QAP, single-machine scheduling, multiple-knapsack, fixed-charge flow, maximum
> independent set, 4×4 Latin square / sudoku, 3-SAT, cardinality / transaction-cost
> portfolios (`test_qap_assignment`, `test_job_scheduling`, `test_multiknapsack`,
> `test_fixed_charge_flow`, `test_max_independent_set`, `test_sudoku_latin`,
> `test_3sat_feasibility`, `test_cardinality_portfolio`, `test_txncost_portfolio`,
> …), plus the mixed-integer LP / quadratic cases (`test_integer_domain_*`,
> `test_mixed_integer`).

Mathilda's SHGO *can* run these (integer coordinates flow through
`nm_local_polish`'s integer descent), but the appropriate engines for
combinatorial problems are `"DifferentialEvolution"` / `"RandomSearch"`, whose
scipy comparison lives in `benchmarks/63-global-optimization`.

---

## Summary

Across the `NMinimize` testbed, SHGO's applicability envelope is sharp and the
comparison confirms it:

- **Inside the envelope** (bounded, ≤~6-D, multimodal or convex-constrained):
  Mathilda's SHGO reaches the same global as scipy's and is typically **much
  faster** (9/11 harness cases ahead, up to 126×), the two exceptions being the
  simplicial-refinement cost in 3–5 D.
- **Constrained problems**: Mathilda is **more robust** — it solves linear,
  equality-constrained, and 15–20-D convex testbed problems that scipy's shgo
  returns wrong, infeasible, or does not finish.
- **Outside the envelope** (high-dim multimodal, sharp needles, combinatorial):
  SHGO is the wrong tool in either implementation; `NMinimize`'s DE / SA /
  RandomSearch are the engines (benchmark 63).

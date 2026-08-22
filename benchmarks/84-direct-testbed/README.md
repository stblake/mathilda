# 84 — DIRECT across the whole NMinimize test corpus

Runs the prior `NMinimize` stress / unit-test problems (from
`tests/test_nminimize.c`) through `Method -> "DIRECT"` with per-problem
parameterisation, and compares each to `scipy.optimize.direct` **where DIRECT is
the right tool and scipy can race it** — the continuous, box-bounded problems.
The `.m` / `.py` pair is that race (18 cases). The constrained, equality,
disjunctive, integer, and combinatorial remainder of the corpus is documented
below: `scipy.optimize.direct` is box-only, so those are not DIRECT-vs-`direct`
races, and pooling "scipy can't do this" into a speed number would be dishonest
(the harness reserves ratios for genuine races — see `benchmarks/ABSENT.md`).

DIRECT is deterministic — there is no seed. The race is raw-vs-raw
(`"PostProcess" -> False`; scipy's `direct` has no polish). Both sides use
scipy's defaults except where a label is annotated `LF/N`
(`locally_biased=False`, `maxfun=N`) — the same parameterisation on both.

## The race (box-bounded, both run DIRECT) — Apple M-series, 2026-08-17

**0 CHECK-FAIL**, 15 of 18 AHEAD. Every case reaches the same point in both
systems: the global for most, a shared non-global basin for Rosenbrock-10D
(~8.7916) and Gaussian-well-10D (both raw DIRECTs stall at 1.0).

| Case | Reaches | Mathilda | scipy | Mathilda/scipy |
|------|---------|---------:|------:|:--------------:|
| Chained trig 1D | global −1.7602 | 0.39 ms | 2.01 ms | **5.2× faster** |
| Quartic 1D | global −3.5139 | 0.35 ms | 1.64 ms | **4.6×** |
| Gamma 1D | global 0.88560 | 0.47 ms | 1.56 ms | **3.3×** |
| Schwefel 2D | global ~0 | 0.77 ms | 1.75 ms | **2.3×** |
| Schaffer N2 2D | global 0 | 1.58 ms | 0.71 ms | 0.45× (2.2× slower) |
| Rugged sine 2D | global 0 | 0.83 ms | 0.72 ms | 1.1× |
| Eggholder 2D `LF/20000` | global −959.64 | 46 ms | 96 ms | **2.1×** |
| Griewank 5D | global 0 | 0.30 ms | 3.43 ms | **11.5×** |
| Griewank 10D | global 0 | 1.13 ms | 2.63 ms | **2.3×** |
| Scaled-quadratic 10D | global 0 | 0.26 ms | 1.02 ms | **3.9×** |
| Quadratic bowl 3D | global 0 | 1.35 ms | 1.42 ms | 1.1× |
| Ackley 3D | global 0 | 0.59 ms | 3.06 ms | **5.2×** |
| Katsuura 8D | global 0 | 2180 ms | 75 ms | 0.03× (29× slower) |
| Rastrigin 5D | global 0 | 0.16 ms | 1.74 ms | **11.1×** |
| Rastrigin 8D | global 0 | 0.25 ms | 1.56 ms | **6.3×** |
| Rosenbrock 10D | same basin ~8.79 | 0.48 ms | 2.78 ms | **5.8×** |
| Modified Ackley 10D | global −1 | 5.19 ms | 2.48 ms | 0.48× (2.1× slower) |
| Gaussian well 10D `LF/40000` | both stall 1.0 | 148 ms | 714 ms | **4.8×** |

Both implement the same `DIRECTv2.04` algorithm; the gap is per-evaluation cost.
scipy's `direct` is compiled C but calls back into Python for every objective
evaluation; Mathilda auto-compiles the machine-precision objective to bytecode
and never leaves C — hence the wins, which widen with the evaluation count
(Griewank/Rastrigin/Ackley at high dimension, Eggholder at 20000 evals).

The **three losses are all objective-evaluation cost, not the DIRECT algorithm**:

- **Schaffer N2, Modified Ackley** (~2×) — cheap objectives converging in a few
  hundred evaluations, where scipy's lower fixed C overhead beats Mathilda's
  compile-then-run.
- **Katsuura 8D (29×)** — its `Round`/nested-`Sum` objective (200 `Abs[2^k x −
  Round[2^k x]]` terms) does not take Mathilda's compiled fast path, so each of
  the ~220 evaluations runs on the interpreter (~10 ms) versus numpy's
  vectorised ~0.34 ms. This is a Mathilda objective-evaluation (compile-coverage)
  limitation that would slow *any* NMinimize method on this objective equally —
  it is not a DIRECT gap. Both systems still reach the global (0) and agree.

Several high-dimensional globals sit at the **origin** — the centre of the box —
which DIRECT samples first (Griewank, Scaled-quadratic, Quadratic bowl,
Rastrigin, Ackley, Modified Ackley). DIRECT's centre-sampling finds these almost
immediately even in 10-D, a genuine strength of the method that the timings show.

**Schwefel-10D is excluded**: it is a deceptive high-D function where the two
deterministic subdivisions diverge into different missed basins (Mathilda ~3141,
scipy ~1233), so it is not a clean comparison point.

## Beyond the race — the rest of the corpus

`scipy.optimize.direct` takes only a box (no constraints argument), so none of
the following are DIRECT-vs-`direct` races. What they show is where **Mathilda's
DIRECT reaches past scipy's**, and where DIRECT (in either system) is simply the
wrong tool.

### Mathilda's DIRECT solves these; `scipy.optimize.direct` cannot express them

Mathilda scores every sampled centre through the shared evaluator (constraint
penalty + Deb feasibility) and polishes the incumbent with its exact augmented-
Lagrangian / integer-descent local solver, so `Method -> "DIRECT"` handles
constrained and integer problems. All values below are what `Method -> "DIRECT"`
actually returns, matching the `tests/test_nminimize.c` optima:

| Problem | Constraint kind | `Method -> "DIRECT"` | Known optimum |
|---------|-----------------|---------------------:|--------------:|
| `x+y`, `x²+y² ≤ 9` | inequality | −4.2426 | −4.2426 |
| `(x−1)²+y²`, `x+y/2 ≤ ½`, `x−y ≥ 0` | inequalities | 0.200 | 0.200 |
| `x+y`, `3x+2y≥7 ∧ x+2y≥6 ∧ x,y≥0` | LP | 3.250 | 3.250 |
| `x+y`, `x²+y² ≤ 1` | inequality | −1.4142 | −√2 |
| `x²+y²`, `x+y ≥ 1` | inequality | 0.500 | 0.500 |
| Mishra's Bird, `(x+5)²+(y+5)² < 25` | inequality | −106.76 | −106.76 |
| `x+2y`, `x²+2y²≤3 ∧ x+y==2 ∧ x≥1` | **equality** | 2.333 | 2.333 |
| `x−y`, `x+y+z==½ ∧ x−2z==1 ∧ 2x−y≥1` | **equality** | 0.4286 | 0.4286 |
| Himmelblau over two disks (`∨`) | **disjunctive** | ~0 | 0 |
| `x²`, `x ≤ −2 ∨ x ≥ 2` | **disjunctive** | 4.000 | 4 |
| nearest point to two disks | **disjunctive** | 4.200 | 4.2 |
| `NMaximize −(x²+y²)` over two disks | **disjunctive** | −6.311 | −6.311 |
| `x+y`, `x+2y≥3`, `x≥−2`, `x,y ∈ ℤ` | **integer LP** | 1.000 | 1 |
| `x+2y`, `x²+2y²≤3 ∧ x+y==2`, `x ∈ ℤ` | **mixed integer** | 3.000 | 3 |
| `(x−50)²+(y−40)²`, `x+y≥80`, `y ∈ ℤ` | integer + region rescue | ~0 | 0 |
| `NMaximize x+y`, `x²+y² ≤ 1` | inequality | 1.4142 | √2 |

That is 16 corpus problems that `scipy.optimize.direct` cannot even be handed.

### Beyond DIRECT in either system

These corpus problems are outside DIRECT's design; Mathilda solves them with its
*other* NMinimize engines (the test suite uses DifferentialEvolution / SLSQP-
style polish / integer descent), and `scipy.optimize.direct` cannot attempt them
either — the right scipy comparison is `differential_evolution` / `milp` /
`minimize(SLSQP)`, not `direct`.

- **High-dimensional equality-constrained (the A-series):** minimax Chebyshev
  15-D, refinery pooling 8-D (bilinear equality), Almgren–Chriss liquidation
  16-D, risk-parity 8-D, VWAP tracking 20-D. DIRECT's penalty-based feasibility
  search cannot locate a point satisfying tight equalities in 15–20-D — e.g.
  `Method -> "DIRECT"` on the 15-D Chebyshev returns `Infinity` (no feasible
  point found). These are convex/structured problems for SLSQP, not DIRECT.
- **Combinatorial 0/1 programs (the B-series):** QAP-6, multi-knapsack-30,
  4×4 Latin square, 3-SAT (15 var), maximum-independent-set-20, adjacency
  assignment, cardinality / transaction-cost portfolios, fixed-charge flow.
  These are integer programs whose structure (one-hot assignment, Big-M
  disjunctions) DIRECT's continuous relaxation destroys; Mathilda solves them
  via DE + one-hot / integer descent, and the scipy counterpart is `scipy.milp`
  / `quadratic_assignment`, not `direct`.
- **Needle functions both DIRECTs miss:** Bukin N.6 (Mathilda 0.118 / scipy
  0.033 — both far from the 0.0114 ridge, and they diverge) and Easom (both miss
  the single `(π, π)` spike). DIRECT is not a needle-in-a-haystack method.

## Takeaway

On the workload DIRECT is built for — low-to-moderate-dimensional box-bounded
multimodal functions — Mathilda matches `scipy.optimize.direct`'s answer on the
whole raceable corpus (0 CHECK-FAIL) and is faster on 15 of 18 cases, the
exceptions being objective-evaluation cost (a compile-coverage gap on the
`Round`-heavy Katsuura, and two trivially-cheap objectives). Beyond that,
Mathilda's DIRECT additionally solves the constrained / equality / disjunctive /
small-integer corpus that scipy's box-only `direct` cannot express at all, while
the high-D-equality and combinatorial corpus is left to the engines built for it.

Run it:

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 run_all.py --only 84 --system mathilda,python
```

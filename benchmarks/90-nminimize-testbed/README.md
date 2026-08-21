# 90 — NMinimize robustness testbed (full corpus)

Companion to [`89-nminimize-nmaximize`](../89-nminimize-nmaximize), mirroring the
`79/80`, `81/82`, `83/84`, `85/86` experiment+testbed pairing.

89 applies a fair-comparison envelope so every ratio compares like with like.
**This file holds the landscapes that envelope excludes**, and asks a different
question: not *how fast*, but *did it solve the problem at all*.

## The check carries solution quality, not the objective value

Each case emits `Boole[Abs[fbest - fstar] < 0.001]` against the published global
optimum — `1` found it, `0` did not. A `CHECK-FAIL` therefore means precisely
**"these systems disagree about whether they solved it"**.

This is the only way the harness can express robustness. `run_all.py:521-523`
discards the timings of any `CHECK-FAIL` row, so a testbed that checked the
objective value would throw away every interesting row and report nothing.
Encoding quality *into* the check turns a discarded row into a measurement.

## Result

9 cases: **3 AHEAD, 6 CHECK-FAIL, 0 INCOMPLETE**. Every CHECK-FAIL is a real
solve/no-solve disagreement, not a tolerance artefact.

| case | Mathilda | Mathematica | scipy | reading |
|---|:--:|:--:|:--:|---|
| T1 schwefel 5d | ✗ | — | ✓ | Mathilda misses (returns 2075.19 vs 0) |
| T2 griewank 5d | ✗ | — | ✓ | Mathilda misses (returns 0.0246444 vs 0) |
| T3 drop-wave 2d | **✓** | — | ✗ | **scipy** misses (-0.9362 vs -1.0) |
| T4 rastrigin 10d | ✗ | ✓ | ✓ | Mathilda misses; solves it fine at 5-D |
| T5 styblinski-tang 5d | ✓ | — | ✓ | both solve — 0.536 ms vs 96.4 ms (**180×**) |
| T6 bukin n6 | ✗ | — | ✗ | both miss (razor ridge) |
| T7 eggholder 2d | ✗ | — | ✗ | both miss |
| F1 ineq feasibility | ✗ | ✓ | ✓ | **Mathilda returns an infeasible point** |
| F2 eq feasibility | ✗ | ✓ | ✓ | **Mathilda returns an infeasible point** |

## What this says

**Mathilda trades robustness for speed, and the trade is steep in both
directions.** On the envelope (experiment 89) it is 7×–350× faster than the best
competitor. On the hard corpus it solves **2 of 7** where scipy solves **4 of 7**.
Those two facts are the same fact: the default search budget is small, which is
why it returns in 0.5 ms where scipy takes 96 ms, and why it stops early on
landscapes that need persistence. `nm_de.c:157-177` breaks out as soon as the
feasible sub-population's objective spread collapses.

**The constraint-feasibility rows are the most actionable finding here.** `F1`
and `F2` are not "found a worse optimum" — they are *the returned point does not
satisfy the constraint*. `NMinimize[{x^2+y^2, x+y >= 2}, {x,y}]` gives
`x -> 0.999971, y -> 0.999929`, summing to `1.99990`. Both competitors return
feasible points to machine precision. The violation is fixed at ~1e-4 and does
not respond to `AccuracyGoal`, `PrecisionGoal` or `MaxIterations`, which points
at the augmented-Lagrangian polish termination in `findmin_penalty.c` rather
than at the global search.

**T3 is why the envelope in 89 is not special pleading.** Mathilda finds
drop-wave's true optimum and scipy does not.

## Published global optima used by the checks

Values are the standard ones from the Surjanovic & Bingham virtual library of
simulation experiments (`sfu.ca/~ssurjano/optimization.html`), which is the
reference `scipy`'s own test suite follows.

| function | domain used (scipy side) | global optimum |
|---|---|---|
| Schwefel 5-D | [-500, 500]^5 | 0 at x_i = 420.9687 |
| Griewank 5-D | [-600, 600]^5 | 0 at the origin |
| Drop-wave 2-D | [-5.12, 5.12]^2 | -1 at the origin |
| Rastrigin 10-D | [-5.12, 5.12]^10 | 0 at the origin |
| Styblinski-Tang 5-D | [-5, 5]^5 | -195.830828518 at x_i = -2.903534 |
| Bukin N.6 | [-15,-5] x [-3,3] | 0 at (-10, 1) |
| Eggholder 2-D | [-512, 512]^2 | -959.6406627 at (512, 404.2319) |

Tolerance is 1e-3 on the objective — the question is "right basin or not", not
digit count.

## A harness limitation worth naming

`T6` and `T7` are classified **AHEAD with a timing** even though *both systems
failed*. They agree on `0`, so the check passes and the race proceeds — but
racing two failures is meaningless, and 0.272 ms vs 415 ms on Bukin N.6 says
only that Mathilda gives up faster. The harness has no way to express "agreed,
but agreed on failure". Read `T6`/`T7` as *no result*, not as a win.

## Note: the canonical report is stale

`benchmarks/REPORT.md` and the last `benchmarks/history.jsonl` record are dated
**2026-08-06**. `NMinimize`/`NMaximize` landed **2026-08-14**
(`docs/spec/changelog/2026-08-10.md:1110`) and the SHGO/DualAnnealing/DIRECT/
BasinHopping methods on **2026-08-17**, so *no* optimization experiment added
since — 63, 79–86, and these two — has ever appeared in the canonical report.
`benchmarks/ABSENT.md:26` still reads `` `NMinimize` | _declared, not yet a row_ ``.
A full `make bench-gap` (all 90 experiments, ~20 min) is required to surface
them; `--only` runs write `REPORT.partial.md` and deliberately leave the
canonical files untouched (`run_all.py:1377`).

## Reproduce

```bash
python3 benchmarks/run_all.py --only 90 --check-labels
```

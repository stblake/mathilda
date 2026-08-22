# 89 — NMinimize / NMaximize, timed race

Measures `src/numerical_calculus/nm_driver.c`, `nm_de.c`, `findmin_nm_common.c`
against **scipy 1.17.1** (`differential_evolution`, `minimize(SLSQP)`) and
**Mathematica** (`wolframscript`, running this identical `.m`).

Companion: [`90-nminimize-testbed`](../90-nminimize-testbed) — the landscapes
excluded from the envelope below.

## Result

18 cases, **18 AHEAD, 0 SLOWER, 0 CHECK-FAIL, 0 INCOMPLETE**.
Mathilda is between **1.4×** and **350×** faster than the best competitor.

| case | Mathilda | Mathematica | scipy | ratio (M / best other) |
|---|---:|---:|---:|---:|
| D1 rastrigin 2d (DE) | 0.114 ms | 213.8 ms | 24.9 ms | **0.0046** (219×) |
| D2 rastrigin 5d (DE) | 2.833 ms | 787.1 ms | 166.1 ms | **0.0171** (59×) |
| D3 ackley 10d (DE) | 10.04 ms | 1188 ms | 742.3 ms | **0.0135** (74×) |
| D4 rosenbrock 5d (DE) | 1.529 ms | 916.6 ms | 536.2 ms | **0.0029** (350×) |
| D5 levy 5d (DE) | 1.112 ms | 740.2 ms | 220.3 ms | **0.0050** (198×) |
| D6 sphere 10d (DE) | 2.359 ms | 1117 ms | 589.5 ms | **0.0040** (250×) |
| A1 branin 2d (auto) | 0.083 ms | — | 16.8 ms | 0.0049 (203×) |
| A2 six-hump camel (auto) | 0.086 ms | — | 12.4 ms | 0.0070 (144×) |
| A3 beale 2d (auto) | 0.094 ms | 4.333 ms | 37.7 ms | 0.0217 (46×) |
| A4 cross-in-tray (auto) | 1.802 ms | — | 13.0 ms | 0.139 (7×) |
| M1 nmaximize styblinski 5d | 0.533 ms | — | 89.4 ms | 0.0060 (168×) |
| M3 wrapper base nminimize | 0.120 ms | 211.4 ms | 26.7 ms | 0.0045 (223×) |
| M4 wrapper same via nmaximize | ~0.14 ms | — | ~17 ms | — |
| C1 ineq constrained | 0.255 ms | 0.941 ms | **0.189 ms** | 1.349 |
| C2 eq constrained | 0.218 ms | 0.414 ms | **0.186 ms** | 1.172 |
| C3 mixed integer | 0.019 ms | 0.646 ms | 5.540 ms | 0.029 (34×) |
| I1 rastrigin 5d explicit vars | 2.933 ms | 740.6 ms | 163.3 ms | 0.0180 (56×) |
| I2 rastrigin 5d indexed vars | **119.7 ms** | 773.4 ms | 162.2 ms | **0.738** (1.4×) |

## Three findings this experiment exists to record

**1. Naming `Method` silently cuts the DE budget 7.5×.** `nm_de.c:77-86` gives
DE `maxgen = 150n` under `Method -> Automatic` with no `MaxIterations`, but only
**100** when `Method` is named explicitly. So the obvious way to build an
engine-matched race — pin `Method -> "DifferentialEvolution"` on both sides —
handicaps Mathilda badly. Under the implicit budget, Rastrigin 5-D returns
`2.98488` on seed 1 and **5 of 6 seeds fail to find the global minimum**; at 750
generations every seed returns `0.0`. The `$DE` binding pins
`MaxIterations -> 1500` for exactly this reason. Without it this file would
report Mathilda as both faster *and* wrong.

**2. Indexed variables cost 41×, and nearly erase the advantage.** `I1` and `I2`
are the same mathematics in two spellings. Explicit variables: 2.9 ms, 56× ahead
of scipy. `Table[v[i], {i,1,5}]`: 119.7 ms, **1.4×** — from a 56× lead to
statistical noise. scipy and Mathematica see one vector function in both rows, so
their near-identical timings are the control. This corroborates the ~1000×
indexed-variable note in [`85-basin-hopping`](../85-basin-hopping/README.md) and
quantifies it on a clean workload.

**3. Constrained problems are the one place Mathilda is not ahead.** `C1`/`C2`
are the only rows where scipy wins (ratios 1.35 and 1.17 — under the 1.5
`SLOWER_AT` threshold, so still classed AHEAD, but the lead is gone). They also
hide a defect: Mathilda returns `1.9998` where the true optimum is `2`, with
`x -> 0.999971, y -> 0.999929`, so `x + y = 1.99990` and the constraint
`x + y >= 2` is **violated by 1e-4**. The value is invariant under
`AccuracyGoal`, `PrecisionGoal` and `MaxIterations`, so it is not a budget
problem. `CHECK` rounds at `10^3` here so the speed comparison survives; the
feasibility itself is measured as `F1`/`F2` in experiment 90 rather than being
rounded away.

## Fair-comparison envelope

Following the convention set by [`81-dual-annealing`](../81-dual-annealing/README.md)
and [`85-basin-hopping`](../85-basin-hopping/README.md): a case is excluded from
the timed race when the systems converge to **different basins**, because a
1e-6 check race would then compare different answers rather than different
speeds. Excluded, all moved to `90-nminimize-testbed`:

| excluded | Mathilda | competitor | who is wrong |
|---|---|---|---|
| Schwefel 5-D | 2075.19 | 0.0 (scipy) | Mathilda |
| Griewank 5-D | 0.0246444 | 0.0 (scipy) | Mathilda |
| Rastrigin 10-D | not solved | solved (scipy, Mathematica) | Mathilda |
| Drop-wave 2-D | **-1.0 (correct)** | -0.9362 (scipy) | **scipy** |
| Bukin N.6, Eggholder | not solved | not solved | both |

The envelope is not a shield: drop-wave is excluded because *scipy* fails it.

## Method notes and known asymmetries

- **Bounds.** scipy's `differential_evolution` requires a box; Mathilda does not,
  and grows unbounded coordinates by powers of 10 (`nm_driver.c:397-454`). Each
  function gets its standard published domain, which is the most favourable
  honest choice **for scipy** — Mathilda searches its default ±10 box with no
  such hint. Where the standard domain is tighter than ±10 this hands scipy an
  advantage; that is deliberate and untuned.
- **`Method -> Automatic` has no scipy analogue.** The `A*` rows use
  `differential_evolution` as the closest one. They answer "what do I get if I
  just call it", not "which DE is faster" — that is what `D*` is for.
- **Missing Mathematica cells** (`—` in A1/A2/A4/M1) are cases where
  `wolframscript` emitted no `BENCH` line. Not investigated; the Python column
  carries those rows.
- Timing is min-of-3 after one untimed warm-up, `AbsoluteTiming` /
  `perf_counter`, per `benchmarks/harness.m:53-90`.

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
python3 benchmarks/run_all.py --only 89 --check-labels
```

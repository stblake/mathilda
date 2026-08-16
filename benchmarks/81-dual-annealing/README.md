# 81 — Dual Annealing (`NMinimize`, `Method -> {"DualAnnealing", …}`)

Races Mathilda's `nm_dual_annealing` (`src/numerical_calculus/findmin.c`) against
[`scipy.optimize.dual_annealing`][da] on nine standard bounded multimodal
benchmarks. Both implement Generalized Simulated Annealing (Tsallis / Xiang): a
heavy-tailed visiting distribution, a generalized Metropolis acceptance rule with
reannealing, and a local search after each Markov chain.

## Method

- **Same algorithm parameters on both sides** — scipy's defaults (`visit 2.62`,
  `accept -5`, `initial_temp 5230`, `maxiter 1000`), and a fixed seed
  (`RandomSeed -> 1` / `seed=1`) so both runs are reproducible. The RNG streams
  differ (Mathilda SplitMix64, scipy PCG64), so the trajectories differ; the race
  is like-for-like on the algorithm, not on the stream.
- **Check** = the objective at the global optimum, `Round[10^6 · f]` on the
  Mathilda side and `int(floor(1e6·f + 0.5))` on the scipy side. All nine cases
  agree to 1e-6 (0 CHECK-FAIL).
- 2-D objectives are written with bare variables so they hit Mathilda's compiled
  machine-precision path — the fair analogue of scipy's numpy objective. The 5-D
  Ackley is built over fresh global symbols and passed through `Evaluate[]` (the
  same idiom as `79-shgo`). An inline `Sum[…x[i]…]` would fall to the interpreter
  and is ~1000× slower here — a Mathilda indexed-variable limitation unrelated to
  the Dual Annealing engine, not a property worth timing.

## Cases

Himmelblau, Booth, Beale, Rosenbrock, six-hump camel, Ackley (2-D and 5-D),
Styblinski-Tang, Rastrigin — all with a global both systems reach and polish to
agreement.

## Result

**9/9 AHEAD, 0 CHECK-FAIL.** Mathilda reaches the identical global on every case
and is ~100×–600× faster per solve (Himmelblau 0.19 ms vs 112 ms; Ackley-5D
2.3 ms vs 441 ms). The gap is per-evaluation cost: dual annealing spends ~10⁴
objective evaluations per run, and Mathilda's is compiled to bytecode where scipy
calls back into Python each time.

## Fair-comparison envelope

These nine functions are ones where both implementations reach the same global, so
the 1e-6 check is unambiguous. A few sharp-needle functions (the Bukin ridge, a
narrow Gaussian well) are missed by both and left out; on the broader deceptive
functions Mathilda actually reaches a *deeper* basin than scipy at a matched seed
(Eggholder −956.9 vs −894.6, drop-wave −1.0 vs −0.936). The full corpus sweep,
with all five regimes, is `benchmarks/82-dual-annealing-testbed/`.

## Run

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py \
    --only 81 --system mathilda,python
```

(writes `REPORT.partial.md` / `results/<date>-partial.json` only — the canonical
weekly artifacts are left untouched; do not commit the `.partial` files.)

[da]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html

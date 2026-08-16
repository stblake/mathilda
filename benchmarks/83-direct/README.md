# 83 — DIRECT (DIviding RECTangles) vs `scipy.optimize.direct`

A head-to-head race of Mathilda's `NMinimize` DIRECT engine
(`src/numerical_calculus/nm_direct.c`) against `scipy.optimize.direct` on nine
standard bounded multimodal benchmarks — Himmelblau, Booth, Beale, Rosenbrock,
six-hump camel, Ackley (2-D and 5-D), Styblinski-Tang, Rastrigin.

DIRECT is deterministic, so there is no seed on either side. Both columns use
scipy's default parameters (`locally_biased=True`, `eps=1e-4`, `maxiter=1000`,
`maxfun=1000·n`) and run **raw** — scipy's `direct` does no local polish, so the
`.m` column sets `"PostProcess" -> False` for an apples-to-apples comparison of
the *same algorithm*. The check is the reported objective rounded to `10³`: two
raw DIRECT runs agree at the basin, not to machine precision, so the rounding is
coarser than the polished-engine benchmarks (81 dual-annealing).

## Result (Apple M-series, 2026-08-17)

All nine cases reach the identical basin (**0 CHECK-FAIL**) and Mathilda is
**AHEAD on all nine, ~1.5×–20× faster per solve**:

| Case | Mathilda | scipy | speed-up |
|------|---------:|------:|---------:|
| Himmelblau 2D | 0.75 ms | 1.12 ms | 1.5× |
| Booth 2D | 0.81 ms | 1.51 ms | 1.9× |
| Beale 2D | 0.91 ms | 2.51 ms | 2.8× |
| Rosenbrock 2D | 0.72 ms | 1.69 ms | 2.4× |
| Six-hump camel 2D | 0.96 ms | 5.11 ms | 5.3× |
| Ackley 2D | 0.40 ms | 3.43 ms | 8.6× |
| Ackley 5D | 1.12 ms | 2.33 ms | 2.1× |
| Styblinski-Tang 2D | 0.96 ms | 19.5 ms | 20.3× |
| Rastrigin 2D | 0.85 ms | 1.85 ms | 2.2× |

Both implement the same `DIRECTv2.04` algorithm; the gap is per-evaluation
overhead. scipy's `direct` is compiled C, but it still calls back into Python for
every objective evaluation, where Mathilda auto-compiles the machine-precision
objective to bytecode and never leaves C. The gap widens with the evaluation
count — Styblinski-Tang and six-hump camel run to the `maxfun` cap (2001 evals),
so the Python-callback cost dominates there.

Run it:

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 run_all.py --only 83 --system mathilda,python
```

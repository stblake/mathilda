---
source: src/datetime.c
---
**Algorithm.** `builtin_repeated_timing` takes `(expr)` or `(expr, t)` where `t` (default 1.0 s) is the target measurement budget. It repeatedly `evaluate`s the expression, recording each per-run CPU time via `clock()`/`CLOCKS_PER_SEC` into a growable array, looping until at least 4 runs have completed *and* the accumulated wall budget `t` is reached. The first run's result is kept; later results are freed. It then `qsort`s the timings and returns a **trimmed mean**: it discards the lowest and highest quartile (`count/4` from each end) and averages the middle, returning `{averageSeconds, firstResult}` as a `List` of `EXPR_REAL` and the evaluated expression. Measures CPU time per run.

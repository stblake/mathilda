---
source: src/ml/pca.c
---
**Algorithm.** `builtin_standardize` reads its argument into a row-major `n × dim`
buffer with `ml_read_data` — a flat list loads as `n` observations of *one* variable
(`was_vector`), a matrix as `n` rows of `dim` columns — and calls `ml_standardize(x, n,
dim, rescale = true)`. That kernel makes two passes over the columns: `ml_column_mean`
accumulates each column mean, then `ml_column_sd` accumulates the **sample** standard
deviation with the `n − 1` divisor (the same divisor `StandardDeviation` uses, so the
result agrees with `(x − Mean[x])/StandardDeviation[x]` written by hand). Each entry is
then shifted by its column mean and divided by its column standard deviation. A column
of zero variance is the one special case: rather than divide by zero and propagate a
`NaN` through every downstream row reduction, the kernel leaves such a column at exactly
`0.0` (`v = (sd[j] > 0.0) ? v / sd[j] : 0.0`).

**Data structures.** One contiguous `double` buffer from `ml_read_data`
(`na_load_matrix` / `na_load_vector`), plus a `dim`-length mean and SD scratch. The
result is rebuilt as a plain `List` through `ml_list_of_reals` / `ml_list_matrix`, so
the evaluator's own packing gate — not the ML code — decides whether it is held as a
buffer.

**Complexity / limits.** `O(n · dim)` time, two passes, `O(dim)` extra space. The mean
buffer is `calloc`-ed only to settle a phantom-uninitialised-read warning on the
degenerate `dim == 0` path; it carries no runtime cost.

# NDArray performance — Mathilda vs NumPy vs Mathematica

Wall-clock time on a dense **10,000,000-element `float64`** array (matmul on two
**1000×1000** matrices), best of several runs, same machine (Intel i9-9880H, 8
cores, AVX2; macOS). NumPy 2.4.4; Mathematica 14.0 (packed arrays). Reproduce
with `comparisons/ndarray_reductions_bench.sh`.

Mathilda's NDArray fast paths are **multithreaded** (the shared
`nd_parallel_for` / `nd_parallel_reduce` machinery) and route matrix multiply
through the platform **BLAS** (`cblas_dgemm`). Timing must be measured as
wall-clock, not `Timing[]` (which reports summed CPU time across threads).

| Operation | Mathilda | NumPy | Mathematica |
|-----------|---------:|------:|------------:|
| `Total`             |     3.1 ms |   4.3 ms |   **2.5 ms** |
| `Mean`              |     3.1 ms |   4.3 ms |   **2.7 ms** |
| `StandardDeviation` | **5.8 ms** |  31.7 ms |  19.6 ms |
| `Max`               | **2.4 ms** |   4.2 ms |   4.0 ms |
| `Median`            |   110 ms   | 111 ms   |  **71 ms**   |
| `Sort`              |   310 ms   | **148 ms** | 745 ms   |
| `arr + 2.5`         | **8.0 ms** |  11.7 ms |  13.6 ms |
| `3 arr`             | **7.7 ms** |  12.3 ms |  14.1 ms |
| `arr + arr`         | **7.7 ms** |  11.5 ms |  14.6 ms |
| `arr^2`             |    15 ms   | **11.5 ms** | 14.9 ms |
| `Sin[arr]`          | **9.7 ms** |  68 ms   |  13.2 ms |
| `Exp[arr]`          | **10 ms**  |  51 ms   |  13.2 ms |
| `A . B` (1000²)     |   9.2 ms   |   8.6 ms |   **6.6 ms** |

(Bold = fastest of the three.)

## Summary

After multithreading + BLAS, Mathilda is the **fastest of the three on 7 of 13**
operations (StandardDeviation, Max, the three elementwise arithmetic forms, and
— by a wide margin — `Sin`/`Exp`, whose NumPy versions aren't SIMD-vectorized),
and within a small constant factor of the best on `Total`, `Mean`, `arr^2`, and
matrix multiply. Two operations still trail the best rival:

- **`Sort`** — NumPy leads with its AVX-2 vectorized quicksort (Intel
  *x86-simd-sort*). Mathilda now uses an 8-pass LSD **radix sort** on the double
  bit-patterns (2.6× faster than the previous introsort, and well ahead of
  Mathematica); a portable scalar radix still trails a bespoke SIMD sort by ~2×.
- **`Median`** — a single sequential quickselect; matches NumPy, ~1.5× behind
  Mathematica's.

## How

- **Reductions** (`Total`/`Mean`/`Variance`/`StandardDeviation`/`Max`/`Min`):
  `nd_parallel_reduce` splits the flat range across cores, each folding a private
  pairwise partial, then combines — reaching memory bandwidth on ~4 cores.
  Columnwise (matrix) reductions parallelize over output columns.
- **Elementwise arithmetic** (`Plus`/`Times`/`Power`) and **`Clip`**: the scalar
  operands are folded once so the hot loop touches only array buffers
  (`c * arr` / `arr + c` vectorize), and the element range is threaded. Small
  integer powers (`x^2`, `x^3`, `1/x`) skip the libm `pow()`.
- **Elementwise functions** (`Sin`/`Exp`/`Log`/…): already threaded; the
  real-axis-escaping tail (`Sqrt`/`Log`/`ArcSin`) is now threaded too.
- **`Dot`**: rank-2 float64 matrix multiply routes to `cblas_dgemm` (Apple
  Accelerate / system BLAS — multithreaded + vectorized), replacing a naive
  O(m·n·k) scalar loop.

All results are exact matches to Mathilda's equivalent `List` computation, and
every summation uses pairwise accumulation, so accuracy matches NumPy.

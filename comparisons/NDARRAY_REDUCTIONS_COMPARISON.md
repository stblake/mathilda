# NDArray reductions & sort — Mathilda vs NumPy vs Mathematica

Time to reduce or sort a dense **10,000,000-element `float64`** array, best of
3–5 runs, same machine (Apple M-series, macOS). Mathilda uses
`First[Timing[…]]` (CPU time; ≈ wall for these serial memory-bound loops); NumPy
uses `min(perf_counter)` (NumPy 2.4.4); Mathematica uses
`Min[Table[First[AbsoluteTiming[…]], …]]` (Mathematica 14, packed array).
Reproduce with `comparisons/ndarray_reductions_bench.sh`.

| Operation | Mathilda | NumPy | Mathematica |
|-----------|---------:|------:|------------:|
| `Total`             |   5.2 ms |   4.7 ms |   **2.4 ms** |
| `Mean`              |   5.2 ms |   4.3 ms |   **2.6 ms** |
| `StandardDeviation` | **11.9 ms** |  33.2 ms |  20.6 ms |
| `Max`               |   7.1 ms |   4.3 ms |   **2.4 ms** |
| `Median`            | 108 ms   | 114 ms   |  **67 ms**   |
| `Sort`              | 858 ms   | **150 ms** | 734 ms   |

(Bold = fastest of the three.)

## Reading the table

- **O(n) reductions (`Total`, `Mean`, `Max`)** — all three are memory-bound.
  Mathematica's packed-array kernels lead (~2.4 ms ≈ single-thread memory
  bandwidth on this machine); Mathilda and NumPy trail by a small constant
  factor. Mathilda's loops use several independent accumulators so the compiler
  vectorizes an associativity-preserving reduction without `-ffast-math`; the
  residual gap to Mathematica is hand-tuned SIMD kernel quality, not algorithm.
- **`StandardDeviation` / `Variance`** — Mathilda is **fastest** here: a two-pass
  (mean, then Σ|x−μ|²) through the same pairwise/unrolled kernel beats both
  NumPy and Mathematica.
- **`Median`** — quickselect (O(n), no full sort). Competitive with NumPy;
  Mathematica's is fastest.
- **`Sort`** — NumPy leads decisively: NumPy 2.x sorts with a vectorized
  AVX-2/AVX-512 quicksort (Intel *x86-simd-sort*). Mathilda's fast path is a
  scalar median-of-three introsort — already faster than the old `qsort`
  behaviour and comparable to Mathematica's `Sort`, but a portable scalar sort
  cannot match a bespoke SIMD sort. This is a deliberate scope boundary.

## Takeaway

Mathilda's NDArray reductions are in the same performance class as NumPy and
Mathematica — matching or beating them on the variance family and holding a small
constant-factor gap on the cheapest memory-bound reductions — while returning an
`NDArray` (or scalar) exactly as they do. All results are exact matches to
Mathilda's equivalent `List` computation (the fast paths fall back to the `List`
path — the oracle — for any dtype/shape/spec they do not handle), and every
summation uses pairwise accumulation, so accuracy matches NumPy.

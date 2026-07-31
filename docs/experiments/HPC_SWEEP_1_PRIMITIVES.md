# Experiment 9 — First and second HPC sweeps: 43 classical kernels

**Date**: 2026-07-30 → 2026-07-31 ·
**Commits**: `6c78653` and follow-ups ·
**Harness**: [`comparisons/hpc_bench.py`](../../comparisons/hpc_bench.py) ·
**Report**: [`docs/design/performance.md`](../design/performance.md) §2–§8

Common method in [`README.md`](README.md).

---

## Hypothesis

Experiments 1–8 each fixed something a *profile* pointed at. A profile only
shows you the code you are already running, which means it cannot find the thing
you are doing 40× too slowly *everywhere*, because that cost is spread evenly
and looks like the baseline.

The hypothesis was that a broad, adversarial benchmark against a mature
competitor would find a different class of problem than profiling does — and
specifically that it would find **fast paths that exist and are never reached**.

## What was built

38 kernels in the first pass, 5 more in the second, across:

dense linear algebra · spectral transforms · STREAM-style array traffic ·
PDE stencils · compiled scalar loops · integer and combinatorial search ·
arbitrary-precision arithmetic · a Krylov solver · direct convolution ·
ODE time integration · an irregular hash-keyed reduction · interpolation

Plus the harness discipline that made the numbers trustworthy — all of it
retained for the later sweeps:

- **`AbsoluteTiming`, never `Timing[]`.** `Timing[A . A]` on a 600×600 reads
  10.5–13.7 ms against a true 1.9 ms, because Accelerate's `dgemm` is threaded
  and `Timing` sums CPU over threads. Any comparison of a threaded path against
  a serial one using `Timing[]` is wrong by roughly the core count.
- **`AbsoluteTiming` did not exist in Mathilda** and was written for this.
- **Min of N reps, max recorded as a tripwire.** A builtin with an internal
  cache answers the second rep in ~0 s; reporting the min alone makes the row
  vacuous.
- **Values must be forced.** Wolfram evaluates `(p + k) q;` on two 10⁶-bit
  integers in 138 µs and `IntegerLength[(p + k) q]` in 2.96 ms — the first never
  materialises the product. Every row was re-checked in both forms.
- **Checked answers.** Each kernel's scalar result is compared across systems.

## What it found

### The 40× and 80× rows: a complete fast path nobody could reach

`Fourier` was the worst result in the document. `fourier.c` had had a full
NDArray implementation all along — and the packing transparency gate
materialised the buffer into 2²⁰ `Expr` nodes before the builtin could see it.

> `Fourier[packedList]` — **986 ms**.
> `Fourier[NDArray[…]]` — **62.8 ms**. *On identical data.*

Nothing was slow. The fast path was simply unreachable, because the head was not
on the aware list. Two more causes on the same row (an identity `b`-gather doing
an integer division and a modulo per element per axis, and a bounce through a
second buffer) took it to within 1.15×/1.30× of Mathematica.

This is the finding the experiment was for, and it recurred: the same shape hit
`Nest`, then `Outer`, then `ListConvolve`. `make check-packed-aware` exists
because of it.

### `Interpolation` was quadratic in its own table

151 µs per point at 10⁴ nodes, 906 µs at 5×10⁴. At 10⁵ nodes the benchmark could
not be timed at all — 228 s.

Four per-point costs, none depending on the point. Three were table-wide facts
recomputed per call, now memoised on the table's identity (sound because
`expr_copy` is a refcount bump, so the memo holds a reference). The fourth was a
genuine O(n²): `build_grid` inserted coordinates one at a time into a sorted
array scanning from the front, so coordinates **arriving in order** — which is
every table from `Table[{x, f[x]}, {x, a, b, dx}]` — walked the whole grid on
every insertion. Sort-and-unique instead: the 160000-node grid build went from
**4.96 s to 21 ms**.

**1.494 s → 7.65 ms (195×)**, and from 74× slower than Wolfram to 2.6× faster.

### `Tally` was *slower* packed than unpacked

No algorithm to fix — `Tally` already hashes. The cost was boxing: it reads
`list->args[i]`, so a packed buffer was materialised into one `Expr` per element
first, and every probe hashed and compared `Expr` nodes. Over 10⁶ packed int64:
**147 ms packed against 96 ms for the identical data as a plain `List`.**

A representation that is a net loss is the failure mode worth having a benchmark
for. It is invisible to profiling, because the profile of the packed run looks
exactly like the profile of the unpacked one, only longer.

### The `PadRight` row is a semantic fact, not an inefficiency

`PadRight[{1., 2., 3.}, 5]` is `{1., 2., 3., 0, 0}` — exact Integer zeros beside
Reals — in **both** systems, so neither can hold it in a uniform buffer and both
fall back to a boxed list. Given an explicit Real fill, Mathilda is the faster of
the two by 1.29×.

Worth stating because it is the kind of row that invites a "fix" that would
change an answer.

## Results

Full tables in [`performance.md`](../design/performance.md) §2–§8. The shape:

| Benchmark | Mathilda | Mathematica 14.0 | NumPy / Python | vs WL | vs NumPy |
|---|---:|---:|---:|---:|---:|
| Matrix multiply, 1000x1000 | 7.42 ms | 6.46 ms | 7.07 ms | 1/1.15x | 1/1.05x |
| Fourier, 2^20 reals | 22.58 ms | 20.99 ms | 40.75 ms | 1/1.08x | 1.80x |
| Sort | 318.55 ms | 714.17 ms | 155.76 ms | 2.24x | 1/2.05x |
| Jacobi 5-point relaxation, 512^2, 100 sweeps | 223.01 ms | 110.71 ms | 68.78 ms | 1/2.01x | 1/3.24x |
| Mandelbrot, 800x800, 100 iterations | 769.23 ms | 1.621 s | 980.20 ms | 2.11x | 1.27x |
| Sieve of Eratosthenes to 10^7 | 664.56 ms | 784.49 ms | 26.44 ms | 1.18x | 1/25.14x |
| 50000! (exact) | 2.62 ms | 2.29 ms | 45.16 ms | 1/1.14x | 17.25x |
| Conjugate gradient, 2D Poisson 256^2, 100 iterations | 147.90 ms | 82.44 ms | 45.09 ms | 1/1.79x | 1/3.28x |
| ListConvolve, 1024^2 image, 5x5 kernel | 35.67 ms | 33.82 ms | 39.71 ms | 1/1.05x | 1.11x |
| NDSolve, Lorenz system to t = 200 | 26.70 ms | 47.78 ms | 962.09 ms | 1.79x | 36.03x |
| Interpolation, 10^4 nodes, 10^4 evaluations | 6.39 ms | 20.00 ms | 798 us | 3.13x | 1/8.00x |
| Tally, 10^7 integers into 10^4 bins | 122.08 ms | 17.03 ms | 106.48 ms | 1/7.17x | 1/1.15x |

> **Caveat on the Jacobi row.** It reads 223 ms above and 128 ms measured alone.
> Bisecting the benchmark prefix found the trigger: running one 1000×1000 `dgemm`
> first makes Mathilda's threaded stencil 1.45× slower, reproducibly. See
> [`HPC_SWEEP_3_NUMPY_GAP.md`](HPC_SWEEP_3_NUMPY_GAP.md) §"A separate finding".
> Both threaded-stencil rows (Jacobi, Game of Life) follow the linear-algebra
> group in the run order and should be read with that in mind.

## What the second sweep added

Five kernels chosen to run where none of the 38 did. Two came back faster than
Wolfram with no work at all — `NDSolve` at 1.52×, which is the RHS compiler
paying off on the workload it was built for.

The important one was **conjugate gradient**, which was slower than Wolfram
despite every one of its parts being at or ahead of parity. A composition can be
slow in ways none of its components are. That single observation is why the
third sweep is built out of *application pipelines* rather than primitives — see
[`HPC_SWEEP_2_APPLICATIONS.md`](HPC_SWEEP_2_APPLICATIONS.md).

## Still open

The gate allowlist in `hpc_bench.py` is a to-do list, not a settled state.
Remaining from these two sweeps: `ListConvolve` (addressed in experiment 11),
`Tally`'s hash probe, and the `Reverse`/`RotateLeft`/`Differences`/`Accumulate`
band of serial buffer work — partly addressed in experiment 11, partly still
wanting SIMD.

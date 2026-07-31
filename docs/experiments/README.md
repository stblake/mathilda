# Execution-speed experiments, 2026-07-26 → 2026-07-31

Eleven experiments, run over six days, each aimed at one hypothesis about why
Mathilda was slower than it should be. One file per experiment. Every file
reports the same three columns — **Mathilda**, **Mathematica 14.0**, and
**Python** (NumPy/SciPy where the kernel is vectorisable, CPython otherwise) —
because the two comparisons answer different questions:

- **Mathematica** says whether Mathilda is behind a *competitor*: a mature CAS
  with the same evaluation semantics, the same automatic compilation, and the
  same automatic packed arrays. A gap here is a gap in engineering.
- **NumPy** says whether Mathilda is behind *the machine*. On this host NumPy
  links the same Apple Accelerate BLAS that Mathilda does, so on the dense rows
  all three call byte-identical kernels and any spread is pure overhead. On the
  array rows NumPy is a thin memory-bandwidth reference.

A row can read acceptably against one column and badly against the other, and
that difference is usually the finding. The sieve is 1.18× *ahead* of
Mathematica and 25× behind NumPy; before the fourth sweep the return-series
kernel was 13.5× behind Mathematica and 66.5× behind NumPy. Only the second
number separates "slower than a competitor" from "slower than a memory copy".

Where NumPy has no equivalent (arbitrary precision, symbolic rule dispatch,
`PrimePi`) the column reads `—`. Where the honest Python answer is a CPython
loop or `mpmath` rather than a library call, the row says so — **numba is not
installed on this host**, and those rows must not be read as library results.

## Method

Common to every experiment; stated once here rather than eleven times.

- **Machine**: Intel Core i9-9880H (8C/16T, AVX2), 16 GB, macOS 15.7.4.
- **Versions**: Mathilda at `-O3` with GMP 6.3, MPFR 4.2.2, FLINT 3.6,
  FFTW 3.3.11, Apple Accelerate. Mathematica 14.0.0. NumPy 2.4.4 / SciPy 1.17.1
  on CPython 3.11, also linked against Accelerate.
- **Wall clock**: `AbsoluteTiming` in both CAS, `time.perf_counter` in Python.
  Never `Timing[]`, which sums CPU time over threads and so overstates a
  threaded path by roughly the core count.
- **Minimum of several repetitions** after one untimed warm-up, with the maximum
  recorded as a caching tripwire.
- **The same source text** runs in Mathilda and Mathematica; the Python column is
  a separate implementation of the same algorithm in the same order.
- **Values are checked, not just timed.** Each benchmark has a cheap scalar
  answer recorded from all three systems and compared. A timing row is
  meaningless until they agree — the check pins the *algorithm*, so the columns
  cannot silently be timing three different computations.

Reproduce with [`comparisons/hpc_bench.py`](../../comparisons/hpc_bench.py).
Within-Mathilda before/after deltas use the kill switches (`MATHILDA_NO_PACK`,
`MATHILDA_NO_NUMLOOP`, `MATHILDA_NO_AUTOCOMPILE`, `$AutoArrayPacking`,
`$AutoCompilation`), so "before" is the same binary with one subsystem disabled
rather than an older build measured on a different day.

## The experiments

| # | Experiment | Dates | Headline |
|---|---|---|---|
| 1 | [`COMPILE_ENGINE.md`](COMPILE_ENGINE.md) | 07-26 → 07-28 | A typed bytecode VM for numeric code: **~234×** over the interpreter |
| 2 | [`COMPILE_ARRAY_FUSION.md`](COMPILE_ARRAY_FUSION.md) | 07-27 | Machine arrays and fused elementwise loops: **3.2–6.6×** at 10⁶ |
| 3 | [`COMPILE_OPTIMISING_CODEGEN.md`](COMPILE_OPTIMISING_CODEGEN.md) | 07-27 | CSE, folding, DCE, LICM, threaded dispatch: **1.5×** on top |
| 4 | [`AUTO_COMPILATION.md`](AUTO_COMPILATION.md) | 07-26 → 07-29 | The compiler runs without `Compile[]`: **6×–353×** |
| 5 | [`MACHINE_INTEGERS.md`](MACHINE_INTEGERS.md) | 07-29 → 07-30 | int64 as a peer of float64, without losing exactness |
| 6 | [`PACKED_ARRAYS.md`](PACKED_ARRAYS.md) | 07-30 | Dense lists become machine buffers invisibly: up to **56000×** |
| 7 | [`BLAS_LAPACK_ROUTING.md`](BLAS_LAPACK_ROUTING.md) | 07-30 | Reaching the vendor kernels, and a rejected scratch pool |
| 8 | [`RANDOM_NUMBER_GENERATION.md`](RANDOM_NUMBER_GENERATION.md) | 07-30 | A machine-precision generator: **53×–70×** |
| 9 | [`HPC_SWEEP_1_PRIMITIVES.md`](HPC_SWEEP_1_PRIMITIVES.md) | 07-30 | 43 classical kernels against Mathematica; what it exposed |
| 10 | [`HPC_SWEEP_2_APPLICATIONS.md`](HPC_SWEEP_2_APPLICATIONS.md) | 07-31 | Eight real-application pipelines, three systems |
| 11 | [`HPC_SWEEP_3_NUMPY_GAP.md`](HPC_SWEEP_3_NUMPY_GAP.md) | 07-31 | Closing on NumPy: structural ops, scans, convolution |

## The one finding that recurs

Six of the seven fixes in experiment 10, and four of the six in experiment 11,
are the same defect wearing different clothes:

> **An operation had a working buffer path *and* a working list path, and quietly
> took the second whenever the two representations met.**

The packing threshold judges each value **in isolation**, but a binary
operation's cost is set by its **largest** operand. A 32-element vector,
correctly left unpacked, forced a 20000×40 `dgemv` through the symbolic
evaluator. The rule that resolves it, now applied at seven sites:

> **Pack the small operand up; never materialise the large one down.**

The second recurring finding is about *tests*, not code — see
[`HPC_SWEEP_3_NUMPY_GAP.md`](HPC_SWEEP_3_NUMPY_GAP.md) §"What the tests missed":
three separate assertions passed while describing behaviour that was wrong,
because each compared numeric *values* where the defect was in element *heads*,
or used data below the packing threshold so the path under test never ran.

## See also

- [`docs/design/performance.md`](../design/performance.md) — the combined
  three-system table across all 51 benchmarks.
- [`docs/design/compile.md`](../design/compile.md) — the compiler's design.
- [`docs/design/packed_arrays.md`](../design/packed_arrays.md) — the packing model.
- [`plans/HPC_IMPROVEMENT_PLAN.md`](../../plans/HPC_IMPROVEMENT_PLAN.md) — what
  each sweep left open.

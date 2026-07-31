# Mathilda vs Mathematica vs NumPy on classical HPC kernels

A like-for-like comparison of Mathilda against Wolfram Language 14.0 — and, from
the third sweep on, against NumPy 2.4.4 — on well-known high-performance-computing
workloads: dense linear algebra, FFT, STREAM-style array traffic, PDE stencils,
compiled scalar kernels, integer and combinatorial search, arbitrary-precision
arithmetic, and eight pipeline-shaped kernels taken from real applications.

Companion to [`compile.md`](compile.md) (the compiler's design) and
[`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) (one worked problem
in depth). This document is the broad sweep: many kernels, shallow each.

Reproduce with [`comparisons/hpc_bench.py`](../../comparisons/hpc_bench.py).

---

## 1. Method

**The machine.** Intel Core i9-9880H (8 cores / 16 threads, AVX2), 16 GB,
macOS 15.7.4. Mathilda at `-O3` with GMP 6.3, MPFR 4.2.2, FLINT 3.6, FFTW 3.3.11
and Apple Accelerate for BLAS/LAPACK. Mathematica 14.0.0. NumPy 2.4.4 with SciPy
1.17.1 on CPython 3.11, also linked against Accelerate. Measured 2026-07-30,
extended 2026-07-31.

**Why three systems.** Mathematica says whether Mathilda is behind a competitor.
NumPy — running the *same* Accelerate BLAS — says whether it is behind the
machine. On the dense rows all three call identical kernels, so any spread is
pure overhead; on the array rows NumPy is a thin memory-bandwidth reference.
Where NumPy has no equivalent (arbitrary precision, symbolic rule dispatch,
`PrimePi`) the column reads `—`, and where the honest Python answer is a scalar
CPython loop or mpmath rather than NumPy, the row says so — numba is not
installed on this host, and those rows must not be read as library comparisons.

**Elapsed wall clock, via `AbsoluteTiming` in both CAS and `time.perf_counter`
in Python.** Not `Timing[]`,
which reports CPU time summed over threads: on this 8-core host
`Timing[A . A]` for a 600×600 matrix reads 10.5–13.7 ms against a true 1.9 ms,
because Accelerate's `dgemm` is threaded. Any comparison of a threaded path
against a serial one using `Timing[]` is off by roughly the core count.

`AbsoluteTiming` did not previously exist in Mathilda and was added for this
document; `AbsoluteTime[]` has one-second resolution and is no use for
benchmarking.

**Minimum of several repetitions**, after one untimed warm-up, with the
*maximum* also recorded. The max is a tripwire: a builtin with an internal cache
answers the second repetition in ~0 s, and reporting the minimum alone would make
the row vacuous. Where a cache could not be defeated by varying the argument, the
first cold call is reported instead and labelled.

**Mathilda and Mathematica run the same source text**; the NumPy column is a
separate implementation of the same algorithm, in the same order. Each benchmark
with a cheap scalar answer has that answer recorded from every system and
compared — a timing row is only meaningful once they agree, and the check pins
the *algorithm*, so the columns cannot silently be timing three different
computations. Every value in the tables below agreed.

**Values must be forced.** A benchmark that discards its result can measure the
discarding. Wolfram evaluates `(p + k) q;` on two 10⁶-bit integers in 138 µs and
`IntegerLength[(p + k) q]` in 2.96 ms — the first does not materialise the
product. Every row was re-checked in both forms; only bignum arithmetic differed,
but the check is the point, not the outcome.

**Idiomatic code, default settings.** Both systems automatically compile numeric
`Table`/`Do` bodies and both automatically store dense numeric lists as machine
buffers, so this is like against like. Where a benchmark is written with an
explicit `Compile[]` it is because that is how the kernel would actually be
written for performance, and both systems get the same treatment.

---

## 2. Dense linear algebra (BLAS / LAPACK)

| Benchmark | Mathilda | Mathematica 14.0 | |
|---|---:|---:|---|
| Matrix multiply, 1000×1000 | 8.09 ms | **6.54 ms** | 1.24× |
| `Det`, 500×500 | **1.45 ms** | 1.58 ms | 1.09× faster |
| `LinearSolve`, 1000×1000 | 15.7 ms | **6.45 ms** | 2.44× |
| `Inverse`, 500×500 | 7.37 ms | **3.03 ms** | 2.43× |
| `SingularValueDecomposition`, 300×300 | 51.5 ms | **8.62 ms** | 5.97× |
| `Eigenvalues`, 300×300 symmetric | 22.5 ms | **3.14 ms** | 7.18× |
| `QRDecomposition`, 500×500 | 59.3 ms | **4.70 ms** | 12.6× |

`Det` is now *faster* than Mathematica and matrix multiply is close, because both
systems reach the same Accelerate kernels and the conversion into them is no
longer element-wise (plan phase 3: `na_load_matrix` converted row-major to
column-major one `ndt_get` at a time; it is now a cache-blocked transpose, or a
`memcpy` where no transpose is needed). The spread widens exactly where Mathilda stops using LAPACK:
`QRDecomposition`, `Eigenvalues` and `SingularValueDecomposition` run in-house
numeric kernels. `Eigenvalues` uses Mathilda's own QR iteration, kept in
preference to LAPACK because the eigenvalue *ordering* convention (|λ| ties broken
by position) cannot be reproduced from LAPACK output without risking parity —
a deliberate trade recorded in
[`NDARRAY_REDUCTIONS_COMPARISON.md`](../../comparisons/NDARRAY_REDUCTIONS_COMPARISON.md).
`QRDecomposition` was 18.7× and is now 14.0×: it now routes to `dgeqrf`+`dorgqr`
(plan phase 1.3), but LAPACK's own factorisation is only a few milliseconds of
that. The rest is the boundary — `na_load_matrix`/`na_build_matrix` convert
element by element (phase 3) — plus ~24 ms materialising the `{q, r}` pair to
plain Lists, which the no-nesting invariant requires and Mathematica does not.
That is the case for doing phase 3 before the remaining decomposition work.

---

## 3. Spectral

| Benchmark | Mathilda | Mathematica 14.0 | | was |
|---|---:|---:|---|---:|
| `Fourier`, 2²⁰ reals | 24.5 ms | **21.3 ms** | 1.15× | 40.4× |
| `Fourier`, 2²⁰ complex | 24.7 ms | **18.9 ms** | 1.30× | 80.2× |

**Fixed 2026-07-30** — this was the worst result in the document at 40× and 80×,
and is now within 1.15×/1.30×. Three causes, in order of size
([plan phase 1.1](../../plans/HPC_IMPROVEMENT_PLAN.md)):

1. **`Fourier` was not packed-aware.** `fourier.c` has had a complete NDArray
   fast path all along, but the transparency gate materialised the buffer into
   2²⁰ `Expr` nodes before the builtin could see it. `Fourier[packedList]`
   measured 986 ms against `Fourier[NDArray[...]]` 62.8 ms on identical data.
2. **The b-gather cost an integer division and a modulo per element per axis.**
   `FourierParameters -> {a, b}` folds a `b`-gather onto the transform, but
   `b = 1` — the default, and so nearly every call — makes that gather the
   *identity*. It is now a scale in place: no index arithmetic and one fewer
   full-size allocation and copy. This was the single hottest symbol in the
   process, with more samples than every FFTW symbol combined.
3. **The FFT bounced through a second buffer.** `fftw_plan_dft` now plans and
   executes in place on the caller's array. Safe because `FFTW_ESTIMATE` does not
   touch the arrays while planning, and the plan is built per call on the actual
   pointer, so a SIMD plan's alignment assumption is by construction satisfied.

Plan *caching* was considered and rejected on measurement: `FFTW_ESTIMATE`
planning for 2²⁰ costs 0.15 ms against a 17.6 ms execute. FFTW's own floor for
this transform is 17.6 ms, so 24.5 ms is ~7 ms of marshalling above the library
itself.

---

## 4. Array and memory primitives (n = 10⁷)

| Benchmark | Mathilda | Mathematica 14.0 | NumPy 2.4.4 | vs WL | vs NumPy |
|---|---:|---:|---:|---:|---:|
| Sort | 318.55 ms | 714.17 ms | 155.76 ms | 2.24x | 1/2.05x |
| Riffle with a Real | 118.87 ms | 225.67 ms | 91.13 ms | 1.90x | 1/1.30x |
| PadRight, Real fill | 18.97 ms | 21.85 ms | 17.22 ms | 1.15x | 1/1.10x |
| STREAM triad, a = b + 3 c | 33.30 ms | 31.82 ms | 29.62 ms | 1/1.05x | 1/1.12x |
| Sin (elementwise) | 13.53 ms | 11.71 ms | 52.68 ms | 1/1.16x | 3.89x |
| Exp (elementwise) | 13.37 ms | 11.50 ms | 54.86 ms | 1/1.16x | 4.10x |
| Total (reduction) | 2.95 ms | 2.29 ms | 4.05 ms | 1/1.29x | 1.37x |
| Partition, length 2 | 18.17 ms | 12.52 ms | 15.91 ms | 1/1.45x | 1/1.14x |
| Join (two 10^7 vectors) | 90.19 ms | 59.03 ms | 91.03 ms | 1/1.53x | 1.01x |
| Reverse | 33.05 ms | 11.29 ms | 15.76 ms | 1/2.93x | 1/2.10x |
| Dot (inner product) | 5.81 ms | 5.22 ms | 6.90 ms | 1/1.11x | 1.19x |
| PadRight, default (exact 0) fill | 1.169 s | 466.27 ms | — | 1/2.51x | — |
| RotateLeft | 38.69 ms | 10.09 ms | 16.54 ms | 1/3.83x | 1/2.34x |
| Differences | 16.28 ms | 14.33 ms | 15.60 ms | 1/1.14x | 1/1.04x |
| Accumulate (prefix scan) | 19.42 ms | 15.13 ms | 28.73 ms | 1/1.28x | 1.48x |

Mathilda leads on `Sort` (an 8-pass LSD radix sort on the double bit patterns,
against Mathematica's comparison sort), on `Riffle`, on `PadRight` with an
explicit Real fill, and — because the elementwise kernels are threaded and
NumPy's ufuncs are not — is **3.9× and 4.1× ahead of NumPy** on `Sin` and `Exp`.

`Differences` and `Accumulate` were in the 2–4× band against Mathematica and are
no longer: both went through the `ndt_get`/`ndt_set` dtype choke point, which is
the right default and two indirect calls per element on data that is one
subtraction away from a plain `double` loop. `Accumulate` at rank 1 also emitted
a full loop prologue around a single add, because its inner extent is a runtime
value that happens to be 1. See §10.

`Reverse` and `RotateLeft` remain: they are *serial* buffer work against
vectorised equivalents, memory-bound, and the fix is SIMD and threading in
`ndstruct.c`, not a better algorithm.

`PadRight` with its **default** fill is the odd row, and it is a semantic fact
rather than an inefficiency: `PadRight[{1., 2., 3.}, 5]` is
`{1., 2., 3., 0, 0}` — exact Integer zeros beside Reals — in *both* systems, so
neither can hold the result in a uniform buffer and both fall back to a boxed
list. Given an explicit Real fill Mathilda is the faster of the two.

---

## 5. Stencils and PDE relaxation

| Benchmark | Mathilda | Mathematica 14.0 | |
|---|---:|---:|---|
| Jacobi 5-point relaxation, 512², 100 sweeps | 147 ms | **114 ms** | 1.29× |
| Game of Life, 256², 100 generations | **260 ms** | 89.9 ms | 2.90× |

The Jacobi row is the classical vectorised stencil,

```mathematica
jac[u_] := (RotateLeft[u,{1,0}] + RotateRight[u,{1,0}]
          + RotateLeft[u,{0,1}] + RotateRight[u,{0,1}])/4.;
Nest[jac, u0, 100]
```

and at 1.20× it is essentially at parity. It reached that from **21.6 s** during
this exercise — see [§9](#9-what-the-first-sweep-changed).

**Game of Life is the second-worst result here, and the cause is one missing
kernel.** The idiomatic vectorised Life step counts neighbours with shifted sums
and then selects with `UnitStep`:

```mathematica
life[g_] := With[{k = nb[g]},
  UnitStep[k-3] UnitStep[3-k] + UnitStep[k-2] UnitStep[2-k] g];
```

`UnitStep` had **no NDArray kernel** and cost ~500 ns/element — 5.0 s for a
single `UnitStep` over 10⁷ elements. It could not be given an ordinary real
kernel, because it answers with *exact Integers* (`UnitStep[{-1., 1.}]` is
`{0, 1}`, not `{0., 1.}`) and a float64-closed kernel would change an element's
head. A **narrowing** category — real in, `NDT_INT64` out — was added
([plan phase 2](../../plans/HPC_IMPROVEMENT_PLAN.md)) and `UnitStep` is now 113×
faster (44 ms), with `Floor`, `Ceiling`, `Round`, `IntegerPart` and `Sign` off
`NOT_AWARE` as well.

**That did not move this row at all** — Life was never `UnitStep`-bound. Chasing
it found three more things, the first of which was a wrong answer:

1. **`RotateLeft[list, i]` with a symbolic `i` returned the list unrotated.**
   `rotate_rec` defaulted its per-axis amount to 0 and never rejected a
   non-integer spec, so it silently rotated by zero where Mathematica — and
   Mathilda's own `RotateRight` — leave the call unevaluated. `Sum`'s closed-form
   stage then saw a body with no `i` dependence, concluded it was constant, and
   returned `9 m` instead of the neighbourhood sum:
   `Sum[RotateLeft[{1,2,3},i],{i,0,2}]` gave `{3,6,9}` where Mathematica gives
   `{6,6,6}`. **This benchmark had been measuring a wrong computation all along**,
   in both the packed and unpacked paths.
2. **`Sum` attempted its closed form on array-valued bodies at fixed cost.** The
   polynomial/geometric stages take the body apart symbolically, so the cost
   scales with the *body*, not the range: a **one-term** `Sum` over a 256² grid
   cost 197 ms. It is now skipped for a short range over an array body, where the
   expansion is bounded and cheap — but kept for long ranges, where it genuinely
   works (`Sum[m, {i,1,10^5}]` for constant `m` answers in 0.4 ms).
3. **The packed-array DownValue exemption covered only float64.** A user helper
   `f[x_] := …` binds opaquely and reads no element, so it was exempted from
   materialising — but `packed_int64_ok` was never extended to match, and Life's
   grid is *integer*. `probe[q_] := NDArrayQ[q]` answered `False` for a packed
   integer argument and `True` for a real one.

Together: **65.7 s → 260 ms, 253×**, now 2.90× off Mathematica, and Mathilda and
Mathematica agree on the answer (both report 320 live cells after 15 generations
of the same 40×40 grid — the check that should have been there from the start).

---

## 6. Scalar kernels via `Compile[]`

| Benchmark | Mathilda | Mathematica 14.0 | Python | vs WL | vs Python |
|---|---:|---:|---:|---:|---:|
| Lennard-Jones energy, 1452 bodies (all pairs) | 133.19 ms | 319.51 ms | 91.25 ms | 2.40x | 1/1.46x |
| Mandelbrot, 800x800, 100 iterations | 769.23 ms | 1.621 s | 980.20 ms | 2.11x | 1.27x |
| Monte Carlo pi, 10^7 samples (vectorized) | 162.36 ms | 250.40 ms | 204.53 ms | 1.54x | 1.26x |
| Logistic map, 10^7 iterations | 183.76 ms | 254.36 ms | 460.57 ms | 1.38x | 2.51x |

**Mathilda's compiled scalar code is faster than Wolfram's on all four compiled
kernels**, by 1.4–2.4×. This is the same result
[`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) reports for a
stencil (1.8–2.1×), now confirmed on four unrelated bodies: a nested all-pairs
`For` loop with indexed array reads (Lennard-Jones), an escape-time loop with a
compound `&&` guard (Mandelbrot), a bare arithmetic recurrence (the logistic
map), and a vectorised Monte Carlo.

The Python column is **vectorised NumPy** for Lennard-Jones (`triu_indices` over
all pairs), Mandelbrot (a grid-wide escape-time iteration) and Monte Carlo π, and
a **CPython loop** for the logistic map, which is a scalar recurrence with no
vectorised form (numba is not installed on this host). So the first three rows
are genuine library comparisons and Mathilda leads two of them; Lennard-Jones is
where NumPy's vectorised all-pairs form wins, at 91 ms against 133 ms. The
logistic-map row is the one that is not a library comparison — read 2.51× as
"what a Python user without a JIT measures".

Monte Carlo π was the `UnitStep` gap of §5, then became RNG-bound at 2.73 s — of
which `RandomReal[{0,1}, 10^7]`, called twice, was 2.4 s. See
[`RANDOM_NUMBER_GENERATION.md`](../experiments/RANDOM_NUMBER_GENERATION.md) for
what that turned out to be and what fixing it cost.

---

## 7. Integer, combinatorial and arbitrary precision

| Benchmark | Mathilda | Mathematica 14.0 | |
|---|---:|---:|---|
| Collatz longest chain below 10⁶ | **4.65 s** | 6.98 s | 1.50× faster |
| Sieve of Eratosthenes to 10⁷ | **666 ms** | 787 ms | 1.18× faster |
| `50000!` (exact) | 2.62 ms | **2.26 ms** | 1.16× |
| `3^1000000` | 3.13 ms | **2.50 ms** | 1.25× |
| Naive recursive `fib(25)` | 368 ms | **134 ms** | 2.74× |
| π to 100,000 digits | 36.6 ms | **12.8 ms** | 2.86× |
| `PrimePi[10^9]` | 14.7 ms | **4.26 ms** | 3.45× |
| Product of two 10⁶-bit integers | 3.37 ms | **2.96 ms** | 1.14× |

Both compiled integer searches (sieve, Collatz) are faster in Mathilda. Bignum
factorial and power are at parity — both systems are calling a well-tuned GMP-class
library.

**Naive recursive `fib`** at 2.74× measures pure rule-dispatch overhead:
`fib[k_] := fib[k-1] + fib[k-2]` is 243k pattern matches and Mathilda's matcher
is slower per call.

The 10⁶-bit multiply row nearly became this document's biggest wrong claim. Timed
as `(p + k) q;` Wolfram reported **138 µs** against Mathilda's 3.37 ms — a 24×
gap. It is an artifact: with the value discarded, Wolfram does not materialise a
bignum product. Forcing it (`IntegerLength[(p + k) q]`) gives **2.96 ms**, and
the two systems are at parity. Independent confirmation: `mpz_mul` on exactly
these operands, in a standalone C program, takes **2.93 ms** — Mathilda adds 5%
over GMP, and 21×-faster-than-GMP was never a plausible number. FLINT's
`fmpz_mul` measures 2.92 ms, i.e. it delegates to GMP at this size, so there is
no faster library already linked to switch to.

The same check on every other row found nothing: for array operations Wolfram's
discarded and forced timings agree within noise, as do `50000!`, `3^1000000` and
`N[Pi, 100000]`. The laziness is specific to bignum arithmetic. **Discarding a
result is not always free, and is sometimes too free.**

Two rows needed care to measure at all. `PrimePi[10^9]` and π to 100k digits both
memoise so hard that no argument variation defeats the cache: Wolfram answered
repeat calls to `PrimePi[10^9 + k]` in 1 µs, which is not a computation. Both
columns above are cold first calls in a fresh process. Mathilda's `PrimePi` is a
Deléglise–Rivat implementation and 14.7 ms for π(10⁹) is a genuine computation;
both systems agree the answer is 50847534.

---

## 8. Second sweep — five more kernels (2026-07-31)

The first 38 benchmarks covered dense linear algebra, spectral transforms,
streaming array traffic, stencils, compiled scalar loops and bignum arithmetic.
Five kernels were added to run somewhere none of those did: a Krylov iterative
solver, direct convolution, ODE time integration, an irregular hash-keyed
reduction, and interpolation.

| Benchmark | Mathilda | Mathematica 14.0 | | was |
|---|---:|---:|---|---:|
| `Interpolation`, 10⁴ nodes, 10⁴ evaluations | **7.65 ms** | 19.9 ms | 2.60× faster | 1.494 s |
| `NDSolve`, Lorenz system to t = 200 | **31.5 ms** | 47.8 ms | 1.52× faster | — |
| Conjugate gradient, 2D Poisson 256², 100 iterations | 155 ms | **84.0 ms** | 1.84× | — |
| `Tally`, 10⁷ integers into 10⁴ bins | 129 ms | **17.9 ms** | 7.21× | 1.230 s |
| `ListConvolve`, 1024² image, 5×5 kernel | 312 ms | **33.7 ms** | 9.27× | 420 ms |

All five are from one full `hpc_bench.py` run, so they share a kernel session with
the 38 above — which is how the `NDSolve` row turned out to be measuring the wrong
thing; see the last subsection.

Two of the five are faster than Wolfram. **`NDSolve` is the one that was faster
without any work here** — adaptive-step integration of the Lorenz system to
t = 200 runs in 31.5 ms against 47.8 ms, which is the RHS compiler paying off on
a workload it was built for.

The other three each found something, and the three "was" figures are what this
section changed.

### `Interpolation` evaluation was quadratic in its own table

Evaluating an `InterpolatingFunction` cost time proportional to the number of
nodes in its table — 151 µs per point at 10⁴ nodes, 906 µs at 5×10⁴ — so
resampling a table at as many points as it has nodes was O(n²). At 10⁵ nodes the
benchmark could not be timed at all: it took 228 s.

Four separate per-point costs, none of which depends on the point being
evaluated. Three were table-wide facts recomputed per call (`fill_values`
rebuilding the whole value tensor, `table_Ksupplied`/`common_scan_inexact`/
`numeric_expr_is_mpfr` each walking every node, and the grid cache testing its
own hit with an O(n) `expr_eq`); all three are now memoised on the table's
identity, which is sound because `expr_copy` is a refcount bump and the memo
holds a reference.

The fourth was a genuine O(n²): `build_grid` inserted coordinates one at a time
into a sorted array, scanning from the front. Coordinates that **arrive in
order** — which is every table from `Table[{x, f[x]}, {x, a, b, dx}]` — walk the
whole grid built so far on every insertion. Sort-and-unique instead: the
160,000-node grid build went from 4.96 s to 21 ms.

**1.494 s → 7.65 ms, 195×**, and from 74× slower than Wolfram to 2.6× faster.

### `Tally` was slower packed than unpacked

`Tally` already hashes, so there was no algorithm to fix. The cost was boxing:
it reads `list->args[i]`, so a packed buffer had to be materialised into one
`Expr` per element before the builtin could begin, and every probe then hashed
and compared `Expr` nodes. Over 10⁶ packed int64 that measured 147 ms against
**96 ms for the identical data held as a plain `List`** — the packed
representation was a net loss, which is the failure mode
[N3](../../plans/HPC_IMPROVEMENT_PLAN.md) exists to make visible.

`ndred_tally` keys on the machine word instead. The instructive part was sizing
the table: for *n* elements the obvious capacity is 2n, and at 10⁷ that is a
2^25-slot table — 268 MB to zero and then random-probe — for 10⁴ live entries.
Sizing it by the *distinct* count instead, discovered by growing, was worth
another 2.7×.

**1.230 s → 129 ms, 9.5×**, closing 69× to 7.2×. What is left is one random
probe per element.

### `ListConvolve` chose the FFT by measuring the wrong thing

The engine choice was a threshold on total work, (output elements) × (kernel
elements) ≥ 4096 ⇒ FFT. That quantity is the cost of the **direct** engine; the
FFT engine's cost does not depend on the kernel size at all. So every large
problem took three full-size complex transforms however small its kernel — and a
small kernel over a large image is precisely what image filtering is. It is now
a comparison between two estimated costs.

Switching to the direct engine made it *worse*, at 1.20 s, which is the useful
part of this result: the direct engine was integer-division bound. Decomposing
the kernel's flat index sat in the innermost loop, four 64-bit divisions per
multiply-add at rank 2 — ~46 ns per MAC where the arithmetic costs ~1 ns. The
decomposition depends only on the kernel position, of which there are 25 rather
than 26 million.

**420 ms → 312 ms**, and no further: at 13 ns per multiply-add what remains is
not arithmetic. `ListConvolve` has no NDArray path, so it materialises 10⁶ input
`Expr` nodes and allocates 10⁶ output ones, and `expr_free` alone was a sixth of
the profile. Making it packed-aware end to end is the outstanding item.

### The harness disagreed with itself about agreement

Every inexact row was reported as a value mismatch on the way out while the table
printed beside it showed none. The table used `values_agree`, which compares to
the shorter of the two printings — Mathilda prints 6 significant figures where
Wolfram prints 17 — and the exit-code check used `vm == vw`. So `0.00182813` and
`0.0018281272300857572` were simultaneously the same answer and a mismatch. The
exit check now uses `values_agree` too.

### A benchmark that measured nothing

The conjugate-gradient row was first written with a smooth source,
`Sin[2 π i/n] Cos[4 π j/n]`. That is an **eigenvector** of the discrete periodic
Laplacian, so CG annihilates it in one iteration: the residual reached 3.6×10⁻²⁴
after 10 iterations and the remaining 90 measured denormal arithmetic. The source
is now a point source with its mean removed — broad spectrum, and the periodic
operator is singular on the constants, so the mean removal is what makes the
problem well-posed at all. Both halves were necessary and only one of them is
obvious.

---

## 9. Third sweep — eight real-world kernels, three systems (2026-07-31)

The first two sweeps could say "Mathilda is N× Mathematica". That answers whether
we are behind a competitor, not whether we are behind *the machine*. This sweep
adds **NumPy 2.4.4** (Python 3.11, SciPy 1.17.1) as a third column. On this host
NumPy links the same Accelerate BLAS Mathilda does, so on the dense rows all
three systems run identical kernels and any spread is pure interpreter overhead;
on the array rows NumPy is a thin, well-optimised memory-bandwidth reference.
That makes it the ruler, and it is the reason most of what follows was found.

The eight kernels are chosen for **depth** rather than breadth: each is a
pipeline out of a real application, not a single builtin. The second sweep's
conjugate-gradient row is the reason — a composition can be slow in ways none of
its parts are, and only a composition finds that.

| id | kernel | field | what it stresses |
|---|---|---|---|
| `bsmc` | Black–Scholes Monte Carlo, 10⁷ paths | finance | RNG at scale, Box–Muller, reduction |
| `tseries` | EMA + rolling vol + max drawdown, 10⁶ | quant finance | sequential scans |
| `logreg` | Logistic regression, 200000×32, 100 GD steps | ML | matrix–vector both ways, loop-carried state |
| `kmeans` | k-means, 100000 points, 8 dims, k = 16, 20 its | ML | masked reductions, argmin, scatter-as-dot |
| `nbody` | All-pairs gravity, 1024 bodies, 10 Verlet steps | astrophysics | `Outer`, megabyte intermediates |
| `heat3d` | 3D 7-point stencil, 128³, 50 steps | CFD | rank-3 shifts, 16 MB working set |
| `psd` | Welch PSD, 2²² samples, 1024 blocks of 4096 | signal processing | many medium FFTs, windowing |
| `imgpipe` | Separable Gaussian blur + Sobel, 1024² | vision | rank-2 `ListCorrelate`, `Transpose` |

Every kernel has a deterministic scalar check that all three systems agree on.
That is a stronger constraint than it sounds: the check pins the *algorithm*, so
the three columns cannot silently be timing three different computations. None
disagreed.

### Results

| Benchmark | Mathilda | Mathematica 14.0 | NumPy 2.4.4 | vs WL | vs NumPy |
|---|---:|---:|---:|---:|---:|
| Black–Scholes Monte Carlo, 10⁷ paths | 565.03 ms | 349.05 ms | 559.82 ms | 1/1.62× | 1/1.01× |
| Return series: EMA + vol + drawdown, 10⁶ | 2.021 s | 149.36 ms | 30.38 ms | 1/13.53× | 1/66.52× |
| Logistic regression, 200000×32, 100 steps | 2.706 s | 2.145 s | 839.44 ms | 1/1.26× | 1/3.22× |
| k-means, 100000×8, k = 16, 20 its | 2.261 s | 7.326 s | 765.45 ms | **3.24×** | 1/2.95× |
| N-body all-pairs gravity, 1024 bodies | 438.89 ms | 6.982 s | 396.36 ms | **15.91×** | 1/1.11× |
| 3D heat equation, 128³, 50 steps | 3.630 s | 2.545 s | 2.250 s | 1/1.43× | 1/1.61× |
| Welch PSD, 2²² samples, 1024 blocks | 149.66 ms | 90.15 ms | 93.89 ms | 1/1.66× | 1/1.59× |
| Gaussian blur + Sobel edges, 1024² | 2.185 s | 96.96 ms | 112.95 ms | 1/22.53× | 1/19.34× |

Two of the eight started out unmeasurable: k-means and logistic regression did
not finish, and N-body took 55.5 s. What that cost, and what fixing it was
worth, is below.

### One defect, seven disguises

Six of the seven fixes this sweep produced are the same bug:

> The packing decision is made about **one** value in isolation, but the cost is
> paid in proportion to **another**.

`PACK_MIN_ELEMENTS` is 250, chosen for blast radius rather than for break-even,
and that is right. What was wrong is that a 32-element vector — correctly left
unpacked — could drag a 6.4-million-element matrix onto the symbolic path with
it. Every site below had a working buffer path and a working List path, and
quietly took the second whenever the two met.

| what | measured | after |
|---|---:|---:|
| `Dot`, one operand a plain List (20000×40 by 40) | 320 ms | **0.31 ms** |
| `Total[m, {1}]` vs the identical `Total[m]` (1000²) | 258 ms | **1.28 ms** |
| `Outer[Subtract, v, v]` (1024², float64) | 808 ms | **4.03 ms** |
| rank-2 ⊖ rank-1 threading (8×100000) | 322 ms | **2.52 ms** |
| `a + b`, one plain operand (10⁶) | 418 ms | **50.5 ms** |
| the N-body step re-applied to its own output | 5.75 s | **43 ms** |
| `MapThread[Min, …]` (16×100000) | 365 ms | **3.87 ms** |

The rule that came out of it, now applied at `Dot`, at `Plus`/`Times`, and at
the evaluator's Listable gate: **pack the small operand up; never materialise the
large one down.** Packing 40 doubles to save 800,000 symbolic multiplies is never
a bad trade, and packing is value-preserving by contract, so there is nothing to
weigh.

The largest single win was not a kernel at all but a **return statement**. A
function that returns several arrays destroyed all of them on the way out —
`{px, py, pz, vx, vy, vz}` came back as six plain Lists — and the caller paid on
the *next* operation, not the one it could see. `List` cannot be made
packed-aware (a List node holding `EXPR_NDARRAY` elements is exactly the
malformed shape the gate exists to prevent), so the fix is to pack the whole
thing: `pack_sniff` already absorbs already-packed rows, so a list of n packed
vectors *is* a rank-2 buffer.

Three of these needed the head added to `pack.c`'s aware list before the fast
path could be reached at all — the gate had already materialised the arguments.
That is now the third time (`Nest`, `Fourier`, and here `Outer`), which is
exactly what `make check-packed-aware` exists to catch and did.

### What NumPy was worth as a control

Three rows would have read as acceptable against Mathematica alone and do not
against NumPy:

- **`Sieve of Eratosthenes`** — 1.20× *ahead* of Mathematica, 26.9× behind
  NumPy. Both CAS run a compiled scalar loop; NumPy runs `s[i*i::i] = False`, a
  strided memset. The gap is an algorithmic idiom neither CAS expresses, not an
  execution gap.
- **`tseries`** — 13.5× behind Mathematica but 66.5× behind NumPy, because
  SciPy's `lfilter` and `np.maximum.accumulate` do a scan in C and neither CAS
  has a vectorized scan at all. `Accumulate` covers the prefix sum; an EMA is a
  general linear recurrence and a drawdown needs a running maximum, and both are
  `FoldList`, i.e. one interpreted call per element.
- **`Fibonacci`/`pi`** — NumPy's columns here are CPython recursion and mpmath,
  and are labelled as such in the harness. They are not NumPy results and should
  not be read as a library comparison.

Conversely, three rows say Mathilda is at the machine's limit and the remaining
gap to Mathematica is small: STREAM triad (1.03× / 1.18×), `Total` (1.15× /
1.49× *ahead*), `Dot` (1.07× / 1.11×).

### Still open

- **`ListConvolve`/`ListCorrelate` have no NDArray path.** The image pipeline is
  22.5× behind Mathematica and 19.3× behind NumPy, and it is all boxing: 10⁶
  input `Expr` nodes materialised and 10⁶ output ones allocated, four times over.
  This is the single largest remaining gap and is HPC plan item 7.1.
- **No vectorized scan.** `FoldList` with a pure function is the only way to
  write an EMA or a running maximum, and it costs one interpreted evaluation per
  element.
- **k-means is still 2.95× behind NumPy** after a 6× improvement, and
  `logreg` 3.22×; both are now composition costs spread thinly rather than one
  hot spot.

---

## 10. Fourth sweep — closing the distance to NumPy (2026-07-31)

The third sweep left two rows badly behind NumPy — the return series at 66× and
the image pipeline at 19×. Neither is a linear-algebra row, so BLAS was not the
answer. This sweep probed **primitives** directly against their NumPy
equivalents rather than running applications, on the hypothesis that the
remaining cost was structural rather than arithmetic.

The elementwise kernels came back at the memory floor (`v v` 1.8× NumPy, `Log`
1.3×), which settled it: nothing above them was about arithmetic. Everything
expensive was either **structural** — an O(1) or memcpy operation walking a
boxed list — or a **scan**.

### The probe, 10⁶ float64, before

| op | Mathilda | NumPy | ratio |
|---|---:|---:|---:|
| `Most[v]` | 270 ms | 0.48 ms | 563× |
| `Rest[v]` | 223 ms | 0.48 ms | 466× |
| `FoldList[a #1+b #2 &, …]` packed | 1246 ms | 3.02 ms (`lfilter`) | 413× |
| `FoldList[Max, …]` | 605 ms | 2.73 ms | 222× |
| `ListConvolve[k, v]` | 263 ms | 1.40 ms | 188× |
| `First[v]` | 123 ms | ~0 | O(1) done in O(n) |
| `Clip[v, {a,b}]` | 356 ms | 4.66 ms | 76× |
| `Differences[v]` | 6.6 ms | 0.41 ms | 16× |
| `ListCorrelate[k5×5, im]` 1024² | 376 ms | 52.5 ms | 7.2× |
| `Accumulate[v]` | 7.9 ms | 2.09 ms | 3.8× |

### What was fixed

| fix | before | after |
|---|---:|---:|
| `First` / `Last` read the buffer instead of materialising it | 123 / 103 ms | **0.01 / 0.00 ms** |
| `Most` / `Rest` as a leading-axis `memcpy` | 270 / 223 ms | **3.42 / 0.71 ms** |
| `Clip` back on the buffer, gated on the *bounds* | 356 ms | **0.93 ms** |
| `ndred_scan` — `Plus`/`Times`/`Max`/`Min`, both spellings | 605 ms | **2.00 ms** |
| `numloop_fold_impl` reads a packed argument | 1246 ms | **22.6 ms** |
| `ListCorrelate` 5×5 over 1024² on the buffer | 376 ms | **46.5 ms** |
| `ListConvolve` 5-tap over 10⁶ | 263 ms | **16.1 ms** |
| `Accumulate` float64 arm, rank-1 specialised | 7.9 ms | **3.5 ms** |
| `Differences` float64 arm | 6.6 ms | **0.97 ms** |
| `Part`, contiguous span as a block copy | 14.0 ms | **7.4 ms** |
| `Outer` with one operand below the threshold | 849 ms | **10.9 ms** |

Four of those are the third sweep's finding again — *the packing decision is made
about one value in isolation, and the cost is paid in proportion to another* —
and the rest are the first sweep's: a fast path that exists and cannot be
reached, or a per-element choke point on data that is one `memcpy` away.

Two are worth singling out.

**`Clip` had been removed from the aware list entirely** because an exact bound
makes the answer non-uniform (`Clip[{-2., 0., 2.}, {-1, 1}]` is `{-1, 0., 1}`).
The exactness problem was real but belonged to the **bounds**, not the head:
Real bounds are uniform and safe, and an exact bound is safe whenever an
in-range scan proves nothing is clipped, because then the answer just *is* the
input.

**The convolution's cost was an integer division.** The output multi-index was
computed as `(o / Lstr[ax]) % Ldims[ax]` — a 64-bit division per axis per
output, with runtime strides the compiler cannot strength-reduce. On a rank-1
filter that was the whole cost: 10⁶ divisions swamping 5×10⁶ multiply-adds, and
the vectorisable interior loop it was hiding never got to show. Replacing it with
an odometer took the rank-2 correlation past SciPy's `correlate2d`.

### `dgemm` slows the next threaded loop by 1.45×

The Jacobi row reads 223 ms in a full run and 128 ms measured alone, while every
other row is stable between the two. Bisecting the benchmark prefix found a
single trigger — running one 1000×1000 `dgemm` first:

| | Jacobi, 512², 100 sweeps |
|---|---:|
| alone | 130.3, 132.2, 132.9 ms |
| after `matmul` | 188.6, 190.2, 190.5 ms |

Reproducible to ±1%. Accelerate's threaded `dgemm` leaves the process in a state
where Mathilda's own `nd_parallel_for` runs 1.45× slower, most likely its worker
threads persisting and competing for cores. This is not confined to the
benchmark: any program that multiplies a matrix and then does array work pays
it. Recorded as an open item rather than fixed — the fix is a thread-pool policy
question that wants its own experiment. **The two threaded-stencil rows (Jacobi,
Game of Life) follow the linear-algebra group in run order and should be read
with this in mind.**

### Result

| Benchmark | before | after | Mathematica 14.0 | NumPy 2.4.4 | |
|---|---:|---:|---:|---:|---|
| Return series, 10⁶ | 2.02 s | **55.6 ms** | 154 ms | 30.4 ms | **2.78× faster than WL** |
| Gaussian blur + Sobel, 1024² | 2.19 s | **103 ms** | 94.5 ms | 114 ms | **1.10× ahead of NumPy** |
| `ListConvolve`, 1024², 5×5 | 312 ms | **35.7 ms** | 33.8 ms | 39.7 ms | **1.11× ahead of NumPy** |
| `Differences`, 10⁷ | 56.6 ms | **16.3 ms** | 14.3 ms | 15.6 ms | |
| `Accumulate`, 10⁷ | 63.8 ms | **19.4 ms** | 15.1 ms | 28.7 ms | **1.48× ahead of NumPy** |

Per-experiment detail:
[`docs/experiments/HPC_SWEEP_3_NUMPY_GAP.md`](../experiments/HPC_SWEEP_3_NUMPY_GAP.md).

---

## 11. What the first sweep changed

Writing the benchmark found more than it measured. Every item below was
diagnosed from a measurement, fixed, and verified by a differential sweep
(the same expression over a packed buffer and over the identical plain list,
requiring byte-identical printed output).

### A crash, under default settings

`LinearSolve[A, b]` on a 90×90 machine-real system **stack-overflowed**. Three
faults compounded:

1. The linear-algebra heads had never been marked packed-aware, so the
   transparency gate materialised every packed matrix *before*
   `linalg_call_has_ndarray` could see it — which then read false and sent a
   machine-real solve down the exact fraction-free (Bareiss) path, one polynomial
   GCD per pivot, recursing.
2. Automatic packing keys on total element count, so a 90×90 system arrives with
   the matrix packed (8100 elements) and the right-hand side **not** (90). The
   NDArray fast path declined a plain-List rhs "to keep behaviour uniform" — and
   that deferral was the crash. `na_load_vector` already accepted an ordinary
   numeric List.
3. `Method -> Automatic` was a plain alias for `DivisionFreeRowReduction` even
   for an inexact matrix, where being fraction-free buys nothing. Now resolved by
   exactness, so the overflow cannot return via `$AutoArrayPacking = False`.

`LinearSolve` at 1000×1000 now runs in 17.5 ms with a verified residual.

### The structural family was on the marshalling path

`RotateLeft`, `RotateRight`, `Join`, `Partition`, `Differences`, `Riffle`,
`PadLeft` and `PadRight` had no buffer walk: each materialised one `Expr` per
element, ran the generic List implementation, and re-packed. These are the
operations array-style numerical code is written out of.

| operation | before | after | |
|---|---:|---:|---|
| `Differences`, 10⁶ reals | 834 ms | 6.2 ms | 135× |
| `Partition[v, 2]`, 10⁶ | 257 ms | 1.3 ms | 192× |
| `Join`, two 10⁶ vectors | 352 ms | 3.3 ms | 106× |
| `Riffle[v, 0.]`, 10⁶ | 206 ms | 6.9 ms | 30× |
| `RotateLeft`, 512×512 matrix | 42.6 ms | 1.1 ms | 39× |
| `PadRight[v, n, 0.]`, 10⁶ | 114 ms | 0.48 ms | 237× |

A rotation permutes contiguous blocks, so all of this is `memcpy` work. The care
went into the *exactness* rule rather than the copying: a head that introduces a
new element (`Riffle`'s separator, `Pad`'s fill) may only use the buffer when that
element is exactly representable at the buffer's dtype with a matching head, which
is why `PadRight[x, n]` with the default exact `0` still declines.

### `Nest` materialised its state every iteration

`Nest`, `NestList`, `NestWhile`, `NestWhileList`, `FixedPoint` and
`FixedPointList` were not packed-aware, although `Fold` and `FoldList` were and
the implementations are the same shape — `nest_impl` copies the value and applies
a function, never inspecting its structure.

The 512² Jacobi sweep body evaluates in **1.10 ms** standalone. `Nest[f, u0, 10]`
took **2.19 s**, against **0.017 s** for the identical `Do[u = f[u], {10}]` loop,
and answered with an unpacked list. Now 0.0185 s — **118×**.

### Defining a helper function cost 100×

Even after that, the Jacobi benchmark ran at 21.6 s, because it is written the way
anyone would write it:

```mathematica
jac[u_] := (RotateLeft[u,{1,0}] + ...)/4.;
```

A user symbol has no `packed_aware` bit, so `jac[packedArray]` materialised —
correctly in general, since the matcher cannot descend a buffer to match
`f[{a_, b_}]`. But the shape that appears in numerical code is a *bare* pattern
variable, bound and substituted whole, and nothing looks inside the value.

The gate now exempts a head whose every DownValue binds opaquely: each top-level
argument of every rule's left-hand side is exactly `Pattern[sym, Blank[]]`, with
no head restriction, no test, no condition and no nested structure. Anything else
still materialises. **21.6 s → 0.137 s, 158×**, and the identical body written as
a pure function — which the gate already exempted — agrees to the last bit.

### Two-argument `If` was outside the compilable subset

`If[test, var = val]` — a conditional store, how anyone writes a running maximum —
had no lowering, and the compilable subset is a *cliff, not a slope*: one head
outside it costs the whole body. Both the Sieve of Eratosthenes and the Collatz
search silently fell back to the interpreter on this single form.

`CompileDiagnostics` named it exactly: `"Compiled" -> False`,
`"Subexpression" -> "If[len > bl, bl = len]"`. It lowers like `While` — guard,
jump over the body, answer 0 — and joins `stmt_valued_head` so that 0 can never
be a program's result. `compile.md` §11 had specified it; it was never emitted.

- Collatz to 10⁶: **240 s → 4.65 s** (52×)
- Sieve to 10⁷: did not compile at all → **666 ms**, answering 664579 = π(10⁷)

### `Norm` of a matrix

`Norm[m]` computed for a packed matrix (LAPACK `dlange`, `gesdd` for the spectral
norm) and returned *unevaluated* for the same value held as a plain list. The
plain path now routes an inexact numeric matrix through the same kernel, matching
Mathematica to the last digit on every test matrix. Exact matrices still stay
symbolic: `Norm[{{1,2},{4,5}}]` is `Sqrt[23 + 2 Sqrt[130]]`, which Mathilda has no
symbolic SVD to produce, and answering 6.76783 to an exact question would be
worse than not answering.

### A pre-existing wrong compiled answer, and one bad benchmark

`Inverse` under `OneStepRowReduction` and `CofactorExpansion` answers
`Inverse[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}]` with exact Integer `1`s among the
Reals, where `DivisionFreeRowReduction` and Mathematica both give `1.0`. Not
introduced here and not on any default path, so it is recorded in
`src/linalg/inv.c` rather than papered over — it is why `Inverse`'s `Automatic`
deliberately still selects `DivisionFreeRowReduction`.

And one failure was the benchmark's own. The Lennard-Jones point set was built
from `Mod[a k, m]`, which cycles — so with 3n coordinates it repeats *points*,
two coincident bodies make r = 0, the 1/r¹² term overflows, and the compiled call
correctly bails to the interpreter on the non-finite result. That read as an 89×
performance cliff between n = 1000 and n = 1200 and answered `Indeterminate`. The
point set is now a cubic lattice. **A benchmark can be wrong in the same ways a
program can.**

---

## 12. Summary

Counts from the full three-system run of 2026-07-31, after the fourth sweep, and
computed from the run's own JSON rather than by eye.

Of 51 benchmarks, Mathilda is **faster than Mathematica on 15**, within 1.5× of
it on a further 18, and more than 1.5× behind on 17 (`PrimePi` has no comparable
Mathematica timing). Against NumPy, on the 49 rows where a NumPy column exists:
**faster on 20**, within 1.5× on 10, and more than 1.5× behind on 19.

So on two thirds of the suite Mathilda is now ahead of, or within 1.5× of, both
systems.

**Where Mathilda is ahead** it is consistently in *compiled* code — all four
compiled scalar kernels, both compiled integer searches — plus `Sort`, the
structural buffer ops, `NDSolve` (1.79× Mathematica and **36× SciPy**),
`Interpolation` (3.13×), and now the whole convolution and scan family. The
compiled-code advantage of 1.4–2.4× matches what `COMPILE_EXAMPLE.md` measures
independently on a stencil.

**Where Mathilda is behind**, the causes are specific and named, not diffuse.
The list is shorter than it was; items struck through were closed by §10.

1. ~~**No narrowing (float→int64) kernel category.**~~ Closed 2026-07-30 — see
   [`MACHINE_INTEGERS.md`](../experiments/MACHINE_INTEGERS.md).
2. **`QRDecomposition`, `Eigenvalues`, `SVD` bypass LAPACK** — 6–14×. Partly
   deliberate: `Eigenvalues`' ordering convention cannot be reproduced from
   LAPACK output without risking parity.
3. **Serial buffer ops** — now just `Reverse` and `RotateLeft`, in a 2.9–3.8×
   band. `Differences`, `Accumulate` and `Dot` left it in §10.
4. **Pattern-match dispatch** — 2.65× on naive recursive `fib`, and the one row
   where CPython beats both CAS.
5. ~~**Heads with no NDArray path pay for boxing** (`ListConvolve`).~~ Closed in
   §10: `ListConvolve` is now 1.05× Mathematica and **1.11× ahead of NumPy**.
   `Tally` remains, at 7.17×.
6. ~~**No vectorized scan.**~~ Closed in §10. The running maximum is now faster
   than `np.maximum.accumulate`; the general linear recurrence is 7.5× `lfilter`,
   down from 413×.
7. **`Transpose` copies where NumPy views.** This is now the largest application
   gap: logistic regression is 3.06× NumPy, and the reason is that
   `Transpose[X]` on a 200000×32 matrix costs 27 ms and the kernel does it 100
   times, because it is loop-invariant but re-evaluated. Closing it needs a
   strided NDArray view that every consumer honours — a design change, not a
   fast path.
8. **No interpreter-level fusion.** k-means at 2.92× NumPy is three passes over
   6.4 MB where one would do. `Compile[]` fuses; ordinary array code does not.
9. **The integer band.** The sieve is 1.18× *ahead* of Mathematica and 25×
   behind NumPy — the clearest illustration of why the third column exists.
10. **`dgemm` slows the next threaded loop by 1.45×** (§10). A thread-pool
    interaction, not an algorithm.

Arbitrary-precision arithmetic is *not* on this list: bignum multiply, factorial
and power are all at parity, because both systems call GMP.

None of these is a design problem except items 7 and 8, which are honest design
changes rather than unfilled paths.

**The methodological point worth keeping.** The largest wins across all four
sweeps — 56000×, 658×, 563×, 466×, 413×, 253×, 195× — were, without exception,
cases where an operation *worked correctly and quietly took the slow path*.
Nothing failed, no test broke, and no amount of reading the code would have found
them. They were found by measuring one operation at a time against the same
operation in another system.

The corollary, which cost time three separate times: **a test that passes can
still be describing the wrong behaviour.** Three assertions in this codebase
enshrined bugs — one compared numeric distance where the defect was in element
*heads*, one asserted a fallback as though it were an invariant, and two used
data below the packing threshold so the path under test never ran. All three
passed for months.

---

## 13. Reproducing

```bash
make -j$(sysctl -n hw.ncpu)
python3 comparisons/hpc_bench.py                       # full run, all three systems
python3 comparisons/hpc_bench.py --only fft,sort       # a subset, by id
python3 comparisons/hpc_bench.py --scale 0.05          # smaller sizes, smoke test
python3 comparisons/hpc_bench.py --system mathilda     # one system
python3 comparisons/hpc_bench.py --system mathilda,numpy
```

The harness prints the markdown table above, writes raw JSON with `--json`, and
exits nonzero if any benchmark's recorded answer differs between systems. Set
`WOLFRAMSCRIPT` if `wolframscript` is not at the default macOS path, and
`HPC_PYTHON` for a Python that has NumPy and SciPy (it defaults to
`/usr/local/bin/python3.11`, and is deliberately not the interpreter running the
harness).

## See also

- [`plans/HPC_IMPROVEMENT_PLAN.md`](../../plans/HPC_IMPROVEMENT_PLAN.md) — the
  plan to close every gap in §10, sequenced by value and risk.
- [`compile.md`](compile.md) — the compiler's design and milestones.
- [`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) — one HPC problem
  end to end, with the bytecode.
- [`packed_arrays.md`](packed_arrays.md) — the automatic packing model that most
  of §10 is about.
- [`NDARRAY_REDUCTIONS_COMPARISON.md`](../../comparisons/NDARRAY_REDUCTIONS_COMPARISON.md)
  — Mathilda vs NumPy vs Mathematica on reductions and elementwise kernels.

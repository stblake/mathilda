# Mathilda vs Mathematica on classical HPC kernels

A like-for-like comparison of Mathilda against Wolfram Language 14.0 on
well-known high-performance-computing workloads: dense linear algebra, FFT,
STREAM-style array traffic, PDE stencils, compiled scalar kernels, integer and
combinatorial search, and arbitrary-precision arithmetic.

Companion to [`compile.md`](compile.md) (the compiler's design) and
[`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) (one worked problem
in depth). This document is the broad sweep: many kernels, shallow each.

Reproduce with [`comparisons/hpc_bench.py`](../../comparisons/hpc_bench.py).

---

## 1. Method

**The machine.** Intel Core i9-9880H (8 cores / 16 threads, AVX2), 16 GB,
macOS 15.7.4. Mathilda at `-O3` with GMP 6.3, MPFR 4.2.2, FLINT 3.6, FFTW 3.3.11
and Apple Accelerate for BLAS/LAPACK. Mathematica 14.0.0. Measured 2026-07-30.

**Elapsed wall clock, via `AbsoluteTiming` in both systems.** Not `Timing[]`,
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

**Both systems run the same source text**, and each benchmark with a cheap scalar
answer has that answer recorded from both systems and compared — a timing row is
only meaningful once the two agree. Every value in the tables below agreed.

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

| Benchmark | Mathilda | Mathematica 14.0 | |
|---|---:|---:|---|
| `Sort` | **314 ms** | 731 ms | 2.33× faster |
| `Riffle[x, 0.]` | **127 ms** | 230 ms | 1.80× faster |
| `PadRight`, Real fill | **18.7 ms** | 24.1 ms | 1.29× faster |
| STREAM triad, `x + 3. y` | **32.3 ms** | 35.2 ms | 1.09× faster |
| `Sin` (elementwise) | 13.2 ms | 12.9 ms | 1.03× |
| `Exp` (elementwise) | 13.4 ms | 12.5 ms | 1.07× |
| `Total` (reduction) | 2.62 ms | 2.34 ms | 1.12× |
| `Partition[x, 2]` | 17.8 ms | 12.6 ms | 1.42× |
| `Join` (two 10⁷ vectors) | 95.9 ms | 65.9 ms | 1.46× |
| `Reverse` | 33.8 ms | 15.1 ms | 2.23× |
| `Dot` (inner product) | **5.84 ms** | 5.08 ms | 1.15× |
| `PadRight`, default (exact `0`) fill | 1.19 s | 483 ms | 2.46× |
| `RotateLeft` | 37.1 ms | 12.9 ms | 2.88× |
| `Differences` | 56.6 ms | 15.1 ms | 3.75× |
| `Accumulate` (prefix scan) | 63.8 ms | 16.3 ms | 3.92× |

Mathilda leads on `Sort` (an 8-pass LSD radix sort on the double bit patterns,
against Mathematica's comparison sort) and is at parity on the threaded
elementwise kernels and reductions. The remaining 2–4× band —
`Reverse`, `RotateLeft`, `Differences`, `Accumulate`, `Dot` — is *serial* buffer
work against Mathematica's vectorised equivalents. These are memory-bound; the
fix is SIMD and threading in `ndstruct.c`/`ndreduce.c`, not a better algorithm.

`PadRight` with its **default** fill is the odd row, and it is a semantic fact
rather than an inefficiency: `PadRight[{1., 2., 3.}, 5]` is
`{1., 2., 3., 0, 0}` — exact Integer zeros beside Reals — in *both* systems, so
neither can hold the result in a uniform buffer and both fall back to a boxed
list. Given an explicit Real fill (`PadRight[x, n, 0.]`) Mathilda is the faster
of the two, by 1.29×.

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
this exercise — see [§8](#8-what-this-exercise-changed).

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

| Benchmark | Mathilda | Mathematica 14.0 | |
|---|---:|---:|---|
| Lennard-Jones energy, 1452 bodies, all pairs | **150 ms** | 326 ms | 2.17× faster |
| Mandelbrot, 800×800, 100 iterations | **785 ms** | 1.67 s | 2.13× faster |
| Logistic map, 10⁷ iterations | **178 ms** | 257 ms | 1.44× faster |
| Monte Carlo π, 10⁷ samples (vectorised) | **2.73 s** | 275 ms | 9.94× |

**Mathilda's compiled scalar code is faster than Wolfram's on all three compiled
kernels**, by 1.4–2.2×. This is the same result
[`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) reports for a
stencil (1.8–2.1×), now confirmed on three unrelated bodies: a nested
all-pairs `For` loop with indexed array reads (Lennard-Jones), an escape-time
loop with a compound `&&` guard (Mandelbrot), and a bare arithmetic recurrence
(logistic map).

Monte Carlo π was the `UnitStep` gap of §5 and is now **RNG-bound**: after the
narrowing kernels it is 2.73 s, of which `RandomReal[{0,1}, 10^7]` — called
twice — accounts for 2.4 s. The arithmetic and reduction together are ~130 ms.
Mathematica runs the whole benchmark in 275 ms, so its generator is roughly 10×
faster than `src/random.c`.

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

## 8. What this exercise changed

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

## 9. Summary

Of 38 benchmarks, Mathilda is faster on 9, within 1.5× on a further 12, and more
than 1.5× behind on 17.

**Where Mathilda is ahead** it is consistently in *compiled* code — all three
compiled scalar kernels, both compiled integer searches — plus `Sort` and the
structural buffer ops that were rewritten for this document. The compiled-code
advantage of 1.4–2.2× matches what `COMPILE_EXAMPLE.md` measures independently
on a stencil.

**Where Mathilda is behind**, the causes are specific and mostly named above,
not diffuse:

1. **No narrowing (float→int64) kernel category** — `UnitStep` and the
   `Floor`/`Ceiling`/`Round`/`IntegerPart`/`Sign`/`Quotient` family. Costs Game of
   Life 658× and Monte Carlo π 27×.
2. **`QRDecomposition`, `Eigenvalues`, `SVD` bypass LAPACK** — 6–19×.
3. **Serial buffer ops** — `Reverse`, `RotateLeft`, `Differences`, `Accumulate`,
   `Dot` sit in a 2–4× band against vectorised equivalents.
4. **Pattern-match dispatch** — 2.74× on naive recursive `fib`.

`Fourier` was the largest gap and is now closed; see §3.

Arbitrary-precision arithmetic is *not* on this list: bignum multiply, factorial
and power are all at parity, because both systems are calling GMP.

None of these is a design problem; each is a specific unfilled path, and the
ratios say which order to fill them in.

The methodological point worth keeping: the largest wins in §8 (658×, 237×,
192×, 158×, 135×, 118×, 52×) were all cases where an operation *worked correctly
and quietly took the slow path*. Nothing failed, no test broke, and no amount of
reading the code would have found them. They were found by measuring one
operation at a time against the same operation in another system, and by asking
`CompileDiagnostics` whether a body it accepted was actually compiled.

---

## 10. Reproducing

```bash
make -j$(sysctl -n hw.ncpu)
python3 comparisons/hpc_bench.py                    # full run, both systems
python3 comparisons/hpc_bench.py --only fft,sort    # a subset, by id
python3 comparisons/hpc_bench.py --scale 0.05       # smaller sizes, smoke test
python3 comparisons/hpc_bench.py --system mathilda  # one system only
```

The harness prints the markdown table above, writes raw JSON with `--json`, and
exits nonzero if any benchmark's recorded answer differs between the two systems.
Set `WOLFRAMSCRIPT` if `wolframscript` is not at the default macOS path.

## See also

- [`plans/HPC_IMPROVEMENT_PLAN.md`](../../plans/HPC_IMPROVEMENT_PLAN.md) — the
  plan to close every gap in §9, sequenced by value and risk.
- [`compile.md`](compile.md) — the compiler's design and milestones.
- [`COMPILE_EXAMPLE.md`](../compile_example/COMPILE_EXAMPLE.md) — one HPC problem
  end to end, with the bytecode.
- [`packed_arrays.md`](packed_arrays.md) — the automatic packing model that most
  of §8 is about.
- [`NDARRAY_REDUCTIONS_COMPARISON.md`](../../comparisons/NDARRAY_REDUCTIONS_COMPARISON.md)
  — Mathilda vs NumPy vs Mathematica on reductions and elementwise kernels.

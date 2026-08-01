# HPC improvement plan — reaching parity with Mathematica

Status: **proposed**. Target: no benchmark in
[`docs/design/performance.md`](../docs/design/performance.md) more than 1.2×
slower than Wolfram Language 14.0 on the same machine.

Baseline measured 2026-07-30, Intel Core i9-9880H (8 cores, AVX2), macOS 15.7.4.
Every number below is a measurement, not an estimate, unless marked *projected*.

---

## 0. Where we are

38 benchmarks: **9 faster than Mathematica, 10 within 1.5×, 19 more than 1.5×
behind.** The 19 gaps, ratio-ordered:

| | benchmark | Mathilda | WL 14.0 | root cause | phase |
|---:|---|---:|---:|---|:--:|
| 658× | Game of Life, 256², 100 gen | 65.7 s | 99.8 ms | no narrowing kernel (`UnitStep`) | 2 |
| 80× | `Fourier`, 2²⁰ complex | 1.54 s | 19.2 ms | never reaches FFTW | 1 |
| 40× | `Fourier`, 2²⁰ real | 855 ms | 21.2 ms | never reaches FFTW | 1 |
| 27× | Monte Carlo π, 10⁷ | 6.75 s | 253 ms | no narrowing kernel (`UnitStep`) | 2 |
| 18.7× | `QRDecomposition`, 500² | 83.0 ms | 4.43 ms | in-house, LAPACK bound but unused | 1 |
| 7.0× | `Eigenvalues`, 300² sym | 21.5 ms | 3.08 ms | in-house QR iteration | 1 |
| 6.5× | `SingularValueDecomposition`, 300² | 54.6 ms | 8.40 ms | in-house | 1 |
| 3.9× | `Accumulate` (prefix scan) | 63.8 ms | 16.3 ms | serial, scalar accessors | 4 |
| 3.8× | `Differences` | 56.6 ms | 15.1 ms | serial, scalar accessors | 4 |
| 3.5× | `PrimePi[10^9]` | 14.7 ms | 4.26 ms | algorithmic | 6 |
| 2.9× | `RotateLeft` | 37.1 ms | 12.9 ms | one redundant pass; serial | 4 |
| 2.9× | π to 100k digits | 36.6 ms | 12.8 ms | MPFR vs Chudnovsky | 6 |
| 2.7× | naive recursive `fib(25)` | 368 ms | 134 ms | rule-dispatch cost | 5 |
| 2.7× | `Inverse`, 500² | 8.18 ms | 2.99 ms | LAPACK boundary marshalling | 3 |
| 2.5× | `LinearSolve`, 1000² | 17.5 ms | 6.89 ms | LAPACK boundary marshalling | 3 |
| 2.5× | `PadRight`, default fill | 1.19 s | 483 ms | boxed-list construction | 5 |
| 2.4× | `Dot` (inner product) | 12.3 ms | 5.19 ms | rank-1 misses `cblas_ddot` | 1 |
| 2.2× | `Reverse` | 33.8 ms | 15.1 ms | serial | 4 |
| 1.6× | `Det`, 500² | 2.99 ms | 1.82 ms | LAPACK boundary marshalling | 3 |

**The shape of this list is the plan's main finding.** Only three items are
algorithmic (`PrimePi`, π digits, `fib`). Everything else is a path that exists
and is not reached, a marshalling cost around a path that *is* reached, or a loop
that has not been threaded. That is why a parity target is realistic.

---

## 1. Non-negotiables

Every item below is subject to these. They are not aspirations; two of them have
already cost this codebase a crash and a wrong answer.

**N1 — Never trade an exact answer for a fast one.** `Total[Range[10^6]]` is the
Integer `500000500000`, `Mean[Range[10]]` is `11/2`, and `Norm[{{1,2},{4,5}}]` is
`Sqrt[23 + 2 Sqrt[130]]`. A fast path that cannot produce the exact answer
declines and the ordinary path runs. This is `pack.h`'s contract and
`compile_internal.h`'s, and it is why `UnitStep` cannot simply be given a float
kernel (§3.2).

**N2 — Representation is never observable.** Packing changes storage and nothing
else: not a value, not an element's head, not a printed form, not an ordering.
Every item that touches a packed path ships with a differential sweep — the same
expression over a buffer and over the identical plain list, requiring
byte-identical printed output.

**N3 — A missing opt-in must not be silently slow.** The 2026-07-30 round found
four subsystems (linear algebra, the `Nest` family, the structural family, user
DownValues) that were correct, quiet, and 30–658× off. Phase 0 adds the audit
that would have caught them.

**N4 — Parity is a gate, not a report.** Every phase lands with its benchmark
row in a checked-in gate that fails the build on regression, on the
`bench_eval.c` normalized-baseline pattern.

---

## Phase 0 — the audit that makes the rest visible  ✅ DONE (2026-07-30)

*Estimated ~1 day; actual ~1 day. Prerequisite for everything.*

**It found a live bug on its first armed run**, which is the argument for having
done it first. `Differences` on an int64 buffer measured **89.3 ms** against
0.83 ms for the same-size float64 buffer. The int64 arm tested
`!ci_sub_i64(...)`, but `ci_sub_i64` is `__builtin_sub_overflow` and returns
**true on overflow** — so the loop abandoned on the first *successful*
subtraction and every call fell back to `delist_repack`. Correct answer, correct
dtype, correct element heads, **1193× too slow**, and invisible to every value
test in the suite. Fixed; the row is now 45 µs.

That is precisely the failure mode Phase 0 exists to make visible, and it was
introduced by the very work that motivated this plan.

**0.1 `tests/bench_pack.c`** — done. 21 workloads on the `bench_eval.c`
normalized-baseline pattern, wired into ctest as `bench_pack`. Every row is an
operation whose slow path is at least an order of magnitude worse, so a
regression cannot hide inside the 2.5× threshold.

Two design points worth keeping: the calibration is a list of exact **Rationals**
(a Rational leaf can never live in a machine buffer, so `pack_sniff` rejects the
first element and the calibration is immune to the thing being measured — the
trap `bench_eval.c` documents), and `main()` **asserts its inputs are actually
packed** before timing anything. Without that second check the file could pass
while measuring the unpacked path and silently re-baseline against it.

**0.2 An allowlist-vs-dispatch audit** — done, as
`tools/check_packed_aware.py` + `make check-packed-aware`, and wired into CI
beside `check-c99`. It scans 630 registered builtins, finds the 42 with an
NDArray dispatch and the 80 with an ND kernel, and diffs against `pack.c`'s
`AWARE` / `NOT_AWARE` / `INT64_OK` lists. Three checks: a head with a fast path
missing from `AWARE`; an `INT64_OK` entry that is not also in `AWARE` (the
narrower claim without the broader one does nothing); and a `NOT_AWARE` entry
with no kernel to clear (a note, not a failure — clearing an unset bit is
harmless but suggests drift).

Deliberate absences live in an `EXEMPT` table **with a reason each**, so a reader
can tell "considered and rejected" from "never noticed".

Verified non-vacuous: removing `Det`, `Inverse`, `LinearSolve` and `Differences`
from `AWARE` makes it name exactly those four, with the file and the marker that
identified each.

**0.3 A "did it take the fast path?" probe** — done. `MATHILDA_PACK_DIAG=gate`
now records every materialisation the transparency gate performs and prints an
exit summary ordered by elements. It catches the half of the problem no static
check can see: a head with **no fast path at all** that real workloads keep
handing buffers to. `Fourier` is exactly that case — 40–80× behind and invisible
to 0.2, because there is no dispatch site to find.

First run over a small corpus named `Fourier`, `UnitStep` and `Floor` — the
Phase 1.1 and Phase 2 targets — alongside `Cases` and `Position`, which are
correct by omission (they must descend the structure).

---

## Phase 1 — routing: paths that exist but are not reached

*Cost: ~1 week. Highest value per unit of risk in the plan: no new numerics, only
wiring, and each item has an existing correct implementation to differential-test
against.*

### 1.1 `Fourier` → FFTW (40× and 80×)  ✅ DONE (2026-07-30)

**Result: 40.4× → 1.15× and 80.2× → 1.30×** (855 ms → 24.5 ms, 1.54 s → 24.7 ms).
FFTW's own floor for this transform is 17.6 ms, so ~7 ms of marshalling remains.

What the plan got right, and what measurement corrected:

- **Right:** the head was not packed-aware, so the gate materialised the buffer
  before `fourier.c`'s existing NDArray fast path could see it. 986 ms against
  62.8 ms for the identical value written `NDArray[...]`.
- **Wrong:** "cache plans". Measured, `FFTW_ESTIMATE` planning for 2²⁰ costs
  **0.15 ms** against a 17.6 ms execute — not worth a cache, and a cache would
  have added an alignment hazard for nothing.
- **Missed entirely, and it was the biggest item after the gate:** the `b`-gather
  ran an integer division *and* a modulo **per element per axis**, when `b = 1`
  (the default) makes that gather the identity. Profiling put
  `machine_transform_buf` above every FFTW symbol combined. Scaling in place
  took 57 ms → 30 ms.
- Also: planning in place instead of through an `fftw_malloc`'d bounce buffer
  (30 ms → 22 ms), and a `hypot` per element in the real-collapse scan replaced
  by `max(|re|,|im|)` — the test is against a tolerance of 16·N·eps, where a
  factor of √2 is not a distinction (65 ms → 57 ms).

The lesson for the remaining phases: **profile before optimising.** Two of the
four wins here were not in the plan, and one item that *was* in the plan was
worth 0.15 ms.

FFTW is linked and used (`src/fourier.c:281–292`), but `fourier.c` has **no
`is_ndarray` fast path at all** and `Fourier` is not in the packed-aware
allowlist. So a packed 2²⁰ buffer is materialised into 2²⁰ `Expr` nodes by the
gate, read back one at a time into an `fftw_complex` buffer, transformed, and
rebuilt as 2²⁰ `Expr`s. 855 ms for 2²⁰ points is ~815 ns/point: that is the
marshalling, not the FFT.

- Mark `Fourier`/`InverseFourier` packed-aware; add a rank-N `is_ndarray` entry
  that hands `data` straight to FFTW and returns a packed `NDT_COMPLEX64` result.
- Real input should use `fftw_plan_dft_r2c` rather than widening to complex first.
- **Cache plans.** `fftw_plan_dft` is currently created and destroyed on every
  call. Key a small plan cache on (rank, dims, sign, kind); FFTW plans are
  reusable across buffers with `fftw_execute_dft`.
- Same treatment for `Fourier`'s DCT/DST siblings at `:988`.

*Projected: 855 ms → ~25 ms; 1.54 s → ~22 ms.* Risk: low. `Fourier`'s existing
numeric path stays as the fallback and is the differential reference.

### 1.2 `Dot` rank-1 → `cblas_ddot` (2.37×)  ✅ DONE (2026-07-30)

**Result: 2.37× → 1.15×** (12.3 ms → 5.84 ms at n = 10^7).

`dot2` tried `nd_blas_matmul` (rank 2 only) then fell to `ndarray_dot2`'s scalar
loop; `cblas_ddot` was bound but reached only by the LAPACK availability probe.
The BLAS path now covers every shape BLAS has: `dgemm`, `dgemv` both ways, and
`ddot`. Nothing needed marshalling — an NDArray is already the contiguous
row-major float64 buffer BLAS wants.

The plan's contingency turned out to be needed: Accelerate's `ddot` is
single-threaded, and one core reached 20 GB/s against Mathematica's 31 GB/s on a
160 MB working set. Chunking it across `nd_parallel_reduce` — the same machinery
`Total` already uses, for the same reason — took 7.86 ms to 5.84 ms. Summation
order therefore reassociates for large n, as it already does for `Total`;
`nd_parallel_reduce` runs one serial chunk below its threading threshold, so
small vectors stay bit-identical.

`dot2` tries `nd_blas_matmul` (rank-2 only) and then falls to `ndarray_dot2`,
Mathilda's own scalar loop. `cblas_ddot` is bound (`blas_bridge.c:73`) but reached
only by the LAPACK probe. Route rank-1 × rank-1 float64 to it, and rank-2 × rank-1
to `cblas_dgemv` (also already bound, `:216`).

*Projected: 12.3 ms → ~6 ms.* Note the target is bandwidth-bound: 10⁷ doubles is
160 MB of traffic, so 5.19 ms is 31 GB/s and near this machine's ceiling. If
Accelerate's `ddot` is single-threaded, this needs `nd_parallel_reduce` over
chunked `cblas_ddot` calls instead.

### 1.3 `QRDecomposition` → LAPACK (18.7×)  ⚠️ PARTIAL (2026-07-30)

**Result: 18.7× → 14.0×** (83.0 ms → 58.4 ms). Correct, verified bit-identical
on 240 matrices — and a long way short of the ~6 ms this plan projected.

The routing itself works: `mat_qr_mathilda` (dgeqrf + dorgqr) is wired in ahead
of the delist, so a packed matrix goes straight to LAPACK. Three gates were
needed and each was found by a test rather than by reading:

- **Inexact matrices only.** `QRDecomposition[{{1,2},{3,4}}]` is exact
  (`1/Sqrt[10]`), and `na_load_matrix` would have flattened it to `-0.316228`,
  with the opposite sign too. Same gate `builtin_norm` needed.
- **Real matrices only.** Mathilda's complex convention is
  `m == ConjugateTranspose[q].r`, so `q` is `conj(Q^T)`, not the plain transpose.
- **Comfortably full rank only.** The in-house Gram-Schmidt truncates its answer
  to the numerical rank (`q` is `rank x n`), which the test suite asserts and
  which neither LAPACK nor Mathematica do. That decision comes out of the
  orthogonalisation residual and cannot be reproduced from `R`, so near-singular
  input falls back.

**Why it stalled, and what it means for 1.4 and 1.5.** LAPACK's own dgeqrf +
dorgqr on a 500x500 is a few milliseconds. The rest is marshalling:

  - `na_load_matrix` converts to column-major **element by element**, through an
    out-of-line `ndt_get` per element with a cache-hostile scattered write —
    this is exactly Phase 3.
  - `na_build_matrix` builds each output the same way.
  - `{q, r}` is a **List of two arrays**, and the no-nesting invariant forbids a
    packed array inside a plain List, so both outputs are materialised: 500k
    `Expr` nodes, measured at ~24 ms of the 58 ms. Mathematica has no such
    restriction and returns packed arrays inside the list.

This plan's own sequencing note said "Phase 3 should land before or alongside
Phase 1.3–1.5 so the new LAPACK paths are not built on the slow boundary."
That was right and I did not follow it. **Phase 3 should be done next**, before
1.4 (SVD) and 1.5 (Eigenvalues), which return two and three arrays respectively
and will hit the same two walls harder.

The output-materialisation half is not a marshalling bug but a structural
consequence of the no-nesting invariant, and deserves its own decision rather
than being absorbed into Phase 3.

The largest linear-algebra gap and the one with no excuse: `mat_lapack_dgeqrf` is
already bound and used elsewhere (`lapack_bridge.c:238`). `builtin_qrdecomposition`
routes an NDArray to `linalg_delist_and_reeval`, i.e. straight back to the
symbolic path.

Add `ndla_qrdecomposition`: `dgeqrf` + `dorgqr`, returning `{Q, R}`. Note that
`{Q, R}` is a **List of arrays**, so both must be materialised — a buffer must
never sit inside a plain `List` (the no-nesting invariant). That costs one
conversion and is still ~15× better than today.

*Projected: 83.0 ms → ~6 ms.*

### 1.4 `SingularValueDecomposition` → LAPACK (6.5×)

`mat_lapack_dgesdd` is bound and already used by `ndla_matrix_norm_direct` for the
spectral norm. Extend to the full `{U, Σ, V}` form with `jobz='A'`.

*Projected: 54.6 ms → ~10 ms.* Risk: **sign convention.** LAPACK's singular
vectors are determined only up to sign; the existing in-house SVD fixes signs some
way, and tests pin it. Normalise (e.g. force the largest-magnitude component of
each column positive) and differential-test against the current output before
switching the default.

### 1.5 `Eigenvalues` / `Eigenvectors` → LAPACK (7.0×)

Deliberately deferred once before, for a stated reason:
[`NDARRAY_REDUCTIONS_COMPARISON.md`](../comparisons/NDARRAY_REDUCTIONS_COMPARISON.md)
records that the eigenvalue **ordering** convention (by |λ| descending, ties broken
by position) cannot be reproduced from LAPACK output without risking parity.

That is a reason to build an ordering shim, not to skip the item:

- Symmetric/Hermitian input → `dsyevd`/`zheevd`; general → `dgeev`/`zgeev`.
- Sort the returned pairs into Mathilda's convention, carrying eigenvectors with
  their values. Ties broken by original position require a **stable** sort keyed
  on the LAPACK output index.
- Differential-test the whole `eigen_tests` corpus against the in-house path, and
  keep the in-house path for any case where the shim cannot reproduce the order
  (degenerate clusters especially).

*Projected: 21.5 ms → ~4 ms.* Risk: medium — the highest-risk item in Phase 1,
which is why it is last in it.

---

## Phase 2 — the narrowing kernel category  ◐ DONE, targets still blocked (2026-07-30)

**The category works.** `NDUnaryKernel` gained `to_int_r` / `to_int_i` / `to_int`
(additive, so the Compile VM and every other kernel consumer are untouched), and
`ndarray_map_unary` writes an `NDT_INT64` result for them. `UnitStep` — which had
no kernel at all — is registered, and `Floor`, `Ceiling`, `Round`, `IntegerPart`
and `Sign` came **off** `NOT_AWARE`.

| 10⁷ elements | before | after | |
|---|---:|---:|---|
| `UnitStep`, real | 5004 ms | 44 ms | **113×** |
| `UnitStep`, int64 | 5004 ms | 27 ms | **185×** |
| `Floor` | materialised | 37 ms | |
| `Sign` | materialised | 76 ms | |
| `Round` | materialised | 101 ms | |

Both arms were needed, not just the narrowing one: Monte Carlo π's `UnitStep`
sees float64, but Game of Life's sees **int64** (its neighbour count is an
integer grid). Verified by a 66-case differential sweep including past-int64
values, infinities, complex, Rationals and MPFR — a real input that will not fit
an int64 abandons the whole array so the List path answers with a bignum.

**Neither headline benchmark reached parity, because neither was actually
`UnitStep`-bound once measured.**

| | before | after | vs WL |
|---|---:|---:|---:|
| Monte Carlo π, 10⁷ | 6.75 s | **2.73 s** | 9.94× (was 26.7×) |
| Game of Life, 256², 100 gen | 65.7 s | 67.9 s | 722× (unchanged) |

Breaking each down:

- **Monte Carlo is now RNG-bound.** `RandomReal[{0,1}, 10^7]` alone is **1.22 s**,
  and the benchmark calls it twice — 2.4 s of the 2.73 s. Everything else totals
  ~0.13 s (`u^2` 13 ms, the subtraction 84 ms, `UnitStep` 27 ms, `Total` 7 ms).
  Mathematica does the *whole* benchmark in 275 ms, so its RNG produces 2×10⁷
  doubles in a fraction of that. **New item: `random.c` throughput**, roughly 10×.

- **Game of Life is `Sum`-bound, and spectacularly.**
  `Sum[RotateLeft[m,{i,j}], {i,-1,1}, {j,-1,1}]` on a packed 256² integer grid
  takes **463 ms** and returns an **unpacked** result, where writing the same nine
  terms as an explicit `Plus` takes **0.81 ms** — a **572× gap** for identical
  arithmetic. A single `RotateLeft` is 57 µs, so the nine shifts plus eight adds
  should be under 1 ms. `Sum` over an array-valued body is not using the buffer
  path at all. That, not `UnitStep`, is the whole of the remaining 65 s.
  **New item: `Sum`/`Product` accumulation over array-valued bodies.**

  **Done 2026-07-30 — 65.7 s → 260 ms (253×), now 2.90×.** And the first thing it
  found was not a performance bug at all: `RotateLeft[list, i]` with a symbolic
  `i` returned the list UNROTATED (`rotate_rec` defaulted its amount to 0 and
  never rejected a non-integer spec), so `Sum`'s closed-form stage saw a body
  with no `i` dependence and returned `9 m`.
  `Sum[RotateLeft[{1,2,3},i],{i,0,2}]` gave `{3,6,9}` against Mathematica's
  `{6,6,6}`. **This benchmark had been measuring a wrong computation from the
  day it was written**, which no timing could have revealed — only an answer
  check. Mathilda and Mathematica now agree on a 40×40 evolution.

  The other two: `Sum` no longer attempts its closed form for a short range over
  an array-valued body (that attempt costs time proportional to the BODY, so a
  one-term Sum over a 256² grid cost 197 ms and then fell through anyway), and
  the DownValue exemption from §"Defining a helper function cost 100×" was
  extended to int64 — it had covered only float64, and Life's grid is integer, so
  every helper call materialised.

Both were invisible before this phase: `UnitStep` at 500 ns/element was large
enough to hide them.

### Original analysis

*Cost: ~3 days. Fixes the single largest gap (658×) and a 27× one, and unblocks
five kernels that are currently forbidden from the fast path.*

**The problem.** `UnitStep[x]` answers with an **exact Integer** (`UnitStep[{-1., 1.}]`
is `{0, 1}`, not `{0., 1.}`). Mathilda's kernel categories are "real-closed"
(float64 → float64) and "escaping" (real → complex). There is no float64 → int64
category, so `UnitStep` has **no kernel at all** and costs ~500 ns/element — a
measured 5.0 s for one `UnitStep` over 10⁷ elements.

That single omission is Game of Life's 658× and Monte Carlo π's 27×, because both
select with `UnitStep` over an array.

The same gap is why `Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign` and `Im`
sit on `pack.c`'s `NOT_AWARE` list: they *do* have real-closed kernels, which keep
the float64 dtype and write `1.0` where the list gives the exact `1` — a wrong
element head, so they must materialise.

**The work.**

- Add a `UK_NARROW(NAME, REXPR)` macro beside `UK_CLOSED`/`UK_ESC` in
  `ndkernels.c`: real input, `NDT_INT64` output buffer, with an overflow check per
  element that abandons the whole result (the `ci_*_i64` contract) so a value
  outside int64 falls to the exact path rather than wrapping.
- Register `UnitStep`, `Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign`, and
  `Boole`.
- The dispatcher must pick the narrowing kernel **only for a real input**:
  `Floor` of an exact Integer list is the identity and must stay exact;
  `Floor[3/2]` is `1`, a Rational input the buffer cannot hold.
- Remove the six heads from `NOT_AWARE`, which is what makes the win reach the
  packed surface.
- `Clip` needs its own treatment (its clamp bounds may themselves be exact) and
  stays on `NOT_AWARE` for now.

*Projected: Game of Life 65.7 s → ~150 ms; Monte Carlo π 6.75 s → ~300 ms; plus
`Floor`/`Round`/`Sign` over a packed list going from "materialise" to "one pass".*

Verification: the existing 134+166-case head/kernel differential sweep that
originally *put* those six heads on `NOT_AWARE` is the exact test that must now
pass with them off it.

---

## Phase 3 — the LAPACK boundary  ◐ PARTIAL (2026-07-30)

**Done: the element-wise marshalling.** `na_load_matrix` / `na_build_matrix`
converted between row-major and LAPACK's column-major **element by element**,
through an out-of-line `ndt_get` per element with a cache-hostile scattered
write. For the float64 real case — which is every measured path — that is not
marshalling at all, only a layout change: same-layout is a `memcpy`, and
row-major to column-major IS a transpose, now 32x32 cache-blocked like
`ndstruct_transpose`.

| | before | after | vs WL |
|---|---:|---:|---:|
| `Det`, 500² | 2.88 ms | **1.45 ms** | **1.09× FASTER** |
| `Inverse`, 500² | 9.29 ms | 7.37 ms | 2.43× (was 3.02×) |
| `LinearSolve`, 1000² | 17.0 ms | 15.7 ms | 2.44× (was 2.73×) |

`Det` now beats Mathematica. `QRDecomposition` did not move (59 ms): its cost is
`dorgqr` plus the output materialisation, not the load.

**Still open, and re-profiling changed what the next item should be.** After the
fix, `LinearSolve` at 1000² samples as roughly:

    BLAS/LAPACK       ~2400 samples   (63%)
    mach_vm_map / mach_vm_deallocate / madvise
                      ~1034 samples   (~27%)
    na_load_matrix     ~376 samples   (10%)

A quarter of the time is the kernel mapping and unmapping memory.

**A scratch pool was built for this, measured, and REJECTED (2026-07-30).** The
reasoning looked sound — a 1000×1000 double matrix is 8 MB, macOS serves
allocations that large straight from `mmap` and returns them with
`munmap`/`madvise`, and `LinearSolve` discards and re-requests exactly that size
every call. A size-keyed pool of parked blocks was implemented in `numarray.c`
(header-prefixed sizes, mutex-guarded, opt-in via `na_load_matrix_scratch` so the
~116 plain `free()` calls elsewhere could not turn into heap bugs) and the three
measured paths converted.

It made no difference at all. Min-of-five, two runs each:

| | with pool | without |
|---|---|---|
| `LinearSolve`, 1000² | 15.45, 15.42 ms | 15.56, 15.12 ms |
| `Inverse`, 500² | 8.24, 8.94 ms | 8.49, 8.02 ms |
| `Det`, 500² | 1.75, 1.54 ms | 1.44, 1.63 ms |

**Why it was the wrong target, and the lesson.** The flat "top of stack" profile
attributes samples to a symbol, not to a culprit. Reading the CALL GRAPH instead
shows those allocations belong to Accelerate:

    30 ??? (in libBLAS.dylib)
     : 29 _szone_free (in libsystem_malloc.dylib)
     : | 13 free_large
     : | + 13 madvise

They are libBLAS's own per-call internal workspace, inside its threaded
implementation — nothing Mathilda allocates, and nothing a pool on our side can
reach. Our 8 MB load buffer is a single `malloc`/`free` pair per call and was
never the problem.

*A flat profile says where the time goes; only the call graph says whose it is.*
The pool was reverted rather than kept "in case it helps later": it is ~150 lines
carrying a class of heap bug (a pooled pointer reaching `free`) for no measured
benefit.

That leaves the transpose-avoidance bullet below as the remaining item here.

The plan's third bullet — **avoid the transpose entirely** — is still the
structurally right answer and is untouched. Solving with a row-major `A` is
solving `Aᵀ`, so `dgetrf` on the row-major buffer followed by `dgetrs` with
`trans='T'` gives `A x = b` with no transposition at all. It needs `dgetrs`
bound, which it is not yet.

---

### Original analysis

*Cost: ~2 days. Fixes three rows at once and everything Phase 1 adds inherits it.*

**The evidence.** `A . B` — row-major, no transposing load — is at **parity**
(8.04 vs 6.71 ms) despite being the largest flop count in the group. Every
operation that goes through a *transposing* load is 1.6–2.7× off: `Det` 1.64×,
`LinearSolve` 2.54×, `Inverse` 2.74×.

**The cause.** `na_load_matrix` (`src/linalg/numarray.c`) converts an NDArray to
LAPACK's column-major layout **element by element**, through a `ndt_get` call per
element and a scattered write:

```c
for (int i = 0; i < r; i++)
    for (int j = 0; j < c; j++) {
        double re, im;
        ndt_get(e->data.ndarray.data, i * c + j, dt, &re, &im);
        na_put(out, na_index(i, j, r, c, colmajor), re, im, want_complex);
    }
```

For a 1000×1000 that is 10⁶ function calls with a cache-hostile write pattern.
`ndstruct_transpose` already does a 32×32 cache-blocked transpose of the same data
in 1.21 ms; this loop is several times that, twice per call (in and out).

**The work.**

- Specialise the common case: float64 source, float64 destination, no complex
  marshalling → reuse the cache-blocked transpose kernel from `ndstruct.c` rather
  than the generic accessor loop.
- Where the source dtype already matches and no transpose is needed
  (`colmajor == false`, float64), `memcpy` the whole buffer.
- Better: **avoid the transpose entirely.** LAPACKE exposes row-major variants
  (`LAPACK_ROW_MAJOR`), and for a symmetric operation like `dgesv` a row-major
  matrix is the transpose of the column-major one — solving `Aᵀx = b` needs only a
  flag. Evaluate this before writing a faster transpose; it may delete the work.
- Same treatment for `na_build_matrix` on the way out.

*Projected: `LinearSolve` 17.5 → ~8 ms, `Inverse` 8.18 → ~3.5 ms, `Det` 2.99 →
~2 ms.* Risk: low. The loads are exercised by every LAPACK test already.

---

## Phase 4 — thread and vectorise the serial buffer ops

*Cost: ~3 days.*

`nd_parallel_for` / `nd_parallel_reduce` already exist and are used by the
reductions and elementwise kernels — which is exactly why `Total`, `Sin` and `Exp`
are at parity. The structural ops written in the 2026-07-30 round are correct but
serial and go through scalar accessors.

| op | now | traffic | achieved | target |
|---|---:|---|---:|---:|
| `Reverse`, 10⁷ | 33.8 ms | 160 MB | 4.7 GB/s | ~15 ms |
| `RotateLeft`, 10⁷ | 37.1 ms | 240 MB | 6.5 GB/s | ~13 ms |
| `Differences`, 10⁷ | 56.6 ms | 160 MB | 2.8 GB/s | ~15 ms |
| `Accumulate`, 10⁷ | 63.8 ms | 160 MB | 2.5 GB/s | ~18 ms |

- **`RotateLeft` has one redundant pass.** `nd_rotate_axes` copies the input into
  a scratch buffer and *then* rotates into a second — three touches of the data
  where two suffice. Rotate directly from the source on the first non-identity
  axis. Free ~30% before any threading.
- **`Differences`** currently reads through `ndt_get`/`ndt_set` per element. Take a
  typed `double*`/`int64_t*` pointer per dtype and thread the range. (The int64 arm
  must keep its `ci_sub_i64` overflow abandonment.)
- **`Accumulate`** is a prefix scan — sequential in the naive form, but the
  standard two-pass parallel scan (per-block partial sums, then offset each block)
  parallelises it exactly. The float64 arm must keep left-to-right summation order
  within a block so the answer does not change; the block offsets are added in
  order, so the result is deterministic but *not* bit-identical to the serial
  scan. **Decide this explicitly** — if bit-identity is required, thread only the
  int64 arm, where it is exact.
- **`Reverse`/`Join`/`Partition`** are pure block moves: thread the memcpy.

Risk: the `Accumulate` reassociation is the only semantic question in this phase,
and it is called out above rather than discovered later.

---

## Phase 5 — allocation and dispatch

*Cost: weeks, not days. The deepest items; lowest value per unit of risk.*

**5.1 `PadRight` with the default fill (2.46×).** Both systems produce a *boxed*
mixed exact/inexact list here — this is not a packing gap, it is `Expr`
construction cost measured directly: 10⁷ nodes in 1.19 s vs Wolfram's 483 ms.
Shares a cause with 5.2.

**5.2 Rule dispatch (2.74× on naive `fib`).** `fib[k_] := fib[k-1] + fib[k-2]` at
`k = 25` is 243k pattern matches, and Mathilda's is ~2.7× more expensive per call.
The `dispatch_arity` index in `symtab.h` already skips rules that cannot match;
the remaining cost is in `MatchEnv` allocation and the substitution walk. Profile
before designing: this is the one item in the plan where the cause is *not* yet
established, and it should not be attempted until it is.

Both are worth doing and neither should be started before Phases 1–4, because
they are the highest-risk, highest-effort, lowest-ratio items on the list.

---

## Phase 6 — number theory and precision

*Cost: ~1 week. Genuinely algorithmic; the only items where matching Wolfram means
implementing a better method rather than reaching an existing one.*

**6.1 π to 100k digits (2.86×).** Mathilda goes through MPFR's π. Wolfram is
~3× faster, which suggests a tuned Chudnovsky with binary splitting. Check first
whether MPFR's `mpfr_const_pi` is the bottleneck or the conversion to a Mathilda
Real is — the latter would be a Phase-3-style marshalling fix, not an algorithm.

**6.2 `PrimePi[10^9]` (3.45×).** Mathilda has Deléglise–Rivat and 14.7 ms is a
genuine computation. Wolfram's 4.26 ms is cold-measured in a fresh process, so it
is real too. Likely a tuning gap (sieve segment size, the μ/φ table cut-off)
rather than a different algorithm. Profile the phases of the existing
implementation before touching it.

**Explicitly not here: arbitrary-precision arithmetic.** Bignum multiply,
factorial and power are all at parity — both systems call GMP. An earlier draft of
`performance.md` reported a 24× multiply gap; it was a measurement artifact
(Wolfram does not materialise a discarded bignum product), and `mpz_mul` on the
same operands in a standalone C program takes 2.93 ms against Mathilda's 3.07 ms.
There is nothing to win here.

---

## Phase 7 — the second sweep's leftovers (added 2026-07-31)

*Five kernels were added to `hpc_bench.py` to probe subsystems the original 38
never touched — a Krylov solver, direct convolution, ODE integration, an
irregular hash-keyed reduction, interpolation. Three of the five found something
and were fixed on the spot (see `docs/design/performance.md` §8); what is below
is what those fixes did not finish.*

Results, after the fixes:

| | benchmark | Mathilda | WL 14.0 | status |
|---:|---|---:|---:|---|
| 2.60× faster | `Interpolation`, 10⁴ nodes / 10⁴ evals | 7.65 ms | 19.9 ms | ✅ was 74× behind |
| 1.52× faster | `NDSolve`, Lorenz to t = 200 | 31.5 ms | 47.8 ms | ✅ no work needed |
| 1.84× | Conjugate gradient, 256², 100 its | 155 ms | 84.0 ms | 7.3 |
| 7.21× | `Tally`, 10⁷ → 10⁴ bins | 129 ms | 17.9 ms | 7.2 (was 69×) |
| 9.27× | `ListConvolve`, 1024², 5×5 | 312 ms | 33.7 ms | 7.1 (was 12.8×) |

**7.1 `ListConvolve` / `ListCorrelate` are not packed-aware.** The engine choice
and the direct engine's inner loop are both fixed, and what is left is boxing:
10⁶ input `Expr` nodes materialised from the buffer, 10⁶ output nodes allocated
and freed, `expr_free` alone a sixth of the profile. At 13 ns per multiply-add
the arithmetic is no longer the cost. Needs the standard treatment — accept
`EXPR_NDARRAY` input, write an NDArray result, join `AWARE`/`INT64_OK`. This is
the single largest remaining ratio in the second sweep and the most mechanical.
*Projected: within 2× of Wolfram, from 9.27×.*

**7.2 `Tally`'s probe.** `ndred_tally` now hashes machine words and the table
grows to the distinct count, which took it from 69× to 7.2×. The remaining
16 ns/element is one random probe per element against a table that fits in L2.
Options, in order of appeal: a direct-mapped counter array when the value range
is small and dense (the common case — histogramming into bins), which removes
hashing entirely; storing the key inline in the slot to remove one dependent
load; a narrower slot type for a smaller footprint. Measure before choosing.
*Projected: 2–3×, from 7.21×.*

**7.3 The conjugate-gradient row (1.84×) is a composition cost, not a kernel
one.** Every operation in the loop already has a buffer path — `RotateLeft`,
elementwise `Plus`/`Times`, `Total` with a level spec — and the Jacobi row built
from the same pieces sits at 1.29×. The difference is that a Krylov iteration
alternates stencil work with *global reductions*, so it has more temporaries and
more scalar-to-array broadcasts per sweep. Profile before proposing anything;
this is a candidate for Phase 5's allocation work rather than a new kernel.

**7.4 `NDSolve` segfaults when a dependent variable already has a value.** Found
while running the second sweep, not caused by it — it reproduces on a clean
checkout:

```mathematica
x = RandomReal[{0, 1}, 100000];
NDSolve[{x'[t] == 10.(y[t] - x[t]), ...}, {x, y, z}, {t, 0, 5}]   (* SIGSEGV *)
```

Size-dependent: 10⁵ elements crashes, 10³ does not, and a scalar `x = 1.` does
not — so it looks like unbounded recursion over the materialised list rather than
a bad pointer. Wolfram solves the (different, and wrong) problem and returns.
This is a crash under default settings from ordinary user input — a symbol left
bound from earlier in a session — and should be fixed ahead of the performance
items above.

**7.5 Re-examine what else has no fast path** — *done, see Phase 8.* `Outer` was
the flagged item and is fixed; the sweep that found it found six more.

---

## Phase 8 — the third sweep's leftovers (added 2026-07-31)

*Eight pipeline-shaped kernels from real applications, run against Mathematica
**and NumPy** (`docs/design/performance.md` §9). Seven fixes came out of it; what
is below is what they did not finish. Items 7.1 and 7.4 above are unchanged and
remain the highest-value entries on this list.*

Results, after the fixes:

| | benchmark | Mathilda | WL 14.0 | NumPy | status |
|---:|---|---:|---:|---:|---|
| 15.9× faster | N-body all-pairs, 1024 bodies | 439 ms | 6.98 s | 396 ms | ✅ was 55.5 s |
| 3.24× faster | k-means, 100000×8, k=16 | 2.26 s | 7.33 s | 765 ms | ✅ was unbounded |
| 1.01× | Black–Scholes MC, 10⁷ paths | 565 ms | 349 ms | 560 ms | ✅ at NumPy parity |
| 1.26× | Logistic regression, 200000×32 | 2.71 s | 2.15 s | 839 ms | 8.3 |
| 1.43× | 3D heat equation, 128³ | 3.63 s | 2.55 s | 2.25 s | — |
| 1.66× | Welch PSD, 1024 blocks | 150 ms | 90.2 ms | 93.9 ms | — |
| 13.5× | Return series (EMA/vol/drawdown) | 2.02 s | 149 ms | 30.4 ms | 8.1 |
| 22.5× | Gaussian blur + Sobel, 1024² | 2.19 s | 97.0 ms | 113 ms | 7.1 |

**8.1 No vectorized scan.** `Accumulate` is the only one, and it only does
`Plus`. An exponential moving average is a general linear recurrence and a
maximum drawdown needs a running maximum; both are written `FoldList[f, x0,
rest]` with a pure function, which is one interpreted evaluation per element.
Mathilda 2.02 s, Mathematica 149 ms, SciPy `lfilter` + `np.maximum.accumulate`
30.4 ms — so this is the one row where Mathilda is far behind *both* other
systems, and the only structural gap the third sweep found that is not a missing
buffer path. Two tractable pieces: `FoldList[Max, …]` / `FoldList[Min, …]` as a
buffer scan (running extremum is common and trivially vectorizable), and a
first-order linear recurrence recogniser for `a #1 + b #2 &`. The general case
needs the `Compile[]` engine, which already exists — routing `FoldList` with a
compilable pure function through it is the principled answer.
*Projected: 5–10×, from 13.5× behind Wolfram.*

**8.2 The remaining ML composition costs.** Logistic regression is 3.22× behind
NumPy and k-means 2.95×, both after large wins, and neither now has a hot spot —
the time is spread across allocation and evaluation overhead on many
medium-sized array temporaries. Same category as 7.3 (conjugate gradient) and
the same recommendation: this is Phase 5 allocation work, not a new kernel, and
nothing should start before a profile attributes the cost to a structure.

**8.3 `Outer` covers only float64 and five heads.** `nd_outer2` handles `Plus`,
`Subtract`, `Times`, `Min`, `Max` on float64 operands and declines everything
else to the tree-building path — so `Outer[Times, intVector, intVector]` is still
~750 ns per element. Extending it to int64 needs the `ci_*_i64` overflow-abandon
contract; extending it to `Divide`/`Power` needs the non-finite cases to agree
with the interpreter's symbolic forms (`ComplexInfinity`, `Indeterminate`), which
is why they were left out rather than because they are hard.

**8.4 The gate materialises rather than packs in the remaining mixed cases.**
`pack_lift_listable_args` lifts a plain List to meet a buffer only when the
shapes are equal, or a prefix under `packed_broadcast_ok` (Plus and Times).
`Power` with a rank mismatch still threads, because
`ndarray_elementwise_power` requires `same_shape`. Giving Power the same
broadcast pre-pass `ndarray_elementwise` has would let it join the list.

**8.5 Packing a large plain List is itself the cost now.** `a + b` with one
plain 10⁶ operand went 418 ms → 50.5 ms, and ~26 ms of what remains is
`pack_force` reading a million boxed `Expr` nodes once. That is inherent to
meeting an unpacked value, but it means the fix pays best when the plain operand
is *small* — which is the case it was written for. Worth knowing before
optimising further: the answer is for the producer to pack, not for the consumer
to keep re-packing.

---

## 7. Sequencing

```
Phase 0  audit + gate            ──┐  (prerequisite: makes regressions visible)
Phase 1  routing to FFTW/LAPACK  ──┼─→ independent, parallelisable
Phase 2  narrowing kernels       ──┤
Phase 3  LAPACK marshalling      ──┘  (do before/with 1.3–1.5: they inherit it)
Phase 4  thread the buffer ops       (independent)
Phase 6  number theory               (independent)
Phase 5  allocation + dispatch       (last: deepest, least certain)
```

Phases 1–4 and 6 touch disjoint files and can proceed in parallel. Phase 3 should
land before or alongside Phase 1.3–1.5 so the new LAPACK paths are not built on
the slow boundary.

---

## 8. Projected end state

Applying the projections above to the 19 gaps:

| after | benchmarks >1.5× behind | worst remaining |
|---|---:|---|
| today | 19 | 658× |
| Phase 1 | 14 | 658× |
| Phase 1+2 | 12 | 3.9× |
| Phase 1+2+3 | 9 | 3.9× |
| Phase 1+2+3+4 | 5 | 3.5× |
| + Phase 6 | 3 | 2.7× |
| + Phase 5 | 0–2 | ~1.5× |

**Phases 1–4 alone take the worst gap from 658× to under 4× and cut the count of
>1.5× gaps from 19 to 5.** They are also the low-risk phases: routing to a library
that is already linked, a new kernel category, a marshalling fix, and threading
loops with machinery that already exists.

Parity ("no row worse than 1.2×") requires Phase 5, and Phase 5.2 in particular is
not yet understood well enough to promise. The honest target is: **everything
except rule-dispatch at or near parity, and rule-dispatch measured and scoped
before it is attempted.**

---

## 9. Explicit non-goals

- **Native codegen for `Compile[]`.** Not needed for this target — Mathilda's
  compiled code is *already* 1.4–2.2× faster than Wolfram's on every compiled
  benchmark, and `COMPILE_EXAMPLE.md` measures Wolfram's own native-C target as
  *slower* than its bytecode VM on a stencil.
- **Beating GMP.** §6 above.
- **GPU offload.** Out of scope for parity with a CPU-bound competitor.
- **Changing the exactness contract to win a benchmark.** N1. `PadRight`'s
  default-fill row stays as it is, because Mathematica produces the same mixed
  list; the only thing to fix there is allocation cost.

---

## 10. Verification

Per item:

1. **Differential sweep** — the operation over a packed buffer against the same
   value as a plain list, printed forms byte-identical (N2). The existing
   `ladiff`/`sdiff`-style sweeps are the template.
2. **Cross-system value agreement** — `hpc_bench.py` already fails if the two
   systems disagree on a benchmark's recorded answer.
3. **The full 395-binary suite**, with only the known-stale `simplify_tests`
   expectation failing.
4. **`bench_pack` / `bench_eval` / `bench_compile` / `bench_assoc` /
   `bench_ndarray_linalg`** all within gate, with the new row's baseline recorded.
5. **Valgrind** definitely-lost unchanged from the macOS start-up baseline.

And per phase: re-run `comparisons/hpc_bench.py` and update
`docs/design/performance.md`'s tables and its §9 counts. The document is the
scoreboard; it should never be stale relative to the code.

---

## Phase 9 — the fourth sweep's leftovers (added 2026-07-31)

Phase 8's items 8.1 (no vectorized scan), 8.3 (`Outer` coverage) and HPC 7.1
(`ListConvolve`/`ListCorrelate` packed-aware) are **DONE** — see
`docs/experiments/11-hpc-sweep-numpy-gap/README.md` and `performance.md` §10. What that
sweep left, in value order:

**9.1 — `Transpose` copies where NumPy views. (largest application gap)**
Logistic regression is 3.06× NumPy for one reason: `Transpose[X]` on a
200000×32 matrix costs 27 ms and the kernel does it 100 times, because it is
loop-invariant but re-evaluated. NumPy's `.T` is a view and costs nothing.

The clean fix is a strided/transposed NDArray view, which **every consumer must
honour** — that is a design change, not a fast path, and it interacts with the
no-nesting invariant. A cheaper 80% is to teach `ndarray_dot2` to accept a
transpose flag and have `Dot` recognise `Transpose[a] . b` *before* the argument
is evaluated, which needs `Dot` to hold its arguments — also not free. Decide
which before building either.

**9.2 — No interpreter-level fusion.** k-means is 2.92× NumPy: three passes over
6.4 MB where one would do. `Compile[]` fuses (experiment 2); ordinary array code
does not. This is the same 3× that shows on every unfused ML-shaped kernel.

**9.3 — `dgemm` slows the next threaded loop by 1.45×.** Reproducible to ±1%
(`performance.md` §10). Accelerate's worker threads persist and compete with
`nd_parallel_for`. Options: yield before our own parallel region, pin, cap our
thread count when Accelerate has been used, or share its pool. **Measure which,
do not guess** — the last thread-pool hypothesis in this plan (the linalg
scratch pool) was wrong and cost 150 lines.

**9.4 — 1-D convolution is 11× NumPy.** The per-output interior test is now the
cost. Splitting the output range into boundary/interior/boundary would let the
interior run as a flat loop with no test at all.

**9.5 — `Part` with a span still materialises a position array per axis.**
14.0 → 7.4 ms came from the block copy; the remaining 7.4 ms is
`build_axis_selector` writing 10⁶ int64 positions before the copy. A
`start/step/n` selector representation makes it a pure memcpy.

**9.6 — `Tally` at 7.17×** and **9.7 — `Reverse`/`RotateLeft` at 2.9–3.8×** are
unchanged from Phase 7.

**9.8 — `ArrayReshape` does not exist.** Reshaping a buffer is a metadata
change; it is currently unavailable at any speed.

**Still the highest priority overall: HPC 7.4, the pre-existing `NDSolve`
segfault.** It is a crash, not a slowdown, and it has now outlived two sweeps.

---

## Phase 10 — the fifth sweep's leftovers (added 2026-07-31)

Eight application domains, experiments 12–19; eleven fixes, tabulated in
`docs/design/performance.md` §12. Phase 9's items **9.1 (Transpose views)** and
**9.2 (no interpreter fusion)** are now the two largest named gaps in the whole
suite and are unchanged. What this sweep adds:

**10.1 — A mixed-dtype tuple return destroys every array in it.** `{a, b}` of two
600² float64 matrices costs 0.55 ms and packs into a rank-3 buffer; `{a, b, m}`
with an `int64` mask costs **53 ms** and does not pack, and the caller then pays
**107.8 ms against 0.90 ms — 120×** on its next operation. This is the third sweep's
return-statement defect one dtype away: `pack_sniff` absorbs *n* packed rows of
the same class, and declines when the classes differ, after which the no-nesting
invariant requires every element to be materialised.

Neither obvious fix works. Making `List` packed-aware would let a plain `List`
hold `EXPR_NDARRAY` elements — the shape the transparency gate exists to prevent
— and would turn every unaware head's O(argc) scan into O(tree). Widening the
integer row to float64 is a value change. The real options are a heterogeneous
packed tuple (a design change) or a diagnostic that names the offending element.
**The workaround is one character** (`1. UnitStep[…]`) and should be documented
before anything is built.

**Not an instance of this: `DiagonalMatrix` (resolved 2026-08-01).** Its `Real`
diagonal was filed here because its exact zeros made the matrix two-headed. The
zeros were simply *wrong* — Mathematica gives all `Real`s — and correcting them
made it one dtype, 320× → 1.09× NumPy. The distinction is worth keeping: 10.1 is
about dtypes the caller genuinely chose and the callee must carry; that was
about a dtype the callee invented and got wrong. Before filing something here,
check the offending element is one the user actually supplied.

**10.2 — The packing threshold is a floor under numeric linear algebra.**
`builtin_inverse` dispatches on `linalg_call_has_ndarray`, so a 6×6 `Real` matrix
— 36 elements, correctly never packed — cannot reach `ndla_inverse` and runs a
fraction-free symbolic Gauss-Jordan at **~800 µs**, against ~1 µs for
`dgetrf`+`dgetri`. A 2×2 costs 26.8 µs; a 6×6 matrix–vector product 12.1 µs. That
is the whole of the 21× gap on the Kalman row.

The fix is to dispatch on "is this a matrix of machine numbers" rather than "is
this already a buffer", at `Inverse`, `Dot`, `LinearSolve` and `Det`. **There is
a trap**: `ndla_inverse` declines by calling `linalg_delist_and_reeval`, which
re-evaluates the call — so a naive lift re-enters `builtin_inverse`, packs again,
declines again, and loops forever, on exactly the singular inputs that are
hardest to test. Each entry point needs a `_try` form that returns `NULL` instead
of degrading. Do that first, then lift.

**10.3 — Mixed `float64 × int64` elementwise is ~25× a pure float multiply.**
`da1 UnitStep[z1]` on 131072 elements costs 1.37 ms against ~0.05 ms for a
float64 product of the same size, because the mixed pair goes through the generic
`ndt_get`/`ndt_set` accessors instead of a widening loop. This is the shape
*every* comparison mask produces, so it is in the inner loop of all branch-free
array code — the MLP's ReLU derivative, the ray tracer's hit mask, the LJ cut-off.
One dtype pair in `ndarray_elementwise`.

**10.4 — Per-call dispatch is now the ceiling on short-body loops.** Three
experiments hit the same wall independently: the spectral solvers (16) spend
~130 µs of a 166 µs time step outside FFTW and the arithmetic; the option pricers
(15) sit at 1.55–3.3× NumPy on 4000–25000 interpreter iterations; the Kalman
filter (18) is nothing but this. No further packing work moves any of them. This
is Phase 5.2, and it now has three independent profiles asking for it rather than
a hypothesis.

**10.5 — No strided views, twice over.** Phase 9.1 asked for a transposed view;
experiment 14 adds the sliding window (`Partition[a, k, 1]` materialises 6 × 10⁶
elements where NumPy's `sliding_window_view` returns a view and copies nothing).
These are the *same* missing feature — a stride vector on `NDArray` that every
consumer honours — and closing either properly closes both. Worth scoping as one
piece of work rather than two.

**10.6 — `RowReduce` of a `Real` matrix still leaks exact `Integer` 1s.**
`RowReduce[{{2., 0.}, {0., 4.}}]` gives `{{1, 0.}, {0., 1}}` where Mathematica
gives all `Real`s. `exact_div_wrapper`'s zero was fixed by this sweep; the
diagonal 1 comes from the divide/eliminate steps and is **pre-existing** — there
is already a comment about it in `src/linalg/inv.c`.

**DONE 2026-08-01, and note what checking cost.** This was first filed as an
instance of 10.1 ("a mixed matrix cannot be packed"), then re-filed hours later
as an exactness bug on the argument that `RowReduce` of a machine matrix "should
be uniformly machine-real, exactly as `DiagonalMatrix` should be" — reasoning by
analogy from a fix that had just worked.

**Mathematica does not agree.** Its `RowReduce[{{2., 4.}, {1., 3.}}]` is
`{{1, 0.}, {0, 1}}`, element heads `{Integer, Real, Integer, Integer}`, and
`Developer`PackedArrayQ` on it is `False`. Mathematica's own RREF of a machine
matrix is two-headed and unpacked. So the second filing was as unchecked as the
first; it happened to point at a change worth making, for a reason that was not
true.

The change landed anyway, as a **stated divergence**: Mathilda follows the
project rule that a routine handed a machine array answers with a machine array,
so the pivot 1 and the eliminated 0 take the input's exactness and the result
packs. Values agree with Mathematica; only the heads of the structural entries
differ, and only for an inexact input. Pinned in `test_packed_list.c` with the
divergence spelled out beside the assertion.

*Twice in one day an exactness claim about Mathematica was written down without
being run. The tool that ends that is `tools/check_array_exactness.py`, whose
EXEMPT table refuses an entry without the output that justifies it.*

**10.6a — The rest of the two-headed surface (done 2026-08-01).** Sweeping the
element heads of 342 results found five more routines inventing an exact element
inside a machine-real one: `VandermondeMatrix`'s `x^0` column, `HankelMatrix`'s
pad zeros, `ToeplitzMatrix`'s un-widened entries, `MatrixPower[m, 0]`'s
identity, and `NullSpace`'s free-variable slot. All fixed, all now pack, all
verified against Mathematica — which agrees on every one of those five.

The eight *legitimate* mixtures are the more useful half of the output, because
they are where the rule stops: `Chop`'s zero is exact on purpose, the
`PadRight`/`Riffle`/`Insert`/`Append`/`ReplacePart` family places a
caller-supplied element verbatim, `Join` mixes what it was given, `Tally`
returns `{value, count}` pairs, and `LUDecomposition`'s exact elements are the
permutation vector in a three-part tuple.

**Still open in this class:** `NullSpace` of an inexact matrix is not
orthonormalized — Mathematica's `NullSpace[{{1., 2.}, {2., 4.}}]` is
`{{-0.894…, 0.447…}}` against our `{{-2., 1.}}`. Same basis, different
normalisation; a larger behavioural change than exactness and not attempted
here.

**10.7 — No `SparseArray`.** Experiment 12 works because a uniform-degree graph
is a dense `n × d` neighbour matrix, which is a real and common case. A general
sparse matrix needs a row-pointer array and a segmented reduction that does not
go through `Partition`.

**10.8 — The gather still writes a position array.** `build_axis_selector`
materialises an `int64*` of the full index length before the copy loop, so a
1.6 × 10⁶-index gather writes 13 MB it then reads once. Reading the index buffer
*inside* the copy loop removes the pass. Same allocation as item 9.5's spans, and
worth most of the remaining 2.7× against NumPy on the graph rows.

**Still the highest priority overall: HPC 7.4, the pre-existing `NDSolve`
segfault.** It is a crash, not a slowdown, and it has now outlived three sweeps.

---

## 11. The coverage sweep's register (experiment 20)

Every item below is **measured**, not estimated: the figures come from
`tools/numeric_sweep.py`, which times 283 probes over the numeric surface
against NumPy/SciPy on the same data. Re-run any row with
`python3 tools/numeric_sweep.py --only <id>`.

These are lettered C.x rather than folded into the numbered phases, because they
did not come from the same place. Phases 1–10 were each derived from a
**workload**: measure a kernel, find what makes it slow. This register comes
from **enumerating all 676 registered builtins** and asking of each whether a
numeric argument reaches machine code — which is why it contains heads no
benchmark in this tree had ever touched.

| # | item | measured | why it is where it is |
|---|---|---|---|
| C.1 | **A boolean array dtype.** `Map[# > 0.5 &, v]` builds 10⁶ `True`/`False` symbols and every consumer walks a boxed list | 432 ms vs 238 µs (**1813×**); `MapThread[And, …]` **4181×** | The largest structural gap on the numeric surface, and the one design change here. `UnitStep[v - 0.5]` is the fast spelling at 3.6 ms, so the *arithmetic* mask already works — it is the *boolean* one that has no representation. The compiled-predicate work (experiment 20) sidesteps it for the seven `Select`-family heads; `Pick`, `Count`, `Position` and boolean algebra still pay. Interacts with the mixed-dtype tuple problem of experiment 13: a bool buffer is a third dtype the no-nesting invariant has to carry |
| C.2 | ~~**`MatrixPower` on a machine matrix** takes an exact/symbolic path~~ **DONE 2026-08-01** | 14.7 s → **827 µs**, **1.08× NumPy** | It was one line, and not the line the note guessed: `builtin_matrixpower` opened with `linalg_delist_and_reeval`, so square-and-multiply ran over boxed matrices. Nothing in it reads an element; it now keeps the buffer and `dot2` reaches `dgemm` |
| C.3 | **The Bessel kernels are registered but do not reach machine code.** `ProductLog` **DONE 2026-08-01** | `BesselK` 44.0 s, `BesselI` 25.5 s, `BesselY` 18.6 s, `BesselJ` 10.3 s per 10⁶. `ProductLog` 28.1 s → **143 ms** | Diagnosed 2026-08-01, and it is **two** causes, neither the guessed one. (a) The four Bessel heads carry **DownValues** from `internal/*.m` (half-integer recurrences, negative-order reflection), and the gate's aware test is `packed_aware && !down_values` — so a registered kernel is unreachable the moment a head acquires a rewrite rule. Every one of those rules is guarded by `Not[NumberQ[z]]` and cannot fire on numeric data, but nothing knows that. (b) The scalar builtin uses `mpfr_jn` (32 µs) while the kernel uses libm `jn` (~0.1 µs), so making the head aware without unifying them would make the packed and plain paths disagree in the last ulp — an N2 violation. **The fix is to give the scalar path the libm kernel**, which makes both consistent AND 300× faster, and needs its own accuracy comparison against MPFR before it lands. `ProductLog` turned out to be a different problem again — see the C.3a row below |
| C.3a | ~~`ProductLog` is 28 µs/element~~ **DONE 2026-08-01, and it was a WRONG ANSWER** | 2.77 s → **14.3 ms** per 10⁵ | The array kernel was failing on ~1 element in 10⁵ and abandoning the whole buffer to MPFR — that was the *symptom*. The *cause* was that `sf_machine_productlog` used the `x → ∞` seed for every `x ≥ 1`, which is `log(0)/0` at `x = 1` and diverges just above it: `ProductLog[1.01]` returned **−338.392** for an answer of 0.5707. Only the fast path was wrong (the interpreter goes to MPFR), so nothing compared them. The kernel now seeds from `log(1+x)` on `1 ≤ x ≤ e` and **verifies `w + log w = log x` before answering** |
| C.4 | ~~**`Extract`** materialises the whole buffer to read three elements~~ **DONE 2026-08-01** | 99.7 ms → **1.0 µs**, **1/1.35× NumPy** | Not the gate: `Extract` was already on `AWARE`. `expr_part` cannot index an NDArray (`is_atomic` is true for one), so it returned NULL and the **post**-gate materialised on the way to rest. Now uses `ndarray_part`, and joins `INT64_OK` |
| C.5 | **The producers that do not build buffers.** `Subdivide`, `Rescale`, `IdentityMatrix`, `DiagonalMatrix`, `UnitVector` **DONE 2026-08-01**; `Table[N[i], …]` and `Array` remain | `Subdivide` 1.87 s → **707 µs** (**1/1.29× NumPy**), `DiagonalMatrix` 70.2 ms → **263 µs** (**1.09×**), `IdentityMatrix` 77 ms → **223 µs** (1.16×), `UnitVector` 130 ms → **720 µs** (3.33×), `Rescale` 2.17 s → **3.1 ms**. `Array[N, 10⁶]` 444 ms and `Table[N[i], …]` 446 ms unchanged | Mostly mechanical, and two things it turned up. (1) `Rescale` was not a buffer problem: it built a fresh `Rescale[element, range]` call *per element*, and the fix is to build the arithmetic once over the whole argument and let Listable `Plus`/`Times` thread. (2) `Subdivide`'s interior points **changed in the last bit**, deliberately — the old value came out of `Times` sorting its `Orderless` factors by value, not from any rounding rule, and the new `min + i·step` is bit-identical to `numpy.linspace`. (3) The `Real` `DiagonalMatrix` row was an **exactness bug**, not a packing limit. It had been filed as an instance of item 10.1 (the mixed-dtype gap) on the grounds that its exact zeros were Mathematica's answer; they are not, and once corrected the matrix is one dtype and packs. 10.1 itself is untouched — a tuple of genuinely different dtypes is still a real problem; `DiagonalMatrix` was never an instance of it. **`Table` and `Array` with an exact iterator are the open half** — see the note below the table |
| C.6 | **`Insert`/`Delete`/`Append`/`Prepend`** are "correct by omission" in `pack.c` and box the whole list | 150–165 ms vs 0.8–1.1 ms (**153–198×**) | NumPy also copies, so this is a real gap rather than a view artifact |
| C.7 | **`MemberQ`/`Count`/`Position`** scan a boxed list | 132 ms / 164 ms (**~570×**) | `Count` and `Position` over a mask are the other half of 12.1 |
| C.8 | **The integer band, again.** `Mod`/`Quotient` have registered binary kernels that decline `int64` | `Mod[iv, 1000]` 572 ms vs 3.9 ms (**146×**); integer group median **70×** | Ninth item on performance.md §15's list, still open. `BitAnd`/`BitOr`/`BitXor`/`BitShiftLeft` do not exist at all |
| C.9 | **Small dense matrices never reach LAPACK** | 6×6 `Det` **122×**, `LinearSolve` **85×**, `Inverse` **75×** | Confirms experiment 18 by direct measurement rather than through an application. Same item as roadmap #3 |
| C.10 | **`N[list]` 367 ms, `Chop` 179 ms, `Catenate` 369 ms, `Sort` on `int64` 390 ms** | 16–315× | Unexamined; each is probably a single missing buffer walk |
| C.11 | ~~`LeastSquares[A500,b500]` and `PseudoInverse[A300]` do not complete in 180 s~~ **DONE 2026-08-01**. `MatrixExp` was a misreading | `LeastSquares` **70.9 ms** (1.69× NumPy), `PseudoInverse` **19.8 ms** (**1/1.07× NumPy**) | Both terminate — the exact pipeline in `inv.c` rationalises a 300×300 machine matrix and row-reduces over ℚ, which is not a hang so much as a wrong choice of algorithm at that size. Both now take one thin `gesdd`. `MatrixExp` does **not** exist (C.12); it returns unevaluated in a microsecond, and the sweep's 180 s was the Python side |
| C.12 | **96 numeric heads do not exist**, including `CholeskyDecomposition`, `KroneckerProduct`, `MatrixExp`, `Diagonal`, `ArrayReshape`, `Ordering`, `BinCounts`, `Quantile`, `Correlation`, `Covariance`, `Standardize`, `SparseArray`, the Chebyshev/Hermite/Laguerre/Jacobi polynomials and the elliptic integrals | — | A **coverage** gap, not a speed gap, and out of scope for an HPC plan — but it is what `tools/numeric_coverage.py --missing` reports and it belongs somewhere. Several (`Diagonal`, `ArrayReshape`, `Ordering`) are one-liners over an existing buffer |

**The open half of C.5: an exact iterator.** `Table[N[i], {i, 10⁶}]` and
`Array[N, 10⁶]` both run one interpreter evaluation per element (~440 ns) and
box the result, and neither is fixed by a buffer alone — the cost is the
evaluation, so the fix is to compile the body. `Table`'s existing auto-compile
is gated on an **inexact iterator** for a good reason: compiling with a real-typed
variable and an exact iterator would answer `Table[i^2, …]` in floats where the
interpreter gives exact Integers.

The sound generalisation is to declare the iterator **`CT_INT`** and require the
compiled program's static result type to be `CT_INT` as well; then the elements
are machine integers and the interpreter would have produced Integers too, for
everything in the compilable subset (`i^2`, `2i+1`, `Mod[i,7]`, `Floor[i/2]`,
`If[i>3,1,2]`). `CT_REAL` out of a `CT_INT` in must **not** be accepted —
`i/2` is a Rational and `Sqrt[i]` is symbolic — which is exactly why
`Table[N[i], …]` needs a separate argument (`N` is *defined* to produce an
inexact result) rather than falling out of the type rule. Not attempted here.

**Where the surface already stands**, so that this list is read in proportion:
the elementwise group's median is **1/1.07× NumPy** and the special-function
group's **1/1.34×** — both at or ahead of parity, across 89 probes. The gaps
above are concentrated in masks, construction, integer arrays and small
matrices, not spread across the numeric surface.

---

## 12. Risks

| risk | mitigation |
|---|---|
| LAPACK eigen/SVD conventions (order, sign) differ from the in-house paths and break parity | Ordering/sign shim, differential-tested against the current implementation across the whole `eigen_tests`/`svd` corpus; keep the in-house path as the fallback for degenerate cases |
| A narrowing kernel writes `1.0` where the list gives `1`, or wraps on overflow | The `NOT_AWARE` sweep that *created* that list is the acceptance test; overflow abandons the whole result via the `ci_*_i64` contract |
| Threaded `Accumulate` reassociates and changes the last ULP | Decided up front (§4), not discovered: either accept a documented deterministic-but-different result, or thread only the exact int64 arm |
| Phase 5.2 (dispatch) turns into an open-ended evaluator rewrite | Profile-first gate: no work starts until the cost is attributed to a specific structure |
| A fast path is added but never reached, and nobody notices | Phase 0.2 and 0.3 exist precisely for this; it has already happened four times |
| FFTW plan cache grows without bound or is not thread-safe | Bounded LRU keyed on shape; FFTW planning is not thread-safe, so guard plan creation with a mutex and use `fftw_execute_dft` (which is) for execution |

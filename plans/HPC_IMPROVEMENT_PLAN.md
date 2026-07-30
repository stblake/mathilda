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

## Phase 2 — the narrowing kernel category

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

**A quarter of the time is the kernel mapping and unmapping memory.** A 1000×1000
double matrix is 8 MB, and macOS's allocator serves allocations that large
straight from `mmap`, returning them with `munmap`/`madvise` — so every call pays
page-table work for buffers that are immediately discarded and re-requested. A
small size-keyed scratch pool in the linalg bridges would remove it. That is now
the largest remaining item here, and it was not in this plan: it only became
visible once the element-wise loop stopped hiding it.

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

## 11. Risks

| risk | mitigation |
|---|---|
| LAPACK eigen/SVD conventions (order, sign) differ from the in-house paths and break parity | Ordering/sign shim, differential-tested against the current implementation across the whole `eigen_tests`/`svd` corpus; keep the in-house path as the fallback for degenerate cases |
| A narrowing kernel writes `1.0` where the list gives `1`, or wraps on overflow | The `NOT_AWARE` sweep that *created* that list is the acceptance test; overflow abandons the whole result via the `ci_*_i64` contract |
| Threaded `Accumulate` reassociates and changes the last ULP | Decided up front (§4), not discovered: either accept a documented deterministic-but-different result, or thread only the exact int64 arm |
| Phase 5.2 (dispatch) turns into an open-ended evaluator rewrite | Profile-first gate: no work starts until the cost is attributed to a specific structure |
| A fast path is added but never reached, and nobody notices | Phase 0.2 and 0.3 exist precisely for this; it has already happened four times |
| FFTW plan cache grows without bound or is not thread-safe | Bounded LRU keyed on shape; FFTW planning is not thread-safe, so guard plan creation with a mutex and use `fftw_execute_dft` (which is) for execution |

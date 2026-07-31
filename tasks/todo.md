# Fourth HPC sweep — closing the distance to NumPy (2026-07-31)

Plus: eleven per-experiment writeups covering 2026-07-26 → 07-31, in
`docs/experiments/`.

Previous sweeps: `plans/PACKED_ARRAYS_TODO.md`, `docs/design/performance.md`
§8–§10.

## Plan

- [x] Probe every primitive the eight application kernels use, against NumPy
- [x] Fix 1 — `First`/`Last`/`Most`/`Rest` read the buffer instead of delisting
- [x] Fix 2 — `Clip` back on the buffer, with an exactness gate that is correct
- [x] Fix 3 — the scan (`Fold`/`FoldList`) reads a buffer and emits a buffer
- [x] Fix 4 — scan kernels for a bare `Max`/`Min`/`Plus`/`Times` operator
- [x] Fix 5 — `ListConvolve`/`ListCorrelate` on the buffer (HPC plan 7.1)
- [x] Fix 6 — `Accumulate`, `Differences`, `Part` spans, `Outer` (re-probe)
- [x] Verify: differential tests, full suite, gates, valgrind
- [x] Re-run the 51-benchmark comparison, three systems
- [x] Write `docs/experiments/<NAME>.md` for every experiment of the last 4 days
- [x] Document: `performance.md` §10, weekly changelog, HPC plan Phase 9

## Results

### The probe that chose the work (10⁶ float64, before)

| op | Mathilda | NumPy | ratio |
|---|---:|---:|---:|
| `Most[v]` | 270 ms | 0.48 ms | 563× |
| `Rest[v]` | 223 ms | 0.48 ms | 466× |
| `FoldList[a #1+b #2 &, …]` packed | 1246 ms | 3.02 ms | 413× |
| `FoldList[Max, …]` | 605 ms | 2.73 ms | 222× |
| `ListConvolve[k, v]` | 263 ms | 1.40 ms | 188× |
| `First[v]` | 123 ms | ~0 | O(1) done in O(n) |
| `Clip[v, {a,b}]` | 356 ms | 4.66 ms | 76× |
| `v v` | 0.81 ms | 0.45 ms | 1.8× ✓ |
| `Log[v]` | 8.8 ms | 6.8 ms | 1.3× ✓ |

The last two rows are the control, and they settled where to look: **the
elementwise kernels were already at the memory floor**, so nothing above them
was about arithmetic. Everything expensive was structural or a scan.

### After

| primitive | before | after |
|---|---:|---:|
| `First` / `Last` | 123 / 103 ms | **0.01 / 0.00 ms** |
| `Most` / `Rest` | 270 / 223 ms | **3.42 / 0.71 ms** |
| `Clip` | 356 ms | **0.93 ms** |
| `FoldList[Max, …]` | 605 ms | **2.00 ms** (NumPy 2.73) |
| EMA `FoldList` | 1246 ms | **22.6 ms** |
| `ListCorrelate` 5×5 / 1024² | 376 ms | **46.5 ms** (SciPy 52.5) |
| `ListConvolve` 5-tap / 10⁶ | 263 ms | **16.1 ms** |
| `Accumulate` | 7.9 ms | **3.5 ms** |
| `Differences` | 6.6 ms | **0.97 ms** |
| `Part` span | 14.0 ms | **7.4 ms** |
| `Outer`, one small operand | 849 ms | **10.9 ms** |

| benchmark | before | after | WL 14.0 | NumPy | |
|---|---:|---:|---:|---:|---|
| Return series, 10⁶ | 2.02 s | **55.6 ms** | 154 ms | 30.4 ms | **2.78× faster than WL** |
| Gaussian blur + Sobel, 1024² | 2.19 s | **103 ms** | 94.5 ms | 114 ms | **1.10× ahead of NumPy** |
| `ListConvolve`, 1024², 5×5 | 312 ms | **35.7 ms** | 33.8 ms | 39.7 ms | **1.11× ahead of NumPy** |
| `Differences`, 10⁷ | 56.6 ms | **16.3 ms** | 14.3 ms | 15.6 ms | |
| `Accumulate`, 10⁷ | 63.8 ms | **19.4 ms** | 15.1 ms | 28.7 ms | **1.48× ahead of NumPy** |

Full run: 51 benchmarks, three systems, no value mismatches. Of 51, Mathilda is
faster than Mathematica on 15 and within 1.5× on 18 more; against NumPy, faster
on 20 and within 1.5× on 10 more.

## Review

### What the probe method was worth

The previous three sweeps started from *applications* and worked inward. This
one started from *primitives measured against NumPy*, and that ordering is why
`First[v]` was found: at the application level a 123 ms element read is
invisible, and in a profile it looks like the function you called. It only
became obvious next to `Drop[v, 250]` at 0.88 ms on the same data.

**Always probe with a control that should cost the same.** Two of the eleven
rows in the probe existed only as controls, and they are what said "stop looking
at arithmetic".

### Two fixes that were not what they looked like

**The convolution's cost was an integer division.** After the buffer path and
the interior hoist, I added an affine-stride inner loop for the rank-1 case —
and it measured *worse*. Not because it was wrong, but because
`(o / Lstr[ax]) % Ldims[ax]` per output was swamping it: 10⁶ 64-bit divisions
against 5×10⁶ multiply-adds, with runtime strides the compiler cannot
strength-reduce. Replacing it with an odometer took 31.5 ms → 16.1 ms and let
the vectorisable loop that was already there finally show.

**`Clip` had been disabled at the wrong granularity.** It was struck from the
aware list entirely because an exact bound makes the result non-uniform. The
exactness problem was real, but it belonged to the **bounds**, not to the head —
and the blanket fix cost 356 ms per call on data where nothing was wrong. Worse,
it left the visible-`NDArray[…]` case *incorrect*, because the gate never
materialises those, so `ndstruct_clip` ran anyway with the exact bounds coerced
to doubles.

### Things that were wrong about my own work

- **I put three NumPy figures into an experiment writeup that I had not
  measured.** Caught and replaced with placeholders before anything else was
  written, then filled from a real run. Numbers in a document that were guessed
  are worse than no numbers.
- **My first summary counts in `performance.md` were wrong** (15/13/22 and
  17/12/20). They are now computed from the run's own JSON: 15/18/17 against
  Mathematica, 20/10/19 against NumPy.
- **My `numloop` change lost the NDArray presentation**, so the same `FoldList`
  answered with head `NDArray` for a `Plus` operator and head `List` for the
  equivalent lambda. Caught by `compile_tests`' cross-spelling parity check, not
  by anything I wrote.

### What the tests missed — three times

1. **A test asserted the fallback, not the invariant.**
   `chk_eq("NDArrayQ[Clip[NDArray[{-2.,0.,2.}]]]", "True")` asserted that `Clip`
   answered `{-1., 0., 1.}` where the List path answers `{-1, 0., 1}` — it was
   enshrining a wrong answer. The `chk_array` checks beside it could not catch
   it: they compare numeric *distance*, and the difference is in element
   **heads**.
2. **A regression visible only across two spellings** (the `FoldList` one above).
3. **Tests below `PACK_MIN_ELEMENTS` test nothing.** Second appearance; every
   new case in `test_fourth_sweep_fast_paths` uses data above the threshold.

### A finding that was not part of the plan

The Jacobi row read 223 ms in a full benchmark run and 128 ms measured alone,
while every other row was stable. Bisecting the benchmark *prefix* found one
trigger: running a single 1000×1000 `dgemm` first makes Mathilda's threaded
stencil **1.45× slower**, reproducibly to ±1%. Accelerate's worker threads
persist and compete with `nd_parallel_for`.

That is a real property of the system, not benchmark noise, and it affects any
program that multiplies a matrix and then does array work. Recorded as HPC plan
item 9.3 rather than fixed — the fix is a thread-pool policy question, and the
last thread-pool hypothesis in this plan (the linalg scratch pool) was wrong.

### Verification

- Full 395-binary suite: clean apart from one **pre-existing** `simplify_tests`
  failure (a symbolic radical case), confirmed pre-existing in the third sweep.
- `test_fourth_sweep_fast_paths` — ~200 differential cases over all six paths,
  each packed form against the identical plain one, every form each path must
  **decline**, and explicit presentation-parity assertions.
- `test_ndarray_reduce.c`'s `Clip` assertions rewritten to be sensitive to
  element heads.
- `make check-c99`, `make check-packed-aware` pass.
- Valgrind unchanged from the macOS start-up baseline.

## The experiment writeups

`docs/experiments/` — eleven files, one per experiment, 2026-07-26 → 07-31, each
comparing Mathilda against **both** Mathematica 14.0 and Python.

The reason for both columns, stated once in the index and borne out repeatedly:
Mathematica says whether we are behind a *competitor*; NumPy — on this host
linking the same Accelerate BLAS — says whether we are behind *the machine*. A
row can read acceptably against one and badly against the other, and that
difference is usually the finding.

## Left open (HPC plan Phase 9)

1. **`Transpose` copies where NumPy views** — the largest application gap
   (logistic regression, 3.06× NumPy). Needs a strided NDArray view, which is a
   design change rather than a fast path.
2. **No interpreter-level fusion** — k-means, 2.92× NumPy.
3. **`dgemm` vs our thread pool**, above.
4. 1-D convolution at 11× NumPy; `Part` span's position array; `Tally`;
   `Reverse`/`RotateLeft`; `ArrayReshape` does not exist.

**Still the highest priority overall: the pre-existing `NDSolve` segfault
(HPC 7.4).** It is a crash, not a slowdown, and it has now outlived two sweeps.

# QRDecomposition Performance Plan

## 1. Motivation — measured baseline (2026-05-21)

`QRDecomposition` currently drives Modified Gram-Schmidt through the
evaluator. Inexact inputs go through the shared
`rationalise → exact MGS → numericalise` round-trip used by
`PseudoInverse` / `Eigenvalues`. The cost dominates everything: every
scalar add, multiply, conjugate and Together travels the head lookup →
attribute resolution → listable thread → builtin-dispatch path.

Side-by-side timings (Mathematica 13.x vs current Mathilda, same
random seed, single thread, M-series Mac):

| Case                | Mathematica | Mathilda    | Slowdown |
|---------------------|-------------|-------------|----------|
| machine real n=5    |    1.7 µs   |    9.4 ms   | ~5,500×  |
| machine real n=8    |     ~5 µs   |   87.9 ms   | ~17,000× |
| machine real n=10   |    3.2 µs   |  > 10 min   | enormous |
| MPFR-30 n=5         |     65 µs   |    7.5 ms   | ~115×    |
| MPFR-30 n=10        |    312 µs   |  > 10 min   | enormous |
| symbolic int n=3    |    453 µs   |    418 µs   | **~1×**  |
| symbolic int n=4    |    696 µs   |    858 µs   | ~1.2×    |

**Symbolic is already at parity** with Mathematica. The fix targets
**only** the inexact paths: machine precision (gap of ~10⁴×) and
arbitrary precision (gap of ~10²–10³×).

Cause of the gap:
1. We have no native double / MPFR kernel — every op crosses the
   evaluator boundary.
2. Inexact inputs are first rationalised, exactly computed, then
   numericalised, which is wholly wasted work for floating-point input.

## 2. Goals

1. **Machine-precision kernel**: 64-bit-float QR via LAPACK (`dgeqp3`,
   `zgeqp3`, `dorgqr`, `zungqr`), bypassing the evaluator entirely.
   Target: within ~5× of Mathematica on n ≤ 100 (Mathematica also calls
   LAPACK; the residual gap is just our packing/unpacking overhead and
   smaller-n call overhead).
2. **Arbitrary-precision kernel**: MPFR-based Householder QR, with
   complex stored as a pair of MPFR arrays (matching the existing
   eigen-module convention — no MPC dependency). Target: within ~3×
   of Mathematica at 30-digit and 100-digit precision.
3. **Dispatch**: classify inputs once with `common_scan_inexact()`,
   route to LAPACK / MPFR / symbolic. No regression on the symbolic
   path.
4. **BLAS / LAPACK integration**: a project-wide capability, not just
   QR — wire `USE_LAPACK` like `USE_ECM` and `USE_MPFR`, with
   Accelerate on macOS and OpenBLAS+LAPACKE on Linux. Other linalg
   builtins (`Inverse`, `LinearSolve`, `Det`, `Eigenvalues`, ...) can
   later opt into the same kernel.

Non-goals:
- Tuning the symbolic path (already at Mathematica parity).
- Adding a third arbitrary-precision library (MPC, FLINT, Arb). MPFR
  pair-of-arrays is the existing convention and proven adequate by the
  eigen kernels.
- GPU / multi-threaded execution (BLAS already threads internally on
  Accelerate / OpenBLAS — we don't need to do anything).

## 3. Architecture

Two new kernels behind a single dispatcher, mirroring `eigen.c`:

```
QRDecomposition[m, opts]              (src/linalg/qrdecomp.c)
        │
        ▼
qr_dispatch(m, opts)
   │
   ├─ common_scan_inexact(m)  →  { has_inexact, min_bits }
   │
   ├─ has_inexact && min_bits ≤ 53        ───► qr_machine    (new)
   │                                              │
   │                                              ▼
   │                                       LAPACK: dgeqp3 / zgeqp3,
   │                                       dorgqr / zungqr
   │
   ├─ has_inexact && min_bits >  53        ───► qr_mpfr      (new)
   │                                              │
   │                                              ▼
   │                                       Householder QR over
   │                                       MPFR arrays
   │
   └─ otherwise (exact / symbolic)         ───► qr_symbolic  (existing,
                                                              renamed)
```

The existing symbolic implementation becomes `qr_symbolic` (renamed
from `mgs_core` + wrappers); the entry point `builtin_qrdecomposition`
is reduced to option parsing + dispatch.

## 4. File layout

New files (mirroring the eigen split):

```
src/linalg/
├── qrdecomp.c           (existing, becomes thin dispatcher + symbolic)
├── qrdecomp.h           (existing)
├── qrdecomp_internal.h  (NEW — shared MatD/MatM facade + kernel decls)
├── qrdecomp_machine.c   (NEW — LAPACK kernel)
└── qrdecomp_mpfr.c      (NEW — MPFR Householder kernel)
```

Where possible we reuse, not copy:

- `MatD`, `MatM` come from `src/linalg/eigen_internal.h` (line 49–67).
  We need to widen them slightly because QR has non-square inputs.
  Add `size_t rows, cols` fields and keep `n` as an alias **or** add
  a separate `MatDRect` / `MatMRect`. Decision: **widen MatD / MatM in
  place** — non-square eigenvalues don't exist, so the field is dead
  code there, but it removes duplication. Coordinate with eigen owner
  (one-line change, all eigen callsites set `rows = cols = n` until they
  also opt-in).
- `matD_load`, `matM_load`, `mpfr_array_alloc/free` already exist in
  `src/linalg/eigen_common.c` and become the shared loader once they
  accept the widened type.
- `common_scan_inexact` from `src/common.h` provides the dispatch
  signal.

## 5. BLAS / LAPACK integration (project-wide)

### 5.1 Distribution strategy — four-tier autodetection

The goal is **`git clone && make` works for ~every user on ~every
platform**, with no manual configure step. The detection ladder runs
top to bottom, taking the first hit:

| Tier | Trigger                                                       | Effect                                                                 |
|------|---------------------------------------------------------------|------------------------------------------------------------------------|
| 1    | `uname -s == Darwin`                                          | Link `-framework Accelerate`. **Zero user install** — shipped with every macOS since 10.4. |
| 2    | `pkg-config --exists lapacke` succeeds                        | Use `pkg-config --libs lapacke` (covers OpenBLAS, MKL, Netlib LAPACKE on Debian/Ubuntu/Fedora/Arch/Nix). |
| 3    | `lapacke.h` present under `/usr/include` or `/usr/local/include` | Link `-llapacke -llapack -lblas` (Debian fallback when pkg-config isn't installed). |
| 4    | Nothing found                                                 | Print one-line warning, set `USE_LAPACK := 0`, build proceeds via existing MPFR/symbolic path. **No build failure.** |

The auto-degrade in tier 4 matches the existing `USE_MPFR=0` /
`USE_ECM=0` policy in the makefile, preserving today's "build always
works" property.

Optional fifth tier — `USE_VENDORED_LAPACK=1` — mirrors `build_ecm.sh`:
`build_openblas.sh` clones a pinned OpenBLAS release into
`src/external/openblas/`, builds a static archive, links it. **Off by
default** because it pulls in a Fortran compiler dependency and adds
minutes to the first build.

README install snippet:

```
macOS:           nothing (Accelerate is built-in)
Ubuntu/Debian:   sudo apt install liblapacke-dev libopenblas-dev
Fedora/RHEL:     sudo dnf install lapack-devel lapacke-devel openblas-devel
Arch:            sudo pacman -S openblas lapacke
NixOS:           use the nix-shell shipped in shell.nix
```

### 5.2 API surface — Fortran ABI behind a thin wrapper

Apple's Accelerate framework ships the **Fortran ABI** LAPACK
(column-major, trailing-underscore symbols) plus CBLAS. It does
**not** ship `LAPACKE_*` row-major wrappers, even with
`ACCELERATE_NEW_LAPACK` defined — those would defeat the
zero-install property by requiring a Homebrew dep.

The portable lowest-common-denominator that works on every
distribution (Accelerate, OpenBLAS, MKL, Netlib) is the **Fortran
ABI** itself: `dgeqp3_`, `dgeqrf_`, `dorgqr_`, with pointers for
every scalar and column-major storage. We hide it once in:

- `src/linalg/lapack.h` — platform-agnostic header. Includes
  `<Accelerate/Accelerate.h>` on macOS, `<cblas.h>` on Linux, and
  declares `extern` Fortran prototypes for the LAPACK routines we
  use (so the project never needs `<lapacke.h>`).
- `src/linalg/lapack.c` — thin row-major C wrappers
  `mat_lapack_dgeqp3(...)`, etc. that handle the column-major
  transpose and pointer-to-int Fortran calling convention. Kernel
  code only touches these wrappers, never the raw Fortran symbols.

This is the same pattern used by NumPy / SciPy / Eigen / BLIS — call
the Fortran ABI, wrap once, hide the rest.

### 5.3 Build wiring

`makefile` — add after the `USE_MPFR` block:

```make
USE_LAPACK ?= 1
ifeq ($(USE_LAPACK), 1)
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    CFLAGS  += -DUSE_LAPACK -DMATHILDA_USE_ACCELERATE
    LDFLAGS += -framework Accelerate
  else ifneq ($(shell pkg-config --exists lapacke && echo y 2>/dev/null),)
    CFLAGS  += -DUSE_LAPACK $(shell pkg-config --cflags lapacke)
    LDFLAGS += $(shell pkg-config --libs lapacke)
  else ifneq ($(wildcard /usr/include/lapacke.h)$(wildcard /usr/local/include/lapacke.h),)
    CFLAGS  += -DUSE_LAPACK
    LDFLAGS += -llapacke -llapack -lblas
  else
    $(warning LAPACK/LAPACKE not detected; building with USE_LAPACK=0)
    override USE_LAPACK := 0
  endif
endif
```

`tests/CMakeLists.txt` — `option(USE_LAPACK ON)` plus
`find_library(ACCELERATE Accelerate)` on Apple,
`find_package(LAPACK)` + `find_path(LAPACKE_INCLUDE_DIR lapacke.h)` on
Linux, same auto-off on miss.

When `USE_LAPACK=0`, the machine kernel falls back to the existing
symbolic/MPFR path with a one-time runtime warning (same pattern as
the existing `matsol_warn_once` helper).

### 5.4 LAPACK routines used by QR

| Routine    | Purpose                                       | Wrapper name              |
|------------|-----------------------------------------------|---------------------------|
| `dgeqp3_`  | real QR with column pivoting                  | `mat_lapack_dgeqp3`       |
| `zgeqp3_`  | complex QR with column pivoting               | `mat_lapack_zgeqp3`       |
| `dgeqrf_`  | real QR (no pivoting; faster path)            | `mat_lapack_dgeqrf`       |
| `zgeqrf_`  | complex QR (no pivoting)                      | `mat_lapack_zgeqrf`       |
| `dorgqr_`  | form Q from real Householder reflectors       | `mat_lapack_dorgqr`       |
| `zungqr_`  | form Q from complex Householder reflectors    | `mat_lapack_zungqr`       |
| `cblas_ddot` (BLAS) | Phase-1 probe used by smoke test    | `mathilda_lapack_probe`   |

### 5.5 Anti-patterns (do not do)

- **Don't vendor reference LAPACK source by default** — ~500K LoC of
  Fortran 77, drags in `gfortran`, minutes-slow first build. Only
  worth it as an opt-in flag for reproducibility-critical builds.
- **Don't hand-roll CBLAS in C** — you'll be 5–10× slower than
  OpenBLAS and own a permanent maintenance liability.
- **Don't hard-fail the build without LAPACK** — today Mathilda
  builds with one `make` on a fresh clone. Tier 4 auto-degrade
  preserves that property.
- **Don't depend on `<lapacke.h>` directly** — missing on Apple
  Accelerate, varies across Linux distros. Use the Fortran ABI.
- **Don't use Apple's deprecated CLAPACK with no-underscore symbols**
  — they were warned-then-removed in macOS 13.3. Use the current
  trailing-underscore Fortran ABI.

### 5.6 Beyond QR — follow-up opportunities

Once `USE_LAPACK` lands, the following can switch to LAPACK fast paths
in follow-up PRs (out of scope for this plan, but the wiring makes
them one-file changes each):

- `Inverse` → `dgetrf_` + `dgetri_`
- `LinearSolve` → `dgesv_`
- `Det` → `dgetrf_` (product of U diag with sign tracking)
- `Eigenvalues` real-machine path → `dgeev_` (replacing the
  hand-rolled Hessenberg+QR in `eigen_direct.c`)
- `LeastSquares` → `dgels_`
- `Norm` 2-norm singular-value path → `dgesvd_`
- `MatrixRank` numeric path → `dgesvd_`

## 6. Implementation phases

### Phase 1 — wiring (no behaviour change)
1. Add `USE_LAPACK` to `makefile` with four-tier detection ladder
   (Apple → pkg-config lapacke → conventional lapacke → auto-off).
2. Add matching `option(USE_LAPACK)` to `tests/CMakeLists.txt`.
3. Add `src/linalg/lapack.h` — platform-papering header. Defines
   `mathilda_lapack_probe()` and `extern` Fortran prototypes for
   the QR routines (declared, used in Phase 3).
4. Add `src/linalg/lapack.c` — `mathilda_lapack_probe()` calls
   `cblas_ddot` to confirm runtime linkage.
5. Add `tests/test_lapack.c` smoke test.
6. Verify build under `USE_LAPACK=1` and `USE_LAPACK=0`.
7. Run existing 70-QR test suite — must still pass unchanged.

### Phase 2 — refactor existing QR
1. Rename `mgs_core` → `qr_symbolic_core`, keep behaviour identical.
2. Extract dispatch from `builtin_qrdecomposition` into `qr_dispatch`.
3. Add `qrdecomp_internal.h` with `MatD` / `MatM` decls (re-exported
   from eigen for now).
4. All existing tests must still pass.

### Phase 3 — machine-precision kernel (`qrdecomp_machine.c`)
1. Load `m` into `MatD` (column-major). Detect real vs complex via
   `is_complex` flag.
2. If `Pivoting -> False`: call `dgeqrf` / `zgeqrf`.
3. If `Pivoting -> True`: call `dgeqp3` / `zgeqp3`.
4. Form Q via `dorgqr` / `zungqr` at rank `min(rows, cols)`.
5. Determine numerical rank from R diagonal (drop rows where
   `|R[i,i]| < tol * |R[0,0]|`) for the "thin" QR contract.
6. Convert back to Mathilda Lists (Q rows = `ConjugateTranspose[Q_lapack]`
   to match the public API `m == ConjugateTranspose[q].r`).
7. Apply `MachinePrecision` to all output reals.

### Phase 4 — MPFR Householder kernel (`qrdecomp_mpfr.c`)
1. Load `m` into `MatM` at working precision = `min_bits`.
2. For complex: pair of MPFR arrays `re`, `im`; arithmetic via the
   helpers already in `eigen_common.c`.
3. Householder reflectors over MPFR — chosen over MGS because it is
   numerically stable enough to publish without iterative
   refinement, matching what Mathematica does at MPFR precision.
4. Column pivoting (full, like `dgeqp3`) when `Pivoting -> True`.
5. Rank detection: `|R[i,i]| < 2^(-bits/2) * |R[0,0]|`.
6. Build output reals as MPFR-precision `Real[..., prec]` nodes.

### Phase 5 — tests (significant coverage required)

Test strategy is **dispatch ↔ kernel ↔ identity**: prove the
dispatcher picks the right kernel, prove each kernel produces a
correct factorisation, prove each kernel honours the public API
contract. Five new test binaries plus regression coverage:

#### 5.1 `tests/test_lapack.c` (Phase 1, smoke)
- `mathilda_lapack_probe()` returns 1 when `USE_LAPACK=1`
- Builds and links under `USE_LAPACK=0` (probe returns 0, no missing
  symbols)
- `cblas_ddot` of orthogonal unit vectors → 0; of parallel → 1

#### 5.2 `tests/test_qrdecomposition_machine.c` (Phase 3)
Coverage matrix, every cell must pass:

| Axis                 | Cases                                                |
|----------------------|------------------------------------------------------|
| Element type         | `MachinePrecision` real, `MachinePrecision` complex  |
| Shape                | 1×1, 2×2, 3×3, 4×4, 5×5, 10×10, 20×20, 50×50, 100×100 |
| Aspect ratio         | square, tall (3×2, 5×3, 20×10), wide (2×3, 3×5, 10×20) |
| Rank                 | full, rank-deficient (one duplicated column), zero matrix |
| Conditioning         | identity, well-conditioned random, Hilbert(6), Vandermonde(10) |
| Pivoting             | `Pivoting -> False`, `Pivoting -> True`               |

For each: assert `||m - q^H·r||_∞ ≤ 16·ε·||m||_∞` and
`||q·q^H − I||_∞ ≤ 16·ε`. For pivoted form additionally assert
`||m·p − q^H·r||_∞ ≤ 16·ε·||m||_∞` and `|r[i,i]|` is non-increasing.

Numerical-stability test: factor `H = Hilbert[6]` and check the
recovered `r[0,0]` matches the known true value to within
`κ(H)·ε`.

Cross-kernel agreement: round-trip an exact rational matrix through
both the machine kernel (after `N[.]`) and the symbolic kernel
(directly), `N[.]` the symbolic result, assert agreement to `16·ε`.

#### 5.3 `tests/test_qrdecomposition_mpfr.c` (Phase 4)
Same coverage matrix as machine, at precisions 30, 60, 100, 200 digits.

Per-precision tolerance: `2^(-bits + 16)` for reconstruction,
`2^(-bits/2)` for orthonormality (loss-of-orthogonality bound for
Householder is known to be O(2^(-bits/2)) on ill-conditioned input).

**Precision-monotonicity** test: same matrix, factor at 30/60/100/200
digits, assert reconstruction error decreases monotonically by ≥ a
factor of 2^(bits_step) between successive precisions.

#### 5.4 `tests/test_qrdecomposition_dispatch.c` (Phase 3 + 4)
- All-`Integer` input → symbolic kernel (no LAPACK call). Verify by
  setting a test-only counter in each kernel and asserting which
  one fired.
- All-`Real` input → machine kernel
- All-MPFR-100 input → MPFR kernel
- Mixed `Real` + MPFR-100 → MPFR kernel at 53 bits (min_bits policy
  from `common_scan_inexact`)
- Mixed `Real` + `Integer` → machine kernel
- Mixed MPFR-30 + MPFR-100 → MPFR kernel at 30 (min)
- `USE_LAPACK=0` build: machine input routes to MPFR kernel with a
  one-time warning, factorisation still correct

#### 5.5 `tests/test_qrdecomposition.c` (existing, regression)
All 70 existing assertions must pass unchanged. The symbolic path
is untouched, so this is a "did we accidentally break it" guard.

#### 5.6 `tests/bench_qr.c` (Phase 5 microbench, not a unit test)
- Print per-(precision × shape × pivoting) wall time.
- Reads `tests/data/qr_mathematica.json` (checked-in Mma baseline,
  same seed) and prints ratio.
- Asserts each ratio is ≤ the success-bar value from §7. **Hard-fail
  on regression**, so a future change that slows QR back to evaluator
  speed is caught in CI.

#### 5.7 Valgrind coverage
Each new test binary runs under `valgrind --leak-check=full
--error-exitcode=1`. Zero Mathilda-code leaks required.

#### 5.8 Fuzz pass
A small property-based check:
```
for i in 0..10000:
   random shape m×n in [1..6], random rank, random entries
   QR -> assert identities
```
Catches off-by-one and aspect-ratio bugs the curated suite misses.

### Phase 6 — valgrind & docs
1. `valgrind --leak-check=full` on focused QR binary — zero leaks
   from Mathilda code.
2. Update `docs/spec/builtins/linear-algebra.md` with the new
   precision-dispatch behaviour.
3. Add changelog entry in `docs/spec/changelog/2026-05.md`.

## 7. Test corpus & success criteria

Concrete pass bar before merging:

| Case               | Mma t  | Mathilda before | Mathilda after target |
|--------------------|--------|-----------------|-----------------------|
| machine real n=20  |   8 µs |       (>10 min) |          ≤ 100 µs     |
| machine real n=100 |   ~ms  |       (intractable) |       ≤ 5 ms      |
| machine cplx n=20  |  18 µs |       (>10 min) |          ≤ 200 µs     |
| MPFR-30 n=20       | 2.7 ms |       (>10 min) |          ≤ 10 ms      |
| MPFR-100 n=20      | 2.3 ms |       (>10 min) |          ≤ 30 ms      |
| symbolic int n=4   | 696 µs |          858 µs |          ≤ 1 ms       |

Functional pass bar:
- Reconstruction error `||m - q^H.r||_∞ / ||m||_∞` ≤ `10·eps` (machine)
  or `10·2^(-bits)` (MPFR).
- Orthonormality `||q.q^H - I||_∞` to the same tolerance.
- Pivoted form: `||m.p - q^H.r||_∞` to the same tolerance, and the
  permutation chooses columns in descending residual-norm order.
- Rank-deficient input: `Length[q] == Length[r] == numerical_rank`.

## 8. Risks & mitigations

| Risk                                          | Mitigation                                                                 |
|-----------------------------------------------|----------------------------------------------------------------------------|
| LAPACK on Linux needs `lapacke.h` (row-major API) which isn't always packaged | Detect at configure; fall back to column-major + manual transpose if `lapacke.h` missing |
| Accelerate's CLAPACK signatures take pointers everywhere (Fortran ABI) | Wrap behind `src/linalg/lapack.c` once, hide the noise from kernels        |
| Widening `MatD`/`MatM` could break eigen tests | One-line audit per eigen callsite (all already pass `n=n`); land widening behind a single PR with eigen tests in CI |
| MPFR Householder rounding vs. Mathematica's choice of algorithm | Both use Householder; tolerance check is "close enough", not bit-equal     |
| `USE_LAPACK=0` builds for environments without LAPACK | Stub falls back to existing symbolic kernel with a one-time warning, like `USE_MPFR=0` |
| Complex via paired MPFR is awkward for inner products | Already proven in `eigen_common.c`; reuse those primitives                 |

## 9. Out of scope (future work)

- Switching `Inverse` / `LinearSolve` / `Det` / `Eigenvalues` /
  `LeastSquares` to LAPACK — done as separate, follow-up PRs once
  `USE_LAPACK` lands.
- SVD-based `MatrixRank` / `PseudoInverse` numerical path.
- Iterative refinement after LAPACK QR (Mathematica does this for
  `WorkingPrecision -> MachinePrecision` rank-deficient input).
- MPC-based complex MPFR kernel (current pair-of-arrays approach is
  cheaper and already idiomatic in this codebase).

## 10. Checklist

- [x] Phase 1 — `USE_LAPACK` build flag (makefile + CMake)
- [x] Phase 1 — `src/linalg/lapack.h` + `lapack.c` wrappers
- [x] Phase 1 — `tests/test_lapack.c` smoke test (8 assertions)
- [x] Phase 1 — `USE_LAPACK=0` graceful-degrade build verified
- [x] Phase 1 — 70-QR + 90 linalg + 75 eigen + 39 matinv + 33 linsolve regression clean
- [x] Phase 2 — `qrdecomp_internal.h` private header (QrOpts, kernel decls)
- [x] Phase 2 — `mgs_core` renamed to `qr_symbolic_core`; `qr_symbolic_dispatch` + `qr_dispatch` extracted
- [x] Phase 2 — `builtin_qrdecomposition` reduced to opts + shape check + dispatch
- [x] Phase 2 — `USE_LAPACK=1` + `USE_LAPACK=0` regression sweep clean (lapack + QR + linalg + matinv + linsolve + eigen + matrank + ConjugateTranspose)
- [x] Phase 3 — Machine kernel (real + complex, pivoting + no-pivot)
- [x] Phase 3 — `mat_lapack_d/zgeqrf` + `d/zgeqp3` + `d/zungqr` Fortran-ABI wrappers in `lapack.{c,h}`
- [x] Phase 3 — `qrdecomp_machine.c` + `qr_machine_dispatch` wired into `qr_dispatch`
- [x] Phase 3 — Reconstruction + orthonormality bounds, rank-deficient + zero handling, pivoting monotonicity, cross-kernel agreement
- [x] Phase 3 — `tests/test_qrdecomposition_machine.c` (18 tests, 32 pass lines); skip-clean under `USE_LAPACK=0`
- [x] Phase 3 — Full linalg regression sweep clean (lapack + QR + linalg + matsol_methods + matinv + matinv_methods + matrank + eigen + nullspace + conjugate_transpose + mateigen_{direct,arnoldi,banded})
- [x] Phase 3 — docs/spec/builtins/linear-algebra.md + docs/spec/changelog/2026-05.md updated
- [x] Phase 4 — MPFR Householder kernel (`qrdecomp_mpfr.c`, real + complex, pivoting)
- [x] Phase 4 — Column-swap mirror update on already-stored R rows (LAPACK gets this for free via in-place upper-triangle storage)
- [x] Phase 4 — `qr_mpfr_dispatch` wired into `qr_dispatch` (`min_bits > 53` route, NULL fall-through to symbolic)
- [x] Phase 4 — `tests/test_qrdecomposition_mpfr.c` (17 tests, precision 30/60/100/200, plus precision-monotonicity)
- [x] Phase 4 — Full linalg regression sweep clean (15 binaries)
- [x] Phase 5 — `tests/bench_qr.c` microbench (per-shape × precision × pivoting wall time, no JSON baseline yet — that's a separate Mma-side plumbing task)
- [x] Phase 6 — Valgrind clean on `qrdecomposition_mpfr_tests` (zero Mathilda-code leaks)
- [x] Phase 6 — `docs/spec/builtins/linear-algebra.md` + `docs/spec/changelog/2026-05.md` updated for the MPFR kernel

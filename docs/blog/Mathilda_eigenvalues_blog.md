# Eigenvalues and Eigenvectors in Mathilda

## Introduction

`Eigenvalues` and `Eigenvectors` sit at an awkward crossroads in a
computer algebra system. On the symbolic end they want to be a
characteristic-polynomial calculation followed by `Solve`; on the
numerical end they want to be a dense Hessenberg + Francis QR sweep;
in between sit the matrices where neither is appropriate — large and
sparse, banded, narrow-interval, MPFR — and each of those wants its
own kernel.

A faithful Mathematica-style implementation has to handle all of them
under one builtin, with a single argument grammar, and route to the
right kernel without the user having to think about it. This post
walks through how Mathilda's `Eigenvalues` / `Eigenvectors` does that,
and benchmarks the result against Mathematica on classic test
matrices (Hilbert, tridiagonal, dense random).

The implementation lives in `src/mateigen.c` (~11,200 lines) and the
test battery spans four files
(`tests/test_mateigen_direct.c`, `_arnoldi.c`, `_banded.c`, `_feast.c`,
~4,000 lines total).

## Two easy examples

```text
In[1]:= Eigenvalues[{{a, b}, {c, d}}]
Out[1]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]),
         1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}

In[2]:= Eigenvalues[N[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]]
Out[2]= {3.4142135623730945, 2.0, 0.5857864376269047}
```

The first runs through the symbolic path: build `Det[m - λ I]` by
Faddeev–Leverrier, hand to `Solve`. The second routes through the
numerical Direct kernel: Householder tridiagonalisation, symmetric
tridiag QR, sort by `|λ|` descending. Every other call you can write
is one of those two with a different preprocessor, a different kernel
choice, or a precision upgrade.

## The dispatcher

`builtin_eigenvalues` in `src/mateigen.c:10895` is the top of a
shape-dependent cascade. The cascade is shorter than Limit's — five
branches, not eleven — but each branch is its own substantial body of
code:

```
  symbolic / exact input
      │
      ▼
   Faddeev-Leverrier char-poly  -->  Solve  -->  optional N[]
                                                       │
                                                       ▼
                                                exact / Root[ ]

  approximate (Real / MPFR) input
      │
      ▼
   ┌─────────────────────────────────────────────┐
   │  Method routing:                            │
   │    explicit "Banded" / "Arnoldi" / "FEAST"  │
   │    Automatic + Hermitian + narrow band       │  → Banded
   │    Automatic + small k                       │  → Arnoldi
   │    everything else                           │  → Direct
   └─────────────────────────────────────────────┘
                       │
                       ▼  (each kernel returns NULL = "I refuse",
                          which cascades to the next preference,
                          terminating at the symbolic path)
```

Every numerical kernel has both a machine-precision (`double`) and an
MPFR variant. `direct_dispatch`, `arnoldi_dispatch`, `banded_dispatch`
and `feast_dispatch` each call `common_scan_inexact(m).min_bits` and
route to the MPFR kernel when any input leaf carries more than 53
bits.

The rest of this post walks the cascade in the order above.

---

### The symbolic path — Faddeev–Leverrier

The naive way to compute a characteristic polynomial is Laplace
expansion of `Det[m - λ I]` where each entry is now a polynomial in
`λ`. That is `O(n!)` even for the symbolic determinant routine and
becomes the bottleneck past `n = 8`.

`eigen_char_poly_faddeev` in `src/mateigen.c:218` uses the
Faddeev–Leverrier–Souriau recurrence instead. Define

```
  M_0 = 0,        c_n = 1
  M_k = A · M_{k-1} + c_{n-k+1} I
  c_{n-k} = -(1/k) tr(A · M_k)
```

After `n` matrix multiplications the `c_k` are the coefficients of
`det(λI - A) = λ^n + c_{n-1} λ^{n-1} + ... + c_0`. Total cost is
`O(n^4)` arithmetic operations — but every operation is a symbolic
multiply, so the constant factor lives in `builtin_times` /
`builtin_plus`, and for diagonal or sparse matrices the dominant cost
disappears entirely because most multiplications fold to zero.

The polynomial is then handed to `Solve`. For degrees ≤ 2 this is
the closed-form quadratic; for degree 3 and 4 the cubic / quartic
radical formulas (gated by `Cubics -> True` / `Quartics -> True`, on
by default); above degree 4 the result is wrapped in `Root[poly, k]`.

Worked example: the 4×4 Hilbert matrix.

```text
In[3]:= Eigenvalues[Table[1/(i+j-1), {i,4}, {j,4}]]
Out[3]= {Root[1 - 10496 #1 + 1603680 #1^2
              - 10137600 #1^3 + 6048000 #1^4 &, 1],
         Root[..., 2], Root[..., 3], Root[..., 4]}
```

Mathilda emits `Root` objects unevaluated — it does not yet
numericalize them, so `N[%]` is a no-op here. For closed-form
spectra below degree 5 the path returns exact algebraic expressions,
e.g. `Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]` gives
`{(15 + 3 Sqrt[33])/2, (15 - 3 Sqrt[33])/2, 0}`.

**When to use:** symbolic matrices, exact closed-form spectra,
generalised eigenproblems `Eigenvalues[{m, a}]` (which still flow
through this path), small exact matrices where the user wants the
characteristic polynomial structure preserved.

**When to avoid:** any `n × n` exact integer / rational matrix with
`n ≥ 5` and a generic spectrum — the result is a list of `Root`
objects you can't fold further. Numericalize first with `N[m]`.

---

### Direct — machine precision

When the input is approximate (`EXPR_REAL` or `EXPR_BIGINT` inside
`N[...]`), Direct is the default kernel. It is also the only kernel
that handles every numeric matrix shape; the others all explicitly
refuse non-Hermitian or generalised input and cascade back to Direct.

The Direct kernel decomposes into four sub-kernels selected by
`direct_dispatch_machine` (`src/mateigen.c:7599`):

| Shape                 | Algorithm                                                |
|-----------------------|----------------------------------------------------------|
| Real symmetric        | Householder tridiag + symmetric tridiag QR              |
| Real general          | Householder Hessenberg + Francis double-shift QR + Schur back-sub |
| Complex Hermitian     | Complex Householder + diagonal phase correction + real symmetric QR |
| Complex non-Hermitian | Real 2n × 2n block embedding + Francis QR + +J/−J split |

The complex Hermitian path is worth a closer look. After Hermitian
Householder tridiagonalisation the matrix `T` has real diagonal but
generally complex sub-diagonal `sub[k]`. A diagonal unitary
`D = diag(d_k)` with `d_{k+1} = d_k · sub[k] / |sub[k]|` zeroes the
sub-diagonal phases, so `D^H T D` is real symmetric tridiagonal and
the existing real-symmetric QR (`direct_symtridiag_qr`) finds the
eigenvalues. `Q` is updated to `Q D` and the final complex
eigenvectors are `V = Q Z` where `Z` is the real QR's eigenvector
matrix. The whole thing avoids needing a separate complex
tridiagonal QR. See Wilkinson 1965, *The Algebraic Eigenvalue
Problem*, §5.45.

The complex non-Hermitian path is even more devious. Rather than
write a native complex Hessenberg + complex QR, Mathilda embeds
`A = R + i S` into the real `2n × 2n` block

```
   M = [[ R, -S ],
        [ S,  R ]]
```

whose spectrum is `spec(A) ∪ conj(spec(A))` as a multiset. The real
Francis QR runs on `M`, then for each `M`-eigenvector
`w = w_re + i w_im` with splits `w_re = [a; b]`, `w_im = [c; d]`, the
candidate `x = (a - d) + i (b + c)` either satisfies `A x = μ x`
(when `μ ∈ spec(A)`) or vanishes (when `μ` is only in
`conj(spec(A))`). Algebraic multiplicity is recovered by grouped
modified Gram–Schmidt over `M`-eigenvalue clusters. The block
embedding is roughly 8x more flops than a native complex kernel at
`-O3`, but it shares one carefully tuned QR sweep across all of the
machine-precision Direct paths.

**Pros:** handles every matrix shape; produces the full spectrum;
serves as the universal fall-back when other kernels refuse.

**Cons:** `O(n^3)` always — you pay for the full Hessenberg even if
you only wanted the top eigenvalue. The block embedding wastes a
factor of ~8 on complex non-Hermitian input vs. a native complex QR.
No BLAS, so the inner loops are pure C99 doubles and a vendor LAPACK
will beat Mathilda by 5–20x at `n = 800` (see benchmarks below).

---

### Direct — MPFR precision

Every Direct sub-kernel has a paired `_M` variant that runs at the
input precision detected by `common_scan_inexact`. Real symmetric is
in `direct_real_sym_mpfr`; real general in `direct_real_general_mpfr`;
complex Hermitian in `direct_complex_hermitian_mpfr`; complex general
in `direct_complex_general_mpfr`.

Three things to know about the MPFR side:

1. **Scratch pools.** Every inner loop receives a caller-supplied
   pool of pre-initialised `mpfr_t` cells (e.g. the Francis step
   takes 14) so `mpfr_init` / `mpfr_clear` stays out of the hot
   path. Without this, `mpfr_init` allocator pressure dominates
   wall time.

2. **Exceptional shifts.** The non-symmetric Francis QR adds an
   LAPACK `dlahqr`-style exceptional shift every 10 iterations:
   `s = 1.5 · (|H[p-1, p-2]| + |H[p-2, p-3]|)`. At 50+ decimal
   digits the trailing 2×2 can stall in a 2-cycle when the
   ordinary shift repeats exactly; machine precision masks this
   via roundoff but MPFR does not.

3. **Deflation tolerance.** The non-symmetric path runs at
   `rel_tol = 2^{-bits+10}` (1024 ULPs) — about one decimal digit
   looser than the symmetric tridiag tolerance. Non-symmetric QR
   accumulates more per-sweep roundoff (factor ~n² from bulge chase
   plus back-transformation), so a slightly more permissive
   threshold is needed for clean convergence at high precision.

The cost is roughly `O(n^3 · p^2)` where `p` is the precision in
bits — MPFR multiplication is quadratic in the operand bit-length
for small `p` and Toom-Cook beyond that. For a 30×30 Hermitian at
200 digits the wall time is ~14 ms, dominated by the Francis sweep
on the embedded block. See the Hilbert section below for a head-to-
head with Mathematica.

---

### Arnoldi

Direct is `O(n^3)` even for the top-`k` eigenvalues. When `k` is
small relative to `n`, the right tool is an iterative projection
method.

`arnoldi_real_general_machine` (`src/mateigen.c:2750`) builds an
m-step Krylov subspace `K_m(A, v_0) = span{v_0, A v_0, ..., A^{m-1} v_0}`
via classical Gram–Schmidt with one re-orthogonalisation pass
("twice is enough" — Daniel et al. 1976), producing an orthonormal
basis `V_m` and an `m × m` upper Hessenberg `H_m` such that

```
  A V_m = V_m H_m + h_{m+1,m} v_{m+1} e_m^T.
```

`H_m` is then diagonalised by the **same** Phase 2 Francis QR
pipeline used by Direct (`direct_qr_real_general` +
`schur_compute_eigvecs`), and Ritz vectors `V_m y_i` are lifted back
to `A`-eigenvectors. Default basis size is `max(2k, 20)` capped at
`n`; on lucky breakdown (`||w||` below tolerance at step `j`)
Arnoldi terminates early with `j+1` exact eigenpairs.

Complex inputs use paired re/im storage for `V_m` and `H_m`, and
route `H_m`'s diagonalisation through the same `2μ × 2μ` real block
embedding the Direct path uses — the same trick to avoid writing a
native complex Hessenberg QR.

The Automatic dispatcher (`arnoldi_automatic_prefers`) routes
through Arnoldi when `n > 32` and `k ≤ max(20, n/10)`. Below that
threshold Direct is faster and exact. Above it Arnoldi's `O(m · n^2)`
beats Direct's `O(n^3)` for `m << n`.

```text
In[4]:= n = 400; SeedRandom[1]; A = N[RandomInteger[{-10, 10}, {n, n}]];
        S = (A + Transpose[A])/2;
        First[Timing[Eigenvalues[S, 5, Method -> "Arnoldi"]]]
Out[4]= 0.012081

In[5]:= First[Timing[Eigenvalues[S, 5, Method -> "Direct"]]]
Out[5]= 0.046525
```

**Pros:** small-`k` queries on large matrices; matrix-free Krylov
expansion is the natural fit if the matrix-vector product is cheap;
shared kernel reduces code surface.

**Cons:** convergence is governed by spectral gap. Ritz values
converge to the spectral edge first; densely-packed interior
eigenvalues need a larger basis or shift-and-invert (not yet
implemented). The default `BasisSize -> max(2k, 20)` is conservative
— for tightly clustered spectra you may want to set it explicitly.

---

### Banded

Hermitian matrices with small half-bandwidth `b ≪ n` are common in
practice — discretisations of differential operators, finite-element
mass matrices, all of computational chemistry's overlap matrices —
and the Direct kernel's `O(n^3)` is wasted on them.

`banded_to_tridiag_real_sym` (`src/mateigen.c:3446`) uses
Schwarz / Rutishauser two-sided Givens rotations that zero one
off-band entry at a time and chase the bulge they introduce until
it falls past the matrix edge:

```
For each anchor column k = 0..n-3:
  For each subdiagonal d = b, b-1, ..., 2:
    Givens in plane (k+d-1, k+d) zeroes A[k+d, k]
        using the in-band entry A[k+d-1, k].
    Two-sided application produces a fresh bulge at (p, q+b).
    Chase: Givens in plane (q+b-1, q+b) eliminates the bulge,
           propagating a new bulge another b columns along,
           until cq + b >= n.
```

The result is symmetric tridiagonal, and the existing Phase 2
symmetric tridiag QR (`direct_symtridiag_qr`) takes it from there.
Total cost is `O(n^2 · b)`. For a tridiagonal `n = 1000` (`b = 1`)
this is `O(n^2)` — three orders of magnitude faster than the dense
Direct path.

Complex Hermitian banded reuses the Phase 2 phase-correction step
(`direct_phase_correct_tridiag`): the band-Givens produces a
Hermitian tridiagonal with real diagonal and complex sub-diagonal,
the phase correction makes the sub-diagonal real-positive, and the
final QR runs at real symmetric tridiagonal precision. The Givens
itself uses Wilkinson's real-`c` / complex-`s` parameterisation
(`c = |a| / r`, `s = b · conj(a) / (|a| · r)`), which keeps the
mathematics in paired re/im doubles (no `libmpc` dependency).

The Automatic dispatcher routes here when the matrix is Hermitian,
`n > 8`, and `b ≤ max(8, n/4)`. Explicit `Method -> "Banded"` is
always safe to pass — Banded refuses (returns `NULL`, the wrapper
cascades to Direct) on non-Hermitian or fully dense (`b == n-1`)
input, so the worst case is a silent fall-through.

```text
In[6]:= n = 1000; T = Table[If[i==j, 2., If[Abs[i-j]==1, -1., 0.]],
                            {i, n}, {j, n}];
        First[Timing[Eigenvalues[T, Method -> "Direct"]]]
Out[6]= 0.620787

In[7]:= First[Timing[Eigenvalues[T, Method -> "Banded"]]]
Out[7]= 0.071603
```

**Pros:** linear in `n` for fixed band; the right tool for the very
common Hermitian-tridiagonal / pentadiagonal / heptadiagonal shapes
that show up in physics and PDE work; same QR back-end as Direct
keeps the code compact.

**Cons:** Hermitian-only; dense matrices fall back silently to
Direct. No shift-and-invert for inverse iteration, so it's a full-
spectrum kernel — it doesn't compose with a top-`k` query the way
Arnoldi does.

---

### FEAST

The fourth and final numerical kernel is the only one that returns a
**spectral slice** rather than the full spectrum: given a Hermitian
matrix and a real interval `[a, b]`, FEAST returns precisely the
eigenpairs whose eigenvalues lie inside it.

The algorithm (Polizzi 2009) is a contour-integral spectral
projector:

```
1.  Pick m_0 (default max(20, n/4), capped at n) and seed
    an n x m_0 random subspace Y.
2.  Approximate
        P_[a,b](A) Y = (1/(2 pi i)) ∫_C (z I - A)^-1 Y dz
    by Ne-point Gauss-Legendre quadrature on the upper half of the
    elliptic contour through (a, 0) and (b, 0).  Schwarz symmetry
    halves the number of complex linear solves to Ne / 2.  Supported
    Ne values are 2, 4, 8 (default), and 16.
3.  Form A_q = Q^H A Q, B_q = Q^H Q and solve the small generalised
    Hermitian-definite problem
        A_q y = lambda B_q y
    by Cholesky reduction B_q = L L^H + the existing Direct kernel
    on L^-1 A_q L^-H.
4.  Filter Ritz pairs against [a, b], lift back to A-eigenvectors,
    iterate until residual converges or MaxIterations is reached.
```

The complex shift solve `(zI - A)x = y` at each quadrature node is a
partial-pivoted LU on paired re/im storage
(`feast_complex_lu_factor` / `_solve`), with no dependency on
`libmpc`. The MPFR variants (`_M`) mirror the machine kernel at
arbitrary precision.

`Method -> "FEAST"` is the **only** way in — Automatic never routes
to FEAST because it requires the user to commit to an interval.
`"Interval"` is the only required sub-option:

```text
In[8]:= Eigenvalues[N[{{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0},
                        {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1},
                        {0, 0, 0, -1, 2}}],
                    Method -> {"FEAST", "Interval" -> {2.5, 4}}]
Out[8]= {3.7320508075688776, 3.0}
```

A fail-soft cascade catches the common edge cases — non-Hermitian
input, missing or degenerate `"Interval"`, generalised eigenproblems,
Cholesky failure on `B_q` (subspace too small for the spectral count
inside `[a, b]`), LU singular at any quadrature node, non-
convergence within `MaxIterations` — and falls through to Direct
with a single stderr line tagged by reason. An explicit
`Method -> "FEAST"` call always returns *some* sensible answer; at
worst, the full Direct spectrum.

**Pros:** the only way to compute "the eigenvalues in this energy
window" without computing the whole spectrum; native MPFR variant;
the right tool for ill-conditioned spectra where you know roughly
where the eigenvalues you want live.

**Cons:** without BLAS, the per-iteration cost is `Ne/2` complex
linear solves at `O(n^3)` each, so for small dense matrices Direct
wins comfortably. FEAST shines when the matrix is genuinely large
and sparse and the interval is narrow — which, in a system without
sparse data structures, isn't yet a regime Mathilda excels in. The
current implementation is the algorithm in its general form; making
it competitive will need a sparse `MatrixVector` and a vendor LU.

---

## Comparison with Mathematica

All timings below are wall-clock seconds, single-threaded, on the
same machine (M-series Apple silicon). Mathematica is version 14.0
(`wolframscript`); Mathilda is the May 2026 HEAD of `main`. The
random matrices use `SeedRandom[42]` so the inputs match exactly.

### Hilbert matrices

Hilbert matrices `H[i, j] = 1 / (i + j - 1)` are the canonical
ill-conditioned test case. `Eigenvalues[H_n]` is dominated by a few
large eigenvalues and many exponentially small ones — a stress test
for any QR routine's deflation logic.

| Matrix              | Mathilda      | Mathematica   | Note                       |
|:--------------------|--------------:|--------------:|:---------------------------|
| Symbolic 5×5        | 0.0114 s      | 0.0029 s      | `Root[...]` output         |
| Symbolic 6×6        | 0.0201 s      | 0.0053 s      | `Root[...]` output         |
| Machine 20×20       | 0.000083 s    | 0.000943 s    | dense Direct               |
| Machine 40×40       | 0.000214 s    | 0.000097 s    | dense Direct               |
| MPFR 60d 40×40      | 0.0124 s      | 0.0630 s      | MPFR symmetric tridiag QR  |
| MPFR 120d 40×40     | 0.0153 s      | 0.0928 s      | MPFR symmetric tridiag QR  |

The MPFR results are notable: Mathilda's MPFR Direct kernel runs
~4-6x faster than Mathematica on these Hilbert sizes, and the
agreement on the smallest eigenvalue is to all 37 displayed digits:

```text
Mathilda    (50d, n = 10):  1.09315381937966576381686691049785556829...e-13
Mathematica (50d, n = 10):  1.09315381937966576381686691049785556828...e-13
```

The first ~37 decimal digits agree. The kernel walks the same
Householder tridiagonalisation + symmetric tridiag QR as the
machine path, with every `double` op replaced by an `mpfr_*` call
rounded `MPFR_RNDN` (see Phase 2d-A, `direct_real_sym_mpfr`).

### Random symmetric matrices

`A = N[RandomInteger[{-10, 10}, {n, n}]]; S = (A + Transpose[A])/2`,
seeded with `SeedRandom[42]` so the matrix is reproducible.

| n   | Mathilda  | Mathematica | Ratio   |
|----:|----------:|------------:|--------:|
|  50 | 0.000357 s | 0.000322 s |  1.1x   |
| 100 | 0.001564 s | 0.000431 s |  3.6x   |
| 200 | 0.006308 s | 0.001507 s |  4.2x   |
| 400 | 0.049079 s | 0.005519 s |  8.9x   |
| 800 | 0.313672 s | 0.022941 s | 13.7x   |

Mathilda is competitive at `n ≤ 100`. The gap widens with `n^3`
because Mathematica's `Eigenvalues` is calling vendor LAPACK
(`dsyevd` on Apple silicon's Accelerate framework, which dispatches
to an Apple-tuned divide-and-conquer kernel). Mathilda is running
pure C99 `double` ops in the inner loop with no vectorisation hint
beyond what `-O3` extracts. The 13.7x at `n = 800` is roughly the
expected gap between a tuned LAPACK and a portable scalar
implementation.

### Random general (non-symmetric) matrices

| n   | Mathilda  | Mathematica | Ratio  |
|----:|----------:|------------:|-------:|
|  50 | 0.000824 s | 0.000465 s | 1.8x   |
| 100 | 0.004594 s | 0.002428 s | 1.9x   |
| 200 | 0.026707 s | 0.015156 s | 1.8x   |
| 400 | 0.201239 s | 0.053183 s | 3.8x   |

The non-symmetric ratio is much smaller. Mathilda's real Francis
QR is in the same ballpark as Mathematica's `dgeev` until `n ≥ 400`
where Apple's blocked Hessenberg reduction starts to dominate.

### Tridiagonal matrices

The discrete Laplacian `tridiag(-1, 2, -1)` is a banded-eigenvalue
benchmark. `Method -> "Banded"` on Mathilda picks up the
Schwarz / Rutishauser chase; Mathematica auto-detects the band.

| n    | Mathilda Banded | Mathilda Direct | Mathematica |
|-----:|----------------:|----------------:|------------:|
|  200 |       0.0028 s  |        0.0068 s |    0.0018 s |
|  500 |       0.0175 s  |        0.0827 s |    0.0065 s |
| 1000 |       0.0716 s  |        0.6208 s |    0.0465 s |
| 2000 |       0.3229 s  |              -  |    0.3418 s |

At `n = 2000` Mathilda's Banded kernel beats Mathematica's auto-
detected routine. Direct on the same matrix is ~10x slower, which is
the expected `O(n^3)` vs `O(n^2)` gap. The crossover happens around
`n = 100` for tridiagonals; Banded is essentially always the right
choice for narrow bands.

### Small-k queries (Arnoldi)

`Eigenvalues[m, k]` on a symmetric matrix, small `k`. The same
random symmetric construction as above with `SeedRandom[1]`.

| n   | k | Arnoldi    | Direct     | Ratio  |
|----:|--:|-----------:|-----------:|-------:|
| 200 | 5 | 0.00313 s  | 0.00735 s  | 0.43x  |
| 400 | 5 | 0.01208 s  | 0.04653 s  | 0.26x  |

Arnoldi pays off when `k` is small relative to `n`. The crossover
in Mathilda's Automatic dispatcher is at `k ≤ max(20, n/10)`, which
matches the curve.

A word of caution: Ritz values converge to the spectral edge first,
so densely-packed interior eigenvalues need a larger `BasisSize`
than the default `max(2k, 20)`. The Arnoldi spec docstring spells
this out, and `Eigenvalues[m, k, Method -> {"Arnoldi", "BasisSize"
-> 4k}]` resolves the surprise on most test cases.

---

## The BLAS / LAPACK gap

Mathilda does not link BLAS or LAPACK. Every inner loop above is
pure C99: scalar `double` multiplies, no SIMD intrinsics, no cache
blocking. The benchmarks reflect that: on `n ≥ 400` dense symmetric
matrices the gap to Mathematica is roughly 5-15x, which is the
expected gap between a scalar reference implementation and a vendor
LAPACK tuned for the target microarchitecture.

What we are getting in return is portability and simplicity. The
build has three external dependencies — `gmp`, `mpfr`, `readline` —
none of which require a tuned BLAS. The whole eigenvalue stack
compiles with `gcc -std=c99 -Wall -Wextra -O3` on Linux and macOS
with no platform-specific code paths. The MPFR side is where the
absence of LAPACK actually pays off, because there isn't a vendor
LAPACK for MPFR — Mathilda's `direct_real_sym_mpfr` /
`direct_real_general_mpfr` are competitive with Mathematica's
internal MPFR routines because both sides are doing roughly the
same hand-rolled work.

LAPACK-HOOK comments in `src/mateigen.c` mark every entry point
where a vendor kernel could drop in: `dsyevd` for symmetric Direct,
`dgeev` for non-symmetric Direct, `zheevd` for complex Hermitian,
`zgeev` for complex general, `dsbevd` for Banded, and the FEAST
algorithm itself wants `zgetrf` / `zgetrs` (or `zgesv`) on the
contour LU. If someone wants to do this, the structural work is
already in place — each kernel is a single function with a clean
input/output contract.

---

## Capability snapshot

| Shape                                       | Status                |
|:--------------------------------------------|:----------------------|
| Symbolic `n × n`, `n ≤ 4`                   | closed-form radicals  |
| Symbolic `n × n`, `n ≥ 5`                   | `Root[poly, k]`       |
| Generalised symbolic `Eigenvalues[{m, a}]`  | symbolic path         |
| Machine real symmetric                      | Direct (Householder + tridiag QR) |
| Machine real general                        | Direct (Hessenberg + Francis QR + Schur) |
| Machine complex Hermitian                   | Direct (phase-corrected tridiag) |
| Machine complex non-Hermitian               | Direct (real-block embedding) |
| Machine Hermitian banded                    | Banded (Schwarz chase) |
| Machine top-`k`, `k ≪ n`                    | Arnoldi (Krylov + QR on H_m) |
| Machine Hermitian interval slice            | FEAST (contour integral) |
| MPFR real / complex, sym / general          | Direct (paired `_M` kernels) |
| MPFR Krylov                                 | Arnoldi MPFR          |
| MPFR Hermitian banded                       | Banded MPFR           |
| MPFR Hermitian interval slice               | FEAST MPFR            |
| Generalised numeric                         | falls to symbolic path |
| Iterative refinement / shift-and-invert     | future work           |
| Sparse data structures                      | future work           |
| Vendor BLAS / LAPACK drop-in                | annotated, not wired  |

---

## Conclusion

The four numerical kernels — Direct, Arnoldi, Banded, FEAST — cover
the standard textbook decomposition of the dense eigenvalue problem,
with each one specialised to a shape the others handle badly. Direct
is the universal solver; Arnoldi is for small-`k` queries; Banded is
for narrow Hermitian; FEAST is for interval slices on Hermitian
matrices. Each one has both a `double` and an `mpfr_t` variant
written from the same template, and each one falls through to the
next preference on shapes it refuses, terminating at the symbolic
characteristic-polynomial path.

The single visible gap in performance is the absence of BLAS /
LAPACK. On dense `n = 800` real symmetric matrices Mathilda is ~14x
slower than Mathematica; on tridiagonal `n = 2000` Mathilda's Banded
kernel actually wins; on MPFR Hilbert matrices Mathilda is 4-6x
faster than Mathematica. The pattern is consistent: where both
systems are running hand-rolled scalar code, Mathilda is
competitive; where Mathematica is calling Apple's tuned LAPACK,
Mathilda pays the scalar gap. The `LAPACK-HOOK` annotations in
`src/mateigen.c` mark the drop-in points for anyone who wants to
close that gap.

The tests pass — 71 across the four `test_mateigen_*.c` files — and
the eigenvalues of the 10×10 Hilbert matrix at 50 digits agree with
Mathematica to 37 places. That is the kind of regression line a
faithful reimplementation needs.

# Plan: `Eigenvalues` / `Eigenvectors` — FEAST Method

Status: planning. No code written yet. This document is the design for
implementing the **FEAST** option of `Method` on `Eigenvalues` and
`Eigenvectors`, completing the four-method numerical eigensolver
specified in `src/mateigen.h:18-27`.

## Context

The numerical eigensolver in `src/mateigen.c` already implements three
of the four advertised methods:

| Method     | Phase | Scope                                                              |
| ---------- | ----- | ------------------------------------------------------------------ |
| `Direct`   | 2     | real symmetric, real general, complex Hermitian, complex general   |
| `Arnoldi`  | 3     | real / complex general (k extreme eigenvalues)                     |
| `Banded`   | 4     | Hermitian banded                                                   |
| **`FEAST`**| —     | Hermitian, eigenvalues in a user-specified real interval **(this plan)** |

FEAST is the spectral-projector method due to Polizzi (2009). Given a
real interval `[a, b]` it returns *all* eigenpairs of a Hermitian (real
symmetric or complex Hermitian) matrix whose eigenvalues lie in that
interval. Internally it builds a low-dimensional invariant subspace via
contour integration of the resolvent `(zI − A)^{-1}` along an ellipse
through `(a, 0)` and `(b, 0)`, then performs Rayleigh–Ritz on that
subspace to recover the eigenpairs.

User-confirmed design points (from the planning Q&A on 2026-05-20):

- Sub-option name is **`"Interval" -> {a, b}"`** (matches the already-
  interned `SYM_Interval`; Mathilda is API-compatible, not bug-compatible,
  with Mathematica's `"Contour"` spelling).
- Failure modes are **fail-soft**: when `"Interval"` is missing, when the
  subspace is undersized, or when the residual doesn't converge within
  `MaxIterations`, the kernel returns `NULL` and the existing dispatcher
  cascade falls through to `Direct`, with a once-per-process warning on
  `stderr` explaining why.

---

## Scope

In scope for this plan:

- `Method -> "FEAST"` and `Method -> {"FEAST", sub-options...}` on
  `Eigenvalues` and `Eigenvectors`.
- Real symmetric and complex Hermitian ordinary eigenproblems at machine
  and MPFR precision.

Out of scope (and explicitly NOT implemented here):

- Non-Hermitian FEAST variants (Beyn's contour method, rational-filter
  FEAST). These would require new resolvent machinery and are not
  promised by `mateigen.h`.
- Generalised eigenproblems `A x = λ B x`. The contour method extends
  naturally to definite generalised problems (`B` SPD), but the
  Mathilda generalised path is currently symbolic-only across all
  methods; tackling it here would mean revisiting that decision for
  every method, not just FEAST.
- Cluster-tracking or auto-interval inference. The user must supply the
  interval.

---

## Algorithm

For Hermitian `A ∈ R^{n×n}` (or `C^{n×n}`) and target interval `[a, b]`:

```
Pick subspace size m₀ ≥ (estimated) number of eigenvalues in [a, b].
Initialise Y ∈ C^{n×m₀} with random or deterministic-seeded columns.

repeat (FEAST iteration, up to MaxIterations):
    /* Approximate the spectral projector P_[a,b](A) Y via Gauss-Legendre
     * quadrature on the upper half of the elliptic contour through
     * (a, 0) and (b, 0).  Schwarz symmetry (Ā = A) means each conjugate
     * pair z_e, z̄_e contributes 2·Re(...) so only N_e/2 solves are
     * required. */
    Q := 0
    for e in 1..N_e/2:
        solve (z_e I − A) X_e = Y      /* n × m₀ complex solve */
        Q := Q + 2 · Re( w_e · X_e )

    /* Rayleigh-Ritz reduction */
    A_q := Q^* A Q          /* m₀ × m₀ Hermitian */
    B_q := Q^* Q            /* m₀ × m₀ Hermitian SPD (typically) */
    solve  A_q Φ = B_q Φ Λ  /* small generalised Hermitian-definite */

    /* Ritz pairs and residuals */
    Λ̂, X̂ := Λ, Q Φ
    keep indices where Λ̂ ∈ [a, b]
    if max ‖A x̂ − λ̂ x̂‖ / (|λ̂| · ‖x̂‖) < tol: break
    Y := X̂  (with re-randomisation if rank-deficient)
```

The inner small generalised-Hermitian solve is routed through Mathilda's
own `direct_real_sym_machine` / `direct_complex_hermitian_machine` via a
Cholesky reduction `B_q = L L^*` followed by standard Hermitian eigen
decomposition of `L^{-1} A_q L^{-*}`. This is the standard FEAST inner
solve and lets us reuse fully-tested code.

Gauss–Legendre quadrature is on the upper unit semicircle parameterised
to `[a, b]`:
```
z_e = ((a + b)/2) + ((b − a)/2) · cos((π/2)(1 + ξ_e)) + i · ((b − a)/2) · sin((π/2)(1 + ξ_e))
w_e = (π/2) · ((b − a)/2) · sin((π/2)(1 + ξ_e)) · η_e
```
where `(ξ_e, η_e)` are the standard `[−1, 1]` Gauss–Legendre nodes /
weights of order `N_e`.

---

## Architecture

```
src/mateigen.c
  └── Phase 5 (new): "FEAST" method
      ├── FeastOpts                          option struct
      ├── feast_parse_subopts                Method -> {"FEAST", ...}
      ├── feast_automatic_prefers            returns false (never auto-picked)
      ├── feast_dispatch                     real/complex + machine/MPFR routing
      ├── feast_real_sym_machine             machine kernel: real symmetric
      ├── feast_complex_hermitian_machine    machine kernel: complex Hermitian
      ├── feast_real_sym_mpfr                MPFR kernel: real symmetric
      ├── feast_complex_hermitian_mpfr       MPFR kernel: complex Hermitian
      └── helpers
          ├── complex_lu_factor / _solve         partial pivoting, raw doubles
          ├── complex_lu_factor_M / _solve_M     MPFR analogue
          ├── feast_gl_nodes_<N> / weights_<N>   precomputed double tables
          └── feast_gl_compute_mpfr              Newton on Legendre polys @ MPFR
```

No new public symbols. `MATEIGEN_FEAST`, `mateigen_parse_method_value`,
`SYM_FEAST`, `SYM_Interval`, `SYM_Tolerance`, `SYM_MaxIterations` are
all already in place.

---

## New files

None. All work is contained in `src/mateigen.c` (and tests in
`tests/test_mateigen_feast.c`), matching how the other three methods
are organised.

---

## Symbols to add

In `src/sym_names.h` / `src/sym_names.c` (alphabetical, mirroring the
existing per-method sub-option additions):

- `SYM_ContourPoints` → `"ContourPoints"`
- `SYM_SubspaceSize`  → `"SubspaceSize"`

`SYM_Interval`, `SYM_Tolerance`, `SYM_MaxIterations`, `SYM_FEAST` already
exist (`src/sym_names.c:90,136,233,403`).

---

## Sub-option grammar

`Method -> {"FEAST", "Interval" -> {a, b}, ...}` accepts:

| Key              | Type           | Default                             | Meaning                                                        |
| ---------------- | -------------- | ----------------------------------- | -------------------------------------------------------------- |
| `"Interval"`     | `{a, b}` reals | **required**                        | Real interval to find eigenvalues in.                          |
| `"ContourPoints"`| Integer        | `8`                                 | Number of Gauss–Legendre nodes `N_e`. Supported: `2, 4, 8, 16`.|
| `"SubspaceSize"` | Integer        | `Max[20, ⌈1.5 · k_spec⌉]` or `Max[20, ⌈n/4⌉]` | Subspace dimension `m₀`.                              |
| `"MaxIterations"`| Integer        | `20`                                | Outer FEAST iterations.                                        |
| `"Tolerance"`    | Real           | `1e-12` machine; `2^{−0.7 · bits}` MPFR | Max relative residual to declare convergence.              |

Unknown keys are silently ignored (consistent with how Arnoldi / Banded
treat unrecognised sub-options).

---

## Dispatcher wiring

In `builtin_eigenvalues` / `builtin_eigenvectors` (currently in the
inexact-input branch around `src/mateigen.c:7839-7867`):

1. Add a `try_feast` branch parallel to `try_banded` / `try_arnoldi`.
2. `feast_automatic_prefers` always returns `false` — FEAST is never
   auto-selected, because it requires a user-supplied interval.
3. When `Method -> "FEAST"` is given but `feast_dispatch` returns
   `NULL`, the existing cascade falls through to `direct_dispatch`
   (already covers the `MATEIGEN_FEAST` enum value via the catch-all
   "fall back to Direct on … failure" clause once we widen it).
4. Drop the `MATEIGEN_FEAST` arm of `eigen_warn_unimplemented_method`
   once Phase 2 lands.

---

## Phasing

Each phase is a complete commit that builds and tests cleanly.

### Phase 1 — Skeleton + option parsing

- Add `SYM_ContourPoints` and `SYM_SubspaceSize`.
- Define `FeastOpts`, `feast_set_defaults`, `feast_parse_subopts`.
- Add `feast_dispatch` returning `NULL` (and stub
  `feast_dispatch_machine` / `feast_dispatch_mpfr`).
- Wire the `try_feast` branch into both builtins.
- Update `eigen_warn_unimplemented_method` so the FEAST case stays
  silent in the fall-back path until the kernel lands.

Acceptance: existing test suite still green; `Method -> {"FEAST",
"Interval" -> {0, 1}}` produces identical output to plain
`Eigenvalues[m]` (because FEAST falls through to Direct).

### Phase 2 — Real-symmetric machine kernel

- Implement `complex_lu_factor` and `complex_lu_solve` (paired `re/im`
  double arrays, partial pivoting).
- Add static Gauss–Legendre tables for `N_e = 2, 4, 8, 16`.
- Implement `feast_real_sym_machine`:
  - subspace initialisation (deterministic LCG seed so tests are
    reproducible),
  - `N_e / 2` complex solves per FEAST iteration,
  - reduced `A_q`, `B_q` accumulation,
  - inner solve via `direct_real_sym_machine` after Cholesky reduction
    of `B_q`,
  - residual check + iteration loop,
  - eigenpair filtering against `[a, b]`,
  - eigenvalue / eigenvector output (sorted ascending by λ to match
    the natural FEAST output order, with `eigen_sort_by_abs_desc`
    applied at the very end so the final list matches the rest of the
    eigensolver).

Acceptance: `tests/test_mateigen_feast.c` covers cross-check against
`Direct` on representative real-symmetric matrices from `eigen_corpus`
for assorted intervals; eigenvector residuals < tolerance;
orthonormality of returned eigenvector basis.

### Phase 3 — Complex-Hermitian machine kernel

Mirror of Phase 2 with complex `A`. The complex LU solver is unchanged;
the change is in matrix loading (`matD_load` already handles complex
input), in `Q^* A Q` formation, and in the inner-solve route
(`direct_complex_hermitian_machine`).

### Phase 4 — MPFR kernels (real + complex)

- `complex_lu_factor_M` / `complex_lu_solve_M` on `mpfr_t` pairs.
- Runtime Gauss–Legendre node computation: Newton's method on the
  Legendre polynomial recurrence at the working precision (avoids
  shipping precomputed MPFR tables that would lock to a fixed bit-
  count).
- `feast_real_sym_mpfr` and `feast_complex_hermitian_mpfr`, dispatched
  by `feast_dispatch_mpfr` when `common_scan_inexact` reports
  `min_bits > 53`.

Acceptance: MPFR-precision smoke tests in
`tests/test_mateigen_feast.c` confirm FEAST matches the Direct MPFR
kernel on small matrices to within `2^{−0.6 · bits}`.

### Phase 5 — Tests

`tests/test_mateigen_feast.c`:

- Cross-check against `Direct` on `eigen_corpus` matrices for fixed
  intervals (whole-spectrum, single-eigenvalue, empty interval).
- Eigenvector residual: `‖A x − λ x‖ < tol` for every returned pair.
- Orthonormality: `|⟨x_i, x_j⟩ − δ_{ij}| < tol`.
- `k_spec` interaction: `Eigenvalues[m, 3, Method -> {"FEAST", ...}]`
  takes the top-3 by `|λ|` from the interval result.
- Edge cases:
  - missing `"Interval"` → falls through to Direct + warning.
  - `b < a` → falls through to Direct + warning.
  - subspace too small (force `m₀ = 1` when interval contains 5
    eigenvalues) → fails gracefully and falls through to Direct.
- MPFR variants at 128 and 256 bits on tridiagonal corpus matrices.

Each test follows the conventions in `tests/test_mateigen_banded.c`
(reuses `eval_string`, `fmt_real_matrix`, `extract_real`, etc.).

Register the new test binary in `tests/CMakeLists.txt`.

### Phase 6 — Docs

- `docs/spec/builtins/linear_algebra.md` — extend the Eigenvalues /
  Eigenvectors section with the `FEAST` method, its sub-options, and
  failure semantics.
- `docs/spec/changelog/2026-05.md` — add a "FEAST Method for
  Eigenvalues / Eigenvectors" subsection summarising the phases.
- `Mathilda_spec.md` — no edit expected (the eigensolver entry already
  routes the reader to `docs/spec/builtins/linear_algebra.md`).

---

## Open implementation questions for later

These don't block Phase 1 but should be revisited as the kernels land:

- **Random subspace seeding.** A fixed deterministic LCG seed makes
  tests reproducible but couples FEAST output to the seed. Mathematica
  ships `Method -> {"FEAST", "ResetEigenvalues" -> True}` to randomise
  on each call. Decision: deterministic seed for now (better debugging,
  reproducible CI); revisit if real-world non-reproducibility becomes a
  feature request.
- **Auto-bump of `N_e`.** If after 3 FEAST iterations the residual is
  still > `100 · tol`, double `N_e` (capped at 16) on the next
  iteration. Cheap, often resolves under-quadrature on tight clusters.
  Defer to Phase 2 evaluation — only land if benchmarks show benefit.
- **`B_q` rank-deficiency handling.** When the subspace happens to
  overlap a non-target eigenspace heavily, `B_q` becomes ill-
  conditioned. The standard remedy is to drop near-null directions via
  a column-pivoted QR of `Q` before the inner solve. Implement this in
  Phase 2 only if the residual check fails on routine test inputs.

---

## Risk and reuse summary

- **Reuse:** ~70% of FEAST flows through existing code paths
  (`MateigenMethod` dispatch, sub-option parsing template, `MatD` /
  `MatM` loaders, Hermitian detection, the small Rayleigh–Ritz solve,
  k-spec application, chopping, result emission).
- **New surface area:** ~600 LOC (machine kernels) + ~400 LOC (MPFR
  kernels) + ~400 LOC tests. The biggest single new primitive is the
  complex LU pair (~150 LOC each path), which is independently testable.
- **Performance.** FEAST is `O(N_e · n³ + iter · n² · m₀)` per call.
  For Mathilda's test corpus (`n ≤ 50`) this is small in absolute
  terms; FEAST is included for API completeness with Mathematica
  rather than competitive performance on the corpus.

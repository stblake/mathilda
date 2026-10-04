---
source: src/numerical_roots/nsolve.c
references:
  - "H. M. Möller and H. J. Stetter, *Multivariate polynomial equations with multiple zeros solved by matrix eigenproblems*, Numer. Math. **70** (1995) 311–329."
  - "D. A. Cox, J. Little and D. O'Shea, *Using Algebraic Geometry*, 2nd ed. (Springer, 2005), ch. 2 — the eigenvalue method for zero-dimensional ideals."
---
**Algorithm.** `builtin_nsolve` reads `NSolve[expr [, vars [, dom [, prec]]]]`,
peels options from the tail, defaults `vars` to the collected non-constant
symbols and `dom` to Complexes, and dispatches:

1. **Univariate polynomial** → `NRoots` directly (NRoots never frees its
   argument), forwarding `PrecisionGoal`/`AccuracyGoal` so its polishing and
   accuracy contract govern the roots; the disjunction is repackaged into
   `{{x -> r1}, …}` with the `Reals` filter and `MaxRoots` cap applied. A huge
   literal exponent is guarded (`NSolve::deg`) before any machinery allocates.
2. **Square zero-dimensional polynomial system** (`nsolve_system.c`) → the
   **eigenvalue / multiplication-matrix (Möller–Stetter) method**: a greVlex
   Gröbner basis (`gb_buchberger`) gives the quotient ring `A = Q[x]/I`; its
   standard-monomial basis is enumerated; rational multiplication matrices
   `M_{x_i}` are built by normal-form reduction; a generic linear form `M_l = Σ
   c_i M_{x_i}` (deterministic seeded coefficients) is formed, and the
   eigenvalues/eigenvectors of `M_l` at MPFR precision
   (`eigen_all_eigenvectors_real_mpfr`) give each coordinate as `x_i(p) =
   (M_{x_i} v)[j]/v[j]`. Every candidate is verified against the original
   residuals. `Method -> "Symbolic"` instead does lexicographic **elimination**
   (solve the univariate generator with `NRoots`, back-substitute, recurse,
   verify); `"Homotopy"` currently routes to the same eigenvalue engine.
3. **Fallback** → symbolic `Solve` then numericalisation, dropping provably
   extraneous roots; a univariate non-polynomial last resort seeds `FindRoot`
   from a real grid (plus `±2i` unless `Reals`), verified and deduplicated.

Results are a list of rule-lists: `{}` no solutions, `{{}}` the universal
solution. Both the univariate path and both system solvers call `builtin_nroots`
directly, so the NRoots engines (Aberth / companion / Jenkins–Traub) are
NSolve's numeric backbone.

**Data structures.** The Gröbner engine works over `Q` (`GBPoly`); the
multiplication matrices are `mpq_t` rationals, the linear-form matrix and its
per-variable companions `mpfr_t`, and the eigen buffers and recovered
coordinates MPFR/`ncpx`. The standard-monomial basis is a flat `int[d·nvar]`.
`want_machine` holds when no precision digit count is given; otherwise the system
runs at `target_bits + max(32, target_bits/2)` bits and emits MPFR values.

**Complexity / limits.** The eigenproblem is `O(d³)` in the quotient-ring
dimension `d`; hard caps `NSYS_MAX_DIM = 256`, `NSYS_MAX_BOX = 200000`, and a
per-generator total-degree gate `NSYS_MAX_TDEG = 60` make a too-large or
positive-dimensional system fall back / stay unevaluated. The univariate degree
guard is `NSOLVE_MAX_POLY_DEGREE = 10000`. Options: `MaxRoots`, `Method`
(`Automatic` | `"EndomorphismMatrix"` | `"Homotopy"` | `"Symbolic"`),
`WorkingPrecision` (also a trailing positional digit count), `AccuracyGoal`
(default `MachinePrecision`, forwarded to NRoots), `PrecisionGoal`,
`VerifySolutions` (default on), `RandomSeeding` (seed for the generic linear
form, default 1234). `NSolve[expr, vars, Reals]` filters to real values; the
`Integers` domain is left to `Solve`. Diagnostics route through `mth_message`.

# Group E benchmarks — 10 advanced-numerical-analysis experiments (Mathilda vs Python 3.11 + numpy/scipy)

## Plan
- [x] Wire a new report group E (53–62) into `benchmarks/run_all.py` `group_of()`
- [x] Add `benchmarks/requirements.txt` (numpy/scipy/mpmath) + group-E table & Python-3.11 note in README
- [x] 53 matrix-decompositions (LU/QR/SVD/PseudoInverse/rank/NullSpace; Cholesky ABSENT)
- [x] 54 eigenproblems (sym/general/generalized/Arnoldi; Eigensystem ABSENT)
- [x] 55 vectorized special functions (Fresnel/Erf/integrals/Beta vectorized; BesselI/K no kernel)
- [x] 56 multidim quadrature (2D/3D/oscillatory/singular/semi-infinite/dependent-bounds)
- [x] 57 stiff ODE + PDE (stiff scalar/Robertson/harmonic/VdP + heat/wave MethodOfLines)
- [x] 58 nonlinear systems (2×2, Burden–Faires 3×3, Broyden tridiagonal N=10/40)
- [x] 59 polynomial roots (NSolve/NRoots dense 20/50/100, roots-of-unity, Wilkinson)
- [x] 60 DCT/DST (types 1/2/4, DST 1/2) + 2D Fourier, normalization reconciled in checks
- [x] 61 regularized least squares (LeastSquares + Fit Tikhonov/ridge; FindFit ABSENT)
- [x] 62 arbitrary precision (N[…,p]/NSum/NIntegrate/FindRoot vs mpmath), deep-digit checks
- [x] Run subset under Python 3.11 via `HPC_PYTHON`, `--check-labels`; produce gap report

## Review
- **Result:** 65 cases, **0 INCOMPLETE, 0 CHECK-FAIL** (every case runs, every check agrees within 1e-6),
  32 SLOWER, 30 AHEAD, 3 ABSENT. Coverage 85.7% (30/35 declared heads). Wall clock 1.8 min.
- Only shared-code edit: one clause in `run_all.py` `group_of()`. No Mathilda source touched.
- Run it: `HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 53,54,55,56,57,58,59,60,61,62 --system mathilda,python`
- Outputs: `benchmarks/REPORT.partial.md`, `ABSENT.partial.md`, `results/2026-08-12-partial.json`.

### Top gaps surfaced (the "drive improvements" queue)
1. `NullSpace` on a float matrix takes a non-machine path — **2428×** (9.6 s vs 4 ms).
2. `NSolve`/`NRoots` high-degree — **60–110×**; Wilkinson-15 **1616×** (symbolic preprocessing of the product form).
3. Generalized `Eigenvalues[{A,B}]` has **no LAPACK path** — symbolic char-poly root-finding; **1136×**, returns Root[] at n≥6, hangs at n≥8.
4. `BesselI`/`BesselK` over arrays have **no SIMD kernel** (scalar threading) — **105×/209×**; the other special functions are vectorized and mostly AHEAD.
5. `LUDecomposition` **31.5×**, `FindRoot` systems scale poorly (Broyden N=40 **19×**), symmetric eig **6.7×**, Fit ridge **8.1×**, DCT-2/4 **~3.6×**.
- Also `FindRoot`/`NSum` under-deliver requested WorkingPrecision (FindRoot ~19 correct digits at WP→100).

### Where Mathilda already wins (regression guards)
- NDSolve (compiled RHS) beats scipy solve_ivp 5–40×; arbitrary precision beats (pure-Python) mpmath up to ~200×
  (N[Gamma[1/3],1000]: 19 ms vs 3.9 s); most vectorized special functions and several NIntegrate cases AHEAD.

### Absences (feature work)
- `CholeskyDecomposition`, `Eigensystem`, `FindFit` (declared, benched on the Python side only).

## Not done (deferred per scope)
- Implementing the kernel/feature fixes above — the user chose "author + run + gap report, then review".

---

# NMinimize / NMaximize — numerical global optimization (2026-08-14)

Plan: `~/.claude/plans/let-s-plan-the-implementation-sequential-cosmos.md`.
Global-optimization driver layered on the existing `FindMinimum` machinery
(`src/findmin.c`), reusing its variable binding, `{f,cons}` constraint parsing,
penalty/BFGS local solvers, MPFR path, and result builders.

## Done
- `builtin_nminimize` + `nm_minimize_driver` + 4 engines (DifferentialEvolution
  default, NelderMead, RandomSearch, SimulatedAnnealing) with Deb feasibility
  rules; `builtin_nmaximize` negate-and-wrap. All in `src/findmin.c`.
- Mixed-integer via `Element[x, Integers]` (lattice search + integer coordinate
  descent; integer results). Empty feasible set → `{Infinity, {x->Indeterminate}}`.
- `SYM_NMinimize`/`SYM_NMaximize` (sym_names), `HoldAll|Protected` registration,
  docstrings (info.c), `Options[]` defaults (options_builtin.c), docs
  (calculus.md) + changelog (2026-08-10.md).
- Deterministic SplitMix64 PRNG (fixed seed; `"RandomSeed"` override) → reproducible.
- Internal solver diagnostics muted during search (`g_fm_quiet`) so successful
  solves are silent.

## Verification
- `tests/test_nminimize.c` (29 cases) — all pass. `findmin_tests`,
  `options_tests`, `core_tests` — no regression. `make check-c99` clean.
  Valgrind identical to `Sin[1.0]` baseline (no NMinimize frames). Spec examples
  reproduced (quartic, disk, linear-program, equality, integer, infeasible, dual).

## Two bugs found + fixed during bring-up
1. `{f, c1, c2, ...}` (>2-element list) silently dropped all constraints — only a
   2-element `{f, cons}` was recognized. Now collects every trailing element into
   an implicit `And`.
2. Contradictory strict bounds (`x>2 && x<1`) were "repaired" by averaging,
   hiding infeasibility. Now detects an empty box → `{Infinity, ...}`.

## Deferred (per agreed scope)
- Vector/matrix/region variables (`Vectors`/`Matrices`), `VectorGreaterEqual`,
  `Or` constraints, non-Integers/Reals domains, general constraints at high
  `WorkingPrecision`. Each abstains with an `NMinimize::nimpl` message.

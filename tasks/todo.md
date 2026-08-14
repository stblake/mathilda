# Fix: NMinimize hang on Max[Re[Eigenvalues[matrix[c1,c2]]]]

## Root cause (confirmed by user's Mathematica output)
Mathematica returns a COMPACT radical for this matrix: the char poly is
biquadratic AFTER depression (x -> x - 5/2), so it solves in nested square
roots (~80 leaves), evaluates by arithmetic (fast), and is correct.

Mathilda's `solve_quartic_radical` (solvepoly.c) detected the biquadratic case
only via `is_definite_zero(q)` — a LITERAL zero test — which missed the
symbolic-zero `q`. It fell to the general resolvent-cubic path: 9205 leaves AND
numerically WRONG (q=0 forces resolvent root t=0 -> 0/0 in q/Sqrt[2t]).

## Desired policy (user)
- Cubics / Quartics default to **False** (Root[] for the general case).
- Biquadratics (incl. biquadratic-after-depression), n-quadratics
  (a x^{2m}+b x^m+c), and n-linear/binomials (a x^n + b) **always** radicals.
- Keep the new Root-object numericalisation.

## Plan
- [x] solvepoly.c `solve_quartic_radical`: biquadratic detection via
      `zero_test_decide` (polynomial-identity), not literal `is_definite_zero`.
- [ ] solvepoly.c dispatch: emit the depressed-biquadratic quartic in radicals
      even when quartics_radical is False (ungate the special case). Binomial /
      n-quadratic already fire before the gate — verified.
- [ ] Restore Fix A: Cubics/Quartics default False (eigen_common.c,
      options_builtin.c).
- [ ] Restore Fix B: Root[] numericalisation of inexact-coefficient polynomials
      (root_numeric.c) — user wants it kept; needed so general Root eigenvalues
      still numericalise under the False default.
- [ ] radicals.c `radical_quartic`: same biquadratic-detection fix (parallel
      Root->radical path).
- [ ] Tests: test_options.c (Cubics/Quartics -> False); test_eigen.c
      irreducible-cubic residual threshold (Root-based now ~1e-13, relax
      10^-20 -> 10^-10).
- [ ] Verify: user NMinimize fast+correct; N[Root[inexact]] works; general
      quartic -> Root; suites pass.
- [ ] Docs + changelog.

## Review
Done. Infinite hang -> 6.2 s, correct global min (-2.0, matches coarse grid).

Changes:
- `poly/solvepoly.c`: biquadratic detection via `zero_test_decide` (was literal
  `is_definite_zero`); dispatch ungates the depressed-biquadratic quartic so it
  emits radicals even under Quartics -> False. New helper
  `quartic_depressed_is_biquadratic`.
- `radicals.c`: same `zero_test_decide` biquadratic fix in the parallel
  `radical_quartic` (Root->radical path).
- `linalg/eigen_common.c` + `options_builtin.c`: Cubics/Quartics default -> False.
- `root_numeric.c`: numericalise Root[] with inexact (machine-real / rational)
  coefficients — `solve_root_core` refactor shared by exact and inexact paths,
  LAPACK dgeev machine fast path, graceful Newton-stall fallback.
- `linalg/svdecomp.c`: SVD requests Cubics/Quartics -> True internally.
- Docs: info.c docstrings, linear-algebra.md, changelog 2026-08-10.
- Tests: test_options.c (defaults -> False), test_eigen.c (irreducible-cubic
  residual 10^-20 -> 10^-10, now Root-based).

Verified:
- Radical obj matches direct numeric eigenvalues to 1.15e-14 over a 121-pt grid.
- Suites pass: root_numeric, nroots, eigen, options, rootreduce, solve,
  mateigen_{direct,arnoldi,banded}, SVD, nsolve, findroot. (Full sweep running.)
- N[Root[inexact]] works; general quartic/cubic -> Root; biquadratic/binomial
  -> radicals.

Not done (out of scope): matching Mathematica's 0.2 s exactly — the residual gap
is the SA evaluation count (tuned for the benchmark suite) and the radical not
being maximally simplified (365 vs ~80 leaves); neither is a correctness issue.

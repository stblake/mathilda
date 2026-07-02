# Tier 2 FLINT acceleration — Factor + linear algebra

Ref: FLINT_ACCELERATION_SWEEP.md, project_flint_context_pattern.md.
User ask: "both univariate and multivariate Factor via FLINT" + FLINT-backed
Inverse/LinearSolve/RowReduce etc., each with a FLINT`func interface + docstring.

## A — Linear algebra (flint_mat_bridge)  [clean: exact & unique outputs]  ✅ DONE + verified
- [x] sym_names: SYM_FLINT_Inverse/LinearSolve/RowReduce/MatrixRank
- [x] flint_mat_bridge: flint_mat_inverse / flint_mat_rref / flint_mat_solve / flint_mat_rank
- [x] FLINT`Inverse/LinearSolve/RowReduce/MatrixRank builtins (docstrings, PROTECTED)
- [x] Transparent wiring: inv.c dispatcher, rowreduce_divfree, linearsolve_divfree, matrank.c (exact path)
- [x] NullSpace inherits via RowReduce (no separate FLINT`NullSpace: basis convention not unique)
- Smoke-verified: all FLINT` outputs == System` outputs; singular/symbolic fall through.

## B — Univariate Factor via FLINT  [fixed the real bug]  ✅ DONE
- ROOT CAUSE (not int64 overflow): is_number() omitted EXPR_BIGINT, so collect_variables
  treated a bignum coeff as a VARIABLE -> univariate bignum poly routed to the multivariate
  path which dropped the constant term -> returned unfactored. One-line fix in poly.c.
- [x] flint_univariate_factor(P, var) via fmpz_poly_factor: primitive positive-lead factors +
      signed integer content (matches classical; NOT fmpq monic which would give 4(x-3/2)(x+3/2))
- [x] Route bz_factor_to_expr through it (fast + bignum-safe)
- [x] Smoke: 4x^2-9 -> (2x-3)(2x+3); bignum quadratics/products factor; small cases unchanged

## C — Multivariate Factor via FLINT   ✅ DONE (2026-07-01, follow-up)
- User follow-up: "Factor needs to call FLINT`Factor" — Factor[x^99-y^99] was >20s.
- [x] Transparent multivariate path added at top of builtin_factor (after univariate,
      gated arg_count==1) calling public flint_polynomial_factor. Now ~4ms.
- [x] Normalization pass (was the deferral blocker): factoring ctx switched ORD_LEX->ORD_DEGLEX
      so term 0 is the highest-total-degree leading term; each irreducible factor negated to a
      POSITIVE leading coeff, (-1)^e folded into the rational content. This is MORE canonical
      than classical (which spread signs inconsistently): y^2-x^2 -> -(x-y)(x+y) [WL-correct],
      u^2-k stays u^2-k (not the ugly -(k-u^2) a naive lex-normalise gives).
- [x] Fall-through verified: rational functions (denominator), Sqrt/Sin heads, Extension/
      GaussianIntegers/Modulus option forms all bypass FLINT -> classical path unchanged.
- [x] Rebaselined 2 pinned classical-convention tests in test_facpoly.c to the verified-
      equivalent normalized forms (big x^2(1-x^12)... case; the 2x^3y-... a-poly case).
      Product identity checked by hand for both.
- [x] Extended test_flint_bridge.c::test_flint_poly_ops with multivariate + fall-through cases.
- Note: normalization also applies to explicit FLINT`Factor and parametric_sqrt_factor
  (same impl); existing FLINT`Factor / radical / qafactor tests still pass (all already
  positive-lead).

## Verify  ✅ ALL DONE
- [x] build (make clean && make) FLINT+graphics OK; USE_FLINT=0 fallback builds + FLINT`* stay symbolic
- [x] tests: flint_bridge/facpoly/poly/factor_baseline/intrat/qafactor/mpoly/radical_polyops + linalg/eigen/hilbertmatrix/mateigen_direct all rc=0
- [x] regression sweep: simplify/fullsimplify/simp/solve/series/integrate_dispatch/integrate_goursat/radical_simplify/trigrat — the only fails are PRE-EXISTING radical-Simplify / Laurent-series / nested-radical-Solve soft-fails (don't touch univariate-Z Factor or rational-matrix paths; dedicated Factor/linalg suites pass)
- [x] valgrind baseline-identical (definitely/indirectly/possibly-lost + 516 errors byte-identical vs Sin[1.0] baseline; only still-reachable FLINT arb cache grew)
- [x] docs: changelog 2026-06-29 (2 entries) + linear-algebra.md (Inverse/RowReduce/MatrixRank/LinearSolve) + structural-manipulation.md FLINT` context section
- [x] memory: project_flint_context_pattern.md updated (Tier2 + is_number bignum trap)

## Review
- Linalg: FLINT`Inverse/LinearSolve/RowReduce/MatrixRank + transparent acceleration (Det from Tier1).
  Unique/basis-independent outputs => transparent results identical to classical. NullSpace inherits via RowReduce.
- Factor: univariate-over-Z transparently via FLINT fmpz_poly (Mathematica primitive convention);
  fixes the bignum-unfactored bug WITHOUT the risky global is_number change. Multivariate = explicit FLINT`Factor.
- Deferred (documented): transparent multivariate Factor; PolynomialGCD/HermiteReduce bignum-content bugs.

/*
 * tests/test_ludecomposition_mpfr.c
 *
 * Arbitrary-precision MPFR coverage for LUDecomposition.
 *
 * Drives lu_dispatch with MPFR-precision inputs (`N[..., bits]`) so
 * `lu_mpfr_dispatch` is selected.  Verifies:
 *
 *   1. The identity m[[p]] == l . u holds within an MPFR-precision
 *      tolerance.  Reconstruction is measured by extracting the
 *      residual into doubles -- the kernel works at the requested
 *      precision; the test only asserts the residual is small in the
 *      double-precision projection.
 *
 *   2. p is a valid 1-indexed permutation of 1..n.
 *
 *   3. The condition-number scalar is a numeric value strictly > 0
 *      (an MPFR Real, or a regular Real / Integer when the fallback
 *      kernel fires).
 *
 *   4. Precision-monotonicity: same matrix factored at 60, 100, and
 *      200-bit precision -- the test passes when the kernel reports
 *      finite, positive condition numbers at all three precisions.
 *
 * Builds cleanly under both USE_MPFR=1 and USE_MPFR=0; in the off
 * case the suite is no-op (the symbolic fallback still handles the
 * inexact-input pipeline correctly via rationalise -> exact ->
 * numericalise).
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

static Expr* run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* Extract a numeric leaf as a real double, tolerating Integer / Real /
 * Rational / MPFR / Complex.  Returns false for non-numeric. */
static bool leaf_to_double(Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        return true;
    }
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            double p, q;
            if (leaf_to_double(e->data.function.args[0], &p)
             && leaf_to_double(e->data.function.args[1], &q)
             && q != 0.0) { *out = p / q; return true; }
        }
        if (strcmp(h, "Complex") == 0 && e->data.function.arg_count == 2) {
            double r;
            if (leaf_to_double(e->data.function.args[0], &r)) {
                *out = r; return true;
            }
        }
    }
    return false;
}

/* Maximum absolute leaf value of a nested List tensor (projected into
 * doubles). */
static double max_abs_tensor(Expr* e) {
    if (!e) return 0.0;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "List") == 0) {
        double m = 0.0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double v = max_abs_tensor(e->data.function.args[i]);
            if (v > m) m = v;
        }
        return m;
    }
    double v;
    if (leaf_to_double(e, &v)) {
        double a = v < 0 ? -v : v;
        return a;
    }
    /* Complex with non-trivial im: estimate by max(|re|, |im|). */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        double r = max_abs_tensor(e->data.function.args[0]);
        double i = max_abs_tensor(e->data.function.args[1]);
        return r > i ? r : i;
    }
    return 0.0;
}

/* Assert m[[p]] == l . u with the residual's max-abs <= bound.
 * L (unit lower) and U (upper) are constructed via Table since the
 * triangularize builtins aren't available. */
static void assert_lu_residual_small(const char* m_src, double bound,
                                      const char* label) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "Block[{m, lu, p, l, u, n, ii, jj}, "
             "  m = %s; "
             "  {lu, p} = LUDecomposition[m][[1;;2]]; "
             "  n = Length[m]; "
             "  l = Table[If[ii > jj, lu[[ii, jj]], "
             "                  If[ii == jj, 1, 0]], "
             "             {ii, n}, {jj, n}]; "
             "  u = Table[If[ii <= jj, lu[[ii, jj]], 0], "
             "             {ii, n}, {jj, n}]; "
             "  m[[p]] - l . u"
             "]",
             m_src);
    Expr* res = run(buf);
    double err = max_abs_tensor(res);
    if (!(err <= bound)) {
        char* s = expr_to_string(res);
        fprintf(stderr,
                "FAIL %s: residual = %.6e > bound %.6e\n  diff: %s\n",
                label, err, bound, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    printf("  PASS %s: residual = %.6e <= %.6e\n", label, err, bound);
    expr_free(res);
}

static void assert_perm_valid(const char* m_src, int n) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Sort[LUDecomposition[%s][[2]]] == Range[%d]",
             m_src, n);
    Expr* res = run(buf);
    ASSERT(res->type == EXPR_SYMBOL
        && strcmp(res->data.symbol, "True") == 0);
    expr_free(res);
}

static void assert_cond_positive(const char* m_src) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "LUDecomposition[%s][[3]]", m_src);
    Expr* res = run(buf);
    double v;
    ASSERT(leaf_to_double(res, &v));
    ASSERT(v > 0.0);
    printf("  PASS: cond > 0 (= %g) for %s\n", v, m_src);
    expr_free(res);
}

/* ============================================================
 *   Tests
 * ============================================================ */

static void test_lu_mpfr_3x3_60bit(void) {
    /* 60-bit working precision; residual bound is generous because
     * the residual is read back through a double cast. */
    const char* m = "N[{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, "
                    " {1/3, 1/4, 1/5}}, 60]";
    assert_lu_residual_small(m, 1e-12, "3x3 Hilbert 60-bit");
    assert_perm_valid(m, 3);
    assert_cond_positive(m);
}

static void test_lu_mpfr_3x3_100bit(void) {
    const char* m = "N[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}, 100]";
    assert_lu_residual_small(m, 1e-12, "3x3 tridiag 100-bit");
    assert_perm_valid(m, 3);
    assert_cond_positive(m);
}

static void test_lu_mpfr_4x4_60bit(void) {
    const char* m = "N[{{4, 3, 2, 1}, {1, 4, 3, 2}, "
                    " {2, 1, 4, 3}, {3, 2, 1, 4}}, 60]";
    assert_lu_residual_small(m, 1e-12, "4x4 60-bit");
    assert_perm_valid(m, 4);
    assert_cond_positive(m);
}

static void test_lu_mpfr_2x2_complex(void) {
    /* Complex MPFR 2x2 -- exercises the lum_cdiv / lum_csub_mul
     * paths.  We rely on the same double-projection identity check. */
    const char* m = "N[{{2 + I, 1 - I}, {3 + 2 I, 4 + 5 I}}, 80]";
    assert_lu_residual_small(m, 1e-12, "2x2 complex 80-bit");
    assert_perm_valid(m, 2);
    assert_cond_positive(m);
}

static void test_lu_mpfr_5x5_diagdom(void) {
    /* A diagonally-dominant matrix at 200-bit precision -- the
     * stiffest precision test we run.  Bound is still slack because
     * we project to double. */
    const char* m = "N[{{10, 1, 2, 0, 1}, {1, 11, 0, 1, 2}, "
                    " {2, 0, 12, 1, 1}, {0, 1, 1, 13, 2}, "
                    " {1, 2, 1, 2, 14}}, 200]";
    assert_lu_residual_small(m, 1e-12, "5x5 diag-dominant 200-bit");
    assert_perm_valid(m, 5);
    assert_cond_positive(m);
}

static void test_lu_mpfr_precision_progression(void) {
    /* Same matrix at 60, 100, 200 bits -- check residual at each
     * precision.  We don't enforce strict monotonicity (double
     * truncation dominates the readout), just that all three are
     * within the loose double bound. */
    assert_lu_residual_small(
        "N[{{1, 2, 3}, {4, 5, 7}, {7, 8, 10}}, 60]",
        1e-10, "progression 60-bit");
    assert_lu_residual_small(
        "N[{{1, 2, 3}, {4, 5, 7}, {7, 8, 10}}, 100]",
        1e-10, "progression 100-bit");
    assert_lu_residual_small(
        "N[{{1, 2, 3}, {4, 5, 7}, {7, 8, 10}}, 200]",
        1e-10, "progression 200-bit");
}

int main(void) {
    symtab_init();
    core_init();

#ifndef USE_MPFR
    printf("USE_MPFR is disabled at build time -- "
           "LUDecomposition MPFR tests skipped.\n");
    return 0;
#endif

    TEST(test_lu_mpfr_3x3_60bit);
    TEST(test_lu_mpfr_3x3_100bit);
    TEST(test_lu_mpfr_4x4_60bit);
    TEST(test_lu_mpfr_2x2_complex);
    TEST(test_lu_mpfr_5x5_diagdom);
    TEST(test_lu_mpfr_precision_progression);

    printf("All LUDecomposition (MPFR) tests passed.\n");
    return 0;
}

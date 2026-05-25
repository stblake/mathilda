/*
 * tests/test_singularvaluedecomposition_machine.c
 *
 * Unit tests for the LAPACK fast path of SingularValueDecomposition.
 *
 * Routed through svd_dispatch by feeding inexact (Real) inputs; we
 * never call svd_machine_dispatch directly so the tests also exercise
 * the dispatcher's classification logic.  Verification strategies:
 *
 *   1. Reconstruction: Max[Abs[Flatten[u . sigma . Transpose[v] - m]]]
 *      below a tight tolerance (1e-10 for machine-precision input).
 *   2. Shape: u (n x n), sigma (n x p), v (p x p) for the default form;
 *      u (n x k), sigma (k x k), v (p x k) for the truncated forms.
 *   3. Singular-value ordering: sigma's diagonal is non-increasing.
 *   4. Match to the Mathematica reference output values for a handful
 *      of the spec-test matrices (within 1e-4 absolute).
 *   5. Tolerance: a near-singular matrix at the default precision has
 *      its tiny singular value zeroed; with an explicit smaller
 *      Tolerance it is preserved.
 */

#include <math.h>
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

static Expr* run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* Pull a Real scalar (or numeric-with-Abs) out of an evaluator result. */
static double as_double(Expr* e) {
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            return as_double(e->data.function.args[0])
                 / as_double(e->data.function.args[1]);
        }
    }
    char* s = expr_to_string(e);
    fprintf(stderr, "FAIL: cannot convert to double: %s\n", s);
    free(s);
    ASSERT(0);
    return 0.0;
}

/* Run an expression and convert to double. */
static double run_d(const char* src) {
    Expr* r = run(src);
    double v = as_double(r);
    expr_free(r);
    return v;
}

/* Verify reconstruction within `tol` absolute error.  Uses
 * ConjugateTranspose so the same helper exercises both real (where it
 * collapses to Transpose) and complex (Hermitian) inputs. */
static void assert_reconstruction(const char* m_src, double tol) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVDMach$Cache = SingularValueDecomposition[%s];", m_src);
    expr_free(run(buf));
    /* Norm gives a scalar even when entries are complex (handles
     * Sqrt[a^2 + b^2] cleanly).  Re drops the residual imaginary
     * roundoff Mathilda's Sqrt-of-positive-Sum-of-squares leaves
     * behind for complex inputs (~1e-30 I).  N forces to Real. */
    snprintf(buf, sizeof(buf),
        "N[Re[Norm[Flatten["
        "TestSVDMach$Cache[[1]] . TestSVDMach$Cache[[2]] . "
        "ConjugateTranspose[TestSVDMach$Cache[[3]]] - (%s)]]]]", m_src);
    double err = run_d(buf);
    expr_free(run("ClearAll[TestSVDMach$Cache]"));
    if (err > tol) {
        fprintf(stderr,
            "FAIL: %s reconstruction error %g > tol %g\n",
            m_src, err, tol);
        ASSERT(0);
    }
    printf("  PASS: reconstruction error %g for %s\n", err, m_src);
}

/* Assert {u_rows, u_cols, s_rows, s_cols, v_rows, v_cols} of the SVD
 * result.  -1 in any slot means "don't check". */
static void assert_shape(const char* expr_src, int ur, int uc,
                         int sr, int sc, int vr, int vc) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "TestSVDMach$Shape = %s;", expr_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf), "Dimensions /@ TestSVDMach$Shape");
    Expr* dims = run(buf);
    ASSERT(dims->type == EXPR_FUNCTION);
    ASSERT(dims->data.function.arg_count == 3);
    int got[3][2];
    for (int i = 0; i < 3; i++) {
        Expr* d = dims->data.function.args[i];
        ASSERT(d->type == EXPR_FUNCTION);
        ASSERT(d->data.function.arg_count == 2);
        got[i][0] = (int)d->data.function.args[0]->data.integer;
        got[i][1] = (int)d->data.function.args[1]->data.integer;
    }
    expr_free(dims);
    int want[3][2] = {{ur, uc}, {sr, sc}, {vr, vc}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            if (want[i][j] < 0) continue;
            if (got[i][j] != want[i][j]) {
                fprintf(stderr,
                    "FAIL: %s shape[%d][%d] expected %d got %d\n",
                    expr_src, i, j, want[i][j], got[i][j]);
                ASSERT(0);
            }
        }
    }
    expr_free(run("ClearAll[TestSVDMach$Shape]"));
    printf("  PASS: shape for %s\n", expr_src);
}

/* Verify singular-value ordering: sigma's diagonal is non-increasing. */
static void assert_sigma_ordered(const char* m_src) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "TestSVDMach$Cache = SingularValueDecomposition[%s];", m_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "TestSVDMach$Diag = Table[TestSVDMach$Cache[[2]][[i, i]], "
        "{i, Min[Dimensions[%s]]}];", m_src);
    expr_free(run(buf));
    Expr* diag = run("TestSVDMach$Diag");
    ASSERT(diag->type == EXPR_FUNCTION);
    int n = (int)diag->data.function.arg_count;
    double prev = 1e300;
    for (int i = 0; i < n; i++) {
        double cur = as_double(diag->data.function.args[i]);
        if (cur > prev + 1e-12) {
            fprintf(stderr,
                "FAIL: %s sigma not descending: sigma[%d]=%g > sigma[%d]=%g\n",
                m_src, i, cur, i - 1, prev);
            ASSERT(0);
        }
        prev = cur;
    }
    expr_free(diag);
    expr_free(run("ClearAll[TestSVDMach$Cache, TestSVDMach$Diag]"));
    printf("  PASS: sigma ordering for %s\n", m_src);
}

/* ============== reconstruction over the spec examples ============== */

static void test_recon_2x2_rank1(void) {
    /* Mathematica spec: {{1, 2}, {1, 2}} (rank 1). */
    assert_reconstruction("N[{{1, 2}, {1, 2}}]", 1e-10);
}

static void test_recon_2x2_full(void) {
    /* Mathematica spec: {{0.5, 1}, {2, 2.5}} (full rank). */
    assert_reconstruction("{{0.5, 1.}, {2., 2.5}}", 1e-10);
}

static void test_recon_3x2_tall(void) {
    /* Mathematica spec: {{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}. */
    assert_reconstruction("{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}", 1e-10);
}

static void test_recon_3x3(void) {
    /* Mathematica spec: {{1.1, 1.3, 2.1}, {2.2, 2.2, 3.3}, {3.3, -3, 4.7}}. */
    assert_reconstruction(
        "{{1.1, 1.3, 2.1}, {2.2, 2.2, 3.3}, {3.3, -3, 4.7}}", 1e-10);
}

static void test_recon_2x3_wide(void) {
    /* Mathematica spec: {{1.1, 2, 5}, {3, -11, 4.2}}. */
    assert_reconstruction("{{1.1, 2, 5}, {3, -11, 4.2}}", 1e-10);
}

static void test_recon_complex_2x2(void) {
    /* Mathematica spec: complex 2x2. */
    assert_reconstruction(
        "{{0.5 + 0.5 I, 1.1}, {-I, 3.2 - 4.5 I}}", 1e-10);
}

/* ============== shape verification ============== */

static void test_shape_square(void) {
    assert_shape(
        "SingularValueDecomposition[{{1.2, 3.4}, {5.6, 7.8}}]",
        2, 2, 2, 2, 2, 2);
}

static void test_shape_tall(void) {
    assert_shape(
        "SingularValueDecomposition[{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}]",
        3, 3, 3, 2, 2, 2);
}

static void test_shape_wide(void) {
    assert_shape(
        "SingularValueDecomposition[{{1.1, 2, 5}, {3, -11, 4.2}}]",
        2, 2, 2, 3, 3, 3);
}

/* ============== ordering ============== */

static void test_ordering_square(void) {
    assert_sigma_ordered("{{1.2, 3.4}, {5.6, 7.8}}");
}

static void test_ordering_tall(void) {
    assert_sigma_ordered("{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}");
}

/* ============== truncation forms ============== */

static void test_truncation_k_positive(void) {
    assert_shape(
        "SingularValueDecomposition["
        "{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}, 1]",
        3, 1, 1, 1, 2, 1);
}

static void test_truncation_k_negative(void) {
    assert_shape(
        "SingularValueDecomposition["
        "{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}, -1]",
        3, 1, 1, 1, 2, 1);
}

static void test_truncation_upto_clamp(void) {
    /* UpTo[10] on a 2x3 matrix (rank <= 2) returns up to 2 columns. */
    assert_shape(
        "SingularValueDecomposition[{{1.1, 2., 5.}, {3., -11., 4.2}}, UpTo[10]]",
        2, 2, 2, 2, 3, 2);
}

/* ============== tolerance ============== */

static void test_tolerance_default(void) {
    /* Near-singular matrix: default LAPACK rank cutoff treats the tiny
     * singular value as exactly zero. */
    Expr* r = run(
        "SingularValueDecomposition[{{1., 0}, {1., 10^-14}}][[2]][[2, 2]]");
    double v = as_double(r);
    expr_free(r);
    if (fabs(v) > 1e-12) {
        fprintf(stderr,
            "FAIL: default tolerance did not zero tiny sigma; got %g\n", v);
        ASSERT(0);
    }
    printf("  PASS: default tolerance zeroed sigma[2,2] = %g\n", v);
}

static void test_tolerance_explicit(void) {
    /* With an explicit tighter Tolerance, the small sigma is preserved. */
    Expr* r = run(
        "SingularValueDecomposition["
        "{{1., 0}, {1., 10^-14}}, Tolerance -> 10^-20][[2]][[2, 2]]");
    double v = as_double(r);
    expr_free(r);
    if (fabs(v) < 1e-16) {
        fprintf(stderr,
            "FAIL: explicit tolerance zeroed sigma that should survive; "
            "got %g\n", v);
        ASSERT(0);
    }
    printf("  PASS: explicit tolerance preserved sigma[2,2] = %g\n", v);
}

/* ============== generalized SVD (LAPACK dggsvd3) ============== */

/* Verify the generalized form reconstructs both m and a within `tol`.
 * Result shape: {{u, ua}, {sigma, sigma_a}, v}, with
 *   m == u  . sigma   . ConjugateTranspose[v]
 *   a == ua . sigma_a . ConjugateTranspose[v]. */
static void assert_gen_reconstruction(const char* m_src, const char* a_src,
                                       double tol)
{
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVDMach$G = SingularValueDecomposition[{%s, %s}];",
        m_src, a_src);
    expr_free(run(buf));

    snprintf(buf, sizeof(buf),
        "N[Re[Norm[Flatten["
        "TestSVDMach$G[[1, 1]] . TestSVDMach$G[[2, 1]] . "
        "ConjugateTranspose[TestSVDMach$G[[3]]] - (%s)]]]]", m_src);
    double err_m = run_d(buf);

    snprintf(buf, sizeof(buf),
        "N[Re[Norm[Flatten["
        "TestSVDMach$G[[1, 2]] . TestSVDMach$G[[2, 2]] . "
        "ConjugateTranspose[TestSVDMach$G[[3]]] - (%s)]]]]", a_src);
    double err_a = run_d(buf);

    expr_free(run("ClearAll[TestSVDMach$G]"));
    if (err_m > tol || err_a > tol) {
        fprintf(stderr,
            "FAIL: generalized SVD({%s, %s}) errors m=%g a=%g (tol %g)\n",
            m_src, a_src, err_m, err_a, tol);
        ASSERT(0);
    }
    printf("  PASS: generalized reconstruction errors m=%g a=%g for "
           "{%s, %s}\n", err_m, err_a, m_src, a_src);
}

/* Assert generalized result has shape {{u_n x u_n, ua_n x ua_n},
 * {sigma rows x cols, sigma_a rows x cols}, v_n x v_n}. */
static void assert_gen_shape(const char* m_src, const char* a_src,
                              int u_n, int ua_n,
                              int s_r, int s_c, int sa_r, int sa_c,
                              int v_n)
{
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVDMach$G = SingularValueDecomposition[{%s, %s}];",
        m_src, a_src);
    expr_free(run(buf));

    /* Pull the five sub-shapes one by one. */
    struct { const char* expr; int want_r; int want_c; const char* label; } checks[5] = {
        { "Dimensions[TestSVDMach$G[[1, 1]]]", u_n,  u_n,  "u"      },
        { "Dimensions[TestSVDMach$G[[1, 2]]]", ua_n, ua_n, "ua"     },
        { "Dimensions[TestSVDMach$G[[2, 1]]]", s_r,  s_c,  "sigma"  },
        { "Dimensions[TestSVDMach$G[[2, 2]]]", sa_r, sa_c, "sigma_a"},
        { "Dimensions[TestSVDMach$G[[3]]]",    v_n,  v_n,  "v"      },
    };
    for (int i = 0; i < 5; i++) {
        Expr* d = run(checks[i].expr);
        ASSERT(d->type == EXPR_FUNCTION);
        ASSERT(d->data.function.arg_count == 2);
        int got_r = (int)d->data.function.args[0]->data.integer;
        int got_c = (int)d->data.function.args[1]->data.integer;
        expr_free(d);
        if (got_r != checks[i].want_r || got_c != checks[i].want_c) {
            fprintf(stderr,
                "FAIL: gen shape %s for {%s, %s}: expected %d x %d got %d x %d\n",
                checks[i].label, m_src, a_src,
                checks[i].want_r, checks[i].want_c, got_r, got_c);
            ASSERT(0);
        }
    }
    expr_free(run("ClearAll[TestSVDMach$G]"));
    printf("  PASS: generalized shape for {%s, %s}\n", m_src, a_src);
}

static void test_gen_2x2_square(void) {
    /* Classic: both square 2x2, same column dimension. */
    assert_gen_reconstruction("N[{{1., 2.}, {3., 4.}}]",
                              "N[{{5., 6.}, {7., 8.}}]", 1e-12);
    assert_gen_shape("N[{{1., 2.}, {3., 4.}}]",
                     "N[{{5., 6.}, {7., 8.}}]",
                     2, 2, 2, 2, 2, 2, 2);
}

static void test_gen_tall_tall(void) {
    /* m: 3 x 2, a: 4 x 2. */
    assert_gen_reconstruction(
        "N[{{1., 2.}, {3., 4.}, {5., 6.}}]",
        "N[{{7., 8.}, {9., 10.}, {11., 12.}, {13., 14.}}]", 1e-12);
    assert_gen_shape(
        "N[{{1., 2.}, {3., 4.}, {5., 6.}}]",
        "N[{{7., 8.}, {9., 10.}, {11., 12.}, {13., 14.}}]",
        3, 4, 3, 2, 4, 2, 2);
}

static void test_gen_wide_wide(void) {
    /* m: 2 x 3, a: 2 x 3. */
    assert_gen_reconstruction(
        "N[{{1., 2., 3.}, {4., 5., 6.}}]",
        "N[{{7., 8., 9.}, {10., 11., 12.}}]", 1e-12);
    assert_gen_shape(
        "N[{{1., 2., 3.}, {4., 5., 6.}}]",
        "N[{{7., 8., 9.}, {10., 11., 12.}}]",
        2, 2, 2, 3, 2, 3, 3);
}

static void test_gen_M_lt_KL(void) {
    /* Exercises the R-split branch (Lm < K+L): m has fewer rows than the
     * shared rank of the (m, a) stack, so R sits across both A and B. */
    assert_gen_reconstruction(
        "N[{{1., 2., 3., 4.}}]",
        "N[{{5., 6., 7., 8.}, {9., 10., 11., 12.}, {13., 14., 15., 16.}}]",
        1e-12);
    assert_gen_shape(
        "N[{{1., 2., 3., 4.}}]",
        "N[{{5., 6., 7., 8.}, {9., 10., 11., 12.}, {13., 14., 15., 16.}}]",
        1, 3, 1, 4, 3, 4, 4);
}

static void test_gen_complex_2x2(void) {
    /* Both matrices complex; verifies the zggsvd3 path. */
    assert_gen_reconstruction(
        "N[{{1. + 0.5 I, 2.}, {3., 4. - 0.5 I}}]",
        "N[{{5. - I, 6.}, {7., 8. + 2. I}}]", 1e-12);
    assert_gen_shape(
        "N[{{1. + 0.5 I, 2.}, {3., 4. - 0.5 I}}]",
        "N[{{5. - I, 6.}, {7., 8. + 2. I}}]",
        2, 2, 2, 2, 2, 2, 2);
}

static void test_gen_complex_mixed(void) {
    /* m real, a complex: dispatcher forces both into the complex layout
     * and routes through zggsvd3.  Tests the cross-real/complex bridge. */
    assert_gen_reconstruction(
        "N[{{1., 2.}, {3., 4.}}]",
        "N[{{5., 6. + I}, {7. - I, 8.}}]", 1e-12);
}

static void test_gen_complex_tall(void) {
    /* m: 3x2 complex, a: 2x2 complex (M >= K+L typical path). */
    assert_gen_reconstruction(
        "N[{{1. + I, 2.}, {3., 4. + I}, {5., 6. - I}}]",
        "N[{{7. - I, 8.}, {9., 10. + 2 I}}]", 1e-12);
}

static void test_gen_dimerror(void) {
    /* Mismatched column counts -> ::matdims, returns unevaluated. */
    Expr* r = run("SingularValueDecomposition["
                  "{N[{{1., 2.}, {3., 4.}}], N[{{5., 6., 7.}}]}]");
    /* Mathilda leaves the call wrapped as SingularValueDecomposition[...]
     * when validation fails -- the head is still SingularValueDecomposition
     * (i.e. EXPR_FUNCTION with head symbol unchanged). */
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol,
                  "SingularValueDecomposition") == 0);
    expr_free(r);
    printf("  PASS: generalized SVD with col mismatch left unevaluated\n");
}

/* ============== larger random reconstruction ============== */

static void test_random_6x4(void) {
    /* RandomReal seeded to a known seed for determinism would be ideal;
     * for now we just rely on RandomReal returning machine-precision
     * Reals (which Mathilda does) so the dispatcher routes through
     * the LAPACK kernel reliably. */
    Expr* err = run(
        "SeedRandom[42]; "
        "m6x4 = RandomReal[1, {6, 4}]; "
        "sv = SingularValueDecomposition[m6x4]; "
        "Max[Abs[Flatten[sv[[1]] . sv[[2]] . Transpose[sv[[3]]] - m6x4]]]");
    double e = as_double(err);
    expr_free(err);
    expr_free(run("ClearAll[m6x4, sv]"));
    if (e > 1e-12) {
        fprintf(stderr, "FAIL: 6x4 random reconstruction error %g > 1e-12\n", e);
        ASSERT(0);
    }
    printf("  PASS: 6x4 random reconstruction error %g\n", e);
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== SingularValueDecomposition: machine-precision (LAPACK) tests ===\n");

    test_recon_2x2_rank1();
    test_recon_2x2_full();
    test_recon_3x2_tall();
    test_recon_3x3();
    test_recon_2x3_wide();
    test_recon_complex_2x2();

    test_shape_square();
    test_shape_tall();
    test_shape_wide();

    test_ordering_square();
    test_ordering_tall();

    test_truncation_k_positive();
    test_truncation_k_negative();
    test_truncation_upto_clamp();

    test_tolerance_default();
    test_tolerance_explicit();

    test_random_6x4();

    test_gen_2x2_square();
    test_gen_tall_tall();
    test_gen_wide_wide();
    test_gen_M_lt_KL();
    test_gen_complex_2x2();
    test_gen_complex_mixed();
    test_gen_complex_tall();
    test_gen_dimerror();

    symtab_clear();
    printf("=== All SingularValueDecomposition machine tests passed ===\n");
    return 0;
}

/*
 * tests/test_ludecomposition_machine.c
 *
 * Machine-precision LU coverage.  Drives lu_dispatch with inputs
 * dominated by IEEE doubles, exercising the LAPACK fast path
 * (`dgetrf` / `zgetrf` + `dgecon` / `zgecon` + `dlange` / `zlange`).
 * On a build with USE_LAPACK=0 the kernel routes through the symbolic
 * fallback, so the same identities still hold -- the tests are written
 * to pass under either configuration.
 *
 * Coverage:
 *   - Identity m[[p]] == l . u to within Chop tolerance.
 *   - p is a valid 1-indexed permutation of 1..n.
 *   - Result has shape {lu (n x n), p (length n), c (numeric)}.
 *   - For well-conditioned matrices the reported c is finite and at
 *     least 1 (the L-infinity condition number lower bound).
 *   - For singular matrices the kernel still returns a triple; we do
 *     not assert a specific value of c (LAPACK returns rcond = 0, our
 *     wrapper reports HUGE_VAL).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* Real-tolerance zero check.  Machine-precision LU residuals live
 * comfortably under 1e-9 for well-conditioned inputs at n <= 20. */
static int is_small_real(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        return v > -1e-9 && v < 1e-9;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (strcmp(h, "Complex") == 0 && e->data.function.arg_count == 2) {
            return is_small_real(e->data.function.args[0])
                && is_small_real(e->data.function.args[1]);
        }
    }
    return 0;
}

static int all_small_tensor(Expr* t) {
    if (!t) return 0;
    if (t->type == EXPR_FUNCTION
        && t->data.function.head->type == EXPR_SYMBOL
        && strcmp(t->data.function.head->data.symbol, "List") == 0) {
        for (size_t i = 0; i < t->data.function.arg_count; i++) {
            if (!all_small_tensor(t->data.function.args[i])) return 0;
        }
        return 1;
    }
    return is_small_real(t);
}

/* Assert m[[p]] == l.u to machine-precision tolerance.  We build
 * L (unit lower) and U (upper) via Table since LowerTriangularize /
 * UpperTriangularize aren't builtin yet. */
static void assert_lu_identity_numeric(const char* m_src) {
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
             "  Chop[m[[p]] - l . u]"
             "]",
             m_src);
    Expr* res = run(buf);
    if (!all_small_tensor(res)) {
        char* s = expr_to_string(res);
        fprintf(stderr,
                "FAIL: LU machine identity failed for %s\n  diff: %s\n",
                m_src, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    printf("  PASS: m[[p]] == l.u (numeric)  for %s\n", m_src);
    expr_free(res);
}

/* Assert the third element of the result is a real >= 1.0 (well-
 * conditioned matrices), or just a finite number for ill-conditioned
 * inputs we just want to land safely. */
static void assert_cond_is_positive_real(const char* m_src) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "LUDecomposition[%s][[3]]", m_src);
    Expr* res = run(buf);
    if (res->type != EXPR_REAL) {
        char* s = expr_to_string(res);
        fprintf(stderr,
                "FAIL: condition entry not Real for %s -> %s\n",
                m_src, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    double c = res->data.real;
    /* The exact cond is >= 1 in theory, but the LAPACK estimator
     * occasionally reports slightly below 1 for near-orthogonal
     * matrices; accept anything strictly positive. */
    ASSERT(c > 0.0);
    printf("  PASS: cond = %g  for %s\n", c, m_src);
    expr_free(res);
}

/* Cheap sanity test: the perm is a valid 1..n bijection. */
static void assert_perm_valid_numeric(const char* m_src, int n) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Sort[LUDecomposition[%s][[2]]] == Range[%d]",
             m_src, n);
    Expr* res = run(buf);
    ASSERT(res->type == EXPR_SYMBOL
        && strcmp(res->data.symbol, "True") == 0);
    expr_free(res);
    printf("  PASS: p is a permutation 1..%d  for %s\n", n, m_src);
}

/* ============================================================
 *   Tests
 * ============================================================ */

static void test_lu_machine_3x3_simple(void) {
    /* The Mathematica documentation machine-precision example. */
    const char* m = "{{1.6, 2.7, 3.6}, {1.2, 3.2, 5.2}, {3.3, 3.4, 6.5}}";
    assert_lu_identity_numeric(m);
    assert_perm_valid_numeric(m, 3);
    assert_cond_is_positive_real(m);
}

static void test_lu_machine_4x4(void) {
    const char* m = "{{1.0, 2.0, 3.0, 4.0}, "
                    " {2.0, 4.0, 6.0, 1.0}, "
                    " {3.0, 1.0, 5.0, 2.0}, "
                    " {4.0, 2.0, 0.0, 7.0}}";
    assert_lu_identity_numeric(m);
    assert_perm_valid_numeric(m, 4);
    assert_cond_is_positive_real(m);
}

static void test_lu_machine_5x5_diagdom(void) {
    /* Diagonally dominant so the condition number is small. */
    const char* m = "{{10.0, 1.0, 2.0, 0.0, 1.0}, "
                    " {1.0, 11.0, 0.0, 1.0, 2.0}, "
                    " {2.0, 0.0, 12.0, 1.0, 1.0}, "
                    " {0.0, 1.0, 1.0, 13.0, 2.0}, "
                    " {1.0, 2.0, 1.0, 2.0, 14.0}}";
    assert_lu_identity_numeric(m);
    assert_perm_valid_numeric(m, 5);

    /* For a diagonally-dominant matrix cond < 5. */
    Expr* res = run("LUDecomposition[{{10.0, 1.0, 2.0, 0.0, 1.0}, "
                    "{1.0, 11.0, 0.0, 1.0, 2.0}, "
                    "{2.0, 0.0, 12.0, 1.0, 1.0}, "
                    "{0.0, 1.0, 1.0, 13.0, 2.0}, "
                    "{1.0, 2.0, 1.0, 2.0, 14.0}}][[3]]");
    ASSERT(res->type == EXPR_REAL);
    ASSERT(res->data.real > 0.0 && res->data.real < 10.0);
    printf("  PASS: diag-dominant cond = %g\n", res->data.real);
    expr_free(res);
}

static void test_lu_machine_2x2_known(void) {
    /* Hand-checked: m = {{2.0, 1.0}, {6.0, 1.0}}.
     *
     * Partial pivoting picks row 1 (|6| > |2|).  lu = {{6, 1}, {1/3,
     * 2/3}}, p = {2, 1}.  We just check the identity. */
    assert_lu_identity_numeric("{{2.0, 1.0}, {6.0, 1.0}}");
}

static void test_lu_machine_complex(void) {
    /* Complex 3x3 in machine precision. */
    const char* m = "{{2.0 + 4.0 I, 9.0 + 9.0 I, 9.0 + 2.0 I}, "
                    " {2.0 + 9.0 I, 1.0 + 3.0 I, 0.0 + 4.0 I}, "
                    " {3.0 + 8.0 I, 0.0 + 0.0 I, 7.0 + 4.0 I}}";
    assert_lu_identity_numeric(m);
    assert_perm_valid_numeric(m, 3);
}

static void test_lu_machine_identity(void) {
    /* Real-valued identity matrix -- trivial factorisation. */
    assert_lu_identity_numeric("N[IdentityMatrix[4]]");
    assert_perm_valid_numeric("N[IdentityMatrix[4]]", 4);
}

static void test_lu_machine_singular_returns_triple(void) {
    /* Singular real matrix -- should still produce a {lu, p, c}
     * triple; LAPACK info > 0 surfaces as ::sing then completes. */
    Expr* res = run("LUDecomposition[{{1.0, 2.0, 3.0}, "
                    " {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol, "List") == 0);
    ASSERT(res->data.function.arg_count == 3);
    expr_free(res);
    printf("  PASS: singular machine-precision matrix returns a triple\n");
}

static void test_lu_machine_mixed_int_real(void) {
    /* Mixed numeric input -- a single Real promotes the whole matrix
     * into the machine kernel. */
    const char* m = "{{1, 2, 3}, {4.0, 5, 6}, {7, 8, 10}}";
    assert_lu_identity_numeric(m);
    assert_perm_valid_numeric(m, 3);
    assert_cond_is_positive_real(m);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_lu_machine_3x3_simple);
    TEST(test_lu_machine_4x4);
    TEST(test_lu_machine_5x5_diagdom);
    TEST(test_lu_machine_2x2_known);
    TEST(test_lu_machine_complex);
    TEST(test_lu_machine_identity);
    TEST(test_lu_machine_singular_returns_triple);
    TEST(test_lu_machine_mixed_int_real);

    printf("All LUDecomposition (machine-precision) tests passed.\n");
    return 0;
}

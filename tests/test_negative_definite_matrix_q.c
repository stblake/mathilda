/* Unit tests for NegativeDefiniteMatrixQ.
 *
 * A matrix m is negative definite iff Re[Conjugate[x] . m . x] < 0 for
 * every nonzero vector x, equivalently iff -m is positive definite,
 * i.e. the negated Hermitian part -(m + ConjugateTranspose[m]) / 2
 * admits a Cholesky factorisation.
 *
 * These tests exercise the *evaluator* path (parse -> evaluate -> print)
 * so attribute application is in scope.  Coverage:
 *
 *   - Real symmetric numeric matrices (small + negated identity).
 *   - Real asymmetric numeric matrices (Hermitian part is what matters).
 *   - Complex Hermitian numeric matrices.
 *   - Diagonal matrices with negative / non-negative entries.
 *   - Cholesky-built -m . Transpose[m] ND matrices of varied sizes.
 *   - Rejection of indefinite, positive-definite, and zero matrices.
 *   - Mixed numeric kinds (Integer / Real / Rational / BigInt) all work.
 *   - Symbolic matrices return False (we do not prove symbolic ND).
 *   - Shape/structure rejections: non-list, non-square, ragged, empty,
 *     vectors, 3-D tensors, scalars.
 *   - Diagnostics: argx for 0 args or >= 2 args (call left unevaluated).
 *   - Attribute introspection (Protected, NOT Listable) and docstring.
 *   - sym_names.c interning of SYM_NegativeDefiniteMatrixQ.
 *   - Cross-check with PositiveDefiniteMatrixQ: NegativeDefiniteMatrixQ[m]
 *     should agree with PositiveDefiniteMatrixQ[-m] on numeric input.
 *   - Repeated evaluation under a tight loop to catch double-frees /
 *     dangling-pointer regressions under valgrind.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

/* --- Small real-symmetric negative-definite matrices ------------- */

static void test_2x2_real_symmetric(void) {
    /* Documented example: {{-5, 1}, {1, -4}} is ND (det = 19, diag < 0). */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-5, 1}, {1, -4}}]",
                   "True", 0);
}

static void test_2x2_real_asymmetric(void) {
    /* Documented example: {{-2.3, -1.2}, {0.6, -3.7}} -- the Hermitian
     * part H = (m + m^T)/2 = {{-2.3, -0.3}, {-0.3, -3.7}} is ND
     * (-H is PD; det(-H) = 8.42 > 0, diag(-H) > 0). */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-2.3, -1.2}, {0.6, -3.7}}]",
                   "True", 0);
}

static void test_negated_identity(void) {
    /* -I_n is ND for all n: -I = -1 * IdentityMatrix[n]. */
    assert_eval_eq("NegativeDefiniteMatrixQ[-IdentityMatrix[1]]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[-IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[-IdentityMatrix[5]]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[-IdentityMatrix[10]]", "True", 0);
    /* Identity matrices are NOT negative definite. */
    assert_eval_eq("NegativeDefiniteMatrixQ[IdentityMatrix[3]]", "False", 0);
}

static void test_1x1_negative(void) {
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1}}]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-42}}]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-0.001}}]", "True", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1/3}}]", "True", 0);
}

static void test_1x1_zero_or_positive(void) {
    /* {{0}} is negative semi-definite but NOT negative definite. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{0}}]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{1}}]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{0.0001}}]", "False", 0);
}

/* --- Complex Hermitian / asymmetric ------------------------------- */

static void test_2x2_complex_hermitian(void) {
    /* Documented example: {{-1, 2I}, {-I, -4}} -- non-Hermitian, but
     * the Hermitian part (m + m^H)/2 = {{-1, 3I/2}, {-3I/2, -4}}.
     * -H = {{1, -3I/2}, {3I/2, 4}} has det = 4 - 9/4 = 7/4 > 0 and
     * positive diag, so -H is PD => H (and hence m) is ND. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, 2*I}, {-I, -4}}]",
                   "True", 0);
}

static void test_2x2_complex_hermitian_strict(void) {
    /* A strictly Hermitian ND example: {{-2, 1-I}, {1+I, -3}}.
     * -m = {{2, -(1-I)}, {-(1+I), 3}} which is Hermitian PD:
     *   minor1 = 2 > 0; det = 6 - |1-I|^2 = 6 - 2 = 4 > 0. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-2, 1 - I}, {1 + I, -3}}]",
                   "True", 0);
}

static void test_complex_mixed_diagonal(void) {
    /* Diagonal m = {{-1, 0}, {0, 2}} is indefinite: one negative, one
     * positive entry => not ND. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, 0}, {0, 2}}]",
                   "False", 0);
    /* Documented: m = {{-1, -2I}, {I, -4}} -- Hermitian part
     * (m + m^H)/2 = {{-1, (-2I - I)/2}, {(I - (-2I))/2 conj ..., -4}}
     * Let's compute: m[0,1] = -2I, m[1,0]^* = (-I)^* = I.  H[0,1] =
     * (-2I + I)/2 = -I/2.  H = {{-1, -I/2}, {I/2, -4}}.  -H is PD iff
     * det(-H) = 4 - 1/4 = 15/4 > 0 and positive diag => True. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, -2*I}, {I, -4}}]",
                   "True", 0);
}

/* --- 3x3 numeric matrices ----------------------------------------- */

static void test_3x3_real_not_nd(void) {
    /* Documented example with Pi, E, ... -- numerically indefinite. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3},"
        " {5, Sqrt[2], 5}}]",
        "False", 0);
}

static void test_3x3_diag_negative(void) {
    /* Documented: diagonal m is ND iff every diagonal entry has
     * negative real part. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1, 0, 0}, {0, -2, 0}, {0, 0, -4}}]",
        "True", 0);
    /* Rational negative diagonal. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1/2, 0, 0}, {0, -1/4, 0}, {0, 0, -7}}]",
        "True", 0);
}

static void test_3x3_diag_zero_or_mixed(void) {
    /* Documented: a zero on the diagonal => not ND (only semi-definite). */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1, 0, 0}, {0, -2, 0}, {0, 0, 0}}]",
        "False", 0);
    /* Mixed sign: positive entry => indefinite => not ND. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1, 0, 0}, {0, 1, 0}, {0, 0, -2}}]",
        "False", 0);
}

/* --- Matrix products -m . Transpose[m] are ND when m has full rank - */

static void test_mtm_nd_real(void) {
    /* For any real m with full row rank, m . m^T is PD, so
     * -(m . m^T) is ND.  Use a small fixed matrix for determinism. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ["
        "  Module[{m = {{1.0, 2.0, 3.0}, {4.0, 5.0, 7.0}, {8.0, 1.0, 2.0}}},"
        "    -(m . Transpose[m])]]",
        "True", 0);
}

static void test_mtm_nd_complex(void) {
    /* -(m . ConjugateTranspose[m]) is Hermitian ND for any full-rank
     * complex m. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ["
        "  Module[{m = {{1 + I, 2, 3 - I}, {0, 1, 2 + I}, {4 - I, 0, 1}}},"
        "    -(m . ConjugateTranspose[m])]]",
        "True", 0);
}

/* --- Negated Hilbert matrix (ND for all n) ------------------------- */

static void test_negated_hilbert_5(void) {
    /* -H[i,j] = -1/(i + j - 1) gives a negative Hilbert matrix, ND for
     * all n.  Cap at n = 8 where double precision still resolves it. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[Table[-1/(i + j - 1), {i, 5}, {j, 5}]]",
        "True", 0);
}

static void test_negated_hilbert_8(void) {
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[Table[-1/(i + j - 1), {i, 8}, {j, 8}]]",
        "True", 0);
}

/* --- BigInt / mixed numeric kinds --------------------------------- */

static void test_bigint_diagonal(void) {
    /* -10^25 exceeds int64 (negated) -- forces BigInt path. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-(10^25), 0}, {0, -(10^25)}}]",
        "True", 0);
    /* Mixed sign with bigint => not ND. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{10^25, 0}, {0, -(10^25)}}]",
        "False", 0);
}

static void test_mixed_numeric_kinds(void) {
    /* Integer / Real / Rational all in the same ND matrix.  Hermitian
     * part = {{-2, -1/4 + 1/4}, ...} = {{-2, 0}, {0, -3.0}} after
     * symmetrising the off-diagonals -- diagonal-dominant negative. */
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-2, 1/2}, {-1/2, -3.0}}]",
        "True", 0);
}

/* --- Symbolic / non-numeric input --------------------------------- */

static void test_symbolic_returns_false(void) {
    /* Documented: {{-1, a}, {b, -2}} -- symbolic entries => return
     * False (we don't prove symbolic ND). */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, a}, {b, -2}}]",
                   "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{a, b}, {c, d}}]",
                   "False", 0);
    /* Mixed numeric / symbolic. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, x}, {x, -2}}]",
                   "False", 0);
}

static void test_pure_symbol_returns_false(void) {
    /* A bare symbol (not a list) -- predicate still returns False. */
    assert_eval_eq("NegativeDefiniteMatrixQ[m]", "False", 0);
    /* Documented: Sqrt[3] is a scalar, not a matrix -- returns False. */
    assert_eval_eq("NegativeDefiniteMatrixQ[Sqrt[3]]", "False", 0);
}

/* --- Shape / structure rejections --------------------------------- */

static void test_non_matrix_inputs(void) {
    assert_eval_eq("NegativeDefiniteMatrixQ[5]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[-5]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[3.14]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[\"foo\"]", "False", 0);
    /* Vectors (1-D lists). */
    assert_eval_eq("NegativeDefiniteMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{-1, -2, -3}]", "False", 0);
    /* Empty list. */
    assert_eval_eq("NegativeDefiniteMatrixQ[{}]", "False", 0);
}

static void test_non_square(void) {
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, -2, -3}}]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[{{-1}, {-2}, {-3}}]", "False", 0);
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1, -2, -3}, {-4, -5, -6}}]", "False", 0);
}

static void test_ragged(void) {
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{-1, -2}, {-3, -4, -5}}]", "False", 0);
}

static void test_three_d_tensor(void) {
    assert_eval_eq(
        "NegativeDefiniteMatrixQ[{{{-1, -2}, {-3, -4}}}]", "False", 0);
}

/* --- Diagnostics (argx) ------------------------------------------- */

static void test_zero_args_emits_argx(void) {
    /* No args -> evaluator should not crash; result keeps the head. */
    const char* in = "NegativeDefiniteMatrixQ[]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "NegativeDefiniteMatrixQ") != NULL,
               "expected unevaluated NegativeDefiniteMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_two_args_emits_argx(void) {
    /* Documented diagnostic: NegativeDefiniteMatrixQ[1,2] left
     * unevaluated. */
    const char* in = "NegativeDefiniteMatrixQ[1, 2]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "NegativeDefiniteMatrixQ") != NULL,
               "expected unevaluated NegativeDefiniteMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_three_args_emits_argx(void) {
    const char* in = "NegativeDefiniteMatrixQ[1, 2, 3]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "NegativeDefiniteMatrixQ") != NULL,
               "expected unevaluated NegativeDefiniteMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Cross-checks against negated PD matrices --------------------- */

static void test_cross_check_nd_battery(void) {
    /* Each test case (m, expected): explicit reasoning in the comment. */
    struct {
        const char* matrix;
        const char* expected;
    } cases[] = {
        /* Strictly diagonally dominant negative symmetric => ND. */
        { "{{-10, 1, 2}, {1, -10, 1}, {2, 1, -10}}",                "True"  },
        /* Off-diagonals dominate => not ND. */
        { "{{-1, 2, 0}, {2, -1, 0}, {0, 0, -1}}",                   "False" },
        /* Negated tridiagonal SPD. */
        { "{{-4, 1, 0}, {1, -4, 1}, {0, 1, -4}}",                   "True"  },
        /* Zero matrix is negative semi-definite, NOT negative definite. */
        { "{{0, 0}, {0, 0}}",                                       "False" },
        /* Standard ND example with rational entries. */
        { "{{-2, 1/2}, {1/2, -2}}",                                 "True"  },
        /* PD example is NOT ND. */
        { "{{5, -1}, {-1, 4}}",                                     "False" },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "NegativeDefiniteMatrixQ[%s]", cases[i].matrix);
        assert_eval_eq(buf, cases[i].expected, 0);
    }
}

/* --- Cross-check: NegativeDefiniteMatrixQ[m] == PositiveDefiniteMatrixQ[-m]
 * for any numeric matrix m.  This is the *definition* of negative
 * definiteness so the two predicates must agree on all numeric inputs.
 * Symbolic inputs collapse to "False == False" trivially. */
static void test_cross_check_against_pd(void) {
    const char* matrices[] = {
        "{{-5, 1}, {1, -4}}",
        "{{5, -1}, {-1, 4}}",
        "{{-1, 0}, {0, -2}}",
        "{{1, 0}, {0, 2}}",
        "{{-1, 2*I}, {-2*I, -3}}",
        "{{-1, 0, 0}, {0, -2, 0}, {0, 0, -4}}",
        "{{0, 0}, {0, 0}}",
        "{{-1, 0, 0}, {0, 0, 0}, {0, 0, -1}}",
        "{{-2, 1/2}, {1/2, -2}}",
    };
    for (size_t i = 0; i < sizeof(matrices)/sizeof(matrices[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "NegativeDefiniteMatrixQ[%s] == PositiveDefiniteMatrixQ[-(%s)]",
                 matrices[i], matrices[i]);
        assert_eval_eq(buf, "True", 0);
    }
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("NegativeDefiniteMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* Should NOT be Listable -- input is a matrix, not a list of matrices. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("NegativeDefiniteMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "negative definite") != NULL);
    ASSERT(strstr(def->docstring, "Cholesky") != NULL);
}

static void test_sym_pointer_interned(void) {
    /* SYM_NegativeDefiniteMatrixQ must be initialised by sym_names_init()
     * and stably point to the interned name. */
    ASSERT(SYM_NegativeDefiniteMatrixQ != NULL);
    ASSERT(strcmp(SYM_NegativeDefiniteMatrixQ,
                  "NegativeDefiniteMatrixQ") == 0);
}

/* --- Repeated evaluation: leak / double-free sanity --------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Run a mix of accept / reject / shape-reject / symbolic cases
     * many times; any allocator misuse (double-free, leak, dangling
     * pointer) will typically surface as a crash, wrong answer, or
     * valgrind diagnostic by the second or third iteration. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("NegativeDefiniteMatrixQ[{{-5, 1}, {1, -4}}]",
                       "True", 0);
        assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, -2*I}, {I, -4}}]",
                       "True", 0);
        assert_eval_eq("NegativeDefiniteMatrixQ[{{-1, -2, -3}, {-4, -5, -6}}]",
                       "False", 0);
        assert_eval_eq("NegativeDefiniteMatrixQ[{{a, b}, {c, d}}]",
                       "False", 0);
        assert_eval_eq("NegativeDefiniteMatrixQ[-IdentityMatrix[4]]",
                       "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_2x2_real_symmetric);
    TEST(test_2x2_real_asymmetric);
    TEST(test_negated_identity);
    TEST(test_1x1_negative);
    TEST(test_1x1_zero_or_positive);

    TEST(test_2x2_complex_hermitian);
    TEST(test_2x2_complex_hermitian_strict);
    TEST(test_complex_mixed_diagonal);

    TEST(test_3x3_real_not_nd);
    TEST(test_3x3_diag_negative);
    TEST(test_3x3_diag_zero_or_mixed);

    TEST(test_mtm_nd_real);
    TEST(test_mtm_nd_complex);

    TEST(test_negated_hilbert_5);
    TEST(test_negated_hilbert_8);

    TEST(test_bigint_diagonal);
    TEST(test_mixed_numeric_kinds);

    TEST(test_symbolic_returns_false);
    TEST(test_pure_symbol_returns_false);

    TEST(test_non_matrix_inputs);
    TEST(test_non_square);
    TEST(test_ragged);
    TEST(test_three_d_tensor);

    TEST(test_zero_args_emits_argx);
    TEST(test_two_args_emits_argx);
    TEST(test_three_args_emits_argx);

    TEST(test_cross_check_nd_battery);
    TEST(test_cross_check_against_pd);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All NegativeDefiniteMatrixQ tests passed!\n");
    return 0;
}

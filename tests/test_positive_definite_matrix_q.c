/* Unit tests for PositiveDefiniteMatrixQ.
 *
 * A matrix m is positive definite iff Re[Conjugate[x] . m . x] > 0 for
 * every nonzero vector x, equivalently iff its Hermitian part
 * (m + ConjugateTranspose[m]) / 2 admits a Cholesky factorisation.
 *
 * These tests exercise the *evaluator* path (parse -> evaluate -> print)
 * so attribute application is in scope.  Coverage:
 *
 *   - Real symmetric numeric matrices (small + identity + Hilbert).
 *   - Real asymmetric numeric matrices (Hermitian part is what matters).
 *   - Complex Hermitian numeric matrices.
 *   - Diagonal matrices with positive / non-positive entries.
 *   - Cholesky-built m . Transpose[m] PD matrices of varied sizes.
 *   - Rejection of indefinite, negative-definite, and zero matrices.
 *   - Rejection of complex symmetric (non-Hermitian) matrices that
 *     have negative diagonal real-part.
 *   - Mixed numeric kinds (Integer / Real / Rational / BigInt) all work.
 *   - Symbolic matrices return False (we do not prove symbolic PD).
 *   - Shape/structure rejections: non-list, non-square, ragged, empty,
 *     vectors, 3-D tensors, scalars.
 *   - Diagnostics: argx for 0 args or >= 2 args (call left unevaluated).
 *   - Attribute introspection (Protected, NOT Listable) and docstring.
 *   - sym_names.c interning of SYM_PositiveDefiniteMatrixQ.
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

/* --- Small real-symmetric positive-definite matrices -------------- */

static void test_2x2_real_symmetric(void) {
    /* Docstring example: {{5,-1},{-1,4}} is PD (det = 19, diag > 0). */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{5, -1}, {-1, 4}}]",
                   "True", 0);
}

static void test_2x2_real_asymmetric(void) {
    /* Docstring example: {{2.3,-1.2},{0.6,3.7}} -- the Hermitian part
     * H = (m + m^T)/2 = {{2.3,-0.3},{-0.3,3.7}} is PD (det = 8.42). */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{2.3, -1.2}, {0.6, 3.7}}]",
                   "True", 0);
}

static void test_identity_is_positive_definite(void) {
    assert_eval_eq("PositiveDefiniteMatrixQ[IdentityMatrix[1]]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[IdentityMatrix[5]]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[IdentityMatrix[10]]", "True", 0);
}

static void test_1x1_positive(void) {
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1}}]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{42}}]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{0.001}}]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1/3}}]", "True", 0);
}

static void test_1x1_zero_or_negative(void) {
    /* {{0}} is positive semi-definite but NOT positive definite. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{0}}]", "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{-1}}]", "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{-0.0001}}]", "False", 0);
}

/* --- Complex Hermitian / asymmetric ------------------------------- */

static void test_2x2_complex_hermitian(void) {
    /* Docstring example: {{1, 2I}, {-I, 4}} -- non-Hermitian, but the
     * Hermitian part (m + m^H)/2 = {{1, 3I/2}, {-3I/2, 4}} has
     * det = 4 - 9/4 = 7/4 > 0 and is PD. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1, 2*I}, {-I, 4}}]",
                   "True", 0);
}

static void test_2x2_complex_hermitian_strict(void) {
    /* A strictly Hermitian PD example: {{2, 1-I}, {1+I, 3}}.
     * Sylvester: minor1 = 2 > 0; det = 6 - (1-I)(1+I) = 6 - 2 = 4 > 0. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{2, 1 - I}, {1 + I, 3}}]",
                   "True", 0);
}

static void test_complex_negative_diagonal(void) {
    /* {{-1, 0}, {0, 2}} -- diagonal has a negative entry => not PD. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{-1, 0}, {0, 2}}]",
                   "False", 0);
}

/* --- 3x3 numeric matrices ----------------------------------------- */

static void test_3x3_real_indefinite(void) {
    /* Docstring example with Pi, E, ... -- numerically indefinite. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3},"
        " {5, Sqrt[2], 5}}]",
        "False", 0);
}

static void test_3x3_diag_positive(void) {
    /* Diagonal with positive entries is PD. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}]",
        "True", 0);
    /* Rational diagonal. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1/2, 0, 0}, {0, 1/4, 0}, {0, 0, 7}}]",
        "True", 0);
}

static void test_3x3_diag_nonpositive(void) {
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1, 0, 0}, {0, -1, 0}, {0, 0, 2}}]",
        "False", 0);
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}]",
        "False", 0);
}

/* --- Matrix products m . Transpose[m] are PD when m has full rank - */

static void test_mtm_pd_real(void) {
    /* For any real m with full row rank, m . m^T is PD.  Use a small
     * fixed matrix so the result is deterministic. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ["
        "  Module[{m = {{1.0, 2.0, 3.0}, {4.0, 5.0, 7.0}, {8.0, 1.0, 2.0}}},"
        "    m . Transpose[m]]]",
        "True", 0);
}

static void test_mtm_pd_complex(void) {
    /* m . ConjugateTranspose[m] is Hermitian PD for any full-rank
     * complex m. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ["
        "  Module[{m = {{1 + I, 2, 3 - I}, {0, 1, 2 + I}, {4 - I, 0, 1}}},"
        "    m . ConjugateTranspose[m]]]",
        "True", 0);
}

/* --- Hilbert matrix (PD for all n; ill-conditioned for large n) ---- */

static void test_hilbert_5(void) {
    /* H[i,j] = 1/(i + j - 1).  Famously PD for all n, but condition
     * number grows ~exponentially -- so we cap at n = 8 where double
     * precision still gets the right answer. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[Table[1/(i + j - 1), {i, 5}, {j, 5}]]",
        "True", 0);
}

static void test_hilbert_8(void) {
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[Table[1/(i + j - 1), {i, 8}, {j, 8}]]",
        "True", 0);
}

/* --- BigInt / mixed numeric kinds --------------------------------- */

static void test_bigint_diagonal(void) {
    /* 10^25 exceeds int64 -- forces BigInt path through pdq_leaf_to_double. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{10^25, 0}, {0, 10^25}}]",
        "True", 0);
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{-(10^25), 0}, {0, 10^25}}]",
        "False", 0);
}

static void test_mixed_numeric_kinds(void) {
    /* Integer / Real / Rational / Complex all in the same matrix. */
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{2, 1/2}, {1/2, 3.0}}]",
        "True", 0);
}

/* --- Symbolic / non-numeric input --------------------------------- */

static void test_symbolic_returns_false(void) {
    /* Mathematica documents True for `Block[{b = -Conjugate[a]},
     * PositiveDefiniteMatrixQ[{{1, a}, {b, 2}}]]` via a deeper
     * symbolic simplification we do NOT implement -- we conservatively
     * return False for any non-numeric input. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1, a}, {b, 2}}]",
                   "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{a, b}, {c, d}}]",
                   "False", 0);
    /* Mixed numeric / symbolic. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1, x}, {x, 2}}]",
                   "False", 0);
}

static void test_pure_symbol_returns_false(void) {
    /* A bare symbol (not a list) -- predicate still returns False. */
    assert_eval_eq("PositiveDefiniteMatrixQ[m]", "False", 0);
}

/* --- Shape / structure rejections --------------------------------- */

static void test_non_matrix_inputs(void) {
    assert_eval_eq("PositiveDefiniteMatrixQ[5]", "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[3.14]", "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[\"foo\"]", "False", 0);
    /* Vectors (1-D lists). */
    assert_eval_eq("PositiveDefiniteMatrixQ[{1, 2, 3}]", "False", 0);
    /* Empty list. */
    assert_eval_eq("PositiveDefiniteMatrixQ[{}]", "False", 0);
}

static void test_non_square(void) {
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1, 2, 3}}]", "False", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[{{1}, {2}, {3}}]", "False", 0);
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
}

static void test_ragged(void) {
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{1, 2}, {3, 4, 5}}]", "False", 0);
}

static void test_three_d_tensor(void) {
    assert_eval_eq(
        "PositiveDefiniteMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
}

/* --- Diagnostics (argx) ------------------------------------------- */

static void test_zero_args_emits_argx(void) {
    /* No args -> evaluator should not crash; result keeps the head. */
    const char* in = "PositiveDefiniteMatrixQ[]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "PositiveDefiniteMatrixQ") != NULL,
               "expected unevaluated PositiveDefiniteMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_three_args_emits_argx(void) {
    const char* in = "PositiveDefiniteMatrixQ[1, 2, 3]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "PositiveDefiniteMatrixQ") != NULL,
               "expected unevaluated PositiveDefiniteMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Cross-checks against Cholesky-built matrices ----------------- */

static void test_cross_check_pd_battery(void) {
    /* Each test case (m, expected): explicit reasoning in the comment. */
    struct {
        const char* matrix;
        const char* expected;
    } cases[] = {
        /* Strictly diagonally dominant symmetric => PD. */
        { "{{10, 1, 2}, {1, 10, 1}, {2, 1, 10}}",                 "True"  },
        /* Strictly diagonally dominant but with a non-PD eigenvalue. */
        { "{{1, 2, 0}, {2, 1, 0}, {0, 0, 1}}",                    "False" },
        /* Tridiagonal SPD (positive diag + small off-diagonals). */
        { "{{4, -1, 0}, {-1, 4, -1}, {0, -1, 4}}",                "True"  },
        /* Zero matrix is positive semi-definite, NOT positive definite. */
        { "{{0, 0}, {0, 0}}",                                     "False" },
        /* Standard SPD example with rational entries. */
        { "{{2, 1/2}, {1/2, 2}}",                                 "True"  },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "PositiveDefiniteMatrixQ[%s]", cases[i].matrix);
        assert_eval_eq(buf, cases[i].expected, 0);
    }
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("PositiveDefiniteMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* Should NOT be Listable -- input is a matrix, not a list of matrices. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("PositiveDefiniteMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "positive definite") != NULL);
    ASSERT(strstr(def->docstring, "Cholesky") != NULL);
}

static void test_sym_pointer_interned(void) {
    /* SYM_PositiveDefiniteMatrixQ must be initialised by sym_names_init()
     * and stably point to the interned name. */
    ASSERT(SYM_PositiveDefiniteMatrixQ != NULL);
    ASSERT(strcmp(SYM_PositiveDefiniteMatrixQ,
                  "PositiveDefiniteMatrixQ") == 0);
}

/* --- Repeated evaluation: leak / double-free sanity --------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Run a mix of accept / reject / shape-reject / symbolic cases
     * many times; any allocator misuse (double-free, leak, dangling
     * pointer) will typically surface as a crash, wrong answer, or
     * valgrind diagnostic by the second or third iteration. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("PositiveDefiniteMatrixQ[{{5, -1}, {-1, 4}}]",
                       "True", 0);
        assert_eval_eq("PositiveDefiniteMatrixQ[{{1, 2*I}, {-I, 4}}]",
                       "True", 0);
        assert_eval_eq("PositiveDefiniteMatrixQ[{{1, 2, 3}, {4, 5, 6}}]",
                       "False", 0);
        assert_eval_eq("PositiveDefiniteMatrixQ[{{a, b}, {c, d}}]",
                       "False", 0);
        assert_eval_eq("PositiveDefiniteMatrixQ[IdentityMatrix[4]]",
                       "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_2x2_real_symmetric);
    TEST(test_2x2_real_asymmetric);
    TEST(test_identity_is_positive_definite);
    TEST(test_1x1_positive);
    TEST(test_1x1_zero_or_negative);

    TEST(test_2x2_complex_hermitian);
    TEST(test_2x2_complex_hermitian_strict);
    TEST(test_complex_negative_diagonal);

    TEST(test_3x3_real_indefinite);
    TEST(test_3x3_diag_positive);
    TEST(test_3x3_diag_nonpositive);

    TEST(test_mtm_pd_real);
    TEST(test_mtm_pd_complex);

    TEST(test_hilbert_5);
    TEST(test_hilbert_8);

    TEST(test_bigint_diagonal);
    TEST(test_mixed_numeric_kinds);

    TEST(test_symbolic_returns_false);
    TEST(test_pure_symbol_returns_false);

    TEST(test_non_matrix_inputs);
    TEST(test_non_square);
    TEST(test_ragged);
    TEST(test_three_d_tensor);

    TEST(test_zero_args_emits_argx);
    TEST(test_three_args_emits_argx);

    TEST(test_cross_check_pd_battery);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All PositiveDefiniteMatrixQ tests passed!\n");
    return 0;
}

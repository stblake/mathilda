/* Unit tests for SquareMatrixQ.
 *
 * A matrix m is square iff Dimensions[m] == {n, n}, i.e. equal number of
 * rows and columns.  SquareMatrixQ is a pure shape test; it does not
 * inspect entries beyond rejecting deeper nesting (anything whose entries
 * are themselves Lists is treated as a higher-rank tensor and rejected).
 *
 * Coverage parallels test_symmetric_matrix_q.c at the shape layer and
 * is intentionally exhaustive on the structural axis (empty / vector /
 * scalar / ragged / rectangular / 3-D tensor / single-element / large
 * symbolic / IdentityMatrix / mixed numeric+symbolic).
 *
 * SquareMatrixQ accepts a single positional argument and supports no
 * options; trailing arguments must therefore leave the call unevaluated.
 *
 *   - Numerical square matrices (small + larger + identity).
 *   - Symbolic square matrices (mixed and fully symbolic).
 *   - Cross-check against Dimensions[m] == {n, n}.
 *   - Cross-check against MatrixQ on the same battery (every square is
 *     a matrix; not every matrix is square).
 *   - Rejection of: non-list, scalar, string, vector, empty list,
 *     {{}} (1 x 0), rectangular (n x m, n != m), ragged rows, 3-D
 *     tensors, single-row, single-column.
 *   - Pure unbound symbol (m, x) returns False rather than leaving the
 *     call unevaluated.
 *   - Extra arguments (no options supported) leave the call unevaluated.
 *   - Attribute introspection (Protected; NOT Listable) and docstring
 *     presence.
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

/* --- Numerical square matrices ------------------------------------ */

static void test_square_2x2_numeric(void) {
    /* Docstring example: {{1,2},{3,4}}. */
    assert_eval_eq("SquareMatrixQ[{{1, 2}, {3, 4}}]", "True", 0);
    assert_eval_eq("SquareMatrixQ[{{1.0, 2.0}, {3.0, 4.0}}]", "True", 0);
    assert_eval_eq("SquareMatrixQ[{{0, 0}, {0, 0}}]", "True", 0);
}

static void test_square_3x3_numeric(void) {
    assert_eval_eq(
        "SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "True", 0);
}

static void test_square_1x1(void) {
    /* The smallest square: a single entry, numeric or symbolic. */
    assert_eval_eq("SquareMatrixQ[{{0}}]", "True", 0);
    assert_eval_eq("SquareMatrixQ[{{1}}]", "True", 0);
    assert_eval_eq("SquareMatrixQ[{{a}}]", "True", 0);
}

static void test_square_identity_matrix(void) {
    assert_eval_eq("SquareMatrixQ[IdentityMatrix[1]]", "True", 0);
    assert_eval_eq("SquareMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("SquareMatrixQ[IdentityMatrix[5]]", "True", 0);
}

static void test_square_larger(void) {
    /* 4x4 — large enough that the row-length loop matters. */
    assert_eval_eq(
        "SquareMatrixQ[{{1, 2, 3, 4}, {5, 6, 7, 8},"
        " {9, 10, 11, 12}, {13, 14, 15, 16}}]",
        "True", 0);
}

/* --- Symbolic square matrices ------------------------------------- */

static void test_square_symbolic_2x2(void) {
    /* All-symbolic. */
    assert_eval_eq("SquareMatrixQ[{{a, b}, {c, d}}]", "True", 0);
    /* Mixed numeric + symbolic. */
    assert_eval_eq("SquareMatrixQ[{{1, a}, {b, 2}}]", "True", 0);
}

static void test_square_symbolic_3x3(void) {
    assert_eval_eq(
        "SquareMatrixQ[{{a, b, c}, {d, e, f}, {g, h, i}}]", "True", 0);
}

static void test_square_unevaluated_expr_entries(void) {
    /* Entries that are themselves function calls (but not Lists) are
     * fine -- SquareMatrixQ is a shape test. */
    assert_eval_eq("SquareMatrixQ[{{Sin[x], Cos[x]}, {y + z, y*z}}]",
                   "True", 0);
}

/* --- Non-matrix scalar / vector inputs ---------------------------- */

static void test_scalar_inputs(void) {
    /* Atoms are never matrices. */
    assert_eval_eq("SquareMatrixQ[5]", "False", 0);
    assert_eval_eq("SquareMatrixQ[3.14]", "False", 0);
    assert_eval_eq("SquareMatrixQ[\"hello\"]", "False", 0);
    /* Pure symbol -- SquareMatrixQ resolves to False, not unevaluated. */
    assert_eval_eq("SquareMatrixQ[m]", "False", 0);
    assert_eval_eq("SquareMatrixQ[x]", "False", 0);
}

static void test_vector_inputs(void) {
    /* 1-D Lists are vectors, not matrices. */
    assert_eval_eq("SquareMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("SquareMatrixQ[{a, b, c}]", "False", 0);
    assert_eval_eq("SquareMatrixQ[{1}]", "False", 0);
}

static void test_empty_list(void) {
    /* {} has dimension {0} -- not a matrix at all. */
    assert_eval_eq("SquareMatrixQ[{}]", "False", 0);
}

static void test_single_empty_row(void) {
    /* {{}} has dimensions {1, 0} -- 1 row, 0 columns, not square. */
    assert_eval_eq("SquareMatrixQ[{{}}]", "False", 0);
    /* {{},{}} has dimensions {2, 0}. */
    assert_eval_eq("SquareMatrixQ[{{}, {}}]", "False", 0);
}

/* --- Non-square matrices ------------------------------------------ */

static void test_single_row_rectangular(void) {
    /* Docstring example: {{1,2,3}} is 1x3. */
    assert_eval_eq("SquareMatrixQ[{{1, 2, 3}}]", "False", 0);
    assert_eval_eq("SquareMatrixQ[{{1, 2}}]", "False", 0);
}

static void test_single_column_rectangular(void) {
    /* {{1},{2},{3}} is 3x1. */
    assert_eval_eq("SquareMatrixQ[{{1}, {2}, {3}}]", "False", 0);
    assert_eval_eq("SquareMatrixQ[{{a}, {b}}]", "False", 0);
}

static void test_rectangular_2x3(void) {
    assert_eval_eq("SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
}

static void test_rectangular_3x2(void) {
    assert_eval_eq(
        "SquareMatrixQ[{{1, 2}, {3, 4}, {5, 6}}]", "False", 0);
}

static void test_rectangular_3x4(void) {
    assert_eval_eq(
        "SquareMatrixQ[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]",
        "False", 0);
}

/* --- Ragged rows -------------------------------------------------- */

static void test_ragged_short_then_long(void) {
    /* Docstring example: {{1},{2,3}}.  Row 0 has 1 col, row 1 has 2. */
    assert_eval_eq("SquareMatrixQ[{{1}, {2, 3}}]", "False", 0);
}

static void test_ragged_long_then_short(void) {
    assert_eval_eq("SquareMatrixQ[{{1, 2}, {3}}]", "False", 0);
}

static void test_ragged_three_rows(void) {
    assert_eval_eq(
        "SquareMatrixQ[{{1, 2, 3}, {4, 5}, {6, 7, 8}}]", "False", 0);
}

/* --- Higher-rank tensors / nesting -------------------------------- */

static void test_three_d_tensor_rejected(void) {
    /* Outer dims may be square but entries are themselves Lists. */
    assert_eval_eq(
        "SquareMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
    assert_eval_eq(
        "SquareMatrixQ[{{{1}, {2}}, {{3}, {4}}}]", "False", 0);
}

static void test_mixed_atom_and_list_entry_rejected(void) {
    /* One entry is a list, breaking the "list of atoms in each row"
     * shape -- even though row lengths match.  SquareMatrixQ treats
     * any list entry as deeper nesting. */
    assert_eval_eq(
        "SquareMatrixQ[{{1, {2, 3}}, {4, 5}}]", "False", 0);
}

/* --- Cross-check against Dimensions ------------------------------- */

static void test_matches_dimensions_eq_nn(void) {
    /* For each matrix, the predicate should agree with
     * Dimensions[m] === {n, n}. */
    struct {
        const char* matrix;
        const char* expected;
    } cases[] = {
        { "{{1, 2}, {3, 4}}",                                  "True"  },
        { "{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}",                 "True"  },
        { "{{a, b, c}, {d, e, f}, {g, h, i}}",                 "True"  },
        { "{{1}}",                                              "True"  },
        { "IdentityMatrix[4]",                                  "True"  },
        { "{{1, 2}, {3, 4}, {5, 6}}",                          "False" },
        { "{{1, 2, 3}, {4, 5, 6}}",                            "False" },
        { "{{1, 2, 3}}",                                        "False" },
        { "{{1}, {2}, {3}}",                                    "False" },
        { "{{1}, {2, 3}}",                                      "False" },
        { "{}",                                                 "False" },
        { "{{}}",                                               "False" },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf), "SquareMatrixQ[%s]", cases[i].matrix);
        assert_eval_eq(buf, cases[i].expected, 0);
    }
}

/* --- Cross-check against MatrixQ ---------------------------------- */

static void test_implies_matrixq(void) {
    /* Every matrix where SquareMatrixQ -> True must also satisfy
     * MatrixQ -> True (square implies matrix). */
    const char* squares[] = {
        "{{1, 2}, {3, 4}}",
        "{{a, b, c}, {d, e, f}, {g, h, i}}",
        "{{1}}",
        "IdentityMatrix[3]",
    };
    for (size_t i = 0; i < sizeof(squares)/sizeof(squares[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf), "SquareMatrixQ[%s]", squares[i]);
        assert_eval_eq(buf, "True", 0);
        snprintf(buf, sizeof(buf), "MatrixQ[%s]", squares[i]);
        assert_eval_eq(buf, "True", 0);
    }
}

static void test_matrixq_not_implies_squarematrixq(void) {
    /* MatrixQ -> True does not imply SquareMatrixQ -> True. */
    assert_eval_eq("MatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "True", 0);
    assert_eval_eq("SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
}

/* --- Argument count / option handling ----------------------------- */

/* Evaluate `in`, assert the printed result still contains the head
 * "SquareMatrixQ" (i.e. the call was left unevaluated), and free
 * everything.  Used for the argx (wrong-arity) cases. */
static void assert_squarematrixq_unevaluated(const char* in) {
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "SquareMatrixQ") != NULL,
               "expected unevaluated SquareMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_zero_arg_stays_unevaluated(void) {
    /* SquareMatrixQ[] -- 0 args -> argx diagnostic + unevaluated. */
    assert_squarematrixq_unevaluated("SquareMatrixQ[]");
}

static void test_three_args_stays_unevaluated(void) {
    /* SquareMatrixQ[1, 2, 3] -- 3 args -> argx diagnostic + unevaluated.
     * Mirrors the surface behaviour the user observed in the REPL. */
    assert_squarematrixq_unevaluated("SquareMatrixQ[1, 2, 3]");
}

static void test_two_args_stays_unevaluated(void) {
    /* SquareMatrixQ takes exactly one argument; extra positional args
     * leave the call unevaluated so user typos surface. */
    assert_squarematrixq_unevaluated(
        "SquareMatrixQ[{{1, 2}, {3, 4}}, foo]");
}

static void test_rule_extra_arg_stays_unevaluated(void) {
    /* No options are defined -- even a well-formed Rule[] is rejected
     * as a wrong-arity call. */
    assert_squarematrixq_unevaluated(
        "SquareMatrixQ[{{1, 2}, {3, 4}}, Tolerance -> 0.1]");
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("SquareMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* SquareMatrixQ is NOT Listable -- input is a single matrix. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("SquareMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "square") != NULL);
    /* The docstring should reference the {n, n} shape contract. */
    ASSERT(strstr(def->docstring, "{n, n}") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_SquareMatrixQ != NULL);
    ASSERT(strcmp(SYM_SquareMatrixQ, "SquareMatrixQ") == 0);
}

/* --- Repeated evaluation: leak sanity ----------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Tight loop over a mix of accepted and rejected shapes to surface
     * any double-free / dangling-pointer regression under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
            "True", 0);
        assert_eval_eq(
            "SquareMatrixQ[{{a, b}, {c, d}}]", "True", 0);
        assert_eval_eq(
            "SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
        assert_eval_eq(
            "SquareMatrixQ[{{1}, {2, 3}}]", "False", 0);
        assert_eval_eq("SquareMatrixQ[{1, 2, 3}]", "False", 0);
        assert_eval_eq("SquareMatrixQ[{}]", "False", 0);
        assert_eval_eq(
            "SquareMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_square_2x2_numeric);
    TEST(test_square_3x3_numeric);
    TEST(test_square_1x1);
    TEST(test_square_identity_matrix);
    TEST(test_square_larger);

    TEST(test_square_symbolic_2x2);
    TEST(test_square_symbolic_3x3);
    TEST(test_square_unevaluated_expr_entries);

    TEST(test_scalar_inputs);
    TEST(test_vector_inputs);
    TEST(test_empty_list);
    TEST(test_single_empty_row);

    TEST(test_single_row_rectangular);
    TEST(test_single_column_rectangular);
    TEST(test_rectangular_2x3);
    TEST(test_rectangular_3x2);
    TEST(test_rectangular_3x4);

    TEST(test_ragged_short_then_long);
    TEST(test_ragged_long_then_short);
    TEST(test_ragged_three_rows);

    TEST(test_three_d_tensor_rejected);
    TEST(test_mixed_atom_and_list_entry_rejected);

    TEST(test_matches_dimensions_eq_nn);
    TEST(test_implies_matrixq);
    TEST(test_matrixq_not_implies_squarematrixq);

    TEST(test_zero_arg_stays_unevaluated);
    TEST(test_three_args_stays_unevaluated);
    TEST(test_two_args_stays_unevaluated);
    TEST(test_rule_extra_arg_stays_unevaluated);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All SquareMatrixQ tests passed!\n");
    return 0;
}

/* Unit tests for ConjugateTranspose.
 *
 * Covers:
 *   - 2-D real matrices (no change beyond Transpose for purely real entries)
 *   - 2-D complex matrices (numeric + Complex atoms)
 *   - 1-D vectors (shape preserved, entries conjugated)
 *   - Symbolic matrices (Conjugate wraps unknown symbols)
 *   - 2-arg form ConjugateTranspose[m, spec]
 *   - 3-D tensors with non-trivial permutations
 *   - Attribute introspection (Protected)
 *   - ConjugateTranspose[m] == Conjugate[Transpose[m]] equivalence
 *   - Idempotence: ConjugateTranspose[ConjugateTranspose[m]] == m
 *   - Edge cases that should stay unevaluated (atoms, wrong arity)
 *
 * Each test parses and evaluates input strings and compares the printed
 * result against the expected output, which is the canonical Mathilda
 * surface form.  We deliberately exercise the *evaluator* path -- not the
 * builtin directly -- so the dispatch, attribute application, and
 * Conjugate-listable threading are all in scope.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include <string.h>
#include <stdlib.h>

static void test_real_matrix_2x2(void) {
    /* For a purely real matrix the conjugate is a no-op, so the result is
     * just the transpose. */
    assert_eval_eq("ConjugateTranspose[{{1, 2}, {3, 4}}]",
                   "{{1, 3}, {2, 4}}", 0);
    assert_eval_eq("ConjugateTranspose[{{1, 2, 3}, {4, 5, 6}}]",
                   "{{1, 4}, {2, 5}, {3, 6}}", 0);
}

static void test_complex_matrix_2x3(void) {
    /* From the Mathematica docs: m = {{1, 2*I, 3}, {3 + 4*I, 5, I}}.
     *
     *   Transpose[m]            = {{1, 3+4I}, {2*I, 5}, {3, I}}
     *   Conjugate[Transpose[m]] = {{1, 3-4I}, {-2*I, 5}, {3, -I}}
     */
    assert_eval_eq("ConjugateTranspose[{{1, 2*I, 3}, {3 + 4*I, 5, I}}]",
                   "{{1, 3 - 4*I}, {-2*I, 5}, {3, -I}}", 0);
}

static void test_complex_matrix_symbolic_entries(void) {
    /* ConjugateTranspose[{{a, b}, {1+2I, 3+4I}}] -- a, b are symbolic. */
    assert_eval_eq("ConjugateTranspose[{{a, b}, {1 + 2*I, 3 + 4*I}}]",
                   "{{Conjugate[a], 1 - 2*I}, {Conjugate[b], 3 - 4*I}}", 0);
}

static void test_row_matrix_becomes_column(void) {
    /* Conjugate-transposing a row matrix (1xn) gives an nx1 column matrix.
     *
     *   r = {{1.5, 2.2*I, 3.1 + 4.4*I}}
     *   ConjugateTranspose[r] = {{1.5}, {-2.2*I}, {3.1 - 4.4*I}}
     *
     * Inexact real entries propagate, so the imaginary zero on the first
     * entry stays printed as `1.5`. */
    assert_eval_eq(
        "ConjugateTranspose[{{1.5, 2.2*I, 3.1 + 4.4*I}}]",
        "{{1.5}, {0.0 - 2.2*I}, {3.1 - 4.4*I}}", 0);
}

static void test_vector_keeps_shape(void) {
    /* On a 1-D vector ConjugateTranspose is element-wise Conjugate.
     * The list shape must not change. */
    assert_eval_eq("ConjugateTranspose[{1, 2*I, 3 + 4*I}]",
                   "{1, -2*I, 3 - 4*I}", 0);
    assert_eval_eq("ConjugateTranspose[{a, b, c}]",
                   "{Conjugate[a], Conjugate[b], Conjugate[c]}", 0);
    assert_eval_eq("ConjugateTranspose[{1.5, 2.2*I, 3.1 + 4.4*I}]",
                   "{1.5, 0.0 - 2.2*I, 3.1 - 4.4*I}", 0);
}

static void test_diagonal_extraction(void) {
    /* Transpose[m, {1, 1}] extracts the diagonal; ConjugateTranspose with
     * the same spec extracts the conjugated diagonal. */
    assert_eval_eq("ConjugateTranspose[{{1, 2}, {3, 4 + 5*I}}, {1, 1}]",
                   "{1, 4 - 5*I}", 0);
}

static void test_two_arg_permutation(void) {
    /* For a 2-D matrix, the spec {2, 1} is the same as the default. */
    assert_eval_eq(
        "ConjugateTranspose[{{1 + I, 2}, {3, 4 - 5*I}}, {2, 1}]",
        "{{1 - I, 3}, {2, 4 + 5*I}}", 0);
}

static void test_three_d_tensor_permutation(void) {
    /* Build a 2x2x2 tensor with complex entries and check that the level
     * permutation {3, 2, 1} reorders correctly (and conjugates). */
    assert_eval_eq(
        "ConjugateTranspose[{{{1 + I, 2}, {3, 4 - I}}, {{5, 6 + 2*I}, {7, 8}}}, {3, 2, 1}]",
        "{{{1 - I, 5}, {3, 7}}, {{2, 6 - 2*I}, {4 + I, 8}}}", 0);
}

static void test_matches_conjugate_of_transpose(void) {
    /* For matrices, ConjugateTranspose[m] is equivalent to
     * Conjugate[Transpose[m]]. */
    const char* inputs[] = {
        "{{1, 2*I, 3}, {3 + 4*I, 5, I}}",
        "{{1 + I, 2 - 3*I}, {4 + 5*I, 6}}",
        "{{a, b}, {c + I d, e}}",
        "{{1, 2}, {3, 4}, {5, 6}}",
    };
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i) {
        char buf1[512];
        char buf2[512];
        snprintf(buf1, sizeof(buf1), "ConjugateTranspose[%s]", inputs[i]);
        snprintf(buf2, sizeof(buf2), "Conjugate[Transpose[%s]]", inputs[i]);

        Expr* p1 = parse_expression(buf1);
        Expr* e1 = evaluate(p1);
        Expr* p2 = parse_expression(buf2);
        Expr* e2 = evaluate(p2);

        char* s1 = expr_to_string(e1);
        char* s2 = expr_to_string(e2);
        if (strcmp(s1, s2) != 0) {
            fprintf(stderr,
                "FAIL: ConjugateTranspose != Conjugate[Transpose] for %s\n"
                "  ConjugateTranspose: %s\n"
                "  Conjugate[Tr]:      %s\n",
                inputs[i], s1, s2);
            ASSERT(0);
        }
        free(s1); free(s2);
        expr_free(p1); expr_free(e1);
        expr_free(p2); expr_free(e2);
    }
}

static void test_idempotent_on_complex_matrix(void) {
    /* ConjugateTranspose is an involution on matrices (it's its own
     * inverse): ConjugateTranspose[ConjugateTranspose[m]] == m. */
    /* Only test on numeric matrices: Mathilda does not (yet) simplify
     * Conjugate[Conjugate[x]] for symbolic x, so the involution check
     * would fail on symbolic entries. */
    const char* inputs[] = {
        "{{1, 2*I, 3}, {3 + 4*I, 5, I}}",
        "{{1 + I, 2 - 3*I}, {4 + 5*I, 6}}",
        "{{1, 2}, {3, 4}}",
    };
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i) {
        char buf_double[512];
        char buf_single[512];
        snprintf(buf_double, sizeof(buf_double),
                 "ConjugateTranspose[ConjugateTranspose[%s]]", inputs[i]);
        snprintf(buf_single, sizeof(buf_single), "%s", inputs[i]);

        Expr* pd = parse_expression(buf_double);
        Expr* ed = evaluate(pd);
        Expr* ps = parse_expression(buf_single);
        Expr* es = evaluate(ps);

        char* sd = expr_to_string(ed);
        char* ss = expr_to_string(es);
        if (strcmp(sd, ss) != 0) {
            fprintf(stderr,
                "FAIL: idempotence broken on %s\n"
                "  double: %s\n"
                "  single: %s\n",
                inputs[i], sd, ss);
            ASSERT(0);
        }
        free(sd); free(ss);
        expr_free(pd); expr_free(ed);
        expr_free(ps); expr_free(es);
    }
}

static void test_stays_unevaluated(void) {
    /* Non-array input: stays unevaluated. */
    assert_eval_eq("ConjugateTranspose[m]", "ConjugateTranspose[m]", 0);
    /* Wrong arity: stays unevaluated. */
    assert_eval_eq("ConjugateTranspose[]", "ConjugateTranspose[]", 0);
    assert_eval_eq(
        "ConjugateTranspose[{{1, 2}, {3, 4}}, {1, 2}, extra]",
        "ConjugateTranspose[{{1, 2}, {3, 4}}, {1, 2}, extra]", 0);
}

static void test_protected_attribute(void) {
    /* ConjugateTranspose must be Protected. */
    SymbolDef* def = symtab_get_def("ConjugateTranspose");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* It is NOT Listable: that would make ConjugateTranspose[{m1, m2}]
     * thread into {ConjugateTranspose[m1], ConjugateTranspose[m2]}, but
     * the Mathematica convention is that the outer argument is one
     * matrix, not a list of matrices. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("ConjugateTranspose");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "conjugate transpose") != NULL ||
           strstr(def->docstring, "Conjugate[Transpose") != NULL);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_real_matrix_2x2);
    TEST(test_complex_matrix_2x3);
    TEST(test_complex_matrix_symbolic_entries);
    TEST(test_row_matrix_becomes_column);
    TEST(test_vector_keeps_shape);
    TEST(test_diagonal_extraction);
    TEST(test_two_arg_permutation);
    TEST(test_three_d_tensor_permutation);
    TEST(test_matches_conjugate_of_transpose);
    TEST(test_idempotent_on_complex_matrix);
    TEST(test_stays_unevaluated);
    TEST(test_protected_attribute);
    TEST(test_docstring_set);

    printf("All ConjugateTranspose tests passed!\n");
    return 0;
}

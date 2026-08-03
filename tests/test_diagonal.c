/* Unit tests for Diagonal.
 *
 *   Diagonal[m]     -- the leading diagonal {m[[1,1]], m[[2,2]], ...}, length
 *                      Min[rows, cols] (so it works for a non-square m).
 *   Diagonal[m, k]  -- the k-th diagonal: k > 0 above the leading diagonal,
 *                      k < 0 below; an out-of-range k gives {}.
 *
 * Coverage:
 *   - Leading, super- (k>0) and sub-diagonals (k<0), symbolic and integer.
 *   - Non-square matrices (wide and tall).
 *   - Machine-real and complex entries; exact-integer entries stay exact
 *     (including past 2^53, which a double round-trip would lose).
 *   - Higher-rank tensors: Diagonal of a rank-3 array is rank-2.
 *   - Out-of-range k gives {}; a vector / non-matrix is left unevaluated.
 *   - The packed / NDArray fast path (ndstruct_diagonal): a visible NDArray in
 *     gives a visible NDArray out with the right values.
 *   - Compile[] lowering (ND_FNS) for Diagonal, and the COMPILE_MISSING.md §1
 *     reductions Tr / Det / MatrixRank / Norm (ND_REDS): each compiled result
 *     equals the interpreter's.
 *   - Attributes (Protected), docstring, and the interned SYM_Diagonal.
 *   - Repeated-evaluation loop to surface double-free / leak regressions under
 *     valgrind.
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

/* Evaluate `input` and assert it is returned unevaluated (printed form equal). */
static void assert_unevaluated(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "expected unevaluated %s, got: %s", expected, s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* Evaluate `input` and assert the printed result CONTAINS `substr`. Used for the
 * CompileDiagnostics association, whose full text is verbose and version-y. */
static void assert_eval_contains(const char* input, const char* substr) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, substr) != NULL,
               "expected result of %s to contain \"%s\", got: %s",
               input, substr, s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Leading diagonal ------------------------------------------------- */

static void test_leading_symbolic(void) {
    assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}]", "{a, e, i}", 0);
}

static void test_leading_integer(void) {
    assert_eval_eq("Diagonal[{{2, 3, 1}, {2, 2, 1}, {3, 1, 2}}]", "{2, 2, 2}", 0);
}

/* --- Super- and sub-diagonals ---------------------------------------- */

static void test_superdiagonal(void) {
    assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, 1]", "{b, f}", 0);
}

static void test_subdiagonal(void) {
    assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, -1]", "{d, h}", 0);
}

static void test_subdiagonal_k2(void) {
    assert_eval_eq(
        "Diagonal[{{a, b, c, d}, {e, f, g, h}, {i, j, k, l}, {m, n, o, p}}, -2]",
        "{i, n}", 0);
}

/* --- Non-square ------------------------------------------------------- */

static void test_nonsquare_wide(void) {
    assert_eval_eq("Diagonal[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]",
                   "{1, 6, 11}", 0);
}

static void test_nonsquare_tall(void) {
    assert_eval_eq("Diagonal[{{3, 2, 2}, {2, 3, -2}, {4, 2, 1}, {3, 7, 9}}]",
                   "{3, 3, 1}", 0);
}

/* --- Machine real / complex ------------------------------------------ */

static void test_machine_real(void) {
    assert_eval_eq(
        "Diagonal[{{1.1, 12.2, 3.23}, {2.3, 42.2, 35.3}, {1.2, 3.1, 2.3}}]",
        "{1.1, 42.2, 2.3}", 0);
}

static void test_complex_superdiagonal(void) {
    /* Symbolic entries (4 Pi, E) keep this on the generic List path. */
    assert_eval_eq(
        "Diagonal[{{1. + I, 2, 3 - 2 I}, {0, 4 Pi, 5 I}, {E, 0, 6}}, 1]",
        "{2, 5*I}", 0);
}

/* --- Exactness -------------------------------------------------------- */

static void test_integer_element_stays_exact(void) {
    assert_eval_eq("Head[Diagonal[{{2, 3}, {4, 5}}][[1]]]", "Integer", 0);
}

static void test_int64_exact_beyond_2_53(void) {
    /* memcpy, not a double round-trip: these differ only in the low bits. */
    assert_eval_eq(
        "Diagonal[{{9007199254740993, 0}, {0, 9007199254740995}}]",
        "{9007199254740993, 9007199254740995}", 0);
}

/* --- Higher rank ------------------------------------------------------ */

static void test_rank3_tensor(void) {
    /* Diagonal of a rank-3 array is rank-2: m[[1,1]] and m[[2,2]]. */
    assert_eval_eq("Diagonal[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
                   "{{1, 2}, {7, 8}}", 0);
}

/* --- Empty / malformed ----------------------------------------------- */

static void test_out_of_range_k(void) {
    assert_eval_eq("Diagonal[{{1, 2}, {3, 4}}, 5]", "{}", 0);
    assert_eval_eq("Diagonal[{{1, 2}, {3, 4}}, -5]", "{}", 0);
}

static void test_nonmatrix_unevaluated(void) {
    assert_unevaluated("Diagonal[{1, 2, 3}]", "Diagonal[{1, 2, 3}]");
    assert_unevaluated("Diagonal[x]", "Diagonal[x]");
    assert_unevaluated("Diagonal[{{1, 2}, {3, 4}}, x]",
                       "Diagonal[{{1, 2}, {3, 4}}, x]");
}

/* --- Packed / NDArray fast path -------------------------------------- */

static void test_ndarray_leading(void) {
    /* Visible NDArray in -> visible NDArray out, correct values. */
    assert_eval_eq("NDArrayQ[Diagonal[NDArray[{{1.5, 2.5}, {3.5, 4.5}}]]]", "True", 0);
    assert_eval_eq("Normal[Diagonal[NDArray[{{1.5, 2.5}, {3.5, 4.5}}]]] === {1.5, 4.5}",
                   "True", 0);
}

static void test_ndarray_super_and_sub(void) {
    assert_eval_eq(
        "Normal[Diagonal[NDArray[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}], 1]] "
        "=== {2., 6.}", "True", 0);
    assert_eval_eq(
        "Normal[Diagonal[NDArray[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}], -1]] "
        "=== {4., 8.}", "True", 0);
}

static void test_ndarray_complex(void) {
    assert_eval_eq(
        "Normal[Diagonal[NDArray[{{1. + 2. I, 3.}, {5., 7. + 8. I}}]]] "
        "=== {1. + 2. I, 7. + 8. I}", "True", 0);
}

static void test_ndarray_agrees_with_list(void) {
    /* The buffer path and the generic List path give the same values. */
    assert_eval_eq(
        "Normal[Diagonal[NDArray[{{1., 2.}, {3., 4.}}]]] "
        "=== Diagonal[{{1., 2.}, {3., 4.}}]", "True", 0);
}

/* --- Compile[] lowering ---------------------------------------------- */

static void test_compile_diagonal(void) {
    assert_eval_contains("CompileDiagnostics[{{m, _Real, 2}}, Diagonal[m]]",
                         "Compiled\" -> True");
    /* Compiled value equals the interpreter's, for k = 0, 1, -1. */
    assert_eval_eq(
        "Compile[{{m, _Real, 2}}, Diagonal[m]][{{1., 2.}, {3., 4.}}] === {1., 4.}",
        "True", 0);
    assert_eval_eq(
        "Compile[{{m, _Real, 2}, {k, _Integer}}, Diagonal[m, k]]"
        "[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}, 1] === {2., 6.}", "True", 0);
    assert_eval_eq(
        "Compile[{{m, _Real, 2}, {k, _Integer}}, Diagonal[m, k]]"
        "[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}, -1] === {4., 8.}", "True", 0);
    /* Fused with a reduction downstream (the cliff scenario). */
    assert_eval_eq(
        "Compile[{{m, _Real, 2}}, Total[Diagonal[m]]][{{1., 2.}, {3., 4.}}] === 5.",
        "True", 0);
}

static void test_compile_reductions(void) {
    assert_eval_contains("CompileDiagnostics[{{m, _Real, 2}}, Tr[m]]", "Compiled\" -> True");
    assert_eval_contains("CompileDiagnostics[{{m, _Real, 2}}, Det[m]]", "Compiled\" -> True");
    assert_eval_contains("CompileDiagnostics[{{m, _Real, 2}}, MatrixRank[m]]", "Compiled\" -> True");
    assert_eval_contains("CompileDiagnostics[{{m, _Real, 2}}, Norm[m]]", "Compiled\" -> True");
    assert_eval_contains("CompileDiagnostics[{{v, _Real, 1}}, Norm[v]]", "Compiled\" -> True");

    assert_eval_eq("Compile[{{m, _Real, 2}}, Tr[m]][{{1., 2.}, {3., 4.}}] === 5.", "True", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, Det[m]][{{1., 2.}, {3., 4.}}] === -2.", "True", 0);
    assert_eval_eq("Compile[{{m, _Real, 2}}, MatrixRank[m]][{{1., 2.}, {2., 4.}}] === 1", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Norm[v]][{3., 4.}] === 5.", "True", 0);
}

/* --- Attributes / docstring / symbol --------------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("Diagonal");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("Diagonal");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "diagonal") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_Diagonal != NULL);
    ASSERT(strcmp(SYM_Diagonal, "Diagonal") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    for (int t = 0; t < 50; t++) {
        assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}]", "{a, e, i}", 0);
        assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, 1]", "{b, f}", 0);
        assert_eval_eq("Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, -1]", "{d, h}", 0);
        assert_eval_eq("Normal[Diagonal[NDArray[{{1., 2.}, {3., 4.}}]]] === {1., 4.}",
                       "True", 0);
        assert_eval_eq("Diagonal[{{1, 2}, {3, 4}}, 5]", "{}", 0);
        assert_unevaluated("Diagonal[{1, 2, 3}]", "Diagonal[{1, 2, 3}]");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_leading_symbolic);
    TEST(test_leading_integer);

    TEST(test_superdiagonal);
    TEST(test_subdiagonal);
    TEST(test_subdiagonal_k2);

    TEST(test_nonsquare_wide);
    TEST(test_nonsquare_tall);

    TEST(test_machine_real);
    TEST(test_complex_superdiagonal);

    TEST(test_integer_element_stays_exact);
    TEST(test_int64_exact_beyond_2_53);

    TEST(test_rank3_tensor);

    TEST(test_out_of_range_k);
    TEST(test_nonmatrix_unevaluated);

    TEST(test_ndarray_leading);
    TEST(test_ndarray_super_and_sub);
    TEST(test_ndarray_complex);
    TEST(test_ndarray_agrees_with_list);

    TEST(test_compile_diagonal);
    TEST(test_compile_reductions);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All Diagonal tests passed!\n");
    symtab_clear();
    return 0;
}

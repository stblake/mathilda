#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "FlattenAt %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- The documented examples ---------- */

static void test_fa_basic(void) {
    /* Flatten the sublist at position 2. */
    run_full("FlattenAt[{a, {b, c}, {d, e}, {f}}, 2]",
             "List[a, b, c, List[d, e], List[f]]");
}

static void test_fa_multi_positions(void) {
    /* Two positions, resolved against the original expression (no index shift):
       position 2 -> b, c and position 4 -> f. */
    run_full("FlattenAt[{a, {b, c}, {d, e}, {f}}, {{2}, {4}}]",
             "List[a, b, c, List[d, e], f]");
}

static void test_fa_nonlist_head(void) {
    /* FlattenAt splices the arguments of any head, not just List. */
    run_full("FlattenAt[f[g[1, 2], g[3, 4], g[5]], {2}]",
             "f[g[1, 2], 3, 4, g[5]]");
}

static void test_fa_removes_head(void) {
    /* The head of the part at the position is removed and its elements spliced. */
    run_full("FlattenAt[{1, {{2}, {3}}, 4}, {2}]",
             "List[1, List[2], List[3], 4]");
}

/* ---------- Position-spec equivalences ---------- */

static void test_fa_single_path_equals_bare(void) {
    /* {2} (a single path) means the same as the bare index 2. */
    run_full("FlattenAt[{a, {b, c}}, {2}]", "List[a, b, c]");
    run_full("FlattenAt[{a, {b, c}}, 2]", "List[a, b, c]");
}

static void test_fa_negative_index(void) {
    /* A negative index counts from the end. */
    run_full("FlattenAt[{a, {b, c}, {d, e}}, -1]",
             "List[a, List[b, c], d, e]");
}

static void test_fa_deep_path(void) {
    /* A single deep position {2, 2} addresses a nested part. */
    run_full("FlattenAt[{a, {b, {c, d}}}, {2, 2}]",
             "List[a, List[b, c, d]]");
}

static void test_fa_multi_nonlist_head(void) {
    /* Several positions on a non-List head, with an untouched element between. */
    run_full("FlattenAt[f[g[1, 2], h[3, 4], g[5]], {{1}, {3}}]",
             "f[1, 2, h[3, 4], 5]");
}

/* ---------- Atoms have nothing to splice ---------- */

static void test_fa_atom_is_noop(void) {
    /* A symbol at the position is atomic -> the expression is unchanged. */
    run_full("FlattenAt[{a, b, c}, 1]", "List[a, b, c]");
}

static void test_fa_rational_not_spliced(void) {
    /* Rational is stored as a function but is atomic: it must not become
       Sequence[1, 2]. */
    run_full("FlattenAt[{1/2, x}, 1]", "List[Rational[1, 2], x]");
}

static void test_fa_complex_not_spliced(void) {
    run_full("FlattenAt[{2 + 3 I, x}, 1]", "List[Complex[2, 3], x]");
}

/* ---------- The result re-evaluates after splicing ---------- */

static void test_fa_reevaluates_after_splice(void) {
    /* Hold is not held by FlattenAt's argument evaluation; once its head is
       removed the spliced 1 + 1 evaluates to 2. */
    run_full("FlattenAt[{Hold[1 + 1], b}, 1]", "List[2, b]");
}

/* ---------- Visible NDArray materializes to a ragged list ---------- */

static void test_fa_ndarray_ragged(void) {
    /* A packed rank-2 array becomes ragged when one row is flattened, so it is
       materialized to a nested List and never repacked. */
    run_full("FlattenAt[NDArray[{{1, 2}, {3, 4}}, DataType -> \"int64\"], 1]",
             "List[1, 2, List[3, 4]]");
}

/* ---------- One-argument form is the identity ---------- */

static void test_fa_one_arg_identity(void) {
    /* No positions given -> nothing is flattened. */
    run_full("FlattenAt[{a, {b, c}}]", "List[a, List[b, c]]");
}

/* ---------- Out-of-range / arity -> unevaluated ---------- */

static void test_fa_out_of_range_unevaluated(void) {
    run_full("FlattenAt[{a, b, c}, 5]",
             "FlattenAt[List[a, b, c], 5]");
    run_full("FlattenAt[{a, b, c}, -5]",
             "FlattenAt[List[a, b, c], -5]");
}

static void test_fa_wrong_arity(void) {
    /* 0 and 3 args emit the argt message (to stderr) and stay unevaluated. */
    run_full("FlattenAt[]", "FlattenAt[]");
    run_full("FlattenAt[a, b, c]", "FlattenAt[a, b, c]");
}

/* ---------- Attributes & documentation ---------- */

static void test_fa_attributes_protected(void) {
    run_full("MemberQ[Attributes[FlattenAt], Protected]", "True");
}

static void test_fa_docstring_present(void) {
    SymbolDef* def = symtab_get_def("FlattenAt");
    ASSERT_MSG(def != NULL && def->docstring != NULL && def->docstring[0] != '\0',
               "FlattenAt should have a non-empty docstring");
}

int main(void) {
    symtab_init();
    core_init();

    /* Documented examples */
    TEST(test_fa_basic);
    TEST(test_fa_multi_positions);
    TEST(test_fa_nonlist_head);
    TEST(test_fa_removes_head);

    /* Position-spec equivalences */
    TEST(test_fa_single_path_equals_bare);
    TEST(test_fa_negative_index);
    TEST(test_fa_deep_path);
    TEST(test_fa_multi_nonlist_head);

    /* Atoms */
    TEST(test_fa_atom_is_noop);
    TEST(test_fa_rational_not_spliced);
    TEST(test_fa_complex_not_spliced);

    /* Re-evaluation */
    TEST(test_fa_reevaluates_after_splice);

    /* NDArray */
    TEST(test_fa_ndarray_ragged);

    /* One-argument identity */
    TEST(test_fa_one_arg_identity);

    /* Unevaluated / arity */
    TEST(test_fa_out_of_range_unevaluated);
    TEST(test_fa_wrong_arity);

    /* Attributes & docs */
    TEST(test_fa_attributes_protected);
    TEST(test_fa_docstring_present);

    printf("All FlattenAt tests passed!\n");
    return 0;
}

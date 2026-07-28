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
    ASSERT_MSG(strcmp(s, expected) == 0, "MapAt %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "MapAt %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* --- Simple integer-position cases --- */

static void test_mapat_basic_integer_position(void) {
    run_full("MapAt[f, {a, b, c, d}, 2]", "List[a, f[b], c, d]");
}

static void test_mapat_first_position(void) {
    run_full("MapAt[f, {a, b, c, d}, 1]", "List[f[a], b, c, d]");
}

static void test_mapat_last_position(void) {
    run_full("MapAt[f, {a, b, c, d}, 4]", "List[a, b, c, f[d]]");
}

static void test_mapat_out_of_range(void) {
    /* A position that does not exist leaves MapAt unevaluated, as Part does
     * for {a,b,c}[[99]]. */
    run_full("MapAt[f, {a, b, c}, 99]", "MapAt[f, List[a, b, c], 99]");
    run_full("MapAt[f, {a, b, c}, -99]", "MapAt[f, List[a, b, c], -99]");
}

/* --- Negative indices --- */

static void test_mapat_negative_index(void) {
    run_full("MapAt[h, {{a, b, c}, {d, e}, f, g}, -3]",
             "List[List[a, b, c], h[List[d, e]], f, g]");
}

static void test_mapat_negative_index_last(void) {
    run_full("MapAt[f, {a, b, c, d}, -1]", "List[a, b, c, f[d]]");
}

static void test_mapat_negative_index_first(void) {
    run_full("MapAt[f, {a, b, c, d}, -4]", "List[f[a], b, c, d]");
}

/* --- Position as a single list (deep path) --- */

static void test_mapat_nested_path_simple(void) {
    run_full("MapAt[f, {{a, b, c}, {d, e}}, {2, 1}]",
             "List[List[a, b, c], List[f[d], e]]");
}

static void test_mapat_nested_path_as_single_list(void) {
    run_full("MapAt[h, {{a, b, c}, {d, e}, f, g}, {2, 1}]",
             "List[List[a, b, c], List[h[d], e], f, g]");
}

static void test_mapat_nested_negative(void) {
    run_full("MapAt[f, {{a, b, c}, {d, e}}, {-1, -1}]",
             "List[List[a, b, c], List[d, f[e]]]");
}

/* --- Multiple positions --- */

static void test_mapat_multiple_top_level(void) {
    run_full("MapAt[f, {a, b, c, d}, {{1}, {4}}]",
             "List[f[a], b, c, f[d]]");
}

static void test_mapat_multiple_list_of_singletons(void) {
    run_full("MapAt[h, {{a, b, c}, {d, e}, f, g}, {{2}, {1}}]",
             "List[h[List[a, b, c]], h[List[d, e]], f, g]");
}

static void test_mapat_multiple_nested(void) {
    run_full("MapAt[h, {{a, b, c}, {d, e}, f, g}, {{1, 1}, {2, 2}, {3}}]",
             "List[List[h[a], b, c], List[d, h[e]], h[f], g]");
}

static void test_mapat_repeated_position(void) {
    /* The same position appearing twice must apply f twice. */
    run_full("MapAt[f, {a, b, c}, {{2}, {2}}]",
             "List[a, f[f[b]], c]");
}

/* --- All --- */

static void test_mapat_all_column(void) {
    run_full("MapAt[f, {{a, b, c}, {d, e}}, {All, 2}]",
             "List[List[a, f[b], c], List[d, f[e]]]");
}

static void test_mapat_all_top_level(void) {
    run_full("MapAt[f, {a, b, c}, {All}]",
             "List[f[a], f[b], f[c]]");
}

/* --- Span --- */

static void test_mapat_span(void) {
    run_full("MapAt[f, {1, 2, 3, 4, 5, 6}, 3;;4]",
             "List[1, 2, f[3], f[4], 5, 6]");
}

static void test_mapat_span_with_step(void) {
    run_full("MapAt[f, {1, 2, 3, 4, 5, 6}, 1;;6;;2]",
             "List[f[1], 2, f[3], 4, f[5], 6]");
}

static void test_mapat_span_negative(void) {
    run_full("MapAt[f, {1, 2, 3, 4, 5, 6}, -2;;-1]",
             "List[1, 2, 3, 4, f[5], f[6]]");
}

/* --- Arbitrary heads --- */

static void test_mapat_plus_head(void) {
    /* Plus is Orderless -- the evaluator canonicalises args after MapAt. */
    run_infix("MapAt[f, a + b + c + d, 2]", "a + c + d + f[b]");
}

static void test_mapat_power_paths(void) {
    run_infix("MapAt[f, x^2 + y^2, {{1, 1}, {2, 1}}]", "f[x]^2 + f[y]^2");
}

/* --- Head (position 0) --- */

static void test_mapat_head_zero(void) {
    run_full("MapAt[f, {a, b, c}, 0]", "f[List][a, b, c]");
}

static void test_mapat_head_zero_path(void) {
    run_full("MapAt[f, {a, b, c}, {0}]", "f[List][a, b, c]");
}

/* --- Identity / edge cases --- */

static void test_mapat_empty_position_list(void) {
    /* {} is an empty list of POSITIONS, so nothing is mapped. The position of
     * the whole expression is the empty path, spelled {{}}. */
    run_full("MapAt[f, h[a, b], {}]", "h[a, b]");
    run_full("MapAt[f, h[a, b], {{}}]", "f[h[a, b]]");
    run_full("MapAt[f, {a, b, c}, {}]", "List[a, b, c]");
}

static void test_mapat_atomic_integer_position(void) {
    /* Cannot descend into an atom: the position does not exist. */
    run_full("MapAt[f, 5, {1}]", "MapAt[f, 5, List[1]]");
    /* Rational and Complex are atoms too, so MapAt must not manufacture a
     * Rational[f[1], 2] out of 1/2. */
    run_full("MapAt[f, 1/2, 1]", "MapAt[f, Rational[1, 2], 1]");
}

/* --- Association positions --- */

static void test_mapat_assoc_does_not_mutate_input(void) {
    /* expr_copy is a refcount bump, so rebuilding an entry in place would
     * corrupt the caller's variable. */
    run_full("(max1 = <|\"x\" -> 1, \"y\" -> 2|>; MapAt[f, max1, \"x\"]; max1)",
             "Association[Rule[\"x\", 1], Rule[\"y\", 2]]");
    run_full("(max2 = <|\"k\" -> {1, 2}|>; MapAt[f, max2, {\"k\", 1}]; max2)",
             "Association[Rule[\"k\", List[1, 2]]]");
}

static void test_mapat_assoc_all(void) {
    run_full("MapAt[f, <|\"a\" -> 1, \"b\" -> 2|>, All]",
             "Association[Rule[\"a\", f[1]], Rule[\"b\", f[2]]]");
}

static void test_mapat_assoc_span(void) {
    run_full("MapAt[f, <|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, 1;;2]",
             "Association[Rule[\"a\", f[1]], Rule[\"b\", f[2]], Rule[\"c\", 3]]");
}

static void test_mapat_assoc_head(void) {
    run_full("MapAt[f, <|\"a\" -> 1|>, 0]", "f[Association][Rule[\"a\", 1]]");
}

static void test_mapat_assoc_absent_key(void) {
    run_full("MapAt[f, <|\"a\" -> 1|>, \"zz\"]",
             "MapAt[f, Association[Rule[\"a\", 1]], \"zz\"]");
}

static void test_mapat_assoc_repeated_key(void) {
    run_full("MapAt[f, <|\"a\" -> 1|>, {{\"a\"}, {\"a\"}}]",
             "Association[Rule[\"a\", f[f[1]]]]");
}

/* --- Operator form --- */

static void test_mapat_operator_form(void) {
    run_full("MapAt[f, 2][{a, b, c}]", "List[a, f[b], c]");
}

static void test_mapat_operator_form_pure_function(void) {
    /* The nested Function must not have its slot captured by the outer one. */
    run_full("MapAt[#^2 &, 2][{2, 3, 5}]", "List[2, 9, 5]");
}

static void test_mapat_operator_form_nested_path(void) {
    run_full("MapAt[f, {2, 1}][{{a, b}, {c, d}}]",
             "List[List[a, b], List[f[c], d]]");
}

/* --- Held arguments --- */

static void test_mapat_does_not_force_evaluation(void) {
    /* f[part] is left for the evaluator, so a surrounding Hold suppresses it. */
    run_full("MapAt[f, Hold[1 + 1], 1]", "Hold[f[Plus[1, 1]]]");
}

/* --- Malformed specs --- */

static void test_mapat_span_zero_step(void) {
    run_full("MapAt[f, {1, 2, 3}, 1;;3;;0]", "MapAt[f, List[1, 2, 3], Span[1, 3, 0]]");
}

static void test_mapat_span_upto(void) {
    /* UpTo clamps to the length, shared with Part's span normaliser. */
    run_full("MapAt[f, {a, b, c}, 1;;UpTo[9]]", "List[f[a], f[b], f[c]]");
}

static void test_mapat_key_on_non_association(void) {
    run_full("MapAt[f, {a, b}, Key[\"x\"]]", "MapAt[f, List[a, b], Key[\"x\"]]");
}

static void test_mapat_leaves_unmatched_untouched(void) {
    /* Verify the other branches are copied, not shared. */
    run_full("MapAt[f, {{a, b}, {c, d}}, {1, 2}]",
             "List[List[a, f[b]], List[c, d]]");
}

/* --- Composition with pure functions --- */

static void test_mapat_pure_function(void) {
    run_full("MapAt[#^2 &, {2, 3, 5}, 2]", "List[2, 9, 5]");
}

static void test_mapat_numeric_f(void) {
    run_infix("MapAt[Sin, {0, Pi/2, Pi}, 2]", "{0, 1, Pi}");
}

/* --- Repeated application verifying no aliasing --- */

static void test_mapat_applies_repeatedly_same_part(void) {
    /* Same path appearing 3 times should wrap f three times. */
    run_full("MapAt[g, {a, b}, {{1}, {1}, {1}}]", "List[g[g[g[a]]], b]");
}

/* --- Numeric / symbolic combos --- */

static void test_mapat_numeric_nested(void) {
    run_full("MapAt[f, {{1, 2}, {3, 4}}, {2, 2}]",
             "List[List[1, 2], List[3, f[4]]]");
}

static void test_mapat_all_at_deep(void) {
    run_full("MapAt[f, {{a, b}, {c, d}, {e, g}}, {All, 1}]",
             "List[List[f[a], b], List[f[c], d], List[f[e], g]]");
}

/* --- Documentation-informed additional tests --- */

static void test_mapat_second_top_level(void) {
    run_full("MapAt[h, {{a, b, c}, {d, e}, f, g}, 2]",
             "List[List[a, b, c], h[List[d, e]], f, g]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_mapat_basic_integer_position);
    TEST(test_mapat_first_position);
    TEST(test_mapat_last_position);
    TEST(test_mapat_out_of_range);
    TEST(test_mapat_negative_index);
    TEST(test_mapat_negative_index_last);
    TEST(test_mapat_negative_index_first);
    TEST(test_mapat_nested_path_simple);
    TEST(test_mapat_nested_path_as_single_list);
    TEST(test_mapat_nested_negative);
    TEST(test_mapat_multiple_top_level);
    TEST(test_mapat_multiple_list_of_singletons);
    TEST(test_mapat_multiple_nested);
    TEST(test_mapat_repeated_position);
    TEST(test_mapat_all_column);
    TEST(test_mapat_all_top_level);
    TEST(test_mapat_span);
    TEST(test_mapat_span_with_step);
    TEST(test_mapat_span_negative);
    TEST(test_mapat_plus_head);
    TEST(test_mapat_power_paths);
    TEST(test_mapat_head_zero);
    TEST(test_mapat_head_zero_path);
    TEST(test_mapat_empty_position_list);
    TEST(test_mapat_atomic_integer_position);
    TEST(test_mapat_assoc_does_not_mutate_input);
    TEST(test_mapat_assoc_all);
    TEST(test_mapat_assoc_span);
    TEST(test_mapat_assoc_head);
    TEST(test_mapat_assoc_absent_key);
    TEST(test_mapat_assoc_repeated_key);
    TEST(test_mapat_operator_form);
    TEST(test_mapat_operator_form_pure_function);
    TEST(test_mapat_operator_form_nested_path);
    TEST(test_mapat_does_not_force_evaluation);
    TEST(test_mapat_span_zero_step);
    TEST(test_mapat_span_upto);
    TEST(test_mapat_key_on_non_association);
    TEST(test_mapat_leaves_unmatched_untouched);
    TEST(test_mapat_pure_function);
    TEST(test_mapat_numeric_f);
    TEST(test_mapat_applies_repeatedly_same_part);
    TEST(test_mapat_numeric_nested);
    TEST(test_mapat_all_at_deep);
    TEST(test_mapat_second_top_level);

    printf("All MapAt tests passed!\n");
    return 0;
}

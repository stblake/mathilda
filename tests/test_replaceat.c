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
    ASSERT_MSG(strcmp(s, expected) == 0, "ReplaceAt %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "ReplaceAt %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ----- All examples from the spec, as documented ----- */

static void test_replace_at_basic_position(void) {
    /* ReplaceAt[{a,a,a,a},a->xx,2] -> {a,xx,a,a} */
    run_full("ReplaceAt[{a,a,a,a}, a -> xx, 2]", "List[a, xx, a, a]");
}

static void test_replace_at_multiple_positions(void) {
    /* ReplaceAt[{a,a,a,a},a->xx,{{1},{4}}] -> {xx,a,a,xx} */
    run_full("ReplaceAt[{a,a,a,a}, a -> xx, {{1},{4}}]", "List[xx, a, a, xx]");
}

static void test_replace_at_nested_part(void) {
    /* ReplaceAt[{{a,a},{a,a}},a->xx,{2,1}] -> {{a,a},{xx,a}} */
    run_full("ReplaceAt[{{a,a},{a,a}}, a -> xx, {2,1}]",
             "List[List[a, a], List[xx, a]]");
}

static void test_replace_at_negative_index(void) {
    /* ReplaceAt[{a,a,a,a},a->xx,-2] -> {a,a,xx,a} */
    run_full("ReplaceAt[{a,a,a,a}, a -> xx, -2]", "List[a, a, xx, a]");
}

static void test_replace_at_negative_nested(void) {
    /* ReplaceAt[{{a,a,a},{a,a,a}},a->xx,{-1,-2}] -> {{a,a,a},{a,xx,a}} */
    run_full("ReplaceAt[{{a,a,a},{a,a,a}}, a -> xx, {-1,-2}]",
             "List[List[a, a, a], List[a, xx, a]]");
}

static void test_replace_at_pattern_delayed(void) {
    /* ReplaceAt[{1,2,3,4}, x_:>2x-1, {{2},{4}}] -> {1,3,3,7} */
    run_full("ReplaceAt[{1,2,3,4}, x_ :> 2 x - 1, {{2},{4}}]",
             "List[1, 3, 3, 7]");
}

static void test_replace_at_multi_rules(void) {
    /* ReplaceAt[{a,b,c,d}, {a->xx, _->yy}, {{1},{2},{4}}] -> {xx,yy,c,yy} */
    run_full("ReplaceAt[{a,b,c,d}, {a -> xx, _ -> yy}, {{1},{2},{4}}]",
             "List[xx, yy, c, yy]");
}

static void test_replace_at_all_in_path(void) {
    /* ReplaceAt[{{a,a},{a,a}}, a->xx, {All,2}] -> {{a,xx},{a,xx}} */
    run_full("ReplaceAt[{{a,a},{a,a}}, a -> xx, {All,2}]",
             "List[List[a, xx], List[a, xx]]");
}

static void test_replace_at_single_2(void) {
    /* ReplaceAt[{{a,b},{c,d},e}, x_:>f[x], 2] -> {{a,b},f[{c,d}],e} */
    run_full("ReplaceAt[{{a,b},{c,d},e}, x_ :> f[x], 2]",
             "List[List[a, b], f[List[c, d]], e]");
}

static void test_replace_at_single_neg1(void) {
    /* ReplaceAt[{{a,b},{c,d},e}, x_:>f[x], -1] -> {{a,b},{c,d},f[e]} */
    run_full("ReplaceAt[{{a,b},{c,d},e}, x_ :> f[x], -1]",
             "List[List[a, b], List[c, d], f[e]]");
}

static void test_replace_at_path_2_1(void) {
    /* ReplaceAt[{{a,b},{c,d},e}, x_:>f[x], {2,1}] -> {{a,b},{f[c],d},e} */
    run_full("ReplaceAt[{{a,b},{c,d},e}, x_ :> f[x], {2,1}]",
             "List[List[a, b], List[f[c], d], e]");
}

static void test_replace_at_several_at_top(void) {
    /* ReplaceAt[{{a,b},{c,d},e}, x_:>f[x], {{1},{3}}] -> {f[{a,b}],{c,d},f[e]} */
    run_full("ReplaceAt[{{a,b},{c,d},e}, x_ :> f[x], {{1},{3}}]",
             "List[f[List[a, b]], List[c, d], f[e]]");
}

static void test_replace_at_several_nested(void) {
    /* ReplaceAt[{{a,b},{c,d},e}, x_:>f[x], {{1,2},{2,2},{3}}] -> {{a,f[b]},{c,f[d]},f[e]} */
    run_full("ReplaceAt[{{a,b},{c,d},e}, x_ :> f[x], {{1,2},{2,2},{3}}]",
             "List[List[a, f[b]], List[c, f[d]], f[e]]");
}

static void test_replace_at_span(void) {
    /* ReplaceAt[{a,a,a,a,a}, a->xx, 2;;4] -> {a,xx,xx,xx,a} */
    run_full("ReplaceAt[{a,a,a,a,a}, a -> xx, 2;;4]",
             "List[a, xx, xx, xx, a]");
}

static void test_replace_at_plus_head(void) {
    /* ReplaceAt[a+b+c+d, _->x, 2] -> a+c+d+x (Plus is Orderless+Flat, re-canonicalises) */
    run_infix("ReplaceAt[a+b+c+d, _ -> x, 2]", "a + c + d + x");
}

static void test_replace_at_nested_in_plus(void) {
    /* ReplaceAt[x^2+y^2, _->z, {{1,1},{2,1}}] -> 2 z^2 */
    run_infix("ReplaceAt[x^2 + y^2, _ -> z, {{1,1},{2,1}}]", "2 z^2");
}

static void test_replace_at_head_index_zero(void) {
    /* ReplaceAt[{a,b,c}, _->f, 0] -> f[a,b,c] */
    run_full("ReplaceAt[{a,b,c}, _ -> f, 0]", "f[a, b, c]");
}

/* ----- Additional edge cases ----- */

static void test_replace_at_no_match_unchanged(void) {
    /* No rule matches => part is left as-is */
    run_full("ReplaceAt[{a,b,c}, q -> z, 2]", "List[a, b, c]");
}

static void test_replace_at_out_of_range(void) {
    /* A position that does not exist leaves ReplaceAt unevaluated, as Part
     * does for {a,b,c}[[99]]. */
    run_full("ReplaceAt[{a,b,c}, _ -> z, 99]", "ReplaceAt[List[a, b, c], Rule[Blank[], z], 99]");
    run_full("ReplaceAt[{a,b,c}, _ -> z, -99]", "ReplaceAt[List[a, b, c], Rule[Blank[], z], -99]");
}

static void test_replace_at_repeated_position(void) {
    /* Repeated positions cause rules to be applied repeatedly to that part */
    run_full("ReplaceAt[{a,b,c}, x_ :> f[x], {{2},{2}}]",
             "List[a, f[f[b]], c]");
}

static void test_replace_at_single_rule_path_neg_in_nested(void) {
    /* Negative index in deep path */
    run_full("ReplaceAt[{{a,b},{c,d}}, x_ :> g[x], {-1,-1}]",
             "List[List[a, b], List[c, g[d]]]");
}

static void test_replace_at_span_with_step(void) {
    /* Span with step */
    run_full("ReplaceAt[{a,b,c,d,e}, _ -> z, 1;;5;;2]",
             "List[z, b, z, d, z]");
}

static void test_replace_at_span_negative_endpoints(void) {
    run_full("ReplaceAt[{a,b,c,d,e}, _ -> z, -3;;-1]",
             "List[a, b, z, z, z]");
}

static void test_replace_at_span_all_form(void) {
    run_full("ReplaceAt[{a,b,c}, _ -> z, 1;;All]",
             "List[z, z, z]");
}

static void test_replace_at_empty_position_list(void) {
    /* {} is an empty list of POSITIONS, so nothing is replaced. The position
     * of the whole expression is the empty path, spelled {{}}. */
    run_full("ReplaceAt[a, a -> b, {}]", "a");
    run_full("ReplaceAt[a, a -> b, {{}}]", "b");
}

static void test_replace_at_atomic_no_descend(void) {
    /* Cannot descend into an atom: the position does not exist. */
    run_full("ReplaceAt[a, _ -> z, 1]", "ReplaceAt[a, Rule[Blank[], z], 1]");
}

/* ----- Association positions (shared walker gave ReplaceAt these) ----- */

static void test_replace_at_assoc_key(void) {
    run_full("ReplaceAt[<|\"a\" -> 1, \"b\" -> 2|>, 1 -> 9, Key[\"a\"]]",
             "Association[Rule[\"a\", 9], Rule[\"b\", 2]]");
}

static void test_replace_at_assoc_string_key(void) {
    run_full("ReplaceAt[<|\"a\" -> 1, \"b\" -> 2|>, 2 -> 9, \"b\"]",
             "Association[Rule[\"a\", 1], Rule[\"b\", 9]]");
}

static void test_replace_at_assoc_positional(void) {
    run_full("ReplaceAt[<|\"a\" -> 1, \"b\" -> 2|>, x_Integer :> 10 x, 2]",
             "Association[Rule[\"a\", 1], Rule[\"b\", 20]]");
}

static void test_replace_at_assoc_all(void) {
    run_full("ReplaceAt[<|\"a\" -> 1, \"b\" -> 2|>, x_Integer :> 10 x, All]",
             "Association[Rule[\"a\", 10], Rule[\"b\", 20]]");
}

static void test_replace_at_assoc_does_not_mutate_input(void) {
    /* expr_copy is a refcount bump, so rebuilding an association entry in
     * place would corrupt the caller's variable. */
    run_full("(rax = <|\"p\" -> 1|>; ReplaceAt[rax, 1 -> 9, Key[\"p\"]]; rax)",
             "Association[Rule[\"p\", 1]]");
}

static void test_replace_at_first_rule_wins(void) {
    /* Rule order matters: the first matching rule wins. */
    run_full("ReplaceAt[{a,b,c}, {_ -> y, _ -> z}, 1]",
             "List[y, b, c]");
}

static void test_replace_at_pure_function(void) {
    /* Replacement may be a pure function expression; bindings substituted. */
    run_full("ReplaceAt[{1,2,3}, x_ :> x^2, 2]",
             "List[1, 4, 3]");
}

static void test_replace_at_immediate_rule(void) {
    /* Rule (->), not RuleDelayed. RHS is evaluated only once, at rule
     * creation; the matched binding is substituted by replace_bindings. */
    run_full("ReplaceAt[{a,b,c}, x_ -> g[x], 2]",
             "List[a, g[b], c]");
}

static void test_replace_at_multipath_plus(void) {
    /* Multiple positions inside a Plus */
    run_infix("ReplaceAt[a+b+c+d, _ -> z, {{1},{3}}]",
             "b + d + 2 z");
}

static void test_replace_at_head_replacement_in_function(void) {
    /* Replace head of an arbitrary head, not just List */
    run_full("ReplaceAt[h[a,b,c], _ -> g, 0]", "g[a, b, c]");
}

static void test_replace_at_head_in_nested(void) {
    /* Head replacement inside a nested expression via path */
    run_full("ReplaceAt[{h[a,b], c}, _ -> k, {1, 0}]",
             "List[k[a, b], c]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_replace_at_basic_position);
    TEST(test_replace_at_multiple_positions);
    TEST(test_replace_at_nested_part);
    TEST(test_replace_at_negative_index);
    TEST(test_replace_at_negative_nested);
    TEST(test_replace_at_pattern_delayed);
    TEST(test_replace_at_multi_rules);
    TEST(test_replace_at_all_in_path);
    TEST(test_replace_at_single_2);
    TEST(test_replace_at_single_neg1);
    TEST(test_replace_at_path_2_1);
    TEST(test_replace_at_several_at_top);
    TEST(test_replace_at_several_nested);
    TEST(test_replace_at_span);
    TEST(test_replace_at_plus_head);
    TEST(test_replace_at_nested_in_plus);
    TEST(test_replace_at_head_index_zero);

    TEST(test_replace_at_no_match_unchanged);
    TEST(test_replace_at_out_of_range);
    TEST(test_replace_at_repeated_position);
    TEST(test_replace_at_single_rule_path_neg_in_nested);
    TEST(test_replace_at_span_with_step);
    TEST(test_replace_at_span_negative_endpoints);
    TEST(test_replace_at_span_all_form);
    TEST(test_replace_at_empty_position_list);
    TEST(test_replace_at_atomic_no_descend);
    TEST(test_replace_at_assoc_key);
    TEST(test_replace_at_assoc_string_key);
    TEST(test_replace_at_assoc_positional);
    TEST(test_replace_at_assoc_all);
    TEST(test_replace_at_assoc_does_not_mutate_input);
    TEST(test_replace_at_first_rule_wins);
    TEST(test_replace_at_pure_function);
    TEST(test_replace_at_immediate_rule);
    TEST(test_replace_at_multipath_plus);
    TEST(test_replace_at_head_replacement_in_function);
    TEST(test_replace_at_head_in_nested);

    printf("All ReplaceAt tests passed!\n");
    return 0;
}

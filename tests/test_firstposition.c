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

/* Compare the FullForm of the evaluated input against expected. FullForm keeps
 * string quotes (Missing["NotFound"], Key["b"], "NoStrings") and renders
 * positions as List[...], so every case is unambiguous. */
static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "FirstPosition %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- The documented examples ---------- */

static void test_fp_basic_first(void) {
    /* First position at which b occurs. */
    run_full("FirstPosition[{a,b,a,a,b,c,b},b]", "List[2]");
}

static void test_fp_nested(void) {
    run_full("FirstPosition[{{a,a,b},{b,a,a},{a,b,a}},b]", "List[1, 3]");
}

static void test_fp_not_found(void) {
    run_full("FirstPosition[{x,y,z},b]", "Missing[\"NotFound\"]");
}

static void test_fp_power_pattern(void) {
    /* First position at which x to any power appears. */
    run_full("FirstPosition[{1+x^2,5,x^4,a+(1+x^2)^2},x^_]", "List[1, 2]");
}

/* ---------- Associations return a Key ---------- */

static void test_fp_assoc_nested(void) {
    run_full("FirstPosition[<|1->1+x^2,2-><|\"a\"->x^2|>,3->x^4,4->a+(1+x^2)^2|>,x^_]",
             "List[Key[1], 2]");
}

static void test_fp_assoc_prime_value(void) {
    /* First position with a prime value; the key "b" holds the first prime (2). */
    run_full("FirstPosition[<|\"a\"->1,\"b\"->2,\"c\"->3,\"d\"->4|>,_Integer?PrimeQ]",
             "List[Key[\"b\"]]");
}

/* ---------- Custom default value ---------- */

static void test_fp_custom_default(void) {
    run_full("FirstPosition[{1,2,3}, _?StringQ, \"NoStrings\"]", "\"NoStrings\"");
}

static void test_fp_custom_default_symbols(void) {
    /* Default returned on a miss over a plain list (no association). */
    run_full("FirstPosition[{a,b,a,a,b,c,b}, z, \"none\"]", "\"none\"");
}

/* ---------- Look for a value anywhere vs. at a fixed level ---------- */

static void test_fp_anywhere(void) {
    run_full("FirstPosition[{{1,2},{2,3},{3,1}},3]", "List[2, 2]");
}

static void test_fp_outer_level_only(void) {
    /* Level {1}: 3 only appears at level 2, so nothing is found -> the default. */
    run_full("FirstPosition[{{1,2},{2,3},{3,1}},3,Missing[\"NotFound\"],{1}]",
             "Missing[\"NotFound\"]");
}

static void test_fp_integer_levelspec_match(void) {
    /* Integer levelspec 2 == levels 1..2; the 3 at {2,2} is found. */
    run_full("FirstPosition[{{1,2},{2,3},{3,1}},3,Missing[\"NotFound\"],2]",
             "List[2, 2]");
}

/* ---------- Heads default True (searches heads), Heads -> False excludes ---------- */

static void test_fp_head_default(void) {
    run_full("FirstPosition[x^2+y^2,Power]", "List[1, 0]");
}

static void test_fp_head_false(void) {
    run_full("FirstPosition[x^2+y^2,Power,Heads->False]", "Missing[\"NotFound\"]");
}

static void test_fp_head_option_with_default_and_levelspec(void) {
    /* Full 4-positional form plus an explicit Heads -> True option. */
    run_full("FirstPosition[x^2+y^2, Power, Missing[\"NotFound\"], {0,Infinity}, Heads->True]",
             "List[1, 0]");
}

/* ---------- Pattern matching is not numerical equality ---------- */

static void test_fp_numeric_pattern_no_match(void) {
    run_full("FirstPosition[Range[-1,1,0.05],0.1]", "Missing[\"NotFound\"]");
}

static void test_fp_numeric_equality_condition(void) {
    run_full("FirstPosition[Range[-1,1,0.05],n_ /; n==0.1]", "List[23]");
}

/* ---------- Level 0 / whole-expression / arity ---------- */

static void test_fp_whole_expression(void) {
    /* The default {0, Infinity} includes level 0; a match of the whole
       expression is position {}. */
    run_full("FirstPosition[a+b, a+b]", "List[]");
}

static void test_fp_no_default_no_match(void) {
    run_full("FirstPosition[{1,2,3}, _?StringQ]", "Missing[\"NotFound\"]");
}

static void test_fp_one_arg_unevaluated(void) {
    /* Fewer than two arguments: FirstPosition stays unevaluated. */
    run_full("FirstPosition[{1,2,3}]", "FirstPosition[List[1, 2, 3]]");
}

/* ---------- default is evaluated only when returned (HoldRest) ---------- */

static void test_fp_default_lazy_on_match(void) {
    /* On a match the (held) default is never evaluated, so its side effect
       does not fire. */
    run_full("fpSideA = 0", "0");
    run_full("FirstPosition[{1, 2, 3}, 2, fpSideA = 99]", "List[2]");
    run_full("fpSideA", "0");
}

static void test_fp_default_evaluated_on_miss(void) {
    /* No match: the default is evaluated (returning 99 and firing its side
       effect). */
    run_full("fpSideB = 0", "0");
    run_full("FirstPosition[{1, 2, 3}, _?StringQ, fpSideB = 99]", "99");
    run_full("fpSideB", "99");
}

/* ---------- Attributes & documentation ---------- */

static void test_fp_attributes(void) {
    run_full("MemberQ[Attributes[FirstPosition], HoldRest]", "True");
    run_full("MemberQ[Attributes[FirstPosition], Protected]", "True");
}

static void test_fp_docstring_present(void) {
    SymbolDef* def = symtab_get_def("FirstPosition");
    ASSERT_MSG(def != NULL && def->docstring != NULL && def->docstring[0] != '\0',
               "FirstPosition should have a non-empty docstring");
}

int main(void) {
    symtab_init();
    core_init();

    /* Documented examples */
    TEST(test_fp_basic_first);
    TEST(test_fp_nested);
    TEST(test_fp_not_found);
    TEST(test_fp_power_pattern);

    /* Associations */
    TEST(test_fp_assoc_nested);
    TEST(test_fp_assoc_prime_value);

    /* Custom defaults */
    TEST(test_fp_custom_default);
    TEST(test_fp_custom_default_symbols);

    /* Levels */
    TEST(test_fp_anywhere);
    TEST(test_fp_outer_level_only);
    TEST(test_fp_integer_levelspec_match);

    /* Heads option */
    TEST(test_fp_head_default);
    TEST(test_fp_head_false);
    TEST(test_fp_head_option_with_default_and_levelspec);

    /* Pattern vs numeric equality */
    TEST(test_fp_numeric_pattern_no_match);
    TEST(test_fp_numeric_equality_condition);

    /* Level 0 / arity */
    TEST(test_fp_whole_expression);
    TEST(test_fp_no_default_no_match);
    TEST(test_fp_one_arg_unevaluated);

    /* Lazy default (HoldRest) */
    TEST(test_fp_default_lazy_on_match);
    TEST(test_fp_default_evaluated_on_miss);

    /* Attributes & docs */
    TEST(test_fp_attributes);
    TEST(test_fp_docstring_present);

    printf("All FirstPosition tests passed!\n");
    return 0;
}

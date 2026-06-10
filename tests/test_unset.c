/* Tests for Unset (`lhs =.`).
 *
 * Two halves:
 *   1. Parsing -- `=.` is a low-precedence postfix operator producing
 *      Unset[lhs]; it must not collide with a real literal on the RHS of
 *      Set (`k =.5` is Set[k, 0.5], not Unset).
 *   2. Semantics -- Unset removes the single rule whose left-hand side is
 *      `lhs` up to renaming of bound pattern variables: a bare symbol's
 *      OwnValue, or one specific DownValue on a head symbol, leaving sibling
 *      rules intact. Protected symbols are never touched; the result is
 *      always Null. Plus HoldFirst behaviour and memory hygiene.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "test_utils.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate a statement for its side effects and discard the result. */
static void run(const char* input) {
    struct Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    expr_free(evaluated);
}

/* Parse-only check: `input` must parse to `expected` in FullForm, with no
 * evaluation (isolates the parser from the Unset builtin). */
static void assert_parse_fullform(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    char* str = expr_to_string_fullform(parsed);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL parse: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
    }
    ASSERT(strcmp(str, expected) == 0);
    free(str);
    expr_free(parsed);
}

/* ------------------------------------------------------------------------
 *  Parsing
 * ---------------------------------------------------------------------- */

static void test_parse_symbol(void) {
    assert_parse_fullform("x =.", "Unset[x]");
}

static void test_parse_function(void) {
    assert_parse_fullform("f[x_] =.",
                          "Unset[f[Pattern[x, Blank[]]]]");
}

static void test_parse_low_precedence_captures_product(void) {
    /* `=.` binds looser than Times, so it wraps the whole product. */
    assert_parse_fullform("a b =.", "Unset[Times[a, b]]");
}

static void test_parse_does_not_swallow_real_literal(void) {
    /* `=.5` is Set with a 0.5 RHS, NOT Unset. */
    assert_parse_fullform("k =.5", "Set[k, 0.5]");
    assert_parse_fullform("k =.25", "Set[k, 0.25]");
}

static void test_parse_inside_hold(void) {
    assert_parse_fullform("Hold[x =.]", "Hold[Unset[x]]");
}

static void test_parse_condition_lhs(void) {
    assert_parse_fullform("g[x_] /; x > 0 =.",
        "Unset[Condition[g[Pattern[x, Blank[]]], Greater[x, 0]]]");
}

/* ------------------------------------------------------------------------
 *  OwnValue removal
 * ---------------------------------------------------------------------- */

static void test_unset_clears_value(void) {
    run("uv = 5");
    assert_eval_eq("uv", "5", 0);
    assert_eval_eq("uv =.", "Null", 0);
    assert_eval_eq("uv", "uv", 0);
    assert_eval_eq("OwnValues[uv]", "{}", 0);
}

static void test_unset_builtin_form(void) {
    /* Unset[x] is identical to `x =.`. */
    run("ub = 7");
    assert_eval_eq("Unset[ub]", "Null", 0);
    assert_eval_eq("ub", "ub", 0);
}

static void test_unset_holdfirst_targets_symbol(void) {
    /* With HoldFirst, Unset[ux] sees the symbol ux, not its value 5: the
     * value is cleared rather than Unset[5] being attempted. */
    run("uh = 5");
    run("uh =.");
    assert_eval_eq("uh", "uh", 0);
}

/* ------------------------------------------------------------------------
 *  DownValue removal
 * ---------------------------------------------------------------------- */

static void test_unset_clears_function(void) {
    run("uf[x_] := x^2");
    assert_eval_eq("uf[3]", "9", 0);
    assert_eval_eq("uf[x_] =.", "Null", 0);
    assert_eval_eq("uf[3]", "uf[3]", 0);
    assert_eval_eq("DownValues[uf]", "{}", 0);
}

static void test_unset_alpha_renaming(void) {
    /* Defined with x_, unset with y_: rules match up to renaming of the
     * bound pattern variable. */
    run("ur[x_] := x^2");
    assert_eval_eq("ur[4]", "16", 0);
    run("ur[y_] =.");
    assert_eval_eq("ur[4]", "ur[4]", 0);
}

static void test_unset_selective_downvalue(void) {
    /* `us[1] =.` removes only the us[1] rule; us[2] survives. */
    run("us[1] = 10");
    run("us[2] = 20");
    run("us[1] =.");
    assert_eval_eq("us[1]", "us[1]", 0);
    assert_eval_eq("us[2]", "20", 0);
}

static void test_unset_specific_keeps_general(void) {
    /* A specific rule can be removed without touching the general one. */
    run("ug[0] = 1");
    run("ug[n_] := n + 1");
    assert_eval_eq("ug[0]", "1", 0);
    run("ug[0] =.");
    /* ug[0] now falls through to the general rule. */
    assert_eval_eq("ug[0]", "1", 0);   /* 0 + 1 */
    assert_eval_eq("ug[5]", "6", 0);
}

static void test_unset_condition_lhs(void) {
    run("uc[x_] /; x > 0 := x + 100");
    assert_eval_eq("uc[5]", "105", 0);
    run("uc[x_] /; x > 0 =.");
    assert_eval_eq("uc[5]", "uc[5]", 0);
}

/* ------------------------------------------------------------------------
 *  No-ops / guards / return value
 * ---------------------------------------------------------------------- */

static void test_unset_nonexistent_is_harmless(void) {
    assert_eval_eq("uNever =.", "Null", 0);
    assert_eval_eq("uNever", "uNever", 0);
}

static void test_unset_returns_null(void) {
    run("urn = 1");
    assert_eval_eq("urn =.", "Null", 0);
}

static void test_unset_protected_is_blocked(void) {
    /* Sin is Protected: `Sin[x] =.` prints Unset::wrsym to stderr, returns
     * Null, and leaves Sin intact. */
    assert_eval_eq("Sin[x] =.", "Null", 0);
    assert_eval_eq("Sin[0]", "0", 0);
}

static void test_unset_does_not_break_builtin(void) {
    /* Plus is Protected; `Plus[a, b] =.` must not remove anything. */
    run("Plus[a, b] =.");
    assert_eval_eq("1 + 2", "3", 0);
}

/* ------------------------------------------------------------------------
 *  Attributes
 * ---------------------------------------------------------------------- */

static void test_unset_attributes(void) {
    assert_eval_eq("Attributes[Unset]", "{HoldFirst, Protected}", 0);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_unset_define_unset_loop(void) {
    for (int i = 0; i < 200; i++) {
        run("loopU[x_] := x^2 + x + 1");
        run("loopV = {1, 2, 3, loopU}");
        run("loopU[x_] =.");
        run("loopV =.");
    }
    assert_eval_eq("loopU[1]", "loopU[1]", 0);
    assert_eval_eq("loopV", "loopV", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Parsing */
    TEST(test_parse_symbol);
    TEST(test_parse_function);
    TEST(test_parse_low_precedence_captures_product);
    TEST(test_parse_does_not_swallow_real_literal);
    TEST(test_parse_inside_hold);
    TEST(test_parse_condition_lhs);

    /* OwnValue removal */
    TEST(test_unset_clears_value);
    TEST(test_unset_builtin_form);
    TEST(test_unset_holdfirst_targets_symbol);

    /* DownValue removal */
    TEST(test_unset_clears_function);
    TEST(test_unset_alpha_renaming);
    TEST(test_unset_selective_downvalue);
    TEST(test_unset_specific_keeps_general);
    TEST(test_unset_condition_lhs);

    /* No-ops / guards / return value */
    TEST(test_unset_nonexistent_is_harmless);
    TEST(test_unset_returns_null);
    TEST(test_unset_protected_is_blocked);
    TEST(test_unset_does_not_break_builtin);

    /* Attributes */
    TEST(test_unset_attributes);

    /* Memory hygiene */
    TEST(test_unset_define_unset_loop);

    printf("All unset_tests passed.\n");
    return 0;
}

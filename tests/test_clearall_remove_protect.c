/* Tests for ClearAll, Remove, Protect, Unprotect.
 *
 * Cover, for each builtin:
 *   - core semantics (values/attributes/definition cleared, removed, or
 *     protection toggled);
 *   - the List-of-specs argument form (e.g. ClearAll[{a, b}]);
 *   - the Locked / Protected guards (ClearAll and Remove never touch a
 *     Protected or Locked symbol; Protect/Unprotect never touch Locked);
 *   - HoldAll (the symbol, not its value, is operated on);
 *   - the WL-faithful return values (Protect/Unprotect return the list of
 *     changed names; ClearAll/Remove return Null);
 *   - that the head symbols themselves carry the right attributes;
 *   - that protection actually blocks redefinition end-to-end;
 *   - basic memory hygiene under a tight define/clear loop.
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

/* ------------------------------------------------------------------------
 *  ClearAll
 * ---------------------------------------------------------------------- */

static void test_clearall_clears_downvalues(void) {
    run("caF[x_]:=x^2");
    assert_eval_eq("caF[3]", "9", 0);
    assert_eval_eq("ClearAll[caF]", "Null", 0);
    /* DownValue gone -> symbol is inert again. */
    assert_eval_eq("caF[3]", "caF[3]", 0);
    assert_eval_eq("DownValues[caF]", "{}", 0);
}

static void test_clearall_clears_ownvalue(void) {
    run("caV = 42");
    assert_eval_eq("caV", "42", 0);
    run("ClearAll[caV]");
    assert_eval_eq("caV", "caV", 0);
}

static void test_clearall_clears_attributes(void) {
    run("caA[x_]:=x");
    run("SetAttributes[caA, Listable]");
    assert_eval_eq("Attributes[caA]", "{Listable}", 0);
    run("ClearAll[caA]");
    /* ClearAll, unlike Clear, also wipes attributes. */
    assert_eval_eq("Attributes[caA]", "{}", 0);
}

static void test_clearall_multiple_args(void) {
    run("caM1 = 1");
    run("caM2 = 2");
    run("ClearAll[caM1, caM2]");
    assert_eval_eq("caM1", "caM1", 0);
    assert_eval_eq("caM2", "caM2", 0);
}

static void test_clearall_list_spec(void) {
    run("caL1 = 1");
    run("caL2 = 2");
    run("ClearAll[{caL1, caL2}]");
    assert_eval_eq("caL1", "caL1", 0);
    assert_eval_eq("caL2", "caL2", 0);
}

static void test_clearall_skips_protected(void) {
    run("caP = 7");
    run("Protect[caP]");
    run("ClearAll[caP]");
    /* Protected -> ClearAll must not touch it. */
    assert_eval_eq("caP", "7", 0);
    run("Unprotect[caP]");
    run("ClearAll[caP]");
    assert_eval_eq("caP", "caP", 0);
}

static void test_clearall_does_not_break_builtin(void) {
    /* Plus is Protected; ClearAll[Plus] must be a no-op. */
    run("ClearAll[Plus]");
    assert_eval_eq("1 + 2", "3", 0);
}

static void test_clearall_returns_null(void) {
    assert_eval_eq("ClearAll[caNeverDefined]", "Null", 0);
}

/* ------------------------------------------------------------------------
 *  Remove
 * ---------------------------------------------------------------------- */

static void test_remove_symbol_value(void) {
    run("rmX = 2");
    assert_eval_eq("rmX", "2", 0);
    assert_eval_eq("Remove[rmX]", "Null", 0);
    assert_eval_eq("rmX", "rmX", 0);
}

static void test_remove_function(void) {
    run("rmF[x_]:=x + 1");
    assert_eval_eq("rmF[4]", "5", 0);
    run("Remove[rmF]");
    assert_eval_eq("rmF[4]", "rmF[4]", 0);
    assert_eval_eq("DownValues[rmF]", "{}", 0);
}

static void test_remove_recreates_clean(void) {
    /* After removal the name is reusable as a fresh, empty symbol. */
    run("rmR = 99");
    run("Remove[rmR]");
    assert_eval_eq("Attributes[rmR]", "{}", 0);
    run("rmR[y_]:=y^2");
    assert_eval_eq("rmR[5]", "25", 0);
}

static void test_remove_list_spec(void) {
    run("rmA = 1");
    run("rmB = 2");
    run("Remove[{rmA, rmB}]");
    assert_eval_eq("rmA", "rmA", 0);
    assert_eval_eq("rmB", "rmB", 0);
}

static void test_remove_skips_protected_builtin(void) {
    /* Times is Protected -> Remove must not delete it. */
    run("Remove[Times]");
    assert_eval_eq("3 * 4", "12", 0);
}

/* ------------------------------------------------------------------------
 *  Protect
 * ---------------------------------------------------------------------- */

static void test_protect_returns_changed_names(void) {
    run("ClearAll[prA, prB]");
    run("Unprotect[prA, prB]");
    assert_eval_eq("Protect[prA, prB]", "{\"prA\", \"prB\"}", 0);
    run("Unprotect[prA, prB]");
}

static void test_protect_already_protected_returns_empty(void) {
    run("Unprotect[prC]");
    run("Protect[prC]");
    /* Re-protecting an already-Protected symbol changes nothing. */
    assert_eval_eq("Protect[prC]", "{}", 0);
    run("Unprotect[prC]");
}

static void test_protect_sets_attribute(void) {
    run("Unprotect[prD]");
    run("Protect[prD]");
    assert_eval_eq("Attributes[prD]", "{Protected}", 0);
    run("Unprotect[prD]");
}

static void test_protect_blocks_redefinition(void) {
    run("prE[x_]:=x^2");
    run("Protect[prE]");
    /* This SetDelayed is blocked (prints prE::wrsym to stderr); the old
     * definition survives. */
    run("prE[x_]:=x^3");
    assert_eval_eq("prE[2]", "4", 0);
    run("Unprotect[prE]");
}

static void test_protect_list_spec(void) {
    run("Unprotect[prL1, prL2]");
    assert_eval_eq("Protect[{prL1, prL2}]", "{\"prL1\", \"prL2\"}", 0);
    run("Unprotect[prL1, prL2]");
}

/* ------------------------------------------------------------------------
 *  Unprotect
 * ---------------------------------------------------------------------- */

static void test_unprotect_returns_changed_names(void) {
    run("Protect[unA]");
    assert_eval_eq("Unprotect[unA]", "{\"unA\"}", 0);
}

static void test_unprotect_not_protected_returns_empty(void) {
    run("ClearAll[unB]");
    assert_eval_eq("Unprotect[unB]", "{}", 0);
}

static void test_unprotect_allows_redefinition(void) {
    run("unC[x_]:=x^2");
    run("Protect[unC]");
    run("unC[x_]:=x^3");     /* blocked */
    assert_eval_eq("unC[2]", "4", 0);
    run("Unprotect[unC]");
    run("unC[x_]:=x^3");     /* same pattern -> replaces in place */
    assert_eval_eq("unC[2]", "8", 0);
}

static void test_unprotect_then_protect_roundtrip(void) {
    run("ClearAll[unD]");
    assert_eval_eq("Protect[unD]", "{\"unD\"}", 0);
    assert_eval_eq("Unprotect[unD]", "{\"unD\"}", 0);
    assert_eval_eq("Attributes[unD]", "{}", 0);
}

/* ------------------------------------------------------------------------
 *  HoldAll behaviour (operate on the symbol, not its value)
 * ---------------------------------------------------------------------- */

static void test_holdall_protect_targets_symbol(void) {
    run("haX = 5");
    /* With HoldAll, Protect sees the symbol haX, not the value 5. */
    assert_eval_eq("Protect[haX]", "{\"haX\"}", 0);
    assert_eval_eq("Attributes[haX]", "{Protected}", 0);
    run("Unprotect[haX]");
    run("ClearAll[haX]");
}

/* ------------------------------------------------------------------------
 *  Head-symbol attributes
 * ---------------------------------------------------------------------- */

static void test_head_attributes(void) {
    assert_eval_eq("Attributes[ClearAll]", "{HoldAll, Protected}", 0);
    assert_eval_eq("Attributes[Remove]", "{HoldAll, Locked, Protected}", 0);
    assert_eval_eq("Attributes[Protect]", "{HoldAll, Protected}", 0);
    assert_eval_eq("Attributes[Unprotect]", "{HoldAll, Protected}", 0);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_define_clear_loop(void) {
    for (int i = 0; i < 200; i++) {
        run("loopSym[x_]:=x^2 + x + 1");
        run("loopSym2 = {1, 2, 3, loopSym}");
        run("ClearAll[loopSym, loopSym2]");
    }
    assert_eval_eq("loopSym[1]", "loopSym[1]", 0);
}

static void test_remove_loop(void) {
    for (int i = 0; i < 200; i++) {
        run("rmLoop = {a, b, c}");
        run("rmLoop2[n_]:=n!");
        run("Remove[rmLoop, rmLoop2]");
    }
    assert_eval_eq("rmLoop", "rmLoop", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* ClearAll */
    TEST(test_clearall_clears_downvalues);
    TEST(test_clearall_clears_ownvalue);
    TEST(test_clearall_clears_attributes);
    TEST(test_clearall_multiple_args);
    TEST(test_clearall_list_spec);
    TEST(test_clearall_skips_protected);
    TEST(test_clearall_does_not_break_builtin);
    TEST(test_clearall_returns_null);

    /* Remove */
    TEST(test_remove_symbol_value);
    TEST(test_remove_function);
    TEST(test_remove_recreates_clean);
    TEST(test_remove_list_spec);
    TEST(test_remove_skips_protected_builtin);

    /* Protect */
    TEST(test_protect_returns_changed_names);
    TEST(test_protect_already_protected_returns_empty);
    TEST(test_protect_sets_attribute);
    TEST(test_protect_blocks_redefinition);
    TEST(test_protect_list_spec);

    /* Unprotect */
    TEST(test_unprotect_returns_changed_names);
    TEST(test_unprotect_not_protected_returns_empty);
    TEST(test_unprotect_allows_redefinition);
    TEST(test_unprotect_then_protect_roundtrip);

    /* HoldAll */
    TEST(test_holdall_protect_targets_symbol);

    /* Head attributes */
    TEST(test_head_attributes);

    /* Memory hygiene */
    TEST(test_define_clear_loop);
    TEST(test_remove_loop);

    printf("All clearall_remove_protect_tests passed.\n");
    return 0;
}

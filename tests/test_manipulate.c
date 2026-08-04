/* Symbolic-construction tests for Manipulate[]. These never call the real
 * Raylib windowing path -- MATHILDA_NO_GRAPHICS_WINDOW forces
 * graphics_manipulate() to no-op regardless of how USE_GRAPHICS resolved
 * for this build, so the suite stays headless everywhere. */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdlib.h>
#include <stdio.h>

void test_manipulate_range_control_returns_null(void) {
    assert_eval_eq("Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]", "Null", 0);
}

void test_manipulate_range_control_with_step_returns_null(void) {
    assert_eval_eq("Manipulate[x + n, {n, 0, 10, 0.5}]", "Null", 0);
}

void test_manipulate_range_control_with_explicit_default_returns_null(void) {
    assert_eval_eq("Manipulate[x + n, {{n, 3}, 0, 10}]", "Null", 0);
}

void test_manipulate_discrete_control_returns_null(void) {
    assert_eval_eq("Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]", "Null", 0);
}

void test_manipulate_discrete_control_with_explicit_default_returns_null(void) {
    assert_eval_eq("Manipulate[f, {{f, Cos[x]}, {Sin[x], Cos[x], x^2}}]", "Null", 0);
}

void test_manipulate_multiple_controls_returns_null(void) {
    assert_eval_eq(
        "Manipulate[a + b, {a, 0, 1}, {b, {1, 2, 3}}]", "Null", 0);
}

void test_manipulate_3d_content_returns_null(void) {
    /* Regression guard for the Graphics2D/Graphics3D content-dispatch
     * split: a Plot3D body must still return Null, not stay unevaluated
     * or hit the text fallback path. */
    assert_eval_eq(
        "Manipulate[Plot3D[Sin[x + n] Cos[y], {x, -3, 3}, {y, -3, 3}], {n, 0, 3}]",
        "Null", 0);
}

void test_manipulate_no_control_args_stays_unevaluated(void) {
    assert_eval_eq("Manipulate[x^2]", "Manipulate[x^2]", 0);
}

void test_manipulate_malformed_control_stays_unevaluated(void) {
    /* Not a symbol for the control variable. */
    assert_eval_eq("Manipulate[x^2, {1, 0, 10}]", "Manipulate[x^2, {1, 0, 10}]", 0);
    /* Empty discrete-value list. */
    assert_eval_eq("Manipulate[x^2, {u, {}}]", "Manipulate[x^2, {u, {}}]", 0);
    /* Wrong arity (single bound). */
    assert_eval_eq("Manipulate[x^2, {u, 0}]", "Manipulate[x^2, {u, 0}]", 0);
}

void test_manipulate_is_holdall(void) {
    /* The body must not be pre-evaluated -- an undefined symbol should
     * reach builtin_manipulate untouched, not trigger evaluation errors. */
    assert_eval_eq("Hold[Manipulate[undefinedSymbolXyz, {u, 0, 1}]]",
                    "Hold[Manipulate[undefinedSymbolXyz, {u, 0, 1}]]", 0);
}

void test_manipulate_attributes_registered(void) {
    assert_eval_eq("MemberQ[Attributes[Manipulate], HoldAll]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Manipulate], Protected]", "True", 0);
}

int main(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    symtab_init();
    core_init();

    TEST(test_manipulate_range_control_returns_null);
    TEST(test_manipulate_range_control_with_step_returns_null);
    TEST(test_manipulate_range_control_with_explicit_default_returns_null);
    TEST(test_manipulate_discrete_control_returns_null);
    TEST(test_manipulate_discrete_control_with_explicit_default_returns_null);
    TEST(test_manipulate_multiple_controls_returns_null);
    TEST(test_manipulate_3d_content_returns_null);
    TEST(test_manipulate_no_control_args_stays_unevaluated);
    TEST(test_manipulate_malformed_control_stays_unevaluated);
    TEST(test_manipulate_is_holdall);
    TEST(test_manipulate_attributes_registered);

    printf("All manipulate tests passed.\n");
    return 0;
}

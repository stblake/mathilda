/* Tests for UnitBox.
 *
 * Cover: interior points, the closed boundary at +-1/2 (both the exact
 * rational and the equivalent Real forms), points outside the box,
 * +-Infinity, an exact symbolic real (Pi) resolved by the same numerical
 * certification UnitStep/Ramp use, symbolic pass-through, the complex
 * (non-real) rejection path and the zero-imaginary-part workaround,
 * Listable threading over a List, argument-count errors, the (Listable,
 * NumericFunction, Protected) attributes -- and NOT Orderless, since UnitBox
 * is unary -- and a memory-hygiene loop.
 *
 * The headline cases mirror the UnitBox docstring (see
 * docs/spec/builtins/elementary-functions.md and
 * docs/spec/changelog/2026-08-17.md).
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Interior points
 * ---------------------------------------------------------------------- */

static void test_unitbox_interior(void) {
    assert_eval_eq("UnitBox[0]", "1", 0);
    assert_eval_eq("UnitBox[0.3]", "1", 0);
    assert_eval_eq("UnitBox[-0.3]", "1", 0);
    assert_eval_eq("UnitBox[1/4]", "1", 0);
    assert_eval_eq("UnitBox[-1/4]", "1", 0);
}

/* ------------------------------------------------------------------------
 *  Closed boundary at +-1/2
 * ---------------------------------------------------------------------- */

static void test_unitbox_boundary_rational(void) {
    assert_eval_eq("UnitBox[1/2]", "1", 0);
    assert_eval_eq("UnitBox[-1/2]", "1", 0);
}

static void test_unitbox_boundary_real(void) {
    assert_eval_eq("UnitBox[0.5]", "1", 0);
    assert_eval_eq("UnitBox[-0.5]", "1", 0);
}

/* ------------------------------------------------------------------------
 *  Outside the box
 * ---------------------------------------------------------------------- */

static void test_unitbox_outside(void) {
    assert_eval_eq("UnitBox[0.6]", "0", 0);
    assert_eval_eq("UnitBox[-0.6]", "0", 0);
    assert_eval_eq("UnitBox[3]", "0", 0);
    assert_eval_eq("UnitBox[-3]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Infinity
 * ---------------------------------------------------------------------- */

static void test_unitbox_infinity(void) {
    assert_eval_eq("UnitBox[Infinity]", "0", 0);
    assert_eval_eq("UnitBox[-Infinity]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Exact symbolic real: certified via the UnitStep/Ramp mechanism
 * ---------------------------------------------------------------------- */

static void test_unitbox_symbolic_exact_real(void) {
    /* Pi ~ 3.14159 is well outside the box; certifiable at low precision. */
    assert_eval_eq("UnitBox[Pi]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Symbolic pass-through
 * ---------------------------------------------------------------------- */

static void test_unitbox_symbolic_passthrough(void) {
    assert_eval_eq("UnitBox[x]", "UnitBox[x]", 0);
}

/* ------------------------------------------------------------------------
 *  Complex rejection and the zero-imaginary-part workaround
 * ---------------------------------------------------------------------- */

static void test_unitbox_rejects_complex(void) {
    assert_eval_eq("UnitBox[1 + I]", "UnitBox[1 + I]", 0);
}

static void test_unitbox_complex_zero_imag(void) {
    /* A Complex whose imaginary part is exactly zero resolves via its real
     * part -- same as UnitStep/Ramp. */
    assert_eval_eq("UnitBox[0.3 + 0 I]", "1", 0);
}

/* ------------------------------------------------------------------------
 *  Listable threading
 * ---------------------------------------------------------------------- */

static void test_unitbox_threads_over_list(void) {
    assert_eval_eq("UnitBox[{-1, -0.5, 0, 0.5, 1}]", "{0, 1, 1, 1, 0}", 0);
}

/* ------------------------------------------------------------------------
 *  Argument-count errors
 * ---------------------------------------------------------------------- */

static void test_unitbox_no_args_unevaluated(void) {
    assert_eval_eq("UnitBox[]", "UnitBox[]", 0);
}

static void test_unitbox_too_many_args_unevaluated(void) {
    assert_eval_eq("UnitBox[1, 2]", "UnitBox[1, 2]", 0);
}

/* ------------------------------------------------------------------------
 *  Attributes
 * ---------------------------------------------------------------------- */

static void test_unitbox_attributes(void) {
    Expr* parsed = parse_expression("Attributes[UnitBox]");
    Expr* evaluated = evaluate(parsed);
    char* str = expr_to_string(evaluated);
    ASSERT_MSG(strstr(str, "Listable") != NULL,
               "expected Listable in attributes, got: %s", str);
    ASSERT_MSG(strstr(str, "NumericFunction") != NULL,
               "expected NumericFunction in attributes, got: %s", str);
    ASSERT_MSG(strstr(str, "Protected") != NULL,
               "expected Protected in attributes, got: %s", str);
    ASSERT_MSG(strstr(str, "Orderless") == NULL,
               "UnitBox is unary and must NOT be Orderless, got: %s", str);
    free(str);
    expr_free(parsed);
    expr_free(evaluated);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_unitbox_memory_loop(void) {
    /* Exercise every shape exercised above repeatedly so valgrind has many
     * chances to catch a leak in any path, including the two-evaluate()
     * construction inside builtin_unitbox. */
    const char* cases[] = {
        "UnitBox[0]",
        "UnitBox[0.3]",
        "UnitBox[1/2]",
        "UnitBox[-1/2]",
        "UnitBox[0.6]",
        "UnitBox[-3]",
        "UnitBox[Infinity]",
        "UnitBox[-Infinity]",
        "UnitBox[Pi]",
        "UnitBox[x]",
        "UnitBox[1 + I]",
        "UnitBox[0.3 + 0 I]",
        "UnitBox[{-1, -0.5, 0, 0.5, 1}]",
        "UnitBox[]",
        "UnitBox[1, 2]",
        NULL
    };
    for (int rep = 0; rep < 20; rep++) {
        for (int i = 0; cases[i]; i++) {
            Expr* p = parse_expression(cases[i]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
    }
}

/* ------------------------------------------------------------------------
 *  Main
 * ---------------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_unitbox_interior);

    TEST(test_unitbox_boundary_rational);
    TEST(test_unitbox_boundary_real);

    TEST(test_unitbox_outside);

    TEST(test_unitbox_infinity);

    TEST(test_unitbox_symbolic_exact_real);
    TEST(test_unitbox_symbolic_passthrough);

    TEST(test_unitbox_rejects_complex);
    TEST(test_unitbox_complex_zero_imag);

    TEST(test_unitbox_threads_over_list);

    TEST(test_unitbox_no_args_unevaluated);
    TEST(test_unitbox_too_many_args_unevaluated);

    TEST(test_unitbox_attributes);

    TEST(test_unitbox_memory_loop);

    printf("All UnitBox tests passed.\n");
    return 0;
}

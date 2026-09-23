/*
 * test_inactive.c — the Inactive / Activate primitive (M59).
 *
 * Inactive[f] is an inert wrapper for a head: Inactive[f][args] evaluates its args
 * but does not fire f's rules, so an integral can be held symbolic without paying
 * (or spinning on) the integration cascade.  D applies the fundamental theorem of
 * calculus to Inactive[Integrate] WITHOUT evaluating the integral, and Activate
 * reactivates it.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}
static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}
#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* Inactive[Integrate][...] stays inert -- even a non-elementary (normally spinning)
 * integrand is held instantly. */
static void t_inactive_inert(void) {
    ASSERT_TRUE("Head[Inactive[Integrate][1/Sqrt[y Log[y] + 3], y]] === Inactive[Integrate]");
    ASSERT_TRUE("FreeQ[Inactive[Integrate][1/p[y], y], Integrate[__]]");   /* not the active head */
}

/* D applies the fundamental theorem of calculus to an inactive integral without
 * evaluating it: D[Inactive[Integrate][f, u], u] == f; a different variable gives 0
 * when the integrand is free of it. */
static void t_inactive_ftc(void) {
    ASSERT_TRUE("D[Inactive[Integrate][1/p[y], y], y] === 1/p[y]");
    ASSERT_TRUE("D[Inactive[Integrate][1/p[y], y], x] === 0");
    ASSERT_TRUE("D[Inactive[Integrate][1/p[y], y] - x, x] === -1");
    /* the implicit-function rule closes: y' = -Gx/Gy = p */
    ASSERT_TRUE("Module[{G = Inactive[Integrate][1/p[y], y] - x}, "
                "Simplify[-D[G, x]/D[G, y]] === p[y]]");
}

/* Activate reactivates and re-evaluates. */
static void t_inactive_activate(void) {
    ASSERT_TRUE("Activate[Inactive[Integrate][2 y, y]] === y^2");
    ASSERT_TRUE("Activate[3 + Inactive[Integrate][2 y, y]] === 3 + y^2");
    ASSERT_TRUE("Activate[Inactive[Plus][1, 2]] === 3");
}

static void t_inactive_protected(void) {
    ASSERT_TRUE("MemberQ[Attributes[Inactive], Protected]");
    ASSERT_TRUE("MemberQ[Attributes[Activate], Protected]");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_inactive_inert);
    TEST(t_inactive_ftc);
    TEST(t_inactive_activate);
    TEST(t_inactive_protected);

    printf("All Inactive/Activate tests passed.\n");
    return 0;
}

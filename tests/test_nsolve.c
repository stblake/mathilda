/* Tests for NSolve — numerical solutions of equations and systems.
 *
 * Cover: univariate polynomial equations (real/complex roots, multiplicity,
 * Reals filter, MaxRoots, WorkingPrecision, positional precision, auto
 * variable, inexact coefficients); empty and universal solution sets; linear
 * systems (unique / underdetermined / inconsistent); zero-dimensional
 * polynomial systems via the eigenvalue method and the elimination
 * ("Symbolic") method, with method agreement; VerifySolutions; the Protected
 * attribute and docstring; and memory hygiene.
 *
 * Solutions are verified *inside* the language: univariate roots are pulled
 * out with Part[res, k, 1, 2]; system solutions are substituted back into the
 * residual list and the largest |f_i| is bounded.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* Largest |f_i| after substituting every solution of `call` into the residual
 * list `resids` is below `tol`. */
static bool system_ok(const char* resids, const char* call, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Total[Map[Function[s, Total[Abs[N[%s /. s]]]], %s]] < %s",
             resids, call, tol);
    bool ok = lang_true(buf);
    if (!ok) fprintf(stderr, "  system_ok FAIL: %s on %s\n", resids, call);
    return ok;
}
#define ASSERT_SYS(resids, call, tol) \
    ASSERT_MSG(system_ok((resids), (call), (tol)), "residual: %s on %s", (resids), (call))

/* ---------------------------------------------------------------------- */

static void test_univariate_real(void) {
    /* x^5 - 2x + 3 == 0 : five roots, one real. */
    ASSERT_TRUE("Length[NSolve[x^5-2x+3==0, x]] == 5");
    ASSERT_TRUE("Length[NSolve[x^5-2x+3==0, x, Reals]] == 1");
    ASSERT_TRUE("N[Abs[Part[NSolve[x^5-2x+3==0,x,Reals],1,1,2] + 1.42361]] < 1e-4");

    /* Bare expression == 0 and a list-of-one variable. */
    ASSERT_TRUE("Length[NSolve[x^3+3x+1==0, {x}]] == 3");

    /* Solutions substitute back to ~0. */
    ASSERT_SYS("{x^5-2x+3}", "NSolve[x^5-2x+3==0, x]", "1e-8");
}

static void test_univariate_multiplicity(void) {
    /* (x^2-1)(x^4-1) : multiplicities -> 6 roots; reals -> 4. */
    ASSERT_TRUE("Length[NSolve[(x^2-1)(x^4-1)==0, x]] == 6");
    ASSERT_TRUE("Length[NSolve[(x^2-1)(x^4-1)==0, x, Reals]] == 4");
    /* (x-1)^2 (x-2)^3 : five solutions repeated by multiplicity. */
    ASSERT_TRUE("Length[NSolve[(x-1)^2 (x-2)^3==0, x]] == 5");
}

static void test_univariate_options(void) {
    /* MaxRoots caps the count. */
    ASSERT_TRUE("Length[NSolve[x^7-2==0, x, MaxRoots->3]] == 3");

    /* WorkingPrecision via option and as a trailing positional argument. */
    ASSERT_TRUE("N[Abs[Part[NSolve[x^2-2==0,x,Reals,WorkingPrecision->30],1,1,2]+Sqrt[2]],40] < 10^-25");
    ASSERT_TRUE("N[Abs[Part[NSolve[x^5-x+2,x,Reals,30],1,1,2]+1.26716830454212431725],40] < 10^-15");

    /* Inexact coefficients. */
    ASSERT_TRUE("Length[NSolve[x^3-1.234x+5.678==0, x]] == 3");
    ASSERT_SYS("{x^3-1.234x+5.678}", "NSolve[x^3-1.234x+5.678==0,x]", "1e-6");

    /* Auto variable: NSolve[expr] collects the variable. */
    ASSERT_TRUE("Length[NSolve[x^2+3x+1==0]] == 2");
}

static void test_degenerate(void) {
    /* Contradiction -> {} ; tautology -> {{}}. */
    ASSERT_TRUE("NSolve[x==1 && x==2, x] === {}");
    ASSERT_TRUE("NSolve[x==x, x] === {{}}");
}

static void test_linear_systems(void) {
    /* Unique solution. */
    ASSERT_TRUE("Length[NSolve[{x+2y+3z==4, 3x+4y+5z==6, 7x+9y+8z==10}, {x,y,z}]] == 1");
    ASSERT_SYS("{x+2y+3z-4, 3x+4y+5z-6, 7x+9y+8z-10}",
               "NSolve[{x+2y+3z==4,3x+4y+5z==6,7x+9y+8z==10},{x,y,z}]", "1e-6");

    /* Inconsistent -> {}. */
    ASSERT_TRUE("NSolve[x+2y+3z==4 && 3x+4y+5z==6 && 6x+7y+8z==0, {x,y,z}] === {}");

    /* Underdetermined: a non-empty (parametric) family via the Solve fallback. */
    ASSERT_TRUE("Length[NSolve[x+2y+3z==4 && 3x+4y+5z==6 && 6x+7y+8z==9, {x,y,z}]] >= 1");
}

static void test_poly_systems(void) {
    /* Six complex solutions; all satisfy both equations. */
    ASSERT_TRUE("Length[NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}]] == 6");
    ASSERT_SYS("{x^2+y^2-1, x^3-y^3-2}",
               "NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}]", "1e-6");

    /* Mixed-degree square system. */
    ASSERT_SYS("{x^2+y^3-1, 2x+3y-4}",
               "NSolve[{x^2+y^3==1, 2x+3y==4}, {x,y}]", "1e-6");
    /* Reals restricts to the single real solution. */
    ASSERT_TRUE("Length[NSolve[{x^2+y^3==1, 2x+3y==4}, {x,y}, Reals]] == 1");

    /* Arbitrary precision system: residual far below machine epsilon. */
    ASSERT_SYS("{x^2+y^2-1, x^3-y^3-2}",
               "NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}, WorkingPrecision->25]", "10^-20");
}

static void test_method_agreement(void) {
    /* Every method returns the same number of solutions, all valid. */
    const char* sys = "{x^2+2y^2==3, x^3-4x y==5}";
    const char* res = "{x^2+2y^2-3, x^3-4x y-5}";
    char buf[1024];

    const char* methods[] = { "", ", Method->\"EndomorphismMatrix\"",
                              ", Method->\"Homotopy\"", ", Method->\"Symbolic\"" };
    for (int i = 0; i < 4; i++) {
        snprintf(buf, sizeof buf, "Length[NSolve[%s, {x,y}%s]] == 6", sys, methods[i]);
        ASSERT_MSG(lang_true(buf), "method count: %s", methods[i]);
        snprintf(buf, sizeof buf, "NSolve[%s, {x,y}%s]", sys, methods[i]);
        ASSERT_SYS(res, buf, "1e-6");
    }
}

static void test_verify_solutions(void) {
    /* VerifySolutions->False still returns the solutions for a clean system. */
    ASSERT_TRUE("Length[NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}, VerifySolutions->False]] >= 6");
    ASSERT_TRUE("Length[NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}, VerifySolutions->True]] == 6");
}

static void test_transcendental_seed(void) {
    /* E^x - x == 7 has two real roots; Solve cannot reduce it, so the
     * FindRoot grid-seeding fallback supplies them. */
    ASSERT_TRUE("Length[NSolve[E^x-x==7, x, Reals]] == 2");
    ASSERT_SYS("{E^x-x-7}", "NSolve[E^x-x==7, x, Reals]", "1e-6");
}

static void test_radical_verification(void) {
    /* Radical substitution (t = x^(1/6)) introduces extraneous complex roots
     * that do not satisfy the original equation; they must be filtered out, so
     * only the single real root x ~ 1.80863 survives. */
    ASSERT_TRUE("Length[NSolve[Sqrt[x] + 3 x^(1/3) == 5, x]] == 1");
    ASSERT_TRUE("Length[NSolve[Sqrt[x] + 3 x^(1/3) == 5, x, Reals]] == 1");
    ASSERT_TRUE("N[Abs[Part[NSolve[Sqrt[x]+3 x^(1/3)==5,x],1,1,2] - 1.80863]] < 1e-4");
    /* Every returned root genuinely satisfies the original equation. */
    ASSERT_SYS("{Sqrt[x] + 3 x^(1/3) - 5}", "NSolve[Sqrt[x]+3 x^(1/3)==5,x]", "1e-6");
}

static void test_degree_and_lcm_guards(void) {
    /* Enormous polynomial degree: leave unevaluated rather than hang. */
    ASSERT_TRUE("Head[NSolve[x^1000000-2 x+3==0, x, MaxRoots->3]] === NSolve");
    /* Absurd radical denominator (x^(p/67890)): the radical path bails on the
     * LCM cap, so NSolve must not hang (here the grid seeder still finds real
     * roots that all satisfy the equation). */
    ASSERT_TRUE("Head[NSolve[2 x^(123451/67890) - x^2 + 4 Sqrt[x] - 4 x - 9/8 == 0, "
                "x, Reals]] === List");
}

static void test_attributes(void) {
    ASSERT_TRUE("MemberQ[Attributes[NSolve], Protected]");
}

static void test_memory_loop(void) {
    /* Exercise the univariate, system, and elimination paths repeatedly. */
    for (int i = 0; i < 30; i++) {
        free(eval_str("NSolve[x^4-3x+1==0, x]"));
        free(eval_str("NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}]"));
        free(eval_str("NSolve[{x^2+2y^2==3, x^3-4x y==5}, {x,y}, Method->\"Symbolic\"]"));
        free(eval_str("NSolve[x==1 && x==2, x]"));
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_univariate_real);
    TEST(test_univariate_multiplicity);
    TEST(test_univariate_options);
    TEST(test_degenerate);
    TEST(test_linear_systems);
    TEST(test_poly_systems);
    TEST(test_method_agreement);
    TEST(test_verify_solutions);
    TEST(test_radical_verification);
    TEST(test_degree_and_lcm_guards);
    TEST(test_transcendental_seed);
    TEST(test_attributes);
    TEST(test_memory_loop);

    printf("All nsolve_tests passed.\n");
    return 0;
}

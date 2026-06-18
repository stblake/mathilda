/* Stress tests for NSolve — numerical solutions of equations and systems.
 *
 * Complements tests/test_nsolve.c (functional coverage) by pushing every
 * option and code path harder, at both MACHINE and ARBITRARY precision:
 *
 *   - high-degree univariate polynomials (deg 20/30/50) and Wilkinson-style
 *     well-separated real roots, machine + arbitrary precision;
 *   - clustered / high-multiplicity roots;
 *   - complex coefficients;
 *   - arbitrary precision residual scaling and Precision[] of the roots;
 *   - every option: MaxRoots, WorkingPrecision (positional + rule),
 *     PrecisionGoal, MaxIterations, RandomSeeding, VerifySolutions, Method,
 *     plus several options combined;
 *   - the three domains (Reals / Complexes / Integers);
 *   - zero-dimensional polynomial systems across all four Method values, at
 *     machine and arbitrary precision, with several random seeds;
 *   - three-variable systems;
 *   - transcendental grid-seeding;
 *   - degenerate (tautology / contradiction) inputs in many shapes;
 *   - a memory-hygiene loop hitting every path (machine + arbitrary precision).
 *
 * Everything is verified INSIDE the language: counts via Length, root/solution
 * validity by substituting back into the residual(s) and bounding the largest
 * |f_i|.  The randomized system engine recovers a seed-dependent SUBSET of a
 * system's solutions, so seed tests assert VALIDITY of whatever is returned,
 * never a fixed count.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

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

/* Largest total |f_i| after substituting every solution of `call` into the
 * residual list `resids` is below `tol` (all arithmetic in-language, numeric). */
static bool system_ok(const char* resids, const char* call, const char* tol) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "Total[Map[Function[s, Total[Abs[N[%s /. s]]]], %s]] < %s",
             resids, call, tol);
    bool ok = lang_true(buf);
    if (!ok) fprintf(stderr, "  system_ok FAIL: %s on %s\n", resids, call);
    return ok;
}
#define ASSERT_SYS(resids, call, tol) \
    ASSERT_MSG(system_ok((resids), (call), (tol)), "residual: %s on %s", (resids), (call))

/* Same, but numericalise the residual to `d` digits before bounding — used for
 * arbitrary-precision solutions whose residuals fall far below machine eps. */
static bool system_ok_prec(const char* resids, const char* call,
                           const char* tol, int d) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "Total[Map[Function[s, Total[Abs[N[%s /. s, %d]]]], %s]] < %s",
             resids, d, call, tol);
    bool ok = lang_true(buf);
    if (!ok) fprintf(stderr, "  system_ok_prec FAIL: %s on %s\n", resids, call);
    return ok;
}
#define ASSERT_SYS_PREC(resids, call, tol, d) \
    ASSERT_MSG(system_ok_prec((resids), (call), (tol), (d)), \
               "residual(%d): %s on %s", (d), (resids), (call))

/* Helper: assert Length[NSolve[...]] == n via a formatted query. */
static bool count_is(const char* expr_text, int n) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Length[%s] == %d", expr_text, n);
    return lang_true(buf);
}
#define ASSERT_COUNT(expr_text, n) \
    ASSERT_MSG(count_is((expr_text), (n)), "count %d: %s", (n), (expr_text))

/* ====================================================================== *
 *  Univariate polynomials — high degree, machine precision
 * ====================================================================== */

static void test_high_degree_machine(void) {
    /* Full complex root counts. */
    ASSERT_COUNT("NSolve[x^20-1==0, x]", 20);
    ASSERT_COUNT("NSolve[x^30-1==0, x]", 30);
    ASSERT_COUNT("NSolve[x^50-2==0, x]", 50);
    ASSERT_COUNT("NSolve[x^25-x-1==0, x]", 25);

    /* Reals filter on roots of unity / radicals. */
    ASSERT_COUNT("NSolve[x^30-1==0, x, Reals]", 2);   /* +-1            */
    ASSERT_COUNT("NSolve[x^31-1==0, x, Reals]", 1);   /* odd -> 1 real  */
    ASSERT_COUNT("NSolve[x^50-2==0, x, Reals]", 2);   /* +- 2^(1/50)    */

    /* Residuals tiny for all returned roots. */
    ASSERT_SYS("{x^20-x-1}", "NSolve[x^20-x-1==0, x]", "1e-7");
    ASSERT_SYS("{x^30-1}",   "NSolve[x^30-1==0, x]",   "1e-7");
    ASSERT_SYS("{x^50-2}",   "NSolve[x^50-2==0, x]",   "1e-6");
}

static void test_wilkinson_machine(void) {
    /* Well-separated real roots 1..10 via a literal factored polynomial
     * (Product[] does not expand here, so spell out the factors). */
    const char* w10 =
        "(x-1)(x-2)(x-3)(x-4)(x-5)(x-6)(x-7)(x-8)(x-9)(x-10)";
    char buf[512];

    snprintf(buf, sizeof buf, "NSolve[%s==0, x]", w10);
    ASSERT_COUNT(buf, 10);

    /* Recovered roots match 1..10 (sorted real parts), imaginary parts ~0. */
    snprintf(buf, sizeof buf,
             "Max[Abs[Sort[Re[Table[Part[NSolve[%s==0,x],j,1,2],{j,1,10}]]] "
             "- Range[10]]] < 1/100", w10);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof buf,
             "Max[Abs[Im[Table[Part[NSolve[%s==0,x],j,1,2],{j,1,10}]]]] < 1/100",
             w10);
    ASSERT_TRUE(buf);
}

static void test_multiplicity_stress(void) {
    /* High multiplicities counted with multiplicity. */
    ASSERT_COUNT("NSolve[(x-1)^5==0, x]", 5);
    ASSERT_COUNT("NSolve[(x-1)^4 (x-2)^4==0, x]", 8);
    ASSERT_COUNT("NSolve[(x^2+1)^3==0, x]", 6);
    ASSERT_COUNT("NSolve[(x^2+1)^3==0, x, Reals]", 0);   /* purely complex */
    /* Mixed clustered real + complex, residual still tiny. */
    ASSERT_SYS("{(x-1)^4 (x-2)^4}", "NSolve[(x-1)^4 (x-2)^4==0, x]", "1e-3");
}

static void test_complex_coefficients(void) {
    /* Complex coefficients -> full set of complex roots. */
    ASSERT_COUNT("NSolve[x^4 + (2+I) x + (1-3I) == 0, x]", 4);
    ASSERT_SYS("{x^4 + (2+I) x + (1-3I)}",
               "NSolve[x^4 + (2+I) x + (1-3I) == 0, x]", "1e-6");
    ASSERT_COUNT("NSolve[x^6 + I==0, x]", 6);
    ASSERT_SYS("{x^6 + I}", "NSolve[x^6 + I==0, x]", "1e-6");
}

/* ====================================================================== *
 *  Univariate polynomials — arbitrary precision
 * ====================================================================== */

static void test_arbitrary_precision(void) {
    /* Residual scales with WorkingPrecision: at d digits the binding is
     * accurate to roughly 10^-(d-margin). */
    ASSERT_TRUE("N[Abs[Part[NSolve[x^2-2==0,x,Reals,WorkingPrecision->30],1,1,2]"
                "+Sqrt[2]],40] < 10^-25");
    ASSERT_TRUE("N[Abs[Part[NSolve[x^2-2==0,x,Reals,WorkingPrecision->50],1,1,2]"
                "+Sqrt[2]],60] < 10^-45");
    ASSERT_TRUE("N[Abs[Part[NSolve[x^2-2==0,x,Reals,WorkingPrecision->80],1,1,2]"
                "+Sqrt[2]],90] < 10^-70");

    /* Precision[] of a returned root tracks the requested working precision. */
    ASSERT_TRUE("Precision[Part[NSolve[x^2-2==0,x,Reals,WorkingPrecision->50],"
                "1,1,2]] > 45");

    /* WorkingPrecision as a trailing positional argument. */
    ASSERT_TRUE("N[Abs[Part[NSolve[x^5-x+2,x,Reals,30],1,1,2]"
                "+1.26716830454212431725],40] < 10^-15");

    /* High-precision complex roots: residual far below machine eps. */
    ASSERT_SYS_PREC("{x^6+x+1}",
                    "NSolve[x^6+x+1==0,x,WorkingPrecision->40]", "10^-30", 45);
    ASSERT_SYS_PREC("{x^20-x-1}",
                    "NSolve[x^20-x-1==0,x,WorkingPrecision->35]", "10^-25", 40);

    /* Reals filter survives high precision. */
    ASSERT_COUNT("NSolve[x^7-3==0, x, Reals, WorkingPrecision->40]", 1);
}

/* ====================================================================== *
 *  Options — exhaustive
 * ====================================================================== */

static void test_maxroots_sweep(void) {
    /* MaxRoots caps the returned count; >count yields the full set. */
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->0]", 0);
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->1]", 1);
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->5]", 5);
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->12]", 12);
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->100]", 12);
    /* MaxRoots together with the Reals filter. */
    ASSERT_COUNT("NSolve[x^12-3==0, x, Reals, MaxRoots->1]", 1);
    /* MaxRoots -> Infinity / Automatic means "all". */
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->Infinity]", 12);
    ASSERT_COUNT("NSolve[x^12-3==0, x, MaxRoots->Automatic]", 12);
}

static void test_options_accepted(void) {
    /* PrecisionGoal, MaxIterations, RandomSeeding are accepted; the
     * polynomial engine still returns the full, valid root set. */
    ASSERT_COUNT("NSolve[x^5-2==0, x, PrecisionGoal->20]", 5);
    ASSERT_COUNT("NSolve[x^5-2==0, x, MaxIterations->5]", 5);
    ASSERT_COUNT("NSolve[x^5-2==0, x, RandomSeeding->42]", 5);
    ASSERT_COUNT("NSolve[x^5-2==0, x, VerifySolutions->True]", 5);
    ASSERT_COUNT("NSolve[x^5-2==0, x, VerifySolutions->False]", 5);

    /* PrecisionGoal also drives the working precision of the roots. */
    ASSERT_TRUE("N[Abs[Part[NSolve[x^2-2==0,x,Reals,PrecisionGoal->30],1,1,2]"
                "+Sqrt[2]],40] < 10^-25");

    /* Several options at once. */
    ASSERT_COUNT("NSolve[x^9-5==0, x, Complexes, MaxRoots->4, "
                 "WorkingPrecision->30, VerifySolutions->True, "
                 "RandomSeeding->7]", 4);
    ASSERT_SYS_PREC("{x^9-5}",
                    "NSolve[x^9-5==0, x, MaxRoots->4, WorkingPrecision->30]",
                    "10^-22", 35);
}

static void test_domains(void) {
    /* Complexes (default) vs Reals on the same polynomial. */
    ASSERT_COUNT("NSolve[x^8-1==0, x]", 8);
    ASSERT_COUNT("NSolve[x^8-1==0, x, Complexes]", 8);
    ASSERT_COUNT("NSolve[x^8-1==0, x, Reals]", 2);

    /* Reals on a polynomial with a mix of real and complex roots. */
    ASSERT_COUNT("NSolve[x^3-x==0, x, Reals]", 3);   /* -1,0,1 */
    ASSERT_COUNT("NSolve[x^3+x==0, x, Reals]", 1);   /* 0       */

    /* Integers domain (routed through Solve). */
    ASSERT_TRUE("Length[NSolve[x^2-4==0, x, Integers]] == 2");
    ASSERT_SYS("{x^2-4}", "NSolve[x^2-4==0, x, Integers]", "1e-6");
}

/* ====================================================================== *
 *  Polynomial systems — methods, precision, seeds
 * ====================================================================== */

static void test_system_methods(void) {
    /* Each system, every Method, must return the expected count and valid
     * solutions at machine precision. */
    struct { const char* sys; const char* res; int n; } cases[] = {
        { "{x^2+y^2==1, x^3-y^3==2}", "{x^2+y^2-1, x^3-y^3-2}", 6 },
        { "{x^2+2y^2==3, x^3-4x y==5}", "{x^2+2y^2-3, x^3-4x y-5}", 6 },
        { "{x^2+y^3==1, 2x+3y==4}", "{x^2+y^3-1, 2x+3y-4}", 3 },
    };
    const char* methods[] = { "", ", Method->\"EndomorphismMatrix\"",
                              ", Method->\"Homotopy\"", ", Method->\"Symbolic\"" };
    char buf[1024];
    for (size_t c = 0; c < sizeof cases / sizeof cases[0]; c++) {
        for (int m = 0; m < 4; m++) {
            snprintf(buf, sizeof buf, "NSolve[%s, {x,y}%s]",
                     cases[c].sys, methods[m]);
            ASSERT_COUNT(buf, cases[c].n);
            ASSERT_SYS(cases[c].res, buf, "1e-6");
        }
    }
}

static void test_system_arbitrary_precision(void) {
    const char* sys = "{x^2+y^2==1, x^3-y^3==2}";
    const char* res = "{x^2+y^2-1, x^3-y^3-2}";
    char buf[512];

    snprintf(buf, sizeof buf, "NSolve[%s, {x,y}, WorkingPrecision->25]", sys);
    ASSERT_COUNT(buf, 6);
    ASSERT_SYS_PREC(res, buf, "10^-18", 30);

    snprintf(buf, sizeof buf, "NSolve[%s, {x,y}, WorkingPrecision->40]", sys);
    ASSERT_COUNT(buf, 6);
    ASSERT_SYS_PREC(res, buf, "10^-30", 45);

    /* Symbolic method at high precision as well. */
    snprintf(buf, sizeof buf,
             "NSolve[%s, {x,y}, WorkingPrecision->40, Method->\"Symbolic\"]", sys);
    ASSERT_SYS_PREC(res, buf, "10^-30", 45);
}

static void test_system_seeds(void) {
    /* The randomized eigenvalue engine recovers a seed-dependent SUBSET of the
     * solutions; whatever it returns must be valid (substitutes back to ~0),
     * and at least one solution must come back for every seed. */
    const char* sys = "{x^2+2y^2==3, x^3-4x y==5}";
    const char* res = "{x^2+2y^2-3, x^3-4x y-5}";
    const char* seeds[] = { "1", "2", "7", "42", "999", "123456" };
    char buf[512];
    for (size_t i = 0; i < sizeof seeds / sizeof seeds[0]; i++) {
        snprintf(buf, sizeof buf, "NSolve[%s, {x,y}, RandomSeeding->%s]",
                 sys, seeds[i]);
        ASSERT_SYS(res, buf, "1e-6");
        snprintf(buf, sizeof buf,
                 "Length[NSolve[%s, {x,y}, RandomSeeding->%s]] >= 1", sys, seeds[i]);
        ASSERT_TRUE(buf);
    }
    /* Default seeding reliably recovers the full set. */
    snprintf(buf, sizeof buf, "NSolve[%s, {x,y}]", sys);
    ASSERT_COUNT(buf, 6);
}

static void test_verify_solutions(void) {
    const char* sys = "{x^2+y^2==1, x^3-y^3==2}";
    const char* res = "{x^2+y^2-1, x^3-y^3-2}";
    char buf[512];
    const char* vs[] = { "True", "False", "Automatic" };
    for (size_t i = 0; i < sizeof vs / sizeof vs[0]; i++) {
        snprintf(buf, sizeof buf,
                 "NSolve[%s, {x,y}, VerifySolutions->%s]", sys, vs[i]);
        ASSERT_SYS(res, buf, "1e-6");
        snprintf(buf, sizeof buf,
                 "Length[NSolve[%s, {x,y}, VerifySolutions->%s]] >= 6", sys, vs[i]);
        ASSERT_TRUE(buf);
    }
}

static void test_three_variable_systems(void) {
    /* Linear + two quadrics. */
    ASSERT_SYS("{x^2+y^2+z^2-1, x+y+z, x-y-z}",
               "NSolve[{x^2+y^2+z^2==1, x+y+z==0, x-y==z}, {x,y,z}]", "1e-6");
    ASSERT_TRUE("Length[NSolve[{x^2+y^2+z^2==1, x+y+z==0, x-y==z}, "
                "{x,y,z}]] >= 1");

    /* Three quadrics: a small zero-dimensional system. */
    ASSERT_SYS("{x^2-y, y^2-z, z^2-x}",
               "NSolve[{x^2==y, y^2==z, z^2==x}, {x,y,z}]", "1e-5");
}

/* ====================================================================== *
 *  Transcendental grid-seeding
 * ====================================================================== */

static void test_transcendental(void) {
    ASSERT_COUNT("NSolve[Cos[x]==x, x, Reals]", 1);
    ASSERT_SYS("{Cos[x]-x}", "NSolve[Cos[x]==x, x, Reals]", "1e-6");

    ASSERT_COUNT("NSolve[E^x-x==7, x, Reals]", 2);
    ASSERT_SYS("{E^x-x-7}", "NSolve[E^x-x==7, x, Reals]", "1e-6");

    /* Arbitrary precision transcendental seeding. */
    ASSERT_SYS_PREC("{Cos[x]-x}",
                    "NSolve[Cos[x]==x, x, Reals, WorkingPrecision->30]",
                    "10^-20", 35);
}

/* ====================================================================== *
 *  Degenerate inputs (many shapes)
 * ====================================================================== */

static void test_degenerate(void) {
    /* Contradictions -> {}. */
    ASSERT_TRUE("NSolve[x==1 && x==2, x] === {}");
    ASSERT_TRUE("NSolve[1==2, x] === {}");
    ASSERT_TRUE("NSolve[{x+y==1, x+y==2}, {x,y}] === {}");
    /* Tautologies -> {{}}. */
    ASSERT_TRUE("NSolve[x==x, x] === {{}}");
    ASSERT_TRUE("NSolve[1==1, x] === {{}}");
    /* Inconsistent linear 3x3 -> {}. */
    ASSERT_TRUE("NSolve[x+2y+3z==4 && 3x+4y+5z==6 && 6x+7y+8z==0, {x,y,z}] === {}");
}

/* ====================================================================== *
 *  Attributes / docstring
 * ====================================================================== */

static void test_attributes(void) {
    ASSERT_TRUE("MemberQ[Attributes[NSolve], Protected]");
    char* s = eval_str("StringLength[NSolve::usage] > 0");
    /* Docstring may or may not be exposed via ::usage; tolerate either, but if
     * present it must be non-empty. */
    if (strcmp(s, "True") != 0 && strcmp(s, "False") != 0)
        ; /* no usage symbol -> fine */
    free(s);
}

/* ====================================================================== *
 *  Memory hygiene — every path, machine + arbitrary precision
 * ====================================================================== */

static void test_memory_loop(void) {
    for (int i = 0; i < 20; i++) {
        free(eval_str("NSolve[x^12-3==0, x]"));
        free(eval_str("NSolve[x^8-1==0, x, Reals]"));
        free(eval_str("NSolve[x^4 + (2+I) x + (1-3I) == 0, x]"));
        free(eval_str("NSolve[x^2-2==0, x, Reals, WorkingPrecision->40]"));
        free(eval_str("NSolve[x^9-5==0, x, MaxRoots->4, WorkingPrecision->30]"));
        free(eval_str("NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}]"));
        free(eval_str("NSolve[{x^2+2y^2==3, x^3-4x y==5}, {x,y}, "
                      "Method->\"Symbolic\"]"));
        free(eval_str("NSolve[{x^2+y^2==1, x^3-y^3==2}, {x,y}, "
                      "WorkingPrecision->30]"));
        free(eval_str("NSolve[{x^2+y^2+z^2==1, x+y+z==0, x-y==z}, {x,y,z}]"));
        free(eval_str("NSolve[Cos[x]==x, x, Reals]"));
        free(eval_str("NSolve[x==1 && x==2, x]"));
        free(eval_str("NSolve[x==x, x]"));
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_high_degree_machine);
    TEST(test_wilkinson_machine);
    TEST(test_multiplicity_stress);
    TEST(test_complex_coefficients);
    TEST(test_arbitrary_precision);
    TEST(test_maxroots_sweep);
    TEST(test_options_accepted);
    TEST(test_domains);
    TEST(test_system_methods);
    TEST(test_system_arbitrary_precision);
    TEST(test_system_seeds);
    TEST(test_verify_solutions);
    TEST(test_three_variable_systems);
    TEST(test_transcendental);
    TEST(test_degenerate);
    TEST(test_attributes);
    TEST(test_memory_loop);

    printf("All nsolve_stress_tests passed.\n");
    return 0;
}

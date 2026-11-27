/* test_radical_field_simplify.c -- the radical function-field normal form
 * (nfrat_radical_normal_form, Simplify's simp_algebraic_field) and the
 * Together / Cancel kernel layer for radicals next to constant algebraic atoms
 * (nfrat_radical_kernels).
 *
 *   1. The reported case: D of a GeneralisedGoursat antiderivative of
 *      (x - 1 - I Sqrt[3])/((x + 1 + I Sqrt[3]) (x^3 + 8)^(1/3)), a sum of
 *      fractions over Q(2^(1/3), (-1)^(1/3), Sqrt[3]) and (8 + x^3)^(k/3),
 *      simplifies back to the integrand, quickly.  The derivative is stored in
 *      data/radical_field/goursat_cbrt_repro_derivative.m (argv[1] is its
 *      path), so the test does not depend on the integrator.
 *   2. Together returns ONE fraction for a sum over radicals of polynomials
 *      and constant radicals (the radicals as independent kernels, as in the
 *      reference); Cancel of that fraction is one fraction too.
 *   3. Several radical classes, radicands equal up to a positive rational,
 *      bare-symbol radicands, several variables.
 *   4. Branch traps: nothing that is false on part of the plane is ever
 *      returned (Sqrt[x] Sqrt[1 + x] vs Sqrt[x + x^2], Sqrt[-x] vs I Sqrt[x],
 *      (x^3)^(1/3) vs x, (-8 - x^3)^(1/3) vs -(8 + x^3)^(1/3)).
 *   5. Simplify never returns a larger result than its input.
 *
 * Expected values are derived from first principles (see the comments) and
 * every simplification is also checked numerically at complex points, off the
 * branch cuts.  Uses ASSERT_MSG (never libc assert: Release defines NDEBUG).
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static char repro_path[1200] = "data/radical_field/goursat_cbrt_repro_derivative.m";

static char* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_inputform_display(res);
    expr_free(res);
    expr_free(e);
    return s;
}

/* Evaluate `input`; assert its InputForm equals `expected`. */
static void run_test(const char* input, const char* expected) {
    char* s = eval_str(input);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
}

static void run_true(const char* input) { run_test(input, "True"); }

/* f and Simplify[f] agree at two complex points (30 digits), and the result
 * is no larger than f. */
static void check_simplify_sound(const char* f) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{ff = %s, ss}, ss = Simplify[ff]; "
             "LeafCount[ss] <= LeafCount[ff] && "
             "Abs[N[(ss - ff) /. x -> N[3/7 + 2 I/5, 40]]] < 10^-25 && "
             "Abs[N[(ss - ff) /. x -> N[-5/3 + I/9, 40]]] < 10^-25]", f);
    run_true(buf);
}

/* ---- 1. the reported case ---- */

static void test_repro(void) {
    char buf[2048];
    /* The integrand is (x - 1 - I Sqrt[3])/((x + 1 + I Sqrt[3]) y),
     * y = (8 + x^3)^(1/3).  With w = (-1)^(1/3) = 1/2 + I Sqrt[3]/2 the
     * constants are 1 + I Sqrt[3] = 2 w, so the integrand is
     * (x - 2 w)/((x + 2 w) y): the normal form Sum c_m y^m has the single
     * term m = 2, c_2 = (x - 2 w)/((x + 2 w)(8 + x^3)), whose radicand factor
     * is absorbed (c_2 y^2 = c_2 (8 + x^3) y^-1). */
    snprintf(buf, sizeof buf,
             "dF = Get[\"%s\"]; f = (x - 1 - Sqrt[3]*I)/(x + 1 + Sqrt[3]*I)/(x^3 + 8)^(1/3); "
             "LeafCount[dF] == 902", repro_path);
    run_true(buf);
    run_test("Simplify[dF]", "(-2*(-1)^(1/3) + x)/((2*(-1)^(1/3) + x)*(8 + x^3)^(1/3))");
    run_test("Simplify[Simplify[dF] - f]", "0");
    run_test("Simplify[dF - f]", "0");
    run_true("Abs[N[(Simplify[dF] - f) /. x -> N[7/5 - I/3, 40]]] < 10^-25");
    run_true("LeafCount[Simplify[dF]] <= LeafCount[f]");
    /* FullSimplify reaches the same form through the same transform */
    run_test("FullSimplify[dF]", "(-2*(-1)^(1/3) + x)/((2*(-1)^(1/3) + x)*(8 + x^3)^(1/3))");
}

/* The reported 3.3 s (7.4 s on a loaded machine) is ~0.15 s now; the bound
 * leaves room for a loaded CI runner. */
static void test_repro_timing(void) {
    Expr* e = parse_expression("First[AbsoluteTiming[Simplify[dF]]]");
    Expr* r = evaluate(e);
    double t = r->type == EXPR_REAL ? r->data.real : 99.0;
    printf("  Simplify[dF]: %.3f s\n", t);
    ASSERT_MSG(t < 2.0, "Simplify[dF] took %.3f s (bound 2 s)", t);
    expr_free(r); expr_free(e);
    e = parse_expression("First[AbsoluteTiming[Together[dF]]]");
    r = evaluate(e);
    t = r->type == EXPR_REAL ? r->data.real : 99.0;
    printf("  Together[dF]: %.3f s\n", t);
    ASSERT_MSG(t < 2.0, "Together[dF] took %.3f s (bound 2 s)", t);
    expr_free(r); expr_free(e);
}

/* ---- 2. Together / Cancel: one fraction ---- */

static void test_together_one_fraction(void) {
    /* The reference: Together of a sum of fractions over radicals of
     * polynomials and constant radicals is ONE fraction (Times[num, den^-1]). */
    run_test("Head[Together[dF]]", "Times");
    run_test("Head[Cancel[Together[dF]]]", "Times");
    run_true("Abs[N[(Together[dF] - f) /. x -> N[7/5 - I/3, 40]]] < 10^-25");
    run_true("Abs[N[(Cancel[Together[dF]] - f) /. x -> N[7/5 - I/3, 40]]] < 10^-25");
    /* the reference: (1 + 2^(1/3))/((2^(1/3) - y) (1 + y)) */
    run_test("Together[1/(1 + (8 + x^3)^(1/3)) + 1/(2^(1/3) - (8 + x^3)^(1/3))]",
             "(-1 - 2^(1/3))/((1 + (8 + x^3)^(1/3))*(-2^(1/3) + (8 + x^3)^(1/3)))");
    run_test("Head[Together[1/(1 + (8 + x^3)^(1/3)) + Sqrt[3]/(1 + (8 + x^3)^(2/3))]]", "Times");
    run_test("Head[Together[(-1)^(1/3)/(x + (1 + x^2)^(1/4)) + 2^(1/3)/(x - Sqrt[1 + x^2])]]", "Times");
    run_test("Head[Together[I/(x + Sqrt[1 - x^2]) + 1/(x + (1 - x)^(1/3))]]", "Times");
    /* a bare-symbol radicand is eliminated, x = (x^(1/2))^2: the common
     * factor Sqrt[x] cancels, as in the reference */
    run_true("Module[{t = Together[1/(x + Sqrt[x]) + Sqrt[2]]}, Head[t] === Times && "
             "Simplify[t - (1 + Sqrt[2] Sqrt[x] + Sqrt[2] x)/(Sqrt[x] (1 + Sqrt[x]))] === 0 && "
             "FreeQ[Denominator[t], x^2]]");
    /* documented semantics unchanged: no rationalisation of a numeric
     * algebraic denominator, no reduction of a Root power */
    run_test("Together[x/(1 + Sqrt[2]) + 1]", "(1 + Sqrt[2] + x)/(1 + Sqrt[2])");
    run_test("Together[Root[#^5 - # - 1 &, 1]^6/(x + 1) + 1/x]",
             "(1 + x + x*Root[#1^5 - #1 - 1 &, 1]^6)/(x*(1 + x))");
    /* the radicals stay in the denominator: Together does not rationalise */
    run_true("!FreeQ[Denominator[Together[1/(1 + (8 + x^3)^(1/3)) + 2^(1/3)/x]], (8 + x^3)^(1/3)]");
    /* Cancel keeps a sum term by term */
    run_test("Head[Cancel[1/(1 + (8 + x^3)^(1/3)) + 2^(1/3)/(1 - (8 + x^3)^(1/3))]]", "Plus");
}

/* ---- 3. several radicals, ratios, variables ---- */

static void test_several_radicals(void) {
    /* (2 x)^(1/3) = 2^(1/3) x^(1/3) for principal values: arg(2 x) = arg(x) */
    run_test("Simplify[(2 x)^(1/3) - 2^(1/3) x^(1/3)]", "0");
    run_test("Simplify[Sqrt[4 x + 4] - 2 Sqrt[1 + x]]", "0");
    /* (8 + x^3)^(4/3) = (8 + x^3) (8 + x^3)^(1/3): an integer power of one radical */
    run_test("Simplify[(8 + x^3)^(4/3) - (8 + x^3) (8 + x^3)^(1/3)]", "0");
    /* two classes and a constant: x Sqrt[1 + x] (1 - x)^(1/3) collected */
    check_simplify_sound("D[x^2 Sqrt[1 + x] (1 - x)^(1/3) + 2^(1/3) Sqrt[1 + x], x]");
    check_simplify_sound("D[(1 + x^3)^(2/3)/x^2 + (-1)^(1/3) (1 + x^3)^(1/3), x]");
    check_simplify_sound("D[Sqrt[3] (x^4 + 1)^(1/4)/x + I (x^4 + 1)^(3/4)/x^3, x]");
    /* a rational radicand: h = (1 + x) Sqrt[(1 - x)/(1 + x)] is integral */
    check_simplify_sound("Sqrt[(1 - x)/(1 + x)] (1 + x)/(1 - x)");
    /* several variables: rationalise over Q(a)(x) */
    run_test("Simplify[1/(Sqrt[x] - Sqrt[y]) - (Sqrt[x] + Sqrt[y])/(x - y)]", "0");
    run_test("Simplify[D[Integrate[1/(x^3 (a + b x)^(1/3)), x], x]]", "1/(x^3*(a + b*x)^(1/3))");
    /* the radical in a denominator rationalised and re-absorbed */
    run_test("Simplify[2/(1 + x^3)^(1/3) - (2 (1 + x^3)^(2/3))/x^3]", "-2/(x^3*(1 + x^3)^(1/3))");
    run_test("Simplify[1/(1 + (8 + x^3)^(1/3)) + 1/(1 - (8 + x^3)^(1/3))]",
             "-2/(-1 + (8 + x^3)^(2/3))");
}

/* ---- 4. branch traps ---- */

static void test_branch_traps(void) {
    /* Sqrt[a] Sqrt[b] = Sqrt[a b] fails for a = b = -1 */
    run_test("Simplify[Sqrt[x] Sqrt[1 + x] - Sqrt[x + x^2]]", "Sqrt[x]*Sqrt[1 + x] - Sqrt[x + x^2]");
    /* Sqrt[-x] = I Sqrt[x] fails for x < 0 */
    run_test("Simplify[Sqrt[-x] - I Sqrt[x]]", "Sqrt[-x] - I*Sqrt[x]");
    /* (x^3)^(1/3) = x fails for x = -1 ((-1)^(1/3) is not -1) */
    run_test("Simplify[(x^3)^(1/3) - x]", "-x + (x^3)^(1/3)");
    run_test("Simplify[Sqrt[x^2] - x]", "-x + Sqrt[x^2]");
    /* (-B)^(1/3) = -B^(1/3) fails for B > 0 (principal cube root of a
     * negative number is not real) */
    run_test("Simplify[(-8 - x^3)^(1/3) + (8 + x^3)^(1/3)]", "(-8 - x^3)^(1/3) + (8 + x^3)^(1/3)");
    /* (1/x)^(1/2) = 1/Sqrt[x] fails on the negative axis */
    run_test("Simplify[Sqrt[1/x] - 1/Sqrt[x]]", "Sqrt[x^(-1)] - 1/Sqrt[x]");
    /* constants: exact principal values ((-2)^(1/3) = (-1)^(1/3) 2^(1/3) holds) */
    run_test("Simplify[(-2)^(1/3) - (-1)^(1/3) 2^(1/3)]", "0");
}

/* ---- 5. never larger ---- */

static void test_never_larger(void) {
    const char* cases[] = {
        "Sqrt[2]/(Sqrt[2] + Sqrt[1 + x^2])",
        "(x - 2^(1/3) y^(1/3))/(x^3 - 2 y)",
        "(a + b)/(Sqrt[a] + Sqrt[b])",
        "1/(1 + (8 + x^3)^(1/3)) + Sqrt[3] x/(2 + (8 + x^3)^(2/3))",
        "(1 + Sqrt[x])^3/(x^(1/3) - (-1)^(1/3))",
    };
    char buf[1024];
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        snprintf(buf, sizeof buf, "LeafCount[Simplify[%s]] <= LeafCount[%s]", cases[i], cases[i]);
        run_true(buf);
        check_simplify_sound(cases[i]);
    }
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    printf("Starting radical_field_simplify_tests\n");
    if (argc > 1) snprintf(repro_path, sizeof repro_path, "%s", argv[1]);
    symtab_init();
    core_init();

    TEST(test_repro);
    TEST(test_repro_timing);
    TEST(test_together_one_fraction);
    TEST(test_several_radicals);
    TEST(test_branch_traps);
    TEST(test_never_larger);

    printf("All radical_field_simplify tests passed!\n");
    return 0;
}

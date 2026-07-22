/* test_ratcanon_spec.c — the canonical output-form SPEC for the
 * Together/Cancel/Simplify rewrite (RATCANON_REWRITE_PLAN.md).
 *
 * This is the CONTRACT every later phase is tested against. Each `spec(...)` row
 * pins the WL-faithful canonical form (§0.4 of the plan) that the real
 * `Together`/`Cancel` builtins must produce — most already do today, so this
 * doubles as the regression baseline for Phases 3–4 (the switch + delete). Each
 * `proto_ok(...)` row proves the Phase-1 prototype's one-front-end + one-
 * reduction pipeline is a correct normalization. `spec_todo(...)` rows document
 * targets not yet met (soft — reported, not asserted; flip to `spec` when met).
 *
 * §0.4 recap: WL-faithful (radicals kept in the denominator, NOT rationalized);
 * top-level num/den expanded; kernel arguments normalized but not force-expanded;
 * one denominator-sign convention; Cancel leaves a sum of fractions uncombined,
 * Together combines.
 */
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int todo_pass = 0, todo_fail = 0;

/* Hard assertion: FullForm[evaluate(input)] == expected. */
static void spec(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Soft (documented) target: reported, not asserted. */
static void spec_todo(const char* input, const char* desired) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, desired) == 0) { todo_pass++; printf("  TODO now-PASSES: %s\n", input); }
    else { todo_fail++; printf("  TODO pending: %s\n    want: %s\n    got:  %s\n", input, desired, got); }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Prototype correctness: RatCanonPrototype[e] is math-equal to e, i.e.
 * Together[RatCanonPrototype[e] - (e)] collapses to 0. */
static void proto_ok(const char* inner) {
    char buf[512];
    snprintf(buf, sizeof buf, "Together[RatCanonPrototype[%s] - (%s)]", inner, inner);
    spec(buf, "0");
}

/* -------------------------------------------------------------------------- */

static void test_plain_q(void) {
    spec("Cancel[(1+x)/(1-x^2)]", "Times[-1, Power[Plus[-1, x], -1]]");
    spec("Cancel[(x^2-1)/(x-1)]", "Plus[1, x]");
    spec("Cancel[(x^3-1)/(x-1)]", "Plus[1, x, Power[x, 2]]");
    spec("Cancel[(x^4-1)/(x^2-1)]", "Plus[1, Power[x, 2]]");
    spec("Cancel[(x^2+2x+1)/(x+1)]", "Plus[1, x]");
    spec("Cancel[(2 x^2+4x)/(x^2-4)]", "Times[2, x, Power[Plus[-2, x], -1]]");
    spec("Cancel[(6x^2-6)/(3x-3)]", "Plus[2, Times[2, x]]");
    spec("Cancel[(a x^2 - a)/(x-1)]", "Plus[a, Times[a, x]]");
    spec("Together[1/x+1/y]", "Times[Power[x, -1], Power[y, -1], Plus[x, y]]");
    spec("Together[x/(x+1)+1/(x+1)]", "1");
    spec("Together[1/(x^2-1)+1/(x+1)]", "Times[x, Power[Plus[-1, Power[x, 2]], -1]]");
    spec("Together[2/(x-1)-2/(x+1)]", "Times[4, Power[Plus[-1, Power[x, 2]], -1]]");
    spec("Together[1/(a+b)+1/(a-b)]",
         "Times[2, a, Power[Plus[Power[a, 2], Times[-1, Power[b, 2]]], -1]]");
    spec("Together[b/(x-a)+b/(x+a)]",
         "Times[-2, b, x, Power[Plus[Power[a, 2], Times[-1, Power[x, 2]]], -1]]");
    spec("Together[a/(x-1)+b/(x-2)]",
         "Times[Power[Plus[2, Times[-3, x], Power[x, 2]], -1], "
         "Plus[Times[-2, a], Times[-1, b], Times[a, x], Times[b, x]]]");
}

static void test_transcendental(void) {
    spec("Together[1/(1+Log[x])+1/Log[x]]",
         "Times[Power[Plus[Log[x], Power[Log[x], 2]], -1], Plus[1, Times[2, Log[x]]]]");
    spec("Cancel[(Log[x]^2-1)/(Log[x]-1)]", "Plus[1, Log[x]]");
    spec("Together[1/Log[x]-1/Log[y]]",
         "Times[Power[Log[x], -1], Power[Log[y], -1], Plus[Times[-1, Log[x]], Log[y]]]");
    spec("Together[1/(1+ArcTan[x])+1/ArcTan[x]]",
         "Times[Power[Plus[ArcTan[x], Power[ArcTan[x], 2]], -1], Plus[1, Times[2, ArcTan[x]]]]");
    spec("Together[(E^(2x)-1)/(E^x-1)]", "Plus[1, Power[E, x]]");
    spec("Cancel[(E^(2x)-1)/(E^x+1)]", "Plus[-1, Power[E, x]]");
    spec("Together[1/(E^x-1)+1/(E^x+1)]",
         "Times[2, Power[E, x], Power[Plus[-1, Power[E, Times[2, x]]], -1]]");
    /* Tan is a transcendental generator reduced via FLINT (one consistent sign
     * convention: leading-coefficient-positive denominator), so this is
     * -2/(-1+Tan^2), math-equal to the classical path's 2/(1-Tan^2). */
    spec("Together[1/(1+Tan[x])+1/(1-Tan[x])]",
         "Times[-2, Power[Plus[-1, Power[Tan[x], 2]], -1]]");
    spec("Together[1/Sinh[x]+1/Cosh[x]]", "Plus[Csch[x], Sech[x]]");
}

static void test_gaussian(void) {
    spec("Together[1/(x-I)+1/(x+I)]", "Times[2, x, Power[Plus[1, Power[x, 2]], -1]]");
    spec("Together[1/(1-I x)+1/(1+I x)]", "Times[2, Power[Plus[1, Power[x, 2]], -1]]");
}

static void test_number_field(void) {
    /* Hard cancellation over Q(Sqrt d) — the field GCD must see (x-Sqrt d). */
    spec("Cancel[(x^2-2)/(x-Sqrt[2])]", "Plus[Power[2, Rational[1, 2]], x]");
    spec("Cancel[(x^2-3)/(x+Sqrt[3])]", "Plus[Times[-1, Power[3, Rational[1, 2]]], x]");
    spec("Cancel[(x-Sqrt[2])(x+Sqrt[2])/(x-Sqrt[2])]", "Plus[Power[2, Rational[1, 2]], x]");
    /* WL-faithful: radical KEPT in the denominator, not rationalized. */
    spec("Cancel[1/(x-Sqrt[2])]", "Power[Plus[Times[-1, Power[2, Rational[1, 2]]], x], -1]");
    /* Conjugate combine clears the radical naturally. */
    spec("Together[1/(x-Sqrt[2])+1/(x+Sqrt[2])]",
         "Times[2, x, Power[Plus[-2, Power[x, 2]], -1]]");
    /* Two distinct radicals: Q(Sqrt3, Sqrt5) tower kept. */
    spec("Together[1/(x-Sqrt[3])+1/(x-Sqrt[5])]",
         "Times[Plus[Times[-1, Power[3, Rational[1, 2]]], Times[-1, Power[5, Rational[1, 2]]], "
         "Times[2, x]], Power[Plus[Power[15, Rational[1, 2]], "
         "Times[-1, Times[Power[3, Rational[1, 2]], x]], "
         "Times[-1, Times[Power[5, Rational[1, 2]], x]], Power[x, 2]], -1]]");
    /* WL keeps a coprime constant fraction (no rationalization for Cancel). */
    spec("Cancel[(1+Sqrt[2])/(1-Sqrt[2])]",
         "Times[Plus[1, Power[2, Rational[1, 2]]], "
         "Power[Plus[1, Times[-1, Power[2, Rational[1, 2]]]], -1]]");
}

static void test_cyclotomic(void) {
    spec("Cancel[x/(x-(-1)^(1/3))]",
         "Times[x, Power[Plus[Times[-1, Power[-1, Rational[1, 3]]], x], -1]]");
    spec("Together[1/(x-(-1)^(1/3))+1/(x+(-1)^(1/3))]",
         "Times[2, x, Power[Plus[Times[-1, Power[-1, Rational[2, 3]]], Power[x, 2]], -1]]");
}

static void test_symbolic_radical(void) {
    /* Q(params)(Sqrt k), k a free symbol — the Sqrt[k] regime. */
    spec("Together[1/(b-a Sqrt[k])+1/(b+a Sqrt[k])]",
         "Times[-2, b, Power[Plus[Times[-1, Power[b, 2]], Times[Power[a, 2], k]], -1]]");
    spec("Cancel[(a^2-b)/(a-Sqrt[b])]", "Plus[a, Power[b, Rational[1, 2]]]");
    /* Cube-root relation cancellation. */
    spec("Cancel[(y-1)/(y^(1/3)-1)]",
         "Plus[1, Power[y, Rational[1, 3]], Power[y, Rational[2, 3]]]");
}

static void test_cancel_vs_together(void) {
    /* Cancel leaves a sum of fractions uncombined; Together combines. */
    spec("Cancel[a/b+c/d]",
         "Plus[Times[a, Power[b, -1]], Times[c, Power[d, -1]]]");
    spec("Together[a/b+c/d]",
         "Times[Power[b, -1], Power[d, -1], Plus[Times[b, c], Times[a, d]]]");
}

static void test_prototype_correct(void) {
    /* Phase-1: one substitution front-end + one FLINT reduction is correct
     * across all four representative regimes. */
    proto_ok("(1+x)/(1-x^2)");                       /* plain Q            */
    proto_ok("1/(1+Log[x])+1/Log[x]");               /* transcendental     */
    proto_ok("1/(x-I)+1/(x+I)");                     /* Q(i)               */
    proto_ok("1/(b-a Sqrt[k])+1/(b+a Sqrt[k])");     /* symbolic Sqrt[k]   */
    proto_ok("(E^(2x)-1)/(E^x-1)");                  /* exp commensurate   */
    proto_ok("1/(x-Sqrt[2])+1/(x+Sqrt[2])");         /* Q(Sqrt d)          */
}

static void test_higher_radicands(void) {
    /* Q(Sqrt d), d >= 5: must NOT crash (TOGETHER_ALGEBRAIC_OVERFLOW.md) and must
     * give the WL-faithful combined form. Currently passes — pinned as a guard. */
    spec("Together[1/(x^2-5)+1/(x-Sqrt[5])]",
         "Times[Power[Plus[-5, Power[x, 2]], -1], Plus[1, Power[5, Rational[1, 2]], x]]");
    spec("Together[1/(x^2-7)+1/(x-Sqrt[7])]",
         "Times[Power[Plus[-7, Power[x, 2]], -1], Plus[1, Power[7, Rational[1, 2]], x]]");
}

int main(void) {
    symtab_init();
    core_init();
    (void)spec_todo;   /* retained for future-phase soft targets */

    printf("Running ratcanon SPEC tests...\n");
    TEST(test_plain_q);
    TEST(test_transcendental);
    TEST(test_gaussian);
    TEST(test_number_field);
    TEST(test_cyclotomic);
    TEST(test_symbolic_radical);
    TEST(test_cancel_vs_together);
    TEST(test_prototype_correct);
    TEST(test_higher_radicands);
    printf("All ratcanon SPEC assertions passed. (TODO: %d passing, %d pending)\n",
           todo_pass, todo_fail);
    return 0;
}

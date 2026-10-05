/* test_integrate_intrep.c
 *
 * Tests for definite integration by special-function integral representation
 * (src/calculus/integrate_intrep.c) and the Euler -> Beta*2F1 extension
 * (src/calculus/integrate_beta.c: integrate_euler_2f1_try):
 *   E^(-p x) BesselJ[nu, q x]          -> Laplace-Bessel
 *   E^(-A Cosh[x]) Cosh[n x]           -> BesselK[n, A]
 *   x^(nu-1) E^(-A x - B/x)            -> 2 (B/A)^(nu/2) BesselK[nu, 2 Sqrt[A B]]
 *   Cos[p x^3 + q x]                   -> Pi (3p)^(-1/3) AiryAi[q (3p)^(-1/3)]
 *   x^(a-1)(1-x)^(b-1)(alpha+beta x)^e -> Euler Beta*2F1 on [0,1]
 *
 * Each closing case is SELF-CERTIFYING: it asserts Simplify[result - reference]
 * === 0, so a decline (which leaves an Integrate head) can never pass -- the
 * difference would not reduce to 0.  The negative controls pin that an
 * out-of-domain integrand stays unevaluated.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Self-certifying numeric check: assert {FreeQ[input, Integrate], Chop[N[(input -
 * reference) /. sample]]} === {True, 0}.  Numeric (not Simplify) because the closed
 * forms are equivalent across representations Simplify will not reconcile to 0
 * (exp vs BesselK, symbolic-exponent Beta forms); the FreeQ guard makes a decline
 * (which leaves an Integrate head) FAIL rather than pass.  `sample` substitutes the
 * free parameters to generic rationals. */
static void check_num(const char* input, const char* reference, const char* sample) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
             "{FreeQ[%s, Integrate], Chop[N[((%s) - (%s)) /. %s], 10^-7]}",
             input, input, reference, sample);
    Expr* p = parse_expression(buf);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, "{True, 0}") != 0)
        fprintf(stderr, "FAIL: %s\n  {FreeQ, Chop[N[diff]]} = %s\n", input, s);
    ASSERT_STR_EQ(s, "{True, 0}");
    free(s);
    expr_free(p);
    expr_free(r);
}

/* Assert the integrand stays unevaluated (result still carries an Integrate). */
static void check_unevaluated(const char* input) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "FreeQ[%s, Integrate]", input);
    Expr* p = parse_expression(buf);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, "False") != 0)
        fprintf(stderr, "FAIL (should stay unevaluated): %s\n  got FreeQ: %s\n", input, s);
    ASSERT_STR_EQ(s, "False");
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------- Laplace-Bessel: E^(-p x) BesselJ[nu, q x] ---------------------- */
static void test_laplace_bessel(void) {
    /* #2: nu=0 -> 1/Sqrt[a^2+c^2]. */
    check_num(
        "Integrate[Exp[-c x] BesselJ[0, a x], {x,0,Infinity}, Assumptions -> c>0 && Element[a,Reals]]",
        "1/Sqrt[a^2 + c^2]", "{a -> 11/10, c -> 23/10}");
    /* nu=1, q>0: (Sqrt[q^2+p^2]-p)/(q Sqrt[q^2+p^2]). */
    check_num(
        "Integrate[Exp[-c x] BesselJ[1, a x], {x,0,Infinity}, Assumptions -> c>0 && a>0]",
        "(Sqrt[a^2+c^2]-c)/(a Sqrt[a^2+c^2])", "{a -> 11/10, c -> 23/10}");
    /* pinned method, concrete. */
    check_num(
        "Integrate[Exp[-2 x] BesselJ[0, 3 x], {x,0,Infinity}, Method -> \"IntegralRepresentation\"]",
        "1/Sqrt[13]", "{}");
}

/* -------- Bessel-K (cosh): E^(-A Cosh[x]) Cosh[n x] ---------------------- */
static void test_besselk_cosh(void) {
    /* #15: -> BesselK[n, a]. */
    check_num(
        "Integrate[Exp[-a Cosh[x]] Cosh[n x], {x,0,Infinity}, Assumptions -> a>0 && Element[n,Reals]]",
        "BesselK[n, a]", "{a -> 3/2, n -> 7/10}");
    /* n=0 (no outer Cosh factor) -> BesselK[0, a]. */
    check_num(
        "Integrate[Exp[-a Cosh[x]], {x,0,Infinity}, Assumptions -> a>0]",
        "BesselK[0, a]", "{a -> 3/2}");
}

/* -------- Bessel-K (exp): x^(nu-1) E^(-A x - B/x) ------------------------ */
static void test_besselk_exp(void) {
    /* #19: nu=1/2 -> elementary Sqrt[Pi/a] e^(-2 Sqrt[a b]). */
    check_num(
        "Integrate[Exp[-a x - b/x]/Sqrt[x], {x,0,Infinity}, Assumptions -> a>0 && b>0]",
        "Sqrt[Pi/a] Exp[-2 Sqrt[a b]]", "{a -> 13/10, b -> 21/10}");
    /* nu=1 (bare exp): 2 (b/a)^(1/2) BesselK[1, 2 Sqrt[a b]]. */
    check_num(
        "Integrate[Exp[-a x - b/x], {x,0,Infinity}, Assumptions -> a>0 && b>0]",
        "2 Sqrt[b/a] BesselK[1, 2 Sqrt[a b]]", "{a -> 13/10, b -> 21/10}");
}

/* -------- Airy: Cos[p x^3 + q x] ----------------------------------------- */
static void test_airy(void) {
    /* #12: Cos[x^3/3 + a x] -> Pi AiryAi[a]. */
    check_num(
        "Integrate[Cos[x^3/3 + a x], {x,0,Infinity}, Assumptions -> Element[a,Reals]]",
        "Pi AiryAi[a]", "{a -> 7/10}");
    /* scaled cubic Cos[x^3] = Cos[1 x^3 + 0 x] -> Pi 3^(-1/3) AiryAi[0]. */
    check_num(
        "Integrate[Cos[x^3], {x,0,Infinity}]",
        "Pi 3^(-1/3) AiryAi[0]", "{}");
}

/* -------- Euler finite-interval -> Beta * 2F1 ---------------------------- */
static void test_euler_2f1(void) {
    /* #18: -> Beta[a,b] c^(-b) (1+c)^(-a). */
    check_num(
        "Integrate[(x^(a-1)(1-x)^(b-1))/(x+c)^(a+b), {x,0,1}, Assumptions -> a>0 && b>0 && c>0]",
        "Beta[a,b] c^(-b) (1+c)^(-a)", "{a -> 2, b -> 3, c -> 4}");
    /* concrete a=2,b=3,c=4 -> 1/19200. */
    check_num(
        "Integrate[(x^(2-1)(1-x)^(3-1))/(x+4)^(2+3), {x,0,1}]",
        "1/19200", "{}");
    /* general (1 - z x)^(-s) -> Beta 2F1. */
    check_num(
        "Integrate[x^(a-1)(1-x)^(b-1)(1 - z x)^(-s), {x,0,1}, Assumptions -> a>0 && b>0 && 0<z<1 && s>0]",
        "Beta[a,b] Hypergeometric2F1[s, a, a+b, z]", "{a -> 2, b -> 3, z -> 2/5, s -> 6/5}");
}

/* -------- Lerch / Hurwitz: x^(s-1) E^(-a x)/(1 - z E^(-c x)) ------------- */
static void test_lerch_hurwitz(void) {
    /* #6: z=1 -> Gamma[s] HurwitzZeta[s, a]. */
    check_num(
        "Integrate[Exp[-a x] x^(s-1)/(1 - Exp[-x]), {x,0,Infinity}, Assumptions -> s>1 && a>0]",
        "Gamma[s] HurwitzZeta[s, a]", "{s -> 5/2, a -> 7/5}");
    /* scaled c=2 -> 2^(-s) Gamma[s] HurwitzZeta[s, a/2] (the separate A>0, c>0 gate). */
    check_num(
        "Integrate[x^(s-1) Exp[-a x]/(1 - Exp[-2 x]), {x,0,Infinity}, Assumptions -> s>1 && a>0]",
        "2^(-s) Gamma[s] HurwitzZeta[s, a/2]", "{s -> 5/2, a -> 7/5}");
    /* fugacity z=1/2 -> Gamma[s] LerchPhi[1/2, s, a]. */
    check_num(
        "Integrate[x^(s-1) Exp[-a x]/(1 - (1/2) Exp[-x]), {x,0,Infinity}, Assumptions -> s>0 && a>0]",
        "Gamma[s] LerchPhi[1/2, s, a]", "{s -> 5/2, a -> 7/5}");
    /* z=-1 (eta-like) -> Gamma[s] LerchPhi[-1, s, a]. */
    check_num(
        "Integrate[x^(s-1) Exp[-a x]/(1 + Exp[-x]), {x,0,Infinity}, Assumptions -> s>0 && a>0]",
        "Gamma[s] LerchPhi[-1, s, a]", "{s -> 5/2, a -> 7/5}");
    /* pinned method, concrete s=2,a=1 -> HurwitzZeta[2,1] = Pi^2/6. */
    check_num(
        "Integrate[x Exp[-x]/(1 - Exp[-x]), {x,0,Infinity}, Method -> \"IntegralRepresentation\"]",
        "Pi^2/6", "{}");
}

/* -------- Negative controls: must stay unevaluated ----------------------- */
static void test_declines(void) {
    /* Divergent: no decay factor (a=0), so HurwitzZeta[s,0] would be a pole. */
    check_unevaluated("Integrate[x^(s-1)/(1 - Exp[-x]), {x,0,Infinity}, Method -> \"IntegralRepresentation\", Assumptions -> s>1]");
    /* s=1 lands on the Hurwitz-zeta pole: must decline (Re s > 1 gate). */
    check_unevaluated("Integrate[Exp[-a x]/(1 - Exp[-x]), {x,0,Infinity}, Method -> \"IntegralRepresentation\", Assumptions -> a>0]");
    /* No assumptions: convergence gate (c>0) cannot be proved. */
    check_unevaluated("Integrate[Exp[-c x] BesselJ[0, a x], {x,0,Infinity}]");
    /* Wrong sign (growth, not decay). */
    check_unevaluated("Integrate[Exp[a x - b/x]/Sqrt[x], {x,0,Infinity}, Assumptions -> a>0 && b>0]");
    /* Not a recognised shape. */
    check_unevaluated("Integrate[Exp[-c x] Sin[a x]^3, {x,0,Infinity}, Method -> \"IntegralRepresentation\"]");
    /* Euler with the extra factor's base not positive on [0,1] (alpha<0). */
    check_unevaluated("Integrate[(x^(a-1)(1-x)^(b-1))/(x-2)^(a+b), {x,0,1}, Method -> \"Beta\", Assumptions -> a>0 && b>0]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_laplace_bessel);
    TEST(test_besselk_cosh);
    TEST(test_besselk_exp);
    TEST(test_airy);
    TEST(test_euler_2f1);
    TEST(test_lerch_hurwitz);
    TEST(test_declines);

    printf("All integrate_intrep tests passed.\n");
    return 0;
}

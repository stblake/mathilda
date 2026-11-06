/* Tests for NRoots — numerical roots of a univariate polynomial equation.
 *
 * Cover: real and complex roots; complex coefficients; multiplicity (identical
 * equations); exact zero roots (trailing x^m factor); canonical ordering;
 * arbitrary precision via PrecisionGoal; degenerate equations (True/False);
 * the explicit Method settings; unevaluated forms; Protected attribute; and
 * memory hygiene.
 *
 * Roots are verified *inside* the language: either compared to a known value
 * (N[Abs[root - expected]] < tol) or substituted back into the polynomial
 * (residual ~ 0).  The disjunction  x==r1 || ... || x==rd  is indexed with
 * Part[res, k, 2] to pull out the k-th root.
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

/* |Part[call, idx, 2] - expected| < tol */
static bool root_close(const char* call, int idx, const char* expected, double tol) {
    char buf[2048];
    snprintf(buf, sizeof buf,
             "N[Abs[Re[Part[%s,%d,2]-(%s)]]+Abs[Im[Part[%s,%d,2]-(%s)]]] < %.17g",
             call, idx, expected, call, idx, expected, tol);
    bool ok = lang_true(buf);
    if (!ok) fprintf(stderr, "  root_close FAIL: %s[[%d,2]] ~= %s\n", call, idx, expected);
    return ok;
}

#define ASSERT_ROOT(call, idx, exp, tol) \
    ASSERT_MSG(root_close((call), (idx), (exp), (tol)), \
               "%s[[%d,2]] ~= %s (tol %g)", (call), (idx), (exp), (tol))

#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* ---------------------------------------------------------------------- */

static void test_real_roots(void) {
    /* x^2 - 2 == 0 : -Sqrt2 then Sqrt2 (reals ascending). */
    ASSERT_TRUE("Length[NRoots[x^2-2==0,x]] == 2");
    ASSERT_ROOT("NRoots[x^2-2==0,x]", 1, "-Sqrt[2]", 1e-10);
    ASSERT_ROOT("NRoots[x^2-2==0,x]", 2, "Sqrt[2]",  1e-10);

    /* (x-1)(x-2)(x-3) = x^3 - 6 x^2 + 11 x - 6. */
    ASSERT_TRUE("Length[NRoots[x^3-6x^2+11x-6==0,x]] == 3");
    ASSERT_ROOT("NRoots[x^3-6x^2+11x-6==0,x]", 1, "1", 1e-9);
    ASSERT_ROOT("NRoots[x^3-6x^2+11x-6==0,x]", 2, "2", 1e-9);
    ASSERT_ROOT("NRoots[x^3-6x^2+11x-6==0,x]", 3, "3", 1e-9);

    /* Linear : single bare equation (not a disjunction). */
    ASSERT_TRUE("N[Abs[Part[NRoots[2x-6==0,x],2] - 3]] < 1e-12");
}

static void test_cubic_example(void) {
    /* The documentation example: 1 + 2x + 3x^2 + 4x^3 == 0. */
    const char* C = "NRoots[1+2x+3x^2+4x^3==0,x]";
    ASSERT_TRUE("Length[NRoots[1+2x+3x^2+4x^3==0,x]] == 3");
    ASSERT_ROOT(C, 1, "-0.605829", 1e-4);
    ASSERT_ROOT(C, 2, "-0.0720854 - 0.638327 I", 1e-4);
    ASSERT_ROOT(C, 3, "-0.0720854 + 0.638327 I", 1e-4);
    /* Residual: every root substituted back gives ~0. */
    ASSERT_TRUE("N[Abs[1+2 r+3 r^2+4 r^3] /. r->Part[NRoots[1+2x+3x^2+4x^3==0,x],1,2]] < 1e-9");
    ASSERT_TRUE("N[Abs[1+2 r+3 r^2+4 r^3] /. r->Part[NRoots[1+2x+3x^2+4x^3==0,x],2,2]] < 1e-9");
}

static void test_complex_roots(void) {
    /* x^2 + 1 : -I then I. */
    ASSERT_ROOT("NRoots[x^2+1==0,x]", 1, "-I", 1e-10);
    ASSERT_ROOT("NRoots[x^2+1==0,x]", 2, "I",  1e-10);

    /* x^4 - 1 : -1, 1, -I, I (reals first). */
    ASSERT_TRUE("Length[NRoots[x^4-1==0,x]] == 4");
    ASSERT_ROOT("NRoots[x^4-1==0,x]", 1, "-1", 1e-10);
    ASSERT_ROOT("NRoots[x^4-1==0,x]", 2, "1",  1e-10);
    ASSERT_ROOT("NRoots[x^4-1==0,x]", 3, "-I", 1e-10);
    ASSERT_ROOT("NRoots[x^4-1==0,x]", 4, "I",  1e-10);

    /* x^3 - 1 : 1, then the primitive cube roots. */
    ASSERT_ROOT("NRoots[x^3-1==0,x]", 1, "1", 1e-10);
    ASSERT_ROOT("NRoots[x^3-1==0,x]", 2, "-1/2 - Sqrt[3]/2 I", 1e-9);
    ASSERT_ROOT("NRoots[x^3-1==0,x]", 3, "-1/2 + Sqrt[3]/2 I", 1e-9);
}

static void test_complex_coefficients(void) {
    /* x^2 = 3 + 4 I  ->  +-(2 + I). */
    ASSERT_TRUE("Length[NRoots[x^2-(3+4 I)==0,x]] == 2");
    ASSERT_ROOT("NRoots[x^2-(3+4 I)==0,x]", 1, "-2 - I", 1e-9);
    ASSERT_ROOT("NRoots[x^2-(3+4 I)==0,x]", 2, "2 + I",  1e-9);

    /* (x - (1+I)) (x - (2-3I)) = x^2 - (3-2I) x + (5 - I). */
    ASSERT_TRUE("N[Abs[Part[NRoots[x^2-(3-2I)x+(5-I)==0,x],1,2] - (1+I)]] < 1e-8 || "
                "N[Abs[Part[NRoots[x^2-(3-2I)x+(5-I)==0,x],2,2] - (1+I)]] < 1e-8");
}

static void test_multiplicity(void) {
    /* (x-1)^3 : three identical equations. */
    ASSERT_TRUE("Length[NRoots[(x-1)^3==0,x]] == 3");
    ASSERT_ROOT("NRoots[(x-1)^3==0,x]", 1, "1", 1e-4);
    ASSERT_ROOT("NRoots[(x-1)^3==0,x]", 3, "1", 1e-4);
    /* Identical equations print identically. */
    ASSERT_TRUE("Part[NRoots[(x-1)^3==0,x],1] == Part[NRoots[(x-1)^3==0,x],3]");

    /* x^2 - 2x + 1 = (x-1)^2. */
    ASSERT_TRUE("Length[NRoots[x^2-2x+1==0,x]] == 2");
    ASSERT_ROOT("NRoots[x^2-2x+1==0,x]", 1, "1", 1e-6);

    /* (x^2+1)^2 : -I, -I, I, I. */
    ASSERT_TRUE("Length[NRoots[(x^2+1)^2==0,x]] == 4");
    ASSERT_ROOT("NRoots[(x^2+1)^2==0,x]", 1, "-I", 1e-4);
    ASSERT_ROOT("NRoots[(x^2+1)^2==0,x]", 4, "I",  1e-4);
}

static void test_zero_roots(void) {
    /* x^3 : triple root at 0. */
    ASSERT_TRUE("Length[NRoots[x^3==0,x]] == 3");
    ASSERT_ROOT("NRoots[x^3==0,x]", 1, "0", 1e-12);
    ASSERT_ROOT("NRoots[x^3==0,x]", 3, "0", 1e-12);

    /* x^3 - x = x(x-1)(x+1). */
    ASSERT_TRUE("Length[NRoots[x^3-x==0,x]] == 3");
    ASSERT_ROOT("NRoots[x^3-x==0,x]", 1, "-1", 1e-10);
    ASSERT_ROOT("NRoots[x^3-x==0,x]", 2, "0",  1e-10);
    ASSERT_ROOT("NRoots[x^3-x==0,x]", 3, "1",  1e-10);
}

static void test_precision(void) {
    /* PrecisionGoal -> 30 : 30-digit Sqrt[2]. */
    ASSERT_TRUE("N[Abs[Part[NRoots[x^2-2==0,x,PrecisionGoal->30],2,2] - Sqrt[2]], 40] < 10^-29");
    /* High-precision cube root of 2. */
    ASSERT_TRUE("N[Abs[Part[NRoots[x^3-2==0,x,PrecisionGoal->40],1,2] - 2^(1/3)], 50] < 10^-38");
}

static void test_degenerate(void) {
    ASSERT_STR_EQ(eval_str("NRoots[1==0,x]"), "False");
    ASSERT_STR_EQ(eval_str("NRoots[1==1,x]"), "True");
    ASSERT_STR_EQ(eval_str("NRoots[2==2,x]"), "True");
}

static void test_method_aberth(void) {
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"Aberth\"]", 1, "-Sqrt[2]", 1e-10);
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"Aberth\"]", 2, "Sqrt[2]",  1e-10);
    ASSERT_ROOT("NRoots[x^3-1==0,x,Method->\"Aberth\"]", 2, "-1/2 - Sqrt[3]/2 I", 1e-9);
    /* MaxIterations is accepted. */
    ASSERT_ROOT("NRoots[x^2-2==0,x,MaxIterations->200]", 2, "Sqrt[2]", 1e-10);
}

static void test_method_companion(void) {
    /* Real coefficients: direct real QR. */
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"CompanionMatrix\"]", 1, "-Sqrt[2]", 1e-10);
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"CompanionMatrix\"]", 2, "Sqrt[2]",  1e-10);
    ASSERT_TRUE("Length[NRoots[1+2x+3x^2+4x^3==0,x,Method->\"CompanionMatrix\"]] == 3");
    ASSERT_ROOT("NRoots[x^3-6x^2+11x-6==0,x,Method->\"CompanionMatrix\"]", 2, "2", 1e-9);
    /* Multiplicity. */
    ASSERT_TRUE("Length[NRoots[(x-1)^3==0,x,Method->\"CompanionMatrix\"]] == 3");
    ASSERT_ROOT("NRoots[(x-1)^3==0,x,Method->\"CompanionMatrix\"]", 1, "1", 1e-4);

    /* Complex coefficients via the 2n real embedding. */
    ASSERT_ROOT("NRoots[x^2-(3+4 I)==0,x,Method->\"CompanionMatrix\"]", 1, "-2 - I", 1e-9);
    ASSERT_ROOT("NRoots[x^2-(3+4 I)==0,x,Method->\"CompanionMatrix\"]", 2, "2 + I",  1e-9);
    /* Complex polynomial with a real root: 1 must appear once, not twice. */
    ASSERT_TRUE("Length[NRoots[x^2-(1+I)x+I==0,x,Method->\"CompanionMatrix\"]] == 2");
    ASSERT_TRUE("(N[Abs[Part[NRoots[x^2-(1+I)x+I==0,x,Method->\"CompanionMatrix\"],1,2]-1]]<1e-8) || "
                "(N[Abs[Part[NRoots[x^2-(1+I)x+I==0,x,Method->\"CompanionMatrix\"],2,2]-1]]<1e-8)");
    /* Complex double root + simple root. */
    ASSERT_TRUE("Length[NRoots[(x-(1+I))^2 (x-2)==0,x,Method->\"CompanionMatrix\"]] == 3");
}

static void test_method_jenkinstraub(void) {
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"JenkinsTraub\"]", 1, "-Sqrt[2]", 1e-10);
    ASSERT_ROOT("NRoots[x^2-2==0,x,Method->\"JenkinsTraub\"]", 2, "Sqrt[2]",  1e-10);
    ASSERT_TRUE("Length[NRoots[1+2x+3x^2+4x^3==0,x,Method->\"JenkinsTraub\"]] == 3");
    ASSERT_ROOT("NRoots[1+2x+3x^2+4x^3==0,x,Method->\"JenkinsTraub\"]", 1, "-0.605829", 1e-4);
    ASSERT_ROOT("NRoots[x^3-6x^2+11x-6==0,x,Method->\"JenkinsTraub\"]", 3, "3", 1e-9);
    /* Complex roots and coefficients. */
    ASSERT_ROOT("NRoots[x^3-1==0,x,Method->\"JenkinsTraub\"]", 2, "-1/2 - Sqrt[3]/2 I", 1e-9);
    ASSERT_ROOT("NRoots[x^2-(3+4 I)==0,x,Method->\"JenkinsTraub\"]", 2, "2 + I", 1e-9);
    /* Multiplicity (found m times by deflation). */
    ASSERT_TRUE("Length[NRoots[(x-1)^3==0,x,Method->\"JenkinsTraub\"]] == 3");
    ASSERT_ROOT("NRoots[(x-1)^3==0,x,Method->\"JenkinsTraub\"]", 1, "1", 1e-4);
    /* High precision. */
    ASSERT_TRUE("N[Abs[Part[NRoots[x^2-2==0,x,Method->\"JenkinsTraub\",PrecisionGoal->30],2,2] - Sqrt[2]], 40] < 10^-29");
}

/* ---------------------------------------------------------------------- *
 *  Cross-method / cross-precision stress battery.
 *
 *  Each polynomial is solved with every Method at both machine and 25-digit
 *  precision, and the canonically-ordered roots are compared element-by-element
 *  against the known exact set.  Covers real and complex coefficients, real and
 *  complex roots, and multiplicity.
 * ---------------------------------------------------------------------- */
/* Solve once into a scratch symbol, then check the count and every root against
 * the known set element-by-element (avoids With/Table over a held body). */
static void check_roots(const char* poly, const char* expected, int deg,
                        const char* method, const char* prec, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf, "nrtmp = NRoots[(%s) == 0, x%s%s]", poly, method, prec);
    free(eval_str(buf));

    snprintf(buf, sizeof buf, "Length[nrtmp] == %d", deg);
    char* s = eval_str(buf);
    ASSERT_MSG(strcmp(s, "True") == 0,
               "NRoots[(%s)==0, x%s%s] count != %d", poly, method, prec, deg);
    free(s);

    for (int k = 1; k <= deg; k++) {
        snprintf(buf, sizeof buf,
            "N[Abs[Re[Part[nrtmp, %d, 2] - Part[%s, %d]]] + "
            "Abs[Im[Part[nrtmp, %d, 2] - Part[%s, %d]]], 45] < %s",
            k, expected, k, k, expected, k, tol);
        s = eval_str(buf);
        ASSERT_MSG(strcmp(s, "True") == 0,
                   "NRoots[(%s)==0, x%s%s] root %d mismatch (=> %s)",
                   poly, method, prec, k, s);
        free(s);
    }
}

static void test_stress_all_methods(void) {
    struct { const char* poly; const char* exp; int deg; } cases[] = {
        /* real coefficients, real roots */
        { "x^2-2", "{-Sqrt[2], Sqrt[2]}", 2 },
        { "x^3-6x^2+11x-6", "{1, 2, 3}", 3 },
        { "x^2-2x+1", "{1, 1}", 2 },
        /* real coefficients, complex roots */
        { "x^2+1", "{-I, I}", 2 },
        { "x^4-1", "{-1, 1, -I, I}", 4 },
        { "x^3-1", "{1, -1/2 - Sqrt[3]/2 I, -1/2 + Sqrt[3]/2 I}", 3 },
        { "x^3-2", "{2^(1/3), -2^(1/3)/2 - 2^(1/3) Sqrt[3]/2 I, -2^(1/3)/2 + 2^(1/3) Sqrt[3]/2 I}", 3 },
        { "x^4+1", "{-1/Sqrt[2] - I/Sqrt[2], -1/Sqrt[2] + I/Sqrt[2], 1/Sqrt[2] - I/Sqrt[2], 1/Sqrt[2] + I/Sqrt[2]}", 4 },
        { "x^6-1", "{-1, 1, -1/2 - Sqrt[3]/2 I, -1/2 + Sqrt[3]/2 I, 1/2 - Sqrt[3]/2 I, 1/2 + Sqrt[3]/2 I}", 6 },
        /* real coefficients, multiplicity (exact squarefree path) */
        { "(x-1)^2 (x-3)", "{1, 1, 3}", 3 },
        { "(x^2-2)^3", "{-Sqrt[2], -Sqrt[2], -Sqrt[2], Sqrt[2], Sqrt[2], Sqrt[2]}", 6 },
        { "(x-1)^4", "{1, 1, 1, 1}", 4 },
        /* complex coefficients, distinct roots */
        { "x^2-(3+4 I)", "{-2 - I, 2 + I}", 2 },
        { "x^3-I", "{-Sqrt[3]/2 + I/2, -I, Sqrt[3]/2 + I/2}", 3 },
        { "(x-(1+I)) (x-(2-I))", "{1 + I, 2 - I}", 2 },
    };
    const char* methods[] = {
        "", ", Method->\"Aberth\"", ", Method->\"CompanionMatrix\"", ", Method->\"JenkinsTraub\""
    };
    struct { const char* opt; const char* tol; } precs[] = {
        { "", "1*10^-6" },
        { ", PrecisionGoal->25", "1*10^-20" },
    };
    int nc = (int)(sizeof(cases)/sizeof(cases[0]));
    for (int i = 0; i < nc; i++)
        for (int m = 0; m < 4; m++)
            for (int pp = 0; pp < 2; pp++)
                check_roots(cases[i].poly, cases[i].exp, cases[i].deg,
                            methods[m], precs[pp].opt, precs[pp].tol);
}

static void test_stress_high_multiplicity(void) {
    /* The reported bug: a high power must stay well-conditioned (squarefree
     * path), giving exactly the repeated roots at full precision and fast.
     * (x^2-2)^30 ⇒ 60 roots, each within 1e-25 of ±Sqrt[2] at 30 digits. */
    free(eval_str("nrbig = NRoots[(x^2-2)^30 == 0, x, PrecisionGoal -> 30]"));
    ASSERT_TRUE("Length[nrbig] == 60");
    for (int k = 1; k <= 60; k++) {
        char buf[512];
        snprintf(buf, sizeof buf,
            "N[Min[Abs[Part[nrbig,%d,2] - Sqrt[2]], Abs[Part[nrbig,%d,2] + Sqrt[2]]], 45] < 10^-25",
            k, k);
        char* s = eval_str(buf);
        ASSERT_MSG(strcmp(s, "True") == 0, "(x^2-2)^30 root %d off ±Sqrt[2] (=> %s)", k, s);
        free(s);
    }
    /* High multiplicity at machine precision via the squarefree path. */
    check_roots("(x-2)^20", "Table[2, {20}]", 20, ", Method->\"JenkinsTraub\"", "", "1*10^-8");
    check_roots("(x^2+1)^10", "Join[Table[-I, {10}], Table[I, {10}]]", 20, ", Method->\"CompanionMatrix\"", "", "1*10^-6");
}

static void test_unevaluated(void) {
    /* Non-polynomial in x -> unevaluated. */
    ASSERT_STR_EQ(eval_str("Head[NRoots[Sin[x]==0,x]]"), "NRoots");
    /* Second argument not a symbol -> unevaluated. */
    ASSERT_STR_EQ(eval_str("Head[NRoots[x^2-2==0,5]]"), "NRoots");
    /* Unknown method -> unevaluated. */
    ASSERT_STR_EQ(eval_str("Head[NRoots[x^2-2==0,x,Method->\"Bogus\"]]"), "NRoots");
}

static void test_attributes(void) {
    ASSERT_TRUE("MemberQ[Attributes[NRoots], Protected]");
}

static void test_memory_loop(void) {
    const char* inputs[] = {
        "NRoots[1+2x+3x^2+4x^3==0,x]",
        "NRoots[x^4-1==0,x]",
        "NRoots[(x-1)^3==0,x]",
        "NRoots[x^2-(3+4 I)==0,x]",
        "NRoots[x^3-x==0,x]",
        "NRoots[x^2-2==0,x,PrecisionGoal->30]",
    };
    for (int i = 0; i < 3; i++)
        for (size_t k = 0; k < sizeof(inputs)/sizeof(inputs[0]); k++) {
            Expr* p = parse_expression(inputs[k]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
}

/* ---------------------------------------------------------------------- *
 *  High degree with a cancelling expansion (the adaptive-precision engine)
 * ---------------------------------------------------------------------- */

/* The monomial expansion of LaguerreL[1000, x] has |c_k x^k| up to ~1e1118
 * where |L_1000| ~ 1e856 (and ~1e470 digits of cancellation mid-spectrum):
 * every root needs its own evaluation precision.  Truth: the largest zero
 * 3943.24739484527095238972810775344564096526363934564746059856857487
 * (computed independently by Newton on the three-term recurrence at 100
 * digits), and the zero count 1000 (L_n has n simple positive zeros). */
static void test_laguerre_1000_wp50(void) {
    ASSERT_TRUE("Module[{r = x /. NSolve[LaguerreL[1000, x] == 0, x, WorkingPrecision -> 50]}, "
                "Length[r] == 1000 && Count[r, _Real] == 1000 && Min[r] > 0 && "
                "Abs[Max[r] - 3943.247394845270952389728107753445640965263639345647460598568574874867`70] "
                "< 10^-44 && Precision[Max[r]] >= 49]");
}

/* Degree-200..300 orthogonal polynomials at 30 digits: every zero real, and
 * sampled zeros bracketed to 10^-27 relative by a sign change of the
 * polynomial evaluated numerically at 90 digits (its own recurrence, not the
 * cancelling monomial expansion the root finder works on). */
static bool polished_ok(const char* fn, int n, const char* idx) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Module[{f = %s, r}, r = Sort[x /. NSolve[f[x] == 0, x, WorkingPrecision -> 30]]; "
        "Length[r] == %d && Count[r, _Real] == %d && "
        "And @@ Table[Sign[f[SetPrecision[r[[k]] (1 - 10^-27), 90]]] "
        "Sign[f[SetPrecision[r[[k]] (1 + 10^-27), 90]]] == -1, {k, %s}]]",
        fn, n, n, idx);
    return lang_true(buf);
}

static void test_orthogonal_high_degree(void) {
    ASSERT_MSG(polished_ok("LaguerreL[200, #] &", 200, "{1, 100, 200}"), "LaguerreL[200]");
    ASSERT_MSG(polished_ok("HermiteH[300, #] &", 300, "{1, 150, 151, 300}"), "HermiteH[300]");
    ASSERT_MSG(polished_ok("LegendreP[256, #] &", 256, "{1, 128, 256}"), "LegendreP[256]");
    /* ChebyshevT zeros are known in closed form: cos((2k-1) Pi / (2n)). */
    ASSERT_TRUE("Module[{r = Sort[x /. NSolve[ChebyshevT[400, x] == 0, x, WorkingPrecision -> 30]]}, "
                "Length[r] == 400 && Max[Abs[r - Sort[N[Table[Cos[(2 k - 1) Pi/800], {k, 400}], 40]]]] < 10^-28]");
}

/* Machine precision on the same polynomials: the coefficients overflow a
 * double (L_1000) or cancel far beyond it, so the machine answer must come
 * from the escalating engine -- correct digits, never garbage. */
static void test_machine_high_degree(void) {
    ASSERT_TRUE("Module[{r = x /. NSolve[LaguerreL[1000, x] == 0, x]}, "
                "Length[r] == 1000 && Count[r, _Real] == 1000 && "
                "Abs[Max[r] - 3943.2473948452709523897] < 10^-9 && Precision[Max[r]] === MachinePrecision]");
    ASSERT_TRUE("Module[{r = Sort[List @@ NRoots[ChebyshevT[300, x] == 0, x][[All, 2]]]}, "
                "Length[r] == 300 && Max[Abs[r - Sort[N[Table[Cos[(2 k - 1) Pi/600], {k, 300}]]]]] < 10^-14]");
}

/* Wilkinson's polynomial: exact integer roots, a classical ill-conditioned
 * expansion (perturbing one coefficient by 2^-23 moves roots into the
 * complex plane), at machine precision and at 30 digits. */
static void test_wilkinson(void) {
    ASSERT_TRUE("Module[{r = Sort[x /. NSolve[Product[x - k, {k, 20}] == 0, x]]}, "
                "Length[r] == 20 && Max[Abs[r - Range[20]]] < 10^-12]");
    ASSERT_TRUE("Module[{r = Sort[x /. NSolve[Product[x - k, {k, 20}] == 0, x, WorkingPrecision -> 30]]}, "
                "Length[r] == 20 && Max[Abs[r - Range[20]]] < 10^-28]");
    /* Close rational roots k/40, degree 40. */
    ASSERT_TRUE("Module[{r = Sort[x /. NSolve[Product[40 x - k, {k, 40}] == 0, x, WorkingPrecision -> 25]]}, "
                "Length[r] == 40 && Max[Abs[r - Range[40]/40]] < 10^-23]");
}

/* Mixed real and complex roots at high degree with a cancelling expansion:
 * the real walk takes the real zeros, Aberth the rest. */
static void test_mixed_high_degree(void) {
    ASSERT_TRUE("Module[{r = x /. NSolve[(x^2 + 1) LaguerreL[60, x] == 0, x, WorkingPrecision -> 30]}, "
                "Length[r] == 62 && Count[r, _Real] == 60 && "
                "Max[Abs[Sort[Select[r, Head[#] === Complex &], Im[#1] < Im[#2] &] - {-I, I}]] < 10^-28]");
    ASSERT_TRUE("Module[{r = x /. NSolve[x^100 - 1 == 0, x]}, Length[r] == 100 && "
                "Max[Abs[r^100 - 1]] < 10^-10 && Count[r, _Real] == 2]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_real_roots);
    TEST(test_cubic_example);
    TEST(test_complex_roots);
    TEST(test_complex_coefficients);
    TEST(test_multiplicity);
    TEST(test_zero_roots);
    TEST(test_precision);
    TEST(test_degenerate);
    TEST(test_method_aberth);
    TEST(test_method_companion);
    TEST(test_method_jenkinstraub);
    TEST(test_stress_all_methods);
    TEST(test_stress_high_multiplicity);
    TEST(test_unevaluated);
    TEST(test_attributes);
    TEST(test_memory_loop);
    TEST(test_laguerre_1000_wp50);
    TEST(test_orthogonal_high_degree);
    TEST(test_machine_high_degree);
    TEST(test_wilkinson);
    TEST(test_mixed_high_degree);

    printf("All nroots_tests passed.\n");
    return 0;
}

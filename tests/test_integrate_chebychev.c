/* test_integrate_chebychev.c
 *
 * Unit tests for Integrate`ChebychevAlgebraic -- the Chebychev binomial
 * differential integrator (x^p (a x^r + b)^q).  Examples are drawn from
 * chebychev_problems.tex: problem 1 (compute), problem 3 (compute the
 * elementary, decline the rest), the worked Type III example, and a couple of
 * problem-4 cases (non-elementary elliptic integrals).
 *
 * Correctness of an antiderivative is asserted two ways together, because a
 * differentiate-back check alone is fooled by an *unevaluated* integral
 * (D[Integrate[f,x],x] - f == 0 trivially):
 *   1. the integral actually closes      -- FreeQ[Integrate[f,x], Integrate],
 *   2. it differentiates back to f       -- PossibleZeroQ[D[..] - f].
 * PossibleZeroQ (numeric sampling over x and the symbolic parameters a, b) is
 * used instead of Simplify, whose radical residues can stall.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Evaluate `input`, assert its printed form equals `expected` (always aborts on
 * mismatch, even under NDEBUG where libc assert() is compiled out). */
static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* Assert that Integrate[f, x] yields a closed, correct antiderivative. */
static void ok(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "FreeQ[Integrate[%s, x], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[D[Integrate[%s, x], x] - (%s)]", f, f);
    check_eq(buf, "True");
}

/* Assert that the forced ChebychevAlgebraic method declines (non-elementary or
 * out of shape) and leaves the integral unevaluated. */
static void declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Head[Integrate[%s, x, Method -> \"ChebychevAlgebraic\"]]", f);
    check_eq(buf, "Integrate");
}

/* ------------------------------------------------------------------ */
/* Problem 1: compute the Chebychev integrals (Types I, II, III).     */
/* ------------------------------------------------------------------ */
static void test_problem1(void) {
    ok("x^(1/3) (a Sqrt[x] + b)^3");          /* (i)   Type I   */
    ok("x^(3/2) (a x^(5/6) + b)^(1/3)");      /* (ii)  Type II  */
    ok("x^(4/3)/(a Sqrt[x] + b)^(2/3)");      /* (iii) Type III */
    ok("(a x^(3/4) + b)^(2/3)/x^3");          /* (iv)  Type III */
}

/* ------------------------------------------------------------------ */
/* Problem 3: compute the elementary ones, decline the rest.          */
/* ------------------------------------------------------------------ */
static void test_problem3_elementary(void) {
    ok("(a x^(1/5) + b)^(3/4)");              /* (i)    Type II  */
    ok("x^(1/3) (a x^2 + b)^(1/3)");          /* (iii)  Type III */
    ok("1/(a x^3 + b)^(1/3)");                /* (v)    Type II  */
    ok("x^4/(a x^4 + b)^(5/4)");              /* (vii)  Type III */
    ok("x^(4/5) (a/x^(4/5) + b)^(5/4)");      /* (viii) Type III, r<0 */
    ok("x^(3/2)/(a/Sqrt[x] + b)^(1/2)");      /* (x)    Type II,  r<0 */
}

static void test_problem3_nonelementary(void) {
    declines("(a x^2 + b)^(1/3)");            /* (ii)  5/6   */
    declines("Sqrt[a x^3 + b]");              /* (iv)  5/6   */
    declines("x/(a x^3 + b)^(1/3)");          /* (vi)  1/3   */
    declines("x^(3/5)/(a x^(3/5) + b)^(3/5)");/* (ix)  31/15 */
}

/* ------------------------------------------------------------------ */
/* The worked Type III example from the chapter text.                 */
/* ------------------------------------------------------------------ */
static void test_worked_example(void) {
    /* Integrate[x^(-5/6) Sqrt[x^(1/3) - 1], x]
     *   = 3 x^(1/6) Sqrt[x^(1/3) - 1] - 3 ArcTanh[Sqrt[x^(1/3)-1]/x^(1/6)].
     * (The text prints ArcTan, but the real antiderivative of a t^2 - 1
     * denominator is ArcTanh; the differentiate-back check is authoritative.) */
    ok("x^(-5/6) Sqrt[x^(1/3) - 1]");
}

/* ------------------------------------------------------------------ */
/* Problem 4(i): 1/Sqrt[x^(2n) +/- 1] is non-elementary for n > 1.    */
/* ------------------------------------------------------------------ */
static void test_problem4_nonelementary(void) {
    declines("1/Sqrt[x^4 + 1]");              /* n = 2, elliptic   */
    declines("1/Sqrt[x^6 - 1]");              /* n = 3, hyperell.  */
}

/* ------------------------------------------------------------------ */
/* Bases that are x^s (a x^r + b) only after expansion: a monomial    */
/* shared by both terms, or constants * x^s * binomial (DSolve        */
/* AutonomousReduction quadratures 2.1.2-267 / 2.1.2-1168).            */
/* ------------------------------------------------------------------ */
static void test_common_monomial_base(void) {
    ok("1/Sqrt[x^(5/2) + x^3]");              /* x^(-5/4)(1+x^(1/2))^(-1/2) */
    ok("-(1/Sqrt[2 a x^(5/2) + 4 b x^3])");   /* 2.1.2-267, symbolic       */
    ok("(-6 x^2)/(a Sqrt[(-12 x^3 (-1 - (2 b x^3)/a))/a])"); /* 2.1.2-1168 */
    ok("Sqrt[x^2 + x^4]");                    /* integer s: W kept         */
    ok("x Sqrt[a x^2 + b x^4]");
    ok("x^(1/3) (x^2 + x^(8/3))^(1/3)");      /* odd-denominator s         */
    /* The x < 0 branch: Sqrt[x^2 + x^4] = -x Sqrt[1 + x^2] there, so the
     * antiderivative must carry the locally constant Sqrt[x^2+x^4]/(x Sqrt[..])
     * rather than x Sqrt[1 + x^2]: its derivative is real and equal at x = -2. */
    check_eq("With[{F = Integrate[Sqrt[x^2 + x^4], x]}, "
             "Chop[N[(D[F, x] - Sqrt[x^2 + x^4]) /. x -> -2]] == 0]", "True");
    /* Still non-elementary after normalisation: x^(9/2)(1+x^2)^(1/2). */
    declines("x^2 Sqrt[x^5 + x^7]");
}

/* ------------------------------------------------------------------ */
/* Non-elementary binomials in the Automatic cascade (v0.6763): the    */
/* 2F1 antiderivative, valid on every branch, or a named form.         */
/* ------------------------------------------------------------------ */
/* The antiderivative of f closes and differentiates back to f at points on
 * both sides of every branch point (|x| < 1 and |x| > 1, x < 0 and x > 0). */
static void special_ok(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "With[{F = Quiet[Integrate[%s, x]]}, FreeQ[F, Integrate] && "
             "AllTrue[{3/10, 7/10, -2/5, 13/10, -2}, "
             "Abs[Quiet[N[(D[F, x] - (%s)) /. {x -> #, a -> 3/7, b -> -5/3}, 20]]] < 10^-12 &]]",
             f, f);
    check_eq(buf, "True");
}

static void test_nonelementary_hypergeometric(void) {
    special_ok("x^3/Sqrt[2 + x^6]");
    special_ok("1/Sqrt[-1 + x^6]");          /* b < 0: (1 + a x^r/b)^(-q) kept */
    special_ok("x^(-1/3)/Sqrt[1 - x]");
    special_ok("(1 + x^3)^(1/3)");
    special_ok("x^2 (b + a x^5)^(1/3)");     /* symbolic a, b */
    /* the named form FunctionExpand finds is kept */
    check_eq("FreeQ[Integrate[1/Sqrt[1 + x^4], x], Hypergeometric2F1 | HypergeometricPFQ]",
             "True");
    special_ok("1/Sqrt[1 + x^4]");
    /* the exact 2F1 form, and no Integrate::nonelem once it closes */
    check_eq("Check[Integrate[x^3/Sqrt[2 + x^6], x], msg, {Integrate::nonelem}]",
             "(x^4 Hypergeometric2F1[1/2, 2/3, 5/3, -1/2 x^6])/(4 Sqrt[2])");
    /* FunctionExpand of a 2F1 with no named form declines promptly: its
     * integration probes are binomials of this class (was ~6.6 s) */
    check_eq("First[AbsoluteTiming[FunctionExpand[Hypergeometric2F1[2/3, 1/2, 5/3, z]]]] < 3",
             "True");
}

/* ------------------------------------------------------------------ */
/* Method plumbing + strictness.                                      */
/* ------------------------------------------------------------------ */
static void test_plumbing(void) {
    /* Direct package head closes and round-trips. */
    check_eq("FreeQ[Integrate`ChebychevAlgebraic[x^(1/3) (a Sqrt[x] + b)^3, x],"
             " Integrate`ChebychevAlgebraic]", "True");
    check_eq("PossibleZeroQ[D[Integrate`ChebychevAlgebraic["
             "x^(1/3) (a Sqrt[x] + b)^3, x], x] - x^(1/3) (a Sqrt[x] + b)^3]",
             "True");

    /* Method option reaches the same routine. */
    check_eq("FreeQ[Integrate[x^(1/3) (a Sqrt[x] + b)^3, x,"
             " Method -> \"ChebychevAlgebraic\"], Integrate]", "True");

    /* Strict: a non-binomial integrand is declined (no fallback). */
    declines("Sin[x]");

    /* Automatic cascade routes a Chebychev binomial that the rational stage
     * cannot handle to this method. */
    check_eq("FreeQ[Integrate[x^(3/2) (a x^(5/6) + b)^(1/3), x], Integrate]",
             "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_problem1);
    TEST(test_problem3_elementary);
    TEST(test_problem3_nonelementary);
    TEST(test_worked_example);
    TEST(test_problem4_nonelementary);
    TEST(test_common_monomial_base);
    TEST(test_nonelementary_hypergeometric);
    TEST(test_plumbing);

    printf("All Integrate ChebychevAlgebraic tests passed!\n");
    return 0;
}

/* test_integrate_goursat.c
 *
 * Unit tests for Integrate`GoursatAlgebraic -- Goursat's pseudo-elliptic
 * algorithm and its cube-/fourth-root generalisations for F(t)/R(t)^p,
 * p in {1/2, 1/3, 2/3, 1/4, 3/4}.  Examples are the worked cases from
 * GoursatAppendix.wl / GoursatExamples.wl and the two preprints.
 *
 * Correctness is asserted two ways together, because a differentiate-back
 * check alone is fooled by an *unevaluated* integral:
 *   1. the integral actually closes  -- FreeQ[Integrate[f,t], Integrate],
 *   2. it differentiates back to f    -- PossibleZeroQ[D[..] - f].
 * (The method itself already runs a differentiate-back guard internally, so a
 * closed result is also a verified one; the test re-checks independently.)
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

/* Assert Integrate[f, t] yields a closed, correct antiderivative. */
static void ok(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "FreeQ[Integrate[%s, t], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[D[Integrate[%s, t], t] - (%s)]", f, f);
    check_eq(buf, "True");
}

/* Assert the forced GoursatAlgebraic method declines (non-elementary /
 * obstructed / non-harmonic / out of shape) and stays unevaluated. */
static void declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Head[Integrate[%s, t, Method -> \"GoursatAlgebraic\"]]", f);
    check_eq(buf, "Integrate");
}

/* Like `ok`, but FORCES Method -> "GoursatAlgebraic" so the result is
 * guaranteed to come from the Goursat integrator (not the Automatic cascade):
 * the integral closes and differentiates back to f. */
static void okg(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "FreeQ[Integrate[%s, t, Method -> \"GoursatAlgebraic\"], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[D[Integrate[%s, t, Method -> \"GoursatAlgebraic\"], t]"
             " - (%s)]", f, f);
    check_eq(buf, "True");
}

/* Forced GoursatAlgebraic, but verify by a NUMERIC differentiate-back at the
 * real point t = `pt` (a string like "17/5").  For closed forms carrying
 * complex Logs / large nested radicals on which PossibleZeroQ's real-only
 * sampling misfires across a branch cut.  Closure already implies the method's
 * internal diff_back_ok guard accepted the result. */
static void okg_num(const char* f, const char* pt) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "FreeQ[Integrate[%s, t, Method -> \"GoursatAlgebraic\"], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof(buf),
             "N[Abs[(D[Integrate[%s, t, Method -> \"GoursatAlgebraic\"], t] - (%s))"
             " /. t -> %s], 20] < 1/100000000000", f, f, pt);
    check_eq(buf, "True");
}

/* ------------------------------------------------------------------ */
/* Graded exercise ladder (easy -> hard), mirroring GOURSAT_EXERCISES.md.  */
/* Every case is FORCED through Method -> "GoursatAlgebraic", so each is     */
/* genuinely solved by the Goursat integrator.  The ladder spans positive   */
/* (radical in the numerator) and negative (denominator) exponents for      */
/* p in {1/2, 1/3, 2/3, 1/4, 3/4} and exercises every involution equation.  */
/* ------------------------------------------------------------------ */
static void test_graded(void) {
    /* -- Tier 1: recognise + fold; antiderivative is a single radical power. */
    okg("t^2/(t^3-1)^(1/3)");      /* neg p=1/3  -> (t^3-1)^(2/3)/2  */
    okg("t^2/(t^3-1)^(2/3)");      /* neg p=2/3  -> (t^3-1)^(1/3)    */
    okg("t^3/(t^4-1)^(3/4)");      /* neg p=3/4  -> (t^4-1)^(1/4)    */
    okg("t^3/(t^4-1)^(1/4)");      /* neg p=1/4  -> (t^4-1)^(3/4)/3  */

    /* -- Tier 2: positive exponent (radical in the numerator). */
    okg("t^2 (t^3-1)^(1/3)");      /* pos p=1/3 -> (t^3-1)^(4/3)/4   */
    okg("t^2 (t^3-1)^(2/3)");      /* pos p=2/3 -> (t^3-1)^(5/3)/5   */
    okg("t^3 (t^4-1)^(1/4)");      /* pos p=1/4 -> (t^4-1)^(5/4)/5   */
    okg("t^3 (t^4-1)^(3/4)");      /* pos p=3/4 -> (t^4-1)^(7/4)/7   */

    /* -- Tier 3: genuine cube-root involution (ArcTan + Log). */
    okg("1/(t^3-1)^(1/3)");        /* p=1/3, H1 == 0                 */
    okg("t/(t^3-1)^(2/3)");        /* p=2/3, H0 == 0 (dual)          */
    okg("1/(t (t^3-1)^(2/3))");    /* p=2/3, H0 == 0, pole cofactor  */
    okg("1/(t^3-8)^(1/3)");        /* p=1/3, shifted radicand        */

    /* -- Tier 4: genuine fourth-root involution (ArcTan + ArcTanh),
     *    harmonic quartic. */
    okg("1/(t^4-1)^(1/4)");        /* p=1/4, H1 == H2 == 0           */
    okg("t^2/(t^4-1)^(3/4)");      /* p=3/4, H0 == H1 == 0 (dual)    */
    okg("1/(t^4-16)^(1/4)");       /* harmonic quartic != t^4-1      */
    okg("t^2/(t^4-16)^(3/4)");     /* harmonic quartic != t^4-1      */

    /* -- Tier 5: Goursat's original square-root V4 (ArcTanh + Log). */
    okg("t/Sqrt[(t^2-1)(t^2-4)]");                /* V4, P0 == 0      */
    okg("(t^2-2)/(t Sqrt[(t^2-1)(t^2-4)])");      /* finite-f.p. S=2/t */
    okg("(t^2-2)/(t Sqrt[(t^2-4)(t^2-9)])");      /* different radicand */
    okg("(t^4 + 2 t^3 - 4)/(t^2 Sqrt[(t^2-1)(t^2-4)])"); /* neg-lc descent */

    /* -- Tier 6: positive square-root exponent (numerator; large closed
     *    form -> numeric differentiate-back). */
    okg_num("t Sqrt[(t^2-1)(t^2-4)]", "7/2");     /* pos p=1/2        */

    /* -- Tier 7 (hardest): period-3 higher symmetry.  Order-3 Mobius S with
     *    F a non-trivial character (F + F(S) + F(S^2) == 0); the V4 trivial
     *    projection is non-zero.  Complex-Log closed form -> numeric check. */
    okg_num("(t-1)/((t+2) Sqrt[t^3-1])", "17/5");

    /* -- Involution gates genuinely fail: these must DECLINE. */
    declines("t/(t^3-1)^(1/3)");          /* cube H1 != 0            */
    declines("1/(t^3-1)^(2/3)");          /* cube H0 != 0            */
    declines("t/(t^4-1)^(1/4)");          /* fourth V1 obstructive   */
    declines("t^2 (t^4-1)^(3/4)");        /* non-elementary (Chebyshev) */
    declines("t^2/Sqrt[(t^2-1)(t^2-4)]"); /* V4-invariant: elliptic  */
    declines("1/Sqrt[t^3-1]");            /* period-3 proj. != 0     */
}

/* ------------------------------------------------------------------ */
/* Square-root case: Goursat's V4 (p = 1/2).                          */
/* ------------------------------------------------------------------ */
static void test_sqrt(void) {
    ok("t/Sqrt[(t^2-1)(t^2-4)]");                 /* Example 5.1     */
    ok("(t^2-2)/(t Sqrt[(t^2-1)(t^2-4)])");       /* finite f.p. S=2/t */
    ok("(t + (t^2-2)/t)/Sqrt[(t^2-1)(t^2-4)]");   /* two components  */

    /* Regression: a component whose S2-involution reduction has a NEGATIVE
     * leading coefficient lc in its descended quadratic Sqrt[lc Q].  The radical
     * must stay Sqrt[lc Q]; factoring Sqrt[lc] into the prefactor is a branch-cut
     * bug (Sqrt[lc] Sqrt[Q] != Sqrt[lc Q] when Q < 0) that previously left this
     * integral unevaluated with a 1/0 message. */
    ok("(t^4 + 2 t^3 - 4)/(t^2 Sqrt[(t^2-1)(t^2-4)])");

    declines("t^2/Sqrt[(t^2-1)(t^2-4)]");         /* V4-invariant: elliptic */
    declines("t/Sqrt[t^3+1]");                    /* cubic radicand: elliptic */
    declines("t/Sqrt[t^4+t+1]");                  /* irreducible quartic: elliptic */
}

/* ------------------------------------------------------------------ */
/* Cube-root case (p = 1/3, 2/3).                                     */
/* ------------------------------------------------------------------ */
static void test_cube(void) {
    ok("1/(t^3-1)^(1/3)");                        /* Example 5.2          */
    ok("t^2/(t^3-1)^(1/3)");                      /* -> (t^3-1)^(2/3)/2   */
    ok("t^2/(t^3-1)^(2/3)");                      /* -> (t^3-1)^(1/3)     */

    declines("t/(t^3-1)^(1/3)");                  /* obstructed (H1 != 0) */
    declines("1/(t^2+t+1)^(1/3)");                /* quadratic radicand   */
}

/* ------------------------------------------------------------------ */
/* Fourth-root case (p = 1/4, 3/4) on the harmonic quartic t^4 - 1.   */
/* ------------------------------------------------------------------ */
static void test_fourth(void) {
    ok("1/(t^4-1)^(1/4)");                        /* V0 elementary        */
    ok("t^3/(t^4-1)^(1/4)");                      /* -> (t^4-1)^(3/4)/3   */
    ok("t^2/(t^4-1)^(3/4)");                      /* V2 elementary (dual) */
    ok("t^3/(t^4-1)^(3/4)");                      /* -> (t^4-1)^(1/4)     */

    declines("t/(t^4-1)^(1/4)");                  /* V1 obstructive       */
    declines("t^2/(t^4-1)^(1/4)");                /* V2 obstructive       */
    declines("1/(t^4-1)^(3/4)");                  /* V0 obstructive (flip)*/
    declines("1/((t-1)(t-2)(t-3)(t-5))^(1/4)");   /* non-harmonic         */
}

/* ------------------------------------------------------------------ */
/* Period-3 higher-symmetry case for sqrt(cubic) (Goursat 1887 Sec.4). */
/* R = t^3 - 1 has an order-3 Mobius S fixing one ramification point    */
/* and cycling the other three; F is a non-trivial period-3 character   */
/* (F(S) = alpha F).  NOT covered by the V4 Theorems 1-2: the V4 trivial */
/* projection is non-zero, but F + F(S) + F(S^2) vanishes.              */
/* ------------------------------------------------------------------ */
static void test_period3(void) {
    /* Goursat's Section-4 worked example:
     *   (1/3) dx/Sqrt[x(1-x)] = (t-1)/(t+2) dt/Sqrt[t^3-1],  x=((t-1)/(t+2))^3,
     * i.e. Integrate[(t-1)/((t+2) Sqrt[t^3-1]), t] = -(2/3) ArcTan[...].
     * The closed antiderivative carries un-simplified nested radicals on which
     * PossibleZeroQ's real-only sampling misfires across the t<1 branch cut, so
     * (unlike `ok`) verify closure plus a numeric differentiate-back at a real
     * point t>1 where the integrand is real.  (Closure already implies the
     * method's internal diff_back_ok guard accepted it.) */
    check_eq("FreeQ[Integrate[(t-1)/((t+2) Sqrt[t^3-1]), t], Integrate]", "True");
    check_eq("FreeQ[Integrate[(t-1)/((t+2) Sqrt[t^3-1]), t,"
             " Method -> \"GoursatAlgebraic\"], Integrate]", "True");
    check_eq("N[Abs[(D[Integrate[(t-1)/((t+2) Sqrt[t^3-1]), t], t]"
             " - (t-1)/((t+2) Sqrt[t^3-1])) /. t -> 17/5], 20] < 1/100000000000",
             "True");

    /* Regression (output normalisation): the antiderivative is now returned over
     * the integrand's OWN radical Sqrt[t^3-1] -- not the Mobius-substituted
     * nested radical Sqrt[w^2-w], w=((t-1)/(t+2))^3 -- so it verifies
     * SYMBOLICALLY via D[..] // Simplify (a single shared radical reduces over
     * the rational field), the user-facing behaviour the form-fix restores. */
    check_eq("Simplify[D[Integrate[(t-1)/((t+2) Sqrt[t^3-1]), t], t]"
             " - (t-1)/((t+2) Sqrt[t^3-1])]", "0");

    /* CYCLOTOMIC-fixed-point character: the order-3 Mobius fixing the cube root
     * of unity alpha = (-1)^(2/3) has fixed points {alpha, -2 alpha} (Goursat's
     * ((t-alpha)/(t+2alpha))^3 substitution).  Previously this hung
     * (ComplexInfinity from Solve+Simplify returning nested radicals
     * Sqrt[-9(-1)^(1/3)] that left roots-of-unity uncollapsed in
     * to_function_of_power); now the clean fixed points keep the reduction in
     * Q(zeta_3) and the descended Sqrt[quadratic] integrates after monic
     * rescaling.  Closure + real-axis differentiate-back (the closed form is a
     * complex Log on which PossibleZeroQ's real sampling misfires). */
    check_eq("FreeQ[Integrate[(t-(-1)^(2/3))/((t+2(-1)^(2/3)) Sqrt[t^3-1]), t,"
             " Method -> \"GoursatAlgebraic\"], Integrate]", "True");
    check_eq("N[Abs[(D[Integrate[(t-(-1)^(2/3))/((t+2(-1)^(2/3)) Sqrt[t^3-1]), t], t]"
             " - (t-(-1)^(2/3))/((t+2(-1)^(2/3)) Sqrt[t^3-1])) /. t -> 3], 20]"
             " < 1/100000000000", "True");

    /* Genuinely elliptic cubic radicands still decline (the period-3 trivial
     * projection does not vanish / no order-3 character). */
    declines("1/Sqrt[t^3-1]");
}

/* ------------------------------------------------------------------ */
/* Method plumbing + strictness.                                      */
/* ------------------------------------------------------------------ */
static void test_plumbing(void) {
    /* Direct package head closes and round-trips. */
    check_eq("FreeQ[Integrate`GoursatAlgebraic[1/(t^3-1)^(1/3), t],"
             " Integrate`GoursatAlgebraic]", "True");
    check_eq("PossibleZeroQ[D[Integrate`GoursatAlgebraic[1/(t^3-1)^(1/3), t], t]"
             " - 1/(t^3-1)^(1/3)]", "True");

    /* Method option reaches the same routine. */
    check_eq("FreeQ[Integrate[t/Sqrt[(t^2-1)(t^2-4)], t,"
             " Method -> \"GoursatAlgebraic\"], Integrate]", "True");

    /* Strict: a non-pseudo-elliptic integrand is declined (no fallback). */
    declines("Sin[t]");

    /* Automatic cascade routes a pseudo-elliptic integrand here. */
    check_eq("FreeQ[Integrate[t^2/(t^3-1)^(2/3), t], Integrate]", "True");
}

/* ------------------------------------------------------------------ */
/* Third-kind cube-root logarithmic part.                             */
/*                                                                    */
/* When the order-3 eigendescent obstruction H1 != 0 but the cofactor */
/* F has a pole at a NON-branch point, F/R^(1/3) is still elementary   */
/* and its antiderivative is a sum of logs of R^(1/3) - omega^j kappa t.*/
/* Closed forms carry complex Logs -> numeric differentiate-back.      */
/* ------------------------------------------------------------------ */
static void test_thirdkind(void) {
    /* Numeric instance (k = 2) of the parametric family. */
    okg_num("(2 - 3 t)/((1 - 3 t) (t (1 - t) (1 - 2 t))^(1/3))", "17/5");
    /* A second numeric member (k = 3). */
    okg_num("(2 - 4 t)/((1 - 4 t) (t (1 - t) (1 - 3 t))^(1/3))", "17/5");

    /* Parametric radicand k: closes and differentiates back (k pinned for the
     * decisive numeric check). */
    check_eq("FreeQ[Integrate[(2 - (k+1) t)/((1 - (k+1) t) (t (1-t)(1-k t))^(1/3)),"
             " t, Method -> \"GoursatAlgebraic\"], Integrate]", "True");
    check_eq("N[Abs[(D[Integrate[(2 - (k+1) t)/((1 - (k+1) t) (t (1-t)(1-k t))^(1/3)),"
             " t, Method -> \"GoursatAlgebraic\"], t]"
             " - (2 - (k+1) t)/((1 - (k+1) t) (t (1-t)(1-k t))^(1/3)))"
             " /. {k -> 3, t -> 17/5}], 20] < 1/100000000000", "True");

    /* The Automatic cascade also lands the parametric integral (no Method). */
    check_eq("FreeQ[Integrate[(2 - (k+1) t)/((1 - (k+1) t) (t (1-t)(1-k t))^(1/3)), t],"
             " Integrate]", "True");

    /* A polynomial cofactor (no non-branch pole) is genuinely second-kind /
     * elliptic and must still DECLINE -- the third-kind ansatz must not
     * over-reach and diff_back must reject a spurious log form. */
    declines("t/(t^3-1)^(1/3)");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_sqrt);
    TEST(test_cube);
    TEST(test_fourth);
    TEST(test_period3);
    TEST(test_thirdkind);
    TEST(test_plumbing);
    TEST(test_graded);

    printf("All Integrate GoursatAlgebraic tests passed!\n");
    return 0;
}

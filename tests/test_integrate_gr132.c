/* test_integrate_gr132.c
 *
 * Regression suite for issue #132: wrong Integrate results on the
 * Gradshteyn & Ryzhik (7th ed.) corpus.  Every row of the issue is pinned
 * here by its reference value at each listed instantiation, so a row cannot
 * silently go wrong again.  The reference values were checked independently
 * (NIntegrate at two instantiations, differentiation back for an
 * antiderivative) before the issue was filed.
 *
 * Two checks:
 *   check_val(integral, rules, ref)  -- the integral evaluates to a closed
 *       form free of Integrate, and N[form /. rules] agrees with ref to a
 *       relative 1e-6 (absolute when |ref| < 1);
 *   check_antideriv(f, x)            -- Integrate[f, x] closes and
 *       D[F, x] - f vanishes at x = 31/100, 77/100, 143/100.
 *
 * The rows are grouped by the root cause each one exposed; the comment above
 * each group names the general rule that was wrong.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failures = 0;

static Expr* eval_str(const char* s) {
    Expr* p = parse_expression(s);
    if (!p) return NULL;
    Expr* r = evaluate(p);
    expr_free(p);
    return r;
}

static bool is_true(Expr* e) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, "True") == 0;
}

/* `integral` closes, and at each instantiation in `rules` (a Mathilda list of
 * rule lists, or "{{}}" for none) its value agrees with the matching entry of
 * `refs` (a Mathilda list of numbers). */
static void check_val(const char* integral, const char* rules, const char* refs) {
    size_t n = strlen(integral) + strlen(rules) + strlen(refs) + 512;
    char* buf = malloc(n);
    snprintf(buf, n,
        "Module[{r = Quiet[%s], rs = %s, vs = %s, v},"
        "  FreeQ[r, Integrate] && And @@ Table["
        "    v = Quiet[N[r /. rs[[i]], 20]];"
        "    NumericQ[v] && Abs[v - vs[[i]]] < 10^-6 Max[1, Abs[vs[[i]]]],"
        "  {i, Length[rs]}]]",
        integral, rules, refs);
    Expr* r = eval_str(buf);
    if (!is_true(r)) {
        Expr* got = eval_str(integral);
        char* gs = got ? expr_to_string(got) : NULL;
        fprintf(stderr, "FAIL: %s\n  at %s\n  expected %s\n  got: %s\n",
                integral, rules, refs, gs ? gs : "(null)");
        free(gs);
        if (got) expr_free(got);
        g_failures++;
    }
    if (r) expr_free(r);
    free(buf);
}

/* Integrate[f, x] closes and differentiates back to f at three points. */
static void check_antideriv(const char* f) {
    size_t n = 2 * strlen(f) + 512;
    char* buf = malloc(n);
    snprintf(buf, n,
        "Module[{F = Quiet[Integrate[%s, x]], g = %s},"
        "  FreeQ[F, Integrate] && And @@ Table["
        "    Abs[Quiet[N[(D[F, x] - g) /. x -> p, 30]]] < 10^-10,"
        "  {p, {31/100, 77/100, 143/100}}]]",
        f, f);
    Expr* r = eval_str(buf);
    if (!is_true(r)) {
        fprintf(stderr, "FAIL: antiderivative of %s\n", f);
        g_failures++;
    }
    if (r) expr_free(r);
    free(buf);
}

/* -------------------------------------------------------------------------
 * Ei off the real axis.  Ei(z) -> I Pi Sign[Im z] as Re z -> -Infinity, not 0
 * (Gruntz's Ei expansion lacked the connection constant; Limit composed
 * Ei[Limit[g]] after the imaginary part of g was absorbed into -Infinity; the
 * Newton-Leibniz boundary substitution Ei[Log[0] - 7 I/3] -> Ei[-Infinity]
 * did the same).  Every Laplace-type integral of 1/(t^2 + a^2) went through
 * one of these.
 * ---------------------------------------------------------------------- */
static void test_ei_connection(void) {
    check_val("Limit[ExpIntegralEi[t - I], t -> -Infinity]", "{{}}", "{-Pi I}");
    check_val("Limit[ExpIntegralEi[t + 7/3 I], t -> -Infinity]", "{{}}", "{Pi I}");
    check_val("Limit[ExpIntegralEi[t - I a], t -> -Infinity, Assumptions -> a > 0]",
              "{{}}", "{-Pi I}");
    check_val("Limit[ExpIntegralEi[Log[x] - I], x -> 0, Direction -> \"FromAbove\"]",
              "{{}}", "{-Pi I}");
    check_val("Integrate[E^((-7*Tan[x])/4), {x, 0, Pi/2}]", "{{}}",
              "{0.4386142036879010}");                                  /* 3.341 */
    check_val("Integrate[(x*Sec[x]^2)/E^((7*Tan[x])/3), {x, 0, Pi/2}]", "{{}}",
              "{0.1525142744519606}");                                  /* 3.963.2 */
    check_val("Integrate[(x*Sec[x]^3*(-Cos[x] + (7*Sin[x])/3))/E^((7*Tan[x])/3), {x, 0, Pi/2}]",
              "{{}}", "{0.1159450517871212}");                          /* 3.964.1 */
    check_val("Integrate[(49/9 + Log[x]^2)^(-1), {x, 0, 1}]", "{{}}",
              "{0.1525142744519606}");                                  /* 4.213.1 */
    check_val("Integrate[Log[x]/(25/36 + Log[x]^2), {x, 0, 1}]", "{{}}",
              "{-0.4170403262489011}");                                 /* 4.213.3 */
    check_val("Integrate[Log[9/64 + x^2]/E^(x/5), {x, 0, Infinity}]", "{{}}",
              "{11.40046498543361}");                                   /* 4.338.1 */
    check_val("Integrate[ArcTan[(8*x)/5]/E^((5*x)/8), {x, 0, Infinity}]", "{{}}",
              "{1.506794487553691}");                                   /* 4.551.3 */
    check_val("Integrate[ArcCot[(8*x)/3]/E^(x/5), {x, 0, Infinity}]", "{{}}",
              "{1.150621463704887}");                                   /* 4.551.4 */
    /* the real-axis limits are unchanged */
    check_val("Limit[ExpIntegralEi[-t], t -> Infinity]", "{{}}", "{0}");
}

/* -------------------------------------------------------------------------
 * Fourier integrals of a rational function with a frequency (or pole) whose
 * sign the assumptions leave open.  The residue families classified the
 * closure at one representative point and returned that half-plane's answer
 * for the whole region; they now split the region at the sign of each such
 * real parameter.
 * ---------------------------------------------------------------------- */
static void test_fourier_sign_split(void) {
    check_val("Integrate[E^((-I)*p*x)/(a^2 + x^2), {x, -Infinity, Infinity}, "
              "Assumptions -> (a != 0) && (Element[a, Reals]) && (Element[p, Reals])]",
              "{{a -> 5/8, p -> 3/8}, {a -> 5/4, p -> 5/6}, {a -> -5/4, p -> -5/6}}",
              "{3.976326945213576, 0.8868491913697098, 0.8868491913697098}"); /* 3.354.5 */
    check_val("Integrate[E^(I*x*xi)/(a^2 + x^2), {x, -Infinity, Infinity}, "
              "Assumptions -> (Re[a] > 0) && (Element[xi, Reals])]",
              "{{a -> 5/6, xi -> 1/5}, {a -> 1/5, xi -> 1/5}, {a -> 1/5, xi -> -1/5}}",
              "{3.191160921977243, 15.09204521843798, 15.09204521843798}");  /* 17.23.14 */
    check_val("Integrate[(E^(I*x*xi)*x)/(a^2 + x^2), {x, -Infinity, Infinity}, "
              "Assumptions -> (Re[a] > 0) && (Element[xi, Reals])]",
              "{{a -> 5/8, xi -> 7/4}, {a -> 5/8, xi -> 5/8}, {a -> 5/8, xi -> -5/8}}",
              "{1.052301726914921 I, 2.125707920271893 I, -2.125707920271893 I}"); /* 17.23.15 */
    /* both closures still agree with the single-sign forms */
    check_val("Integrate[Cos[xi x]/(1 + x^2), {x, -Infinity, Infinity}, "
              "Assumptions -> Element[xi, Reals]]",
              "{{xi -> 2}, {xi -> -2}, {xi -> 0}}",
              "{0.4251682622789167, 0.4251682622789167, 3.141592653589793}");
}

/* -------------------------------------------------------------------------
 * Branch cuts of the antiderivative.  A smooth integrand antidifferentiated to
 * F with a two-argument ArcTan[u, v], which jumps by 2 Pi where v changes sign
 * with u < 0; F(b) - F(a) telescoped across the jump, and Limit read the
 * principal value on the cut from either side.  Now Limit takes one-sided
 * limits at a cut by the side's branch, and Newton-Leibniz splits at F's cut
 * crossings.
 * ---------------------------------------------------------------------- */
static void test_branch_cut_jumps(void) {
    check_val("Limit[ArcTan[-1, y], y -> 0, Direction -> \"FromBelow\"]", "{{}}", "{-Pi}");
    check_val("Limit[ArcTan[-1, y], y -> 0, Direction -> \"FromAbove\"]", "{{}}", "{Pi}");
    check_val("Limit[Log[-1 + I y], y -> 0, Direction -> \"FromBelow\"]", "{{}}", "{-I Pi}");
    check_val("Limit[Sqrt[-1 + I y], y -> 0, Direction -> \"FromBelow\"]", "{{}}", "{-I}");
    check_val("Limit[ArcTan[Cos[2x], -Sin[2x]], x -> Pi/2, Direction -> \"FromBelow\"]",
              "{{}}", "{-Pi}");
    check_val("Limit[Log[x - 2], x -> 1]", "{{}}", "{I Pi}");   /* real argument: on the cut */
    check_val("Integrate[(Sin[x]*Sin[11*x])/(26/25 - (2*Cos[2*x])/5), {x, 0, Pi}]",
              "{{}}", "{0.0004188790204786}");                          /* 3.613.8a */
    check_val("Integrate[Log[61/36 - (5*Cos[2*x])/3]*Sin[x]*Sin[9*x], {x, 0, Pi}]",
              "{{}}", "{-0.0631267813944708}");                         /* 4.397.11 */
    check_val("Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]", "{{}}", "{3.627598728468436}");
}

/* -------------------------------------------------------------------------
 * Rows already correct when this suite was written (fixed by earlier work);
 * pinned so they stay that way.
 * ---------------------------------------------------------------------- */
static void test_already_fixed(void) {
    check_antideriv("1/Sqrt[-Cos[(3*x)/4]]");                           /* 2.619.1 */
    check_val("Integrate[E^(-(-7/(4*x) + x)^2), {x, -Infinity, Infinity}]", "{{}}",
              "{1.772453850905516}");                                   /* 3.324.2 */
    check_val("Integrate[E^(-(b/x^2) - a*x^2), {x, 0, Infinity}, "
              "Assumptions -> (a > 0) && (b > 0)]",
              "{{a -> 1/5, b -> 7/4}, {a -> 5/8, b -> 3/8}}",
              "{0.6069679815244005, 0.4256975082080158}");              /* 3.325 S */
    check_val("Integrate[E^(-7/(4*x^2) - x^2/5), {x, 0, Infinity}]", "{{}}",
              "{0.6069679815244005}");                                  /* 3.325 N */
    check_val("Integrate[E^((-2*x^2)/5 - (7*x^4)/3), {x, 0, Infinity}]", "{{}}",
              "{0.6742393862812157}");                                  /* 3.469.1 N */
    check_val("Integrate[(-1 + E^(-(a/x^2)))/E^(mu*x^2), {x, 0, Infinity}, "
              "Assumptions -> (Re[mu] > 0) && (Re[a] > 0)]",
              "{{a -> 7/4, mu -> 7/4}, {a -> 5/8, mu -> 3/8}}",
              "{-0.6496946161125394, -0.8976293891793217}");            /* 3.472.1 S */
    check_val("Integrate[(-1 + E^(-7/(4*x^2)))/E^((7*x^2)/4), {x, 0, Infinity}]", "{{}}",
              "{-0.6496946161125394}");                                 /* 3.472.1 N */
    check_val("Integrate[E^(-5/(8*x^2) - (5*x^2)/8)*x^2, {x, 0, Infinity}]", "{{}}",
              "{0.5781084731486058}");                                  /* 3.472.2 N */
    check_val("Integrate[Sec[x]*(-Csc[x]^mu + Sin[x]^mu), {x, 0, Pi/2}, "
              "Assumptions -> (Abs[Re[mu]] < 1)]",
              "{{mu -> 3/8}, {mu -> 5/8}}",
              "{-1.049572550086451, -2.350862834655024}");              /* 3.685.3a */
    check_val("Integrate[((15*(1 - 2*x^2))/8 - (35*(1 - 2*x^2)^3)/4 + "
              "(63*(1 - 2*x^2)^5)/8)*Cos[(3*x)/8], {x, 0, 1}]", "{{}}",
              "{0.0909637787288092}");                                  /* 7.244.2 */
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_ei_connection);
    TEST(test_fourier_sign_split);
    TEST(test_branch_cut_jumps);
    TEST(test_already_fixed);

    if (g_failures) {
        printf("%d integrate_gr132 check(s) FAILED.\n", g_failures);
        return 1;
    }
    printf("All integrate_gr132 tests passed.\n");
    return 0;
}

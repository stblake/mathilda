/* test_risch_rde_tower.c
 *
 * Gap 1a — the recursive Risch differential equation over the transcendental
 * TOWER (Bronstein Symbolic Integration I, Chapter 6, lifted from the base field
 * C(x) to K_m = C(x)(t_0..t_m)).  Exercises, per component:
 *
 *   Risch`RischDE[f, g, x]            — the driver, base AND tower field
 *   Risch`SPDE[a, b, c, x, n]         — Rothstein's degree-reducing box (§6.4)
 *   Risch`PolyRischDENoCancel[b,c,x,n]— non-cancellation polynomial solve (§6.5)
 *
 * plus the integrator end-to-end on integrands whose Laurent-coefficient Risch DE
 * now routes through the literal Bronstein stack (rde_tower) instead of the bounded
 * ansatz.  Every "solve" assertion re-checks the EXACT identity D[y] + f y == g
 * (so a wrong reconstruction is caught regardless of form); every "decline"
 * asserts the call stays unevaluated (the recursive Risch procedure is never wrong,
 * only silent).  1a covers the PRIMITIVE (RT_LOG) top, non-cancellation branch; the
 * exponential top and the cancellation / antidifferentiation branches (which need
 * later increments / LimitedIntegrate) are asserted to decline cleanly.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "print.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* Evaluate `src`, return the FullForm string (caller frees). */
static char* eval_fullform(const char* src) {
    Expr* e = parse_expression(src);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    expr_free(r);
    expr_free(e);
    return s;
}

/* Risch`RischDE[f, g, x] with g = D[y]+f y must SOLVE (not stay unevaluated) and
 * the returned solution must satisfy D[sol]+f sol == g exactly. */
static void assert_rde_y(const char* f, const char* y) {
    char buf[4096];
    snprintf(buf, sizeof(buf),
        "With[{gg = Together[D[%s, x] + (%s)*(%s)]}, "
        "With[{sol = Risch`RischDE[%s, gg, x]}, "
        "If[Head[sol] === Risch`RischDE, $Unsolved, "
        "Together[D[sol, x] + (%s)*sol - gg]]]]",
        y, f, y, f, f);
    char* s = eval_fullform(buf);
    if (strcmp(s, "0") != 0)
        printf("FAIL: RischDE(f=%s, y=%s) -> %s (expected solved, residual 0)\n", f, y, s);
    ASSERT_MSG(strcmp(s, "0") == 0,
        "RischDE(f=%s, y=%s): expected solved with residual 0, got %s", f, y, s);
    free(s);
}

/* Risch`RischDE[f, g, x] with g given DIRECTLY (not g = D[y]+f y) must SOLVE and
 * the solution must satisfy D[sol]+f sol == g.  Used for antidifferentiation
 * (D h = g) where forming g = D[y] would drop the kernel that defines the field. */
static void assert_rde_solves(const char* f, const char* g) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "With[{sol = Risch`RischDE[%s, %s, x]}, "
        "If[Head[sol] === Risch`RischDE, $Unsolved, "
        "Together[Expand[D[sol, x] + (%s)*sol - (%s)]]]]",
        f, g, f, g);
    char* s = eval_fullform(buf);
    if (strcmp(s, "0") != 0)
        printf("FAIL: RischDE(f=%s, g=%s) -> %s (expected solved, residual 0)\n", f, g, s);
    ASSERT_MSG(strcmp(s, "0") == 0,
        "RischDE(f=%s, g=%s): expected solved with residual 0, got %s", f, g, s);
    free(s);
}

/* Risch`RischDE[f, g, x] must stay UNEVALUATED (no rational solution in scope). */
static void assert_rde_nosol(const char* f, const char* g) {
    char buf[2048];
    snprintf(buf, sizeof(buf), "Risch`RischDE[%s, %s, x]", f, g);
    char* s = eval_fullform(buf);
    bool unev = (strncmp(s, "Risch`RischDE[", 14) == 0);
    if (!unev) printf("FAIL: RischDE(f=%s, g=%s) -> %s (expected unevaluated)\n", f, g, s);
    ASSERT_MSG(unev, "RischDE(f=%s, g=%s): expected unevaluated, got %s", f, g, s);
    free(s);
}

/* Risch`SPDE[a,b,c,x,n] must return exactly `expected` (FullForm). */
static void assert_spde(const char* a, const char* b, const char* c, int n,
                        const char* expected) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Risch`SPDE[%s, %s, %s, x, %d]", a, b, c, n);
    char* s = eval_fullform(buf);
    if (strcmp(s, expected) != 0)
        printf("FAIL: SPDE[%s,%s,%s,%d] -> %s (expected %s)\n", a, b, c, n, s, expected);
    ASSERT_MSG(strcmp(s, expected) == 0,
        "SPDE[%s,%s,%s,%d]: expected %s, got %s", a, b, c, n, expected, s);
    free(s);
}

/* Risch`PolyRischDENoCancel[b, D[q]+b q, x, n] must return a poly solution whose
 * identity D[sol]+b sol == D[q]+b q holds exactly. */
static void assert_prde_q(const char* b, const char* q, int n) {
    char buf[3072];
    snprintf(buf, sizeof(buf),
        "With[{cc = Expand[D[%s, x] + (%s)*(%s)]}, "
        "With[{sol = Risch`PolyRischDENoCancel[%s, cc, x, %d]}, "
        "If[sol === $Failed, $Unsolved, Expand[D[sol, x] + (%s)*sol - cc]]]]",
        q, b, q, b, n, b);
    char* s = eval_fullform(buf);
    if (strcmp(s, "0") != 0)
        printf("FAIL: PolyRischDENoCancel(b=%s, q=%s) -> %s\n", b, q, s);
    ASSERT_MSG(strcmp(s, "0") == 0,
        "PolyRischDENoCancel(b=%s, q=%s): expected residual 0, got %s", b, q, s);
    free(s);
}

/* Risch`PolyRischDENoCancel[b,c,x,n] must return $Failed (no bounded solution). */
static void assert_prde_fail(const char* b, const char* c, int n) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Risch`PolyRischDENoCancel[%s, %s, x, %d]", b, c, n);
    char* s = eval_fullform(buf);
    if (strcmp(s, "$Failed") != 0)
        printf("FAIL: PolyRischDENoCancel[%s,%s,%d] -> %s (expected $Failed)\n", b, c, n, s);
    ASSERT_MSG(strcmp(s, "$Failed") == 0,
        "PolyRischDENoCancel[%s,%s,%d]: expected $Failed, got %s", b, c, n, s);
    free(s);
}

/* Integrate[g, x] must be elementary and diff back to g exactly. */
static void assert_integrates(const char* g) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "With[{ii = Integrate[%s, x]}, "
        "If[!FreeQ[ii, Integrate], $Unintegrated, "
        "Simplify[D[ii, x] - (%s)]]]", g, g);
    char* s = eval_fullform(buf);
    if (strcmp(s, "0") != 0)
        printf("FAIL: Integrate(%s) -> %s (expected diff-back 0)\n", g, s);
    ASSERT_MSG(strcmp(s, "0") == 0,
        "Integrate(%s): expected elementary, diff-back 0, got %s", g, s);
    free(s);
}

/* ------------------------------------------------------------------ */
/* Base field C(x) — regression on the pre-existing rde_base path.     */
/* ------------------------------------------------------------------ */
static void test_base_field(void) {
    assert_rde_y("1", "x");                 /* y'+y = 1+x  -> x            */
    assert_rde_y("1", "x^2 - 2 x + 2");     /* y'+y polynomial            */
    assert_rde_y("2 x", "1");               /* y'+2x y = 2x -> 1 (const)  */
    assert_rde_y("1/x", "x^2");             /* rational coefficient       */
    assert_rde_nosol("2 x", "1");           /* Erf: q'+2x q=1 no rational  */
    assert_rde_nosol("1", "1/x");           /* Ei: q'+q=1/x no rational    */
}

/* ------------------------------------------------------------------ */
/* Risch`SPDE — Rothstein's box (base field), direct against the algebra. */
/* ------------------------------------------------------------------ */
static void test_spde_box(void) {
    /* a = 1 (already reduced): returns (b/a, c/a, n, 1, 0). */
    assert_spde("1", "x", "x^2", 3, "List[x, Power[x, 2], 3, 1, 0]");
    /* a in k* nonunit: divides through. */
    assert_spde("2", "x", "4", 2, "List[Times[Rational[1, 2], x], 2, 2, 1, 0]");
    /* deg(a) > 0: one reduction consumes everything into beta (alpha 0). */
    assert_spde("x", "1", "Plus[x^2, x]", 4,
                "List[0, 0, 0, 0, Plus[Times[Rational[1, 2], x], Times[Rational[1, 3], Power[x, 2]]]]");
    assert_spde("Plus[x, 1]", "2", "3", 0, "List[0, 0, 0, 0, Rational[3, 2]]");
    /* n < 0 with c != 0: no solution. */
    assert_spde("1", "1", "1", -1, "$Failed");
}

/* ------------------------------------------------------------------ */
/* Risch`PolyRischDENoCancel — non-cancellation polynomial solve.       */
/* ------------------------------------------------------------------ */
static void test_polyrischde_nocancel(void) {
    assert_prde_q("1", "1", 3);             /* q'+q = 1 -> 1             */
    assert_prde_q("x", "x^2 + 1", 5);       /* q'+x q                    */
    assert_prde_q("1", "x^3 - 2 x", 6);     /* higher degree             */
    assert_prde_q("x^2", "3 x^4 + x", 8);   /* deg(b)=2                  */
    assert_prde_fail("x", "x^2 + 2 x", 5);  /* leading + constant clash  */
}

/* ------------------------------------------------------------------ */
/* Tower field — depth-1 primitive (Log[x]) top, NON-cancellation.      */
/* deg_tau(f) >= 1 so the leading terms never cancel.                   */
/* ------------------------------------------------------------------ */
static void test_tower_depth1_log(void) {
    /* Polynomial-in-tau solutions. */
    assert_rde_y("Log[x]/x", "Log[x]");                 /* y = tau       */
    assert_rde_y("Log[x]/x", "Log[x]^2");               /* y = tau^2     */
    assert_rde_y("Log[x]/x", "Log[x]^3 + 2 Log[x]");    /* y = tau^3+..  */
    assert_rde_y("Log[x]/x", "3 Log[x]^4 - Log[x]");    /* higher degree */
    /* Rational-in-tau solutions (a genuine tau-pole; normal-denom reduction). */
    assert_rde_y("2 Log[x]/x", "1/Log[x]");             /* y = 1/tau     */
    assert_rde_y("Log[x]/x", "1/(1 + Log[x])");         /* non-monomial denom */
    assert_rde_y("2 Log[x]/x", "Log[x] + 1/Log[x]");    /* mixed poly+pole    */
    /* x-dependent (k = C(x)) lower-field coefficients. */
    assert_rde_y("Log[x]/x", "x Log[x]^2");
    assert_rde_y("Log[x]/x", "Log[x]/(x + 1)");
}

/* ------------------------------------------------------------------ */
/* Tower field — depth-2 and depth-3 nested-Log towers.                 */
/* The derivation abstraction must generalise (k itself carries logs).  */
/* ------------------------------------------------------------------ */
static void test_tower_deep_log(void) {
    /* depth 2: tau = Log[Log[x]], k = C(x, Log[x]). */
    assert_rde_y("Log[Log[x]]/(x Log[x])", "Log[Log[x]]");
    assert_rde_y("Log[Log[x]]/(x Log[x])", "Log[Log[x]]^2");
    assert_rde_y("2 Log[Log[x]]/(x Log[x])", "1/Log[Log[x]]");
    assert_rde_y("Log[Log[x]]/(x Log[x])", "1/(1 + Log[Log[x]])");
    /* depth 3: tau = Log[Log[Log[x]]]. */
    assert_rde_y("Log[Log[Log[x]]]/(x Log[x] Log[Log[x]])", "Log[Log[Log[x]]]");
    assert_rde_y("Log[Log[Log[x]]]/(x Log[x] Log[Log[x]])", "Log[Log[Log[x]]]^2");
}

/* ------------------------------------------------------------------ */
/* Tower field — exponential (RT_EXP) top, NON-cancellation (Gap 1b).   */
/* tau hyperexponential (delta=1, deriv-preserving); the special part of */
/* the denominator (tau-power poles) is cleared by RdeSpecialDenomExp.   */
/* ------------------------------------------------------------------ */
static void test_tower_exp(void) {
    /* tau = E^x, eta = Dtau/tau = 1.  Non-cancellation: deg_tau(f) >= 1. */
    assert_rde_y("E^x", "1/(1 + E^x)");
    assert_rde_y("E^x", "E^x/(2 + E^x)");
    assert_rde_y("E^x", "1/(1 + E^x)^2");           /* repeated normal pole */
    assert_rde_y("2 E^x", "1/(3 + E^x)");
    assert_rde_y("E^x", "E^x");                     /* q_bar polynomial in tau */
    assert_rde_y("E^x", "E^(2 x)/(1 + E^x)");
    /* depth-2 exp-over-exp: tau = E^(E^x), eta = E^x — the derivation
     * abstraction must generalise (k = C(x, E^x) itself carries an exponential). */
    assert_rde_y("E^(E^x)", "1/(1 + E^(E^x))");
    assert_rde_y("E^(E^x)", "E^(E^x)/(2 + E^(E^x))");
}

/* ------------------------------------------------------------------ */
/* Tower field — the CANCELLATION case (Gap 1c): b in k* (deg_tau(b)=0),  */
/* PolyRischDECancel{Prim,Exp} (§6.6), which recurses one tower level     */
/* down (a lower-field Risch DE per coefficient).                        */
/* ------------------------------------------------------------------ */
static void test_tower_cancellation(void) {
    /* Primitive top, b = const in k: PolyRischDECancelPrim -> RischDE over C(x). */
    assert_rde_y("1", "Log[x]");
    assert_rde_y("1", "Log[x]^2");
    assert_rde_y("2", "3 Log[x] + 1");
    /* Exponential top, b in k: PolyRischDECancelExp (coefficient b + m*eta). */
    assert_rde_y("2", "E^(-x)");
    assert_rde_y("3", "E^(-2 x)");
    assert_rde_y("1", "E^(E^x)");                   /* q_bar poly; cancel -> nocancel */
    /* Deep: the per-coefficient lower Risch DE is itself an exponential
     * non-cancellation solve over C(x, E^x). */
    assert_rde_y("1", "1/(1 + E^(E^x))");
}

/* ------------------------------------------------------------------ */
/* Declines — never wrong, only silent.  Genuine no-solution + the      */
/* branches deferred to later increments / LimitedIntegrate (Gap 2).    */
/* ------------------------------------------------------------------ */
static void test_tower_declines(void) {
    /* Non-elementary RHS: no rational solution in the tower field. */
    assert_rde_nosol("Log[x]/x", "1/x");
    assert_rde_nosol("Log[x]/x", "Log[x]");
}

/* ------------------------------------------------------------------ */
/* Tower field — the ANTIDIFFERENTIATION branch (Gap 2): sp.b == 0, so   */
/* D h = c, solved by integrating c within the field K_m via the         */
/* integrator (the Bronstein RDE<->integrator mutual recursion).         */
/* ------------------------------------------------------------------ */
static void test_tower_antidiff(void) {
    /* Primitive top: D h = g integrates within k[tau] (poly antiderivative).
     * g is passed directly (it carries the Log kernel that defines the field). */
    assert_rde_solves("0", "Log[x]");               /* x Log[x] - x */
    assert_rde_solves("0", "Log[x]^2");
    assert_rde_solves("0", "x Log[x]");
    /* Exponential top: the Laurent-coefficient antiderivative, incl. a
     * non-constant coefficient (x E^x -> (x-1) E^x — the gate-Expand case). */
    assert_rde_solves("0", "E^x");
    assert_rde_solves("0", "x E^x");
    assert_rde_solves("0", "x^2 E^x");
    assert_rde_solves("0", "E^(2 x)");
    /* No rational (elementary-in-field) antiderivative -> decline: the integral
     * of 1/(1+Log[x]) is ~ExpIntegralEi (non-elementary). */
    assert_rde_nosol("0", "1/(1 + Log[x])");
}

/* ------------------------------------------------------------------ */
/* Integrator end-to-end — Laurent-coefficient RDEs now via rde_tower.  */
/* ------------------------------------------------------------------ */
static void test_integrator_endtoend(void) {
    /* E^(Log[x]^2): the i=1 Laurent coefficient RDE q'+(2Log/x)q = 2/x-1/(x Log^2)
     * -> q = 1/Log[x], solved literally by the tower RDE. */
    assert_integrates("(2/x - 1/(x Log[x]^2)) E^(Log[x]^2)");
    /* Same family, polynomial-in-tau coefficient. */
    assert_integrates("(2 Log[x]/x) E^(Log[x]^2)");
    /* Exp-over-exp: the E^(E^x) Laurent coefficient RDE q'+E^x q = p is solved by
     * the literal 1b exponential tower stack. */
    assert_integrates("D[E^(E^x)/(1 + E^x), x]");
    /* Antidifferentiation (Gap 2): the E^(E^x) i=0 Laurent coefficient integrates
     * within the field, x E^x + 1 -> x, giving x E^(E^x). */
    assert_integrates("E^(E^x) (x E^x + 1)");
    /* The non-elementary sibling still declines (guard against over-eager solve). */
    {
        char* s = eval_fullform("FreeQ[Integrate[E^(Log[x]^2), x], Integrate]");
        ASSERT_MSG(strcmp(s, "False") == 0,
            "Integrate[E^(Log[x]^2)] must stay unintegrated, got FreeQ=%s", s);
        free(s);
    }
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== Risch DE over the tower (Gap 1a, Bronstein Ch.6) ===\n");
    TEST(test_base_field);
    TEST(test_spde_box);
    TEST(test_polyrischde_nocancel);
    TEST(test_tower_depth1_log);
    TEST(test_tower_deep_log);
    TEST(test_tower_exp);
    TEST(test_tower_cancellation);
    TEST(test_tower_declines);
    TEST(test_tower_antidiff);
    TEST(test_integrator_endtoend);
    printf("All Risch DE tower tests passed.\n");

    symtab_clear();
    return 0;
}

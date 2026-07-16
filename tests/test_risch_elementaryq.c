/* test_risch_elementaryq.c
 *
 * P3 decision half — Risch`ElementaryIntegralQ[f, x], the Bronstein
 * elementary-integrability decision predicate.  Returns:
 *   True  — f has an elementary antiderivative (exhibited by the integrator);
 *   False — the recursive Risch decision procedure PROVES none exists, via
 *           §5.6 residue criterion (Thm 5.6.1(ii): a non-constant residue), a
 *           Ch.6 Risch DE with no rational solution, or the §5.8 Dc!=0 certificate;
 *   unevaluated — the verdict is outside the single-tower field scope.
 *
 * Sound by construction: a Boolean is emitted only behind an exact certificate,
 * never a bounded-ansatz "gave up".  The classic non-elementary integrands
 * E^x/x (ExpIntegralEi), E^(x^2) (Erf) and 1/Log[x] (LogIntegral) — which the
 * integrator answers with special functions, i.e. NON-elementary — must decide
 * False, while their elementary siblings decide True.
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

/* Evaluate Risch`ElementaryIntegralQ[f, x] and assert the FullForm result. */
static void assert_decision(const char* f, const char* expected) {
    char buf[512];
    snprintf(buf, sizeof(buf), "Risch`ElementaryIntegralQ[%s, x]", f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: ElementaryIntegralQ[%s] expected %s, got %s\n", f, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0,
        "ElementaryIntegralQ[%s]: expected %s, got %s", f, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

static void assert_true(const char* f)  { assert_decision(f, "True"); }
static void assert_false(const char* f) { assert_decision(f, "False"); }
/* Undecided: the call stays unevaluated (head unchanged). */
static void assert_undecided(const char* f) {
    char buf[512];
    snprintf(buf, sizeof(buf), "Risch`ElementaryIntegralQ[%s, x]", f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    bool undec = (strncmp(s, "Risch`ElementaryIntegralQ[", 26) == 0);
    if (!undec) printf("FAIL: ElementaryIntegralQ[%s] expected unevaluated, got %s\n", f, s);
    ASSERT_MSG(undec, "ElementaryIntegralQ[%s]: expected unevaluated, got %s", f, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* ---- False: PROVABLY non-elementary via the Risch DE certificate (Ch.6) ---- */
static void test_false_risch_de(void) {
    assert_false("E^x/x");                 /* ExpIntegralEi — RDE q'+q=1/x no soln  */
    assert_false("E^(x^2)");               /* Erf — RDE q'+2x q=1 no polynomial soln */
    assert_false("E^(-x^2)");              /* Erf sibling                            */
    assert_false("E^(E^x)/(1 + E^(E^x))"); /* non-elementary depth-2 tower           */
}

/* ---- False: PROVABLY non-elementary via the residue criterion (Thm 5.6.1) --- */
static void test_false_residue(void) {
    assert_false("1/Log[x]");              /* LogIntegral — residue x is non-constant */
}

/* ---- False: the polynomial-in-Log dilogarithm obstruction (Thm 5.6.1(ii)) ----
 * A (theta^{>=1})-level residual that integrates to an x-dependent ArcTan (from an
 * irreducible-quadratic log argument) has no elementary integral — the
 * dilogarithm obstruction, e.g. Integrate[Log[1+x^2]/(1+x^2)] = (dilog). */
static void test_false_poly_log_dilog(void) {
    assert_false("(3 x + 1) Log[x^2 + 1]^5"); /* reported In[7]: NOT elementary   */
    assert_false("Log[1 + x^2]/(1 + x^2)");   /* the depth-1 dilog obstruction    */
    assert_false("Log[1 + x^2]^2");           /* cross-term dilog (Catalan value) */
}

/* ---- False: trig / inverse-trig integrands, decided through the SAME Gaussian
 * exponential tower the integrator uses.  rt_decide_field exponentializes with
 * TrigToExp (+ rt_powers_to_exp) so a bare trig kernel becomes E^(I ...), then
 * reads the residue / Risch-DE certificate off the tower.  Each of these is a
 * classic special-function integral (non-elementary), so must decide False. */
static void test_false_trig_frontend(void) {
    assert_false("Sin[x]/x");        /* SinIntegral — E^(I x)/x RDE no soln     */
    assert_false("Sin[x]/x^2");      /* SinIntegral sibling                     */
    assert_false("Sin[x^2]");        /* FresnelS — E^(I x^2) Erf-type no soln   */
    assert_false("Cos[x^2]");        /* FresnelC sibling                        */
    assert_false("Cos[E^x]");        /* CosIntegral of E^x (depth-2 tower)      */
    assert_false("Cos[x Log[x]]");   /* trig of a logarithm: non-elementary     */
}

/* ---- True: an elementary antiderivative is exhibited ----------------------- */
static void test_true_elementary(void) {
    assert_true("E^x");
    assert_true("x E^x");
    assert_true("E^x/(1 + E^x)");          /* Log[1+E^x]        */
    assert_true("1/(x Log[x])");           /* Log[Log[x]]       */
    assert_true("Log[x]/x");               /* Log[x]^2/2        */
    assert_true("1/(1 + x^2)");            /* ArcTan[x] (rational base case) */
    assert_true("Tan[x]");                 /* -Log[Cos[x]] (hypertangent)    */
    assert_true("1/(x (Log[x]^2 + 1))");   /* elementary complex-residue ArcTan */
    /* Polynomial * Log: degree-1 in Log is ALWAYS elementary (one IBP -> rational
     * integral), even with an irreducible-quadratic argument whose bottom-level
     * ArcTan is a legitimate part of the answer. */
    assert_true("(x^5 - 1) Log[x^2 - x + 1]"); /* reported In[12] */
    assert_true("x Log[x^2 + 1]");
    assert_true("Log[x]^3");               /* linear argument: no dilog obstruction */
    /* Across the dispatch: fractional log/exp, nested tower, mixed tower. */
    assert_true("E^x/(E^(2 x) + 1)");      /* ArcTan[E^x] (algebraic residues)   */
    assert_true("1/(x Log[x] Log[Log[x]])"); /* Log[Log[Log[x]]] (nested log tower) */
    assert_true("1/(3 + Tan[x]^2)");       /* real hypertangent (§5.10)          */
    /* Trig integrands with an elementary antiderivative, exhibited before the
     * field decision is ever consulted (True dominates). */
    assert_true("Sin[x]");                 /* -Cos[x]                            */
    assert_true("Sec[x]^2");               /* Tan[x]                             */
    assert_true("Sin[Log[x]]");            /* trig of a log, now integrable      */
    assert_true("E^x Sin[x]");             /* damped oscillation                 */
}

/* ---- Undecided: outside the transcendental-tower field scope --------------- */
static void test_undecided_out_of_scope(void) {
    assert_undecided("1/Sqrt[1 + x^3]");   /* algebraic (elliptic) — no log/exp tower */
    /* Algebraic function of x: an exp/trig of a radical is elementary via an
     * ALGEBRAIC substitution the transcendental algorithm does not perform, so
     * the field decision must decline to UNKNOWN, never a spurious False.  (These
     * would blow the pure-transcendental scope if exponentialized naively.) */
    assert_undecided("E^Sqrt[x]");         /* = 2 E^Sqrt[x] (Sqrt[x] - 1): elementary */
    assert_undecided("Cos[Sqrt[x]]");      /* elementary via u = Sqrt[x]              */
    assert_undecided("Sin[x^(1/3)]");      /* elementary via u = x^(1/3)              */
    assert_undecided("Sqrt[Sin[x]]");      /* algebraic over the trig kernel          */
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== Risch`ElementaryIntegralQ decision predicate (P3) ===\n");
    TEST(test_false_risch_de);
    TEST(test_false_residue);
    TEST(test_false_poly_log_dilog);
    TEST(test_false_trig_frontend);
    TEST(test_true_elementary);
    TEST(test_undecided_out_of_scope);
    printf("All ElementaryIntegralQ tests passed.\n");

    symtab_clear();
    return 0;
}

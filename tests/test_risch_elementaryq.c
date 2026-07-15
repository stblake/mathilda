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
}

/* ---- Undecided: outside the transcendental-tower field scope --------------- */
static void test_undecided_out_of_scope(void) {
    assert_undecided("1/Sqrt[1 + x^3]");   /* algebraic (elliptic) — no log/exp tower */
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== Risch`ElementaryIntegralQ decision predicate (P3) ===\n");
    TEST(test_false_risch_de);
    TEST(test_false_residue);
    TEST(test_false_poly_log_dilog);
    TEST(test_true_elementary);
    TEST(test_undecided_out_of_scope);
    printf("All ElementaryIntegralQ tests passed.\n");

    symtab_clear();
    return 0;
}

/* test_risch_residue_split.c
 *
 * Phase-2 tests for the A2-2 PARTIAL log part (Bronstein, Symbolic Integration I,
 * Thm 5.6.1 / the residue-criterion κ_D split).  When the Rothstein-Trager
 * residue resultant of a proper rational function over a single transcendental
 * monomial mixes constant and non-constant residues (r = r_s · r_n), the
 * integrator now returns the elementary constant-residue logs from r_s PLUS an
 * unevaluated Integrate[remainder, x] for the non-constant part r_n, instead of
 * declining the whole integral.
 *
 * Also covers the Fundamental Theorem of Calculus derivative rule
 * D[Integrate[f, x], x] = f, which (a) is independently useful and (b) makes the
 * partial antiderivative verify by differentiation: D[logs + Integrate[rem,x]] =
 * D[logs] + rem = integrand.
 *
 * Correctness is asserted by the universal predicate Simplify[D[∫f] − f] === 0
 * (robust to surface form), plus structural checks (a partial result contains an
 * Integrate[] term; a fully-elementary result does not; a fully non-elementary
 * input declines).
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

/* Assert the FullForm of `input` equals `expected`. */
static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert Simplify[D[head[f, x], x] − (f)] === 0. */
static void assert_diff_zero_via(const char* head, const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "Simplify[D[%s[%s, x], x] - (%s)]", head, f, f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, "0") != 0)
        printf("FAIL diff-back (%s): %s\n  -> %s\n", head, f, s);
    ASSERT_MSG(strcmp(s, "0") == 0, "diff-back (%s) nonzero for %s: got %s", head, f, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert FreeQ[<head[f,x] evaluated>, Integrate] === expected_free.
 * expected_free=false pins a PARTIAL result (carries Integrate[remainder,x]);
 * expected_free=true pins a fully-elementary result (no unintegrated part). */
static void assert_integrate_presence(const char* head, const char* f, bool expected_free) {
    char buf[1400];
    snprintf(buf, sizeof(buf), "FreeQ[%s[%s, x], Integrate]", head, f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    const char* want = expected_free ? "True" : "False";
    if (strcmp(s, want) != 0)
        printf("FAIL presence: %s[%s]\n  FreeQ[.,Integrate] expected %s got %s\n",
               head, f, want, s);
    ASSERT_MSG(strcmp(s, want) == 0, "%s[%s]: FreeQ expected %s got %s", head, f, want, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert head[f,x] stays unevaluated (declines with the same head). */
static void assert_head_unevaluated(const char* head, const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf), "%s[%s, x]", head, f);
    Expr* parsed = parse_expression(buf);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol.name, head) == 0,
        "expected unevaluated %s[...] for: %s", head, f);
    expr_free(result);
}

/* ---- Fundamental theorem of calculus derivative rule --------------------- */
static void test_ftc_derivative(void) {
    run_test("D[Integrate[f[x], x], x]", "f[x]");
    run_test("D[Integrate[Sin[x^2], x], x]", "Sin[Power[x, 2]]");
    run_test("D[Integrate[a x^3 + b, x], x]", "Plus[b, Times[a, Power[x, 3]]]");
    /* Different integration variable / definite forms are NOT the FTC: f[y] is
     * free of x, so its x-derivative is 0. */
    run_test("D[Integrate[f[y], y], x]", "0");
}

/* ---- A2-2 partial log part (mixed residues) ------------------------------ */

/* The canonical mixed case over a single logarithm: 1/(x Log[x]) has a constant
 * residue (→ Log[Log[x]]); 1/(x (Log[x]^2 − x)) has non-constant residues
 * (Log[x]^2 = x depends on x) → non-elementary.  Combined, the residue resultant
 * is mixed and the integrator returns the elementary log plus Integrate[rem]. */
#define MIXED_LOG "1/(x Log[x]) + 1/(x (Log[x]^2 - x))"

static void test_partial_log_mixed(void) {
    /* Correct by differentiation (FTC closes the Integrate[rem] term). */
    assert_diff_zero_via("Integrate", MIXED_LOG);
    assert_diff_zero_via("Integrate`RischTranscendental", MIXED_LOG);
    /* It is a genuine PARTIAL result: an Integrate[remainder, x] survives. */
    assert_integrate_presence("Integrate`RischTranscendental", MIXED_LOG, false);
    /* The elementary constant-residue part is exactly Log[Log[x]]. */
    run_test("Simplify[First[Integrate`RischTranscendental[" MIXED_LOG ", x]] - Log[Log[x]]]", "0");
}

/* Scope boundary: the EXPONENTIAL (hyperexponential) monomial does not yet
 * surface a partial log part — the §5.9 coupled Laurent reconciliation of the
 * exponential proper-part pipeline cannot carry a proper simple remainder, so a
 * mixed exp resultant declines (refuses the partial) rather than emitting a
 * wrong or half-reconciled form.  Pinned so a future extension flips it
 * deliberately.  (Via linearity the user-facing Integrate still splits an
 * explicit SUM of an elementary and a non-elementary exp term.) */
static void test_exp_partial_out_of_scope(void) {
    assert_head_unevaluated("Integrate`RischTranscendental",
                            "(1 + x + 2 E^x)/((1 + E^x) (x + E^x))");
}

/* ---- Non-regression: all-constant residues stay fully elementary --------- */
static void test_all_constant_unchanged(void) {
    /* Fully elementary — NO Integrate[] term must appear. */
    assert_integrate_presence("Integrate`RischTranscendental", "1/(x (Log[x]^2 + Log[x]))", true);
    assert_integrate_presence("Integrate`RischTranscendental", "1/(x Log[x]^2)", true);
    assert_diff_zero_via("Integrate`RischTranscendental", "1/(x (Log[x]^2 + Log[x]))");
    /* Fully elementary complex-residue case (ArcTan) also stays whole. */
    assert_integrate_presence("Integrate`RischTranscendental", "1/(x (Log[x]^2 + 1))", true);
}

/* ---- Non-regression: all-non-constant residues decline (no spurious 0+∫,
 *      no infinite Integrate[remainder] regress) ---------------------------- */
static void test_all_nonconstant_declines(void) {
    assert_head_unevaluated("Integrate`RischTranscendental", "1/(x (Log[x]^2 - x))");
    assert_head_unevaluated("Integrate`RischTranscendental", "1/(x (Log[x]^3 - x))");
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== residue-criterion partial log part (Phase 2) ===\n");
    test_ftc_derivative();
    test_partial_log_mixed();
    test_exp_partial_out_of_scope();
    test_all_constant_unchanged();
    test_all_nonconstant_declines();
    printf("All residue-split tests passed.\n");

    symtab_clear();
    return 0;
}

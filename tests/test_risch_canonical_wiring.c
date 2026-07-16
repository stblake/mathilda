/* test_risch_canonical_wiring.c
 *
 * Phase-1 wiring tests for the canonical-representation SPINE of the recursive
 * transcendental Risch integrator.  As of the P0 wiring, the RT_LOG branch of
 * rt_field_integrate (integrate_risch_transcendental.c) dispatches on the
 * Bronstein §3.5 canonical split F = f_p + f_s + f_n (via
 * risch_canonical_representation) instead of the former ad-hoc
 * PolynomialQuotient/Remainder gate.  This suite pins two things:
 *
 *   (A) the split's PRIMITIVE-monomial invariant the spine relies on — a
 *       logarithmic (primitive) monomial has no non-trivial special polynomial,
 *       so f_s ≡ 0 — plus exact reconstruction f_p + f_s + f_n = F.  Verified
 *       through the Risch`CanonicalRepresentation builtin.  A contrasting
 *       exponential case shows f_s ≠ 0 (the split really does isolate the
 *       special part), also reconstructing exactly.
 *
 *   (B) that primitive-tower integrands still integrate correctly THROUGH the
 *       new spine — pure polynomial part, proper normal part, mixed, rational
 *       lower-field coefficients, and depth-2/3 log towers — asserted by the
 *       universal diff-back predicate Simplify[D[∫f] − f] === 0 (surface-form
 *       independent).  These are behaviour-preservation guards: the closed forms
 *       must remain correct after the spine replaced the ad-hoc split.
 *
 * Equality of split parts is checked with Expand[Together[a − b]] == 0 (robust
 * to cosmetic ordering), matching test_risch_canonical.c.
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

/* Assert that the rational expression `diff` is identically zero. */
static void run_zero(const char* diff) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Expand[Together[%s]]", diff);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, "0") != 0)
        printf("FAIL: %s\n  expected: 0\n  got:      %s\n", diff, s);
    ASSERT_MSG(strcmp(s, "0") == 0, "%s: expected 0, got %s", diff, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert that the FullForm of `input` equals `expected`. */
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

/* Assert Simplify[D[∫[f], x] − (f)] === 0 for a chosen integrator head.  `head`
 * is either "Integrate`RischTranscendental" (drives the recursive tower
 * integrator — and thus the canonical spine — directly) or "Integrate" (the
 * user surface, which may reach the spine via linearity/normalization). */
static void assert_diff_zero_via(const char* head, const char* f) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "Simplify[D[%s[%s, x], x] - (%s)]", head, f, f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, "0") != 0)
        printf("FAIL diff-back (%s): %s\n  D[∫]-f simplified to: %s\n", head, f, s);
    ASSERT_MSG(strcmp(s, "0") == 0, "diff-back nonzero (%s) for %s: got %s", head, f, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Direct package-head diff-back (recursive tower integrator / canonical spine). */
static void assert_rm_diff_zero(const char* f) {
    assert_diff_zero_via("Integrate`RischTranscendental", f);
}

/* ---- (A) canonical split: primitive f_s ≡ 0 + exact reconstruction -------- */

/* Single logarithmic monomial t = Log[x]:  deriv = {x -> 1, t -> 1/x}. */
#define LOGDERIV "{x -> 1, t -> 1/x}"
#define LOGCR(F) "Risch`CanonicalRepresentation[" F ", t, " LOGDERIV "]"

static void test_primitive_special_part_vanishes(void) {
    /* For every primitive-monomial input the special part (2nd component) is 0. */
    run_test("Part[" LOGCR("t^2 + t + 1/x") ", 2]", "0");
    run_test("Part[" LOGCR("1/t^2") ", 2]", "0");
    run_test("Part[" LOGCR("(t + 1)/(t^2 (t + 1/x))") ", 2]", "0");
    run_test("Part[" LOGCR("(x t^3 + t + 1)/(t^2 + 1)") ", 2]", "0");

    /* Reconstruction f_p + f_s + f_n = F for the same inputs. */
    run_zero("Total[" LOGCR("t^2 + t + 1/x") "] - (t^2 + t + 1/x)");
    run_zero("Total[" LOGCR("1/t^2") "] - (1/t^2)");
    run_zero("Total[" LOGCR("(t + 1)/(t^2 (t + 1/x))") "] - ((t + 1)/(t^2 (t + 1/x)))");
    run_zero("Total[" LOGCR("(x t^3 + t + 1)/(t^2 + 1)") "] - ((x t^3 + t + 1)/(t^2 + 1))");
}

static void test_exponential_special_part_isolated(void) {
    /* Contrast: an exponential monomial t = E^x (deriv {x -> 1, t -> t}) DOES
     * carry a special polynomial (t itself), so the special part is nonzero and
     * captures the t = 0 pole; reconstruction still exact. */
    #define EXPCR(F) "Risch`CanonicalRepresentation[" F ", t, {x -> 1, t -> t}]"
    /* 1/(t^2 (1 + t)) = -1/t + 1/t^2 + 1/(1+t): the t = 0 poles (denominator a
     * power of t) are SPECIAL for an exponential, 1/(1+t) is normal. */
    run_zero("Part[" EXPCR("1/(t^2 (1 + t))") ", 2] - (1/t^2 - 1/t)"); /* f_s special */
    run_zero("Part[" EXPCR("1/(t^2 (1 + t))") ", 3] - (1/(1 + t))");    /* f_n = 1/(1+t) */
    run_zero("Total[" EXPCR("1/(t^2 (1 + t))") "] - (1/(t^2 (1 + t)))");
    #undef EXPCR
}

/* ---- (B) integration through the spine (diff-back == 0) -------------------- */

static void test_spine_polynomial_part(void) {
    /* Pure polynomial part f_p in t = Log[x] (§5.8 primitive-polynomial). */
    assert_rm_diff_zero("Log[x]");
    assert_rm_diff_zero("Log[x]^2");
    assert_rm_diff_zero("Log[x]^3");
    /* Rational lower-field coefficient on the polynomial part. */
    assert_rm_diff_zero("Log[x]/x^2");
    assert_rm_diff_zero("(2 + 3 Log[x])/x");
}

static void test_spine_normal_part(void) {
    /* Proper normal part f_n: Hermite (repeated pole) + residue log part. */
    assert_rm_diff_zero("1/(x Log[x]^2)");            /* repeated pole -> -1/Log[x] */
    assert_rm_diff_zero("1/(x Log[x])");              /* residue log  -> Log[Log[x]] */
    assert_rm_diff_zero("1/(x (Log[x]^2 + Log[x]))"); /* split residue -> ArcTanh   */
    assert_rm_diff_zero("1/(x Log[x]^3)");            /* higher repeated pole       */
}

static void test_spine_mixed_part(void) {
    /* Proper combined normal part: a repeated pole (Hermite) AND residue factors
     * in ONE denominator — the canonical f_n exercised end to end through the
     * recursive tower head. */
    assert_rm_diff_zero("1/(x Log[x]^2 (Log[x] + 1))");        /* Hermite + residue    */
    assert_rm_diff_zero("1/(x Log[x] (Log[x] + 1))");          /* two residues         */
    assert_rm_diff_zero("(Log[x] + 2)/(x Log[x]^2 (Log[x] + 3))"); /* residue num + pole */
}

static void test_spine_improper_via_integrate(void) {
    /* Improper polynomial-plus-proper single-log sums.  The recursive tower head
     * declines the COMBINED improper rational (a pre-existing limitation of the
     * package head — confirmed identical before the canonical spine landed); the
     * user-facing Integrate reaches the spine per summand via linearity and
     * closes them.  These guard that the spine keeps producing correct closed
     * forms at the user surface. */
    assert_diff_zero_via("Integrate", "Log[x] + 1/(x Log[x]^2)");
    assert_diff_zero_via("Integrate", "Log[x]^2 + 1/(x Log[x])");
    assert_diff_zero_via("Integrate", "(x Log[x]^3 + 1)/(x Log[x]^2)");
}

static void test_spine_nested_towers(void) {
    /* Depth-2 / depth-3 primitive (log) towers exercised through the recursive
     * field integrator (the spine at each level). */
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]])");           /* -> Log[Log[Log[x]]] */
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]]^2)");
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]] Log[Log[Log[x]]])");
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== canonical-representation wiring (Phase 1) ===\n");
    test_primitive_special_part_vanishes();
    test_exponential_special_part_isolated();
    test_spine_polynomial_part();
    test_spine_normal_part();
    test_spine_mixed_part();
    test_spine_improper_via_integrate();
    test_spine_nested_towers();
    printf("All canonical-wiring tests passed.\n");

    symtab_clear();
    return 0;
}

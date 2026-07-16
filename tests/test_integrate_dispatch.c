/* test_integrate_dispatch.c
 *
 * Smoke tests for the three-stage Integrate cascade and the Method
 * option added 2026-05-15.  Covers:
 *  - Cascade routing: rational integrand goes through BronsteinRational, an
 *    elementary one through RischNorman, an explicit CRCTable call
 *    survives without infinite-looping on the formerly-divergent
 *    inputs (Formula 49 family).
 *  - Method option: strict passthrough; unknown method bubbles back
 *    with Integrate::method.
 *  - Termination: pathological inputs that would have looped under
 *    the pre-2026-05-15 CRC table all return within a small budget.
 *
 * What we deliberately DO NOT test here:
 *  - Numerical correctness of every CRC formula.  Most rules don't
 *    fire today because Mathilda's matcher does not fully support
 *    /;-guarded multi-arg patterns (a separate work item); the
 *    cascade is in place for when that lands.
 *  - End-to-end Risch correctness (covered by intrischnorman_tests).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static double elapsed_seconds(const char* input) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    if (result) expr_free(result);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    return (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
}

static void test_cascade_rational(void) {
    /* Default cascade — first stage closes. */
    assert_eval_eq("Integrate[x^2, x]", "1/3 x^3", 0);
    assert_eval_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]", 0);
}

static void test_cascade_risch(void) {
    /* Default cascade — Risch closes once Rational gives up. */
    Expr* parsed = parse_expression("Integrate[Sin[x], x]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    /* Any non-NULL, non-Integrate-headed result is success — the
     * specific form differs across Risch-Norman builds. */
    ASSERT(result != NULL);
    ASSERT(!(result->type == EXPR_FUNCTION
             && result->data.function.head
             && result->data.function.head->type == EXPR_SYMBOL
             && strcmp(result->data.function.head->data.symbol.name, "Integrate") == 0));
    expr_free(result);
}

static void test_method_strict_rational(void) {
    /* Method -> "BronsteinRational" closes the polynomial case. */
    assert_eval_eq("Integrate[x^3, x, Method -> \"BronsteinRational\"]", "1/4 x^4", 0);

    /* Method -> "BronsteinRational" on a non-rational integrand bubbles back
     * with no Risch / CRC fallback. */
    Expr* parsed = parse_expression("Integrate[Sin[x], x, Method -> \"BronsteinRational\"]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL
        && result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol.name, "Integrate") == 0);
    expr_free(result);
}

static void test_method_invalid(void) {
    /* Method -> "Bogus" emits Integrate::method (to stderr) and bubbles
     * back unevaluated with no crash. */
    Expr* parsed = parse_expression("Integrate[x, y, Method -> \"Bogus\"]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol.name, "Integrate") == 0);
    expr_free(result);
}

static void test_crctable_termination(void) {
    /* These integrands match the LHS of Formula 49 (line 119 in the
     * CRC .m file) with non-positive integer exponents — the
     * pre-2026-05-15 rule would have spun forever.  With the
     * IntegerQ + bound guard added in this change, each call must
     * return promptly. */
    double t1 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-3), x, Method -> \"CRCTable\"]");
    double t2 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-5), x, Method -> \"CRCTable\"]");
    double t3 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-7), x, Method -> \"CRCTable\"]");

    /* A correct termination guard returns within a few hundred
     * milliseconds even at REPL warm-up.  Allow 5s wall-clock as a
     * very generous ceiling — Formula 49's old form would lock the
     * REPL indefinitely. */
    ASSERT_MSG(t1 < 5.0, "Formula 49 (n=-3) took %.2fs (>5s) — likely diverging.", t1);
    ASSERT_MSG(t2 < 5.0, "Formula 49 (n=-5) took %.2fs (>5s) — likely diverging.", t2);
    ASSERT_MSG(t3 < 5.0, "Formula 49 (n=-7) took %.2fs (>5s) — likely diverging.", t3);
}

static void test_crctable_reciprocal_sqrt_branch(void) {
    /* Formulas 398-401: 1/Sqrt[1 +- Cos/Sin[x]].  The classical
     * Sqrt[2] Log[Tan[...]] primitives are real on only every other
     * inter-pole interval of the integrand and sign-wrong elsewhere (the
     * integrand carries an absolute value).  The branch-correct forms keep
     * the literal radical so D[result, x] - integrand vanishes on every
     * pole-bounded interval, including the one the old form got wrong.
     *
     * Each check picks a point on an interval where the OLD form failed:
     *   1/Sqrt[1-Sin] at x=0.7 (old form -> complex; this was the corpus
     *   DIFF-NONZERO regression), the Cos/+Sin siblings at x=8 (second
     *   inter-pole interval).  We assert Abs[diff] < 1e-6 numerically. */
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[1 - Sin[x]], x, Method -> \"CRCTable\"], x]"
        " - 1/Sqrt[1 - Sin[x]]) /. x -> 0.7]] < 0.000001", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[1 - Sin[x]], x, Method -> \"CRCTable\"], x]"
        " - 1/Sqrt[1 - Sin[x]]) /. x -> 2.6]] < 0.000001", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[1 - Cos[x]], x, Method -> \"CRCTable\"], x]"
        " - 1/Sqrt[1 - Cos[x]]) /. x -> 8.0]] < 0.000001", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[1 + Cos[x]], x, Method -> \"CRCTable\"], x]"
        " - 1/Sqrt[1 + Cos[x]]) /. x -> 8.0]] < 0.000001", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[1 + Sin[x]], x, Method -> \"CRCTable\"], x]"
        " - 1/Sqrt[1 + Sin[x]]) /. x -> 8.0]] < 0.000001", "True", 0);
}

static void test_arctanh_real_branch(void) {
    /* Radical antiderivatives that a substitution integrator would express
     * with ArcTanh[g], |g| > 1 on the real axis, must be repaired to the
     * derivative-identical, real-valued ArcCoth[g] (intsimp_finalize).  Before
     * the repair, Integrate[Sqrt[x + Sqrt[x]], x] evaluated at x = 2 returned
     * 2.68952 - I Pi/8 because ArcTanh[Sqrt[1 + Sqrt[x]]/x^(1/4)] sits on its
     * branch cut for every x > 0.
     *
     * We assert (a) the antiderivative is real on x > 0, and (b) the
     * derivative still matches the integrand.  The derivative check guards
     * correctness independent of the surd form; the imaginary-part check
     * guards the real-branch selection. */
    assert_eval_eq(
        "Abs[Im[N[Integrate[Sqrt[x + Sqrt[x]], x] /. x -> 2]]] < 0.000001",
        "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[Sqrt[x + Sqrt[x]], x], x] - Sqrt[x + Sqrt[x]])"
        " /. x -> 2.3]] < 0.000001", "True", 0);
    /* Sibling case: 1/(x Sqrt[x + 1]) -> -2 ArcCoth[Sqrt[1 + x]] (real). */
    assert_eval_eq(
        "Abs[Im[N[Integrate[1/(x Sqrt[x + 1]), x] /. x -> 3]]] < 0.000001",
        "True", 0);
    /* Non-radical ArcTanh must be left intact (gate is radical-only). */
    assert_eval_eq("Integrate[1/(1 - x^2), x]", "ArcTanh[x]", 0);
}

static void test_root_kernel_keeps_radicals(void) {
    /* Derivative-divides on a unit-root kernel (u = Sqrt[x]) must close the
     * substitution with x -> u^2 and PowerExpand so the antiderivative stays in
     * the integrand's OWN radical Sqrt[x + Sqrt[x]] rather than falling through
     * to the Eliminate/Solve search's x^(1/4) / Sqrt[1 + Sqrt[x]] form.
     *
     * The decisive, form-independent check is that the derivative simplifies
     * back to exactly the integrand (canonical printed form Sqrt[Sqrt[x] + x]);
     * we also pin a numeric derivative-back for the x^2-weighted sibling, which
     * exercises the (u^2)^(5/2) -> u^5 PowerExpand path. */
    assert_eval_eq(
        "D[Integrate[Sqrt[x + Sqrt[x]], x], x] // Simplify",
        "Sqrt[Sqrt[x] + x]", 0);
    assert_eval_eq(
        "D[Integrate[x^2 Sqrt[x + Sqrt[x]], x], x] // Simplify",
        "x^2 Sqrt[Sqrt[x] + x]", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[x^2 Sqrt[x + Sqrt[x]], x], x]"
        " - x^2 Sqrt[x + Sqrt[x]]) /. x -> 1.7]] < 0.000001", "True", 0);
    /* The result must not carry the foreign x^(1/4) generator. */
    assert_eval_eq(
        "FreeQ[Integrate[Sqrt[x + Sqrt[x]], x], x^(1/4)]", "True", 0);
}

static void test_crctable_secant_cosecant_powers(void) {
    /* Sec[a x]^n and Csc[a x]^n reduction rules (added 2026-07-16).  The CRC
     * table only carried the reciprocal-power forms 1/Cos[a x]^m and
     * 1/Sin[a x]^m, which never fire on the Sec / Csc heads — so
     * Integrate[Csc[x]^4, x, Method -> "CRCTable"] returned unevaluated even
     * though the hyperbolic Csch[x]^4 already reduced.
     *
     * The Csc recursion carries a + on its integral term (Csc^2 = 1 + Cot^2),
     * unlike the Csch analogue's - (Csch^2 = Coth^2 - 1); getting that sign
     * wrong leaves an extra -4/3 Csc^2 in the derivative.  We assert the
     * antiderivative differentiates back to the integrand for both even and
     * odd powers and for a nontrivial linear argument a x. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sec[x]^4, x, Method -> \"CRCTable\"], x] - Sec[x]^4]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Sec[x]^5, x, Method -> \"CRCTable\"], x] - Sec[x]^5]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Sec[2 x]^4, x, Method -> \"CRCTable\"], x] - Sec[2 x]^4]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Csc[x]^4, x, Method -> \"CRCTable\"], x] - Csc[x]^4]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Csc[x]^7, x, Method -> \"CRCTable\"], x] - Csc[x]^7]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Csc[3 x]^5, x, Method -> \"CRCTable\"], x] - Csc[3 x]^5]",
        "0", 0);
}

static void test_crctable_trig_power_family(void) {
    /* Regression guard for the full trig / hyperbolic integer-power reduction
     * family (six circular + six hyperbolic).  Each must differentiate back to
     * the integrand across even and odd powers — a single {0,0,0,0,0} list per
     * head locks the whole family against future rule removals. */
    const char* heads[] = {"Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
                           "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch"};
    char buf[256];
    for (size_t i = 0; i < sizeof(heads) / sizeof(heads[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Table[Simplify[D[Integrate[%s[x]^n, x, Method -> \"CRCTable\"], x] "
            "- %s[x]^n], {n, 2, 6}]", heads[i], heads[i]);
        assert_eval_eq(buf, "{0, 0, 0, 0, 0}", 0);
    }
}

static void test_crctable_hyperbolic_products(void) {
    /* Sinh^m Cosh^n products (added 2026-07-16).  The CRC table's Formula 547
     * reduction bottomed out on single-Cosh / single-Sinh factors that had no
     * base rule, and only bound a bare argument x — so Sinh[x]^4 Cosh[x]^3 and
     * every a != 1 case returned unevaluated.  Verified by diff-back for even
     * and odd Cosh powers and for a nontrivial argument. */
    const char* cases[] = {
        "Sinh[x]^4 Cosh[x]^3", "Sinh[x]^3 Cosh[x]^4", "Sinh[x] Cosh[x]",
        "Sinh[x]^5 Cosh[x]", "Sinh[x] Cosh[x]^5", "Sinh[2 x]^4 Cosh[2 x]^3",
    };
    char buf[256];
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[%s, x, Method -> \"CRCTable\"], x] - (%s)]",
            cases[i], cases[i]);
        assert_eval_eq(buf, "0", 0);
    }
}

static void test_crctable_reciprocal_head_products(void) {
    /* Reciprocal-head products.  Mathilda canonicalises 1/Sin -> Csc,
     * 1/Cos -> Sec, 1/Sinh -> Csch, 1/Cosh -> Sech, so the CRC
     * 1/(Sin^m Cos^n) and 1/(Sinh^m Cosh^n) rules never matched an evaluated
     * integrand (Csc[x]^4 Sec[x]^3 failed in BOTH circular and hyperbolic).
     * The head-keyed reductions close them; verified by diff-back. */
    const char* cases[] = {
        "Csc[x]^4 Sec[x]^3", "Csc[x] Sec[x]", "Csc[x] Sec[x]^4",
        "Csc[2 x]^3 Sec[2 x]^2",
        "Csch[x]^4 Sech[x]^3", "Csch[x] Sech[x]", "Csch[x] Sech[x]^4",
        "Csch[2 x]^4 Sech[2 x]^3",
    };
    char buf[256];
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[%s, x, Method -> \"CRCTable\"], x] - (%s)]",
            cases[i], cases[i]);
        assert_eval_eq(buf, "0", 0);
    }
}

static void test_crctable_exp_hyperbolic(void) {
    /* E^(a x) Sinh^n / Cosh^n (added 2026-07-16).  Unlike the circular
     * E^(a x) Sin^n / Cos^n reductions (denominator a^2 + n^2 b^2, never zero),
     * the hyperbolic denominator a^2 - n^2 b^2 vanishes at the resonance
     * a = n b — exactly the E^x Cosh[x] case (a = b = 1), which carries a
     * secular x/2 term.  The Cosh recursion also flips sign relative to Sinh
     * (cosh^2 = 1 + sinh^2 vs sinh^2 = cosh^2 - 1).  Diff-back over n = 1..5,
     * covering the resonant a = b base and non-resonant a != b. */
    char buf[256];
    for (int n = 1; n <= 5; n++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[Exp[x] Cosh[x]^%d, x, Method -> \"CRCTable\"], x] "
            "- Exp[x] Cosh[x]^%d]", n, n);
        assert_eval_eq(buf, "0", 0);
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[Exp[x] Sinh[x]^%d, x, Method -> \"CRCTable\"], x] "
            "- Exp[x] Sinh[x]^%d]", n, n);
        assert_eval_eq(buf, "0", 0);
    }
    /* Non-resonant argument (a=2, b=3): denominator never hits zero for n<=5. */
    assert_eval_eq(
        "Simplify[D[Integrate[Exp[2 x] Cosh[3 x]^3, x, Method -> \"CRCTable\"], x] "
        "- Exp[2 x] Cosh[3 x]^3]", "0", 0);
}

static void test_crctable_hyperbolic_argument_generality(void) {
    /* The single-function hyperbolic bases were generalised from bare x to
     * argument a x; without this every a != 1 reduction that bottomed out on a
     * base (e.g. Sech[2 x]^3 -> ... -> Integrate[Sech[2 x]]) failed to close. */
    const char* heads[] = {"Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch"};
    char buf[256];
    for (size_t i = 0; i < sizeof(heads) / sizeof(heads[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[%s[2 x]^3, x, Method -> \"CRCTable\"], x] "
            "- %s[2 x]^3]", heads[i], heads[i]);
        assert_eval_eq(buf, "0", 0);
    }
}

static void test_crctable_mixed_tan_sec_families(void) {
    /* Tan^m Sec^n / Cot^m Csc^n and hyperbolic Tanh^m Sech^n / Coth^m Csch^n
     * (added 2026-07-16).  These are the canonical forms Mathilda folds the
     * Formula 307/308 quotients into (Cos^2/Sin^3 -> Cot^2 Csc, Sin^2/Cos^3 ->
     * Sec Tan^2), so the quotient rules never matched.  The n = 1 (single
     * Sec/Csc) case — including Sec Tan and Csc Cot, which the circular CRC
     * tables omit — is covered via an optional exponent.  Full 4x4 diff-back. */
    const char* fams[][2] = {
        {"Tan", "Sec"}, {"Cot", "Csc"}, {"Tanh", "Sech"}, {"Coth", "Csch"},
    };
    char buf[256];
    for (size_t f = 0; f < sizeof(fams) / sizeof(fams[0]); f++) {
        snprintf(buf, sizeof(buf),
            "Flatten[Table[Simplify[D[Integrate[%s[x]^m %s[x]^n, x, "
            "Method -> \"CRCTable\"], x] - %s[x]^m %s[x]^n], {m, 1, 4}, {n, 1, 4}]]",
            fams[f][0], fams[f][1], fams[f][0], fams[f][1]);
        assert_eval_eq(buf,
            "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}", 0);
    }
    /* The quotient input forms (which canonicalise into the above) and Sec Tan. */
    assert_eval_eq(
        "Simplify[D[Integrate[Cos[x]^2/Sin[x]^3, x, Method -> \"CRCTable\"], x] "
        "- Cos[x]^2/Sin[x]^3]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Sinh[x]^2/Cosh[x]^3, x, Method -> \"CRCTable\"], x] "
        "- Sinh[x]^2/Cosh[x]^3]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Sec[x] Tan[x], x, Method -> \"CRCTable\"], x] "
        "- Sec[x] Tan[x]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[Csc[2 x] Cot[2 x], x, Method -> \"CRCTable\"], x] "
        "- Csc[2 x] Cot[2 x]]", "0", 0);
}

static void test_crctable_polynomial_hyperbolic(void) {
    /* Polynomial x hyperbolic (added 2026-07-16): x^n Sinh/Cosh (generalised
     * from bare x to a x), x Sinh^m / x Cosh^m (Formula 381-386 analogues), and
     * x Csch^n / x Sech^n = x/Sinh^n, x/Cosh^n (Formula 410-413 analogues, keyed
     * on the canonical Csch/Sech heads).  All diff-back to the integrand.
     *
     * (Sinh[a x]/x^m -> SinhIntegral is intentionally NOT added: SinhIntegral /
     * CoshIntegral carry no derivative rule yet, so such results are neither
     * verifiable nor differentiable — a separate special-function task.) */
    char buf[256];
    for (int k = 1; k <= 4; k++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[x^%d Sinh[2 x], x, Method -> \"CRCTable\"], x] "
            "- x^%d Sinh[2 x]]", k, k);
        assert_eval_eq(buf, "0", 0);
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[x^%d Cosh[3 x], x, Method -> \"CRCTable\"], x] "
            "- x^%d Cosh[3 x]]", k, k);
        assert_eval_eq(buf, "0", 0);
    }
    const char* powcases[] = {
        "x Sinh[a x]^2", "x^2 Sinh[a x]^2", "x Sinh[a x]^3",
        "x Cosh[a x]^2", "x^2 Cosh[a x]^2", "x Cosh[a x]^3",
    };
    for (size_t i = 0; i < sizeof(powcases) / sizeof(powcases[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Simplify[D[Integrate[%s, x, Method -> \"CRCTable\"], x] - (%s)]",
            powcases[i], powcases[i]);
        assert_eval_eq(buf, "0", 0);
    }
    /* x/hyperbolic^n and x/trig^n (even n; keyed on canonical Csch/Sech/Csc/Sec). */
    const char* recips[] = {"Csch", "Sech", "Csc", "Sec"};
    for (size_t i = 0; i < sizeof(recips) / sizeof(recips[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Flatten[Table[Simplify[D[Integrate[x %s[x]^n, x, Method -> "
            "\"CRCTable\"], x] - x %s[x]^n], {n, 2, 6, 2}]]", recips[i], recips[i]);
        assert_eval_eq(buf, "{0, 0, 0}", 0);
    }
    /* The reciprocal input forms that canonicalise into the above. */
    assert_eval_eq(
        "Simplify[D[Integrate[x/Sinh[x]^2, x, Method -> \"CRCTable\"], x] "
        "- x/Sinh[x]^2]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[x/Cosh[2 x]^4, x, Method -> \"CRCTable\"], x] "
        "- x/Cosh[2 x]^4]", "0", 0);
}

void test_integrate_dispatch(void) {
    symtab_init();
    core_init();

    TEST(test_cascade_rational);
    TEST(test_cascade_risch);
    TEST(test_method_strict_rational);
    TEST(test_method_invalid);
    TEST(test_crctable_termination);
    TEST(test_crctable_reciprocal_sqrt_branch);
    TEST(test_crctable_secant_cosecant_powers);
    TEST(test_crctable_trig_power_family);
    TEST(test_crctable_hyperbolic_products);
    TEST(test_crctable_reciprocal_head_products);
    TEST(test_crctable_exp_hyperbolic);
    TEST(test_crctable_hyperbolic_argument_generality);
    TEST(test_crctable_mixed_tan_sec_families);
    TEST(test_crctable_polynomial_hyperbolic);
    TEST(test_arctanh_real_branch);
    TEST(test_root_kernel_keeps_radicals);

    printf("All Integrate dispatch tests passed!\n");
}

int main(void) {
    test_integrate_dispatch();
    return 0;
}

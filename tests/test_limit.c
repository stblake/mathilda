/*
 * test_limit.c -- Correctness tests for the native Limit builtin.
 *
 * Exercises the public surface of src/limit.c across all four return
 * categories called out in plans/limit_candidate_spec.md:
 *   - definite numeric / symbolic limits,
 *   - directional limits (including the Mathematica sign flip),
 *   - bounded-oscillation Interval returns,
 *   - unevaluated fall-through for unknown heads.
 *
 * Organised so that each new test can be added with a single assertion
 * line; the helpers do the plumbing.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Compare the InputForm (default) printed result of `input` to `expected`.
 * Mismatches print both strings for easy debugging. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Limit mismatch for %s:\n    expected: %s\n    got:      %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* Structural-equality check: accept any spelling that the printer might
 * emit by testing whether (result) - (expected) simplifies to 0.
 *
 * Mathematica's Infinity arithmetic refuses to fold Infinity - Infinity
 * (it is genuinely Indeterminate), so for results that are exactly an
 * infinity flavour we fall back to a literal printed-form comparison
 * rather than the subtractive equivalence check. */
static void check_equiv(const char* input, const char* expected) {
    /* First try direct printed-form equality -- this catches Infinity,
     * -Infinity, ComplexInfinity, Indeterminate, and any other case the
     * subtractive check cannot distinguish. */
    Expr* lhs = parse_expression(input);
    ASSERT_MSG(lhs != NULL, "Failed to parse: %s", input);
    Expr* lhs_v = evaluate(lhs);
    char* lhs_s = expr_to_string(lhs_v);
    Expr* rhs = parse_expression(expected);
    Expr* rhs_v = rhs ? evaluate(rhs) : NULL;
    char* rhs_s = rhs_v ? expr_to_string(rhs_v) : NULL;
    bool literal_match = rhs_s && strcmp(lhs_s, rhs_s) == 0;
    free(lhs_s);
    if (rhs_s) free(rhs_s);
    expr_free(lhs_v);
    expr_free(lhs);
    if (rhs_v) expr_free(rhs_v);
    if (rhs) expr_free(rhs);
    if (literal_match) return;

    char buf[2048];
    snprintf(buf, sizeof(buf), "Expand[Together[(%s) - (%s)]]", input, expected);
    Expr* e = parse_expression(buf);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", buf);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, "0") == 0,
               "Limit not equivalent: %s\n    expected ~ %s\n    diff:       %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* ----------------------------------------------------------------- */
/* Layer 1 -- fast paths                                              */
/* ----------------------------------------------------------------- */
static void test_fast_paths(void) {
    /* Free of the variable. */
    check("Limit[7, x -> 0]", "7");
    check("Limit[a, x -> 5]", "a");
    /* The variable itself. */
    check("Limit[x, x -> 5]", "5");
    /* Continuous polynomial / trig. */
    check("Limit[x^2, x -> 2]", "4");
    check("Limit[a x^2 + 4 x + 7, x -> 1]", "11 + a");
    check("Limit[Sin[Pi/2], x -> 0]", "1");
}

/* ----------------------------------------------------------------- */
/* Layer 2 -- series-based                                            */
/* ----------------------------------------------------------------- */
static void test_series(void) {
    check("Limit[Sin[x]/x, x -> 0]", "1");
    check("Limit[(Cos[x]-1)/(Exp[x^2]-1), x -> 0]", "-1/2");
    check("Limit[(Tan[x]-x)/x^3, x -> 0]", "1/3");
    check("Limit[Sin[2 x]/Sin[x], x -> Pi]", "-2");
    check("Limit[(x^a - 1)/a, a -> 0]", "Log[x]");
}

/* ----------------------------------------------------------------- */
/* Layer 3 -- rational functions                                      */
/* ----------------------------------------------------------------- */
static void test_rational(void) {
    check("Limit[x/(x + 1), x -> Infinity]", "1");
    check("Limit[3 x^2/(x^2 - 2), x -> -Infinity]", "3");
    check("Limit[x/(x - 2 a), x -> a]", "-1");
}

/* ----------------------------------------------------------------- */
/* Layer 5 -- log reduction                                           */
/* ----------------------------------------------------------------- */
static void test_log_reduction(void) {
    /* Classical 1^Infinity forms. */
    check("Limit[(1 + A x)^(1/x), x -> 0]", "E^A");
    check("Limit[(1 + a/x)^(b x), x -> Infinity]", "E^(a b)");
}

/* ----------------------------------------------------------------- */
/* Layer 4-ish -- RP-form families via Series                         */
/* ----------------------------------------------------------------- */
static void test_rp_forms(void) {
    check_equiv("Limit[Sqrt[x^2 + a x] - Sqrt[x^2 + b x], x -> Infinity]",
                "(a - b)/2");
}

/* ----------------------------------------------------------------- */
/* Directional limits + sign-convention tests                         */
/* ----------------------------------------------------------------- */
static void test_directions(void) {
    check("Limit[1/x, x -> 0, Direction -> \"FromAbove\"]", "Infinity");
    check("Limit[1/x, x -> 0, Direction -> \"FromBelow\"]", "-Infinity");
    /* Mathematica sign convention: -1 == FromAbove, +1 == FromBelow. */
    check("Limit[1/x, x -> 0, Direction -> -1]",           "Infinity");
    check("Limit[1/x, x -> 0, Direction -> 1]",            "-Infinity");
}

/* ----------------------------------------------------------------- */
/* Layer 6 -- bounded oscillation (Interval)                          */
/* ----------------------------------------------------------------- */
static void test_interval_returns(void) {
    /* Unbounded-oscillation limits return Indeterminate, matching
     * Mathematica. Previously we returned Interval[{-1, 1}] which is a
     * correct bound but not a limit value. */
    check("Limit[Sin[1/x], x -> 0]", "Indeterminate");
    check("Limit[Sin[x], x -> Infinity]", "Indeterminate");
}

/* ----------------------------------------------------------------- */
/* Iterated / multivariate forms                                      */
/* ----------------------------------------------------------------- */
static void test_iterated(void) {
    /* Iterated right-to-left: first y -> 0 gives x/x = 1, then x -> 0
     * of 1 is 1. */
    check("Limit[x/(x + y), {x -> 0, y -> 0}]", "1");
}

/* ----------------------------------------------------------------- */
/* Calling-form error handling                                        */
/* ----------------------------------------------------------------- */
static void test_unevaluated(void) {
    /* Non-rule second argument: Limit[f, a] should stay unevaluated. */
    check("Limit[x, x]", "Limit[x, x]");
    /* Unknown direction option: unevaluated. */
    check("Limit[x, x -> 0, Direction -> foo]",
          "Limit[x, x -> 0, Direction -> foo]");
}

/* ----------------------------------------------------------------- */
/* Reciprocal trig at singular points                                 */
/*                                                                    */
/* Every one of these collapses the 0 * ComplexInfinity trap unless   */
/* Csc/Sec/Cot are rewritten in terms of Sin/Cos before the layered   */
/* dispatcher runs. Covers the bulk of the new regression suite.     */
/* ----------------------------------------------------------------- */
static void test_reciprocal_trig(void) {
    check_equiv("Limit[x Csc[x], x -> 0]",              "1");
    check_equiv("Limit[x Csc[2 x], x -> 0]",            "1/2");
    check_equiv("Limit[Csc[7 x] Tan[2 x], x -> 0]",     "2/7");
    check_equiv("Limit[Csc[a x] Tan[b x], x -> 0]",     "b/a");
    check_equiv("Limit[Cot[4 x] Tan[3 x], x -> 0]",     "3/4");
    check_equiv("Limit[Cot[a x] Tan[b x], x -> 0]",     "b/a");
    check_equiv("Limit[x Cot[a x], x -> 0]",            "1/a");
    check_equiv("Limit[x Cot[x], x -> 0]",              "1");
    check_equiv("Limit[x^2 Csc[x], x -> 0]",            "0");
    check_equiv("Limit[(-1 + E^x) Cot[a x], x -> 0]",   "1/a");
    check_equiv("Limit[Csc[x] (x + Tan[x]), x -> 0]",   "2");
    check_equiv("Limit[(x - ArcSin[x]) Csc[x]^3, x -> 0]", "-1/6");
    check_equiv("Limit[Csc[Pi x] Log[x], x -> 1]",      "-1/Pi");
    check_equiv("Limit[Csch[a x] Sin[b x], x -> 0]",    "b/a");
    check_equiv("Limit[Csc[x]^2 Log[Cos[x]], x -> 0]",  "-1/2");
    check_equiv("Limit[Csc[2 x] Sin[3 x], x -> Pi]",    "-3/2");
    check_equiv("Limit[Sec[2 x] (1 - Tan[x]), x -> Pi/4]", "1");
    check_equiv("Limit[(-Pi/2 + x) Tan[3 x], x -> Pi/2]", "-1/3");
}

/* ----------------------------------------------------------------- */
/* Bounded envelope at +Infinity                                      */
/*                                                                    */
/* Families where |f| is dominated by something that decays to zero;  */
/* the squeeze layer short-circuits before Series trips on Sin[t^2]   */
/* (which has no Taylor expansion at infinity).                       */
/* ----------------------------------------------------------------- */
static void test_bounded_envelope(void) {
    check_equiv("Limit[Sin[t^2]/t^2, t -> Infinity]",         "0");
    check_equiv("Limit[(1 - Cos[x])/x, x -> Infinity]",       "0");
    check_equiv("Limit[(1 + Sin[x])/x, x -> Infinity]",       "0");
    check_equiv("Limit[(x Sin[x])/(5 + x^2), x -> Infinity]", "0");
    /* Interior bounded oscillation at 0 -- Sin[1/x] has no series at 0,
     * so the squeeze is the cleanest route.                              */
    check_equiv("Limit[x^2 Sin[1/x], x -> 0]",                "0");

    /* Bounded base raised to a divergent positive power. The base stays in
     * [-1/2, 1/2] while 1/x^2 -> +Infinity, so |f| <= (1/2)^(1/x^2) -> 0.  */
    check_equiv("Limit[(Sin[1/x]/2)^(1/x^2), x -> 0]",        "0");
    check_equiv("Limit[(Cos[1/x]/3)^(1/x^2), x -> 0]",        "0");
    /* At Infinity the exponent -> 0, so f -> 1 (envelope must not fire).  */
    check_equiv("Limit[(Sin[1/x]/2)^(1/x^2), x -> Infinity]", "1");
    /* Base magnitude bound >= 1: (2 Sin[1/x])^(1/x^2) is unbounded, so no
     * false 0 -- the shape is left unevaluated.                          */
    check("Limit[(2 Sin[1/x])^(1/x^2), x -> 0]",
          "Limit[(2 Sin[1/x])^(1/x^2), x -> 0]");
    /* Negative divergent exponent: the bound inequality flips and a
     * vanishing base blows up, so this too stays unevaluated.            */
    check("Limit[(Sin[1/x]/2)^(-1/x^2), x -> 0]",
          "Limit[(Sin[1/x]/2)^(-1/x^2), x -> 0]");

    /* Shrinking (x-dependent) magnitude bound: |x Sin[1/x]/2| <= Abs[x]/2,
     * whose limit is 0 < 1, so the squeeze fires even though the bound is
     * not a constant.                                                    */
    check_equiv("Limit[(x Sin[1/x]/2)^(1/x^2), x -> 0]",     "0");

    /* Base bounded below by a constant > 1 (band [3/2, 5/2], positive):
     * base^(1/x^2) -> +Infinity as the exponent -> +Infinity.           */
    check("Limit[(2 + Sin[1/x]/2)^(1/x^2), x -> 0]",         "Infinity");

    /* Bounded base with a *linearly divergent* exponent at Infinity. The
     * general log-reduction/Series path used to recurse without progress
     * (magnitude bound keeps the oscillator inside Abs[Log[..]]) and hang;
     * the squeeze now settles both quickly.                              */
    check_equiv("Limit[(Cos[x]/2)^x, x -> Infinity]",        "0");
    check_equiv("Limit[(Sin[x]/2)^x, x -> Infinity]",        "0");

    /* Genuine 1^Infinity with an unbounded-oscillation log (x Sin[x] has no
     * limit): must terminate and stay unevaluated, not hang.            */
    check("Limit[(1 + Sin[x]/x)^(x^2), x -> Infinity]",
          "Limit[(1 + Sin[x]/x)^x^2, x -> Infinity]");

    /* Convergent 1^Infinity through a bounded head (Cos[1/x] -> 1): once
     * the squeeze stops the hang, Series recovers the true value E^(-1/2). */
    check_equiv("Limit[Cos[1/x]^(x^2), x -> Infinity]",      "1/Sqrt[E]");
}

/* ----------------------------------------------------------------- */
/* ArcTan / ArcCot with divergent inner argument                     */
/* ----------------------------------------------------------------- */
static void test_arctan_infinity(void) {
    check_equiv("Limit[ArcTan[x^2 - x^4], x -> Infinity]", "-Pi/2");
}

/* ----------------------------------------------------------------- */
/* Asymptotic values of the builtin math functions at +/-Infinity.   */
/*                                                                    */
/* layer_compose_at_infinity applies f to the (infinite) inner limit  */
/* and lets the builtin fold it. Every function that self-evaluates   */
/* at Infinity is therefore covered by a single mechanism. Oscillatory */
/* heads (Sin, Cos) stay Indeterminate.                               */
/* ----------------------------------------------------------------- */
static void test_asymptotic_at_infinity(void) {
    /* Error functions -- the motivating gap. */
    check("Limit[Erf[x], x -> Infinity]",   "1");
    check("Limit[Erf[x], x -> -Infinity]",  "-1");
    check("Limit[Erfc[x], x -> Infinity]",  "0");
    check("Limit[Erfc[x], x -> -Infinity]", "2");
    check("Limit[Erfi[x], x -> Infinity]",  "Infinity");
    check("Limit[Erfi[x], x -> -Infinity]", "-Infinity");
    /* Composed inner argument: inner limit is still +/-Infinity. */
    check("Limit[Erf[3 x + 1], x -> Infinity]",  "1");
    check("Limit[Erf[1 - x^3], x -> Infinity]",  "-1");

    /* Exponential / logarithm. */
    check("Limit[Exp[x], x -> Infinity]",   "Infinity");
    check("Limit[Exp[x], x -> -Infinity]",  "0");
    check("Limit[Log[x], x -> Infinity]",   "Infinity");

    /* Hyperbolic + inverse hyperbolic. */
    check("Limit[Sinh[x], x -> Infinity]",  "Infinity");
    check("Limit[Cosh[x], x -> -Infinity]", "Infinity");
    check("Limit[Tanh[x], x -> Infinity]",  "1");
    check("Limit[Tanh[x], x -> -Infinity]", "-1");
    check("Limit[Coth[x], x -> Infinity]",  "1");
    check("Limit[Sech[x], x -> Infinity]",  "0");
    check("Limit[Csch[x], x -> Infinity]",  "0");
    check("Limit[ArcSinh[x], x -> Infinity]", "Infinity");
    check("Limit[ArcCoth[x], x -> Infinity]", "0");

    /* Inverse trig (values now self-evaluate in the builtins). */
    check_equiv("Limit[ArcTan[x], x -> Infinity]",  "Pi/2");
    check_equiv("Limit[ArcTan[x], x -> -Infinity]", "-Pi/2");
    check("Limit[ArcCot[x], x -> Infinity]",  "0");
    check("Limit[ArcCot[x], x -> -Infinity]", "0");   /* odd: 0 (not Pi) */
    check_equiv("Limit[ArcSec[x], x -> Infinity]",  "Pi/2");
    check("Limit[ArcCsc[x], x -> Infinity]",  "0");

    /* Special functions. */
    check("Limit[Gamma[x], x -> Infinity]",     "Infinity");
    check("Limit[Zeta[x], x -> Infinity]",      "1");
    check("Limit[ProductLog[x], x -> Infinity]", "Infinity");

    /* Oscillatory heads have no limit at Infinity. */
    check("Limit[Sin[x], x -> Infinity]", "Indeterminate");
    check("Limit[Cos[x], x -> Infinity]", "Indeterminate");
}

/* ----------------------------------------------------------------- */
/* Log-reduction of Power-in-x-exponent indeterminate forms          */
/*                                                                    */
/* Covers both the single-Power case (handled by the old code) and   */
/* the product-of-Powers generalization.                              */
/* ----------------------------------------------------------------- */
static void test_exp_indeterminate(void) {
    /* Classical 1^Infinity / 0^0 shapes that continuous substitution
     * cannot touch. */
    check_equiv("Limit[(-1 + x)^(-1 + x), x -> 1]",                 "1");
    check_equiv("Limit[Cos[x]^Cot[x], x -> 0]",                     "1");
    check_equiv("Limit[Sec[x]^Cot[x], x -> 0]",                     "1");
    check_equiv("Limit[(1 - 1/x^2)^Cot[1/x], x -> Infinity]",       "1");
    check_equiv("Limit[(1 - 1/x)^(x^2), x -> Infinity]",            "0");
    check_equiv("Limit[(1 + Sin[a x])^Cot[b x], x -> 0]",           "E^(a/b)");
    check_equiv("Limit[x^Tan[(Pi x)/2], x -> 1]",                   "E^(-2/Pi)");
    check_equiv("Limit[((-3 + 2*x)/(5 + 2*x))^(1 + 2*x), x -> Infinity]",
                "1/E^8");
    /* Polynomial-dominates-exponential decay (product-of-Powers log-
     * reduction discovers this via Log[E^(-x) x^20] = -x + 20 Log[x] -> -Inf). */
    check_equiv("Limit[E^(-x) x^20, x -> Infinity]", "0");
}

/* ----------------------------------------------------------------- */
/* Continuous-substitution residuals                                  */
/*                                                                    */
/* Power[0, positive] appears naturally when substituting directly:   */
/* Sqrt[x-1]/x at x=1 evaluates to Sqrt[0]/1 which Mathilda's Power    */
/* evaluator leaves un-folded. The limit module has a small folder   */
/* so the answer printed is the plain 0.                              */
/* ----------------------------------------------------------------- */
static void test_power_zero_folding(void) {
    check_equiv("Limit[Sqrt[x - 1]/x, x -> 1]", "0");
}

/* ----------------------------------------------------------------- */
/* Miscellaneous exact cases                                          */
/* ----------------------------------------------------------------- */
static void test_misc(void) {
    check_equiv("Limit[x (-Sqrt[x] + Sqrt[2 + x]), x -> Infinity]", "Infinity");
    check_equiv("Limit[E^(9/x)*(-2 + x) - x, x -> Infinity]", "7");
    check_equiv("Limit[(Tan[x] - x)/x^3, x -> 0]", "1/3");
}

/* Regression set added 2026-04-20 (batch 2).
 * Includes the Factor[(1+t)^(-1/2)] segfault case and eight challenging
 * limits reported from the REPL, covering Csc * Log, (1+sin)^cot forms,
 * ArcTan at infinity, Csc^3 singularity, Log merge at infinity, and
 * Sin/Cos removable singularities. */
static void test_regression_batch2(void) {
    check_equiv("Limit[1/(t*Sqrt[t + 1]) - 1/t, t -> 0]", "-1/2");
    check_equiv("Limit[Csc[Pi*x]*Log[x], x -> 1]", "-1/Pi");
    check_equiv("Limit[(1 + Sin[a x])^Cot[b x], x -> 0]", "E^(a/b)");
    check_equiv("Limit[ArcTan[x^2 - x^4], x -> Infinity]", "-Pi/2");
    check_equiv("Limit[(x - ArcSin[x]) Csc[x]^3, x -> 0]", "-1/6");
    check_equiv("Limit[Cos[x]^Cot[x], x -> 0]", "1");
    check_equiv("Limit[-x + Log[2 + E^x], x -> Infinity]", "0");
    check_equiv("Limit[(Sin[x] - Cos[x])/Cos[2x], x -> Pi/4]", "-1/Sqrt[2]");
    check_equiv("Limit[(Sqrt[x^2 + 4] - Sqrt[8])/(x^2 - 4), x -> 2]", "Sqrt[2]/8");
    /* Radical-ratio post-processing: Sqrt[6]/Sqrt[2] -> Sqrt[3], and the
     * 0/0 Sqrt-denominator path must not emit a spurious Power::infy. */
    check("Limit[Sqrt[x^2 - 9]/Sqrt[2 x - 6], x -> 3]", "Sqrt[3]");
}

/* Regression set added 2026-04-21.
 * Covers a batch of failures reported from the REPL:
 *   - Display normalisation (-1 Infinity -> -Infinity)
 *   - Unknown / discontinuous head refusal (stay unevaluated)
 *   - Even-order pole returning +Infinity rather than ComplexInfinity
 *   - Oscillatory limit returning Indeterminate (not Interval[{-1,1}])
 *   - Direction -> I parsed as complex (stays unevaluated, doesn't error)
 *   - Leaked limit variable in directional-infinity coefficient
 */
static void test_regression_batch3(void) {
    /* -Infinity display: formerly printed as "-1 Infinity". */
    check("Limit[1/x, x -> 0, Direction -> \"FromBelow\"]", "-Infinity");
    check("Limit[-1/(x-2)^2, x -> 2]", "-Infinity");

    /* Even-order pole: two-sided limit is +Infinity, not ComplexInfinity. */
    check("Limit[1/(x-2)^2, x -> 2]", "Infinity");
    check("Limit[1/x^2, x -> 0]", "Infinity");
    /* Odd-order pole keeps the old ComplexInfinity behaviour. */
    check("Limit[1/(x-2)^3, x -> 2]", "ComplexInfinity");

    /* Unknown head stays unevaluated instead of silently substituting. */
    check("Limit[f[x], x -> a]", "Limit[f[x], x -> a]");

    /* Discontinuous heads (Floor, Ceiling, FractionalPart, ...) are treated
     * as opaque: Limit refuses to substitute rather than pick one side. */
    check("Limit[Ceiling[5 - x^2], x -> 1]",
          "Limit[Ceiling[5 - x^2], x -> 1]");
    check("Limit[FractionalPart[x^2] Sin[x], x -> 2, Direction -> \"FromBelow\"]",
          "Limit[FractionalPart[x^2] Sin[x], x -> 2, Direction -> \"FromBelow\"]");

    /* Unbounded oscillation -> Indeterminate, not a bounded interval. */
    check("Limit[Sin[1/x], x -> 0]", "Indeterminate");
    check("Limit[Sin[x], x -> Infinity]", "Indeterminate");
    check("Limit[Cos[x], x -> Infinity]", "Indeterminate");

    /* Sin[x]/x at infinity: squeeze to 0. */
    check("Limit[Sin[x]/x, x -> Infinity]", "0");

    /* Series left the limit variable inside a DirectedInfinity coefficient
     * (Log[x]/x at 0). Bail to unevaluated instead. */
    check("Limit[Log[x]/x, x -> 0]", "Limit[Log[x]/x, x -> 0]");

    /* Direction -> I now threads and flips the imaginary part (covered
     * by test_wp9_branch_cuts below). The previous batch-3 assertion
     * pinned the unevaluated fall-through; WP-9 replaces it. */
}

/* ----------------------------------------------------------------- */
/* WP-8 -- Finite-point direct-substitution fast path                 */
/*                                                                    */
/* Numeric limit points where the expression is analytic must resolve */
/* by plain substitution, even when the expression is large enough    */
/* that the Series or L'Hospital layers would spin on it. The two     */
/* tests below each time-bound their own execution indirectly: if     */
/* substitution didn't fire we'd block the test runner indefinitely.  */
/* ----------------------------------------------------------------- */
static void test_wp8_numeric_point_bypass(void) {
    /* Previously hung: Log[1-(Log[Exp[z]/z-1]+Log[z])/z]/z at z=100.
     * Now returns the analytic value via direct substitution. We assert
     * structural equivalence rather than a single canonical string, so
     * printer-formatting changes don't break the test. */
    check_equiv(
        "Limit[Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z, z -> 100]",
        "Log[1 + (-Log[100] - Log[-1 + E^100/100])/100]/100");

    /* Trivial analytic-at-numeric-point cases that exercise the new
     * fast path; they used to resolve through the Together path too,
     * but we pin them so a future refactor can't silently drop the
     * short-circuit. */
    check("Limit[Exp[z] + z^3, z -> 2]", "8 + E^2");
    check("Limit[Log[z] + 1/z, z -> 1]", "1");

    /* Sanity: the fast path must NOT short-cut when the denominator
     * actually vanishes at the point. Sin[x]/x at 0 must still give 1
     * through the Series / L'Hospital layers, not 0 via substitution. */
    check("Limit[Sin[x]/x, x -> 0]", "1");
    check("Limit[(1 - Cos[x])/x^2, x -> 0]", "1/2");
}

/* ----------------------------------------------------------------- */
/* WP-2 -- b^x kernel: series inversion with symbolic-Log coefs       */
/*                                                                    */
/* Previously hung inside so_inv because each iteration multiplied    */
/* growing polynomials in Log[a], Log[b], ... and Mathilda's evaluator */
/* did not collapse them to a canonical form. The fix caps the so_inv */
/* iteration count based on total leaf count of input coefficients:   */
/* large symbolic coefficients use a smaller N, so we still get the   */
/* leading Laurent term (enough for Limit) in bounded time.           */
/* ----------------------------------------------------------------- */
static void test_wp2_bx_kernel(void) {
    /* Used to hang; now resolves via Series + so_inv cap. Matches
     * Mathematica's raw output (no auto-Log[25] -> 2 Log[5] fold). */
    check_equiv("Limit[(25^x - 5^x)/(4^x - 2^x), x -> 0]",
                "(Log[25] - Log[5])/(Log[4] - Log[2])");

    /* Classic derivative-definition limit. Unchanged by WP-2 but a
     * useful control that the cap didn't regress the common case. */
    check_equiv("Limit[(a^x - 1)/x, x -> 0]", "Log[a]");
    check_equiv("Limit[(2^x - 1)/x, x -> 0]", "Log[2]");

    /* Two distinct bases in both numerator and denominator was the
     * original hang. */
    check_equiv("Limit[(3^x - 2^x)/(5^x - 2^x), x -> 0]",
                "(Log[3] - Log[2])/(Log[5] - Log[2])");

    /* Pin Series on the problematic shapes so a future tighter cap
     * doesn't silently drop the leading term. */
    check_equiv("Normal[Series[1/(2^x - 3^x), {x, 0, 0}]]",
                "1/((Log[2] - Log[3]) x) - (Log[2]^2 - Log[3]^2)"
                "/(2 (Log[2] - Log[3])^2)");
}

/* ----------------------------------------------------------------- */
/* WP-4 -- Pole-with-sign-disagreement classifier                     */
/*                                                                    */
/* Mathematica's Direction -> Reals specifically asks for the real-   */
/* line answer at a pole. If the two one-sided limits disagree in     */
/* sign (odd-order real pole, Sqrt-type branch singularity, Tan at    */
/* Pi/2) the real-line limit does not exist -> Indeterminate.         */
/*                                                                    */
/* Direction -> Complexes asks for the radial / all-direction answer; */
/* any isolated pole is ComplexInfinity regardless of parity.         */
/*                                                                    */
/* The implicit two-sided default (no Direction option) keeps the old */
/* ComplexInfinity fall-back -- it is the older Mathematica convention */
/* for unqualified two-sided limits on rational functions and our     */
/* internal tests pin it.                                             */
/* ----------------------------------------------------------------- */
static void test_wp4_pole_sign_disagreement(void) {
    /* Tan[x] at Pi/2: two-sided disagreement (+inf vs -inf). */
    check("Limit[Tan[x], x -> Pi/2]",                       "ComplexInfinity");
    check("Limit[Tan[x], x -> Pi/2, Direction -> Reals]",   "Indeterminate");
    check("Limit[Tan[x], x -> Pi/2, Direction -> Complexes]", "ComplexInfinity");
    /* One-sided still goes to the signed-infinity answer. */
    check("Limit[Tan[x], x -> Pi/2, Direction -> \"FromBelow\"]", "Infinity");
    check("Limit[Tan[x], x -> Pi/2, Direction -> \"FromAbove\"]", "-Infinity");

    /* Sqrt branch singularity: Sqrt[x] - 3 at x = 9 hits 0 but only
     * approaches from above on the principal branch. */
    check("Limit[(x - 3)/(Sqrt[x] - 3), x -> 9, Direction -> Reals]",
          "Indeterminate");

    /* Sin[Tan[x]] at Pi/2 Reals: Tan diverges both ways with sign
     * disagreement, so the composed limit is Indeterminate too.
     * Previously returned Interval[{-1, 1}] then Indeterminate via WP-3
     * batch; this test pins the Direction -> Reals path specifically. */
    check("Limit[Sin[Tan[x]], x -> Pi/2, Direction -> Reals]",
          "Indeterminate");

    /* Odd-order rational-function pole: default direction keeps the
     * ComplexInfinity convention, explicit Reals asks for Indeterminate,
     * Complexes says ComplexInfinity. */
    check("Limit[1/(x-2)^3, x -> 2]",                       "ComplexInfinity");
    check("Limit[1/(x-2)^3, x -> 2, Direction -> Reals]",   "Indeterminate");
    check("Limit[1/(x-2)^3, x -> 2, Direction -> Complexes]", "ComplexInfinity");

    /* Even-order pole: both sides agree on +Infinity regardless of the
     * direction-mode interpretation. */
    check("Limit[1/(x-2)^2, x -> 2]",                        "Infinity");
    check("Limit[1/(x-2)^2, x -> 2, Direction -> Reals]",    "Infinity");
    check("Limit[1/(x-2)^2, x -> 2, Direction -> Complexes]", "ComplexInfinity");
}

/* ----------------------------------------------------------------- */
/* WP-1 -- Multivariate path-dependence analyzer                      */
/*                                                                    */
/* Joint limits {x1,...,xn} -> {a1,...,an} run a polar / spherical    */
/* substitution plus a direction-sampling cross-check. The result:    */
/*   - polar r-limit free of angles -> that value;                    */
/*   - polar r-limit depends on angles -> sample axis and diagonal    */
/*     paths; if they all agree return the common value, else         */
/*     Indeterminate.                                                 */
/*                                                                    */
/* The origin-case fast path also scans for any reciprocal subterm    */
/* whose base vanishes at the joint point, which catches 0/0 inside   */
/* ArcTan/Sin/etc. that Mathilda would otherwise fold to 0.            */
/* ----------------------------------------------------------------- */
static void test_wp1_multivariate(void) {
    /* Single-direction collapse: u = x y, so Tan[u]/u -> 1 as u -> 0. */
    check("Limit[Tan[x y]/(x y), {x, y} -> {0, 0}]", "1");

    /* Classic path-dependent origin limit. */
    check("Limit[y/(x + y), {x, y} -> {0, 0}]", "Indeterminate");

    /* Polar r-limit independent of angle. */
    check("Limit[(x^3 + y^3)/(x^2 + y^2), {x, y} -> {0, 0}]", "0");

    /* 0/0 inside ArcTan: without the inner-divide-by-zero scan,
     * Mathilda's ArcTan[0/0] -> ArcTan[0] = 0 fold would produce a
     * spurious finite answer. Sampling gives disagreement. */
    check("Limit[ArcTan[y^2/(x^2 + x^3)], {x, y} -> {0, 0}]",
          "Indeterminate");

    /* 3D spherical. Depends on both theta and phi, diverges between
     * path samples. */
    check("Limit[(x z)/(x^2 + y^2 + z^2), {x, y, z} -> {0, 0, 0}]",
          "Indeterminate");

    /* 2D reduction of the 3D case with y pinned at 0. */
    check("Limit[(x z)/(x^2 + y^2 + z^2) /. y -> 0, {x, z} -> {0, 0}]",
          "Indeterminate");

    /* Joint at infinity: ArcTan[y/x] is path-dependent along the
     * positive orthant (different y/x slopes give different angles). */
    check("Limit[ArcTan[y/x], {x, y} -> {Infinity, Infinity}]",
          "Indeterminate");
    check("Limit[x/y, {x, y} -> {Infinity, Infinity}]",
          "Indeterminate");

    /* Polar limit that IS continuous at origin; should pass through. */
    check("Limit[x y/(x^2 + y^2 + 1), {x, y} -> {0, 0}]", "0");
}

/* ----------------------------------------------------------------- */
/* WP-3 -- Series of x^a at a nonzero finite point                    */
/*                                                                    */
/* Previously `try_factor_power_prefactor` pulled `x^a` out of any    */
/* top-level Times regardless of x0; at x0 = 1 that left the other    */
/* factors to be expanded in isolation, producing 0. Two fixes combine */
/* here:                                                              */
/*   - The prefactor factoring is now gated to x0 = 0 or Infinity.    */
/*   - Apart preprocessing is skipped when the expression has a       */
/*     non-rational power (e.g. x^a), because Mathilda's Apart         */
/*     collapses such inputs to 0.                                    */
/* ----------------------------------------------------------------- */
static void test_wp3_xpow_a_at_nonzero(void) {
    /* The motivating regression. 3 (x^a - a x + a - 1) has a second-order
     * zero at x = 1 with leading coefficient 3 * a(a-1)/2; dividing by
     * (x-1)^2 gives the stated limit. */
    check_equiv("Limit[3 (x^a - a x + a - 1)/(x-1)^2, x -> 1]",
                "3/2 a (a - 1)");

    /* Pin Series at the nonzero point so future regressions don't
     * silently collapse to 0 again. */
    check_equiv("Normal[Series[x^a, {x, 1, 2}]]",
                "1 + a (x - 1) + (1/2) a (a - 1) (x - 1)^2");
    check_equiv("Normal[Series[x^a/(x-1)^2, {x, 1, 0}]]",
                "1/(x - 1)^2 + a/(x - 1) + (1/2) a (a - 1)");
    check_equiv("Normal[Series[x^a (x-1), {x, 1, 2}]]",
                "(x - 1) + a (x - 1)^2");
}

/* ----------------------------------------------------------------- */
/* WP-6 -- Sinh/Cosh/Tanh exponentialisation at infinity              */
/*                                                                    */
/* Rewrites Sinh[z] -> (E^z - E^-z)/2, Cosh[z] -> (E^z + E^-z)/2,     */
/* Tanh[z] -> ratio, then Expand[]s the whole expression. After this  */
/* the term-wise Plus layer catches cancelling Exp[kx] summands that  */
/* the series kernels on Sinh/Cosh themselves couldn't see.           */
/* ----------------------------------------------------------------- */
static void test_wp6_hyperbolic_at_infinity(void) {
    /* The motivating failure: (1 + Sinh[x])/Exp[x] -> 1/2. After the
     * rewrite this expands to 1/2 + E^-x - E^-2x/2, which is a
     * term-wise sum with a finite limit. */
    check("Limit[(1 + Sinh[x])/Exp[x], x -> Infinity]", "1/2");
    check("Limit[Cosh[x]/Exp[x], x -> Infinity]", "1/2");

    /* -Infinity case: Sinh[-x] rewrites to (E^-x - E^x)/2; Expand
     * still sums term-wise to the right value. */
    check("Limit[Cosh[x]/Exp[-x], x -> -Infinity]", "1/2");

    /* Limits at finite points must keep their old behaviour -- we
     * only rewrite at Infinity. Sinh[x] at 0 is 0, Cosh[x] at 0 is 1
     * via Taylor expansion, not the Exp rewrite. */
    check("Limit[Sinh[x], x -> 0]", "0");
    check("Limit[Cosh[x], x -> 0]", "1");
    check("Limit[Sinh[x]/x, x -> 0]", "1");
}

/* ----------------------------------------------------------------- */
/* WP-5 -- Sign-aware envelope / dominant-term at infinity            */
/*                                                                    */
/* A Plus at +/-Infinity where one term strictly dominates (higher    */
/* polynomial growth than everything else) evaluates to that term's   */
/* limit. Bounded oscillators (Sin/Cos/Tanh/ArcTan) count as growth   */
/* exponent 0, so they are absorbed by any positive-degree dominant   */
/* polynomial term.                                                   */
/* ----------------------------------------------------------------- */
static void test_wp5_dominant_term_at_infinity(void) {
    /* Dominant x^2 absorbs x Sin[x^2] (growth 1 * 0 = 1). */
    check("Limit[x^2 + x Sin[x^2], x -> Infinity]", "Infinity");

    /* Dominant x absorbs Sin[x] (bounded). */
    check("Limit[x + Sin[x], x -> Infinity]", "Infinity");
    check("Limit[x - Cos[x], x -> Infinity]", "Infinity");

    /* At -Infinity: dominant x -> -Infinity, Sin[x] is still bounded. */
    check("Limit[x + Sin[x], x -> -Infinity]", "-Infinity");

    /* Two bounded terms with no dominant -> unresolved (safer than a
     * wrong finite answer). */
    check("Limit[Sin[x^2] + Cos[x], x -> Infinity]",
          "Limit[Sin[x^2] + Cos[x], x -> Infinity]");
}

/* ----------------------------------------------------------------- */
/* WP-9 -- Complex-direction branch-cut handling                      */
/*                                                                    */
/* At a branch point on the negative real axis (Sqrt, Log, etc.):     */
/*   Direction -> Reals   keeps the principal-branch value            */
/*   Direction -> I       flips the sign of the imaginary part        */
/*   Direction -> Complexes returns Indeterminate for a branch point  */
/*     whose radial approach gives different values on either side    */
/*     of the cut.                                                    */
/* Limit also threads over a top-level List in its first argument.    */
/* ----------------------------------------------------------------- */
static void test_wp9_branch_cuts(void) {
    /* Sqrt at z = -1. */
    check("Limit[Sqrt[z], z -> -1, Direction -> Reals]",      "I");
    check("Limit[Sqrt[z], z -> -1, Direction -> I]",          "-I");
    check("Limit[Sqrt[z], z -> -1, Direction -> Complexes]",  "Indeterminate");

    /* Log at z = -1. */
    check("Limit[Log[z], z -> -1, Direction -> Reals]",      "I Pi");
    check("Limit[Log[z], z -> -1, Direction -> I]",          "-I Pi");
    check("Limit[Log[z], z -> -1, Direction -> Complexes]",  "Indeterminate");

    /* List threading in the first argument. */
    check("Limit[{Sqrt[z], Log[z]}, z -> -1, Direction -> I]",
          "{-I, -I Pi}");
    check("Limit[{Sqrt[z], Log[z]}, z -> -1, Direction -> Complexes]",
          "{Indeterminate, Indeterminate}");

    /* Real-valued poles ignore the I/Complexes post-pass -- still a
     * signed infinity or ComplexInfinity. */
    check("Limit[1/(x-2), x -> 2, Direction -> I]", "ComplexInfinity");
}

/* ----------------------------------------------------------------- */
/* WP-7 -- Gruntz-lite: Log[sum] with dominant summand                */
/*                                                                    */
/* Partial Gruntz implementation: rewrites Log[dom + rest...] at      */
/* +Infinity into Log[dom] + Log[1 + rest/dom] when one summand       */
/* strictly dominates (higher growth_exponent). This unblocks a       */
/* subset of iterated-log limits. Full Hardy-field comparative        */
/* asymptotics (stacked Exp, nested logs with multiple dominance      */
/* levels, Sin[x] + Log[x-a]/Log[E^x-E^a], etc.) remain future work.  */
/* ----------------------------------------------------------------- */
static void test_wp7_gruntz_log_sum(void) {
    /* A direct test of the rewrite: Log[x + Log[x]] at infinity.
     * dom = x, rest = Log[x]. The ratio Log[x]/x -> 0 so the limit of
     * Log[x + Log[x]] - Log[x] = Log[1 + Log[x]/x] -> Log[1] = 0. */
    check_equiv("Limit[Log[x + Log[x]] - Log[x], x -> Infinity]", "0");

    /* A finite non-zero case: Log[x^2 + x] - 2 Log[x]. Rewrite: Log[x]
     * doesn't dominate x^2, so the Log[sum] rewrites with dom = x^2,
     * rest = x, ratio x/x^2 = 1/x -> 0. Log[x^2 + x] -> 2 Log[x], so
     * the difference -> 0. */
    check_equiv("Limit[Log[x^2 + x] - 2 Log[x], x -> Infinity]", "0");
}

/* ----------------------------------------------------------------- */
/* Abs[g(x)] in directional limits                                    */
/*                                                                    */
/* Abs is continuous everywhere but non-differentiable at the kink.   */
/* L'Hospital and Series both fail on Abs[g(x)] near a zero of g and  */
/* would otherwise emit Derivative[1][Abs][a] -- a non-numeric        */
/* symbolic value -- as a "clean" answer. The Abs-rewrite layer       */
/* resolves these by direction:                                       */
/*                                                                    */
/*   FromAbove (x -> a^+):  Abs[g(x)] -> sign(g near a^+) g(x)        */
/*   FromBelow (x -> a^-):  Abs[g(x)] -> sign(g near a^-) g(x)        */
/*                                                                    */
/* For two-sided limits where the one-sided values disagree, the      */
/* genuine answer is Indeterminate. The continuity identity           */
/* Abs[Limit[g]] = Limit[Abs[g]] handles divergent inner limits        */
/* (Abs[Tan[x]] near Pi/2, Abs[1/x] near 0, Abs[Log[x]] near 0+, ...).*/
/* ----------------------------------------------------------------- */
static void test_abs_directional(void) {
    /* The canonical kink: Abs[x]/x at x = 0. Direction -> -1 and +1 use
     * Mathematica's sign convention (-1 = approach from above). */
    check_equiv("Limit[Abs[x]/x, x -> 0, Direction -> -1]", "1");
    check_equiv("Limit[Abs[x]/x, x -> 0, Direction -> 1]", "-1");
    check_equiv("Limit[Abs[x]/x, x -> 0, Direction -> \"FromAbove\"]", "1");
    check_equiv("Limit[Abs[x]/x, x -> 0, Direction -> \"FromBelow\"]", "-1");
    check_equiv("Limit[Abs[x]/x, x -> 0]", "Indeterminate");

    /* Polynomial argument with a zero at the limit point: factor
     * (x-1)(x+1) so the sign of x^2 - 1 changes through x = 1. */
    check_equiv("Limit[Abs[x^2 - 1]/(x - 1), x -> 1, Direction -> -1]", "2");
    check_equiv("Limit[Abs[x^2 - 1]/(x - 1), x -> 1, Direction -> 1]", "-2");
    check_equiv("Limit[Abs[x^2 - 1]/(x - 1), x -> 1]", "Indeterminate");

    /* Linear argument vanishing at a non-zero point. */
    check_equiv("Limit[Abs[x + 1]/(x + 1), x -> -1, Direction -> -1]", "1");
    check_equiv("Limit[Abs[x + 1]/(x + 1), x -> -1, Direction -> 1]", "-1");
    check_equiv("Limit[Abs[x + 1]/(x + 1), x -> -1]", "Indeterminate");

    /* Even-order vanishing: x^2 >= 0 globally so both sides agree. */
    check_equiv("Limit[Abs[x^2]/x^2, x -> 0]", "1");
    check_equiv("Limit[Abs[x^2]/x, x -> 0]", "0");

    /* Cubic kink (odd order vanishing): sides disagree. */
    check_equiv("Limit[Abs[x^3]/x^3, x -> 0, Direction -> -1]", "1");
    check_equiv("Limit[Abs[x^3]/x^3, x -> 0, Direction -> 1]", "-1");
    check_equiv("Limit[Abs[x^3]/x^3, x -> 0]", "Indeterminate");

    /* Trig argument vanishing at the limit point: sign of Sin[x] near 0
     * matches the side. */
    check_equiv("Limit[Abs[Sin[x]]/x, x -> 0, Direction -> -1]", "1");
    check_equiv("Limit[Abs[Sin[x]]/x, x -> 0, Direction -> 1]", "-1");
    check_equiv("Limit[Abs[Sin[x]]/x, x -> 0]", "Indeterminate");

    /* Non-trivial Taylor leading term: Sin[x] - x = -x^3/6 + ... */
    check_equiv("Limit[Abs[Sin[x] - x]/x^3, x -> 0, Direction -> -1]", "1/6");
    check_equiv("Limit[Abs[Sin[x] - x]/x^3, x -> 0, Direction -> 1]", "-1/6");
    check_equiv("Limit[Abs[Sin[x] - x]/x^3, x -> 0]", "Indeterminate");

    /* Even-leading-order argument: Cos[x] - 1 = -x^2/2 + ... is negative
     * on both sides of 0, so Abs[Cos[x] - 1] = 1 - Cos[x] and the limit
     * is unambiguous (no direction split needed). */
    check_equiv("Limit[Abs[Cos[x] - 1]/x^2, x -> 0]", "1/2");
    check_equiv("Limit[Abs[1 - Cos[x]]/x^2, x -> 0]", "1/2");

    /* Exponential-minus-one: leading term is x. */
    check_equiv("Limit[Abs[E^x - 1]/x, x -> 0, Direction -> -1]", "1");
    check_equiv("Limit[Abs[E^x - 1]/x, x -> 0, Direction -> 1]", "-1");
    check_equiv("Limit[Abs[E^(2 x) - 1]/x, x -> 0, Direction -> -1]", "2");
    check_equiv("Limit[Abs[E^(2 x) - 1]/x, x -> 0, Direction -> 1]", "-2");
    check_equiv("Limit[Abs[E^(2 x) - 1]/x, x -> 0]", "Indeterminate");

    /* Sign-cancellation in a sum: Abs[x] + x equals 2x for x>0, 0 for x<0. */
    check_equiv("Limit[(Abs[x] + x)/x, x -> 0, Direction -> -1]", "2");
    check_equiv("Limit[(Abs[x] + x)/x, x -> 0, Direction -> 1]", "0");

    /* x/Abs[x] is the sign function -- one-sided gives ±1, two-sided is
     * the canonical disagreement. */
    check_equiv("Limit[x/Abs[x], x -> 0, Direction -> -1]", "1");
    check_equiv("Limit[x/Abs[x], x -> 0, Direction -> 1]", "-1");
    check_equiv("Limit[x/Abs[x], x -> 0]", "Indeterminate");

    /* Continuity identity Abs[Limit[g]] = Limit[Abs[g]] for divergent
     * inner limits. */
    check_equiv("Limit[Abs[1/x], x -> 0, Direction -> -1]", "Infinity");
    check_equiv("Limit[Abs[1/x], x -> 0, Direction -> 1]", "Infinity");
    check_equiv("Limit[Abs[1/x], x -> 0]", "Infinity");
    check_equiv("Limit[Abs[Tan[x]], x -> Pi/2, Direction -> -1]", "Infinity");
    check_equiv("Limit[Abs[Tan[x]], x -> Pi/2, Direction -> 1]", "Infinity");
    check_equiv("Limit[Abs[Tan[x]], x -> Pi/2]", "Infinity");
    check_equiv("Limit[Abs[Log[x]], x -> 0, Direction -> -1]", "Infinity");

    /* Abs of the limit variable approaching a finite point or infinity
     * stays continuous. */
    check_equiv("Limit[Abs[x], x -> 0]", "0");
    check_equiv("Limit[x Abs[x], x -> 0]", "0");
    check_equiv("Limit[Abs[x]/x, x -> Infinity]", "1");
    check_equiv("Limit[Abs[x]/x, x -> -Infinity]", "-1");

    /* Nested / product Abs: each subterm resolves independently. */
    check_equiv("Limit[Abs[x - 1] Abs[x + 1], x -> 0]", "1");

    /* Bounded oscillator times Abs[x] -- |Sin[1/x]| <= 1 squeezes to 0. */
    check_equiv("Limit[Abs[x] Sin[1/x], x -> 0]", "0");
}

/* ----------------------------------------------------------------- */
/* Special functions with a logarithmic leading-order term            */
/*                                                                    */
/* ExpIntegralEi[x] at x -> 0 has no Taylor series: its Series leading */
/* (x^0) coefficient is EulerGamma + Log[x], which itself diverges.   */
/* The series layer reads that coefficient and recurses on its limit  */
/* (every higher-order term vanishes), giving -Infinity. LogIntegral  */
/* takes the finite value 0 at the origin and must stay 0.            */
/* ----------------------------------------------------------------- */
static void test_special_log_leading(void) {
    check_equiv("Limit[ExpIntegralEi[x], x -> 0]", "-Infinity");
    check_equiv("Limit[LogIntegral[x], x -> 0]",   "0");
    /* The leading coefficient itself, resolved via finite-point
     * single-dominant termwise reduction. */
    check_equiv("Limit[EulerGamma + Log[x], x -> 0]", "-Infinity");
    /* Plain Log baseline (handled by layer_log_of_finite) still works. */
    check_equiv("Limit[Log[x], x -> 0]", "-Infinity");
}

/* ----------------------------------------------------------------- */
/* Method option                                                      */
/*                                                                    */
/* Method -> Automatic (default) runs the full cascade; a named       */
/* method restricts the top-level call to a single strategy group and */
/* leaves Limit unevaluated when that group does not apply.           */
/* ----------------------------------------------------------------- */
static void test_method(void) {
    /* Options[Limit] advertises Method -> Automatic. */
    check("MemberQ[Options[Limit], Method -> Automatic]", "True");

    /* Default is unchanged by the Method plumbing. */
    check("Limit[Sin[x]/x, x -> 0]", "1");
    check("Limit[Sin[x]/x, x -> 0, Method -> Automatic]", "1");

    /* Each named method resolves the shape it is meant for. */
    check("Limit[Sin[x]/x, x -> 0, Method -> \"Series\"]", "1");
    check("Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> \"RationalFunction\"]", "2");
    check("Limit[(x^2 - 1)/(x - 1), x -> 1, Method -> \"Substitution\"]", "2");
    check("Limit[(Exp[x] - 1)/x, x -> 0, Method -> \"LHospital\"]", "1");
    check("Limit[(1 + 1/x)^x, x -> Infinity, Method -> \"Asymptotic\"]", "E");
    check("Limit[Sin[x^2]/x, x -> Infinity, Method -> \"Bounded\"]", "0");

    /* A named method that does not apply to this shape leaves the whole
     * expression unevaluated (the RationalFunction shortcut cannot see a
     * transcendental numerator). */
    check("Limit[Sin[x]/x, x -> 0, Method -> \"RationalFunction\"]",
          "Limit[Sin[x]/x, x -> 0, Method -> \"RationalFunction\"]");

    /* An unrecognised method value is left unevaluated. */
    check("Limit[Sin[x]/x, x -> 0, Method -> \"Bogus\"]",
          "Limit[Sin[x]/x, x -> 0, Method -> \"Bogus\"]");

    /* The restriction is confined to the outermost call: the two-sided
     * probe recurses into one-sided sub-limits that run the full cascade,
     * so Series still resolves the pole at Pi/2. */
    check("Limit[Tan[x], x -> Pi/2, Method -> \"Series\"]", "ComplexInfinity");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_fast_paths);
    TEST(test_series);
    TEST(test_rational);
    TEST(test_log_reduction);
    TEST(test_rp_forms);
    TEST(test_directions);
    TEST(test_interval_returns);
    TEST(test_iterated);
    TEST(test_unevaluated);

    TEST(test_reciprocal_trig);
    TEST(test_bounded_envelope);
    TEST(test_arctan_infinity);
    TEST(test_asymptotic_at_infinity);
    TEST(test_exp_indeterminate);
    TEST(test_power_zero_folding);
    TEST(test_misc);
    TEST(test_regression_batch2);
    TEST(test_regression_batch3);
    TEST(test_wp8_numeric_point_bypass);
    TEST(test_wp2_bx_kernel);
    TEST(test_wp4_pole_sign_disagreement);
    TEST(test_wp1_multivariate);
    TEST(test_wp3_xpow_a_at_nonzero);
    TEST(test_wp6_hyperbolic_at_infinity);
    TEST(test_wp5_dominant_term_at_infinity);
    TEST(test_wp9_branch_cuts);
    TEST(test_wp7_gruntz_log_sum);
    TEST(test_abs_directional);
    TEST(test_special_log_leading);
    TEST(test_method);

    printf("All limit tests passed.\n");
    return 0;
}

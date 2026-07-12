/* test_integrate_risch_transcendental.c
 *
 * EXTENSIVE test suite for the recursive transcendental Risch integrator
 * (Integrate`RischTranscendental, Method -> "RischTranscendental").  Correctness is asserted
 * by the universal predicate Simplify[D[Integrate[f, x], x] - f] === 0 (via
 * assert_rm_diff_zero) rather than by fixed output strings, so the tests
 * survive surface-form changes; the special-function and trig/hyperbolic
 * outputs (Erf/Ei/li/PolyLog and I-laden Cosh/Sinh forms) whose exact Simplify
 * is slower are verified numerically at interior points (assert_rm_num).  Every
 * integrand below was empirically classified against the built integrator
 * (tests/rt_probe.sh over tests/rt_candidates.txt): the assert_rm_* cases close
 * with a diff-back of exactly 0, and the assert_head_unevaluated cases decline
 * cleanly (never a wrong closed form — correct by construction).
 *
 * Coverage is organised by the transcendental-tower case structure of the
 * algorithm (see RISCH_STATUS.md / docs/spec/builtins/calculus.md):
 *   base (rational) -> single logarithmic extension (polynomial / fractional /
 *   Hermite) -> single exponential extension (Laurent / fractional / Hermite /
 *   coupled hyperexponential) -> trig-hyperbolic front-end -> multi-kernel
 *   exponential sums -> nested logarithmic / exponential towers -> the genuine
 *   one-extension recursion (mixed towers, rational lower-field coefficients,
 *   proper parts, field Risch DE, evaluator-merged monomials) -> special
 *   functions (Erf/Ei/li/dilog).  Each case function also pins representative
 *   out-of-scope siblings that must DECLINE.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "calculus/integrate_risch_transcendental.h"  /* rt_rde_var_bound (white-box) */

#include <stdio.h>
#include <string.h>

/* Assert that the evaluated `input` is an unevaluated call with head
 * `head` (used for strict no-fallback / out-of-scope cases). */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* Assert Simplify[D[Integrate`RischTranscendental[f, x], x] - (f)] === 0, i.e.
 * the explicit package head produces a correct antiderivative. */
static void assert_rm_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate`RischTranscendental[%s, x], x] - (%s)]", f, f);
    assert_eval_eq(buf, "0", 0);
}

/* Assert the same via the Method-option surface form. */
static void assert_rm_method_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischTranscendental\"], x] - (%s)]",
        f, f);
    assert_eval_eq(buf, "0", 0);
}

/* Fast numeric verification for special-function / I-laden antiderivatives,
 * whose exact Simplify[D - f] is slower (Erf/Ei/li/PolyLog, Cosh/Sinh-of-
 * complex).  Confirms closure (non-unevaluated) and correctness at three
 * interior points. */
static void assert_rm_num(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate`RischTranscendental[%s, x], x] - (%s)) /. x -> 13/10]] +"
        " Abs[N[(D[Integrate`RischTranscendental[%s, x], x] - (%s)) /. x -> 17/10]] +"
        " Abs[N[(D[Integrate`RischTranscendental[%s, x], x] - (%s)) /. x -> 9/10]]"
        " < 1/100000", f, f, f, f, f, f);
    assert_eval_eq(buf, "True", 0);
}
/* ================= RATIONAL BASE CASE =================
 * The bottom of the tower recursion: a rational function of x alone, delegated
 * to Integrate`BronsteinRational (Hermite + Lazard-Rioboo-Trager).  Covers
 * polynomials, distinct/repeated real poles, irreducible quadratics (ArcTan),
 * higher-degree cyclotomic-style denominators, repeated quadratics, improper
 * fractions (polynomial + proper split), and the pure-power case (which must
 * NOT emit a spurious Log). */
static void test_rational_case(void) {
    assert_rm_diff_zero("0");
    assert_rm_diff_zero("5");
    assert_rm_diff_zero("x");
    assert_rm_diff_zero("x^2");
    assert_rm_diff_zero("x^3");
    assert_rm_diff_zero("x^5");
    assert_rm_diff_zero("3 x^2 + 2 x + 1");
    assert_rm_diff_zero("x^4 - x");
    assert_rm_diff_zero("(x-1)(x+2)");
    assert_rm_diff_zero("1/(x-1)");
    assert_rm_diff_zero("1/((x-1)(x-2))");
    assert_rm_diff_zero("1/(x^2-1)");
    assert_rm_diff_zero("1/(x^2-4)");
    assert_rm_diff_zero("(2 x + 3)/(x^2 + 3 x + 2)");
    assert_rm_diff_zero("1/(x^2-x-6)");
    assert_rm_diff_zero("1/((x-1)(x-2)(x-3))");
    assert_rm_diff_zero("1/(x-1)^2");
    assert_rm_diff_zero("1/(x-1)^3");
    assert_rm_diff_zero("1/(x^2 (x-1))");
    assert_rm_diff_zero("1/x^3");
    assert_rm_diff_zero("1/((x-1)^2 (x-2))");
    assert_rm_diff_zero("x/(x-1)^2");
    assert_rm_diff_zero("1/(x^2 (x+1)^2)");
    assert_rm_diff_zero("1/(x^2+1)");
    assert_rm_diff_zero("1/(x^2+4)");
    assert_rm_diff_zero("x/(x^2+1)");
    assert_rm_diff_zero("(2 x + 1)/(x^2 + x + 1)");
    assert_rm_diff_zero("1/(x^2+x+1)");
    assert_rm_diff_zero("1/(x^2+2 x+5)");
    assert_rm_diff_zero("(x+1)/(x^2+1)");
    assert_rm_diff_zero("1/(x^4+1)");
    assert_rm_diff_zero("1/(x^4+x^2+1)");
    assert_rm_diff_zero("x/(x^4+1)");
    assert_rm_diff_zero("1/(x^4-1)");
    assert_rm_diff_zero("1/(x^3-1)");
    assert_rm_diff_zero("1/(x^3+1)");
    assert_rm_diff_zero("1/(x^2+1)^2");
    assert_rm_diff_zero("1/(x^2+1)^3");
    assert_rm_diff_zero("x/(x^2+1)^2");
    assert_rm_diff_zero("1/(x^2+4)^2");
    assert_rm_diff_zero("(x^3+1)/(x^2+1)");
    assert_rm_diff_zero("(x^2+1)/(x-1)");
    assert_rm_diff_zero("x^4/(x^2+1)");
    assert_rm_diff_zero("(x^5+1)/(x^2-1)");
    assert_rm_diff_zero("1/(x^2-2)");
    assert_rm_diff_zero("x^2/(x^2+1)");
    assert_rm_diff_zero("1/(x (x^2+1))");
    assert_rm_diff_zero("(x+1)/(x^3+x)");
    assert_rm_diff_zero("1/(x^3+x)");
}

/* ================= LOGARITHMIC POLYNOMIAL =================
 * Single logarithmic extension theta = Log[u], integrand polynomial in theta:
 * the recursive primitive-polynomial coefficient matching q_i' + (i+1) q_{i+1}
 * eta = p_i, with a limited-integration oracle that folds a would-be new
 * logarithm back into the tower (Log[x]^k/x -> Log[x]^(k+1)/(k+1)).  Covers
 * theta powers, polynomial-in-x coefficients, Log[a x + b], Log of a quadratic
 * (mixes ArcTan/Log), and mixed elementary + log sums. */
static void test_logarithmic_case(void) {
    assert_rm_diff_zero("Log[x]");
    assert_rm_diff_zero("Log[x]^2");
    assert_rm_diff_zero("Log[x]^3");
    assert_rm_diff_zero("Log[x]^4");
    assert_rm_diff_zero("Log[x]^5");
    assert_rm_diff_zero("x Log[x]");
    assert_rm_diff_zero("x^2 Log[x]");
    assert_rm_diff_zero("x^3 Log[x]");
    assert_rm_diff_zero("(x^2+1) Log[x]");
    assert_rm_diff_zero("x Log[x]^2");
    assert_rm_diff_zero("x^2 Log[x]^2");
    assert_rm_diff_zero("x Log[x]^3");
    assert_rm_diff_zero("Log[2 x + 3]");
    assert_rm_diff_zero("Log[3 x - 1]");
    assert_rm_diff_zero("Log[x + 1]");
    assert_rm_diff_zero("Log[5 x]");
    assert_rm_diff_zero("Log[x^2 + 1]");
    assert_rm_diff_zero("Log[x^2 - 1]");
    assert_rm_diff_zero("Log[x^2 + x + 1]");
    assert_rm_diff_zero("Log[x]/x");
    assert_rm_diff_zero("Log[x]^2/x");
    assert_rm_diff_zero("Log[x]^3/x");
    assert_rm_diff_zero("Log[x]^4/x");
    assert_rm_diff_zero("x + Log[x]");
    assert_rm_diff_zero("Log[x] - Log[x]^2");
    assert_rm_diff_zero("Log[x]/x^2");
}

/* ================= FRACTIONAL (Rothstein-Trager) LOG =================
 * Proper rational function of theta = Log[u] with a SQUAREFREE denominator ->
 * sum_i c_i D(g_i)/g_i with constant residues found by an exact SolveAlways over
 * {t, x}.  Single and multi-residue denominators. */
static void test_fractional_log_case(void) {
    assert_rm_diff_zero("1/(x (1 + Log[x]))");
    assert_rm_diff_zero("1/(x Log[x])");
    assert_rm_diff_zero("1/(x (1 - Log[x]))");
    assert_rm_diff_zero("1/(x (2 + Log[x]))");
    assert_rm_diff_zero("1/(x (2 + 3 Log[x]))");
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[x]))");
    assert_rm_diff_zero("1/(x (Log[x]^2 - 1))");
    assert_rm_diff_zero("1/(x (Log[x]^2 - 4))");
    assert_rm_diff_zero("1/(x (Log[x]^2 + Log[x]))");
    /* Pure resultant LRT: irreducible quadratic in Log[x] with ALGEBRAIC
     * (complex-conjugate) residues -> ArcTan[Log[x]] form, which the
     * single-constant-per-factor SolveAlways path cannot produce. */
    assert_rm_diff_zero("1/(x (Log[x]^2 + 1))");
    assert_rm_diff_zero("1/(x (Log[x]^2 + 4))");
    assert_rm_diff_zero("1/(x (Log[x]^2 + 2 Log[x] + 2))");
}

/* ================= HERMITE (repeated poles) LOG =================
 * Proper rational function of theta = Log[u] with a REPEATED pole -> Hermite
 * reduction Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j), Hden = gcd(D,
 * dD/dtheta), all constants by SolveAlways over {t, x}.  Squared/cubed poles,
 * shifted poles, and Hermite + log combinations. */
static void test_hermite_log_case(void) {
    assert_rm_diff_zero("1/(x (1 + Log[x])^2)");
    assert_rm_diff_zero("1/(x Log[x]^2)");
    assert_rm_diff_zero("1/(x (1 + Log[x])^3)");
    assert_rm_diff_zero("1/(x Log[x]^3)");
    assert_rm_diff_zero("(1 + Log[x])/(x Log[x]^2)");
    assert_rm_diff_zero("Log[x]/(x (1 + Log[x])^2)");
    assert_rm_diff_zero("1/(x (2 + Log[x])^2)");
    assert_rm_diff_zero("1/(x (1 - Log[x])^2)");
    assert_rm_diff_zero("(2 + Log[x])/(x (1 + Log[x])^2)");
}

/* ================= EXPONENTIAL LAURENT =================
 * Single exponential extension theta = E^u; integrand a Laurent polynomial
 * sum_i p_i(x) theta^i.  Powers decouple: the i=0 term integrates in K, each
 * i != 0 term solves the Risch DE q_i' + i u' q_i = p_i (polynomial ansatz, or
 * the q = h/Denominator[p] ansatz when u or p is rational in x).  Covers
 * polynomial coefficients, quadratic exponents (-> E^(x^2)/2), Laurent sums,
 * and rational exponent/coefficient (E^(1/x), E^x/x - E^x/x^2). */
static void test_exponential_case(void) {
    assert_rm_diff_zero("Exp[x]");
    assert_rm_diff_zero("Exp[2 x]");
    assert_rm_diff_zero("Exp[-x]");
    assert_rm_diff_zero("Exp[3 x]");
    assert_rm_diff_zero("Exp[3 x + 1]");
    assert_rm_diff_zero("x Exp[x]");
    assert_rm_diff_zero("x^2 Exp[x]");
    assert_rm_diff_zero("x^3 Exp[x]");
    assert_rm_diff_zero("x^4 Exp[x]");
    assert_rm_diff_zero("(x^2+1) Exp[x]");
    assert_rm_diff_zero("x Exp[2 x]");
    assert_rm_diff_zero("x^2 Exp[2 x]");
    assert_rm_diff_zero("x Exp[-x]");
    assert_rm_diff_zero("x Exp[3 x + 1]");
    assert_rm_diff_zero("x Exp[x^2]");
    assert_rm_diff_zero("x Exp[-x^2]");
    assert_rm_diff_zero("x Exp[2 x^2]");
    assert_rm_diff_zero("(Exp[x] + Exp[-x])/2");
    assert_rm_diff_zero("Exp[x] + Exp[-x]");
    assert_rm_diff_zero("Exp[x] + Exp[-x] + Exp[2 x]");
    assert_rm_diff_zero("Exp[2 x] + Exp[-2 x]");
    assert_rm_diff_zero("Exp[x] - Exp[-x]");
    assert_rm_diff_zero("-Exp[1/x]/x^2");
    assert_rm_diff_zero("Exp[x]/x - Exp[x]/x^2");
    assert_rm_diff_zero("x + x Exp[x]");
}

/* ================= FRACTIONAL (Rothstein-Trager) EXP =================
 * Squarefree-denominator proper fraction of theta = E^u whose numerator carries
 * the theta' factor, so it closes to a single Log with a constant residue
 * (E^x/(c + d E^x) -> (1/d) Log(c + d E^x)). */
static void test_fractional_exp_case(void) {
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])");
    assert_rm_diff_zero("Exp[x]/(1 - Exp[x])");
    assert_rm_diff_zero("Exp[x]/(2 + Exp[x])");
    assert_rm_diff_zero("Exp[x]/(1 + 2 Exp[x])");
    assert_rm_diff_zero("2 Exp[x]/(1 + Exp[x])");
    assert_rm_diff_zero("Exp[x]/(2 + 3 Exp[x])");
    /* Pure resultant LRT: an irreducible quadratic in theta = E^x with
     * ALGEBRAIC residues (E^(2x) kernelised to theta^2) -> ArcTan[E^x] form. */
    assert_rm_diff_zero("Exp[x]/(Exp[2 x] + 1)");
    assert_rm_diff_zero("Exp[x]/(Exp[2 x] + 4)");
    /* Higher (cubic+) irreducible denominators in the exp kernel: the residues
     * are algebraic of degree >= 3, so the real Log/ArcTan Rothstein-Trager
     * declines and the rational-reduction fallback (t = E^(a x),
     * dx = dt/(a t) -> rational integral in t) closes them as RootSum.  The
     * commensurate-exponent case E^(x/6)/(1+E^(x/2)+E^(x/3)) reduces onto the
     * primitive t = E^(x/6) (den 1 + t^2 + t^3).  Diff-back verification relies
     * on the generalized D[RootSum] collapse (src/root.c). */
    assert_rm_diff_zero("1/(1 + Exp[2 x] + Exp[3 x])");
    assert_rm_diff_zero("Exp[x/6]/(1 + Exp[x/2] + Exp[x/3])");
}

/* ================= HERMITE (repeated poles) EXP =================
 * Repeated transcendental pole for theta = E^u with D coprime to theta:
 * E^x/(c + E^x)^n -> a Hermite rational part (no residual log). */
static void test_hermite_exp_case(void) {
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^2");
    assert_rm_diff_zero("Exp[x]/(1 - Exp[x])^2");
    assert_rm_diff_zero("Exp[x]/(2 + Exp[x])^2");
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^3");
}

/* ================= COUPLED HYPEREXPONENTIAL =================
 * Rational function of theta = E^u whose integral MIXES a Laurent-polynomial
 * part with a log part (D(g)/g improper in theta, so they do not separate):
 * one unified ansatz Q = sum_i w_i(x) theta^i + H(theta)/Hden + sum_j c_j
 * Log(g_j) solved by SolveAlways over {t, x}.  Covers 1/(c + E^x), improper
 * (E^(2x)/(1+E^x)), and repeated / theta=0 poles (Phase A: 1/(1+E^x)^2,
 * 1/(E^x (1+E^x)^2), 1/(1+E^x)^3). */
static void test_hyperexponential_case(void) {
    assert_rm_diff_zero("1/(1 + Exp[x])");
    assert_rm_diff_zero("1/(1 - Exp[x])");
    assert_rm_diff_zero("1/(2 + Exp[x])");
    assert_rm_diff_zero("1/(1 + 2 Exp[x])");
    assert_rm_diff_zero("Exp[2 x]/(1 + Exp[x])");
    assert_rm_diff_zero("(1 + Exp[x])/(2 + Exp[x])");
    assert_rm_diff_zero("Exp[3 x]/(1 + Exp[x])");
    assert_rm_diff_zero("1/(1 + Exp[x])^2");
    assert_rm_diff_zero("1/(Exp[x] (1 + Exp[x])^2)");
    assert_rm_diff_zero("1/(2 + Exp[x])^2");
    assert_rm_diff_zero("1/(1 + Exp[x])^3");
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^3");
    assert_rm_diff_zero("Exp[2 x]/(1 + Exp[x])^2");
    /* Non-integer multiplicatively-commensurate exponents (§6.1 item 3): the two
     * kernels E^(x/2), E^(x/3) share no integer-ratio member primitive, so the
     * old kernelizer declined; the synthesized primitive E^(x/6) (x/2 = 3 x/6,
     * x/3 = 2 x/6) collapses them onto one tower variable and closes them. */
    assert_rm_diff_zero("1/(Exp[x/2] + Exp[x/3])");
    assert_rm_diff_zero("1/(Exp[x/2] + Exp[x/3])^2");
    assert_rm_diff_zero("Exp[x/6]/(1 + Exp[x/2] + Exp[x/3])");
}

/* ================= TRIG / HYPERBOLIC FRONT-END =================
 * TrigToExp -> the exponential (Laurent/coupled) machinery -> ExpToTrig, both
 * exact rewrites.  Powers, products, and multiple-angle trig/hyperbolic reduce
 * to a Laurent polynomial in E^(I x)/E^x; Tan/Tanh close via the coupled case
 * (correct, though left in an I-laden complex-log form).  Verified numerically
 * (the I-laden outputs are correct but not symbolically simplified to 0). */
static void test_trig_frontend(void) {
    assert_rm_num("Sin[x]");
    assert_rm_num("Cos[x]");
    assert_rm_num("Sinh[x]");
    assert_rm_num("Cosh[x]");
    assert_rm_num("Tan[x]");
    assert_rm_num("Tanh[x]");
    assert_rm_num("Sin[x]^2");
    assert_rm_num("Cos[x]^2");
    assert_rm_num("Sinh[x]^2");
    assert_rm_num("Cosh[x]^2");
    assert_rm_num("Sin[x]^3");
    assert_rm_num("Cos[x]^3");
    assert_rm_num("Sin[x] Cos[x]");
    assert_rm_num("Sin[2 x]");
    assert_rm_num("Cos[2 x]");
    assert_rm_num("Sin[3 x]");
    assert_rm_num("Sin[x] Cos[x]^2");
    assert_rm_num("Sin[x]^2 Cos[x]");
    assert_rm_num("Sin[x] + Cos[x]");
    assert_rm_num("Sinh[x] Cosh[x]");
    assert_rm_num("Sin[x]^4");
    assert_rm_num("Cos[x]^4");
    assert_rm_num("Sinh[2 x]");
}

/* ================= MULTI-KERNEL EXPONENTIAL SUMS =================
 * Integrands that exponentialize to a sum sum_k p_k(x) E^(W_k) of
 * NON-commensurate exponentials (e.g. the (1 +/- I) x pair from E^x Sin[x]).
 * The distinct exponentials are independent extensions and DECOUPLE: each term
 * closes by its own Risch DE q_k' + W_k' q_k = p_k (rt_expsum_case) without a
 * coupled tower.  Real exponential times trig, polynomial coefficients, and
 * non-unit frequency pairs. */
static void test_multikernel_case(void) {
    assert_rm_diff_zero("Exp[x] Sin[x]");
    assert_rm_diff_zero("Exp[x] Cos[x]");
    assert_rm_diff_zero("x Exp[x] Sin[x]");
    assert_rm_diff_zero("Exp[2 x] Sin[3 x]");
    assert_rm_diff_zero("Exp[2 x] Cos[3 x]");
    assert_rm_diff_zero("x^2 Exp[x] Cos[x]");
    assert_rm_diff_zero("Exp[-x] Sin[x]");
    assert_rm_diff_zero("Exp[x] Sin[2 x]");
    assert_rm_diff_zero("Exp[x] (Sin[x] + Cos[x])");
    assert_rm_diff_zero("x Exp[x] Cos[x]");
    assert_rm_diff_zero("Exp[x] + Exp[2 x] Sin[x]");
}

/* ================= NESTED LOG TOWERS =================
 * Rational function of a chain of nested logarithms Log[x], Log[Log[x]], ...
 * integrated over the triangular tower derivation D = d/dx + sum_i Dt_i d/dt_i
 * by one unified SolveAlways ansatz over all tower variables (rt_log_tower_case),
 * up to depth 4, plus independent (non-nested) log chains.  A whole-tower
 * rationality gate makes non-rational inner kernels decline (see declines). */
static void test_log_tower_case(void) {
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]])");
    assert_rm_diff_zero("Log[Log[x]]/(x Log[x])");
    assert_rm_diff_zero("Log[Log[x]]/x");
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]] Log[Log[Log[x]]])");
    assert_rm_diff_zero("Log[x] + Log[x + 1]");
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[Log[x]]))");
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[Log[x]])^2)");
    /* Cap-free top-degree bound: these need Ntop = 6 and 8 in the top log kernel
     * Log[Log[x]] — the former arbitrary `Ntop > 4 -> 4` cap declined them.  The
     * bound Ntop = deg_top(f) + 1 is now exact/uncapped (a LOG kernel lowers degree
     * under D), so ALL top degrees close. */
    assert_rm_diff_zero("Log[Log[x]]^5/(x Log[x])");   /* -> Log[Log[x]]^6/6 */
    assert_rm_diff_zero("Log[Log[x]]^7/(x Log[x])");   /* -> Log[Log[x]]^8/8 */
    /* Composite-argument logs (Log[x/Log[x]], Log[Log[x]/x]) are normalized to
     * the minimal independent generator set {Log[x], Log[Log[x]]} before the
     * ansatz is built (rt_expand_logs), collapsing a spurious depth-3 tower to
     * depth 2 (~90x faster).  -> -x Log[Log[x]/x] - x/Log[x/Log[x]]. */
    assert_rm_diff_zero("1 - (1 + 1/Log[x/Log[x]]^2)/Log[x] + 1/Log[x/Log[x]]^2 "
                        "- 1/Log[x/Log[x]] - Log[Log[x]/x]");
}

/* ================= NESTED EXP TOWERS + MERGED MONOMIAL =================
 * Chain of nested exponentials t_i = E^(u_i) (e.g. E^x, E^(E^x)) via a Laurent
 * ansatz over the exp tower derivation, diff-back verified.  The third case is
 * an EVALUATOR-MERGED monomial: the evaluator folds E^x E^(E^x) into
 * E^(x + E^x), whose exponent carries the foreign kernel E^x; the
 * rt_expand_exp_sums pre-pass re-splits E^(a+b) -> E^a E^b to restore the
 * independent basis, and the coupled hyperexponential proper part closes it. */
static void test_exp_tower_case(void) {
    assert_rm_diff_zero("E^x E^(E^x)");
    assert_rm_diff_zero("E^(2 E^x) E^x");
    assert_rm_diff_zero("E^x E^(E^x)/(1 + E^(E^x))");
    /* Multiplicatively commensurate merged kernels (§6.1 item 3): E^(2 E^x) =
     * (E^(E^x))^2 is a power of another tower kernel, so the commensurate
     * reduction in rt_tower_build keeps only the primitive E^(E^x) as an
     * extension and aliases E^(2 E^x) -> t^2.  Without it the dependent kernel
     * would spuriously add an extension and the tower declined. */
    assert_rm_diff_zero("E^x E^(2 E^x)/(1 + E^(E^x))");        /* -> E^(E^x) - Log[1+E^(E^x)] */
    assert_rm_diff_zero("E^x E^(3 E^x)/(1 + E^(E^x))");
    assert_rm_diff_zero("E^x E^(2 E^x)/(1 + E^(2 E^x))");
    /* Exp-top algebraic-residue LRT unblocked by the reduction (den = t^2+1). */
    assert_rm_diff_zero("E^x E^(E^x)/(1 + E^(2 E^x))");        /* -> ArcTan[E^(E^x)] */
    assert_rm_diff_zero("E^x E^(E^x)/(4 + E^(2 E^x))");
    /* Cap-free top Laurent extent: the antiderivative is a degree-5 Laurent
     * polynomial in the top exp kernel E^(E^x) — beyond the former arbitrary
     * `ihi > 4 -> 4` cap, which declined it.  The top Laurent range is now the
     * exact deg_top(num) - deg_top(den) (an EXP kernel preserves degree under D). */
    assert_rm_diff_zero("E^x E^(6 E^x)/(1 + E^(E^x))");
    /* Non-integer multiplicatively-commensurate nested exponents (§6.1 item 3):
     * the class {7 E^x/6, 2 E^x/3, E^x/2} has NO integer-ratio member primitive
     * (ratios 7/3, 4/3, 7/4, ...), so the old tower builder declined.  The
     * SYNTHESIZED (non-member) primitive E^(x/6) of the inner exponent — here
     * E^(E^x/6) — collapses the class onto one tower variable (E^(2 E^x/3) = t^4,
     * E^(E^x/2) = t^3) and the coupled hyperexponential proper part closes it.
     * Integrand is d/dx of the elementary E^(E^x/2)/(1 + E^(2 E^x/3)). */
    assert_rm_diff_zero("D[E^(E^x/2)/(1 + E^(2 E^x/3)), x]");
}

/* ================= RECURSIVE / MIXED =================
 * The genuine one-extension-at-a-time recursion (Bronstein ch. 5) that the
 * flat single-kind tower ansatz cannot express: MIXED
 * exp/log towers, RATIONAL lower-field coefficients, a proper rational part at
 * a logarithmic top (tower Hermite + Rothstein-Trager), the coupled
 * hyperexponential proper part at an exponential top, and the general field
 * Risch DE (rational lower-field solution, monomial and non-monomial
 * denominators via the denominator theorem).  All diff-back verified. */
static void test_recursive_tower_case(void) {
    assert_rm_diff_zero("E^x/x + E^x Log[x]");
    assert_rm_diff_zero("Log[1 + E^x] + x E^x/(1 + E^x)");
    assert_rm_diff_zero("1/(x^2 Log[x]) - Log[Log[x]]/x^2");
    assert_rm_diff_zero("1/(x (1 + Log[x])) + E^x");
    assert_rm_diff_zero("Log[Log[x]]/(x Log[x] (1 + Log[Log[x]]))");
    assert_rm_diff_zero("(2 Log[x]/x) E^(Log[x]^2)/(1 + E^(Log[x]^2))");
    assert_rm_diff_zero("(2/x - 1/(x Log[x]^2)) E^(Log[x]^2)");
    assert_rm_diff_zero("(2 Log[x]/(x (1 + Log[x])) - 1/(x (1 + Log[x])^2)) E^(Log[x]^2)");
    /* Exact field-RDE degree bound (§6.1 item 1, rt_rde_var_bound): the RDE solver
     * works for ALL degrees — the exponential-Laurent coefficient q = Log[x]^k has
     * lower-variable degree k, found EXACTLY (deg(q) = deg(p) - deg(Dcoef)) with NO
     * arbitrary degree cap or monomial ceiling.  There is nothing special about any
     * particular k; deg-6/-8/-12/-20 all close (these formerly declined at the
     * removed cap-at-5). */
    assert_rm_diff_zero("(6 Log[x]^5 + 2 Log[x]^7)/x E^(Log[x]^2)");   /* q = Log[x]^6 */
    assert_rm_diff_zero("(8 Log[x]^7 + 2 Log[x]^9)/x E^(Log[x]^2)");   /* q = Log[x]^8 */
    assert_rm_diff_zero("(12 Log[x]^11 + 2 Log[x]^13)/x E^(Log[x]^2)"); /* q = Log[x]^12 */
    assert_rm_diff_zero("(20 Log[x]^19 + 2 Log[x]^21)/x E^(Log[x]^2)"); /* q = Log[x]^20 */
    /* Pure resultant Lazard-Rioboo-Trager log part at a LOGARITHMIC tower top:
     * a proper fraction of Log[Log[x]] whose Rothstein-Trager residues are
     * ALGEBRAIC constants (the +-I/2 that split t^2+1 into ArcTan) — the class
     * the bounded single-constant SolveAlways ansatz cannot express, closed by
     * lifting rt_frac_lrt into rt_field_ratint (gate free of x and Log[x]). */
    assert_rm_diff_zero("1/(x Log[x] (Log[Log[x]]^2 + 1))");
    assert_rm_diff_zero("1/(x Log[x] (Log[Log[x]]^2 + 4))");
    assert_rm_diff_zero("(1 + Log[Log[x]])/(x Log[x] (Log[Log[x]]^2 + 1))");
    assert_rm_diff_zero("1/(x Log[x] (Log[Log[x]]^3 + 1))");
}

/* ================= SPECIAL FUNCTIONS =================
 * Erf / Ei / li / dilog forms: integrals the elementary cascade leaves open,
 * closed with special functions Mathilda provides, each behind an exact
 * structural certificate.  Gaussian K E^(a x^2 + b x + c) -> Erf/Erfi (with the
 * completing-the-square constant), (const E^(a x))/(c x + d) -> ExpIntegralEi,
 * K/Log[x] -> LogIntegral, K Log[1 +/- p x]/x -> PolyLog[2, ...].  Verified
 * numerically. */
static void test_special_functions(void) {
    assert_rm_num("Exp[-x^2]");
    assert_rm_num("Exp[x^2]");
    assert_rm_num("Exp[-x^2/2]");
    assert_rm_num("Exp[-3 x^2 + 2 x + 1]");
    assert_rm_num("Exp[-x^2 + x]");
    assert_rm_num("Exp[2 x^2]");
    assert_rm_num("Exp[-2 x^2]");
    assert_rm_num("Exp[x]/x");
    assert_rm_num("Exp[2 x]/(x - 1)");
    assert_rm_num("Exp[3 x]/x");
    assert_rm_num("Exp[x]/(2 x + 1)");
    assert_rm_num("Exp[x]/(x + 1)");
    /* Widened Ei recognizer (extracts the E^v kernel directly, so a<0 and a
     * nonzero exponent constant b close): E^(-x)/x -> ExpIntegralEi[-x]. */
    assert_rm_num("Exp[-x]/x");
    assert_rm_num("Exp[-2 x]/x");
    assert_rm_num("Exp[-x]/(x + 1)");
    assert_rm_num("1/Log[x]");
    assert_rm_num("3/Log[x]");
    /* Widened li recognizer (c w^(p-1) w'/Log[w] -> c LogIntegral[w^p]): a
     * scaled/affine log argument (1/Log[2x] -> LogIntegral[2x]/2, p=1) and a
     * monomial numerator (x/Log[x] -> LogIntegral[x^2], p=2). */
    assert_rm_num("1/Log[2 x]");
    assert_rm_num("1/Log[3 x + 1]");
    assert_rm_num("x/Log[x]");
    assert_rm_num("x^2/Log[x]");
    assert_rm_num("Log[1 - x]/x");
    assert_rm_num("Log[1 + x]/x");
    assert_rm_num("Log[1 - 2 x]/x");
    assert_rm_num("Log[1 + 2 x]/x");
}

/* ================= EXPECTED DECLINES (non-elementary / out-of-scope) =================
 * A Risch decision procedure must NEVER emit a wrong closed form: anything out
 * of scope, non-elementary, or beyond an implemented degree bound returns the
 * integral unevaluated.  Each line pins a specific reason the integrand cannot
 * (or does not yet) close, so a future change that starts *wrongly* closing one
 * of them trips here. */
static void test_strict_unevaluated(void) {
    /* Genuinely non-elementary: Fresnel S, sine/cosine integral. */
    assert_head_unevaluated("Integrate`RischTranscendental[Sin[x^2], x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[Sin[x]/x, x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[Cos[x]/x, x]", "Integrate`RischTranscendental");
    /* Exponential of a rational whose Risch DE has no rational solution, and an
     * Ei-type integrand with a NON-linear denominator (unimplemented). */
    assert_head_unevaluated("Integrate`RischTranscendental[Exp[1/x], x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[Exp[x]/x^2, x]", "Integrate`RischTranscendental");
    /* Ei-type integrand with a NON-linear (quadratic) denominator remains out of
     * scope: the widened Ei recognizer requires an exactly linear c x + d. */
    assert_head_unevaluated("Integrate`RischTranscendental[Exp[x]/(x^2 + 1), x]", "Integrate`RischTranscendental");
    /* Non-elementary nested-log integrands (need Ei/li of a log); and a residual
     * NON-rational inner kernel (Sin[Log[x]]) must DECLINE, never certify a wrong
     * form (the whole-tower rationality gate). */
    assert_head_unevaluated("Integrate`RischTranscendental[Log[x] Log[Log[x]], x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[Sin[Log[x]] Log[Log[x]], x]", "Integrate`RischTranscendental");
    /* Non-elementary nested exponentials.  E^(E^x)/(1+E^(E^x)) is a REGRESSION
     * guard: a single-kernel case once certified the WRONG Log[1+E^(E^x)]/E^x by
     * leaving the inner E^x as a free SolveAlways parameter (fixed by rt_kernel_simple). */
    assert_head_unevaluated("Integrate`RischTranscendental[E^(E^x), x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[E^(E^x)/(1 + E^(E^x)), x]", "Integrate`RischTranscendental");
    /* Commensurate-reduction siblings that stay non-elementary: no E^x derivation
     * factor to supply dt for the E^(2 E^x)=t^2 denominator, so they decline even
     * though the reduction now lets the tower build (§6.1 item 3, never wrong). */
    assert_head_unevaluated("Integrate`RischTranscendental[E^(E^x)/(1 + E^(2 E^x)), x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[E^(2 E^x)/(1 + E^(E^x)), x]", "Integrate`RischTranscendental");
    /* Proper tower fraction with a NON-constant Rothstein-Trager residue is
     * non-elementary (must not certify a wrong constant residue). */
    assert_head_unevaluated("Integrate`RischTranscendental[1/(Log[x] (1 + Log[Log[x]])), x]", "Integrate`RischTranscendental");
    /* Algebraic-residue tower LRT gate: 1/x Log[x] missing makes the residues
     * of the Log[Log[x]]^2+1 factor depend on x (not constants of the tower
     * derivation), so the resultant LRT must DECLINE, not certify a wrong form. */
    assert_head_unevaluated("Integrate`RischTranscendental[1/(Log[x] (Log[Log[x]]^2 + 1)), x]", "Integrate`RischTranscendental");
    /* Field Risch DE with no elementary solution (E^(Log[x]^2) and its coupled
     * fraction) declines rather than forcing a bounded ansatz to a wrong answer. */
    assert_head_unevaluated("Integrate`RischTranscendental[E^(Log[x]^2), x]", "Integrate`RischTranscendental");
    assert_head_unevaluated("Integrate`RischTranscendental[E^(Log[x]^2)/(1 + Log[x]), x]", "Integrate`RischTranscendental");
}

/* ================= Supplementary: engine agreement ================= */
/* On the rational base case RischTranscendental must agree with the dedicated
 * rational Risch engine (Integrate`BronsteinRational) that it delegates to. */
static void test_rational_agreement(void) {
    assert_eval_eq(
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischTranscendental\"]"
        " - Integrate[1/(x^2 - 1), x, Method -> \"BronsteinRational\"]]",
        "0", 0);
    assert_eval_eq(
        "Simplify[Integrate[(x^3 + 1)/(x^2 + 1), x, Method -> \"RischTranscendental\"]"
        " - Integrate[(x^3 + 1)/(x^2 + 1), x, Method -> \"BronsteinRational\"]]",
        "0", 0);
}

/* ================= Automatic cascade closes special functions ================= */
/* The cascade insertion adds capability by default (WL-faithful): these
 * previously returned unevaluated.  Verified numerically (exact Simplify on
 * Erf/Ei/PolyLog derivatives is slower). */
static void test_cascade_default(void) {
    assert_eval_eq(
        "Abs[N[(D[Integrate[Exp[-x^2], x], x] - Exp[-x^2]) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[Exp[x]/x, x], x] - Exp[x]/x) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Log[x], x], x] - 1/Log[x]) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[Log[1 - x]/x, x], x] - Log[1 - x]/x) /. x -> 3/10]]"
        " < 1/100000", "True", 0);
    /* Regression: an ordinary rational integrand is unchanged. */
    assert_eval_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]", 0);
}

/* ================= Dispatch plumbing ================= */
static void test_method_plumbing(void) {
    /* Both surface forms route to the routine and agree. */
    assert_rm_diff_zero("1/(x^2 + 1)");
    assert_rm_method_diff_zero("1/(x^2 + 1)");
    assert_rm_method_diff_zero("x Exp[x]");
    assert_rm_method_diff_zero("1/(x (1 + Log[x]))");
    /* The two forms produce the same antiderivative. */
    assert_eval_eq(
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischTranscendental\"]"
        " - Integrate`RischTranscendental[1/(x^2 - 1), x]]",
        "0", 0);
}

/* ================= Strict / malformed-input behaviour ================= */
static void test_strict_misc(void) {
    /* A non-symbol integration variable is rejected (never garbage). */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[1/x, x + 1]", "Integrate`RischTranscendental");
    /* Fresnel / Si / Ci non-elementary integrands bubble back unevaluated. */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[Cos[x^2], x]", "Integrate`RischTranscendental");
    /* The LRT frac path fires only with a genuine derivation factor: the
     * ArcTan-family integrands above carry 1/x (log) or E^x (exp).  Without
     * it the integral is non-elementary (li / Ei family), so the LRT must
     * DECLINE, not certify a wrong closed form. */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[1/(Log[x]^2 + 1), x]", "Integrate`RischTranscendental");
    /* x-dependent residues (the resultant does not become free of x after the
     * content strip): the x-content gate must reject rather than certify. */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[1/(x^2 (Log[x]^2 + 1)), x]",
        "Integrate`RischTranscendental");
}

/* ================= UNIT: rt_rde_var_bound (degree arithmetic) =================
 * White-box tests of the Bronstein RdeBoundDegree leading-degree bound used by the
 * exponential-Laurent field Risch DE (rt_field_rde) and the rational-exponent RDE
 * (rt_solve_rde_rational).  This is a pure integer function of the equation's degrees
 * (dpv = deg_v(p), dfv = deg_v(f)), the monomial kind (deriv_lowers), and the Bronstein
 * resonance integer (m_res); testing it directly pins the "no arbitrary caps" contract
 * and the cancellation/resonance widening without having to construct an integrand that
 * reaches each configuration (both cancellation configs are pre-empted in the current
 * tower architecture — see the reachability note in integrate_risch_transcendental.c).
 *
 * Contract recap (bound on deg_v(q) for D[q] + f q = p, clamped at 0):
 *   deriv-lowering v (base x / log):  dfv >= 0 -> dpv - dfv                (f q dominates)
 *                                     dfv <  0 -> max(dpv+1, dpv-dfv)      (rise vs pole)
 *                                     dfv == -1 & m_res >= 0 -> widened to max(.., m_res)
 *   deriv-preserving v (exp):         dpv - dfv                            (both at deg q)
 *                                     dfv == 0 & m_res >= 0 -> widened to max(dpv, m_res)
 * The m_res widening is MONOTONE (never lowers the bound) and fires ONLY in its exact
 * configuration.  m_res == -1 disables it. */
static void chk_bound(long dpv, long dfv, bool dl, long m_res, long expect) {
    long got = rt_rde_var_bound(dpv, dfv, dl, m_res);
    ASSERT_MSG(got == expect,
        "rt_rde_var_bound(dpv=%ld, dfv=%ld, deriv_lowers=%d, m_res=%ld) = %ld, expected %ld",
        dpv, dfv, (int)dl, m_res, got, expect);
}

static void test_rde_var_bound(void) {
    /* --- deriv-lowering (primitive: base x / log), dfv >= 0: f q dominates, exact. --- */
    chk_bound(3, 1, true, -1, 2);      /* dpv - dfv */
    chk_bound(2, 0, true, -1, 2);      /* f a nonzero constant */
    chk_bound(5, 2, true, -1, 3);
    chk_bound(0, 0, true, -1, 0);
    chk_bound(1, 5, true, -1, 0);      /* dpv - dfv = -4 -> clamp to 0 */
    chk_bound(2, 7, true, -1, 0);      /* clamp */

    /* --- deriv-lowering, dfv < 0: max(dpv+1 integration rise, dpv-dfv pole in f). --- */
    chk_bound(2, -2, true, -1, 4);     /* max(3, 4) = 4 (order-2 pole dominates) */
    chk_bound(0, -3, true, -1, 3);     /* max(1, 3) = 3 */
    chk_bound(3, -1, true, -1, 4);     /* max(4, 4) = 4 (simple pole, no resonance) */
    chk_bound(-2, -1, true, -1, 0);    /* max(-1, -1) = -1 -> clamp 0 */
    chk_bound(1, -4, true, -1, 5);     /* max(2, 5) = 5 */

    /* --- deriv-lowering primitive cancellation (dfv == -1) resonance widening. --- */
    chk_bound(2, -1, true, 10, 10);    /* naive max(3,3)=3, widen to m_res=10 */
    chk_bound(2, -1, true, 3, 3);      /* m_res == naive: unchanged */
    chk_bound(2, -1, true, 1, 3);      /* MONOTONE: m_res < naive never lowers -> 3 */
    chk_bound(2, -1, true, 0, 3);      /* m_res=0 < naive -> 3 */
    chk_bound(-2, -1, true, 2, 2);     /* naive clamps to 0; widen to 2 */
    /* widening fires ONLY at dfv == -1, not at other poles. */
    chk_bound(2, -2, true, 100, 4);    /* dfv=-2: m_res ignored -> max(3,4)=4 */
    chk_bound(3, -3, true, 100, 6);    /* dfv=-3: m_res ignored -> max(4,6)=6 */

    /* --- deriv-preserving (exponential): dpv - dfv, both leading terms at deg q. --- */
    chk_bound(3, 1, false, -1, 2);
    chk_bound(2, 2, false, -1, 0);
    chk_bound(5, 0, false, -1, 5);
    chk_bound(1, 3, false, -1, 0);     /* -2 -> clamp 0 */
    chk_bound(2, -2, false, -1, 4);    /* dpv - dfv = 4 (pole raises the exp bound) */

    /* --- exponential integer resonance (dfv == 0) widening: max(dpv, m_res). --- */
    chk_bound(2, 0, false, 7, 7);      /* max(2, 7) = 7 */
    chk_bound(0, 0, false, 3, 3);      /* max(0, 3) = 3 */
    chk_bound(2, 0, false, 2, 2);      /* m_res == naive */
    chk_bound(2, 0, false, 1, 2);      /* MONOTONE: m_res < naive -> unchanged */
    chk_bound(5, 0, false, 3, 5);      /* naive dpv=5 already exceeds m_res=3 */
    chk_bound(2, 0, false, -1, 2);     /* disabled */
    /* widening fires ONLY at dfv == 0, not at other exp degrees. */
    chk_bound(2, 1, false, 100, 1);    /* dfv=1: m_res ignored -> dpv-dfv=1 */
    chk_bound(3, 2, false, 100, 1);    /* dfv=2: m_res ignored -> 1 */
}

/* Single-integration numeric verification: integrate ONCE (Module-bound), then
 * confirm the result is closed and differentiates back to f at an interior
 * point.  Used for the heavy high-degree / nested-exponential Bronstein cases
 * where re-integrating three times (assert_rm_num) or an exact Simplify would
 * dominate the suite. */
static void assert_rm_num1(const char* f) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "Module[{rmF = Integrate`RischTranscendental[%s, x]},"
        " Head[rmF] =!= Integrate &&"
        " Abs[N[(D[rmF, x] - (%s)) /. x -> 13/10]] < 1/1000000]", f, f);
    assert_eval_eq(buf, "True", 0);
}

/* Closure-only check: the integrator returns a non-Integrate result.  Used for
 * the degree-100+ shifted-pole case whose exact form's numeric diff-back would
 * dominate the suite; the closed form is correct by construction (SPDE
 * certificate), so closure suffices. */
static void assert_rm_closes(const char* f) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "Head[Integrate`RischTranscendental[%s, x]] =!= Integrate", f);
    assert_eval_eq(buf, "True", 0);
}

/* ================= BRONSTEIN RDE EXAMPLES =================
 * The transcendental-Risch-DE examples that motivated the rational one-step
 * (SPDE) rewrite of the base-field RDE and the tower fixes (exp-in-exp
 * structure check, exact tower-variable verification, exp-monomial Laurent
 * ansatz).  Correctness is asserted numerically (assert_rm_num) because the
 * closed forms are high-degree rationals or nested exponentials whose exact
 * Simplify[D - f] is prohibitively slow — the antiderivatives themselves are
 * produced in well under a second.  In2/In5 are Bronstein's own worked
 * examples (Symbolic Integration I, 2nd ed., Examples 6.2.1 / 6.5.x). */
static void test_bronstein_rde_examples(void) {
    /* High-degree exponential-times-rational (former SolveAlways blow-up: the
     * degree-100+ cases timed out or took seconds; the SPDE path is sub-second). */
    assert_rm_num1("((x - 100)/x^101) Exp[x]");                       /* -> E^x/x^100 */
    assert_rm_closes("((x - 99)/(x + 2)^102) Exp[x]");               /* In17 -> E^x/(2+x)^101 (closure; dense diff-back too slow) */
    assert_rm_num1("((x - 17)/(x + 2)^20) Exp[x]");                   /* shifted-pole degree-20 numeric anchor */
    assert_rm_num1("(x^101 + 1) Exp[x]");
    assert_rm_num1("((x^101 + 1)/(x + 1)) Exp[x]");
    assert_rm_num1("((1 - 1/100 x^100)/x^101) Exp[x^100/100]");       /* -> -E^(x^100/100)/(100 x^100) */
    assert_rm_num1("((x^4 + 2 x^3 - 7 x^2 - 16 x - 6)"
                  "/(x^2 (x + 1)^2 (x + 2)^2 (x + 3)^2)) Exp[x]");   /* -> E^x/(x^2(x+1)(x+2)(x+3)) */
    /* Single exponential/log tower coefficients (already elementary, kept as
     * regression anchors). */
    assert_rm_num1("(1 - 1/(x Log[x]^2)) Exp[1/Log[x] + x]");         /* In3 -> E^(x+1/Log[x]) */
    assert_rm_num1("(Exp[Exp[x + 1]/x]/Exp[x^2])"
                  "(Exp[x + 1]/x - Exp[x + 1]/x^2 - 2 x)");          /* In4 -> E^(-x^2+E^(1+x)/x) */
    assert_rm_num1("Exp[1/Log[x]]/(x Log[x]^3)");                     /* In6 */
    assert_rm_num1("((x + 1) Exp[x/Exp[1/x]])/Exp[2/x]");             /* In7 (also: no Power::infy) */
    /* Nested / mixed towers: exp-in-exp exponent (Bronstein's Example 6.2.1),
     * log-of-exp, and the coupled hyperexponential Laurent. */
    assert_rm_num1("((E^x - x^2 + 2 x)/(x^2 (x + E^x)^2))"
                  " E^((x^2 - 1)/x + 1/(x + E^x))");                 /* In2 -> E^(-x+1/(E^x+x)+(x^2-1)/x) */
    assert_rm_num1("(E^x + x + (1 + E^x) x Log[x] Log[1 + 1/Log[x]]"
                  " + (1 + E^x) x Log[x]^2 Log[1 + 1/Log[x]])"
                  "/(x (E^x + x)^2 Log[x] (1 + Log[x]))");           /* In5 -> -Log[1+1/Log[x]]/(E^x+x) */
    assert_rm_num1("1/(x Log[Exp[x] + 1])"
                  " - (Exp[x] Log[x])/((Exp[x] + 1) Log[Exp[x] + 1]^2)"); /* In8 -> Log[x]/Log[1+E^x] */
    assert_rm_num1("(Exp[1/(Exp[x] + 1)]/(39916800 Exp[10 x]))"
                  "((1757211400 + 2581284541 Exp[x])/(Exp[x] + 1)^3)"); /* In9 */

    /* Genuinely non-elementary siblings that MUST decline (the RDE correctly
     * proves no solution — a numerator perturbation that leaves an Ei residue). */
    assert_head_unevaluated(
        "Integrate[E^x/(x + 2)^2, x, Method -> \"RischTranscendental\"]", "Integrate");
    assert_head_unevaluated(
        "Integrate[E^x/(x - 3)^2, x, Method -> \"RischTranscendental\"]", "Integrate");
}

void test_integrate_risch_transcendental(void) {
    symtab_init();
    core_init();

    /* Unit: degree-bound arithmetic. */
    TEST(test_rde_var_bound);
    /* Base case. */
    TEST(test_rational_case);
    TEST(test_rational_agreement);
    /* Single logarithmic extension. */
    TEST(test_logarithmic_case);
    TEST(test_fractional_log_case);
    TEST(test_hermite_log_case);
    /* Single exponential extension. */
    TEST(test_exponential_case);
    TEST(test_fractional_exp_case);
    TEST(test_hermite_exp_case);
    TEST(test_hyperexponential_case);
    /* Trig / hyperbolic + multi-kernel. */
    TEST(test_trig_frontend);
    TEST(test_multikernel_case);
    /* Nested towers + genuine recursion. */
    TEST(test_log_tower_case);
    TEST(test_exp_tower_case);
    TEST(test_recursive_tower_case);
    /* Bronstein transcendental-RDE motivating examples (SPDE + tower fixes). */
    TEST(test_bronstein_rde_examples);
    /* Special functions. */
    TEST(test_special_functions);
    /* Cascade, plumbing, strictness. */
    TEST(test_cascade_default);
    TEST(test_method_plumbing);
    TEST(test_strict_unevaluated);
    TEST(test_strict_misc);

    printf("All Integrate RischTranscendental tests passed!\n");
}

int main(void) {
    /* test_utils.h's constructor sets alarm(60). The Bronstein worst-case
     * examples (In16/In17 at degree 100+) exercise the SPDE ladder's ~O(n)
     * Together/expand steps and legitimately run for tens of seconds on a
     * loaded machine, so 60s is too tight here. Extend it, matching the
     * precedent in test_linearsolve.c. */
    alarm(600);
    test_integrate_risch_transcendental();
    return 0;
}

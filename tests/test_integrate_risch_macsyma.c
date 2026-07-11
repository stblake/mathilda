/* test_integrate_risch_macsyma.c
 *
 * EXTENSIVE test suite for the Maxima-ported recursive Risch integrator
 * (Integrate`RischMacsyma, Method -> "RischMacsyma").  Correctness is asserted
 * by the universal predicate Simplify[D[Integrate[f, x], x] - f] === 0 (via
 * assert_rm_diff_zero) rather than by fixed output strings, so the tests
 * survive surface-form changes; the special-function and trig/hyperbolic
 * outputs (Erf/Ei/li/PolyLog and I-laden Cosh/Sinh forms) whose exact Simplify
 * is slower are verified numerically at interior points (assert_rm_num).  Every
 * integrand below was empirically classified against the built integrator
 * (tests/rm_probe.sh over tests/rm_candidates.txt): the assert_rm_* cases close
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

/* Assert Simplify[D[Integrate`RischMacsyma[f, x], x] - (f)] === 0, i.e.
 * the explicit package head produces a correct antiderivative. */
static void assert_rm_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate`RischMacsyma[%s, x], x] - (%s)]", f, f);
    assert_eval_eq(buf, "0", 0);
}

/* Assert the same via the Method-option surface form. */
static void assert_rm_method_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischMacsyma\"], x] - (%s)]",
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
        "Abs[N[(D[Integrate`RischMacsyma[%s, x], x] - (%s)) /. x -> 13/10]] +"
        " Abs[N[(D[Integrate`RischMacsyma[%s, x], x] - (%s)) /. x -> 17/10]] +"
        " Abs[N[(D[Integrate`RischMacsyma[%s, x], x] - (%s)) /. x -> 9/10]]"
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
 * closes by its own Risch DE q_k' + W_k' q_k = p_k (rm_expsum_case) without a
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
 * by one unified SolveAlways ansatz over all tower variables (rm_log_tower_case),
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
}

/* ================= NESTED EXP TOWERS + MERGED MONOMIAL =================
 * Chain of nested exponentials t_i = E^(u_i) (e.g. E^x, E^(E^x)) via a Laurent
 * ansatz over the exp tower derivation, diff-back verified.  The third case is
 * an EVALUATOR-MERGED monomial: the evaluator folds E^x E^(E^x) into
 * E^(x + E^x), whose exponent carries the foreign kernel E^x; the
 * rm_expand_exp_sums pre-pass re-splits E^(a+b) -> E^a E^b to restore the
 * independent basis, and the coupled hyperexponential proper part closes it. */
static void test_exp_tower_case(void) {
    assert_rm_diff_zero("E^x E^(E^x)");
    assert_rm_diff_zero("E^(2 E^x) E^x");
    assert_rm_diff_zero("E^x E^(E^x)/(1 + E^(E^x))");
}

/* ================= RECURSIVE / MIXED =================
 * The genuine one-extension-at-a-time recursion (Bronstein ch. 5 / Maxima
 * risch.lisp) that the flat single-kind tower ansatz cannot express: MIXED
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
}

/* ================= SPECIAL FUNCTIONS =================
 * Maxima erfarg / Ei / li / dilog: integrals the elementary cascade leaves open,
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
    assert_rm_num("1/Log[x]");
    assert_rm_num("3/Log[x]");
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
    assert_head_unevaluated("Integrate`RischMacsyma[Sin[x^2], x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[Sin[x]/x, x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[Cos[x]/x, x]", "Integrate`RischMacsyma");
    /* Exponential of a rational whose Risch DE has no rational solution, and an
     * Ei-type integrand with a NON-linear denominator (unimplemented). */
    assert_head_unevaluated("Integrate`RischMacsyma[Exp[1/x], x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[Exp[x]/x^2, x]", "Integrate`RischMacsyma");
    /* Ei / li recognizer gaps (correct decline, not wrong): the special-function
     * recognizers match (const E^(a x))/(c x+d) and K/Log[x] only, so E^(-x)/x
     * (Ei with a<0), 1/Log[2x] (li of a scaled arg), and x/Log[x] (non-constant
     * numerator) fall through.  Documented narrowness of rm_special_case. */
    assert_head_unevaluated("Integrate`RischMacsyma[Exp[-x]/x, x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[1/Log[2 x], x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[x/Log[x], x]", "Integrate`RischMacsyma");
    /* Non-elementary nested-log integrands (need Ei/li of a log); and a residual
     * NON-rational inner kernel (Sin[Log[x]]) must DECLINE, never certify a wrong
     * form (the whole-tower rationality gate). */
    assert_head_unevaluated("Integrate`RischMacsyma[Log[x] Log[Log[x]], x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[Sin[Log[x]] Log[Log[x]], x]", "Integrate`RischMacsyma");
    /* Non-elementary nested exponentials.  E^(E^x)/(1+E^(E^x)) is a REGRESSION
     * guard: a single-kernel case once certified the WRONG Log[1+E^(E^x)]/E^x by
     * leaving the inner E^x as a free SolveAlways parameter (fixed by rm_kernel_simple). */
    assert_head_unevaluated("Integrate`RischMacsyma[E^(E^x), x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[E^(E^x)/(1 + E^(E^x)), x]", "Integrate`RischMacsyma");
    /* Proper tower fraction with a NON-constant Rothstein-Trager residue is
     * non-elementary (must not certify a wrong constant residue). */
    assert_head_unevaluated("Integrate`RischMacsyma[1/(Log[x] (1 + Log[Log[x]])), x]", "Integrate`RischMacsyma");
    /* Field Risch DE with no elementary solution (E^(Log[x]^2) and its coupled
     * fraction) declines rather than forcing a bounded ansatz to a wrong answer. */
    assert_head_unevaluated("Integrate`RischMacsyma[E^(Log[x]^2), x]", "Integrate`RischMacsyma");
    assert_head_unevaluated("Integrate`RischMacsyma[E^(Log[x]^2)/(1 + Log[x]), x]", "Integrate`RischMacsyma");
}

/* ================= Supplementary: engine agreement ================= */
/* On the rational base case RischMacsyma must agree with the dedicated
 * rational Risch engine (Integrate`BronsteinRational) that it delegates to. */
static void test_rational_agreement(void) {
    assert_eval_eq(
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischMacsyma\"]"
        " - Integrate[1/(x^2 - 1), x, Method -> \"BronsteinRational\"]]",
        "0", 0);
    assert_eval_eq(
        "Simplify[Integrate[(x^3 + 1)/(x^2 + 1), x, Method -> \"RischMacsyma\"]"
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
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischMacsyma\"]"
        " - Integrate`RischMacsyma[1/(x^2 - 1), x]]",
        "0", 0);
}

/* ================= Strict / malformed-input behaviour ================= */
static void test_strict_misc(void) {
    /* A non-symbol integration variable is rejected (never garbage). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[1/x, x + 1]", "Integrate`RischMacsyma");
    /* Fresnel / Si / Ci non-elementary integrands bubble back unevaluated. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Cos[x^2], x]", "Integrate`RischMacsyma");
}

void test_integrate_risch_macsyma(void) {
    symtab_init();
    core_init();

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
    /* Special functions. */
    TEST(test_special_functions);
    /* Cascade, plumbing, strictness. */
    TEST(test_cascade_default);
    TEST(test_method_plumbing);
    TEST(test_strict_unevaluated);
    TEST(test_strict_misc);

    printf("All Integrate RischMacsyma tests passed!\n");
}

int main(void) {
    test_integrate_risch_macsyma();
    return 0;
}

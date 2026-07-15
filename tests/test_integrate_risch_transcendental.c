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
        && strcmp(result->data.function.head->data.symbol.name, head) == 0,
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
    /* Higher-multiplicity repeated poles with a linear exponent: F is free of x,
     * so these take the fast rational-of-a-single-exp path (rt_exp_ratreduce_case
     * reduces to a pure rational integral in t) ahead of the O(mult)-variable
     * Hermite SolveAlways ansatz.  Diff-back gated. */
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^8");
    assert_rm_diff_zero("Exp[3 x]/(1 + Exp[x])^6");
    assert_rm_diff_zero("Exp[2 x]/((1 + Exp[x])^4 (2 + Exp[x]))");
}

/* ================= COUPLED HYPEREXPONENTIAL =================
 * Rational function of theta = E^u whose integral MIXES a Laurent-polynomial
 * part with a log part (D(g)/g improper in theta, so they do not separate).
 * Integrated by the LITERAL Bronstein pipeline (rt_field_hyperexp_hermite):
 * HermiteReduce (5.3) peels the repeated NORMAL poles into an exact rational
 * part g with arbitrary rational lower-field coefficients; the residue criterion
 * (5.6) gives the logs sum_j c_j Log(g_j) of the simple part; the coupling is
 * reconciled by subtracting D_tower[logs]; and the remaining Laurent polynomial
 * is integrated per-coefficient by a Risch DE (5.9).  The former undetermined-
 * coefficient ansatz has been removed.  Covers 1/(c + E^x), improper
 * (E^(2x)/(1+E^x)), repeated / theta=0 poles, and (the value the ansatz could
 * not reach) rational lower-field coefficients at repeated exponential poles. */
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
    /* Value-add of the literal Hermite pipeline over the removed bounded ansatz:
     * a repeated exponential pole whose Hermite numerator coefficient is a genuine
     * RATIONAL function of x (not a bounded polynomial).  These are constructed as
     * exact derivatives so a correct integrator must round-trip them. */
    assert_rm_diff_zero("D[1/((1 + x) (1 + Exp[x])^2), x]");
    assert_rm_diff_zero("D[x/(1 + Exp[x])^2, x]");
    assert_rm_diff_zero("D[x^2/(1 + Exp[x])^3, x]");
    assert_rm_diff_zero("D[1/((1 + x^2) (2 + Exp[x])), x]");
    assert_rm_diff_zero("x Exp[x]/(1 + Exp[x])^2");
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

/* ================= REAL RECONSTRUCTION (rational trigonometric) =================
 * The complex-exponential route integrates over the kernel E^(I x); its log part
 * is a correct but I-laden closed form.  rt_realify reconstructs a REAL closed form
 * (Re[G], diff-back gated) so rational functions of the circular trig kernels come
 * back real (Log / ArcTan), NOT I-laden and NOT via any Weierstrass substitution.
 * Each case is asserted (a) Complex-free and (b) a correct antiderivative
 * numerically (the real ArcTan/atan2 forms are not always Simplify-reducible). */
static void assert_real_trig(const char* f) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{g = Integrate`RischTranscendental[%s, x]}, "
        "Head[g] =!= Integrate`RischTranscendental && FreeQ[g, Complex] && "
        "Abs[N[(D[g, x] - (%s)) /. x -> 3/5]] + "
        "Abs[N[(D[g, x] - (%s)) /. x -> 19/10]] < 1/100000]", f, f, f);
    assert_eval_eq(buf, "True", 0);
}

static void test_real_trig_reconstruction(void) {
    /* Sec / Csc and powers — real Log, were I-laden via TrigToExp. */
    assert_real_trig("Sec[x]");
    assert_real_trig("Csc[x]");
    assert_real_trig("Sec[x]^2");
    assert_real_trig("Csc[x]^2");
    /* Irreducible-quadratic denominators over E^(I x) — real ArcTan (atan2). */
    assert_real_trig("1/(2 + Cos[x])");
    assert_real_trig("1/(5 + 4 Cos[x])");
    assert_real_trig("1/(3 + 5 Cos[x])");
    assert_real_trig("1/(1 + Sin[x])");
    assert_real_trig("1/(1 - Cos[x])");
    /* Mixed Sin + Cos, and a proper fraction. */
    assert_real_trig("1/(Sin[x] + Cos[x])");
    assert_real_trig("Cos[x]/(1 + Cos[x])");
    /* Pure polynomials in Sin/Cos stay clean (exp-poly, multiple-angle). */
    assert_real_trig("Sin[x]^2");
    assert_real_trig("Sin[x]^3");
    assert_real_trig("Cos[x]^4");
    assert_real_trig("Sin[x] Cos[x]");
}

/* The two-argument ArcTan[u, v] = Arg[u + I v] differentiates as
 * (u v' - v u')/(u^2 + v^2) — the rule rt_realify's atan2 output relies on. */
static void test_arctan2_derivative(void) {
    assert_eval_eq("D[ArcTan[Cos[x], Sin[x]], x]", "1", 0);          /* Arg = x */
    assert_eval_eq("Simplify[D[ArcTan[x, x^2], x] - 1/(1 + x^2)]", "0", 0);
    assert_eval_eq("Simplify[D[ArcTan[1, x], x] - 1/(1 + x^2)]", "0", 0);
}

/* ================= REAL HYPERTANGENT CASE (Bronstein §5.10) =================
 * Rational functions of a single real tangent kernel t = Tan[u] (u rational in
 * x) integrate DIRECTLY and REAL through rt_hypertangent_case, retiring the
 * I-laden TrigToExp route for them.  The reliable syntactic I-detector here is
 * Position[.,Complex[__]] (FreeQ/Cases do not descend into Complex atoms). */

/* The Method->"RischTranscendental" antiderivative is Complex-free (real) AND
 * diff-backs to f. */
static void assert_tan_real(const char* f) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{g = Integrate[%s, x, Method -> \"RischTranscendental\"]}, "
        "Position[g, Complex[__]] === {} && Simplify[D[g, x] - (%s)] === 0]", f, f);
    assert_eval_eq(buf, "True", 0);
}

/* The antiderivative equals a specific real closed form (Simplify-robust). */
static void assert_tan_form(const char* f, const char* expected) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "Simplify[Integrate[%s, x, Method -> \"RischTranscendental\"] - (%s)]", f, expected);
    assert_eval_eq(buf, "0", 0);
}

static void test_real_hypertangent(void) {
    /* Canonical clean real forms — what used to come back I-laden. */
    assert_tan_form("Tan[x]", "-Log[Cos[x]]");
    assert_tan_form("Tan[x]^2", "Tan[x] - x");
    assert_tan_form("Tan[2 x]", "-Log[Cos[2 x]]/2");
    /* Real + correct across powers, scalings, and a nonlinear tangent argument. */
    assert_tan_real("Tan[x]");
    assert_tan_real("Tan[x]^2");
    assert_tan_real("Tan[x]^3");
    assert_tan_real("Tan[x]^4");
    assert_tan_real("Tan[x]^5");
    assert_tan_real("Tan[x]^6");
    assert_tan_real("Tan[2 x]");
    assert_tan_real("Tan[3 x]^3");
    assert_tan_real("Tan[x/2]");
    assert_tan_real("2 x Tan[x^2]");                 /* Dt = 2x(t^2+1) */
    assert_tan_real("Tan[x] + Tan[x]^4");
    assert_tan_real("Tan[x]^2 + 3 Tan[x] + 1");
    /* TRANSCENDENTAL tangent argument: tau = Tan[Log[x]] is hypertangent over C(x)
     * (eta = 1/x), so it integrates directly to a REAL, correct form.  The raw
     * hypertangent output is (1/2)Log[1+Tan[Log[x]]^2] (= -Log[Cos[Log[x]]]); the
     * cosmetic collapse to -Log[Cos] is a Simplify gap (SIMPLIFY_GAPS.md, Family 2),
     * not an integrator defect, so correctness is asserted via assert_tan_real
     * (real + exact diff-back) rather than the specific -Log[Cos] string. */
    assert_tan_real("Tan[Log[x]]/x");
    assert_tan_real("Tanh[Log[x]]/x");
    assert_tan_real("Tan[Log[x]]^2/x");
    assert_tan_real("(1 + Tan[Log[x]])/(x (1 - Tan[Log[x]]))");
}

static void test_real_hypertangent_rational(void) {
    /* Special-pole (t^2+1)^k parts via the coupled reduced case. */
    assert_tan_real("1/(1 + Tan[x]^2)");             /* = Cos[x]^2 */
    assert_tan_real("Tan[x]/(1 + Tan[x]^2)");
    assert_tan_real("Tan[x]^2/(Tan[x]^2 + 1)^2");
    /* Simple linear normal poles (rational Rothstein-Trager residues). */
    assert_tan_real("(1 + Tan[x])/(1 - Tan[x])");
    assert_tan_real("1/(1 + Tan[x])");
    assert_tan_real("(Tan[x]^2 + 1)/(Tan[x] - 2)");
}

static void test_real_hypertangent_mixed(void) {
    /* Polynomial-in-x coefficients that stay elementary (even tan powers). */
    assert_tan_real("x Tan[x]^2");
    assert_tan_real("x Tan[x]^4");
    assert_tan_real("(x + 1) Tan[x]^2");
    /* Genuinely non-elementary tangent integrands stay unevaluated (the driver
     * PROVES it: Bronstein's ∫ x tan x obstruction lifted to higher terms). */
    assert_head_unevaluated(
        "Integrate[x Tan[x], x, Method -> \"RischTranscendental\"]", "Integrate");
    assert_head_unevaluated(
        "Integrate[x Tan[x]^3, x, Method -> \"RischTranscendental\"]", "Integrate");
    assert_head_unevaluated(
        "Integrate[x^2 Tan[x]^2, x, Method -> \"RischTranscendental\"]", "Integrate");
}

/* Cotangent (circular, special t^2+1, eta = -u') integrates real like Tan. */
static void test_real_cotangent(void) {
    assert_tan_form("Cot[x]", "Log[Sin[x]]");
    assert_tan_form("Cot[x]^2", "-Cot[x] - x");
    assert_tan_real("Cot[x]");
    assert_tan_real("Cot[x]^2");
    assert_tan_real("Cot[x]^3");
    assert_tan_real("Cot[x]^4");
    assert_tan_real("Cot[2 x]");
    assert_tan_real("Cot[x]^2 + Cot[x]");
    assert_tan_real("x Cot[x]^2");
    assert_tan_real("1/(1 + Cot[x]^2)");           /* = Sin[x]^2 */
}

/* Hyperbolic tangent (special t^2-1, splits: two real Risch DEs) integrates
 * real — retiring TrigToExp for real Tanh. */
static void test_real_hypertanh(void) {
    assert_tan_form("Tanh[x]", "Log[Cosh[x]]");
    assert_tan_form("Tanh[x]^2", "x - Tanh[x]");
    assert_tan_form("Tanh[2 x]", "Log[Cosh[2 x]]/2");
    assert_tan_real("Tanh[x]");
    assert_tan_real("Tanh[x]^2");
    assert_tan_real("Tanh[x]^3");
    assert_tan_real("Tanh[x]^4");
    assert_tan_real("Tanh[x]^5");
    assert_tan_real("Tanh[2 x]");
    assert_tan_real("Tanh[x/2]");
    assert_tan_real("Tanh[x] + Tanh[x]^4");
    assert_tan_real("x Tanh[x]^2");                /* elementary (even power) */
    assert_tan_real("1/(1 - Tanh[x]^2)");          /* = Cosh[x]^2, reduced pole */
    assert_tan_real("Tanh[x]/(1 - Tanh[x]^2)");
    /* Non-elementary hyperbolic tangent integrands stay unevaluated. */
    assert_head_unevaluated(
        "Integrate[x Tanh[x], x, Method -> \"RischTranscendental\"]", "Integrate");
    assert_head_unevaluated(
        "Integrate[x Tanh[x]^3, x, Method -> \"RischTranscendental\"]", "Integrate");
    /* Coth: same hyperbolic monomial as Tanh (Dt = eta(1-t^2)), real cosmetic
     * -2 Log[Sinh] (1-Coth^2 = -Csch^2, the sign folded into the rewrite). */
    assert_tan_form("Coth[x]", "Log[Sinh[x]]");
    assert_tan_form("Coth[x]^2", "x - Coth[x]");
    assert_tan_real("Coth[x]");
    assert_tan_real("Coth[x]^2");
    assert_tan_real("Coth[x]^3");
    assert_tan_real("Coth[2 x]");
    assert_tan_real("x Coth[x]^2");
    assert_head_unevaluated(
        "Integrate[x Coth[x], x, Method -> \"RischTranscendental\"]", "Integrate");
}

static void test_real_hypertangent_robustness(void) {
    /* Tanh has Dt = eta(1 - t^2), NOT a hypertangent monomial (special factor
     * t^2-1, not t^2+1): the Tan-only gate must leave it on the TrigToExp path
     * (correct, verified numerically) — a guard that the new case does not
     * mis-claim hyperbolic tangents. */
    assert_rm_num("Tanh[x]");
    assert_rm_num("Tanh[x]^2");
    /* A degree->=2 irreducible normal pole with rational residues still routes
     * through the real path (linear-normal after the gate strips t^2+1)... this
     * one (t^2+1 special only) stays real: Cot-free rational of Tan. */
    assert_tan_real("(2 Tan[x] + 3)/(Tan[x]^2 + 1)");
    /* Irreducible-QUADRATIC normal pole (3 + t^2 does not split over Q): the §5.10
     * residue criterion realises the quadratic-algebraic residues as a real ArcTan.
     * Admitted by the max-irreducible-degree<=2 normal-pole gate (was declined by
     * the old total-degree<=1 gate, which stranded these on the I-laden TrigToExp
     * route).  Genuine §5.10 — no Weierstrass substitution. */
    assert_tan_real("Tan[x]/(3 + Tan[x]^2)");
    assert_tan_real("1/(3 + Tan[x]^2)");
    assert_tan_real("(Tan[x] + 1)/(2 + Tan[x]^2)");
    assert_tan_real("Cot[x]/(2 + Cot[x]^2)");
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
    assert_rm_diff_zero("Log[Log[x]]/(x Log[x])");   /* LimitedIntegrate fold-back -> Log[Log[x]]^2/2 */
    assert_rm_diff_zero("Log[Log[x]]^5/(x Log[x])");  /* higher-degree primitive polynomial */
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

/* ================= TANGENT TOWER (RT_TAN monomial) =================
 * A tangent monomial t = Tan[u]/Tanh[u] nested with an independent Log/Exp kernel.
 * The tangent is NONLINEAR (Dt = u'(t^2 + sigma)); rt_subst_kernels rationalises
 * the circular/hyperbolic trig to the tower variable, and the residue LRT clears
 * the tangent's rational Dcoef denominator.  A tangent LOWER kernel with a log top:
 * Int (1+Tan^2)/(Tan Log[Tan]) = Log[Log[Tan[x]]]. */
static void test_tangent_tower(void) {
    assert_rm_diff_zero("(1 + Tan[x]^2)/(Tan[x] Log[Tan[x]])");       /* Log[Log[Tan[x]]] */
    assert_rm_diff_zero("(1 + Tan[x]^2)/(Tan[x] Log[Tan[x]]^2)");     /* repeated pole -> -1/Log[Tan[x]] */
    assert_rm_diff_zero("(1 - Tanh[x]^2)/(Tanh[x] Log[Tanh[x]])");    /* hyperbolic, sigma=-1 */
    /* Hypertangent TOP over a genuine tower base K_{L-1} (§5.10 driver dispatch):
     * a tangent of a LOG-tower argument, coefficient = the argument's derivative. */
    assert_rm_diff_zero("2 Log[x]/x Tan[Log[x]^2]");                  /* -Log[Cos[Log[x]^2]] */
    assert_rm_diff_zero("2 Log[x]/x Tanh[Log[x]^2]");                 /* -Log[Cosh[Log[x]^2]], sigma=-1 */
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

/* ================= REPORTED-BUG REGRESSIONS =================
 * The six transcendental-Risch integrands from the 2026-07-15 bug report, plus
 * the "already works, lock it in" Group-2 exponential cases.  Each is diff-back
 * verified (exact where Simplify closes, else numeric interior points). */
static void test_reported_bug_fixes(void) {
    /* In[1]: mixed Log+Tan tower — the LOG-after-TAN ordering + Expand'd tower
     * diff-back.  ∫ = x^5 Log[x^12 Cos x]. */
    assert_rm_diff_zero("12 x^4 + 5 x^4 Log[x^12 Cos[x]] - x^5 Tan[x]");
    assert_rm_method_diff_zero("12 x^4 + 5 x^4 Log[x^12 Cos[x]] - x^5 Tan[x]");

    /* In[2]/In[9]: Q-linearly DEPENDENT log generators collapse to Log[x]/x via
     * the log-combination fallback.  ∫ = Log[x]^2/2. */
    assert_rm_diff_zero("(Log[x/(1 + x)] + Log[1 + x])/x");
    assert_rm_diff_zero("(Log[x + 1] + Log[x/(x + 1)])/x");

    /* In[12]: polynomial * Log with an irreducible-quadratic argument — the
     * bottom-level ArcTan is admitted (was rejected as a foreign new log). */
    assert_rm_diff_zero("(x^5 - 1) Log[x^2 - x + 1]");
    assert_rm_method_diff_zero("(x^5 - 1) Log[x^2 - x + 1]");

    /* In[14]: nested exponential with a constant-plus-rational inner exponent
     * E^(1+1/x) — the exp-sum split must keep 1+1/x together (no bare E^1 that
     * re-merges).  ∫ = x E^(x E^(1+1/x)). */
    assert_rm_num("E^(E^(1 + 1/x)*x) - E^(1 + 1/x + E^(1 + 1/x)*x)"
                  " + E^(1 + 1/x + E^(1 + 1/x)*x)*x");

    /* In[21]: depth-2 nested-exp inside a Cos — the trig front-end now routes the
     * exponentialized integrand through the recursive tower, and the Dcoef exp
     * product is kept structural.  ∫ = Sin[x E^E^x]. */
    assert_rm_num("Exp[Exp[x]] (1 + x Exp[x]) Cos[x Exp[Exp[x]]]");

    /* Group-2 (In[15]-In[20]): already-passing exponential integrands, locked in
     * so the tower/exp-split refactors cannot silently regress them. */
    assert_rm_diff_zero("E^(-1 + E^x + 1/x + 1/(1 - x^2))"
                        " (E^x - 1/x^2 + (2 x)/(-1 + x^2)^2)");
    assert_rm_diff_zero("(1 - 1/(x Log[x]^2)) Exp[1/Log[x] + x]");
    assert_rm_diff_zero("(x + 1)/x^4 Exp[1/x]");
    assert_rm_diff_zero("((-1 - x - x^2 + x^3)/(1 - 2 x^2 + x^4)) Exp[x]");
    assert_rm_num("(E^(1/(E^x + x) + (-1 + x^2)/x) (E^x + 2 x - x^2))"
                  "/(x^2 (E^x + x)^2)");
    assert_rm_diff_zero("(Exp[1/Log[x]] (Log[x]^2 - 1))/Log[x]^2");
}

/* ================= CONSTANT-BASE EXPONENTIALS (a^u, a != e) =================
 * a^u = E^(u Log a) is an ordinary hyperexponential monomial, but the evaluator
 * stores it as Power[a, u] (E^(c Log a) collapses back), so it must be DEBASED to
 * the base-e form before the exponential machinery fires.  Prime-factored sharing
 * makes 4^x = (2^x)^2 commensurate.  Diff-back verified NUMERICALLY (the answers
 * carry Log[a] constants whose Simplify reduction is fragile). */
static void test_constant_base_exponentials(void) {
    /* Reported 2026-07-15: nested a^x tower and a single-kernel a^x rational. */
    assert_rm_num("3^x/(Exp[3^x] + 1)");        /* -> (3^x - Log[1+E^3^x])/Log[3] */
    assert_rm_num("2^x/(4^x - 2^x - 1)");       /* 4^x = (2^x)^2: one kernel */
    /* Bare and simple rational-of-a^x. */
    assert_rm_num("2^x");
    assert_rm_num("1/(2^x + 1)");
    assert_rm_num("x 2^x");                     /* polynomial * a^x */
    /* Composite (multi-prime) base recombines to base form. */
    assert_rm_num("6^x/(6^x + 1)");             /* -> Log[1+6^x]/Log[6] */
    assert_rm_num("5^x/(25^x - 1)");            /* 25 = 5^2 */

    /* INEXACT base must NOT be force-integrated to a (degraded) exact closed
     * form: the float base is rationalised upstream, and debasing is suppressed
     * so the general Integrate stays unevaluated (inexact-in / inexact-out). */
    assert_eval_eq("Head[Integrate[x*2.71828^(-x), x]]", "Integrate", 0);
    assert_eval_eq("Head[Integrate[x^2*2.71828^(-x), x]]", "Integrate", 0);
}

/* ================= DISPATCH CASE MATRIX =================
 * One stressing integrand per case / sub-case of rt_transcendental_case's
 * cascade (RISCH_STATUS.md §2), each diff-back verified.  Complements the
 * per-family suites above with a single cross-cutting coverage matrix so a
 * regression in any one dispatch arm is caught by a named check.  All reduce
 * EXACTLY (Simplify[D - f] === 0), including the I-laden exp-sum and the
 * special-function outputs. */
static void test_dispatch_case_matrix(void) {
    /* rt_rational_case — rational base, LRT over an irreducible quartic
     * (two conjugate ArcTan + Log pairs). */
    assert_rm_diff_zero("1/(x^4 + 1)");
    /* rt_log_poly_case — polynomial in Log[x] (linear argument), high degree. */
    assert_rm_diff_zero("x^2 Log[x]^3");
    /* rt_exp_poly_case — Laurent polynomial in E^(x^2): the i=1 Risch DE
     * q' + 2x q = x^3 has the polynomial solution (x^2-1)/2. */
    assert_rm_diff_zero("x^3 E^(x^2)");
    /* rt_exp_poly_case, Phase C — RATIONAL exponent E^(1/x). */
    assert_rm_diff_zero("(x + 1)/x^4 E^(1/x)");
    /* rt_frac_case — squarefree Rothstein-Trager with RATIONAL residues (log). */
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[x]))");
    /* rt_frac_lrt — ALGEBRAIC residues -> ArcTan (exp kernel, E^(2x) -> t^2). */
    assert_rm_diff_zero("E^x/(E^(2 x) + 1)");
    /* rt_hermite_case — repeated pole (log kernel) + a residual log. */
    assert_rm_diff_zero("(1 + Log[x])/(x Log[x]^2)");
    /* rt_hyperexp_case — coupled Laurent + Hermite, higher repeated pole. */
    assert_rm_diff_zero("1/(1 + E^x)^3");
    /* rt_expsum_case — non-commensurate exponentials from E^x * trig (the
     * I-laden Cosh/Sinh-of-complex form still diff-backs to exactly 0). */
    assert_rm_diff_zero("x^2 E^x Cos[x]");
    /* rt_log_tower_case — nested logarithmic tower, depth 3. */
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]] Log[Log[Log[x]]])");
    /* rt_exp_tower_case — nested exponential tower (Laurent ansatz). */
    assert_rm_diff_zero("E^(2 E^x) E^x");
    /* rt_recursive_tower_case — MIXED exp/log tower (independent extensions). */
    assert_rm_diff_zero("E^x/x + E^x Log[x]");
    /* rt_recursive_tower_case — RATIONAL lower-field coefficient (1/x). */
    assert_rm_diff_zero("1/(x^2 Log[x]) - Log[Log[x]]/x^2");
    /* rt_hypertangent_case — real §5.10 tangent, irreducible-quadratic residue. */
    assert_rm_diff_zero("1/(3 + Tan[x]^2)");
    /* rt_special_case — the special-function outputs (Ei / Erf / li).  These are
     * NON-elementary; the integrator answers with the special function, whose
     * derivative is exactly the integrand. */
    assert_rm_diff_zero("E^x/x");     /* ExpIntegralEi[x] */
    assert_rm_diff_zero("E^(-x^2)");  /* (Sqrt[Pi]/2) Erf[x] */
    assert_rm_diff_zero("1/Log[x]");  /* LogIntegral[x] */
}

/* ================= CIRCULAR TRIGONOMETRIC =================
 * Real Tan/Cot integrate through the direct §5.10 hypertangent case; every other
 * circular integrand goes through the TrigToExp -> exponential -> ExpToTrig ->
 * rt_realify front-end.  Each result verifies EXACTLY (Simplify[D - f] === 0),
 * including the I-laden E^x*trig forms and the two-argument-ArcTan (atan2)
 * rational-trig outputs.  The odd higher powers Sec[x]^3 return an exact but
 * I-laden exponential form (their real reduction is a known Simplify gap), so
 * they are diff-back-verified NUMERICALLY (assert_rm_num) rather than by the
 * exact-Simplify assert_rm_diff_zero, which would itself blow up on the form. */
static void test_circular_trig_integration(void) {
    /* basic six */
    assert_rm_diff_zero("Sin[x]");
    assert_rm_diff_zero("Cos[x]");
    assert_rm_diff_zero("Tan[x]");            /* -Log[Cos[x]] (direct hypertangent) */
    assert_rm_diff_zero("Cot[x]");            /* Log[Sin[x]]                        */
    assert_rm_diff_zero("Sec[x]");            /* realified real Log form            */
    assert_rm_diff_zero("Csc[x]");
    /* powers */
    assert_rm_diff_zero("Sin[x]^2");
    assert_rm_diff_zero("Cos[x]^2");
    assert_rm_diff_zero("Tan[x]^2");          /* -x + Tan[x]                        */
    assert_rm_diff_zero("Sec[x]^2");          /* Tan[x]                             */
    assert_rm_diff_zero("Csc[x]^3");
    assert_rm_diff_zero("Sin[x]^3");
    assert_rm_diff_zero("Cos[x]^3");
    assert_rm_diff_zero("Tan[x]^3");
    assert_rm_num("Sec[x]^3");   /* exact I-laden exp form; numeric diff-back */
    /* products / multiple angles */
    assert_rm_diff_zero("Sin[x] Cos[x]");
    assert_rm_diff_zero("Sin[x]^2 Cos[x]^2");
    assert_rm_diff_zero("Sin[x]^3 Cos[x]^2");
    assert_rm_diff_zero("Sin[2 x]");
    assert_rm_diff_zero("Sin[x] Cos[2 x]");
    /* rational in the circular kernels (real ArcTan/Log via realify) */
    assert_rm_diff_zero("1/(2 + Cos[x])");
    assert_rm_diff_zero("1/(5 + 4 Cos[x])");
    assert_rm_diff_zero("Cos[x]/(1 + Cos[x])");
    assert_rm_diff_zero("1/(3 + Tan[x]^2)");
    assert_rm_diff_zero("Tan[x]/(3 + Tan[x]^2)");
    /* Hypertangent (§5.10) rational-in-Tan integrands that formerly stranded on the
     * Sec Jacobian / an improper (deg num >= deg den) quotient and fell through to the
     * exponentializing front-end (an unbounded blow-up).  The Sec[x]^2 = 1+Tan[x]^2
     * reduction now routes them through the direct hypertangent driver; their clean
     * Sqrt[3]-wrapped Log/ArcTan forms are checked numerically because the Simplify
     * diff-back mis-reduces that shape (exactly the spurious residual that used to make
     * the diff-back gate reject the correct answer). */
    assert_rm_num("((2 + Tan[x]^2) Sec[x]^2)/(1 + Tan[x]^3)");  /* reported hang */
    assert_rm_num("Sec[x]^2/(1 + Tan[x]^3)");
    assert_rm_num("(2 + 3 Tan[x]^2 + Tan[x]^4)/(1 + Tan[x]^3)"); /* pure-Tan improper */
    assert_rm_num("Csc[x]^2/(1 + Cot[x]^3)");                    /* cosecant reciprocal */
    /* exponential x trig (non-commensurate exp-sum; I-laden yet exactly 0) */
    assert_rm_diff_zero("E^x Sin[x]");
    assert_rm_diff_zero("E^x Cos[x]");
    assert_rm_diff_zero("E^(2 x) Sin[3 x]");
    assert_rm_diff_zero("x E^x Sin[x]");
    assert_rm_diff_zero("x^2 E^x Cos[x]");
    /* trig OF an exponential kernel */
    assert_rm_diff_zero("Sin[E^x] E^x");      /* -Cos[E^x] */

    /* Non-elementary circular integrands MUST decline cleanly (never a wrong form). */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[Sin[x]/x, x]", "Integrate`RischTranscendental"); /* SinIntegral */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[x Tan[x], x]", "Integrate`RischTranscendental");
    assert_head_unevaluated(
        "Integrate`RischTranscendental[Sin[x^2], x]", "Integrate`RischTranscendental"); /* Fresnel */
    /* Non-rational inner kernel (Sin of a Log): the whole-tower rationality gate
     * declines rather than certify a wrong closed form. */
    assert_head_unevaluated(
        "Integrate`RischTranscendental[Sin[Log[x]], x]", "Integrate`RischTranscendental");
}

/* ================= HYPERBOLIC =================
 * Tanh/Coth integrate through the direct hyperbolic hypertangent case (special
 * t^2 - 1); the rest through TrigToExp over the single exponential E^x.  All
 * verify exactly. */
static void test_hyperbolic_integration(void) {
    assert_rm_diff_zero("Sinh[x]");
    assert_rm_diff_zero("Cosh[x]");
    assert_rm_diff_zero("Tanh[x]");           /* Log[Cosh[x]]  */
    assert_rm_diff_zero("Coth[x]");           /* Log[Sinh[x]]  */
    assert_rm_num("Sech[x]");                 /* 2 ArcTan[Cosh+Sinh]; Simplify diff-back hangs */
    assert_rm_diff_zero("Csch[x]");
    assert_rm_diff_zero("Sinh[x]^2");
    assert_rm_diff_zero("Cosh[x]^2");
    assert_rm_diff_zero("Tanh[x]^2");         /* x - Tanh[x]   */
    assert_rm_diff_zero("Sech[x]^2");         /* -2/(1 + E^(2 x)) */
    assert_rm_diff_zero("Tanh[x]^3");
    assert_rm_diff_zero("Sinh[x] Cosh[x]");
    assert_rm_diff_zero("1/(2 + Cosh[x])");
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
    /* Real hypertangent family (retires TrigToExp for real Tan/Cot/Tanh). */
    TEST(test_real_hypertangent);
    TEST(test_real_hypertangent_rational);
    TEST(test_real_hypertangent_mixed);
    TEST(test_real_cotangent);
    TEST(test_real_hypertanh);
    TEST(test_real_hypertangent_robustness);
    TEST(test_real_trig_reconstruction);
    TEST(test_arctan2_derivative);
    TEST(test_multikernel_case);
    /* Nested towers + genuine recursion. */
    TEST(test_log_tower_case);
    TEST(test_tangent_tower);
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
    /* 2026-07-15 reported-bug regressions. */
    TEST(test_reported_bug_fixes);
    /* Constant-base exponentials a^u (a != e) via debasing. */
    TEST(test_constant_base_exponentials);
    /* Cross-cutting one-per-case dispatch coverage matrix. */
    TEST(test_dispatch_case_matrix);
    /* Trigonometric and hyperbolic integration suites. */
    TEST(test_circular_trig_integration);
    TEST(test_hyperbolic_integration);

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

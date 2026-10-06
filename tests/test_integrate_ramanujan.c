/* test_integrate_ramanujan.c
 *
 * Tests for half-line definite integration by the Mellin-transform / Ramanujan
 * Master Theorem method (Integrate`RamanujanMasterTheorem, Method ->
 * "RamanujanMasterTheorem" / "Mellin").  Correctness is asserted by
 * Simplify[<integral> - <expected>, <assumptions>] === 0 -- the check is
 * symbolic and lives in the test, not the method (the method performs no
 * NIntegrate crosscheck, per project rule).
 *
 * The pinned Method is used throughout so the tests exercise this method in
 * isolation, independent of the Automatic cascade's residue / Newton-Leibniz
 * stages (which run first).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert Simplify[<integral> - <expected>, <assum>] === 0. */
static void assert_closes(const char* integral, const char* expected,
                          const char* assum) {
    char buf[1024];
    if (assum && assum[0])
        snprintf(buf, sizeof(buf), "Simplify[(%s) - (%s), %s]", integral, expected, assum);
    else
        snprintf(buf, sizeof(buf), "Simplify[(%s) - (%s)]", integral, expected);
    assert_eval_eq(buf, "0", 0);
}

/* Strip a ConditionalExpression[value, strip] wrapper and assert
 * Simplify[value - expected, assum] === 0.  Used for the symbolic-parameter
 * cases, whose strip is carried as a ConditionalExpression. */
static void assert_cond_closes(const char* integral, const char* expected,
                               const char* assum) {
    char buf[2048];
    const char* a = (assum && assum[0]) ? assum : "True";
    snprintf(buf, sizeof(buf),
        "Simplify[((%s) /. ConditionalExpression[cv_, cc_] :> cv) - (%s), %s]",
        integral, expected, a);
    assert_eval_eq(buf, "0", 0);
}

/* Strip the ConditionalExpression and assert the value equals `expected`
 * numerically at the substitution `subst` (a strip-interior parameter point).
 * Used where the closed forms differ only by trig/Gamma-reflection identities
 * that Simplify does not close symbolically. */
static void assert_cond_num(const char* integral, const char* expected,
                            const char* subst) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Chop[N[(((%s) /. ConditionalExpression[cv_, cc_] :> cv) /. %s) "
        "- ((%s) /. %s), 20]]",
        integral, subst, expected, subst);
    assert_eval_eq(buf, "0", 0);
}

/* Assert that the evaluated `input` is an unevaluated call with head `head`. */
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

/* Gaussian moments Integrate[x^n Exp[-x^2], {x,0,Inf}] = (1/2) Gamma[(n+1)/2]. */
static void test_gaussian_moments(void) {
    assert_closes("Integrate[Exp[-x^2], {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "Sqrt[Pi]/2", "");
    assert_closes("Integrate[x^2 Exp[-x^2], {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "Sqrt[Pi]/4", "");
    /* Scaling lam: Integrate[Exp[-b x^2]] = (1/2) Sqrt[Pi/b]. */
    assert_closes("Integrate[Exp[-3 x^2], {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "Sqrt[Pi]/(2 Sqrt[3])", "");
}

/* Exponential / Gamma: Integrate[x^(s-1) Exp[-x], {x,0,Inf}] = Gamma[s]. */
static void test_exponential_gamma(void) {
    assert_closes("Integrate[x^(s-1) Exp[-x], {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\", Assumptions -> s > 0]",
                  "Gamma[s]", "s > 0");
    /* Simple decay rate: Integrate[Exp[-a x], {x,0,Inf}] = 1/a. */
    assert_closes("Integrate[Exp[-a x], {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\", Assumptions -> a > 0]",
                  "1/a", "a > 0");
}

/* Algebraic binomial (p + q x^m)^(-a): Beta-integral family. */
static void test_algebraic_binomial(void) {
    assert_closes("Integrate[1/(1+x^2), {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "Pi/2", "");
    assert_closes("Integrate[1/(1+x)^2, {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "1", "");
    /* Branch power: Integrate[Sqrt[x]/(1+x)^2] = Pi/2. */
    assert_closes("Integrate[Sqrt[x]/(1+x)^2, {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]",
                  "Pi/2", "");
}

/* Trig transforms, including the removable s=0 (Sin[x]/x) Limit fallback. */
static void test_trig(void) {
    /* Dirichlet: Integrate[Sin[x]/x] = Pi/2 (s = 0, resolved by Limit). */
    assert_closes("Integrate[Sin[x]/x, {x, 0, Infinity}, Method -> \"RamanujanMasterTheorem\"]",
                  "Pi/2", "");
    /* Fresnel-type cosine: Integrate[Cos[x]/Sqrt[x]] = Sqrt[Pi/2]. */
    assert_closes("Integrate[Cos[x]/Sqrt[x], {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]",
                  "Sqrt[Pi/2]", "");
    /* Scaling with an assumed-positive rate. */
    assert_closes("Integrate[Sin[a x]/x, {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\", Assumptions -> a > 0]",
                  "Pi/2", "a > 0");
}

/* Bessel-Mellin: Integrate[BesselJ[0, x], {x,0,Inf}] = 1. */
static void test_bessel(void) {
    assert_closes("Integrate[BesselJ[0, x], {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]",
                  "1", "");
}

/* Linearity: a sum of individually-convergent terms is summed term by term. */
static void test_linearity(void) {
    assert_closes("Integrate[(1 + x^2) Exp[-x^2], {x, 0, Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]",
                  "3 Sqrt[Pi]/4", "");
}

/* Reachable three ways: Method string, the "Mellin" alias, the direct builtin. */
static void test_routing(void) {
    assert_closes("Integrate[Exp[-x^2], {x, 0, Infinity}, Method -> \"Mellin\"]",
                  "Sqrt[Pi]/2", "");
    assert_closes("Integrate`RamanujanMasterTheorem[x^2 Exp[-x^2], {x, 0, Infinity}]",
                  "Sqrt[Pi]/4", "");
    /* Also picked up by the Automatic cascade: Gaussian moments are closed by
     * this method alone (residue declines, FTC has no Erf antiderivative). */
    assert_closes("Integrate[x^2 Exp[-x^2], {x, 0, Infinity}]", "Sqrt[Pi]/4", "");
}

/* Tier A -- the strip is carried as a ConditionalExpression when the assumptions
 * do not prove it, and collapses to the bare value when they do. */
static void test_conditional_strip(void) {
    /* No assumptions: honest ConditionalExpression. */
    assert_eval_eq("Head[Integrate[x^(s-1) Exp[-x], {x,0,Infinity}, Method -> \"Mellin\"]]",
                   "ConditionalExpression", 0);
    /* Assumptions prove the strip -> bare Gamma[s] (existing exponential test
     * also covers this).  Here confirm the collapse explicitly. */
    assert_eval_eq("Integrate[x^(s-1) Exp[-x], {x,0,Infinity}, "
                   "Method -> \"Mellin\", Assumptions -> s > 0]", "Gamma[s]", 0);
}

/* Tier B -- monomial internal substitution g(x^k) (Sin[Sqrt x], BesselJ[.,2Sqrt x]). */
static void test_monomial_substitution(void) {
    assert_cond_num(
        "Integrate[x^(s-1) Sin[Sqrt[x]]/Sqrt[x], {x,0,Infinity}, Method -> \"Mellin\"]",
        "Pi/(Sin[Pi s] Gamma[2 - 2 s])", "{s -> 1/2}");
    /* Ramanujan's canonical example. */
    assert_cond_closes(
        "Integrate[x^(s-1) BesselJ[nu, 2 Sqrt[x]]/x^(nu/2), {x,0,Infinity}, Method -> \"Mellin\"]",
        "Gamma[s]/Gamma[1 + nu - s]", "");
}

/* Tier C -- Log, ArcTan, and the pFq master kernel (1F1 / 2F1). */
static void test_log_arctan_pfq(void) {
    assert_cond_num(
        "Integrate[x^(s-1) Log[1+x]/x, {x,0,Infinity}, Method -> \"Mellin\"]",
        "Pi/((1 - s) Sin[Pi s])", "{s -> 1/2}");
    assert_cond_num(
        "Integrate[x^(s-1) ArcTan[Sqrt[x]]/Sqrt[x], {x,0,Infinity}, Method -> \"Mellin\"]",
        "Pi/((1 - 2 s) Sin[Pi s])", "{s -> 1/4}");
    assert_cond_closes(
        "Integrate[x^(s-1) Hypergeometric1F1[a,c,-x], {x,0,Infinity}, Method -> \"Mellin\"]",
        "(Gamma[c] Gamma[s] Gamma[a-s])/(Gamma[a] Gamma[c-s])", "");
    assert_cond_closes(
        "Integrate[x^(s-1) Hypergeometric2F1[a,b,c,-x], {x,0,Infinity}, Method -> \"Mellin\"]",
        "(Gamma[c] Gamma[s] Gamma[a-s] Gamma[b-s])/(Gamma[a] Gamma[b] Gamma[c-s])", "");
}

/* Tier D/G -- hypergeometric reductions: lower incomplete gamma, Erf, and the
 * product BesselJ^2 (a Mellin convolution closed via the J^2 -> 1F2 identity). */
static void test_reductions(void) {
    /* Cancellation kernel Gamma[a]-Gamma[a,x]: must be handled before Expand. */
    assert_cond_closes(
        "Integrate[x^(s-1) (Gamma[a]-Gamma[a,x])/x^a, {x,0,Infinity}, Method -> \"Mellin\"]",
        "Gamma[s]/(a - s)", "");
    assert_cond_closes(
        "Integrate[x^(s-1) Erf[Sqrt[x]]/Sqrt[x], {x,0,Infinity}, Method -> \"Mellin\"]",
        "(Gamma[s] Gamma[1/2 - s])/(Sqrt[Pi] Gamma[3/2 - s])", "");
    /* BesselJ[nu,2Sqrt x]^2/x^nu: numeric spot-check (nu=2, s=1 -> 1/3). */
    assert_eval_eq(
        "Chop[N[((Integrate[x^(s-1) BesselJ[nu,2 Sqrt[x]]^2/x^nu, {x,0,Infinity}, "
        "Method -> \"Mellin\"] /. ConditionalExpression[cv_,cc_] :> cv) "
        "/. {nu -> 2, s -> 1}) - 1/3]]", "0", 0);
}

/* Tier E -- PolyLog / Dirichlet-series kernel, with the nu=1 = Log consistency. */
static void test_polylog(void) {
    assert_cond_num(
        "Integrate[x^(s-1) (-PolyLog[nu,-x]/x), {x,0,Infinity}, Method -> \"Mellin\"]",
        "Pi/((1 - s)^nu Sin[Pi s])", "{s -> 1/2, nu -> 2}");
    /* nu = 1: -PolyLog[1,-x] = Log[1+x], so this must equal the Log transform. */
    assert_cond_num(
        "Integrate[x^(s-1) (-PolyLog[1,-x]/x), {x,0,Infinity}, Method -> \"Mellin\"]",
        "Pi/((1 - s) Sin[Pi s])", "{s -> 1/2}");
}

/* Tier F -- parametric differentiation: Log-power / Log*(1+x)^(-w0) kernels. */
static void test_parametric_differentiation(void) {
    /* Log[1+x]/(x(1+x)):  n=1, w0=1. */
    assert_cond_num(
        "Integrate[x^(s-1) Log[1+x]/(x(1+x)), {x,0,Infinity}, Method -> \"Mellin\"]",
        "(Pi (EulerGamma + PolyGamma[0, 2 - s]))/Sin[Pi s]", "{s -> 1/2}");
    /* Log[1+x]^2/x^2:  n=2, w0=0. */
    assert_cond_num(
        "Integrate[x^(s-1) Log[1+x]^2/x^2, {x,0,Infinity}, Method -> \"Mellin\"]",
        "(2 Pi (EulerGamma + PolyGamma[0, 2 - s]))/((2 - s) Sin[Pi s])", "{s -> 1/2}");
}

/* Exponential-geometric kernel 1/(e^(c x) + gamma): Bose-Einstein (Gamma Zeta),
 * Fermi-Dirac (Gamma eta), and general fugacity (Gamma PolyLog). */
static void test_bose_fermi(void) {
    /* Bose-Einstein / Debye: Gamma(4) Zeta(4) = Pi^4/15. */
    assert_closes("Integrate[x^3/(Exp[x]-1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Pi^4/15", NULL);
    /* Bose s = 2. */
    assert_closes("Integrate[x/(Exp[x]-1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Pi^2/6", NULL);
    /* Fermi-Dirac: Gamma(4) eta(4) = 7 Pi^4/120. */
    assert_closes("Integrate[x^3/(Exp[x]+1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "7 Pi^4/120", NULL);
    /* Fermi s = 1: the -PolyLog(s,-1) form stays finite where (1-2^(1-s))Zeta(s)
     * would be 0*Infinity. */
    assert_closes("Integrate[1/(Exp[x]+1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Log[2]", NULL);
    /* Scaled rate c = 2: c^(-s) factor -> Pi^4/240. */
    assert_closes("Integrate[x^3/(Exp[2 x]-1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Pi^4/240", NULL);
    /* Symbolic order (Bose): Gamma(s) Zeta(s), strip s > 1 discharged. */
    assert_cond_closes("Integrate[x^(s-1)/(Exp[x]-1), {x,0,Infinity}, "
                       "Method -> \"RamanujanMasterTheorem\", Assumptions -> s > 1]",
                       "Gamma[s] Zeta[s]", "s > 1");
    /* Monomial substitution: 1/(e^Sqrt[x]-1) -> 2 Zeta(2) = Pi^2/3. */
    assert_closes("Integrate[1/(Exp[Sqrt[x]]-1), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Pi^2/3", NULL);
    /* Divergent Bose side gamma' = -2 < -1 (interior pole at x = Log[2]) ->
     * unevaluated (never a wrong value).  The relaxed gate admits gamma' > 1
     * (high-fugacity Fermi, see test_high_fugacity_fermi) but not gamma' < -1. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[1/(Exp[x]-2), {x,0,Infinity}]",
        "Integrate`RamanujanMasterTheorem");
    /* Divergent Bose at s = 1 (strip s > 1 excludes it) -> unevaluated. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[1/(Exp[x]-1), {x,0,Infinity}]",
        "Integrate`RamanujanMasterTheorem");
}

/* Symbolic fugacity: the general Bose integral 1/(z^-1 e^x - 1) with a SYMBOLIC
 * z admitted by the interval-bound gate when the Assumptions confine gamma' to
 * (-1, 1].  Int_0^Inf x^(s-1)/(z^-1 e^x - 1) dx = Gamma(s) PolyLog(s, z). */
static void test_symbolic_fugacity(void) {
    /* The reported case: fugacity form 1/(z^-1 e^x - 1) -> Gamma(s) PolyLog(s, z),
     * under Re[s] > 1 && 0 < z < 1 (gamma' = -z is proved into (-1, 0) by the
     * interval gate; the assumption engine cannot discharge -1 <= -z <= 1). */
    assert_cond_closes(
        "Integrate[x^(s-1)/(z^-1 Exp[x]-1), {x,0,Infinity}, "
        "Method -> \"RamanujanMasterTheorem\", Assumptions -> Re[s] > 1 && 0 < z < 1]",
        "Gamma[s] PolyLog[s, z]", "Re[s] > 1 && 0 < z < 1");
    /* Equivalent denominator form e^x - z (A_eff = 1 constant): the extra 1/z. */
    assert_cond_closes(
        "Integrate[x^(s-1)/(Exp[x]-z), {x,0,Infinity}, "
        "Method -> \"RamanujanMasterTheorem\", Assumptions -> Re[s] > 1 && 0 < z < 1]",
        "Gamma[s] PolyLog[s, z]/z", "Re[s] > 1 && 0 < z < 1");
    /* Concrete fugacity, symbolic order: 2 Gamma(s) PolyLog(s, 1/2). */
    assert_cond_closes(
        "Integrate[x^(s-1)/(Exp[x]-1/2), {x,0,Infinity}, "
        "Method -> \"RamanujanMasterTheorem\", Assumptions -> Re[s] > 1]",
        "2 Gamma[s] PolyLog[s, 1/2]", "Re[s] > 1");
    /* Soundness of the interval gate: with NO bound on z it cannot prove the
     * fugacity is admissible -> declines (never a wrong value). */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[x^(s-1)/(z^-1 Exp[x]-1), {x,0,Infinity}]",
        "Integrate`RamanujanMasterTheorem");
    /* Provably out-of-range fugacity (z > 2 => |gamma'| > 1, divergent) -> declines. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[x^(s-1)/(Exp[x]-z), {x,0,Infinity}, "
        "Assumptions -> Re[s] > 1 && z > 2]",
        "Integrate`RamanujanMasterTheorem");
}

/* High-fugacity / positive-chemical-potential Fermi-Dirac (gamma' > 1): the
 * DEGENERATE Fermi gas (electrons in metals, white-dwarf matter).  The
 * denominator e^(c x) + gamma' with gamma' > 1 has no interior pole, so the
 * integral converges for Re s > 0 and equals -Gamma(s) c^(-s) PolyLog(s,-gamma')/gamma'
 * by the Fermi-Dirac analytic continuation -- the geometric series diverges here,
 * the integral does not.  These were declined before the fugacity gate was relaxed
 * from the series bound |gamma'| <= 1 to the true convergence bound gamma' >= -1. */
static void test_high_fugacity_fermi(void) {
    /* Concrete, s = 1: Int_0^Inf dx/(e^x + a) = (1/a) Log(1 + a).  a = 2 -> Log[3]/2. */
    assert_closes("Integrate[1/(Exp[x]+2), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Log[3]/2", NULL);
    /* Concrete, s = 2: -PolyLog[2,-2]/2 (matches NIntegrate ~= 0.718373). */
    assert_closes("Integrate[x/(Exp[x]+2), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "-PolyLog[2,-2]/2", NULL);
    /* Larger fugacity, s = 1: a = 5 -> Log[6]/5. */
    assert_closes("Integrate[1/(Exp[x]+5), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Log[6]/5", NULL);
    /* Symbolic degenerate Fermi gas: 1/(e^x + z), z > 0 admitted by the interval
     * gate (gamma' = z > 0 > -1) -> -Gamma(s) PolyLog(s,-z)/z, strip Re s > 0. */
    assert_cond_closes("Integrate[x^(s-1)/(Exp[x]+z), {x,0,Infinity}, "
                       "Method -> \"RamanujanMasterTheorem\", Assumptions -> s > 0 && z > 0]",
                       "-Gamma[s] PolyLog[s,-z]/z", "s > 0 && z > 0");
    /* The divergent Bose side (gamma' < -1) must STILL decline: the relaxation
     * lifted only the upper bound.  x^(s-1)/(e^x - z), z > 2 has gamma' = -z < -1
     * (interior pole at x = Log[z]) -> unevaluated. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[x^(s-1)/(Exp[x]-z), {x,0,Infinity}, "
        "Assumptions -> Re[s] > 1 && z > 2]",
        "Integrate`RamanujanMasterTheorem");
}

/* Frullani integrals: Integrate[(f(a x) - f(b x))/x] = (f(0) - f(Inf)) Log(b/a). */
static void test_frullani(void) {
    /* Classic exp Frullani. */
    assert_closes("Integrate[(Exp[-2 x]-Exp[-5 x])/x, {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Log[5/2]", NULL);
    /* Symbolic scales (output form is Log[b] - Log[a]). */
    assert_closes("Integrate[(Exp[-a x]-Exp[-b x])/x, {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\", Assumptions -> a>0 && b>0]",
                  "Log[b] - Log[a]", "a>0 && b>0");
    /* ArcTan Frullani: f(0)=0, f(Inf)=Pi/2 -> (Pi/2) Log(5/2).  Verified
     * numerically: Simplify does not combine Log[2/5] and Log[5/2] under a Pi
     * coefficient, though the value is exact. */
    assert_eval_eq(
        "Chop[N[Integrate[(ArcTan[5 x]-ArcTan[2 x])/x, {x,0,Infinity}, "
        "Method -> \"RamanujanMasterTheorem\"]] - N[(Pi/2) Log[5/2]]]", "0", 0);
}

/* Log-weighted Mellin: Integrate[x^(s-1) Log[x]^k R(x)] = d^k/ds^k M_R(s). */
static void test_log_weighted(void) {
    /* Odd Log weight over a transform symmetric about its centre -> 0. */
    assert_closes("Integrate[Log[x]/(1+x^2), {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "0", NULL);
    assert_closes("Integrate[Log[x]/(1+x)^2, {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "0", NULL);
    /* Nonzero: Integrate[x Log[x] e^(-x)] = Gamma'(2) = 1 - EulerGamma. */
    assert_closes("Integrate[x Log[x] Exp[-x], {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "1 - EulerGamma", NULL);
    /* Log^2 weight = Pi^3/8 (exact closed form 1/4 Pi PolyGamma[1,1/2]; Mathilda
     * does not reduce PolyGamma[1,1/2] = Pi^2/2, so verify numerically). */
    assert_eval_eq(
        "Chop[N[Integrate[Log[x]^2/(1+x^2), {x,0,Infinity}, "
        "Method -> \"RamanujanMasterTheorem\"]] - N[Pi^3/8]]", "0", 0);
    /* Product of two logs: Integrate[Log[x] Log[1+x^2]/x^2] = Pi.  Here the base
     * kernel is itself a Log (not rational), so the Mellin weight Log[x] pulls a
     * d/ds down onto the Log[1+x^2]/x^2 transform -- a strictly stronger case
     * than the rational-R(x) log-weighted integrals above.  Verified to 22
     * digits against the exact Pi. */
    assert_closes("Integrate[(Log[x] Log[x^2+1])/x^2, {x,0,Infinity}, "
                  "Method -> \"RamanujanMasterTheorem\"]", "Pi", NULL);
    /* The exact reported input, through the full Automatic cascade. */
    assert_closes("Integrate[(Log[x] Log[x^2+1])/x^2, {x,0,Infinity}]", "Pi", NULL);
}

/* Safety: out-of-scope / non-half-line inputs return unevaluated, never wrong,
 * and never hang. */
static void test_declines_cleanly(void) {
    /* Pinned, but not a half-line integral -> strict, no fallback. */
    assert_head_unevaluated(
        "Integrate[Exp[-x], {x, 0, 1}, Method -> \"RamanujanMasterTheorem\"]", "Integrate");
    /* Product of two transcendental kernels = Mellin convolution, out of scope. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[Exp[-x] Cos[x], {x, 0, Infinity}]",
        "Integrate`RamanujanMasterTheorem");
    /* Positivity of the scaling rate unknown (no assumption) -> declines. */
    assert_head_unevaluated(
        "Integrate`RamanujanMasterTheorem[Sin[a x], {x, 0, Infinity}]",
        "Integrate`RamanujanMasterTheorem");
}

/* Assumptions handling: the x-aware undefined-function gate lets a List-bearing
 * HypergeometricPFQ integrand reach the engine (In[6-8]); the Refine-augmented
 * strip discharge collapses a ConditionalExpression once Re[...]/sign
 * assumptions pin the strip; an x-dependent arbitrary function still skips the
 * improper methods and stays unevaluated. */
static void test_mellin_assumptions(void) {
    /* 2F1 / 1F1 / 1F2 now close to the exact Gamma-ratio (gate fix). */
    assert_cond_closes(
        "Integrate[x^(s-1) Hypergeometric2F1[a,b,c,-x], {x,0,Infinity}, "
        "Assumptions -> 0 < Re[s] < Min[Re[a],Re[b]], Method -> \"Mellin\"]",
        "(Gamma[c] Gamma[s] Gamma[a-s] Gamma[b-s])/(Gamma[a] Gamma[b] Gamma[c-s])", "");
    assert_cond_closes(
        "Integrate[x^(s-1) HypergeometricPFQ[{a},{b,c},-x], {x,0,Infinity}, "
        "Assumptions -> 0 < Re[s] < Re[a], Method -> \"Mellin\"]",
        "(Gamma[b] Gamma[c] Gamma[s] Gamma[a-s])/(Gamma[a] Gamma[b-s] Gamma[c-s])", "");
    /* Refine discharge: Re[s]>0 && Re[a]>0 collapses the strip to the bare value. */
    assert_eval_eq(
        "Head[Integrate[x^(s-1) Exp[-a x^2], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, Re[a]>0}, Method -> \"Mellin\"]]", "Times", 0);
    assert_closes(
        "Integrate[x^(s-1) Exp[-a x^2], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, Re[a]>0}, Method -> \"Mellin\"]",
        "(1/2) a^(-s/2) Gamma[s/2]", "Re[s]>0 && Re[a]>0");
    /* Protective: an x-dependent arbitrary function keeps the integral inert. */
    assert_head_unevaluated(
        "Integrate[ff[x] Exp[-x^2], {x,0,Infinity}]", "Integrate");
}

/* The special-function transforms closed by the new recognizers (BesselK, Airy),
 * the Route-A pFq reductions (SinIntegral, StruveH), and the Route-B operational
 * calculus (Erfc, CosIntegral, ExpIntegralE).  Reflection-form results (Si/Ci/
 * StruveH) are checked numerically at a strip-interior point. */
static void test_mellin_special_functions(void) {
    /* BesselK (dedicated recognizer); anchor 2^(s-2) a^-s Gamma[s/2]^2. */
    assert_cond_closes(
        "Integrate[x^(s-1) BesselK[0,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, a>0}, Method -> \"Mellin\"]",
        "2^(s-2) a^(-s) Gamma[s/2]^2", "a>0");
    /* AiryAi (dedicated recognizer); anchor Integrate[AiryAi[x]] = 1/3 at s=1. */
    assert_cond_closes(
        "Integrate[x^(s-1) AiryAi[a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, a>0}, Method -> \"Mellin\"]",
        "a^(-s) Gamma[s]/(3^((s+2)/3) Gamma[(s+2)/3])", "a>0");
    assert_eval_eq(
        "Integrate[x^(s-1) AiryAi[x], {x,0,Infinity}, Assumptions -> Re[s]>0, "
        "Method -> \"Mellin\"] /. s -> 1", "1/3", 0);
    /* Erfc (Route B: K' = Gaussian); anchor Integrate[Erfc[x]] = 1/Sqrt[Pi]. */
    assert_cond_closes(
        "Integrate[x^(s-1) Erfc[a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, a>0}, Method -> \"Mellin\"]",
        "a^(-s) Gamma[(s+1)/2]/(Sqrt[Pi] s)", "a>0");
    /* ExpIntegralE (Route B: K' = -Exp/x); anchor Integrate[ExpIntegralE[1,x]] = 1. */
    assert_cond_closes(
        "Integrate[x^(s-1) ExpIntegralE[1,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, a>0}, Method -> \"Mellin\"]",
        "a^(-s) Gamma[s]/s", "a>0");
    assert_eval_eq(
        "Integrate[x^(s-1) ExpIntegralE[1,x], {x,0,Infinity}, Assumptions -> Re[s]>0, "
        "Method -> \"Mellin\"] /. s -> 1", "1", 0);
    /* SinIntegral (Route A: 1F2), reflection form -> numeric at s = -1/2. */
    assert_cond_num(
        "Integrate[x^(s-1) SinIntegral[a x], {x,0,Infinity}, "
        "Assumptions -> {-1<Re[s]<0, a>0}, Method -> \"Mellin\"]",
        "-a^(-s) Gamma[s] Sin[Pi s/2]/s", "{s -> -1/2, a -> 2}");
    /* CosIntegral (Route B: cos/x), reflection form -> numeric at s = 1/2. */
    assert_cond_num(
        "Integrate[x^(s-1) CosIntegral[a x], {x,0,Infinity}, "
        "Assumptions -> {0<Re[s]<1, a>0}, Method -> \"Mellin\"]",
        "-a^(-s) Gamma[s] Cos[Pi s/2]/s", "{s -> 1/2, a -> 2}");
    /* StruveH (Route A: 1F2), reflection form -> numeric at s = 1/4. */
    assert_cond_num(
        "Integrate[x^(s-1) StruveH[0,a x], {x,0,Infinity}, "
        "Assumptions -> {-1<Re[s]<1/2, a>0}, Method -> \"Mellin\"]",
        "Pi 2^(s-1) a^(-s)/(Cos[Pi s/2] Gamma[1-s/2]^2)", "{s -> 1/4, a -> 2}");
}

/* Single-function extensions (G1/G2): EllipticK, BesselY, upper incomplete
 * Gamma[nu,z], Csch, Sech.  Numeric spot-checks for the reflection/HurwitzZeta
 * forms; symbolic closure where the Gamma ratios coincide. */
static void test_mellin_single_extensions(void) {
    /* EllipticK -> 2F1 (Route A); 0 < Re s < 1/2. */
    assert_cond_closes(
        "Integrate[x^(s-1) EllipticK[-a x], {x,0,Infinity}, "
        "Assumptions -> {0<Re[s]<1/2, a>0}, Method -> \"Mellin\"]",
        "(Gamma[s] Gamma[1/2-s]^2 a^(-s))/(2 Gamma[1-s])", "0<Re[s]<1/2 && a>0");
    /* BesselY -> cos/sin combination of BesselJ (Route A); numeric at interior. */
    assert_cond_num(
        "Integrate[x^(s-1) BesselY[nu,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>Abs[Re[nu]], Re[s]<3/2, a>0}, Method -> \"Mellin\"]",
        "2^(s-1) a^(-s) (Cos[nu Pi] Gamma[(nu+s)/2]/Gamma[1+(nu-s)/2] "
        "- Gamma[(-nu+s)/2]/Gamma[1+(-nu-s)/2])/Sin[nu Pi]",
        "{s -> 1, nu -> 1/3, a -> 2}");
    /* Upper incomplete Gamma[nu, a x] (dedicated recognizer, IBP identity). */
    assert_cond_closes(
        "Integrate[x^(s-1) Gamma[nu, a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, Re[s+nu]>0, a>0}, Method -> \"Mellin\"]",
        "a^(-s) Gamma[s+nu]/s", "Re[s]>0 && Re[s+nu]>0 && a>0");
    /* Csch -> (1-2^-s) Zeta (hyperbolic Dirichlet recognizer); Re s > 1. */
    assert_cond_closes(
        "Integrate[x^(s-1) Csch[a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>1, a>0}, Method -> \"Mellin\"]",
        "2 a^(-s) Gamma[s] (1-2^(-s)) Zeta[s]", "Re[s]>1 && a>0");
    /* Sech -> LerchPhi/beta (hyperbolic Dirichlet); numeric at s = 2. */
    assert_cond_num(
        "Integrate[x^(s-1) Sech[a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>0, a>0}, Method -> \"Mellin\"]",
        "2^(1-s) a^(-s) Gamma[s] LerchPhi[-1,s,1/2]", "{s -> 2, a -> 1}");
}

/* Mellin convolution -- products of two transcendental kernels (G3): the
 * Bessel-product families (Weber-Schafheitlin / Barnes / Kummer) and the
 * exp/Gaussian x BesselJ families (2F1 / 1F1).  Numeric spot-checks at a point
 * strictly interior to the convergence strip. */
static void test_mellin_convolution(void) {
    /* J.J, equal scale (Weber-Schafheitlin -> Gamma ratio via 2F3). */
    assert_cond_num(
        "Integrate[x^(s-1) BesselJ[mu,a x] BesselJ[nu,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s+mu+nu]>0, Re[s]<1, a>0}, Method -> \"Mellin\"]",
        "2^(s-1) a^(-s) Gamma[1-s] Gamma[(mu+nu+s)/2]"
        "/(Gamma[(mu-nu-s)/2+1] Gamma[(nu-mu-s)/2+1] Gamma[(mu+nu-s)/2+1])",
        "{s -> 0, mu -> 1, nu -> 2, a -> 1}");
    /* K.K, equal scale (Barnes' first lemma -> Gamma ratio). */
    assert_cond_num(
        "Integrate[x^(s-1) BesselK[mu,a x] BesselK[nu,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s]>Abs[Re[mu]]+Abs[Re[nu]], a>0}, Method -> \"Mellin\"]",
        "2^(s-3) a^(-s) Gamma[(s+mu+nu)/2] Gamma[(s+mu-nu)/2] "
        "Gamma[(s-mu+nu)/2] Gamma[(s-mu-nu)/2]/Gamma[s]",
        "{s -> 2, mu -> 1/3, nu -> 1/4, a -> 1}");
    /* J.K, equal order and scale (Kummer's theorem -> Gamma ratio). */
    assert_cond_num(
        "Integrate[x^(s-1) BesselJ[nu,a x] BesselK[nu,a x], {x,0,Infinity}, "
        "Assumptions -> {Re[s+2 nu]>0, a>0}, Method -> \"Mellin\"]",
        "2^(s-2) a^(-s) Gamma[s/2] Gamma[1+nu/2+s/4]/((nu+s/2) Gamma[1+nu/2-s/4])",
        "{s -> 3/2, nu -> 1, a -> 2}");
    /* Exp[-a x] BesselJ[nu, b x], distinct scales (-> 2F1 in -b^2/a^2). */
    assert_cond_num(
        "Integrate[x^(s-1) Exp[-a x] BesselJ[nu,b x], {x,0,Infinity}, "
        "Assumptions -> {Re[s+nu]>0, Re[a]>0, b>0}, Method -> \"Mellin\"]",
        "(b/2)^nu Gamma[s+nu]/(a^(s+nu) Gamma[nu+1]) "
        "Hypergeometric2F1[(s+nu)/2,(s+nu+1)/2,nu+1,-b^2/a^2]",
        "{s -> 2, nu -> 1, a -> 3, b -> 2}");
    /* Exp[-a^2 x^2] BesselJ[nu, b x] (Weber's 2nd exp integral -> 1F1). */
    assert_cond_num(
        "Integrate[x^(s-1) Exp[-a^2 x^2] BesselJ[nu,b x], {x,0,Infinity}, "
        "Assumptions -> {Re[s+nu]>0, Re[a^2]>0, b>0}, Method -> \"Mellin\"]",
        "b^nu Gamma[(s+nu)/2]/(2^(nu+1) a^(s+nu) Gamma[nu+1]) "
        "Hypergeometric1F1[(s+nu)/2,nu+1,-b^2/(4 a^2)]",
        "{s -> 1, nu -> 1/2, a -> 2, b -> 3}");
}

void test_integrate_ramanujan(void) {
    symtab_init();
    core_init();

    TEST(test_gaussian_moments);
    TEST(test_exponential_gamma);
    TEST(test_algebraic_binomial);
    TEST(test_trig);
    TEST(test_bessel);
    TEST(test_linearity);
    TEST(test_routing);
    TEST(test_conditional_strip);
    TEST(test_monomial_substitution);
    TEST(test_log_arctan_pfq);
    TEST(test_mellin_assumptions);
    TEST(test_mellin_special_functions);
    TEST(test_mellin_single_extensions);
    TEST(test_mellin_convolution);
    TEST(test_reductions);
    TEST(test_polylog);
    TEST(test_parametric_differentiation);
    TEST(test_bose_fermi);
    TEST(test_symbolic_fugacity);
    TEST(test_high_fugacity_fermi);
    TEST(test_frullani);
    TEST(test_log_weighted);
    TEST(test_declines_cleanly);

    printf("All Integrate Ramanujan Master Theorem tests passed!\n");
}

int main(void) {
    test_integrate_ramanujan();
    return 0;
}

/* test_integrate_roach.c
 *
 * Integrate[..., Method -> "Roach"] (src/internal/RoachIntegrate.m): definite
 * integrals by Meijer G functions -- every factor a G function (PBM vol. 3
 * sec. 8.4), the one- and two-G integration theorems under the PBM 2.24.1
 * conditions, G -> named functions by Roach's Formula algorithm.
 *
 * Cases are the integrals of Kelly Roach, "Solving Integrals with the Quantum
 * Computer Algebra System" (2008), and the ISSAC '97 gallery. Each answer is
 * checked numerically at a sample point against a split-interval NIntegrate,
 * and some are also pinned to Roach's printed form.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* split NIntegrate reference, defined once */
static void define_ref(void) {
    Expr* e = parse_expression(
        "RoachTest`splitN[f_, {x_, a_, b_}] := Module[{p = If[b === Infinity,"
        " Join[{a}, a + {1, 4, 16, 64, 256}, {Infinity}], Subdivide[a, b, 4]]},"
        " Total[Table[NIntegrate[f, {x, p[[i]], p[[i + 1]]}], {i, Length[p] - 1}]]]");
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
    /* Cauchy principal value about the poles ps, by symmetric excision:
     * int_0^d (f(p + u) + f(p - u)) du on [p - d, p + d] (bounded integrand),
     * the rest of [a, b] ordinary */
    e = parse_expression(
        "RoachTest`pvN[f_, {x_, a_, b_}, ps_List] := Module[{d = 1/4, pts, tot = 0, u},"
        " pts = Join[{a}, Flatten[{# - d, # + d} & /@ ps], {b}];"
        " Do[tot += If[pts[[2 i]] === Infinity,"
        "   Total[Table[NIntegrate[f, {x, #[[j]], #[[j + 1]]}],"
        "     {j, Length[#] - 1}]] &[Join[{pts[[2 i - 1]]}, pts[[2 i - 1]] + {1, 4, 16, 64, 256}, {Infinity}]],"
        "   NIntegrate[f, {x, pts[[2 i - 1]], pts[[2 i]]}]], {i, Length[pts]/2}];"
        " Do[tot += NIntegrate[(f /. x -> p + u) + (f /. x -> p - u), {u, 0, d}],"
        "   {p, ps}];"
        " tot]");
    r = evaluate(e);
    expr_free(e); expr_free(r);
}

/* Integrate[f, rng, PrincipalValue -> True, Assumptions -> as(, Method -> "Roach")]
 * closes and agrees at `sub` with the symmetric-excision PV quadrature about
 * the poles `ps` to 1e-8. */
static void check_pv(const char* f, const char* rng, const char* as, const char* sub,
                     const char* ps, int pinned) {
    char buf[4096], opt[512];
    /* no Assumptions option for "True": the PV Newton-Leibniz stage under an
     * explicit Assumptions -> True is slow on ArcTan[x]/(x (1 - x^2)) (Limit) */
    if (strcmp(as, "True") == 0) opt[0] = '\0';
    else snprintf(opt, sizeof opt, ", Assumptions -> %s", as);
    snprintf(buf, sizeof buf,
        "With[{r = Integrate[%s, %s, PrincipalValue -> True%s%s]},"
        " FreeQ[r, Integrate] && With[{v = N[(r /. ConditionalExpression[w_, _] :> w) /. %s],"
        " ref = RoachTest`pvN[Evaluate[(%s) /. %s], Evaluate[%s /. %s], %s]},"
        " Abs[v - ref] < 10^-8 (1 + Abs[ref])]]",
        f, rng, opt, pinned ? ", Method -> \"Roach\"" : "", sub, f, sub, rng, sub, ps);
    assert_eval_eq(buf, "True", 0);
}

/* Cauchy principal values by Roach's PV theorem: the kernel 1/(1 - x) under
 * PV is Pi G^{1,1}_{2,2}(x | 0, 1/2; 0, 1/2) (Mellin transform Pi Cot[Pi w])
 * inside the one- and two-G theorems */
static void test_principal_value(void) {
    /* Roach 2008 #20; the reference value is -Catalan */
    assert_eval_eq("Integrate[ArcTan[x]/(1 - x^2), {x, 0, Infinity}, PrincipalValue -> True,"
                   " Method -> \"Roach\"]", "-Catalan", 0);
    assert_eval_eq("Integrate[ArcTan[x]/(1 - x^2), {x, 0, Infinity}, PrincipalValue -> True]",
                   "-Catalan", 0);
    check_pv("ArcTan[x]/(1 - x^2)", "{x, 0, Infinity}", "True", "{}", "{1}", 1);
    /* the default cascade closes it without a stray Integrate::idiv from the
     * Newton-Leibniz stage that declined before the Roach fallback */
    assert_eval_eq("Check[Integrate[ArcTan[x]/(1 - x^2), {x, 0, Infinity},"
                   " PrincipalValue -> True], failed]", "-Catalan", 0);
    /* one G: the Mellin transform of the kernel */
    assert_eval_eq("Integrate[x^(s - 1)/(1 - x), {x, 0, Infinity}, PrincipalValue -> True,"
                   " Assumptions -> 0 < s < 1, Method -> \"Roach\"] === Pi Cot[Pi s]", "True", 0);
    check_pv("x^(s - 1)/(1 - x)", "{x, 0, Infinity}", "0 < s < 1", "{s -> 1/3}", "{1}", 1);
    check_pv("x^(s - 1)/(1 - x^2)", "{x, 0, Infinity}", "0 < s < 2", "{s -> 3/5}", "{1}", 1);
    check_pv("x^(s - 1)/(a - x)", "{x, 0, Infinity}", "0 < s < 1 && a > 0",
             "{s -> 2/7, a -> 3/2}", "{3/2}", 1);
    /* Log[x] by differentiation in the exponent */
    assert_eval_eq("Integrate[Log[x]/(1 - x^2), {x, 0, Infinity}, PrincipalValue -> True,"
                   " Method -> \"Roach\"] === -Pi^2/4", "True", 0);
    check_pv("Log[x]/(1 - x^2)", "{x, 0, Infinity}", "True", "{}", "{1}", 1);
    /* two G: Bessel and ArcTan against the kernel, through the default cascade */
    check_pv("BesselJ[0, x]/(1 - x^2)", "{x, 0, Infinity}", "True", "{}", "{1}", 0);
    check_pv("ArcTan[x]/(x (1 - x^2))", "{x, 0, Infinity}", "True", "{}", "{1}", 0);
    /* declines: a double pole has no Cauchy PV, ArcTan[x]/(1 - x) diverges at
     * Infinity */
    assert_eval_eq("Head[Integrate[x^(s - 1)/(1 - x)^2, {x, 0, Infinity}, PrincipalValue -> True,"
                   " Assumptions -> 0 < s < 1, Method -> \"Roach\"]]", "Integrate", 0);
    assert_eval_eq("Head[Quiet[Integrate[ArcTan[x]/(1 - x), {x, 0, Infinity},"
                   " PrincipalValue -> True]]]", "Integrate", 0);
}

/* Integrate[f, {x, a, b}, Assumptions -> as, Method -> "Roach"] closes (no
 * Integrate left) and agrees with the split NIntegrate at `sub` to 1e-7. */
static void check_roach(const char* f, const char* rng, const char* as, const char* sub) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{r = Integrate[%s, %s, Assumptions -> %s, Method -> \"Roach\"]},"
        " FreeQ[r, Integrate] && With[{v = N[(r /. ConditionalExpression[w_, _] :> w) /. %s],"
        " ref = RoachTest`splitN[Evaluate[(%s) /. %s], Evaluate[%s /. %s]]},"
        " Abs[v - ref] < 10^-7 (1 + Abs[ref])]]",
        f, rng, as, sub, f, sub, rng, sub);
    assert_eval_eq(buf, "True", 0);
}

/* Roach 2008, p. 3-4 */
static void test_paper(void) {
    check_roach("Exp[-al x] BesselJ[0, be x]^2", "{x, 0, Infinity}", "al > 0 && be > 0", "{al -> 2, be -> 1/3}");
    check_roach("(a^2 - x^2)^(n - 1/2)", "{x, 0, a}", "a > 0 && n > 0", "{a -> 2, n -> 1/3}");
    check_roach("Exp[-b x - x^2/(4 a)]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x Exp[-x n]/(Exp[x] + 1)", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
    check_roach("ArcSin[x] x^(2 n - 1)", "{x, 0, 1}", "n > 0", "{n -> 2/3}");
    check_roach("x^(a - 1)/Sqrt[x^b + 1]", "{x, 0, Infinity}", "b > 0 && a > 0 && a < b/2", "{a -> 1/3, b -> 3}");
    check_roach("x Exp[-2 x n]/(Exp[x] + 1)", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
    check_roach("x Exp[-3 x]/(Exp[-x] + 1)", "{x, 0, Infinity}", "True", "{}");
    check_roach("x^(b - 1) (a x + 1)^(-c)", "{x, 0, Infinity}", "a > 0 && b > 0 && c > b", "{a -> 2, b -> 1/3, c -> 7/5}");
    check_roach("(SinIntegral[a x] - Pi/2) Exp[-b x]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x^(a - 1) (b + x)^(-n)", "{x, 0, Infinity}", "b > 0 && a > 0 && n > a", "{a -> 1/3, b -> 2, n -> 7/5}");
    check_roach("BesselJ[1, a x] Exp[-x^2 b]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("Exp[-a x] x^n/(b + x)", "{x, 0, Infinity}", "a > 0 && b > 0 && n > -1", "{a -> 2, b -> 1/3, n -> 2/3}");
    check_roach("Log[x + 1] Cos[a x]/x", "{x, 0, Infinity}", "a > 0", "{a -> 2}");
    check_roach("x Erfc[a x] Exp[-x^2 b]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x/(Exp[3 x] - 1)^(1/3)", "{x, 0, Infinity}", "True", "{}");
    check_roach("x Exp[-x]/Sqrt[Exp[x] - 1]", "{x, 0, Infinity}", "True", "{}");
    check_roach("x^2 Csch[x] Exp[-2 x n]", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
}

/* t = E^-x with Log[t]^k by differentiation in the exponent, and the Lerch
 * term for 1/(1 -+ t^m): Roach's printed values */
static void test_exponential(void) {
    assert_eval_eq("Integrate[x/(Exp[3 x] - 1)^(1/3), {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Pi (Pi + 3 Sqrt[3] Log[3])/27", "True", 0);
    assert_eval_eq("Integrate[x Exp[-x]/Sqrt[Exp[x] - 1], {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Pi Log[2] - Pi/2", "True", 0);
    assert_eval_eq("Integrate[x^2 Csch[x] Exp[-2 x n], {x, 0, Infinity}, Assumptions -> n > 0,"
                   " Method -> \"Roach\"] === Zeta[3, n + 1/2]/2", "True", 0);
}

/* Roach's printed forms */
static void test_forms(void) {
    assert_eval_eq("Integrate[(SinIntegral[a x] - Pi/2) Exp[-b x], {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0, Method -> \"Roach\"] === (ArcTan[a/b] - Pi/2)/b", "True", 0);
    assert_eval_eq("Integrate[BesselJ[1, a x] Exp[-x^2 b], {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0, Method -> \"Roach\"] === (1 - E^(-a^2/(4 b)))/a", "True", 0);
    assert_eval_eq("Integrate[x^(b - 1) (a x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0 && c > b, Method -> \"Roach\"] ==="
                   " Gamma[b] Gamma[c - b]/(a^b Gamma[c])", "True", 0);
    /* #9 by t = Tan[x]: (1 - t^2)^a on [0, 1].  The exponents PowerExpand
     * leaves are expanded, so (1 + t^2)^(-1 - a + (2 + 2 a)/2) is 1 and not a
     * spurious second G (whose convolution carried Gamma[0] and was rejected
     * by the numeric check). */
    assert_eval_eq("Integrate[Cos[x]^(-2 a - 2) Cos[2 x]^a, {x, 0, Pi/4},"
                   " Assumptions -> a > 0, Method -> \"Roach\"] ==="
                   " Sqrt[Pi] Gamma[1 + a]/(2 Gamma[3/2 + a])", "True", 0);
}

/* a symbolic finite upper limit z: the rescale x -> z y gives Log[z y], and a
 * G-form carrying theta(1 - z y) must not stand in for theta(1 - y) -- it
 * returned the value at z = 1 (Mathar M6460), which the numeric check passed
 * because FindInstance chose z -> 1 */
static void test_finite_symbolic_limit(void) {
    check_roach("Log[t]/(1 + a^2 t^2)", "{t, 0, z}", "a > 0 && z > 0", "{a -> 7/3, z -> 7/4}");
    check_roach("Log[t]/(1 + a^2 t^2)", "{t, 0, z}", "a > 0 && z > 0", "{a -> 7/4, z -> 5/4}");
    assert_eval_eq("FreeQ[Integrate[Log[t]/(1 + a^2 t^2), {t, 0, z},"
                   " Assumptions -> a > 0 && z > 0, Method -> \"Roach\"], z]", "False", 0);
}

/* Roach 2008 #18: the convolution is G^{3,4}_{4,4} at 1 in the log case,
 * reduced by the residue sum of its trigonometric Mellin-Barnes integrand
 * (Hrep`unitG). #24: a growing exponential against Ei(-x^2), by analytic
 * continuation from the decaying case (Roach`int2grow). */
static void test_unit_and_growing(void) {
    check_roach("ArcCot[x]/(x + 1)", "{x, 0, Infinity}", "True", "{}");
    check_roach("ArcTan[x]/(x (x + 1))", "{x, 0, Infinity}", "True", "{}");
    check_roach("x ArcCot[x]/((1 + x) (1 + x^2))", "{x, 0, Infinity}", "True", "{}");
    assert_eval_eq("Integrate[ArcCot[x]/(x + 1), {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Catalan + Pi Log[2]/4", "True", 0);
    assert_eval_eq("Integrate[ArcCot[x]/(x + 1), {x, 0, Infinity}] === Catalan + Pi Log[2]/4", "True", 0);
    check_roach("ExpIntegralEi[-x^2] Exp[a^2 x^2]", "{x, 0, Infinity}", "0 < a < 1", "{a -> 1/3}");
    check_roach("ExpIntegralEi[-x] Exp[b x] x", "{x, 0, Infinity}", "0 < b < 1", "{b -> 1/2}");
    check_roach("Erfc[x] Exp[b x^2]", "{x, 0, Infinity}", "0 < b < 1", "{b -> 1/2}");
    assert_eval_eq("Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> 0 < a < 1, Method -> \"Roach\"] === -Sqrt[Pi] ArcSin[a]/a", "True", 0);
    assert_eval_eq("Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> 0 < a < 1] === -Sqrt[Pi] ArcSin[a]/a", "True", 0);
    /* divergent: the exponential grows faster than Ei(-x^2) decays */
    assert_eval_eq("Head[Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> a > 1, Method -> \"Roach\"]]", "Integrate", 0);
}

/* Adamchik and Marichev, "The algorithm for calculating integrals of
 * hypergeometric type functions and its realization in REDUCE system"
 * (ISSAC '90). Section 5's demonstration, its p = v + 2 special case, the
 * indefinite / variable-upper-limit formula (26) with the paper's eq. (27),
 * and the table entries of eq. (17) and PBM 8.4.23. */
static void check_antiderivative(const char* f, const char* sub) {
    /* closes, is free of MeijerG, and differentiates back to f at `sub` */
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{F = Integrate[%s, x, Method -> \"Roach\"]},"
        " FreeQ[F, Integrate | MeijerG] &&"
        " Abs[N[(D[F, x] - (%s)) /. %s, 30]] < 10^-15]", f, f, sub);
    assert_eval_eq(buf, "True", 0);
}
static void test_adamchik_marichev(void) {
    /* sec. 5 demo, under the DEFAULT cascade (it timed out there: Newton-
     * Leibniz and DiffUnderInt ground on it; the early Roach slot closes it) */
    const char* demo = "Integrate[x^(p - 1) (1 - x^2)^(b - 1) BesselJ[v, c x], {x, 0, 1},"
                       " Assumptions -> c > 0 && b > 0 && p + v > 0]";
    char buf[2048];
    snprintf(buf, sizeof buf,
        "TimeConstrained[%s, 30, $Aborted] === (2^(-1 - v) c^v Gamma[b] Gamma[(p + v)/2]"
        " HypergeometricPFQ[{(p + v)/2}, {(2 b + p + v)/2, 1 + v}, -c^2/4])/"
        "(Gamma[1 + v] Gamma[(2 b + p + v)/2])", demo);
    assert_eval_eq(buf, "True", 0);
    check_roach("x^(p - 1) (1 - x^2)^(b - 1) BesselJ[v, c x]", "{x, 0, 1}",
                "c > 0 && b > 0 && p + v > 0", "{p -> 5/3, b -> 3/2, v -> 1/3, c -> 7/4}");
    /* p = v + 2: the paper's 2^(b-1) Gamma(b) J_(b+v)(c)/c^b */
    assert_eval_eq("Integrate[x^(v + 1) (1 - x^2)^(b - 1) BesselJ[v, c x], {x, 0, 1},"
                   " Assumptions -> c > 0 && b > 0 && v > -1] === 2^(b - 1) Gamma[b] BesselJ[b + v, c]/c^b",
                   "True", 0);
    /* eq. 27: a variable upper limit, conditioned on z >= 0 (the G forms are
     * x > 0 forms: Erf[t] = G(t^2)/Sqrt[Pi] is even in t) */
    check_roach("Sin[x]/Sqrt[x]", "{x, 0, z}", "z > 0", "{z -> 7/4}");
    check_roach("x^2 Erf[x]", "{x, 0, z}", "z > 0", "{z -> 7/4}");
    check_roach("x^(1/3) BesselJ[1, x]", "{x, 1/2, z}", "z > 1/2", "{z -> 9/4}");
    assert_eval_eq("Head[Integrate[x^2 Erf[x], {x, 0, z}, Method -> \"Roach\"]]",
                   "ConditionalExpression", 0);
    /* eq. 26 as an antiderivative */
    check_antiderivative("Sin[x]/Sqrt[x]", "{x -> 7/4}");
    check_antiderivative("x^a BesselJ[v, x]", "{x -> 7/4, a -> 1/3, v -> 2/5}");
    check_antiderivative("x^(1/3) Erfc[x]", "{x -> 7/4}");
    check_antiderivative("x BesselJ[a, x] BesselJ[b, x]", "{x -> 7/4, a -> 1/3, b -> 2/5}");
    check_antiderivative("x^2 Exp[-x^3/2]", "{x -> 7/4}");
    /* the antiderivative of an odd-even pair must hold for x < 0 too */
    assert_eval_eq("With[{F = Integrate[x^2 Erf[x], x, Method -> \"Roach\"]},"
                   " FreeQ[F, Integrate] && Abs[N[(D[F, x] - x^2 Erf[x]) /. x -> -7/4, 30]] < 10^-15]",
                   "True", 0);
    /* dx/x: the G is undefined (a Gamma pole); the power rule gives Log */
    assert_eval_eq("Integrate[1/x, x, Method -> \"Roach\"]", "Log[x]", 0);
    /* squares and products of one argument (PBM 8.4.23): two G functions */
    assert_eval_eq("Integrate[BesselK[0, x]^2, {x, 0, Infinity}] === Pi^2/4", "True", 0);
    check_roach("BesselK[0, x]^2 x^s", "{x, 0, Infinity}", "s > 0", "{s -> 2/3}");
    check_roach("x^2 Exp[-x] BesselK[0, x]^2", "{x, 0, Infinity}", "True", "{}");
    check_roach("Exp[-x] BesselJ[0, x] BesselJ[1, x]", "{x, 0, Infinity}", "True", "{}");
    check_roach("Erfc[x]^2", "{x, 0, Infinity}", "True", "{}");
    /* eq. 17 and the confluent functions as integrands */
    check_roach("x^(s - 1) Hypergeometric2F1[a, b, c, -x]", "{x, 0, Infinity}",
                "0 < s < a && s < b && c > 0", "{s -> 1/3, a -> 3/2, b -> 7/4, c -> 5/2}");
    check_roach("Exp[-x] Hypergeometric1F1[a, c, -x]", "{x, 0, Infinity}",
                "a > 0 && c > 0", "{a -> 1/3, c -> 5/2}");
    check_roach("HypergeometricU[a, b, x] x^(s - 1)", "{x, 0, Infinity}",
                "s > 0 && s < a && s > b - 1 && a > 0", "{s -> 1/2, a -> 3/2, b -> 1/3}");
    /* a MeijerG integrand (it evaluates to U and Exp before Integrate sees it) */
    check_roach("MeijerG[{{1/2}, {}}, {{0, 0}, {}}, x] MeijerG[{{}, {}}, {{1/3}, {}}, x^(3/2)]",
                "{x, 0, Infinity}", "True", "{}");
}

/* the stage declines (stays unevaluated) rather than guessing */
static void test_declines(void) {
    assert_eval_eq("Head[Integrate[Sin[x^x], {x, 0, 1}, Method -> \"Roach\"]]", "Integrate", 0);
    /* divergent: x^(b-1) (1 + x)^-c with c < b */
    assert_eval_eq("Head[Integrate[x^(b - 1) (x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> b > 0 && c > 0 && c < b, Method -> \"Roach\"]]", "Integrate", 0);
}

/* the stage leaves no global state behind (pattern capture, bindings) */
static void test_hygiene(void) {
    assert_eval_eq("Integrate[Exp[-al x] BesselJ[0, be x]^2, {x, 0, Infinity},"
                   " Assumptions -> al > 0 && be > 0, Method -> \"Roach\"]; {ValueQ[x], ValueQ[al], ValueQ[be]}",
                   "{False, False, False}", 0);
    assert_eval_eq("Check[Integrate[x^(b - 1) (a x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0 && c > b, Method -> \"Roach\"]; ok, failed]", "ok", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    define_ref();
    TEST(test_paper);
    TEST(test_forms);
    TEST(test_exponential);
    TEST(test_finite_symbolic_limit);
    TEST(test_unit_and_growing);
    TEST(test_declines);
    TEST(test_hygiene);
    TEST(test_principal_value);
    TEST(test_adamchik_marichev);
    printf("All Integrate Roach tests passed!\n");
    return 0;
}

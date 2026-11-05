/* test_hrep.c
 *
 * Roach's Formula algorithm for pFq (src/special_functions/hrep_*.c):
 * FunctionExpand[HypergeometricPFQ[...]] and
 * HypergeometricPFQ`Representation.
 *
 * Correctness is asserted numerically. Every closed form must agree with
 * the defining series to 25 digits at a sample point, so the tests survive
 * surface-form changes. A few printed forms, Roach's own, are pinned
 * exactly as well.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* The closed form of `f` (no HypergeometricPFQ left) agrees with f at
 * z -> z0 under the parameter rules `sub` to 25 digits. */
static void check_closed(const char* f, const char* sub, const char* z0) {
    char buf[2048];
    snprintf(buf, sizeof buf,
        "With[{r = FunctionExpand[%s]}, FreeQ[r, HypergeometricPFQ] && "
        "Abs[N[(r - (%s)) /. %s /. z -> %s, 40]] < 10^-25]", f, f, sub, z0);
    assert_eval_eq(buf, "True", 0);
}

/* Roach '96 sec. 2: 2F1(-3/2, -1/2; 1/2; z) walked from 2F1(1/2, 1/2; 3/2; z). */
static void test_lead_example(void) {
    assert_eval_eq("FunctionExpand[Hypergeometric2F1[-3/2, -1/2, 1/2, z]]",
        "1/2 Sqrt[1 - z] (2 + z) + 3/2 Sqrt[z] ArcSin[Sqrt[z]]", 0);
    check_closed("Hypergeometric2F1[-3/2, -1/2, 1/2, z]", "{}", "3/10");
    check_closed("Hypergeometric2F1[-3/2, -1/2, 1/2, z]", "{}", "-7/3");
}

/* The ISSAC '96 gallery entries within the Phase 1 certificate set. */
static void test_gallery(void) {
    /* F(1/2; 9/2; z): exp and erfi */
    check_closed("Hypergeometric1F1[1/2, 9/2, z]", "{}", "3/10");
    check_closed("Hypergeometric1F1[1/2, 9/2, z]", "{}", "-5/2");
    /* F(-3/2; -1/2, 1/2; z) = (1 + 2z) cosh 2 sqrt z + sqrt z sinh 2 sqrt z - 4 z^(3/2) Shi */
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{-3/2}, {-1/2, 1/2}, z]]",
        "(1 + 2 z) Cosh[2 Sqrt[z]] + Sqrt[z] Sinh[2 Sqrt[z]] - 4 z^(3/2) SinhIntegral[2 Sqrt[z]]", 0);
    check_closed("HypergeometricPFQ[{-3/2}, {-1/2, 1/2}, z]", "{}", "3/10");
    /* F(-3/2, -1/2; 2; z) in K and E */
    check_closed("Hypergeometric2F1[-3/2, -1/2, 2, z]", "{}", "3/10");
    check_closed("Hypergeometric2F1[-3/2, -1/2, 2, z]", "{}", "-4");
}

/* Walks in every direction: raise and lower, upper and lower parameters. */
static void test_walk_directions(void) {
    check_closed("Hypergeometric2F1[-1/2, 1/2, 3/2, z]", "{}", "3/10");   /* a down */
    check_closed("Hypergeometric2F1[5/2, 1/2, 7/2, z]", "{}", "3/10");    /* a up, b up */
    check_closed("Hypergeometric2F1[1/2, 1/2, -5/2, z]", "{}", "3/10");   /* b down */
    check_closed("Hypergeometric1F1[3, 5, z]", "{}", "7/3");
    check_closed("Hypergeometric2F1[1, 1, 7/2, z]", "{}", "3/10");
    check_closed("HypergeometricPFQ[{}, {7/2}, z]", "{}", "3/10");
    check_closed("HypergeometricPFQ[{2, 3}, {1}, z]", "{}", "3/10");
    /* long walks */
    check_closed("Hypergeometric2F1[1/2, 1/2, 41/2, z]", "{}", "3/10");
    check_closed("Hypergeometric1F1[1/2, 61/2, z]", "{}", "3/10");
}

/* Derivative formula: an upper parameter exceeds a lower one by an
 * integer, so the order drops. */
static void test_derivative_formula(void) {
    check_closed("HypergeometricPFQ[{3/2, 2}, {1, 5/2}, z]", "{}", "3/10");
    check_closed("HypergeometricPFQ[{5/2, 1/2, 4}, {7/2, 2}, z]", "{}", "1/5");
}

/* Symbolic parameters: family certificates. */
static void test_symbolic(void) {
    check_closed("Hypergeometric2F1[a, a + 1/2, 5/2, z]", "a -> 1/7", "3/10");
    check_closed("Hypergeometric1F1[n, n + 3, z]", "n -> 2/9", "-3/10");
    check_closed("HypergeometricPFQ[{}, {b + 2}, z]", "b -> 1/3", "3/10");
    check_closed("Hypergeometric2F1[1, b, b + 3, z]", "b -> 2/5", "3/10");
    check_closed("Hypergeometric2F1[a, -a + 2, 1/2, z]", "a -> 1/3", "3/10");
}

/* Every lookup certificate against its defining series: family parameters
 * at sample rationals, w at 1/7 and -1/5, to 25 digits. This is where a
 * transcription error in a table entry is caught. */
static void test_certificates(void) {
    assert_eval_eq(
        "With[{rows = HypergeometricPFQ`Certificates[]},"
        " Select[rows, Function[r, Module[{s = {Hrep`al -> 3/11, Hrep`be -> 5/13}, f, ok = True},"
        "  If[ReleaseHold[r[[4]] /. s] === False, False,"
        "   f = Hrep`Inamed[ReleaseHold[r[[3]] /. s]];"
        "   Do[ok = ok && Abs[N[(f - HypergeometricPFQ[r[[1]] /. s, r[[2]] /. s, Hrep`w]) /. Hrep`w -> w0, 40]] < 10^-25,"
        "    {w0, {1/7, -1/5}}]; !ok]]]]]", "{}", 0);
}

/* Generated certificates (integral formula), the Lerch Phi formula, the
 * transformations, and the ISSAC '96 gallery entries they reach. */
static void test_generated(void) {
    /* integral formula: F(-1/2, 1, 2; 3, 4; z), Roach's gallery, elementary + log */
    check_closed("HypergeometricPFQ[{-1/2, 1, 2}, {3, 4}, z]", "{}", "3/10");
    check_closed("HypergeometricPFQ[{-1/2, 1, 2}, {3, 4}, z]", "{}", "-2/7");
    /* Lerch Phi: dilogarithms of sqrt z */
    check_closed("HypergeometricPFQ[{-1/2, 1/2, 1}, {3/2, 5/2}, z]", "{}", "3/10");
    check_closed("HypergeometricPFQ[{-1/2, 1/2, 1}, {3/2, 5/2}, z]", "{}", "-2/7");
    check_closed("HypergeometricPFQ[{1, 1, 1}, {2, 2}, z]", "{}", "3/10");
    /* Clausen certificate: ArcSin^2 */
    check_closed("HypergeometricPFQ[{1, 2, 3}, {1/2, 4}, z]", "{}", "3/10");
    /* 1F2(1; 2, b) family and the Bessel kernel: J0, J1 of 2 sqrt z */
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{1, 3}, {2, 4, 5}, -z]] ==="
        " (6 (-12 + z^2))/z^3 - (72 BesselJ[0, 2 Sqrt[z]])/z^3 + (144 BesselJ[1, 2 Sqrt[z]])/z^(7/2)", "True", 0);
    check_closed("HypergeometricPFQ[{1, 3}, {2, 4, 5}, -z]", "{}", "3/10");
    /* 2F2 by the derivative formula to 1F1(a; 2a) */
    check_closed("HypergeometricPFQ[{-3/2, -1/2}, {-5/2, 1}, z]", "{}", "3/10");
    /* Pfaff image into the LerchPhi family */
    check_closed("Hypergeometric2F1[1, 1, 2/3, z]", "{}", "3/10");
    check_closed("Hypergeometric2F1[1, 1, 2/3, z]", "{}", "-2/7");
}

/* Roach 2008 p. 6: values at z = 1, -1, 9/10, 8/9. */
static void test_values(void) {
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{-15/4, 1, 1}, {-7/2, 2}, 1]] === (2252 + 315 Pi + 630 Log[2])/665", "True", 0);
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{1, 1, 3/2}, {1/2, 9/2}, -1]] === 1841/15 - 98 Sqrt[2] ArcSinh[1]", "True", 0);
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{-11/4, -5/2, -5/2}, {-7/2, 1/4}, 1]] === (355 Gamma[1/4]^2)/(56 Sqrt[2] Sqrt[Pi])", "True", 0);
    assert_eval_eq("Abs[N[FunctionExpand[HypergeometricPFQ[{-1/2, -1/2, 5/2}, {3/2, 7/2}, -1]]"
                   " - (805/(9 2^(13/2)) - 15 Log[Sqrt[2] + 1]/128), 30]] < 10^-25", "True", 0);
    /* Kummer at -1 against the series summed directly (alternating, converges) */
    assert_eval_eq("Abs[N[FunctionExpand[Hypergeometric2F1[-26/7, -26/7, 1, -1]] -"
                   " 95 Gamma[1/7] Gamma[5/14] Cos[3 Pi/14]/(39 2^(2/7) Pi^(3/2)), 30]] < 10^-25", "True", 0);
    assert_eval_eq("Abs[N[FunctionExpand[Hypergeometric2F1[1, 1, 2/3, 8/9]] -"
                   " (3^(3/2) Pi + 9 - 3 Log[3]), 30]] < 10^-25", "True", 0);
    assert_eval_eq("FunctionExpand[Hypergeometric2F1[1/2, 1/2, 1, 9/10]] === (2 EllipticK[9/10])/Pi", "True", 0);
    /* never Indeterminate at a removable singularity */
    assert_eval_eq("FreeQ[FunctionExpand[HypergeometricPFQ[{-11/4, -5/2, -5/2}, {-7/2, 1/4}, 1]],"
                   " Indeterminate]", "True", 0);
}

/* The engine is a probe region: its internal Integrate/Limit messages do not
 * reach the user or an enclosing Check. */
static void test_quiet(void) {
    assert_eval_eq("Check[FunctionExpand[HypergeometricPFQ[{-1/2, 1/2, 1}, {3/2, 5/2}, z]]; ok, failed]",
                   "ok", 0);
}

/* The G walk (Roach '97), including the log case the Slater route cannot do. */
static void check_g(const char* g) {
    char buf[1024];
    snprintf(buf, sizeof buf,
        "With[{r = FunctionExpand[%s]}, FreeQ[r, MeijerG | HypergeometricPFQ] &&"
        " Abs[N[(r - (%s)) /. z -> 7/10, 30]] < 10^-20 && Abs[N[(r - (%s)) /. z -> 3, 30]] < 10^-20]",
        g, g, g);
    assert_eval_eq(buf, "True", 0);
}
static void test_meijer(void) {
    assert_eval_eq("FunctionExpand[MeijerG[{{}, {1}}, {{0, 0}, {}}, z]] === Gamma[0, z]", "True", 0);
    assert_eval_eq("FunctionExpand[MeijerG[{{}, {-1/2}}, {{0, 0}, {-1/2}}, z]] === BesselY[0, 2 Sqrt[z]]", "True", 0);
    check_g("MeijerG[{{}, {2}}, {{1, 0}, {}}, z]");
    check_g("MeijerG[{{}, {1}}, {{0, 0}, {1/2}}, z]");
    check_g("MeijerG[{{}, {2}}, {{1, 1}, {1/2}}, z]");
    check_g("MeijerG[{{}, {1}}, {{1/2, 0}, {}}, z]");
    check_g("MeijerG[{{}, {2}}, {{3/2, 0}, {}}, z]");
}

/* G^{m,n}_{p,p}(1) by residue summation: the Mellin-Barnes integrand is
 * trigonometric times rational, the residues a rational function of k
 * summed as Lerch Phi at +-1. */
static void test_unit_g(void) {
    /* the convolution of ArcCot[x] with 1/(1 + x) (Roach 2008, #18): log case */
    assert_eval_eq("FunctionExpand[MeijerG[{{0, 0, 1/2, 1/2}, {}}, {{0, 0, 1/2}, {-1/2}}, 1]] ==="
                   " Pi (4 Catalan + Pi Log[2])", "True", 0);
    assert_eval_eq("FunctionExpand[MeijerG[{{0, 1/2, 1/2}, {1}}, {{0, 0, 1/2, 1/2}, {}}, 1]] ==="
                   " Pi (4 Catalan + Pi Log[2])", "True", 0);
    /* ArcCot's own G at 1: 2 ArcCot[1] */
    assert_eval_eq("FunctionExpand[MeijerG[{{1/2}, {1}}, {{1/2, 0}, {}}, 1]] === Pi/2", "True", 0);
    assert_eval_eq("With[{g = MeijerG[{{1/3}, {2/3}}, {{0, 1/3}, {}}, 1]},"
                   " Abs[N[FunctionExpand[g] - g, 20]] < 10^-15]", "True", 0);
    /* the Slater terms diverge at 1 (2F1(2/3, 1/3; 2/3; 1)): no ComplexInfinity */
    assert_eval_eq("FreeQ[FunctionExpand[MeijerG[{{1/3, 2/3}, {}}, {{0, 1/3}, {}}, 1]], ComplexInfinity]",
                   "True", 0);
    /* not trigonometric (an unpaired Gamma): declines */
    assert_eval_eq("Head[FunctionExpand[MeijerG[{{1/2, 1/2}, {1}}, {{1/2, 0, 0}, {}}, 1]]]", "MeijerG", 0);
}

/* Struve L certificates and integer-order reduction: Roach's gallery
 * 1F2(3/2; 5/2, 5; z) in I0, I1, L0, L1 */
static void test_struve(void) {
    check_closed("HypergeometricPFQ[{3/2}, {5/2, 5}, z]", "{}", "7/10");
    check_closed("HypergeometricPFQ[{3/2}, {5/2, 5}, z]", "{}", "-3");
    assert_eval_eq("FreeQ[FunctionExpand[HypergeometricPFQ[{3/2}, {5/2, 5}, z]],"
                   " (BesselI | StruveL)[n_ /; n > 1 || n < 0, _]]", "True", 0);
    check_closed("HypergeometricPFQ[{1}, {3/2, 7/3}, z]", "{}", "7/10");
}

/* Declines stay unevaluated. */
static void test_declines(void) {
    assert_eval_eq("FunctionExpand[HypergeometricPFQ[{a}, {c}, z]]",
                   "HypergeometricPFQ[{a}, {c}, z]", 0);
    assert_eval_eq("HypergeometricPFQ`Representation[{1/3}, {1/5}, z]",
                   "HypergeometricPFQ`Representation[{1/3}, {1/5}, z]", 0);
    /* automatic evaluation is unchanged */
    assert_eval_eq("Hypergeometric2F1[-3/2, -1/2, 1/2, z]",
                   "HypergeometricPFQ[{-3/2, -1/2}, {1/2}, z]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_lead_example);
    TEST(test_gallery);
    TEST(test_walk_directions);
    TEST(test_derivative_formula);
    TEST(test_symbolic);
    TEST(test_certificates);
    TEST(test_generated);
    TEST(test_values);
    TEST(test_quiet);
    TEST(test_meijer);
    TEST(test_unit_g);
    TEST(test_struve);
    TEST(test_declines);
    printf("All hrep tests passed!\n");
    return 0;
}

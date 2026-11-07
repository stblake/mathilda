/* Tests for RiemannR[x], Riemann's prime-counting function (Gram series).
 *
 * Reference values: mpmath riemannr / an independent mpmath Gram-series sum at
 * a few hundred digits (the 0 < x < 1 cases cancel through e^|Log x| and need
 * that), and the SIAM 100-Digit Challenge problem D4 truth (the largest zero
 * of R, t* = Log x* = -34142.12818460646499263...).
 *
 * Covers: exact arguments (RiemannR[1] = 1, others unevaluated), machine real
 * on both sides of 1 including the double-double accuracy at x = 10^300,
 * machine complex and negative x, arbitrary precision in the cancelling
 * region, the derivative (symbolic inert, numeric value), the four surfaces
 * (scalar, packed, NDArray, Compile) agreeing bit for bit, attributes, and the
 * D4 root by FindRoot. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        printf("RiemannR test failed: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    }
    free(s);
    expr_free(e);
    expr_free(r);
}

static void test_exact(void) {
    check("RiemannR[1]", "1");
    check("RiemannR[2]", "RiemannR[2]");
    check("RiemannR[0]", "RiemannR[0]");
    check("RiemannR[3/2]", "RiemannR[3/2]");
    check("RiemannR[x]", "RiemannR[x]");
    check("RiemannR[{1, y}]", "{1, RiemannR[y]}");
}

/* Machine reals: relative error at the ulp level (the double-double series). */
static void test_machine_real(void) {
    check("Abs[RiemannR[2.] - 1.541009016187131883288503786627546543531] < 4*^-16", "True");
    check("Abs[RiemannR[0.5] - 0.6635262381124574212296788759118420698842] < 2*^-16", "True");
    check("Abs[RiemannR[0.3] - 0.497534732335077013642737508401322083607] < 2*^-16", "True");
    check("Abs[RiemannR[1.*^6]/78527.39942912770485887029214095925103488 - 1] < 4*^-16", "True");
    /* a plain log(x) would cost ~300 ulps here: dR/R ~ dLog x */
    check("Abs[RiemannR[1.*^300]/1.449750052669336365059039839805422974402*^297 - 1] < 4*^-16", "True");
    /* below 1/e: the MPFR path at 53 bits, rounded */
    check("Abs[RiemannR[1.*^-5] - 0.01190373103082508357472495811321459040983] < 1*^-17", "True");
    check("RiemannR[0.]", "0.0");
    check("Head[RiemannR[2.]]", "Real");
}

static void test_complex(void) {
    check("Abs[RiemannR[-1.] - (-0.2678760993850765965240829186555242505945 + "
          "0.7428692097276307001108319448622983625228 I)] < 1*^-15", "True");
    check("Abs[RiemannR[1. + I] - (1.076837233762951963425840851534240332915 + "
          "0.5746065575833807851708724266904059824519 I)] < 1*^-15", "True");
    check("Abs[RiemannR[N[1 + I, 40]] - (1.076837233762951963425840851534240332915 + "
          "0.5746065575833807851708724266904059824519 I)] < 10^-38", "True");
}

/* Arbitrary precision, including deep cancellation at 0 < x < 1. */
static void test_mpfr(void) {
    check("Abs[RiemannR[N[100, 40]] - 25.66163326692418259322679794035569814997`45] < 10^-37", "True");
    check("Precision[RiemannR[N[100, 40]]]", "40.0");
    check("Abs[RiemannR[N[1/3, 40]] - 0.5274463980337391768005332626706584002093`45] < 10^-39", "True");
    /* Log x = -1000 and -3000: the terms reach e^1000 and e^3000 before
     * cancelling to 10^-8 and 10^-10 (independent Gram sum at 1400 digits) */
    check("Abs[RiemannR[N[Exp[-1000], 30]]/2.189126073708485811597206`25*^-8 - 1] < 10^-24", "True");
    check("Abs[RiemannR[N[Exp[-3000], 30]]/8.127624007363882862396705`25*^-10 - 1] < 10^-24", "True");
    check("Abs[RiemannR[N[-1, 40]] - (-0.2678760993850765965240829186555242505945 + "
          "0.7428692097276307001108319448622983625228 I)] < 10^-38", "True");
}

static void test_derivative(void) {
    check("D[RiemannR[x], x]", "Derivative[1][RiemannR][x]");
    /* R'(2) = 0.492873059153156128131756166775 (mpmath diff of riemannr) */
    check("Abs[Derivative[1][RiemannR][2.] - 0.492873059153156128131756166775] < 1*^-15", "True");
    check("Abs[Derivative[1][RiemannR][N[2, 30]] - 0.492873059153156128131756166775] < 10^-28", "True");
    check("Abs[D[RiemannR[x], x] /. x -> 0.1] > 0", "True");
}

/* scalar == unpacked List == packed List == NDArray == Compile, bit for bit. */
static void test_surfaces(void) {
    check("xs = {0.5, 2., 10., 1.*^6, 1.*^300, 0.37, 7.3}; "
          "cf = Compile[{{v, _Real, 1}}, RiemannR[v]]; cs = Compile[{{x, _Real}}, RiemannR[x]]; "
          "RiemannR /@ xs == RiemannR[xs] == RiemannR[Developer`ToPackedArray[xs]] == "
          "Normal[RiemannR[NDArray[xs]]] == cf[xs] == (cs /@ xs)", "True");
    check("Developer`PackedArrayQ[RiemannR[Developer`ToPackedArray[{0.5, 2., 10.}]]]", "True");
    check("CompileDiagnostics[{{v, _Real, 1}}, RiemannR[v]][[1]]", "\"Compiled\" -> True");
    check("CompileDiagnostics[{{x, _Real}}, RiemannR[x]][[1]]", "\"Compiled\" -> True");
    /* a declined element (x < 1/e, x < 0) still answers through the List path */
    check("RiemannR[Developer`ToPackedArray[{0.1, 2.}]] == {RiemannR[0.1], RiemannR[2.]}", "True");
    check("Compile[{{x, _Real}}, RiemannR[x]][-1.] == RiemannR[-1.]", "True");
}

static void test_attributes(void) {
    check("Attributes[RiemannR]", "{Listable, NumericFunction, Protected}");
}

/* SIAM 100-Digit Challenge D4: the largest zero of R, x* = e^t*, needs the
 * series at ~49000 bits (terms to 10^14820 cancelling to 10^-13).  R changes
 * sign across t* = -34142.12818460646499263004606798354053561706...: one
 * evaluation on each side, 10^-30 away, at 40 digits.  (The full FindRoot at
 * WorkingPrecision -> 60, 50 correct digits, is benchmark 179's d4.wp50-builtin
 * row; it takes ~17 s, too long for a unit suite.) */
static void test_d4_root(void) {
    check("tz = -34142.12818460646499263004606798354053561706110688218996483979`70; "
          "{Sign[RiemannR[Exp[N[tz - 10^-30, 40]]]], Sign[RiemannR[Exp[N[tz + 10^-30, 40]]]]}",
          "{-1, 1}");   /* R > 0 exactly for x > x* */
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_exact);
    TEST(test_machine_real);
    TEST(test_complex);
    TEST(test_mpfr);
    TEST(test_derivative);
    TEST(test_surfaces);
    TEST(test_attributes);
    TEST(test_d4_root);
    printf("All RiemannR tests passed!\n");
    return 0;
}

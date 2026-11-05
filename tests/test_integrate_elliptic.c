/* test_integrate_elliptic.c
 *
 * The elliptic definite-integration stage Integrate`Elliptic
 * (src/calculus/integrate_elliptic.c): Int_a^b R(x, Sqrt[P(x)]) dx with P a
 * squarefree cubic or quartic, in closed form with EllipticF / EllipticE /
 * EllipticPi / EllipticK and elementary terms.
 *
 * The seven worked examples of Roach (2008) carry 30-digit reference values
 * computed independently (50-digit tanh-sinh quadrature in mpmath, with the
 * singular endpoints and interior split points made explicit); the root-
 * configuration sweep references were computed the same way.  NIntegrate is
 * not used as a reference: at WorkingPrecision 40 it returns 5.0111 for
 * Int_-1^3 x^2/Sqrt[(x+1)(3-x)(x^2-2x+2)] = 2 Sqrt[5] EllipticE[4/5] = 5.2704.
 * Every case must come back in CLOSED form (no Integrate left).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* The stage returns a closed form (no Integrate`Elliptic / Integrate head
 * left) whose 30-digit value agrees with `ref` to 10^-20. */
static void check_value(const char* f, const char* spec, const char* ref) {
    char buf[2048];
    snprintf(buf, sizeof buf,
             "With[{r = Integrate`Elliptic[%s, %s]}, "
             "FreeQ[r, Integrate`Elliptic | Integrate] && "
             "Abs[N[r - (%s), 40]] < 10^-20 (1 + Abs[N[%s, 40]])]",
             f, spec, ref, ref);
    check_eq(buf, "True");
}

static void check_declines(const char* f, const char* spec) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Head[Integrate`Elliptic[%s, %s]]", f, spec);
    check_eq(buf, "Integrate`Elliptic");
}

/* Roach 2008 */
static void test_roach(void) {
    check_value("1/((z - 5) Sqrt[(z^2 + 1) (z^2 - 4)])", "{z, -Infinity, -3}",
                "-0.030157804163465903564055244573393`30");
    check_value("1/((z + 14) Sqrt[(z - 1) ((z - 13)^2 + 81)])", "{z, 2, Infinity}",
                "0.035522982787892870969054399101842`30");
    check_value("1/((18 - z) Sqrt[(1 - z) ((16 - z)^2 + 64)])", "{z, -Infinity, 0}",
                "0.016074534639182151420953829232576`30");
    check_value("Sqrt[(z^2 - 1)/(9 - z^2)]", "{z, 1, 2}",
                "0.42935345715334728773589768427085`30");
    check_value("1/((z + 5) Sqrt[(z - 1) (5 - z) ((z - 1)^2 + 9)])", "{z, 4, 5}",
                "0.023040613456465724040052272451206`30");
    check_value("Sqrt[z/((z - 1) (z - 2))]", "{z, 2, 3}",
                "2.6696301067066588095129015821`28");
    check_value("1/((6 - z) Sqrt[(5 - z) (z - 2) (z - 3) z^3])", "{z, 3, 4}",
                "0.085529431421027252717085752561018`30");
    /* the paper's closed form, Sqrt[6] + 2^(3/2) (F - E) at phi = Pi/4, m = 1/2 */
    check_value("Sqrt[z/((z - 1) (z - 2))]", "{z, 2, 3}",
                "Sqrt[6] + 2^(3/2) EllipticF[Pi/4, 1/2] - 2^(3/2) EllipticE[Pi/4, 1/2]");
    /* the structure: incomplete F, E, Pi with exact arguments */
    check_eq("FreeQ[Integrate`Elliptic[1/((6 - z) Sqrt[(5 - z) (z - 2) (z - 3) z^3]), {z, 3, 4}], "
             "Sqrt[(5 - z) (z - 2) (z - 3) z^3]] && !FreeQ[Integrate`Elliptic["
             "1/((6 - z) Sqrt[(5 - z) (z - 2) (z - 3) z^3]), {z, 3, 4}], EllipticPi]", "True");
}

/* Root configurations, with and without poles, finite and infinite limits */
static void test_configurations(void) {
    /* four real roots */
    check_value("1/Sqrt[(x + 1) x (x - 2) (x - 5)]", "{x, 1/2, 3/2}",
                "0.37306480336959975931453446856132`30");
    check_value("x^2/Sqrt[(x + 1) x (x - 2) (x - 5)]", "{x, 0, 2}",
                "1.5951791086804783601035359613547`30");
    check_value("(x^2 + 2 x - 1)/((x + 3)^2 Sqrt[(x + 1) x (x - 2) (x - 5)])", "{x, 1/2, 3/2}",
                "0.044575238637510329792541361527227`30");
    check_value("1/((x + 3) Sqrt[(x + 1) x (x - 2) (x - 5)])", "{x, 5, Infinity}",
                "0.040963337198107672964181630263862`30");
    check_value("1/((x - 3) Sqrt[(x + 1) x (x - 2) (x - 5)])", "{x, -Infinity, -1}",
                "-0.1309998484893109777663555770166`30");
    /* two real roots and a complex pair */
    check_value("(x + 2)/((x + 3)^2 Sqrt[(x + 1) (3 - x) (x^2 - 2 x + 2)])", "{x, -1/2, 5/2}",
                "0.25137026374841169928150165837103`30");
    check_value("x^2/Sqrt[(x + 1) (3 - x) (x^2 - 2 x + 2)]", "{x, -1, 3}",
                "5.2703671631912601846706579443711`30");
    /* two complex pairs (x^4 + 1: irreducible over Q, roots in radicals) */
    check_value("1/((x + 3) Sqrt[(x^2 + 1) (x^2 + 2 x + 5)])", "{x, -2, 7}",
                "0.53224667618730211002440917044603`30");
    check_value("1/Sqrt[(x^2 + 1) (x^2 + 2 x + 5)]", "{x, -Infinity, Infinity}",
                "2.0992762520697304203200007956436`30");
    check_value("(x^2 + 3)/((x - 1)^2 Sqrt[(x^2 + 1) (x^2 + 2 x + 5)])", "{x, 2, Infinity}",
                "0.87125962533680626994219213793805`30");
    check_value("1/Sqrt[x^4 + 1]", "{x, 0, 1}",
                "0.92703733865068595921692517359763`30");
    check_value("(1 + x^3)/((2 + x)^2 Sqrt[x^4 + 1])", "{x, -1, 3}",
                "0.75848329238020567301613385350555`30");
    /* cubics: three real roots, one real root, a Cardano root */
    check_value("1/((x + 1/2) Sqrt[x (x + 1) (x - 2)])", "{x, 2, Infinity}",
                "0.40247161339850652109103157513306`30");
    check_value("Sqrt[x (x + 1) (2 - x)]", "{x, 0, 2}",
                "2.2033457318247437718068935054437`30");
    check_value("1/((x - 3)^2 Sqrt[-x (x + 1) (x - 2)])", "{x, -1, 0}",
                "(-0.17103117762932451390913584397047`30) I");
    check_value("x/Sqrt[(x - 2) (x^2 + 1)]", "{x, 2, 6}",
                "3.775588485625273887155485858947`30");
    check_value("1/((x + 1) Sqrt[(x - 2) (x^2 + 1)])", "{x, 3, Infinity}",
                "0.13844638307668577631434429100896`30");
    check_value("1/Sqrt[x^3 - 2]", "{x, 2, Infinity}",
                "1.4424544993425006729814685499474`30");
    /* second kind; negative radicand (an I factor); rational radicand */
    check_value("Sqrt[(1 - x^2) (4 - x^2)]", "{x, 0, 1}",
                "1.5200399881728984324569102410382`30");
    check_value("1/Sqrt[-(x^2 + 1) (x^2 + 4)]", "{x, 0, 1}",
                "(-0.42561187453559270452847425168529`30) I");
    check_value("Sqrt[(x + 1)/(x (x - 2))]", "{x, 3, 4}",
                "0.94089643997179028424783196839771`30");
    /* a rational part beside the radical one; a reversed interval;
     * a radical in the denominator; a half-line ending at a root */
    check_value("1/x + 1/Sqrt[x (x + 1) (x - 2)]", "{x, 3, 4}",
                "0.49911316113262529258765859357438`30");
    check_value("1/Sqrt[x (x + 1) (x - 2)]", "{x, 4, 3}",
                "-0.21143108868084436514843958758055`30");
    check_value("(1 + Sqrt[x (x^2 + 3)])/(2 x + Sqrt[x (x^2 + 3)])", "{x, 3/2, 2}",
                "0.31567309883448934623645344311946`30");
    check_value("1/Sqrt[x (x + 1) (x - 2)]", "{x, 2, Infinity}",
                "2.0021547609122124721593191722035`30");
    /* both limits at roots: complete integrals */
    check_value("1/Sqrt[(1 - x^2) (4 - x^2)]", "{x, -1, 1}", "EllipticK[1/4]");
}

static void test_declines(void) {
    check_declines("1/((x - 1) Sqrt[x (x + 1) (x - 2)])", "{x, 3/2, 3}");   /* root inside */
    check_declines("1/((x - 3) Sqrt[x (x + 1) (x - 2)])", "{x, 2, 4}");     /* pole inside */
    check_declines("1/Sqrt[1 - x^2]", "{x, 0, 1/2}");                       /* not elliptic */
    check_declines("1/Sqrt[x^5 + 1]", "{x, 0, 1}");                         /* hyperelliptic */
    check_declines("1/Sqrt[x^4 + a]", "{x, 0, 1}");                         /* symbolic */
    check_declines("1/(Sqrt[x^3 + 1] Sqrt[x^3 + 2])", "{x, 0, 1}");         /* two kernels */
    check_declines("x^2/Sqrt[x (x + 1) (x - 2)]", "{x, 3, Infinity}");      /* divergent */
    /* the internal probes leave no message behind */
    check_eq("Check[Integrate`Elliptic[1/Sqrt[x (x + 1) (x - 2)], {x, 2, Infinity}]; ok, failed]", "ok");
}

/* Roach 2008 #29 through the DEFAULT cascade.  Newton-Leibniz runs first, and
 * its indefinite search used to spend 13-16 s (55-90 s under load) in
 * DerivativeDivides' Eliminate/Solve pass: each shift kernel z - 1, z - 2,
 * 1/(z - 1), 1/(z - 2) re-expressed the same elliptic integrand in u and paid
 * a full cascade descent per Solve branch.  The elliptic-kernel gate sends it
 * down the direct-quotient pass only, as for a polynomial cubic radicand. */
static void test_default_cascade(void) {
    check_eq("With[{t = AbsoluteTiming[Integrate[Sqrt[z/((z - 1) (z - 2))], {z, 2, 3}]]}, "
             "{t[[1]] < 12, FreeQ[t[[2]], Integrate], "
             "Abs[N[t[[2]] - 2.6696301067066588095129015821`28, 30]] < 10^-20}]",
             "{True, True, True}");
    /* the indefinite form declines (non-elementary) promptly */
    check_eq("With[{t = AbsoluteTiming[Quiet[Integrate[Sqrt[z/((z - 1) (z - 2))], z]]]}, "
             "{t[[1]] < 8, Head[t[[2]]]}]", "{True, Integrate}");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_default_cascade);
    TEST(test_roach);
    TEST(test_configurations);
    TEST(test_declines);

    printf("All Integrate`Elliptic tests passed!\n");
    return 0;
}

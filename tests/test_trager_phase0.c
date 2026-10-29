/* test_trager_phase0.c
 *
 * Unit tests for the front end of the Trager integrator (src/internal/trager):
 *   Common.m   -- the failure protocol (throwNonElementary / throwDecline /
 *                 throwInternal, catchTrager, tragerFailQ, failToList) and the
 *                 predicates on Q (rationalNumberQ, rationalPolyQ,
 *                 rationalFunctionQ, rationalInXYQ, rdeg, tragerOption);
 *   Validate.m -- radical detection and the canonical form
 *                 f = Sum a_i(x) yh^i (radicalPowers, detectRadical with every
 *                 decline reason, substituteRadical, reduceYDegree,
 *                 powerBasisCoords, powerBasisInverse, canonicalForm,
 *                 radicalExpression);
 *   Reduce.m   -- the per-component curve (componentCurve in its three
 *                 regimes: d = 1 exact, d > 1 with a new radical, everything
 *                 absorbed m = 1; reducedCurveQ);
 *   Genus.m    -- curveGenus over Q and modulo p.
 *
 * Why these matter: canonicalForm is the only place where the integrand is
 * turned into algebra, so a wrong coordinate is a wrong integrand and hence a
 * wrong (yet self-consistently "verified") answer; the genus is the good-
 * reduction test of TorsionModP.m, so a wrong genus is a false certificate.
 *
 * How expectations are derived:
 *   - canonicalForm is checked by RECONSTRUCTION (tRecon): f is re-expressed
 *     in yh with substituteRadical and compared with Sum a_i yh^i modulo
 *     yh^n - g; and numerically (tNumCF) at a point where every radicand is
 *     positive, using radicalExpression for yh -- plus hand-computed
 *     coordinates for small cases;
 *   - componentCurve values are worked by hand from i e_j = n q_j + r_j and
 *     re-checked by the identity (yh^i)^n = s^n Y^n (exact case);
 *   - a table of 31 genera computed BY HAND from Riemann-Hurwitz
 *     (2g - 2 = -2m + Sum deg p (m - gcd(m, e_p)) + (m - gcd(m, deg G))), with
 *     the hyperelliptic (genus g for degree 2g+1 / 2g+2), Fermat
 *     ((n-1)(n-2)/2) and superelliptic-with-repeated-factor families, and
 *     reductions mod p including two of bad reduction (the genus drops).
 *
 * Symbols made by Unique (Y123) are never asserted in printed form.  Top-level
 * names carry a `t` prefix so they cannot collide with Trager`Private`
 * symbols, which are on $ContextPath.
 */
#include "trager_test_utils.h"

static void setup_phase0(void) {
    char* s = tr_eval(
        /* f equals Sum a_i yh^i modulo yh^n - g (exact). */
        "tRecon[tf_] := Module[{tcf = canonicalForm[tf, x], tn, tY, tfy, tdf},"
        "  tn = tcf[\"n\"]; tY = tcf[\"y\"];"
        "  If[tn === 1, Return[Together[tf - tcf[\"coeffs\"][[1]]] === 0]];"
        "  tfy = substituteRadical[tf, x, tcf[\"b\"], tn, tyy] /. tyy -> tY/tcf[\"q\"];"
        "  tdf = Together[tfy - tcf[\"coeffs\"] . tY^Range[0, tn - 1]];"
        "  reduceYDegree[Numerator[tdf], tY, tn, tcf[\"g\"]] === 0];"
        /* f equals Sum a_i (q b^(1/n))^i numerically at x = x0 (radicands > 0). */
        "tNumCF[tf_, tx0_] := Module[{tcf = canonicalForm[tf, x], tR},"
        "  tR = radicalExpression[tcf] /. x -> tx0;"
        "  Abs[N[(tf /. x -> tx0) - (tcf[\"coeffs\"] /. x -> tx0) . tR^Range[0, tcf[\"n\"] - 1], 40]] < 10^-30];"
        /* g = q^(n-1) Numerator[b] and the radical identity yh^n = g numerically */
        "tGCheck[tf_] := Module[{tcf = canonicalForm[tf, x]},"
        "  Expand[tcf[\"g\"] - tcf[\"q\"]^(tcf[\"n\"] - 1) Numerator[Together[tcf[\"b\"]]]] === 0 &&"
        "  Abs[N[(radicalExpression[tcf]^tcf[\"n\"] - tcf[\"g\"]) /. x -> 17/5, 40]] < 10^-30];"
        "True");
    free(s);
}

/* ---- Common.m: the failure protocol -------------------------------------- */
static void test_failure_protocol(void) {
    TR_TRUE("catchTrager[throwDecline[\"foo\"]] === tragerFail[\"Decline\", \"foo\", <||>]");
    TR_TRUE("catchTrager[throwNonElementary[\"torsion\", <|\"order\" -> 0|>]] ==="
            " tragerFail[\"NonElementary\", \"torsion\", <|\"order\" -> 0|>]");
    TR_TRUE("catchTrager[throwInternal[\"inv\", <|\"a\" -> 1|>]] === tragerFail[\"Internal\", \"inv\", <|\"a\" -> 1|>]");
    TR_TRUE("catchTrager[throwInternal[\"inv\"]] === tragerFail[\"Internal\", \"inv\", <||>]");
    /* a value passes through untouched */
    TR_EQ("catchTrager[1 + 1]", "2");
    TR_TRUE("catchTrager[{1, \"a\"}] === {1, \"a\"}");
    /* HoldFirst: the throw stops the rest of the compound expression */
    TR_TRUE("tSide = 0; catchTrager[throwDecline[\"a\"]; tSide = 1]; tSide === 0");
    TR_TRUE("tSide = 0; catchTrager[tSide = 1; 5] === 5 && tSide === 1");
    /* a throw deep inside a function is caught */
    TR_TRUE("tThrower[k_] := If[k > 3, throwDecline[\"deep\", <|\"k\" -> k|>], tThrower[k + 1]];"
            " catchTrager[tThrower[0]] === tragerFail[\"Decline\", \"deep\", <|\"k\" -> 4|>]");
    /* a throw with another tag is NOT intercepted */
    TR_EQ("Catch[catchTrager[Throw[5, \"Other\"]], \"Other\"]", "5");
    /* nested: the innermost catchTrager catches */
    TR_TRUE("catchTrager[{catchTrager[throwDecline[\"in\"]], 2}] === {tragerFail[\"Decline\", \"in\", <||>], 2}");
    /* a non-string reason does not match the definition (stays unevaluated) */
    TR_TRUE("catchTrager[throwDecline[3]] === throwDecline[3]");

    TR_TRUE("tragerFailQ[tragerFail[\"Decline\", \"x\", <||>]]");
    TR_TRUE("tragerFailQ[catchTrager[throwNonElementary[\"y\"]]]");
    TR_FALSE("tragerFailQ[{\"failed\", \"x\"}]");
    TR_FALSE("tragerFailQ[3]");
    TR_FALSE("tragerFailQ[catchTrager[7]]");

    TR_TRUE("failToList[tragerFail[\"NonElementary\", \"r\", <||>]] === {\"not elementary\", \"r\"}");
    TR_TRUE("failToList[tragerFail[\"Decline\", \"r\", <|\"z\" -> 1|>]] === {\"failed\", \"r\"}");
    TR_TRUE("failToList[tragerFail[\"Internal\", \"inv\", <||>]] === {\"failed\", \"internal: inv\"}");
    TR_TRUE("failToList[catchTrager[throwInternal[\"kInvZero\"]]] === {\"failed\", \"internal: kInvZero\"}");
}

/* ---- Common.m: predicates on Q ------------------------------------------- */
static void test_predicates(void) {
    TR_TRUE("rationalNumberQ[3]");
    TR_TRUE("rationalNumberQ[-2/5]");
    TR_TRUE("rationalNumberQ[0]");
    TR_TRUE("rationalNumberQ[10^40 + 1]");
    TR_FALSE("rationalNumberQ[1.5]");
    TR_FALSE("rationalNumberQ[Sqrt[2]]");
    TR_FALSE("rationalNumberQ[x]");
    TR_FALSE("rationalNumberQ[I]");
    TR_FALSE("rationalNumberQ[1/2 + I]");

    TR_TRUE("rationalPolyQ[x^2/3 + 1, x]");
    TR_TRUE("rationalPolyQ[(x + 1)^3, x]");
    TR_TRUE("rationalPolyQ[0, x]");
    TR_TRUE("rationalPolyQ[7/2, x]");
    TR_FALSE("rationalPolyQ[Sqrt[2] x, x]");
    TR_FALSE("rationalPolyQ[x + a, x]");
    TR_FALSE("rationalPolyQ[1/x, x]");
    TR_FALSE("rationalPolyQ[1.5 x, x]");
    TR_FALSE("rationalPolyQ[Sqrt[x], x]");

    TR_TRUE("rationalFunctionQ[1/(x + 1), x]");
    TR_TRUE("rationalFunctionQ[(x + 1/2)/(3 x), x]");
    TR_TRUE("rationalFunctionQ[x^2 + 1/x, x]");
    TR_TRUE("rationalFunctionQ[5, x]");
    TR_FALSE("rationalFunctionQ[Sqrt[x], x]");
    TR_FALSE("rationalFunctionQ[a/x, x]");
    TR_FALSE("rationalFunctionQ[1/(x + I), x]");

    TR_TRUE("rationalInXYQ[(x + y^2/3)/(1 + x y), x, y]");
    TR_TRUE("rationalInXYQ[y^5 + x/y, x, y]");
    TR_TRUE("rationalInXYQ[2, x, y]");
    TR_FALSE("rationalInXYQ[Sqrt[x] y, x, y]");
    TR_FALSE("rationalInXYQ[2.0 x + y, x, y]");
    TR_FALSE("rationalInXYQ[I y, x, y]");
    TR_FALSE("rationalInXYQ[x + y + a, x, y]");

    TR_EQ("rdeg[0, x]", "-Infinity");
    TR_EQ("rdeg[(x^3 + 1)/(x + 1)^5, x]", "-2");
    TR_EQ("rdeg[5, x]", "0");
    TR_EQ("rdeg[x^2/(x - 1), x]", "1");
    TR_EQ("rdeg[(x^2 - 1)/(x - 1), x]", "1");
    TR_EQ("rdeg[1/x - 1/x, x]", "-Infinity");
    TR_EQ("rdeg[x^4 + 1/x, x]", "4");

    TR_EQ("{tragerOption[\"ShiftBound\"], tragerOption[\"HermiteMaxIter\"], tragerOption[\"MaxFieldDegree\"],"
          " tragerOption[\"MaxTorsionSearch\"], tragerOption[\"MaxPrimes\"], tragerOption[\"KauersMaxPower\"]}",
          "{64, 200, 32, 400, 200, 16}");
    TR_TRUE("AllTrue[{\"ShiftBound\", \"HermiteMaxIter\", \"MaxFieldDegree\", \"MaxTorsionSearch\","
            " \"MaxPrimes\", \"KauersMaxPower\"}, IntegerQ[tragerOption[#]] && tragerOption[#] > 0 &]");
}

/* ---- Validate.m: radical detection --------------------------------------- */
#define DECLINE_REASON(f, reason) \
    TR_TRUE("Module[{e = catchTrager[detectRadical[" f ", x]]}," \
            " MatchQ[e, tragerFail[\"Decline\", _, _]] && e[[2]] === \"" reason "\"]")

static void test_detect(void) {
    TR_TRUE("Sort[radicalPowers[x Sqrt[x + 1] + 1/(x + 1)^(1/3) + 2^(1/2), x]] ==="
            " Sort[{(1 + x)^(-1/3), Sqrt[1 + x]}]");
    TR_EQ("radicalPowers[x^2 + 1/x, x]", "{}");
    /* duplicates are removed */
    TR_EQ("Length[radicalPowers[Sqrt[x] + 1/(1 + Sqrt[x]), x]]", "1");

    TR_TRUE("detectRadical[Sqrt[x^2 + 1]/x, x] === <|\"b\" -> 1 + x^2, \"n\" -> 2|>");
    TR_EQ("detectRadical[x^2, x]", "None");
    TR_EQ("detectRadical[(x + 1)/(x - 3), x]", "None");
    /* b^(1/2) and b^(1/3) are powers of b^(1/6) */
    TR_TRUE("detectRadical[Sqrt[x] + x^(1/3), x] === <|\"b\" -> x, \"n\" -> 6|>");
    TR_TRUE("detectRadical[(x + 1)^(1/3) + (x + 1)^(-1/2), x] === <|\"b\" -> 1 + x, \"n\" -> 6|>");
    /* the reciprocal radicand 1/b is merged with b */
    TR_TRUE("detectRadical[Sqrt[1/(x + 1)] + Sqrt[x + 1], x] === <|\"b\" -> 1/(1 + x), \"n\" -> 2|>");
    TR_TRUE("detectRadical[((x - 1)/(x + 1))^(1/3) + ((x + 1)/(x - 1))^(1/3), x][\"n\"] === 3");
    TR_TRUE("detectRadical[Sqrt[(x + 1)/(x - 1)] + 1/x, x] === <|\"b\" -> (1 + x)/(-1 + x), \"n\" -> 2|>");
    TR_TRUE("detectRadical[x^(5/2), x] === <|\"b\" -> x, \"n\" -> 2|>");
    TR_TRUE("detectRadical[(x^2 + 1)^(-3/4), x] === <|\"b\" -> 1 + x^2, \"n\" -> 4|>");
    TR_TRUE("detectRadical[Sqrt[x]/3 + 1/2, x] === <|\"b\" -> x, \"n\" -> 2|>");
    TR_TRUE("detectRadical[x^(2/5) + x^(3/10), x][\"n\"] === 10");

    DECLINE_REASON("Sin[x] Sqrt[x]", "transcendental function in the integrand");
    DECLINE_REASON("Exp[x] Sqrt[x]", "transcendental function in the integrand");
    DECLINE_REASON("E^x Sqrt[x]", "transcendental function in the integrand");
    DECLINE_REASON("Sqrt[x] Log[x]", "transcendental function in the integrand");
    DECLINE_REASON("ArcTan[x] + Sqrt[x]", "transcendental function in the integrand");
    DECLINE_REASON("Cosh[x]/x", "transcendental function in the integrand");
    DECLINE_REASON("1.5 Sqrt[x]", "coefficients outside Q");
    DECLINE_REASON("I Sqrt[x]", "coefficients outside Q");
    DECLINE_REASON("Sqrt[x + I]", "coefficients outside Q");
    DECLINE_REASON("Sqrt[2] Sqrt[x]", "algebraic constant in the integrand (coefficients must be rational)");
    DECLINE_REASON("x + 3^(1/3)", "algebraic constant in the integrand (coefficients must be rational)");
    DECLINE_REASON("Sqrt[1 + Sqrt[x]]", "nested radicals");
    DECLINE_REASON("(x + x^(1/3))^(1/2)", "nested radicals");
    DECLINE_REASON("Sqrt[x] + Sqrt[x + 1]", "more than one radicand (a radical compositum)");
    DECLINE_REASON("x^(1/3) + (x^2 + 1)^(1/3)", "more than one radicand (a radical compositum)");
    DECLINE_REASON("Sqrt[x + a]", "radicand is not a rational function over Q");
    /* the compositum decline carries both radicands */
    TR_TRUE("Sort[catchTrager[detectRadical[Sqrt[x] + Sqrt[x + 1], x]][[3]][\"Radicands\"]] === Sort[{x, 1 + x}]");
}

/* ---- Validate.m: substitution and the power basis ------------------------ */
static void test_power_basis(void) {
    /* x^(-3/2) = y^-3, Sqrt[1/x] = y^-1 (b = x, n = 2) */
    TR_TRUE("substituteRadical[Sqrt[x] + x^(-3/2), x, x, 2, tyy] === tyy + tyy^-3");
    TR_TRUE("Together[substituteRadical[Sqrt[x] + Sqrt[1/x], x, x, 2, tyy] - (tyy + 1/tyy)] === 0");
    /* b^(1/2), b^(1/3) with n = 6 -> y^3, y^2 */
    TR_TRUE("substituteRadical[Sqrt[x + 1] + (x + 1)^(1/3), x, 1 + x, 6, tyy] === tyy^2 + tyy^3");
    /* a reciprocal radicand gives the negated power */
    TR_TRUE("substituteRadical[Sqrt[1/(x + 1)] x, x, 1 + x, 2, tyy] === x/tyy");
    /* a constant radical is not touched */
    TR_TRUE("substituteRadical[x^2 + 1, x, x, 2, tyy] === x^2 + 1");

    TR_TRUE("reduceYDegree[tyy^5 + tyy, tyy, 2, x] === tyy + x^2 tyy");
    TR_TRUE("reduceYDegree[tyy + 3, tyy, 2, x] === 3 + tyy");
    TR_TRUE("reduceYDegree[tyy^3, tyy, 3, x^2 + 1] === 1 + x^2");
    TR_TRUE("reduceYDegree[(tyy + 1)^4, tyy, 3, x] === Expand[1 + 4 x + (4 + x) tyy + 6 tyy^2]");
    TR_EQ("reduceYDegree[0, tyy, 4, x]", "0");

    TR_TRUE("powerBasisCoords[tyy^4, tyy, 3, x + 1] === {0, 1 + x, 0}");
    TR_EQ("powerBasisCoords[0, tyy, 3, x]", "{0, 0, 0}");
    TR_EQ("powerBasisCoords[3, tyy, 2, x]", "{3, 0}");
    /* (y + 1)^2 = y^2 + 2y + 1 = (x + 1) + 2 y modulo y^2 - x */
    TR_TRUE("powerBasisCoords[(tyy + 1)^2, tyy, 2, x] === {1 + x, 2}");
    TR_TRUE("powerBasisCoords[tyy^7, tyy, 4, x^2] === {0, 0, 0, x^2}");

    /* 1/(1 + y) = (1 - y)/(1 - x) on y^2 = x */
    TR_TRUE("Together[powerBasisInverse[1 + tyy, tyy, 2, x, x] - {1/(1 - x), -1/(1 - x)}] === {0, 0}");
    /* 1/y = y^2/x on y^3 = x */
    TR_TRUE("powerBasisInverse[tyy, tyy, 3, x, x] === {0, 0, 1/x}");
    /* inverse property: q * q^-1 = 1 in Q(x)[y]/(y^n - g) */
    TR_TRUE("Module[{tq = 1 + x tyy + tyy^3, tinv}, tinv = powerBasisInverse[tq, tyy, 4, x^3 + 2, x];"
            " Together[powerBasisCoords[Expand[tq tinv . tyy^Range[0, 3]], tyy, 4, x^3 + 2]] === {1, 0, 0, 0}]");
    TR_TRUE("Module[{tq = x - tyy^2, tinv}, tinv = powerBasisInverse[tq, tyy, 3, x^2 + 1, x];"
            " Together[powerBasisCoords[Expand[tq tinv . tyy^Range[0, 2]], tyy, 3, x^2 + 1]] === {1, 0, 0}]");
    /* a zero divisor (y^2 - x^2 is reducible) declines */
    TR_TRUE("catchTrager[powerBasisInverse[tyy - x, tyy, 2, x^2, x]] ==="
            " tragerFail[\"Decline\", \"the denominator is a zero divisor: the radical is reducible\", <||>]");
}

/* ---- Validate.m: canonical form ------------------------------------------- */
static void test_canonical_form(void) {
    /* hand values */
    TR_TRUE("Module[{tcf = canonicalForm[Sqrt[x^2 + 1]/x, x]}, {tcf[\"b\"], tcf[\"n\"], tcf[\"q\"], tcf[\"g\"], tcf[\"coeffs\"]} ==="
            " {1 + x^2, 2, 1, 1 + x^2, {0, 1/x}} && Head[tcf[\"y\"]] === Symbol && tcf[\"x\"] === x]");
    /* 1/Sqrt[x] = y/x */
    TR_TRUE("canonicalForm[1/Sqrt[x], x][\"coeffs\"] === {0, 1/x}");
    TR_TRUE("canonicalForm[x^(1/3) + x^(2/3) + 1, x][\"coeffs\"] === {1, 1, 1}");
    /* 1/(1 + y) = (1 - y)/(1 - x) */
    TR_TRUE("Together[canonicalForm[1/(1 + Sqrt[x]), x][\"coeffs\"] - {1/(1 - x), -1/(1 - x)}] === {0, 0}");
    /* b = x/(x+1): q = x + 1, g = x(x + 1), yh = (x+1) y, so y = yh/(x+1) */
    TR_TRUE("Module[{tcf = canonicalForm[Sqrt[x/(x + 1)], x]}, tcf[\"q\"] === 1 + x &&"
            " Expand[tcf[\"g\"]] === x + x^2 && Together[tcf[\"coeffs\"] - {0, 1/(1 + x)}] === {0, 0}]");
    /* (x^2+1)^(-3/4) = y^-3 = y/(x^2 + 1) */
    TR_TRUE("Together[canonicalForm[(x^2 + 1)^(-3/4), x][\"coeffs\"] - {0, 1/(1 + x^2), 0, 0}] === {0, 0, 0, 0}");
    /* Sqrt[x] + x^(1/3) = y^3 + y^2 with y = x^(1/6) */
    TR_TRUE("canonicalForm[Sqrt[x] + x^(1/3), x][\"coeffs\"] === {0, 0, 1, 1, 0, 0}");
    /* reciprocal radicand: b = 1/(1+x), q = 1+x, g = 1+x; Sqrt[1/(x+1)] = yh/(1+x),
       Sqrt[x+1] = yh, so f = yh (2 + x)/(1 + x) */
    TR_TRUE("Module[{tcf = canonicalForm[Sqrt[1/(x + 1)] + Sqrt[x + 1], x]}, Expand[tcf[\"g\"]] === 1 + x &&"
            " Together[tcf[\"coeffs\"] - {0, (2 + x)/(1 + x)}] === {0, 0}]");
    /* (x+1)^(5/2) = (1 + x)^2 y */
    TR_TRUE("Together[canonicalForm[(x + 1)^(5/2), x][\"coeffs\"] - {0, (1 + x)^2}] === {0, 0}");
    /* no radical: n = 1 */
    TR_TRUE("Module[{tcf = canonicalForm[x/(x^2 + 1), x]}, tcf[\"n\"] === 1 && tcf[\"y\"] === None &&"
            " tcf[\"coeffs\"] === {x/(1 + x^2)}]");
    TR_TRUE("canonicalForm[(x^2 - 1)/(x - 1), x][\"coeffs\"] === {1 + x}");

    /* reconstruction (exact) */
    TR_TRUE("tRecon[Sqrt[x^2 + 1]/x]");
    TR_TRUE("tRecon[1/(1 + (x/(x + 1))^(1/3))]");
    TR_TRUE("tRecon[(x + Sqrt[x^3 + 1])/(x - Sqrt[x^3 + 1])]");
    TR_TRUE("tRecon[1/(x^(1/3) + x^(1/2) + 1)]");
    TR_TRUE("tRecon[(x^2 + 1)^(2/3)/(x^3 + (x^2 + 1)^(1/3))]");
    TR_TRUE("tRecon[Sqrt[(x - 1)/(x + 2)]/(1 + x Sqrt[(x + 2)/(x - 1)])]");
    TR_TRUE("tRecon[x^2/(x^4 + 1)^(1/4) + (x^4 + 1)^(3/4)]");
    TR_TRUE("tRecon[1/(2 - (x^5 + x)^(1/5))]");
    TR_TRUE("tRecon[(3 x + 1)/(x^2 - 7)]");

    /* numeric agreement at a point with positive radicands */
    TR_TRUE("tNumCF[Sqrt[x^2 + 1]/x, 3/7]");
    TR_TRUE("tNumCF[1/(1 + (x/(x + 1))^(1/3)), 5/2]");
    TR_TRUE("tNumCF[(x + Sqrt[x^3 + 1])/(x - Sqrt[x^3 + 1]), 2/3]");
    TR_TRUE("tNumCF[1/(x^(1/3) + x^(1/2) + 1), 7/3]");
    TR_TRUE("tNumCF[Sqrt[(x - 1)/(x + 2)]/(1 + x Sqrt[(x + 2)/(x - 1)]), 4]");
    TR_TRUE("tNumCF[Sqrt[1/(x + 1)] + Sqrt[x + 1], 1/3]");
    TR_TRUE("tNumCF[1/(2 - (x^5 + x)^(1/5)), 1/2]");

    /* g = q^(n-1) num(b) and yh^n = g */
    TR_TRUE("tGCheck[Sqrt[x/(x + 1)]]");
    TR_TRUE("tGCheck[((x^2 + 1)/(x^3 - 2))^(1/4)]");
    TR_TRUE("tGCheck[1/(1 + (x/(x + 1))^(1/3))]");
    TR_TRUE("Module[{tcf = canonicalForm[((x^2 + 1)/(x^3 - 2))^(1/4), x]},"
            " Expand[tcf[\"g\"]] === Expand[(x^3 - 2)^3 (x^2 + 1)] && Exponent[tcf[\"g\"], x] === 11]");

    /* declines */
    TR_TRUE("catchTrager[canonicalForm[1/(Sqrt[x^2] - x), x]][[2]] === \"the denominator is a zero divisor: the radical is reducible\"");
    TR_TRUE("catchTrager[canonicalForm[x + a, x]][[2]] === \"integrand is not algebraic over Q(x)\"");
    TR_TRUE("catchTrager[canonicalForm[Sqrt[x] + a, x]][[2]] === \"integrand is not a rational function of x and the radical over Q\"");
    TR_TRUE("catchTrager[canonicalForm[Sqrt[x] + Sqrt[x + 1], x]][[2]] === \"more than one radicand (a radical compositum)\"");
    TR_TRUE("catchTrager[canonicalForm[Sin[x], x]][[1]] === \"Decline\"");

    /* radicalExpression = q b^(1/n) */
    TR_TRUE("radicalExpression[canonicalForm[Sqrt[x/(x + 1)], x]] === (1 + x) Sqrt[x/(1 + x)]");
    TR_TRUE("radicalExpression[canonicalForm[(x^3 + 1)^(1/3)/x, x]] === (1 + x^3)^(1/3)");
}

/* ---- Reduce.m: component curves ------------------------------------------ */
static void test_component_curve(void) {
    /* d = 1, exact */
    TR_TRUE("componentCurve[2, x^2 + 1, 1, x] === <|\"m\" -> 2, \"G\" -> 1 + x^2, \"s\" -> 1, \"exact\" -> True|>");
    /* constant factor kept in the exact case: 3 (x + 1) */
    TR_TRUE("componentCurve[2, 3 (x + 1), 1, x] === <|\"m\" -> 2, \"G\" -> 3 + 3 x, \"s\" -> 1, \"exact\" -> True|>");
    /* x^5 = x^3 x^2: s = x, G = x^2 */
    TR_TRUE("componentCurve[3, x^5, 1, x] === <|\"m\" -> 3, \"G\" -> x^2, \"s\" -> x, \"exact\" -> True|>");
    TR_TRUE("componentCurve[6, x^3 (x + 1)^2, 1, x] === <|\"m\" -> 6, \"G\" -> Expand[x^3 (x + 1)^2], \"s\" -> 1, \"exact\" -> True|>");
    /* i = 5 on y^6 = x^3 (x+1)^2: 15 = 6*2 + 3, 10 = 6 + 4; gcd(6, 3, 4) = 1:
       s = x^2 (x + 1), G = x^3 (x + 1)^4 */
    TR_TRUE("componentCurve[6, x^3 (x + 1)^2, 5, x] === <|\"m\" -> 6, \"G\" -> Expand[x^3 (x + 1)^4],"
            " \"s\" -> Expand[x^2 (x + 1)], \"exact\" -> True|>");
    /* exact case with a constant: i = 2 on y^3 = 2x -> (2x)^2 = 4 x^2 */
    TR_TRUE("componentCurve[3, 2 x, 2, x] === <|\"m\" -> 3, \"G\" -> 4 x^2, \"s\" -> 1, \"exact\" -> True|>");
    /* d > 1: i = 2 on y^6 = x^3 (x+1)^2: 6 = 6 + 0, 4: d = gcd(6, 4) = 2,
       m = 3, H = (x + 1)^2, s = x */
    TR_TRUE("componentCurve[6, x^3 (x + 1)^2, 2, x] === <|\"m\" -> 3, \"G\" -> Expand[(1 + x)^2], \"s\" -> x, \"exact\" -> False|>");
    /* i = 3: 9 = 6 + 3, 6 = 6 + 0: d = 3, m = 2, H = x, s = x (x + 1) */
    TR_TRUE("componentCurve[6, x^3 (x + 1)^2, 3, x] === <|\"m\" -> 2, \"G\" -> x, \"s\" -> x + x^2, \"exact\" -> False|>");
    /* the constant is absorbed by rho when d > 1 */
    TR_TRUE("componentCurve[6, 4 x^3 (x + 1)^2, 3, x] === componentCurve[6, x^3 (x + 1)^2, 3, x]");
    /* i = 2 on y^4 = x^2 + 1 (e = 1): 2 = 0*4 + 2, d = 2, m = 2, H = x^2 + 1 */
    TR_TRUE("componentCurve[4, x^2 + 1, 2, x] === <|\"m\" -> 2, \"G\" -> 1 + x^2, \"s\" -> 1, \"exact\" -> False|>");
    /* everything absorbed: i = 2 on y^4 = x^2 (x+1)^2 -> m = 1, G = 1, s = x (x + 1) */
    TR_TRUE("componentCurve[4, x^2 (x + 1)^2, 2, x] === <|\"m\" -> 1, \"G\" -> 1, \"s\" -> x + x^2, \"exact\" -> False|>");
    TR_TRUE("componentCurve[2, 5 x^2, 1, x] === <|\"m\" -> 1, \"G\" -> 1, \"s\" -> x, \"exact\" -> False|>");
    /* i = 0 is the rational component */
    TR_TRUE("componentCurve[3, x^2 + 1, 0, x][\"m\"] === 1");

    /* exact case identity: (yh^i)^n = s^n G, i.e. g^i = s^n G */
    TR_TRUE("AllTrue[{{6, x^3 (x + 1)^2, 1}, {6, x^3 (x + 1)^2, 5}, {3, x^5, 1}, {3, 2 x, 2}, {5, x^7 (x^2 + 1)^3, 3}},"
            " Module[{tc = componentCurve[#[[1]], #[[2]], #[[3]], x]}, tc[\"exact\"] &&"
            " Expand[#[[2]]^#[[3]] - tc[\"s\"]^#[[1]] tc[\"G\"]] === 0] &]");
    /* inexact case: g^i / s^n is a constant times G^d (d = n/m) */
    TR_TRUE("AllTrue[{{6, x^3 (x + 1)^2, 2}, {6, 4 x^3 (x + 1)^2, 3}, {4, x^2 + 1, 2}, {4, 3 x^2 (x + 1)^2, 2}},"
            " Module[{tc = componentCurve[#[[1]], #[[2]], #[[3]], x], tq}, tq = Together[#[[2]]^#[[3]]/(tc[\"s\"]^#[[1]] tc[\"G\"]^(#[[1]]/tc[\"m\"]))];"
            " !tc[\"exact\"] && FreeQ[tq, x]] &]");
    /* every component curve with m >= 2 is reduced */
    TR_TRUE("AllTrue[{{6, x^3 (x + 1)^2, 1}, {6, x^3 (x + 1)^2, 2}, {6, x^3 (x + 1)^2, 3}, {6, x^3 (x + 1)^2, 5},"
            " {4, x^2 + 1, 2}, {3, x^5, 1}, {5, x^7 (x^2 + 1)^3, 3}},"
            " Module[{tc = componentCurve[#[[1]], #[[2]], #[[3]], x]}, reducedCurveQ[tc[\"m\"], tc[\"G\"], x]] &]");

    TR_TRUE("reducedCurveQ[2, x^2 + 1, x]");
    TR_TRUE("reducedCurveQ[6, x^3 (x + 1)^2, x]");
    TR_TRUE("reducedCurveQ[3, 2 x^2, x]");
    TR_FALSE("reducedCurveQ[4, x^2 (x + 1)^2, x]");     /* gcd 2 */
    TR_FALSE("reducedCurveQ[2, x^3, x]");               /* multiplicity >= m */
    TR_FALSE("reducedCurveQ[1, x, x]");                 /* m < 2 */
    TR_FALSE("reducedCurveQ[2, 5, x]");                 /* constant G */
    TR_FALSE("reducedCurveQ[2, Sqrt[2] x + 1, x]");     /* not over Q */
    TR_FALSE("reducedCurveQ[2, x + a, x]");
    TR_FALSE("reducedCurveQ[6, x^2 (x + 1)^4, x]");     /* gcd(6, 2, 4) = 2 */
}

/* ---- Genus.m --------------------------------------------------------------- */
/* Each value computed by hand from 2g - 2 = -2m + Sum deg p (m - gcd(m, e))
   + (m - gcd(m, deg G)). */
static void test_genus(void) {
    /* hyperelliptic y^2 = squarefree of degree 2g+1 / 2g+2 */
    TR_EQ("curveGenus[2, x, x]", "0");
    TR_EQ("curveGenus[2, x^2 + 1, x]", "0");
    TR_EQ("curveGenus[2, x^2 - 2, x]", "0");
    TR_EQ("curveGenus[2, x^3 + 1, x]", "1");
    TR_EQ("curveGenus[2, 3 x^3 - 5, x]", "1");
    TR_EQ("curveGenus[2, x^4 + 1, x]", "1");
    TR_EQ("curveGenus[2, x^5 - x + 1, x]", "2");
    TR_EQ("curveGenus[2, x^6 + 1, x]", "2");
    TR_EQ("curveGenus[2, x^7 - 1, x]", "3");
    TR_EQ("curveGenus[2, x^8 - 1, x]", "3");
    TR_EQ("curveGenus[2, x^9 - x, x]", "4");
    /* Fermat y^n = x^n + 1: (n - 1)(n - 2)/2 */
    TR_EQ("curveGenus[3, x^3 + 1, x]", "1");
    TR_EQ("curveGenus[4, x^4 + 1, x]", "3");
    TR_EQ("curveGenus[5, x^5 + 1, x]", "6");
    TR_EQ("curveGenus[6, x^6 + 1, x]", "10");
    TR_EQ("curveGenus[7, x^7 + 1, x]", "15");
    /* superelliptic, simple factors */
    TR_EQ("curveGenus[3, x, x]", "0");
    TR_EQ("curveGenus[3, x^2 + 1, x]", "1");
    TR_EQ("curveGenus[3, x^4 + 1, x]", "3");
    TR_EQ("curveGenus[3, x (x - 1) (x + 1) (x - 2), x]", "3");
    TR_EQ("curveGenus[5, x^2 + 1, x]", "2");
    TR_EQ("curveGenus[6, x (x - 1), x]", "2");
    TR_EQ("curveGenus[7, x (x - 1), x]", "3");
    /* repeated factors */
    TR_EQ("curveGenus[3, x^2 (x + 1), x]", "0");
    TR_EQ("curveGenus[4, x^2 (x + 1), x]", "1");
    TR_EQ("curveGenus[6, x^3 (x + 1)^2, x]", "1");
    TR_EQ("curveGenus[4, x (x^2 + 1), x]", "3");
    TR_EQ("curveGenus[4, x^2 (x^2 + 1), x]", "1");
    TR_EQ("curveGenus[5, x^2 (x - 1)^3, x]", "0");
    TR_EQ("curveGenus[4, x^3, x]", "0");
    TR_EQ("curveGenus[3, x^2 (x^2 + 1)^2 (x - 3), x]", "3");
    /* the hyperelliptic family (x^k - 2 is Eisenstein-irreducible): g = Floor[(k - 1)/2] */
    TR_TRUE("AllTrue[Range[1, 12], curveGenus[2, x^# - 2, x] === Floor[(# - 1)/2] &]");
    TR_TRUE("AllTrue[Range[2, 9], curveGenus[#, x^# + 1, x] === (# - 1) (# - 2)/2 &]");

    /* modulo p (p not dividing m) */
    TR_EQ("curveGenus[2, x^5 - x + 1, x, 5]", "2");
    TR_EQ("curveGenus[2, x^3 - x, x, 3]", "1");
    TR_EQ("curveGenus[2, x^4 + 1, x, 3]", "1");
    TR_EQ("curveGenus[2, x^6 + 1, x, 5]", "2");
    /* rational coefficient reduced first: x^3 + 1/2 = x^3 + 4 (mod 7) */
    TR_EQ("curveGenus[2, x^3 + 1/2, x, 7]", "1");
    TR_EQ("curveGenus[3, x^4 + 1, x, 5]", "3");
    /* bad reduction: x^3 + 1 = (x + 1)^3 mod 3, x^4 + 1 = (x + 1)^4 mod 2 */
    TR_EQ("curveGenus[2, x^3 + 1, x, 3]", "0");
    TR_EQ("curveGenus[3, x^4 + 1, x, 2]", "0");
    /* x^2 + 1 = (x + 2)(x + 3) mod 5 vs irreducible over Q: same genus */
    TR_EQ("curveGenus[3, x^2 + 1, x, 5]", "1");
    /* the leading coefficient vanishing mod p drops the degree: 5 x^3 + x + 1 */
    TR_EQ("curveGenus[2, 5 x^3 + x + 1, x, 5]", "0");
    TR_EQ("curveGenus[2, 5 x^3 + x + 1, x]", "1");

    /* curveGenusFromFactors directly */
    TR_EQ("curveGenusFromFactors[2, {{x, 1}, {x - 1, 1}, {x + 1, 1}}, 3, x]", "1");
    TR_EQ("curveGenusFromFactors[4, {{1 + x^2, 1}, {x, 2}}, 4, x]", "1");
}

int main(void) {
    tr_init();
    setup_phase0();
    test_failure_protocol();
    test_predicates();
    test_detect();
    test_power_basis();
    test_canonical_form();
    test_component_curve();
    test_genus();
    printf("trager_phase0: %d assertions passed\n", tr_assertions);
    return 0;
}

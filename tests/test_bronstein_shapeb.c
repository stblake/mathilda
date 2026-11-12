/* test_bronstein_shapeb.c -- Shape B, a radical y^n = b(x) at the BOTTOM of
 * a Log/Exp tower (BRONSTEIN_PLAN.md "Shape B"):
 *
 *   - the front end: radical-argument kernels (shape "BCurve"), the genus-0
 *     pre-pass (Moebius / Euler substitution), the canonical presentation;
 *   - K0.m, arithmetic in k0 = k(Y) and k0[theta];
 *   - BaseAlgebraic.m, the component method: every routine directly, the
 *     closes (each verified by PossibleZeroQ), the certificates (each with its
 *     exact reason, each a theorem: the residue criterion over k0, a certified
 *     "no solution" of a Risch DE over k, or Trager's own), and a seeded
 *     property test (random D[F] integrands close).
 */
#include "bronstein_test_utils.h"

static void test_front_end(void) {
    BR_TRUE("brRadicalOverXQ[Log[x + Sqrt[1 + x^3]], x]");
    BR_FALSE("brRadicalOverXQ[Sqrt[Log[x]] + x, x]");
    /* one radical-argument kernel is admitted (the curve-monomial case) ... */
    BR_TRUE("Module[{r = brKernelsToSymbols[Log[Sqrt[x] + 1], x]}, r[[1]] === r[[2, 1]][\"sym\"] && TrueQ[r[[2, 1]][\"algebraic\"]]]");
    BR_TRUE("brTower[Log[x + Sqrt[1 + x^3]] x, x][\"shape\"] === \"BCurve\"");
    BR_TRUE("brTower[Log[x] Sqrt[1 + x^3], x][\"shape\"] === \"B\"");
    /* ... but not over or under another monomial, nor with a radical over the tower */
    BR_TRUE("MemberQ[{{\"Decline\", \"radical inside a logarithm or exponential\"},"
            " {\"Decline\", \"a further monomial over a kernel with a radical argument (Shape B, curve-monomial case, takes one kernel)\"}},"
            " brFailTag[brCatch[brKernelsToSymbols[Log[x] + Log[Sqrt[x] + 1], x]]]]");
    BR_THROWS("brKernelsToSymbols[Log[Sqrt[Log[x]] + 1], x]", "Decline", "radical inside a logarithm or exponential");
}

static void test_genus0(void) {
    BR_TRUE("brConicPoint[x^2 + 1, x] === {\"Infinity\", 1}");
    BR_TRUE("brConicPoint[2 x^2 + 1, x] === {0, 1}");
    BR_TRUE("brConicPoint[2 x^2 + 7, x] === {1, 3}");
    BR_TRUE("brConicPoint[3 x^2 + 2, x] === None");
    BR_TRUE("Together[brEulerFlatten[Sqrt[1 + x^2], x, 1 + x^2, u][\"x\"] - (u^2 - 1)/(2 u)] === 0");
    BR_TRUE("Together[brEulerFlatten[Sqrt[1 + 2 x^2], x, 1 + 2 x^2, u][\"x\"] + 2 u/(u^2 - 2)] === 0");
    BR_TRUE("brEulerFlatten[1, x, 1 + x^3, u] === None");
    /* the parametrisation is an identity: x(u(x)) = x and y(u(x)) = Sqrt[b] */
    BR_TRUE("Module[{r = brEulerFlatten[Sqrt[2 x^2 + 7], x, 2 x^2 + 7, u]},"
            " PossibleZeroQ[(r[\"x\"] /. r[\"back\"]) - x]]");
    BR_TRUE("brGenus0Flatten[Log[x + Sqrt[1 + x^2]], x, u][\"back\"] === (u -> x + Sqrt[1 + x^2])");
    BR_TRUE("brGenus0Flatten[Log[x] Sqrt[1 + x^3], x, u] === None");
    BR_TRUE("brGenus0Flatten[Sqrt[1 + x^2], x, u] === None");
    BR_TRUE("brCurveCanon[Log[(x^2 + x Sqrt[1 + x^2])/(x + Sqrt[1 + x^2])], x, 1 + x^2, 2] === Log[x]");
    BR_TRUE("brCurveCanon[1/(x + Sqrt[1 + x^2]), x, 1 + x^2, 2] === Sqrt[1 + x^2] - x");
    BR_TRUE("(brCanonDone[Log, x] /. brCanonDone[h_, a_] :> h[a]) === Log[x]");
    /* kernels containing the radical, decided after the substitution */
    BR_ANTIDERIV("D[Log[x] Log[x + Sqrt[x^2 + 2]], x]", "Bronstein");
    BR_ANTIDERIV("Log[x] Log[x + Sqrt[x^2 + 2]]", "Bronstein");
    BR_ANTIDERIV("D[x Log[1 + Sqrt[x + 3]], x]", "Bronstein");
    BR_TRUE("BronsteinIntegrate[D[x Log[1 + Sqrt[x + 3]], x], x] === x Log[1 + Sqrt[3 + x]]");
}

static void test_k0(void) {
    br_run("Bk = curveBasis[2, 1 + x^3, x]; Tk = brTower[Log[x] Sqrt[1 + x^3], x];"
           " Fk = brK0Field[Tk]; th = Tk[\"syms\"][[1]]; True");
    BR_TRUE("brK0Field[Tk] === fieldFn[Tk[\"syms\"]]");
    BR_TRUE("brK0Scalar[Bk, x] === {x, 0}");
    BR_TRUE("brK0YMultiple[Bk, x] === {0, x}");
    BR_TRUE("brK0FromY[fieldFn[{}], Bk, 2 + x Y, Y] === {2, x}");
    BR_TRUE("brK0ToY[Bk, {2, x}, Y] === 2 + x Y");
    BR_TRUE("brK0Times[Fk, Bk, {th, 1}, {th, -1}] === {th^2 - 1 - x^3, 0}");
    BR_TRUE("brK0Times[Fk, Bk, {th, 1}, brK0Inverse[Fk, Bk, {th, 1}]] === {1, 0}");
    BR_THROWS("brK0Inverse[Fk, Bk, {0, 0}]", "Internal", "brK0InverseZero");
    /* D(theta + Y) = 1/x + s Y */
    BR_TRUE("brK0D[Tk, Fk, Bk, {th, 1}] === {1/x, 3 x^2/(2 (1 + x^3))}");
    BR_TRUE("brK0Degree[{th^2 + x, th}, th] === 2");
    BR_TRUE("brK0Degree[{0, 0}, th] === -Infinity");
    BR_THROWS("brK0Degree[{1/th, 0}, th]", "Internal", "brK0NotPolynomial");
    BR_TRUE("brK0Coefficient[Fk, {th^2 + x th, 3 th}, th, 1] === {x, 3}");
    BR_TRUE("brK0Monic[Fk, Bk, {2 th, 2}, th] === {th, 1}");
    /* theta^2 - Y^2 = (theta - Y)(theta + Y): factors over k0, not over k */
    BR_TRUE("brK0PolyRem[Fk, Bk, {th^2 - 1 - x^3, 0}, {th, 1}, th] === {0, 0}");
    BR_TRUE("brK0PolyGCD[Fk, Bk, {th^2 - 1 - x^3, 0}, {th, -1}, th] === {th, -1}");
    BR_TRUE("brK0PolyGCD[Fk, Bk, {th, 0}, {th + 1, 0}, th] === {1, 0}");
    BR_TRUE("brK0PolyGCD[Fk, Bk, {0, 0}, {0, 0}, th] === {0, 0}");
    BR_TRUE("brK0Trace[brK0Field[Tk, {tt, tt^2 - 2, Sqrt[2]}], tt^2 x + tt] === 4 x");
    /* s a + t b = g: (theta - Y) and (theta + Y) are coprime in k0[theta] */
    BR_TRUE("Module[{e = brK0PolyExtGCD[Fk, Bk, {th, -1}, {th, 1}, th]},"
            " e[[1]] === {1, 0} && afZeroQ[Fk, afMinus[Fk, afPlus[Fk, brK0Times[Fk, Bk, e[[2]], {th, -1}],"
            "   brK0Times[Fk, Bk, e[[3]], {th, 1}]], {1, 0}]]]");
    BR_TRUE("brK0PolyExtGCD[Fk, Bk, {th^2 - 1 - x^3, 0}, {th, 1}, th][[1]] === {th, 1}");
}

static void test_component_routines(void) {
    br_run("Tc = brTower[Log[x] Sqrt[1 + x^3], x]; thc = Tc[\"syms\"][[1]]; cvc = brShapeBCurve[2, 1 + x^3, x];"
           " Te = brTower[E^x Sqrt[1 + x^3], x]; True");
    BR_TRUE("Together[cvc[\"s\"] - 3 x^2/(2 (1 + x^3))] === 0");
    BR_TRUE("cvc[\"B\"][\"d\"] === {1, 1} && cvc[\"m\"] === 2 && cvc[\"x\"] === x");
    BR_TRUE("brShapeBTop[Tc][\"theta\"] === thc && brShapeBTop[Tc][\"kind\"] === \"Log\"");
    BR_TRUE("brShapeBRest[Tc][\"syms\"] === {} && brShapeBRest[Tc][\"shape\"] === \"B\"");
    BR_TRUE("brShapeBResult[1, {}, 0, Y] === <|\"Alg\" -> 1, \"Logs\" -> {}, \"Extra\" -> 0, \"Y\" -> Y|>");
    BR_TRUE("brShapeBExpression[brShapeBJoin[brShapeBResult[x Y, {{1, x}}, 0, Y], brShapeBResult[0, {{2, x + Y}}, x, Y]]]"
            " === x + x Y + Log[x] + 2 Log[x + Y]");
    BR_THROWS("brShapeBJoin[brShapeBResult[0, {}, 0, Y], brShapeBResult[0, {}, 0, Z]]", "Internal", "brShapeBJoinY");
    /* the base: Trager on the curve */
    BR_TRUE("brShapeBExpression[brShapeBBase[cvc, 3 x^2/(2 (1 + x^3)), Y]] === Y");
    BR_THROWS("brShapeBBase[cvc, 1/(1 + x^3), Y]", "NonElementary",
              "Trager (algebraic base): the Hermite remainder is a nonzero holomorphic differential (no residues)");
    /* the remainder split: Laurent part and the proper normal part */
    BR_TRUE("brShapeBSplitRemainder[brShapeBTop[Te], (2 + brT1^2)/(brT1 (1 + brT1))] ==="
            " <|\"L\" -> {{-1, 2}, {0, 1}}, \"a\" -> -3, \"d\" -> 1 + brT1|>");
    BR_TRUE("brShapeBResidueRoots[z^3 - 2 z, z] === {{z^2 - 2, 2}}");
    BR_TRUE("brShapeBLaurent[brShapeBTop[Tc], cvc, {x thc, thc^2 + 1/thc}] === {{{-1, 1}, {2, 1}}, x thc}");
    /* Risch DEs over the base: (beta Y)' = beta' + s beta */
    BR_TRUE("brShapeBRDE[brShapeBRest[Tc], cvc[\"s\"], cvc[\"s\"], \"r\"] === 1");
    BR_THROWS("brShapeBRDE[brShapeBRest[Tc], cvc[\"s\"], 1, \"no rde\"]", "NonElementary", "no rde");
    /* the polynomial part and the component: D(theta Y) = (1/x + s theta) Y */
    BR_TRUE("brShapeBExpression[brShapeBPolyPart[brShapeBTop[Tc], cvc, {{1, cvc[\"s\"]}, {0, 1/x}}, Y]] === thc Y");
    BR_TRUE("brShapeBExpression[brShapeBComponent[Tc, cvc, 1/x + cvc[\"s\"] thc, Y]] === thc Y");
    /* the log part over k0: residues +-Sqrt[2], logands theta -+ Sqrt[2] Y */
    BR_TRUE("Module[{lp = brShapeBLogPart[brShapeBTop[Tc], cvc, 4 (1/x - cvc[\"s\"] thc), thc^2 - 2 (1 + x^3), Y]},"
            " afZeroQ[brK0Field[Tc], lp[\"Rem\"]] && Length[lp[\"Logs\"]] === 1 && lp[\"Logs\"][[1, 3]] === \"RootSum\"]");
    BR_THROWS("brShapeBLogPart[brShapeBTop[Tc], cvc, 1, thc^2 - 2 (1 + x^3), Y]", "NonElementary",
              "non-constant residues (resultant criterion over the algebraic base, [B] Prop 2.9)");
    BR_TRUE("Module[{f = D[Log[x] Sqrt[1 + x^3], x], T, r}, T = brTower[f, x]; r = brShapeB[T, f, x];"
            " r[\"verified\"] === False && PossibleZeroQ[D[r[\"F\"], x] - f]]");
    BR_THROWS("brShapeB[brTower[Log[x] Sqrt[1 + x^3], x], Log[x] Sqrt[1 + x^3], x]", "NonElementary",
              "no algebraic part: the coefficient Risch DE D b + s b = c has no solution over the base");
}

static const char* const closes[] = {
    "D[Sqrt[1 + x^2] Log[x], x]",
    "D[E^x Sqrt[1 + x^2], x]",
    "D[E^x Sqrt[1 + 2 x^2], x]",
    "D[Log[x] Sqrt[1 + x^3], x]",
    "D[Log[x]^2 (1 + x^3)^(1/3), x]",
    "D[E^x Log[x] Sqrt[1 + x^3], x]",                       /* two monomials */
    "D[Log[Log[x] + Sqrt[1 + x^3]], x]",                    /* log part over k0 */
    "D[Log[E^x + Sqrt[1 + x^3]], x]",                       /* exponential theta, log part over k0 */
    "Together[4 (1/x - 3 x^2 Log[x]/(2 (1 + x^3))) Sqrt[1 + x^3]/(Log[x]^2 - 2 (1 + x^3))]", /* residues +-Sqrt[2] */
    "D[Sqrt[x] E^x Log[x], x]",
    "D[Log[x] Sqrt[1 + x^3] + Log[1 + x], x]",
    "D[(1 + x^2)^(1/3) E^(2 x)/x, x]",
    "D[Log[Log[x]] Sqrt[1 + x^4], x]",
    NULL
};

static const char* const certificates[][2] = {
    {"E^x/Sqrt[1 + x^2]", "no algebraic part: the coefficient Risch DE D b + (s + j eta') b = c has no solution over the base"},
    {"E^x/Sqrt[1 + x^3]", "no algebraic part: the coefficient Risch DE D b + (s + j eta') b = c has no solution over the base"},
    {"Log[x] Sqrt[1 + x^3]", "no algebraic part: the coefficient Risch DE D b + s b = c has no solution over the base"},
    {"Log[x]/Sqrt[1 + x^3]", "no algebraic part: the coefficient Risch DE D b + s b = c has no solution over the base"},
    {"Sqrt[x^2 - 3] Log[x]", "no algebraic part: the coefficient Risch DE D b + s b = c has no solution over the base"},
    {NULL, NULL}
};

static void test_integrals(void) {
    char buf[1024];
    for (int i = 0; closes[i]; i++) BR_ANTIDERIV(closes[i], "Bronstein");
    BR_TRUE("BronsteinIntegrate[D[E^x Sqrt[1 + x^2], x], x] === E^x Sqrt[1 + x^2]");
    BR_TRUE("BronsteinIntegrate[D[Log[x]^2 (1 + x^3)^(1/3), x], x] === (1 + x^3)^(1/3) Log[x]^2");
    for (int i = 0; certificates[i][0]; i++) {
        snprintf(buf, sizeof buf, "BronsteinIntegrate[%s, x] === {\"not elementary\", \"%s\"}",
                 certificates[i][0], certificates[i][1]);
        BR_TRUE(buf);
        snprintf(buf, sizeof buf, "Head[Integrate[%s, x, Method -> \"Bronstein\"]]", certificates[i][0]);
        BR_EQ(buf, "Integrate");
    }
    /* the certificate is announced, Mathematica-style */
    BR_EQ("Check[Integrate[E^x/Sqrt[1 + x^2], x, Method -> \"Bronstein\"]; \"no msg\", \"msg\", {Integrate::nonelem}]", "\"msg\"");
}

static void test_properties(void) {
    /* random D[F], F = (p0 + p1 theta) Y with theta Log or Exp and a genus >= 1 or n = 3 curve */
    br_run("SeedRandom[20261112]; brSR[] := RandomInteger[{-3, 3}];"
           " brBF[] := Module[{th = RandomChoice[{Log[x], E^x, E^(2 x), Log[1 + x]}],"
           "    y = RandomChoice[{Sqrt[1 + x^3], (1 + x^2)^(1/3), Sqrt[x^4 + 2], (x^3 - 2)^(1/3)}]},"
           "   ((brSR[] + x)/(x + 2 + brSR[]^2) + (1 + brSR[] x) th) y]; True");
    BR_TRUE("And @@ Table[Module[{F = brBF[], f, r}, f = Together[D[F, x]];"
            "   r = BronsteinIntegrate[f, x]; !ListQ[r] && PossibleZeroQ[D[r, x] - f]], {16}]");
}

int main(void) {
    br_init();
    test_front_end();
    test_genus0();
    test_k0();
    test_component_routines();
    test_integrals();
    test_properties();
    BR_DONE("bronstein_shapeb_tests");
    return 0;
}

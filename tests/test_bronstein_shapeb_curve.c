/* test_bronstein_shapeb_curve.c -- Shape B, the curve-monomial case: one
 * Log/Exp kernel whose argument contains the radical, on a curve of any genus
 * (src/internal/bronstein/BaseCurve.m).  Every routine directly; Bronstein's
 * ISSAC'91 integral; closes (each verified by PossibleZeroQ); certificates with
 * their exact reasons (each a theorem: the Risch DE on the curve, limited
 * integration over k0, Trager's algebraic-base certificate); the decline over a
 * larger tower; and a seeded property test.
 */
#include "bronstein_test_utils.h"

static void test_context(void) {
    br_run("Tb = brTower[E^(x Sqrt[x^3 + 1])/x, x]; Cb = brBCContext[Tb]; thb = Tb[\"syms\"][[1]]; Yb = Cb[\"Y\"];"
           " Tl = brTower[x Log[x + Sqrt[x^3 + 1]], x]; Cl = brBCContext[Tl]; thl = Tl[\"syms\"][[1]]; Yl = Cl[\"Y\"]; True");
    BR_TRUE("Cb[\"kind\"] === \"Exp\" && Cb[\"m\"] === 2 && Cb[\"G\"] === 1 + x^3");
    BR_TRUE("Cb[\"h\"] === x Yb");
    BR_TRUE("Together[Cb[\"Dh\"] - Yb (2 + 5 x^3)/(2 (1 + x^3))] === 0");
    BR_TRUE("Together[Cb[\"Dth\"] - thb Cb[\"Dh\"]] === 0");
    BR_TRUE("Cl[\"kind\"] === \"Log\" && Together[Cl[\"h\"] - (x + Yl)] === 0");
    BR_TRUE("brBCNorm[Cb, Yb^3] === Yb (1 + x^3)");
    BR_TRUE("brBCList[Cb, Yb + x] === {x, 1}");
    BR_TRUE("brBCExpr[Cb, {x, 1}] === Yb + x");
    BR_TRUE("Together[brBCD[Cb, Yb] - 3 x^2 Yb/(2 (1 + x^3))] === 0");
    BR_TRUE("brBCSurface[Cb, Yb thb] === E^(x Sqrt[1 + x^3]) Sqrt[1 + x^3]");
    BR_TRUE("brBCSplitDen[Cb, Yb/(thb^2 (x + 1))] === {Yb/(1 + x), thb^2}");
    BR_TRUE("Together[brBCFromSurface[Cl, x Sqrt[1 + x^3] + 1/Sqrt[1 + x^3]] - Yl (1 + x + x^4)/(1 + x^3)] === 0");
    BR_TRUE("brBCLogFree[x + Log[x] + 2 ArcTan[x] + RootSum[#^2 + 1 &, Log[x - #] &]] === x");
}

static void test_stages(void) {
    BR_TRUE("Together[brBCInverseMod[Cl, Yl, thl^2 + 1] - Yl/(1 + x^3)] === 0");
    BR_TRUE("brBCRemTheta[Cl, thl^3 Yl, thl^2 + 1] === -thl Yl");
    BR_TRUE("Module[{h = brBCHermite[Cl, brBCD[Cl, Yl/thl]]}, Together[h[\"H\"] - Yl/thl] === 0 && h[\"R\"] === 0]");
    BR_TRUE("brBCLaurent[Cb, Yb/thb + x thb^2] === {{-1, Yb}, {2, x}}");
    BR_TRUE("brBCSplit[Cb, (1 + thb)/(thb (thb - 1))] === <|\"L\" -> {{-1, -1}}, \"a\" -> 2, \"d\" -> thb - 1|>");
    BR_TRUE("Module[{lp = brBCLogPart[Cl, brBCD[Cl, thl^2 - 2], thl^2 - 2]}, lp[\"Rem\"] === 0 && Length[lp[\"Logs\"]] === 1]");
    BR_THROWS("brBCLogPart[Cl, 1, thl^2 - 2]", "NonElementary",
              "non-constant residues (resultant criterion over the algebraic base, [B] Prop 2.9)");
    BR_TRUE("brBCIntegrateK0[Cl, 3 x^2/(2 Yl)] === Sqrt[1 + x^3]");
    BR_THROWS("brBCIntegrateK0[Cl, thl]", "Internal", "brBCIntegrateK0Theta");
    BR_TRUE("brBCLimited[Cl, 1] === {x, 0}");
    BR_TRUE("brBCLimited[Cl, Cl[\"Dth\"]] === {0, 1}");
    BR_THROWS("brBCLimited[Cl, Yl]", "NonElementary",
              "algebraic base (Trager): the Hermite remainder is a nonzero holomorphic differential (no residues)");
    BR_TRUE("brBCPolyExp[Cb, {{1, brBCD[Cb, x Yb]}}] === {thb, 0}");
    BR_THROWS("brBCPolyExp[Cb, {{1, 1}}]", "NonElementary",
              "no algebraic part: the Risch DE D b + j h' b = c has no solution over the algebraic base ([B] Part 3)");
    BR_TRUE("brBCPolyLog[Cl, {{1, 1}, {0, brBCNorm[Cl, x Cl[\"Dth\"]]}}] === {thl x, 0}");
    BR_TRUE("Module[{f = D[E^(x Sqrt[x^3 + 1])/x, x], r}, r = brShapeBMono[brTower[f, x], f, x];"
            " r[\"verified\"] === False && PossibleZeroQ[D[r[\"F\"], x] - f]]");
}

static const char* const closes[] = {
    "(((5 x^4 + 2 x - 2)/x^2) (1 + 1/Sqrt[x^3 + 1]) + x/Sqrt[x^3 + 1]) E^(x Sqrt[x^3 + 1])",  /* ISSAC'91 */
    "D[x Log[x + Sqrt[x^3 + 1]], x]",
    "D[Log[x + Sqrt[1 + x^3]]^2, x]",
    "D[Sqrt[1 + x^3] E^Sqrt[1 + x^3], x]",
    "D[E^(x Sqrt[1 + x^3])/x, x]",
    "D[(x + Sqrt[1 + x^4]) Log[x^2 + Sqrt[1 + x^4]], x]",
    NULL
};

static const char* const certificates[][2] = {
    {"E^Sqrt[1 + x^3]", "no algebraic part: the Risch DE D b + j h' b = c has no solution over the algebraic base ([B] Part 3)"},
    {"Log[x + Sqrt[1 + x^3]]", "algebraic base (Trager): the Hermite remainder has a pole at infinity"},
    {"Sqrt[2 + 3 x^2] Log[x + Sqrt[2 + 3 x^2]]", "no algebraic part: the limited integration problem has no solution over the algebraic base"},
    {NULL, NULL}
};

static void test_integrals(void) {
    char buf[1024];
    for (int i = 0; closes[i]; i++) BR_ANTIDERIV(closes[i], "Bronstein");
    BR_TRUE("Together[BronsteinIntegrate[(((5 x^4 + 2 x - 2)/x^2) (1 + 1/Sqrt[x^3 + 1]) + x/Sqrt[x^3 + 1]) E^(x Sqrt[x^3 + 1]), x]"
            " - (2 + 2 Sqrt[x^3 + 1]) E^(x Sqrt[x^3 + 1])/x] === 0");
    for (int i = 0; certificates[i][0]; i++) {
        snprintf(buf, sizeof buf, "BronsteinIntegrate[%s, x] === {\"not elementary\", \"%s\"}",
                 certificates[i][0], certificates[i][1]);
        BR_TRUE(buf);
    }
    /* scope: a radical-argument kernel over or under another monomial declines */
    BR_TRUE("First[BronsteinIntegrate[Log[x] E^(x Sqrt[x^3 + 1]), x]] === \"implementation incomplete\"");
}

static void test_properties(void) {
    br_run("SeedRandom[20261113]; rq[] := RandomInteger[{-2, 2}];"
           " brCF[] := Module[{y = RandomChoice[{Sqrt[1 + x^3], Sqrt[x^4 + 2], Sqrt[x^3 - x + 3]}], th},"
           "   th = RandomChoice[{E^(x y), Log[x + y], E^y, Log[x^2 + y]}];"
           "   ((rq[] + x)/(x + 2 + rq[]^2) + y (1 + rq[] x)) th]; True");
    BR_TRUE("And @@ Table[Module[{F = brCF[], f, r}, f = Together[D[F, x]];"
            "   r = BronsteinIntegrate[f, x]; !ListQ[r] && PossibleZeroQ[D[r, x] - f]], {10}]");
}

int main(void) {
    br_init();
    test_context();
    test_stages();
    test_integrals();
    test_properties();
    BR_DONE("bronstein_shapeb_curve_tests");
    return 0;
}

/* test_bronstein_phase5.c -- Phase 5, milestone M1 of the Bronstein campaign
 * (Campaign D3; docs/notes/bronstein_extensions.tex, sections E2 and E1, "M1"):
 * two monomials over the curve k0 = Q(x)(Y), src/internal/bronstein/BaseTower.m.
 *
 *  - admission over the curve (E2.2): the structure-theorem span on normal
 *    forms (brDFSpan), the zero-offset rule (A3) for dependent kernels;
 *  - the field record and its levels (brDTContext, brDTLevel);
 *  - one level of integration (brDTInt) and the level-0 oracles: Trager
 *    (brDTIntegrateK0), its Liouville data, LD_0/PLD_0 (brDTPLD0), "is a
 *    derivative in k0" (brDTIsDerivK0), RDE_0 (brDTRDE0);
 *  - LI_1 (brDTLimited) and RDE_1 over k1 ([SI] Ch. 6: brDTRDE and its steps);
 *  - the note's examples: Example 1 (f1) and 2 (f2) decided, g3 declined, the
 *    E2 admission example (Log[Sqrt[1 + x^3]] beside Log[1 + x^3]) eliminated;
 *  - every M1 decline class with its exact reason;
 *  - two soundness fixes found on the way: brBCLogPart's residues for a d whose
 *    leading coefficient is not 1, and brKummerSplit over a radical kernel.
 */
#include "bronstein_test_utils.h"

static const char* const Y3 = "Sqrt[1 + x^3]";

static void setup(void) {
    /* Example 2 of E2.5: f2 = E^y + D(y + Log[x])/(y + Log[x]); the tower orders it {E^y, Log[x]} */
    br_run("f1 = (2 E^(x Sqrt[1 + x^3]) x + 5 E^(x Sqrt[1 + x^3]) x^4 - 2 E^(x Sqrt[1 + x^3]) Sqrt[1 + x^3]"
           " + 2 x^2 Sqrt[1 + x^3] Log[x + Sqrt[1 + x^3]])/(2 x^2 Sqrt[1 + x^3]);"
           " f2 = (2 E^Sqrt[1 + x^3] x + 3 x^3 + 2 E^Sqrt[1 + x^3] x^4 + 2 Sqrt[1 + x^3]"
           " + 2 E^Sqrt[1 + x^3] x Sqrt[1 + x^3] Log[x])/(2 x + 2 x^4 + 2 x Sqrt[1 + x^3] Log[x]);"
           " g3 = D[Log[x + Sqrt[1 + x^3]], x] Log[x + Sqrt[1 + x^3]]/(Log[2 x + 2 Sqrt[1 + x^3]] - Log[x + Sqrt[1 + x^3]]); True");
    /* Example 1's tower: t1 = E^(x y) (exponential), t2 = Log[x + y] (primitive) */
    br_run("T1 = brTower[Together[f1], x]; C1x = brDTContext[T1]; s1 = C1x[\"syms\"]; Y1 = C1x[\"Y\"]; True");
    /* a primitive bottom: t1 = Log[x], t2 = E^(t1^2 + y) (exponential top: RDE_1 over k0(Log[x])) */
    br_run("T2 = brTower[E^(Log[x]^2 + Sqrt[1 + x^3]) + Log[x], x];"
           " C2x = brDTContext[T2]; s2 = C2x[\"syms\"]; Y2 = C2x[\"Y\"]; True");
    (void)Y3;
}

static void test_curve_and_admission(void) {
    br_run("cv = brDTCurve[1 + x^3, 2, x]; Yc = cv[\"Y\"]; True");
    BR_TRUE("cv[\"m\"] === 2 && cv[\"G\"] === 1 + x^3 && cv[\"phi\"] === 1");
    BR_TRUE("Together[cv[\"sY\"] - 3 x^2/(2 (1 + x^3))] === 0");
    BR_TRUE("brDTCurveOf[{x + Sqrt[1 + x^3], 1/Sqrt[1 + x^3]}, x][\"G\"] === 1 + x^3");
    BR_TRUE("brDTCurveOf[{(1 + x^3)^(1/3)}, x][\"m\"] === 3");
    BR_THROWS("brDTCurve[x^2, 2, x]", "Unimplemented",
              "tower over the curve: Y^n - b is not absolutely irreducible over Q(x) (a constant-field extension, OL2)");
    BR_THROWS("brDTCurveOf[{x + Sqrt[1/(1 + x^3)], Sqrt[1 + x^3]}, x]", "Unimplemented",
              "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)");
    BR_TRUE("brDTToY[cv, (1 + x^3)^(3/2) + x/Sqrt[1 + x^3]] === Yc^3 + x/Yc");
    BR_THROWS("brDTToY[cv, Sqrt[2 + x^3]]", "Unimplemented",
              "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)");
    BR_TRUE("brDTRationalQ[x/2 + 1/(1 + x), {x}]");
    BR_FALSE("brDTRationalQ[Sqrt[2] x, {x}]");
    BR_FALSE("brDTRationalQ[Pi x, {x}]");
    BR_FALSE("brDTRationalQ[I x, {x}]");
    BR_FALSE("brDTRationalQ[a x, {x}]");
    /* the span: D Log[1 + x^3] = 2 D Log[y]; D(x y) is independent of D Log[x] */
    BR_TRUE("brDFSpan[cv, 3 x^2/(1 + x^3), {3 x^2/(2 (1 + x^3))}, {}] === {2}");
    BR_FALSE("brDFSpan[cv, brBCNorm[cv, D[x, x] Yc + x 3 x^2/(2 (1 + x^3)) Yc], {1/x}, {}]");
    BR_FALSE("brDFSpan[cv, 1/x, {}, {}]");
    BR_THROWS("brDFSpan[cv, 0, {1/x}, {}]", "Internal", "brDFSpanZero");
    /* branch differences: Log[1 + x^3] - 2 Log[Sqrt[1 + x^3]] vanishes everywhere; Log[x^2] - 2 Log[x] does not */
    BR_TRUE("brDTZeroOffsetQ[Log[1 + x^3] - 2 Log[Sqrt[1 + x^3]], x]");
    BR_FALSE("brDTZeroOffsetQ[Log[x^2] - 2 Log[x], x]");
    /* monomial data over the curve: Example 2's tower {E^y, Log[x]} */
    br_run("Tf2 = brTower[f2, x]; md = brDTMonoData[cv, Tf2[\"monos\"]]; True");
    BR_TRUE("Tf2[\"shape\"] === \"BTower\" && Tf2[\"kinds\"] === {\"Exp\", \"Log\"}");
    BR_TRUE("md[[1]][\"argY\"] === Yc && Together[md[[1]][\"eta\"] - 3 x^2 Yc/(2 (1 + x^3))] === 0");
    BR_TRUE("md[[2]][\"DY\"] === 1/x && md[[2]][\"eta\"] === 1/x");
    /* brDTAdmit directly: a new monomial, a dependent logarithm (r = 2, kappa = 1, offset 0), a fractional exponential */
    br_run("m1 = brTower[Log[Sqrt[1 + x^3]], x][\"monos\"]; me = brTower[E^(x Sqrt[1 + x^3]), x][\"monos\"]; True");
    BR_TRUE("First[brDTAdmit[x, m1, \"Log\", x, Log[x], cv]] === \"new\"");
    BR_TRUE("brDTAdmit[x, m1, \"Log\", 1 + x^3, Log[1 + x^3], cv] === {\"rewrite\", 2 m1[[1]][\"sym\"]}");
    BR_THROWS("brDTAdmit[x, me, \"Exp\", x Sqrt[1 + x^3]/2, E^(x Sqrt[1 + x^3]/2), cv]", "Unimplemented",
              "a dependent exponential over the curve with a fractional coefficient (an algebraic element; A3, OL6)");
    BR_TRUE("brDTAdmit[x, me, \"Exp\", 2 x Sqrt[1 + x^3], E^(2 x Sqrt[1 + x^3]), cv] === {\"rewrite\", me[[1]][\"sym\"]^2}");
    BR_THROWS("brDTAdmit[x, me, \"Exp\", x Sqrt[1 + x^3] + 1, E^(x Sqrt[1 + x^3] + 1), cv]", "Unimplemented",
              "a dependent exponential over the curve with a nonzero constant offset (a factor E^c: the constant field would leave Q; A3, OL2)");
}

static void test_declines(void) {
    /* g3 (E2, Example 3): Log[2 x + 2 z] - Log[x + z] = Log[2], kappa = 2 */
    BR_THROWS("brTower[g3, x]", "Unimplemented",
              "a dependent logarithm over the curve whose constant offset Log[kappa] is nonzero (kappa != 1: the constant field would leave Q; A3, OL2)");
    BR_TRUE("!MatchQ[BronsteinIntegrate[g3, x], {\"not elementary\", _}]");
    /* Log[x - z] beside Log[x + z], z = Sqrt[1 + x^2]: kappa = -1 (offset I Pi) */
    BR_THROWS("brTower[Log[x - Sqrt[1 + x^2]] Log[x + Sqrt[1 + x^2]], x]", "Unimplemented",
              "a dependent logarithm over the curve whose constant offset Log[kappa] is nonzero (kappa != 1: the constant field would leave Q; A3, OL2)");
    /* E^(x z + 1) beside E^(x z): a factor E */
    BR_THROWS("brTower[E^(x Sqrt[1 + x^3] + 1) + E^(x Sqrt[1 + x^3]) Log[x], x]", "Unimplemented",
              "a dependent exponential over the curve with a nonzero constant offset (a factor E^c: the constant field would leave Q; A3, OL2)");
    /* a fractional exponential coefficient: E^(x y/2) after E^(x y) (E^(x y/2) first is a rewrite: E^(x y) = t^2) */
    BR_TRUE("Module[{r = brCatch[brTower[E^(x Sqrt[1 + x^3]/2) Log[x] + E^(x Sqrt[1 + x^3]), x]]},"
            " brFailTag[r] === {\"Unimplemented\", \"a dependent exponential over the curve with a fractional coefficient (an algebraic element; A3, OL6)\"}"
            " || (!brFailQ[r] && Length[r[\"monos\"]] === 2)]");
    /* three monomials */
    BR_THROWS("brTower[Log[x] Log[x + Sqrt[1 + x^3]] E^Sqrt[1 + x^3], x]", "Unimplemented",
              "three or more monomials over the curve (Phase 5 M1 takes two)");
    /* several radicals / another spelling inside the kernels */
    BR_THROWS("brTower[Log[x + Sqrt[1 + x^3]] Log[x + Sqrt[2 + x^3]], x]", "Unimplemented",
              "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)");
    BR_THROWS("brTower[Log[x + Sqrt[1/(1 + x^3)]] Log[x + Sqrt[1 + x^3]], x]", "Unimplemented",
              "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)");
    /* a radicand over a monomial: Log[1 + Sqrt[Log[x]]] (G-mid-1) */
    BR_THROWS("brTower[Log[1 + Sqrt[Log[x]]]/x, x]", "Unimplemented",
              "a radical whose radicand is not rational in x inside a logarithm or exponential");
    /* a constant outside Q (H_Q): Gaussian, then a parameter */
    BR_THROWS("brDTContext[brTower[E^(I x) Log[x + Sqrt[1 + x^3]], x]]", "Unimplemented",
              "a constant other than a rational number in a tower over the curve (hypothesis H_Q; OL2)");
    BR_THROWS("Block[{tragerParams = {a}}, brDTContext[brTower[E^(a x) Log[x + Sqrt[1 + x^3]], x]]]", "Unimplemented",
              "a constant other than a rational number in a tower over the curve (hypothesis H_Q; OL2)");
    BR_THROWS("brTower[Log[1 + Sqrt[Log[x]^2 + x]] Log[x + Sqrt[1 + x^3]], x]", "Unimplemented",
              "a radical whose radicand is not rational in x inside a logarithm or exponential");
}

static void test_context_levels(void) {
    BR_TRUE("T1[\"shape\"] === \"BTower\" && T1[\"kinds\"] === {\"Exp\", \"Log\"}");
    BR_TRUE("C1x[\"argsY\"] === {x Y1, x + Y1}");
    BR_TRUE("Together[C1x[\"etas\"][[1]] - (2 + 5 x^3) Y1/(2 (1 + x^3))] === 0");
    br_run("L1 = brDTLevel[C1x, 1]; L2 = brDTLevel[C1x, 2]; L0 = brDTLevel[C1x, 0]; True");
    BR_TRUE("L1[\"th\"] === s1[[1]] && L1[\"kind\"] === \"Exp\" && L1[\"lsyms\"] === {s1[[1]]}");
    BR_TRUE("L2[\"th\"] === s1[[2]] && L2[\"lsyms\"] === s1 && L0[\"lsyms\"] === {} && L0[\"th\"] === None");
    BR_TRUE("Together[L1[\"Dh\"] - C1x[\"etas\"][[1]]] === 0");
    BR_TRUE("brDTSurface[C1x, Y1 s1[[1]] + s1[[2]]] === E^(x Sqrt[1 + x^3]) Sqrt[1 + x^3] + Log[x + Sqrt[1 + x^3]]");
    /* brBCD on a level record differentiates every monomial of the level */
    BR_TRUE("brBCNorm[L2, brBCD[L2, s1[[2]]] - C1x[\"Ds\"][[2]]] === 0 && brBCNorm[L2, brBCD[L2, s1[[1]]] - s1[[1]] C1x[\"etas\"][[1]]] === 0");
    BR_TRUE("brBCSurface[L2, s1[[1]]] === E^(x Sqrt[1 + x^3])");
}

static void test_level0_oracles(void) {
    /* INT_0 and its Liouville data */
    BR_TRUE("brDTIntegrateK0[C1x, 3 x^2/(2 Y1)] === Sqrt[1 + x^3]");
    BR_THROWS("brDTIntegrateK0[C1x, s1[[1]]]", "Internal", "brDTIntegrateK0Tower");
    BR_TRUE("brDTLiouville0[C1x, x + Log[x]/2] === {x, {{1/2, x}}}");
    BR_TRUE("brDTLiouville0[C1x, Sqrt[1 + x^3] + 3 Log[x + Sqrt[1 + x^3]]] === {Y1, {{3, x + Y1}}}");
    BR_THROWS("brDTLiouville0[C1x, Sqrt[2] Log[x]]", "Unimplemented",
              "an irrational logarithm coefficient in the integral over k0 inside LD/PLD/LI (gap I2)");
    BR_THROWS("brDTLiouville0[C1x, ArcTan[x]]", "Unimplemented",
              "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)");
    /* LD_0 / PLD_0 (Prop E1ld): E1.5's LD_0(-Dy) has no solution; 1/x = D x/x; 1/(2x) needs n0 = 2 */
    BR_TRUE("brDTPLD0[C1x, -3 x^2 Y1/(2 (1 + x^3)), 0] === {\"none\"}");
    BR_TRUE("brDTPLD0[C1x, 1/x, 0] === {\"radical\", 1, 0, {{1, x}}}");
    BR_TRUE("brDTPLD0[C1x, 1/(2 x), 0] === {\"radical\", 2, 0, {{1, x}}}");
    BR_TRUE("brDTPLD0[C1x, 0, 1/x] === {\"radical\", 1, 0, {}}");
    BR_TRUE("Module[{e = C1x[\"etas\"][[1]]}, brDTPLD0[C1x, brBCNorm[C1x, e + 1/x], e] === {\"radical\", 1, 1, {{1, x}}}]");
    BR_TRUE("brDTPLD0[C1x, Y1, 0] === {\"none\"}");                                  /* Int y is elliptic */
    /* is a derivative in k0 */
    BR_TRUE("brDTIsDerivK0[C1x, 3 x^2/(2 Y1)] === {\"yes\", Y1}");
    BR_TRUE("brDTIsDerivK0[C1x, 1/x] === {\"no\"}");
    BR_TRUE("brDTIsDerivK0[C1x, Y1] === {\"no\"}");
    BR_TRUE("brDTIsDerivK0[C1x, 0] === {\"yes\", 0}");
    /* RDE_0: E1.5's two equations */
    BR_TRUE("brDTRDE0[C1x, 3 x^2/(2 Y1), 1] === {\"NoSolution\"}");
    BR_TRUE("brDTRDE0[C1x, C1x[\"etas\"][[1]], brBCNorm[C1x, (2 + 5 x^3)/(2 x Y1) - 1/x^2]] === {\"Solution\", 1/x}");
    BR_THROWS("brDTRDE0[C1x, s1[[1]], 1]", "Internal", "brDTRDE0Tower");
}

static void test_int_levels(void) {
    /* INT_1 and INT_2 on Example 1's tower */
    BR_TRUE("Module[{r = brDTInt[C1x, 1, brBCNorm[C1x, brBCD[L1, s1[[1]]/x]]]}, Together[r[\"alg\"] - s1[[1]]/x] === 0 && r[\"logs\"] === {} && r[\"surf\"] === 0]");
    BR_TRUE("Module[{r = brDTInt[C1x, 2, brBCNorm[C1x, brBCD[L2, x s1[[2]]]]]}, Together[brBCNorm[L2, brBCD[L2, r[\"alg\"]] + D[r[\"surf\"], x] - brBCD[L2, x s1[[2]]]]] === 0]");
    BR_TRUE("Module[{r = brDTInt[C1x, 2, brBCNorm[C1x, brBCD[L2, s1[[2]] + Y1]/(s1[[2]] + Y1)]]}, r[\"logs\"] =!= {} && r[\"rs\"] === {}]");
    BR_TRUE("brDTInt[C1x, 0, 0] === <|\"alg\" -> 0, \"logs\" -> {}, \"rs\" -> {}, \"surf\" -> 0|>");
    /* f1 over ONE denominator (no D1 split possible): certified, chain (d) then (e) */
    BR_THROWS("brShapeBTower[T1, Together[f1], x]", "NonElementary",
              "algebraic base (Trager): the Hermite remainder has a pole at infinity");
    BR_TRUE("Module[{f = D[E^(x Sqrt[1 + x^3])/x + x Log[x + Sqrt[1 + x^3]], x], r},"
            " r = brShapeBTower[brTower[f, x], f, x]; r[\"verified\"] === False && PossibleZeroQ[D[r[\"F\"], x] - f]]");
}

static void test_li1(void) {
    /* LI_1(p; D t2), t2 = Log[x + y] over k1 = k0(E^(x y)) */
    BR_TRUE("brDTLimited[C1x, 1] === {x, 0}");
    BR_TRUE("brDTLimited[C1x, C1x[\"Ds\"][[2]]] === {0, 1}");
    BR_TRUE("brDTLimited[C1x, 0] === {0, 0}");
    BR_TRUE("Module[{w = brDTLimited[C1x, brBCNorm[L1, brBCD[L1, s1[[1]]/x] + 2 C1x[\"Ds\"][[2]]]]},"
            " Together[w[[1]] - s1[[1]]/x] === 0 && w[[2]] === 2]");
    BR_THROWS("brDTLimited[C1x, Y1]", "NonElementary",
              "no algebraic part: the limited integration problem has no solution over the tower on the curve (its integrand is not elementary over k1; [SI] Thm 5.8.1)");
    BR_THROWS("brDTLimited[C1x, 1/x]", "NonElementary",
              "no algebraic part: the limited integration problem has no solution over the tower on the curve (the logarithmic part of the integral over k1 is not in the span of the structure generators; [SI] Thm 5.8.1, Cor 9.3.1)");
}

static void test_rde1(void) {
    /* RDE_1 over k1 = k0(Log[x]) (primitive t1) and over k0(E^(x y)) (exponential t1) */
    br_run("P1 = brDTLevel[C2x, 1]; t = s2[[1]]; True");
    BR_TRUE("C2x[\"kinds\"] === {\"Log\", \"Exp\"}");
    BR_TRUE("brDTDeg[P1, t^2 Y2 + x] === 2 && brDTDeg[P1, 0] === -Infinity");
    BR_THROWS("brDTDeg[P1, 1/t]", "Internal", "brDTNotPolynomial");
    BR_TRUE("brDTCoeff[P1, t^2 Y2 + x t + 1, 1] === x");
    BR_TRUE("brDTDivide[P1, t^2 + Y2, t + 1] === {t - 1, 1 + Y2}");
    BR_TRUE("brDTGCD[P1, (t + Y2) (t - 1), (t + Y2) (t + x)] === t + Y2");
    BR_TRUE("brDTNormalDen[P1, 1/((t + x)^2 (t - Y2))] === Expand[(t + x)^2 (t^2 - 1 - x^3)]");
    BR_TRUE("brDTNormalDen[L1, 1/(s1[[1]]^2 (s1[[1]] + x))] === s1[[1]] + x");   /* the special t^2 dropped */
    /* degree bounds: no cancellation, LD_0 cancellation candidate (E1.5: LD_0(-Dy) has none) */
    BR_TRUE("brDTDegreeBound[C2x, 1, brBCNorm[P1, 3 x^2/(2 Y2)], 1] === 1");
    BR_TRUE("brDTDegreeBound[C2x, 1, t, t^3] === 2");
    BR_TRUE("brDTDegreeBound[C1x, 1, C1x[\"etas\"][[1]], s1[[1]]] === 1");
    /* SPDE: a of degree 0 returns at once; gcd(a, b) not dividing c is a proof */
    BR_TRUE("brDTSPDE[C2x, 1, t, t^2, 3] === {\"Reduced\", t, t^2, 3, 1, 0}");
    BR_TRUE("brDTSPDE[C2x, t, t, 1, 3] === {\"NoSolution\"}");
    BR_TRUE("brDTSPDE[C2x, t, 1, 0, 3] === {\"Zero\", 1, 0}");
    BR_TRUE("Module[{r = brDTSPDE[C2x, t + 1, 1, brBCNorm[P1, (t + 1) brBCD[P1, t^2] + t^2], 4]},"
            " r[[1]] === \"Reduced\" || r[[1]] === \"Zero\"]");
    /* D q = c in k0[t]: t = Log[x], D(t^2/2) = t/x; E1.5's c = Dy t is D(y t) - y/x, not a derivative */
    BR_TRUE("brDTPolyDeriv[C2x, t/x, 3] === {\"Solution\", t^2/2}");
    BR_TRUE("brDTPolyDeriv[C2x, t/x, 1] === {\"NoSolution\"}");
    BR_TRUE("brDTPolyDeriv[C2x, brBCNorm[P1, 3 x^2 t/(2 Y2)], 3] === {\"NoSolution\"}");
    BR_TRUE("brDTPolyDeriv[C1x, brBCNorm[L1, brBCD[L1, s1[[1]]/x + Y1]], 2] === {\"Solution\", brBCNorm[L1, s1[[1]]/x + Y1]}");
    /* after SPDE: the three cases */
    BR_TRUE("brDTPolyRDE[C2x, t, t^2 + 1/x, 3] === {\"Solution\", t}");                 /* b large */
    BR_TRUE("brDTPolyRDE[C2x, t, 1, 3] === {\"NoSolution\"}");
    BR_TRUE("brDTPolyRDE[C2x, brBCNorm[P1, 3 x^2/(2 Y2)], 1, 3] === {\"NoSolution\"}");   /* E1.5: cancellation, RDE_0(Dy, 1) */
    BR_TRUE("Module[{e = brBCNorm[P1, 3 x^2/(2 Y2)]}, brDTPolyRDE[C2x, e, brBCNorm[P1, brBCD[P1, t Y2] + e t Y2], 3] === {\"Solution\", t Y2}]");
    BR_THROWS("brDTPolyRDE[C2x, -1/x, 1, 3]", "Unimplemented",
              "the cancellation case of the Risch DE over k1 with a homogeneous solution (or one that needs ROOT_0)");
    /* RDE_1 end to end: E1.5's RDE_1(Dy, 1) over k0(Log[x]) has no solution; solvable ones verify */
    BR_TRUE("brDTRDE[C2x, brBCNorm[P1, 3 x^2/(2 Y2)], 1] === {\"NoSolution\"}");
    BR_TRUE("Module[{e = C2x[\"etas\"][[2]], q = t/x + Y2}, brDTRDE[C2x, e, brBCNorm[P1, brBCD[P1, q] + e q]] === {\"Solution\", brBCNorm[P1, q]}]");
    BR_TRUE("Module[{e = C2x[\"etas\"][[2]], q = 1/(t + x)}, brDTRDE[C2x, e, brBCNorm[P1, brBCD[P1, q] + e q]] === {\"Solution\", brBCNorm[P1, q]}]");
    BR_TRUE("brDTRDE[C2x, C2x[\"etas\"][[2]], 0] === {\"Solution\", 0}");
    /* exponential bottom: special denominators (q = 1/t1) */
    BR_TRUE("Module[{e = brBCNorm[L1, 2 C1x[\"etas\"][[1]]], q = brBCNorm[L1, x/s1[[1]] + Y1]},"
            " brDTRDE[C1x, e, brBCNorm[L1, brBCD[L1, q] + e q]] === {\"Solution\", q}]");
}

static void test_examples(void) {
    /* Example 2: f2 not elementary (RDE over the tower: Int E^y has none) */
    BR_EQ("BronsteinIntegrate[f2, x][[1]]", "\"not elementary\"");
    BR_THROWS("brIntegrateBody[f2, x]", "NonElementary",
              "no algebraic part: the Risch DE D b + i w' b = p_i has no solution over the tower on the curve ([SI] Thm 5.9.1, E1)");
    /* Example 1 over one denominator: decided without the D1 split */
    BR_THROWS("brIntegrateBody[Together[f1], x]", "NonElementary",
              "algebraic base (Trager): the Hermite remainder has a pole at infinity");
    /* E2's admission example: Log[Sqrt[1 + x^3]] beside Log[1 + x^3] (r = 2, offset 0) is eliminated, so
     * 1/(1 + Log[1 + x^3] - 2 Log[Sqrt[1 + x^3]]) is the integrand 1 */
    BR_ANTIDERIV("1/(1 + Log[1 + x^3] - 2 Log[Sqrt[1 + x^3]])", "Bronstein");
    BR_ANTIDERIV("D[x Log[Sqrt[1 + x^3]] Log[1 + x^3], x]", "Bronstein");
    /* closes over two monomials, each verified */
    BR_ANTIDERIV("D[E^(x Sqrt[1 + x^3]) Log[x], x]", "Bronstein");
    BR_ANTIDERIV("D[Log[x] Log[x + Sqrt[1 + x^3]] + Sqrt[1 + x^3] Log[x], x]", "Bronstein");
    BR_ANTIDERIV("D[E^Sqrt[1 + x^3]/(Log[x] + Sqrt[1 + x^3]), x]", "Bronstein");
    BR_ANTIDERIV("D[Log[Sqrt[1 + x^4] + E^Sqrt[1 + x^4] Log[1 + x Sqrt[1 + x^4]]], x]", "Bronstein");
    BR_ANTIDERIV("D[x Log[1 + (1 + x^3)^(1/3)] E^x, x]", "Bronstein");
    BR_ANTIDERIV("D[E^(x Sqrt[1 + x^2]) Log[x], x]", "Bronstein");
    /* a certificate whose Laurent coefficient is a non-elementary RDE over k1 */
    BR_EQ("BronsteinIntegrate[E^(x Sqrt[1 + x^3]) Log[x], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Log[x] Log[x + Sqrt[1 + x^3]], x][[1]]", "\"not elementary\"");
    /* 256 bgenusp-2-9 is f2 */
    BR_EQ("BronsteinIntegrate[(2*E^Sqrt[1 + x^3]*x + 3*x^3 + 2*E^Sqrt[1 + x^3]*x^4 + 2*Sqrt[1 + x^3] + 2*E^Sqrt[1 + x^3]*x*Sqrt[1 + x^3]*Log[x])"
          "/(2*x + 2*x^4 + 2*x*Sqrt[1 + x^3]*Log[x]), x][[1]]", "\"not elementary\"");
}

static void test_soundness_fixes(void) {
    /* brBCLogPart: a denominator whose leading coefficient in theta is x (residues were scaled by it) */
    BR_ANTIDERIV("D[Log[x Log[x + Sqrt[1 + x^3]] + 1], x]", "Bronstein");
    BR_ANTIDERIV("D[Log[x E^(x Sqrt[1 + x^3]) + 1], x]", "Bronstein");
    BR_ANTIDERIV("D[3 Log[Sqrt[1 + x^4] + E^Sqrt[1 + x^4] Log[1 + x Sqrt[1 + x^4]]], x]", "Bronstein");
    /* brKummerSplit: a kernel with a radical argument is not fixed by the sign changes */
    BR_TRUE("brKummerSplit[Sqrt[1 + x] Log[1 + Sqrt[1 + x^3]] + Sqrt[1 + x^3] Log[x], x] === None");
    BR_ANTIDERIV("-(1/(x*Sqrt[1 + x])) + Log[x]/(2*(1 + x)^(3/2)) + (9*x^2*Log[x])/(2*(1 + x)*(1 + Sqrt[1 + x^3]))"
                 " + (3*Sqrt[1 + x^3]*Log[1 + Sqrt[1 + x^3]])/(x*(1 + x)) + (9*x^2*Log[x]*Log[1 + Sqrt[1 + x^3]])/(2*(1 + x)*Sqrt[1 + x^3])"
                 " - (3*Sqrt[1 + x^3]*Log[x]*Log[1 + Sqrt[1 + x^3]])/(1 + x)^2", "Bronstein");
}

/* False certificates of the Galois split (v0.6645-v0.6650): a kernel whose
 * argument involves a radical is not fixed by that radical's sign change, so
 * brKummerSplit may not treat it as a scalar.  Each D[F] is elementary. */
static void test_kummer_split_kernel_radicals(void) {
    BR_TRUE("!MatchQ[BronsteinIntegrate[D[Sqrt[1 + x^2] Log[x + Sqrt[x]], x], x], {\"not elementary\", _}]");
    BR_TRUE("!MatchQ[BronsteinIntegrate[D[Sqrt[x] Log[1 + Sqrt[1 + x^3]], x], x], {\"not elementary\", _}]");
    BR_TRUE("!MatchQ[BronsteinIntegrate[D[Sqrt[1 + x] Log[1 + Sqrt[x]] + Sqrt[x] Log[x], x], x], {\"not elementary\", _}]");
    BR_TRUE("!MatchQ[BronsteinIntegrate[D[x E^(x Sqrt[1 + x^2])/Sqrt[x], x], x], {\"not elementary\", _}]");
    BR_TRUE("brKummerSplit[Sqrt[x] Log[1 + Sqrt[1 + x^3]], x] === None");
    /* brBCLogPart: the residues of the monic denominator (was scaled by lc) */
    BR_TRUE("With[{f = D[Log[x Log[x + Sqrt[1 + x^3]] + 1], x]}, With[{F = BronsteinIntegrate[f, x]}, "
            "!ListQ[F] && PossibleZeroQ[D[F, x] - f]]]");
}

int main(void) {
    br_init();
    setup();
    test_curve_and_admission();
    test_declines();
    test_context_levels();
    test_level0_oracles();
    test_int_levels();
    test_li1();
    test_rde1();
    test_examples();
    test_soundness_fixes();
    test_kummer_split_kernel_radicals();
    BR_DONE("bronstein_phase5_tests");
    return 0;
}

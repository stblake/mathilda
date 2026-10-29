/* test_trager_logterms.c -- the logarithmic part of a component:
 *   TragerLogTerms.m  tragerLogTerms, logPartResidual            ([T] Ch 5-6)
 *   KauersLogTerms.m  kauersLogTerms, groebnerLex, leadingCoefficientLex,
 *                     kauersPrincipal, kauersResidual, kauersTrace,
 *                     kauersTraceRational, kauersTraceNF, companionMatrix,
 *                     kauersLogExpressions                        ([K] 2008)
 *   LogToArcTan.m     tidyLog, dropNumericContent, expandRootSums,
 *                     complexSplit, realFormTerm, realForm
 *
 * Fixtures are component contexts built exactly as Pipeline.m's
 * integrateComponent builds them (helper tlCtx[c, m, G]: the component is
 * c(x) Y dx on Y^m = G, so c Y is the integrand), plus the Hermite remainder
 * under the key "rem".
 *
 * Independent ground truth:
 *   * Bronstein Ex. 6 ([B] p.20): Int dx/(x Sqrt[1-x^3]) = 1/3 Log[(y-1)/(y+1)];
 *     Ex. 4 ([B] p.19): Int dx/(x Sqrt[x^8+1]) = 1/4 Log[x^4/(1+y)].  The
 *     package's logands are compared with these up to a constant factor.
 *   * Traces: kauersTrace (the trace over Q[Z]/(q) of Z p'/p) against the
 *     direct numeric sum over the roots of q at a point of the curve, for
 *     hand-made p (not the heuristic's own).
 *   * Companion matrix: Tr C^j against Newton's identities.
 *   * Real form: D[realForm[F] - F] == 0 for complex-form antiderivatives F
 *     of real integrands, derivatives taken by hand.
 */
#include "trager_test_utils.h"

static void define_helpers(void) {
    char* r = tr_eval(
        "tlCtx[cc_, mm_Integer, GG_] := Module[{F0 = fieldQ[], B0, f0, sh, z0, w0, B1, her, r0, cr, num, den, rs},"
        "  B0 = curveBasis[mm, GG, x]; f0 = PadRight[{0, Together[cc]}, mm];"
        "  If[regularAtInfinityQ[F0, B0, f0], w0 = x; B1 = B0,"
        "    z0 = Unique[\"z\"]; sh = shiftCurve[cc, mm, GG, x, z0]; w0 = z0;"
        "    B1 = curveBasis[mm, sh[\"G\"], z0]; f0 = PadRight[{0, sh[\"c\"]}, mm]];"
        "  her = hermiteReduce[B1, f0]; r0 = her[\"Remainder\"]; cr = Together[r0[[2]]];"
        "  num = Numerator[cr]; den = Denominator[cr];"
        "  rs = componentResidues[num, den, mm, B1[\"G\"], w0];"
        "  Join[rs, <|\"B\" -> B1, \"genus\" -> curveGenus[mm, GG, x], \"N\" -> num, \"E\" -> den,"
        "    \"Ep\" -> Expand[D[den, w0]], \"Ysym\" -> Unique[\"Yk\"], \"rem\" -> r0|>]]; True");
    free(r);
    /* Numeric trace check: Sum_{q(g)=0} g p'(g)/p(g) at (w0, Y0) vs the exact
     * kauersTrace, p'(x, Y) = p_x + p_Y G' Y/(m G) (derivation on the curve). */
    r = tr_eval(
        "tlTraceOK[ctx_, q_, p_, w0_] := Module[{B = ctx[\"B\"], w, Yk = ctx[\"Ysym\"], Z = ctx[\"Z\"], G, m, Y0, rts, dp, num, sym},"
        "  w = B[\"x\"]; G = B[\"G\"]; m = B[\"m\"];"
        "  Y0 = N[(G /. w -> w0)^(1/m), 60];"
        "  rts = N[Z /. Solve[q == 0, Z], 60];"
        "  dp = D[p, w] + D[p, Yk] D[G, w] Yk/(m G);"
        "  num = Total[Function[g0, (g0 dp/p) /. Z -> g0] /@ rts] /. {w -> w0, Yk -> Y0};"
        "  sym = afToY[B, kauersTrace[ctx, q, p], Yk] /. {w -> w0, Yk -> Y0};"
        "  TrueQ[Abs[N[num - sym, 40]] < 10^-30]]; True");
    free(r);
    r = tr_eval(
        "tlEx6 = tlCtx[1/(x (1 - x^3)), 2, 1 - x^3];"
        "tlEx4 = tlCtx[1/(x (x^8 + 1)), 2, x^8 + 1];"
        "tlS6 = tlCtx[1/(x^2 - 2), 2, x^2 + 1];"          /* Sqrt[x^2+1]/(x^2-2): Q(Sqrt[6]) */
        "tlG0 = tlCtx[1/(x (x^2 + 1)), 2, x^2 + 1];"
        "tlQ4 = tlCtx[1/(x^4 + 1), 4, (x^4 + 1)^3];"      /* (x^4+1)^(-1/4) */
        "tlP35 = tlCtx[5/(x - 3), 2, x^3 - 2];"           /* infinite order */
        "tlP23 = tlCtx[1/((x - 2) (x^3 + 1)), 2, x^3 + 1];" /* 1/((x-2)Sqrt[x^3+1]) */
        "tlR8 = tlCtx[1/(x^2 + 1), 2, x^2 - 4 x - 1];"    /* residues in Q(Sqrt[1+-2I]) */
        "True");
    free(r);
}

/* ----------------------------------------------------------- tragerLogTerms */
static void test_trager_log_terms(void) {
    /* Bronstein Ex. 6: one term, coefficient 1/3, logand c (y-1)/(y+1) */
    TR_TRUE("Module[{lg = tragerLogTerms[tlEx6]}, Length[lg] == 1 && lg[[1, 1]] === 1/3]");
    TR_TRUE("Module[{lg = tragerLogTerms[tlEx6], u},"
            " u = afToY[tlEx6[\"B\"], lg[[1, 2]], Y] /. Y -> Sqrt[1 - x^3];"
            " TrueQ[PossibleZeroQ[D[u (Sqrt[1 - x^3] + 1)/(Sqrt[1 - x^3] - 1), x]]]]");
    /* Bronstein Ex. 4: coefficient 1/4, logand c x^4/(1+y) */
    TR_TRUE("Module[{lg = tragerLogTerms[tlEx4]}, Length[lg] == 1 && lg[[1, 1]] === 1/4]");
    TR_TRUE("Module[{lg = tragerLogTerms[tlEx4], u},"
            " u = afToY[tlEx4[\"B\"], lg[[1, 2]], Y] /. Y -> Sqrt[1 + x^8];"
            " TrueQ[PossibleZeroQ[D[u (1 + Sqrt[1 + x^8])/x^4, x]]]]");
    /* one term per Q-basis direction; coefficients are basis/N */
    TR_EQ("Length[tragerLogTerms[tlS6]]", "2");
    TR_TRUE("Module[{lg = tragerLogTerms[tlS6]}, lg[[All, 1]] === tlS6[\"basis\"]]");   /* genus 0: N = 1 */
    TR_EQ("Length[tragerLogTerms[tlQ4]]", "2");
    TR_TRUE("Module[{lg = tragerLogTerms[tlQ4]}, lg[[All, 1]] === (kc[tlQ4[\"F\"], #/4] & /@ tlQ4[\"basis\"])]");
    TR_TRUE("Length[tragerLogTerms[tlR8]] == Length[tlR8[\"basis\"]] && Length[tlR8[\"basis\"]] == 3");
    /* genus 0 always closes */
    TR_TRUE("ListQ[tragerLogTerms[tlG0]] && Length[tragerLogTerms[tlG0]] == 1");
    /* certificates: infinite order (orders mod p disagree) */
    TR_TRUE("MatchQ[catchTrager[tragerLogTerms[tlP35]],"
            " tragerFail[\"NonElementary\", s_String /; StringContainsQ[s, \"infinite order\"], _]]");
    /* certificate: (1, 3) on y^2 = x^3 + 8 (Nagell-Lutz: infinite order) */
    TR_TRUE("MatchQ[catchTrager[tragerLogTerms[tlCtx[1/((x - 1) (x^3 + 8)), 2, x^3 + 8]]],"
            " tragerFail[\"NonElementary\", _String, _]]");
    /* a torsion divisor (order 3) still produces a log term: the remainder minus
     * it is a nonzero holomorphic differential (checked below) */
    TR_TRUE("Module[{lg = tragerLogTerms[tlP23]}, Length[lg] == 1 && lg[[1, 1]] === tlP23[\"basis\"][[1]]/3]");
    /* every logand is a unit of the right kind: its divisor is N_j delta_j,
     * so its log-derivative has residues = coefficient times b: the
     * logPartResidual below is zero exactly. */
}

/* ----------------------------------------------------------- logPartResidual */
static void test_log_part_residual(void) {
    TR_TRUE("afZeroQ[tlEx6[\"F\"], logPartResidual[tlEx6[\"F\"], tlEx6[\"B\"], tlEx6[\"rem\"], tragerLogTerms[tlEx6]]]");
    TR_TRUE("afZeroQ[tlEx4[\"F\"], logPartResidual[tlEx4[\"F\"], tlEx4[\"B\"], tlEx4[\"rem\"], tragerLogTerms[tlEx4]]]");
    TR_TRUE("afZeroQ[tlS6[\"F\"], logPartResidual[tlS6[\"F\"], tlS6[\"B\"], tlS6[\"rem\"], tragerLogTerms[tlS6]]]");
    TR_TRUE("afZeroQ[tlG0[\"F\"], logPartResidual[tlG0[\"F\"], tlG0[\"B\"], tlG0[\"rem\"], tragerLogTerms[tlG0]]]");
    TR_TRUE("afZeroQ[tlQ4[\"F\"], logPartResidual[tlQ4[\"F\"], tlQ4[\"B\"], tlQ4[\"rem\"], tragerLogTerms[tlQ4]]]");
    TR_TRUE("afZeroQ[tlR8[\"F\"], logPartResidual[tlR8[\"F\"], tlR8[\"B\"], tlR8[\"rem\"], tragerLogTerms[tlR8]]]");
    /* perturbed coefficient -> nonzero */
    TR_FALSE("Module[{lg = tragerLogTerms[tlEx6]}, afZeroQ[tlEx6[\"F\"],"
             " logPartResidual[tlEx6[\"F\"], tlEx6[\"B\"], tlEx6[\"rem\"], {{1/2, lg[[1, 2]]}}]]]");
    TR_FALSE("Module[{lg = tragerLogTerms[tlEx4]}, afZeroQ[tlEx4[\"F\"],"
             " logPartResidual[tlEx4[\"F\"], tlEx4[\"B\"], tlEx4[\"rem\"], {{1/8, lg[[1, 2]]}}]]]");
    TR_FALSE("Module[{lg = tragerLogTerms[tlS6]}, afZeroQ[tlS6[\"F\"],"
             " logPartResidual[tlS6[\"F\"], tlS6[\"B\"], tlS6[\"rem\"], ReplacePart[lg, {2, 1} -> 2 lg[[2, 1]]]]]]");
    /* perturbed logand (u -> u + 1) -> nonzero */
    TR_FALSE("Module[{lg = tragerLogTerms[tlEx6]}, afZeroQ[tlEx6[\"F\"],"
             " logPartResidual[tlEx6[\"F\"], tlEx6[\"B\"], tlEx6[\"rem\"], {{1/3, afPlus[tlEx6[\"F\"], lg[[1, 2]], afOne[tlEx6[\"B\"]]]}}]]]");
    /* dropping a term -> nonzero */
    TR_FALSE("Module[{lg = tragerLogTerms[tlS6]}, afZeroQ[tlS6[\"F\"],"
             " logPartResidual[tlS6[\"F\"], tlS6[\"B\"], tlS6[\"rem\"], Most[lg]]]]");
    /* with no logs it is the remainder itself */
    TR_TRUE("logPartResidual[tlEx6[\"F\"], tlEx6[\"B\"], tlEx6[\"rem\"], {}] === tlEx6[\"rem\"]");
    /* the remainder of 1/((x-2) Sqrt[x^3+1]) minus its torsion log part is a
     * NONZERO holomorphic differential (FriCAS also proves it non-elementary) */
    TR_FALSE("afZeroQ[tlP23[\"F\"], logPartResidual[tlP23[\"F\"], tlP23[\"B\"], tlP23[\"rem\"], tragerLogTerms[tlP23]]]");
    /* the identity it encodes, by hand on Ex. 6: P r - c P' with P = u */
    TR_TRUE("Module[{lg = tragerLogTerms[tlEx6], F = tlEx6[\"F\"], B = tlEx6[\"B\"], u},"
            " u = lg[[1, 2]];"
            " logPartResidual[F, B, tlEx6[\"rem\"], lg] === afMinus[F, afTimes[F, B, tlEx6[\"rem\"], u], afScale[F, 1/3, afD[F, B, u]]]]");
}

/* ----------------------------------------------------------- groebnerLex, LC */
static void test_groebner_helpers(void) {
    TR_EQ("groebnerLex[{x^2 + y^2 - 1, x - y}, {x, y}]", "{x - y, -1/2 + y^2}");
    TR_EQ("groebnerLex[{2 x y - 4, 3 y^2 - 3}, {x, y}]", "{x - 2 y, -1 + y^2}");
    TR_EQ("groebnerLex[{3 x - 6}, {x}]", "{-2 + x}");
    TR_EQ("groebnerLex[{1, x}, {x}]", "{1}");
    /* canonical: the same ideal from different generators compares === */
    TR_TRUE("groebnerLex[{x^2 - 1, x y - 1}, {x, y}] === groebnerLex[{y^2 - 1, x - y, 5 (x^2 - 1)}, {x, y}]");
    /* every element monic in its lex leading term */
    TR_TRUE("AllTrue[groebnerLex[{x^3 - 2 x y, x^2 y - 2 y^2 + x}, {x, y}], leadingCoefficientLex[#, {x, y}] === 1 &]");
    TR_EQ("leadingCoefficientLex[3 x^2 y + 5 x y^3 + 7, {x, y}]", "3");
    TR_EQ("leadingCoefficientLex[3 x^2 y + 5 x y^3 + 7, {y, x}]", "5");
    TR_EQ("leadingCoefficientLex[-2 y + 4, {x, y}]", "-2");
    TR_EQ("leadingCoefficientLex[0, {x}]", "1");
    TR_EQ("leadingCoefficientLex[7, {x, y}]", "7");
    TR_EQ("leadingCoefficientLex[(2 x + 1)^2, {x}]", "4");
    TR_EQ("leadingCoefficientLex[x y z - 9 x^2 + z^5, {x, y, z}]", "-9");
}

/* ----------------------------------------------------------- kauersPrincipal */
static void test_kauers_principal(void) {
    /* on y^2 = x (genus 0) the point (1, 1) is cut out by y - 1 */
    TR_EQ("kauersPrincipal[groebnerLex[{y^2 - x, x - 1, y - 1}, {x, y}], {y^2 - x}, {x, y}]", "-1 + y");
    /* on y^2 = x^3 + 1 the flex (0, 1) is not principal (affine) but 3(0,1) is:
     * the tangent y = 1 meets the curve only there */
    TR_EQ("kauersPrincipal[groebnerLex[{y^2 - x^3 - 1, x, y - 1}, {x, y}], {y^2 - x^3 - 1}, {x, y}]", "None");
    TR_EQ("Module[{J = groebnerLex[{y^2 - x^3 - 1, x, y - 1}, {x, y}], J3},"
          " J3 = groebnerLex[Join[Flatten[Outer[Times, J, J, J]], {y^2 - x^3 - 1}], {x, y}];"
          " kauersPrincipal[J3, {y^2 - x^3 - 1}, {x, y}]]", "-1 + y");
    /* the ideal of the base itself has no candidate */
    TR_EQ("kauersPrincipal[groebnerLex[{y^2 - x}, {x, y}], {y^2 - x}, {x, y}]", "None");
    /* (2, 3) on y^2 = x^3 + 1: the tangent y = 2x - 1 meets at 2P + (0, -1), so
     * the 2nd power is not principal on the affine curve */
    TR_EQ("Module[{J = groebnerLex[{y^2 - x^3 - 1, x - 2, y - 3}, {x, y}], J2},"
          " J2 = groebnerLex[Join[Flatten[Outer[Times, J, J]], {y^2 - x^3 - 1}], {x, y}];"
          " kauersPrincipal[J2, {y^2 - x^3 - 1}, {x, y}]]", "None");
}

/* ----------------------------------------------------------- kauersLogTerms */
static void test_kauers_log_terms(void) {
    /* Ex. 6: factors Z -+ 1, power 3, logands Y -+ 1 (Bronstein's (y-1)/(y+1)) */
    TR_TRUE("Module[{kl = kauersLogTerms[tlEx6], Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]},"
            " Sort[kl] === Sort[{{-1 + Z, 3, -1 + Yk}, {1 + Z, 3, 1 + Yk}}]]");
    TR_TRUE("Module[{kl = kauersLogTerms[tlEx4], Z = tlEx4[\"Z\"], Yk = tlEx4[\"Ysym\"]},"
            " Sort[kl] === Sort[{{-1 + Z, 8, -1 + Yk}, {1 + Z, 8, 1 + Yk}}]]");
    /* every q is an irreducible factor of R, all factors covered */
    TR_TRUE("Module[{kl = kauersLogTerms[tlS6], Z = tlS6[\"Z\"]},"
            " Total[Exponent[#, Z] & /@ kl[[All, 1]]] == Exponent[tlS6[\"R\"], Z] &&"
            " PolynomialRemainder[tlS6[\"R\"], Expand[Times @@ kl[[All, 1]]], Z] === 0]");
    TR_TRUE("Module[{kl = kauersLogTerms[tlR8]}, Length[kl] == 3 && Sort[Exponent[#, tlR8[\"Z\"]] & /@ kl[[All, 1]]] === {1, 1, 4}]");
    /* no residues -> no terms */
    TR_EQ("kauersLogTerms[tlCtx[1/(x^4 + 1), 4, x^4 + 1]]", "{}");
    /* declines (never a certificate) where the heuristic does not close */
    TR_TRUE("MatchQ[catchTrager[kauersLogTerms[tlQ4]], tragerFail[\"Decline\", _String, _]]");
    TR_TRUE("MatchQ[catchTrager[kauersLogTerms[tlP35]], tragerFail[\"Decline\", _String, _]]");
    TR_TRUE("MatchQ[catchTrager[kauersLogTerms[tlP23]], tragerFail[\"Decline\", _String, _]]");
    /* KauersMaxPower bound: Ex. 6 needs k = 3 */
    TR_TRUE("Module[{old = tragerOption[\"KauersMaxPower\"], r0},"
            " tragerOption[\"KauersMaxPower\"] = 2;"
            " r0 = catchTrager[kauersLogTerms[tlEx6]];"
            " tragerOption[\"KauersMaxPower\"] = old;"
            " MatchQ[r0, tragerFail[\"Decline\", _, _]]]");
    TR_EQ("tragerOption[\"KauersMaxPower\"]", "16");
}

/* ----------------------------------------------------------- kauersResidual */
static void test_kauers_residual(void) {
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[tlEx6, tlEx6[\"rem\"], kauersLogTerms[tlEx6]]]");
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[tlEx4, tlEx4[\"rem\"], kauersLogTerms[tlEx4]]]");
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[tlS6, tlS6[\"rem\"], kauersLogTerms[tlS6]]]");
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[tlG0, tlG0[\"rem\"], kauersLogTerms[tlG0]]]");
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[tlR8, tlR8[\"rem\"], kauersLogTerms[tlR8]]]");
    /* wrong power k -> nonzero */
    TR_FALSE("Module[{kl = kauersLogTerms[tlEx6]}, afZeroQ[fieldQ[],"
             " kauersResidual[tlEx6, tlEx6[\"rem\"], ReplacePart[kl, {1, 2} -> 2]]]]");
    /* a dropped factor -> nonzero */
    TR_FALSE("Module[{kl = kauersLogTerms[tlS6]}, afZeroQ[fieldQ[], kauersResidual[tlS6, tlS6[\"rem\"], Most[kl]]]]");
    /* a wrong logand (p + 1) on the degree-2 factor -> nonzero */
    TR_FALSE("Module[{kl = kauersLogTerms[tlS6], i},"
             " i = First[Select[Range[Length[kl]], Exponent[kl[[#, 1]], tlS6[\"Z\"]] == 2 &]];"
             " afZeroQ[fieldQ[], kauersResidual[tlS6, tlS6[\"rem\"], ReplacePart[kl, {i, 3} -> kl[[i, 3]] + 1]]]]");
    /* no logs: the remainder */
    TR_TRUE("kauersResidual[tlEx6, tlEx6[\"rem\"], {}] === tlEx6[\"rem\"]");
}

/* ----------------------------------------------------------- traces */
static void test_kauers_trace(void) {
    /* linear q: the value at the root, = Z p'/p at Z = 1 */
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]}, tlTraceOK[tlEx6, Z - 1, Yk - 1, 1/3]]");
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]}, tlTraceOK[tlEx6, 2 Z + 3, x Yk + Z, 1/2]]");
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]},"
            " kauersTrace[tlEx6, Z - 1, Yk - 1] === kauersTraceRational[tlEx6, Z - 1, Yk - 1]]");
    /* quadratic and quartic q, hand-made p */
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]}, tlTraceOK[tlEx6, Z^2 - 3, Yk + Z x + 1, 1/2]]");
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]}, tlTraceOK[tlEx6, Z^2 + Z + 1, Yk - Z, 2/7]]");
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]}, tlTraceOK[tlEx6, Z^3 - 2, Z^2 Yk + x - Z, 1/5]]");
    TR_TRUE("Module[{Z = tlEx4[\"Z\"], Yk = tlEx4[\"Ysym\"]}, tlTraceOK[tlEx4, Z^4 - 4 Z^2 + 2, (x - Z) Yk + Z^3, 3/4]]");
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]},"
            " kauersTrace[tlEx6, Z^2 - 3, Yk + Z x + 1] === kauersTraceNF[tlEx6, Z^2 - 3, Yk + Z x + 1]]");
    /* non-monic q gives the same trace as its monic associate */
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], Yk = tlEx6[\"Ysym\"]},"
            " Together[kauersTrace[tlEx6, 5 Z^2 - 15, Yk + Z] - kauersTrace[tlEx6, Z^2 - 3, Yk + Z]] === {0, 0}]");
    /* the heuristic's own factors (shifted curve, w = z) */
    TR_TRUE("Module[{kl = kauersLogTerms[tlS6]}, AllTrue[kl, tlTraceOK[tlS6, #[[1]], #[[3]], 3/7] &]]");
    TR_TRUE("Module[{kl = kauersLogTerms[tlR8]}, AllTrue[kl, tlTraceOK[tlR8, #[[1]], #[[3]], 2/9] &]]");
    /* a trace of a p free of Z and Y: Sum gamma * p'/p = (Sum gamma) p'/p */
    TR_TRUE("Module[{Z = tlEx6[\"Z\"], tr}, tr = kauersTrace[tlEx6, Z^2 - 2 Z - 5, x^2 + 1];"
            " Together[tr - {2 (2 x)/(x^2 + 1), 0}] === {0, 0}]");
}

/* ----------------------------------------------------------- companionMatrix */
static void test_companion(void) {
    TR_EQ("companionMatrix[Z^2 - 3, Z]", "{{0, 3}, {1, 0}}");
    TR_EQ("companionMatrix[Z - 5, Z]", "{{5}}");
    /* char poly = mu */
    TR_TRUE("Expand[Det[Z IdentityMatrix[3] - companionMatrix[Z^3 - 2 Z + 5, Z]] - (Z^3 - 2 Z + 5)] === 0");
    TR_TRUE("Expand[Det[Z IdentityMatrix[5] - companionMatrix[Z^5 + 3 Z^4 - Z + 7/2, Z]] - (Z^5 + 3 Z^4 - Z + 7/2)] === 0");
    /* power sums by Newton: Z^3 - 2Z + 5 has e1 = 0, e2 = -2, e3 = -5:
     * p0 = 3, p1 = 0, p2 = e1^2 - 2 e2 = 4, p3 = e1 p2 - e2 p1 + 3 e3 = -15,
     * p4 = e1 p3 - e2 p2 + e3 p1 = 8 */
    TR_EQ("Table[Tr[MatrixPower[companionMatrix[Z^3 - 2 Z + 5, Z], j]], {j, 0, 4}]", "{3, 0, 4, -15, 8}");
    /* Z^2 + 1: roots +-I, power sums 2, 0, -2, 0, 2 */
    TR_EQ("Table[Tr[MatrixPower[companionMatrix[Z^2 + 1, Z], j]], {j, 0, 4}]", "{2, 0, -2, 0, 2}");
    /* (Z - 1)(Z - 2)(Z - 4): p_j = 1 + 2^j + 4^j */
    TR_TRUE("Table[Tr[MatrixPower[companionMatrix[Expand[(Z - 1) (Z - 2) (Z - 4)], Z], j]], {j, 0, 6}] ==="
            " Table[1 + 2^j + 4^j, {j, 0, 6}]");
}

/* ----------------------------------------------------------- kauersLogExpressions */
static void test_kauers_log_expressions(void) {
    /* linear factors: {gamma/k, p(gamma)} with Y substituted */
    TR_TRUE("Module[{e = kauersLogExpressions[tlEx6, kauersLogTerms[tlEx6], Y]},"
            " Sort[e] === Sort[{{1/3, -1 + Y}, {-1/3, 1 + Y}}]]");
    TR_TRUE("Module[{e = kauersLogExpressions[tlEx4, kauersLogTerms[tlEx4], Y]},"
            " Sort[e] === Sort[{{1/8, -1 + Y}, {-1/8, 1 + Y}}]]");
    /* a degree-2 factor becomes a RootSum triple {1/k, RootSum[..], \"RootSum\"} */
    TR_TRUE("Module[{e = kauersLogExpressions[tlS6, kauersLogTerms[tlS6], Y], rs},"
            " rs = Select[e, Length[#] == 3 &];"
            " Length[rs] == 1 && rs[[1, 3]] === \"RootSum\" && Head[rs[[1, 2]]] === RootSum && FreeQ[e, tlS6[\"Ysym\"]]]");
    /* the expressions integrate the remainder: Ex. 6 in x with Y = Sqrt[1-x^3]:
     * 1/3 Log[Y-1] - 1/3 Log[Y+1] has derivative 1/(x Sqrt[1-x^3]) */
    TR_TRUE("Module[{e = kauersLogExpressions[tlEx6, kauersLogTerms[tlEx6], Y], F0},"
            " F0 = Total[#[[1]] Log[#[[2]]] & /@ e] /. Y -> Sqrt[1 - x^3];"
            " TrueQ[PossibleZeroQ[D[F0, x] - 1/(x Sqrt[1 - x^3])]]]");
    TR_EQ("kauersLogExpressions[tlEx6, {}, Y]", "{}");
}

/* ----------------------------------------------------------- LogToArcTan */
static void test_log_to_arctan(void) {
    /* dropNumericContent */
    TR_EQ("dropNumericContent[6 x + 4]", "2 + 3 x");
    TR_EQ("dropNumericContent[-6 x - 4]", "2 + 3 x");
    TR_EQ("dropNumericContent[0]", "0");
    TR_EQ("dropNumericContent[7]", "1");
    TR_EQ("dropNumericContent[x/3 + 1/6]", "1 + 2 x");
    TR_EQ("dropNumericContent[Sqrt[2] x + Sqrt[2]]", "Sqrt[2] + Sqrt[2] x");   /* irrational content kept */
    TR_TRUE("FreeQ[Together[dropNumericContent[2 I x + 4 I]/(2 I x + 4 I)], x]");   /* a constant multiple */
    /* tidyLog: Log[c u] - Log[u] constant */
    TR_EQ("tidyLog[(6 x + 4)/(2 x^2 - 8), x]", "(2 + 3 x)/(-4 + x^2)");
    TR_EQ("tidyLog[3 (x + Sqrt[x^2 + 1])/(9 x), x]", "(x + Sqrt[1 + x^2])/x");
    TR_EQ("tidyLog[5, x]", "1");
    TR_TRUE("Module[{u = (12 x^3 - 4 + 8 Sqrt[x^3 + 1])/(6 x^2 + 18)},"
            " TrueQ[PossibleZeroQ[D[Log[tidyLog[u, x]] - Log[u], x]]]]");
    /* expandRootSums */
    TR_EQ("expandRootSums[RootSum[#^2 - 2 &, # Log[x - #] &]]", "Sqrt[2] Log[-Sqrt[2] + x] - Sqrt[2] Log[Sqrt[2] + x]");
    TR_EQ("expandRootSums[RootSum[# + 3 &, Log[x - #] &]]", "Log[3 + x]");
    TR_EQ("expandRootSums[tqa + RootSum[#^3 - 2 &, # Log[x - #] &]]", "tqa + RootSum[#1^3 - 2 &, #1 Log[x - #1] &]");
    TR_TRUE("Module[{e = RootSum[#^2 + # + 1 &, # Log[x - #] &]},"
            " FreeQ[expandRootSums[e], RootSum] && TrueQ[PossibleZeroQ[D[expandRootSums[e], x] + (x + 2)/(x^2 + x + 1)]]]");
    TR_EQ("expandRootSums[x^2]", "x^2");
    /* complexSplit: e = re + I im, both free of I */
    TR_EQ("complexSplit[(1 + 2 I)/(3 - I x)]", "{(3 - 2 x)/(9 + x^2), (6 + x)/(9 + x^2)}");
    TR_EQ("complexSplit[x + I tqy]", "{x, tqy}");
    TR_EQ("complexSplit[I]", "{0, 1}");
    TR_EQ("complexSplit[5]", "{5, 0}");
    TR_EQ("complexSplit[1/(1 + I)]", "{1/2, -1/2}");
    TR_TRUE("Module[{e = (x + 2 I Sqrt[x^2 + 1])/(x^2 - I x + 3), cs}, cs = complexSplit[e];"
            " FreeQ[cs, Complex] && TrueQ[PossibleZeroQ[cs[[1]] + I cs[[2]] - e]]]");
    TR_TRUE("Module[{e = (2 - I)^3 Sqrt[x]/(x + I)^2, cs}, cs = complexSplit[e];"
            " FreeQ[cs, Complex] && TrueQ[PossibleZeroQ[cs[[1]] + I cs[[2]] - e]]]");
    /* realFormTerm */
    /* Im part: -(1/2) ArcTan[B/A] with A = (x^2-1)/(x^2+1), B = -2x/(x^2+1).
     * arcTanTerm picks the orientation whose quotient denominator changes sign
     * fewer times: 2x (once) beats -1 + x^2 (twice), giving (1/2) ArcTan[A/B],
     * which differs from the other form by a piecewise constant. */
    TR_EQ("realFormTerm[I/2, (x - I)/(x + I), x]", "1/2 ArcTan[(1 - x^2)/(2 x)]");
    TR_TRUE("PossibleZeroQ[D[1/2 ArcTan[(1 - x^2)/(2 x)] - 1/2 ArcTan[(2 x)/(-1 + x^2)], x]]");
    TR_EQ("realFormTerm[1, x + 1, x]", "Log[1 + x]");
    TR_EQ("realFormTerm[1/2, 2 I x, x]", "1/2 Log[x]");          /* A = 0 */
    TR_EQ("realFormTerm[3, 2 x + 4, x]", "3 Log[2 + x]");        /* I-free: tidy */
    TR_TRUE("Module[{t = realFormTerm[1 + I, x + I, x]}, FreeQ[t, Complex] &&"
            " TrueQ[PossibleZeroQ[D[t, x] - ComplexExpand[Re[D[(1 + I) Log[x + I], x]]]]]]");
    /* realForm: D[realForm[F] - F] == 0 for complex antiderivatives of real
     * integrands */
    TR_TRUE("Module[{F0 = I/2 Log[(x - I)/(x + I)] + x^2, G0}, G0 = realForm[F0, x];"
            " FreeQ[G0, Complex] && TrueQ[PossibleZeroQ[D[G0 - F0, x]]]]");
    TR_TRUE("Module[{F0 = RootSum[4 #^2 + 1 &, # Log[x - 2 #] &], G0}, G0 = realForm[F0, x];"
            " FreeQ[G0, Complex] && FreeQ[G0, RootSum] && TrueQ[PossibleZeroQ[D[G0, x] + 1/(x^2 + 1)]]]");
    /* 1/Sqrt[1 - x^2] = -I Log[x (1 - I Y/x)] in complex form (integrateComponent) */
    TR_TRUE("Module[{r0 = integrateComponent[1/(1 - x^2), x, 2, 1 - x^2, \"Trager\"], F0, G0},"
            " F0 = r0[\"Algebraic\"] + Total[#[[1]] Log[#[[2]]] & /@ r0[\"Logs\"]] /. r0[\"Y\"] -> Sqrt[1 - x^2];"
            " G0 = realForm[F0, x];"
            " !FreeQ[F0, Complex] && FreeQ[G0, Complex] && TrueQ[PossibleZeroQ[D[G0, x] - 1/Sqrt[1 - x^2]]]]");
    /* (1 + x^2)^(-1) ArcTan-type with algebraic terms kept */
    TR_TRUE("Module[{F0 = Sqrt[x^2 + 1] + (1/2 - I/2) Log[x + I] + (1/2 + I/2) Log[x - I], G0},"
            " G0 = realForm[F0, x]; FreeQ[G0, Complex] && TrueQ[PossibleZeroQ[D[G0 - F0, x]]]]");
    /* an I-free answer is only tidied */
    TR_EQ("realForm[3 Log[(2 x + 2)/4], x]", "3 Log[1 + x]");
    TR_EQ("realForm[x^2/2, x]", "1/2 x^2");
}

/* The per-orbit construction (one Galois orbit of residues over Q[Z]/(q))
 * and the presentation helpers added with it. */
static void test_orbit_construction(void) {
    /* The remainder Y/((x^2 + 2)(x^2 + 1)) dx on Y^2 = x^2 + 1 (regular at
     * infinity, genus 0): residues at the places over x = +-I Sqrt[2], where
     * Y = +-I, are N Y/E' = +-1/(2 Sqrt[2]) -- the irreducible q = Z^2 - 1/8. */
    char* r = tr_eval(
        "toc = 1/((x^2 + 2) (x^2 + 1)); toB = curveBasis[2, 1 + x^2, x];"
        "toR = componentResiduePolynomial[1, (x^2 + 2) (x^2 + 1), 2, 1 + x^2, x];"
        "toX = Join[toR, <|\"B\" -> toB, \"genus\" -> 0, \"N\" -> 1, \"E\" -> Expand[(x^2 + 2) (x^2 + 1)],"
        " \"Ep\" -> Expand[D[(x^2 + 2) (x^2 + 1), x]], \"Ysym\" -> toY|>]; True");
    free(r);
    TR_TRUE("regularAtInfinityQ[fieldQ[], toB, {0, toc}]");
    TR_TRUE("PossibleZeroQ[toR[\"R\"] - (toR[\"Z\"]^2 - 1/8)]");
    TR_TRUE("KeyExistsQ[toR, \"Bp\"] && !KeyExistsQ[toR, \"F\"]");
    TR_TRUE("PossibleZeroQ[toR[\"Bp\"] - (x^2 + 2)]");
    /* reference point: ordinary (G and E nonzero there) */
    TR_EQ("referencePoint[toX]", "0");
    /* idealDegree = number of places with that residue: over each of the two
     * roots of x^2 + 2 exactly one of the two places carries the residue Z
     * (the other carries -Z), so the divisor has degree 2 */
    TR_EQ("idealDegree[fieldNF[toR[\"Z\"], toR[\"R\"], None], toB, residueIdeal[fieldNF[toR[\"Z\"], toR[\"R\"], None], toB, toR[\"Z\"], 1, toX[\"Ep\"], toX[\"Bp\"]]]", "2");
    /* the orbit construction closes (genus 0), its exact trace check passes,
     * k = m N = 2, and a perturbed coefficient fails the check */
    r = tr_eval("toO = tragerOrbitLogTerms[toX]; True");
    free(r);
    TR_TRUE("ListQ[toO] && Length[toO] === 1 && toO[[1, 2]] === 2");
    TR_TRUE("afZeroQ[fieldQ[], kauersResidual[toX, {0, toc}, toO]]");
    TR_FALSE("afZeroQ[fieldQ[], kauersResidual[toX, {0, toc}, ReplacePart[toO, {1, 2} -> 4]]]");
    /* both constructions give verified log parts for the same remainder */
    TR_TRUE("ListQ[tragerOrbitLogPart[toX, {0, toc}, toY]] && ListQ[tragerBasisLogPart[toX, {0, toc}, toY]]");
    /* logDerivativeNF: u'/u for u = 1 + Y Z over Q(i) on Y^2 = z^2 - 1 */
    TR_TRUE("Module[{K = fieldNF[Zl, Zl^2 + 1, None], Bl = curveBasis[2, z^2 - 1, z], ld},"
            " ld = logDerivativeNF[K, Bl, 1 + Yl Zl, Yl];"
            " PossibleZeroQ[(afToY[Bl, ld, Yl] - (Zl z Yl/(z^2 - 1))/(1 + Yl Zl)) /. {Yl -> Sqrt[z^2 - 1], Zl -> I}]]");
    TR_TRUE("logDerivativeNF[fieldNF[Zl, Zl^2 + 1, None], curveBasis[2, z^2 - 1, z], 7, Yl] === {0, 0}");
}

static void test_presentation_helpers(void) {
    /* arcTanTerm chooses the quotient whose denominator changes sign less:
     * -ArcTan[Sqrt[1 - x^2]/x] jumps at 0, ArcTan[x/Sqrt[1 - x^2]] does not */
    TR_EQ("arcTanTerm[1, Sqrt[1 - x^2], -x, x]", "ArcTan[x/Sqrt[1 - x^2]]");
    TR_TRUE("PossibleZeroQ[D[arcTanTerm[1, Sqrt[1 - x^2], -x, x] + ArcTan[-x/Sqrt[1 - x^2]], x]]");
    TR_EQ("signChanges[x, x]", "1");
    TR_EQ("signChanges[x^2 - 4, x]", "2");
    TR_EQ("signChanges[1 + x^2, x]", "0");
    TR_EQ("signChanges[Sqrt[1 - x^2], x]", "0");      /* complex values skipped */
    /* tidyRootSum applies the tidy function inside the held RootSum body */
    TR_TRUE("Module[{rs = RootSum[#^2 - 2 &, # Log[x^2 + 2 x + 1 - #] &]},"
            " tidyRootSum[rs, Factor] =!= rs && FreeQ[tidyRootSum[rs, Factor], Factor]]");
    TR_EQ("tidyRootSum[Log[x], Factor]", "Log[x]");
    /* evaluateRootSumBodies simplifies substituted radicals inside bodies */
    TR_EQ("evaluateRootSumBodies[RootSum[#^2 + 1 &, Function[Log[((1 + x)^(1/3))^2 - #] #]]]",
          "RootSum[#1^2 + 1 &, Log[(1 + x)^(2/3) - #1] #1 &]");
}

int main(void) {
    tr_init();
    define_helpers();
    test_trager_log_terms();
    test_log_part_residual();
    test_groebner_helpers();
    test_kauers_principal();
    test_kauers_log_terms();
    test_kauers_residual();
    test_kauers_trace();
    test_companion();
    test_kauers_log_expressions();
    test_log_to_arctan();
    test_orbit_construction();
    test_presentation_helpers();
    printf("trager_logterms: %d assertions passed\n", tr_assertions);
    return 0;
}

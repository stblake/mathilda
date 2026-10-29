/* test_trager_torsion.c -- src/internal/trager/TorsionModP.m: the order of a
 * residue divisor class by reduction modulo primes ([T] Ch 6 Sec 2-3
 * pp.66-70; [B] p.19).
 *
 * Functions covered: pPart, denominatorsPrimeToQ, nonzeroModQ, goodPrime,
 * reduceResidue, orderModP, goodPrimes, torsionOrder.
 *
 * Fixtures.  A component context is built exactly the way Pipeline.m's
 * integrateComponent builds it (shift at infinity if needed, Hermite
 * reduction, componentResidues, then B / genus / N / Ep / Ysym), by the
 * helper tqCtx[c, m, G] defined below: the component is c(x) Y dx on
 * Y^m = G(x).  Note that c Y is the integrand, so 1/(x Sqrt[1 - x^3]) is
 * tqCtx[1/(x (1 - x^3)), 2, 1 - x^3].
 *
 * Independent ground truth.
 *   * Orders modulo p on genus-1 curves Y^2 = x^3 + b are checked against
 *     brute-force elliptic-curve arithmetic over GF(p) written here (tqEcOrd2):
 *     for an integrand c Y with residues +-r at the two places (x0, +-y0),
 *     the residue divisor P - (-P) has class 2P in the group (O at infinity),
 *     so its order mod p is the order of 2P in E(GF(p)).
 *   * Residues are computed by hand: the residue of c(x) Y dx at a simple
 *     pole x0 of c with Y = y0 is y0 Res_{x0} c (invariant under the shift).
 *   * Torsion orders: Bronstein Ex. 6 (order 3: (0,1) has order 3 on
 *     y^2 = x^3 + 1 after x -> -x, and the divisor class is 2P), Bronstein
 *     Ex. 4 (order 4, [B] p.19), (x^4+1)^(-1/4) (order 4), (2,3) on
 *     y^2 = x^3 + 1 (order 6, so 2P has order 3).
 *   * Infinite order: (3, 5) on y^2 = x^3 - 2.  2P = (129/100, ...) is not
 *     integral, so by Nagell-Lutz P (and 2P) have infinite order; likewise
 *     (1, 3) on y^2 = x^3 + 8 (2P has x = -7/4).
 *   * Good primes: Legendre symbols ((2/p) = 1 iff p = +-1 mod 8), and
 *     discriminants computed by hand.
 */
#include "trager_test_utils.h"

/* Component-context builder (mirrors Pipeline.m integrateComponent) and a
 * brute-force elliptic-curve order over GF(p) for y^2 = x^3 + a x + b. */
static void define_helpers(void) {
    char* r = tr_eval(
        "tqCtx[cc_, mm_Integer, GG_] := Module[{F0 = fieldQ[], B0, f0, sh, z0, w0, B1, her, r0, cr, num, den, rs},"
        "  B0 = curveBasis[mm, GG, x]; f0 = PadRight[{0, Together[cc]}, mm];"
        "  If[regularAtInfinityQ[F0, B0, f0], w0 = x; B1 = B0,"
        "    z0 = Unique[\"z\"]; sh = shiftCurve[cc, mm, GG, x, z0]; w0 = z0;"
        "    B1 = curveBasis[mm, sh[\"G\"], z0]; f0 = PadRight[{0, sh[\"c\"]}, mm]];"
        "  her = hermiteReduce[B1, f0]; r0 = her[\"Remainder\"]; cr = Together[r0[[2]]];"
        "  num = Numerator[cr]; den = Denominator[cr];"
        "  rs = componentResidues[num, den, mm, B1[\"G\"], w0];"
        "  Join[rs, <|\"B\" -> B1, \"genus\" -> curveGenus[mm, GG, x], \"N\" -> num, \"E\" -> den,"
        "    \"Ep\" -> Expand[D[den, w0]], \"Ysym\" -> Unique[\"Yk\"]|>]]; True");
    free(r);
    r = tr_eval(
        "tqEcAdd[P_, Q_, a_, p_] := Which[P === \"O\", Q, Q === \"O\", P,"
        "  Mod[P[[1]] - Q[[1]], p] == 0 && Mod[P[[2]] + Q[[2]], p] == 0, \"O\","
        "  True, Module[{l, x3},"
        "    l = If[Mod[P[[1]] - Q[[1]], p] == 0,"
        "      Mod[(3 P[[1]]^2 + a) PowerMod[2 P[[2]], -1, p], p],"
        "      Mod[(Q[[2]] - P[[2]]) PowerMod[Q[[1]] - P[[1]], -1, p], p]];"
        "    x3 = Mod[l^2 - P[[1]] - Q[[1]], p]; {x3, Mod[l (P[[1]] - x3) - P[[2]], p]}]];"
        "tqEcOrd[P_, a_, p_] := Module[{Q = P, k = 1}, While[Q =!= \"O\", Q = tqEcAdd[Q, P, a, p]; k++]; k];"
        "tqEcOrd2[P_, a_, p_] := tqEcOrd[tqEcAdd[P, P, a, p], a, p];"
        "tqSqrtMod[c_, p_] := First[Select[Range[0, p - 1], Mod[#^2 - c, p] == 0 &]]; True");
    free(r);
    /* Shared fixtures (names chosen not to collide with package symbols). */
    r = tr_eval(
        "tqEx6 = tqCtx[1/(x (1 - x^3)), 2, 1 - x^3];"            /* Bronstein Ex. 6 */
        "tqEx4 = tqCtx[1/(x (x^8 + 1)), 2, x^8 + 1];"            /* Bronstein Ex. 4 */
        "tqP23 = tqCtx[1/(x - 2), 2, x^3 + 1];"                  /* (2, 3), order 6 */
        "tqP35 = tqCtx[5/(x - 3), 2, x^3 - 2];"                  /* (3, 5), infinite */
        "tqSq2 = tqCtx[1/(x - 1), 2, x^3 + 1];"                  /* (1, Sqrt[2]) */
        "tqG0 = tqCtx[1/(x (x^2 + 1)), 2, x^2 + 1];"             /* genus 0 */
        "tqQ4 = tqCtx[1/(x^4 + 1), 4, (x^4 + 1)^3];"             /* (x^4+1)^(-1/4) */
        "True");
    free(r);
}

/* ----------------------------------------------------------- pPart */
static void test_ppart(void) {
    TR_EQ("pPart[24, 2]", "8");
    TR_EQ("pPart[24, 3]", "3");
    TR_EQ("pPart[24, 5]", "1");
    TR_EQ("pPart[1, 7]", "1");
    TR_EQ("pPart[7, 7]", "7");
    TR_EQ("pPart[2^10 3^4 7, 2]", "1024");
    TR_EQ("pPart[2^10 3^4 7, 3]", "81");
    TR_EQ("pPart[2^10 3^4 7, 7]", "7");
    TR_EQ("pPart[2^10 3^4 7, 11]", "1");
    TR_EQ("pPart[5^20, 5]", "95367431640625");
    /* n = pPart[n, p] * (prime-to-p part) for a sweep */
    TR_TRUE("AllTrue[Range[1, 300], Function[n, AllTrue[{2, 3, 5, 7}, Mod[n/pPart[n, #], #] =!= 0 &]]]");
    TR_TRUE("AllTrue[Range[1, 300], Function[n, Times @@ (pPart[n, #] & /@ {2, 3, 5, 7, 11, 13, 17}) == "
            "Times @@ (#[[1]]^#[[2]] & /@ Select[FactorInteger[n], #[[1]] <= 17 &])]]");
    /* unevaluated on non-integers (pattern restricted) */
    TR_EQ("Head[pPart[3/2, 2]]", "pPart");
}

/* ----------------------------------------------------------- denominatorsPrimeToQ */
static void test_denominators(void) {
    TR_TRUE("denominatorsPrimeToQ[{1, 2, x^2 + 3}, 2]");
    TR_TRUE("denominatorsPrimeToQ[{1/3, x/5}, 7]");
    TR_FALSE("denominatorsPrimeToQ[{1/3, x/5}, 5]");
    TR_FALSE("denominatorsPrimeToQ[{1/3, x/5}, 3]");
    TR_FALSE("denominatorsPrimeToQ[x^2 + 1/6 x, 2]");
    TR_FALSE("denominatorsPrimeToQ[x^2 + 1/6 x, 3]");
    TR_TRUE("denominatorsPrimeToQ[x^2 + 1/6 x, 5]");
    TR_TRUE("denominatorsPrimeToQ[7/3, 7]");                    /* numerator divisible is fine */
    TR_FALSE("denominatorsPrimeToQ[{{1, {2, {3/14}}}}, 7]");    /* deep nesting */
    TR_TRUE("denominatorsPrimeToQ[{{1, {2, {3/14}}}}, 5]");
    TR_FALSE("denominatorsPrimeToQ[-5/22 t^3 + x, 11]");
    TR_TRUE("denominatorsPrimeToQ[{}, 3]");
    TR_TRUE("denominatorsPrimeToQ[0, 3]");
    /* the denominator 9 = 3^2 */
    TR_FALSE("denominatorsPrimeToQ[{x/9}, 3]");
}

/* ----------------------------------------------------------- nonzeroModQ */
static void test_nonzeromod(void) {
    TR_FALSE("nonzeroModQ[0, 5]");
    TR_FALSE("nonzeroModQ[10, 5]");
    TR_FALSE("nonzeroModQ[-35, 7]");
    TR_TRUE("nonzeroModQ[1, 2]");
    TR_TRUE("nonzeroModQ[11, 5]");
    TR_TRUE("nonzeroModQ[-1, 3]");
    TR_FALSE("nonzeroModQ[7/3, 7]");        /* 7 * 3^-1 = 0 mod 7 */
    TR_TRUE("nonzeroModQ[7/3, 5]");
    TR_TRUE("nonzeroModQ[1/2, 3]");
    TR_FALSE("nonzeroModQ[26/5, 13]");
    /* the discriminant of x^3 + 1 is -27 */
    TR_FALSE("nonzeroModQ[Discriminant[x^3 + 1, x], 3]");
    TR_TRUE("nonzeroModQ[Discriminant[x^3 + 1, x], 5]");
}

/* ----------------------------------------------------------- reduceResidue */
static void test_reduce_residue(void) {
    TR_EQ("reduceResidue[fieldQ[], 1/3, 7, None]", "5");         /* 3 * 5 = 15 = 1 */
    TR_EQ("reduceResidue[fieldQ[], -1, 7, None]", "6");
    TR_EQ("reduceResidue[fieldQ[], 21, 5, None]", "1");
    TR_EQ("reduceResidue[fieldQ[], -21, 11, None]", "1");
    TR_EQ("reduceResidue[fieldQ[], 0, 11, None]", "0");
    /* Q(Sqrt[2]) = Q[t]/(t^2 - 2): rho = 3 is a square root of 2 mod 7 */
    TR_EQ("Module[{tqF = fieldNF[tqt, tqt^2 - 2, Sqrt[2]]}, reduceResidue[tqF, 3 + tqt, 7, 3]]", "6");
    TR_EQ("Module[{tqF = fieldNF[tqt, tqt^2 - 2, Sqrt[2]]}, reduceResidue[tqF, tqt^2, 7, 3]]", "2");
    TR_EQ("Module[{tqF = fieldNF[tqt, tqt^2 - 2, Sqrt[2]]}, reduceResidue[tqF, tqt/2, 7, 3]]", "5");
    /* the image of t is a root of mu mod p: t^2 - 2 -> 0 */
    TR_EQ("Module[{tqF = fieldNF[tqt, tqt^2 - 2, Sqrt[2]]}, reduceResidue[tqF, tqt^2 - 2, 17, 6]]", "0");
    /* reduction is a ring homomorphism: (a b) -> a' b' */
    TR_TRUE("Module[{tqF = fieldNF[tqt, tqt^3 - tqt - 1, Root[#^3 - # - 1 &, 1]], tqa = 2 + tqt^2, tqb = 3 tqt - 1/2, p0 = 5, rh},"
            " rh = First[Select[Range[0, 4], Mod[#^3 - # - 1, 5] == 0 &]];"
            " Mod[reduceResidue[tqF, kc[tqF, tqa tqb], p0, rh] - reduceResidue[tqF, tqa, p0, rh] reduceResidue[tqF, tqb, p0, rh], p0] == 0]");
}

/* ----------------------------------------------------------- goodPrime */
static void test_good_prime(void) {
    /* sanity of the fixtures (residues by hand: Res = y0 Res_{x0} c) */
    TR_EQ("tqEx6[\"residues\"]", "{1, -1}");
    TR_EQ("Sort[tqP23[\"residues\"]]", "{-3, 3}");       /* y0 = 3 at x = 2 */
    TR_EQ("Sort[tqP35[\"residues\"]]", "{-25, 25}");     /* 5 * 5 */
    TR_EQ("tqEx6[\"genus\"]", "1");
    TR_EQ("tqEx4[\"genus\"]", "3");
    TR_EQ("tqQ4[\"genus\"]", "3");

    /* p | m: m = 2 rejects 2; m = 4 rejects 2 */
    TR_EQ("goodPrime[tqEx6, 2]", "None");
    TR_EQ("goodPrime[tqQ4, 2]", "None");
    /* bad reduction of the curve: 1 - x^3 = (1 - x)^3 mod 3 */
    TR_EQ("goodPrime[tqEx6, 3]", "None");
    /* good primes accepted, with rho = None over Q, and a frame centre a: a
     * root of 1 - x^3 mod p (multiplicity 1, prime to m = 2) */
    TR_TRUE("AllTrue[{5, 7, 11, 13, 17, 19, 23}, Module[{gp = goodPrime[tqEx6, #]},"
            " AssociationQ[gp] && gp[\"p\"] === # && gp[\"rho\"] === None &&"
            " Mod[1 - gp[\"a\"]^3, #] === 0] &]");

    /* p dividing a denominator: residues +-1/5 (otherwise identical to Ex. 6) */
    TR_TRUE("Module[{tqE = tqCtx[1/(5 x (1 - x^3)), 2, 1 - x^3]},"
            " Sort[tqE[\"residues\"]] === {-1/5, 1/5} && goodPrime[tqE, 5] === None &&"
            " goodPrime[tqE, 7] =!= None && goodPrime[tqEx6, 5] =!= None]");

    /* bad reduction of G: disc(x^3 + x + 3) = -4 - 243 = -247 = -13 * 19.
     * The integrand 1/((x+1) y) has residues +-1 at (-1, +-1). */
    TR_TRUE("Module[{tqB = tqCtx[1/((x + 1) (x^3 + x + 3)), 2, x^3 + x + 3]},"
            " Sort[tqB[\"residues\"]] === {-1, 1} &&"
            " goodPrime[tqB, 13] === None && goodPrime[tqB, 19] === None &&"
            /* the other small primes are good exactly when x^3 + x + 3 has a
             * root mod p (the reduced class arithmetic centres its frame there) */
            " AllTrue[{3, 5, 7, 11, 17, 23, 29},"
            "  (goodPrime[tqB, #] =!= None) === (Select[Range[0, # - 1], Function[r, Mod[r^3 + r + 3, #] == 0]] =!= {}) &]]");

    /* residue zero mod p: residues +-25 -> 5 rejected; disc(x^3 - 2) = -108
     * rejects 2 and 3; 7 .. 37 good */
    TR_EQ("goodPrime[tqP35, 5]", "None");
    TR_EQ("goodPrime[tqP35, 3]", "None");
    TR_TRUE("AllTrue[{7, 11, 13, 17, 19, 23, 29, 31, 37}, goodPrime[tqP35, #] =!= None &]");

    /* residues colliding mod p: 2 (3x+1) y/(x (x - 2)) on y^2 = x^3 + 1 has
     * residues +-1 (at x = 0: 2*1/(-2) y0) and +-21 (at x = 2: 2*7/2 * 3).
     * 21 = 1 mod 5, 21 = -1 mod 11, 21 = 0 mod 3, 7.  13, 17, 19 good. */
    TR_TRUE("Module[{tqA = tqCtx[2 (3 x + 1)/(x (x - 2)), 2, x^3 + 1]},"
            " Sort[tqA[\"residues\"]] === {-21, -1, 1, 21} &&"
            " goodPrime[tqA, 5] === None && goodPrime[tqA, 11] === None &&"
            " goodPrime[tqA, 7] === None && goodPrime[tqA, 3] === None &&"
            " AllTrue[{13, 17, 19, 23}, goodPrime[tqA, #] =!= None &]]");

    /* Bp losing squarefreeness: Bp = x^2 + x + 2 has discriminant -7, a double
     * root mod 7 (the two places, and so their residues, collide). */
    TR_TRUE("Module[{tqD = tqCtx[(2 x + 1)/(x^2 + x + 2), 2, x^3 + 2]},"
            " goodPrime[tqD, 7] === None]");

    /* residue-field minimal polynomial without a root mod p: residues
     * +-Sqrt[2]; 2 is a square mod p iff p = +-1 mod 8. */
    TR_TRUE("tqSq2[\"F\"][\"Type\"] === \"NF\" && Expand[tqSq2[\"F\"][\"mu\"] /. tqSq2[\"F\"][\"t\"] -> Sqrt[2]] === 0");
    TR_TRUE("AllTrue[{11, 13, 19, 29, 37, 43}, goodPrime[tqSq2, #] === None &]");
    TR_TRUE("AllTrue[{7, 17, 23, 31, 41}, Module[{gp = goodPrime[tqSq2, #]},"
            " AssociationQ[gp] && gp[\"p\"] === # && IntegerQ[gp[\"rho\"]] && Mod[gp[\"rho\"]^2 - 2, #] == 0] &]");
    /* the Jacobi-symbol rule over the first 60 odd primes > 5 (3 and 5 are
     * rejected for other reasons: disc(x^3+1) = -27; 5 kills the shift data
     * or not -- only the implication "accepted => (2/p) = 1" is asserted) */
    TR_TRUE("AllTrue[Prime[Range[4, 60]], goodPrime[tqSq2, #] === None || JacobiSymbol[2, #] === 1 &]");
    /* genus 0 contexts: every odd prime not dividing the data is good */
    /* accepted exactly when G has a root mod p for the reduced frame */
    TR_TRUE("Module[{Gw = tqG0[\"B\"][\"G\"], w = tqG0[\"B\"][\"x\"]}, AllTrue[{3, 5, 7, 11, 13},"
            " (goodPrime[tqG0, #] =!= None) === (Select[Range[0, # - 1], Function[r, Mod[Gw /. w -> r, #] == 0]] =!= {}) &]]");
}

/* ----------------------------------------------------------- goodPrimes */
static void test_good_primes(void) {
    TR_EQ("#[\"p\"] & /@ goodPrimes[tqEx6, 4]", "{5, 7, 11, 13}");
    TR_EQ("#[\"p\"] & /@ goodPrimes[tqP35, 3]", "{7, 11, 13}");
    TR_EQ("#[\"p\"] & /@ goodPrimes[tqSq2, 5]", "{7, 17, 23, 31, 41}");
    TR_EQ("Length[goodPrimes[tqEx6, 1]]", "1");
    /* every returned record is goodPrime's own, distinct, ascending */
    TR_TRUE("Module[{gs = goodPrimes[tqSq2, 6]}, AllTrue[gs, goodPrime[tqSq2, #[\"p\"]] === # &] &&"
            " OrderedQ[#[\"p\"] & /@ gs] && DuplicateFreeQ[#[\"p\"] & /@ gs]]");
    /* and no good prime was skipped: all primes below the largest one that
     * are not returned are bad */
    TR_TRUE("Module[{gs = #[\"p\"] & /@ goodPrimes[tqSq2, 5]},"
            " AllTrue[Complement[Prime[Range[PrimePi[Last[gs]]]], gs], goodPrime[tqSq2, #] === None &]]");
    /* a degree-8 residue field (residues of (2x+1) y/(x^2+x+2) on y^2 = x^3+2)
     * still finds two primes (split primes have density 1/8) */
    TR_TRUE("Module[{tqD = tqCtx[(2 x + 1)/(x^2 + x + 2), 2, x^3 + 2], gs},"
            " gs = catchTrager[goodPrimes[tqD, 2]];"
            " ListQ[gs] && Length[gs] == 2 && AllTrue[gs, IntegerQ[#[\"rho\"]] &]]");
    /* Decline when MaxPrimes is exhausted (temporarily lower the bound) */
    TR_TRUE("Module[{old = tragerOption[\"MaxPrimes\"], r0},"
            " tragerOption[\"MaxPrimes\"] = 2;"
            " r0 = catchTrager[goodPrimes[tqEx6, 2]];"
            " tragerOption[\"MaxPrimes\"] = old;"
            " MatchQ[r0, tragerFail[\"Decline\", \"no primes of good reduction found\", _]]]");
    TR_EQ("tragerOption[\"MaxPrimes\"]", "200");   /* restored */
    /* 3 tried (2 is skipped), 5 accepted: MaxPrimes = 2 suffices for one */
    TR_TRUE("Module[{old = tragerOption[\"MaxPrimes\"], r0},"
            " tragerOption[\"MaxPrimes\"] = 2;"
            " r0 = catchTrager[goodPrimes[tqEx6, 1]];"
            " tragerOption[\"MaxPrimes\"] = old;"
            " (#[\"p\"] & /@ r0) === {5}]");
}

/* ----------------------------------------------------------- orderModP */
static void test_order_mod_p(void) {
    /* Bronstein Ex. 6: order 3 at every good prime */
    TR_TRUE("AllTrue[goodPrimes[tqEx6, 6], orderModP[tqEx6, {1, -1}, #] === 3 &]");
    /* (2, 3) on y^2 = x^3 + 1: the class 2P, by brute-force EC arithmetic */
    TR_TRUE("AllTrue[goodPrimes[tqP23, 6], orderModP[tqP23, tqP23[\"coeffs\"][[All, 1]], #] ==="
            " tqEcOrd2[{2, 3}, 0, #[\"p\"]] &]");
    /* (3, 5) on y^2 = x^3 - 2 (infinite order over Q): orders 7, 6, 19, 9, 3 */
    TR_EQ("orderModP[tqP35, tqP35[\"coeffs\"][[All, 1]], #] & /@ goodPrimes[tqP35, 5]", "{7, 6, 19, 9, 3}");
    TR_TRUE("AllTrue[goodPrimes[tqP35, 8], orderModP[tqP35, tqP35[\"coeffs\"][[All, 1]], #] ==="
            " tqEcOrd2[{3, 5}, 0, #[\"p\"]] &]");
    /* residue field Q(Sqrt[2]): the point (1, rho) over GF(p) */
    TR_TRUE("AllTrue[goodPrimes[tqSq2, 5], orderModP[tqSq2, tqSq2[\"coeffs\"][[All, 1]], #] ==="
            " tqEcOrd2[{1, #[\"rho\"]}, 0, #[\"p\"]] &]");
    TR_EQ("orderModP[tqSq2, tqSq2[\"coeffs\"][[All, 1]], #] & /@ goodPrimes[tqSq2, 5]", "{3, 9, 6, 3, 21}");
    /* group structure: the order of k delta is o / gcd(o, k) */
    TR_TRUE("AllTrue[goodPrimes[tqP35, 5], Function[gp, Module[{o = orderModP[tqP35, {1, -1}, gp]},"
            " AllTrue[{2, 3, -1, 6}, orderModP[tqP35, # {1, -1}, gp] === o/GCD[o, #] &]]]]");
    /* the zero divisor has order 1; so does a principal one (genus 0) */
    TR_TRUE("AllTrue[goodPrimes[tqP35, 3], orderModP[tqP35, {0, 0}, #] === 1 &]");
    TR_TRUE("AllTrue[goodPrimes[tqG0, 3], orderModP[tqG0, tqG0[\"coeffs\"][[All, 1]], #] === 1 &]");
    /* Bronstein Ex. 4 (genus 3): order 4 */
    TR_TRUE("AllTrue[goodPrimes[tqEx4, 3], orderModP[tqEx4, {1, -1}, #] === 4 &]");
    /* bounded by the Weil bound (1 + Sqrt[p])^(2 g) */
    TR_TRUE("AllTrue[goodPrimes[tqP35, 5], orderModP[tqP35, {1, -1}, #] <= N[(1 + Sqrt[#[\"p\"]])^2] &]");
    /* MaxTorsionSearch exhausted -> Decline (order 3 > cap 2) */
    TR_TRUE("Module[{old = tragerOption[\"MaxTorsionSearch\"], r0},"
            " tragerOption[\"MaxTorsionSearch\"] = 1;"
            " r0 = catchTrager[orderModP[tqEx6, {1, -1}, First[goodPrimes[tqEx6, 1]]]];"
            " tragerOption[\"MaxTorsionSearch\"] = old;"
            " MatchQ[r0, tragerFail[\"Decline\", _String, _Association]]]");
    TR_EQ("tragerOption[\"MaxTorsionSearch\"]", "400");
}

/* ----------------------------------------------------------- torsionOrder */
static void test_torsion_order(void) {
    /* genus 0 -> 1 without any prime */
    TR_EQ("torsionOrder[tqG0, tqG0[\"coeffs\"][[All, 1]]]", "1");
    TR_EQ("torsionOrder[tqCtx[1/(x^2 - 2), 2, x^2 + 1], {1, -1, 0, 0}]", "1");
    /* Bronstein Ex. 6: 3; Ex. 4: 4 */
    TR_EQ("torsionOrder[tqEx6, {1, -1}]", "3");
    TR_EQ("torsionOrder[tqEx4, {1, -1}]", "4");
    /* (x^4+1)^(-1/4): component Y = yh^3 on Y^4 = (x^4+1)^3, residues
     * {1, -1, -I, I} (Q(I)), two basis directions, each of order 4 */
    TR_TRUE("tqQ4[\"F\"][\"Type\"] === \"NF\" && Length[tqQ4[\"residues\"]] == 4 && Length[tqQ4[\"basis\"]] == 2");
    TR_EQ("torsionOrder[tqQ4, #] & /@ Transpose[tqQ4[\"coeffs\"]]", "{4, 4}");
    /* (2, 3) on y^2 = x^3 + 1 has order 6, so 2P has order 3 */
    TR_EQ("tqEcOrd[{2, 3}, 0, 101]", "6");    /* sanity of the brute force: 6 | #E(F_101) */
    TR_EQ("torsionOrder[tqP23, {1, -1}]", "3");
    /* the candidate order of a multiple: 3 (P - (-P)) is principal -> 1 */
    TR_EQ("torsionOrder[tqP23, {3, -3}]", "1");
    TR_EQ("torsionOrder[tqEx4, {2, -2}]", "2");
    /* infinite order: (3, 5) on y^2 = x^3 - 2 (Nagell-Lutz) -> certificate */
    TR_TRUE("MatchQ[catchTrager[torsionOrder[tqP35, {1, -1}]], tragerFail[\"NonElementary\", _String, _Association]]");
    TR_TRUE("Module[{r0 = catchTrager[torsionOrder[tqP35, {1, -1}]]},"
            " r0[[3]][\"Primes\"] === {7, 11} && r0[[3]][\"OrderModP\"] === 7]");
    /* (1, 3) on y^2 = x^3 + 8 is of infinite order too: either the prime-to-pq
     * parts differ (certificate here) or the candidate N delta is not
     * principal (certificate in TragerLogTerms) -- never an order */
    TR_TRUE("Module[{tqH = tqCtx[1/(x - 1), 2, x^3 + 8], r0},"
            " r0 = catchTrager[torsionOrder[tqH, {1, -1}]];"
            " MatchQ[r0, tragerFail[\"NonElementary\", _, _]] ||"
            " (IntegerQ[r0] && !principalQ[tqH[\"F\"], tqH[\"B\"],"
            "   idealPower[tqH[\"F\"], tqH[\"B\"], idealTimes[tqH[\"F\"], tqH[\"B\"],"
            "     residueIdeal[tqH[\"F\"], tqH[\"B\"], tqH[\"residues\"][[1]], tqH[\"N\"], tqH[\"Ep\"], tqH[\"Bp\"]],"
            "     idealInverse[tqH[\"F\"], tqH[\"B\"], residueIdeal[tqH[\"F\"], tqH[\"B\"], tqH[\"residues\"][[2]], tqH[\"N\"], tqH[\"Ep\"], tqH[\"Bp\"]]]], r0]])]");
    /* residue field Q(Sqrt[2]): (1, Sqrt[2]) on y^2 = x^3 + 1 has infinite
     * order (2P has x = -7/8, whose 2-adic valuation -3 exceeds the torsion
     * bound) -> certificate */
    TR_TRUE("MatchQ[catchTrager[torsionOrder[tqSq2, tqSq2[\"coeffs\"][[All, 1]]]], tragerFail[\"NonElementary\", _, _]]");
    /* the returned N really makes N delta principal (Ex. 6, Ex. 4) */
    TR_TRUE("Module[{tqI = idealTimes[tqEx6[\"F\"], tqEx6[\"B\"],"
            "  residueIdeal[tqEx6[\"F\"], tqEx6[\"B\"], 1, tqEx6[\"N\"], tqEx6[\"Ep\"], tqEx6[\"Bp\"]],"
            "  idealInverse[tqEx6[\"F\"], tqEx6[\"B\"], residueIdeal[tqEx6[\"F\"], tqEx6[\"B\"], -1, tqEx6[\"N\"], tqEx6[\"Ep\"], tqEx6[\"Bp\"]]]]},"
            " principalQ[tqEx6[\"F\"], tqEx6[\"B\"], idealPower[tqEx6[\"F\"], tqEx6[\"B\"], tqI, 3]] &&"
            " !principalQ[tqEx6[\"F\"], tqEx6[\"B\"], tqI] &&"
            " !principalQ[tqEx6[\"F\"], tqEx6[\"B\"], idealPower[tqEx6[\"F\"], tqEx6[\"B\"], tqI, 2]]]");
    /* a torsion point of order 5 on 11a3: y^2 + y = x^3 - x^2, i.e.
     * Y^2 = 4x^3 - 4x^2 + 1 with Y = 2y + 1 and P = (0, 1).  The integrand
     * Y/x has residues +-1 at (0, +-1); 2P has order 5 as well. */
    TR_EQ("torsionOrder[tqCtx[1/x, 2, 4 x^3 - 4 x^2 + 1], {1, -1}]", "5");
    /* a point of order 4 on y^2 = x^3 + 4x (the 2-torsion of P = (2, 4) is
     * (0, 0)): 2P = (0, 0) has order 2, so P - (-P) has order 2 */
    TR_EQ("tqEcOrd[{2, 4}, 4, 1009]", "4");
    TR_EQ("torsionOrder[tqCtx[1/(x - 2), 2, x^3 + 4 x], {1, -1}]", "2");
}

int main(void) {
    tr_init();
    define_helpers();
    test_ppart();
    test_denominators();
    test_nonzeromod();
    test_reduce_residue();
    test_good_prime();
    test_good_primes();
    test_order_mod_p();
    test_torsion_order();
    printf("trager_torsion: %d assertions passed\n", tr_assertions);
    return 0;
}

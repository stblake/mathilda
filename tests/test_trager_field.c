/* test_trager_field.c
 *
 * Unit tests for src/internal/trager/Field.m -- the coefficient fields of the
 * Trager integrator (Q, the residue field Q(theta) = Q[t]/(mu), and GF(p)) and
 * the univariate polynomial / rational-function / linear-algebra routines
 * built on them.  Every later stage (HNF, normalisation at infinity, the
 * principality test, torsion mod p) runs field-generically through these
 * routines, so an error here becomes a wrong divisor or a false certificate
 * far downstream.  The suite therefore exercises every function on all three
 * field types, with small primes (GF(2), GF(3), GF(5), GF(7)) and number
 * fields of degree 2 (t^2 = 2) and 3 (t^3 = 2).
 *
 * Expected values are derived independently:
 *   - exact values computed by hand (e.g. 1/(1 + 2^(1/3)) = (1 - t + t^2)/3,
 *     the norm of x + t^2 over Q(2^(1/3)) is x^3 + 4, 1/3 = 2 in GF(5));
 *   - invariants checked against the plain builtins: a = q b + r with
 *     deg r < deg b, the Bezout identity s a + t b = g, T.A = R for the row
 *     echelon transform, c.A = 0 / c.A = b for the null and solve routines,
 *     kNorm against Resultant, rc against Together after clearing mu.
 *
 * Fields used throughout:  tQ = Q,  tN2 = Q(Sqrt[2]),  tN3 = Q(2^(1/3)),
 * tG2, tG3, tG5, tG7 = GF(2), GF(3), GF(5), GF(7).  Top-level names carry a
 * `t` prefix (tA1, ...) so they cannot collide with Trager`Private` symbols,
 * which are on $ContextPath.  Calls that may throw are wrapped in catchTrager.
 *
 * Bugs found while writing the suite are kept with their CORRECT expectation
 * in known_failures_field(), run only when TRAGER_KNOWN_FAIL is set.
 */
#include "trager_test_utils.h"

static void setup_fields(void) {
    char* s = tr_eval(
        "tQ = fieldQ[]; tN2 = fieldNF[t, t^2 - 2, Sqrt[2]];"
        " tN3 = fieldNF[t, t^3 - 2, 2^(1/3)];"
        " tG2 = fieldGF[2]; tG3 = fieldGF[3]; tG5 = fieldGF[5]; tG7 = fieldGF[7]; True");
    free(s);
}

/* ---- constructors and fieldDegree --------------------------------------- */
static void test_constructors(void) {
    TR_EQ("fieldQ[]", "<|\"Type\" -> \"Q\"|>");
    TR_EQ("fieldGF[5]", "<|\"Type\" -> \"GF\", \"p\" -> 5|>");
    TR_TRUE("tN2[\"Type\"] === \"NF\"");
    TR_EQ("tN2[\"mu\"]", "-2 + t^2");
    /* mu is stored expanded */
    TR_TRUE("fieldNF[t, (t - 1)(t + 1) - 1, Sqrt[2]][\"mu\"] === -2 + t^2");
    TR_EQ("tN3[\"theta\"]", "2^(1/3)");
    TR_EQ("{fieldDegree[tQ], fieldDegree[tG7], fieldDegree[tN2], fieldDegree[tN3]}",
          "{1, 1, 2, 3}");
}

/* ---- gfMod and kc --------------------------------------------------------- */
static void test_canonical(void) {
    /* 1/3 = 2 (mod 5): x^2/3 + 5x + 7 -> 2x^2 + 2 */
    TR_TRUE("gfMod[x^2/3 + 5 x + 7, 5] === 2 + 2 x^2");
    TR_EQ("gfMod[-1, 7]", "6");
    TR_EQ("gfMod[3 x + 3, 3]", "0");
    /* -1/2 = 3 (mod 7), since 2*3 = 6 = -1 */
    TR_EQ("gfMod[-1/2, 7]", "3");
    TR_TRUE("gfMod[x^3 - x, 3] === x^3 + 2 x");
    TR_TRUE("gfMod[(x + 1)^2, 2] === x^2 + 1");

    TR_TRUE("kc[tQ, (x + 1)^2] === 1 + 2 x + x^2");
    TR_TRUE("kc[tN2, t^3 + x t^2] === 2 t + 2 x");
    TR_TRUE("kc[tN2, t^4 x] === 4 x");
    /* (1 + t)^3 = 1 + 3t + 3t^2 + 2 over Q(2^(1/3)) */
    TR_TRUE("kc[tN3, (1 + t)^3] === 3 + 3 t + 3 t^2");
    TR_TRUE("kc[tN3, t^5] === 2 t^2");
    TR_TRUE("kc[tG5, 7 x^2 + x/2] === 3 x + 2 x^2");
    /* kc is idempotent on every field */
    TR_TRUE("AllTrue[{tQ, tN2, tN3, tG2, tG5}, kc[#, kc[#, (t + 2 x + 3)^4/7]] === kc[#, (t + 2 x + 3)^4/7] &]");

    TR_TRUE("kZeroQ[tN2, t^2 - 2]");
    TR_FALSE("kZeroQ[tQ, t^2 - 2]");
    TR_TRUE("kZeroQ[tG5, 5 x^3 - 10]");
    TR_FALSE("kZeroQ[tG5, 5 x^3 - 11]");
    TR_TRUE("kZeroQ[tN3, x (t^3 - 2)^2]");
    TR_TRUE("kZeroQ[tQ, 0]");
}

/* ---- scalar inverse ------------------------------------------------------- */
static void test_kinv(void) {
    TR_EQ("kInv[tQ, 3/4]", "4/3");
    TR_EQ("kInv[tG5, 2]", "3");
    TR_EQ("kInv[tG7, 3]", "5");
    TR_EQ("kInv[tG2, 1]", "1");
    /* a rational scalar over GF(7): (1/2)^-1 = 2 */
    TR_EQ("kInv[tG7, 1/2]", "2");
    TR_TRUE("kInv[tN2, 1 + t] === -1 + t");
    TR_TRUE("kInv[tN2, t] === t/2");
    TR_EQ("kInv[tN2, 2]", "1/2");
    TR_TRUE("kInv[tN3, t] === t^2/2");
    /* 1/(1 + 2^(1/3)) = (1 - t + t^2)/(1 + t^3) = (1 - t + t^2)/3 */
    TR_TRUE("kInv[tN3, 1 + t] === Expand[(1 - t + t^2)/3]");
    /* non-reduced input is reduced first: t^3 + t = 2 + t */
    TR_TRUE("kc[tN3, (t^3 + t) kInv[tN3, t^3 + t]] === 1");
    TR_TRUE("AllTrue[{1 + t + t^2, 3 - t^2, 5 t^2 - 7}, kc[tN3, # kInv[tN3, #]] === 1 &]");
    /* zero is an internal error on every field */
    TR_EQ("catchTrager[kInv[tQ, 0]]", "tragerFail[\"Internal\", \"kInvZero\", <||>]");
    TR_EQ("catchTrager[kInv[tG5, 10]]", "tragerFail[\"Internal\", \"kInvZero\", <||>]");
    TR_EQ("catchTrager[kInv[tN2, t^2 - 2]]", "tragerFail[\"Internal\", \"kInvZero\", <||>]");
    /* a zero divisor of a reducible "mu" is detected, not silently inverted */
    TR_EQ("catchTrager[kInv[fieldNF[t, t^2 - 1, 1], t - 1]][[2]] === \"kInvNotUnit\"", "True");
}

/* ---- degree, leading coefficient, monic ---------------------------------- */
static void test_deg_lc_monic(void) {
    TR_EQ("kDeg[tQ, 0, x]", "-Infinity");
    TR_EQ("kDeg[tQ, 5, x]", "0");
    TR_EQ("kDeg[tN2, (t^2 - 2) x^3 + x, x]", "1");
    TR_EQ("kDeg[tG5, 5 x^3 + x^2, x]", "2");
    TR_EQ("kDeg[tG3, 3 x + 6, x]", "-Infinity");

    TR_EQ("kLC[tQ, 0, x]", "0");
    TR_EQ("kLC[tQ, 3 x^2 - x, x]", "3");
    TR_TRUE("kLC[tN2, (t + t^2) x^2 + x, x] === 2 + t");
    TR_EQ("kLC[tN2, (t^2 - 2) x^2 + 7 x, x]", "7");
    TR_EQ("kLC[tG3, 4 x^2 + x, x]", "1");

    TR_TRUE("kMonic[tQ, 2 x^2 + 4, x] === 2 + x^2");
    TR_TRUE("kMonic[tN2, t x + 1, x] === t/2 + x");
    /* 1/3 = 2 (mod 5): 3x + 1 -> x + 2 */
    TR_TRUE("kMonic[tG5, 3 x + 1, x] === 2 + x");
    TR_EQ("kMonic[tG7, 0, x]", "0");
    TR_EQ("kMonic[tN3, 0, x]", "0");
    TR_EQ("kMonic[tN3, 3 t, x]", "1");
    /* monic means leading coefficient 1, on every field */
    TR_TRUE("AllTrue[{tN2, tN3}, kLC[#, kMonic[#, (2 + t) x^3 - 5 x + 1, x], x] === 1 &]");
    TR_TRUE("AllTrue[{tQ, tG3, tG7}, kLC[#, kMonic[#, 2 x^3 - 5 x + 1, x], x] === 1 &]");
}

/* ---- division with remainder --------------------------------------------- */
/* Assert a = q b + r, deg r < deg b over F (a property of kDivRem). */
#define DIVREM_OK(F, a, b) \
    TR_TRUE("Module[{qr = kDivRem[" F ", " a ", " b ", x]}," \
            " kZeroQ[" F ", (" a ") - qr[[1]] (" b ") - qr[[2]]] &&" \
            " kDeg[" F ", qr[[2]], x] < kDeg[" F ", " b ", x] &&" \
            " qr === {kc[" F ", qr[[1]]], kc[" F ", qr[[2]]]}]")

static void test_divrem(void) {
    TR_TRUE("kDivRem[tQ, x^3 + 1, x + 1, x] === {1 - x + x^2, 0}");
    TR_TRUE("kDivRem[tQ, x^3, x^2 + 1, x] === {x, -x}");
    TR_TRUE("kDivRem[tQ, 1, x, x] === {0, 1}");
    TR_TRUE("kDivRem[tQ, x^2, 2, x] === {x^2/2, 0}");
    /* GF(2): x^3 + x + 1 = (x + 1)(x^2 + x) + 1 */
    TR_TRUE("kDivRem[tG2, x^3 + x + 1, x + 1, x] === {x + x^2, 1}");
    /* GF(3): x^2 + 1 is irreducible, x^4 - 1 = (x^2 + 1)(x^2 - 1) */
    TR_TRUE("kDivRem[tG3, x^4 - 1, x^2 + 1, x] === {2 + x^2, 0}");
    /* GF(5): rational input coefficients are reduced first */
    TR_TRUE("kDivRem[tG5, x^2/2, 2 x, x] === {4 x, 0}");
    /* NF: (x^2 - 2) = (x - t)(x + t), x^3 - 2 = (x - t)(x^2 + t x + t^2) */
    TR_TRUE("kDivRem[tN2, x^2 - 2, x - t, x] === {t + x, 0}");
    TR_TRUE("kDivRem[tN3, x^3 - 2, x - t, x] === {t^2 + t x + x^2, 0}");
    TR_TRUE("kDivRemNF[tN3, x^3 - 2, x - t, x] === kDivRem[tN3, x^3 - 2, x - t, x]");
    /* division by a scalar of the field: (x^2 + 1)/(2t) = t/4 (x^2 + 1) */
    TR_TRUE("kDivRem[tN2, x^2 + 1, 2 t, x] === {t/4 + t x^2/4, 0}");
    /* x^2 = (x - t)(x + t) + t^2 = ... + 2 */
    TR_TRUE("kDivRem[tN2, x^2, x - t, x] === {t + x, 2}");
    TR_EQ("catchTrager[kDivRem[tN2, x, t^2 - 2, x]][[2]] === \"kDivRemByZero\"", "True");
    /* zero dividend */
    TR_TRUE("kDivRem[tN3, 0, x + t, x] === {0, 0}");
    TR_TRUE("kDivRem[tG7, 0, x + 1, x] === {0, 0}");

    DIVREM_OK("tQ", "3 x^5 - x^2/7 + 2", "2 x^2 + x - 1");
    DIVREM_OK("tG2", "x^7 + x^3 + 1", "x^3 + x + 1");
    DIVREM_OK("tG3", "2 x^6 + x^4 + 2 x + 1", "2 x^2 + 1");
    DIVREM_OK("tG7", "5 x^5 + 3 x + 4", "3 x^3 + 6");
    DIVREM_OK("tN2", "(1 + t) x^4 + t x + 3", "t x^2 - 1");
    DIVREM_OK("tN3", "t^2 x^5 + (1 - t) x^2 + t", "(1 + t) x^2 + t^2 x + 1");

    TR_TRUE("kQuo[tN3, x^3 - 2, x - t, x] === t^2 + t x + x^2");
    TR_TRUE("kRem[tN2, x^2, x - t, x] === 2");
    TR_TRUE("kQuo[tG2, x^3 + x + 1, x + 1, x] === x + x^2");
    TR_EQ("kRem[tG2, x^3 + x + 1, x + 1, x]", "1");
    TR_TRUE("kQuo[tQ, x^3, x^2 + 1, x] === x");
    TR_TRUE("kRem[tQ, x^3, x^2 + 1, x] === -x");

    TR_TRUE("kExactQuo[tQ, x^2 - 1, x + 1, x] === -1 + x");
    TR_TRUE("kExactQuo[tN2, x^2 - 2, x + t, x] === -t + x");
    TR_TRUE("kExactQuo[tG3, x^3 - x, x + 2, x] === x + x^2");
    TR_EQ("catchTrager[kExactQuo[tQ, x^2, x + 1, x]][[2]] === \"kExactQuoRemainder\"", "True");
    TR_EQ("catchTrager[kExactQuo[tN2, x^2, x - t, x]][[1]] === \"Internal\"", "True");
}

/* ---- extended gcd and gcd ------------------------------------------------- */
/* Assert the Bezout identity s a + t b = g and that g is monic (or 0). */
#define BEZOUT_OK(F, a, b) \
    TR_TRUE("Module[{eg = kExtGCD[" F ", " a ", " b ", x]}," \
            " kZeroQ[" F ", eg[[2]] (" a ") + eg[[3]] (" b ") - eg[[1]]] &&" \
            " (eg[[1]] === 0 || kLC[" F ", eg[[1]], x] === 1) &&" \
            " kZeroQ[" F ", kRem[" F ", " a ", eg[[1]], x]] &&" \
            " kZeroQ[" F ", kRem[" F ", " b ", eg[[1]], x]]]")

static void test_extgcd(void) {
    TR_EQ("kExtGCD[tQ, 0, 0, x]", "{0, 0, 0}");
    TR_EQ("kExtGCD[tN2, 0, t^2 - 2, x]", "{0, 0, 0}");
    TR_TRUE("kExtGCD[tQ, 2 x^2 + 2, 0, x] === {1 + x^2, 1/2, 0}");
    /* 3^-1 = 2 (mod 5): gcd(0, 3x + 1) = x + 2 with cofactor 2 */
    TR_TRUE("kExtGCD[tG5, 0, 3 x + 1, x] === {2 + x, 0, 2}");
    TR_TRUE("kExtGCD[tN2, t x, 0, x] === {x, t/2, 0}");
    /* gcd over Q(Sqrt[2]) of x^2 - 2 and (x - t)^2 is x - t, cofactors t/4, -t/4 */
    TR_TRUE("kExtGCD[tN2, x^2 - 2, x^2 - 2 t x + 2, x] === {-t + x, t/4, -t/4}");
    /* coprime over Q: s (x^2 + 1) + t (x + 1) = 1 with s = 1/2, t = (1 - x)/2 */
    TR_TRUE("kExtGCD[tQ, x^2 + 1, x + 1, x] === {1, 1/2, 1/2 - x/2}");

    BEZOUT_OK("tQ", "x^4 - 1", "x^3 + x^2 - x - 1");
    BEZOUT_OK("tQ", "3 x^3 + 2", "5 x^2 - 1");
    BEZOUT_OK("tG2", "x^4 + x + 1", "x^2 + 1");
    BEZOUT_OK("tG3", "x^5 + 2 x + 1", "x^3 + 2 x^2 + x");
    BEZOUT_OK("tG7", "3 x^4 + 6", "2 x^2 + 5 x");
    BEZOUT_OK("tN2", "x^3 - 2 x", "(x - t)^2 (x + 1)");
    BEZOUT_OK("tN3", "x^3 - 2", "x^2 - t^2");
    BEZOUT_OK("tN3", "(1 + t) x^3 + x", "t^2 x^2 + t");

    /* kExtGCDEuclid agrees with kExtGCD on Q and GF (where kExtGCD uses
       the builtin PolynomialExtendedGCD); the monic gcd is unique, the
       cofactors are checked by the Bezout identity. */
    TR_TRUE("kExtGCDEuclid[tQ, x^4 - 1, x^3 + x^2 - x - 1, x][[1]] === kExtGCD[tQ, x^4 - 1, x^3 + x^2 - x - 1, x][[1]]");
    TR_TRUE("Module[{eg = kExtGCDEuclid[tG3, x^5 + 2 x + 1, x^3 + 2 x^2 + x, x]},"
            " eg[[1]] === kExtGCD[tG3, x^5 + 2 x + 1, x^3 + 2 x^2 + x, x][[1]] &&"
            " kZeroQ[tG3, eg[[2]] (x^5 + 2 x + 1) + eg[[3]] (x^3 + 2 x^2 + x) - eg[[1]]]]");
    TR_TRUE("kExtGCDEuclid[tN3, x^3 - 2, x^2 - t^2, x][[1]] === -t + x");

    TR_TRUE("kGCD[tQ, x^2 - 1, x^2 + 2 x + 1, x] === 1 + x");
    TR_TRUE("kGCD[tQ, 2 x + 2, 4 x + 4, x] === 1 + x");
    /* over GF(2): x^2 + 1 = (x + 1)^2; over GF(3) x^2 + 1 is irreducible */
    TR_TRUE("kGCD[tG2, x^2 + 1, x + 1, x] === 1 + x");
    TR_EQ("kGCD[tG3, x^2 + 1, x + 1, x]", "1");
    /* over GF(5): x^2 + 1 = (x - 2)(x + 2) */
    /* x^2 + x + 3 vanishes at x = 3 (9 + 3 + 3 = 15) but not at x = 2, so the
       common factor is x - 3 = x + 2 */
    TR_TRUE("kGCD[tG5, x^2 + 1, x^2 + x + 3, x] === 2 + x");
    TR_TRUE("kGCD[tN3, x^3 - 2, x^2 - t^2, x] === -t + x");
    TR_EQ("kGCD[tN2, x^2 + 1, x - t, x]", "1");
    TR_EQ("kGCD[tQ, 0, 0, x]", "0");
}

/* ---- norm ------------------------------------------------------------------ */
static void test_norm(void) {
    TR_TRUE("kNorm[tN2, x - t] === -2 + x^2");
    /* Resultant[x - t, t^3 - 2, t] = (-1)^3 Prod (x - alpha) = 2 - x^3 */
    TR_TRUE("kNorm[tN3, x - t] === 2 - x^3");
    /* t^2 has minimal polynomial z^3 - 4: Prod (x + alpha^2) = x^3 + 4 */
    TR_TRUE("kNorm[tN3, x + t^2] === 4 + x^3");
    /* the norm of Sqrt[2] is -2, of 1 + Sqrt[2] is -1 */
    TR_EQ("kNorm[tN2, t]", "-2");
    TR_EQ("kNorm[tN2, 1 + t]", "-1");
    /* reduced first: t^3 - 2 + x is x, a constant in t: Res = x^3 */
    TR_TRUE("kNorm[tN3, t^3 - 2 + x] === x^3");
    TR_TRUE("kNorm[tQ, (x + 1)^2] === 1 + 2 x + x^2");
    TR_TRUE("kNorm[tG3, 4 x + 5] === 2 + x");
    /* the norm is the resultant against mu and is multiplicative */
    TR_TRUE("kNorm[tN3, (1 + t) x^2 + t^2] === Expand[Resultant[(1 + t) x^2 + t^2, t^3 - 2, t]]");
    /* As a resultant against the REDUCED representative, kNorm is the true
       norm Prod sigma(p) times (-1)^(deg_t p deg mu): exactly multiplicative
       when deg mu is even, multiplicative up to sign when it is odd (the
       callers only use it up to a unit: Ideal.m takes kMonic, rcNF divides). */
    TR_TRUE("Expand[kNorm[tN2, (x - t) (x + 1 + t)] - kNorm[tN2, x - t] kNorm[tN2, x + 1 + t]] === 0");
    TR_TRUE("Module[{u = kNorm[tN3, (x - t) (x + t^2)], v = kNorm[tN3, x - t] kNorm[tN3, x + t^2]},"
            " Expand[u - v] === 0 || Expand[u + v] === 0]");
    TR_TRUE("kMonic[tQ, kNorm[tN3, (x - t) (x + t^2)], x] === Expand[(x^3 - 2) (x^3 + 4)]");
    /* the norm is free of t */
    TR_TRUE("FreeQ[kNorm[tN3, t x^2 + (1 - t^2) x + 3 t], t]");
}

/* ---- rational functions -------------------------------------------------- */
static void test_rc(void) {
    TR_TRUE("rc[tQ, (x^2 - 1)/(x - 1)] === 1 + x");
    TR_EQ("rc[tQ, 0]", "0");
    TR_TRUE("rc[tQ, 1/x + 1/(x + 1)] === Together[1/x + 1/(x + 1)]");
    /* GF */
    TR_TRUE("rc[tG5, (x^2 - 1)/(x - 1)] === 1 + x");
    TR_EQ("rc[tG5, 5 x/(x + 1)]", "0");
    /* (3x + 1)/(2x + 4) = 3(x + 2)/(2(x + 2)) = 3/2 = 4 (mod 5) */
    TR_EQ("rc[tG5, (3 x + 1)/(2 x + 4)]", "4");
    /* denominator made monic: 1/(2x + 1) = 3/(x + 3) */
    TR_TRUE("rc[tG5, 1/(2 x + 1)] === 3/(3 + x)");
    /* constant denominator: (x + 1)/2 = 3x + 3 (mod 5) */
    TR_TRUE("rc[tG5, (x + 1)/2] === 3 + 3 x");
    /* GF(2): (x^2 + 1)/(x + 1) = x + 1 */
    TR_TRUE("rc[tG2, (x^2 + 1)/(x + 1)] === 1 + x");
    TR_TRUE("rc[tG3, (x + 2)/(x^2 + 1)] === (2 + x)/(1 + x^2)");
    /* NF: rationalised by the norm */
    TR_TRUE("rc[tN2, 1/(x - t)] === (t + x)/(-2 + x^2)");
    TR_TRUE("rc[tN2, 1/(1 + t)] === -1 + t");
    TR_EQ("rc[tN2, (t^2 - 2)/(x + 1)]", "0");
    TR_TRUE("rc[tN2, (t^3 + x)/x] === (2 t + x)/x");
    /* 1/(x^2 - t) over Q(2^(1/3)): (x^2 - t)(x^4 + t x^2 + t^2) = x^6 - 2 */
    TR_TRUE("Module[{e = rc[tN3, 1/(x^2 - t)]}, FreeQ[Denominator[e], t] &&"
            " PolynomialRemainder[Expand[Numerator[e] (x^2 - t) - Denominator[e]], t^3 - 2, t] === 0]");
    TR_TRUE("rc[tN3, 1/(x^2 - t)] === (t^2 + t x^2 + x^4)/(-2 + x^6)");
    /* a mixed denominator in t and x, with a numerator needing reduction */
    TR_TRUE("Module[{e = rc[tN3, (t^4 x + 1)/((1 + t) x + t^2)]}, FreeQ[Denominator[e], t] &&"
            " PolynomialRemainder[Expand[Numerator[e] ((1 + t) x + t^2) - (t^4 x + 1) Denominator[e]], t^3 - 2, t] === 0 &&"
            " Exponent[Numerator[e], t] < 3]");
    TR_TRUE("rcNF[tN2, 1/(x - t)] === rc[tN2, 1/(x - t)]");
    TR_TRUE("rcNF[tN3, (x + t)/(x - t)] === rc[tN3, (x + t)/(x - t)]");

    TR_TRUE("rcGFMonic[1, 2 x + 1, 5] === 3/(3 + x)");
    /* 3^-1 = 5 (mod 7) */
    TR_TRUE("rcGFMonic[x, 3, 7] === 5 x");
    TR_TRUE("rcGFMonic[x + 1, 4 x^2 + 1, 7] === (2 + 2 x)/(2 + x^2)");

    TR_TRUE("rcZeroQ[tN2, (t^2 - 2)/x]");
    TR_TRUE("rcZeroQ[tG3, 3 x/(x + 1)]");
    TR_TRUE("rcZeroQ[tQ, x/x - 1]");
    TR_FALSE("rcZeroQ[tQ, 1/x]");
    TR_FALSE("rcZeroQ[tN3, t/(x - t)]");
    /* exactness: 1/(x - t) - (x + t)/(x^2 - 2) is zero over Q(Sqrt[2]) */
    TR_TRUE("rcZeroQ[tN2, 1/(x - t) - (x + t)/(x^2 - 2)]");

    TR_TRUE("rcNumDen[tQ, (x^2 - 1)/(x^2 + 2 x + 1)] === {-1 + x, 1 + x}");
    TR_TRUE("rcNumDen[tN2, 1/(x - t)] === {t + x, -2 + x^2}");
    TR_TRUE("rcNumDen[tG5, 1/(2 x + 1)] === {3, 3 + x}");
    TR_TRUE("rcNumDen[tQ, 0] === {0, 1}");

    TR_EQ("rcDeg[tQ, 0, x]", "-Infinity");
    TR_EQ("rcDeg[tQ, (x^3 + 1)/(x + 2), x]", "2");
    TR_EQ("rcDeg[tQ, 1/x^2, x]", "-2");
    TR_EQ("rcDeg[tG5, (5 x^3 + x)/x^2, x]", "-1");
    /* the x^4 term vanishes over Q(2^(1/3)), leaving 1/x^4 */
    TR_EQ("rcDeg[tN3, ((t^3 - 2) x^4 + x)/x^5, x]", "-4");
    TR_EQ("rcDeg[tN2, 7, x]", "0");

    TR_TRUE("rcLeading[tN2, (2 t x^2 + 1)/(3 x^2 + x), x] === 2 t/3");
    TR_EQ("rcLeading[tQ, x^2/(2 x^2 + 1), x]", "1/2");
    /* GF(5): (2x + 1)/(3x) -> 2 * 3^-1 = 2 * 2 = 4 */
    TR_EQ("rcLeading[tG5, (2 x + 1)/(3 x), x]", "4");
    TR_EQ("rcLeading[tQ, 0, x]", "0");
    TR_TRUE("rcLeading[tN2, 1/(x - t), x] === 1");
    TR_TRUE("rcLeading[tN3, t^4 x/(x + 1), x] === 2 t");
}

/* ---- linear algebra -------------------------------------------------------- */
/* T.A = R over F, R is in reduced row echelon form with the given pivots. */
#define ECH_OK(F, A) \
    TR_TRUE("Module[{e = kRowEchelon[" F ", " A "]}," \
            " Map[kc[" F ", #] &, e[[3]] . (" A "), {2}] === e[[1]] &&" \
            " AllTrue[Range[Length[e[[2]]]], e[[1, #, e[[2, #]]]] === 1 &] &&" \
            " AllTrue[Range[Length[e[[2]]]], Count[e[[1, All, e[[2, #]]]], 0] === Length[e[[1]]] - 1 &] &&" \
            " !kZeroQ[" F ", Det[e[[3]]]]]")

static void test_linalg(void) {
    TR_EQ("kRowEchelon[tQ, {{1, 2}, {2, 4}}]", "{{{1, 2}, {0, 0}}, {1}, {{1, 0}, {-2, 1}}}");
    TR_EQ("kRowEchelon[tQ, {{0, 1}, {1, 0}}]", "{{{1, 0}, {0, 1}}, {1, 2}, {{0, 1}, {1, 0}}}");
    TR_EQ("kRowEchelon[tQ, {{2, 4, 6}}]", "{{{1, 2, 3}}, {1}, {{1/2}}}");
    TR_EQ("kRowEchelon[tQ, {{0, 0}, {0, 0}}][[2]]", "{}");
    ECH_OK("tQ", "{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}");
    ECH_OK("tQ", "{{0, 0, 1}, {0, 2, 3}, {0, 4, 6}}");
    ECH_OK("tG2", "{{1, 1, 0}, {0, 1, 1}, {1, 0, 1}}");
    ECH_OK("tG3", "{{1, 2, 0, 1}, {2, 1, 1, 0}, {0, 0, 2, 2}}");
    ECH_OK("tG7", "{{3, 5}, {6, 3}, {1, 1}}");
    ECH_OK("tN2", "{{1, t}, {t, 2}}");
    ECH_OK("tN3", "{{t, 1, 0}, {1, t^2, 1}, {t^2, 0, t}}");

    /* rank: det {{1,1,0},{0,1,1},{1,0,1}} = 2, so rank 3 over Q, 2 over GF(2) */
    TR_EQ("kRank[tQ, {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}}]", "3");
    TR_EQ("kRank[tG2, {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}}]", "2");
    /* det {{1, t}, {t, 2}} = 2 - t^2 = 0 over Q(Sqrt[2]) */
    TR_EQ("kRank[tN2, {{1, t}, {t, 2}}]", "1");
    TR_EQ("kRank[tN3, {{1, t}, {t, 2}}]", "2");
    TR_EQ("kRank[tQ, {}]", "0");
    TR_EQ("kRank[tQ, {{0, 0}, {0, 0}}]", "0");
    TR_EQ("kRank[tG5, {{5, 10}, {1, 2}}]", "1");
    TR_TRUE("kRank[tQ, {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}] === MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]");

    TR_EQ("kLeftNullVector[tQ, {{1, 2}, {2, 4}}]", "{-2, 1}");
    TR_EQ("kLeftNullVector[tQ, {{1, 0}, {0, 1}}]", "None");
    TR_TRUE("Module[{v = kLeftNullVector[tG2, {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}}]},"
            " v =!= None && AllTrue[v . {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}}, kZeroQ[tG2, #] &] && !AllTrue[v, kZeroQ[tG2, #] &]]");
    TR_TRUE("Module[{v = kLeftNullVector[tN2, {{1, t}, {t, 2}}]},"
            " v =!= None && AllTrue[v . {{1, t}, {t, 2}}, kZeroQ[tN2, #] &] && !AllTrue[v, kZeroQ[tN2, #] &]]");
    TR_EQ("kLeftNullVector[tN3, {{1, t}, {t, 2}}]", "None");
    /* more rows than columns always has a null vector */
    TR_TRUE("Module[{v = kLeftNullVector[tG7, {{3, 5}, {6, 3}, {1, 1}}]},"
            " v =!= None && AllTrue[v . {{3, 5}, {6, 3}, {1, 1}}, kZeroQ[tG7, #] &]]");

    TR_EQ("kSolveVector[tQ, {{1, 0}, {0, 1}}, {3, 4}]", "{3, 4}");
    TR_EQ("kSolveVector[tQ, {{1, 2}, {2, 4}}, {1, 0}]", "None");
    TR_TRUE("Module[{c = kSolveVector[tQ, {{1, 2}, {3, 4}}, {5, 6}]}, c . {{1, 2}, {3, 4}} === {5, 6}]");
    TR_TRUE("kSolveVector[tQ, {{1, 2}, {3, 4}}, {5, 6}] === LinearSolve[Transpose[{{1, 2}, {3, 4}}], {5, 6}]");
    TR_TRUE("Module[{c = kSolveVector[tG3, {{1, 2, 0}, {0, 1, 1}}, {1, 1, 2}]},"
            " c =!= None && Map[kc[tG3, #] &, c . {{1, 2, 0}, {0, 1, 1}}] === {1, 1, 2} && c === {1, 2}]");
    TR_EQ("kSolveVector[tG3, {{1, 2, 0}, {0, 1, 1}}, {0, 0, 1}]", "None");
    TR_TRUE("Module[{c = kSolveVector[tN2, {{1, t}, {t, 3}}, {1, 0}]},"
            " c =!= None && Map[kc[tN2, #] &, c . {{1, t}, {t, 3}}] === {1, 0}]");
    /* b = t * row 1 is in the span of the dependent rows {1, t}, {t, 2} */
    TR_TRUE("Module[{c = kSolveVector[tN2, {{1, t}, {t, 2}}, {t, 2}]},"
            " c =!= None && Map[kc[tN2, #] &, c . {{1, t}, {t, 2}}] === {t, 2}]");
    TR_EQ("kSolveVector[tN2, {{1, t}, {t, 2}}, {1, 0}]", "None");
    /* zero target: the zero combination (or any null combination) */
    TR_TRUE("Module[{c = kSolveVector[tQ, {{1, 0}, {0, 1}}, {0, 0}]}, c === {0, 0}]");
}

/* Bugs: kept with the correct expectation, run only under TRAGER_KNOWN_FAIL. */
/* Formerly known failures, now fixed: kRowEchelon of the empty matrix. rc over
 * Q is NOT canonical by design (it is Together's output; its docstring says to
 * compare by rcZeroQ of the difference), so equality is tested that way. */
static void test_fixed_known_failures(void) {
    TR_TRUE("kRowEchelon[tQ, {}] === {{}, {}, {}}");
    TR_EQ("kRank[tQ, {}]", "0");
    TR_TRUE("rcZeroQ[tQ, rc[tQ, x (1 + x)] - rc[tQ, x + x^2]]");
    TR_TRUE("rcZeroQ[tQ, rc[tQ, x (x + 1)/(x - 1)] - rc[tQ, (x^2 + x)/(x - 1)]]");
}

int main(void) {
    tr_init();
    setup_fields();
    test_constructors();
    test_canonical();
    test_kinv();
    test_deg_lc_monic();
    test_divrem();
    test_extgcd();
    test_norm();
    test_rc();
    test_linalg();
    test_fixed_known_failures();
    printf("trager_field: %d assertions passed\n", tr_assertions);
    return 0;
}

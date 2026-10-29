/* test_flint_polymat.c -- the FLINT` polynomial-matrix kernels
 * (src/poly/flint_polymat.c): PolynomialHermiteForm (with "DeterminantMultiple",
 * "Denominator", "Products"), PolynomialDetAdjugate, PolynomialWeakPopovForm,
 * PolynomialQuotientRemainder and PolynomialExtendedGCD over Q, GF(p) and
 * Q(theta).
 *
 * Two independent checks for every kernel, on seeded random inputs:
 *   - agreement with the Trager package's own routines (hnf, hnfModDet,
 *     idealReduce, idealTimes, reducedRows, elementInverse, kDivRemNF,
 *     kExtGCDEuclid) run with tragerOption["NativeKernels"] = False -- the
 *     reference implementation the kernels replace;
 *   - the defining invariants, computed without either: the Hermite form is
 *     upper triangular with monic pivots and reduced entries above them and
 *     every input row reduces to zero modulo it (same module); adj . M = det Id;
 *     a = q b + r with deg r < deg b; s a + t b = g, g monic dividing a and b.
 * Plus the out-of-scope contract: the call stays unevaluated.
 */
#include "trager_test_utils.h"

static void setup(void) {
    char* r = tr_eval(
        /* the reference: package routines with the native kernels switched off */
        "SetAttributes[pmRef, HoldAll];"
        "pmRef[e_] := Module[{v}, tragerOption[\"NativeKernels\"] = False;"
        "  v = Catch[e, \"TragerFail\"]; tragerOption[\"NativeKernels\"] = True; v];"
        /* fields */
        "pmQ = fieldQ[];"
        "pmNF2 = fieldNF[t, t^2 - 2, Sqrt[2]];"
        "pmNF3 = fieldNF[t, t^3 - t - 1, Root[#^3 - # - 1 &, 1]];"
        "pmNF6 = fieldNF[t, t^6 - 3, 3^(1/6)];"
        "pmFields = {pmQ, fieldGF[2], fieldGF[3], fieldGF[101], fieldGF[1000003], pmNF2, pmNF3, pmNF6};"
        /* the package reference is slow over the sextic field: compared against
         * it on these, the sextic field checked by the invariants */
        "pmRefFields = {pmQ, fieldGF[2], fieldGF[3], fieldGF[101], fieldGF[1000003], pmNF2, pmNF3};"
        /* random field elements and polynomials (seeded by the caller) */
        "pmCoef[F_] := Switch[F[\"Type\"],"
        "  \"Q\", RandomInteger[{-9, 9}]/RandomInteger[{1, 4}],"
        "  \"GF\", RandomInteger[{0, F[\"p\"] - 1}],"
        "  \"NF\", Sum[RandomInteger[{-5, 5}]/RandomInteger[{1, 3}] F[\"t\"]^k, {k, 0, fieldDegree[F] - 1}]];"
        "pmPoly[F_, d_] := kc[F, Sum[If[RandomInteger[{0, 3}] === 0, 0, pmCoef[F]] x^k, {k, 0, d}]];"
        "pmMat[F_, r_, c_, d_] := Table[pmPoly[F, RandomInteger[{0, d}]], {r}, {c}];"
        /* invariants of a Hermite form H of the module of `rows` */
        "pmHermiteQ[F_, H_] := Module[{n = Length[H]},"
        "  Length[First[H]] === n &&"
        "  AllTrue[Range[n], Function[j, kLC[F, H[[j, j]], x] === 1 &&"
        "    AllTrue[Range[j + 1, n], H[[#, j]] === 0 &] &&"
        "    AllTrue[Range[j - 1], kDeg[F, H[[#, j]], x] < kDeg[F, H[[j, j]], x] &]]]];"
        /* row reduces to zero modulo the upper-triangular H: same module */
        "pmInModuleQ[F_, H_, row0_] := Module[{row = kc[F, #] & /@ row0, q},"
        "  Do[If[row[[j]] =!= 0,"
        "      q = kDivRem[F, row[[j]], H[[j, j]], x];"
        "      If[q[[2]] =!= 0, Return[False, Module]];"
        "      row = MapThread[kc[F, #1 - q[[1]] #2] &, {row, H[[j]]}]],"
        "    {j, Length[H]}];"
        "  AllTrue[row, # === 0 &]];"
        "pmSameModuleQ[F_, H_, rows_] := pmHermiteQ[F, H] && AllTrue[rows, pmInModuleQ[F, H, #] &];"
        "True");
    free(r);
}

/* ---------------------------------------------------- Hermite normal form */
static void test_hermite_random(void) {
    /* For every field, 12 random matrices of several shapes: the native form
     * satisfies the invariants and equals the package's hnf. */
    TR_TRUE("SeedRandom[11];"
            "AllTrue[pmRefFields, Function[F, AllTrue[Range[6], Function[k,"
            "  Module[{r = {3, 4, 5, 6, 4, 3}[[Mod[k, 6] + 1]], c, M, H, ref},"
            "    c = Min[r, {2, 3, 3, 4, 2, 1}[[Mod[k, 6] + 1]]];"
            "    M = pmMat[F, r, c, 3];"
            "    H = FLINT`PolynomialHermiteForm[M, x, Sequence @@ nativeOpts[F]];"
            "    ref = pmRef[hnf[F, M, x]];"
            "    If[Head[ref] === tragerFail, Head[H] === FLINT`PolynomialHermiteForm,"
            "      H === ref && pmSameModuleQ[F, H, M]]]]]]]");
    /* the sextic field: invariants only (same module, Hermite shape) */
    TR_TRUE("SeedRandom[13];"
            "AllTrue[Range[4], Function[k, Module[{M = pmMat[pmNF6, 4, 2, 2], H},"
            "  H = FLINT`PolynomialHermiteForm[M, x, Sequence @@ nativeOpts[pmNF6]];"
            "  Head[H] === FLINT`PolynomialHermiteForm || pmSameModuleQ[pmNF6, H, M]]]]");
    /* the determinant of the module is the monic product of the pivots */
    TR_TRUE("SeedRandom[12];"
            "AllTrue[{pmQ, fieldGF[101], pmNF3}, Function[F, AllTrue[Range[6], Function[k,"
            "  Module[{M = pmMat[F, 3, 3, 3], H},"
            "    H = FLINT`PolynomialHermiteForm[M, x, Sequence @@ nativeOpts[F]];"
            "    kc[F, Times @@ Diagonal[H]] === kMonic[F, kc[F, Det[M]], x]]]]]]");
    /* worked example over Q and Q(Sqrt[2]) */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x^2 + 1, x}, {x, 1 + x}}, x]", "{{1, -x^2}, {0, 1 + x + x^3}}");
    TR_EQ("FLINT`PolynomialHermiteForm[{{x^2 + t, t x}, {x, 1 + t}}, x, Extension -> {t, t^2 - 2}]",
          "{{1, -1/2 t x}, {0, 2 + t + x^2}}");
    TR_EQ("FLINT`PolynomialHermiteForm[{{2 x + 4, 6}, {3, x}}, x, Modulus -> 5]",
          "{{1, 2 x}, {0, 1 + 2 x + x^2}}");
}

static void test_hermite_moddet(void) {
    /* "DeterminantMultiple" -> D (a multiple of det): same form as plain hnf */
    TR_TRUE("SeedRandom[21];"
            "AllTrue[pmRefFields, Function[F, AllTrue[Range[4], Function[k,"
            "  Module[{M = pmMat[F, 4, 3, 3], ref, D, H},"
            "    ref = pmRef[hnf[F, M, x]];"
            "    If[Head[ref] === tragerFail, True,"
            "      D = kc[F, (Times @@ Diagonal[ref]) pmPoly[F, 2]];"
            "      If[D === 0, True,"
            "        H = FLINT`PolynomialHermiteForm[M, x, \"DeterminantMultiple\" -> D, Sequence @@ nativeOpts[F]];"
            "        H === ref && H === pmRef[hnfModDet[F, M, x, D]]]]]]]]]");
    /* deg D <= 0: ignored */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, \"DeterminantMultiple\" -> 3]",
          "{{1, x}, {0, -1 + x^2}}");
}

static void test_hermite_denominator(void) {
    /* "Denominator": common prime-field factors of den and every entry go */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x (x - 1), 0}, {0, x (x + 1)}}, x, \"Denominator\" -> x^2]",
          "{{{-1 + x, 0}, {0, 1 + x}}, x}");
    TR_EQ("FLINT`PolynomialHermiteForm[{{x^2, 0}, {0, x^2}}, x, \"Denominator\" -> x^3]",
          "{{{1, 0}, {0, 1}}, x}");                        /* a square: two passes */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x + 1, 0}, {0, 1}}, x, \"Denominator\" -> 2 x + 2]",
          "{{{1 + x, 0}, {0, 1}}, 1 + x}");                  /* no common factor; den made monic */
    /* over Q(Sqrt[2]): x^2 - 2 = (x - t)(x + t) splits, but only a factor over
     * Q may be removed -- x - t divides every entry, x^2 - 2 does not */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x - t, 0}, {0, x - t}}, x, \"Denominator\" -> x^2 - 2,"
          " Extension -> {t, t^2 - 2}]",
          "{{{-t + x, 0}, {0, -t + x}}, -2 + x^2}");
    TR_EQ("FLINT`PolynomialHermiteForm[{{x^2 - 2, 0}, {0, t x^2 - 2 t}}, x, \"Denominator\" -> x^3 - 2 x,"
          " Extension -> {t, t^2 - 2}]",
          "{{{1, 0}, {0, 1}}, x}");
    /* agreement with the package's idealReduce over every field */
    TR_TRUE("SeedRandom[31];"
            "AllTrue[pmFields, Function[F, AllTrue[Range[6], Function[k,"
            "  Module[{g = kMonic[F, Expand[x^2 + RandomInteger[{1, 5}] x + 1] /. t -> 0, x],"
            "          M0 = pmMat[F, 3, 3, 2], M, den, ref, nat},"
            "    ref = pmRef[hnf[F, M0, x]];"
            "    If[Head[ref] === tragerFail || kDeg[F, g, x] < 1, True,"
            "      M = Map[kc[F, # g] &, ref, {2}];"
            "      den = kc[F, g^2 (x + 7)];"
            "      nat = FLINT`PolynomialHermiteForm[M, x, \"Denominator\" -> den, Sequence @@ nativeOpts[F]];"
            "      ref = pmRef[idealReduce[F, <|\"x\" -> x|>, <|\"M\" -> M, \"den\" -> den|>]];"
            "      nat === {ref[\"M\"], ref[\"den\"]}]]]]]]");
}

static void test_hermite_products(void) {
    /* "Products" -> {gens, table}: ideal products on real curve bases agree with
     * the package's idealTimes / idealTimesGenerators (products formed in .m) */
    TR_TRUE("Module[{ok = True},"
            "  Do[Module[{B = cb[[1]], F = cb[[2]], I1, I2, nat, ref},"
            "    I1 = residueIdeal[F, B, 1, 1, kc[F, D[B[\"G\"], x]], kc[F, B[\"G\"] + 7]];"
            "    I2 = idealFromGenerators[F, B, {PadRight[{kc[F, x - 3]}, B[\"m\"]], PadRight[{0, 1}, B[\"m\"]]}];"
            "    nat = idealTimes[F, B, I1, I2];"
            "    ref = pmRef[idealTimes[F, B, I1, I2]];"
            "    ok = ok && nat === ref],"
            "   {cb, {{curveBasis[2, x^3 + 1, x], fieldQ[]}, {curveBasis[3, x^2 + 1, x], fieldQ[]},"
            "         {basisMod[curveBasis[4, x^3 - 2 x, x], 101], fieldGF[101]},"
            "         {curveBasis[3, x^2 (x - 1), x], fieldQ[]}}}];"
            "  ok]");
}

/* -------------------------------------------------- determinant / adjugate */
static void test_det_adjugate(void) {
    TR_TRUE("SeedRandom[41];"
            "AllTrue[pmFields, Function[F, AllTrue[Range[8], Function[k,"
            "  Module[{n = 1 + Mod[k, 4], M, r},"
            "    M = pmMat[F, n, n, 3];"
            "    r = FLINT`PolynomialDetAdjugate[M, x, Sequence @@ nativeOpts[F]];"
            "    r[[1]] === kc[F, Det[M]] &&"
            "    Map[kc[F, #] &, r[[2]] . M, {2}] === Map[kc[F, #] &, r[[1]] IdentityMatrix[n], {2}] &&"
            "    r[[2]] === pmRef[polyAdjugate[F, M, x]]]]]]]");
    TR_EQ("FLINT`PolynomialDetAdjugate[{{x^2 + 1, x}, {x, 1 + x}}, x]",
          "{1 + x + x^3, {{1 + x, -x}, {-x, 1 + x^2}}}");
    TR_EQ("FLINT`PolynomialDetAdjugate[{{x + 3}}, x, Modulus -> 5]", "{3 + x, {{1}}}");
    /* "Products": the multiplication matrix of h, as elementInverse uses it */
    TR_TRUE("Module[{B = curveBasis[3, x^2 + 1, x], h = {x, 1, x + 2}, Mh},"
            "  Mh = Table[afTimesPoly[fieldQ[], B, afUnit[B, k], h], {k, 0, 2}];"
            "  FLINT`PolynomialDetAdjugate[IdentityMatrix[3], x, \"Products\" -> {{h}, B[\"mult\"]}] ==="
            "  FLINT`PolynomialDetAdjugate[Mh, x]]");
    TR_TRUE("Module[{B = curveBasis[2, x^3 - 2, x], h = {x - 3, 1}},"
            "  elementInverse[fieldQ[], B, h] === pmRef[elementInverse[fieldQ[], B, h]]]");
}

/* ------------------------------------------------------- weak Popov form */
static void test_weak_popov(void) {
    /* the package's reducedRows on frame bases, native against reference */
    TR_TRUE("SeedRandom[51];"
            "Module[{ok = True},"
            "  Do[Module[{Bt = curveBasisOver[F, m, G, x], rows, nat, ref},"
            "    Do["
            "      rows = Table[pmPoly[F, RandomInteger[{0, 3}]], {m}, {m}];"
            "      If[kc[F, Det[rows]] =!= 0,"
            "        nat = reducedRows[F, Bt, rows];"
            "        ref = pmRef[reducedRows[F, Bt, rows]];"
            "        ok = ok && nat === ref && DuplicateFreeQ[leadingColumn[F, Bt, #] & /@ nat]],"
            "      {6}]],"
            "   {F, {fieldQ[], fieldGF[101]}}, {mG, {{2, x^3 + 2}, {3, x^2 + 1}, {4, x^3 - 2 x}}},"
            "   {m, {mG[[1]]}}, {G, {mG[[2]]}}];"
            "  ok]");
    /* shifted degrees: leading columns distinct afterwards */
    TR_EQ("FLINT`PolynomialWeakPopovForm[{{x^2, x}, {x, 1}}, x, {0, 0}]", "{{0, 0}, {x, 1}}");
    /* det 1: (1, 0) after one step, then (x, 1) - x (1, 0) = (0, 1) */
    TR_EQ("FLINT`PolynomialWeakPopovForm[{{x^2 + 1, x}, {x, 1}}, x, {0, 0}]", "{{1, 0}, {0, 1}}");
    TR_EQ("FLINT`PolynomialWeakPopovForm[{{x, 1}, {1, x}}, x, {0, 5}, \"DegreeScale\" -> 2]",
          "{{x, 1}, {1 - x^2, 0}}");
}

/* ------------------------------------------------ division, extended gcd */
static void test_division_gcd(void) {
    TR_TRUE("SeedRandom[61];"
            "AllTrue[pmFields, Function[F, AllTrue[Range[10], Function[k,"
            "  Module[{a = pmPoly[F, 6], b = pmPoly[F, 3], qr, eg},"
            "    If[b === 0, True,"
            "      qr = FLINT`PolynomialQuotientRemainder[a, b, x, Sequence @@ nativeOpts[F]];"
            "      eg = FLINT`PolynomialExtendedGCD[a, b, x, Sequence @@ nativeOpts[F]];"
            "      kc[F, qr[[1]] b + qr[[2]]] === kc[F, a] && kDeg[F, qr[[2]], x] < kDeg[F, b, x] &&"
            "      kc[F, eg[[2, 1]] a + eg[[2, 2]] b] === eg[[1]] && kLC[F, eg[[1]], x] === 1 &&"
            "      kRem[F, a, eg[[1]], x] === 0 && kRem[F, b, eg[[1]], x] === 0 &&"
            "      If[F[\"Type\"] === \"NF\","
            "        qr === pmRef[kDivRemNF[F, a, b, x]] &&"
            "        eg[[1]] === pmRef[kExtGCDEuclid[F, a, b, x]][[1]], True]]]]]]]");
    TR_EQ("FLINT`PolynomialExtendedGCD[x^2 - 2, x - t, x, Extension -> {t, t^2 - 2}]", "{-t + x, {0, 1}}");
    TR_EQ("FLINT`PolynomialExtendedGCD[0, 0, x]", "{0, {0, 0}}");
    /* the modular gcd over Q(theta): equals the Euclidean one, including a
     * non-integral mu (theta scaled to an algebraic integer internally) */
    TR_TRUE("SeedRandom[62];"
            "AllTrue[{pmNF2, pmNF3, pmNF6, fieldNF[t, t^2 - 1/3, None], fieldNF[t, t^3 - t/4 - 1/8, None]},"
            " Function[F, AllTrue[Range[5], Function[k,"
            "  Module[{g = kMonic[F, pmPoly[F, 2] + x^3, x], a, b, r},"
            "    a = kc[F, g pmPoly[F, 3]]; b = kc[F, g pmPoly[F, 2]];"
            "    r = FLINT`PolynomialFieldGCD[a, b, x, Sequence @@ nativeOpts[F]];"
            "    r === pmRef[First[kExtGCD[F, a, b, x]]] && kRem[F, r, g, x] === 0]]]]]");
    TR_EQ("FLINT`PolynomialFieldGCD[x^2 - 2, x^2 - t x, x, Extension -> {t, t^2 - 2}]", "-t + x");
    TR_EQ("FLINT`PolynomialFieldGCD[x^2 - 1, x^2 + 2 x + 1, x, Modulus -> 7]", "1 + x");
    TR_EQ("FLINT`PolynomialFieldGCD[x^2 - 1, 0, x]", "-1 + x^2");
    TR_TRUE("Module[{qr = FLINT`PolynomialQuotientRemainder[x^3 + 1, 2 x + 1, x, Modulus -> 7]},"
            " gfMod[qr[[1]] (2 x + 1) + qr[[2]] - (x^3 + 1), 7] === 0 && Exponent[qr[[2]], x] <= 0]");
    /* t of degree >= deg mu is reduced on input */
    TR_EQ("FLINT`PolynomialQuotientRemainder[t^2 x, 1, x, Extension -> {t, t^2 - 2}]", "{2 x, 0}");
    /* rationals over GF(p): 1/2 = 4 mod 7 */
    TR_EQ("FLINT`PolynomialQuotientRemainder[x/2, 1, x, Modulus -> 7]", "{4 x, 0}");
}

/* ------------------------------------------------------------ out of scope */
static void test_out_of_scope(void) {
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {2 x, 2}}, x]] === FLINT`PolynomialHermiteForm");   /* rank 1 */
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, y}, {1, x}}, x]] === FLINT`PolynomialHermiteForm");     /* symbol y */
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{Sin[x], 1}, {1, x}}, x]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, Modulus -> 8]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, Extension -> {t, t^2 - 1}]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, Modulus -> 7, Extension -> {t, t^2 - 2}]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, \"Bogus\" -> 1]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1}}, x]] === FLINT`PolynomialHermiteForm");          /* ragged */
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{}, x]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, 3]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{1/x, 1}, {1, x}}, x]] === FLINT`PolynomialHermiteForm");
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, \"Denominator\" -> x - t,"
            " Extension -> {t, t^2 - 2}]] === FLINT`PolynomialHermiteForm");                                   /* den not over Q */
    /* a constant of K is a unit of K[x]: M/(t x) is M/x */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, \"Denominator\" -> t x, Extension -> {t, t^2 - 2}]",
          "{{{1, x}, {0, -1 + x^2}}, x}");
    TR_TRUE("Head[FLINT`PolynomialDetAdjugate[{{x, 1}}, x]] === FLINT`PolynomialDetAdjugate");               /* not square */
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[x, 0, x]] === FLINT`PolynomialQuotientRemainder");
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[x/2, 1, x, Modulus -> 2]] === FLINT`PolynomialQuotientRemainder");
    TR_TRUE("Head[FLINT`PolynomialExtendedGCD[x, 1.5, x]] === FLINT`PolynomialExtendedGCD");
    TR_TRUE("MemberQ[Attributes[FLINT`PolynomialHermiteForm], Protected]");
    TR_TRUE("AllTrue[{FLINT`PolynomialDetAdjugate, FLINT`PolynomialWeakPopovForm, FLINT`PolynomialQuotientRemainder,"
            " FLINT`PolynomialExtendedGCD, FLINT`PolynomialFieldGCD}, MemberQ[Attributes[#], Protected] &]");
    /* large word-size prime */
    TR_EQ("FLINT`PolynomialHermiteForm[{{x, 1}, {1, x}}, x, Modulus -> 2305843009213693951]",
          "{{1, x}, {0, 2305843009213693950 + x^2}}");
}

#include <time.h>
/* per-section wall time on stderr, to keep the suite on the push tier */
#define PM_TIMED(fn) do { clock_t c0 = clock(); fn(); \
    fprintf(stderr, "  %-26s %6.2f s\n", #fn, (double)(clock() - c0) / CLOCKS_PER_SEC); } while (0)

int main(void) {
    tr_init();
    setup();
    PM_TIMED(test_hermite_random);
    PM_TIMED(test_hermite_moddet);
    PM_TIMED(test_hermite_denominator);
    PM_TIMED(test_hermite_products);
    PM_TIMED(test_det_adjugate);
    PM_TIMED(test_weak_popov);
    PM_TIMED(test_division_gcd);
    PM_TIMED(test_out_of_scope);
    printf("flint_polymat: %d assertions passed\n", tr_assertions);
    return 0;
}

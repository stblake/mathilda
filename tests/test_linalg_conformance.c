/* test_linalg_conformance.c -- dense linear-algebra conformance against the
 * reference, one section per defect the differential oracle
 * (tools/sparse_linalg_oracle.py) found in the dense paths.
 *
 * Expected values come from the reference system (wolframscript) wherever the
 * answer is unique; a backend-dependent numeric answer is checked through an
 * invariant (residual, orthogonality) instead.  Failures are counted rather
 * than asserted (Release builds pass -DNDEBUG). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* `input` must emit a message with this tag (Check returns "msg") */
static void chk_msg(const char* input, const char* tag) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Quiet[Check[%s; \"none\", \"msg\", %s]]", input, tag);
    chk(buf, "\"msg\"");
}

/* `input` must emit no message at all */
static void chk_nomsg(const char* input) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Quiet[Check[%s; \"none\", \"msg\"]]", input);
    chk(buf, "\"none\"");
}

static void setup(void) {
    chk("prec[e_] := Union[Flatten[{Precision /@ Flatten[{e}]}]]; "
        "resid[a_, x_, b_] := Max[Abs[Flatten[{a . x - b}]]]; "
        "Null", "Null");
}

/* ---- 1. MPFR SVD / LeastSquares of a wide matrix ------------------------- */

static void test_mpfr_svd_wide(void) {
    /* used to overrun U (one Jacobi column per input column, p > n) */
    chk_true("a = N[{{1, 2, 3}, {4, 5, 6}}, 30]; {u, w, v} = SingularValueDecomposition[a]; "
             "Dimensions /@ {u, w, v} == {{2, 2}, {2, 3}, {3, 3}}");
    chk_true("Max[Abs[Flatten[u . w . ConjugateTranspose[v] - a]]] < 10^-28");
    chk_true("Max[Abs[Flatten[ConjugateTranspose[v] . v - IdentityMatrix[3]]]] < 10^-28");
    chk_true("prec[{u, v}] == {prec[a][[1]]}");
    /* reference: singular values 9.5080320006957241865, 0.7728696356734842916 */
    chk_true("Abs[w[[1, 1]] - 95080320006957241865/10^19] < 10^-18 && "
             "Abs[w[[2, 2]] - 7728696356734842916/10^19] < 10^-18");
    /* reference U (sign included) */
    chk_true("Abs[u[[1, 1]] - 3863177031186115398/10^19] < 10^-18 && "
             "Abs[u[[2, 1]] - 9223657800770583205/10^19] < 10^-18");
    /* complex wide */
    chk_true("c = N[{{1 + I, 2, 3}, {4, 5 - 2 I, 6}}, 30]; {u, w, v} = SingularValueDecomposition[c]; "
             "Max[Abs[Flatten[u . w . ConjugateTranspose[v] - c]]] < 10^-28");
    chk_true("Abs[w[[1, 1]] - 97292513889380419778/10^19] < 10^-18 && "
             "Abs[w[[2, 2]] - 11583036781635380033/10^19] < 10^-18");
    /* truncated form */
    chk("Dimensions /@ SingularValueDecomposition[N[{{1, 2, 3}, {4, 5, 6}}, 30], 1]",
        "{{2, 1}, {1, 1}, {3, 1}}");
    /* LeastSquares through the pseudo-inverse: reference {-1/18, 1/9, 5/18} */
    chk_true("ls = LeastSquares[N[{{1, 2, 3}, {4, 5, 6}}, 30], N[{1, 2}, 30]]; "
             "Max[Abs[ls - {-1/18, 1/9, 5/18}]] < 10^-27");
    chk_true("prec[ls][[1]] > 29");
}

/* ---- 4. SchurDecomposition refuses exact input --------------------------- */

static void test_schur_exact(void) {
    chk_msg("SchurDecomposition[{{1, 2}, {3, 4}}]", "SchurDecomposition::schurf");
    chk("Quiet[SchurDecomposition[{{1, 2}, {3, 4}}]]", "SchurDecomposition[{{1, 2}, {3, 4}}]");
    chk_msg("SchurDecomposition[{{1/2, 2}, {3, 4}}]", "SchurDecomposition::schurf");
    chk_msg("SchurDecomposition[{{1 + I, 2}, {3, 4}}]", "SchurDecomposition::schurf");
    chk_msg("SchurDecomposition[{{1, 2}, {3, Pi}}]", "SchurDecomposition::schurf");
    chk_msg("SchurDecomposition[{{x, 2}, {3, 4}}]", "SchurDecomposition::schurn");
    chk_msg("SchurDecomposition[{{x, 2.}, {3, 4}}]", "SchurDecomposition::schurn");
    chk_msg("SchurDecomposition[{{1, 2}, {3, 4}}, Pivoting -> True]", "SchurDecomposition::schurf");
    chk_msg("SchurDecomposition[{{1, 2, 3}, {4, 5, 6}}]", "SchurDecomposition::matsq");
    chk_msg("SchurDecomposition[{{{1, 2}, {3, 4}}, {{1, 0}, {0, 1}}}]", "SchurDecomposition::exnum");
    chk_msg("SchurDecomposition[{{{1., 2}, {3, 4}}, {{x, 0}, {0, 1}}}]", "SchurDecomposition::exnum");
    /* one inexact entry makes the matrix numeric (contagion) */
    chk("Length[SchurDecomposition[{{1, 2}, {3, 4.}}]]", "2");
    chk("Length[SchurDecomposition[{{{1., 2}, {3, 4}}, {{1, 0}, {0, 1}}}]]", "4");
    /* reference: T[[1, 1]] = -0.3722813232690143 */
    chk_true("{q, t} = SchurDecomposition[{{1., 2}, {3, 4}}]; Abs[t[[1, 1]] + 0.3722813232690143] < 10^-12");
}

/* ---- 7. ::luc on a badly conditioned machine matrix ------------------------ */

static void test_luc(void) {
    /* threshold: rcond < 10^6 eps (2.215e-10 warns, 2.23e-10 does not) */
    chk_nomsg("Inverse[{{1., 0.}, {0., 2.23*^-10}}]");
    chk_msg("Inverse[{{1., 0.}, {0., 2.215*^-10}}]", "Inverse::luc");
    chk_nomsg("Inverse[{{1., 1.}, {1., 1. + 1*^-9}}]");
    chk_msg("Inverse[{{1., 1.}, {1., 1. + 8*^-10}}]", "Inverse::luc");
    chk_msg("Det[{{1., 1.}, {1., 1. + 8*^-10}}]", "Det::luc");
    chk_msg("LinearSolve[{{1., 1.}, {1., 1. + 8*^-10}}, {1., 2.}]", "LinearSolve::luc");
    chk_msg("LUDecomposition[{{1., 1.}, {1., 1. + 8*^-10}}]", "LUDecomposition::luc");
    chk_nomsg("Det[{{1., 1.}, {1., 1. + 1*^-9}}]");
    chk_nomsg("LinearSolve[{{1., 1.}, {1., 1. + 1*^-9}}, {1., 2.}]");
    chk_nomsg("LUDecomposition[{{1., 1.}, {1., 1. + 1*^-9}}]");
    /* Inverse uses the 1-norm rcond, the others the infinity norm */
    chk_msg("Inverse[{{1., 0, 1}, {0, 1, 1}, {0, 0, 6*^-10}}]", "Inverse::luc");
    chk_nomsg("Det[{{1., 0, 1}, {0, 1, 1}, {0, 0, 6*^-10}}]");
    chk_nomsg("Inverse[{{1., 0, 0}, {0, 1, 0}, {1, 1, 6*^-10}}]");
    chk_msg("Det[{{1., 0, 0}, {0, 1, 0}, {1, 1, 6*^-10}}]", "Det::luc");
    chk_msg("LinearSolve[{{1., 0, 0}, {0, 1, 0}, {1, 1, 6*^-10}}, {1., 1, 1}]", "LinearSolve::luc");
    chk_msg("LUDecomposition[{{1., 0, 0}, {0, 1, 0}, {1, 1, 6*^-10}}]", "LUDecomposition::luc");
    /* complex */
    chk_msg("Inverse[{{1. + I, 0.}, {0., 1*^-12}}]", "Inverse::luc");
    chk_msg("Det[{{1. + I, 0.}, {0., 1*^-12}}]", "Det::luc");
    /* every call warns (no once-per-process latch) */
    chk_msg("LUDecomposition[{{1., 1.}, {1., 1. + 8*^-10}}]", "LUDecomposition::luc");
    /* an exactly singular matrix: Det answers 0. without a warning */
    chk("Quiet[Check[Det[{{1., 1.}, {1., 1.}}], \"msg\"]]", "0.0");
    /* MatrixPower[m, -k]: Inverse::luc, and MatrixPower::sing below 100 eps */
    chk_msg("MatrixPower[{{1., 0.}, {0., 2.3*^-14}}, -1]", "Inverse::luc");
    chk_msg("MatrixPower[{{1., 0.}, {0., 2.22*^-14}}, -1]", "MatrixPower::sing");
    chk("Quiet[MatrixPower[{{1., 0.}, {0., 2.22*^-14}}, -1]]",
        "MatrixPower[{{1.0, 0.0}, {0.0, 2.22e-14}}, -1]");
    chk_msg("MatrixPower[{{1, 2}, {2, 4}}, -1]", "MatrixPower::sing");
    chk_msg("MatrixPower[{{1., 2}, {2, 4}}, -2]", "MatrixPower::sing");
    chk_nomsg("MatrixPower[{{1., 0.}, {0., 1*^-9}}, -1]");
}

/* ---- 2. arbitrary precision kept through the Inverse family ------------- */

static void test_mpfr_precision(void) {
    /* reference: Inverse[{{N[9, 30]}}] is {{0.111...`30.}} */
    chk_true("prec[Inverse[{{N[9, 30]}}]] == prec[N[9, 30]]");
    chk_true("m = N[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}, 30]; c = N[{{2 + I, 1, 0}, {1, 3 - I, 1}, {0, 1, 4}}, 30]; "
             "p30 = prec[m][[1]]; "
             "And @@ (prec[#] == {p30} &) /@ {Inverse[m], Inverse[c], MatrixPower[m, -1], "
             "MatrixPower[c, -2], LinearSolve[m, N[{1, 2, 3}, 30]], LinearSolve[c, N[{1, 2, 3}, 30]], "
             "Det[m], Det[c], PseudoInverse[m], First[LUDecomposition[c]], Flatten[QRDecomposition[m]], "
             "Flatten[QRDecomposition[c]], Flatten[SchurDecomposition[m]], Flatten[SchurDecomposition[c]], "
             "Flatten[SchurDecomposition[m, RealBlockDiagonalForm -> False]]}");
    chk_true("{q, t} = SchurDecomposition[c]; Max[Abs[Flatten[q . t . ConjugateTranspose[q] - c]]] < 10^-28 && "
             "Max[Abs[Flatten[q . ConjugateTranspose[q] - IdentityMatrix[3]]]] < 10^-28 && "
             "Max[Abs[Flatten[LowerTriangularize[t, -1]]]] == 0");
    chk_true("Max[Abs[Flatten[Inverse[c] . c - IdentityMatrix[3]]]] < 10^-28");
    chk_true("Max[Abs[Flatten[MatrixPower[c, -2] . c . c - IdentityMatrix[3]]]] < 10^-27");
    /* singular: Inverse::sing, MatrixPower::sing */
    chk_msg("Inverse[N[{{1, 2}, {2, 4}}, 30]]", "Inverse::sing");
    chk_msg("MatrixPower[N[{{1, 2}, {2, 4}}, 30], -2]", "MatrixPower::sing");
    /* a 54 x 54 complex Hermitian system solves at the input precision, fast */
    chk_true("h = Table[N[If[i == j, 60 + i, (Mod[i j, 7] - 3) + I (Mod[i + j, 5] - 2) Sign[j - i]], 30], {i, 54}, {j, 54}]; "
             "bb = N[Range[54], 30]; {tt, xx} = AbsoluteTiming[LinearSolve[h, bb]]; "
             "tt < 2 && Max[Abs[h . xx - bb]] < 10^-25 && prec[xx] == {prec[h][[1]]}");
}

/* ---- 3. RowReduce / NullSpace / LinearSolve of inexact singular matrices --- */

static void test_numeric_rref(void) {
    /* reference forms: exact 1 at the pivots, exact 0 below and in zero rows,
     * an inexact 0. above a machine pivot, exact 0 above an MPFR one */
    chk("RowReduce[{{1., 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "{{1, 0.0, -1.0}, {0, 1, 2.0}, {0, 0, 0}}");
    chk("RowReduce[{{1., 2}, {3, 4}}]", "{{1, 0.0}, {0, 1}}");
    chk("RowReduce[N[{{1, 2}, {3, 4}}, 30]]", "{{1, 0}, {0, 1}}");
    chk("RowReduce[N[{{1 + I, 2}, {3, 4}}, 30]]", "{{1, 0}, {0, 1}}");
    chk("RowReduce[{{1. + I, 2}, {3, 4}}]", "{{1, 0.0}, {0, 1}}");
    /* used to be the identity (the 10^-30 residue taken as a pivot) */
    chk_true("r = RowReduce[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 30]]; "
             "r[[3]] === {0, 0, 0} && r[[1, 1]] === 1 && r[[2, 2]] === 1 && "
             "Abs[r[[1, 3]] + 1] < 10^-28 && Abs[r[[2, 3]] - 2] < 10^-28");
    chk("RowReduce[N[{{1 + I, 2}, {2 + 2 I, 4}}, 30]][[2]]", "{0, 0}");
    chk_true("r = RowReduce[N[{{1 + I, 2}, {2 + 2 I, 4}}, 30]]; Abs[r[[1, 2]] - (1 - I)] < 10^-28");
    /* the machine zero tolerance: max(r, c) eps max|a| (reference) */
    chk("{MatrixRank[{{1., 0.}, {0., 4*^-16}}], RowReduce[{{1., 0.}, {0., 4*^-16}}][[2, 2]], "
        "RowReduce[{{1., 0.}, {0., 4.5*^-16}}][[2, 2]], RowReduce[{{1., 0, 0}, {0, 1., 0}, {0, 0, 6*^-16}}][[3, 3]], "
        "RowReduce[{{1., 0, 0}, {0, 1., 0}, {0, 0, 7*^-16}}][[3, 3]]}", "{1, 0, 1, 0, 1}");
    /* RowReduce::luc: a pivot below 10^6 eps of the largest entry */
    chk_msg("RowReduce[{{1., 0.}, {0., 1*^-10}}]", "RowReduce::luc");
    chk_nomsg("RowReduce[{{1., 0.}, {0., 3*^-10}}]");
    chk_msg("RowReduce[{{1., 1., 1.}, {1., 1. + 1*^-12, 2.}}]", "RowReduce::luc");
    chk_nomsg("RowReduce[{{1., 1.}, {1., 1. + 1*^-12}, {1, 2}}]");
    chk_nomsg("RowReduce[{{1., 1.}, {1., 1. + 3*^-16}}]");
    /* a singular machine complex 11 x 11: the free columns come back */
    chk_true("cm = Table[N[(Mod[i j, 5] - 2) + I (Mod[i + 2 j, 3] - 1)], {i, 11}, {j, 11}]; "
             "cm[[11]] = cm[[1]] + (1 + I) cm[[2]]; rr = RowReduce[cm]; "
             "rr[[11]] === Table[0, {11}] && Max[Abs[Flatten[rr . NullSpace[cm][[1]] ]]] < 10^-10");
    /* NullSpace at arbitrary precision: orthonormal, annihilated */
    chk_true("ns = NullSpace[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 30]]; Length[ns] == 1 && "
             "Abs[Norm[ns[[1]]] - 1] < 10^-28 && "
             "Max[Abs[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 30] . ns[[1]]]] < 10^-28");
    /* reference: {{0.4082482904638630164, -0.8164965809277260327, 0.4082482904638630164}} */
    chk_true("Max[Abs[ns[[1]] - {4082482904638630164, -8164965809277260327, 4082482904638630164}/10^19]] < 10^-18");
    chk_true("ns = NullSpace[N[{{1, 2, 3}, {2, 4, 6}, {3, 6, 9}}, 30]]; Length[ns] == 2 && "
             "Max[Abs[Flatten[ns . ConjugateTranspose[ns] - IdentityMatrix[2]]]] < 10^-28");
    chk_true("ns = NullSpace[N[{{1 + I, 2}, {2 + 2 I, 4}}, 30]]; Length[ns] == 1 && "
             "Max[Abs[N[{{1 + I, 2}, {2 + 2 I, 4}}, 30] . ns[[1]]]] < 10^-28");
    /* LinearSolve at arbitrary precision: every free variable 0 (reference
     * {1., 0``30}, {-1/3, 2/3, 0``30}, {0``30, 0.5}) */
    chk_true("xs = LinearSolve[N[{{1, 2}, {2, 4}}, 30], N[{1, 2}, 30]]; Abs[xs[[1]] - 1] < 10^-28 && xs[[2]] == 0");
    chk_true("xs = LinearSolve[N[{{1, 2, 3}, {4, 5, 6}}, 30], N[{1, 2}, 30]]; "
             "Max[Abs[xs - {-1/3, 2/3, 0}]] < 10^-28");
    chk_true("xs = LinearSolve[N[{{1, 2}, {3, 4}, {5, 6}}, 30], N[{1, 2, 3}, 30]]; "
             "Max[Abs[xs - {0, 1/2}]] < 10^-28");
    chk_msg("LinearSolve[N[{{1, 2}, {2, 4}}, 30], N[{1, 3}, 30]]", "LinearSolve::nosol");
    chk_true("mm = N[{{1, 2}, {3, 4}}, 30]; xs = LinearSolve[mm, N[{{1, 0}, {0, 1}}, 30]]; "
             "Max[Abs[Flatten[mm . xs - IdentityMatrix[2]]]] < 10^-28");
}

/* ---- 5. QRDecomposition: exact radicals fast, machine residual exact ------ */

static void test_qr(void) {
    /* reference: {{{1/Sqrt[10], 3/Sqrt[10]}, {3/Sqrt[10], -(1/Sqrt[10])}},
     *             {{Sqrt[10], 7*Sqrt[2/5]}, {0, Sqrt[2/5]}}} */
    chk("QRDecomposition[{{1, 2}, {3, 4}}]",
        "{{{1/Sqrt[10], 3/Sqrt[10]}, {3/Sqrt[10], -1/Sqrt[10]}}, {{Sqrt[10], 7 Sqrt[2/5]}, {0, Sqrt[2/5]}}}");
    chk_true("m5 = {{1,2,3,4,5},{6,7,8,9,1},{2,3,5,7,11},{1,0,1,0,1},{3,1,4,1,5}}; "
             "{tt, {q5, r5}} = AbsoluteTiming[QRDecomposition[m5]]; tt < 1");
    /* reference column / row, in its radical form */
    chk("q5[[1]]", "{1/Sqrt[51], 2 Sqrt[3/17], 2/Sqrt[51], 1/Sqrt[51], Sqrt[3/17]}");
    chk("r5[[1]]", "{Sqrt[51], 53/Sqrt[51], 74/Sqrt[51], 25 Sqrt[3/17], 49/Sqrt[51]}");
    chk("Simplify[ConjugateTranspose[q5] . r5 - m5]", "{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}");
    chk("Simplify[q5 . ConjugateTranspose[q5]]", "{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}");
    chk_true("m6 = Table[Mod[3 i + 5 j^2, 11] - 5, {i, 6}, {j, 6}]; "
             "{tt, {q6, r6}} = AbsoluteTiming[QRDecomposition[m6]]; tt < 2 && "
             "Simplify[ConjugateTranspose[q6] . r6 - m6] == ConstantArray[0, {6, 6}]");
    /* rank deficient / rectangular: dependent columns dropped */
    chk("QRDecomposition[{{1, 2}, {2, 4}}]", "{{{1/Sqrt[5], 2/Sqrt[5]}}, {{Sqrt[5], 2 Sqrt[5]}}}");
    chk("Dimensions /@ QRDecomposition[{{1, 2, 3}, {4, 5, 6}}]", "{{2, 2}, {2, 3}}");
    chk("Dimensions /@ QRDecomposition[{{1, 2}, {3, 4}, {5, 6}}]", "{{2, 3}, {2, 2}}");
    chk("QRDecomposition[{{1/2, 1}, {0, 1/3}}]", "{{{1, 0}, {0, 1}}, {{1/2, 1}, {0, 1/3}}}");
    /* reference: the exact zero matrix gives one row */
    chk("QRDecomposition[ConstantArray[0, {3, 4}]]", "{{{1, 0, 0}}, {{0, 0, 0, 0}}}");
    chk("QRDecomposition[{{0, 0}, {0, 1}}]", "{{{0, 1}}, {{0, 1}}}");
    /* machine wide: a dependent early column used to drop the later rows of r */
    chk_true("w = Table[N[Mod[7 i + 3 j, 13] - 6 + If[j == 2, 0, 0]], {i, 65}, {j, 75}]; "
             "w[[All, 2]] = 2 w[[All, 1]]; {qw, rw} = QRDecomposition[w]; "
             "Max[Abs[Flatten[ConjugateTranspose[qw] . rw - w]]] < 10^-10");
    chk_true("w = Table[N[Sin[i j]], {i, 65}, {j, 75}]; {qw, rw} = QRDecomposition[w]; "
             "Max[Abs[Flatten[ConjugateTranspose[qw] . rw - w]]] < 10^-10");
}

/* ---- 6. machine complex: mixed Real/Complex matrices take the buffer path -- */

static void test_complex_perf(void) {
    /* Dot of two 166 x 166 matrices mixing machine Reals and Complexes ran the
     * boxed loop for 20 s */
    chk_true("SeedRandom[1]; nn = 166; ma = Table[If[EvenQ[i + j], RandomReal[], RandomComplex[]], {i, nn}, {j, nn}]; "
             "mb = Table[If[OddQ[i + j], RandomReal[], RandomComplex[]], {i, nn}, {j, nn}]; "
             "{tt, pr} = AbsoluteTiming[ma . mb]; tt < 2");
    chk_true("Abs[pr[[7, 9]] - Sum[ma[[7, k]] mb[[k, 9]], {k, nn}]] < 10^-10");
    /* exact beside machine: contagion, the answer is machine */
    chk("{1, 2} . {3., 4.}", "11.0");
    chk("{{1, 2}, {3, 4}} . {1., I}", "{1.0 + 2.0*I, 3.0 + 4.0*I}");
    chk("{{1, 2}, {3, 4}} . {{1, 1/2}, {1, 1}}", "{{3, 5/2}, {7, 11/2}}");
    chk("{p1, p2} . {1., 2}", "1.0 p1 + 2 p2");
    /* the complex factors and their invariants, fast */
    chk_true("c91 = RandomComplex[{-1 - I, 1 + I}, {91, 91}]; {tt, {q, t}} = AbsoluteTiming[SchurDecomposition[c91]]; "
             "{t2, rr} = AbsoluteTiming[q . t . ConjugateTranspose[q]]; tt + t2 < 2 && Max[Abs[Flatten[rr - c91]]] < 10^-12");
    chk_true("h59 = Table[If[i == j, 70. + i, N[(Mod[i j, 7] - 3) + I (Mod[i + j, 5] - 2)]], {i, 59}, {j, 59}]; "
             "h59 = (h59 + ConjugateTranspose[h59])/2; {tt, rr} = AbsoluteTiming[RowReduce[h59]]; "
             "tt < 2 && rr == IdentityMatrix[59]");
    chk_true("{tt, {lu, pp, cc}} = AbsoluteTiming[LUDecomposition[c91]]; tt < 1");
    /* PseudoInverse / LeastSquares of a machine complex matrix ran the
     * rationalised exact pipeline for minutes */
    chk_true("c20 = RandomComplex[{-1 - I, 1 + I}, {20, 20}]; {tt, pc} = AbsoluteTiming[PseudoInverse[c20]]; "
             "tt < 1 && Max[Abs[Flatten[pc . c20 - IdentityMatrix[20]]]] < 10^-10");
    chk_true("w35 = RandomComplex[{-1 - I, 1 + I}, {3, 5}]; pw = PseudoInverse[w35]; "
             "Max[Abs[Flatten[{w35 . pw . w35 - w35, pw . w35 . pw - pw, "
             "ConjugateTranspose[w35 . pw] - w35 . pw, ConjugateTranspose[pw . w35] - pw . w35}]]] < 10^-12");
    chk_true("ls = LeastSquares[{{1. + I, 2}, {3, 4}, {5, 6 I}}, {1, 2, 3}]; "
             "Norm[ls - PseudoInverse[{{1. + I, 2}, {3, 4}, {5, 6 I}}] . {1, 2, 3}] < 10^-12");
    chk_true("Norm[LeastSquares[NDArray[{{1. + I, 2}, {3, 4}, {5, 6 I}}], {1, 2, 3}] - Normal[ls]] < 10^-12");
}

/* ---- 7b. Minors: each minor warns Det::luc, as the reference's Det does ---- */

static void test_minors_luc(void) {
    chk_msg("Minors[{{1., 1., 0}, {1., 1. + 1*^-12, 0}, {0, 0, 1}}]", "Det::luc");
    chk_nomsg("Minors[{{1., 1., 0}, {1., 1. + 1*^-12, 0}, {0, 0, 1}}, 1]");
    chk_nomsg("Minors[{{1., 1.}, {1., 1.}}, 2]");
    chk_nomsg("Minors[{{1., 2, 0}, {3, 4, 0}, {0, 0, 1}}]");
    chk("Minors[{{1., 2, 0}, {3, 4, 0}, {0, 0, 1}}]", "{{-2.0, 0.0, 0.0}, {0.0, 1.0, 2.0}, {0.0, 3.0, 4.0}}");
}

/* ---- 9. matrix functions of normal matrices: no Schur--Parlett blocks ----- */

static void test_matfunc_normal(void) {
    chk_true("SeedRandom[3]; n2 = 200; a2 = RandomReal[{-1, 1}, {n2, n2}]; s2 = (a2 . Transpose[a2])/n2 + IdentityMatrix[n2]; "
             "{tt, fs} = AbsoluteTiming[MatrixFunction[Sin, s2]]; {vals, vecs} = Eigensystem[s2]; "
             "tt < 1 && Max[Abs[Flatten[fs - Transpose[vecs] . DiagonalMatrix[Sin[vals]] . vecs]]] < 10^-12");
    chk_true("c2 = RandomComplex[{-1 - I, 1 + I}, {n2, n2}]; h2 = (c2 . ConjugateTranspose[c2])/n2 + IdentityMatrix[n2]; "
             "{tt, fh} = AbsoluteTiming[MatrixFunction[Sin, h2]]; "
             "tt < 2 && Max[Abs[Flatten[fh - (MatrixExp[I h2] - MatrixExp[-I h2])/(2 I)]]] < 10^-11");
    chk_true("{tt, sq} = AbsoluteTiming[MatrixPower[h2, 1/2]]; tt < 2 && Max[Abs[Flatten[sq . sq - h2]]] < 10^-11");
    chk_true("{tt, sq} = AbsoluteTiming[MatrixPower[s2, 1/2]]; tt < 1 && Max[Abs[Flatten[sq . sq - s2]]] < 10^-11");
    /* arbitrary precision: Parlett on the complex Schur form, not the exact
     * path (a 6 x 6 30-digit matrix ran past any time limit) */
    chk_true("am = N[{{9, 1, 0, 0, 0, 0}, {1, 93/10, -3/2, 0, 0, 4/5}, {0, -3/2, 13/2, 0, 0, 0}, "
             "{0, 0, 0, 8, 0, 0}, {0, 0, 0, 0, 6, 0}, {0, 4/5, 0, 0, 0, 24/5}}, 30]; "
             "{tt, fm} = AbsoluteTiming[MatrixFunction[Sin, am]]; tt < 1 && "
             "Max[Abs[Flatten[fm - (MatrixExp[I am] - MatrixExp[-I am])/(2 I)]]] < 10^-27 && "
             "prec[fm] == prec[am] && Head[fm[[1, 1]]] === Real");
    chk_true("cm3 = N[{{1 + I, 2, 0}, {1, 3 - 2 I, 1}, {0, 1/2, 2}}, 30]; fc = MatrixFunction[Sin, cm3]; "
             "Max[Abs[Flatten[fc - (MatrixExp[I cm3] - MatrixExp[-I cm3])/(2 I)]]] < 10^-27 && prec[fc] == prec[cm3]");
    chk_true("tu = N[{{1, 1, 2}, {0, 2, 1}, {0, 0, 3}}, 30]; "
             "Max[Abs[Flatten[MatrixFunction[Sin, tu] - (MatrixExp[I tu] - MatrixExp[-I tu])/(2 I)]]] < 10^-27");
    /* a defective matrix declines to the exact path: the derivative appears */
    chk_true("jb = MatrixFunction[Sin, N[{{1, 1}, {0, 1}}, 30]]; Abs[jb[[1, 2]] - Cos[1]] < 10^-27");
    chk_true("Max[Abs[Flatten[MatrixFunction[Sin, N[IdentityMatrix[3], 30]] - Sin[1] IdentityMatrix[3]]]] < 10^-28");
    /* a non-normal matrix keeps Schur--Parlett */
    chk_true("{tt, fa} = AbsoluteTiming[MatrixFunction[Sin, a2]]; "
             "tt < 2 && Max[Abs[Flatten[fa - (MatrixExp[I a2] - MatrixExp[-I a2])/(2 I)]]] < 10^-10");
}

/* ---- 8. Norm: the exact and arbitrary-precision 2-norm --------------------- */

static void test_norm(void) {
    /* reference answers */
    chk("Norm[{{4, 1}, {0, 5}}]", "Sqrt[21 + Sqrt[41]]");
    chk_true("Norm[{{1, 2, 3}, {4, 5, 6}}] == Sqrt[(91 + Sqrt[8065])/2]");
    chk_true("Norm[{{1, I}, {2, 3}}] == Sqrt[(15 + Sqrt[173])/2]");
    chk("Norm[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]", "Sqrt[Root[-9 + 244 #1 - 304 #1^2 + #1^3 &, 3]]");
    chk("Norm[{{0, 0}, {0, 0}}]", "0");
    chk("Norm[{{3}}]", "3");
    chk_true("Abs[N[Norm[{{4, 1}, {0, 5}}]] - Norm[{{4., 1.}, {0., 5.}}]] < 10^-12");
    chk_true("Abs[N[Norm[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]] - Norm[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]]] < 10^-12");
    chk_true("Head[Norm[{{p1, p2}, {p3, p4}}]] === Max");
    chk_true("Abs[(Norm[{{p1, p2}, {p3, p4}}] /. {p1 -> 1, p2 -> 2, p3 -> 3, p4 -> 4}) - Norm[N[{{1, 2}, {3, 4}}]]] < 10^-12");
    /* arbitrary precision keeps it (reference 9.5080320006957241865`30.) */
    chk_true("nn = Norm[N[{{1, 2, 3}, {4, 5, 6}}, 30]]; prec[nn][[1]] > 29 && "
             "Abs[nn - 95080320006957241865/10^19] < 10^-18");
    chk_true("Abs[Norm[N[{{1, 2}, {3, 4}, {5, 6}}, 30]] - Norm[N[{{1, 2}, {3, 4}, {5, 6}}]]] < 10^-12");
}

/* ---- 10. ArrayQ[expr, patt, test] ------------------------------------------ */

static void test_arrayq(void) {
    /* every value probed against the reference */
    chk("f = FreeQ[#, _SparseArray] &; "
        "{ArrayQ[{{{1,2},{3,4}}, {1,2}, 3.}], ArrayQ[{{{1,2},{3,4}}, {1,2}, 3.}, _, f], "
        "ArrayQ[{{1,2},{3,4}}, _, f], ArrayQ[{{1,2},{3,4}}, 1, f], ArrayQ[{{1,2},{3,4}}, 2, f], "
        "ArrayQ[{{1,2},{3}}, _, f], ArrayQ[{{1,2},{3}}, 1, f], ArrayQ[{{1,2},{3}}, 2, f]}",
        "{False, True, True, True, True, True, True, False}");
    chk("{ArrayQ[{{1,2},{3,4}}, _, NumberQ], ArrayQ[{{1,2},{3,4}}, 1, NumberQ], "
        "ArrayQ[{{1,2},{3,4}}, _, ListQ], ArrayQ[{{1,2},{3,4}}, 1, ListQ], "
        "ArrayQ[{{1,2},{3}}, _, ListQ], ArrayQ[{{1,2},{3}}, _, NumberQ], "
        "ArrayQ[{{1,2},{3,x}}, 2, NumberQ], ArrayQ[{}, _, NumberQ], ArrayQ[{{},{}}, _, NumberQ], "
        "ArrayQ[3, _, NumberQ], ArrayQ[{{1,2},{3,4}}, 0, ListQ], ArrayQ[{1,{2}}, _, (True&)], "
        "ArrayQ[{{{1}},{{2}}}, 1|2, (True&)]}",
        "{True, False, False, True, True, False, False, True, True, False, False, True, False}");
    chk("{ArrayQ[{{{1}},{{2}}}, 2, ListQ], ArrayQ[{{{1}},{{2}}}, 2], ArrayQ[{{1,2},{3,4}}, 1], "
        "ArrayQ[{{1,2},{3}}, 1], ArrayQ[{{{1}},{{2}}}, 2|3, ListQ], "
        "ArrayQ[{{{1}},{{2}}}, x_ /; x < 3, ListQ], ArrayQ[{{1,2},{3,4}}, 1.0, ListQ], "
        "ArrayQ[{{1,2},{3,4}}, 3, NumberQ], ArrayQ[{{1,2},{3,4}}, -1, NumberQ]}",
        "{True, False, False, False, False, False, False, False, False}");
    chk("{ArrayQ[SparseArray[{{1,2},{3,4}}], 1, ListQ], ArrayQ[SparseArray[{{1,2},{3,4}}], _, NumberQ]}",
        "{True, True}");
    chk("ArrayQ[NDArray[{{1., 2.}, {3., 4.}}], 1, VectorQ]", "True");
}

/* ---- 11. LinearSolve of a singular / rectangular machine system ------------- */

static void test_singular_linearsolve(void) {
    /* reference answers: the minimum-norm solution */
    chk_true("Max[Abs[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 2.}] - {0.2, 0.4}]] < 10^-14");
    chk_true("Max[Abs[LinearSolve[{{1., 2., 3.}, {4., 5., 6.}}, {1., 2.}] - "
             "{-0.05555555555555536, 0.11111111111111109, 0.27777777777777773}]] < 10^-14");
    chk_true("Max[Abs[LinearSolve[{{1., 2., 3.}, {2., 4., 6.}}, {1., 2.}] - "
             "{0.0714285714285714, 0.14285714285714285, 0.21428571428571425}]] < 10^-14");
    chk_true("Max[Abs[LinearSolve[{{1., 2.}, {3., 4.}, {5., 6.}}, {1., 2., 3.}] - {0., 0.5}]] < 10^-14");
    chk_true("Max[Abs[Flatten[LinearSolve[{{1., 2.}, {2., 4.}}, {{1., 1.}, {2., 2.}}] - "
             "{{0.2, 0.2}, {0.4, 0.4}}]]] < 10^-14");
    chk("LinearSolve[{{2., 0., 0.}, {0., 0., 0.}, {0., 0., 1.}}, {1., 0., 1.}]", "{0.5, 0.0, 1.0}");
    chk_true("xs = LinearSolve[{{1. + I, 2.}, {2. + 2 I, 4.}}, {1. + I, 2. + 2. I}]; "
             "resid[{{1. + I, 2.}, {2. + 2 I, 4.}}, xs, {1. + I, 2. + 2. I}] < 10^-14 && "
             "Norm[xs - PseudoInverse[{{1. + I, 2.}, {2. + 2 I, 4.}}] . {1. + I, 2. + 2. I}] < 10^-14");
    /* inconsistent: nosol, unevaluated (the 1e-12 relative residual edge) */
    chk_msg("LinearSolve[{{1., 2.}, {2., 4.}}, {1., 3.}]", "LinearSolve::nosol");
    chk("Quiet[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 3.}]]", "LinearSolve[{{1.0, 2.0}, {2.0, 4.0}}, {1.0, 3.0}]");
    chk_msg("LinearSolve[{{1., 2.}, {3., 4.}, {5., 6.}}, {1., 2., 3.1}]", "LinearSolve::nosol");
    chk("Head[Quiet[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 2. + 5*^-12}]]]", "LinearSolve");
    chk("Head[Quiet[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 2. + 3*^-12}]]]", "List");
    /* the reference refuses a zero matrix even with a zero right-hand side */
    chk_msg("LinearSolve[{{0., 0.}, {0., 0.}}, {0., 0.}]", "LinearSolve::nosol");
    chk("LinearSolve[{{1., 2.}, {2., 4.}}, {0., 0.}]", "{0.0, 0.0}");
    /* packed and NDArray agree with the List */
    chk_true("LinearSolve[NDArray[{{1., 2.}, {2., 4.}}], NDArray[{1., 2.}]] == "
             "NDArray[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 2.}]]");
}

/* A rank-deficient 30-digit real matrix: its zero-eigenvalue cluster makes
 * 2x2 blocks with a near-double eigenvalue, whose split rotation and
 * eigenvalues the MPFR Francis QR computed through tr^2 - 4 det (half the
 * digits lost): the Schur residual was 8e-16 instead of ~1e-28. */
static void test_schur_mpfr_rank_deficient(void) {
    chk_true("SeedRandom[3]; d = N[RandomInteger[{-9, 9}, {12, 12}]/7, 30]; "
             "d[[All, 12]] = d[[All, 1]] + d[[All, 2]]; {q, t} = SchurDecomposition[d]; "
             "Max[Abs[q . t . Transpose[q] - d]] < 10^-25 && Max[Abs[q . Transpose[q] - IdentityMatrix[12]]] < 10^-25");
    chk_true("SeedRandom[5]; d = N[RandomInteger[{-5, 5}, {9, 9}], 50]; d[[All, 9]] = 2 d[[All, 3]]; "
             "d[[All, 8]] = d[[All, 1]] - d[[All, 4]]; {q, t} = SchurDecomposition[d]; "
             "Max[Abs[q . t . Transpose[q] - d]] < 10^-45");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    setup();
    test_mpfr_svd_wide();
    test_mpfr_precision();
    test_numeric_rref();
    test_schur_exact();
    test_qr();
    test_luc();
    test_norm();
    test_complex_perf();
    test_minors_luc();
    test_matfunc_normal();
    test_arrayq();
    test_singular_linearsolve();
    test_schur_mpfr_rank_deficient();
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("linalg_conformance_tests: all passed\n");
    return 0;
}

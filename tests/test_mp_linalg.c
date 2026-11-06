/*
 * tests/test_mp_linalg.c
 *
 * Arbitrary-precision dense linear algebra on ill-conditioned matrices, and
 * the scalar Precision / Accuracy of a matrix.
 *
 *   - LinearSolve[m, b], LinearSolve[m][b] and MatrixRank of N[Hilbert[n], p]:
 *     the condition number of the Hilbert matrix grows like e^(3.5 n), far
 *     beyond 2^53 but below 10^p, so the matrix is regular at p digits and the
 *     solution must agree with the EXACT rational solution (computed by the
 *     exact solver from the exact matrix) to about p - log10(cond) digits;
 *   - LUBackSubstitution against hand-derived solutions and against
 *     LinearSolve, at exact, machine and MPFR precision, real and complex,
 *     vector and matrix right-hand sides, and its messages;
 *   - Precision / Accuracy of lists and matrices are ONE number, the minimum;
 *   - Dot of MPFR vectors and matrices against the exact product;
 *   - Eigenvalues of N[Hilbert[n], p] against the roots of the exact
 *     characteristic polynomial, and of NONSYMMETRIC matrices (a diagonal
 *     similarity of a Hilbert matrix, a matrix with complex eigenvalues, a
 *     numerically low-rank Nystrom kernel matrix) the same way;
 *   - Dot with an exact operand holding zeros (DiagonalMatrix);
 *   - Outer[f, ...] of numeric vectors through compiled code, against the
 *     interpreter's cell-by-cell values.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_check_true(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, "True") != 0) {
        printf("FAIL (expected True): %s\n  got: %.300s\n", input, s);
        ASSERT(0);
    } else {
        printf("PASS: %s\n", input);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

static void run_silent(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(r);
    expr_free(e);
}

static void define_helpers(void) {
    run_silent("hilb[n_] := Table[1/(i + j - 1), {i, n}, {j, n}]");
    /* relative max-norm error of an approximate tensor against an exact one */
    run_silent("relerr[a_, b_] := N[Max[Abs[Flatten[{a - b}]]] / Max[Abs[Flatten[{b}]]], 20]");
}

/* ---- LinearSolve / MatrixRank at arbitrary precision ------------------- */

static void test_linearsolve_hilbert_mpfr(void) {
    /* cond(H_24) ~ 1e35: at 100 digits about 65 correct digits remain */
    run_silent("he = hilb[24]; be = Table[(-1)^k k, {k, 24}]; xe = LinearSolve[he, be]");
    run_check_true("relerr[LinearSolve[N[he, 100], N[be, 100]], xe] < 10^-55");
    run_check_true("relerr[LinearSolve[N[he, 100], be], xe] < 10^-55");
    /* the factor object: no LinearSolve::sing1, same answer, all MPFR */
    run_check_true("Check[f = LinearSolve[N[he, 100]]; True, False]");
    run_check_true("relerr[f[N[be, 100]], xe] < 10^-55");
    run_check_true("relerr[f[be], xe] < 10^-55");
    run_check_true("Abs[Precision[f[N[be, 100]]] - 100] < 1");
    /* a matrix right-hand side: two solution columns */
    run_check_true("relerr[f[Transpose[{be, 2 be}]], Transpose[{xe, 2 xe}]] < 10^-55");
    /* applied repeatedly (the factor is cached, not recomputed) */
    run_check_true("relerr[Nest[f, N[be, 100], 2], LinearSolve[he, xe]] < 10^-20");
    /* a genuinely singular MPFR matrix is still singular */
    run_check_true("Quiet[Check[LinearSolve[N[{{1, 2}, {2, 4}}, 30]], \"sing\", "
                   "LinearSolve::sing1], LinearSolve::sing1] === \"sing\"");
    /* a machine right-hand side solves at machine precision */
    run_check_true("Precision[f[N[be]]] === MachinePrecision");
}

static void test_matrixrank_mpfr(void) {
    run_check_true("MatrixRank[N[hilb[40], 100]] === 40");
    run_check_true("MatrixRank[N[hilb[12], 50]] === 12");
    run_check_true("MatrixRank[N[{{1, 2}, {2, 4}}, 30]] === 1");
    run_check_true("MatrixRank[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 40]] === 2");
}

/* ---- LUBackSubstitution ------------------------------------------------ */

static void test_lubacksubstitution_exact(void) {
    /* {{1, 2}, {3, 4}} x = {5, 6}: x = {-4, 9/2} (Cramer: det -2) */
    run_check_true("LUBackSubstitution[LUDecomposition[{{1, 2}, {3, 4}}], {5, 6}] === {-4, 9/2}");
    /* second column {1, 0}: x = {-2, 3/2} */
    run_check_true("LUBackSubstitution[LUDecomposition[{{1, 2}, {3, 4}}], {{5, 1}, {6, 0}}] === "
                   "{{-4, -2}, {9/2, 3/2}}");
    /* {{1, 2 I}, {3, 4}} x = {1, 1}: det = 4 - 6 I, x1 = (4 - 2 I)/(4 - 6 I) */
    run_check_true("LUBackSubstitution[LUDecomposition[{{1, 2 I}, {3, 4}}], {1, 1}] === "
                   "{7/13 + 4 I/13, -2/13 - 3 I/13}");
    /* symbolic entries */
    run_check_true("Simplify[LUBackSubstitution[LUDecomposition[{{a, b}, {c, d}}], {p, q}] - "
                   "{(d p - b q)/(a d - b c), (a q - c p)/(a d - b c)}] === {0, 0}");
    /* agrees with LinearSolve on an exact 6 x 6 */
    run_check_true("LUBackSubstitution[LUDecomposition[hilb[6]], Range[6]] === LinearSolve[hilb[6], Range[6]]");
}

static void test_lubacksubstitution_machine(void) {
    run_check_true("Max[Abs[LUBackSubstitution[LUDecomposition[N[{{1, 2}, {3, 4}}]], {5, 6}] - {-4, 9/2}]] < 10^-14");
    run_check_true("Union[Head /@ LUBackSubstitution[LUDecomposition[N[{{1, 2}, {3, 4}}]], {5, 6}]] === {Real}");
    run_check_true("Max[Abs[LUBackSubstitution[LUDecomposition[{{1, 2}, {3, 4}}], {5., 6}] - {-4, 9/2}]] < 10^-14");
    run_check_true("Max[Abs[LUBackSubstitution[LUDecomposition[N[{{1, 2 I}, {3, 4}}]], {1, 1}] - "
                   "{7/13 + 4 I/13, -2/13 - 3 I/13}]] < 10^-14");
    /* a 200 x 200 packed system: residual small, answer packed */
    run_silent("SeedRandom[7]; mm = RandomReal[{-1, 1}, {200, 200}] + 10 IdentityMatrix[200]; "
               "bb = RandomReal[1, 200]; xx = LUBackSubstitution[LUDecomposition[mm], bb]");
    run_check_true("Norm[mm . xx - bb] < 10^-10");
    run_check_true("Developer`PackedArrayQ[xx]");
    run_check_true("Norm[xx - LinearSolve[mm, bb]] < 10^-10");
    run_check_true("Norm[LUBackSubstitution[LUDecomposition[mm], Transpose[{bb, 2 bb}]] - "
                   "Transpose[{xx, 2 xx}]] < 10^-10");
}

static void test_lubacksubstitution_mpfr(void) {
    run_check_true("relerr[LUBackSubstitution[LUDecomposition[N[{{1, 2}, {3, 4}}, 30]], {5, 6}], {-4, 9/2}] < 10^-28");
    run_check_true("Abs[Precision[LUBackSubstitution[LUDecomposition[N[{{1, 2}, {3, 4}}, 30]], {5, 6}]] - 30] < 1");
    /* the smaller precision of lu and b */
    run_check_true("Abs[Precision[LUBackSubstitution[LUDecomposition[N[{{1, 2}, {3, 4}}, 40]], N[{5, 6}, 25]]] - 25] < 1");
    run_check_true("relerr[LUBackSubstitution[LUDecomposition[N[{{1, 2 I}, {3, 4}}, 30]], {1, 1}], "
                   "{7/13 + 4 I/13, -2/13 - 3 I/13}] < 10^-28");
    /* the ill-conditioned case: H_24 at 100 digits against the exact solution */
    run_check_true("relerr[LUBackSubstitution[LUDecomposition[N[he, 100]], N[be, 100]], xe] < 10^-55");
    run_check_true("relerr[LUBackSubstitution[LUDecomposition[N[he, 100]], N[be, 100]], "
                   "LinearSolve[N[he, 100], N[be, 100]]] < 10^-55");
}

static void test_lubacksubstitution_messages(void) {
    run_check_true("Quiet[Check[LUBackSubstitution[Quiet[LUDecomposition[{{1, 2}, {2, 4}}]], {1, 2}], "
                   "\"s\", LUBackSubstitution::singlu], LUBackSubstitution::singlu] === \"s\"");
    run_check_true("Quiet[Check[LUBackSubstitution[{{{1, 2}, {3, 4}}, {1, 2}, 0}, {1, 2, 3}], "
                   "\"l\", LUBackSubstitution::lslc], LUBackSubstitution::lslc] === \"l\"");
    run_check_true("Quiet[Check[LUBackSubstitution[LUDecomposition[{{1, 2, 3}, {4, 5, 6}}], {1, 1}], "
                   "\"m\", LUBackSubstitution::matsq], LUBackSubstitution::matsq] === \"m\"");
    run_check_true("Head[LUBackSubstitution[{{{1, 2}, {3, 4}}, {1, 2}}, {1, 2}]] === LUBackSubstitution");
    run_check_true("Head[LUBackSubstitution[{{{1, 2}, {3, 4}}, {1, 1}, 0}, {1, 2}]] === LUBackSubstitution");
    run_check_true("MemberQ[Attributes[LUBackSubstitution], Protected]");
}

/* ---- Precision / Accuracy are scalars ---------------------------------- */

static void test_precision_accuracy_scalar(void) {
    run_check_true("Abs[Precision[N[hilb[4], 30]] - 30] < 0.5");
    run_check_true("NumberQ[Precision[N[hilb[4], 30]]]");
    run_check_true("Abs[Precision[{1, 2.5`20, N[1, 30]}] - 20] < 0.5");
    run_check_true("Precision[{1, 2.5, N[1, 30]}] === MachinePrecision");
    run_check_true("Precision[{{1.5}}] === MachinePrecision");
    run_check_true("Precision[{1, 2}] === Infinity");
    run_check_true("Precision[{x, 2.5`20}] == Precision[2.5`20]");
    /* Accuracy[1.5] = MachinePrecision - Log10[1.5]; 100.`30 has 28 */
    run_check_true("Abs[Accuracy[{1.5, 100.`30}] - (Log10[2.^53] - Log10[1.5])] < 10^-6");
    run_check_true("Accuracy[{1, 2}] === Infinity");
    /* packed arrays: an integer buffer is exact, a real one machine */
    run_check_true("Precision[Range[1000]] === Infinity");
    run_check_true("Accuracy[Range[1000]] === Infinity");
    run_check_true("Precision[N[Range[1000]]] === MachinePrecision");
    run_check_true("Abs[Accuracy[N[Range[1000]]] - (Log10[2.^53] - 3)] < 10^-6");
    run_check_true("Precision[SparseArray[{1 -> 1.5}, 4]] === MachinePrecision");
    run_check_true("Precision[SparseArray[{1 -> 2}, 4]] === Infinity");
    run_check_true("!MemberQ[Attributes[Precision], Listable] && !MemberQ[Attributes[Accuracy], Listable]");
}

/* ---- Dot at arbitrary precision ----------------------------------------- */

static void test_dot_mpfr(void) {
    run_silent("de = hilb[30]; dv = Table[(-1)^i i/7, {i, 30}]; dh = N[de, 60]; dw = N[dv, 60]");
    run_check_true("relerr[dh . dw, de . dv] < 10^-58");
    run_check_true("relerr[dw . dh, dv . de] < 10^-58");
    run_check_true("relerr[dw . dw, dv . dv] < 10^-58");
    run_check_true("relerr[dh . dh, de . de] < 10^-58");
    run_check_true("relerr[dh . dv, de . dv] < 10^-58");
    run_check_true("relerr[de . dw, de . dv] < 10^-58");
    run_check_true("Abs[Precision[dh . dw] - 60] < 1");
    /* the smaller precision wins, as for Times / Plus */
    run_check_true("Abs[Precision[N[de, 40] . dw] - 40] < 1");
    /* machine contagion, exact zeros and symbols keep the generic path */
    run_check_true("Precision[{N[1/3, 20], 1.5} . {3, 3}] === MachinePrecision");
    run_check_true("Coefficient[{N[1/3, 20], x} . {3, 3}, x] === 3");
    run_check_true("Head[{N[1, 20], N[1, 20]} . {N[1, 20], N[-1, 20]}] === Real");
}

/* An exact operand with zeros: Times[0, x] of an MPFR x is an MPFR zero at
 * x's precision, so m . DiagonalMatrix[w] is the column scaling m_ij w_j,
 * every entry MPFR at the smaller precision. */
static void test_dot_mpfr_exact_zeros(void) {
    run_silent("zm = N[Table[1/(i + 2 j), {i, 12}, {j, 12}], 30]; zw = Table[j/5, {j, 12}]; "
               "zd = DiagonalMatrix[zw]; zx = Table[1/(i + 2 j), {i, 12}, {j, 12}]");
    run_check_true("relerr[zm . zd, zx . zd] < 10^-29");
    run_check_true("relerr[zd . zm, zd . zx] < 10^-29");
    run_check_true("Abs[Precision[zm . zd] - 30] < 1");
    run_check_true("MatrixQ[zm . zd, Head[#] === Real &]");
    /* an MPFR diagonal at a lower precision sets the result's */
    run_check_true("Abs[Precision[zm . DiagonalMatrix[N[zw, 20]]] - 20] < 1");
    run_check_true("relerr[zm . DiagonalMatrix[N[zw, 20]], zx . zd] < 10^-19");
    /* a row of zeros gives MPFR zeros, as the generic Times / Plus do */
    run_check_true("Module[{r = {{0, 0}, {1, 2}} . N[{{1/3, 1/5}, {1/7, 1/9}}, 20]}, "
                   "Head[r[[1, 1]]] === Real && r[[1, 1]] == 0 && "
                   "Abs[r[[2, 1]] - 13/21] < 10^-19]");
    run_check_true("relerr[{0, 1, 0, 2} . N[hilb[4], 25], {0, 1, 0, 2} . hilb[4]] < 10^-24");
}

/* ---- Eigenvalues of N[Hilbert[n], p] ------------------------------------ */

static void test_eigenvalues_hilbert_mpfr(void) {
    /* The eigenvalues are the roots of the exact characteristic polynomial.
     * For each computed eigenvalue e, the Newton step cp(r)/cp'(r) taken in
     * EXACT rational arithmetic at r = the exact value of e is, to first
     * order, its error -- independent of the eigen solver. */
    run_silent("cp6 = CharacteristicPolynomial[hilb[6], t]; dcp6 = D[cp6, t]; "
               "ev6 = Eigenvalues[N[hilb[6], 50]]; "
               "nerr[p_, dp_, e_] := Module[{r = SetPrecision[e, Infinity]}, "
               "N[Abs[(p /. t -> r)/((dp /. t -> r) r)], 20]]");
    run_check_true("Length[ev6] == 6 && Min[ev6] > 0");
    /* every eigenvalue, including the smallest (~1.1e-7), to relative accuracy */
    run_check_true("Max[nerr[cp6, dcp6, #] & /@ ev6] < 10^-40");
    /* kappa(H_6) = lambda_max / lambda_min = 1.4951059e7 */
    run_check_true("Abs[Max[ev6]/Min[ev6] - 1.4951059*^7] < 10");
    /* n = 10 at 60 digits: the smallest eigenvalue is ~1.1e-13 */
    run_silent("cp10 = CharacteristicPolynomial[hilb[10], t]; dcp10 = D[cp10, t]; "
               "ev10 = Eigenvalues[N[hilb[10], 60]]");
    run_check_true("Max[nerr[cp10, dcp10, #] & /@ ev10] < 10^-40");
    run_check_true("Abs[Max[ev10]/Min[ev10] - 1.6026*^13] < 10^10");
}

/* ---- Eigenvalues of nonsymmetric MPFR matrices ------------------------- */

static void test_eigenvalues_nonsymmetric_mpfr(void) {
    /* D H D^-1 with D = diag(1..n): nonsymmetric, the Hilbert spectrum
     * (smallest eigenvalue ~1.1e-13 for n = 10), each eigenvalue checked
     * by an exact Newton step on the exact characteristic polynomial. */
    run_silent("ns10 = Table[i/j/(i + j - 1), {i, 10}, {j, 10}]; "
               "nsev = Eigenvalues[N[ns10, 60]]");
    run_check_true("Length[nsev] == 10 && FreeQ[nsev, Complex] && Min[nsev] > 0");
    run_check_true("Max[nerr[cp10, dcp10, #] & /@ nsev] < 10^-40");
    /* complex eigenvalues: a real 5 x 5 with two conjugate pairs */
    run_silent("cm5 = {{1, 2, 3, 0, 1}, {-4, 5, 6, 2, 0}, {7, -8, 9, 1, 3}, "
               "{0, 1, -2, 3, 5}, {2, 0, 1, -6, 4}}; "
               "cpc = CharacteristicPolynomial[cm5, t]; dcpc = D[cpc, t]; "
               "evc = Eigenvalues[N[cm5, 40]]");
    run_check_true("Length[evc] == 5 && Count[evc, _Complex] == 4");
    run_check_true("Max[nerr[cpc, dcpc, #] & /@ evc] < 10^-35");
    /* conjugate pairs come out as exact conjugates */
    run_check_true("Sort[Im[evc]] == -Reverse[Sort[Im[evc]]]");
    /* Eigensystem: A v == lambda v at the working precision */
    run_silent("{esv, esx} = Eigensystem[N[cm5, 40]]");
    run_check_true("Max[Table[Norm[N[cm5, 40] . esx[[k]] - esv[[k]] esx[[k]]], {k, 5}]] < 10^-35");
    /* a 3 x 3 with exact eigenvalues 1, 2, 3 (P diag(1,2,3) P^-1) */
    run_silent("pp = {{1, 2, 0}, {0, 1, 1}, {1, 0, 1}}; "
               "ev3 = Eigenvalues[N[pp . DiagonalMatrix[{1, 2, 3}] . Inverse[pp], 50]]");
    run_check_true("Max[Abs[ev3 - {3, 2, 1}]] < 10^-48");
}

/* The Nystrom matrix of the smooth kernel E^(x + y + x^2 + y^2 + x y + x^2 y^2)
 * on 32 Gauss-Legendre nodes of (0, 1) (and 64, on two halves): numerically
 * low rank, so most of its eigenvalues are rounding noise.  The largest
 * eigenvalue of the integral operator is 37.52914556033537842820658... (the
 * 40-digit, 64-node value, converged in the node count). */
static void test_eigenvalues_nystrom_mpfr(void) {
    run_silent("nyK[x_, y_] := E^(x + y + x^2 + y^2 + x y + x^2 y^2); "
               "nyGL[n_, prec_] := Module[{J, ev, evec, b}, b[k_] := k/Sqrt[4 k^2 - 1]; "
               "J = N[Table[Which[j == i + 1, b[i], i == j + 1, b[j], True, 0], {i, n}, {j, n}], prec]; "
               "{ev, evec} = Eigensystem[J]; "
               "Transpose[Sort[Transpose[{(1 + ev)/2, evec[[All, 1]]^2}]]]]; "
               "nyMat[prec_, c_] := Module[{xs, ws, pts, wts}, {xs, ws} = nyGL[32, prec]; "
               "pts = Flatten[Table[xs + s, {s, 0, c - 1}]/c]; wts = Flatten[Table[ws, {c}]/c]; "
               "Outer[nyK, pts, pts] . DiagonalMatrix[wts]]; "
               "nyTop[prec_, c_] := Max[Re[Eigenvalues[nyMat[prec, c]]]]");
    /* 32 nodes: converged to ~5e-17 in the node count */
    run_check_true("Abs[nyTop[20, 1] - 37.529145560335378428206580483543`40] < 10^-15");
    /* 20 digits (67 bits): rounding at ~1e-18 relative */
    run_check_true("Abs[nyTop[20, 2] - 37.529145560335378428206580483543`40] < 10^-16");
    run_check_true("Abs[nyTop[40, 2] - 37.529145560335378428206580483543`40] < 10^-30");
    run_check_true("Abs[nyTop[MachinePrecision, 2] - 37.529145560335] < 10^-11");
    /* every eigenvalue above the rounding level agrees with the same matrix
     * solved at 60 digits (its entries promoted exactly) */
    run_silent("nyM = nyMat[30, 1]; nyE = Eigenvalues[nyM]; "
               "nyR = Eigenvalues[SetPrecision[nyM, 60]]");
    run_check_true("Length[nyE] == 32");
    run_check_true("Max[Abs[Take[nyE, 8] - Take[nyR, 8]]] < 10^-27");
}

/* ---- Outer through compiled code ---------------------------------------- */

static void test_outer_compiled(void) {
    run_silent("og[x_, y_] := Sin[x] y - Log[1 + x^2]/(2 + y); "
               "oh[x_, y_, z_] := x y - z^2; "
               "oc = 0; ocount[x_, y_] := (oc++; x + y)");
    for (int k = 0; k < 3; k++) {
        const char* p = k == 0 ? "MachinePrecision" : k == 1 ? "20" : "45";
        char buf[512];
        snprintf(buf, sizeof buf, "op = N[Range[24]/25, %s]; oq = N[Range[17]/7 - 1, %s]", p, p);
        run_silent(buf);
        /* bit-identical to the interpreter evaluating each cell */
        run_check_true("Outer[nyK, op, oq] === Table[nyK[u, v], {u, op}, {v, oq}]");
        run_check_true("Outer[og, op, oq] === Table[og[u, v], {u, op}, {v, oq}]");
        run_check_true("Outer[Function[{s, t}, s t + 1], op, oq] === Table[u v + 1, {u, op}, {v, oq}]");
        run_check_true("Outer[oh, op, oq, Take[op, 3]] === "
                       "Table[oh[u, v, w], {u, op}, {v, oq}, {w, Take[op, 3]}]");
        /* a cell the compiled code declines (Log of a negative number) is the
         * interpreter's complex value */
        run_check_true("Outer[Function[{s, t}, Log[s - t]], op, oq] === "
                       "Table[Log[u - v], {u, op}, {v, oq}]");
    }
    /* a side-effecting definition still runs once per cell */
    run_check_true("oc = 0; Outer[ocount, N[Range[6]], N[Range[5]]]; oc == 30");
    /* Orderless: the arguments are sorted before the rule binds them */
    run_silent("SetAttributes[oo, Orderless]; oo[x_, y_] := x - y");
    run_check_true("Outer[oo, N[Range[5]], N[Range[5]]] === Table[oo[u, v], {u, N[Range[5]]}, {v, N[Range[5]]}]");
    /* two definitions: the general path */
    run_silent("o2[x_, y_] := x + y; o2[1., y_] := 0");
    run_check_true("Outer[o2, N[Range[5]], N[Range[5]]][[1]] === {0, 0, 0, 0, 0}");
}

/* ---- Conjugate of a huge MPFR number ------------------------------------ */

static void test_conjugate_huge_mpfr(void) {
    run_check_true("Conjugate[N[10^400, 30]] === N[10^400, 30]");
    run_check_true("Conjugate[N[-10^3000, 30]] === N[-10^3000, 30]");
    run_check_true("Conjugate[10^400] === 10^400");
}

int main(void) {
    alarm(600);
    symtab_init();
    core_init();
    define_helpers();
    printf("Running arbitrary-precision linear algebra tests...\n");
    TEST(test_linearsolve_hilbert_mpfr);
    TEST(test_matrixrank_mpfr);
    TEST(test_lubacksubstitution_exact);
    TEST(test_lubacksubstitution_machine);
    TEST(test_lubacksubstitution_mpfr);
    TEST(test_lubacksubstitution_messages);
    TEST(test_precision_accuracy_scalar);
    TEST(test_dot_mpfr);
    TEST(test_dot_mpfr_exact_zeros);
    TEST(test_eigenvalues_hilbert_mpfr);
    TEST(test_eigenvalues_nonsymmetric_mpfr);
    TEST(test_eigenvalues_nystrom_mpfr);
    TEST(test_outer_compiled);
    TEST(test_conjugate_huge_mpfr);
    printf("All arbitrary-precision linear algebra tests passed!\n");
    symtab_clear();
    return 0;
}

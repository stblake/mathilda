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
 *     characteristic polynomial.
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
    TEST(test_eigenvalues_hilbert_mpfr);
    TEST(test_conjugate_huge_mpfr);
    printf("All arbitrary-precision linear algebra tests passed!\n");
    symtab_clear();
    return 0;
}

/* End-to-end tests for the NDArray[...] fast paths across the linear-algebra
 * builtins (src/linalg/ndlinalg.c and the per-builtin guards).
 *
 * Coverage goals:
 *   - every consumer routine accepts an NDArray argument and returns the right
 *     value (real float64, float32, and complex64 dtypes);
 *   - the fast path agrees with the ordinary List path on the same data;
 *   - matrix-valued heavy ops return a closed-system NDArray (Head == NDArray)
 *     with the expected dtype;
 *   - error / edge shapes (non-square, singular, rectangular) behave sanely;
 *   - the delist-backed routines (Eigen*, QR, LU, SVD, LeastSquares, NullSpace,
 *     RowReduce, PseudoInverse, PositiveDefiniteMatrixQ, constructors, LLL)
 *     evaluate correctly on NDArray input.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------------- Det ---------------- */
static void test_det_real(void) {
    assert_eval_eq("Det[NDArray[{{1., 2.}, {3., 4.}}]]", "-2.0", 0);
    assert_eval_eq("Det[NDArray[{{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}}]]", "18.0", 0);
    /* Fast path agrees with the ordinary List path. */
    assert_eval_eq("Det[NDArray[{{1., 2.}, {3., 4.}}]] == Det[{{1., 2.}, {3., 4.}}]", "True", 0);
    /* float32 and complex64 dtypes take the fast path too. */
    assert_eval_eq("Det[NDArray[{{1, 2}, {3, 4}}, DataType -> \"float32\"]]", "-2.0", 0);
    assert_eval_eq("Det[NDArray[{{1, 2}, {3, 4}}, DataType -> \"complex64\"]]", "-2.0", 0);
}

static void test_det_edge(void) {
    /* Non-square: delist -> Det::matsq, call left unevaluated on the List form. */
    assert_eval_eq("Det[NDArray[{{1., 2., 3.}, {4., 5., 6.}}]]",
                   "Det[{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}]", 0);
    /* Singular matrix: determinant is exactly 0. */
    assert_eval_eq("Det[NDArray[{{1., 2.}, {2., 4.}}]]", "0.0", 0);
}

/* ---------------- Inverse ---------------- */
static void test_inverse_real(void) {
    assert_eval_eq("Inverse[NDArray[{{4., 3.}, {6., 3.}}]]",
                   "NDArray[{{-0.5, 0.5}, {1.0, -0.666667}}]", 0);
    /* Result is a closed-system NDArray of float64. */
    assert_eval_eq("Head[Inverse[NDArray[{{4., 3.}, {6., 3.}}]]]", "NDArray", 0);
    assert_eval_eq("DataType[Inverse[NDArray[{{4., 3.}, {6., 3.}}]]]", "\"float64\"", 0);
    /* A . A^-1 == I (Chop kills rounding noise). */
    assert_eval_eq("Chop[Dot[Inverse[NDArray[{{4., 3.}, {6., 3.}}]], "
                   "NDArray[{{4., 3.}, {6., 3.}}]] - NDArray[{{1., 0.}, {0., 1.}}]]",
                   "NDArray[{{0.0, 0.0}, {0.0, 0.0}}]", 0);
}

static void test_inverse_edge(void) {
    /* Singular: delist -> Inverse::sing, left unevaluated on the List form. */
    assert_eval_eq("Inverse[NDArray[{{1., 2.}, {2., 4.}}]]",
                   "Inverse[{{1.0, 2.0}, {2.0, 4.0}}]", 0);
}

static void test_inverse_complex(void) {
    assert_eval_eq("DataType[Inverse[NDArray[{{4, 3}, {6, 3}}, DataType -> \"complex64\"]]]",
                   "\"complex64\"", 0);
}

/* ---------------- LinearSolve ---------------- */
static void test_linearsolve(void) {
    assert_eval_eq("LinearSolve[NDArray[{{2., 0.}, {0., 4.}}], NDArray[{2., 8.}]]",
                   "NDArray[{1.0, 2.0}]", 0);
    /* Matrix rhs. */
    assert_eval_eq("LinearSolve[NDArray[{{2., 0.}, {0., 4.}}], NDArray[{{2., 4.}, {8., 4.}}]]",
                   "NDArray[{{1.0, 2.0}, {2.0, 1.0}}]", 0);
    /* Complex system. */
    assert_eval_eq("LinearSolve[NDArray[{{1, 0}, {0, 2}}, DataType -> \"complex64\"], "
                   "NDArray[{Complex[1, 1], 4}, DataType -> \"complex64\"]]",
                   "NDArray[{1.0 + 1.0*I, 2.0 + 0.0*I}]", 0);
    /* Singular: delist, left unevaluated. */
    assert_eval_startswith("LinearSolve[NDArray[{{1., 2.}, {2., 4.}}], NDArray[{1., 1.}]]",
                           "LinearSolve[");
}

/* ---------------- MatrixRank ---------------- */
static void test_matrixrank(void) {
    assert_eval_eq("MatrixRank[NDArray[{{1., 0.}, {0., 1.}}]]", "2", 0);
    assert_eval_eq("MatrixRank[NDArray[{{1., 2.}, {2., 4.}}]]", "1", 0);
    assert_eval_eq("MatrixRank[NDArray[{{1., 2., 3.}, {2., 4., 6.}, {1., 1., 1.}}]]", "2", 0);
    /* Agrees with the List path. */
    assert_eval_eq("MatrixRank[NDArray[{{1., 2.}, {3., 4.}}]] == MatrixRank[{{1., 2.}, {3., 4.}}]",
                   "True", 0);
}

/* ---------------- Tr ---------------- */
static void test_tr(void) {
    assert_eval_eq("Tr[NDArray[{{1., 2.}, {3., 4.}}]]", "5.0", 0);
    /* Rectangular: sum over the min-dimension diagonal. */
    assert_eval_eq("Tr[NDArray[{{1., 2., 3.}, {4., 5., 6.}}]]", "6.0", 0);
    /* Complex diagonal collapses an exactly-zero imaginary part. */
    assert_eval_eq("Tr[NDArray[{{Complex[1, 1], 2}, {3, Complex[4, -1]}}, DataType -> \"complex64\"]]",
                   "5.0", 0);
}

/* ---------------- Norm / Normalize ---------------- */
static void test_norm(void) {
    assert_eval_eq("Norm[NDArray[{3., 4.}]]", "5.0", 0);
    assert_eval_eq("Norm[NDArray[{1., 2., 2.}], 1]", "5.0", 0);
    assert_eval_eq("Norm[NDArray[{1., -2., 2.}], Infinity]", "2.0", 0);
    /* Complex vector 2-norm uses the modulus. */
    assert_eval_eq("Norm[NDArray[{Complex[3, 4], 0}, DataType -> \"complex64\"]]", "5.0", 0);
}

static void test_normalize(void) {
    assert_eval_eq("Normalize[NDArray[{3., 4.}]]", "NDArray[{0.6, 0.8}]", 0);
    assert_eval_eq("Head[Normalize[NDArray[{3., 4.}]]]", "NDArray", 0);
    /* Zero vector normalises to itself. */
    assert_eval_eq("Normalize[NDArray[{0., 0.}]]", "NDArray[{0.0, 0.0}]", 0);
}

/* ---------------- Cross ---------------- */
static void test_cross(void) {
    assert_eval_eq("Cross[NDArray[{1., 0., 0.}], NDArray[{0., 1., 0.}]]",
                   "NDArray[{0.0, 0.0, 1.0}]", 0);
    assert_eval_eq("Cross[NDArray[{1., 2., 3.}], NDArray[{4., 5., 6.}]]",
                   "NDArray[{-3.0, 6.0, -3.0}]", 0);
}

/* ---------------- Delist-backed factorisations / predicates ---------------- */
static void test_eigen(void) {
    /* Previously crashed on NDArray; now delists to the numeric kernel. */
    assert_eval_eq("Eigenvalues[NDArray[{{2., 0.}, {0., 3.}}]]", "{3.0, 2.0}", 0);
    assert_eval_eq("Eigenvalues[NDArray[{{2., 0.}, {0., 3.}}]] == "
                   "Eigenvalues[{{2., 0.}, {0., 3.}}]", "True", 0);
}

static void test_factorisations_evaluate(void) {
    /* Each must produce a numeric result (not stay unevaluated) on NDArray. */
    assert_eval_eq("MatrixQ[QRDecomposition[NDArray[{{1., 2.}, {3., 4.}}]][[1]]]", "True", 0);
    assert_eval_eq("MatrixQ[LUDecomposition[NDArray[{{2., 1.}, {1., 3.}}]][[1]]]", "True", 0);
    assert_eval_eq("Length[SingularValueDecomposition[NDArray[{{1., 0.}, {0., 2.}}]]]", "3", 0);
}

/* ---------------- The SVD pair: PseudoInverse / LeastSquares ---------------- */
static void test_svd_pair(void) {
    /* Both used to delist and run the exact rationalised pipeline in inv.c,
     * which does not terminate on a machine matrix of any size worth timing.
     * They now go through one thin gesdd, so a visible NDArray operand answers
     * with a visible NDArray -- the same convention Normalize, Cross and
     * LinearSolve already follow, and the reason the assertion below is exact
     * where it used to only check for a leading brace. */
    assert_eval_eq("LeastSquares[NDArray[{{1., 0.}, {0., 2.}}], NDArray[{2., 6.}]]",
                   "NDArray[{2.0, 3.0}]", 0);
    /* LeastSquares on a square nonsingular system reproduces LinearSolve. */
    assert_eval_eq("LeastSquares[NDArray[{{1., 0.}, {0., 2.}}], NDArray[{2., 6.}]] == "
                   "LinearSolve[NDArray[{{1., 0.}, {0., 2.}}], NDArray[{2., 6.}]]", "True", 0);
    assert_eval_eq("PseudoInverse[NDArray[{{2., 0.}, {0., 4.}}]]",
                   "NDArray[{{0.5, 0.0}, {0.0, 0.25}}]", 0);
    /* The Moore-Penrose identity A A^+ A == A, on a rank-deficient matrix the
     * exact pipeline could not have been asked for at this size. */
    assert_eval_eq("Module[{a = NDArray[{{1., 2.}, {2., 4.}}]}, "
                   "Max[Abs[Flatten[a . PseudoInverse[a] . a - a]]] < 1.*^-12]", "True", 0);
    /* An EXACT matrix keeps the exact pipeline: no buffer, no SVD. */
    assert_eval_eq("PseudoInverse[{{1, 2}, {3, 4}}]", "{{-2, 1}, {3/2, -1/2}}", 0);
}

static void test_predicates_and_constructors(void) {
    assert_eval_eq("PositiveDefiniteMatrixQ[NDArray[{{2., 0.}, {0., 3.}}]]", "True", 0);
    assert_eval_eq("PositiveDefiniteMatrixQ[NDArray[{{-2., 0.}, {0., 3.}}]]", "False", 0);
    assert_eval_eq("NegativeDefiniteMatrixQ[NDArray[{{-2., 0.}, {0., -3.}}]]", "True", 0);
    /* Constructor accepts an NDArray vector. The zeros are REAL, not exact:
     * a machine Real on the diagonal makes the whole matrix machine-real, which
     * is Mathematica's answer and what lets the result pack. */
    assert_eval_eq("DiagonalMatrix[NDArray[{1., 2.}]]", "{{1.0, 0.0}, {0.0, 2.0}}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_det_real);
    TEST(test_det_edge);
    TEST(test_inverse_real);
    TEST(test_inverse_edge);
    TEST(test_inverse_complex);
    TEST(test_linearsolve);
    TEST(test_matrixrank);
    TEST(test_tr);
    TEST(test_norm);
    TEST(test_normalize);
    TEST(test_cross);
    TEST(test_eigen);
    TEST(test_factorisations_evaluate);
    TEST(test_svd_pair);
    TEST(test_predicates_and_constructors);

    printf("All NDArray linalg fast-path tests passed.\n");
    return 0;
}

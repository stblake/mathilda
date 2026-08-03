/* Unit tests for the COMPILE_MISSING.md §4 / §5 Compile[] lowerings:
 *
 *   §4  Transforms (real in, complex or real out), ND_FNS / A_NDFN:
 *         Fourier, InverseFourier  (real/int -> complex, complex_result flag),
 *         FourierDCT, FourierDST   (real -> real, delegated to their builtins).
 *   §5  Matrix producers (rank 1 -> rank 2), ND_FNS / A_NDFN with rank_rule 5:
 *         DiagonalMatrix, HankelMatrix, ToeplitzMatrix, VandermondeMatrix.
 *
 * Each head gets, through the eval path (as in test_compile_linalg.c):
 *   1. "it lowers":  CompileDiagnostics[spec, H[v]] contains Compiled -> True.
 *   2. "compiled == interpreted":  Normal[compiled] equals the interpreter's OWN
 *      NDArray-path result, Normal[H[NDArray[...]]] — bit-identical by
 *      construction (same delegate), Normal[] removing the packed-vs-plain
 *      difference.  NEVER H[plainList] as the reference (that hits the exact
 *      path, e.g. an int Vandermonde would stay exact where the machine path
 *      packs).  Fourier inputs are chosen so the interpreter also answers
 *      complex (asymmetric) or exactly real (symmetric), so the two agree — the
 *      compiled form commits to a complex register, the interpreter's
 *      data-dependent real-collapse only differs for a tiny nonzero imaginary
 *      part, which no test pins.
 *   3. "cliff": a COMPOSED body (Total[Abs[Fourier[v]]], Tr[VandermondeMatrix[v]])
 *      compiles whole — exercises the shared nd_fn_result via the infer_type
 *      branch.  Total[Abs[Fourier[v]]] also proves the register is genuinely
 *      COMPLEX (Abs reads the imaginary part; a mis-typed real register would
 *      give a different sum).
 *   4. "clean decline": a gated operand (a rank-2 DiagonalMatrix operand, an
 *      int FourierDCT operand) reports Compiled -> False.
 *   5. Part C — the interpreter's direct rank-2 buffer producers: numeric forms
 *      pack (NDArrayQ), a mixed int+real Vandermonde stays unpacked, and an
 *      int64-overflowing Vandermonde falls back to the exact bignum answer.
 *   6. A repeated-evaluation loop to surface double-free / leak regressions.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Evaluate `input` and assert the printed result CONTAINS `substr` (the
 * CompileDiagnostics association's full text is verbose and version-y). */
static void assert_eval_contains(const char* input, const char* substr) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, substr) != NULL,
               "expected result of %s to contain \"%s\", got: %s",
               input, substr, s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

static void assert_lowers(const char* diag_call) {
    assert_eval_contains(diag_call, "Compiled\" -> True");
}

/* ------------------------------------------------------------------ *
 *  §4 transforms                                                     *
 * ------------------------------------------------------------------ */

static void test_cf_fourier(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Fourier[v]]");
    /* asymmetric input -> the interpreter also answers complex, so the compiled
     * complex result is bit-identical. */
    assert_true("Normal[Compile[{{v, _Real, 1}}, Fourier[v]][{1., 2., 3., 4.}]]"
                " === Normal[Fourier[NDArray[{1., 2., 3., 4.}]]]");
    /* symmetric input -> imaginary parts are exactly 0; per-element boxing
     * collapses Complex[x, 0.] -> x, so the two still agree. */
    assert_true("Normal[Compile[{{v, _Real, 1}}, Fourier[v]][{1., 1.}]]"
                " === Normal[Fourier[NDArray[{1., 1.}]]]");
    /* int operand: real/int in -> complex out (nd_gather widens int to real). */
    assert_lowers("CompileDiagnostics[{{v, _Integer, 1}}, Fourier[v]]");
    assert_true("Normal[Compile[{{v, _Integer, 1}}, Fourier[v]][{1, 2, 3, 4}]]"
                " === Normal[Fourier[NDArray[{1, 2, 3, 4}]]]");
    /* cliff + proof the register is COMPLEX: Abs reads the imaginary part. */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Total[Abs[Fourier[v]]]]");
    assert_true("Compile[{{v, _Real, 1}}, Total[Abs[Fourier[v]]]][{1., 2., 3., 4.}]"
                " === Total[Abs[Fourier[NDArray[{1., 2., 3., 4.}]]]]");
}

static void test_cf_inverse_fourier(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, InverseFourier[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, InverseFourier[v]][{1., 2., 3., 4.}]]"
                " === Normal[InverseFourier[NDArray[{1., 2., 3., 4.}]]]");
    /* Fourier followed by InverseFourier round-trips to the input. */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, InverseFourier[Fourier[v]]]");
}

static void test_cf_fourier_dct(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, FourierDCT[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, FourierDCT[v]][{1., 2., 3., 4.}]]"
                " === Normal[FourierDCT[NDArray[{1., 2., 3., 4.}]]]");
    /* real -> real, so an int operand is gated out (would mis-promise CT_INT). */
    assert_true("Lookup[CompileDiagnostics[{{v, _Integer, 1}}, FourierDCT[v]], \"Compiled\"]"
                " === False");
}

static void test_cf_fourier_dst(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, FourierDST[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, FourierDST[v]][{1., 2., 3., 4.}]]"
                " === Normal[FourierDST[NDArray[{1., 2., 3., 4.}]]]");
}

/* ------------------------------------------------------------------ *
 *  §5 matrix producers (rank 1 -> rank 2)                            *
 * ------------------------------------------------------------------ */

static void test_cf_diagonalmatrix(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, DiagonalMatrix[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, DiagonalMatrix[v]][{1., 2., 3.}]]"
                " === Normal[DiagonalMatrix[NDArray[{1., 2., 3.}]]]");
    /* integer operand: preserves int64.  The reference MUST spell DataType ->
     * "int64" — a bare NDArray[{4,5,6}] defaults to float64, which would not
     * match the int64 the _Integer register produces. */
    assert_true("Normal[Compile[{{v, _Integer, 1}}, DiagonalMatrix[v]][{4, 5, 6}]]"
                " === Normal[DiagonalMatrix[NDArray[{4, 5, 6}, DataType -> \"int64\"]]]");
    /* rank_rule 5 requires a rank-1 operand: a rank-2 operand does NOT lower. */
    assert_true("Lookup[CompileDiagnostics[{{m, _Real, 2}}, DiagonalMatrix[m]], \"Compiled\"]"
                " === False");
    /* cliff: Total[Flatten[...]] over the produced matrix compiles whole. */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Total[Flatten[DiagonalMatrix[v]]]]");
    assert_true("Compile[{{v, _Real, 1}}, Total[Flatten[DiagonalMatrix[v]]]][{1., 2., 3.}]"
                " === 6.");
}

static void test_cf_hankelmatrix(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, HankelMatrix[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, HankelMatrix[v]][{1., 2., 3.}]]"
                " === Normal[HankelMatrix[NDArray[{1., 2., 3.}]]]");
    assert_true("Normal[Compile[{{v, _Integer, 1}}, HankelMatrix[v]][{1, 2, 3, 4}]]"
                " === Normal[HankelMatrix[NDArray[{1, 2, 3, 4}, DataType -> \"int64\"]]]");
}

static void test_cf_toeplitzmatrix(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, ToeplitzMatrix[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, ToeplitzMatrix[v]][{1., 2., 3.}]]"
                " === Normal[ToeplitzMatrix[NDArray[{1., 2., 3.}]]]");
    assert_true("Normal[Compile[{{v, _Integer, 1}}, ToeplitzMatrix[v]][{5, 6, 7}]]"
                " === Normal[ToeplitzMatrix[NDArray[{5, 6, 7}, DataType -> \"int64\"]]]");
}

static void test_cf_vandermondematrix(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, VandermondeMatrix[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, VandermondeMatrix[v]][{1., 2., 3.}]]"
                " === Normal[VandermondeMatrix[NDArray[{1., 2., 3.}]]]");
    assert_true("Normal[Compile[{{v, _Integer, 1}}, VandermondeMatrix[v]][{2, 3, 4}]]"
                " === Normal[VandermondeMatrix[NDArray[{2, 3, 4}, DataType -> \"int64\"]]]");
    /* cliff: the trace of the produced matrix compiles whole. */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Tr[VandermondeMatrix[v]]]");
    assert_true("Compile[{{v, _Real, 1}}, Tr[VandermondeMatrix[v]]][{1., 2., 3.}]"
                " === Tr[VandermondeMatrix[NDArray[{1., 2., 3.}]]]");
}

/* ------------------------------------------------------------------ *
 *  Part C — direct rank-2 buffer producers (interpreter)            *
 * ------------------------------------------------------------------ */

static void test_buffer_producers(void) {
    /* numeric single-vector forms pack (NDArrayQ True). */
    assert_true("NDArrayQ[HankelMatrix[{1, 2, 3}]]");
    assert_true("NDArrayQ[ToeplitzMatrix[{5, 6, 7}]]");
    assert_true("NDArrayQ[VandermondeMatrix[{1, 2, 3}]]");
    assert_true("NDArrayQ[HankelMatrix[{1., 2., 3.}]]");
    assert_true("NDArrayQ[VandermondeMatrix[{1., 2., 3.}]]");
    /* Hankel/Toeplitz coerce a mixed int+real argument to all-real (contagion),
     * so it packs; the values match the boxed result. */
    assert_true("Normal[HankelMatrix[{1, 2, 3.}]] === {{1., 2., 3.}, {2., 3., 0.}, {3., 0., 0.}}");
    assert_true("Normal[ToeplitzMatrix[{1, 2, 3.}]] === {{1., 2., 3.}, {2., 1., 2.}, {3., 2., 1.}}");
    /* Vandermonde does NOT coerce (Power[node, e] keeps the node's head), so a
     * mixed list stays a two-head, UNPACKED result — unchanged by the buffer
     * path (its guard requires a uniform node head). */
    assert_true("NDArrayQ[VandermondeMatrix[{1, 2, 3.}]] === False");
    assert_true("FullForm[VandermondeMatrix[{1, 2, 3.}]]"
                " === FullForm[{{1., 1, 1}, {1., 2, 4}, {1., 3., 9.}}]");
    /* integer Vandermonde whose powers overflow int64: the checked-power path
     * abandons the buffer and the exact bignum answer is produced. */
    assert_true("VandermondeMatrix[{2, 100}, 20][[2, 20]] === 100^19");
    assert_true("NDArrayQ[VandermondeMatrix[{2, 100}, 20]] === False");
    /* the integer-n and two-vector forms are unchanged. */
    assert_true("Normal[HankelMatrix[4]] === {{1,2,3,4},{2,3,4,0},{3,4,0,0},{4,0,0,0}}");
    assert_true("Normal[ToeplitzMatrix[{1,2,3},{1,4,5}]] === {{1,4,5},{2,1,4},{3,2,1}}");
    /* symbolic argument still flows through as a nested List. */
    assert_true("HankelMatrix[{a, b}] === {{a, b}, {b, 0}}");
    assert_true("VandermondeMatrix[{a, b}] === {{1, a}, {1, b}}");
}

/* ------------------------------------------------------------------ *
 *  Memory-safety stress loop                                         *
 * ------------------------------------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    for (int t = 0; t < 50; t++) {
        assert_true("Normal[Compile[{{v, _Real, 1}}, Fourier[v]][{1., 2., 3., 4.}]]"
                    " === Normal[Fourier[NDArray[{1., 2., 3., 4.}]]]");
        assert_true("Normal[Compile[{{v, _Real, 1}}, FourierDCT[v]][{1., 2., 3.}]]"
                    " === Normal[FourierDCT[NDArray[{1., 2., 3.}]]]");
        assert_true("Normal[Compile[{{v, _Real, 1}}, DiagonalMatrix[v]][{1., 2., 3.}]]"
                    " === Normal[DiagonalMatrix[NDArray[{1., 2., 3.}]]]");
        assert_true("Normal[Compile[{{v, _Real, 1}}, VandermondeMatrix[v]][{1., 2., 3.}]]"
                    " === Normal[VandermondeMatrix[NDArray[{1., 2., 3.}]]]");
        assert_true("Compile[{{v, _Real, 1}}, Total[Abs[Fourier[v]]]][{1., 2., 3., 4.}]"
                    " === Total[Abs[Fourier[NDArray[{1., 2., 3., 4.}]]]]");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cf_fourier);
    TEST(test_cf_inverse_fourier);
    TEST(test_cf_fourier_dct);
    TEST(test_cf_fourier_dst);

    TEST(test_cf_diagonalmatrix);
    TEST(test_cf_hankelmatrix);
    TEST(test_cf_toeplitzmatrix);
    TEST(test_cf_vandermondematrix);

    TEST(test_buffer_producers);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All Compile transforms tests passed!\n");
    symtab_clear();
    return 0;
}

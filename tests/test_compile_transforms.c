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

/* Two-vector forms HankelMatrix[c, r] / ToeplitzMatrix[c, r]: rank 1 x rank 1
 * -> rank 2, via ND_FN2S / A_NDFN2 with R2_MATRIX + NDF2_SAME.  Corners are kept
 * matching (c_last == r_first for Hankel, c_first == r_first for Toeplitz) so no
 * ::crs warning fires. */
static void test_cf_hankelmatrix_2vec(void) {
    assert_lowers("CompileDiagnostics[{{c, _Real, 1}, {r, _Real, 1}}, HankelMatrix[c, r]]");
    assert_true("Normal[Compile[{{c, _Real, 1}, {r, _Real, 1}}, HankelMatrix[c, r]]"
                "[{1., 2., 3.}, {3., 4., 5.}]]"
                " === Normal[HankelMatrix[NDArray[{1., 2., 3.}], NDArray[{3., 4., 5.}]]]");
    /* rectangular: a length-5 column and a length-3 row -> a 5x3 matrix. */
    assert_true("Normal[Compile[{{c, _Real, 1}, {r, _Real, 1}}, HankelMatrix[c, r]]"
                "[{1., 2., 3., 4., 5.}, {5., 6., 7.}]]"
                " === Normal[HankelMatrix[NDArray[{1., 2., 3., 4., 5.}], NDArray[{5., 6., 7.}]]]");
    /* int + int -> int64 (NDF2_SAME preserves; the reference must be int64). */
    assert_true("Normal[Compile[{{c, _Integer, 1}, {r, _Integer, 1}}, HankelMatrix[c, r]]"
                "[{1, 2, 3}, {3, 4, 5}]]"
                " === Normal[HankelMatrix[NDArray[{1, 2, 3}, DataType -> \"int64\"],"
                " NDArray[{3, 4, 5}, DataType -> \"int64\"]]]");
    /* cliff: Tr of the produced matrix compiles whole. */
    assert_lowers("CompileDiagnostics[{{c, _Real, 1}, {r, _Real, 1}}, Tr[HankelMatrix[c, r]]]");
    /* mixed int/real declines (SAME needs equal dtypes) -> interpreter coerces. */
    assert_true("Lookup[CompileDiagnostics[{{c, _Integer, 1}, {r, _Real, 1}}, HankelMatrix[c, r]],"
                " \"Compiled\"] === False");
    /* a rank-2 first operand is not the two-vector form -> declines. */
    assert_true("Lookup[CompileDiagnostics[{{c, _Real, 2}, {r, _Real, 1}}, HankelMatrix[c, r]],"
                " \"Compiled\"] === False");
}

static void test_cf_toeplitzmatrix_2vec(void) {
    assert_lowers("CompileDiagnostics[{{c, _Real, 1}, {r, _Real, 1}}, ToeplitzMatrix[c, r]]");
    assert_true("Normal[Compile[{{c, _Real, 1}, {r, _Real, 1}}, ToeplitzMatrix[c, r]]"
                "[{1., 2., 3.}, {1., 4., 5.}]]"
                " === Normal[ToeplitzMatrix[NDArray[{1., 2., 3.}], NDArray[{1., 4., 5.}]]]");
    /* rectangular: a length-5 column and a length-3 row -> a 5x3 matrix. */
    assert_true("Normal[Compile[{{c, _Real, 1}, {r, _Real, 1}}, ToeplitzMatrix[c, r]]"
                "[{1., 2., 3., 4., 5.}, {1., 6., 7.}]]"
                " === Normal[ToeplitzMatrix[NDArray[{1., 2., 3., 4., 5.}], NDArray[{1., 6., 7.}]]]");
    assert_true("Normal[Compile[{{c, _Integer, 1}, {r, _Integer, 1}}, ToeplitzMatrix[c, r]]"
                "[{1, 2, 3}, {1, 4, 5}]]"
                " === Normal[ToeplitzMatrix[NDArray[{1, 2, 3}, DataType -> \"int64\"],"
                " NDArray[{1, 4, 5}, DataType -> \"int64\"]]]");
    assert_true("Lookup[CompileDiagnostics[{{c, _Complex, 1}, {r, _Complex, 1}}, ToeplitzMatrix[c, r]],"
                " \"Compiled\"] === False");
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
        assert_true("Normal[Compile[{{c, _Real, 1}, {r, _Real, 1}}, HankelMatrix[c, r]]"
                    "[{1., 2., 3.}, {3., 4., 5.}]]"
                    " === Normal[HankelMatrix[NDArray[{1., 2., 3.}], NDArray[{3., 4., 5.}]]]");
        assert_true("Compile[{{v, _Real, 1}}, Total[Abs[Fourier[v]]]][{1., 2., 3., 4.}]"
                    " === Total[Abs[Fourier[NDArray[{1., 2., 3., 4.}]]]]");
    }
}

/* ------------------------------------------------------------------ *
 *  Norm[v, p] / Norm[m, method]  (two-arg forms + inline-vector expansion)
 * ------------------------------------------------------------------ *
 * Two mechanisms, by operand kind:
 *   - A DECLARED-ARRAY operand delegates to ndla_norm through the V_NORM opcode
 *     with the literal p baked into imm.r; the compiled answer is bit-identical
 *     to Norm[NDArray[...], p] by construction (same delegate).
 *   - A LIST-LITERAL operand is expanded to scalar arithmetic (norm_try_expand),
 *     which is what lets a Norm[{f(x), ...}] body auto-compile inside Table/Plot.
 *     Its value equals the interpreter's within machine tolerance (a different
 *     summation order, exactly as the other autocompiled bodies). */
static void test_cf_norm_array_vector(void) {
    /* every reported two-arg vector form now lowers */
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Norm[v, 1]]");
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Norm[v, Infinity]]");
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Norm[v, 3]]");
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Norm[v, 2.5]]");
    /* bit-identical to the NDArray delegate */
    assert_true("Compile[{{v, _Real, 1}}, Norm[v, 1]][{3., -4., 5.}]"
                " === Norm[NDArray[{3., -4., 5.}], 1]");
    assert_true("Compile[{{v, _Real, 1}}, Norm[v, Infinity]][{3., -4., 5.}]"
                " === Norm[NDArray[{3., -4., 5.}], Infinity]");
    assert_true("Compile[{{v, _Real, 1}}, Norm[v, 3]][{3., -4., 5.}]"
                " === Norm[NDArray[{3., -4., 5.}], 3]");
    /* na==1 (the pre-existing reduction) is untouched */
    assert_true("Compile[{{v, _Real, 1}}, Norm[v]][{3., 4.}]"
                " === Norm[NDArray[{3., 4.}]]");
    /* a rank-1 operand does not accept the string method (ndla declines) */
    assert_true("Lookup[CompileDiagnostics[{{v, _Real, 1}}, Norm[v, \"Frobenius\"]],"
                " \"Compiled\"] === False");
}

static void test_cf_norm_array_matrix(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, Norm[m, \"Frobenius\"]]");
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, Norm[m, 1]]");
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, Norm[m, Infinity]]");
    assert_true("Compile[{{m, _Real, 2}}, Norm[m, \"Frobenius\"]][{{3., 4.}, {0., 0.}}]"
                " === Norm[NDArray[{{3., 4.}, {0., 0.}}], \"Frobenius\"]");
    assert_true("Compile[{{m, _Real, 2}}, Norm[m, 1]][{{1., -2.}, {3., 4.}}]"
                " === Norm[NDArray[{{1., -2.}, {3., 4.}}], 1]");
    assert_true("Compile[{{m, _Real, 2}}, Norm[m, Infinity]][{{1., -2.}, {3., 4.}}]"
                " === Norm[NDArray[{{1., -2.}, {3., 4.}}], Infinity]");
    /* the induced matrix p-norm (p != 1, 2, Inf) has no LAPACK path -> declines */
    assert_true("Lookup[CompileDiagnostics[{{m, _Real, 2}}, Norm[m, 3]],"
                " \"Compiled\"] === False");
}

static void test_cf_norm_list_literal(void) {
    /* the inline-vector case that blocked auto-compilation now lowers */
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Norm[{Sin[x], Cos[x], x}]]");
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Norm[{Sin[x], x}, 1]]");
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Norm[{Sin[x], x}, Infinity]]");
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Norm[{Sin[x], x}, 3]]");
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Norm[{{x, 1.}, {2., x}}, \"Frobenius\"]]");
    /* value equals the interpreter within machine tolerance */
    assert_true("Abs[Compile[{{x, _Real}}, Norm[{Sin[x], Cos[x], x}]][0.7]"
                " - Norm[{Sin[0.7], Cos[0.7], 0.7}]] < 10^-10");
    assert_true("Abs[Compile[{{x, _Real}}, Norm[{Sin[x], x}, Infinity]][0.7]"
                " - Norm[{Sin[0.7], 0.7}, Infinity]] < 10^-10");
    assert_true("Abs[Compile[{{x, _Real}}, Norm[{Sin[x], x}, 3]][0.7]"
                " - Norm[{Sin[0.7], 0.7}, 3]] < 10^-10");
    /* a rank-2 list under an induced (non-Frobenius) norm is not expandable */
    assert_true("Lookup[CompileDiagnostics[{{x, _Real}}, Norm[{{x, 1.}, {2., x}}]],"
                " \"Compiled\"] === False");
}

/* ------------------------------------------------------------------ *
 *  Chop / Clip scalar lowerings (branchless mask-multiply; default bounds)
 * ------------------------------------------------------------------ *
 * Chop[x] = x * (Abs[x] >= delta); the chopped value is a machine 0. (a
 * CompiledFunction returns machine types, the compiled counterpart of the
 * interpreter's exact Integer 0). Clip[x] fills the default-bounds na==1 gap
 * as Min[Max[x, -1.], 1.]. */
static void test_cf_chop(void) {
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Chop[x]]");
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Chop[x, 0.01]]");
    assert_lowers("CompileDiagnostics[{{n, _Integer}}, Chop[n]]");   /* int never chops */
    assert_true("Compile[{{x, _Real}}, Chop[x]][1.0*10^-15] === 0.");
    assert_true("Compile[{{x, _Real}}, Chop[x]][3.5] === 3.5");
    assert_true("Compile[{{x, _Real}}, Chop[x, 0.01]][0.005] === 0.");
    assert_true("Compile[{{x, _Real}}, Chop[x, 0.01]][0.5] === 0.5");
    assert_true("Compile[{{n, _Integer}}, Chop[n]][5] === 5");
    /* a complex operand has no scalar lowering (chop can drop to real) */
    assert_true("Lookup[CompileDiagnostics[{{z, _Complex}}, Chop[z]], \"Compiled\"] === False");
}

static void test_cf_clip_default(void) {
    assert_lowers("CompileDiagnostics[{{x, _Real}}, Clip[x]]");
    assert_true("Compile[{{x, _Real}}, Clip[x]][1.5] === 1.");
    assert_true("Compile[{{x, _Real}}, Clip[x]][-2.] === -1.");
    assert_true("Compile[{{x, _Real}}, Clip[x]][0.5] === 0.5");
    /* the explicit-bounds scalar form still lowers */
    assert_true("Compile[{{x, _Real}}, Clip[x, {0., 2.}]][3.0] === 2.");
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
    TEST(test_cf_hankelmatrix_2vec);
    TEST(test_cf_toeplitzmatrix_2vec);
    TEST(test_cf_vandermondematrix);

    TEST(test_buffer_producers);

    TEST(test_cf_norm_array_vector);
    TEST(test_cf_norm_array_matrix);
    TEST(test_cf_norm_list_literal);

    TEST(test_cf_chop);
    TEST(test_cf_clip_default);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All Compile transforms tests passed!\n");
    symtab_clear();
    return 0;
}

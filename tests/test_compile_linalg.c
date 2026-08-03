/* Unit tests for the COMPILE_MISSING.md §2 / §3 Compile[] lowerings:
 * linear-algebra and structural heads that delegate to the interpreter's own
 * NDArray fast path.
 *
 *   §2 (single array -> array, ND_FNS / A_NDFN):
 *       Inverse, Normalize, MatrixPower, ReverseSort, ConjugateTranspose,
 *       PseudoInverse.
 *   §3 (two arrays -> array or scalar, ND_FN2S / A_NDFN2, V_NDFN2):
 *       Dot (matrix.matrix, matrix.vector, and the vector.vector SCALAR inner
 *       product), LinearSolve, Cross, LeastSquares, ListConvolve, ListCorrelate,
 *       Join.
 *
 * Each head gets, through the eval path:
 *   1. "it lowers":  CompileDiagnostics[spec, H[...]] contains Compiled -> True.
 *   2. "compiled == interpreted":  the compiled value equals the interpreter's
 *      OWN NDArray-path result for the same input.  The reference wraps the
 *      inputs in NDArray[...] so the interpreter takes the SAME ndla_* / BLAS
 *      delegate the compiled code does — bit-identical by construction — and
 *      Normal[] on both sides removes the packed-vs-plain representation
 *      difference (a plain-List reference would instead hit the interpreter's
 *      exact/mixed path: LinearSolve of a plain List gives {1, 2.} where the
 *      machine path gives {1., 2.}).  This tests the LOWERING (the VM rebuilds
 *      the call and round-trips the result), not the kernel, which has its own
 *      tests.
 *   3. "cliff / fusion": a COMPOSED body (Total[...], a.b + 1.0, Tr[Inverse[m]])
 *      compiles whole, which exercises the infer_type branch — without it a
 *      composed body would bail even though the standalone head lowers.
 *   4. "clean decline, no wrong answer": an operand dtype barred by the element
 *      gate reports Compiled -> False and still returns the interpreter's exact
 *      value.
 *   5. A repeated-evaluation loop to surface double-free / leak regressions
 *      under valgrind.
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

/* Evaluate `input` and assert the printed result CONTAINS `substr` (used for the
 * CompileDiagnostics association, whose full text is verbose and version-y). */
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

/* Shorthand: assert a boolean Mathilda expression evaluates to True. */
static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* Shorthand: assert `H[...]` lowers under Compile[]. */
static void assert_lowers(const char* diag_call) {
    assert_eval_contains(diag_call, "Compiled\" -> True");
}

/* ------------------------------------------------------------------ *
 *  §2 single-array heads                                             *
 * ------------------------------------------------------------------ */

static void test_compile_inverse(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, Inverse[m]]");
    assert_true("Normal[Compile[{{m, _Real, 2}}, Inverse[m]][{{4., 7.}, {2., 6.}}]]"
                " === Normal[Inverse[NDArray[{{4., 7.}, {2., 6.}}]]]");
    /* complex operand (gate allows real + complex) */
    assert_lowers("CompileDiagnostics[{{m, _Complex, 2}}, Inverse[m]]");
    assert_true("Normal[Compile[{{m, _Complex, 2}}, Inverse[m]][{{1. + I, 0.}, {0., 2.}}]]"
                " === Normal[Inverse[NDArray[{{1. + I, 0.}, {0., 2.}}]]]");
    /* fused (cliff): Tr of the inverse compiles whole */
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, Tr[Inverse[m]]]");
    assert_true("Compile[{{m, _Real, 2}}, Tr[Inverse[m]]][{{4., 7.}, {2., 6.}}]"
                " === Tr[Inverse[NDArray[{{4., 7.}, {2., 6.}}]]]");
    /* int operand: exact Rationals no machine slot holds -> declines, exact value */
    assert_true("Lookup[CompileDiagnostics[{{m, _Integer, 2}}, Inverse[m]], \"Compiled\"]"
                " === False");
    assert_true("Compile[{{m, _Integer, 2}}, Inverse[m]][{{1, 2}, {3, 4}}]"
                " === Inverse[{{1, 2}, {3, 4}}]");
}

static void test_compile_normalize(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, Normalize[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, Normalize[v]][{3., 4.}]]"
                " === Normal[Normalize[NDArray[{3., 4.}]]]");
    /* rank-2 operand: rank_rule 4 requires rank 1 -> does NOT lower */
    assert_true("Lookup[CompileDiagnostics[{{m, _Real, 2}}, Normalize[m]], \"Compiled\"]"
                " === False");
}

static void test_compile_matrixpower(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}, {k, _Integer}}, MatrixPower[m, k]]");
    assert_true("Normal[Compile[{{m, _Real, 2}}, MatrixPower[m, 3]][{{1., 1.}, {0., 1.}}]]"
                " === Normal[MatrixPower[NDArray[{{1., 1.}, {0., 1.}}], 3]]");
    assert_true("Normal[Compile[{{m, _Real, 2}}, MatrixPower[m, 0]][{{2., 5.}, {1., 3.}}]]"
                " === Normal[MatrixPower[NDArray[{{2., 5.}, {1., 3.}}], 0]]");
}

static void test_compile_reversesort(void) {
    assert_lowers("CompileDiagnostics[{{v, _Real, 1}}, ReverseSort[v]]");
    assert_true("Normal[Compile[{{v, _Real, 1}}, ReverseSort[v]][{3., 1., 2.}]]"
                " === Normal[ReverseSort[NDArray[{3., 1., 2.}]]]");
    /* integer vector preserves dtype (gate int + real) */
    assert_lowers("CompileDiagnostics[{{v, _Integer, 1}}, ReverseSort[v]]");
    assert_true("Normal[Compile[{{v, _Integer, 1}}, ReverseSort[v]][{3, 1, 2}]]"
                " === {3, 2, 1}");
}

static void test_compile_conjugatetranspose(void) {
    assert_lowers("CompileDiagnostics[{{m, _Complex, 2}}, ConjugateTranspose[m]]");
    assert_true("Normal[Compile[{{m, _Complex, 2}}, ConjugateTranspose[m]]"
                "[{{1. + I, 2.}, {3., 4. - I}}]]"
                " === Normal[ConjugateTranspose[NDArray[{{1. + I, 2.}, {3., 4. - I}}]]]");
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, ConjugateTranspose[m]]");
    assert_true("Normal[Compile[{{m, _Real, 2}}, ConjugateTranspose[m]][{{1., 2.}, {3., 4.}}]]"
                " === Normal[ConjugateTranspose[NDArray[{{1., 2.}, {3., 4.}}]]]");
}

static void test_compile_pseudoinverse(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}}, PseudoInverse[m]]");
    assert_true("Normal[Compile[{{m, _Real, 2}}, PseudoInverse[m]]"
                "[{{1., 0.}, {0., 1.}, {1., 1.}}]]"
                " === Normal[PseudoInverse[NDArray[{{1., 0.}, {0., 1.}, {1., 1.}}]]]");
    /* complex is real-only-gated out -> declines */
    assert_true("Lookup[CompileDiagnostics[{{m, _Complex, 2}}, PseudoInverse[m]], \"Compiled\"]"
                " === False");
}

/* ------------------------------------------------------------------ *
 *  §3 two-array heads                                                *
 * ------------------------------------------------------------------ */

static void test_compile_dot_matrix(void) {
    /* matrix . matrix -> A_NDFN2 */
    assert_lowers("CompileDiagnostics[{{a, _Real, 2}, {b, _Real, 2}}, a.b]");
    assert_true("Normal[Compile[{{a, _Real, 2}, {b, _Real, 2}}, a.b]"
                "[{{1., 2.}, {3., 4.}}, {{5., 6.}, {7., 8.}}]]"
                " === Normal[NDArray[{{1., 2.}, {3., 4.}}].NDArray[{{5., 6.}, {7., 8.}}]]");
    /* matrix . vector -> A_NDFN2 (the commonest numeric kernel) */
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}, {v, _Real, 1}}, m.v]");
    assert_true("Compile[{{m, _Real, 2}, {v, _Real, 1}}, m.v]"
                "[{{1., 2.}, {3., 4.}}, {1., 1.}] === {3., 7.}");
    /* complex matrix . vector */
    assert_true("Normal[Compile[{{m, _Complex, 2}, {v, _Complex, 1}}, m.v]"
                "[{{1. + I, 2.}, {3., 4.}}, {1., 1.}]]"
                " === Normal[NDArray[{{1. + I, 2.}, {3., 4.}}].NDArray[{1., 1.}]]");
}

static void test_compile_dot_scalar(void) {
    /* vector . vector -> V_NDFN2 (scalar inner product) */
    assert_lowers("CompileDiagnostics[{{a, _Real, 1}, {b, _Real, 1}}, a.b]");
    assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, a.b][{1., 2., 3.}, {4., 5., 6.}]"
                " === 32.");
    /* fused (cliff): a body containing an inner product must compile WHOLE */
    assert_lowers("CompileDiagnostics[{{a, _Real, 1}, {b, _Real, 1}}, a.b + 1.0]");
    assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, a.b + 1.0]"
                "[{1., 2., 3.}, {4., 5., 6.}] === 33.");
    /* int inner product: exact past 2^53 -> declines, exact value */
    assert_true("Lookup[CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, v.w], \"Compiled\"]"
                " === False");
    assert_true("Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, v.w][{2, 3}, {4, 5}] === 23");
}

static void test_compile_linearsolve(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}, {v, _Real, 1}}, LinearSolve[m, v]]");
    assert_true("Normal[Compile[{{m, _Real, 2}, {v, _Real, 1}}, LinearSolve[m, v]]"
                "[{{2., 0.}, {0., 4.}}, {2., 8.}]]"
                " === Normal[LinearSolve[NDArray[{{2., 0.}, {0., 4.}}], NDArray[{2., 8.}]]]");
    /* matrix rhs -> matrix result (result rank follows rhs rank) */
    assert_true("Normal[Compile[{{m, _Real, 2}, {b, _Real, 2}}, LinearSolve[m, b]]"
                "[{{2., 0.}, {0., 4.}}, {{2., 4.}, {8., 4.}}]]"
                " === Normal[LinearSolve[NDArray[{{2., 0.}, {0., 4.}}],"
                " NDArray[{{2., 4.}, {8., 4.}}]]]");
}

static void test_compile_cross(void) {
    assert_lowers("CompileDiagnostics[{{a, _Real, 1}, {b, _Real, 1}}, Cross[a, b]]");
    assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, Cross[a, b]]"
                "[{1., 0., 0.}, {0., 1., 0.}] === {0., 0., 1.}");
}

static void test_compile_leastsquares(void) {
    assert_lowers("CompileDiagnostics[{{m, _Real, 2}, {v, _Real, 1}}, LeastSquares[m, v]]");
    assert_true("Normal[Compile[{{m, _Real, 2}, {v, _Real, 1}}, LeastSquares[m, v]]"
                "[{{1., 1.}, {1., 2.}, {1., 3.}}, {6., 7., 10.}]]"
                " === Normal[LeastSquares[NDArray[{{1., 1.}, {1., 2.}, {1., 3.}}],"
                " NDArray[{6., 7., 10.}]]]");
}

static void test_compile_listconvolve(void) {
    assert_lowers("CompileDiagnostics[{{k, _Real, 1}, {v, _Real, 1}}, ListConvolve[k, v]]");
    assert_true("Compile[{{k, _Real, 1}, {v, _Real, 1}}, ListConvolve[k, v]]"
                "[{1., 1.}, {1., 2., 3., 4.}] === {3., 5., 7.}");
    assert_lowers("CompileDiagnostics[{{k, _Real, 1}, {v, _Real, 1}}, ListCorrelate[k, v]]");
    assert_true("Compile[{{k, _Real, 1}, {v, _Real, 1}}, ListCorrelate[k, v]]"
                "[{1., 2.}, {1., 2., 3., 4.}] === {5., 8., 11.}");
}

static void test_compile_join(void) {
    assert_lowers("CompileDiagnostics[{{a, _Real, 1}, {b, _Real, 1}}, Join[a, b]]");
    assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, Join[a, b]][{1., 2.}, {3., 4.}]"
                " === {1., 2., 3., 4.}");
    /* integer arrays: Join preserves dtype (SAME rule, gate int too) */
    assert_true("Compile[{{a, _Integer, 1}, {b, _Integer, 1}}, Join[a, b]][{1, 2}, {3, 4}]"
                " === {1, 2, 3, 4}");
    /* rank-2 Join -> rank-2 */
    assert_true("Normal[Compile[{{a, _Real, 2}, {b, _Real, 2}}, Join[a, b]]"
                "[{{1., 2.}}, {{3., 4.}}]] === {{1., 2.}, {3., 4.}}");
}

/* ------------------------------------------------------------------ *
 *  Memory-safety stress loop                                         *
 * ------------------------------------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    for (int t = 0; t < 50; t++) {
        assert_true("Compile[{{m, _Real, 2}, {v, _Real, 1}}, m.v]"
                    "[{{1., 2.}, {3., 4.}}, {1., 1.}] === {3., 7.}");
        assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, a.b][{1., 2., 3.}, {4., 5., 6.}]"
                    " === 32.");
        assert_true("Normal[Compile[{{m, _Real, 2}}, Inverse[m]][{{4., 7.}, {2., 6.}}]]"
                    " === Normal[Inverse[NDArray[{{4., 7.}, {2., 6.}}]]]");
        assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, Join[a, b]][{1., 2.}, {3., 4.}]"
                    " === {1., 2., 3., 4.}");
        assert_true("Compile[{{a, _Real, 1}, {b, _Real, 1}}, Cross[a, b]]"
                    "[{1., 0., 0.}, {0., 1., 0.}] === {0., 0., 1.}");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_compile_inverse);
    TEST(test_compile_normalize);
    TEST(test_compile_matrixpower);
    TEST(test_compile_reversesort);
    TEST(test_compile_conjugatetranspose);
    TEST(test_compile_pseudoinverse);

    TEST(test_compile_dot_matrix);
    TEST(test_compile_dot_scalar);
    TEST(test_compile_linearsolve);
    TEST(test_compile_cross);
    TEST(test_compile_leastsquares);
    TEST(test_compile_listconvolve);
    TEST(test_compile_join);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All Compile linalg tests passed!\n");
    symtab_clear();
    return 0;
}

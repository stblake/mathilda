/* test_dense_struct_fixes.c -- dense structural builtins:
 *
 *   1. KroneckerProduct for List / packed / NDArray operands of any rank
 *      (src/linalg/kronecker.c), next to the native SparseArray kernel;
 *   2. Upper/LowerTriangularize for dense matrices (src/linalg/triangularize.c);
 *   3. Transpose of a vector, and the level-specification messages
 *      (perm1, perm2, newdims, tperm, permspec, nmtx);
 *   4. Reverse never reaching inside a Rational or a Complex;
 *   5. IdentityMatrix::dims;
 *   6. SparseArray with a zero dimension is the empty List;
 *   7. Flatten::flev for a negative level, dense and sparse (and Infinity);
 *   8. Transpose::diagnl for a non-square {1, 1}, dense and sparse.
 *
 * Every expected value was checked against an independent implementation.
 * A message is checked through Check[expr, $Failed] (the funnel's Check
 * surface); its absence through the same Check returning the value.
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
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

/* `input` emits a message (Check fires) and stays unevaluated as `rest`. */
static void chk_msg(const char* input, const char* rest) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Check[%s, $Failed]", input);
    chk(buf, "$Failed");
    snprintf(buf, sizeof buf, "Quiet[%s]", input);
    chk(buf, rest);
}

static void run(const char* input) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    expr_free(e);
}

/* ------------------------------------------------------------ 1. Kronecker */

static void test_kronecker(void) {
    chk("Attributes[KroneckerProduct]", "{Protected}");
    /* two vectors: a column times a row */
    chk("KroneckerProduct[{1, 2}, {3, 4}]", "{{3, 4}, {6, 8}}");
    chk("KroneckerProduct[{{1, 2}, {3, 4}}, {{0, 5}, {6, 7}}]",
        "{{0, 5, 0, 10}, {6, 7, 12, 14}, {0, 15, 0, 20}, {18, 21, 24, 28}}");
    /* a leading vector is a column, a trailing one a row */
    chk("KroneckerProduct[{a, b}, {{1, 2}, {3, 4}}]",
        "{{a, 2 a}, {3 a, 4 a}, {b, 2 b}, {3 b, 4 b}}");
    chk("KroneckerProduct[{{1, 2}, {3, 4}}, {a, b}]", "{{a, b, 2 a, 2 b}, {3 a, 3 b, 4 a, 4 b}}");
    chk("KroneckerProduct[{1, 2}, {3, 4}, {5, 6}]", "{{15, 18, 20, 24}, {30, 36, 40, 48}}");
    chk("KroneckerProduct[{1, 2}, {3, 4}, {5, 6}, {7}]",
        "{{105, 126, 140, 168}, {210, 252, 280, 336}}");
    chk("KroneckerProduct[{{1, 2, 3}}, {{1}, {2}}]", "{{1, 2, 3}, {2, 4, 6}}");
    chk("KroneckerProduct[{1, 2}, {{1, 2}, {3, 4}}]", "{{1, 2}, {3, 4}, {2, 4}, {6, 8}}");
    chk("KroneckerProduct[{1, 2}, {x}]", "{{x}, {2 x}}");
    /* higher ranks */
    chk("KroneckerProduct[{{{1, 2}}}, {{3}, {4}}]", "{{{3}, {4}}, {{6}, {8}}}");
    chk("KroneckerProduct[{{1}}, {{{1, 2}}}]", "{{{1, 2}}}");
    chk("KroneckerProduct[{{1, 2}}, Array[b, {1, 1, 2}]]",
        "{{{b[1, 1, 1], b[1, 1, 2]}, {2 b[1, 1, 1], 2 b[1, 1, 2]}}}");
    chk("Dimensions[KroneckerProduct[Array[a, {2, 3, 4}], Array[b, {5, 6}]]]", "{8, 15, 6}");
    chk("Dimensions[KroneckerProduct[Array[a, {2, 3}], Array[b, {4, 5, 6}]]]", "{8, 15, 6}");
    chk("Dimensions[KroneckerProduct[Array[a, {2, 3, 4}], Array[b, {5, 6, 7}]]]", "{10, 18, 28}");
    chk("Dimensions[KroneckerProduct[Array[a, {2}], Array[b, {3, 4, 5}]]]", "{8, 15}");
    chk("Dimensions[KroneckerProduct[Array[a, {3, 4, 5}], Array[b, {2}]]]", "{15, 8}");
    chk("Dimensions[KroneckerProduct[Array[a, {2}], Array[b, {2, 3, 4, 5}]]]", "{2, 2, 3, 4, 5}");
    chk("Dimensions[KroneckerProduct[Array[a, {2, 3, 4, 5}], Array[b, {6}]]]", "{2, 3, 4, 5, 6}");
    chk("Dimensions[KroneckerProduct[Array[a, {2, 3}], Array[b, {3}]]]", "{2, 9}");
    /* the element order of the generic rule (primes make every product unique) */
    chk("Flatten[KroneckerProduct[ArrayReshape[Prime[Range[24]], {2, 3, 4}], "
        "ArrayReshape[Prime[Range[25, 26]], {2}]]][[{1, 2, 3, 10, 30, 48}]]",
        "{194, 202, 1067, 1313, 7373, 8989}");
    /* three or more: only the last may have rank > 2 */
    chk("Dimensions[KroneckerProduct[{{1}}, {{1}}, Array[c, {2, 2, 2}]]]", "{2, 2, 2}");
    chk("KroneckerProduct[Array[c, {1, 1, 1}], {{1}}, {{1}}]",
        "KroneckerProduct[{{{c[1, 1, 1]}}}, {{1}}, {{1}}]");
    /* empty operands */
    chk("KroneckerProduct[{1, 2}, {}]", "{{}, {}}");
    chk("KroneckerProduct[{{1, 2}}, {}]", "{{}}");
    chk("KroneckerProduct[{}, {}]", "{}");
    /* values: a machine number makes an all-numeric product machine */
    chk("KroneckerProduct[{{1., 2.}}, {{3, 4}}]", "{{3.0, 4.0, 6.0, 8.0}}");
    chk("KroneckerProduct[{1/2, 2}, {3., 1}]", "{{1.5, 0.5}, {6.0, 2.0}}");
    chk("KroneckerProduct[{1, 2.}, {x}]", "{{x}, {2.0 x}}");
    chk("KroneckerProduct[{1, x}, {2.5, 3}]", "{{2.5, 3}, {2.5 x, 3 x}}");
    chk("KroneckerProduct[{{1, 2}}, {{1 + I, 2}}]", "{{1 + I, 2, 2 + 2*I, 4}}");
    chk("KroneckerProduct[{1, 2}, {1.5, I}]",
        "{{1.5 + 0.0*I, 0.0 + 1.0*I}, {3.0 + 0.0*I, 0.0 + 2.0*I}}");
    /* integer products past int64 stay exact */
    chk("KroneckerProduct[{{1, 2}}, {{3, 2^70}}]",
        "{{3, 1180591620717411303424, 6, 2361183241434822606848}}");
    chk("KroneckerProduct[{{2^62, 2}}, {{4, 1}}]",
        "{{18446744073709551616, 4611686018427387904, 8, 2}}");
    /* the packed fast path writes a packed answer */
    run("pm = Partition[Range[4.], 2]; pi = Range[6];");
    chk_true("NDArrayQ[pm] && NDArrayQ[KroneckerProduct[pm, pm]]");
    chk_true("KroneckerProduct[pm, pm] === KroneckerProduct[Normal[pm], Normal[pm]]");
    chk_true("KroneckerProduct[pm, pm] == KroneckerProduct[{{1., 2.}, {3., 4.}}, {{1., 2.}, {3., 4.}}]");
    chk("KroneckerProduct[pi, {1, 2}][[6]]", "{6, 12}");
    chk_true("NDArrayQ[KroneckerProduct[pi, Range[5]]] && "
             "KroneckerProduct[pi, Range[5]] === Outer[Times, Range[6], Range[5]]");
    chk_true("Head[KroneckerProduct[NDArray[{1., 2.}], NDArray[{3., 4.}]]] === NDArray");
    /* the native sparse kernel is unaffected */
    chk("Normal[KroneckerProduct[SparseArray[{1, 2}], {3, 4}]]", "{{3, 4}, {6, 8}}");
    chk("Head[KroneckerProduct[{{1, 2}}, SparseArray[{{3, 4}}]]]", "SparseArray");
    /* messages and unevaluated forms */
    chk_msg("KroneckerProduct[{1, 2}]", "KroneckerProduct[{1, 2}]");
    chk_msg("KroneckerProduct[]", "KroneckerProduct[]");
    chk_msg("KroneckerProduct[{1, 2}, 3]", "KroneckerProduct[{1, 2}, 3]");
    chk_msg("KroneckerProduct[{{1, 2}, {3}}, {1}]", "KroneckerProduct[{{1, 2}, {3}}, {1}]");
    chk("Check[KroneckerProduct[x, y], $Failed]", "KroneckerProduct[x, y]");
    chk("Check[KroneckerProduct[{1, 2}, {3, 4}, x], $Failed]", "KroneckerProduct[{1, 2}, {3, 4}, x]");
    chk("Check[KroneckerProduct[{{1, 2}}, {{1, {2}}}], $Failed]", "KroneckerProduct[{{1, 2}}, {{1, {2}}}]");
}

/* ------------------------------------------------------ 2. triangularize */

static void test_triangularize(void) {
    run("m3 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; r24 = {{1, 2, 3, 4}, {5, 6, 7, 8}};");
    chk("Attributes[UpperTriangularize]", "{Protected}");
    chk("Attributes[LowerTriangularize]", "{Protected}");
    chk("UpperTriangularize[m3]", "{{1, 2, 3}, {0, 5, 6}, {0, 0, 9}}");
    chk("UpperTriangularize[m3, 1]", "{{0, 2, 3}, {0, 0, 6}, {0, 0, 0}}");
    chk("UpperTriangularize[m3, -1]", "{{1, 2, 3}, {4, 5, 6}, {0, 8, 9}}");
    chk("LowerTriangularize[m3]", "{{1, 0, 0}, {4, 5, 0}, {7, 8, 9}}");
    chk("LowerTriangularize[m3, 1]", "{{1, 2, 0}, {4, 5, 6}, {7, 8, 9}}");
    chk("LowerTriangularize[m3, -1]", "{{0, 0, 0}, {4, 0, 0}, {7, 8, 0}}");
    chk("UpperTriangularize[r24]", "{{1, 2, 3, 4}, {0, 6, 7, 8}}");
    chk("LowerTriangularize[r24]", "{{1, 0, 0, 0}, {5, 6, 0, 0}}");
    chk("UpperTriangularize[r24, 2]", "{{0, 0, 3, 4}, {0, 0, 0, 8}}");
    chk("LowerTriangularize[r24, -5]", "{{0, 0, 0, 0}, {0, 0, 0, 0}}");
    chk("UpperTriangularize[m3, 10]", "{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}");
    chk("LowerTriangularize[m3, 10]", "{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}");
    chk("UpperTriangularize[{{1, 2}, {3, 4}}, -2^62]", "{{1, 2}, {3, 4}}");
    chk("UpperTriangularize[{{a, b}, {c, d}}]", "{{a, b}, {0, d}}");
    chk("UpperTriangularize[{{}}]", "{{}}");
    /* the zero: 0. when the matrix holds a machine number anywhere */
    chk("UpperTriangularize[{{1, 2}, {3, 4.}}]", "{{1, 2}, {0.0, 4.0}}");
    chk("UpperTriangularize[{{1.5, I}, {2, 3}}]", "{{1.5, I}, {0.0, 3}}");
    chk("UpperTriangularize[{{f[1.5], 2}, {3, 4}}]", "{{f[1.5], 2}, {0.0, 4}}");
    chk("LowerTriangularize[{{1, x}, {y, 2.5}}]", "{{1, 0.0}, {y, 2.5}}");
    chk("UpperTriangularize[{{1/2, 2}, {3, 4}}]", "{{1/2, 2}, {0, 4}}");
    chk_true("Head[UpperTriangularize[{{1, N[2, 30]}, {3, 4}}][[2, 1]]] === Integer");
    /* buffers: same dtype, packed in and out */
    run("pr = Partition[Range[9.], 3]; pz = Partition[Range[12], 4];");
    chk_true("NDArrayQ[pr] && NDArrayQ[UpperTriangularize[pr]]");
    chk("UpperTriangularize[pr]", "{{1.0, 2.0, 3.0}, {0.0, 5.0, 6.0}, {0.0, 0.0, 9.0}}");
    chk("LowerTriangularize[pr, 1]", "{{1.0, 2.0, 0.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}");
    chk_true("NDArrayQ[pz] && UpperTriangularize[pz, 1] === UpperTriangularize[Normal[pz], 1]");
    chk("LowerTriangularize[pz, -1]", "{{0, 0, 0, 0}, {5, 0, 0, 0}, {9, 10, 0, 0}}");
    chk_true("Head[UpperTriangularize[NDArray[{{1., 2.}, {3., 4.}}]]] === NDArray");
    /* the sparse kernel still answers a SparseArray */
    chk("Head[UpperTriangularize[SparseArray[m3]]]", "SparseArray");
    /* messages */
    chk_msg("UpperTriangularize[{1, 2, 3}]", "UpperTriangularize[{1, 2, 3}]");
    chk_msg("UpperTriangularize[x]", "UpperTriangularize[x]");
    chk_msg("UpperTriangularize[{{1, 2}, {3}}]", "UpperTriangularize[{{1, 2}, {3}}]");
    chk_msg("UpperTriangularize[{{1, {2}}, {3, 4}}]", "UpperTriangularize[{{1, {2}}, {3, 4}}]");
    chk_msg("UpperTriangularize[m3, 1.5]", "UpperTriangularize[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1.5]");
    chk_msg("LowerTriangularize[{{1, 2}, {3, 4}}, 2^70]",
            "LowerTriangularize[{{1, 2}, {3, 4}}, 1180591620717411303424]");
    chk_msg("UpperTriangularize[m3, 1, 2]", "UpperTriangularize[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1, 2]");
    chk_msg("UpperTriangularize[]", "UpperTriangularize[]");
    chk_msg("UpperTriangularize[Range[5.]]", "UpperTriangularize[{1.0, 2.0, 3.0, 4.0, 5.0}]");
    chk("Check[UpperTriangularize[{}], $Failed]", "UpperTriangularize[{}]");
    chk("Check[UpperTriangularize[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}], $Failed]",
        "UpperTriangularize[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]");
}

/* ---------------------------------------------------- 3. vector transpose */

static void test_transpose(void) {
    chk("Transpose[{1, 2, 3}]", "{1, 2, 3}");
    chk("Transpose[{1, 2, 3}, {1}]", "{1, 2, 3}");
    chk("Transpose[{}]", "{}");
    chk("Transpose[{{}, {}}]", "{}");
    chk_true("NDArrayQ[Transpose[Range[10.]]] && Transpose[Range[10.]] === Range[10.]");
    chk("Transpose[{{1, 2, 3}, {4, 5, 6}}, {1}]", "{{1, 2, 3}, {4, 5, 6}}");
    chk("Transpose[{{1, 2, 3}, {4, 5, 6}}, {1, 2, 3}]", "{{1, 2, 3}, {4, 5, 6}}");
    chk("Transpose[{{1, 2}, {4, 5}}, {1, 1}]", "{1, 5}");
    chk("Transpose[Array[a, {2, 3, 2}], {1, 2, 1}]",
        "{{a[1, 1, 1], a[1, 2, 1], a[1, 3, 1]}, {a[2, 1, 2], a[2, 2, 2], a[2, 3, 2]}}");
    chk_msg("Transpose[{1, 2, 3}, {2}]", "Transpose[{1, 2, 3}, {2}]");
    chk_msg("Transpose[{1, 2, 3}, 1]", "Transpose[{1, 2, 3}, 1]");
    chk_msg("Transpose[{1, 2, 3}, {1, 1}]", "Transpose[{1, 2, 3}, {1, 1}]");
    chk_msg("Transpose[{a, {b}}]", "Transpose[{a, {b}}]");
    chk_msg("Transpose[{1, 2}, {1.5}]", "Transpose[{1, 2}, {1.5}]");
    chk_msg("Transpose[{{1, 2, 3}, {4, 5, 6}}, {2, 2}]", "Transpose[{{1, 2, 3}, {4, 5, 6}}, {2, 2}]");
    chk_msg("Transpose[{{1, 2, 3}, {4, 5, 6}}, {3, 1}]", "Transpose[{{1, 2, 3}, {4, 5, 6}}, {3, 1}]");
    chk_msg("Transpose[{{1, 2, 3}, {4, 5, 6}}, {0, 1}]", "Transpose[{{1, 2, 3}, {4, 5, 6}}, {0, 1}]");
    chk_msg("Transpose[{{1, 2}, {3, 4}}, {1, 1, 2}]", "Transpose[{{1, 2}, {3, 4}}, {1, 1, 2}]");
    chk("Check[Transpose[x], $Failed]", "Transpose[x]");
}

/* ------------------------------------------------------------- 4. Reverse */

static void test_reverse(void) {
    chk("Reverse[{5/4, 7/2}, 2]", "{5/4, 7/2}");
    chk("Reverse[{5/4, 7/2}, {1, 2}]", "{7/2, 5/4}");
    chk("Reverse[{1 + 2 I, 3}, 2]", "{1 + 2*I, 3}");
    chk("Reverse[{{5/4, 7/2}, {1, 2}}, 2]", "{{7/2, 5/4}, {2, 1}}");
    chk("Reverse[{{5/4, 7/2}, {1, 2}}, 3]", "{{5/4, 7/2}, {1, 2}}");
    chk("Reverse[{{5/4, 7/2}, {1, 2}}, {1, 3}]", "{{1, 2}, {5/4, 7/2}}");
    chk("Reverse[{f[1, 2], g[3, 4]}, 2]", "{f[2, 1], g[4, 3]}");
    chk("Reverse[f[1/2, 3], 2]", "f[1/2, 3]");
    chk("Reverse[{{1/2}}, {2, 3}]", "{{1/2}}");
    chk("Reverse[{1.5, 2}, 2]", "{1.5, 2}");
}

/* ------------------------------------------------------- 5. IdentityMatrix */

static void test_identity(void) {
    chk_msg("IdentityMatrix[0]", "IdentityMatrix[0]");
    chk_msg("IdentityMatrix[2.5]", "IdentityMatrix[2.5]");
    chk_msg("IdentityMatrix[-1]", "IdentityMatrix[-1]");
    chk_msg("IdentityMatrix[5/2]", "IdentityMatrix[5/2]");
    chk_msg("IdentityMatrix[2^70]", "IdentityMatrix[1180591620717411303424]");
    chk_msg("IdentityMatrix[1 + I]", "IdentityMatrix[1 + I]");
    chk("Check[IdentityMatrix[x], $Failed]", "IdentityMatrix[x]");
    chk("Check[IdentityMatrix[{0, 3}], $Failed]", "IdentityMatrix[{0, 3}]");
    chk("Check[IdentityMatrix[{2, 0}], $Failed]", "IdentityMatrix[{2, 0}]");
    chk("Check[IdentityMatrix[{2, -1}], $Failed]", "IdentityMatrix[{2, -1}]");
    chk("IdentityMatrix[{2, 3}]", "{{1, 0, 0}, {0, 1, 0}}");
    chk("IdentityMatrix[2]", "{{1, 0}, {0, 1}}");
}

/* ------------------------------------------------- 6. zero-size SparseArray */

static void test_sparse_empty(void) {
    chk("Check[SparseArray[{}, {0}], $Failed]", "{}");
    chk("Head[SparseArray[{}, {0}]]", "List");
    chk("SparseArray[{}, {0, 3}]", "{}");
    chk("SparseArray[{}, {3, 0}]", "{{}, {}, {}}");
    chk("SparseArray[{}, {2, 0, 3}]", "{{}, {}}");
    chk("SparseArray[{}, {0}, 5]", "{}");
    chk("SparseArray[{}, 0]", "{}");
    chk("Check[SparseArray[{1 -> 1}, {0}], $Failed]", "{}");
    chk("Head[SparseArray[{}, {3}]]", "SparseArray");
    chk_msg("SparseArray[{}]", "SparseArray[{}]");
}

/* ----------------------------------------------------------- 7. Flatten */

static void test_flatten(void) {
    chk_msg("Flatten[{{1, 2}, {3}}, -1]", "Flatten[{{1, 2}, {3}}, -1]");
    chk_msg("Flatten[{{1, 2}, {3}}, -2]", "Flatten[{{1, 2}, {3}}, -2]");
    chk_msg("Flatten[Range[6], -1]", "Flatten[{1, 2, 3, 4, 5, 6}, -1]");
    run("sf = SparseArray[{{1, 2}, {3, 4}}];");
    chk("Check[Flatten[sf, -1], $Failed]", "$Failed");
    chk("Head[Quiet[Flatten[sf, -1]]]", "Flatten");
    chk("Head[Quiet[Flatten[sf, -1]][[1]]]", "SparseArray");
    chk("Flatten[{{1, 2}, {3}}, Infinity]", "{1, 2, 3}");
    chk("Flatten[{{1, {2}}, {3}}, Infinity]", "{1, 2, 3}");
    chk("Flatten[{{1, 2}, {3}}, 0]", "{{1, 2}, {3}}");
    /* a symbolic level is Flatten::flev, unevaluated (audit Flatten#4) */
    chk("Check[Flatten[{{1, 2}, {3}}, x], $Failed]", "$Failed");
    chk("Quiet[Flatten[{{1, 2}, {3}}, x]]", "Flatten[{{1, 2}, {3}}, x]");
    chk("Normal[Flatten[sf, Infinity]]", "{1, 2, 3, 4}");
}

/* ------------------------------------------------ 8. Transpose::diagnl */

static void test_diagnl(void) {
    chk_msg("Transpose[{{1, 2, 3}, {4, 5, 6}}, {1, 1}]", "Transpose[{{1, 2, 3}, {4, 5, 6}}, {1, 1}]");
    chk_msg("Transpose[N[{{1, 2, 3}, {4, 5, 6}}], {1, 1}]",
            "Transpose[{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, {1, 1}]");
    run("st = SparseArray[{{1, 2, 3}, {4, 5, 6}}];");
    chk("Check[Transpose[st, {1, 1}], $Failed]", "$Failed");
    chk("Head[Quiet[Transpose[st, {1, 1}]][[1]]]", "SparseArray");
    chk("Normal[Transpose[SparseArray[{{1, 2}, {4, 5}}], {1, 1}]]", "{1, 5}");
}

int main(void) {
    symtab_init();
    core_init();
    test_kronecker();
    test_triangularize();
    test_transpose();
    test_reverse();
    test_identity();
    test_sparse_empty();
    test_flatten();
    test_diagnl();
    if (failures) {
        fprintf(stderr, "\n%d dense structural test(s) FAILED\n", failures);
        return 1;
    }
    printf("All dense structural tests passed.\n");
    return 0;
}

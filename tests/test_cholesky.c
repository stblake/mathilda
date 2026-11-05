/* test_cholesky.c -- CholeskyDecomposition[m] (src/linalg/cholesky.c).
 *
 * Exact answers are pinned to the reference's printed form; inexact ones are
 * checked by the invariant ConjugateTranspose[c] . c == m and by structure
 * (upper triangular, positive real diagonal).  Every precision class (exact,
 * symbolic, machine real/complex, MPFR real/complex), every storage (List,
 * packed, NDArray, SparseArray), the messages, and the Compile lowering with
 * bit-identical results.
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_dense_heads.h"

/* Residual of the factorisation, relative to the matrix: tolerance `tol`. */
static void chk_factor(const char* m, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{mm = %s, c}, c = CholeskyDecomposition[mm]; "
             "Max[Abs[ConjugateTranspose[c] . c - mm]] <= %s * Max[Abs[mm]] && "
             "Max[Abs[LowerTriangularize[c, -1]]] == 0]", m, tol);
    chk_true(buf);
}

/* reference answers (wolframscript) for exact input */
static void test_exact(void) {
    chk_if("CholeskyDecomposition[{{4,2},{2,3}}]", "{{2, 1}, {0, Sqrt[2]}}");
    chk_if("CholeskyDecomposition[{{2,1},{1,2}}]", "{{Sqrt[2], 1/Sqrt[2]}, {0, Sqrt[3/2]}}");
    chk_if("CholeskyDecomposition[{{2, I},{-I, 2}}]", "{{Sqrt[2], I/Sqrt[2]}, {0, Sqrt[3/2]}}");
    chk_if("CholeskyDecomposition[{{5}}]", "{{Sqrt[5]}}");
    chk_if("CholeskyDecomposition[{{4,2},{2,3}}/3]", "{{2/Sqrt[3], 1/Sqrt[3]}, {0, Sqrt[2/3]}}");
    chk_if("CholeskyDecomposition[{{4,2,1},{2,3,0},{1,0,5}}]",
           "{{2, 1, 1/2}, {0, Sqrt[2], -1/2*1/Sqrt[2]}, {0, 0, Sqrt[37/2]/2}}");
    chk_true("With[{m = {{4,2,1},{2,3,0},{1,0,5}}}, "
             "Simplify[ConjugateTranspose[CholeskyDecomposition[m]] . CholeskyDecomposition[m] - m] "
             "== {{0,0,0},{0,0,0},{0,0,0}}]");
    chk_true("With[{m = {{2, I, 0},{-I, 3, 1 + I},{0, 1 - I, 4}}}, "
             "Simplify[ConjugateTranspose[CholeskyDecomposition[m]] . CholeskyDecomposition[m] - m] "
             "== ConstantArray[0, {3, 3}]]");
    /* a Hilbert matrix: exact Sqrt of rationals throughout */
    chk_true("With[{m = HilbertMatrix[5]}, "
             "Simplify[Transpose[CholeskyDecomposition[m]] . CholeskyDecomposition[m] - m] "
             "== ConstantArray[0, {5, 5}]]");
}

static void test_symbolic(void) {
    chk_if("CholeskyDecomposition[{{a,b},{b,c}}]",
           "{{Sqrt[a], b/Sqrt[a]}, {0, Sqrt[c - (b*Conjugate[b/Sqrt[a]])/Sqrt[a]]}}");
    /* not Hermitian (Conjugate[a] != a), and a negative numeric pivot is
     * silently unevaluated for a symbolic matrix, as in the reference */
    chk_msg("CholeskyDecomposition[{{a,b},{Conjugate[b],c}}]", "CholeskyDecomposition::herm");
    chk_if("Check[CholeskyDecomposition[{{-1,0},{0,a}}], $Failed]", "CholeskyDecomposition[{{-1, 0}, {0, a}}]");
    chk_if("Check[CholeskyDecomposition[{{a,0},{0,-1}}], $Failed]", "CholeskyDecomposition[{{a, 0}, {0, -1}}]");
}

static void test_machine(void) {
    chk_if("CholeskyDecomposition[{{4.,2.},{2.,3.}}]", "{{2.0, 1.0}, {0.0, 1.4142135623730951}}");
    chk_if("CholeskyDecomposition[{{4,2.},{2,3}}]", "{{2.0, 1.0}, {0.0, 1.4142135623730951}}");
    chk_true("Chop[CholeskyDecomposition[{{2., I},{-I, 2.}}] - "
             "{{1.4142135623730951, 0.7071067811865475 I}, {0, 1.224744871391589}}] == "
             "{{0, 0}, {0, 0}}");
    chk_factor("{{4., 2., 1.}, {2., 3., 0.5}, {1., 0.5, 5.}}", "10.^-15");
    chk_factor("{{4., 1. + I, 0.}, {1. - I, 3., 2. I}, {0., -2. I, 6.}}", "10.^-15");
    chk_factor("Module[{r = RandomReal[{-1, 1}, {60, 60}]}, r . Transpose[r] + 60 IdentityMatrix[60]]",
               "10.^-13");
    chk_factor("Module[{r = RandomComplex[{-1 - I, 1 + I}, {40, 40}]}, "
               "r . ConjugateTranspose[r] + 40 IdentityMatrix[40]]", "10.^-13");
    chk_factor("{{2.5}}", "10.^-15");
    /* the factor of a machine matrix is packed */
    chk_true("NDArrayQ[CholeskyDecomposition[{{4.,2.},{2.,3.}}]]");
}

static void test_mpfr(void) {
    chk_true("With[{c = CholeskyDecomposition[N[{{4,2},{2,3}}, 30]]}, "
             "Precision[c[[2,2]]] > 29 && c[[2,1]] === 0 && "
             "Abs[c[[2,2]] - Sqrt[2]] < 10^-28]");
    chk_true("With[{c = CholeskyDecomposition[N[{{2, I},{-I, 2}}, 30]]}, "
             "Precision[c[[1,2]]] > 29 && Abs[c[[1,2]] - I/Sqrt[2]] < 10^-28 && "
             "Abs[c[[2,2]] - Sqrt[3/2]] < 10^-28]");
    chk_true("With[{m = N[HilbertMatrix[6], 40]}, With[{c = CholeskyDecomposition[m]}, "
             "Max[Abs[Transpose[c] . c - m]] < 10^-35 && Precision[c[[6,6]]] > 39]]");
    chk_msg("CholeskyDecomposition[N[{{1,2},{2,1}}, 30]]", "CholeskyDecomposition::posdef");
}

static void test_storage(void) {
    chk("a = {{4.,2.,1.},{2.,3.,0.5},{1.,0.5,5.}}; Null", "Null");
    chk_true("Head[CholeskyDecomposition[NDArray[a]]] === NDArray");
    chk_true("Normal[CholeskyDecomposition[NDArray[a]]] === CholeskyDecomposition[a]");
    chk_true("CholeskyDecomposition[Developer`ToPackedArray[a]] === CholeskyDecomposition[a]");
    chk_if("CholeskyDecomposition[SparseArray[{{4,2},{2,3}}]]", "{{2, 1}, {0, Sqrt[2]}}");
    chk_true("CholeskyDecomposition[SparseArray[a]] === CholeskyDecomposition[a]");
}

static void test_messages(void) {
    chk_msg("CholeskyDecomposition[{{1,2},{2,1}}]", "CholeskyDecomposition::npdef");
    chk_msg("CholeskyDecomposition[{{0,0},{0,0}}]", "CholeskyDecomposition::npdef");
    chk_msg("CholeskyDecomposition[{{1.,2.},{2.,1.}}]", "CholeskyDecomposition::posdef");
    chk_msg("CholeskyDecomposition[{{1.,0.},{0.,-1.}}]", "CholeskyDecomposition::posdef");
    chk_msg("CholeskyDecomposition[NDArray[{{1.,2.},{2.,1.}}]]", "CholeskyDecomposition::posdef");
    chk_msg("CholeskyDecomposition[{{1,2},{3,4}}]", "CholeskyDecomposition::herm");
    chk_msg("CholeskyDecomposition[{{2, I},{I, 2}}]", "CholeskyDecomposition::herm");
    chk_msg("CholeskyDecomposition[{{2., I},{I, 2.}}]", "CholeskyDecomposition::herm");
    chk_msg("CholeskyDecomposition[{{1.,2.},{2.000001,5.}}]", "CholeskyDecomposition::herm");
    chk_msg("CholeskyDecomposition[{{1,2,3},{4,5,6}}]", "CholeskyDecomposition::matsq");
    chk_msg("CholeskyDecomposition[{}]", "CholeskyDecomposition::matsq");
    chk_msg("CholeskyDecomposition[{{}}]", "CholeskyDecomposition::matsq");
    chk_msg("CholeskyDecomposition[{{4,2},{2,3}}, Method -> Automatic]", "CholeskyDecomposition::optx");
    /* a vector or an atom: silently unevaluated */
    chk_if("Check[CholeskyDecomposition[{1,2}], $Failed]", "CholeskyDecomposition[{1, 2}]");
    chk_if("Check[CholeskyDecomposition[x], $Failed]", "CholeskyDecomposition[x]");
    chk_if("CholeskyDecomposition[{{4,2},{2,3}}, TargetStructure -> \"Dense\"]", "{{2, 1}, {0, Sqrt[2]}}");
    chk_if("Options[CholeskyDecomposition]", "{TargetStructure -> Automatic}");
    chk_if("Attributes[CholeskyDecomposition]", "{Protected}");
}

static void test_compile(void) {
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, CholeskyDecomposition[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, CholeskyDecomposition[m]]");
    chk("cf = Compile[{{m, _Real, 2}}, CholeskyDecomposition[m]]; "
        "cz = Compile[{{m, _Complex, 2}}, CholeskyDecomposition[m]]; Null", "Null");
    chk("a = {{4.,2.,1.},{2.,3.,0.5},{1.,0.5,5.}}; "
        "z = {{4., 1. + I, 0.}, {1. - I, 3., 2. I}, {0., -2. I, 6.}}; Null", "Null");
    /* bit-identical: List == packed == NDArray == Compile */
    chk_true("Normal[cf[a]] === CholeskyDecomposition[a]");
    chk_true("Normal[cf[a]] === Normal[CholeskyDecomposition[NDArray[a]]]");
    chk_true("Normal[cz[z]] === CholeskyDecomposition[z]");
    /* a complex factor is a packed complex array, as the reference's (v0.6221) */
    chk_true("Developer`PackedArrayQ[CholeskyDecomposition[Developer`ToPackedArray[z]]]");
    /* a non-positive-definite operand degrades to the interpreter */
    chk_msg("cf[{{1.,2.},{2.,1.}}]", "CholeskyDecomposition::posdef");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_exact();
    test_symbolic();
    test_machine();
    test_mpfr();
    test_storage();
    test_messages();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d CholeskyDecomposition test(s) FAILED\n", failures);
        return 1;
    }
    printf("All CholeskyDecomposition tests passed.\n");
    return 0;
}

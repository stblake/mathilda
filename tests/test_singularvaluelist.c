/* test_singularvaluelist.c -- SingularValueList (src/linalg/singularvaluelist.c).
 *
 * Exact answers pinned to the reference's printed form; inexact ones checked
 * against the reference digits to a tolerance and against
 * SingularValueDecomposition.  Covers every precision class, the sequence
 * specifications, Tolerance, the generalized form, storage (List, packed,
 * NDArray, SparseArray), messages and the Compile lowering. */

#include "test_dense_heads.h"

static void test_exact(void) {
    chk_if("SingularValueList[{{1,2},{3,4}}]", "{Sqrt[15 + Sqrt[221]], Sqrt[15 - Sqrt[221]]}");
    chk_if("SingularValueList[{{1,2},{2,4}}]", "{5}");
    chk_if("SingularValueList[{{1,2,3},{4,5,6}}]",
           "{Sqrt[(91 + Sqrt[8065])/2], Sqrt[(91 - Sqrt[8065])/2]}");
    chk_if("SingularValueList[{{1, I},{3, 4}}]",
           "{Sqrt[(27 + Sqrt[629])/2], Sqrt[(27 - Sqrt[629])/2]}");
    chk_if("SingularValueList[{{0,0},{0,0}}]", "{0, 0}");
    chk_if("SingularValueList[{{5}}]", "{5}");
    chk_if("SingularValueList[{{1,1},{1,1}}/3]", "{2/3}");
    chk_if("SingularValueList[{{1,2},{3,4}}, 1]", "{Sqrt[15 + Sqrt[221]]}");
    chk_if("SingularValueList[{{1,2},{2,4}}, 2]", "{5, 0}");
    chk_if("SingularValueList[{{1,2},{2,4}}, UpTo[2]]", "{5, 0}");
    chk_if("SingularValueList[{{1.,2.},{2.,4.}}, UpTo[3]]", "{5.0, 0.0}");
    chk_if("SingularValueList[{{1,2},{3,4}}, Tolerance -> 0.5]", "{Sqrt[15 + Sqrt[221]]}");
    chk_if("SingularValueList[{{1,2},{3,4}}, Cubics -> True]", "{Sqrt[15 + Sqrt[221]], Sqrt[15 - Sqrt[221]]}");
    chk_if("SingularValueList[{{a,0},{0,b}}]", "{Sqrt[a*Conjugate[a]], Sqrt[b*Conjugate[b]]}");
}

static void test_machine(void) {
    chk_true("Max[Abs[SingularValueList[{{1.,2.},{3.,4.}}] - {5.464985704219042, 0.3659661906262577}]] < 10.^-14");
    chk_if("SingularValueList[{{1.,2.},{2.,4.}}]", "{5.0}");
    chk_true("Max[Abs[SingularValueList[{{1.,2.,3.},{4.,5.,6.}}] - {9.508032000695719, 0.7728696356734839}]] < 10.^-13");
    chk_true("Max[Abs[SingularValueList[{{1.,I},{3.,4.}}] - {5.102934077957939, 0.9798284523402799}]] < 10.^-13");
    chk_if("SingularValueList[{{0.,0.},{0.,0.}}]", "{0.0, 0.0}");
    chk_if("SingularValueList[{{1.,0.},{0.,1.*^-20}}]", "{1.0}");
    chk_if("SingularValueList[{{1.,0.},{0.,1.*^-3}}, Tolerance -> 0.01]", "{1.0}");
    chk_if("SingularValueList[{{1.,0.},{0.,1.*^-3}}, Tolerance -> 0]", "{1.0, 0.001}");
    chk_if("SingularValueList[{{10.,0},{0,1.}}, Tolerance -> 0.09]", "{10.0, 1.0}");
    chk_if("SingularValueList[{{10.,0},{0,1.}}, Tolerance -> 0.11]", "{10.0}");
    chk_if("SingularValueList[{{10.,0,0},{0,1.,0},{0,0,0}}]", "{10.0, 1.0}");
    chk_if("SingularValueList[{{10.,0,0},{0,1.,0},{0,0,0}}, 3]", "{10.0, 1.0, 0.0}");
    chk_if("SingularValueList[{{10.,0,0},{0,1.,0},{0,0,0}}, -1]", "{0.0}");
    chk_if("SingularValueList[{{1.,2.},{3.,4.}}, 0]", "{}");
    chk_true("SingularValueList[{{1.,2.},{3.,4.}}, {1}] == Take[SingularValueList[{{1.,2.},{3.,4.}}], 1]");
    chk_true("SingularValueList[{{1.,2.},{3.,4.}}, {1, 2}] == SingularValueList[{{1.,2.},{3.,4.}}]");
    chk_true("SingularValueList[{{1.,2.},{3.,4.}}, UpTo[3]] == SingularValueList[{{1.,2.},{3.,4.}}]");
    chk_true("SingularValueList[{{1.,2.},{3.,4.}}, -1] == {Last[SingularValueList[{{1.,2.},{3.,4.}}]]}");
    /* agreement with the full decomposition, real and complex, tall and wide */
    chk_true("Module[{r = RandomReal[{-1, 1}, {30, 20}]}, "
             "Max[Abs[SingularValueList[r] - Diagonal[SingularValueDecomposition[r][[2]]]]] < 10.^-12]");
    chk_true("Module[{r = RandomComplex[{-1 - I, 1 + I}, {12, 25}]}, "
             "Max[Abs[SingularValueList[r] - Diagonal[SingularValueDecomposition[r][[2]]]]] < 10.^-12]");
}

static void test_mpfr(void) {
    chk_true("With[{s = SingularValueList[N[{{1,2},{3,4}}, 30]]}, Precision[s[[1]]] > 29 && "
             "Abs[s[[1]] - Sqrt[15 + Sqrt[221]]] < 10^-28 && Abs[s[[2]] - Sqrt[15 - Sqrt[221]]] < 10^-28]");
    chk_true("With[{s = SingularValueList[N[{{1,I},{3,4}}, 30]]}, Precision[s[[2]]] > 29 && "
             "Abs[s[[2]] - Sqrt[(27 - Sqrt[629])/2]] < 10^-28]");
    chk_true("Length[SingularValueList[N[{{1,2},{2,4}}, 30]]] == 1");
    chk_true("With[{s = SingularValueList[N[{{1,2,3},{4,5,6}}, 40]]}, "
             "Abs[s[[1]] - Sqrt[(91 + Sqrt[8065])/2]] < 10^-38]");
}

/* The real MPFR kernel: Householder bidiagonalisation + bidiagonal QR, with
 * the working precision widened by the spread of the values.  Every expected
 * value here is exact or an independently known closed form. */
static void test_mpfr_bidiag(void) {
    /* Graded, exactly representable: an 8 x 4 matrix with orthonormal
     * columns (I - J/4 is orthogonal) times diag(1, 2^-30, 2^-80, 2^-130)
     * times a 4 x 4 orthogonal Hadamard matrix / 2.  Every entry is dyadic
     * and fits in 60 digits, so the singular values of the MPFR matrix are
     * exactly the diagonal; each must come out to full relative precision
     * although the spread is 2^130. */
    chk_true("Module[{q = IdentityMatrix[8] - ConstantArray[1/4, {8, 8}], "
             "h = {{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}}/2, "
             "d = {1, 2^-30, 2^-80, 2^-130}, s}, "
             "s = SingularValueList[N[q[[All, 1 ;; 4]] . DiagonalMatrix[d] . h, 60]]; "
             "Length[s] == 4 && Max[Abs[s/d - 1]] < 10^-55 && Min[Precision /@ s] >= 59]");
    /* the same, wide (transposed) */
    chk_true("Module[{q = IdentityMatrix[8] - ConstantArray[1/4, {8, 8}], "
             "h = {{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}}/2, "
             "d = {1, 2^-30, 2^-80, 2^-130}, s}, "
             "s = SingularValueList[N[Transpose[q[[All, 1 ;; 4]] . DiagonalMatrix[d] . h], 60]]; "
             "Max[Abs[s/d - 1]] < 10^-55]");
    /* the n x n bidiagonal of ones: sigma_k = 2 Cos[k Pi/(2 n + 1)] */
    chk_true("Module[{n = 30, s}, "
             "s = SingularValueList[N[Table[If[j == i || j == i + 1, 1, 0], {i, n}, {j, n}], 50]]; "
             "Max[Abs[s - Table[2 Cos[k Pi/(2 n + 1)], {k, n}]]] < 10^-47]");
    /* the 50 x 50 Vandermonde matrix on (i-1)/49: smallest value known to
     * 28 digits (two independent algorithms at 400 and 600 digits) */
    chk_true("Module[{s = SingularValueList[N[Table[If[j == 1, 1, ((i - 1)/49)^(j - 1)], "
             "{i, 50}, {j, 50}], 100]]}, "
             "Length[s] == 50 && Abs[Last[s]/1.214052974061843949823397917`30*^-42 - 1] < 10^-26]");
    /* Kahan matrix, n = 120, theta = 6/5 */
    chk_true("Module[{s = SingularValueList[N[Table[Which[i == j, Sin[6/5]^(i - 1), "
             "j > i, -Cos[6/5] Sin[6/5]^(i - 1), True, 0], {i, 120}, {j, 120}], 60]]}, "
             "Abs[Last[s]/4.488013531293049259832894189`30*^-20 - 1] < 10^-26]");
    /* rank deficiency and zero diagonal entries of the bidiagonal */
    chk_true("With[{s = SingularValueList[N[{{1,2},{2,4},{3,6}}, 30]]}, "
             "Length[s] == 1 && Abs[s[[1]] - Sqrt[70]] < 10^-28]");
    chk_true("With[{s = SingularValueList[N[{{0,1,0},{0,0,1},{0,0,0}}, 30], 3]}, "
             "Abs[s[[1]] - 1] < 10^-29 && Abs[s[[2]] - 1] < 10^-29 && s[[3]] == 0]");
    chk_true("With[{s = SingularValueList[N[{{1,1,0},{0,0,1},{0,0,1}}, 30]]}, "
             "Max[Abs[s - {Sqrt[2], Sqrt[2]}]] < 10^-28]");
    chk_if("SingularValueList[N[{{0,0},{0,0},{0,0}}, 30]]", "{0, 0}");
    /* agreement with the eigenvalues of the Gram matrix at twice the digits */
    chk_true("Module[{r = RandomReal[{-1, 1}, {9, 6}, WorkingPrecision -> 40], s, g}, "
             "s = SingularValueList[r]; "
             "g = Sqrt[Reverse[Sort[Eigenvalues[SetPrecision[Transpose[r] . r, 90]]]]]; "
             "Max[Abs[s - g]] < 10^-37]");
}

static void test_generalized(void) {
    chk_true("Max[Abs[SingularValueList[{{{1.,2.},{3.,4.}}, {{2.,0.},{0.,1.}}}] - "
             "{4.738720018687269, 0.21102744961856199}]] < 10.^-13");
    chk_true("Max[Abs[SingularValueList[{{{1.,2.},{3.,4.}}, {{2.,0.},{0.,1.}}}, 1] - {4.738720018687269}]] < 10.^-13");
    chk_true("Max[Abs[SingularValueList[{{{1.,2.,3.},{3.,4.,1.}}, {{2.,0.,1.},{0.,1.,1.},{1.,1.,1.}}}] - "
             "{11.026867565231802, 2.099569407959401}]] < 10.^-12");
    chk_true("Max[Abs[SingularValueList[{{{1.,2.},{3.,4.}}, {{1.,0.},{0.,0.}}}] - {0.4472135954999587}]] < 10.^-13");
    chk_msg("SingularValueList[{{{1,2},{3,4}}, {{2,0},{0,1}}}]", "SingularValueList::gfargs");
}

static void test_storage(void) {
    chk("a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Head[SingularValueList[NDArray[a]]] === NDArray");
    chk_true("Normal[SingularValueList[NDArray[a]]] === SingularValueList[a]");
    chk_true("SingularValueList[Developer`ToPackedArray[a]] === SingularValueList[a]");
    chk_true("SingularValueList[SparseArray[a]] === SingularValueList[a]");
    chk_if("SingularValueList[SparseArray[{{1,2},{3,4}}]]", "{Sqrt[15 + Sqrt[221]], Sqrt[15 - Sqrt[221]]}");
}

static void test_messages(void) {
    chk_msg("SingularValueList[{{1.,2.},{3.,4.}}, 3]", "SingularValueList::take");
    chk_msg("SingularValueList[{{1.,2.},{3.,4.}}, -3]", "SingularValueList::take");
    chk_msg("SingularValueList[{{1.,2.},{3.,4.}}, x]", "SingularValueList::seqs");
    chk_msg("SingularValueList[{1,2}]", "SingularValueList::matrix");
    chk_msg("SingularValueList[{}]", "SingularValueList::matrix");
    chk_msg("SingularValueList[{{1,2},{3,4}}, Method -> \"Foo\"]", "SingularValueList::emeth");
    chk_msg("SingularValueList[{{1.,2.,3.},{4.,5.,6.}}, Tolerance -> x]", "SingularValueList::tol");
    chk_true("SingularValueList[{{1.,2.},{3.,4.}}, Method -> \"Arnoldi\"] == SingularValueList[{{1.,2.},{3.,4.}}]");
    chk_if("Options[SingularValueList]",
           "{Cubics -> False, Method -> Automatic, Quartics -> False, Tolerance -> Automatic}");
    chk_if("Attributes[SingularValueList]", "{Protected}");
}

static void test_compile(void) {
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, SingularValueList[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, SingularValueList[m]]");
    chk("cf = Compile[{{m, _Real, 2}}, SingularValueList[m]]; "
        "cz = Compile[{{m, _Complex, 2}}, SingularValueList[m]]; "
        "a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; z = {{1.,I},{3.,4.}}; Null", "Null");
    chk_true("Normal[cf[a]] === SingularValueList[a]");
    chk_true("Normal[cf[a]] === Normal[SingularValueList[NDArray[a]]]");
    chk_true("Normal[cz[z]] === SingularValueList[z]");
    chk_true("cf[{{0.,0.},{0.,0.}}] == {0.0, 0.0}");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_exact();
    test_machine();
    test_mpfr();
    test_mpfr_bidiag();
    test_generalized();
    test_storage();
    test_messages();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d SingularValueList test(s) FAILED\n", failures);
        return 1;
    }
    printf("All SingularValueList tests passed.\n");
    return 0;
}

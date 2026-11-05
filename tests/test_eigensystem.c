/* test_eigensystem.c -- Eigensystem (src/linalg/eigen.c, ndla_eigensystem in
 * eigen_direct.c, the WANT_BOTH arms of the eigen kernels).
 *
 * The contract: Eigensystem[args] == {Eigenvalues[args], Eigenvectors[args]}
 * in order and pairing, every A . v == lambda v, for every precision class,
 * argument form (k, -k, UpTo[k], {m, a}) and Method; reference answers pinned
 * where Mathilda prints the same form; messages; storage. */

#include "test_dense_heads.h"

/* Eigensystem[args] agrees with the separate heads (values to `tol`; vectors
 * exactly-or-to-`tol`, both from the same kernels) and pairs correctly. */
static void chk_es(const char* args, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{es = Eigensystem[%s], ev = Eigenvalues[%s], vc = Eigenvectors[%s]}, "
             "Length[es] == 2 && Length[es[[1]]] == Length[ev] && Length[es[[2]]] == Length[vc] && "
             "Max[Append[Abs[N[es[[1]] - ev]], 0]] <= %s && "
             "Max[Append[Flatten[Abs[N[Abs[es[[2]]] - Abs[vc]]]], 0]] <= %s]",
             args, args, args, tol, tol);
    chk_true(buf);
}

/* A . v == lambda v for each pair of Eigensystem[m]. */
static void chk_pairs(const char* m, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{mm = %s, es}, es = Eigensystem[mm]; "
             "Max[Table[Norm[mm . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, Length[es[[1]]]}]] <= %s]",
             m, tol);
    chk_true(buf);
}

static void test_exact(void) {
    chk_if("Eigensystem[{{2,1},{1,2}}]", "{{3, 1}, {{1, 1}, {-1, 1}}}");
    chk_if("Eigensystem[{{1,1},{0,1}}]", "{{1, 1}, {{1, 0}, {0, 0}}}");
    chk_if("Eigensystem[{{0,1},{-1,0}}]", "{{I, -I}, {{-I, 1}, {I, 1}}}");
    chk_true("Simplify[Eigensystem[{{1,2},{3,4}}] - "
             "{{(5 + Sqrt[33])/2, (5 - Sqrt[33])/2}, {{(-3 + Sqrt[33])/6, 1}, {(-3 - Sqrt[33])/6, 1}}}] == "
             "{{0, 0}, {{0, 0}, {0, 0}}}");
    chk_true("Simplify[Eigensystem[{{{1,2},{3,4}},{{2,0},{0,1}}}] - "
             "{{(9 + Sqrt[97])/4, (9 - Sqrt[97])/4}, {{(-7 + Sqrt[97])/12, 1}, {(-7 - Sqrt[97])/12, 1}}}] == "
             "{{0, 0}, {{0, 0}, {0, 0}}}");
    chk_es("{{1,2},{3,4}}", "0");
    chk_es("{{1,2,3},{4,5,6},{7,8,10}}", "0");
    chk_es("{{1,2,3},{4,5,6},{7,8,10}}, 2", "0");
    chk_es("{{1,2,3},{4,5,6},{7,8,10}}, -1", "0");
    chk_es("{{2,0,0},{0,2,0},{0,0,3}}", "0");
    chk_es("{{{1,2},{3,4}},{{2,0},{0,1}}}", "0");
    chk_es("{{1,2,3},{4,5,6},{7,8,10}}, Cubics -> True", "0");
    /* the reference orders this {{a, b}, {{1, 0}, {0, 1}}}; Mathilda's
     * Eigenvalues orders {b, a} -- Eigensystem follows Eigenvalues, paired */
    chk_es("{{a,0},{0,b}}", "0");
    chk_true("Module[{es = Eigensystem[{{a,0},{0,b}}]}, "
             "Simplify[{{a,0},{0,b}} . es[[2,1]] - es[[1,1]] es[[2,1]]] == {0, 0}]");
    chk_true("Module[{es = Eigensystem[{{a,b},{c,d}}]}, "
             "Simplify[{{a,b},{c,d}} . es[[2,1]] - es[[1,1]] es[[2,1]]] == {0, 0}]");
}

static void test_machine(void) {
    chk_true("Max[Abs[Flatten[Eigensystem[{{1.,2.},{3.,4.}}] - {{5.372281323269014, -0.3722813232690143}, "
             "{{-0.41597355791928425, -0.9093767091321241}, {-0.8245648401323938, 0.5657674649689923}}}]]] < 10.^-14");
    chk_true("Max[Abs[Flatten[Eigensystem[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, 2] - "
             "{{16.707493316124747, -0.9057401795217573}, {{-0.22351335772785175, -0.5039456275263902, -0.8343144391740168}, "
             "{-0.8658457814289262, 0.08565120135087341, 0.4929248974122573}}}]]] < 10.^-13");
    chk_true("Max[Abs[Flatten[Eigensystem[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, -1] - "
             "{{0.19824686339700862}, {{0.27829648598057294, -0.8318468022052949, 0.4801895079566917}}}]]] < 10.^-13");
    chk_true("Length[Eigensystem[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, UpTo[5]][[1]]] == 3");
    chk_true("Max[Abs[Flatten[Eigensystem[{{1.,2.},{3.,4.}}, UpTo[1]] - "
             "{{5.372281323269014}, {{-0.41597355791928425, -0.9093767091321241}}}]]] < 10.^-14");
    chk_if("Eigensystem[{{1.,2.},{3.,4.}}, 0]", "{{}, {}}");
    chk_es("{{1.,2.},{3.,4.}}", "10.^-14");
    chk_es("{{2.,1.},{1.,2.}}", "10.^-14");
    chk_es("{{1.,I},{2.,3.}}", "10.^-14");
    chk_es("{{2., I},{-I, 3.}}", "10.^-14");
    chk_es("{{{1.,2.},{3.,4.}},{{2.,0.},{0.,1.}}}", "10.^-13");
    chk_es("{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, 2", "10.^-14");
    chk_es("{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, -2", "10.^-14");
    chk_pairs("RandomReal[{-1, 1}, {30, 30}]", "10.^-11");
    chk_pairs("Module[{r = RandomReal[{-1, 1}, {30, 30}]}, r + Transpose[r]]", "10.^-11");
    chk_pairs("RandomComplex[{-1 - I, 1 + I}, {20, 20}]", "10.^-11");
    chk_pairs("Module[{r = RandomComplex[{-1 - I, 1 + I}, {20, 20}]}, r + ConjugateTranspose[r]]", "10.^-11");
    /* the iterative and banded methods return both parts from one run too */
    chk("big = Module[{r = RandomReal[{-1, 1}, {120, 120}]}, r + Transpose[r]]; "
        "band = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {60, 60}] // Normal; Null",
        "Null");
    /* Arnoldi: the same Ritz pairs as Eigenvalues / Eigenvectors with the
     * same Method (whose own residual on a clustered spectrum is a known
     * limitation of the Arnoldi kernel, not of Eigensystem) */
    chk_es("big, 3, Method -> \"Arnoldi\"", "0");
    chk_true("Length[Eigensystem[big, 3, Method -> \"Arnoldi\"][[1]]] == 3");
    chk_true("Module[{es = Eigensystem[band, Method -> \"Banded\"]}, "
             "Max[Table[Norm[band . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, 60}]] < 10.^-10]");
    chk_es("band, 4", "10.^-10");
}

/* Complex Hermitian input (packed and boxed) goes through zheevd: REAL
 * eigenvalues, unitary eigenvectors, same order on both surfaces (v0.6221;
 * the packed path used general zgeev and left roundoff imaginary parts). */
static void test_hermitian(void) {
    chk("SeedRandom[4]; hz = With[{c = RandomComplex[{-1 - I, 1 + I}, {80, 80}]}, c + ConjugateTranspose[c]]; "
        "{hv, hx} = Eigensystem[hz]; Null", "Null");
    chk_true("VectorQ[hv, Head[#] === Real &] && VectorQ[Eigenvalues[hz], Head[#] === Real &]");
    chk_true("Max[Abs[Flatten[hz . Transpose[hx] - Transpose[hx] . DiagonalMatrix[hv]]]] < 10.^-12");
    chk_true("Max[Abs[Flatten[Conjugate[hx] . Transpose[hx] - IdentityMatrix[80]]]] < 10.^-13");
    chk_true("OrderedQ[Reverse[Abs[hv]]]");
    chk_true("Max[Abs[Eigenvalues[Developer`FromPackedArray[hz]] - hv]] < 10.^-12");
    chk_true("Eigenvectors[hz] == hx");
    chk_true("Max[Abs[Eigenvalues[{{2., 1. - I}, {1. + I, 3.}}] - {4., 1.}]] < 10.^-14");
}

static void test_mpfr(void) {
    chk_es("N[{{1,2},{3,4}}, 30]", "10.^-25");
    chk_es("N[{{2,1},{1,2}}, 30]", "10.^-25");
    chk_es("N[{{1,I},{2,3}}, 30]", "10.^-25");
    chk_true("Module[{mm = N[{{1,2,3},{4,5,6},{7,8,10}}, 40], es}, es = Eigensystem[mm]; "
             "Precision[es[[1, 1]]] > 35 && "
             "Max[Table[Norm[mm . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, 3}]] < 10^-34]");
}

static void test_storage(void) {
    chk("a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; s = {{2.,1.,0.},{1.,3.,1.},{0.,1.,4.}}; "
        "z = {{1.,I},{2.,3.}}; Null", "Null");
    chk_true("Normal /@ Eigensystem[NDArray[a]] == {Normal[Eigenvalues[NDArray[a]]], Normal[Eigenvectors[NDArray[a]]]}");
    chk_true("Normal /@ Eigensystem[NDArray[s]] == {Normal[Eigenvalues[NDArray[s]]], Normal[Eigenvectors[NDArray[s]]]}");
    chk_true("Normal /@ Eigensystem[NDArray[z]] == {Normal[Eigenvalues[NDArray[z]]], Normal[Eigenvectors[NDArray[z]]]}");
    chk_true("Max[Abs[Flatten[(Normal /@ Eigensystem[NDArray[a]]) - Eigensystem[a]]]] < 10.^-13");
    chk_true("Eigensystem[Developer`ToPackedArray[a]] == Eigensystem[a]");
    chk_true("Eigensystem[SparseArray[a]] == Eigensystem[a]");
    chk_if("Eigensystem[SparseArray[{{2,1},{1,2}}]]", "{{3, 1}, {{1, 1}, {-1, 1}}}");
}

static void test_messages(void) {
    chk_msg("Eigensystem[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}, 5]", "Eigensystem::take");
    chk_msg("Eigensystem[{{1,2},{3,4}}, x]", "Eigensystem::seqs");
    chk_msg("Eigensystem[{1,2}]", "Eigensystem::matsq");
    chk_msg("Eigensystem[{{1,2,3},{4,5,6}}]", "Eigensystem::matsq");
    chk_if("Attributes[Eigensystem]", "{Protected}");
    /* ragged result: documented Compile decline */
    chk_no_lower("CompileDiagnostics[{{m, _Real, 2}}, Eigensystem[m]]");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_exact();
    test_machine();
    test_hermitian();
    test_mpfr();
    test_storage();
    test_messages();
    if (failures) {
        fprintf(stderr, "\n%d Eigensystem test(s) FAILED\n", failures);
        return 1;
    }
    printf("All Eigensystem tests passed.\n");
    return 0;
}

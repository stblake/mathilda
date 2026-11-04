/* test_hessenberg.c -- HessenbergDecomposition (src/linalg/hessenberg.c).
 *
 * Machine answers pinned to the reference's digits (LAPACK reflector signs);
 * every result checked by the invariants p . h . ConjugateTranspose[p] == m,
 * p unitary, h upper Hessenberg.  MPFR against the machine answer, the
 * messages, storage, and the rank-3 Compile lowering. */

#include "test_dense_heads.h"

static void chk_hess(const char* m, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{mm = %s, p, h, n}, {p, h} = HessenbergDecomposition[mm]; n = Length[mm]; "
             "Max[Abs[p . h . ConjugateTranspose[p] - mm]] <= %s Max[1, Max[Abs[mm]]] && "
             "Max[Abs[p . ConjugateTranspose[p] - IdentityMatrix[n]]] <= %s && "
             "Max[Abs[LowerTriangularize[h, -2]]] == 0]", m, tol, tol);
    chk_true(buf);
}

static void test_machine(void) {
    chk("m3 = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Max[Abs[HessenbergDecomposition[m3] - "
             "{{{1., 0., 0.}, {0., -0.49613893835683376, -0.8682431421244592}, {0., -0.8682431421244592, 0.4961389383568338}}, "
             "{{1., -3.597007303087045, -0.24806946917841688}, {-8.06225774829855, 14.799999999999997, 2.3999999999999995}, "
             "{0., 0.399999999999997, 0.1999999999999986}}}]] < 10.^-13");
    chk_true("Max[Abs[HessenbergDecomposition[{{1.,2.,3.},{4.,I,6.},{7.,8.,10.}}] - "
             "{{{1., 0., 0.}, {0., -0.49613893835683376, 0.843182538700414 - 0.20709746564571568 I}, "
             "{0., -0.8682431421244592, -0.4818185935430938 + 0.11834140894040897 I}}, "
             "{{1., -3.597007303087045, 0.24090929677154685 - 0.05917070447020448 I}, "
             "{-8.06225774829855, 13.569230769230769 + 0.24615384615384606 I, -0.3417980822229101 - 0.3596219085964969 I}, "
             "{0., 1.805972995794362, -3.569230769230769 + 0.7538461538461534 I}}}]] < 10.^-13");
    chk_true("Max[Abs[HessenbergDecomposition[{{1.,2.,3.,4.},{5.,6.,7.,8.},{9.,1.,2.,3.},{4.,5.,6.,8.}}][[2]] - "
             "{{1., -4.798394540253481, 1.7038563377804858, 1.7527930329251211}, "
             "{-11.045361017187261, 11.344262295081963, -4.1799792179391595, -3.3619488147103667}, "
             "{0., -10.50781948633639, 3.96389649221709, 1.9548936442962153}, "
             "{0., 0., 0.3252502155308817, 0.6918412127009405}}]] < 10.^-12");
    chk_if("HessenbergDecomposition[{{1.,2.},{3.,4.}}]", "{{{1.0, 0.0}, {0.0, 1.0}}, {{1.0, 2.0}, {3.0, 4.0}}}");
    chk_if("HessenbergDecomposition[{{1, 2.},{3, 4}}]", "{{{1.0, 0.0}, {0.0, 1.0}}, {{1.0, 2.0}, {3.0, 4.0}}}");
    chk_if("HessenbergDecomposition[{{5.}}]", "{{{1.0}}, {{5.0}}}");
    chk_hess("RandomReal[{-1, 1}, {40, 40}]", "10.^-13");
    chk_hess("RandomComplex[{-1 - I, 1 + I}, {25, 25}]", "10.^-13");
    /* Pivoting: a third factor, the balancing matrix */
    chk_true("With[{r = HessenbergDecomposition[m3, Pivoting -> True]}, Length[r] == 3 && "
             "Max[Abs[r[[3]] - IdentityMatrix[3]]] == 0]");
    chk_true("Module[{mm = {{1., 100000., 0.}, {0.00001, 2., 3.}, {0., 4., 5.}}, p, h, d}, "
             "{p, h, d} = HessenbergDecomposition[mm, Pivoting -> True]; "
             "Max[Abs[mm . d - d . p . h . ConjugateTranspose[p]]] < 10.^-9]");
}

static void test_mpfr(void) {
    chk("mp = N[{{1,2,3},{4,5,6},{7,8,10}}, 30]; Null", "Null");
    chk_true("With[{r = HessenbergDecomposition[mp]}, r[[1, 1, 2]] === 0 && r[[2, 3, 1]] === 0 && "
             "Precision[r[[2, 2, 2]]] > 29 && Max[Abs[N[r] - HessenbergDecomposition[N[mp]]]] < 10.^-13]");
    chk_hess("N[{{1,2,3},{4,5,6},{7,8,10}}, 30]", "10^-27");
    chk_hess("N[{{1,2,3},{4,I,6},{7,8,10}}, 30]", "10^-27");
    chk_hess("SetPrecision[RandomComplex[{-1 - I, 1 + I}, {6, 6}], 40]", "10^-36");
    chk_true("Module[{mm = SetPrecision[RandomReal[{-1, 1}, {6, 6}], 30]}, "
             "Max[Abs[N[HessenbergDecomposition[mm]] - HessenbergDecomposition[N[mm]]]] < 10.^-12]");
    chk_true("Module[{mm = SetPrecision[RandomComplex[{-1 - I, 1 + I}, {5, 5}], 30]}, "
             "Max[Abs[N[HessenbergDecomposition[mm]] - HessenbergDecomposition[N[mm]]]] < 10.^-12]");
}

static void test_messages(void) {
    chk_msg("HessenbergDecomposition[{{1,2,3},{4,5,6},{7,8,10}}]", "HessenbergDecomposition::schurf");
    chk_msg("HessenbergDecomposition[{{a,b},{c,d}}]", "HessenbergDecomposition::schurn");
    chk_msg("HessenbergDecomposition[{{1.,2.,3.},{4.,5.,6.}}]", "HessenbergDecomposition::matsq");
    chk_msg("HessenbergDecomposition[{1.,2.}]", "HessenbergDecomposition::matsq");
    chk_msg("HessenbergDecomposition[{{1.,2.},{3.,4.}}, Method -> Automatic]", "HessenbergDecomposition::optx");
    chk_if("Options[HessenbergDecomposition]", "{Pivoting -> False, TargetStructure -> Automatic}");
    chk_if("Attributes[HessenbergDecomposition]", "{Protected}");
}

static void test_storage(void) {
    chk("m3 = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Head[First[HessenbergDecomposition[NDArray[m3]]]] === NDArray");
    chk_true("Normal /@ HessenbergDecomposition[NDArray[m3]] === HessenbergDecomposition[m3]");
    chk_true("HessenbergDecomposition[Developer`ToPackedArray[m3]] === HessenbergDecomposition[m3]");
    chk_true("HessenbergDecomposition[SparseArray[m3]] === HessenbergDecomposition[m3]");
}

static void test_compile(void) {
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, HessenbergDecomposition[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, HessenbergDecomposition[m]]");
    chk("ch = Compile[{{m, _Real, 2}}, HessenbergDecomposition[m]]; "
        "cz = Compile[{{m, _Complex, 2}}, HessenbergDecomposition[m]]; "
        "m3 = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; z3 = {{1.,2.,3.},{4.,I,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Dimensions[ch[m3]] == {2, 3, 3}");
    chk_true("Normal[ch[m3]] === HessenbergDecomposition[m3]");
    chk_true("Normal[cz[z3]] === HessenbergDecomposition[z3]");
    chk_true("Normal[ch[m3]] === (Normal /@ HessenbergDecomposition[NDArray[m3]])");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_machine();
    test_mpfr();
    test_messages();
    test_storage();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d HessenbergDecomposition test(s) FAILED\n", failures);
        return 1;
    }
    printf("All HessenbergDecomposition tests passed.\n");
    return 0;
}

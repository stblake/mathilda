/* test_orthogonalize.c -- Orthogonalize (src/linalg/orthogonalize.c).
 *
 * Exact and symbolic answers pinned to the reference where Mathilda's
 * arithmetic prints the same canonical form; otherwise checked by value and by
 * orthonormality.  Methods, Tolerance, user inner products, messages, every
 * precision class and storage, and the Compile lowering. */

#include "test_dense_heads.h"

/* rows of `o` orthonormal, or zero */
static void chk_orthonormal(const char* vs) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{o = Orthogonalize[%s], g}, g = o . ConjugateTranspose[o]; "
             "Max[Abs[g - DiagonalMatrix[Map[If[Norm[#] > 1/2, 1, 0] &, o]]]] < 10.^-12]", vs);
    chk_true(buf);
}

static void test_exact(void) {
    chk_if("Orthogonalize[{{1,2},{3,4}}]", "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}}");
    chk_if("Orthogonalize[{{1,0,0},{1,1,0},{1,1,1},{2,3,4}}]",
           "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}");
    chk_if("Orthogonalize[{{1, I},{2, 3}}]",
           "{{1/Sqrt[2], I/Sqrt[2]}, {(2 + 3*I)/Sqrt[26], (3 - 2*I)/Sqrt[26]}}");
    chk_if("Orthogonalize[{{0,0},{1,2}}]", "{{0, 0}, {1/Sqrt[5], 2/Sqrt[5]}}");
    chk_if("Orthogonalize[{{1,2},{3,4},{5,6}}]",
           "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}, {0, 0}}");
    chk_if("Orthogonalize[{}]", "{}");
    chk_if("Orthogonalize[{{1,2}}]", "{{1/Sqrt[5], 2/Sqrt[5]}}");
    /* the reference prints {Sqrt[13/14], -Sqrt[2/91], -3/Sqrt[182]} for the
     * third row: the same numbers */
    chk_true("With[{o = Orthogonalize[{{1,2,3},{2,4,6},{1,0,0}}]}, o[[2]] === {0, 0, 0} && "
             "Simplify[o[[3]] - {Sqrt[13/14], -Sqrt[2/91], -3/Sqrt[182]}] == {0, 0, 0}]");
    chk_if("Orthogonalize[{{1,2},{3,4}}, Dot]", "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}}");
    chk_if("Orthogonalize[{{1,2},{3,4}}, (#1 . #2) &]", "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}}");
    chk_if("Orthogonalize[{{1,2},{3,4}}, Method -> \"Householder\"]",
           "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}}");
    chk_true("With[{o = Orthogonalize[{1, x, x^2}, Integrate[#1 #2, {x, -1, 1}] &]}, "
             "Simplify[o - {1/Sqrt[2], Sqrt[3/2] x, (3 Sqrt[5/2] (-1/3 + x^2))/2}] == {0, 0, 0}]");
}

static void test_symbolic(void) {
    chk_if("Orthogonalize[{{a,b},{c,d}}][[1]]",
           "{a/Sqrt[Abs[a]^2 + Abs[b]^2], b/Sqrt[Abs[a]^2 + Abs[b]^2]}");
    chk_if("Orthogonalize[{{a,b},{c,d}}][[2, 1]]",
           "(c - (a*(c*Conjugate[a] + d*Conjugate[b]))/(Abs[a]^2 + Abs[b]^2))/"
           "Sqrt[Abs[c - (a*(c*Conjugate[a] + d*Conjugate[b]))/(Abs[a]^2 + Abs[b]^2)]^2 + "
           "Abs[d - (b*(c*Conjugate[a] + d*Conjugate[b]))/(Abs[a]^2 + Abs[b]^2)]^2]");
}

static void test_machine(void) {
    chk_if("Orthogonalize[{{1.,2.},{3.,4.}}]",
           "{{0.4472135954999579, 0.8944271909999159}, {0.8944271909999159, -0.4472135954999579}}");
    chk_true("Orthogonalize[{{1.,2.,3.},{2.,4.,6.},{1.,0.,0.}}][[2]] === {0.0, 0.0, 0.0}");
    chk_if("Orthogonalize[{{0.,0.},{1.,2.}}]", "{{0.0, 0.0}, {0.4472135954999579, 0.8944271909999159}}");
    chk_if("Orthogonalize[{{1.,2.},{3.,4.}}, Tolerance -> 0.5]",
           "{{0.4472135954999579, 0.8944271909999159}, {0.0, 0.0}}");
    chk_true("Last[Orthogonalize[{{1.,1.},{1.,1.0000000001}}]] != {0., 0.}");
    chk_true("Max[Abs[Orthogonalize[{{1., I},{2., 3.}}] - {{0.7071067811865476, 0.7071067811865476 I}, "
             "{0.3922322702763681 + 0.5883484054145521 I, 0.5883484054145521 - 0.3922322702763681 I}}]] < 10.^-15");
    chk_true("Max[Abs[Orthogonalize[{{1.,2.},{3.,4.}}, Method -> \"Householder\"] - "
             "{{-0.44721359549995787, -0.8944271909999157}, {-0.8944271909999157, 0.4472135954999581}}]] < 10.^-15");
    chk_true("Max[Abs[Orthogonalize[{{1.,2.},{3.,4.}}, Method -> \"GramSchmidt\"] - "
             "Orthogonalize[{{1,2},{3,4}}]]] < 10.^-14");
    chk_true("Max[Abs[Orthogonalize[{{1.,2.},{3.,4.}}, Method -> \"ModifiedGramSchmidt\"] - "
             "Orthogonalize[{{1,2},{3,4}}]]] < 10.^-14");
    chk_true("Orthogonalize[{{1.,2.},{3.,4.}}, Method -> \"Reorthogonalization\"] === "
             "Orthogonalize[{{1.,2.},{3.,4.}}]");
    chk_orthonormal("RandomReal[{-1, 1}, {8, 12}]");
    chk_orthonormal("RandomComplex[{-1 - I, 1 + I}, {6, 6}]");
    chk_orthonormal("Join[#, {#[[1]] + #[[2]]}] &[RandomReal[{-1, 1}, {4, 7}]]");
    chk_true("Last[Orthogonalize[Join[#, {#[[1]] + #[[2]]}] &[RandomReal[{-1, 1}, {4, 7}]]]] == ConstantArray[0., 7]");
    /* a Hilbert basis: MGS with reorthogonalisation stays orthonormal */
    chk_orthonormal("N[HilbertMatrix[8]]");
}

static void test_mpfr(void) {
    chk_true("With[{o = Orthogonalize[N[{{1,2},{3,4}}, 30]]}, Precision[o[[2,2]]] > 29 && "
             "Max[Abs[o - {{1, 2}, {2, -1}}/Sqrt[5]]] < 10^-28]");
    chk_true("With[{o = Orthogonalize[N[{{1,I},{2,3}}, 30]]}, "
             "Max[Abs[o - {{1/Sqrt[2], I/Sqrt[2]}, {(2 + 3 I)/Sqrt[26], (3 - 2 I)/Sqrt[26]}}]] < 10^-28]");
    chk_true("Orthogonalize[N[{{1,2,3},{2,4,6},{1,0,0}}, 30]][[2]] == {0, 0, 0}");
}

static void test_messages(void) {
    chk_msg("Orthogonalize[{{1,2},{3,4}}, Integrate[#1 #2, {x, 0, 1}] &]", "Orthogonalize::ornfa");
    chk_msg("Orthogonalize[{{1.,2.,3.},{2.,4.,6.},{1.,0.,0.}}, Method -> \"Householder\"]",
            "Orthogonalize::ornhsl");
    chk_msg("Orthogonalize[{{1,2},{3,4}}, Method -> \"Foo\"]", "Orthogonalize::orbop");
    chk_msg("Orthogonalize[{{1,2},{3,4}}, Dot, Method -> \"Householder\"]", "Orthogonalize::ornhip");
    chk_msg("Orthogonalize[{1,2}]", "Orthogonalize::orwdt");
    chk_msg("Orthogonalize[{{1,2},{3}}]", "Orthogonalize::ornt");
    chk_if("Options[Orthogonalize]", "{Method -> Automatic, Tolerance -> Automatic}");
    chk_if("Attributes[Orthogonalize]", "{Protected}");
}

static void test_storage(void) {
    chk("a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Head[Orthogonalize[NDArray[a]]] === NDArray");
    chk_true("Normal[Orthogonalize[NDArray[a]]] === Orthogonalize[a]");
    chk_true("Orthogonalize[Developer`ToPackedArray[a]] === Orthogonalize[a]");
    chk_true("Orthogonalize[SparseArray[a]] === Orthogonalize[a]");
    chk_if("Orthogonalize[SparseArray[{{1,2},{3,4}}]]", "{{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}}");
}

static void test_compile(void) {
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, Orthogonalize[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, Orthogonalize[m]]");
    chk("co = Compile[{{m, _Real, 2}}, Orthogonalize[m]]; cz = Compile[{{m, _Complex, 2}}, Orthogonalize[m]]; "
        "a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; z = {{1., I},{2., 3.}}; Null", "Null");
    chk_true("Normal[co[a]] === Orthogonalize[a]");
    chk_true("Normal[co[a]] === Normal[Orthogonalize[NDArray[a]]]");
    chk_true("Normal[cz[z]] === Orthogonalize[z]");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_exact();
    test_symbolic();
    test_machine();
    test_mpfr();
    test_messages();
    test_storage();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d Orthogonalize test(s) FAILED\n", failures);
        return 1;
    }
    printf("All Orthogonalize tests passed.\n");
    return 0;
}

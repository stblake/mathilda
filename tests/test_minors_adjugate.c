/* test_minors_adjugate.c -- Minors and Adjugate (src/linalg/minors.c).
 *
 * Exact and symbolic answers pinned to the reference's printed form
 * (lexicographic subset order); inexact ones by tolerance and by the
 * invariant Adjugate[m] . m == Det[m] I, including singular matrices. */

#include "test_dense_heads.h"

static void test_minors(void) {
    chk_if("Minors[{{1,2,3},{4,5,6},{7,8,10}}]", "{{-3, -6, -3}, {-6, -11, -4}, {-3, -2, 2}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}]",
           "{{-(b*d) + a*e, -(c*d) + a*f, -(c*e) + b*f}, {-(b*g) + a*h, -(c*g) + a*i, -(c*h) + b*i}, "
           "{-(e*g) + d*h, -(f*g) + d*i, -(f*h) + e*i}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 1]", "{{a, b, c}, {d, e, f}, {g, h, i}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 3]",
           "{{-(c*e*g) + b*f*g + c*d*h - a*f*h - b*d*i + a*e*i}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 0]", "{{1}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 4]", "{}");
    chk_if("Minors[{{1,2},{3,4}}, -1]", "{}");
    chk_if("Minors[{{a,b,c},{d,e,f}}]", "{{a, b, c}, {d, e, f}}");
    chk_if("Minors[{{a,b,c},{d,e,f}}, 2]", "{{-(b*d) + a*e, -(c*d) + a*f, -(c*e) + b*f}}");
    chk_if("Minors[{{a,b},{c,d},{e,f}}]", "{{a, b}, {c, d}, {e, f}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 2, ff]",
           "{{ff[{{a, b}, {d, e}}], ff[{{a, c}, {d, f}}], ff[{{b, c}, {e, f}}]}, "
           "{ff[{{a, b}, {g, h}}], ff[{{a, c}, {g, i}}], ff[{{b, c}, {h, i}}]}, "
           "{ff[{{d, e}, {g, h}}], ff[{{d, f}, {g, i}}], ff[{{e, f}, {h, i}}]}}");
    chk_if("Minors[{{a,b,c},{d,e,f},{g,h,i}}, 2, Tr]",
           "{{a + e, a + f, b + f}, {a + h, a + i, b + i}, {d + h, d + i, e + i}}");
    chk_if("Minors[{{5}}]", "{{1}}");
    chk_if("Minors[{{1/2,2},{3,4}}]", "{{1/2, 2}, {3, 4}}");
    chk_if("Minors[{{1,2},{3,4}}, x]", "Minors[{{1, 2}, {3, 4}}, x]");
    /* machine real / complex / MPFR against the exact minors */
    chk_true("Max[Abs[Minors[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}] - {{-3, -6, -3}, {-6, -11, -4}, {-3, -2, 2}}]] < 10.^-13");
    chk_true("Max[Abs[Minors[{{1.,2.,3.},{4.,I,6.},{7.,8.,10.}}] - Minors[{{1,2,3},{4,I,6},{7,8,10}}]]] < 10.^-13");
    chk_true("With[{mm = Minors[N[{{1,2,3},{4,5,6},{7,8,10}}, 30]]}, "
             "Max[Abs[mm - {{-3, -6, -3}, {-6, -11, -4}, {-3, -2, 2}}]] < 10^-27 && Precision[mm[[1,1]]] > 25]");
    chk_true("Module[{r = RandomReal[{-1, 1}, {5, 6}]}, "
             "Max[Abs[Minors[r, 3] - Minors[Rationalize[r, 0], 3]]] < 10.^-12]");
    chk_true("Dimensions[Minors[RandomReal[1, {5, 6}], 3]] == {10, 20}");
    chk_msg("Minors[{1,2}]", "Minors::matrix");
    chk_msg("Minors[{}]", "Minors::matrix");
    chk_if("Attributes[Minors]", "{Protected}");
}

static void test_adjugate(void) {
    chk_if("Adjugate[{{a,b},{c,d}}]", "{{d, -b}, {-c, a}}");
    chk_if("Adjugate[{{1,2,3},{4,5,6},{7,8,10}}]", "{{2, 4, -3}, {2, -11, 6}, {-3, 6, -3}}");
    chk_if("Adjugate[{{1,2,3},{4,5,6},{7,8,9}}]", "{{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}}");
    chk_if("Adjugate[{{5}}]", "{{1}}");
    chk_if("Adjugate[{{5.}}]", "{{1.0}}");
    chk_if("Adjugate[{{0,0},{0,0}}]", "{{0, 0}, {0, 0}}");
    chk_if("Adjugate[{{1/2,2},{3,4}}]", "{{4, -2}, {-3, 1/2}}");
    chk_true("Adjugate[HilbertMatrix[5]] == Det[HilbertMatrix[5]] Inverse[HilbertMatrix[5]]");
    chk_true("With[{m = {{1,2,3,4},{2,4,6,8},{1,0,1,0},{0,1,1,2}}}, Adjugate[m] . m == Det[m] IdentityMatrix[4]]");
    /* machine: singular and nonsingular, real and complex */
    chk_true("Max[Abs[Adjugate[{{1.,2.,3.},{4.,5.,6.},{7.,8.,9.}}] - {{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}}]] < 10.^-12");
    chk_true("Max[Abs[Adjugate[{{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}] - {{2, 4, -3}, {2, -11, 6}, {-3, 6, -3}}]] < 10.^-12");
    chk_true("Max[Abs[Adjugate[{{1.,2.},{3.,I}}] - {{I, -2}, {-3, 1}}]] < 10.^-14");
    chk_true("Module[{r = RandomComplex[{-1 - I, 1 + I}, {6, 6}]}, "
             "Max[Abs[Adjugate[r] . r - Det[r] IdentityMatrix[6]]] < 10.^-11]");
    chk_true("Module[{r = RandomReal[{-1, 1}, {7, 5}]}, With[{s = r . Transpose[r]}, "
             "Max[Abs[Adjugate[s] . s]] < 10.^-10]]");     /* rank 5 of 7: adj . s == 0 */
    chk_true("With[{a = Adjugate[N[{{1,2,3},{4,5,6},{7,8,10}}, 30]]}, "
             "Max[Abs[a - {{2, 4, -3}, {2, -11, 6}, {-3, 6, -3}}]] < 10^-27 && Precision[a[[1,1]]] > 25]");
    chk_msg("Adjugate[{{1,2,3},{4,5,6}}]", "Adjugate::matsq");
    chk_msg("Adjugate[{1,2}]", "Adjugate::matsq");
    chk_msg("Adjugate[{{1,2},{3,4}}, Modulus -> 5]", "Adjugate::intnm");
    chk_if("Attributes[Adjugate]", "{Protected}");
}

static void test_storage(void) {
    chk("a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; Null", "Null");
    chk_true("Head[Minors[NDArray[a]]] === NDArray && Head[Adjugate[NDArray[a]]] === NDArray");
    chk_true("Normal[Minors[NDArray[a]]] === Minors[a]");
    chk_true("Normal[Adjugate[NDArray[a]]] === Adjugate[a]");
    chk_true("Minors[Developer`ToPackedArray[a]] === Minors[a]");
    chk_true("Adjugate[Developer`ToPackedArray[a]] === Adjugate[a]");
    chk_if("Minors[SparseArray[{{1,2},{3,4}}]]", "{{1, 2}, {3, 4}}");
    chk_if("Adjugate[SparseArray[{{1,2},{3,4}}]]", "{{4, -2}, {-3, 1}}");
}

static void test_compile(void) {
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, Minors[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, Minors[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Real, 2}}, Adjugate[m]]");
    chk_lowers("CompileDiagnostics[{{m, _Complex, 2}}, Adjugate[m]]");
    chk("cm = Compile[{{m, _Real, 2}}, Minors[m]]; ca = Compile[{{m, _Real, 2}}, Adjugate[m]]; "
        "cz = Compile[{{m, _Complex, 2}}, Adjugate[m]]; "
        "a = {{1.,2.,3.},{4.,5.,6.},{7.,8.,10.}}; z = {{1.,2.},{3.,I}}; Null", "Null");
    chk_true("Normal[cm[a]] === Minors[a]");
    chk_true("Normal[ca[a]] === Adjugate[a]");
    chk_true("Normal[ca[a]] === Normal[Adjugate[NDArray[a]]]");
    chk_true("Normal[cz[z]] === Adjugate[z]");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_minors();
    test_adjugate();
    test_storage();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d Minors/Adjugate test(s) FAILED\n", failures);
        return 1;
    }
    printf("All Minors/Adjugate tests passed.\n");
    return 0;
}

/* Tests for auto-compilation of the numeric builtins (Plot/Plot3D now;
 * NIntegrate/FindRoot/Table as they are wired). These drive the real builtins
 * end-to-end and check that the compiled fast path agrees with the interpreter
 * to machine precision, that uncompilable bodies fall back cleanly, and that the
 * non-real-exclusion / fallback semantics are preserved. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* Every sampled Plot point lies on the true curve (compiled path == interpreter). */
void test_plot_parity(void) {
    assert_eval_eq(
        "With[{p = Cases[Plot[Sin[x] + Cos[2 x], {x, -3, 3}], Line[q_] :> q, Infinity][[1]]}, "
        "Max[Table[Abs[(Sin[p[[k,1]]] + Cos[2 p[[k,1]]]) - p[[k,2]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
    /* a rational + exp body */
    assert_eval_eq(
        "With[{p = Cases[Plot[Exp[-x^2] (1 + x)/(2 + x^2), {x, -2, 2}], Line[q_] :> q, Infinity][[1]]}, "
        "Max[Table[Abs[Exp[-p[[k,1]]^2] (1 + p[[k,1]])/(2 + p[[k,1]]^2) - p[[k,2]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
}

/* Body outside the compilable subset (Zeta has no machine kernel) still plots. */
void test_plot_fallback(void) {
    assert_eval_eq("Head[Plot[Zeta[x], {x, 2, 5}]]", "Graphics", 0);
    assert_eval_eq("Head[Plot[Sin[x], {x, 0, Pi}]]", "Graphics", 0);
}

/* Where the body goes complex (Sqrt of a negative), the point is excluded —
 * exactly as the interpreter path does. */
void test_plot_complex_excluded(void) {
    assert_eval_eq(
        "With[{p = Cases[Plot[Sqrt[x], {x, -1, 1}], Line[q_] :> q, Infinity][[1]]}, "
        "Min[Table[p[[k,1]], {k, 1, Length[p]}]] >= 0]",
        "True", 0);
}

/* Plot3D: every polygon vertex lies on the surface. */
void test_plot3d_parity(void) {
    assert_eval_eq(
        "With[{p = Flatten[Cases[Plot3D[Sin[x] Cos[y], {x, 0, 2}, {y, 0, 2}], Polygon[q_] :> q, Infinity], 1]}, "
        "Max[Table[Abs[Sin[p[[k,1]]] Cos[p[[k,2]]] - p[[k,3]]], {k, 1, Length[p]}]] < 10^-9]",
        "True", 0);
    assert_eval_eq("Head[Plot3D[Zeta[x + y], {x, 2, 3}, {y, 2, 3}]]", "Graphics3D", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_plot_parity);
    TEST(test_plot_fallback);
    TEST(test_plot_complex_excluded);
    TEST(test_plot3d_parity);

    printf("All auto-compile tests passed!\n");
    return 0;
}

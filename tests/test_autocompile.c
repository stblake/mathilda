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

/* Table: exact/symbolic iterators are UNTOUCHED; only inexact (machine-real)
 * iterators take the compiled path. */
void test_table_exact_untouched(void) {
    assert_eval_eq("Table[i^2, {i, 1, 5}]", "{1, 4, 9, 16, 25}", 0);          /* exact Integers */
    assert_eval_eq("Table[2^i, {i, 62, 64}]",
                   "{4611686018427387904, 9223372036854775808, 18446744073709551616}", 0); /* BigInt, no int64 overflow */
    assert_eval_eq("Table[Sin[i], {i, 1, 3}]", "{Sin[1], Sin[2], Sin[3]}", 0); /* symbolic */
    assert_eval_eq("Table[1/i, {i, 1, 4}]", "{1, 1/2, 1/3, 1/4}", 0);          /* exact Rationals */
}

/* Real iterator: compiled result matches the interpreter to machine precision. */
void test_table_real_parity(void) {
    assert_eval_eq(
        "Max[Table[Abs[Table[Sin[x] + x^2, {x, 0., 10., 0.1}][[k]] - (Sin[(k-1) 0.1] + ((k-1) 0.1)^2)], "
        "{k, 1, 101}]] < 10^-11",
        "True", 0);
    assert_eval_eq("Table[x^2, {x, 1., 4., 1.}]", "{1.0, 4.0, 9.0, 16.0}", 0);
}

/* Real iterator where the body goes complex → per-element interpreter fallback. */
void test_table_real_complex_fallback(void) {
    /* Sqrt[-1.] must still come back complex (I), not NaN. */
    assert_eval_eq("Chop[Table[Sqrt[x], {x, -1., 1., 1.}] - {I, 0., 1.}] == {0, 0, 0}", "True", 0);
}

/* NIntegrate: finite / half-line / whole-line agree with the interpreter; the
 * oscillatory sub-method (which re-bodies a copied context) stays correct; a
 * body that goes complex falls back per-sample; uncompilable bodies still work. */
void test_nintegrate_parity(void) {
    assert_eval_eq("Abs[NIntegrate[Sin[x]^2, {x, 0, Pi}] - Pi/2] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-x], {x, 0, Infinity}] - 1] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}] - Sqrt[Pi]] < 10^-7", "True", 0);
}

void test_nintegrate_oscillatory(void) {
    /* Regression: the amplitude/phase decomposition copies the sampler context
     * and swaps the body — the compiled program must not carry over. */
    assert_eval_eq("Abs[NIntegrate[Cos[100000 x], {x, 0, 1}] - Sin[100000]/100000] < 10^-9", "True", 0);
}

void test_nintegrate_complex_fallback(void) {
    /* Sqrt[x-1] is imaginary on [0,1], real on [1,2]; the compiled real program
     * returns NaN on the imaginary part and the interpreter supplies it. */
    assert_eval_eq("Chop[NIntegrate[Sqrt[x - 1], {x, 0, 2}] - (2/3 + 2 I/3)] == 0", "True", 0);
}

/* Multi-dimensional NIntegrate (cubature / Monte-Carlo, ni_mc_sample). */
void test_nintegrate_multidim(void) {
    assert_eval_eq("Abs[NIntegrate[x y, {x, 0, 1}, {y, 0, 1}] - 1/4] < 10^-8", "True", 0);
    assert_eval_eq("Abs[NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}] - 1/8] < 10^-7", "True", 0);
    assert_eval_eq("Abs[NIntegrate[Exp[-(x^2 + y^2)], {x, -2, 2}, {y, -2, 2}] "
                   "- (NIntegrate[Exp[-x^2], {x, -2, 2}])^2] < 10^-6", "True", 0);
    /* Monte-Carlo (low default accuracy) */
    assert_eval_eq("Abs[NIntegrate[x^2 + y^2, {x, 0, 1}, {y, 0, 1}, Method -> \"MonteCarlo\"] - 2/3] < 10^-2", "True", 0);
    /* complex-going integrand → per-point interpreter fallback (nonzero imaginary part) */
    assert_eval_eq("Im[Chop[NIntegrate[Sqrt[x + y - 1], {x, 0, 1}, {y, 0, 1}]]] > 0", "True", 0);
}

void test_nintegrate_uncompilable(void) {
    /* Zeta has no machine kernel → the integrand stays on the interpreter path;
     * cross-check the machine result against the arbitrary-precision one. */
    assert_eval_eq("Abs[NIntegrate[Zeta[x]/x^3, {x, 2, 4}] "
                   "- NIntegrate[Zeta[x]/x^3, {x, 2, 4}, WorkingPrecision -> 30]] < 10^-6", "True", 0);
}

/* FindRoot: the compiled machine-real path converges to the same root as the
 * interpreter (checked by residual and by parity against an uncompilable
 * perturbation). Complex/MPFR/inert-derivative paths keep working. */
void test_findroot_parity(void) {
    /* residual f(root) ~ 0 for Newton / secant / bracket */
    assert_eval_eq("With[{r = x /. FindRoot[Cos[x] - x, {x, 0.5}]}, Abs[Cos[r] - r] < 10^-8]", "True", 0);
    assert_eval_eq("With[{r = x /. FindRoot[x^3 - x - 2, {x, 1, 2}]}, Abs[r^3 - r - 2] < 10^-8]", "True", 0);
    assert_eval_eq("With[{r = x /. FindRoot[x^2 - 2, {x, 1.5, 0, 3}]}, Abs[r^2 - 2] < 10^-8]", "True", 0);
    /* compiled result equals the interpreter (uncompilable perturbation) result */
    assert_eval_eq("Abs[(x /. FindRoot[Cos[x] - x, {x, 0.5}]) "
                   "- (x /. FindRoot[Cos[x] - x + 10^-290 Zeta[x + 100], {x, 0.5}])] < 10^-11", "True", 0);
}

void test_findroot_complex_and_inert(void) {
    assert_eval_eq("Chop[(x /. FindRoot[x^2 + 1, {x, I}]) - I] == 0", "True", 0);       /* complex root */
    assert_eval_eq("With[{r = x /. FindRoot[Zeta[x] - 2, {x, 1.5}]}, Abs[Zeta[r] - 2] < 10^-8]", "True", 0); /* inert deriv → FD */
}

/* FindRoot systems: each component f_i and Jacobian entry is compiled as a
 * function of all variables; the linear-solve Newton is otherwise unchanged. */
void test_findroot_system(void) {
    assert_eval_eq("{x, y} /. FindRoot[{x + y == 3, x - y == 1}, {{x, 0}, {y, 0}}]", "{2.0, 1.0}", 0);
    assert_eval_eq("With[{s = {x, y} /. FindRoot[{x^2 + y^2 == 4, x - y == 0}, {{x, 1}, {y, 1}}]}, "
                   "Max[Abs[s - {Sqrt[2], Sqrt[2]}]] < 10^-8]", "True", 0);
    assert_eval_eq("With[{s = {x, y, z} /. FindRoot[{x + y + z == 6, x^2 - y == 2, z - x == 1}, "
                   "{{x, 1}, {y, 1}, {z, 1}}]}, "
                   "Max[Abs[{s[[1]]+s[[2]]+s[[3]]-6, s[[1]]^2-s[[2]]-2, s[[3]]-s[[1]]-1}]] < 10^-8]", "True", 0);
    /* one component uncompilable (LogGamma) but with an evaluable derivative →
     * that component and its Jacobian row use the interpreter, the rest compile. */
    assert_eval_eq("With[{s = {x, y} /. FindRoot[{LogGamma[x] - y == 0, x - 3 == 0}, {{x, 2.5}, {y, 0}}]}, "
                   "Abs[s[[2]] - LogGamma[3]] < 10^-7]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_plot_parity);
    TEST(test_plot_fallback);
    TEST(test_plot_complex_excluded);
    TEST(test_plot3d_parity);
    TEST(test_table_exact_untouched);
    TEST(test_table_real_parity);
    TEST(test_table_real_complex_fallback);
    TEST(test_nintegrate_parity);
    TEST(test_nintegrate_multidim);
    TEST(test_nintegrate_oscillatory);
    TEST(test_nintegrate_complex_fallback);
    TEST(test_nintegrate_uncompilable);
    TEST(test_findroot_parity);
    TEST(test_findroot_complex_and_inert);
    TEST(test_findroot_system);

    printf("All auto-compile tests passed!\n");
    return 0;
}

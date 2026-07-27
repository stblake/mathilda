/* Tests for the user-facing Compile[] / CompiledFunction object (M1b).
 *
 * These drive the whole pipeline (parse → evaluate → apply), unlike
 * test_compile.c which unit-tests the compile_expr engine directly.  Float
 * results are checked via a Mathilda equality/tolerance expression that
 * evaluates to True, so no fragile decimal-string comparisons are needed;
 * integer and symbolic results are compared by exact printed form. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* Integer-typed compiled functions return exact Integers. */
void test_cf_integer(void) {
    assert_eval_eq("Compile[{{n, _Integer}}, n^2 + 1][5]", "26", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, a b - b][6, 7]", "35", 0);
}

/* Real path, default (bare-symbol) type, and a machine-precision identity. */
void test_cf_real(void) {
    assert_eval_eq("Compile[{{x, _Real}}, x^2 + 1][3.0] == 10", "True", 0);
    assert_eval_eq("Compile[{x}, x + x][2.5] == 5", "True", 0);          /* default _Real */
    assert_eval_eq("Compile[{{x, _Real}}, x^2][3] == 9", "True", 0);      /* Integer arg → Real */
    assert_eval_eq("Abs[Compile[{{x, _Real}}, Sin[x]^2 + Cos[x]^2][1.234] - 1] < 10^-12", "True", 0);
}

/* Complex arguments and results. */
void test_cf_complex(void) {
    assert_eval_eq("Chop[Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I] - (-3 + 4 I)] == 0", "True", 0);
}

/* Symbolic argument → interpreter fallback (still produces the right value). */
void test_cf_symbolic_fallback(void) {
    assert_eval_eq("Compile[{{x, _Real}}, x^2 + 1][a]", "1 + a^2", 0);
}

/* Body outside the compilable subset (Zeta has no machine kernel) → the object
 * is still built and application falls back to the interpreter. */
void test_cf_uncompilable_fallback(void) {
    assert_eval_eq("Compile[{{x, _Real}}, Zeta[x]][2] == Zeta[2]", "True", 0);
}

/* Procedural / control-flow bodies (M2c) reachable through Compile. */
void test_cf_procedural(void) {
    assert_eval_eq("Abs[Compile[{{n, _Integer}}, Module[{s = 0.}, Do[s = s + i, {i, 1, n}]; s]][10] - 55] < 10^-9", "True", 0);
    assert_eval_eq("Compile[{{x, _Real}}, Nest[Function[u, u/2], x, 3]][8.0] == 1", "True", 0);
}

/* The object prints as CompiledFunction[...]; a wrong-arity application stays
 * unevaluated (the object is not consumed). */
void test_cf_object_and_arity(void) {
    assert_eval_startswith("Compile[{{x, _Real}}, x^2]", "CompiledFunction[{x}");
    assert_eval_startswith("Compile[{{x, _Real}}, x][1.0, 2.0]", "CompiledFunction[{x}");
}

/* Array argspec: `{v, _Real, r}` declares a rank-r machine array, the same
 * spelling the Wolfram Language uses.  Until this landed the whole array engine
 * — fusion, any-rank, the reductions — was reachable only from C. */
void test_cf_array_argspec(void) {
    /* A List argument is packed at the boundary and comes back as a List: the
     * result KIND must follow the argument kind, or the compiled path would
     * answer with a different head from the interpreter fallback on the same
     * input.  The two lines below are exactly that pair — the second has a
     * symbolic element, so it falls back — and they must agree in shape. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, v^2 + 2 v + 1][{1., 2., 3.}]",
                   "{4.0, 9.0, 16.0}", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][{1., 2.}]]", "List", 0);
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][{1., 2., x}]]", "List", 0);

    /* An NDArray argument stays an NDArray — borrowed in, new one out. */
    assert_eval_eq("Head[Compile[{{v, _Real, 1}}, 2 v][NDArray[{1., 2.}]]]", "NDArray", 0);

    /* Any rank. */
    assert_eval_eq("Compile[{{m, _Real, 2}}, m + 1][{{1., 2.}, {3., 4.}}]",
                   "{{2.0, 3.0}, {4.0, 5.0}}", 0);
    assert_eval_eq("Compile[{{t, _Real, 3}}, 2 t][{{{1., 2.}}, {{3., 4.}}}]",
                   "{{{2.0, 4.0}}, {{6.0, 8.0}}}", 0);

    /* Reduction to a scalar, and several array parameters. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, Total[v^2]][{1., 2., 3.}] == 14", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}, {u, _Real, 1}}, v u + v][{1., 2.}, {3., 4.}]",
                   "{4.0, 10.0}", 0);

    /* Complex elements. */
    assert_eval_eq("Chop[Total[Compile[{{v, _Complex, 1}}, v v][{1. + 2. I, 3.}]"
                   " - {-3 + 4 I, 9}]] == 0", "True", 0);

    /* Parity with the interpreter over a libm chain — same body, same data. */
    assert_eval_eq("Max[Abs[Compile[{{v, _Real, 1}}, Sin[v] Exp[-v] + Sqrt[v]][{0.5, 1.5, 2.5}]"
                   " - (Sin[#] Exp[-#] + Sqrt[#] & /@ {0.5, 1.5, 2.5})]] < 10^-12", "True", 0);

    /* Shapes that must NOT take the fast path: a ragged list cannot be packed,
     * a rank mismatch is not this function's signature, and a symbolic element
     * is not a machine number.  All three fall back and stay correct. */
    assert_eval_eq("Compile[{{v, _Real, 1}}, 2 v][{{1., 2.}, {3.}}]",
                   "{{2.0, 4.0}, {6.0}}", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, 2 v][{{1., 2.}, {3., 4.}}]",
                   "{{2.0, 4.0}, {6.0, 8.0}}", 0);

    /* Malformed argspecs must leave Compile[] unevaluated rather than guess. */
    assert_eval_eq("Head[Compile[{{v, _Real, 0.5}}, v]]", "Compile", 0);
    assert_eval_eq("Head[Compile[{{v, _Integer, 2}}, v]]", "Compile", 0);  /* no integer dtype */
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cf_integer);
    TEST(test_cf_real);
    TEST(test_cf_complex);
    TEST(test_cf_symbolic_fallback);
    TEST(test_cf_uncompilable_fallback);
    TEST(test_cf_procedural);
    TEST(test_cf_object_and_arity);
    TEST(test_cf_array_argspec);

    printf("All CompiledFunction tests passed!\n");
    return 0;
}

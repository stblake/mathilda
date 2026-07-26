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

    printf("All CompiledFunction tests passed!\n");
    return 0;
}

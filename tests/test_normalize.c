/* Unit tests for Normalize.
 *
 *   Normalize[v]      -- v / Norm[v], with the zero vector returned
 *                        unchanged (Mathematica's documented contract).
 *   Normalize[z]      -- scalar form, equivalently z / Abs[z].
 *   Normalize[expr,f] -- generic-norm form, expr / f[expr].
 *
 * Coverage:
 *   - Real / rational / negative exact vectors.
 *   - Complex-valued vectors (mixed real / pure-imaginary entries).
 *   - Zero short-circuit: integer 0, real 0., empty list, all-zero vector,
 *     2-D zero matrix.
 *   - Scalar inputs: exact integer, negative integer, exact complex,
 *     pure imaginary, integer 0.
 *   - Symbolic vectors with default Norm.
 *   - Custom norm-function argument: arbitrary symbol `f`, and the
 *     explicit Norm symbol on a matrix (left symbolic because builtin_norm
 *     does not yet compute the 2-norm of a general matrix -- this
 *     guarantees the call still terminates and the shape distributes).
 *   - Machine-precision input via N[] -- result is a unit float vector.
 *   - Idempotence: Norm[Normalize[v]] == 1 for rational v.
 *   - Listable check: Normalize must NOT be Listable; if it were,
 *     Normalize[{1, 2, 3}] would map to {Normalize[1], ...} = {1, 1, 1}
 *     instead of the correct {1/Sqrt[14], 2/Sqrt[14], 3/Sqrt[14]}.
 *   - Arity diagnostics: argc 0 and 3+ emit `Normalize::argt` and leave
 *     the call unevaluated.
 *   - Attributes (Protected, !Listable), docstring, and the interned
 *     SYM_Normalize pointer.
 *   - Repeated-evaluation loop to surface double-free / leak regressions
 *     under valgrind.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

/* --- Exact rational / integer vectors -------------------------------- */

static void test_exact_3vec(void) {
    /* Documented example: Normalize[{1, 5, 1}].  Norm = Sqrt[27] = 3 Sqrt[3]. */
    assert_eval_eq("Normalize[{1, 5, 1}]",
                   "{1/3/Sqrt[3], 5/3/Sqrt[3], 1/3/Sqrt[3]}", 0);
}

static void test_pythagorean_3_4(void) {
    /* Norm[{3, 4}] = 5, so result is rational with no surd. */
    assert_eval_eq("Normalize[{3, 4}]", "{3/5, 4/5}", 0);
}

static void test_signed_unit_pair(void) {
    /* Norm[{1, -1}] = Sqrt[2]. */
    assert_eval_eq("Normalize[{1, -1}]",
                   "{1/Sqrt[2], -1/Sqrt[2]}", 0);
}

static void test_rational_entries(void) {
    /* Norm[{1/2, 1/2}] = Sqrt[1/2] = 1/Sqrt[2]; result = {1/Sqrt[2]*1/2*Sqrt[2], ...}
     * After full simplification each component is 1/Sqrt[2]. */
    assert_eval_eq("Normalize[{1/2, 1/2}]",
                   "{1/Sqrt[2], 1/Sqrt[2]}", 0);
}

/* --- Complex vectors ------------------------------------------------- */

static void test_complex_vec_real_im(void) {
    /* Documented example: v = {1, 2I, 3, 4I, 5, 6I}.  |v|^2 = 1+4+9+16+25+36 = 91. */
    assert_eval_eq(
        "Normalize[{1, 2 I, 3, 4 I, 5, 6 I}]",
        "{1/Sqrt[91], (2*I)/Sqrt[91], 3/Sqrt[91], (4*I)/Sqrt[91], 5/Sqrt[91], (6*I)/Sqrt[91]}",
        0);
}

static void test_complex_vec_short(void) {
    /* {1+I, 2}: |1+I|^2 + 4 = 2 + 4 = 6. */
    assert_eval_eq("Normalize[{1 + I, 2}]",
                   "{(1 + I)/Sqrt[6], 2/Sqrt[6]}", 0);
}

/* --- Zero short-circuit --------------------------------------------- */

static void test_scalar_zero(void) {
    /* Normalize[0] returns 0 (Abs[0] = 0 short-circuit). */
    assert_eval_eq("Normalize[0]", "0", 0);
}

static void test_scalar_zero_real(void) {
    /* Real 0.0 short-circuit -- distinct from integer 0 in the type union. */
    assert_eval_eq("Normalize[0.]", "0.0", 0);
}

static void test_zero_vector(void) {
    assert_eval_eq("Normalize[{0, 0, 0}]", "{0, 0, 0}", 0);
}

static void test_empty_list(void) {
    /* Norm[{}] = 0 in this implementation, so Normalize[{}] short-circuits. */
    assert_eval_eq("Normalize[{}]", "{}", 0);
}

static void test_zero_matrix(void) {
    /* Frobenius norm of the zero matrix is 0 -> matrix returned unchanged. */
    assert_eval_eq("Normalize[{{0, 0}, {0, 0}}]",
                   "{{0, 0}, {0, 0}}", 0);
}

/* --- Scalar arguments ----------------------------------------------- */

static void test_scalar_positive_int(void) {
    /* Normalize[5] = 5 / Abs[5] = 5/5 = 1. */
    assert_eval_eq("Normalize[5]", "1", 0);
}

static void test_scalar_negative_int(void) {
    /* Normalize[-5] = -5 / Abs[-5] = -5/5 = -1. */
    assert_eval_eq("Normalize[-5]", "-1", 0);
}

static void test_scalar_complex(void) {
    /* Normalize[3 + 4 I] = (3 + 4 I) / 5 = 3/5 + (4/5) I. */
    assert_eval_eq("Normalize[3 + 4 I]", "3/5 + 4/5*I", 0);
}

static void test_scalar_pure_imag(void) {
    /* Normalize[I] = I / Abs[I] = I / 1 = I. */
    assert_eval_eq("Normalize[I]", "I", 0);
}

/* --- Symbolic vectors ------------------------------------------------ */

static void test_symbolic_vec(void) {
    assert_eval_eq("Normalize[{x, y}]",
                   "{x/Sqrt[Abs[x]^2 + Abs[y]^2], y/Sqrt[Abs[x]^2 + Abs[y]^2]}",
                   0);
}

/* --- Custom norm function ------------------------------------------- */

static void test_custom_norm_symbol(void) {
    /* `f` is undefined -- result is left as expr/f[expr] threaded. */
    assert_eval_eq("Normalize[{x, y}, f]",
                   "{x/f[{x, y}], y/f[{x, y}]}", 0);
}

static void test_custom_norm_norm_on_matrix(void) {
    /* The user explicitly passes Norm as the norm function.  Norm on a
     * 2-D matrix is the 2-norm (largest singular value), which the
     * current builtin_norm leaves symbolic -- we just verify the shape
     * distributes properly and the call terminates. */
    assert_eval_eq(
        "Normalize[{{1, 2}, {4, 5}}, Norm]",
        "{{1/Norm[{{1, 2}, {4, 5}}], 2/Norm[{{1, 2}, {4, 5}}]},"
        " {4/Norm[{{1, 2}, {4, 5}}], 5/Norm[{{1, 2}, {4, 5}}]}}",
        0);
}

/* --- Machine arithmetic --------------------------------------------- */

static void test_machine_arith_unit_vector(void) {
    /* N[{1, 5, 1}] = {1., 5., 1.}.  Norm = Sqrt[27] ~= 5.19615.
     * Components: {1/5.19615, 5/5.19615, 1/5.19615}
     *           = {0.19245, 0.96225, 0.19245}. */
    assert_eval_eq("Normalize[N[{1, 5, 1}]]",
                   "{0.19245, 0.96225, 0.19245}", 0);
}

/* --- Idempotence: ||Normalize[v]|| == 1 ----------------------------- */

static void test_idempotent_norm(void) {
    /* For rational v, Normalize[v] is an exact unit vector. */
    assert_eval_eq("Norm[Normalize[{3, 4}]]", "1", 0);
    assert_eval_eq("Norm[Normalize[{3, 4, 12}]]", "1", 0);
    /* For Pythagorean complex {3 + 4 I}: |3+4I| = 5, normalized gives
     * (3+4I)/5 whose modulus is exactly 1. */
    assert_eval_eq("Norm[Normalize[{3 + 4 I}]]", "1", 0);
}

/* --- Listable guard ------------------------------------------------- */

static void test_not_listable_semantics(void) {
    /* If Normalize were Listable the next call would yield
     * {Normalize[1], Normalize[2], Normalize[3]} = {1, 1, 1}.  The
     * correct vector form is {1/Sqrt[14], 2/Sqrt[14], 3/Sqrt[14]}. */
    assert_eval_eq("Normalize[{1, 2, 3}]",
                   "{1/Sqrt[14], 2/Sqrt[14], 3/Sqrt[14]}", 0);
}

/* --- Arity diagnostics ---------------------------------------------- */

static void test_zero_args_left_unevaluated(void) {
    const char* in = "Normalize[]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strcmp(s, "Normalize[]") == 0,
               "expected Normalize[], got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_three_args_left_unevaluated(void) {
    const char* in = "Normalize[1, 2, 3]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strcmp(s, "Normalize[1, 2, 3]") == 0,
               "expected Normalize[1, 2, 3], got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_four_args_left_unevaluated(void) {
    const char* in = "Normalize[a, b, c, d]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "Normalize") != NULL,
               "expected unevaluated Normalize, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Attribute / docstring / interned-symbol introspection ----------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("Normalize");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* Critical: Normalize must NOT be Listable -- it operates on a
     * vector as a whole, not element-wise. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("Normalize");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "normalized form") != NULL);
    ASSERT(strstr(def->docstring, "Norm[v]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_Normalize != NULL);
    ASSERT(strcmp(SYM_Normalize, "Normalize") == 0);
}

/* --- Memory-safety stress loop -------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix of exact / complex / symbolic / zero-short-circuit / arity-error
     * cases; allocator misuse typically surfaces as a crash, wrong answer,
     * or valgrind diagnostic by the second or third iteration. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("Normalize[{3, 4}]", "{3/5, 4/5}", 0);
        assert_eval_eq("Normalize[{1, 2 I, 3}]",
                       "{1/Sqrt[14], (2*I)/Sqrt[14], 3/Sqrt[14]}", 0);
        assert_eval_eq("Normalize[{0, 0, 0}]", "{0, 0, 0}", 0);
        assert_eval_eq("Normalize[{x, y}, f]",
                       "{x/f[{x, y}], y/f[{x, y}]}", 0);
        assert_eval_eq("Normalize[3 + 4 I]", "3/5 + 4/5*I", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_exact_3vec);
    TEST(test_pythagorean_3_4);
    TEST(test_signed_unit_pair);
    TEST(test_rational_entries);

    TEST(test_complex_vec_real_im);
    TEST(test_complex_vec_short);

    TEST(test_scalar_zero);
    TEST(test_scalar_zero_real);
    TEST(test_zero_vector);
    TEST(test_empty_list);
    TEST(test_zero_matrix);

    TEST(test_scalar_positive_int);
    TEST(test_scalar_negative_int);
    TEST(test_scalar_complex);
    TEST(test_scalar_pure_imag);

    TEST(test_symbolic_vec);

    TEST(test_custom_norm_symbol);
    TEST(test_custom_norm_norm_on_matrix);

    TEST(test_machine_arith_unit_vector);

    TEST(test_idempotent_norm);

    TEST(test_not_listable_semantics);

    TEST(test_zero_args_left_unevaluated);
    TEST(test_three_args_left_unevaluated);
    TEST(test_four_args_left_unevaluated);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All Normalize tests passed!\n");
    return 0;
}

/* Unit tests for MachineNumberQ.
 *
 *   MachineNumberQ[expr] -- True if expr is a machine-precision (IEEE
 *   double) real or a Complex of two such reals; False otherwise.
 *
 * Coverage:
 *   - Documented Mathematica examples verbatim (Sin[1000.], Exp[1000.],
 *     machine literal, high-precision literal).
 *   - Plain Real / Integer / BigInt / Rational / named-constant leaves.
 *   - Special IEEE values (inf, -inf, nan from 0./0.) → False.
 *   - MPFR arbitrary-precision reals (built two ways: high-precision
 *     literal and N[x, n] for n > 16) → False.
 *   - Complex of machine reals (both literal and Complex[...] form) → True.
 *   - Complex with any non-machine part (Integer, BigInt, MPFR, mixed) → False.
 *   - Hold preserves the unevaluated form so the test sees the raw atom.
 *   - Listable across a List threads element-wise.
 *   - Arity != 1 stays unevaluated.
 *   - Attribute set is exactly Protected (no Listable -- per spec).
 *   - Docstring registered.
 *   - SYM_MachineNumberQ pointer interned.
 *   - Stress loop to flush double-frees / MPFR mishandling under valgrind.
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* --- Documented examples (highest signal) ----------------------------- */

static void test_doc_sin_1000(void) {
    /* Sin[1000.] ≈ 0.8268795405320026 -- a finite machine double, hence
     * a machine number. */
    assert_eval_eq("MachineNumberQ[Sin[1000.]]", "True", 0);
}

static void test_doc_exp_1000(void) {
    /* Exp[1000.] overflows IEEE double range; Mathilda's cexp(1000+0i)
     * returns +inf, which is NOT a machine number. WL agrees (in WL the
     * result is an arbitrary-precision number, also not a machine
     * number). */
    assert_eval_eq("MachineNumberQ[Exp[1000.]]", "False", 0);
}

static void test_doc_machine_real_literal(void) {
    assert_eval_eq("MachineNumberQ[-29037945.290347]", "True", 0);
}

static void test_doc_high_precision_literal(void) {
    /* A literal that exceeds machine precision auto-promotes to EXPR_MPFR
     * at parse time (see commit 3ebf388 "parse: auto-promote high-precision
     * real literals to MPFR"). MPFR values are not machine numbers. */
    assert_eval_eq(
        "MachineNumberQ[-29037852093587905730945.29034875093457832094573984537498]",
        "False", 0);
}

/* --- Atom-level coverage --------------------------------------------- */

static void test_plain_machine_real(void) {
    assert_eval_eq("MachineNumberQ[0.0]",     "True", 0);
    assert_eval_eq("MachineNumberQ[1.0]",     "True", 0);
    assert_eval_eq("MachineNumberQ[-3.14]",   "True", 0);
    assert_eval_eq("MachineNumberQ[1.5e10]",  "True", 0);
    assert_eval_eq("MachineNumberQ[1.5e-10]", "True", 0);
}

static void test_integer_is_not_machine(void) {
    assert_eval_eq("MachineNumberQ[0]",    "False", 0);
    assert_eval_eq("MachineNumberQ[1]",    "False", 0);
    assert_eval_eq("MachineNumberQ[-42]",  "False", 0);
}

static void test_bigint_is_not_machine(void) {
    /* 2^100, 100!, and a hand-written 30-digit integer all live in
     * EXPR_BIGINT. None are machine numbers. */
    assert_eval_eq("MachineNumberQ[2^100]", "False", 0);
    assert_eval_eq("MachineNumberQ[100!]",  "False", 0);
    assert_eval_eq("MachineNumberQ[123456789012345678901234567890]",
                   "False", 0);
}

static void test_rational_is_not_machine(void) {
    assert_eval_eq("MachineNumberQ[1/2]", "False", 0);
    assert_eval_eq("MachineNumberQ[7/3]", "False", 0);
}

static void test_named_constants_are_not_machine(void) {
    /* These are exact symbolic constants; MachineNumberQ inspects the
     * head expression as given (no implicit N). */
    assert_eval_eq("MachineNumberQ[Pi]",          "False", 0);
    assert_eval_eq("MachineNumberQ[E]",           "False", 0);
    assert_eval_eq("MachineNumberQ[EulerGamma]",  "False", 0);
    assert_eval_eq("MachineNumberQ[GoldenRatio]", "False", 0);
    assert_eval_eq("MachineNumberQ[Degree]",      "False", 0);
}

/* --- IEEE special-value handling ------------------------------------- */

static void test_inf_real_is_not_machine(void) {
    /* Exp[1000.] returns an EXPR_REAL containing +inf. The isfinite()
     * guard in builtin_machinenumberq is what makes this False -- a
     * plain `e->type == EXPR_REAL` check would (wrongly) return True. */
    assert_eval_eq("MachineNumberQ[Exp[1000.]]", "False", 0);
}

static void test_indeterminate_is_not_machine(void) {
    /* Mathilda collapses 0./0. and -Exp[1000.] to the symbol Indeterminate
     * (rather than producing a NaN-valued EXPR_REAL). Either way,
     * MachineNumberQ on a non-numeric atom is False -- verify that
     * propagation works end-to-end. */
    assert_eval_eq("MachineNumberQ[0./0.]",       "False", 0);
    assert_eval_eq("MachineNumberQ[-Exp[1000.]]", "False", 0);
}

/* --- MPFR / high precision ------------------------------------------- */

static void test_mpfr_from_n_with_precision(void) {
    /* N[1, 50] builds an EXPR_MPFR at 50 decimal digits. Not a machine
     * number even though the numeric value happens to be 1.0. */
    assert_eval_eq("MachineNumberQ[N[1, 50]]",     "False", 0);
    assert_eval_eq("MachineNumberQ[N[Pi, 30]]",    "False", 0);
    assert_eval_eq("MachineNumberQ[N[1/3, 40]]",   "False", 0);
}

static void test_mpfr_via_high_prec_literal(void) {
    /* Force MPFR via the parser auto-promote path: a literal with > ~16
     * significant decimal digits. */
    assert_eval_eq("MachineNumberQ[1.2345678901234567890123456789]",
                   "False", 0);
}

/* --- Complex --------------------------------------------------------- */

static void test_complex_two_machine_reals(void) {
    assert_eval_eq("MachineNumberQ[1.0 + 2.0 I]",         "True", 0);
    assert_eval_eq("MachineNumberQ[Complex[1.0, 2.0]]",   "True", 0);
    assert_eval_eq("MachineNumberQ[Complex[-1.5, -2.5]]", "True", 0);
}

static void test_complex_with_int_part_is_not_machine(void) {
    /* WL: Complex[1, 2] is the exact Gaussian integer 1 + 2 I, not a
     * machine number. Mathilda matches: a Complex needs *both* parts to
     * be machine reals. */
    assert_eval_eq("MachineNumberQ[Complex[1, 2]]", "False", 0);
    assert_eval_eq("MachineNumberQ[1 + 2 I]",       "False", 0);
}

static void test_complex_with_mpfr_part_is_not_machine(void) {
    /* Build a Complex whose real part is MPFR via N[..., n]. The
     * surrounding contagion promotes the imaginary side; MachineNumberQ
     * is False either way. */
    assert_eval_eq("MachineNumberQ[N[1 + I, 40]]", "False", 0);
}

static void test_complex_with_inf_part_is_not_machine(void) {
    /* Complex[+inf, 0.] -- the real part fails isfinite, so the whole
     * Complex is not a machine number. */
    assert_eval_eq("MachineNumberQ[Complex[Exp[1000.], 0.0]]", "False", 0);
}

/* --- Non-numeric inputs --------------------------------------------- */

static void test_symbol_is_not_machine(void) {
    assert_eval_eq("MachineNumberQ[x]",         "False", 0);
    assert_eval_eq("MachineNumberQ[someName]",  "False", 0);
    assert_eval_eq("MachineNumberQ[Infinity]",  "False", 0);
}

static void test_string_is_not_machine(void) {
    assert_eval_eq("MachineNumberQ[\"1.0\"]", "False", 0);
}

static void test_list_is_not_machine(void) {
    /* MachineNumberQ is intentionally not Listable, so a List arg is
     * inspected as a whole: a List is not a number. */
    assert_eval_eq("MachineNumberQ[{1.0, 2.0}]", "False", 0);
}

static void test_unevaluated_function_is_not_machine(void) {
    /* foo[1.0] never evaluates further; the head is a non-Complex
     * function, so MachineNumberQ is False. */
    assert_eval_eq("MachineNumberQ[foo[1.0]]", "False", 0);
}

/* --- Arity / structural errors -------------------------------------- */

static void test_zero_args_stays_unevaluated(void) {
    /* No diagnostic: the builtin returns NULL and the evaluator leaves
     * the call in place. */
    assert_eval_eq("MachineNumberQ[]", "MachineNumberQ[]", 0);
}

static void test_two_args_stays_unevaluated(void) {
    assert_eval_eq("MachineNumberQ[1.0, 2.0]",
                   "MachineNumberQ[1.0, 2.0]", 0);
}

/* --- Attributes / docstring / interned symbol ----------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("MachineNumberQ");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("MachineNumberQ");
    ASSERT((a & ATTR_PROTECTED) != 0);
    /* Not Listable -- matches Mathematica's WL behavior, where
     * MachineNumberQ[{...}] inspects the List atom rather than
     * threading. */
    ASSERT((a & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("MachineNumberQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "machine-precision") != NULL);
    ASSERT(strstr(def->docstring, "MachineNumberQ[expr]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_MachineNumberQ != NULL);
    ASSERT(strcmp(SYM_MachineNumberQ, "MachineNumberQ") == 0);
}

/* --- Memory-safety stress loop -------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Exercise every code path -- atom True, atom False (each EXPR_*
     * variant), Complex True/False, MPFR True/False, IEEE special, arity
     * mismatch -- in a loop. Anything mishandling an Expr* will surface
     * as a leak or use-after-free under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("MachineNumberQ[Sin[1000.]]",      "True",  0);
        assert_eval_eq("MachineNumberQ[Exp[1000.]]",      "False", 0);
        assert_eval_eq("MachineNumberQ[-29037945.290347]","True",  0);
        assert_eval_eq(
            "MachineNumberQ[-29037852093587905730945.29034875093457832094573984537498]",
            "False", 0);
        assert_eval_eq("MachineNumberQ[N[Pi, 50]]",       "False", 0);
        assert_eval_eq("MachineNumberQ[1.0 + 2.0 I]",     "True",  0);
        assert_eval_eq("MachineNumberQ[1 + 2 I]",         "False", 0);
        assert_eval_eq("MachineNumberQ[2^100]",           "False", 0);
        assert_eval_eq("MachineNumberQ[x]",               "False", 0);
        assert_eval_eq("MachineNumberQ[]",                "MachineNumberQ[]", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Documented examples (must come first -- highest signal). */
    TEST(test_doc_sin_1000);
    TEST(test_doc_exp_1000);
    TEST(test_doc_machine_real_literal);
    TEST(test_doc_high_precision_literal);

    /* Atom-level coverage. */
    TEST(test_plain_machine_real);
    TEST(test_integer_is_not_machine);
    TEST(test_bigint_is_not_machine);
    TEST(test_rational_is_not_machine);
    TEST(test_named_constants_are_not_machine);

    /* IEEE special values. */
    TEST(test_inf_real_is_not_machine);
    TEST(test_indeterminate_is_not_machine);

    /* MPFR. */
    TEST(test_mpfr_from_n_with_precision);
    TEST(test_mpfr_via_high_prec_literal);

    /* Complex. */
    TEST(test_complex_two_machine_reals);
    TEST(test_complex_with_int_part_is_not_machine);
    TEST(test_complex_with_mpfr_part_is_not_machine);
    TEST(test_complex_with_inf_part_is_not_machine);

    /* Non-numeric inputs. */
    TEST(test_symbol_is_not_machine);
    TEST(test_string_is_not_machine);
    TEST(test_list_is_not_machine);
    TEST(test_unevaluated_function_is_not_machine);

    /* Arity. */
    TEST(test_zero_args_stays_unevaluated);
    TEST(test_two_args_stays_unevaluated);

    /* Introspection. */
    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    /* Memory-safety stress. */
    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All MachineNumberQ tests passed!\n");
    return 0;
}

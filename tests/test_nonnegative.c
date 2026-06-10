/* Tests for NonNegative[x].
 *
 * NonNegative[x] gives True when x is a real number that is positive or zero,
 * False when x is a manifestly negative real or a non-real complex value, and
 * is left unevaluated for non-numeric arguments. Covers exact
 * integers/rationals/bigints, machine and arbitrary-precision (MPFR) reals,
 * symbolic constants, numeric functions, complex values, Listable threading,
 * attributes, and arity diagnostics. Reference values cross-checked against the
 * running binary and the task specification. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- exact numbers --------------------------------------------------- */

void test_nonnegative_exact_integers() {
    assert_eval_eq("NonNegative[1]", "True", 0);
    assert_eval_eq("NonNegative[5]", "True", 0);
    assert_eval_eq("NonNegative[-5]", "False", 0);
    assert_eval_eq("NonNegative[0]", "True", 0);
    /* BigInt path (well beyond int64). */
    assert_eval_eq("NonNegative[10^40]", "True", 0);
    assert_eval_eq("NonNegative[-10^40]", "False", 0);
}

void test_nonnegative_rationals() {
    assert_eval_eq("NonNegative[3/4]", "True", 0);
    assert_eval_eq("NonNegative[-3/4]", "False", 0);
    assert_eval_eq("NonNegative[-1/1000000000000]", "False", 0);
    assert_eval_eq("NonNegative[7/3]", "True", 0);
}

/* ---- machine reals --------------------------------------------------- */

void test_nonnegative_machine_reals() {
    assert_eval_eq("NonNegative[1.6]", "True", 0);
    assert_eval_eq("NonNegative[-2.5]", "False", 0);
    assert_eval_eq("NonNegative[0.0]", "True", 0);
    assert_eval_eq("NonNegative[-1.0*10^-12]", "False", 0);
}

/* ---- arbitrary precision (MPFR) -------------------------------------- */

void test_nonnegative_arbitrary_precision() {
    assert_eval_eq("NonNegative[N[Pi, 50]]", "True", 0);
    assert_eval_eq("NonNegative[-N[Pi, 50]]", "False", 0);
    assert_eval_eq("NonNegative[-N[E, 100]]", "False", 0);
    /* Difference that is a tiny but genuinely negative high-precision real. */
    assert_eval_eq("NonNegative[3 - N[Pi, 50]]", "False", 0);
}

/* ---- symbolic constants & numeric functions -------------------------- */

void test_nonnegative_constants_and_functions() {
    assert_eval_eq("NonNegative[Pi]", "True", 0);
    assert_eval_eq("NonNegative[E]", "True", 0);
    assert_eval_eq("NonNegative[Sqrt[2]]", "True", 0);
    assert_eval_eq("NonNegative[E^2]", "True", 0);
    assert_eval_eq("NonNegative[Log[2]]", "True", 0);
    assert_eval_eq("NonNegative[Log[1/2]]", "False", 0);   /* log(1/2) ~ -0.693 */
    /* Spec example: Sin[10^5] ~ 0.0357 > 0, so NonNegative is True. */
    assert_eval_eq("NonNegative[Sin[10^5]]", "True", 0);
    assert_eval_eq("NonNegative[Sin[4]]", "False", 0);     /* sin(4) ~ -0.7568 */
    /* Exact numeric combinations. */
    assert_eval_eq("NonNegative[1 - Pi]", "False", 0);
    assert_eval_eq("NonNegative[Pi - 3]", "True", 0);
}

/* ---- complex values are not non-negative ----------------------------- */

void test_nonnegative_complex() {
    assert_eval_eq("NonNegative[I]", "False", 0);
    assert_eval_eq("NonNegative[1 + I]", "False", 0);
    assert_eval_eq("NonNegative[2 - 3 I]", "False", 0);
    assert_eval_eq("NonNegative[Sqrt[-2]]", "False", 0);
    /* A Complex literal whose imaginary part collapses to zero is real. */
    assert_eval_eq("NonNegative[-2 + 0 I]", "False", 0);
    assert_eval_eq("NonNegative[2 + 0 I]", "True", 0);
}

/* ---- non-numeric arguments stay unevaluated -------------------------- */

void test_nonnegative_symbolic_unevaluated() {
    assert_eval_eq("NonNegative[x]", "NonNegative[x]", 0);
    assert_eval_eq("NonNegative[Sin[y]]", "NonNegative[Sin[y]]", 0);
    assert_eval_eq("NonNegative[a + b]", "NonNegative[a + b]", 0);
    /* FullForm confirms the head/arg survive intact. */
    assert_eval_eq("NonNegative[x]", "NonNegative[x]", 1);
}

/* ---- Listable threading (spec example) ------------------------------- */

void test_nonnegative_listable() {
    assert_eval_eq("NonNegative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]",
                   "{True, True, True, True, False, False, True}", 0);
    assert_eval_eq("NonNegative[{x, Sin[y]}]",
                   "{NonNegative[x], NonNegative[Sin[y]]}", 0);
    assert_eval_eq("NonNegative[{}]", "{}", 0);
    /* Nested lists thread element by element. */
    assert_eval_eq("NonNegative[{{1, -1}, {0, -2}}]",
                   "{{True, False}, {True, False}}", 0);
}

/* ---- attributes ------------------------------------------------------ */

void test_nonnegative_attributes() {
    SymbolDef* def = symtab_get_def("NonNegative");
    ASSERT(def != NULL);
    ASSERT_MSG(def->attributes & ATTR_LISTABLE,
               "NonNegative must be Listable");
    ASSERT_MSG(def->attributes & ATTR_PROTECTED,
               "NonNegative must be Protected");
    /* Attributes[NonNegative] should report exactly these two. */
    assert_eval_eq("Attributes[NonNegative]", "{Listable, Protected}", 0);
}

/* ---- docstring ------------------------------------------------------- */

void test_nonnegative_docstring() {
    const char* doc = symtab_get_docstring("NonNegative");
    ASSERT_MSG(doc != NULL && doc[0] != '\0',
               "NonNegative must have a non-empty docstring");
}

/* ---- arity: argx diagnostic + unevaluated result --------------------- */

void test_nonnegative_arity() {
    /* Wrong argument counts emit NonNegative::argx on stderr and leave the
     * expression unevaluated. We only check the returned (unevaluated) form
     * here; the diagnostic itself goes to stderr. */
    assert_eval_eq("NonNegative[]", "NonNegative[]", 0);
    assert_eval_eq("NonNegative[1, 2, 3]", "NonNegative[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_nonnegative_exact_integers);
    TEST(test_nonnegative_rationals);
    TEST(test_nonnegative_machine_reals);
    TEST(test_nonnegative_arbitrary_precision);
    TEST(test_nonnegative_constants_and_functions);
    TEST(test_nonnegative_complex);
    TEST(test_nonnegative_symbolic_unevaluated);
    TEST(test_nonnegative_listable);
    TEST(test_nonnegative_attributes);
    TEST(test_nonnegative_docstring);
    TEST(test_nonnegative_arity);

    symtab_clear();
    printf("All NonNegative tests passed.\n");
    return 0;
}

/* Tests for NonPositive[x].
 *
 * NonPositive[x] gives True when x is a real number that is negative or zero,
 * False when x is a manifestly positive real or a non-real complex value, and
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

void test_nonpositive_exact_integers() {
    assert_eval_eq("NonPositive[1]", "False", 0);
    assert_eval_eq("NonPositive[5]", "False", 0);
    assert_eval_eq("NonPositive[-5]", "True", 0);
    assert_eval_eq("NonPositive[0]", "True", 0);
    /* BigInt path (well beyond int64). */
    assert_eval_eq("NonPositive[10^40]", "False", 0);
    assert_eval_eq("NonPositive[-10^40]", "True", 0);
}

void test_nonpositive_rationals() {
    assert_eval_eq("NonPositive[3/4]", "False", 0);
    assert_eval_eq("NonPositive[-3/4]", "True", 0);
    assert_eval_eq("NonPositive[-1/1000000000000]", "True", 0);
    assert_eval_eq("NonPositive[7/3]", "False", 0);
}

/* ---- machine reals --------------------------------------------------- */

void test_nonpositive_machine_reals() {
    assert_eval_eq("NonPositive[1.6]", "False", 0);
    assert_eval_eq("NonPositive[-2.5]", "True", 0);
    assert_eval_eq("NonPositive[0.0]", "True", 0);
    assert_eval_eq("NonPositive[-1.0*10^-12]", "True", 0);
}

/* ---- arbitrary precision (MPFR) -------------------------------------- */

void test_nonpositive_arbitrary_precision() {
    assert_eval_eq("NonPositive[N[Pi, 50]]", "False", 0);
    assert_eval_eq("NonPositive[-N[Pi, 50]]", "True", 0);
    assert_eval_eq("NonPositive[-N[E, 100]]", "True", 0);
    /* Difference that is a tiny but genuinely negative high-precision real. */
    assert_eval_eq("NonPositive[3 - N[Pi, 50]]", "True", 0);
}

/* ---- symbolic constants & numeric functions -------------------------- */

void test_nonpositive_constants_and_functions() {
    assert_eval_eq("NonPositive[Pi]", "False", 0);
    assert_eval_eq("NonPositive[E]", "False", 0);
    assert_eval_eq("NonPositive[Sqrt[2]]", "False", 0);
    assert_eval_eq("NonPositive[E^2]", "False", 0);
    assert_eval_eq("NonPositive[Log[2]]", "False", 0);
    assert_eval_eq("NonPositive[Log[1/2]]", "True", 0);    /* log(1/2) ~ -0.693 */
    /* Spec example: Sin[10^5] ~ 0.0357 > 0, so NonPositive is False. */
    assert_eval_eq("NonPositive[Sin[10^5]]", "False", 0);
    assert_eval_eq("NonPositive[Sin[4]]", "True", 0);      /* sin(4) ~ -0.7568 */
    /* Exact numeric combinations. */
    assert_eval_eq("NonPositive[1 - Pi]", "True", 0);
    assert_eval_eq("NonPositive[Pi - 3]", "False", 0);
}

/* ---- complex values are not non-positive ----------------------------- */

void test_nonpositive_complex() {
    assert_eval_eq("NonPositive[I]", "False", 0);
    assert_eval_eq("NonPositive[1 + I]", "False", 0);
    assert_eval_eq("NonPositive[2 - 3 I]", "False", 0);
    assert_eval_eq("NonPositive[Sqrt[-2]]", "False", 0);
    /* A Complex literal whose imaginary part collapses to zero is real. */
    assert_eval_eq("NonPositive[-2 + 0 I]", "True", 0);
    assert_eval_eq("NonPositive[2 + 0 I]", "False", 0);
}

/* ---- non-numeric arguments stay unevaluated -------------------------- */

void test_nonpositive_symbolic_unevaluated() {
    assert_eval_eq("NonPositive[x]", "NonPositive[x]", 0);
    assert_eval_eq("NonPositive[Sin[y]]", "NonPositive[Sin[y]]", 0);
    assert_eval_eq("NonPositive[a + b]", "NonPositive[a + b]", 0);
    /* FullForm confirms the head/arg survive intact. */
    assert_eval_eq("NonPositive[x]", "NonPositive[x]", 1);
}

/* ---- Listable threading (spec example) ------------------------------- */

void test_nonpositive_listable() {
    assert_eval_eq("NonPositive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]",
                   "{False, False, False, True, True, False, False}", 0);
    assert_eval_eq("NonPositive[{x, Sin[y]}]",
                   "{NonPositive[x], NonPositive[Sin[y]]}", 0);
    assert_eval_eq("NonPositive[{}]", "{}", 0);
    /* Nested lists thread element by element. */
    assert_eval_eq("NonPositive[{{1, -1}, {0, -2}}]",
                   "{{False, True}, {True, True}}", 0);
}

/* ---- attributes ------------------------------------------------------ */

void test_nonpositive_attributes() {
    SymbolDef* def = symtab_get_def("NonPositive");
    ASSERT(def != NULL);
    ASSERT_MSG(def->attributes & ATTR_LISTABLE,
               "NonPositive must be Listable");
    ASSERT_MSG(def->attributes & ATTR_PROTECTED,
               "NonPositive must be Protected");
    /* Attributes[NonPositive] should report exactly these two. */
    assert_eval_eq("Attributes[NonPositive]", "{Listable, Protected}", 0);
}

/* ---- docstring ------------------------------------------------------- */

void test_nonpositive_docstring() {
    const char* doc = symtab_get_docstring("NonPositive");
    ASSERT_MSG(doc != NULL && doc[0] != '\0',
               "NonPositive must have a non-empty docstring");
}

/* ---- arity: argx diagnostic + unevaluated result --------------------- */

void test_nonpositive_arity() {
    /* Wrong argument counts emit NonPositive::argx on stderr and leave the
     * expression unevaluated. We only check the returned (unevaluated) form
     * here; the diagnostic itself goes to stderr. */
    assert_eval_eq("NonPositive[]", "NonPositive[]", 0);
    assert_eval_eq("NonPositive[1, 2, 3]", "NonPositive[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_nonpositive_exact_integers);
    TEST(test_nonpositive_rationals);
    TEST(test_nonpositive_machine_reals);
    TEST(test_nonpositive_arbitrary_precision);
    TEST(test_nonpositive_constants_and_functions);
    TEST(test_nonpositive_complex);
    TEST(test_nonpositive_symbolic_unevaluated);
    TEST(test_nonpositive_listable);
    TEST(test_nonpositive_attributes);
    TEST(test_nonpositive_docstring);
    TEST(test_nonpositive_arity);

    symtab_clear();
    printf("All NonPositive tests passed.\n");
    return 0;
}

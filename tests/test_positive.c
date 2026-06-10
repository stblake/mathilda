/* Tests for Positive[x].
 *
 * Positive[x] gives True when x is a positive real number, False when x is a
 * manifestly non-positive numeric quantity (negative real, zero, or a non-real
 * complex value), and is left unevaluated for non-numeric arguments. Covers
 * exact integers/rationals/bigints, machine and arbitrary-precision (MPFR)
 * reals, symbolic constants, numeric functions, complex values, Listable
 * threading, attributes, and arity diagnostics. Reference values cross-checked
 * against the running binary and the task specification. */

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

void test_positive_exact_integers() {
    assert_eval_eq("Positive[1]", "True", 0);
    assert_eval_eq("Positive[5]", "True", 0);
    assert_eval_eq("Positive[-5]", "False", 0);
    assert_eval_eq("Positive[0]", "False", 0);
    /* BigInt path (well beyond int64). */
    assert_eval_eq("Positive[10^40]", "True", 0);
    assert_eval_eq("Positive[-10^40]", "False", 0);
}

void test_positive_rationals() {
    assert_eval_eq("Positive[3/4]", "True", 0);
    assert_eval_eq("Positive[-3/4]", "False", 0);
    assert_eval_eq("Positive[-1/1000000000000]", "False", 0);
    assert_eval_eq("Positive[7/3]", "True", 0);
}

/* ---- machine reals --------------------------------------------------- */

void test_positive_machine_reals() {
    assert_eval_eq("Positive[1.6]", "True", 0);
    assert_eval_eq("Positive[-2.5]", "False", 0);
    assert_eval_eq("Positive[0.0]", "False", 0);
    assert_eval_eq("Positive[1.0*10^-12]", "True", 0);
}

/* ---- arbitrary precision (MPFR) -------------------------------------- */

void test_positive_arbitrary_precision() {
    assert_eval_eq("Positive[N[Pi, 50]]", "True", 0);
    assert_eval_eq("Positive[-N[Pi, 50]]", "False", 0);
    assert_eval_eq("Positive[N[E, 100]]", "True", 0);
    /* Difference that is a tiny but genuinely positive high-precision real. */
    assert_eval_eq("Positive[N[Pi, 50] - 3]", "True", 0);
}

/* ---- symbolic constants & numeric functions -------------------------- */

void test_positive_constants_and_functions() {
    assert_eval_eq("Positive[Pi]", "True", 0);
    assert_eval_eq("Positive[E]", "True", 0);
    assert_eval_eq("Positive[Sqrt[2]]", "True", 0);
    assert_eval_eq("Positive[E^2]", "True", 0);
    assert_eval_eq("Positive[Log[2]]", "True", 0);
    assert_eval_eq("Positive[Log[1/2]]", "False", 0);
    /* Spec example: Sin[10^5] ~ 0.0357 > 0. */
    assert_eval_eq("Positive[Sin[10^5]]", "True", 0);
    assert_eval_eq("Positive[Sin[4]]", "False", 0);   /* sin(4) ~ -0.7568 */
    /* An exact numeric combination that is negative. */
    assert_eval_eq("Positive[1 - Pi]", "False", 0);
    assert_eval_eq("Positive[Pi - 3]", "True", 0);
}

/* ---- complex values are not positive --------------------------------- */

void test_positive_complex() {
    assert_eval_eq("Positive[I]", "False", 0);
    assert_eval_eq("Positive[1 + I]", "False", 0);
    assert_eval_eq("Positive[2 - 3 I]", "False", 0);
    assert_eval_eq("Positive[Sqrt[-2]]", "False", 0);
    /* A Complex literal whose imaginary part collapses to zero is real. */
    assert_eval_eq("Positive[2 + 0 I]", "True", 0);
    assert_eval_eq("Positive[-2 + 0 I]", "False", 0);
}

/* ---- non-numeric arguments stay unevaluated -------------------------- */

void test_positive_symbolic_unevaluated() {
    assert_eval_eq("Positive[x]", "Positive[x]", 0);
    assert_eval_eq("Positive[Sin[y]]", "Positive[Sin[y]]", 0);
    assert_eval_eq("Positive[a + b]", "Positive[a + b]", 0);
    /* FullForm confirms the head/arg survive intact. */
    assert_eval_eq("Positive[x]", "Positive[x]", 1);
}

/* ---- Listable threading (spec example) ------------------------------- */

void test_positive_listable() {
    assert_eval_eq("Positive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]",
                   "{True, True, True, False, False, False, True}", 0);
    assert_eval_eq("Positive[{x, Sin[y]}]",
                   "{Positive[x], Positive[Sin[y]]}", 0);
    assert_eval_eq("Positive[{}]", "{}", 0);
    /* Nested lists thread element by element. */
    assert_eval_eq("Positive[{{1, -1}, {0, 2}}]",
                   "{{True, False}, {False, True}}", 0);
}

/* ---- attributes ------------------------------------------------------ */

void test_positive_attributes() {
    SymbolDef* def = symtab_get_def("Positive");
    ASSERT(def != NULL);
    ASSERT_MSG(def->attributes & ATTR_LISTABLE,
               "Positive must be Listable");
    ASSERT_MSG(def->attributes & ATTR_PROTECTED,
               "Positive must be Protected");
    /* Attributes[Positive] should report exactly these two. */
    assert_eval_eq("Attributes[Positive]", "{Listable, Protected}", 0);
}

/* ---- docstring ------------------------------------------------------- */

void test_positive_docstring() {
    const char* doc = symtab_get_docstring("Positive");
    ASSERT_MSG(doc != NULL && doc[0] != '\0',
               "Positive must have a non-empty docstring");
}

/* ---- arity: argx diagnostic + unevaluated result --------------------- */

void test_positive_arity() {
    /* Wrong argument counts emit Positive::argx on stderr and leave the
     * expression unevaluated. We only check the returned (unevaluated) form
     * here; the diagnostic itself goes to stderr. */
    assert_eval_eq("Positive[]", "Positive[]", 0);
    assert_eval_eq("Positive[1, 2, 3]", "Positive[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_positive_exact_integers);
    TEST(test_positive_rationals);
    TEST(test_positive_machine_reals);
    TEST(test_positive_arbitrary_precision);
    TEST(test_positive_constants_and_functions);
    TEST(test_positive_complex);
    TEST(test_positive_symbolic_unevaluated);
    TEST(test_positive_listable);
    TEST(test_positive_attributes);
    TEST(test_positive_docstring);
    TEST(test_positive_arity);

    symtab_clear();
    printf("All Positive tests passed.\n");
    return 0;
}

/* Tests for Negative[x].
 *
 * Negative[x] gives True when x is a negative real number, False when x is a
 * manifestly non-negative numeric quantity (positive real, zero, or a non-real
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

void test_negative_exact_integers() {
    assert_eval_eq("Negative[1]", "False", 0);
    assert_eval_eq("Negative[5]", "False", 0);
    assert_eval_eq("Negative[-5]", "True", 0);
    assert_eval_eq("Negative[0]", "False", 0);
    /* BigInt path (well beyond int64). */
    assert_eval_eq("Negative[10^40]", "False", 0);
    assert_eval_eq("Negative[-10^40]", "True", 0);
}

void test_negative_rationals() {
    assert_eval_eq("Negative[3/4]", "False", 0);
    assert_eval_eq("Negative[-3/4]", "True", 0);
    assert_eval_eq("Negative[-1/1000000000000]", "True", 0);
    assert_eval_eq("Negative[7/3]", "False", 0);
}

/* ---- machine reals --------------------------------------------------- */

void test_negative_machine_reals() {
    assert_eval_eq("Negative[1.6]", "False", 0);
    assert_eval_eq("Negative[-2.5]", "True", 0);
    assert_eval_eq("Negative[0.0]", "False", 0);
    assert_eval_eq("Negative[-1.0*10^-12]", "True", 0);
}

/* ---- arbitrary precision (MPFR) -------------------------------------- */

void test_negative_arbitrary_precision() {
    assert_eval_eq("Negative[N[Pi, 50]]", "False", 0);
    assert_eval_eq("Negative[-N[Pi, 50]]", "True", 0);
    assert_eval_eq("Negative[-N[E, 100]]", "True", 0);
    /* Difference that is a tiny but genuinely negative high-precision real. */
    assert_eval_eq("Negative[3 - N[Pi, 50]]", "True", 0);
}

/* ---- symbolic constants & numeric functions -------------------------- */

void test_negative_constants_and_functions() {
    assert_eval_eq("Negative[Pi]", "False", 0);
    assert_eval_eq("Negative[E]", "False", 0);
    assert_eval_eq("Negative[Sqrt[2]]", "False", 0);
    assert_eval_eq("Negative[E^2]", "False", 0);
    assert_eval_eq("Negative[Log[2]]", "False", 0);
    assert_eval_eq("Negative[Log[1/2]]", "True", 0);   /* log(1/2) ~ -0.693 */
    /* Spec example: Sin[10^5] ~ 0.0357 > 0, so Negative is False. */
    assert_eval_eq("Negative[Sin[10^5]]", "False", 0);
    assert_eval_eq("Negative[Sin[4]]", "True", 0);     /* sin(4) ~ -0.7568 */
    /* Exact numeric combinations. */
    assert_eval_eq("Negative[1 - Pi]", "True", 0);
    assert_eval_eq("Negative[Pi - 3]", "False", 0);
}

/* ---- complex values are not negative --------------------------------- */

void test_negative_complex() {
    assert_eval_eq("Negative[I]", "False", 0);
    assert_eval_eq("Negative[1 + I]", "False", 0);
    assert_eval_eq("Negative[2 - 3 I]", "False", 0);
    assert_eval_eq("Negative[Sqrt[-2]]", "False", 0);
    /* A Complex literal whose imaginary part collapses to zero is real. */
    assert_eval_eq("Negative[-2 + 0 I]", "True", 0);
    assert_eval_eq("Negative[2 + 0 I]", "False", 0);
}

/* ---- non-numeric arguments stay unevaluated -------------------------- */

void test_negative_symbolic_unevaluated() {
    assert_eval_eq("Negative[x]", "Negative[x]", 0);
    assert_eval_eq("Negative[Sin[y]]", "Negative[Sin[y]]", 0);
    assert_eval_eq("Negative[a + b]", "Negative[a + b]", 0);
    /* FullForm confirms the head/arg survive intact. */
    assert_eval_eq("Negative[x]", "Negative[x]", 1);
}

/* ---- Listable threading (spec example) ------------------------------- */

void test_negative_listable() {
    assert_eval_eq("Negative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]",
                   "{False, False, False, False, True, False, False}", 0);
    assert_eval_eq("Negative[{x, Sin[y]}]",
                   "{Negative[x], Negative[Sin[y]]}", 0);
    assert_eval_eq("Negative[{}]", "{}", 0);
    /* Nested lists thread element by element. */
    assert_eval_eq("Negative[{{1, -1}, {0, -2}}]",
                   "{{False, True}, {False, True}}", 0);
}

/* ---- attributes ------------------------------------------------------ */

void test_negative_attributes() {
    SymbolDef* def = symtab_get_def("Negative");
    ASSERT(def != NULL);
    ASSERT_MSG(def->attributes & ATTR_LISTABLE,
               "Negative must be Listable");
    ASSERT_MSG(def->attributes & ATTR_PROTECTED,
               "Negative must be Protected");
    /* Attributes[Negative] should report exactly these two. */
    assert_eval_eq("Attributes[Negative]", "{Listable, Protected}", 0);
}

/* ---- docstring ------------------------------------------------------- */

void test_negative_docstring() {
    const char* doc = symtab_get_docstring("Negative");
    ASSERT_MSG(doc != NULL && doc[0] != '\0',
               "Negative must have a non-empty docstring");
}

/* ---- arity: argx diagnostic + unevaluated result --------------------- */

void test_negative_arity() {
    /* Wrong argument counts emit Negative::argx on stderr and leave the
     * expression unevaluated. We only check the returned (unevaluated) form
     * here; the diagnostic itself goes to stderr. */
    assert_eval_eq("Negative[]", "Negative[]", 0);
    assert_eval_eq("Negative[1, 2, 3]", "Negative[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_negative_exact_integers);
    TEST(test_negative_rationals);
    TEST(test_negative_machine_reals);
    TEST(test_negative_arbitrary_precision);
    TEST(test_negative_constants_and_functions);
    TEST(test_negative_complex);
    TEST(test_negative_symbolic_unevaluated);
    TEST(test_negative_listable);
    TEST(test_negative_attributes);
    TEST(test_negative_docstring);
    TEST(test_negative_arity);

    symtab_clear();
    printf("All Negative tests passed.\n");
    return 0;
}

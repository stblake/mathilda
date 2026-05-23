#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for Element[x, dom]. The builtin reads the current $Assumptions
 * via the OwnValue (NOT via evaluate, which would re-fire builtin_element
 * recursively) and decides:
 *   - True   if dom membership follows from the assumption set or a
 *            structural / numeric literal.
 *   - False  if the value is structurally outside the domain.
 *   - stays unevaluated otherwise.
 *
 * The tests cover (a) numeric and structural literal decisions, (b)
 * $Assumptions-driven decisions visible inside Assuming[...] bodies,
 * (c) the lattice Integer => Rational => Algebraic => Real => Complex,
 * (d) list-threading in the first argument, and (e) the no-decision
 * case where the result remains unevaluated.
 */

/* ---- Numeric / structural literals ---- */

void test_element_int_in_integers(void)        { assert_eval_eq("Element[5, Integers]", "True", 0); }
void test_element_int_in_reals(void)           { assert_eval_eq("Element[-3, Reals]", "True", 0); }
void test_element_int_in_complexes(void)       { assert_eval_eq("Element[42, Complexes]", "True", 0); }
void test_element_rat_not_in_integers(void)    { assert_eval_eq("Element[5/2, Integers]", "False", 0); }
void test_element_rat_in_rationals(void)       { assert_eval_eq("Element[5/2, Rationals]", "True", 0); }
void test_element_rat_in_reals(void)           { assert_eval_eq("Element[5/2, Reals]", "True", 0); }
void test_element_real_in_reals(void)          { assert_eval_eq("Element[2.5, Reals]", "True", 0); }
void test_element_real_not_in_integers(void)   { assert_eval_eq("Element[2.5, Integers]", "False", 0); }
void test_element_complex_in_complexes(void)   { assert_eval_eq("Element[1 + I, Complexes]", "True", 0); }
void test_element_complex_not_in_reals(void)   { assert_eval_eq("Element[1 + I, Reals]", "False", 0); }
void test_element_complex_not_in_integers(void){ assert_eval_eq("Element[1 + I, Integers]", "False", 0); }

void test_element_true_in_booleans(void)       { assert_eval_eq("Element[True, Booleans]", "True", 0); }
void test_element_false_in_booleans(void)      { assert_eval_eq("Element[False, Booleans]", "True", 0); }

/* PrimeQ-driven literals. */
void test_element_prime_int(void)              { assert_eval_eq("Element[7, Primes]", "True", 0); }
void test_element_composite_int(void)          { assert_eval_eq("Element[15, Primes]", "False", 0); }
void test_element_composite_in_composites(void){ assert_eval_eq("Element[15, Composites]", "True", 0); }

/* ---- Symbolic without assumptions: stays unevaluated ---- */

void test_element_sym_no_assumption(void) {
    assert_eval_eq("Element[x, Reals]", "Element[x, Reals]", 0);
}

/* ---- $Assumptions integration via Assuming ---- */

void test_element_under_reals(void) {
    assert_eval_eq("Assuming[Element[x, Reals], Element[x, Reals]]", "True", 0);
}

void test_element_lattice_integer_to_real(void) {
    /* Integer => Real */
    assert_eval_eq("Assuming[Element[x, Integers], Element[x, Reals]]", "True", 0);
}

void test_element_lattice_real_to_complex(void) {
    /* Real => Complex */
    assert_eval_eq("Assuming[Element[x, Reals], Element[x, Complexes]]", "True", 0);
}

void test_element_lattice_integer_to_rational(void) {
    assert_eval_eq("Assuming[Element[x, Integers], Element[x, Rationals]]", "True", 0);
}

void test_element_lattice_integer_to_algebraic(void) {
    assert_eval_eq("Assuming[Element[x, Integers], Element[x, Algebraics]]", "True", 0);
}

void test_element_inequality_implies_real(void) {
    /* x > 0 implies x in Reals (the prov_re lattice goes
     * positive => non-negative => real). */
    assert_eval_eq("Assuming[x > 0, Element[x, Reals]]", "True", 0);
}

void test_element_dollar_assumptions_restored(void) {
    /* After Assuming exits, the global $Assumptions returns to True and
     * symbolic Element queries drop back to unevaluated form. */
    assert_eval_eq("Assuming[Element[x, Reals], Element[x, Reals]]", "True", 0);
    assert_eval_eq("$Assumptions", "True", 0);
    assert_eval_eq("Element[x, Reals]", "Element[x, Reals]", 0);
}

/* ---- List shorthand: Element[{x1,...,xN}, dom] is a conjunction ---- */

void test_element_threads_over_list(void) {
    /* Element[{x1,...,xN}, dom] is shorthand for the conjunction of
     * Element[xi, dom]; reduces to a single True/False only when every
     * component decides, otherwise stays unchanged so downstream
     * consumers (e.g. Simplify's AssumeCtx) can treat it as a joint
     * assumption on x1..xN. */
    assert_eval_eq("Element[{1, 5/2, 1 + I}, Integers]", "False", 0);
    assert_eval_eq("Element[{1, 5/2, 1 + I}, Reals]",    "False", 0);
    assert_eval_eq("Element[{1, 5/2, 1 + I}, Complexes]", "True", 0);
    /* Undecided components: leave the expression in place. */
    assert_eval_eq("Element[{x, y}, Reals]",
                   "Element[{x, y}, Reals]", 0);
    /* Simplify treats Element[{x,y}, Reals] as Element[x,Reals] &&
     * Element[y,Reals], so Log[E^(x+y)] collapses to x+y. */
    assert_eval_eq("Simplify[Log[E^(x+y)], Element[{x,y}, Reals]]",
                   "x + y", 0);
}

/* ---- Direct fact lookup, even for compound expressions ---- */

void test_element_direct_compound_fact(void) {
    /* If the assumption set names a non-symbol expression, a verbatim
     * Element[..., dom] query should still hit the direct-fact lookup. */
    assert_eval_eq("Assuming[Element[f[x], Reals], Element[f[x], Reals]]",
                   "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_element_int_in_integers);
    TEST(test_element_int_in_reals);
    TEST(test_element_int_in_complexes);
    TEST(test_element_rat_not_in_integers);
    TEST(test_element_rat_in_rationals);
    TEST(test_element_rat_in_reals);
    TEST(test_element_real_in_reals);
    TEST(test_element_real_not_in_integers);
    TEST(test_element_complex_in_complexes);
    TEST(test_element_complex_not_in_reals);
    TEST(test_element_complex_not_in_integers);
    TEST(test_element_true_in_booleans);
    TEST(test_element_false_in_booleans);
    TEST(test_element_prime_int);
    TEST(test_element_composite_int);
    TEST(test_element_composite_in_composites);
    TEST(test_element_sym_no_assumption);
    TEST(test_element_under_reals);
    TEST(test_element_lattice_integer_to_real);
    TEST(test_element_lattice_real_to_complex);
    TEST(test_element_lattice_integer_to_rational);
    TEST(test_element_lattice_integer_to_algebraic);
    TEST(test_element_inequality_implies_real);
    TEST(test_element_dollar_assumptions_restored);
    TEST(test_element_threads_over_list);
    TEST(test_element_direct_compound_fact);

    printf("All Element tests passed!\n");
    return 0;
}

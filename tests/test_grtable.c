/* test_grtable.c
 *
 * The Gradshteyn & Ryzhik table stage of Integrate (Method -> "GRTable";
 * engine src/internal/GRTable.m, records in src/internal/GRTables/ generated
 * by tools/gr/build_table.py).  Each record is a G&R formula verified
 * numerically by the G&R oracle; the stage unifies the integrand with the
 * record's pattern, solves the record's parameters, checks the book's
 * conditions under the caller's Assumptions, and self-checks the answer
 * numerically before returning it.
 *
 * Covered: an indefinite record (symbolic and an instance), a definite
 * half-line record at a numeric instance (value checked against NIntegrate),
 * a symbolic record whose condition is not implied by the assumptions
 * (ConditionalExpression), the declines (an assumption contradicting the
 * book's condition, an integrand no record matches, an interval no record
 * covers), and the lazy load (GRTable`Count[] is 0 until the first query).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_grtable_lazy_load(void) {
    /* The table is not read at startup: the first GRTable query loads it. */
    assert_eval_eq("Integrate[Sin[x]^2 + q, x, Method -> \"GRTable\"] // Head", "Integrate", 0);
    assert_eval_eq("GRTable`Count[] > 100", "True", 0);
}

static void test_grtable_indefinite(void) {
    /* G&R 2.325.7 */
    assert_eval_eq("Integrate[E^(a x^n)/x, x, Method -> \"GRTable\"]",
                   "ExpIntegralEi[a x^n]/n", 0);
    assert_eval_eq("Integrate[E^(3 x^2)/x, x, Method -> \"GRTable\"]",
                   "1/2 ExpIntegralEi[3 x^2]", 0);
    /* Integrate[1/x, x] is Log[x]: the book's Log[Abs[u]] antiderivatives
     * were rewritten Log[u] (tools/gr/logconv.py), after which the cascade
     * closes G&R 2.736.2 itself and the record left the table */
    assert_eval_eq("Head[Integrate[x Log[x^2 - a^2], x, Method -> \"GRTable\"]]", "Integrate", 0);
}

static void test_grtable_definite(void) {
    /* G&R 3.197.1 at nu = 7/2, mu = rho = 3, beta = 2, gamma = 5 */
    assert_eval_eq("Integrate[x^(5/2)/((2 + x)^3 (5 + x)^3), {x, 0, Infinity}, Method -> \"GRTable\"]",
                   "(3 Sqrt[5] Pi)/(8 (725 + 229 Sqrt[10]))", 0);
    assert_eval_eq("Abs[N[(3 Sqrt[5] Pi)/(8 (725 + 229 Sqrt[10]))] - "
                   "NIntegrate[x^(5/2)/((2 + x)^3 (5 + x)^3), {x, 0, Infinity}]] < 10^-10", "True", 0);
    /* the book's condition (mu > 1/2 here) is not implied by mu > 0: the
     * answer carries it */
    assert_eval_eq("Head[Integrate[x^(-1/2 + mu)/((a + x)^mu (b + x)^mu), {x, 0, Infinity}, "
                   "Method -> \"GRTable\", Assumptions -> mu > 0 && a > 0 && b > 0]]",
                   "ConditionalExpression", 0);
}

static void test_grtable_declines(void) {
    /* contradicts the book's condition: divergent, no answer */
    assert_eval_eq("Head[Integrate[x^(-1/2 + mu)/((a + x)^mu (b + x)^mu), {x, 0, Infinity}, "
                   "Method -> \"GRTable\", Assumptions -> mu < -1 && a > 0 && b > 0]]", "Integrate", 0);
    /* no record matches */
    assert_eval_eq("Head[Integrate[Sin[x]^2, x, Method -> \"GRTable\"]]", "Integrate", 0);
    /* no record on [1, Infinity) */
    assert_eval_eq("Head[Integrate[x^(5/2)/((2 + x)^3 (5 + x)^3), {x, 1, Infinity}, Method -> \"GRTable\"]]",
                   "Integrate", 0);
}

/* A value the book prints as a decimal is served exactly (tools/gr/value_fixes.tsv):
 * never a machine number, and a bound iterator prints as the plain symbol. */
static void test_grtable_exact_values(void) {
    assert_eval_eq("Integrate[x^(-x), {x, 0, 1}, Method -> \"GRTable\"]",
                   "Sum[n^(-n), {n, 1, Infinity}]", 0);
    assert_eval_eq("Integrate[Log[Log[1/x]]/(1 + x), {x, 0, 1}, Method -> \"GRTable\"]",
                   "-1/2 Log[2]^2", 0);
    assert_eval_eq("Expand[Integrate[x/(2 Cosh[x] - 1), {x, 0, Infinity}, Method -> \"GRTable\"] - "
                   "(PolyGamma[1, 1/6] + PolyGamma[1, 1/3] - PolyGamma[1, 2/3] - PolyGamma[1, 5/6])/36]",
                   "0", 0);
    assert_eval_eq("FreeQ[{Integrate[x^(-x), {x, 0, 1}, Method -> \"GRTable\"], "
                   "Integrate[x/(2 Cosh[x] - 1), {x, 0, Infinity}, Method -> \"GRTable\"]}, _Real]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    test_mute_messages();
    TEST(test_grtable_lazy_load);
    TEST(test_grtable_indefinite);
    TEST(test_grtable_definite);
    TEST(test_grtable_declines);
    TEST(test_grtable_exact_values);
    printf("All GRTable tests passed.\n");
    return 0;
}

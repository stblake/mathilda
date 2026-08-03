#include "sort.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e); expr_free(res);
}

void test_sort_canonical() {
    run_test("Sort[{d, b, c, a}]", "List[a, b, c, d]");
    run_test("Sort[{\"cat\", \"fish\", \"catfish\", \"Cat\"}]", "List[\"cat\", \"Cat\", \"catfish\", \"fish\"]");
    run_test("Sort[{y, x^2, x + y, y^3}]", "List[Power[x, 2], y, Power[y, 3], Plus[x, y]]");
    run_test("Sort[f[b, a, c]]", "f[a, b, c]");
    // Numbers
    run_test("Sort[{3, 1, 2}]", "List[1, 2, 3]");
    run_test("Sort[{1.5, 1, 2}]", "List[1, 1.5, 2]");
}

void test_sort_custom() {
    // Sort by numerical value with Less
    // Note: Pi and E are symbols, Sqrt[2] is Power[2, 1/2]
    // Order: 1, Sqrt[2], 2, E, 3, Pi
    // Numerically: 1, 1.414, 2, 2.718, 3, 3.14159
    // Wait, Less currently works on numbers. If they are symbols it might not evaluate to True/False.
    // I need to make sure Pi and E have numerical values if I want Less to work, 
    // or use N[Pi] etc. But Sort[list, Less] in Mathematica works if they are numeric.
    
    // For now test with explicit numbers
    run_test("Sort[{3, 1, 2, 1.5}, Less]", "List[1, 1.5, 2, 3]");
    run_test("Sort[{3, 1, 2, 1.5}, Greater]", "List[3, 2, 1.5, 1]");
    
    // Sort by comparing second part
    run_test("Sort[{{a, 2}, {c, 1}, {d, 3}}, #1[[2]] < #2[[2]] &]", 
             "List[List[c, 1], List[a, 2], List[d, 3]]");
}


void test_orderedq() {
    run_test("OrderedQ[{e, e}]", "True");
    run_test("OrderedQ[{1, 4, 2}]", "False");
    run_test("OrderedQ[{\"cat\", \"catfish\", \"fish\"}]", "True");
    run_test("OrderedQ[{1, Sqrt[2], 2, E, 3, Pi}]", "False");
    run_test("OrderedQ[{1, Sqrt[2], 2, E, 3, Pi}, Less]", "True");
    run_test("OrderedQ[{{a, 2}, {c, 1}, {d, 3}}, #1[[2]] < #2[[2]] &]", "False");
    run_test("OrderedQ[{x, y, x + y}]", "True");
    run_test("OrderedQ[f[b, a, c]]", "False");
    run_test("OrderedQ[{4, 3, 2, 1}, Greater]", "True");
    run_test("OrderedQ[{4, 3, 3, 1}, Greater]", "False");
    run_test("OrderedQ[{4, 3, 3, 1}, GreaterEqual]", "True");
}

void test_order() {
    /* Basics: the defining triple. */
    run_test("Order[a, a]", "0");
    run_test("Order[a, b]", "1");
    run_test("Order[b, a]", "-1");
    run_test("{Order[a, a], Order[a, b], Order[b, a]}", "List[0, 1, -1]");

    /* Structural, NOT numerical: Integer sorts before the symbol Pi, but two
     * numeric atoms compare by value. */
    run_test("Order[6, Pi]", "1");
    run_test("Order[6, N[Pi]]", "-1");

    /* The Order@@@ example from the spec, over the 9 tuples of {0,1,2}.
     * Uses the explicit tuple list (what Tuples[{0,1,2},2] yields) so the sort
     * test binary's leak surface does not pull in an unrelated pre-existing
     * Tuples allocation; the literal `Order @@@ Tuples[{0,1,2},2]` form is
     * exercised in the REPL smoke test. */
    run_test("Order @@@ {{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,1},{2,2}}",
             "List[0, 1, 1, -1, 0, 1, -1, -1, 0]");

    /* Numeric tier: by value, mixed Integer/Real. */
    run_test("Order[1, 2]", "1");
    run_test("Order[2, 1]", "-1");
    run_test("Order[1, 1.5]", "1");
    run_test("Order[2, 1.5]", "-1");
    run_test("Order[3, 3]", "0");

    /* Cross-tier ordering (numeric < string < symbol < compound), matching
     * Sort's canonical order. */
    run_test("Order[5, \"a\"]", "1");        /* number before string */
    run_test("Order[\"a\", x]", "1");        /* string before symbol */
    run_test("Order[\"a\", \"b\"]", "1");    /* strings lexicographic */
    run_test("Order[x, f[x]]", "1");         /* bare symbol before compound */
    run_test("Order[x, x^2]", "1");          /* x before x^2 (polynomial order) */
    run_test("Order[x^2, x]", "-1");

    /* Lists compare element-wise; identical compounds give 0. */
    run_test("Order[{1, 2}, {1, 3}]", "1");
    run_test("Order[{1, 3}, {1, 2}]", "-1");
    run_test("Order[f[a, b], f[a, b]]", "0");
    run_test("Order[f[a], f[a, b]]", "1");   /* shorter before longer */

    /* Arity: wrong argument count leaves Order unevaluated. */
    run_test("Order[a]", "Order[a]");
    run_test("Order[a, b, c]", "Order[a, b, c]");
}

void test_ordering() {
    /* Basic permutation and the defining identity list[[Ordering[list]]] == Sort[list]. */
    run_test("Ordering[{c, a, b}]", "List[2, 3, 1]");
    run_test("{c, a, b}[[Ordering[{c, a, b}]]]", "List[a, b, c]");
    run_test("Ordering[{2, 6, 1, 9, 1, 2, 3}]", "List[3, 5, 1, 6, 7, 2, 4]");
    run_test("{2, 6, 1, 9, 1, 2, 3}[[Ordering[{2, 6, 1, 9, 1, 2, 3}]]] === Sort[{2, 6, 1, 9, 1, 2, 3}]",
             "True");

    /* n smallest / n largest, and the empty selection. */
    run_test("Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]", "List[3, 5, 1, 6]");
    run_test("Ordering[{2, 6, 1, 9, 1, 2, 3}, -1]", "List[4]");
    run_test("Ordering[{4, 5, -2, 3}, -1]", "List[2]");
    run_test("Ordering[{2, 6, 1, 9, 1, 2, 3}, 0]", "List[]");

    /* seq is a full Take spec: a {m, n} span, UpTo (clamped), and All. */
    run_test("Ordering[{2, 6, 1, 9, 1, 2, 3}, {4, -1}]", "List[6, 7, 2, 4]");
    run_test("Ordering[{2, 6, 1, 9, 2}, UpTo[6]]", "List[3, 1, 5, 2, 4]");
    run_test("Ordering[{3, 1, 2}, All]", "List[2, 3, 1]");

    /* Custom ordering function. With no ties the permutation is deterministic. */
    run_test("Ordering[{3, 1, 2}, All, Greater]", "List[1, 3, 2]");
    run_test("Ordering[{3, 1, 2}, 1, Greater]", "List[1]");
    /* With ties under a custom (strict) comparator, only the defining invariant
     * list[[ord]] == Sort[list, p] is guaranteed -- the tie order is qsort-defined,
     * exactly as it is for Sort[list, p]. */
    run_test("{2, 6, 1, 9, 1, 2, 3}[[Ordering[{2, 6, 1, 9, 1, 2, 3}, All, Greater]]] "
             "=== Sort[{2, 6, 1, 9, 1, 2, 3}, Greater]", "True");

    /* Inverse permutation. */
    run_test("Ordering[{4, 5, 1, 2, 3}]", "List[3, 4, 5, 1, 2]");

    /* Any head, and an Association (ordered by its values). */
    run_test("Ordering[f[3, 1, 2]]", "List[2, 3, 1]");
    run_test("Ordering[<|1 -> c, 2 -> a, 3 -> b|>]", "List[2, 3, 1]");

    /* Stability: ties break by original position (documented for the default order). */
    run_test("Ordering[{2, 2, 1}]", "List[3, 1, 2]");
    run_test("Ordering[{1, 1, 1}]", "List[1, 2, 3]");
    run_test("Ordering[{5, 5, 5, 5}, 1]", "List[1]");

    /* Reals and mixed exact/inexact. */
    run_test("Ordering[{1.5, 1, 2}]", "List[2, 1, 3]");

    /* Edge cases: empty and singleton. */
    run_test("Ordering[{}]", "List[]");
    run_test("Ordering[{5}]", "List[1]");

    /* Atom -> unevaluated; too many arguments -> unevaluated. */
    run_test("Ordering[5]", "Ordering[5]");
    run_test("Ordering[a, b, c, d]", "Ordering[a, b, c, d]");
}

int main() {
    symtab_init();
    core_init();

    // Define Pi and E for numerical tests if needed
    // symtab_add_own_value("Pi", expr_new_symbol("Pi"), expr_new_real(3.141592653589793));
    // symtab_get_def("Pi")->attributes |= ATTR_NUMERICFUNCTION;

    TEST(test_sort_canonical);
    TEST(test_sort_custom);
    TEST(test_orderedq);
    TEST(test_order);
    TEST(test_ordering);

    printf("All sort tests passed!\n");
    return 0;
}

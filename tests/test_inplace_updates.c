/* Regression tests for the in-place update paths of the data-structure program:
 * AssociateTo / AppendTo / Part assignment mutate a value in place when the
 * variable is its only owner (symtab_value_exclusively_owned), and must still
 * behave exactly as value-semantics demand when anything else holds the value.
 *
 * Every test pairs an in-place-eligible sequence with a check that an aliased
 * copy, a Block-saved value, or a cached derived expression did NOT change. */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }
static void run(const char* src) { Expr* p = parse_expression(src); Expr* r = evaluate(p); expr_free(r); expr_free(p); }

static void test_assoc_inplace_semantics(void) {
    run("a1 = <||>; Do[AssociateTo[a1, i -> i^2], {i, 5}]");
    assert_true("a1 === <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16, 5 -> 25|>");
    run("AssociateTo[a1, 3 -> x]");
    assert_true("a1 === <|1 -> 1, 2 -> 4, 3 -> x, 4 -> 16, 5 -> 25|>");   /* keeps position */
    run("b1 = a1; AssociateTo[a1, 9 -> 1]");
    assert_true("Length[a1] == 6 && Length[b1] == 5 && !KeyExistsQ[b1, 9]");
    run("c1 = <|x -> 1|>; f1 = Length[c1]; AssociateTo[c1, y -> 2]");
    assert_true("Length[c1] == 2 && f1 == 1");                               /* cached read invalidated */
    assert_true("Block[{c1 = c1}, AssociateTo[c1, z -> 3]; Length[c1]] == 3 && Length[c1] == 2");
}
static void test_assoc_appendto_order(void) {
    /* the reference: an existing key moves to the END */
    run("d1 = <|1 -> a, 2 -> b|>; AppendTo[d1, 1 -> z]");
    assert_true("d1 === <|2 -> b, 1 -> z|>");
    run("AppendTo[d1, 3 -> w]");
    assert_true("d1 === <|2 -> b, 1 -> z, 3 -> w|>");
}
static void test_assoc_part_set(void) {
    run("e1 = <|1 -> 1|>; e1[[Key[2]]] = 5; e1[[Key[1]]] = 9; e1[\"s\"] = 3");
    assert_true("e1 === <|1 -> 9, 2 -> 5, \"s\" -> 3|>");
    run("g1 = e1; e1[[Key[7]]] = 0");
    assert_true("Length[e1] == 4 && Length[g1] == 3");
    run("Do[e1[[Key[k]]] = k, {k, 2000}]");
    assert_true("Length[e1] == 2001 && e1[[Key[1999]]] == 1999 && Length[Union[Keys[e1]]] == Length[e1]");
}
static void test_packed_part_set(void) {
    run("p1 = Range[10]; q1 = p1; p1[[3]] = 99");
    assert_true("p1[[3]] == 99 && q1[[3]] == 3");
    run("r1 = Range[5] * 1.5; r1[[-1]] = 2.5");
    assert_true("r1 === {1.5, 3., 4.5, 6., 2.5}");
    run("m1 = {{1, 2}, {3, 4}}; m1[[2, 1]] = 9");
    assert_true("m1 === {{1, 2}, {9, 4}}");
    run("v1 = Range[8]; w1 = v1; Do[v1[[i]] = 0, {i, 8}]");
    assert_true("Total[v1] == 0 && Total[w1] == 36");
    run("t1 = Range[4]");
    assert_true("Block[{t1 = t1}, t1[[1]] = 0; t1] === {0, 2, 3, 4} && t1 === {1, 2, 3, 4}");
    /* a value the buffer cannot hold exactly must still give the List answer */
    run("u1 = Range[4]; u1[[2]] = 2.5");
    assert_true("u1 === {1, 2.5, 3, 4}");
    run("x1 = Range[3]; x1[[2]] = y");
    assert_true("x1 === {1, y, 3}");
}
static void test_list_appendto(void) {
    run("l1 = {}; Do[AppendTo[l1, i], {i, 10}]");
    assert_true("l1 === Range[10]");
    run("m2 = l1; AppendTo[l1, 11]");
    assert_true("Length[l1] == 11 && Length[m2] == 10");
    run("r2 = {1.5, 2.5, 3.5, 4.5}; AppendTo[r2, 5.5]; r3 = r2; AppendTo[r2, 6.5]");
    assert_true("r2 === {1.5, 2.5, 3.5, 4.5, 5.5, 6.5} && r3 === {1.5, 2.5, 3.5, 4.5, 5.5}");
    run("AppendTo[r2, \"s\"]");                                      /* leaves the buffer */
    assert_true("r2 === {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, \"s\"}");
    run("s1 = {\"a\"}; s2 = s1; Do[AppendTo[s1, \"b\"], {3}]");
    assert_true("s1 === {\"a\", \"b\", \"b\", \"b\"} && s2 === {\"a\"}");
    run("w2 = Range[5]; AppendTo[w2, 6]; k2 = Total[w2]; AppendTo[w2, 7]");
    assert_true("k2 == 21 && Total[w2] == 28");                       /* cache not stale */
    assert_true("Block[{w2 = w2}, AppendTo[w2, 100]; Length[w2]] == 8 && Length[w2] == 7");
    run("big = {}; Do[AppendTo[big, i], {i, 20000}]");
    assert_true("big === Range[20000]");
}
static void test_compile_append_order(void) {
    assert_true("Compile[{{p, _Association, _Real}, {x, _Real}}, Append[p, \"b\" -> x]]"
                "[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>, 9.] === Append[<|\"a\" -> 1., \"b\" -> 2., \"c\" -> 3.|>, \"b\" -> 9.]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_assoc_inplace_semantics);
    TEST(test_assoc_appendto_order);
    TEST(test_assoc_part_set);
    TEST(test_packed_part_set);
    TEST(test_list_appendto);
    TEST(test_compile_append_order);
    printf("All in-place update tests passed!\n");
    return 0;
}

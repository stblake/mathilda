#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parse `input`, evaluate, and assert that the printed (non-FullForm)
 * result equals `expected`. */
static void check(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT_MSG(parsed != NULL, "parse failed: %s", input);
    Expr* res = evaluate(parsed);
    char* s = expr_to_string(res);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Thread %s: expected `%s`, got `%s`", input, expected, s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

/* ---------- 1. Basic forms (default head = List) ---------- */

static void test_thread_basic_single_list(void) {
    check("Thread[f[{a, b, c}]]", "{f[a], f[b], f[c]}");
}

static void test_thread_broadcast_scalar(void) {
    check("Thread[f[{a, b, c}, x]]", "{f[a, x], f[b, x], f[c, x]}");
}

static void test_thread_two_lists(void) {
    check("Thread[f[{a, b, c}, {x, y, z}]]", "{f[a, x], f[b, y], f[c, z]}");
}

static void test_thread_scalar_between_lists(void) {
    check("Thread[f[{a, b, c}, h, {x, y, z}]]", "{f[a, h, x], f[b, h, y], f[c, h, z]}");
}

static void test_thread_atom_passthrough(void) {
    /* Atom input: nothing to thread over -- return the atom. */
    check("Thread[5]", "5");
    check("Thread[x]", "x");
}

static void test_thread_no_threadable_args(void) {
    /* f[a, b] has no List argument -- result is f[a, b] unchanged. */
    check("Thread[f[a, b]]", "f[a, b]");
}

static void test_thread_empty_list(void) {
    /* Threading over a single empty list produces an empty wrapper. */
    check("Thread[f[{}]]", "{}");
}

/* ---------- 2. Custom head h ---------- */

static void test_thread_plus_head(void) {
    /* Default does not thread over Plus. */
    check("Thread[f[a + b, c + d]]", "f[a + b, c + d]");
    /* With head=Plus, it does. */
    check("Thread[f[a + b, c + d], Plus]", "f[a, c] + f[b, d]");
}

static void test_thread_equation(void) {
    /* a==b is Equal[a,b]; threading turns it into a list of equalities. */
    check("Thread[{a, b, c} == {x, y, z}]", "{a == x, b == y, c == z}");
}

static void test_thread_apply_to_equation(void) {
    /* Log on an equation, threading over Equal: Log[x]==Log[y]. */
    check("Thread[Log[x == y], Equal]", "Log[x] == Log[y]");
}

static void test_thread_custom_head_explicit_all(void) {
    /* Explicit n=All should match the default. */
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List]",
          "{f[a, r, u, x], f[b, s, v, y]}");
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, All]",
          "{f[a, r, u, x], f[b, s, v, y]}");
}

static void test_thread_unmatched_head_passthrough(void) {
    /* Head h=Times nowhere in args -- nothing threadable, return as-is. */
    check("Thread[f[a + b, c + d], Times]", "f[a + b, c + d]");
}

/* ---------- 3. Position specs ---------- */

static void test_thread_position_none(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, None]",
          "f[{a, b}, {r, s}, {u, v}, {x, y}]");
}

static void test_thread_position_first_n(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, 2]",
          "{f[a, r, {u, v}, {x, y}], f[b, s, {u, v}, {x, y}]}");
}

static void test_thread_position_last_n(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, -2]",
          "{f[{a, b}, {r, s}, u, x], f[{a, b}, {r, s}, v, y]}");
}

static void test_thread_position_single_index(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {2}]",
          "{f[{a, b}, r, {u, v}, {x, y}], f[{a, b}, s, {u, v}, {x, y}]}");
}

static void test_thread_position_range(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {2, 4}]",
          "{f[{a, b}, r, u, x], f[{a, b}, s, v, y]}");
}

static void test_thread_position_step(void) {
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {1, -1, 2}]",
          "{f[a, {r, s}, u, {x, y}], f[b, {r, s}, v, {x, y}]}");
}

static void test_thread_position_negative_step(void) {
    /* {4,1,-1} should select all four positions (4,3,2,1). */
    check("Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {4, 1, -1}]",
          "{f[a, r, u, x], f[b, s, v, y]}");
}

static void test_thread_position_single_negative(void) {
    /* {-1} = only the last argument. */
    check("Thread[f[a, b, {x, y}], List, {-1}]",
          "{f[a, b, x], f[a, b, y]}");
}

static void test_thread_position_first_n_clamped(void) {
    /* n larger than arg count is clamped to all args. */
    check("Thread[f[{a, b}, x], List, 5]",
          "{f[a, x], f[b, x]}");
}

/* ---------- 4. Mismatched length ---------- */

static void test_thread_mismatch_lengths(void) {
    /* Lengths 3 and 2 of List args -- Mathilda returns the expression
     * unchanged (we don't issue the message, but the expression must
     * remain well-formed and not crash). */
    check("Thread[f[{a, b, c}, {x, y}]]", "f[{a, b, c}, {x, y}]");
}

/* ---------- 5. Threading on built-ins ---------- */

static void test_thread_evaluates_result(void) {
    /* Plus is Listable+Orderless; once threaded, each Plus[a, x] is
     * still just a+x, so Thread[Plus[{1,2}, {10,20}]] computes. */
    check("Thread[Plus[{1, 2}, {10, 20}]]", "{11, 22}");
}

static void test_thread_with_orderless_head(void) {
    /* Threading respects Listable on the inner function:
     * Thread[Sin[{0, x}]] -> {Sin[0], Sin[x]} -> {0, Sin[x]}. */
    check("Thread[Sin[{0, x}]]", "{0, Sin[x]}");
}

/* ---------- 6. Unevaluated / shape ---------- */

static void test_thread_returns_unevaluated_on_bad_spec(void) {
    /* Bad spec (string, non-integer) -> Thread should remain
     * unevaluated. We verify the call survives without crashing. */
    Expr* parsed = parse_expression("Thread[f[{a, b}], List, \"bad\"]");
    Expr* res = evaluate(parsed);
    ASSERT(res != NULL);
    expr_free(parsed);
    expr_free(res);
}

static void test_thread_attribute_protected(void) {
    /* Thread is Protected: assignment to Thread should fail without crash. */
    Expr* parsed = parse_expression("Attributes[Thread]");
    Expr* res = evaluate(parsed);
    char* s = expr_to_string(res);
    /* Must contain "Protected" somewhere in the attribute list. */
    ASSERT_MSG(strstr(s, "Protected") != NULL,
               "Thread should have Protected attribute, got %s", s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

/* ---------- 7. Stress: repeated threading does not leak ---------- */

static void test_thread_repeated_calls(void) {
    /* Run the same threading 200 times to give valgrind a clear leak
     * signal if we forget to free anything. */
    for (int i = 0; i < 200; i++) {
        Expr* parsed = parse_expression("Thread[f[{a, b, c}, {x, y, z}]]");
        Expr* res = evaluate(parsed);
        expr_free(parsed);
        expr_free(res);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_thread_basic_single_list);
    TEST(test_thread_broadcast_scalar);
    TEST(test_thread_two_lists);
    TEST(test_thread_scalar_between_lists);
    TEST(test_thread_atom_passthrough);
    TEST(test_thread_no_threadable_args);
    TEST(test_thread_empty_list);

    TEST(test_thread_plus_head);
    TEST(test_thread_equation);
    TEST(test_thread_apply_to_equation);
    TEST(test_thread_custom_head_explicit_all);
    TEST(test_thread_unmatched_head_passthrough);

    TEST(test_thread_position_none);
    TEST(test_thread_position_first_n);
    TEST(test_thread_position_last_n);
    TEST(test_thread_position_single_index);
    TEST(test_thread_position_range);
    TEST(test_thread_position_step);
    TEST(test_thread_position_negative_step);
    TEST(test_thread_position_single_negative);
    TEST(test_thread_position_first_n_clamped);

    TEST(test_thread_mismatch_lengths);
    TEST(test_thread_evaluates_result);
    TEST(test_thread_with_orderless_head);

    TEST(test_thread_returns_unevaluated_on_bad_spec);
    TEST(test_thread_attribute_protected);
    TEST(test_thread_repeated_calls);

    printf("All Thread tests passed!\n");
    return 0;
}

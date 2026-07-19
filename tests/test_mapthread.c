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
               "MapThread %s: expected `%s`, got `%s`", input, expected, s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

/* ---------- 1. Basic level-1 threading ---------- */

static void test_mapthread_two_lists(void) {
    check("MapThread[f, {{a, b, c}, {x, y, z}}]",
          "{f[a, x], f[b, y], f[c, z]}");
}

static void test_mapthread_single_list(void) {
    /* One argument list -- degenerates to Map. */
    check("MapThread[f, {{a, b, c}}]", "{f[a], f[b], f[c]}");
}

static void test_mapthread_four_lists(void) {
    check("MapThread[f, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]",
          "{f[1, 4, 7, 10], f[2, 5, 8, 11], f[3, 6, 9, 12]}");
}

static void test_mapthread_length_one(void) {
    check("MapThread[f, {{a}, {x}}]", "{f[a, x]}");
}

/* ---------- 2. Level specification ---------- */

static void test_mapthread_level_2_matrix(void) {
    check("MapThread[f, {{{a, b}, {c, d}}, {{u, v}, {s, t}}}, 2]",
          "{{f[a, u], f[b, v]}, {f[c, s], f[d, t]}}");
}

static void test_mapthread_level_1_default_nested(void) {
    /* Default level 1 threads over the outer list only; the sublists are
     * passed whole to f. */
    check("MapThread[f, {{{a, b}, {c, d}}, {{u, v}, {w, x}}}]",
          "{f[{a, b}, {u, v}], f[{c, d}, {w, x}]}");
}

static void test_mapthread_level_0(void) {
    /* Level 0 applies f to the whole expressions with no descent. */
    check("MapThread[f, {a, b}, 0]", "f[a, b]");
    check("MapThread[f, {{a, b}, {c, d}}, 0]", "f[{a, b}, {c, d}]");
}

static void test_mapthread_level_2_three_lists(void) {
    check("MapThread[f, {{{1, 2}}, {{3, 4}}, {{5, 6}}}, 2]",
          "{{f[1, 3, 5], f[2, 4, 6]}}");
}

/* ---------- 3. Attribute-driven evaluation of f ---------- */

static void test_mapthread_plus(void) {
    /* Plus adds columns of the matrix (like Total). */
    check("MapThread[Plus, {{a, b, c}, {u, v, w}, {x, y, z}}]",
          "{a + u + x, b + v + y, c + w + z}");
}

static void test_mapthread_list_head(void) {
    /* List gathers columns (like Transpose). */
    check("MapThread[List, {{a, b, c}, {u, v, w}, {x, y, z}}]",
          "{{a, u, x}, {b, v, y}, {c, w, z}}");
}

static void test_mapthread_derivative(void) {
    /* MapThread[D, ...] differentiates each expr wrt its variable. */
    check("MapThread[D, {{x, x y, x z}, {x, y, z}}]", "{1, x, x}");
}

static void test_mapthread_numeric_plus(void) {
    check("MapThread[Plus, {{1, 2}, {10, 20}}]", "{11, 22}");
}

/* ---------- 4. Equivalence with Thread ---------- */

static void test_mapthread_matches_thread(void) {
    /* MapThread[f, args] == Thread[f[args...]] for corresponding lists. */
    check("MapThread[f, {{a, b, c}, {x, y, z}}]",
          "{f[a, x], f[b, y], f[c, z]}");
    check("Thread[f[{a, b, c}, {x, y, z}]]",
          "{f[a, x], f[b, y], f[c, z]}");
}

/* ---------- 5. Associations ---------- */

static void test_mapthread_assoc_single_key(void) {
    check("MapThread[f, {<|a -> 1|>, <|a -> 2|>}]", "<|a -> f[1, 2]|>");
}

static void test_mapthread_assoc_multi_key(void) {
    check("MapThread[f, {<|a -> 1, b -> 2|>, <|a -> 3, b -> 4|>}]",
          "<|a -> f[1, 3], b -> f[2, 4]|>");
}

static void test_mapthread_assoc_plus(void) {
    check("MapThread[Plus, {<|a -> 1, b -> 2|>, <|a -> 10, b -> 20|>}]",
          "<|a -> 11, b -> 22|>");
}

static void test_mapthread_assoc_mismatched_keys_unevaluated(void) {
    /* Different key sets -> structural mismatch -> left unevaluated. */
    Expr* parsed = parse_expression("MapThread[f, {<|a -> 1|>, <|b -> 2|>}]");
    Expr* res = evaluate(parsed);
    char* s = expr_to_string(res);
    ASSERT_MSG(strncmp(s, "MapThread[", 10) == 0,
               "mismatched assoc keys should stay unevaluated, got `%s`", s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

/* ---------- 6. Edge cases / unevaluated forms ---------- */

static void test_mapthread_empty(void) {
    check("MapThread[f, {}]", "{}");
}

static void test_mapthread_empty_sublists(void) {
    check("MapThread[f, {{}, {}}]", "{}");
}

static void test_mapthread_unequal_lengths(void) {
    /* Unequal list lengths -> left unchanged (no crash, no message text). */
    check("MapThread[f, {{a, b, c}, {u, v}}]",
          "MapThread[f, {{a, b, c}, {u, v}}]");
}

static void test_mapthread_nonlist_element(void) {
    /* A non-list among the lists -> left unchanged. */
    check("MapThread[f, {{a, b, c}, t, {u, v, w}}]",
          "MapThread[f, {{a, b, c}, t, {u, v, w}}]");
}

static void test_mapthread_outer_not_list(void) {
    /* Second argument is not a List -> left unchanged. */
    check("MapThread[f, g[{a, b}, {c, d}]]",
          "MapThread[f, g[{a, b}, {c, d}]]");
}

static void test_mapthread_negative_level_unevaluated(void) {
    /* Negative level is invalid -> left unchanged. */
    check("MapThread[f, {{a, b}, {c, d}}, -1]",
          "MapThread[f, {{a, b}, {c, d}}, -1]");
}

static void test_mapthread_noninteger_level_unevaluated(void) {
    check("MapThread[f, {{a, b}, {c, d}}, x]",
          "MapThread[f, {{a, b}, {c, d}}, x]");
}

/* ---------- 7. Arity / attribute ---------- */

static void test_mapthread_arity_zero(void) {
    /* MapThread[] -> unchanged (prints MapThread::argt on stderr). */
    Expr* parsed = parse_expression("MapThread[]");
    Expr* res = evaluate(parsed);
    char* s = expr_to_string(res);
    ASSERT_MSG(strcmp(s, "MapThread[]") == 0,
               "MapThread[] should stay unchanged, got `%s`", s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

static void test_mapthread_arity_one(void) {
    check("MapThread[f]", "MapThread[f]");
}

static void test_mapthread_arity_four(void) {
    check("MapThread[f, {{a}}, 1, extra]", "MapThread[f, {{a}}, 1, extra]");
}

static void test_mapthread_protected(void) {
    Expr* parsed = parse_expression("Attributes[MapThread]");
    Expr* res = evaluate(parsed);
    char* s = expr_to_string(res);
    ASSERT_MSG(strstr(s, "Protected") != NULL,
               "MapThread should have Protected attribute, got %s", s);
    free(s);
    expr_free(parsed);
    expr_free(res);
}

/* ---------- 8. Function-in-list, pure functions ---------- */

static void test_mapthread_pure_function(void) {
    /* Apply the functions in a list to corresponding arguments. */
    check("MapThread[#1[#2] &, {{f, g, h}, {x, y, z}}]",
          "{f[x], g[y], h[z]}");
}

/* ---------- 9. Stress: repeated calls (valgrind leak signal) ---------- */

static void test_mapthread_repeated_success(void) {
    for (int i = 0; i < 200; i++) {
        Expr* parsed = parse_expression(
            "MapThread[f, {{a, b, c}, {x, y, z}, {p, q, r}}]");
        Expr* res = evaluate(parsed);
        expr_free(parsed);
        expr_free(res);
    }
}

static void test_mapthread_repeated_failure(void) {
    /* Exercise the NULL / cleanup paths repeatedly (mismatch + non-list). */
    for (int i = 0; i < 200; i++) {
        Expr* p1 = parse_expression("MapThread[f, {{a, b, c}, {u, v}}]");
        Expr* r1 = evaluate(p1);
        expr_free(p1); expr_free(r1);

        Expr* p2 = parse_expression("MapThread[f, {{a, b}, t}, 2]");
        Expr* r2 = evaluate(p2);
        expr_free(p2); expr_free(r2);

        Expr* p3 = parse_expression("MapThread[f, {<|a -> 1|>, <|b -> 2|>}]");
        Expr* r3 = evaluate(p3);
        expr_free(p3); expr_free(r3);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_mapthread_two_lists);
    TEST(test_mapthread_single_list);
    TEST(test_mapthread_four_lists);
    TEST(test_mapthread_length_one);

    TEST(test_mapthread_level_2_matrix);
    TEST(test_mapthread_level_1_default_nested);
    TEST(test_mapthread_level_0);
    TEST(test_mapthread_level_2_three_lists);

    TEST(test_mapthread_plus);
    TEST(test_mapthread_list_head);
    TEST(test_mapthread_derivative);
    TEST(test_mapthread_numeric_plus);

    TEST(test_mapthread_matches_thread);

    TEST(test_mapthread_assoc_single_key);
    TEST(test_mapthread_assoc_multi_key);
    TEST(test_mapthread_assoc_plus);
    TEST(test_mapthread_assoc_mismatched_keys_unevaluated);

    TEST(test_mapthread_empty);
    TEST(test_mapthread_empty_sublists);
    TEST(test_mapthread_unequal_lengths);
    TEST(test_mapthread_nonlist_element);
    TEST(test_mapthread_outer_not_list);
    TEST(test_mapthread_negative_level_unevaluated);
    TEST(test_mapthread_noninteger_level_unevaluated);

    TEST(test_mapthread_arity_zero);
    TEST(test_mapthread_arity_one);
    TEST(test_mapthread_arity_four);
    TEST(test_mapthread_protected);

    TEST(test_mapthread_pure_function);

    TEST(test_mapthread_repeated_success);
    TEST(test_mapthread_repeated_failure);

    printf("All MapThread tests passed!\n");
    return 0;
}

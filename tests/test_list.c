#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void test_min() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Min[9, 2]", "2"},
        {"Min[{4, 1, 7, 2}]", "1"},
        {"Min[{{-1, 0, 1, 2}, {0, 2, 4, 6}, {-3, -2, -1, 0}}]", "-3"},
        {"Min[Infinity, 5]", "5"},
        {"Min[-1 * Infinity, 5]", "-Infinity"},
        {"Min[{a, b}, {c, d}]", "Min[a, b, c, d]"},
        {"Min[]", "Infinity"},
        {"Min[x, 3, 5]", "Min[3, x]"},
        {"Min[x, x]", "x"},
        /* BigInt and MPFR must compare alongside Integer / Real. */
        {"Min[10^50, 3]", "3"},
        {"Min[10^50, 10^60, 10^40, 5]", "5"},
        {"Min[10^50, 10^60, 10^40]", "10000000000000000000000000000000000000000"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Min test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_max() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Max[9, 2]", "9"},
        {"Max[{4, 1, 7, 2}]", "7"},
        {"Max[Infinity, 5]", "Infinity"},
        {"Max[-1 * Infinity, 5]", "5"},
        {"Max[]", "-Infinity"},
        {"Max[x, 3, 5]", "Max[5, x]"},
        {"Max[x, x]", "x"},
        /* BigInt and MPFR must compare alongside Integer / Real. */
        {"Max[10^50, 1.5]", "100000000000000000000000000000000000000000000000000"},
        {"Max[10^50, 10^60, 10^40, 5]", "1000000000000000000000000000000000000000000000000000000000000"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Max test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_listq() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"ListQ[{a, b, c}]", "True"},
        {"ListQ[a]", "False"},
        {"ListQ[f[a]]", "False"},
        {"ListQ[{}]", "True"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("ListQ test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_vectorq() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"VectorQ[{a, b, c}]", "True"},
        {"VectorQ[{{1}, {2}}]", "False"},
        {"VectorQ[{{1}, {2, 3}}, ListQ]", "True"},
        {"VectorQ[{a, 1.2}, NumberQ]", "False"},
        {"VectorQ[Range[10], IntegerQ]", "True"},
        {"VectorQ[a]", "False"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("VectorQ test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_matrixq() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"MatrixQ[{{a, b}, {3, 4}}]", "True"},
        {"MatrixQ[{{1}, {2, 3}}]", "False"},
        {"MatrixQ[Array[a, {2, 2, 2}]]", "False"},
        {"MatrixQ[Array[a, {2, 2, 2}], ListQ]", "True"},
        {"MatrixQ[{}]", "False"},
        {"MatrixQ[{{}}]", "True"},
        {"MatrixQ[{{1, 2}, {3, 4}}, NumberQ]", "True"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("MatrixQ test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_table_n() {
    Expr* t = parse_expression("Table[x, 3]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_SYMBOL);
        ASSERT_STR_EQ(res->data.function.args[i]->data.symbol.name, "x");
    }
    expr_free(t); expr_free(res);
}

void test_table_imax() {
    Expr* t = parse_expression("Table[i, {i, 3}]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(res->data.function.args[i]->data.integer == i + 1);
    }
    expr_free(t); expr_free(res);
}

void test_table_imin_imax() {
    Expr* t = parse_expression("Table[i, {i, 2, 4}]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(res->data.function.args[i]->data.integer == i + 2);
    }
    expr_free(t); expr_free(res);
}

void test_table_imin_imax_di() {
    Expr* t = parse_expression("Table[i, {i, 1, 5, 2}]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    ASSERT(res->data.function.args[0]->data.integer == 1);
    ASSERT(res->data.function.args[1]->data.integer == 3);
    ASSERT(res->data.function.args[2]->data.integer == 5);
    expr_free(t); expr_free(res);
}

void test_table_list() {
    Expr* t = parse_expression("Table[i, {i, {a, b, c}}]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    ASSERT_STR_EQ(res->data.function.args[0]->data.symbol.name, "a");
    ASSERT_STR_EQ(res->data.function.args[1]->data.symbol.name, "b");
    ASSERT_STR_EQ(res->data.function.args[2]->data.symbol.name, "c");
    expr_free(t); expr_free(res);
}

void test_table_nested() {
    Expr* t = parse_expression("Table[i + j, {i, 2}, {j, 3}]");
    /* 6 elements packs at PACK_MIN_ELEMENTS = 4; the value is unchanged and a
     * packed list is a List, but this walks args[] directly. */
    Expr* res = test_delist(evaluate(t));
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 2);
    
    // row1 = {2, 3, 4}
    Expr* row1 = res->data.function.args[0];
    ASSERT(row1->data.function.arg_count == 3);
    ASSERT(row1->data.function.args[0]->data.integer == 2);
    ASSERT(row1->data.function.args[1]->data.integer == 3);
    ASSERT(row1->data.function.args[2]->data.integer == 4);
    
    // row2 = {3, 4, 5}
    Expr* row2 = res->data.function.args[1];
    ASSERT(row2->data.function.arg_count == 3);
    ASSERT(row2->data.function.args[0]->data.integer == 3);
    ASSERT(row2->data.function.args[1]->data.integer == 4);
    ASSERT(row2->data.function.args[2]->data.integer == 5);
    
    expr_free(t); expr_free(res);
}

void test_range_imax() {
    Expr* t = parse_expression("Range[3]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(res->data.function.args[i]->data.integer == i + 1);
    }
    expr_free(t); expr_free(res);
}

void test_range_imin_imax() {
    Expr* t = parse_expression("Range[2, 4]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(res->data.function.args[i]->data.integer == i + 2);
    }
    expr_free(t); expr_free(res);
}

void test_range_imin_imax_di() {
    Expr* t = parse_expression("Range[1, 5, 2]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    ASSERT(res->data.function.args[0]->data.integer == 1);
    ASSERT(res->data.function.args[1]->data.integer == 3);
    ASSERT(res->data.function.args[2]->data.integer == 5);
    expr_free(t); expr_free(res);
}

void test_range_real() {
    Expr* t = parse_expression("Range[1.5, 3.5]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    ASSERT(res->data.function.args[0]->type == EXPR_REAL);
    ASSERT(res->data.function.args[0]->data.real == 1.5);
    ASSERT(res->data.function.args[1]->data.real == 2.5);
    ASSERT(res->data.function.args[2]->data.real == 3.5);
    expr_free(t); expr_free(res);
}

/* Range used to TRUNCATE SILENTLY at 10^6 elements: `Range[2000000]` answered
 * with 1000001 of them and said nothing, so every downstream Length, Total and
 * Plot was quietly computed on half the data.  A resource limit is legitimate;
 * one that changes the answer instead of declining is not, and one that does so
 * without a message cannot even be noticed.  Found on 2026-08-02 while sizing a
 * benchmark vector -- by accident, downstream of the damage, which is how a
 * silent truncation always gets found.
 *
 * The ceiling is now on BYTES (2 GiB, ~268 million int64 elements), and
 * reaching it leaves the call UNEVALUATED with a Range::toobig message. */
void test_range_no_silent_truncation() {
    /* The exact length, past the old cap, on the packed and the plain path. */
    assert_eval_eq("Length[Range[2000000]]", "2000000", 0);
    assert_eval_eq("Last[Range[2000000]]", "2000000", 0);
    assert_eval_eq("Length[Range[5000000]]", "5000000", 0);
    assert_eval_eq("Total[Range[2000000]] == 2000000*2000001/2", "True", 0);
    /* And the two-argument and stepped forms, which count the same way. */
    assert_eval_eq("Length[Range[1000000, 3000000]]", "2000001", 0);
    assert_eval_eq("Length[Range[1, 4000000, 2]]", "2000000", 0);
    /* An inexact range past the old cap. */
    assert_eval_eq("Length[Range[0., 2000000.]]", "2000001", 0);

    /* Past the resource ceiling the call DECLINES -- head Range, not a
     * plausible short List.  This is the whole point: a caller can test for it. */
    assert_eval_eq("Head[Range[10^12]]", "Range", 0);
    assert_eval_eq("Head[Range[0., 1., 10^-12]]", "Range", 0);

    /* Everything small is untouched, including the exact branch past 2^53 that
     * the count-up-front rewrite fixed. */
    assert_eval_eq("Range[10]", "{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}", 0);
    assert_eval_eq("Range[2, 10, 3]", "{2, 5, 8}", 0);
    assert_eval_eq("Range[10, 1]", "{}", 0);
    assert_eval_eq("Range[0., 1., 0.25]", "{0.0, 0.25, 0.5, 0.75, 1.0}", 0);
    assert_eval_eq("Range[10^18, 10^18 + 3]",
                   "{1000000000000000000, 1000000000000000001, "
                   "1000000000000000002, 1000000000000000003}", 0);
}

void test_array_n() {
    Expr* t = parse_expression("Array[f, 3]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_FUNCTION);
        ASSERT_STR_EQ(res->data.function.args[i]->data.function.head->data.symbol.name, "f");
        ASSERT(res->data.function.args[i]->data.function.arg_count == 1);
        ASSERT(res->data.function.args[i]->data.function.args[0]->data.integer == i + 1);
    }
    expr_free(t); expr_free(res);
}

void test_array_n_r() {
    Expr* t = parse_expression("Array[f, 3, 0]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 3);
    for (int i = 0; i < 3; i++) {
        ASSERT(res->data.function.args[i]->type == EXPR_FUNCTION);
        ASSERT_STR_EQ(res->data.function.args[i]->data.function.head->data.symbol.name, "f");
        ASSERT(res->data.function.args[i]->data.function.arg_count == 1);
        ASSERT(res->data.function.args[i]->data.function.args[0]->data.integer == i);
    }
    expr_free(t); expr_free(res);
}

void test_array_nested() {
    Expr* t = parse_expression("Array[f, {2, 3}]");
    Expr* res = evaluate(t);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "List");
    ASSERT(res->data.function.arg_count == 2);
    
    Expr* row1 = res->data.function.args[0];
    ASSERT(row1->data.function.arg_count == 3);
    ASSERT(row1->data.function.args[0]->data.function.args[0]->data.integer == 1);
    ASSERT(row1->data.function.args[0]->data.function.args[1]->data.integer == 1);
    
    Expr* row2 = res->data.function.args[1];
    ASSERT(row2->data.function.arg_count == 3);
    ASSERT(row2->data.function.args[2]->data.function.args[0]->data.integer == 2);
    ASSERT(row2->data.function.args[2]->data.function.args[1]->data.integer == 3);
    
    expr_free(t); expr_free(res);
}

void test_take() {
    Expr* t1 = parse_expression("Take[{a, b, c, d}, 2]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->type == EXPR_FUNCTION);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT_STR_EQ(res1->data.function.args[0]->data.symbol.name, "a");
    ASSERT_STR_EQ(res1->data.function.args[1]->data.symbol.name, "b");
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("Take[{a, b, c, d}, -2]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->type == EXPR_FUNCTION);
    ASSERT(res2->data.function.arg_count == 2);
    ASSERT_STR_EQ(res2->data.function.args[0]->data.symbol.name, "c");
    ASSERT_STR_EQ(res2->data.function.args[1]->data.symbol.name, "d");
    expr_free(t2); expr_free(res2);

    Expr* t3 = parse_expression("Take[{a, b, c, d}, {2, 3}]");
    Expr* res3 = evaluate(t3);
    ASSERT(res3->type == EXPR_FUNCTION);
    ASSERT(res3->data.function.arg_count == 2);
    ASSERT_STR_EQ(res3->data.function.args[0]->data.symbol.name, "b");
    ASSERT_STR_EQ(res3->data.function.args[1]->data.symbol.name, "c");
    expr_free(t3); expr_free(res3);
}

void test_drop() {
    Expr* t1 = parse_expression("Drop[{a, b, c, d}, 2]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->type == EXPR_FUNCTION);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT_STR_EQ(res1->data.function.args[0]->data.symbol.name, "c");
    ASSERT_STR_EQ(res1->data.function.args[1]->data.symbol.name, "d");
    expr_free(t1); expr_free(res1);
}

void test_flatten() {
    Expr* t1 = parse_expression("Flatten[{{a, b}, {c, {d, e}}}]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->data.function.arg_count == 5);
    ASSERT_STR_EQ(res1->data.function.args[3]->data.symbol.name, "d");
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("Flatten[{{a, b}, {c, {d, e}}}, 1]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->data.function.arg_count == 4);
    ASSERT_STR_EQ(res2->data.function.args[3]->data.function.head->data.symbol.name, "List");
    expr_free(t2); expr_free(res2);
}

void test_partition() {
    Expr* t1 = parse_expression("Partition[{a, b, c, d, e}, 2]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->data.function.arg_count == 2);
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("Partition[{a, b, c, d, e}, 2, 1]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->data.function.arg_count == 4);
    expr_free(t2); expr_free(res2);
}

void test_pick() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* Two-argument form: literal True selects. */
        {"Pick[{a, b, c, d}, {True, False, True, False}]", "{a, c}"},
        {"Pick[{a, b}, {True, True}]", "{a, b}"},
        {"Pick[{a, b}, {False, False}]", "{}"},
        {"Pick[{}, {}]", "{}"},
        /* Only literal True selects; other selector values are dropped. */
        {"Pick[{a, b, c}, {1, True, xyz}]", "{b}"},

        /* Three-argument form: the selector must match the pattern. */
        {"Pick[{a, b, c, d}, {1, 2, 3, 4}, 3]", "{c}"},
        {"Pick[{a, b, c, d}, {1, 2, 3, 4}, _Integer]", "{a, b, c, d}"},
        {"Pick[{a, b, c, d}, {1, 2, 3, 4}, _String]", "{}"},

        /* All levels; sel mirrors the structure of expr. */
        {"Pick[{{a, b}, {c, d}}, {{1, 0}, {0, 1}}, 1]", "{{a}, {d}}"},
        {"Pick[{{a, b}, {c, d}}, {{True, False}, {False, True}}]", "{{a}, {d}}"},
        /* A matching selector at an outer level keeps that element whole. */
        {"Pick[{{a, b}, {c, d}}, {True, False}]", "{{a, b}}"},

        /* The head comes from expr, not from sel. */
        {"Pick[f[a, b, c], {True, False, True}]", "f[a, c]"},

        /* Structure mismatch -> unevaluated, never a partial result. */
        {"Pick[{a, b, c}, {True, False}]", "Pick[{a, b, c}, {True, False}]"},
        {"Pick[{a, b}, {{1, 0}, {0, 1}}, 1]", "Pick[{a, b}, {{1, 0}, {0, 1}}, 1]"},
        {"Pick[{{a, b}, {c}}, {{1, 0}, {0, 1}}, 1]",
         "Pick[{{a, b}, {c}}, {{1, 0}, {0, 1}}, 1]"},
        {"Pick[a, {True}]", "Pick[a, {True}]"},
        {"Pick[{a, b}, sel]", "Pick[{a, b}, sel]"},

        {"Attributes[Pick]", "{Protected}"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Pick test failed: %s expected %s, got %s\n",
                   tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_rotate() {
    Expr* t1 = parse_expression("RotateLeft[{a, b, c}, 1]");
    Expr* res1 = evaluate(t1);
    ASSERT_STR_EQ(res1->data.function.args[0]->data.symbol.name, "b");
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("RotateRight[{a, b, c}, 1]");
    Expr* res2 = evaluate(t2);
    ASSERT_STR_EQ(res2->data.function.args[0]->data.symbol.name, "c");
    expr_free(t2); expr_free(res2);
}

void test_reverse() {
    Expr* t1 = parse_expression("Reverse[{a, b, c}]");
    Expr* res1 = evaluate(t1);
    ASSERT_STR_EQ(res1->data.function.args[0]->data.symbol.name, "c");
    expr_free(t1); expr_free(res1);
}

void test_transpose() {
    Expr* t1 = parse_expression("Transpose[{{a, b}, {c, d}}]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT_STR_EQ(res1->data.function.args[0]->data.function.args[1]->data.symbol.name, "c");
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("Transpose[{{a, b}, {c, d}}, {1, 1}]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->data.function.arg_count == 2);
    ASSERT_STR_EQ(res2->data.function.args[0]->data.symbol.name, "a");
    ASSERT_STR_EQ(res2->data.function.args[1]->data.symbol.name, "d");
    expr_free(t2); expr_free(res2);
}

void test_total() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Total[{a, b, c, d}]", "a + b + c + d"},
        {"Total[{1, 2, 3}]", "6"},
        {"Total[{{1, 2}, {3, 4}}]", "{4, 6}"},
        {"Total[{{1, 2}, {3, 4}}, 2]", "10"},
        {"Total[{{1, 2}, {3, 4}}, {2}]", "{3, 7}"},
        {"Total[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, 2]", "{16, 20}"},
        {"Total[{{1, 2}, {3}}, 2]", "6"},
        {"Total[{{1, 2}, {3}}, {1, 2}]", "6"},
        {"Total[{1, 2, 3}, {1, 2}]", "6"},
        {"Total[{{1, 2}, {3, 4}}, {-1}]", "{3, 7}"},
        {"Total[{{1, 2}, {3, 4}}, Infinity]", "10"},
        {"Total[1]", "1"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Total test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_commonest() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Commonest[{b, a, c, 2, a, b, 1, 2}]", "{b, a, 2}"},
        {"Commonest[{b, a, c, 2, a, b, 1, 2}, 4]", "{b, a, c, 2}"},
        {"Commonest[{b, a, c, 2, a, b, 1, 2}, UpTo[6]]", "{b, a, c, 2, 1}"},
        {"Commonest[{1, 2, 2, 3, 3, 3, 4}]", "{3}"},
        {"Commonest[{a, E, Sin[y], E, a, 7}]", "{a, E}"},
        {"Commonest[{1., 2., 2., 3., 3., 3., 4.}]", "{3.0}"},
        {"Commonest[{a, E, Sin[y], E, a, 1.5, 3}, 10]", "{a, E, Sin[y], 1.5, 3}"},
        {"Commonest[{}]", "{}"},
        {"Commonest[{}, 2]", "{}"},
        /* A count of zero or fewer selects nothing. -1 is the case worth
         * pinning: it used to be the internal "no count given" sentinel, so it
         * answered {1} where every other negative count answered {}. */
        {"Commonest[{1, 1, 2}, 0]", "{}"},
        {"Commonest[{1, 1, 2}, -1]", "{}"},
        {"Commonest[{1, 1, 2}, -2]", "{}"},
        /* The buffer path (src/ndreduce.c). Same answers, and the packed List
         * prints as a List -- see test_packed_list.c for the differential. */
        {"Commonest[ToNDArray[{5, 1, 5, 3, 1, 9}]]", "{5, 1}"},
        {"Commonest[ToNDArray[{5, 1, 5, 3, 1, 9}], 2]", "{5, 1}"},
        {"Commonest[ToNDArray[{1., 2., 2., 3., 3., 3., 4.}]]", "{3.0}"},
        /* A VISIBLE NDArray keeps its presentation, as every ndred_* path does. */
        {"Commonest[NDArray[{1, 2, 2}, DataType -> \"int64\"]]", "NDArray[{2}]"},
        /* Rank 2 has no machine-word key -- the rows are the elements -- so the
         * call is handed back and the List path answers. */
        {"Commonest[ToNDArray[{{1, 2}, {1, 2}, {3, 4}}]]", "{{1, 2}}"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Commonest test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_splitby() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* Core semantics: runs of consecutive elements with equal f[e]. */
        {"SplitBy[{1, 3, 2, 4, 5}, EvenQ]", "{{1, 3}, {2, 4}, {5}}"},
        {"SplitBy[{1, 2, 3, 4, 5, 6}, EvenQ]", "{{1}, {2}, {3}, {4}, {5}, {6}}"},
        {"SplitBy[{1, 1, 2, 2, 3}, Identity]", "{{1, 1}, {2, 2}, {3}}"},

        /* Only adjacent elements group -- this is not GatherBy. */
        {"SplitBy[{2, 1, 4}, EvenQ]", "{{2}, {1}, {4}}"},

        /* Empty list. */
        {"SplitBy[{}, EvenQ]", "{}"},
        {"SplitBy[{}, Identity]", "{}"},

        /* Single element. */
        {"SplitBy[{a}, Identity]", "{{a}}"},
        {"SplitBy[{7}, EvenQ]", "{{7}}"},

        /* All elements share a key: one run holding everything. */
        {"SplitBy[{2, 4, 6, 8}, EvenQ]", "{{2, 4, 6, 8}}"},
        {"SplitBy[{5, 5, 5}, Identity]", "{{5, 5, 5}}"},

        /* All keys distinct: one singleton run per element. */
        {"SplitBy[{1, 2, 3}, Identity]", "{{1}, {2}, {3}}"},

        /* An f that stays unevaluated still groups adjacent elements whose
         * keys are structurally identical (keyfn[x] == keyfn[x]). */
        {"SplitBy[{x, x, y}, keyfn]", "{{x, x}, {y}}"},
        {"SplitBy[{x, y, x}, keyfn]", "{{x}, {y}, {x}}"},

        /* Single-function list form matches the scalar form. */
        {"SplitBy[{1, 3, 2, 4}, {EvenQ}]", "{{1, 3}, {2, 4}}"},

        /* Multiple functions nest one level deeper per function. */
        {"SplitBy[{1, 1, 3, 2, 4, 4}, {EvenQ, Identity}]",
         "{{{1, 1}, {3}}, {{2}, {4, 4}}}"},
        {"SplitBy[{1, 3, 2, 4}, {EvenQ, Identity}]", "{{{1}, {3}}, {{2}, {4}}}"},

        /* Empty function list: nothing to split by, left unevaluated. */
        {"SplitBy[{1, 2, 3}, {}]", "SplitBy[{1, 2, 3}, {}]"},

        /* Atoms have no elements to split; left unevaluated. */
        {"SplitBy[x, EvenQ]", "SplitBy[x, EvenQ]"},

        /* Wrong arity is left unevaluated. */
        {"SplitBy[{1, 2}]", "SplitBy[{1, 2}]"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("SplitBy test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_join_basic() {
    /* Basic concatenation of lists */
    assert_eval_eq("Join[{a, b, c}, {x, y}, {u, v, w}]",
                   "{a, b, c, x, y, u, v, w}", 0);
}

void test_join_two_lists() {
    assert_eval_eq("Join[{1, 2}, {3, 4}]", "{1, 2, 3, 4}", 0);
}

void test_join_single_list() {
    assert_eval_eq("Join[{a, b}]", "{a, b}", 0);
}

void test_join_empty_lists() {
    assert_eval_eq("Join[{}, {a, b}]", "{a, b}", 0);
    assert_eval_eq("Join[{a, b}, {}]", "{a, b}", 0);
    assert_eval_eq("Join[{}, {}]", "{}", 0);
}

void test_join_non_list_head() {
    /* Join works on any head, not just List */
    assert_eval_eq("Join[f[a, b], f[c, d]]", "f[a, b, c, d]", 0);
}

void test_join_mismatched_heads() {
    /* Mismatched heads: should remain unevaluated */
    assert_eval_eq("Join[{a, b}, f[c, d]]", "Join[{a, b}, f[c, d]]", 0);
}

void test_join_level2_matrices() {
    /* Join columns of two matrices */
    assert_eval_eq("Join[{{a, b}, {c, d}}, {{1, 2}, {3, 4}}, 2]",
                   "{{a, b, 1, 2}, {c, d, 3, 4}}", 0);
}

void test_join_level2_ragged() {
    /* Ragged arrays: successive elements at level 2 are concatenated */
    assert_eval_eq("Join[{{1}, {5, 6}}, {{2, 3}, {7}}, {{4}, {8}}, 2]",
                   "{{1, 2, 3, 4}, {5, 6, 7, 8}}", 0);
}

void test_join_level2_ragged_unequal_lengths() {
    /* When one list has fewer rows, extra rows pass through */
    assert_eval_eq("Join[{{x}}, {{1, 2}, {3, 4}}, 2]",
                   "{{x, 1, 2}, {3, 4}}", 0);
}

/* Evaluate `input` and return the printed result (caller frees). */
static char* eval_to_string(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    expr_free(evaluated);
    return str;
}

void test_gather() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* Groups are ordered by FIRST occurrence of their element (1, 7, 3, 2,
         * 9 here), not by sorted order; within a group, input order is kept. */
        {"Gather[{1, 7, 3, 7, 2, 3, 9}]", "{{1}, {7, 7}, {3, 3}, {2}, {9}}"},
        {"Gather[{}]", "{}"},
        {"Gather[{5}]", "{{5}}"},
        {"Gather[{2, 2, 2}]", "{{2, 2, 2}}"},
        {"Gather[{1, 2, 3}]", "{{1}, {2}, {3}}"},
        /* Equal elements are collected from anywhere in the list, so the two
         * non-adjacent a's share a group and that group leads because a occurs
         * first. This is Gather, not Split. */
        {"Gather[{a, b, a}]", "{{a, a}, {b}}"},
        {"Gather[{3, 1, 3, 1, 2}]", "{{3, 3}, {1, 1}, {2}}"},
        {"Gather[{x, y, x, x}]", "{{x, x, x}, {y}}"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }

    /* Gather[list] must agree with GatherBy[list, Identity] on every case
     * above — they share one grouping engine, and this pins that down. */
    const char* args[] = {
        "{1, 7, 3, 7, 2, 3, 9}", "{}", "{5}", "{2, 2, 2}",
        "{1, 2, 3}", "{a, b, a}", "{3, 1, 3, 1, 2}", "{x, y, x, x}",
    };
    for (int i = 0; i < (int)(sizeof(args) / sizeof(args[0])); i++) {
        char gather_in[128], gatherby_in[128];
        snprintf(gather_in,   sizeof(gather_in),   "Gather[%s]", args[i]);
        snprintf(gatherby_in, sizeof(gatherby_in), "GatherBy[%s, Identity]", args[i]);
        char* g  = eval_to_string(gather_in);
        char* gb = eval_to_string(gatherby_in);
        if (strcmp(g, gb) != 0) {
            printf("Gather/GatherBy mismatch for %s: Gather gave %s, "
                   "GatherBy[..., Identity] gave %s\n", args[i], g, gb);
            ASSERT(0);
        }
        free(g);
        free(gb);
    }

    /* Non-list arguments stay unevaluated; wrong arity does too. */
    assert_eval_eq("Gather[x]", "Gather[x]", 0);
    assert_eval_eq("Gather[{1, 1}, foo]", "Gather[{1, 1}, foo]", 0);
}

void test_subsets() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* Power set: increasing length, lexicographic by position within a
         * length. Empty subset first, full set last. */
        {"Subsets[{a, b, c}]",
         "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, c}, {a, b, c}}"},
        {"Subsets[{a, b, c, d}]",
         "{{}, {a}, {b}, {c}, {d}, {a, b}, {a, c}, {a, d}, {b, c}, {b, d}, "
         "{c, d}, {a, b, c}, {a, b, d}, {a, c, d}, {b, c, d}, {a, b, c, d}}"},
        {"Subsets[{}]", "{{}}"},
        {"Subsets[{a}]", "{{}, {a}}"},

        /* Subsets[list, All] behaves as Subsets[list]. */
        {"Subsets[{a, b, c}, All]",
         "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, c}, {a, b, c}}"},

        /* Subsets[list, n] — lengths 0..n inclusive. */
        {"Subsets[{a, b, c}, 2]", "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, c}}"},
        {"Subsets[{a, b, c}, 0]", "{{}}"},
        /* n past the end clamps to the full power set. */
        {"Subsets[{a, b, c}, 7]",
         "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, c}, {a, b, c}}"},
        /* Negative n selects nothing. */
        {"Subsets[{a, b, c}, -1]", "{}"},

        /* Subsets[list, {n}] — exactly length n. */
        {"Subsets[{a, b, c}, {2}]", "{{a, b}, {a, c}, {b, c}}"},
        {"Subsets[{a, b, c}, {0}]", "{{}}"},
        {"Subsets[{a, b, c}, {3}]", "{{a, b, c}}"},
        /* Exactly-n does NOT clamp: no subsets are that long. */
        {"Subsets[{a, b, c}, {5}]", "{}"},

        /* Subsets[list, {nmin, nmax}] — inclusive length range. */
        {"Subsets[{a, b, c}, {1, 2}]", "{{a}, {b}, {c}, {a, b}, {a, c}, {b, c}}"},
        {"Subsets[{a, b, c}, {2, 3}]", "{{a, b}, {a, c}, {b, c}, {a, b, c}}"},
        /* An over-long upper bound clamps. */
        {"Subsets[{a, b, c}, {2, 9}]", "{{a, b}, {a, c}, {b, c}, {a, b, c}}"},
        {"Subsets[{a, b, c}, {1, Infinity}]",
         "{{a}, {b}, {c}, {a, b}, {a, c}, {b, c}, {a, b, c}}"},
        /* nmin > nmax selects nothing. */
        {"Subsets[{a, b, c}, {3, 1}]", "{}"},
        {"Subsets[{a, b, c}, {-2, -1}]", "{}"},
        /* Three-element spec adds a length step. */
        {"Subsets[{a, b, c}, {0, 3, 2}]", "{{}, {a, b}, {a, c}, {b, c}}"},

        /* Subsets[list, spec, s] — first s subsets, same order. */
        {"Subsets[{a, b, c}, All, 3]", "{{}, {a}, {b}}"},
        {"Subsets[{a, b, c}, {2}, 2]", "{{a, b}, {a, c}}"},
        {"Subsets[{a, b, c}, All, 0]", "{}"},
        /* s past the number of available subsets returns all of them. */
        {"Subsets[{a, b}, All, 99]", "{{}, {a}, {b}, {a, b}}"},
        {"Subsets[{a, b}, All, All]", "{{}, {a}, {b}, {a, b}}"},

        /* The head of the input is preserved on the inner subsets; the outer
         * wrapper is always a List. */
        {"Subsets[f[a, b]]", "{f[], f[a], f[b], f[a, b]}"},
        {"Subsets[f[a, b, c], {2}]", "{f[a, b], f[a, c], f[b, c]}"},

        /* Duplicates are distinct by position — no dedup. */
        {"Subsets[{a, a}]", "{{}, {a}, {a}, {a, a}}"},
        {"Subsets[{1, 1, 1}, {2}]", "{{1, 1}, {1, 1}, {1, 1}}"},

        /* Elements are evaluated normally before being distributed. */
        {"Subsets[{1 + 1, 2 * 3}, {1}]", "{{2}, {6}}"},

        /* Not a length spec at all: the call stays unevaluated. */
        {"Subsets[{a, b}, x]", "Subsets[{a, b}, x]"},
        /* An atom has no sublists. */
        {"Subsets[5]", "Subsets[5]"},

        /* PERFORMANCE: the 3-argument form must generate lazily. A 40-element
         * list has 2^40 subsets; materializing them would never finish, so
         * this returning promptly is itself the assertion. */
        {"Subsets[Range[40], All, 5]",
         "{{}, {1}, {2}, {3}, {4}}"},
        {"Subsets[Range[40], {20}, 3]",
         "{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}, "
         "{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21}, "
         "{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22}}"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }
}

/* Riffle — every row here is one row of the acceptance table on md-czh.1. */
void test_riffle() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* A non-list separator goes in every gap. n elements, n - 1 gaps. */
        {"Riffle[{1, 2, 3}, 0]", "{1, 0, 2, 0, 3}"},
        {"Riffle[{1, 2}, 0]", "{1, 0, 2}"},
        /* No gaps to fill, so the separator is ignored entirely. */
        {"Riffle[{1}, 0]", "{1}"},
        {"Riffle[{}, 0]", "{}"},
        {"Riffle[{a, b, c}, x]", "{a, x, b, x, c}"},

        /* A List separator is consumed cyclically, left to right: 3 gaps take
         * x, y, x. */
        {"Riffle[{a, b, c, d}, {x, y}]", "{a, x, b, y, c, x, d}"},
        /* Only 2 gaps here, so z is never used. */
        {"Riffle[{a, b, c}, {x, y, z}]", "{a, x, b, y, c}"},
        {"Riffle[{1, 2, 3}, {x}]", "{1, x, 2, x, 3}"},
        {"Riffle[{a}, {x, y}]", "{a}"},
        {"Riffle[{1, 2, 3}, {0, 0}]", "{1, 0, 2, 0, 3}"},

        /* An empty separator list supplies nothing (the k == 0 case) — must
         * not divide by zero computing the cycling index. */
        {"Riffle[{a, b}, {}]", "{a, b}"},
        /* n == 0 reached with a list separator: the no-gaps check has to come
         * before any 2n - 1 sizing, which would underflow size_t. */
        {"Riffle[{}, {x, y}]", "{}"},
        /* More separators than gaps: both z and w go unused. */
        {"Riffle[{a, b, c}, {x, y, z, w}]", "{a, x, b, y, c}"},
        /* Exactly one gap, which takes only the first separator. */
        {"Riffle[{a, b}, {x, y}]", "{a, x, b}"},
        /* The head of the first argument is preserved, not forced to List. */
        {"Riffle[f[a, b], x]", "f[a, x, b]"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }
}

void test_subdivide() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* Subdivide[n] — n + 1 points spanning 0 to 1. n counts the parts the
         * interval is cut into, not the points produced. */
        {"Subdivide[1]", "{0, 1}"},
        {"Subdivide[2]", "{0, 1/2, 1}"},
        {"Subdivide[3]", "{0, 1/3, 2/3, 1}"},
        /* Rationals reduce: the middle point is 1/2, not 2/4. */
        {"Subdivide[4]", "{0, 1/4, 1/2, 3/4, 1}"},

        /* Subdivide[max, n] — spans 0 to max. */
        {"Subdivide[10, 5]", "{0, 2, 4, 6, 8, 10}"},
        /* Step 5/2: points landing on whole numbers print as integers,
         * mixed with rationals in the same list. */
        {"Subdivide[10, 4]", "{0, 5/2, 5, 15/2, 10}"},

        /* Subdivide[min, max, n] — the lower endpoint is min, not 0. */
        {"Subdivide[1, 3, 4]", "{1, 3/2, 2, 5/2, 3}"},
        {"Subdivide[2, 8, 3]", "{2, 4, 6, 8}"},

        /* DESCENDING INTERVALS. No special case: min + i (max - min)/n is
         * applied verbatim, so max < min simply gives a negative step and
         * still yields n + 1 points. */
        {"Subdivide[3, 1, 4]", "{3, 5/2, 2, 3/2, 1}"},
        /* The two-argument form descends when max is negative. */
        {"Subdivide[-10, 5]", "{0, -2, -4, -6, -8, -10}"},
        /* An interval straddling zero descends through it. */
        {"Subdivide[1, -1, 4]", "{1, 1/2, 0, -1/2, -1}"},
        /* Degenerate interval: step 0 repeats the endpoint. n is a positive
         * integer, which is the only validity condition, so this evaluates. */
        {"Subdivide[5, 5, 2]", "{5, 5, 5}"},

        /* n must be a positive integer; otherwise the call stays
         * unevaluated and prints back exactly as written. */
        {"Subdivide[0]", "Subdivide[0]"},
        {"Subdivide[-1]", "Subdivide[-1]"},
        {"Subdivide[5, 0]", "Subdivide[5, 0]"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }
}

void test_nearest() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        /* ALL tied elements are returned, not a single minimum. Both 1 and 5
         * sit at distance 2 from 3. This is the row the whole design turns on:
         * a quickselect-style tiebreak (RankedMin's ranked_cmp, sort.c:932)
         * exists to make ties impossible and would answer {1}. */
        {"Nearest[{1, 5, 10}, 3]", "{1, 5}"},
        {"Nearest[{4}, 100]", "{4}"},

        /* Ties come back in INPUT order, not sorted order. */
        {"Nearest[{5, 1, 10}, 3]", "{5, 1}"},
        {"Nearest[{-5, 5}, 0]", "{-5, 5}"},
        /* Exact rational tie: both distances are exactly 1/6. A float-keyed
         * comparison could miss this. */
        {"Nearest[{1/3, 2/3}, 1/2]", "{1/3, 2/3}"},

        /* MIXED-TYPE TIES. Distances equal in VALUE but of different ExprType
         * must still tie. These are the rows that fail if the comparison is
         * expr_compare: its canonical order breaks a value tie on the type
         * enum (sort.c:376), so Integer beats Real beats Rational and only one
         * of the tied elements survives. */
        {"Nearest[{0, 2.0}, 1]", "{0, 2.0}"},        /* dists 1 and 1.0   */
        {"Nearest[{3, 1.0}, 2]", "{3, 1.0}"},        /* dists 1 and 1.0   */
        {"Nearest[{-1, 1.0}, 0]", "{-1, 1.0}"},      /* dists 1 and 1.0   */
        {"Nearest[{1.5, 5/2}, 2]", "{1.5, 5/2}"},    /* dists 0.5 and 1/2 */

        /* The opposite direction: distances that differ below double
         * resolution must NOT tie. expr_compare falls back to comparing
         * get_numeric_value() doubles for atoms that are not both
         * integer-like (sort.c:372-377), which would report both as nearest. */
        {"Nearest[{1/3, 1/3 + 1/10^18}, 0]", "{1/3}"},

        /* Bigint elements order exactly. */
        {"Nearest[{10^25, 1}, 0]", "{1}"},
        /* Duplicates are distinct by position, as in Subsets. */
        {"Nearest[{1, 1, 2}, 1]", "{1, 1}"},
        /* An exact hit does not short-circuit: a later element could still tie
         * at distance 0, so the collect pass always runs to the end. */
        {"Nearest[{2, 4, 6}, 4]", "{4}"},

        /* Unique nearest. */
        {"Nearest[{1, 2}, 1.4]", "{1}"},
        {"Nearest[{10, 20, 30}, 100]", "{30}"},
        /* Abs of a complex difference is its modulus, so the complex case
         * falls out of the composition: 5 versus 1. */
        {"Nearest[{3 + 4 I, 1}, 0]", "{1}"},

        /* Empty in, empty out -- checked before the numeric gate, so a
         * symbolic target on an empty list is still {}. */
        {"Nearest[{}, 3]", "{}"},
        {"Nearest[{}, a]", "{}"},

        /* THE GATE. A non-real distance means no definite answer, so the whole
         * call stays unevaluated. Note MinimalBy[{1, a, 3}, Abs[# - 2] &] gives
         * {1, 3} -- expr_compare orders symbols after all numbers, so it drops
         * the symbolic element and answers anyway. That plausible wrong answer
         * is what this row exists to prevent. */
        {"Nearest[{1, a, 3}, 2]", "Nearest[{1, a, 3}, 2]"},
        {"Nearest[{1, 2, 3}, a]", "Nearest[{1, 2, 3}, a]"},
        /* A symbolic REAL is rejected too: Abs[Pi - 3] stays as Abs[-3 + Pi].
         * Numericalizing it (as RankedMin's ranked_numeric_key would) is a
         * deliberate follow-up, not current behaviour. */
        {"Nearest[{Pi, 4}, 3]", "Nearest[{Pi, 4}, 3]"},
        /* A rational with a BIGINT component declines, for a reason that is
         * not Nearest's: builtin_abs does not evaluate one, so the distance
         * comes back as an unevaluated Abs[...] and the numeric gate rejects
         * it. Abs[1/1000] is 1/1000 but Abs[1/10^25] is Abs[1/10^25], and
         * Sign has the same gap. Pinned here so this row flips the day
         * builtin_abs is fixed, rather than the limitation going unnoticed. */
        {"Nearest[{1/10^25, 1}, 0]",
         "Nearest[{1/10000000000000000000000000, 1}, 0]"},

        /* MIXED EXACT / INEXACT DISTANCES are compared as exact rationals,
         * not by subtracting. Subtraction widens the pair to a double and
         * loses the exact operand, which made the comparator INTRANSITIVE
         * with plain int64 values: for a = 2^60, b = 2.0^60, c = 2^60 + 1,
         * a - b and c - b both gave 0.0 while a - c gave -1, so a == b,
         * c == b and a < c all held at once. An intransitive comparator
         * feeding a sort yields order-dependent output.
         *
         * A double is exactly a rational, so lifting both sides with mpq is
         * exact and total. The two rows below pin the semantics against
         * Mathematica, which agrees on both: 1 and 1.0 tie (1.0 lifts to
         * exactly 1), 0.1 and 1/10 do not, and Nearest[{0.9, 11/10}, 1] is
         * {0.9} there too. */
        {"Nearest[{2^60, 2.0^60}, 2^60]", "{1152921504606846976, 1.15292e+18}"},
        {"Nearest[{2^60 + 1, 2.0^60}, 2^60]", "{1.15292e+18}"},
        {"Nearest[{0.9, 11/10}, 1]", "{0.9}"},
        {"Nearest[{1/10, 0.1}, 0]", "{1/10}"},

        /* Arity. The 3-argument n-nearest form is a follow-up and is inert. */
        {"Nearest[{1, 5, 10}]", "Nearest[{1, 5, 10}]"},
        {"Nearest[{1, 5}, 3, 2]", "Nearest[{1, 5}, 3, 2]"},
        {"Nearest[3, 1]", "Nearest[3, 1]"},
        /* Non-List head follows RankedMin (sort.c:1014), not MinimalBy, which
         * accepts and preserves any head. */
        {"Nearest[f[1, 5], 3]", "Nearest[f[1, 5], 3]"},
        /* A visible NDArray is never materialised by the transparency gate, so
         * is_listq is the only guard against a silently truncated answer.
         * Unevaluated is the correct conservative result. */
        {"Nearest[NDArray[{1., 5., 10.}], 3.]",
         "Nearest[NDArray[{1.0, 5.0, 10.0}], 3.0]"},

        /* A PACKED list, by contrast, is materialised on the way in because
         * Nearest is not on pack.c's AWARE list. */
        {"Nearest[Range[5], 3]", "{3}"},

        {"Attributes[Nearest]", "{Protected}"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }
}

/* FindClusters. Every expected value below was produced by the built binary and
 * pasted in -- none is hand-predicted. The Nearest table was written the other
 * way for its mixed-type rows, and that is exactly where a high-severity bug
 * survived 22 passing tests. */
void test_find_clusters() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {

        /* COUNT MODES. The exact/bounded distinction: n forces exactly n
         * clusters, UpTo[n] gives the natural count when that is already at or
         * below n. Both are capped at the distinct-value count -- no method can
         * separate two equal elements. */
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}]", "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 4]",
          "{{1, 2, 3, 1}, {10}, {12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[4]]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[2]]",
          "{{1, 2, 10, 12, 3, 1, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 2]", "{{1, 2, 10, 12, 3, 1, 13}, {25}}"},
        {"FindClusters[{1, 2, 3}, 5]", "{{1}, {2}, {3}}"},
        {"FindClusters[{1, 2, 3}, UpTo[5]]", "{{1, 2, 3}}"},
        {"FindClusters[{7, 7, 7, 7}, 3]", "{{7, 7, 7, 7}}"},
        {"FindClusters[{1, 2, 3}, 1]", "{{1, 2, 3}}"},
        {"FindClusters[{1, 2, 3, 5, 8, 9, 10}, 2]", "{{1, 2, 3, 5}, {8, 9, 10}}"},
        {"FindClusters[{1, 2, 3, 5, 8, 9, 10}, 3]", "{{1, 2, 3}, {5}, {8, 9, 10}}"},
        {"FindClusters[{1, 2, 3, 5, 8, 9, 10}, 4]", "{{1}, {2, 3}, {5}, {8, 9, 10}}"},

        /* A count spec this builtin does not understand leaves the whole call
         * alone rather than guessing. */
        {"FindClusters[{1, 2, 3}, 0]", "FindClusters[{1, 2, 3}, 0]"},
        {"FindClusters[{1, 2, 3}, -1]", "FindClusters[{1, 2, 3}, -1]"},
        {"FindClusters[{1, 2, 3}, x]", "FindClusters[{1, 2, 3}, x]"},
        {"FindClusters[{1, 2, 3}, UpTo[0]]", "FindClusters[{1, 2, 3}, UpTo[0]]"},
        {"FindClusters[{1, 2, 3}, 2.5]", "FindClusters[{1, 2, 3}, 2.5]"},

        /* THE 10 x 3 CAPABILITY MATRIX, one row per cell. KMeans and KMedoids
         * need a count; the five density methods need Automatic; Spectral takes
         * Automatic and UpTo[n] but not a bare n. Transcribed from the
         * allowed-lists in Mathematica's three error messages, which contradict
         * its own documentation bullets on Spectral -- the runtime wins. */
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"Agglomerate\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"Agglomerate\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"Agglomerate\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"SpanningTree\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"SpanningTree\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"SpanningTree\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"KMeans\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"KMeans\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"KMeans\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"KMeans\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"KMedoids\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"KMedoids\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"KMedoids\"]",
          "{{1, 2, 3, 1}, {10}, {12, 13, 25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"KMedoids\"]",
          "{{1, 2, 3, 1}, {10}, {12, 13, 25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"Spectral\"]",
          "{{1, 2, 3, 1}, {10}, {12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"Spectral\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"Spectral\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"Spectral\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"DBSCAN\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"DBSCAN\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"DBSCAN\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"DBSCAN\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"DBSCAN\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"GaussianMixture\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"GaussianMixture\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"GaussianMixture\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"GaussianMixture\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"GaussianMixture\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"JarvisPatrick\"]",
          "{{1, 2, 10, 12, 3, 1, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"JarvisPatrick\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"JarvisPatrick\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"JarvisPatrick\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"JarvisPatrick\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"MeanShift\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"MeanShift\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"MeanShift\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"MeanShift\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"MeanShift\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> \"NeighborhoodContraction\"]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"NeighborhoodContraction\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[3], Method -> \"NeighborhoodContraction\"]"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"NeighborhoodContraction\"]",
          "FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> \"NeighborhoodContraction\"]"},

        /* Method plumbing: Automatic, an unknown name, suboptions, and the
         * options accepted for compatibility. */
        {"FindClusters[{1, 2, 3}, Method -> \"Nonsense\"]",
          "FindClusters[{1, 2, 3}, Method -> \"Nonsense\"]"},
        {"FindClusters[{1, 2, 3}, Method -> Automatic]", "{{1, 2, 3}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> {\"DBSCAN\", \"NeighborhoodRadius\" -> 0.5}]",
          "{{1, 1}, {2}, {10}, {12}, {3}, {13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> {\"DBSCAN\", \"MinPoints\" -> 3}]",
          "{{1, 2, 3, 1}, {10, 12, 13}, {25}}"},
        {"FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, Method -> {\"JarvisPatrick\", \"NeighborCount\" -> 2}]",
          "{{1, 2, 1}, {10, 12, 13}, {3}, {25}}"},
        {"FindClusters[{1, 2, 10, 12}, Method -> {\"DBSCAN\", \"Bogus\" -> 2}]",
          "FindClusters[{1, 2, 10, 12}, Method -> {\"DBSCAN\", \"Bogus\" -> 2}]"},
        {"FindClusters[{1, 2, 10, 12}, DistanceFunction -> EuclideanDistance]",
          "{{1, 2}, {10, 12}}"},
        {"FindClusters[{1, 2, 10, 12}, DistanceFunction -> ManhattanDistance]",
          "{{1, 2}, {10, 12}}"},
        {"FindClusters[{1, 2, 10, 12}, DistanceFunction -> Sin]",
          "FindClusters[{1, 2, 10, 12}, DistanceFunction -> Sin]"},
        {"FindClusters[{1, 2, 10, 12}, CriterionFunction -> \"Silhouette\"]",
          "{{1, 2}, {10, 12}}"},
        {"FindClusters[{1, 2, 10, 12}, PerformanceGoal -> \"Speed\"]", "{{1, 2}, {10, 12}}"},

        /* The Automatic cluster count -- gaps wider than FC_GAP_FACTOR times
         * the median gap are cut. Uniform spacing and identical elements must
         * NOT be split. */
        {"FindClusters[{1, 2, 3, 100, 101, 102}]", "{{1, 2, 3}, {100, 101, 102}}"},
        {"FindClusters[{1, 2, 3, 20, 21, 22, 500}]", "{{1, 2, 3}, {20, 21, 22}, {500}}"},
        {"FindClusters[{1, 2, 3, 4, 5, 6, 7, 8}]", "{{1, 2, 3, 4, 5, 6, 7, 8}}"},
        {"FindClusters[{7, 7, 7, 7}]", "{{7, 7, 7, 7}}"},
        {"FindClusters[{5}]", "{{5}}"},
        {"FindClusters[{-10, -9, 5, 6}]", "{{-10, -9}, {5, 6}}"},

        /* DELIBERATE DIVERGENCES FROM MATHEMATICA, pinned so they surface as a
         * diff if anyone changes them. Mathematica's Automatic count is an
         * unpublished internal index, not a gap rule: it splits {1,4,9,16,25,36}
         * three ways and {5,6} into two singletons. It also clusters symbolic
         * elements as nominal features; we are numeric-only and decline. */
        {"FindClusters[{1, 4, 9, 16, 25, 36}]", "{{1, 4, 9, 16, 25, 36}}"},
        {"FindClusters[{5, 6}]", "{{5, 6}}"},
        {"FindClusters[{1, a, 3}, 2]", "FindClusters[{1, a, 3}, 2]"},

        /* EXACT ARITHMETIC. Ordering and fixed-count gap selection run on the
         * elements themselves via list_numeric_cmp, never on a double
         * projection, so bigint-component rationals work -- Nearest cannot do
         * this, because it routes its distance through Abs, which declines on
         * one (complex.c:418-421). The sub-double row fails if the comparison
         * ever degrades to doubles; the 1/10^25 row fails if a gap is ever
         * rewritten as Abs[b - a]. */
        {"FindClusters[{1/3, 2/3, 10, 31/3}, 2]", "{{1/3, 2/3}, {10, 31/3}}"},
        {"FindClusters[{1/10^25, 2/10^25, 1}, 2]",
          "{{1/10000000000000000000000000, 1/5000000000000000000000000}, {1}}"},
        {"FindClusters[{1/3, 1/3 + 1/10^18, 5}, 2]",
          "{{1/3, 1000000000000000003/3000000000000000000}, {5}}"},
        {"FindClusters[{1, 2.0, 10, 11.0}, 2]", "{{1, 2.0}, {10, 11.0}}"},
        {"FindClusters[{10^25, 1, 2}, 2]", "{{10000000000000000000000000}, {1, 2}}"},

        /* EXACTLY-n WITH HEAVY TIES. Quantile seeding over raw sorted
         * positions puts two centroids on the same value when the data is
         * tie-heavy, and a centroid that never wins a point leaves an empty
         * cluster -- so KMeans and KMedoids returned FEWER clusters than
         * asked for, on data whose distinct count was not the limit.
         * Measured before the fix: {1,1,1,1,2,3} with n=3 gave two clusters.
         * The partition invariant below cannot catch this -- no element is
         * lost, only the promised count is wrong -- so these rows and the
         * exactly-n loop that follows are the guard. */
        {"FindClusters[{1, 1, 1, 1, 2, 3}, 3, Method -> \"KMeans\"]",
          "{{1, 1, 1, 1}, {2}, {3}}"},
        {"FindClusters[{1, 1, 1, 1, 2, 3}, 3, Method -> \"KMedoids\"]",
          "{{1, 1, 1, 1}, {2}, {3}}"},
        {"FindClusters[{1, 1, 1, 1, 2, 3}, 3]", "{{1, 1, 1, 1}, {2}, {3}}"},
        {"FindClusters[{5, 5, 5, 6, 7}, 3, Method -> \"KMeans\"]", "{{5, 5, 5}, {6}, {7}}"},
        {"FindClusters[{1, 1, 2, 2, 3, 3}, 3, Method -> \"KMeans\"]",
          "{{1, 1}, {2, 2}, {3, 3}}"},
        {"FindClusters[{2, 2, 2, 2, 2, 9}, 2, Method -> \"KMedoids\"]",
          "{{2, 2, 2, 2, 2}, {9}}"},

        /* THE ORDER CONTRACT: clusters by first occurrence, elements in input
         * order. These fail if anyone sorts the output. */
        {"FindClusters[{10, 1, 11, 2}]", "{{10, 11}, {1, 2}}"},
        {"FindClusters[{3, 1, 2}, 1]", "{{3, 1, 2}}"},
        {"FindClusters[{1, 1, 2}, 2]", "{{1, 1}, {2}}"},

        /* REGRESSION ROWS FROM CODE REVIEW. Each of these returned a wrong
         * answer before the fix named beside it; none was reachable from the
         * rows above, which is why the review found them and the table did
         * not.
         *
         * Equal elements must never be separated -- JarvisPatrick returned
         * four clusters for nine copies of 7, and DBSCAN split a duplicated
         * pair whenever MinPoints exceeded its multiplicity. Now enforced
         * centrally on exact zero gaps, so no method can reintroduce it.
         *
         * Exactly-n must hold at scales the double projection cannot see:
         * 2^60 and 2^60+1 are distinct but project to one double, and the
         * repair loop hunted boundaries in double space, so KMeans returned
         * ONE cluster for three distinct values.
         *
         * The Automatic threshold took the mean of the two middle gaps for
         * an even-length list, which makes a split arithmetically impossible
         * for two or three elements: the largest gap is itself one of the
         * values averaged. {0, 1, 10^12} was a single cluster. The lower
         * median fixes it, and computing it on the gap Exprs removes the
         * double projection -- and with it the scale cliff where 10^25
         * worked and 10^400 did not. The 10^25 and 10^400 rows now agree,
         * which is the point: one cluster there is the median rule, not an
         * overflow.
         *
         * MeanShift fragmented evenly spaced data into near-singletons
         * (Range[40] gave 34 clusters) because every interior point of a
         * flat density is stationary and the merge tolerance was tied to the
         * bandwidth. UpTo[n] was identical to a bare n for KMeans/KMedoids,
         * collapsing the three-mode design to two for four methods. And the
         * two quadratic methods now decline rather than appearing to hang. */
        {"FindClusters[{7, 7, 7, 7, 7, 7, 7, 7, 7}, Method -> \"JarvisPatrick\"]",
          "{{7, 7, 7, 7, 7, 7, 7, 7, 7}}"},
        {"FindClusters[{7, 7, 100}, Method -> {\"JarvisPatrick\", \"NeighborCount\" -> 1}]",
          "{{7, 7}, {100}}"},
        {"FindClusters[{7, 7}, Method -> {\"DBSCAN\", \"NeighborhoodRadius\" -> 10.0, \"MinPoints\" -> 5}]",
          "{{7, 7}}"},
        {"FindClusters[{1, 1, 100, 101, 102, 103}, Method -> {\"DBSCAN\", \"MinPoints\" -> 3}]",
          "{{1, 1}, {100, 101, 102, 103}}"},
        {"FindClusters[{2^60, 2^60 + 1, 2^60 + 2}, 3, Method -> \"KMeans\"]",
          "{{1152921504606846976}, {1152921504606846977}, {1152921504606846978}}"},
        {"FindClusters[{2^60, 2^60 + 1, 2^60 + 2}, 3, Method -> \"KMedoids\"]",
          "{{1152921504606846976}, {1152921504606846977}, {1152921504606846978}}"},
        {"FindClusters[{1/3, 1/3 + 1/10^18, 5}, 3, Method -> \"KMeans\"]",
          "{{1/3}, {1000000000000000003/3000000000000000000}, {5}}"},
        {"FindClusters[{0, 1, 10^12}]", "{{0, 1}, {1000000000000}}"},
        {"FindClusters[{1, 2, 1000000}]", "{{1, 2}, {1000000}}"},
        {"FindClusters[{0, 1, 10^25, 2*10^25}]",
          "{{0, 1, 10000000000000000000000000, 20000000000000000000000000}}"},
        {"Length[FindClusters[{0, 1, 10^400, 2*10^400}]]", "1"},
        {"Length[FindClusters[Range[40], Method -> \"MeanShift\"]]", "1"},
        {"FindClusters[{1, 1, 1, 1, 100}, Method -> \"MeanShift\"]", "{{1, 1, 1, 1}, {100}}"},
        {"FindClusters[{1, 1, 1, 1, 100}, Method -> \"NeighborhoodContraction\"]",
          "{{1, 1, 1, 1}, {100}}"},
        {"FindClusters[{1, 2, 3, 4, 5, 6, 7, 8}, UpTo[3], Method -> \"KMeans\"]",
          "{{1, 2, 3, 4, 5, 6, 7, 8}}"},
        {"FindClusters[{1, 2, 3}, UpTo[5], Method -> \"KMeans\"]", "{{1, 2, 3}}"},
        {"Options[FindClusters]",
          "{Method -> Automatic, DistanceFunction -> Automatic, CriterionFunction -> Automatic, PerformanceGoal -> Automatic}"},
        {"Head[FindClusters[Range[5000], Method -> \"MeanShift\"]]", "FindClusters"},
        {"Head[FindClusters[Range[5000], Method -> \"NeighborhoodContraction\"]]",
          "FindClusters"},
        {"Head[FindClusters[Range[3000], Method -> \"Spectral\"]]", "FindClusters"},
        {"Length[FindClusters[Range[1000], Method -> \"MeanShift\"]]", "1"},

        /* Shape guards. A visible NDArray is not a List and is never
         * materialised by the transparency gate, so it declines rather than
         * being silently truncated; a packed list is materialised on the way in. */
        {"FindClusters[{}]", "FindClusters[{}]"},
        {"FindClusters[5]", "FindClusters[5]"},
        {"FindClusters[f[1, 2, 3]]", "FindClusters[f[1, 2, 3]]"},
        {"FindClusters[NDArray[{1., 2., 10.}]]", "FindClusters[NDArray[{1.0, 2.0, 10.0}]]"},
        {"FindClusters[Range[10]]", "{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}}"},
        {"FindClusters[]", "FindClusters[]"},

        /* --- n-dimensional numeric vectors ---------------------------
         * A list of equal-length numeric vectors clusters in any
         * dimension. The structure is a real minimum spanning tree over
         * exact SQUARED distances; in one dimension that tree is the
         * sorted adjacency chain, which is why every scalar row above
         * is unchanged. Ragged rows, a scalar mixed among vectors,
         * depth over 2 and a symbolic component all decline. The
         * methods that still read the sorted 1-D projection (KMeans,
         * DBSCAN, ...) decline on vectors rather than running on
         * meaningless data. Cluster ORDER follows first occurrence in
         * input order, so it differs from Mathematica while the
         * partition agrees. */
        {"FindClusters[{{1, 2}, {3, 4}}]", "{{{1, 2}, {3, 4}}}"},
        {"FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}}]", "{{{1, 1}, {1, 2}}, {{9, 9}, {9, 8}}}"},
        {"FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}}, 2]", "{{{1, 1}, {1, 2}}, {{9, 9}, {9, 8}}}"},
        {"FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}}, 3]", "{{{1, 1}}, {{1, 2}}, {{9, 9}, {9, 8}}}"},
        {"FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}}, 4]", "{{{1, 1}}, {{1, 2}}, {{9, 9}}, {{9, 8}}}"},
        {"FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}}, UpTo[3]]", "{{{1, 1}, {1, 2}}, {{9, 9}, {9, 8}}}"},
        {"FindClusters[{{2.5, 3.1}, {5.9, 3.4}, {2.6, 3.0}, {6.1, 3.5}}]", "{{{2.5, 3.1}, {2.6, 3.0}}, {{5.9, 3.4}, {6.1, 3.5}}}"},
        {"FindClusters[{{1, 1, 1}, {1, 1, 2}, {9, 9, 9}}]", "{{{1, 1, 1}, {1, 1, 2}}, {{9, 9, 9}}}"},
        {"FindClusters[{{1}, {2}, {100}}]", "{{{1}, {2}}, {{100}}}"},
        {"FindClusters[{{1}, {2}, {100}}, 3]", "{{{1}}, {{2}}, {{100}}}"},
        {"FindClusters[{{1, 1}, {1, 1}, {9, 9}}, 3]", "{{{1, 1}, {1, 1}}, {{9, 9}}}"},
        {"FindClusters[{{1, 1}, {1, 1}, {1, 1}}]", "{{{1, 1}, {1, 1}, {1, 1}}}"},
        {"FindClusters[{{1, 2}, {1, 2}, {9, 9}, {9, 9}}, 4]", "{{{1, 2}, {1, 2}}, {{9, 9}, {9, 9}}}"},
        {"FindClusters[{{1/3, 1/7}, {1/3, 1/7 + 1/10^20}, {5, 5}}]", "{{{1/3, 1/7}, {1/3, 100000000000000000007/700000000000000000000}}, {{5, 5}}}"},
        {"FindClusters[{{0, 0}, {0, 1}, {1, 0}, {1, 1}, {50, 50}}]", "{{{0, 0}, {0, 1}, {1, 0}, {1, 1}}, {{50, 50}}}"},
        {"FindClusters[{{1, 1}, {2, 2}, {3, 3}}, Method -> \"SpanningTree\"]", "{{{1, 1}, {2, 2}, {3, 3}}}"},
        {"FindClusters[{{1, 1}, {2, 2}, {30, 30}}, Method -> \"Agglomerate\"]", "{{{1, 1}, {2, 2}}, {{30, 30}}}"},
        {"FindClusters[{{1, 2}, {3}}]", "FindClusters[{{1, 2}, {3}}]"},
        {"FindClusters[{{1, 2}, 3}]", "FindClusters[{{1, 2}, 3}]"},
        {"FindClusters[{3, {1, 2}}]", "FindClusters[{3, {1, 2}}]"},
        {"FindClusters[{{{1}}, {{2}}}]", "FindClusters[{{{1}}, {{2}}}]"},
        {"FindClusters[{{1, a}, {2, 3}}]", "FindClusters[{{1, a}, {2, 3}}]"},
        {"FindClusters[{{}, {}}]", "FindClusters[{{}, {}}]"},
        {"FindClusters[{{1, 1}, {9, 9}}, 2, Method -> \"KMeans\"]", "FindClusters[{{1, 1}, {9, 9}}, 2, Method -> \"KMeans\"]"},
        {"FindClusters[{{1, 1}, {9, 9}}, Method -> \"DBSCAN\"]", "FindClusters[{{1, 1}, {9, 9}}, Method -> \"DBSCAN\"]"},


        /* --- strings and colours -------------------------------------
         * A colour is a compound expression carrying numeric arguments,
         * i.e. the same shape as a vector, so RGBColor and GrayLevel are
         * points in their own space with no colour-specific code.
         * Strings have no coordinates at all and cluster on exact
         * integer edit distances -- the gap methods only ever ask for
         * pairwise distances, never positions.
         *
         * The {1/2, 1/3, 10} row is a guard, not a feature: Rational is
         * itself a compound expression, so reading any compound head as
         * a point would cluster 1/2 by numerator and denominator. Heads
         * are an explicit list and scalars are tested first. Complex is
         * declined for the same reason, unchanged from before. */
        {"FindClusters[{\"GGTTT\", \"GGGGT\", \"AACCC\", \"AAACC\"}, 2]", "{{\"GGTTT\", \"GGGGT\"}, {\"AACCC\", \"AAACC\"}}"},
        {"FindClusters[{\"cat\", \"cot\", \"dog\", \"dig\"}, 2]", "{{\"cat\", \"cot\"}, {\"dog\", \"dig\"}}"},
        {"FindClusters[{\"abc\", \"abc\", \"xyz\"}, 3]", "{{\"abc\", \"abc\"}, {\"xyz\"}}"},
        {"FindClusters[{RGBColor[1, 0, 0], RGBColor[9/10, 1/10, 0], RGBColor[0, 0, 1]}, 2]", "{{RGBColor[1, 0, 0], RGBColor[9/10, 1/10, 0]}, {RGBColor[0, 0, 1]}}"},
        {"FindClusters[{GrayLevel[0], GrayLevel[1/10], GrayLevel[1]}, 2]", "{{GrayLevel[0], GrayLevel[1/10]}, {GrayLevel[1]}}"},
        {"FindClusters[{1/2, 1/3, 10}]", "{{1/2, 1/3}, {10}}"},
        {"FindClusters[{\"a\", 5}]", "FindClusters[{\"a\", 5}]"},
        {"FindClusters[{RGBColor[1, 0, 0], GrayLevel[0]}]", "FindClusters[{RGBColor[1, 0, 0], GrayLevel[0]}]"},
        {"FindClusters[{RGBColor[1, 0, 0], RGBColor[0, 1, 0]}, Method -> \"KMeans\"]", "FindClusters[{RGBColor[1, 0, 0], RGBColor[0, 1, 0]}, Method -> \"KMeans\"]"},
        {"FindClusters[{2 + 3 I, 1 + I}]", "FindClusters[{2 + 3*I, 1 + I}]"},

        {"Attributes[FindClusters]", "{Protected}"},
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        assert_eval_eq(tests[i].input, tests[i].expected, 0);
    }

    /* The partition invariant, as a loop rather than a table: for every method
     * in a mode it accepts, the clusters must together be a permutation of the
     * input -- no element lost, none duplicated. One assertion covering a whole
     * class of bug across ten independent implementations, and the check that
     * would catch a scatter/remap error no single expected-string row would. */
    static const char* auto_methods[] = {
        "Agglomerate", "SpanningTree", "Spectral", "DBSCAN",
        "GaussianMixture", "JarvisPatrick", "MeanShift", "NeighborhoodContraction"
    };
    static const char* counted_methods[] = {
        "Agglomerate", "SpanningTree", "KMeans", "KMedoids"
    };
    char buf[256];
    for (int i = 0; i < (int)(sizeof(auto_methods) / sizeof(auto_methods[0])); i++) {
        snprintf(buf, sizeof buf,
                 "Sort[Flatten[FindClusters[{1,2,10,12,3,1,13,25}, Method -> \"%s\"]]] "
                 "=== Sort[{1,2,10,12,3,1,13,25}]", auto_methods[i]);
        assert_eval_eq(buf, "True", 0);
    }
    for (int i = 0; i < (int)(sizeof(counted_methods) / sizeof(counted_methods[0])); i++) {
        snprintf(buf, sizeof buf,
                 "Sort[Flatten[FindClusters[{1,2,10,12,3,1,13,25}, 3, Method -> \"%s\"]]] "
                 "=== Sort[{1,2,10,12,3,1,13,25}]", counted_methods[i]);
        assert_eval_eq(buf, "True", 0);
    }
    /* The exactly-n contract, as a loop: for a fixed count, every method that
     * accepts one must return exactly Min[n, distinct] clusters. Tie-heavy data
     * is the case that breaks it, and no single expected-string row generalises
     * across methods and counts the way this does. */
    static const char* fixed_methods[] = {
        "Agglomerate", "SpanningTree", "KMeans", "KMedoids"
    };
    static const char* tie_data[] = {
        "{1,1,1,1,2,3}", "{5,5,5,6,7}", "{1,1,2,2,3,3}", "{2,2,2,2,2,9}", "{1,1,1,2}"
    };
    for (int m = 0; m < (int)(sizeof(fixed_methods) / sizeof(fixed_methods[0])); m++) {
        for (int dsi = 0; dsi < (int)(sizeof(tie_data) / sizeof(tie_data[0])); dsi++) {
            for (int nn = 1; nn <= 4; nn++) {
                snprintf(buf, sizeof buf,
                         "Length[FindClusters[%s, %d, Method -> \"%s\"]] == "
                         "Min[%d, Length[Union[%s]]]",
                         tie_data[dsi], nn, fixed_methods[m], nn, tie_data[dsi]);
                assert_eval_eq(buf, "True", 0);
            }
        }
    }
    /* The equal-elements invariant, as a loop. Two methods violated it and a
     * third could have; it is now enforced centrally, and this is what keeps it
     * enforced. No method may return more clusters than the input has distinct
     * values -- which is exactly the statement that identical elements are
     * never separated. */
    static const char* tie_probe[] = {
        "{7,7,7,7,7,7,7,7,7}", "{7,7}", "{1,1,100,101,102,103}",
        "{5,5,5,6,7}", "{2,2,2,2,2,9}"
    };
    for (int m = 0; m < (int)(sizeof(auto_methods) / sizeof(auto_methods[0])); m++) {
        for (int dsi = 0; dsi < (int)(sizeof(tie_probe) / sizeof(tie_probe[0])); dsi++) {
            snprintf(buf, sizeof buf,
                     "Length[FindClusters[%s, Method -> \"%s\"]] <= "
                     "Length[Union[%s]]",
                     tie_probe[dsi], auto_methods[m], tie_probe[dsi]);
            assert_eval_eq(buf, "True", 0);
        }
    }
}

/* The distance functions FindClusters ranks with, and which had no
 * implementation at all before -- Names["*Distance*"] returned only
 * GraphDistance, while FindClusters' DistanceFunction option accepted the names
 * as inert symbols. Every expected value below was cross-checked against
 * Mathematica via wolframscript, including the two conventions that are not
 * derivable: CosineDistance of a zero vector is 0 rather than Indeterminate,
 * and the Abs-then-square definition (not square-then-sum) is what makes a
 * complex component contribute its modulus. */
void test_distance_functions() {
    struct { const char* in; const char* out; } cases[] = {
        {"EuclideanDistance[{1, 2}, {4, 6}]", "5"},
        {"SquaredEuclideanDistance[{1, 2}, {4, 6}]", "25"},
        {"SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]", "58/441"},
        {"ManhattanDistance[{1, 2}, {4, 6}]", "7"},
        {"EuclideanDistance[{0, 0}, {1, 1}]", "Sqrt[2]"},
        {"EuclideanDistance[3, 7]", "4"},
        {"ManhattanDistance[{a}, {b}]", "Abs[a - b]"},
        {"SquaredEuclideanDistance[{3 + 4 I}, {0}]", "25"},
        {"EuclideanDistance[{1, 2}, {1, 2, 3}]", "EuclideanDistance[{1, 2}, {1, 2, 3}]"},
        {"EuclideanDistance[{{1, 2}}, {{3, 4}}]", "EuclideanDistance[{{1, 2}}, {{3, 4}}]"},
        {"CosineDistance[{1, 0}, {0, 1}]", "1"},
        {"CosineDistance[{1, 1}, {1, 0}]", "1 - 1/Sqrt[2]"},
        {"CosineDistance[{1, 2}, {2, 4}]", "0"},
        {"CosineDistance[{1, 0}, {-1, 0}]", "2"},
        {"CosineDistance[{0, 0}, {1, 2}]", "0"},
        {"CosineDistance[3, 4]", "0"},
        {"Attributes[EuclideanDistance]", "{Protected}"},
        {"Attributes[CosineDistance]", "{Protected}"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
        assert_eval_eq(cases[i].in, cases[i].out, 0);

    /* Exact at a scale where Abs itself declines: Abs[-1/10^20] does not
     * evaluate (nor does Sign of it), so composing internal_abs blindly left
     * these unevaluated and FindClusters refused exact high-precision vectors.
     * The distances take the sign through the exact comparator instead. */
    assert_eval_eq("SquaredEuclideanDistance[{1/3, 1/7}, {1/3, 1/7 + 1/10^20}]",
                   "1/10000000000000000000000000000000000000000", 0);
    assert_eval_eq("ManhattanDistance[{0}, {1/10^20}]", "1/100000000000000000000", 0);
    assert_eval_eq("ManhattanDistance[{1/10^20}, {0}]", "1/100000000000000000000", 0);

    /* Sequence distances: strings and lists, cross-checked against
     * Mathematica. HammingDistance requires equal lengths there
     * (::idim) so declining is the faithful behaviour. */
    struct { const char* in; const char* out; } seqs[] = {
        {"EditDistance[\"GGTTT\", \"GGGGT\"]", "2"},
        {"EditDistance[\"kitten\", \"sitting\"]", "3"},
        {"EditDistance[\"\", \"abc\"]", "3"},
        {"EditDistance[\"abc\", \"abc\"]", "0"},
        {"EditDistance[{1, 2, 3}, {1, 3}]", "1"},
        {"HammingDistance[\"GGTTT\", \"GGGGT\"]", "2"},
        {"HammingDistance[{1, 2, 3}, {1, 5, 3}]", "1"},
        {"HammingDistance[\"abc\", \"abcd\"]", "HammingDistance[\"abc\", \"abcd\"]"},
        {"EditDistance[\"abc\", 5]", "EditDistance[\"abc\", 5]"},
        {"Attributes[EditDistance]", "{Protected}"},
    };
    for (size_t i = 0; i < sizeof(seqs) / sizeof(seqs[0]); i++)
        assert_eval_eq(seqs[i].in, seqs[i].out, 0);

    /* Squared Euclidean is monotone in Euclidean -- the property the n-D
     * clustering relies on to rank exactly without taking a root. */
    assert_eval_eq("SquaredEuclideanDistance[{0, 0}, {1, 1}] < "
                   "SquaredEuclideanDistance[{0, 0}, {2, 2}]", "True", 0);
    assert_eval_eq("EuclideanDistance[{0, 0}, {1, 1}] < "
                   "EuclideanDistance[{0, 0}, {2, 2}]", "True", 0);
}

/* The two spanning-tree builders must agree.
 *
 * Machine-precision points take a double Prim (over two orders of magnitude
 * faster: 2000 2-D points went from 1.49 s to 6.8 ms) while exact points keep the
 * Expr-arithmetic builder, which is the only one that can order a Rational or a
 * bigint correctly. Two code paths for one definition is a standing risk, so the
 * same points are clustered both ways -- as integers, which are machine, and
 * divided by a constant, which makes them exact -- and the partitions must
 * match. */
void test_find_clusters_builder_agreement() {
    assert_eval_eq(
        "FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}, {50, 50}, {2, 1}}, 3] === "
        "3 * FindClusters[{{1, 1}, {1, 2}, {9, 9}, {9, 8}, {50, 50}, {2, 1}}/3, 3]",
        "True", 0);
    assert_eval_eq(
        "Module[{p = Table[{RandomInteger[1000], RandomInteger[1000]}, {400}]},"
        " FindClusters[p, 5] === 2 * FindClusters[p/2, 5]]", "True", 0);

    /* Exactness is not lost on either path: two vectors differing by 1/10^20
     * stay together, which a double projection of the ELEMENTS would collapse.
     * (The distances may be computed in double; distinctness never is.) */
    assert_eval_eq("FindClusters[{{1/3, 1/7}, {1/3, 1/7 + 1/10^20}, {5, 5}}]",
                   "{{{1/3, 1/7}, {1/3, 100000000000000000007/700000000000000000000}}, {{5, 5}}}", 0);

    /* Each builder has its own ceiling, because their constants differ. */
    assert_eval_eq("Head[FindClusters[Partition[RandomReal[{0, 1}, 2 * 20001], 2]]]",
                   "FindClusters", 0);
    assert_eval_eq("Head[FindClusters[Table[{RandomInteger[10^6]/7, RandomInteger[10^6]/11}, {2001}]]]",
                   "FindClusters", 0);
}

int main() {
    symtab_init();
    core_init();
    extern void trig_init(void);
    trig_init();
    
    TEST(test_min);
    TEST(test_max);
    TEST(test_listq);
    TEST(test_vectorq);
    TEST(test_matrixq);
    TEST(test_table_n);
    TEST(test_table_imax);
    TEST(test_table_imin_imax);
    TEST(test_table_imin_imax_di);
    TEST(test_table_list);
    TEST(test_table_nested);
    
    TEST(test_range_imax);
    TEST(test_range_imin_imax);
    TEST(test_range_imin_imax_di);
    TEST(test_range_real);
    TEST(test_range_no_silent_truncation);
    
    TEST(test_array_n);
    TEST(test_array_n_r);
    TEST(test_array_nested);
    
    TEST(test_take);
    TEST(test_drop);
    TEST(test_flatten);
    TEST(test_partition);
    TEST(test_pick);
    TEST(test_rotate);
    TEST(test_reverse);
    TEST(test_transpose);
    TEST(test_total);
    TEST(test_commonest);
    TEST(test_splitby);

    TEST(test_join_basic);
    TEST(test_join_two_lists);
    TEST(test_join_single_list);
    TEST(test_join_empty_lists);
    TEST(test_join_non_list_head);
    TEST(test_join_mismatched_heads);
    TEST(test_join_level2_matrices);
    TEST(test_join_level2_ragged);
    TEST(test_join_level2_ragged_unequal_lengths);

    TEST(test_subsets);
    TEST(test_riffle);
    TEST(test_gather);
    TEST(test_subdivide);
    TEST(test_nearest);
    TEST(test_find_clusters);
    TEST(test_distance_functions);
    TEST(test_find_clusters_builder_agreement);

    printf("All list tests passed!\n");
    return 0;
}


#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
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
    Expr* res = evaluate(t);
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
        {"Commonest[{a, E, Sin[y], E, a, 1.5, 3}, 10]", "{a, E, Sin[y], 1.5, 3}"}
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
    
    TEST(test_array_n);
    TEST(test_array_n_r);
    TEST(test_array_nested);
    
    TEST(test_take);
    TEST(test_drop);
    TEST(test_flatten);
    TEST(test_partition);
    TEST(test_rotate);
    TEST(test_reverse);
    TEST(test_transpose);
    TEST(test_total);
    TEST(test_commonest);

    TEST(test_join_basic);
    TEST(test_join_two_lists);
    TEST(test_join_single_list);
    TEST(test_join_empty_lists);
    TEST(test_join_non_list_head);
    TEST(test_join_mismatched_heads);
    TEST(test_join_level2_matrices);
    TEST(test_join_level2_ragged);
    TEST(test_join_level2_ragged_unequal_lengths);

    printf("All list tests passed!\n");
    return 0;
}


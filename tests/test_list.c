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

    printf("All list tests passed!\n");
    return 0;
}


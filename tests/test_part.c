
#include "part.h"
#include "expr.h"
#include "print.h"
#include "test_utils.h"
#include "core.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"

void test_part_single_index() {
    // Setup: {a, b, c}
    Expr* elems[] = {
        expr_new_symbol("a"),
        expr_new_symbol("b"),
        expr_new_symbol("c")
    };
    Expr* list = expr_new_function(expr_new_symbol("List"), elems, 3);
    
    // Test {a,b,c}[[2]] → b
    Expr* index = expr_new_integer(2);
    Expr* result = expr_part(list, &index, 1);
    
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(result->data.symbol.name, "b");
    
    expr_free(list);
    expr_free(index);
    expr_free(result);
}

void test_part_multiple_indices() {
    // Setup: f[a, b, c, d, e]
    Expr* args[] = {
        expr_new_symbol("a"),
        expr_new_symbol("b"),
        expr_new_symbol("c"),
        expr_new_symbol("d"),
        expr_new_symbol("e")
    };
    Expr* f = expr_new_function(expr_new_symbol("f"), args, 5);
    
    // Create index list {1, 3, 4}
    Expr* indices[] = {
        expr_new_integer(1),
        expr_new_integer(3),
        expr_new_integer(4)
    };
    Expr* index_list = expr_new_function(expr_new_symbol("List"), indices, 3);
    
    // Test f[[{1,3,4}]] → f[a,c,d]
    Expr* part_args[] = {index_list};
    Expr* result = expr_part(f, part_args, 1);
    
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(result->data.function.head->data.symbol.name, "f");
    ASSERT(result->data.function.arg_count == 3);
    ASSERT_STR_EQ(result->data.function.args[0]->data.symbol.name, "a");
    ASSERT_STR_EQ(result->data.function.args[1]->data.symbol.name, "c");
    ASSERT_STR_EQ(result->data.function.args[2]->data.symbol.name, "d");
    
    expr_free(f);
    expr_free(index_list);
    expr_free(result);
}

void test_part_nested_extraction() {
    // Setup: {{a,b,c}, {d,e,f}, {g,h,i}}
    Expr* row1[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c")};
    Expr* row2[] = {expr_new_symbol("d"), expr_new_symbol("e"), expr_new_symbol("f")};
    Expr* row3[] = {expr_new_symbol("g"), expr_new_symbol("h"), expr_new_symbol("i")};
    Expr* rows[] = {
      expr_new_function(expr_new_symbol("List"), row1, 3),
      expr_new_function(expr_new_symbol("List"), row2, 3),
      expr_new_function(expr_new_symbol("List"), row3, 3)
    };
    Expr* mat = expr_new_function(expr_new_symbol("List"), rows, 3);
    
    // Test mat[[2,3]] → f
    Expr* indices1[] = {expr_new_integer(2), expr_new_integer(3)};
    Expr* result1 = expr_part(mat, indices1, 2);
    ASSERT(result1 != NULL);
    ASSERT(result1->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(result1->data.symbol.name, "f");
    
    // Test mat[[All,2]] → {b,e,h}
    Expr* all = expr_new_symbol("All");
    Expr* indices2[] = {all, expr_new_integer(2)};
    Expr* result2 = expr_part(mat, indices2, 2);
    ASSERT(result2 != NULL);
    ASSERT(result2->type == EXPR_FUNCTION);
    ASSERT(result2->data.function.arg_count == 3);
    ASSERT_STR_EQ(result2->data.function.args[0]->data.symbol.name, "b");
    ASSERT_STR_EQ(result2->data.function.args[1]->data.symbol.name, "e");
    ASSERT_STR_EQ(result2->data.function.args[2]->data.symbol.name, "h");
    
    expr_free(mat);
    expr_free(result1);
    expr_free(result2);
    expr_free(all);
    expr_free(indices1[0]);
    expr_free(indices1[1]);
    expr_free(indices2[1]);
}


void test_part_all_specifier() {
    /* Create {{1, 2}, {3, 4}, {5, 6}} */
    Expr* row1[] = {expr_new_integer(1), expr_new_integer(2)};
    Expr* row2[] = {expr_new_integer(3), expr_new_integer(4)};
    Expr* row3[] = {expr_new_integer(5), expr_new_integer(6)};
    Expr* rows[] = {
        expr_new_function(expr_new_symbol("List"), row1, 2),
        expr_new_function(expr_new_symbol("List"), row2, 2),
        expr_new_function(expr_new_symbol("List"), row3, 2)
    };
    Expr* matrix = expr_new_function(expr_new_symbol("List"), rows, 3);

    /* Test matrix[[All, 2]] → {2, 4, 6} */
    Expr* all = expr_new_symbol("All");
    Expr* col2 = expr_new_integer(2);
    Expr* indices[] = {all, col2};
    Expr* result = expr_part(matrix, indices, 2);

    /* Verify result */
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_FUNCTION);
    ASSERT(strcmp(result->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(result->data.function.arg_count == 3);
    
    Expr* expected[] = {expr_new_integer(2), expr_new_integer(4), expr_new_integer(6)};
    for (int i = 0; i < 3; i++) {
        ASSERT(result->data.function.args[i]->type == EXPR_INTEGER);
        ASSERT(result->data.function.args[i]->data.integer == expected[i]->data.integer);
        expr_free(expected[i]);
    }

    /* Cleanup */
    expr_free(matrix);
    expr_free(all);
    expr_free(col2);
    expr_free(result);
}



void test_part_negative_index() {
    // Setup: {a, b, c}
    Expr* elems[] = {
        expr_new_symbol("a"),
        expr_new_symbol("b"),
        expr_new_symbol("c")
    };
    Expr* list = expr_new_function(expr_new_symbol("List"), elems, 3);
    
    // Test {a,b,c}[[-1]] → c
    Expr* index = expr_new_integer(-1);
    Expr* result = expr_part(list, &index, 1);
    
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(result->data.symbol.name, "c");
    
    expr_free(list);
    expr_free(index);
    expr_free(result);
}

void test_part_head_extraction() {
    // Setup: f[a, b]
    Expr* elems[] = {
        expr_new_symbol("a"),
        expr_new_symbol("b")
    };
    Expr* f = expr_new_function(expr_new_symbol("f"), elems, 2);
    
    // Test f[a,b][[0]] → f
    Expr* index0 = expr_new_integer(0);
    Expr* result = expr_part(f, &index0, 1);
    
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(result->data.symbol.name, "f");
    
    // Test 42[[0]] → Integer
    Expr* num = expr_new_integer(42);
    Expr* result2 = expr_part(num, &index0, 1);
    ASSERT(result2 != NULL);
    ASSERT(result2->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(result2->data.symbol.name, "Integer");

    expr_free(f);
    expr_free(index0);
    expr_free(result);
    expr_free(num);
    expr_free(result2);
}

void test_part_atomic() {
    Expr* t1 = parse_expression("Complex[2, 3][[1]]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->type == EXPR_FUNCTION); // Should remain unevaluated since part 1 doesn't exist
    ASSERT_STR_EQ(res1->data.function.head->data.symbol.name, "Part");
    expr_free(t1); expr_free(res1);

    Expr* t2 = parse_expression("Rational[2, 3][[1]]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res2->data.function.head->data.symbol.name, "Part");
    expr_free(t2); expr_free(res2);

    Expr* t3 = parse_expression("Complex[2, 3][[0]]");
    Expr* res3 = evaluate(t3);
    ASSERT(res3->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res3->data.symbol.name, "Complex");
    expr_free(t3); expr_free(res3);
}

void test_head() {
    Expr* f_head = expr_new_symbol("f");
    Expr* f_args[] = {expr_new_symbol("a"), expr_new_symbol("b")};
    Expr* f = expr_new_function(f_head, f_args, 2);
    
    // Head[f[a,b]] -> f
    Expr* res1 = expr_head(f);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "f") == 0);
    expr_free(res1);

    // Head[1] -> Integer
    Expr* i = expr_new_integer(1);
    Expr* res2 = expr_head(i);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "Integer") == 0);
    expr_free(res2);

    // Head[1.1] -> Real
    Expr* r = expr_new_real(1.1);
    Expr* res3 = expr_head(r);
    ASSERT(res3 && res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol.name, "Real") == 0);
    expr_free(res3);

    // Head[x] -> Symbol
    Expr* s = expr_new_symbol("x");
    Expr* res4 = expr_head(s);
    ASSERT(res4 && res4->type == EXPR_SYMBOL && strcmp(res4->data.symbol.name, "Symbol") == 0);
    expr_free(res4);

    // Head["s"] -> String
    Expr* str = expr_new_string("s");
    Expr* res5 = expr_head(str);
    ASSERT(res5 && res5->type == EXPR_SYMBOL && strcmp(res5->data.symbol.name, "String") == 0);
    expr_free(res5);

    expr_free(f);
    expr_free(i);
    expr_free(r);
    expr_free(s);
    expr_free(str);
}

void test_first_last_most_rest() {
    // Setup: f[a, b, c]
    Expr* args[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c")};
    Expr* f = expr_new_function(expr_new_symbol("f"), args, 3);
    
    // First[f[a,b,c]] -> a
    Expr* copy1 = expr_copy(f);
    Expr* wrap_f1 = expr_new_function(expr_new_symbol("First"), &copy1, 1);
    Expr* res1 = builtin_first(wrap_f1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "a") == 0);
    expr_free(res1); expr_free(wrap_f1);

    // Last[f[a,b,c]] -> c
    Expr* copy2 = expr_copy(f);
    Expr* wrap_f2 = expr_new_function(expr_new_symbol("Last"), &copy2, 1);
    Expr* res2 = builtin_last(wrap_f2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "c") == 0);
    expr_free(res2); expr_free(wrap_f2);

    // Most[f[a,b,c]] -> f[a, b]
    Expr* copy3 = expr_copy(f);
    Expr* wrap_f3 = expr_new_function(expr_new_symbol("Most"), &copy3, 1);
    Expr* res3 = builtin_most(wrap_f3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && res3->data.function.arg_count == 2);
    ASSERT(strcmp(res3->data.function.args[0]->data.symbol.name, "a") == 0);
    ASSERT(strcmp(res3->data.function.args[1]->data.symbol.name, "b") == 0);
    expr_free(res3); expr_free(wrap_f3);

    // Rest[f[a,b,c]] -> f[b, c]
    Expr* copy4 = expr_copy(f);
    Expr* wrap_f4 = expr_new_function(expr_new_symbol("Rest"), &copy4, 1);
    Expr* res4 = builtin_rest(wrap_f4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && res4->data.function.arg_count == 2);
    ASSERT(strcmp(res4->data.function.args[0]->data.symbol.name, "b") == 0);
    ASSERT(strcmp(res4->data.function.args[1]->data.symbol.name, "c") == 0);
    expr_free(res4); expr_free(wrap_f4);

    expr_free(f);
}


void test_insert() {
    // Setup: {a, b, c}
    Expr* args[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c")};
    Expr* list = expr_new_function(expr_new_symbol("List"), args, 3);
    Expr* x = expr_new_symbol("x");

    // Insert[{a, b, c}, x, 2] -> {a, x, b, c}
    Expr* n2 = expr_new_integer(2);
    Expr* res1 = expr_insert(list, x, n2);
    ASSERT(res1 && res1->data.function.arg_count == 4);
    ASSERT_STR_EQ(res1->data.function.args[1]->data.symbol.name, "x");
    expr_free(res1); expr_free(n2);

    // Insert[{a, b, c}, x, -1] -> {a, b, c, x}
    Expr* nm1 = expr_new_integer(-1);
    Expr* res2 = expr_insert(list, x, nm1);
    ASSERT(res2 && res2->data.function.arg_count == 4);
    ASSERT_STR_EQ(res2->data.function.args[3]->data.symbol.name, "x");
    expr_free(res2); expr_free(nm1);

    // Insert[{{a, b}, {c, d}}, x, {1, 2}] -> {{a, x, b}, {c, d}}
    Expr* sub1_args[] = {expr_new_symbol("a"), expr_new_symbol("b")};
    Expr* sub1 = expr_new_function(expr_new_symbol("List"), sub1_args, 2);
    Expr* sub2_args[] = {expr_new_symbol("c"), expr_new_symbol("d")};
    Expr* sub2 = expr_new_function(expr_new_symbol("List"), sub2_args, 2);
    Expr* mat = expr_new_function(expr_new_symbol("List"), (Expr*[]){sub1, sub2}, 2);
    
    Expr* path_args[] = {expr_new_integer(1), expr_new_integer(2)};
    Expr* path = expr_new_function(expr_new_symbol("List"), path_args, 2);
    Expr* res3 = expr_insert(mat, x, path);
    ASSERT(res3 && res3->data.function.args[0]->data.function.arg_count == 3);
    ASSERT_STR_EQ(res3->data.function.args[0]->data.function.args[1]->data.symbol.name, "x");
    expr_free(res3); expr_free(path); expr_free(mat);

    // Insert[{a, b, c, d}, x, {{1}, {3}}] -> {x, a, b, x, c, d}
    Expr* long_args[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c"), expr_new_symbol("d")};
    Expr* list_long = expr_new_function(expr_new_symbol("List"), long_args, 4);
    Expr* p1 = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_new_integer(1)}, 1);
    Expr* p3 = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_new_integer(3)}, 1);
    Expr* multi_pos = expr_new_function(expr_new_symbol("List"), (Expr*[]){p1, p3}, 2);
    Expr* res4 = expr_insert(list_long, x, multi_pos);
    ASSERT(res4 && res4->data.function.arg_count == 6);
    ASSERT_STR_EQ(res4->data.function.args[0]->data.symbol.name, "x");
    ASSERT_STR_EQ(res4->data.function.args[3]->data.symbol.name, "x");
    expr_free(res4); expr_free(multi_pos); expr_free(list_long);

    expr_free(list);
    expr_free(x);
}

void test_delete() {
    Expr* args[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c"), expr_new_symbol("d")};
    Expr* list = expr_new_function(expr_new_symbol("List"), args, 4);

    // Delete[{a, b, c, d}, 2] -> {a, c, d}
    Expr* n2 = expr_new_integer(2);
    Expr* res1 = expr_delete(list, n2);
    ASSERT(res1 && res1->data.function.arg_count == 3);
    ASSERT_STR_EQ(res1->data.function.args[1]->data.symbol.name, "c");
    expr_free(res1); expr_free(n2);

    // Delete[{a, b, c, d}, -1] -> {a, b, c}
    Expr* nm1 = expr_new_integer(-1);
    Expr* res2 = expr_delete(list, nm1);
    ASSERT(res2 && res2->data.function.arg_count == 3);
    ASSERT_STR_EQ(res2->data.function.args[2]->data.symbol.name, "c");
    expr_free(res2); expr_free(nm1);

    // Delete[{{a, b}, {c, d}}, {1, 2}] -> {{a}, {c, d}}
    Expr* sub1_args[] = {expr_new_symbol("a"), expr_new_symbol("b")};
    Expr* sub1 = expr_new_function(expr_new_symbol("List"), sub1_args, 2);
    Expr* sub2_args[] = {expr_new_symbol("c"), expr_new_symbol("d")};
    Expr* sub2 = expr_new_function(expr_new_symbol("List"), sub2_args, 2);
    Expr* mat = expr_new_function(expr_new_symbol("List"), (Expr*[]){sub1, sub2}, 2);
    
    Expr* path_args[] = {expr_new_integer(1), expr_new_integer(2)};
    Expr* path = expr_new_function(expr_new_symbol("List"), path_args, 2);
    Expr* res3 = expr_delete(mat, path);
    ASSERT(res3 && res3->data.function.args[0]->data.function.arg_count == 1);
    ASSERT_STR_EQ(res3->data.function.args[0]->data.function.args[0]->data.symbol.name, "a");
    expr_free(res3); expr_free(path); expr_free(mat);

    // Delete[{a, b, c, d}, {{1}, {3}}] -> {b, d}
    Expr* long_args[] = {expr_new_symbol("a"), expr_new_symbol("b"), expr_new_symbol("c"), expr_new_symbol("d")};
    Expr* list_long = expr_new_function(expr_new_symbol("List"), long_args, 4);
    Expr* p1 = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_new_integer(1)}, 1);
    Expr* p3 = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_new_integer(3)}, 1);
    Expr* multi_pos = expr_new_function(expr_new_symbol("List"), (Expr*[]){p1, p3}, 2);
    Expr* res4 = expr_delete(list_long, multi_pos);
    ASSERT(res4 && res4->data.function.arg_count == 2);
    ASSERT_STR_EQ(res4->data.function.args[0]->data.symbol.name, "b");
    ASSERT_STR_EQ(res4->data.function.args[1]->data.symbol.name, "d");
    expr_free(res4); expr_free(multi_pos); expr_free(list_long);

    expr_free(list);
}

void test_replace_part() {
    // Basic test: ReplacePart[{a, b, c}, 2 -> x]
    Expr* t1 = parse_expression("ReplacePart[{a, b, c}, 2 -> x]");
    Expr* res1 = evaluate(t1);
    ASSERT(res1->type == EXPR_FUNCTION);
    ASSERT(res1->data.function.arg_count == 3);
    ASSERT_STR_EQ(res1->data.function.args[1]->data.symbol.name, "x");
    expr_free(t1); expr_free(res1);

    // Negative index: ReplacePart[{a, b, c}, -1 -> x]
    Expr* t2 = parse_expression("ReplacePart[{a, b, c}, -1 -> x]");
    Expr* res2 = evaluate(t2);
    ASSERT(res2->type == EXPR_FUNCTION);
    ASSERT(res2->data.function.arg_count == 3);
    ASSERT_STR_EQ(res2->data.function.args[2]->data.symbol.name, "x");
    expr_free(t2); expr_free(res2);

    // Multiple paths: ReplacePart[{{a, b}, {c, d}}, {{1, 2}, {2, 1}} -> x]
    Expr* t3 = parse_expression("ReplacePart[{{a, b}, {c, d}}, {{1, 2}, {2, 1}} -> x]");
    Expr* res3 = evaluate(t3);
    ASSERT(res3->type == EXPR_FUNCTION);
    ASSERT(res3->data.function.arg_count == 2);
    ASSERT_STR_EQ(res3->data.function.args[0]->data.function.args[1]->data.symbol.name, "x");
    ASSERT_STR_EQ(res3->data.function.args[1]->data.function.args[0]->data.symbol.name, "x");
    expr_free(t3); expr_free(res3);

    // Multiple rules: ReplacePart[{a, b, c}, {1 -> x, 3 -> y}]
    Expr* t4 = parse_expression("ReplacePart[{a, b, c}, {1 -> x, 3 -> y}]");
    Expr* res4 = evaluate(t4);
    ASSERT(res4->type == EXPR_FUNCTION);
    ASSERT(res4->data.function.arg_count == 3);
    ASSERT_STR_EQ(res4->data.function.args[0]->data.symbol.name, "x");
    ASSERT_STR_EQ(res4->data.function.args[2]->data.symbol.name, "y");
    expr_free(t4); expr_free(res4);

    // Replace nested head: ReplacePart[f[a], 0 -> g]
    Expr* t5 = parse_expression("ReplacePart[f[a], 0 -> g]");
    Expr* res5 = evaluate(t5);
    ASSERT(res5->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res5->data.function.head->data.symbol.name, "g");
    expr_free(t5); expr_free(res5);
}

void test_span_parsing() {
    Expr* e1 = parse_expression("1;;2");
    char* s1 = expr_to_string_fullform(e1);
    ASSERT_STR_EQ(s1, "Span[1, 2]");
    free(s1); expr_free(e1);

    Expr* e2 = parse_expression("1;;");
    char* s2 = expr_to_string_fullform(e2);
    ASSERT_STR_EQ(s2, "Span[1, All]");
    free(s2); expr_free(e2);

    Expr* e3 = parse_expression(";;2");
    char* s3 = expr_to_string_fullform(e3);
    ASSERT_STR_EQ(s3, "Span[1, 2]");
    free(s3); expr_free(e3);

    Expr* e4 = parse_expression(";;");
    char* s4 = expr_to_string_fullform(e4);
    ASSERT_STR_EQ(s4, "Span[1, All]");
    free(s4); expr_free(e4);

    Expr* e5 = parse_expression("1;;2;;3");
    char* s5 = expr_to_string_fullform(e5);
    ASSERT_STR_EQ(s5, "Span[1, 2, 3]");
    free(s5); expr_free(e5);

    Expr* e6 = parse_expression(";;;;3");
    char* s6 = expr_to_string_fullform(e6);
    ASSERT_STR_EQ(s6, "Span[1, All, 3]");
    free(s6); expr_free(e6);

    Expr* e7 = parse_expression("1;;;;3");
    char* s7 = expr_to_string_fullform(e7);
    ASSERT_STR_EQ(s7, "Span[1, All, 3]");
    free(s7); expr_free(e7);
    
    Expr* e8 = parse_expression(";;2;;3");
    char* s8 = expr_to_string_fullform(e8);
    ASSERT_STR_EQ(s8, "Span[1, 2, 3]");
    free(s8); expr_free(e8);
    
    Expr* e9 = parse_expression("a;;b");
    char* s9 = expr_to_string_fullform(e9);
    ASSERT_STR_EQ(s9, "Span[a, b]");
    free(s9); expr_free(e9);
    
    Expr* e10 = parse_expression("1;;-1;;2");
    char* s10 = expr_to_string_fullform(e10);
    ASSERT_STR_EQ(s10, "Span[1, -1, 2]");
    free(s10); expr_free(e10);
}

void test_span_assignment() {
    Expr* st1 = parse_expression("t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = x; t");
    Expr* rt1 = evaluate(st1);
    char* s_rt1 = expr_to_string_fullform(rt1);
    ASSERT_STR_EQ(s_rt1, "List[a, x, x, x, x, f, g, h]");
    free(s_rt1); expr_free(rt1); expr_free(st1);

    Expr* st2 = parse_expression("t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = {p, q, r, s}; t");
    Expr* rt2 = evaluate(st2);
    char* s_rt2 = expr_to_string_fullform(rt2);
    ASSERT_STR_EQ(s_rt2, "List[a, p, q, r, s, f, g, h]");
    free(s_rt2); expr_free(rt2); expr_free(st2);
}

void test_extract() {
    assert_eval_eq("Extract[{1,2,3,4},{2}]", "2", 0);
    assert_eval_eq("Extract[{2}][{a,b,c,d}]", "b", 0);
    assert_eval_eq("Extract[f[g[1,2],h[3]],{1,2}]", "2", 0);
    assert_eval_eq("Extract[f[g[1,2],h[x^2]],{{1,2},{2,1,1}}]", "{2, x}", 0);
    assert_eval_eq("Extract[{{a,b,c},{d,e,f},{g,h,i}}, {All,2}]", "{b, e, h}", 0);
    assert_eval_eq("e = f[g[1,2],{h[3]}]; p = Position[e,_Integer]; Extract[e,p]", "{1, 2, 3}", 0);
    assert_eval_eq("Clear[e, p]", "Null", 0);
    assert_eval_eq("Extract[{a,b,c,d,e}, {3}]", "c", 0);
    assert_eval_eq("Extract[{a,b,c,d,e}, {{1},{4},{3}}]", "{a, d, c}", 0);
    assert_eval_eq("mat=Array[a,{3,3}]; Extract[mat,{1,3}]", "a[1, 3]", 0);
    assert_eval_eq("mat=Array[a,{3,3}]; Extract[mat,{{1},{3}}]", "{{a[1, 1], a[1, 2], a[1, 3]}, {a[3, 1], a[3, 2], a[3, 3]}}", 0);
    assert_eval_eq("mat=Array[a,{3,3}]; Extract[mat,{{1,3}}]", "{a[1, 3]}", 0);
    assert_eval_eq("mat=Array[a,{3,3}]; Extract[mat,{{1,2},{3,3},{2,1}}]", "{a[1, 2], a[3, 3], a[2, 1]}", 0);
    assert_eval_eq("Extract[{a,b,c,d,e}, {{1},{2},{5}},Hold]", "{Hold[a], Hold[b], Hold[e]}", 0);
    assert_eval_eq("rules={{a->1,b->2},{c->3,d->4,e->5},{f->6}}; Extract[rules,{All,1}]", "{a -> 1, c -> 3, f -> 6}", 0);
    assert_eval_eq("rules={{a->1,b->2},{c->3,d->4,e->5},{f->6}}; Extract[rules,{2,All,2}]", "{3, 4, 5}", 0);
    assert_eval_eq("rules={{a->1,b->2},{c->3,d->4,e->5},{f->6}}; Extract[rules,{All,All,1}]", "{{a, b}, {c, d, e}, {f}}", 0);
    assert_eval_eq("Extract[e[f[1,2,3],g[4,5,6],h[7,8,9]],{All,1}]", "e[1, 4, 7]", 0);
    assert_eval_eq("Extract[Range[10],{3;;5}]", "{3, 4, 5}", 0);
    assert_eval_eq("Extract[{a,b,c,d,e,f},{1;;-1;;2}]", "{a, c, e}", 0);
    assert_eval_eq("Extract[f[1,2,3,4,5,6,7,8,9,10],{-1;;1;;-3}]", "f[10, 7, 4, 1]", 0);
    assert_eval_eq("Extract[{{1,2},{3,4,5}},{All,1;;2}]", "{{1, 2}, {3, 4}}", 0);
    assert_eval_eq("Extract[{{1,2},{3,4,5}},{All,{1,-1}}]", "{{1, 2}, {3, 5}}", 0);
    assert_eval_eq("Extract[{a:>1^2, b:>2^2,c:>3^3},{1;;2,2},Hold]", "Hold[{1^2, 2^2}]", 0);
    /* Scalar position is shorthand for a length-1 path. */
    assert_eval_eq("Extract[{1,2,3}, 0, f]", "f[List]", 0);
    assert_eval_eq("Extract[{1,2,3}, 0]", "List", 0);
    assert_eval_eq("Extract[{1,2,3}, 1]", "1", 0);
    assert_eval_eq("Extract[{1,2,3}, -1]", "3", 0);
    assert_eval_eq("Extract[{1,2,3}, 1, f]", "f[1]", 0);
}

// Part-specific test runner
void run_part_tests() {
    TEST(test_extract);
    TEST(test_span_parsing);
    TEST(test_span_assignment);
    TEST(test_part_single_index);
    TEST(test_part_multiple_indices);
    TEST(test_part_nested_extraction);
    TEST(test_part_all_specifier);
    TEST(test_part_negative_index);
    TEST(test_part_head_extraction);
    TEST(test_part_atomic);
    TEST(test_head);
    TEST(test_first_last_most_rest);
    TEST(test_insert);
    TEST(test_delete);
    TEST(test_replace_part);
}

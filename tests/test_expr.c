
#include "expr.h"
#include "print.h"
#include "test_utils.h"
#include "test_part.h"
#include <string.h>

void test_integer_expr() {
    Expr* e = expr_new_integer(42);
    ASSERT(e != NULL);
    ASSERT(e->type == EXPR_INTEGER);
    ASSERT(e->data.integer == 42);
    expr_free(e);
}

void test_real_expr() {
    Expr* e = expr_new_real(3.14159);
    ASSERT(e != NULL);
    ASSERT(e->type == EXPR_REAL);
    ASSERT(e->data.real == 3.14159);
    expr_free(e);
}

void test_symbol_expr() {
    Expr* e = expr_new_symbol("x");
    ASSERT(e != NULL);
    ASSERT(e->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(e->data.symbol.name, "x");
    expr_free(e);
}

void test_string_expr() {
    Expr* e = expr_new_string("hello");
    ASSERT(e != NULL);
    ASSERT(e->type == EXPR_STRING);
    ASSERT_STR_EQ(e->data.string, "hello");
    expr_free(e);
}

void test_function_expr() {
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {x};
    Expr* e = expr_new_function(expr_new_symbol("Sin"), args, 1);
    
    ASSERT(e != NULL);
    ASSERT(e->type == EXPR_FUNCTION);
    ASSERT(e->data.function.arg_count == 1);
    ASSERT(e->data.function.args[0] == x);
    ASSERT_STR_EQ(e->data.function.head->data.symbol.name, "Sin");
    
    expr_free(e); // Recursively frees x
}

void test_nested_free() {
    // Build Sin[Cos[x]]
    Expr* x = expr_new_symbol("x");
    Expr* cos_args[] = {x};
    Expr* cos = expr_new_function(expr_new_symbol("Cos"), cos_args, 1);
    Expr* sin_args[] = {cos};
    Expr* e = expr_new_function(expr_new_symbol("Sin"), sin_args, 1);
    
    expr_free(e); // Should free all nested expressions
    // Test passes if no memory leaks
}

void test_memory_cleanup() {
    // Verify strings are copied
    char buffer[] = "temp";
    Expr* e = expr_new_string(buffer);
    buffer[0] = 'X'; // Modify original
    ASSERT_STR_EQ(e->data.string, "temp"); // Should be unchanged
    expr_free(e);
}

void test_print_basic() {
    Expr* f = expr_new_function(
        expr_new_symbol("f"),
        (Expr*[]){expr_new_symbol("a"), expr_new_symbol("b")},
        2
    );
    
    printf("Expected: f[a, b]\nActual:   ");
    expr_print(f);
    expr_free(f);
}

void test_print_nested() {
    Expr* nested = expr_new_function(
        expr_new_symbol("g"),
        (Expr*[]){
            expr_new_integer(1),
            expr_new_function(
                expr_new_symbol("h"),
                (Expr*[]){expr_new_symbol("x")},
                1
            )
        },
        2
    );
    
    printf("Expected: g[1, h[x]]\nActual:   ");
    expr_print(nested);
    expr_free(nested);
}

#include "core.h"
#include "symtab.h"

int main() {
    symtab_init();
    core_init();
    TEST(test_integer_expr);
    TEST(test_real_expr);
    TEST(test_symbol_expr);
    TEST(test_string_expr);
    TEST(test_function_expr);
    TEST(test_nested_free);
    TEST(test_memory_cleanup);
    TEST(test_print_basic);
    TEST(test_print_nested);
    
    run_part_tests();
    
    printf("All tests passed!\n");
    return 0;
}

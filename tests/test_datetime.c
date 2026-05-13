#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

void run_test(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        assert(false);
    }
    Expr* res = evaluate(e);
    if (!res) {
        printf("Failed to evaluate: %s\n", input);
        expr_free(e);
        assert(false);
    }
    
    char* s = expr_to_string_fullform(res);
    printf("EVAL: %s -> %s\n", input, s);
    free(s);
    expr_free(e);
    expr_free(res);
}

void test_timing() {
    // Timing should return a list: {time, result}
    // E.g. Timing[Plus[2, 3]] -> {time, 5}
    Expr* p = parse_expression("Timing[2 + 3]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol, "List") == 0);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL); // time
    assert(e->data.function.args[1]->type == EXPR_INTEGER);
    assert(e->data.function.args[1]->data.integer == 5);
    expr_free(e);
}

void test_repeated_timing() {
    Expr* p = parse_expression("RepeatedTiming[2 + 3]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_FUNCTION);
    assert(strcmp(e->data.function.head->data.symbol, "List") == 0);
    assert(e->data.function.arg_count == 2);
    assert(e->data.function.args[0]->type == EXPR_REAL); // time
    assert(e->data.function.args[1]->type == EXPR_INTEGER);
    assert(e->data.function.args[1]->data.integer == 5);
    expr_free(e);
}

static int64_t eval_to_int(const char* input) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_INTEGER);
    int64_t v = e->data.integer;
    expr_free(e);
    return v;
}

static double eval_to_real(const char* input) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    expr_free(p);
    double v;
    if (e->type == EXPR_REAL)         v = e->data.real;
    else if (e->type == EXPR_INTEGER) v = (double)e->data.integer;
    else { assert(0); v = 0.0; }
    expr_free(e);
    return v;
}

void test_absolute_time_datelist_int() {
    /* Exact Mathematica reference values. */
    assert(eval_to_int("AbsoluteTime[{2022,1,1,0,0,0}]") == 3849984000LL);
    assert(eval_to_int("AbsoluteTime[{2022,1}]")          == 3849984000LL);
    assert(eval_to_int("AbsoluteTime[{2022}]")            == 3849984000LL);
}

void test_absolute_time_normalization() {
    /* {2022, 2, 31} normalizes to {2022, 3, 3}. */
    int64_t a = eval_to_int("AbsoluteTime[{2022,2,31}]");
    int64_t b = eval_to_int("AbsoluteTime[{2022,3,3}]");
    assert(a == 3855254400LL);
    assert(a == b);
}

void test_absolute_time_fractional_day() {
    /* {2022, 3, 15.5} = {2022, 3, 15} + half a day (43200 s). */
    double v = eval_to_real("AbsoluteTime[{2022,3,15.5}]");
    int64_t base = eval_to_int("AbsoluteTime[{2022,3,15,0,0,0}]");
    assert(v == (double)base + 43200.0);
}

void test_absolute_time_fractional_hour() {
    /* {2022, 3, 15, 12.3} = {2022, 3, 15} + 12.3 hours. */
    double v = eval_to_real("AbsoluteTime[{2022,3,15,12.3}]");
    int64_t base = eval_to_int("AbsoluteTime[{2022,3,15,0,0,0}]");
    double expected = (double)base + 12.3 * 3600.0;
    assert(v == expected);
}

void test_absolute_time_passthrough() {
    /* AbsoluteTime[t] for numeric t returns t unchanged. */
    assert(eval_to_int("AbsoluteTime[3849984000]") == 3849984000LL);
}

void test_absolute_time_now() {
    /* AbsoluteTime[] returns a real and stays in a plausible 21st-century range. */
    Expr* p = parse_expression("AbsoluteTime[]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_REAL);
    /* 2020-01-01 .. 2200-01-01 in AbsoluteTime seconds. */
    assert(e->data.real > 3.78e9);
    assert(e->data.real < 9.46e9);
    expr_free(e);
}

void test_absolute_time_attributes() {
    /* Must be Protected per the spec. */
    Expr* p = parse_expression("MemberQ[Attributes[AbsoluteTime], Protected]");
    Expr* e = evaluate(p);
    expr_free(p);
    assert(e->type == EXPR_SYMBOL);
    assert(strcmp(e->data.symbol, "True") == 0);
    expr_free(e);
}

int main() {
    symtab_init();
    core_init();

    printf("Running datetime tests...\n");
    test_timing();
    test_repeated_timing();
    test_absolute_time_datelist_int();
    test_absolute_time_normalization();
    test_absolute_time_fractional_day();
    test_absolute_time_fractional_hour();
    test_absolute_time_passthrough();
    test_absolute_time_now();
    test_absolute_time_attributes();
    printf("All datetime tests passed!\n");
    symtab_clear();
    return 0;
}

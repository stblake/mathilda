#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "print_latex.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void test_print_basic() {
    // We can't easily capture stdout here without pipe/dup, 
    // but we can check if it returns Null.
    Expr* e = parse_expression("Print[1, \" \", x + y]");
    Expr* res = evaluate(e);
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.symbol.name, "Null") == 0);
    expr_free(e);
    expr_free(res);
}

void test_fullform_wrapper() {
    Expr* e = parse_expression("Print[FullForm[x + y]]");
    Expr* res = evaluate(e);
    ASSERT(res != NULL);
    ASSERT(strcmp(res->data.symbol.name, "Null") == 0);
    expr_free(e);
    expr_free(res);
}

void test_inputform_wrapper() {
    Expr* e = parse_expression("Print[InputForm[x + y]]");
    Expr* res = evaluate(e);
    ASSERT(res != NULL);
    ASSERT(strcmp(res->data.symbol.name, "Null") == 0);
    expr_free(e);
    expr_free(res);
}


/* Regression: a negative bigint term in a Plus must print as " - <abs>"
 * rather than " + -<abs>". This previously affected only bignum coefficients
 * (int64 terms were already handled). */
void test_negative_bigint_in_plus() {
    /* Top-level negative bigint term: x - Fibonacci[100] y. */
    Expr* e = parse_expression("x - Fibonacci[100] y");
    Expr* res = evaluate(e);
    char* str = expr_to_string(res);
    ASSERT(strstr(str, "+ -") == NULL);
    ASSERT(strstr(str, "- 354224848179261915075 y") != NULL);
    free(str);
    expr_free(e);
    expr_free(res);

    /* Bare negative bigint plus a symbol. */
    Expr* e2 = parse_expression("-Fibonacci[100] + x");
    Expr* res2 = evaluate(e2);
    char* str2 = expr_to_string(res2);
    ASSERT(strstr(str2, "+ -") == NULL);
    ASSERT(strncmp(str2, "-354224848179261915075", 22) == 0);
    free(str2);
    expr_free(e2);
    expr_free(res2);
}

void test_holdform() {
    Expr* e = parse_expression("HoldForm[1 + 1]");
    Expr* res = evaluate(e);
    
    char* str = expr_to_string(res);
    ASSERT(strcmp(str, "1 + 1") == 0);
    free(str);
    
    char* full = expr_to_string_fullform(res);
    ASSERT(strcmp(full, "HoldForm[Plus[1, 1]]") == 0);
    free(full);

    /* HoldForm must be transparent in LaTeX too (notebook rendering), not leak
     * "HoldForm[...]". This is what makes Trace[]'s HoldForm-wrapped steps
     * render cleanly in the notebook. */
    char* tex = expr_to_latex(res);
    ASSERT(strcmp(tex, "1+1") == 0);
    free(tex);

    expr_free(e);
    expr_free(res);
}

void test_series_latex() {
    /* SeriesData must render as the series it represents in LaTeX (notebook),
     * not as the raw SeriesData[...] container — matching the plain printer.
     * Regression for gh #22. Also checks negative rational coefficients show
     * as subtraction (x - 1/6 x^3), not "+ -1/6". */
    Expr* e = parse_expression("Series[Sin[x], {x, 0, 5}]");
    Expr* res = evaluate(e);

    char* str = expr_to_string(res);
    ASSERT(strcmp(str, "x - 1/6 x^3 + 1/120 x^5 + O[x]^6") == 0);
    free(str);

    char* tex = expr_to_latex(res);
    ASSERT(strcmp(tex, "x-\\frac{1}{6}\\,x^{3}+\\frac{1}{120}\\,x^{5}+O[x]^{6}") == 0);
    free(tex);

    expr_free(e);
    expr_free(res);
}

void test_operator_latex() {
    /* Every infix head used to fall through the LaTeX printer's generic
     * `Head[a, b]` arm, so the notebook typeset a DSolve answer as
     * `Rule[y[x], ...]` -- FullForm dressed up as mathematics, and the most
     * visible output in the application, since every Solve/DSolve/Reduce result
     * is built from these heads.
     *
     * Parsed and NOT evaluated: the printer is what is under test, and
     * evaluation would answer the comparisons instead of printing them. */
    static const struct { const char* in; const char* tex; } cases[] = {
        {"a -> b",        "a\\to b"},
        /* A delayed rule is a different object; one glyph for both would make
           `a -> b` and `a :> b` indistinguishable once typeset. */
        {"a :> b",        "a:\\to b"},
        {"x == 1",        "x=1"},
        {"x == -1",       "x=-1"},                 /* an atom stays bare: not (-1) */
        {"a != b",        "a\\neq b"},
        {"a <= b",        "a\\leq b"},
        {"a >= b",        "a\\geq b"},
        {"a && b",        "a\\land b"},
        {"a && b || c",   "a\\land b\\lor c"},      /* && binds tighter: no parens */
        {"(a || b) && c", "\\left(a\\lor b\\right)\\land c"},
        {"!a",            "\\neg a"},
        {"1 < x < 2",     "1<x<2"},                 /* the Inequality chain */
        {"{a -> 1}",      "\\{a\\to 1\\}"},
        {"(a -> b)^2",    "(a\\to b)^{2}"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        char* tex = expr_to_latex(e);
        if (!tex || strcmp(tex, cases[i].tex) != 0)
            printf("  %s -> \"%s\", expected \"%s\"\n",
                   cases[i].in, tex ? tex : "(null)", cases[i].tex);
        ASSERT(tex && strcmp(tex, cases[i].tex) == 0);
        free(tex);
        expr_free(e);
    }

    /* True/False/Null are words: math mode would set `False` as a product of
     * five italic variables. */
    Expr* f = parse_expression("False");
    char* tex = expr_to_latex(f);
    ASSERT(tex && strcmp(tex, "\\text{False}") == 0);
    free(tex);
    expr_free(f);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_print_basic);
    TEST(test_fullform_wrapper);
    TEST(test_inputform_wrapper);
    TEST(test_negative_bigint_in_plus);
    TEST(test_holdform);
    TEST(test_series_latex);
    TEST(test_operator_latex);

    printf("All print tests passed!\n");
    return 0;
}

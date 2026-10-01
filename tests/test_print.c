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

void test_printer_directive_has_no_latex() {
    /* InputForm, FullForm, TeXForm and NumberForm ask the printer for a
     * notation that is NOT StandardForm, and typesetting has only StandardForm
     * to offer. expr_to_latex must therefore answer with the empty string, the
     * signal both the notebook and the FFI read as "use the text payload".
     *
     * Before this, the wrapper leaked into the LaTeX -- `InputForm[\frac{1}
     * {2}]` -- and because the notebook prefers `latex` over the payload
     * whenever it is non-empty, `expr // InputForm` typeset the StandardForm
     * the reader had just asked not to see. The directive was, in effect,
     * ignored everywhere except the terminal REPL.
     *
     * Each case is checked BOTH ways: the plain printer still honours the
     * directive (that half always worked and must stay working), and the LaTeX
     * writer declines. The nested cases are the reason the check is a scan of
     * the whole tree rather than a look at the outermost head. */
    static const struct { const char* in; const char* text; } cases[] = {
        {"InputForm[1/2]",        "1/2"},
        {"FullForm[a + b]",       "Plus[a, b]"},
        {"TeXForm[a/b]",          "\\frac{a}{b}"},
        {"NumberForm[1.5, 2]",    "1.5"},
        {"Hold[InputForm[x]]",    "Hold[x]"},      /* nested, not at the top */
        {"{InputForm[1/2], 3}",   "{1/2, 3}"},     /* one directive in a list  */
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e   = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        Expr* res = evaluate(e);

        char* str = expr_to_string(res);
        if (!str || strcmp(str, cases[i].text) != 0)
            printf("  %s -> text \"%s\", expected \"%s\"\n",
                   cases[i].in, str ? str : "(null)", cases[i].text);
        ASSERT(str && strcmp(str, cases[i].text) == 0);
        free(str);

        char* tex = expr_to_latex(res);
        if (!tex || tex[0] != '\0')
            printf("  %s -> latex \"%s\", expected \"\"\n",
                   cases[i].in, tex ? tex : "(null)");
        ASSERT(tex && tex[0] == '\0');
        free(tex);

        expr_free(e);
        expr_free(res);
    }

    /* The scan must not cost an ordinary expression its typesetting: a result
     * with no directive anywhere still renders. */
    Expr* plain = parse_expression("{1/2, Sqrt[x]}");
    char* ptex = expr_to_latex(plain);
    ASSERT(ptex && ptex[0] != '\0');
    free(ptex);
    expr_free(plain);
}

/* A factor inside a \frac slot must keep its parentheses whenever the slot
 * holds more than one factor: they are juxtaposed by implicit multiplication,
 * so a Plus among them that loses its brackets makes the typeset output state a
 * DIFFERENT expression.  Every factor was being rendered at precedence 0
 * ("never parenthesise"), which is right only for a slot's sole occupant, so
 * `a/(b (c + d))` came out as \frac{a}{b c+d}, i.e. (bc+d)/… — a wrong answer
 * in the typeset form, and the shape every third-kind elliptic answer has. */
void test_texform_fraction_parenthesisation() {
    static const struct { const char* in; const char* tex; } cases[] = {
        /* the bug: two or more factors in a slot */
        {"TeXForm[a/(b (c + d))]",   "\\frac{a}{b \\left(c+d\\right)}"},
        {"TeXForm[1/(x (1 + x))]",   "\\frac{1}{x \\left(1+x\\right)}"},
        {"TeXForm[x/((a + b) (c + d))]",
         "\\frac{x}{\\left(a+b\\right) \\left(c+d\\right)}"},
        {"TeXForm[(x + 1) (x + 2)/((x + 3) (x + 4))]",
         "\\frac{\\left(1+x\\right) \\left(2+x\\right)}"
         "{\\left(3+x\\right) \\left(4+x\\right)}"},
        /* a lone factor must NOT gain gratuitous brackets */
        {"TeXForm[1/(1 + x)]",       "\\frac{1}{1+x}"},
        {"TeXForm[(a + b)/(c + d)]", "\\frac{a+b}{c+d}"},
        {"TeXForm[a/b]",             "\\frac{a}{b}"},
        /* several atoms in a slot need none either */
        {"TeXForm[1/(x y z)]",       "\\frac{1}{x y z}"},
        {"TeXForm[(a + b)/(x y)]",   "\\frac{a+b}{x y}"},
        /* a parenthesised base already carried its own brackets */
        {"TeXForm[1/(x^2 (1 + x)^3)]",
         "\\frac{1}{x^{2} \\left(1+x\\right)^{3}}"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        Expr* res = evaluate(e);
        char* str = expr_to_string(res);
        if (!str || strcmp(str, cases[i].tex) != 0)
            printf("  %s\n    got      %s\n    expected %s\n",
                   cases[i].in, str ? str : "(null)", cases[i].tex);
        ASSERT(str && strcmp(str, cases[i].tex) == 0);
        free(str);
        expr_free(res);
    }
    printf("PASS: TeXForm parenthesises juxtaposed factors in a fraction\n");
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
    TEST(test_printer_directive_has_no_latex);
    TEST(test_texform_fraction_parenthesisation);

    printf("All print tests passed!\n");
    return 0;
}

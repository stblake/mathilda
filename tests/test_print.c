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

/* A coefficient of exactly INT64_MIN: its negation does not fit an int64,
 * and the printer used to emit "--9223372036854775808". */
void test_int64_min_coefficient() {
    const char* in[3] = { "-9223372036854775808 x", "x - 9223372036854775808 y",
                          "-9223372036854775808 (1 + x)" };
    const char* out[3] = { "-9223372036854775808 x", "x - 9223372036854775808 y",
                           "-9223372036854775808 (1 + x)" };
    for (int i = 0; i < 3; i++) {
        Expr* e = parse_expression(in[i]);
        Expr* res = evaluate(e);
        char* str = expr_to_string(res);
        ASSERT_STR_EQ(out[i], str);
        free(str);
        expr_free(e);
        expr_free(res);
    }
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

void test_standardform_latex() {
    /* StandardForm is the notebook's default display: Mathilda's own names
     * typeset (Sin[x], not \sin(x)), upright symbols, the language's operator
     * spellings (==, &&), and E/I as the double-struck ⅇ/ⅈ. Fractions, powers
     * and Sqrt stay two-dimensional. TraditionalForm (expr_to_latex) must be
     * untouched by any of it -- test_operator_latex pins that side.
     *
     * Parsed, not evaluated, except where noted, so the printer is what is
     * under test. */
    static const struct { const char* in; const char* tex; } cases[] = {
        {"Sin[x]",        "\\mathrm{Sin}{\\left[\\mathrm{x}\\right]}"},
        {"Log[2, z]",     "\\mathrm{Log}{\\left[2, \\mathrm{z}\\right]}"},
        {"Abs[q]",        "\\mathrm{Abs}{\\left[\\mathrm{q}\\right]}"},
        {"E^x",           "\xe2\x85\x87^{\\mathrm{x}}"},
        {"Sqrt[y]",       "\\sqrt{\\mathrm{y}}"},
        {"x^(1/3)",       "\\mathrm{x}^{\\frac{1}{3}}"},   /* no surd */
        {"x == 1",        "\\mathrm{x}\\mathrel{==}1"},
        {"a && b || c",   "\\mathrm{a}\\mathrel{\\&\\&}\\mathrm{b}\\mathrel{||}\\mathrm{c}"},
        {"!a",            "!\\mathrm{a}"},
        {"a === b",       "\\mathrm{a}\\mathrel{===}\\mathrm{b}"},
        {"a -> b",        "\\mathrm{a}\\to \\mathrm{b}"},
        {"Alpha",         "\\mathrm{Alpha}"},     /* a NAME, not the letter */
        {"Pi",            "\\pi"},
        {"C[1]",          "\\mathrm{C}{\\left[1\\right]}"},
        {"$Version",      "\\mathrm{\\$Version}"},
        {"Binomial[n, k]","\\mathrm{Binomial}{\\left[\\mathrm{n}, \\mathrm{k}\\right]}"},
        /* The display wrappers switch notation for their subtree only. */
        {"{Sin[x], TraditionalForm[Cos[x]]}",
         "\\{\\mathrm{Sin}{\\left[\\mathrm{x}\\right]}, \\cos\\left(x\\right)\\}"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        char* tex = expr_to_latex_form(e, LATEX_FORM_STANDARD);
        if (!tex || strcmp(tex, cases[i].tex) != 0)
            printf("  %s -> \"%s\", expected \"%s\"\n",
                   cases[i].in, tex ? tex : "(null)", cases[i].tex);
        ASSERT(tex && strcmp(tex, cases[i].tex) == 0);
        free(tex);
        expr_free(e);
    }

    /* Complex numbers carry ⅈ (evaluated: Complex[3, 4] is an atom-like head). */
    Expr* c_in = parse_expression("3 + 4 I");
    Expr* c = evaluate(c_in);   /* evaluate borrows its argument */
    expr_free(c_in);
    char* ctex = expr_to_latex_form(c, LATEX_FORM_STANDARD);
    ASSERT(ctex && strcmp(ctex, "3+4\\,\xe2\x85\x88") == 0);
    free(ctex);
    expr_free(c);

    /* An outermost TraditionalForm asks for that notation outright: there is
     * no StandardForm rendering ("" -- the consumer shows the Traditional one),
     * while the TraditionalForm LaTeX is the wrapped expression's. */
    Expr* t = parse_expression("TraditionalForm[Sin[x]]");
    char* st = expr_to_latex_form(t, LATEX_FORM_STANDARD);
    char* tt = expr_to_latex(t);
    ASSERT(st && st[0] == '\0');
    ASSERT(tt && strcmp(tt, "\\sin\\left(x\\right)") == 0);
    free(st); free(tt);
    expr_free(t);

    /* ...and StandardForm inside TraditionalForm switches back. */
    Expr* s = parse_expression("StandardForm[Sin[x]]");
    char* ts = expr_to_latex(s);
    ASSERT(ts && strcmp(ts, "\\mathrm{Sin}{\\left[\\mathrm{x}\\right]}") == 0);
    free(ts);
    expr_free(s);

    /* A printer directive still has no typeset form in either notation. */
    Expr* d = parse_expression("InputForm[Sin[x]]");
    char* dt = expr_to_latex_form(d, LATEX_FORM_STANDARD);
    ASSERT(dt && dt[0] == '\0');
    free(dt);
    expr_free(d);

    /* As text the display wrappers are invisible -- except under InputForm. */
    Expr* p_in = parse_expression("TraditionalForm[Sin[x]^2]");
    Expr* p = evaluate(p_in);   /* evaluate borrows its argument */
    expr_free(p_in);
    char* ps = expr_to_string(p);
    ASSERT(ps && strcmp(ps, "Sin[x]^2") == 0);
    free(ps);
    expr_free(p);
    Expr* q_in = parse_expression("InputForm[Hold[StandardForm[x]]]");
    Expr* q = evaluate(q_in);   /* evaluate borrows its argument */
    expr_free(q_in);
    char* qs = expr_to_string(q);
    ASSERT(qs && strcmp(qs, "Hold[StandardForm[x]]") == 0);
    free(qs);
    expr_free(q);
}

void test_selectable_standardform() {
    /* The notebook's StandardForm is SELECTABLE: every rendered subexpression
     * is wrapped in \htmlData{mp=K}{...} and parts.text[K] is its InputForm,
     * which is what a reader copies after selecting part of a result. */
    Expr* in = parse_expression("(1 + x)/Sqrt[2 + y] + Sin[a b]");
    Expr* e = evaluate(in);              /* evaluate borrows its argument */
    expr_free(in);
    LatexParts parts;
    char* tex = expr_to_latex_selectable(e, &parts);
    ASSERT(tex != NULL);
    ASSERT(parts.count > 0);

    /* Part 0 is the whole result, and equals the plain text payload. */
    char* whole = expr_to_string(e);
    ASSERT(strcmp(parts.text[0], whole) == 0);
    free(whole);

    /* Every subexpression is a part -- including the reciprocal in the
     * denominator, which is drawn as one radical and so must select as one. */
    static const char* want[] = { "Sin[a b]", "a b", "(1 + x)/Sqrt[2 + y]",
                                  "1 + x", "Sqrt[2 + y]", "2 + y", "x", "y" };
    for (size_t w = 0; w < sizeof(want) / sizeof(want[0]); w++) {
        int found = 0;
        for (size_t k = 0; k < parts.count; k++)
            if (strcmp(parts.text[k], want[w]) == 0) { found = 1; break; }
        if (!found) printf("  missing part \"%s\"\n", want[w]);
        ASSERT(found);
    }

    /* One tag per part, numbered 0..count-1. */
    for (size_t k = 0; k < parts.count; k++) {
        char tag[48];
        snprintf(tag, sizeof(tag), "\\htmlData{mp=%zu}{", k);
        ASSERT(strstr(tex, tag) != NULL);
    }
    ASSERT(strstr(tex, "\\htmlData{mp=" ) == tex);   /* the whole is outermost */
    free(tex);
    latex_parts_free(&parts);
    ASSERT(parts.count == 0 && parts.text == NULL);

    /* An outer TraditionalForm has no StandardForm, so nothing is selectable. */
    Expr* t = parse_expression("TraditionalForm[Sin[x]]");
    char* tt = expr_to_latex_selectable(t, &parts);
    ASSERT(tt && tt[0] == '\0' && parts.count == 0);
    free(tt);
    expr_free(t);

    /* Over budget: plain untagged StandardForm, no parts. */
    Expr* big_in = parse_expression("Range[10000]");
    Expr* big = evaluate(big_in);
    expr_free(big_in);
    char* bt = expr_to_latex_selectable(big, &parts);
    ASSERT(bt && !strstr(bt, "htmlData") && parts.count == 0);
    free(bt);
    expr_free(big);
    expr_free(e);
}

void test_radical_body_unbracketed() {
    /* A radical groups its body: Sqrt[1 + x] is drawn as the radical over
     * 1+x, never over a bracketed (1+x). Both notations. */
    Expr* in = parse_expression("{Sqrt[1 + x], (1 + x)^(1/3)}");
    Expr* e = evaluate(in);
    expr_free(in);
    char* tt = expr_to_latex(e);
    ASSERT(tt && strcmp(tt, "\\{\\sqrt{1+x}, \\sqrt[3]{1+x}\\}") == 0);
    char* st = expr_to_latex_form(e, LATEX_FORM_STANDARD);
    ASSERT(st && strstr(st, "\\sqrt{1+\\mathrm{x}}") != NULL);
    free(tt); free(st);
    expr_free(e);
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

/* A denominator is the right operand of a left-associative `/`: a factor
 * that itself prints as a quotient (a negative Power, a Rational) must be
 * parenthesised, or the text re-reads as a different value -- "1/1/x" is
 * (1/1)/x = 1/x, not x. */
/* InputForm / ToString[_, InputForm] print machine reals with the shortest
 * digits that read back to the same double (a 17-digit mantissa carries the
 * machine mark `), and *^ literals parse with one correct rounding, so a
 * printed expression re-reads to an identical one. */
void test_inputform_real_roundtrip() {
    static const struct { const char* in; const char* out; } cases[] = {
        {"ToString[{-1.8954942670339809, 2.123456789123, 0.1, 3.}, InputForm]",
         "{-1.8954942670339809, 2.123456789123, 0.1, 3.0}"},
        {"ToString[{N[Pi], 1/3., 1.5*^-5, 6.02*^23}, InputForm]",
         "{3.141592653589793, 0.3333333333333333, 1.5*^-5, 6.02*^23}"},
        {"ToString[1.5000000000000002`*^-5, InputForm]", "1.5000000000000002`*^-5"},
        {"With[{v = {N[Pi], 1/3., 0.1, 1.2345678901234567*^-5, 2.^70, N[E]^-40}}, "
         "ToExpression[ToString[v, InputForm]] === v]", "True"},
        {"1.5*^-5 == 0.000015 && 6.02*^23 === 6.02`*^23 && 6.02*^23 == 6.02*10^23", "True"},
        {"ToString[{2.5, \"a\\\"b\"}, InputForm]", "{2.5, \"a\\\"b\"}"},
        {"ToString[{1.5, N[Pi]}]", "{1.5, 3.14159}"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        Expr* res = evaluate(e);
        char* str = (res->type == EXPR_STRING) ? strdup(res->data.string) : expr_to_string(res);
        if (!str || strcmp(str, cases[i].out) != 0)
            printf("  %s\n    got      %s\n    expected %s\n",
                   cases[i].in, str ? str : "(null)", cases[i].out);
        ASSERT(str && strcmp(str, cases[i].out) == 0);
        free(str);
        expr_free(e);
        expr_free(res);
    }
    printf("PASS: InputForm machine reals round-trip\n");
}

void test_quotient_denominator_parenthesised() {
    static const struct { const char* in; const char* out; } cases[] = {
        {"HoldForm[Power[Power[x, -1], -1]]",          "1/(1/x)"},
        {"HoldForm[Times[a, Power[Power[x, -1], -1]]]", "a/(1/x)"},
        {"HoldForm[Power[Rational[2, 3], -1]]",         "1/(2/3)"},
        {"HoldForm[Power[Power[x, -1], -2]]",           "1/(1/x)^2"},
        {"1/(2 x)",                                     "1/(2 x)"},
        {"a/b^2",                                       "a/b^2"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        Expr* res = evaluate(e);
        char* str = expr_to_string(res);
        if (!str || strcmp(str, cases[i].out) != 0)
            printf("  %s\n    got      %s\n    expected %s\n",
                   cases[i].in, str ? str : "(null)", cases[i].out);
        ASSERT(str && strcmp(str, cases[i].out) == 0);
        free(str);
        expr_free(e);
        expr_free(res);
    }
    printf("PASS: quotient-valued denominators are parenthesised\n");
}

/* Dot prints infix, `a . b`, at the parser's Dot level: tighter than Times
 * and Divide, looser than Power. Parenthesisation follows the reference
 * system's InputForm (Times renders as juxtaposition here). Every form must
 * re-parse to the same expression. */
static void test_dot_infix(void) {
    static const struct { const char* in; const char* out; } cases[] = {
        { "a . b",                 "a . b" },
        { "Dot[a, b, c]",          "a . b . c" },
        { "(a + b) . c",           "(a + b) . c" },
        { "Dot[a^2, b]",           "a^2 . b" },
        { "Dot[a, b]^2",           "(a . b)^2" },
        { "-a . b",                "-a . b" },
        { "Dot[-a, b]",            "(-a) . b" },
        { "Dot[a, -b]",            "a . (-b)" },
        { "Dot[a b, c]",           "(a b) . c" },
        { "Dot[1/2, a]",           "(1/2) . a" },
        { "f[a . b]",              "f[a . b]" },
        { "Dot[a, b] + c",         "c + a . b" },
        { "Dot[a, b] -> c",        "a . b -> c" },
        { "Dot[a, b][x]",          "(a . b)[x]" },
        { "Dot[a]",                "a" },
        { "Hold[Dot[a]]",          "Hold[Dot[a]]" },
        { "{1, 2} . {c, d}",       "c + 2 d" },
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Expr* e = parse_expression(cases[i].in);
        ASSERT(e != NULL);
        Expr* res = evaluate(e);
        char* str = expr_to_string(res);
        if (!str || strcmp(str, cases[i].out) != 0)
            printf("  %s\n    got      %s\n    expected %s\n",
                   cases[i].in, str ? str : "(null)", cases[i].out);
        ASSERT(str && strcmp(str, cases[i].out) == 0);
        /* Round trip: the printed form parses back to the same tree. */
        Expr* back = parse_expression(str);
        ASSERT(back != NULL);
        Expr* back_ev = evaluate(back);
        ASSERT(expr_eq(back_ev, res));
        expr_free(back);
        expr_free(back_ev);
        free(str);
        expr_free(e);
        expr_free(res);
    }
    assert_eval_eq("Dot[a, b]", "Dot[a, b]", 1);   /* FullForm unchanged */
    printf("PASS: Dot prints infix with precedence-aware parentheses\n");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_print_basic);
    TEST(test_fullform_wrapper);
    TEST(test_inputform_wrapper);
    TEST(test_negative_bigint_in_plus);
    TEST(test_int64_min_coefficient);
    TEST(test_holdform);
    TEST(test_series_latex);
    TEST(test_operator_latex);
    TEST(test_printer_directive_has_no_latex);
    TEST(test_standardform_latex);
    TEST(test_selectable_standardform);
    TEST(test_radical_body_unbracketed);
    TEST(test_texform_fraction_parenthesisation);
    TEST(test_quotient_denominator_parenthesised);
    TEST(test_inputform_real_roundtrip);
    TEST(test_dot_infix);

    printf("All print tests passed!\n");
    return 0;
}

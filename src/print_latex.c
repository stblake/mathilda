/*
 * print_latex.c — StandardForm-style LaTeX serialisation of Mathilda Expr* trees.
 *
 * Converts an expression to a KaTeX-compatible LaTeX string.
 * Handles: fractions (via Times/Power/-1 and Rational), roots, integer/real
 * powers, Greek symbols, trig/log/exp functions, sums, products, integrals,
 * limits, lists, and complex numbers.
 *
 * Design:
 *   - expr_to_latex() is the public entry point.
 *   - to_latex_prec(buf, e, ctx_prec) is the recursive workhorse; ctx_prec
 *     controls when to add parentheses.
 *   - A growing-string LBuf avoids repeated reallocation.
 */

#include "print_latex.h"
#include "print.h"          /* expr_to_string, for fallback */
#include "sym_names.h"
#include "expr.h"
#include "ndarray.h"     /* packed-list and NDArray[...] rendering */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

/* =========================================================================
 * Growing string buffer
 * ======================================================================== */

typedef struct { char* s; size_t len, cap; } LBuf;

static void lb_init(LBuf* b) {
    b->cap = 512; b->len = 0;
    b->s = malloc(b->cap);
    if (b->s) b->s[0] = '\0';
}

static void lb_ensure(LBuf* b, size_t extra) {
    if (!b->s) return;
    while (b->len + extra + 1 > b->cap) {
        b->cap *= 2;
        b->s = realloc(b->s, b->cap);
        if (!b->s) return;
    }
}

static void lb_cat(LBuf* b, const char* t) {
    if (!b->s || !t) return;
    size_t n = strlen(t);
    lb_ensure(b, n);
    if (!b->s) return;
    memcpy(b->s + b->len, t, n + 1);
    b->len += n;
}

static void lb_catf(LBuf* b, const char* fmt, ...) {
    char tmp[256];
    va_list ap; va_start(ap, fmt);
    vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    lb_cat(b, tmp);
}

/* =========================================================================
 * Operator precedence for parenthesisation
 * ======================================================================== */

/* Below PREC_ADD, in the parser's own order (docs/spec/operators.md): the
 * assignment/rule arrows bind loosest, then Or, And, Not, then the relations,
 * and only then arithmetic. Numeric values are otherwise arbitrary — only the
 * ordering is load-bearing, since it is what decides parenthesisation. */
#define PREC_SET    1   /* Set, SetDelayed */
#define PREC_RULE   2   /* Rule, RuleDelayed */
#define PREC_OR     4   /* Or */
#define PREC_AND    5   /* And */
#define PREC_NOT    6   /* Not */
#define PREC_REL    8   /* Equal, Less, ... and the Inequality chain */
#define PREC_ADD   10   /* Plus */
#define PREC_MUL   20   /* Times */
#define PREC_NEG   15   /* unary minus (between add and mul) */
#define PREC_POW   30   /* Power */
#define PREC_ATOM  99   /* numbers, symbols — never need parens */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* =========================================================================
 * Heads that print INFIX.
 *
 * Without this table every one of them fell through to the generic
 * `Head[a, b]` arm, so the notebook typeset a rule as `Rule[y[x], ...]` —
 * FullForm dressed up as mathematics. That is what every Solve and DSolve
 * result is made of, so it was the most visible output in the application.
 * The CLI's own TeX renderer (print.c) already had the same list; keeping the
 * two in step is why the operators and spellings here match it exactly.
 * ======================================================================== */

typedef struct { const char* sym; const char* tex; int prec; } InfixTeX;

static const InfixTeX INFIX_MAP[] = {
    {"Set",          "=",             PREC_SET},
    {"SetDelayed",   ":=",            PREC_SET},
    {"Rule",         "\\to ",         PREC_RULE},
    /* `:\to`, not `\to`: a delayed rule is a different object from an immediate
     * one, and rendering both the same way makes two different results
     * indistinguishable on screen. Mirrors SetDelayed's `:=`. */
    {"RuleDelayed",  ":\\to ",        PREC_RULE},
    {"Or",           "\\lor ",        PREC_OR},
    {"And",          "\\land ",       PREC_AND},
    {"Equal",        "=",             PREC_REL},
    {"Unequal",      "\\neq ",        PREC_REL},
    {"Less",         "<",             PREC_REL},
    {"Greater",      ">",             PREC_REL},
    {"LessEqual",    "\\leq ",        PREC_REL},
    {"GreaterEqual", "\\geq ",        PREC_REL},
    {"SameQ",        "\\equiv ",      PREC_REL},
    {"UnsameQ",      "\\not\\equiv ", PREC_REL},
    {NULL, NULL, 0}
};

static const InfixTeX* find_infix(const char* name) {
    for (int i = 0; INFIX_MAP[i].sym; i++)
        if (strcmp(name, INFIX_MAP[i].sym) == 0) return &INFIX_MAP[i];
    return NULL;
}

static int expr_prec(const Expr* e) {
    if (!e) return PREC_ATOM;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT  || e->type == EXPR_SYMBOL ||
        e->type == EXPR_STRING)  return PREC_ATOM;
    if (head_is(e, SYM_Plus))    return PREC_ADD;
    if (head_is(e, SYM_Times))   return PREC_MUL;
    if (head_is(e, SYM_Power))   return PREC_POW;
    /* The infix heads must report their level too, or a nested one never gets
     * parenthesised: `Or[And[a, b], c]` would print as `a \land b \lor c`. */
    if (e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Not && e->data.function.arg_count == 1) return PREC_NOT;
        if (h == SYM_Inequality) return PREC_REL;
        const InfixTeX* op = find_infix(h);
        if (op && e->data.function.arg_count >= 2) return op->prec;
    }
    return PREC_ATOM;
}

/* =========================================================================
 * Symbol → LaTeX mapping
 * ======================================================================== */

typedef struct { const char* sym; const char* latex; } SymTeX;

static const SymTeX SYM_MAP[] = {
    /* constants */
    {"Pi",          "\\pi"},
    {"E",           "e"},
    {"I",           "i"},
    {"Infinity",    "\\infty"},
    {"ComplexInfinity","\\tilde{\\infty}"},
    {"EulerGamma",  "\\gamma"},
    {"GoldenRatio", "\\varphi"},
    {"Catalan",     "G"},          /* Catalan's constant */
    /* Upright, not italic: these are words, and math mode would set `False` as
     * a product of five variables. Matches the CLI TeX renderer (print.c). */
    {"True",          "\\text{True}"},
    {"False",         "\\text{False}"},
    {"Null",          "\\text{Null}"},
    {"Indeterminate", "\\text{Indeterminate}"},
    /* Greek uppercase */
    {"Alpha","\\alpha"},{"Beta","\\beta"},{"Gamma","\\Gamma"},
    {"Delta","\\Delta"},{"Epsilon","\\epsilon"},{"Zeta","\\zeta"},
    {"Eta","\\eta"},{"Theta","\\Theta"},{"Iota","\\iota"},
    {"Kappa","\\kappa"},{"Lambda","\\Lambda"},{"Mu","\\mu"},
    {"Nu","\\nu"},{"Xi","\\Xi"},{"Rho","\\rho"},
    {"Sigma","\\Sigma"},{"Tau","\\tau"},{"Upsilon","\\upsilon"},
    {"Phi","\\Phi"},{"Chi","\\chi"},{"Psi","\\Psi"},{"Omega","\\Omega"},
    {NULL, NULL}
};

static const char* sym_latex(const char* name) {
    for (int i = 0; SYM_MAP[i].sym; i++)
        if (strcmp(name, SYM_MAP[i].sym) == 0) return SYM_MAP[i].latex;
    return name;
}

/* =========================================================================
 * Forward declaration
 * ======================================================================== */

static void to_latex_prec(LBuf* b, const Expr* e, int ctx_prec);

/* Wrap e in braces if its precedence < ctx_prec */
static void to_latex_maybe_paren(LBuf* b, const Expr* e, int ctx_prec) {
    int p = expr_prec(e);
    int need = (p < ctx_prec);
    if (need) lb_cat(b, "(");
    to_latex_prec(b, e, p);
    if (need) lb_cat(b, ")");
}

/* An operand of one of the loose infix operators (=, ->, &&, <, ...).
 *
 * Not to_latex_maybe_paren: that one hands the operand its OWN precedence as the
 * context, which for an atom is PREC_ATOM, and a negative integer parenthesises
 * itself above PREC_ADD -- right for `(-2)^n`, wrong for `x = (-1)`. Here the
 * inner context is PREC_ADD, so a bare `-1` stays bare and only a genuinely
 * looser operator gets brackets. */
static void to_latex_operand(LBuf* b, const Expr* e, int need_prec) {
    int paren = (expr_prec(e) < need_prec);
    if (paren) lb_cat(b, "\\left(");
    to_latex_prec(b, e, PREC_ADD);
    if (paren) lb_cat(b, "\\right)");
}

/* =========================================================================
 * Helpers: detect numeric factors
 * ======================================================================== */

/* True if e is the integer -1 */
static int is_neg_one(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == -1;
}

/* Get the numeric value of an atom; returns 0 and sets ok=0 on failure */
static double atom_value(const Expr* e, int* ok) {
    *ok = 1;
    if (!e) { *ok = 0; return 0; }
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_BIGINT)  return mpz_get_d(e->data.bigint);
    *ok = 0; return 0;
}

/* =========================================================================
 * Fraction detection: split Times[...] into numerator / denominator factors.
 * Returns 1 if the expression has a non-trivial denominator.
 * ======================================================================== */

/* Factors that belong in the denominator (Power[x, n] with n < 0) */
static int is_denom_factor(const Expr* e) {
    if (!head_is(e, SYM_Power) || e->data.function.arg_count < 2) return 0;
    const Expr* exp = e->data.function.args[1];
    int ok;
    double v = atom_value(exp, &ok);
    if (ok && v < 0) return 1;
    /* Rational exponent: Power[x, Rational[-p, q]] */
    if (head_is(exp, SYM_Rational) && exp->data.function.arg_count == 2) {
        double p = atom_value(exp->data.function.args[0], &ok);
        if (ok && p < 0) return 1;
    }
    return 0;
}

/* Collect numerator and denominator factors from a Times expression.
 * caller passes fixed-size arrays; returns 0 if not a fraction. */
#define MAX_FACTORS 32

/* Render a denominator factor: Power[x, -n] → x^{n} (with positive exponent) */
static void render_denom_factor(LBuf* b, const Expr* pow_expr) {
    /* pow_expr is Power[base, neg_exp] */
    const Expr* base = pow_expr->data.function.args[0];
    const Expr* nexp = pow_expr->data.function.args[1];
    to_latex_maybe_paren(b, base, PREC_POW);
    /* Negate the exponent */
    if (nexp->type == EXPR_INTEGER && nexp->data.integer == -1) {
        return; /* just the base */
    }
    lb_cat(b, "^{");
    /* Render |exponent| */
    if (nexp->type == EXPR_INTEGER) {
        lb_catf(b, "%lld", (long long)-nexp->data.integer);
    } else if (head_is(nexp, SYM_Rational)) {
        /* Render as fraction: Rational[-p, q] → p/q */
        int ok; double p = atom_value(nexp->data.function.args[0], &ok);
        double q = atom_value(nexp->data.function.args[1], &ok);
        if (q == 1) { lb_catf(b, "%g", -p); }
        else        { lb_catf(b, "\\frac{%g}{%g}", -p, q); }
    } else {
        to_latex_prec(b, nexp, PREC_ATOM); /* fallback */
    }
    lb_cat(b, "}");
}

/* =========================================================================
 * Trig / named function → LaTeX command mapping
 * ======================================================================== */

typedef struct { const char* sym; const char* latex_cmd; int paren; } FuncTeX;

/* paren=1 → wrap arg in parentheses, paren=0 → use braces (e.g. sqrt) */
static const FuncTeX FUNC_MAP[] = {
    {"Sin","\\sin",1}, {"Cos","\\cos",1}, {"Tan","\\tan",1},
    {"Csc","\\csc",1}, {"Sec","\\sec",1}, {"Cot","\\cot",1},
    {"ArcSin","\\arcsin",1}, {"ArcCos","\\arccos",1}, {"ArcTan","\\arctan",1},
    {"Sinh","\\sinh",1}, {"Cosh","\\cosh",1}, {"Tanh","\\tanh",1},
    {"Log","\\ln",1},
    {"Exp","\\exp",1},
    {"Abs","\\left|%s\\right|",0},  /* special: inline arg */
    {"Sqrt","\\sqrt",0},
    {"Floor","\\lfloor %s \\rfloor",0},
    {"Ceiling","\\lceil %s \\rceil",0},
    {"GCD","\\gcd",1},
    {"LCM","\\mathrm{lcm}",1},
    {"Det","\\det",1},
    {"Tr","\\mathrm{tr}",1},
    {"Max","\\max",1},
    {"Min","\\min",1},
    {NULL,NULL,0}
};

static const FuncTeX* find_func(const char* sym) {
    for (int i = 0; FUNC_MAP[i].sym; i++)
        if (strcmp(sym, FUNC_MAP[i].sym) == 0) return &FUNC_MAP[i];
    return NULL;
}

/* =========================================================================
 * Core recursive renderer
 * ======================================================================== */

static void render_times(LBuf* b, const Expr* e, int ctx_prec);
static void render_plus (LBuf* b, const Expr* e, int ctx_prec);

/* One level of a dense buffer as {...}, recursing on rank. `idx` walks the
 * row-major buffer; ndarray_buffer_element_to_expr is the single place that
 * decides an element's HEAD, so an int64 buffer typesets as integers. */
static void to_latex_packed(LBuf* b, const Expr* a, int level, size_t* idx) {
    const NDArrayData* nd = &a->data.ndarray;
    if (level == nd->rank) {
        Expr* leaf = ndarray_buffer_element_to_expr(nd->data, (*idx)++, nd->dtype);
        to_latex_prec(b, leaf, PREC_ADD);
        expr_free(leaf);
        return;
    }
    lb_cat(b, "\\{");
    for (int64_t i = 0; i < nd->dims[level]; i++) {
        if (i) lb_cat(b, ", ");
        to_latex_packed(b, a, level + 1, idx);
    }
    lb_cat(b, "\\}");
}

static void to_latex_prec(LBuf* b, const Expr* e, int ctx_prec) {
    if (!b->s || !e) return;

    /* ---- Atoms ---- */
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 0 && ctx_prec > PREC_ADD) {
            lb_cat(b, "(");
            lb_catf(b, "%lld", (long long)e->data.integer);
            lb_cat(b, ")");
        } else {
            lb_catf(b, "%lld", (long long)e->data.integer);
        }
        return;
    }
    if (e->type == EXPR_REAL) {
        /* Use up to 6 significant digits; strip trailing zeros */
        char buf[64];
        snprintf(buf, sizeof(buf), "%g", e->data.real);
        lb_cat(b, buf);
        return;
    }
    if (e->type == EXPR_BIGINT) {
        char* s = mpz_get_str(NULL, 10, e->data.bigint);
        lb_cat(b, s);
        free(s);
        return;
    }
    if (e->type == EXPR_SYMBOL) {
        lb_cat(b, sym_latex(e->data.symbol.name));
        return;
    }
    if (e->type == EXPR_STRING) {
        /* Use a directional quote pair so the opening quote is a left double
         * quote: a straight " renders as a closing quote in the math font, so
         * "apples" would show as ”apples”. U+201C / U+201D give “apples”. */
        lb_cat(b, "\\text{\xe2\x80\x9c");
        lb_cat(b, e->data.string);
        lb_cat(b, "\xe2\x80\x9d}");
        return;
    }

    /* ---- A dense buffer: packed List, or an explicit NDArray[...] ----
     * Streamed from the buffer one element at a time, the same way print.c does
     * it, so a 10^6-element value does not materialise 10^6 Expr nodes just to
     * be typeset. Without this arm the LaTeX came back EMPTY for every packed
     * result (and always had for a visible NDArray[...]), because the
     * EXPR_FUNCTION test below is the printer's only gate. */
    if (e->type == EXPR_NDARRAY) {
        bool visible = (e->data.ndarray.present_as == NDA_HEAD_NDARRAY);
        if (visible) lb_cat(b, "\\text{NDArray}\\left[");
        size_t idx = 0;
        to_latex_packed(b, e, 0, &idx);
        if (visible) lb_cat(b, "\\right]");
        return;
    }

    if (e->type != EXPR_FUNCTION) return;

    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    const Expr** args = (const Expr**)e->data.function.args;

    if (!head || head->type != EXPR_SYMBOL) goto fallback;
    const char* hname = head->data.symbol.name;

    /* ---- HoldForm[x] → x (transparent) ----
     * HoldForm suppresses evaluation but is invisible when printed; it must
     * render its argument at the enclosing precedence, exactly as the plain
     * printer does (print.c). Trace[expr] wraps each step in HoldForm, so
     * without this the notebook LaTeX would leak "HoldForm[...]". */
    if (hname == SYM_HoldForm && argc == 1) {
        to_latex_prec(b, args[0], ctx_prec);
        return;
    }

    /* ---- SeriesData[...] → a0 + a1 (x-x0) + ... + O[x-x0]^n ----
     * SeriesData is an inert head that must display as the sum it represents,
     * not as the literal 6-argument container. The plain printer (print.c)
     * builds an equivalent Plus/Times/Power/O tree and delegates; we call the
     * same shared builder so the notebook LaTeX matches the CLI exactly (gh
     * #22). A NULL result (unrenderable shape) falls through to generic. */
    if (hname == SYM_SeriesData && argc == 6) {
        Expr* disp = series_data_to_display_expr((Expr*)e);
        if (disp) {
            to_latex_prec(b, disp, ctx_prec);
            expr_free(disp);
            return;
        }
    }

    /* ---- Rational[p, q] ---- */
    if (hname == SYM_Rational && argc == 2) {
        lb_cat(b, "\\frac{");
        int ok; double p = atom_value(args[0], &ok);
        if (p < 0) { lb_catf(b, "-"); p = -p; }
        lb_catf(b, "%g}{", p);
        to_latex_prec(b, args[1], PREC_ATOM);
        lb_cat(b, "}");
        return;
    }

    /* ---- Complex[a, b] → a + b i ---- */
    if (hname == SYM_Complex && argc == 2) {
        to_latex_prec(b, args[0], PREC_ADD);
        lb_cat(b, "+");
        to_latex_prec(b, args[1], PREC_MUL);
        lb_cat(b, "i");
        return;
    }

    /* ---- C[k] → c_k : the DSolve/Reduce/Integrate generated constant of
     * integration, matching Mathematica's TeXForm (c_1, c_2, ...) and the
     * CLI TeX renderer (print.c). A single-character subscript stays bare;
     * anything longer is braced so LaTeX groups the whole subscript. */
    if (strcmp(hname, "C") == 0 && argc == 1) {
        const Expr* sub = args[0];
        int bare = (sub->type == EXPR_INTEGER
                    && sub->data.integer >= 0 && sub->data.integer <= 9)
                || (sub->type == EXPR_SYMBOL && strlen(sub->data.symbol.name) == 1);
        lb_cat(b, "c_");
        if (!bare) lb_cat(b, "{");
        to_latex_prec(b, sub, PREC_ATOM);
        if (!bare) lb_cat(b, "}");
        return;
    }

    /* ---- The loose infix operators: a -> b, x == 1, a && b, ... ----
     * Above Plus so the whole ladder reads loosest-first. Each operand is
     * bracketed only when it is looser than the operator holding it. */
    {
        const InfixTeX* op = find_infix(hname);
        if (op && argc >= 2) {
            for (size_t i = 0; i < argc; i++) {
                if (i) lb_cat(b, op->tex);
                to_latex_operand(b, args[i], op->prec + 1);
            }
            return;
        }
    }

    /* ---- Not[a] → \neg a ---- */
    if (hname == SYM_Not && argc == 1) {
        lb_cat(b, "\\neg ");
        to_latex_operand(b, args[0], PREC_NOT + 1);
        return;
    }

    /* ---- Inequality[v0, op0, v1, op1, v2, ...] → a < b <= c ----
     * The chained form `1 < x < 2` parses to this rather than to nested Less,
     * and Reduce returns it, so without this arm an interval printed as
     * `Inequality[1, Less, x, Less, 2]`. The odd positions are operator SYMBOLS,
     * looked up in the same table as the two-argument heads. */
    if (hname == SYM_Inequality && argc >= 3 && (argc % 2) == 1) {
        for (size_t i = 0; i < argc; i++) {
            if (i % 2 == 0) {
                to_latex_operand(b, args[i], PREC_REL + 1);
            } else {
                const InfixTeX* rel = args[i]->type == EXPR_SYMBOL
                    ? find_infix(args[i]->data.symbol.name) : NULL;
                /* An unrecognised relation is still better spelled out than
                   dropped, which would invert the meaning of the chain. */
                lb_cat(b, rel ? rel->tex : "\\,?\\,");
            }
        }
        return;
    }

    /* ---- Plus[...] ---- */
    if (hname == SYM_Plus && argc >= 1) {
        render_plus(b, e, ctx_prec);
        return;
    }

    /* ---- Times[...] ---- */
    if (hname == SYM_Times && argc >= 1) {
        render_times(b, e, ctx_prec);
        return;
    }

    /* ---- Power[base, exp] ---- */
    if (hname == SYM_Power && argc == 2) {
        const Expr* base = args[0];
        const Expr* exp  = args[1];
        /* Power[x, -1] → \frac{1}{x} */
        if (exp->type == EXPR_INTEGER && exp->data.integer == -1) {
            lb_cat(b, "\\frac{1}{");
            to_latex_prec(b, base, PREC_ATOM);
            lb_cat(b, "}");
            return;
        }

        /* Power[x, Rational[1, 2]] → \sqrt{x} */
        if (head_is(exp, SYM_Rational) && exp->data.function.arg_count == 2) {
            const Expr* p = exp->data.function.args[0];
            const Expr* q = exp->data.function.args[1];
            if (p->type == EXPR_INTEGER && p->data.integer == 1) {
                if (q->type == EXPR_INTEGER && q->data.integer == 2) {
                    lb_cat(b, "\\sqrt{");
                    to_latex_prec(b, base, PREC_ATOM);
                    lb_cat(b, "}");
                    return;
                }
                /* n-th root */
                lb_cat(b, "\\sqrt[");
                to_latex_prec(b, q, PREC_ATOM);
                lb_cat(b, "]{");
                to_latex_prec(b, base, PREC_ATOM);
                lb_cat(b, "}");
                return;
            }
        }

        /* Power[E, x] → e^{x} */
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E) {
            lb_cat(b, "e^{");
            to_latex_prec(b, exp, PREC_ATOM);
            lb_cat(b, "}");
            return;
        }

        /* General Power[base, exp] → base^{exp} */
        to_latex_maybe_paren(b, base, PREC_POW + 1); /* base needs parens if lower prec */
        lb_cat(b, "^{");
        to_latex_prec(b, exp, PREC_ATOM);
        lb_cat(b, "}");
        return;
    }

    /* ---- Sqrt[x] ---- */
    if (hname == SYM_Sqrt && argc == 1) {
        lb_cat(b, "\\sqrt{");
        to_latex_prec(b, args[0], PREC_ATOM);
        lb_cat(b, "}");
        return;
    }

    /* ---- List[a, b, ...] → {a, b, ...} ---- */
    if (hname == SYM_List) {
        lb_cat(b, "\\{");
        for (size_t i = 0; i < argc; i++) {
            if (i) lb_cat(b, ", ");
            to_latex_prec(b, args[i], PREC_ADD);
        }
        lb_cat(b, "\\}");
        return;
    }

    /* ---- Association[k->v, ...] → ⟨| k → v, ... |⟩ ---- */
    if (hname == SYM_Association) {
        bool all_rules = true;
        for (size_t i = 0; i < argc; i++) {
            const Expr* r = args[i];
            if (!((head_is(r, SYM_Rule) || head_is(r, SYM_RuleDelayed)) &&
                  r->data.function.arg_count == 2)) { all_rules = false; break; }
        }
      if (all_rules) {
        lb_cat(b, "\\left\\langle\\!\\left|\\, ");
        for (size_t i = 0; i < argc; i++) {
            if (i) lb_cat(b, ",\\; ");
            const Expr* rule = args[i];
            if ((head_is(rule, SYM_Rule) || head_is(rule, SYM_RuleDelayed)) &&
                rule->data.function.arg_count == 2) {
                to_latex_prec(b, rule->data.function.args[0], PREC_ATOM);
                lb_cat(b, " \\to ");
                to_latex_prec(b, rule->data.function.args[1], PREC_ADD);
            } else {
                to_latex_prec(b, rule, PREC_ADD);
            }
        }
        lb_cat(b, "\\, \\right|\\!\\right\\rangle");
        return;
      }
      /* else fall through to generic head[args] printing below */
    }

    /* ---- Factorial[n] → n! ---- */
    if (hname == SYM_Factorial && argc == 1) {
        to_latex_maybe_paren(b, args[0], PREC_POW);
        lb_cat(b, "!");
        return;
    }

    /* ---- Binomial[n, k] ---- */
    if (strcmp(hname, "Binomial") == 0 && argc == 2) {
        lb_cat(b, "\\binom{");
        to_latex_prec(b, args[0], PREC_ATOM);
        lb_cat(b, "}{");
        to_latex_prec(b, args[1], PREC_ATOM);
        lb_cat(b, "}");
        return;
    }

    /* ---- Sum[f, {n, a, b}] ---- */
    if (hname == SYM_Sum && argc == 2 && head_is(args[1], SYM_List)
            && args[1]->data.function.arg_count == 3) {
        const Expr* iter = args[1];
        lb_cat(b, "\\sum_{");
        to_latex_prec(b, iter->data.function.args[0], PREC_ATOM);
        lb_cat(b, "=");
        to_latex_prec(b, iter->data.function.args[1], PREC_ATOM);
        lb_cat(b, "}^{");
        to_latex_prec(b, iter->data.function.args[2], PREC_ATOM);
        lb_cat(b, "} ");
        to_latex_prec(b, args[0], PREC_ADD);
        return;
    }

    /* ---- Product[f, {n, a, b}] ---- */
    if (hname == SYM_Product && argc == 2 && head_is(args[1], SYM_List)
            && args[1]->data.function.arg_count == 3) {
        const Expr* iter = args[1];
        lb_cat(b, "\\prod_{");
        to_latex_prec(b, iter->data.function.args[0], PREC_ATOM);
        lb_cat(b, "=");
        to_latex_prec(b, iter->data.function.args[1], PREC_ATOM);
        lb_cat(b, "}^{");
        to_latex_prec(b, iter->data.function.args[2], PREC_ATOM);
        lb_cat(b, "} ");
        to_latex_prec(b, args[0], PREC_ADD);
        return;
    }

    /* ---- Integrate[f, {x, a, b}] ---- */
    if (hname == SYM_Integrate && argc == 2) {
        if (head_is(args[1], SYM_List) && args[1]->data.function.arg_count == 3) {
            const Expr* iter = args[1];
            lb_cat(b, "\\int_{");
            to_latex_prec(b, iter->data.function.args[1], PREC_ATOM);
            lb_cat(b, "}^{");
            to_latex_prec(b, iter->data.function.args[2], PREC_ATOM);
            lb_cat(b, "} ");
            to_latex_prec(b, args[0], PREC_MUL);
            lb_cat(b, "\\,d");
            to_latex_prec(b, iter->data.function.args[0], PREC_ATOM);
        } else {
            lb_cat(b, "\\int ");
            to_latex_prec(b, args[0], PREC_MUL);
            lb_cat(b, "\\,d");
            to_latex_prec(b, args[1], PREC_ATOM);
        }
        return;
    }

    /* ---- D[f, x] → f' or \frac{d}{dx} f ---- */
    if (hname == SYM_D && argc == 2) {
        lb_cat(b, "\\frac{d}{d");
        to_latex_prec(b, args[1], PREC_ATOM);
        lb_cat(b, "}\\left(");
        to_latex_prec(b, args[0], PREC_ADD);
        lb_cat(b, "\\right)");
        return;
    }

    /* ---- Limit[f, x->a] ---- */
    if (strcmp(hname, "Limit") == 0 && argc >= 2) {
        const Expr* rule = args[1];
        lb_cat(b, "\\lim_{");
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            to_latex_prec(b, rule->data.function.args[0], PREC_ATOM);
            lb_cat(b, "\\to ");
            to_latex_prec(b, rule->data.function.args[1], PREC_ATOM);
        } else {
            to_latex_prec(b, rule, PREC_ATOM);
        }
        lb_cat(b, "} ");
        to_latex_prec(b, args[0], PREC_ADD);
        return;
    }

    /* ---- Log[b, x] → \log_b x ---- */
    if (hname == SYM_Log && argc == 2) {
        lb_cat(b, "\\log_{");
        to_latex_prec(b, args[0], PREC_ATOM);
        lb_cat(b, "}\\left(");
        to_latex_prec(b, args[1], PREC_ADD);
        lb_cat(b, "\\right)");
        return;
    }

    /* ---- Named functions: sin, cos, sqrt, ... ---- */
    if (head->type == EXPR_SYMBOL) {
        const FuncTeX* ft = find_func(hname);
        if (ft) {
            if (strcmp(ft->latex_cmd, "\\sqrt") == 0 && argc == 1) {
                lb_cat(b, "\\sqrt{");
                to_latex_prec(b, args[0], PREC_ATOM);
                lb_cat(b, "}");
                return;
            }
            if (ft->paren && argc >= 1) {
                lb_cat(b, ft->latex_cmd);
                lb_cat(b, "\\left(");
                for (size_t i = 0; i < argc; i++) {
                    if (i) lb_cat(b, ", ");
                    to_latex_prec(b, args[i], PREC_ADD);
                }
                lb_cat(b, "\\right)");
                return;
            }
        }
    }

    /* ---- Generic function: head[a, b, ...] ---- */
    fallback: {
        /* For a symbol-headed function, render Name[arg, ...] with each argument
         * in LaTeX, so nested strings (curly quotes via \text), fractions, etc.
         * render properly rather than dumping the plain re-parseable form (which
         * would show straight quotes, e.g. Key["a"]). */
        if (e->type == EXPR_FUNCTION && head && head->type == EXPR_SYMBOL) {
            lb_cat(b, sym_latex(head->data.symbol.name));
            lb_cat(b, "[");
            for (size_t i = 0; i < argc; i++) {
                if (i) lb_cat(b, ", ");
                to_latex_prec(b, args[i], 0);
            }
            lb_cat(b, "]");
            return;
        }
        /* Non-symbol head or atom: fall back to the plain string form. */
        char* s = expr_to_string((Expr*)e);  /* const cast: expr_to_string doesn't modify */
        if (s) { lb_cat(b, s); free(s); }
    }
}

/* =========================================================================
 * Plus renderer — handles subtraction (negative terms)
 * ======================================================================== */

/* True if e is Rational[p, q] with a negative integer numerator p. Such a
 * coefficient (e.g. -1/6) should display as subtraction, matching the plain
 * printer, rather than "+ -1/6" (surfaced by series output, gh #22). */
static int is_neg_rational(const Expr* e) {
    return head_is(e, SYM_Rational) && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer < 0;
}

/* Render |Rational[-p, q]| = \frac{p}{q} (numerator already known negative). */
static void render_rational_abs(LBuf* b, const Expr* e) {
    lb_cat(b, "\\frac{");
    lb_catf(b, "%lld", (long long)(-e->data.function.args[0]->data.integer));
    lb_cat(b, "}{");
    to_latex_prec(b, e->data.function.args[1], PREC_ATOM);
    lb_cat(b, "}");
}

/* True if e is a negative term: Times[-1, ...] or negative integer/real,
 * a negative rational, or a Times led by one of those. */
static int is_negative_term(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0;
    if (is_neg_rational(e))      return 1;
    if (head_is(e, SYM_Times) && e->data.function.arg_count >= 1) {
        const Expr* first = e->data.function.args[0];
        if (first->type == EXPR_INTEGER && first->data.integer < 0) return 1;
        if (first->type == EXPR_REAL    && first->data.real    < 0) return 1;
        if (is_neg_rational(first)) return 1;
    }
    return 0;
}

/* Return -e (negate a negative term) for display as subtraction */
static void render_negate(LBuf* b, const Expr* e) {
    if (e->type == EXPR_INTEGER) { lb_catf(b, "%lld", (long long)-e->data.integer); return; }
    if (e->type == EXPR_REAL)    { lb_catf(b, "%g",             -e->data.real);     return; }
    if (is_neg_rational(e))      { render_rational_abs(b, e); return; }
    if (head_is(e, SYM_Times) && e->data.function.arg_count >= 1) {
        const Expr* first = e->data.function.args[0];
        if (is_neg_rational(first)) {
            /* Times[-p/q, rest...] → (p/q) rest */
            render_rational_abs(b, first);
            for (size_t i = 1; i < e->data.function.arg_count; i++) {
                lb_cat(b, "\\,");
                to_latex_maybe_paren(b, e->data.function.args[i], PREC_MUL + 1);
            }
            return;
        }
        int ok; double v = atom_value(first, &ok);
        if (ok && v == -1 && e->data.function.arg_count == 2) {
            /* Times[-1, x] → just render x */
            to_latex_prec(b, e->data.function.args[1], PREC_MUL);
            return;
        }
        if (ok && v < 0) {
            /* Times[-n, ...] → render n * rest */
            lb_catf(b, "%g", -v);
            for (size_t i = 1; i < e->data.function.arg_count; i++) {
                lb_cat(b, " ");
                to_latex_maybe_paren(b, e->data.function.args[i], PREC_MUL + 1);
            }
            return;
        }
    }
    to_latex_prec(b, e, PREC_MUL);
}

static void render_plus(LBuf* b, const Expr* e, int ctx_prec) {
    size_t argc = e->data.function.arg_count;
    const Expr** args = (const Expr**)e->data.function.args;
    int need_paren = ctx_prec > PREC_ADD;
    if (need_paren) lb_cat(b, "\\left(");
    for (size_t i = 0; i < argc; i++) {
        if (i == 0) {
            to_latex_prec(b, args[i], PREC_ADD);
        } else if (is_negative_term(args[i])) {
            lb_cat(b, "-");
            render_negate(b, args[i]);
        } else {
            lb_cat(b, "+");
            to_latex_prec(b, args[i], PREC_ADD);
        }
    }
    if (need_paren) lb_cat(b, "\\right)");
}

/* =========================================================================
 * Times renderer — handles fractions and implicit multiplication
 * ======================================================================== */

static void render_times(LBuf* b, const Expr* e, int ctx_prec) {
    size_t argc = e->data.function.arg_count;
    const Expr** args = (const Expr**)e->data.function.args;

    /* Collect numerator and denominator factors */
    const Expr* num_f[MAX_FACTORS];
    const Expr* den_f[MAX_FACTORS];
    int nnum = 0, nden = 0;

    /* Separate sign/coefficient from positive factors */
    int sign = 1;

    for (size_t i = 0; i < argc && i < MAX_FACTORS; i++) {
        const Expr* f = args[i];
        if (is_denom_factor(f)) {
            den_f[nden++] = f;
        } else if (is_neg_one(f)) {
            sign = -sign;
        } else {
            num_f[nnum++] = f;
        }
    }

    int is_frac = nden > 0;
    int need_paren = ctx_prec > PREC_MUL && !is_frac;
    if (need_paren) lb_cat(b, "\\left(");
    if (sign < 0) lb_cat(b, "-");

    if (is_frac) {
        lb_cat(b, "\\frac{");
        if (nnum == 0) {
            lb_cat(b, "1");
        } else {
            for (int i = 0; i < nnum; i++) {
                if (i > 0) lb_cat(b, " ");
                to_latex_maybe_paren(b, num_f[i], PREC_MUL + 1);
            }
        }
        lb_cat(b, "}{");
        for (int i = 0; i < nden; i++) {
            if (i > 0) lb_cat(b, " ");
            render_denom_factor(b, den_f[i]);
        }
        lb_cat(b, "}");
    } else {
        for (int i = 0; i < nnum; i++) {
            if (i > 0) {
                /* No \cdot between number and symbol, but yes between two symbols */
                const Expr* prev = num_f[i-1];
                int prev_num = (prev->type == EXPR_INTEGER || prev->type == EXPR_REAL ||
                                prev->type == EXPR_BIGINT  || head_is(prev, SYM_Rational));
                (void)prev_num;
                /* Use thin space between terms */
                lb_cat(b, "\\,");
            }
            to_latex_maybe_paren(b, num_f[i], PREC_MUL + 1);
        }
    }
    if (need_paren) lb_cat(b, "\\right)");
}

/* =========================================================================
 * Printer directives: the expressions that have no typeset form
 *
 * InputForm, FullForm, TeXForm and NumberForm are instructions to the
 * PRINTER, not mathematics. The plain printer (print.c) consumes each one and
 * renders its argument in the notation asked for, so `expr // InputForm` comes
 * out of expr_to_string already in exactly the form the user requested.
 *
 * The LaTeX writer has no such notation to offer — typesetting *is*
 * StandardForm — so whatever it produced for these heads was wrong twice over:
 * the wrapper leaked into the output (`InputForm[\frac{1}{2}]`, which KaTeX
 * renders as an upright product) and what it wrapped was the StandardForm the
 * user had just asked NOT to see. Because both the notebook and the FFI prefer
 * the `latex` field whenever it is non-empty, that leak was the entire visible
 * result: `D[Log[1 - Sqrt[x]] Sqrt[x], x] // InputForm` typeset the derivative
 * and ignored the directive.
 *
 * Emitting no LaTeX is the answer rather than some fallback rendering: it hands
 * the consumer back to the plain-text payload, which is precisely the form the
 * directive asked for. The scan covers the whole tree rather than the outermost
 * head alone, so a directive nested anywhere (`Hold[InputForm[x]]`,
 * `{InputForm[1/2], 3}`) cannot leak either — there is no partial typesetting
 * that would be more faithful than the text the printer has already produced.
 * HoldForm is absent on purpose: it is transparent to both printers, and
 * to_latex_prec renders through it above.
 * ======================================================================== */
static int has_print_directive(const Expr* e) {
    /* Atoms carry none, and an EXPR_NDARRAY holds machine numbers, not heads. */
    if (!e || e->type != EXPR_FUNCTION) return 0;

    const Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol.name;   /* interned: compare by pointer */
        if (h == SYM_InputForm || h == SYM_FullForm
            || h == SYM_TeXForm || h == SYM_NumberForm)
            return 1;
    }
    if (has_print_directive(head)) return 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_print_directive(e->data.function.args[i])) return 1;
    return 0;
}

/* =========================================================================
 * Public entry point
 * ======================================================================== */

char* expr_to_latex(const Expr* e) {
    if (has_print_directive(e)) {
        char* none = malloc(1);          /* "": no typeset form — use the text */
        if (none) none[0] = '\0';
        return none;
    }

    LBuf b;
    lb_init(&b);
    if (!b.s) return NULL;
    to_latex_prec(&b, e, 0);
    if (!b.s) return NULL;
    return b.s;  /* caller frees */
}

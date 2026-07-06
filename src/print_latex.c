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

#define PREC_ADD   10   /* Plus */
#define PREC_MUL   20   /* Times */
#define PREC_NEG   15   /* unary minus (between add and mul) */
#define PREC_POW   30   /* Power */
#define PREC_ATOM  99   /* numbers, symbols — never need parens */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == sym;
}

static int expr_prec(const Expr* e) {
    if (!e) return PREC_ATOM;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT  || e->type == EXPR_SYMBOL ||
        e->type == EXPR_STRING)  return PREC_ATOM;
    if (head_is(e, SYM_Plus))    return PREC_ADD;
    if (head_is(e, SYM_Times))   return PREC_MUL;
    if (head_is(e, SYM_Power))   return PREC_POW;
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
        lb_cat(b, sym_latex(e->data.symbol));
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

    if (e->type != EXPR_FUNCTION) return;

    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    const Expr** args = (const Expr**)e->data.function.args;

    if (!head || head->type != EXPR_SYMBOL) goto fallback;
    const char* hname = head->data.symbol;

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
        if (base->type == EXPR_SYMBOL && base->data.symbol == SYM_E) {
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
            lb_cat(b, sym_latex(head->data.symbol));
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

/* True if e is a negative term: Times[-1, ...] or negative integer/real */
static int is_negative_term(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0;
    if (head_is(e, SYM_Times) && e->data.function.arg_count >= 1) {
        const Expr* first = e->data.function.args[0];
        if (first->type == EXPR_INTEGER && first->data.integer < 0) return 1;
        if (first->type == EXPR_REAL    && first->data.real    < 0) return 1;
    }
    return 0;
}

/* Return -e (negate a negative term) for display as subtraction */
static void render_negate(LBuf* b, const Expr* e) {
    if (e->type == EXPR_INTEGER) { lb_catf(b, "%lld", (long long)-e->data.integer); return; }
    if (e->type == EXPR_REAL)    { lb_catf(b, "%g",             -e->data.real);     return; }
    if (head_is(e, SYM_Times) && e->data.function.arg_count >= 1) {
        const Expr* first = e->data.function.args[0];
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
 * Public entry point
 * ======================================================================== */

char* expr_to_latex(const Expr* e) {
    LBuf b;
    lb_init(&b);
    if (!b.s) return NULL;
    to_latex_prec(&b, e, 0);
    if (!b.s) return NULL;
    return b.s;  /* caller frees */
}

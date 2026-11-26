#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// Print

#include "ds_ops.h"
#include "pvec.h"
#include "print.h"
#include "linkedlist.h"
#include "print_fortran.h"
#include "print_lean.h"
#include "expr.h"
#include "arithmetic.h"
#include "context.h"
#include "attr.h"           /* Definition: Attributes[f] line */
#include "sym_names.h"
#include "numberform.h"
#include "ndarray.h"
#include "sparse/sparse.h"
#include "compile/compiled_function.h"
#include "graph.h"   /* graph_is_list, for the Graph[...] summary form */
#include "graph_hyper.h"   /* hypergraph_is_valid, for Hypergraph[...] */
#include "symtab.h"        /* FormatValues */
#include "match.h"
#include "eval.h"
#include "message.h"
#include "linear_boxes.h"   /* box control characters print as escapes */
#include "utf8.h"
#include "common.h"        /* builtin_arg_error */
#include "numeric.h"       /* numeric_bits_to_digits: precision of an MPFR width */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>

#include <inttypes.h>
#include <math.h>

/* -v as an Integer: v = INT64_MIN has no int64 negation, so it goes through
 * GMP (it printed as --9223372036854775808). */
/* A negative integer exponent (machine or BigInt), negated (owned), or NULL.
 * Power printing treats x^(-n) the same way whatever the width of n: the
 * InputForm parentheses and the OutputForm denominator. */
static Expr* print_neg_int_exp(const Expr* exp);

static Expr* print_neg_int(int64_t v) {
    if (v != INT64_MIN) return expr_new_integer(-v);
    mpz_t z;
    mth_mpz_init_set_i64(z, v);
    mpz_neg(z, z);
    Expr* e = expr_new_bigint_from_mpz(z);
    mpz_clear(z);
    return e;
}

static Expr* print_neg_int_exp(const Expr* exp) {
    if (exp->type == EXPR_INTEGER && exp->data.integer < 0)
        return print_neg_int(exp->data.integer);
    if (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) < 0) {
        mpz_t z;
        mpz_init(z);
        mpz_neg(z, exp->data.bigint);
        Expr* e = expr_new_bigint_from_mpz(z);
        mpz_clear(z);
        return e;
    }
    return NULL;
}

/* An exponent Times[c, rest...] with c a negative number, negated (owned):
 * the reference puts such a power in the denominator of a product, as it
 * does x^(-2) -- c x^(-n) is c/x^n, a E^(-x) is a/E^x, 3^(x/2) 2^(-x) is
 * 3^(x/2)/2^x.  NULL for any other exponent.  A lone power keeps x^(-n). */
static Expr* print_neg_symbolic_exp(const Expr* x) {
    if (x->type != EXPR_FUNCTION || x->data.function.head->type != EXPR_SYMBOL ||
        x->data.function.head->data.symbol.name != SYM_Times || x->data.function.arg_count < 2)
        return NULL;
    const Expr* c = x->data.function.args[0];
    Expr* nc = NULL;
    int64_t n, d;
    if (c->type == EXPR_INTEGER && c->data.integer < 0) nc = print_neg_int(c->data.integer);
    else if (c->type == EXPR_BIGINT && mpz_sgn(c->data.bigint) < 0) {
        mpz_t z; mpz_init(z); mpz_neg(z, c->data.bigint);
        nc = expr_new_bigint_from_mpz(z); mpz_clear(z);
    }
    else if (is_rational((Expr*)c, &n, &d) && n < 0) nc = make_rational_neg(n, d);
    else if (c->type == EXPR_REAL && c->data.real < 0.0) nc = expr_new_real(-c->data.real);
    if (!nc) return NULL;
    size_t k = x->data.function.arg_count;
    bool unit = nc->type == EXPR_INTEGER && nc->data.integer == 1;
    if (unit && k == 2) { expr_free(nc); return expr_copy(x->data.function.args[1]); }
    size_t m = unit ? k - 1 : k;
    Expr** a = malloc(m * sizeof(Expr*));
    size_t j = 0;
    if (!unit) a[j++] = nc; else expr_free(nc);
    for (size_t i = 1; i < k; i++) a[j++] = expr_copy(x->data.function.args[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), a, m);
    free(a);
    return r;
}

/* ---------------------------------------------------------------------------
 * Output sink.  Every output primitive in this file writes through g_print_sink
 * (NULL => the real stdout).  The capture helpers below install an in-memory
 * stream here so expr_to_string & friends collect printed output without
 * reassigning `stdout`, which is a read-only (const) pointer on musl libc.
 * `printf`/`putchar` are redefined as file-local macros so the ~190 existing
 * call sites need no edit; fputs(…, stdout) sites use MTH_POUT explicitly.
 * numberform.c, the one other printer in the print_standard path, writes to
 * mth_out() (same sink) for the same reason.
 * ------------------------------------------------------------------------- */
static FILE* g_print_sink = NULL;          /* NULL => stdout */
#define MTH_POUT (g_print_sink ? g_print_sink : stdout)

FILE* mth_out(void) { return MTH_POUT; }
FILE* mth_out_push(FILE* f) { FILE* prev = g_print_sink; g_print_sink = f; return prev; }
void  mth_out_pop(FILE* prev) { g_print_sink = prev; }

#define printf(...) fprintf(MTH_POUT, __VA_ARGS__)
#define putchar(c)  putc((c), MTH_POUT)

static void print_standard(Expr* e, int parent_prec);
static void print_series_data(Expr* e, int parent_prec);
static bool print_plus_spliced_series(Expr* e);
static void print_tex(Expr* e, int parent_prec);

/* True when every argument of `e` is a two-argument Rule/RuleDelayed, i.e. a
 * well-formed association that should print in <| ... |> form. An Association
 * node with any non-rule argument (e.g. an unevaluated Association[x]) falls
 * back to ordinary head[args] printing, matching Wolfram behaviour. */
static bool assoc_prints_special(const Expr* e) {
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* a = e->data.function.args[i];
        if (!(a->type == EXPR_FUNCTION &&
              a->data.function.head->type == EXPR_SYMBOL &&
              (a->data.function.head->data.symbol.name == SYM_Rule ||
               a->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
              a->data.function.arg_count == 2))
            return false;
    }
    return true;
}

/* Returns a fresh expression holding -b, where b is a bigint node. The result
 * is normalised back to a machine int when it fits. Used by the Plus printers
 * to pull the sign of a negative bigint term onto the operator. */
static Expr* negate_bigint(Expr* b) {
    mpz_t n;
    mpz_init(n);
    mpz_neg(n, b->data.bigint);
    Expr* r = expr_new_bigint_from_mpz(n);
    mpz_clear(n);
    return expr_bigint_normalize(r);
}

/*
 * When non-zero, print_standard renders heads like SeriesData using their
 * literal head[args...] form rather than the pretty mathematical form.
 * Toggled while printing inside an InputForm[...] wrapper.
 */
static int g_inputform_depth = 0;

/* The active NumberForm[...] formatting context, or NULL. Installed by the
 * SYM_NumberForm branch of print_standard while it recurses into the wrapped
 * expression; every numeric leaf then routes through numberform_render_number
 * instead of the default real/integer formatter. Saved/restored so nested
 * NumberForm wrappers compose. */
static const NumberFormCtx* g_numberform_ctx = NULL;

/* True for the atomic numeric types NumberForm knows how to reformat. */
static bool nf_reformattable_atom(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Blank[] / Blank[h] / BlankSequence[...] / BlankNullSequence[...]: what
 * Pattern[x, b] prints as the atomic x_ / x_h / x__ / x___. */
static bool is_plain_blank(const Expr* b) {
    if (b->type != EXPR_FUNCTION || b->data.function.head->type != EXPR_SYMBOL ||
        b->data.function.arg_count > 1) return false;
    const char* h = b->data.function.head->data.symbol.name;
    return h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence;
}
/* A Complex node is printed as re + im I only when both parts are numbers;
 * Complex[a, b] with symbolic parts is not a number and prints as itself
 * (audit Complex#2). */
static bool print_complex_is_numeric(const Expr* e) {
    return e->data.function.arg_count == 2 &&
           expr_is_numeric_like(e->data.function.args[0]) &&
           expr_is_numeric_like(e->data.function.args[1]);
}

/* Rational[p, q] with integer parts: a rational number, printed p/q. */
static bool print_rational_is_numeric(const Expr* e) {
    return e->data.function.arg_count == 2 &&
           expr_is_integer_like(e->data.function.args[0]) &&
           expr_is_integer_like(e->data.function.args[1]);
}

static int get_expr_prec(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 9500;
    if (e->data.function.head->type != EXPR_SYMBOL) return 9500;
    
    const char* head = e->data.function.head->data.symbol.name;

    /* The binary infix operators print infix only with exactly two operands;
     * any other count is shown in functional form, Rule[a, b, c] (audit
     * Rule#2), whose precedence is that of application. */
    if ((head == SYM_Set || head == SYM_SetDelayed || head == SYM_Rule ||
         head == SYM_RuleDelayed || head == SYM_Condition) &&
        e->data.function.arg_count != 2) return 9500;
    if (head == SYM_Set || head == SYM_SetDelayed) return 500;
#if UP_VALUES
    if (head == SYM_UpSet || head == SYM_UpSetDelayed ||
        head == SYM_TagSet || head == SYM_TagSetDelayed || head == SYM_TagUnset) return 500;
#endif
    if (head == SYM_MessageName) return 8600;
    if (head == SYM_Rule || head == SYM_RuleDelayed) return 1500;
    if (head == SYM_DirectedEdge || head == SYM_UndirectedEdge) return 1500;
    if (head == SYM_Condition) return 1700;
    /* p?test, at the parser's PatternTest level: looser than application
     * and the atomic blank forms (x_, _h), tighter than every other infix. */
    if (head == SYM_PatternTest && e->data.function.arg_count == 2) return 7900;
    if (head == SYM_Alternatives) return 2300;
    if (head == SYM_StringExpression) return 2100;   /* a ~~ b, below `|` */
    if ((head == SYM_Repeated || head == SYM_RepeatedNull) &&
        e->data.function.arg_count == 1) return 2500;
    /* x:p (Pattern whose body is not a plain blank; x_ / x__ print atomically). */
    if (head == SYM_Pattern && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        !is_plain_blank(e->data.function.args[1])) return 2200;
    if (head == SYM_And) return 2800;
    if (head == SYM_Or) return 2700;
    /* !a, the parser's prefix Not level: looser than the comparisons
     * (!a == b is !(a == b)), tighter than And. */
    if (head == SYM_Not && e->data.function.arg_count == 1) return 3000;   /* looser than And, so Or nested in And parenthesises */
    if (head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_Inequality) return 3200;
    /* Plus[] / Times[] print as Plus[] / Times[]: atomic, never parenthesised */
    if ((head == SYM_Plus || head == SYM_Times) && e->data.function.arg_count == 0) return 9500;
    if (head == SYM_Plus) return 3500;

    if (head == SYM_Times) return 4500;
    /* NonCommutativeMultiply prints infix only with >= 2 args; the one-arg
     * bracketed form keeps the default 9500 so it is never spuriously wrapped. */
    if (head == SYM_NonCommutativeMultiply && e->data.function.arg_count >= 2) return 5900;
    /* a . b, the parser's Dot level: tighter than Times and Divide (so
     * (a b) . c and (1/2) . a parenthesise, and c a . b does not), looser than
     * Power (a^2 . b, (a . b)^2). Infix only with >= 2 args, like `**`. */
    if (head == SYM_Dot && e->data.function.arg_count >= 2) return 5300;
    if (head == SYM_Divide) return 5000;
    if (head == SYM_Power) {
        if (e->data.function.arg_count == 2) {
            Expr* exp = e->data.function.args[1];
            /* InputForm writes x^(-k), a power (x^(-1) as the reference) */
            if ((exp->type == EXPR_INTEGER && exp->data.integer < 0)
                || (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) < 0))
                return g_inputform_depth > 0 ? 6500 : 5000;
            int64_t n, d;
            if (is_rational(exp, &n, &d) && n < 0) return 5000;
            if (exp->type == EXPR_REAL && exp->data.real < 0.0) return 5000;
        }
        return 6500;
    }
    if (head == SYM_Rational && print_rational_is_numeric(e)) return 5000;
    if (head == SYM_Complex && print_complex_is_numeric(e)) {
        /* Only the bare imaginary unit `I` (Complex[0, 1]) is an atom.
         * `-I` (Complex[0, -1]) prints as a unary-minus form and so
         * takes Plus-level precedence — otherwise Power[-I, n] would
         * print as "-I^n" (i.e. -(I^n)) instead of "(-I)^n". General
         * Complex[a, b] also prints at Plus level. */
        if (e->data.function.arg_count == 2) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0)
                        || (re->type == EXPR_REAL && re->data.real == 0.0);
            if (re_zero && im->type == EXPR_INTEGER && im->data.integer == 1) {
                return 9500;
            }
        }
        return 3500;
    }
    /* Function[body] is the slot-based pure function rendered "body &".
     * Mathematica places "&" just above Postfix; on our scale that is 1000
     * (Postfix 800). */
    if (head == SYM_Function && e->data.function.arg_count == 1) return 1000;
    return 9500;
}

/* Print a string literal with surrounding quotes, re-escaping everything the
 * lexer decodes -- the control characters `\n`, `\t`, `\r` and the quote and
 * backslash themselves -- so the output is a faithful, re-parseable rendering
 * (Mathematica's InputForm: ToString["a\"b", InputForm] is "\"a\\\"b\""). */
/* Print uses OutputForm, in which a String is rendered as its CHARACTERS and
 * not as a quoted literal -- and recursively, so Print[{1, "a"}] is `{1, a}`.
 * Every other consumer of the printer (the REPL's Out[], ToString, FullForm,
 * InputForm, the LaTeX writer) keeps the quoted form, so this is a flag around
 * Print's own traversal rather than a change to the literal printer's contract.
 *
 * Without it every Print of a string in a script differed from Mathematica's,
 * which makes a side-by-side run of the same .m file in the two systems
 * needlessly hard to compare -- see docs/experiments/. */
static bool g_print_output_form = false;

/* Display rules (FormatValues, from Format[...] := ...) apply to what a user
 * sees -- ToString's OutputForm and the REPL's Out[] -- and never to the
 * internal renderings expr_to_string feeds back into evaluation, nor inside
 * InputForm or FullForm, nor in a script's Print (raw, as in wolframscript).
 * Set by those entry points only. */
static bool g_apply_formats = false;
static int g_format_depth = 0;
#define FORMAT_DEPTH_LIMIT 256

/* The display rewrite of e under its symbol's FormatValues, its right-hand
 * side evaluated; NULL when no rule applies.
 *
 * form == NULL (OutputForm / StandardForm): the first plain rule,
 * Format[f[...]] := rhs, stored as f[...] :> rhs.
 * form == InputForm / TeXForm: the first rule for that form,
 * Format[f[...], form] := rhs, stored as MakeBoxes[f[...], form] :>
 * Format[rhs, form]; the Format wrapper is stripped. (The generic
 * MakeBoxes[f[...], FormatType_] twin of a plain rule is never used here:
 * its FormatType is a pattern, not the form.) */
static Expr* format_value_rewrite(Expr* e, const char* form) {
    const Expr* k = e;
    while (k->type == EXPR_FUNCTION) k = k->data.function.head;
    if (k->type != EXPR_SYMBOL) return NULL;
    for (Rule* r = symtab_get_format_values(k->data.symbol.name); r; r = r->next) {
        const Expr* p = r->pattern;
        bool is_mb = p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
                     p->data.function.head->data.symbol.name == SYM_MakeBoxes;
        Expr* pat = r->pattern;
        if (!form) {
            if (is_mb) continue;
        } else {
            if (!is_mb || p->data.function.arg_count != 2) continue;
            const Expr* f = p->data.function.args[1];
            if (f->type != EXPR_SYMBOL || f->data.symbol.name != form) continue;
            pat = p->data.function.args[0];
        }
        MatchEnv* env = env_new();
        Expr* out = NULL;
        if (match(e, pat, env)) {
            Expr* bound = replace_bindings(r->replacement, env);
            /* A delayed rule's rhs /; test (also inside the Format[rhs, form]
             * wrapper a MakeBoxes rule carries) applies only when the test is
             * True: Format[e2[x_]] := "pos" /; x > 0 left e2[-1] showing
             * "pos /; -1 > 0" (audit Format#1). */
            Expr** slot = &bound;
            if (bound->type == EXPR_FUNCTION && bound->data.function.head->type == EXPR_SYMBOL &&
                bound->data.function.head->data.symbol.name == SYM_Format &&
                bound->data.function.arg_count >= 1)
                slot = &bound->data.function.args[0];
            bool guard_ok = true;
            while (r->delayed && guard_ok && (*slot)->type == EXPR_FUNCTION &&
                   (*slot)->data.function.head->type == EXPR_SYMBOL &&
                   (*slot)->data.function.head->data.symbol.name == SYM_Condition &&
                   (*slot)->data.function.arg_count == 2) {
                Expr* t = evaluate((*slot)->data.function.args[1]);
                guard_ok = t->type == EXPR_SYMBOL && t->data.symbol.name == SYM_True;
                expr_free(t);
                if (guard_ok) {
                    /* bound is fresh from replace_bindings, but its nodes may
                     * be shared: rebuild rather than write through. */
                    Expr* v = expr_copy((*slot)->data.function.args[0]);
                    if (slot == &bound) { expr_free(bound); bound = v; slot = &bound; }
                    else {
                        Expr* fa[2];
                        size_t na = bound->data.function.arg_count;
                        fa[0] = v;
                        fa[1] = na >= 2 ? expr_copy(bound->data.function.args[1]) : NULL;
                        Expr* nb = expr_new_function(expr_copy(bound->data.function.head), fa, na >= 2 ? 2 : 1);
                        expr_free(bound);
                        bound = nb;
                        slot = &bound->data.function.args[0];
                    }
                }
            }
            if (guard_ok) out = evaluate(bound);
            expr_free(bound);
        }
        env_free(env);
        if (out && form && out->type == EXPR_FUNCTION &&
            out->data.function.head->type == EXPR_SYMBOL &&
            out->data.function.head->data.symbol.name == SYM_Format &&
            out->data.function.arg_count >= 1) {
            Expr* inner = expr_copy(out->data.function.args[0]);
            expr_free(out);
            out = inner;
        }
        if (out) return out;
    }
    return NULL;
}

/* MakeBoxes[e, form]: a box rule the user gave for e's tag -- the tagged
 * h /: MakeBoxes[h[x_], fmt_] := boxes lands in h's FormatValues, which
 * evaluating MakeBoxes never consulted, so the call stayed unevaluated (audit
 * MakeBoxes#1).  Only rules written as MakeBoxes rules are used here: the
 * MakeBoxes[lhs, FormatType_] :> Format[rhs, FormatType] rule a Format
 * assignment installs is the printer's, not a box value.  NULL leaves the call. */
Expr* builtin_makeboxes(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    if (!symtab_format_value_count) return NULL;
    const Expr* k = res->data.function.args[0];
    while (k->type == EXPR_FUNCTION) k = k->data.function.head;
    if (k->type != EXPR_SYMBOL) return NULL;
    for (Rule* r = symtab_get_format_values(k->data.symbol.name); r; r = r->next) {
        const Expr* p = r->pattern;
        if (!(p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
              p->data.function.head->data.symbol.name == SYM_MakeBoxes)) continue;
        const Expr* rp = r->replacement;
        if (rp->type == EXPR_FUNCTION && rp->data.function.head->type == EXPR_SYMBOL &&
            rp->data.function.head->data.symbol.name == SYM_Format) continue;
        MatchEnv* env = env_new();
        Expr* out = NULL;
        if (match(res, r->pattern, env)) {
            Expr* bound = replace_bindings(r->replacement, env);
            out = r->delayed ? evaluate(bound) : expr_copy(bound);
            expr_free(bound);
        }
        env_free(env);
        if (out) return out;
    }
    return NULL;
}

/* Format[f[...], form] := rhs rewrite of e for a named form (FortranForm,
 * CForm, ...): the rule's right-hand side evaluated, or NULL.  Used by the
 * Fortran/C renderer (print_fortran.c). */
Expr* print_form_format_rewrite(Expr* e, const char* form) {
    if (!symtab_format_value_count || !e) return NULL;
    if (e->type != EXPR_SYMBOL && e->type != EXPR_FUNCTION) return NULL;
    return format_value_rewrite(e, form);
}

/* The box control character (\!, \(, \), \* ... as private-use code
 * points) starting at p, written as its escape -- the reference shows them
 * that way in every form: Print["a\(b\)c"] is a\(b\)c. Returns the bytes
 * consumed, 0 when p does not start one. */
static size_t print_box_escape(const char* p) {
    if ((unsigned char)p[0] != 0xEF) return 0;
    size_t n = 0, i = 0;
    while (n < 3 && p[n]) n++;
    const char* esc = lb_cp_escape(u8_decode(p, n, &i));
    if (!esc) return 0;
    fputs(esc, MTH_POUT);
    return i;
}

static void print_string_literal(const char* s) {
    /* An explicit InputForm[...] inside Print keeps the quoted literal. */
    if (g_print_output_form && g_inputform_depth == 0) {
        if (!strstr(s, "\xEF")) { fputs(s, MTH_POUT); return; }
        for (const char* p = s; *p; ) {
            size_t k = print_box_escape(p);
            if (k) { p += k; continue; }
            putchar((unsigned char)*p++);
        }
        return;
    }
    putchar('"');
    for (const char* p = s; *p; p++) {
        size_t k = print_box_escape(p);
        if (k) { p += k - 1; continue; }
        switch (*p) {
            case '\n': fputs("\\n", MTH_POUT); break;
            case '\t': fputs("\\t", MTH_POUT); break;
            case '\r': fputs("\\r", MTH_POUT); break;
            case '"':  fputs("\\\"", MTH_POUT); break;
            case '\\': fputs("\\\\", MTH_POUT); break;
            default:   putchar((unsigned char)*p); break;
        }
    }
    putchar('"');
}

/* Print a packed list as the nested List it stands for, STREAMING one element
 * at a time off the buffer. Materialising the whole tree first (as the visible
 * NDArray[...] path does) would allocate an Expr per element just to print --
 * fine for a value the user explicitly wrapped, unacceptable for an ordinary
 * List that happens to be packed and may hold 10^6 entries. One leaf is built
 * at a time so element formatting stays byte-identical to the List path
 * (Real formatting, Complex[a,b] -> "a + b I", FullForm quoting) rather than
 * being reimplemented here and drifting.
 *
 * `fullform` selects List[...] over {...}; the leaves follow suit. */
static void print_packed_level(Expr* a, int level, size_t* idx, bool fullform) {
    const NDArrayData* nd = &a->data.ndarray;
    if (level == nd->rank) {
        Expr* leaf = ndarray_buffer_element_to_expr(nd->data, (*idx)++, nd->dtype);
        if (fullform) expr_print_fullform(leaf);
        else          print_standard(leaf, 0);
        expr_free(leaf);
        return;
    }
    printf(fullform ? "List[" : "{");
    for (int64_t i = 0; i < nd->dims[level]; i++) {
        if (i) printf(", ");
        print_packed_level(a, level + 1, idx, fullform);
    }
    printf(fullform ? "]" : "}");
}

/* The directed infinities are symbols and products here (Infinity,
 * ComplexInfinity, Times[-1, Infinity]); FullForm shows the reference's
 * DirectedInfinity form, which reads back to the same value:
 * DirectedInfinity[1], DirectedInfinity[], DirectedInfinity[-1]. A product
 * qualifies when Infinity is one factor and every other factor is a numeric
 * quantity (the direction). */
/* > 0 while printing in FullForm proper (FullForm[...], ToString[e,
 * FullForm]); expr_print_fullform is also the atom printer of every other
 * form, where Infinity is Infinity. */
static int g_fullform_depth = 0;

static bool fullform_directed_infinity(Expr* e) {
    if (g_fullform_depth == 0) return false;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_Infinity) { printf("DirectedInfinity[1]"); return true; }
        if (e->data.symbol.name == SYM_ComplexInfinity) { printf("DirectedInfinity[]"); return true; }
        return false;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Times ||
        e->data.function.arg_count < 2) return false;
    size_t n = e->data.function.arg_count, inf = n;
    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Infinity) {
            if (inf != n) return false;
            inf = i;
        } else if (!expr_is_numeric_quantity(a)) {
            return false;
        }
    }
    if (inf == n) return false;
    printf("DirectedInfinity[");
    if (n == 2) {
        expr_print_fullform(e->data.function.args[inf == 0 ? 1 : 0]);
    } else {
        printf("Times[");
        bool first = true;
        for (size_t i = 0; i < n; i++) {
            if (i == inf) continue;
            if (!first) printf(", ");
            first = false;
            expr_print_fullform(e->data.function.args[i]);
        }
        printf("]");
    }
    printf("]");
    return true;
}

static void print_function_fullform(Expr* e) {
    if (fullform_directed_infinity(e)) return;
    expr_print_fullform(e->data.function.head);
    printf("[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i > 0) printf(", ");
        expr_print_fullform(e->data.function.args[i]);
    }    
    printf("]");
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/*  Arbitrary-precision reals                                          */
/* ------------------------------------------------------------------ */
/*
 * An MPFR real carries only a binary precision. N[x, n] (and the literal
 * x`n) stores ceil(n log2 10) bits; its decimal precision p is n when a
 * whole n maps to those bits (floor(bits log10 2) recovers it), else
 * bits log10 2 -- the value Precision[] reports (15.954589770191003 for the
 * 53 bits of a machine computation that left the double range).
 *
 * The digits shown are those above the value's uncertainty |x| 10^-p: with
 * L = log10|x|, the last digit sits at 10^(floor(L - p) + 1), so
 *     n = floor(L) - floor(L - p)
 * significant digits. A whole p gives n = p for every x; a fractional p
 * depends on the leading digits: 1.`16.5 shows 17 (1.0000000000000000),
 * 0.33333333333333333`16.5 shows 16, a 53-bit 10^400/3 shows 16
 * (3.333333333333333e+399), a 53-bit 0.95 shows 15 (0.950000000000000).
 * Every one of them is significant, so trailing zeros are printed --
 * N[1, 30] is 1.00000000000000000000000000000, not 1.0, which would claim
 * two digits.
 *
 * Display (OutputForm, ToString, Print, TeXForm), for |x| = d.ddd 10^e:
 *   - positional when -5 <= e < 6 and e < n (n = significant digits); the
 *     integer part is padded with zeros when the n digits do not reach the
 *     decimal point (99.99999`2 -> 100.);
 *   - otherwise d.ddd e+XX (TeXForm: d.ddd\times 10^{XX}), the mantissa
 *     still showing all n digits;
 *   e is decided before rounding to n digits, so a value that rounds up
 *   across a power of ten keeps its positional layout (9.99999`3 -> 10.0).
 *   A zero shows its accuracy a the same way: 0. for a <= 5, else 0.e-a.
 *
 * InputForm / FullForm: digits`p (or digits`p*^e for e < -5 or e >= 6)
 * with enough digits to name the binary value exactly (trailing zeros
 * trimmed: 1.`30.) and a precision mark p that the parser maps back to the
 * same number of bits, so the text re-reads to the identical number. A zero
 * is written by its accuracy, 0``a.
 */
#define MTH_LOG10_2 0.30102999566398119521

/* Whole decimal digits (>= 1) carried by a `bits`-bit value: the accuracy
 * shown for an inexact zero. `bits` is the mantissa width; the precision it
 * carries excludes the guard bits (numeric.h). */
static long mpfr_known_digits(mpfr_prec_t bits) {
    long d = (long)floor(numeric_bits_to_digits((long)bits) + 1e-9);
    return d < 1 ? 1 : d;
}

/* Significant digits (>= 1) displayed for a nonzero finite `bits`-bit value
 * whose leading digits, truncated, are D17 (at least one digit): see the
 * block comment above. With f = log10 of the mantissa d.ddd (0 <= f < 1),
 * floor(L) - floor(L - p) = -floor(f - p), which is p itself for a whole
 * precision. p is the precision of the width, not the width (guard bits). */
static long mpfr_display_digits(mpfr_prec_t bits, const char* D17) {
    double p = numeric_bits_to_digits((long)bits);
    if (p == floor(p) && p >= 1.0) return (long)p;      /* whole precision: n = p */
    double m = 0.0, scale = 1.0;
    for (const char* q = D17; *q; q++, scale /= 10.0) m += (*q - '0') * scale;
    double f = m > 1.0 ? log10(m) : 0.0;
    if (f >= 1.0) f = nextafter(1.0, 0.0);
    long n = -(long)floor(f - p);
    return n < 1 ? 1 : n;
}

/* The first n significant digits of |x| (x finite, nonzero), rounded with
 * `rnd`, as a malloc'd NUL-terminated digit string; *exp10 receives e with
 * |x| ~ d.ddd 10^e for those rounded digits. NULL on failure. */
static char* mpfr_sig_digits(const mpfr_t x, long n, mpfr_rnd_t rnd, long* exp10) {
    char* buf = NULL;
    int len = rnd == MPFR_RNDZ ? mpfr_asprintf(&buf, "%.*RZe", (int)(n - 1), x)
                               : mpfr_asprintf(&buf, "%.*RNe", (int)(n - 1), x);
    if (len < 0 || !buf) return NULL;
    char* out = malloc((size_t)len + 1);
    if (!out) { mpfr_free_str(buf); return NULL; }
    size_t k = 0;
    const char* p = buf;
    if (*p == '-' || *p == '+') p++;
    for (; *p && *p != 'e' && *p != 'E'; p++)
        if (isdigit((unsigned char)*p)) out[k++] = *p;
    out[k] = '\0';
    *exp10 = (*p == 'e' || *p == 'E') ? strtol(p + 1, NULL, 10) : 0;
    mpfr_free_str(buf);
    return out;
}

/* Lay digits D (nd of them) with decimal exponent e out positionally. */
static void mpfr_put_positional(const char* D, long nd, long e) {
    if (e >= 0) {
        for (long i = 0; i <= e; i++) putchar(i < nd ? D[i] : '0');
        putchar('.');
        for (long i = e + 1; i < nd; i++) putchar(D[i]);
    } else {
        printf("0.");
        for (long i = 0; i < -e - 1; i++) putchar('0');
        fputs(D, MTH_POUT);
    }
}

/* Mantissa d.ddd of a scientific layout. */
static void mpfr_put_mantissa(const char* D) {
    putchar(D[0]);
    putchar('.');
    fputs(D + 1, MTH_POUT);
}

/* Exponent suffix of the display form: *^50 / *^-5, as a machine real
 * prints, or the TeX \times 10^{e}. */
static void mpfr_put_display_exponent(long e, bool tex) {
    if (tex) {
        char t[32];
        snprintf(t, sizeof t, "%ld", e);
        if (strlen(t) == 1) printf("\\times 10^%s", t);
        else printf("\\times 10^{%s}", t);
    } else {
        /* *^e, the input syntax, so the displayed number re-reads
         * (numeric-tower review C9a); it was e+50 / e-05 */
        printf("*^%ld", e);
    }
}

/* A 53-bit MPFR is a machine-precision number with MPFR's exponent range
 * (the boundary policy's Exp[1000.], 1.*^400), and prints as a machine real
 * does (C9b): six digits for display, the shortest round-trip digits and no
 * precision mark for input -- parse.c reads such a literal back as the same
 * 53-bit value, since it lies outside the double range. %g's layout: an
 * exponent below -4 or at least the digit count is scientific. */
static void print_mpfr_machine(const mpfr_t x, bool input, bool tex) {
    long e = 0;
    char* D = NULL;
    if (input) {
        for (long n = 15; n <= 17; n++) {
            free(D);
            D = mpfr_sig_digits(x, n, MPFR_RNDN, &e);
            if (!D) break;
            /* does this digit string name x at 53 bits? */
            char t[64];
            snprintf(t, sizeof t, "%s%c.%se%ld", mpfr_sgn(x) < 0 ? "-" : "", D[0], D + 1, e);
            mpfr_t y; mpfr_init2(y, 53);
            mpfr_set_str(y, t, 10, MPFR_RNDN);
            bool same = mpfr_equal_p(x, y) != 0;
            mpfr_clear(y);
            if (same) break;
        }
    } else {
        D = mpfr_sig_digits(x, 6, MPFR_RNDN, &e);
    }
    if (!D) { printf("<mpfr print error>"); return; }
    long nd = (long)strlen(D);
    while (nd > 1 && D[nd - 1] == '0') D[--nd] = '\0';
    if (mpfr_sgn(x) < 0) putchar('-');
    long P = input ? 17 : 6;
    if (e < -4 || e >= P) {
        putchar(D[0]); putchar('.');
        fputs(nd > 1 ? D + 1 : "0", MTH_POUT);
        if (input) printf("*^%ld", e); else mpfr_put_display_exponent(e, tex);
    } else {
        mpfr_put_positional(D, nd, e);
        if (e >= nd - 1) putchar('0');          /* 100000. -> 100000.0 */
    }
    free(D);
}

/* Display form of an MPFR real (see the block comment above). */
static void print_mpfr_display(const mpfr_t x, bool tex) {
    if (mpfr_get_prec(x) == 53 && mpfr_number_p(x) && !mpfr_zero_p(x)) {
        print_mpfr_machine(x, false, tex);
        return;
    }
    if (!mpfr_number_p(x)) {            /* NaN / infinity: no digits to show */
        char* buf = NULL;
        if (mpfr_asprintf(&buf, "%Rg", x) >= 0 && buf) { fputs(buf, MTH_POUT); mpfr_free_str(buf); }
        return;
    }
    if (mpfr_zero_p(x)) {
        /* An inexact zero is known to absolute accuracy 10^-n. It has no
         * sign: an MPFR -0 prints as 0. too. */
        long n = mpfr_known_digits(mpfr_get_prec(x));
        printf("0.");
        if (n > 5) mpfr_put_display_exponent(-n, tex);
        return;
    }
    long eu = 0, er = 0;
    /* Truncated leading digits: the exponent before rounding (truncation
     * never carries into the next power of ten) and the mantissa that
     * decides how many digits a fractional precision shows. */
    char* T = mpfr_sig_digits(x, 17, MPFR_RNDZ, &eu);
    if (!T) { printf("<mpfr print error>"); return; }
    long n = mpfr_display_digits(mpfr_get_prec(x), T);
    char* D = mpfr_sig_digits(x, n, MPFR_RNDN, &er);
    if (!D) { free(T); printf("<mpfr print error>"); return; }
    if (mpfr_sgn(x) < 0) putchar('-');
    if (eu < -5 || eu >= 6 || eu >= n) {
        mpfr_put_mantissa(D);
        mpfr_put_display_exponent(er, tex);
    } else {
        mpfr_put_positional(D, (long)strlen(D), er);
    }
    free(T);
    free(D);
}

/* The precision (accuracy, for a zero) mark that parse.c reads back to
 * exactly `bits` bits: bits = numeric_digits_to_bits(mark). The whole number
 * of digits when it does (N[x, 30] -> 30.), else the shortest decimal that
 * does, checked through the same text -> double the parser performs. */
static void mpfr_put_precision_mark(mpfr_prec_t bits) {
    double x = numeric_bits_to_digits((long)bits);
    char t[48];
    if (x == floor(x) && x >= 1.0) {
        printf("%ld.", (long)x);
        return;
    }
    for (int k = 1; k <= 17; k++) {
        double scale = pow(10.0, k);
        snprintf(t, sizeof t, "%.*f", k, floor(x * scale) / scale);
        double m = strtod(t, NULL);
        if (m >= 1.0 && numeric_digits_to_bits(m) == (long)bits) {
            fputs(t, MTH_POUT);
            return;
        }
    }
    snprintf(t, sizeof t, "%.17g", x);
    fputs(t, MTH_POUT);
}

/* InputForm / FullForm text of an MPFR real, which re-reads to the same
 * value at the same precision. */
static void print_mpfr_input(const mpfr_t x) {
    mpfr_prec_t bits = mpfr_get_prec(x);
    if (!mpfr_number_p(x)) { print_mpfr_display(x, false); return; }
    if (bits == 53 && !mpfr_zero_p(x)) { print_mpfr_machine(x, true, false); return; }
    if (mpfr_zero_p(x)) {
        printf("0``");
        mpfr_put_precision_mark(bits);
        return;
    }
    /* 1 + ceil(bits log10 2) digits name every bits-bit value uniquely. */
    long n = 1 + (long)ceil((double)bits * MTH_LOG10_2);
    long e = 0;
    char* D = mpfr_sig_digits(x, n, MPFR_RNDN, &e);
    if (!D) { printf("<mpfr print error>"); return; }
    long nd = (long)strlen(D);
    while (nd > 1 && D[nd - 1] == '0') D[--nd] = '\0';
    if (mpfr_sgn(x) < 0) putchar('-');
    bool sci = e < -5 || e >= 6;
    if (sci) mpfr_put_mantissa(D);
    else mpfr_put_positional(D, nd, e);
    putchar('`');
    mpfr_put_precision_mark(bits);
    if (sci) printf("*^%ld", e);
    free(D);
}
#endif

/* InputForm/FullForm text of a machine real: the shortest decimal that reads
 * back to the same double, so the printed form re-parses to the same value
 * (OutputForm's %g keeps six digits and does not).  At most 17 significant
 * digits, which the parser reads as a machine real.  Exponents use the
 * input syntax *^. */
static void print_real_roundtrip(double d) {
    char buf[64];
    if (d == 0.0) d = 0.0;   /* a negative zero prints as 0.0 (it compares equal) */
    int digits = 15;
    for (; digits <= 17; digits++) {
        snprintf(buf, sizeof(buf), "%.*g", digits, d);
        if (strtod(buf, NULL) == d) break;
    }
    if (digits > 17) digits = 17;
    char* epos = strpbrk(buf, "eE");
    char mant[64], expo[16] = "";
    if (epos) {
        size_t n = (size_t)(epos - buf);
        memcpy(mant, buf, n); mant[n] = '\0';
        int ev = atoi(epos + 1);
        snprintf(expo, sizeof(expo), "*^%d", ev);
    } else {
        snprintf(mant, sizeof(mant), "%s", buf);
    }
    if (!strchr(mant, '.')) strcat(mant, ".0");   /* Mathilda writes 3.0, not 3. */
    /* No precision mark: a literal of up to 17 significant digits reads
     * back as a machine real (parse.c), as in the reference. */
    printf("%s%s", mant, expo);
}

/* Display (OutputForm) text of a machine real: six significant digits, as
 * %g picks them, ALWAYS with a decimal point and with the input syntax *^ for
 * an exponent, so the text re-reads as a Real (numeric-tower review C9a:
 * 1e+15 read back as the Integer-looking 1e+15 symbol product, 2.22045e-16 as
 * 2.22045 e - 16). 1.*^15 prints 1.0*^15, 100000. prints 100000.0. A negative
 * zero prints as 0.0 (it is SameQ to 0.0). */
static void print_real_display(double d) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%g", d == 0.0 ? 0.0 : d);
    char* ep = strpbrk(buf, "eE");
    if (!isfinite(d)) { fputs(buf, MTH_POUT); return; }
    if (ep && atoi(ep + 1) == -5) {
        /* %g turns scientific below 10^-4; the display (as the reference's,
         * and the MPFR display's) stays positional down to 10^-5:
         * 0.000012, not 1.2*^-5. The same six significant digits. */
        *ep = '\0';
        const char* m = buf;
        if (*m == '-') { fputc('-', MTH_POUT); m++; }
        fputs("0.0000", MTH_POUT);
        for (; *m; m++) if (*m != '.') fputc(*m, MTH_POUT);
    } else if (ep) {
        int ev = atoi(ep + 1);
        *ep = '\0';
        fputs(buf, MTH_POUT);
        if (!strchr(buf, '.')) fputs(".0", MTH_POUT);
        printf("*^%d", ev);
    } else {
        fputs(buf, MTH_POUT);
        if (!strchr(buf, '.')) fputs(".0", MTH_POUT);
    }
}

void expr_print_fullform(Expr* e) {
    if (!e) { printf("Null"); return; }
    switch (e->type) {
        case EXPR_INTEGER: printf("%" PRId64, e->data.integer); break;
        case EXPR_REAL: {
            if (g_inputform_depth > 0 && isfinite(e->data.real)) {
                print_real_roundtrip(e->data.real);
                break;
            }
            /* -0.0 prints as 0.0: it is SameQ to 0.0, and printing the sign
             * made exact-zero products like 0 (-1.2) look negative */
            print_real_display(e->data.real);
            break;
        }
        case EXPR_SYMBOL:
            if (fullform_directed_infinity(e)) break;
            printf("%s", context_display_name(e->data.symbol.name));
            break;
        case EXPR_STRING: print_string_literal(e->data.string); break;
        case EXPR_FUNCTION: print_function_fullform(e); break;
        case EXPR_NDARRAY: {
            /* A packed list must be indistinguishable from the nested List it
             * stands for, in every form including FullForm -- so List[...],
             * with no hint that the storage is a buffer. */
            if (is_packed_list(e)) {
                size_t idx = 0;
                print_packed_level(e, 0, &idx, true);
                break;
            }
            /* A visible NDArray[...] is always wrapped so it can never be
             * mistaken for a bare nested List, even in FullForm/InputForm. The
             * inner braces use the standard {..} rendering (not a further
             * List[...] FullForm expansion) so the payload stays readable. */
            printf("NDArray[");
            Expr* nested = ndarray_to_nested_list(e);
            print_standard(nested, 0);
            expr_free(nested);
            printf("]");
            break;
        }
        case EXPR_SPARSE: {
            /* FullForm (and InputForm) print the canonical constructor form,
             * which reads back to the same array. */
            Expr* canon = sp_canonical_expr(e->data.sparse);
            expr_print_fullform(canon);
            expr_free(canon);
            break;
        }
        case EXPR_PVEC:
        case EXPR_LINKEDLIST: {
            /* Payload kinds with no literal that reads back print the summary
             * Name[<n>], the name and length coming from the ds_ops row. */
            const DsOps* o = ds_ops(e->type);
            printf("%s[<%" PRId64 ">]", o->name, o->length(e));
            break;
        }
        case EXPR_COMPILED: {
            const CompiledFunction* cf = e->data.compiled;
            size_t n = compiled_function_num_args(cf);
            const char* const* names = compiled_function_arg_names(cf);
            printf("CompiledFunction[{");
            for (size_t i = 0; i < n; i++) { if (i) printf(", "); printf("%s", names[i] ? names[i] : "?"); }
            printf("}, ");
            print_standard((Expr*)compiled_function_body(cf), 0);
            printf("]");
            break;
        }
        case EXPR_BIGINT: {
            char* str = mpz_get_str(NULL, 10, e->data.bigint);
            printf("%s", str);
            free(str);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR:
            /* Display: every known digit, trailing zeros included;
             * InputForm/FullForm: the re-readable digits`prec spelling. */
            if (g_inputform_depth > 0) print_mpfr_input(e->data.mpfr);
            else print_mpfr_display(e->data.mpfr, false);
            break;
#endif
    }
}

#if UP_VALUES
/* Render the own/up/down values of `name` as assignment lines, one per rule,
 * the way Definition[name] displays them (in that order). Each rule records
 * whether a delayed assignment made it, so it prints with `=` or `:=`. An
 * upvalue prints as `lhs ^= rhs` / `lhs ^:= rhs` when UpSet would have tagged
 * `name` alone, else as `name /: lhs = rhs` / `name /: lhs := rhs`.
 * Empty -> "Null". */
/* One rule as `[tag /: ]lhs op rhs`, a guard written on the right of the
 * assignment shown back on the right (symtab_rule_display_parts). */
static void print_rule_parts(const Rule* r, const char* op, const char* tag) {
    Expr *lhs, *rhs;
    symtab_rule_display_parts(r, &lhs, &rhs);
    if (tag) printf("%s /: ", tag);
    print_standard(lhs, 0);
    printf("%s", op);
    print_standard(rhs, 0);
    expr_free(lhs);
    expr_free(rhs);
}

/* The N, Format and Default values of `name`, as Definition shows them (audit
 * Definition#3): N[f] = 2.5 (stored as N[f, {MachinePrecision,
 * MachinePrecision}]), Format[f[x_]] := "F" and Format[f[x_], TeXForm] :=
 * ... (stored as the display rule, MakeBoxes twins of a plain Format being
 * implementation detail), Default[f] = 0 (stored on Default).  `any` counts
 * the lines printed so far. */
static void print_definition_extra(const char* name, int* any) {
    for (Rule* r = symtab_get_n_values(name); r; r = r->next) {
        Expr *lhs, *rhs;
        symtab_rule_display_parts(r, &lhs, &rhs);
        Expr* shown = lhs;
        if (lhs->type == EXPR_FUNCTION && lhs->data.function.arg_count == 2) {
            Expr* p = lhs->data.function.args[1];
            if (head_is(p, SYM_List) && p->data.function.arg_count == 2 &&
                p->data.function.args[0]->type == EXPR_SYMBOL &&
                p->data.function.args[0]->data.symbol.name == SYM_MachinePrecision &&
                p->data.function.args[1]->type == EXPR_SYMBOL &&
                p->data.function.args[1]->data.symbol.name == SYM_MachinePrecision) {
                Expr* a0 = expr_copy(lhs->data.function.args[0]);
                shown = expr_new_function(expr_copy(lhs->data.function.head), &a0, 1);
            }
        }
        if ((*any)++) printf("\n");
        print_standard(shown, 0);
        printf("%s", r->delayed ? " := " : " = ");
        print_standard(rhs, 0);
        if (shown != lhs) expr_free(shown);
        expr_free(lhs);
        expr_free(rhs);
    }
    for (Rule* r = symtab_get_format_values(name); r; r = r->next) {
        Expr *lhs, *rhs;
        symtab_rule_display_parts(r, &lhs, &rhs);
        Expr* target = lhs;
        Expr* form = NULL;
        bool skip = false;
        if (head_is(target, SYM_HoldPattern) && target->data.function.arg_count == 1)
            target = target->data.function.args[0];
        if (head_is(target, SYM_MakeBoxes) && target->data.function.arg_count == 2) {
            form = target->data.function.args[1];
            target = target->data.function.args[0];
            /* MakeBoxes[lhs, FormatType_]: the twin of a plain Format. */
            if (form->type != EXPR_SYMBOL) skip = true;
            /* its rhs is Format[v, form]: show v */
            if (head_is(rhs, SYM_Format) && rhs->data.function.arg_count == 2) {
                Expr* v = expr_copy(rhs->data.function.args[0]);
                expr_free(rhs);
                rhs = v;
            }
        }
        if (!skip) {
            Expr* fa[2] = { expr_copy(target), form ? expr_copy(form) : NULL };
            Expr* f = expr_new_function(expr_new_symbol(SYM_Format), fa, form ? 2 : 1);
            if ((*any)++) printf("\n");
            /* the rule itself, not its effect: no display rules while printing it */
            bool saved_af = g_apply_formats;
            g_apply_formats = false;
            print_standard(f, 0);
            printf("%s", r->delayed ? " := " : " = ");
            print_standard(rhs, 0);
            g_apply_formats = saved_af;
            expr_free(f);
        }
        expr_free(lhs);
        expr_free(rhs);
    }
    for (Rule* r = symtab_get_down_values(SYM_Default); r; r = r->next) {
        if (!symtab_default_value_rule(r, name)) continue;
        if ((*any)++) printf("\n");
        print_rule_parts(r, r->delayed ? " := " : " = ", NULL);
    }
}

static bool has_definition_extra(const char* name) {
    if (symtab_get_n_values(name) || symtab_get_format_values(name)) return true;
    for (Rule* r = symtab_get_down_values(SYM_Default); r; r = r->next)
        if (symtab_default_value_rule(r, name)) return true;
    return false;
}

static void print_definition_body(const char* name) {
    Rule* ov = symtab_get_own_values(name);
    Rule* dv = symtab_get_down_values(name);
    Rule* uv = symtab_get_up_values(name);
    int any = 0;
    /* The attributes lead and the options close the definition, each set off
     * by a blank line, as the reference shows them (audit Definition#1,
     * Information#1): Attributes[f] = {Listable} ... Options[f] = {...}. */
    if (get_attributes(name) != ATTR_NONE) {
        Expr* a[1] = { expr_new_symbol(name) };
        Expr* call = expr_new_function(expr_new_symbol("Attributes"), a, 1);
        Expr* lst = builtin_attributes(call);
        expr_free(call);
        if (lst) {
            printf("Attributes[%s] = ", context_display_name(name));
            print_standard(lst, 0);
            expr_free(lst);
            any++;
        }
    }
    bool rules = ov || dv || uv || symtab_get_sub_values(name) || has_definition_extra(name);
    if (any && rules) printf("\n\n");
    int nattr = any;
    any = 0;
    for (Rule* r = ov; r; r = r->next) {
        if (any++) printf("\n");
        print_rule_parts(r, r->delayed ? " := " : " = ", NULL);
    }
    for (Rule* r = uv; r; r = r->next) {
        if (any++) printf("\n");
        if (upvalue_renders_as_upset(name, r->pattern))
            print_rule_parts(r, r->delayed ? " ^:= " : " ^= ", NULL);
        else
            print_rule_parts(r, r->delayed ? " := " : " = ", name);
    }
    for (Rule* r = dv; r; r = r->next) {
        if (any++) printf("\n");
        print_rule_parts(r, r->delayed ? " := " : " = ", NULL);
    }
    for (Rule* r = symtab_get_sub_values(name); r; r = r->next) {   /* f[x_][y_] := ... */
        if (any++) printf("\n");
        print_rule_parts(r, r->delayed ? " := " : " = ", NULL);
    }
    print_definition_extra(name, &any);
    Expr* opts = symtab_get_options(name);            /* borrowed */
    if (opts && opts->type == EXPR_FUNCTION && opts->data.function.arg_count > 0) {
        if (any || nattr) printf("\n\n");
        printf("Options[%s] = ", context_display_name(name));
        print_standard(opts, 0);
        any++;
    }
    if (!any && !nattr) printf("Null");
}
#endif

/* A rational coefficient p/q in a product that already has a denominator
 * folds q into that denominator, as Mathematica does: Times[-1/2, a, 1/h]
 * prints "-a/(2 h)", not "-(1/2 a)/h", and Times[1/2, 1/h] prints
 * "1/(2 h)", not "1/2/h". Shared by the plain and TeX Times printers, which
 * both strip the sign first, so the coefficient (if any) is num[0] and is
 * positive. num/den own their entries and have room for one more each. */
static void fold_rational_coeff_into_den(Expr** num, size_t* nnum,
                                         Expr** den, size_t* nden) {
    if (*nden == 0 || *nnum == 0) return;
    Expr* r = num[0];
    if (r->type != EXPR_FUNCTION || r->data.function.head->type != EXPR_SYMBOL ||
        r->data.function.head->data.symbol.name != SYM_Rational ||
        r->data.function.arg_count != 2) return;
    Expr* p = r->data.function.args[0];
    memmove(den + 1, den, sizeof(Expr*) * *nden);
    den[0] = expr_copy(r->data.function.args[1]);
    (*nden)++;
    if (p->type == EXPR_INTEGER && p->data.integer == 1 && *nnum > 1) {
        memmove(num, num + 1, sizeof(Expr*) * (*nnum - 1));
        (*nnum)--;
    } else {
        num[0] = expr_copy(p);
    }
    expr_free(r);
}

/* ------------------------------------------------------------------ */
/*  Times under InputForm                                              */
/* ------------------------------------------------------------------ */

/* The reference's InputForm of a product (ToString[1/100 a + 5 x, InputForm]
 * is "a/100 + 5*x"):
 *   - factors joined by an explicit "*";
 *   - a rational coefficient p/q goes over the denominator, the numerator
 *     parenthesised when it is a product or signed: x/2, (3*x)/4, (-3*x)/4,
 *     Pi/4, (a + b)/2, (a*b)/(c*d);
 *   - a lone reciprocal with an integer exponent reads as a power: x^(-1),
 *     -x^(-1), x^(-2) (2/x, a/b and 1/Sqrt[x] keep the quotient);
 *   - a negated product is -(a*b), a negated single factor -x, and a signed
 *     coefficient -2*x. */
static void print_factor_list(Expr** f, size_t n, bool paren_if_many) {
    bool paren = paren_if_many && n > 1;
    if (paren) printf("(");
    for (size_t i = 0; i < n; i++) {
        if (i > 0) printf("*");
        print_standard(f[i], 4501);
    }
    if (paren) printf(")");
}

/* The factors after the coefficient, printed as InputForm prints a product
 * with no coefficient: a*b, x/y, 1/(x*y), x^(-1).  `wrap` parenthesises a
 * product or quotient of more than one factor (after a unary minus). */
static void print_mpz(const mpz_t z) {
    char* t = mpz_get_str(NULL, 10, z);
    printf("%s", t);
    free(t);
}

/* Is e a (held) product that prints as a quotient -- a factor with a
 * negative numeric exponent?  As the lone numerator of a quotient it needs no
 * parentheses, division being left-associative: HoldForm[(a/b)/c] is a/b/c. */
static bool times_is_quotient(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Times) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* f = e->data.function.args[i];
        if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
            f->data.function.head->data.symbol.name == SYM_Power && f->data.function.arg_count == 2) {
            const Expr* x = f->data.function.args[1];
            int64_t n, d;
            if ((x->type == EXPR_INTEGER && x->data.integer < 0) || (is_rational((Expr*)x, &n, &d) && n < 0))
                return true;
        }
    }
    return false;
}

static void print_inputform_rest(Expr** num, size_t nn, Expr** den, size_t nd,
                                 bool den_int_power, int wrap) {
    if (nn == 0 && nd == 1 && den_int_power && wrap != 2) {   /* x^(-k) */
        Expr* b = den[0];
        Expr* base = b; Expr* ex = NULL;
        if (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL
            && b->data.function.head->data.symbol.name == SYM_Power) {
            base = b->data.function.args[0]; ex = b->data.function.args[1];
        }
        print_standard(base, 6501);
        printf("^(-");
        /* The exponent sits under a unary minus, so a sum must keep its
         * parentheses: (q^(n + 1))^(-1) is q^(-(n + 1)).  At precedence 0 it
         * printed q^(-n + 1), which re-reads as a different expression. */
        if (ex) print_standard(ex, 4501); else printf("1");
        printf(")");
        return;
    }
    if (nd == 0) {
        bool paren = wrap && nn > 1;
        if (paren) printf("(");
        print_factor_list(num, nn, false);
        if (paren) printf(")");
        return;
    }
    /* wrap == 1: around every quotient (after "-"); wrap == 2: only a
     * product (after "-1/q*", where a quotient stays bare) */
    bool paren = wrap == 1;
    if (paren) printf("(");
    if (nn == 0) printf("1");
    else if (nn == 1 && times_is_quotient(num[0])) print_standard(num[0], 4500);   /* a/b/c */
    else print_factor_list(num, nn, true);
    printf("/");
    if (nd > 1) print_factor_list(den, nd, true);
    else print_standard(den[0], 5001);
    if (paren) printf(")");
}

static void print_times_inputform(Expr* e) {
    size_t count = e->data.function.arg_count;
    /* A held, non-canonical product -- a literal factor 1, or two numeric
     * factors (HoldForm[1*x], HoldForm[2*3]) -- prints literally, as written:
     * folding its coefficients would show HoldForm[1*x] as x. An evaluated
     * Times never has either shape. */
    {
        size_t nnum = 0;
        bool has_one = false;
        for (size_t i = 0; i < count; i++) {
            Expr* a = e->data.function.args[i];
            if (a->type == EXPR_INTEGER || a->type == EXPR_BIGINT || a->type == EXPR_REAL) nnum++;
            if (a->type == EXPR_INTEGER && a->data.integer == 1) has_one = true;
        }
        if (count >= 2 && (has_one || nnum >= 2)) {
            for (size_t i = 0; i < count; i++) {
                if (i) printf("*");
                print_standard(e->data.function.args[i], 4501);
            }
            return;
        }
    }
    mpz_t pn, pd;
    mpz_init_set_si(pn, 1);
    mpz_init_set_si(pd, 1);
    Expr* coef_other = NULL;                 /* a Real or Complex coefficient */
    Expr** num = malloc(sizeof(Expr*) * (count + 1));
    Expr** den = malloc(sizeof(Expr*) * (count + 1));
    size_t nn = 0, nd = 0;
    bool den_int_power = true;               /* every denominator exponent an integer */
    for (size_t i = 0; i < count; i++) {
        Expr* a = e->data.function.args[i];
        int64_t rn, rd;
        if (i == 0 && a->type == EXPR_INTEGER) { mth_mpz_set_i64(pn, a->data.integer); continue; }
        if (i == 0 && a->type == EXPR_BIGINT) { mpz_set(pn, a->data.bigint); continue; }
        if (i == 0 && is_rational(a, &rn, &rd)) { mth_mpz_set_i64(pn, rn); mth_mpz_set_i64(pd, rd); continue; }
        if (i == 0 && (a->type == EXPR_REAL
                       || (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
                           && a->data.function.head->data.symbol.name == SYM_Complex))) {
            coef_other = a; continue;
        }
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol.name == SYM_Power
            && a->data.function.arg_count == 2) {
            Expr* x = a->data.function.args[1];
            int64_t xn, xd;
            if (x->type == EXPR_INTEGER && x->data.integer < 0) {
                den[nd++] = x->data.integer == -1 ? expr_copy(a->data.function.args[0])
                    : expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy(a->data.function.args[0]),
                                     print_neg_int(x->data.integer) }, 2);
                continue;
            }
            if (is_rational(x, &xn, &xd) && xn < 0) {
                den[nd++] = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(a->data.function.args[0]), make_rational_neg(xn, xd) }, 2);
                den_int_power = false;
                continue;
            }
            Expr* nx = print_neg_symbolic_exp(x);
            if (nx) {
                den[nd++] = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(a->data.function.args[0]), nx }, 2);
                continue;
            }
        }
        num[nn++] = expr_copy(a);
    }
    bool neg = mpz_sgn(pn) < 0;
    if (neg) mpz_neg(pn, pn);
    bool p_one = mpz_cmp_ui(pn, 1) == 0, q_one = mpz_cmp_ui(pd, 1) == 0;

    if (coef_other && nd > 0) {
        /* a real or complex coefficient is part of the numerator: (1.5*x)/y,
         * 2.5/x, ((2*I)*x)/y, I/y */
        bool np = nn > 0;
        if (np) printf("(");
        print_standard(coef_other, 4500);
        for (size_t i = 0; i < nn; i++) { printf("*"); print_standard(num[i], 4501); }
        if (np) printf(")");
        printf("/");
        if (nd > 1) print_factor_list(den, nd, true);
        else print_standard(den[0], 5001);
    } else if (coef_other) {
        /* (2*I)*Pi, (I/2)*x, 1.5*x */
        print_standard(coef_other, 4500);
        if (nn + nd > 0) {
            printf("*");
            /* a product after the coefficient is not parenthesised
             * ((2*I)*a*b, 0.1*d*x), a quotient keeps its form */
            if (nd == 0) print_factor_list(num, nn, false);
            else print_inputform_rest(num, nn, den, nd, den_int_power, 2);
        }
    } else if (p_one && q_one) {
        /* coefficient +-1: x^(-1), a*b, -x, -(a*b), -(1/(x*y)) */
        bool single = (nd == 0 && nn == 1) || (nn == 0 && nd == 1 && den_int_power);
        if (neg) printf("-");
        print_inputform_rest(num, nn, den, nd, den_int_power, neg && !single ? 1 : 0);
    } else if (neg && p_one) {
        /* -1/q*rest: -1/2*a, -1/2*(a*b), -1/2*x/y, -1/2*1/x^2 */
        printf("-1/");
        print_mpz(pd);
        printf("*");
        print_inputform_rest(num, nn, den, nd, den_int_power, 2);
    } else {
        /* (p*num)/(q*den) with the sign inside the numerator */
        size_t k = (p_one && nn > 0) ? 0 : 1;
        bool np = (k + nn) > 1 && (nd > 0 || !q_one);
        if (np) printf("(");
        if (neg) printf("-");
        if (k) print_mpz(pn);
        for (size_t i = 0; i < nn; i++) {
            if (k || i > 0) printf("*");
            print_standard(num[i], 4501);
        }
        if (np) printf(")");
        if (!q_one || nd > 0) {
            printf("/");
            size_t m = (q_one ? 0 : 1) + nd;
            if (m > 1) printf("(");
            if (!q_one) print_mpz(pd);
            for (size_t i = 0; i < nd; i++) {
                if (!q_one || i > 0) printf("*");
                print_standard(den[i], m > 1 ? 4501 : 5001);
            }
            if (m > 1) printf(")");
        }
    }
    for (size_t i = 0; i < nn; i++) expr_free(num[i]);
    for (size_t i = 0; i < nd; i++) expr_free(den[i]);
    free(num); free(den);
    mpz_clear(pn); mpz_clear(pd);
}

/* InputForm of a Complex whose imaginary part is a non-integer rational p/q,
 * as the reference writes it: I/2, (2*I)/3, -1/2*I, (-2*I)/3, and with a real
 * part 1/2 + I/3, 1/2 - I/2.  False (nothing printed) for other shapes. */
static bool print_complex_rational_im_inputform(Expr* e) {
    Expr* re = e->data.function.args[0];
    Expr* im = e->data.function.args[1];
    int64_t p, q;
    if (!is_rational(im, &p, &q) || q == 1) return false;
    bool re_zero = re->type == EXPR_INTEGER && re->data.integer == 0;
    if (!re_zero && re->type != EXPR_INTEGER && !is_rational(re, &(int64_t){0}, &(int64_t){0}))
        return false;
    bool neg = p < 0;
    int64_t ap = neg ? -p : p;
    if (!re_zero) {
        print_standard(re, 3500);
        printf(neg ? " - " : " + ");
        neg = false;
    }
    if (neg && ap == 1) { printf("-1/%lld*I", (long long)q); return true; }
    if (ap == 1) printf("I/%lld", (long long)q);
    else printf("(%s%lld*I)/%lld", neg ? "-" : "", (long long)ap, (long long)q);
    return true;
}

static void print_standard(Expr* e, int parent_prec) {
    if (!e) { printf("Null"); return; }
    /* Format[...] display rules. A rewrite is itself displayed (and so
     * formatted again: Format[a[x_]] := b[x], Format[b[x_]] := 9 shows 9); a
     * rule that keeps reproducing a match stops at FORMAT_DEPTH_LIMIT, as the
     * reference stops at $RecursionLimit with Format::forml. */
    /* Inside InputForm only Format[..., InputForm] rules apply. */
    if (g_apply_formats && symtab_format_value_count &&
        (e->type == EXPR_SYMBOL || e->type == EXPR_FUNCTION)) {
        if (g_format_depth >= FORMAT_DEPTH_LIMIT) {
            if (g_format_depth == FORMAT_DEPTH_LIMIT)
                mth_message("Format", "forml",
                            "Recursion depth limit $RecursionLimit exceeded while formatting output.");
            g_format_depth++;   /* report once per print */
        } else {
            Expr* shown = format_value_rewrite(e, g_inputform_depth > 0 ? SYM_InputForm : NULL);
            if (shown) {
                g_format_depth++;
                print_standard(shown, parent_prec);
                g_format_depth--;
                expr_free(shown);
                return;
            }
        }
    }
    /* OutputForm (Print, ToString): Style[expr, opts...] is presentation only
     * -- with no front end to style it, it prints as expr. */
    if (g_print_output_form && e->type == EXPR_FUNCTION && e->data.function.arg_count >= 1 &&
        e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol.name == SYM_Style) {
        print_standard(e->data.function.args[0], parent_prec);
        return;
    }
    /* Intercept before the delegation below: a packed list is not an
     * EXPR_FUNCTION, so it would otherwise fall through to
     * expr_print_fullform and come out as List[...] instead of {...}. */
    if (is_packed_list(e)) {
        size_t idx = 0;
        print_packed_level(e, 0, &idx, false);
        return;
    }
    if (e->type == EXPR_SPARSE) {
        /* InputForm: the canonical SparseArray[Automatic, dims, bg, {...}]
         * form, which reads back. Otherwise the summary SparseArray[<n>, dims]
         * (n explicit entries), as the array is shown in output. */
        const SparseArray* sa = e->data.sparse;
        if (g_inputform_depth > 0) {
            Expr* canon = sp_canonical_expr(sa);
            print_standard(canon, parent_prec);
            expr_free(canon);
            return;
        }
        printf("SparseArray[<%" PRId64 ">, {", sp_nnz(sa));
        for (int k = 0; k < sa->rank; k++) printf(k ? ", %" PRId64 : "%" PRId64, sa->dims[k]);
        printf("}]");
        return;
    }
    if (e->type != EXPR_FUNCTION) {
        /* Inside a NumberForm[...] wrapper, reformat every numeric leaf; fall
         * through to the default printer for anything the engine declines
         * (non-finite reals, non-numeric atoms). */
        if (g_numberform_ctx && nf_reformattable_atom(e)
            && numberform_render_number(e, g_numberform_ctx))
            return;
        expr_print_fullform(e);
        return;
    }

    int my_prec = get_expr_prec(e);
    bool need_parens = (my_prec < parent_prec);

    if (need_parens) printf("(");

    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol.name;
        
        /* Plus[] / Times[] have no operator spelling: printing nothing made
         * Hold[{Times[], 2}] read back as Hold[{, 2}] (audit InputForm#3). */
        if ((head == SYM_Plus || head == SYM_Times) && e->data.function.arg_count == 0) {
            printf("%s[]", head);
        }
        else if (head == SYM_SequenceForm) {
            /* SequenceForm[a, b, ...] shows its arguments side by side: strings
             * as their characters, except in InputForm, where each argument
             * keeps its input form ("a"1 for SequenceForm["a", 1]). */
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* a = e->data.function.args[i];
                if (a->type == EXPR_STRING && g_inputform_depth == 0) fputs(a->data.string, MTH_POUT);
                else print_standard(a, 0);
            }
        }
        else if (head == SYM_FullForm && e->data.function.arg_count == 1) {
            /* FullForm keeps string quotes even inside Print / ToString. */
            g_inputform_depth++;
            g_fullform_depth++;
            expr_print_fullform(e->data.function.args[0]);
            g_fullform_depth--;
            g_inputform_depth--;
        }
        else if (head == SYM_TeXForm && e->data.function.arg_count == 1) {
            print_tex(e->data.function.args[0], 0);
        }
        else if ((head == SYM_FortranForm || head == SYM_CForm) &&
                 e->data.function.arg_count == 1) {
            char* txt = fortran_form_string(e->data.function.args[0],
                                            head == SYM_CForm ? FORM_LANG_C : FORM_LANG_FORTRAN, 0);
            printf("%s", txt);
            free(txt);
        }
        else if (head == SYM_LeanForm && e->data.function.arg_count == 1) {
            char* txt = lean_form_string(e->data.function.args[0]);
            printf("%s", txt);
            free(txt);
        }
        else if (head == SYM_InputForm && e->data.function.arg_count == 1) {
            g_inputform_depth++;
            print_standard(e->data.function.args[0], parent_prec);
            g_inputform_depth--;
        }
        else if (head == SYM_SeriesData && e->data.function.arg_count == 6 && g_inputform_depth == 0) {
            print_series_data(e, parent_prec);
        }
        else if (head == SYM_HoldForm && e->data.function.arg_count == 1) {
            /* Invisible in display, but InputForm keeps the wrapper so the
             * text re-reads as the held expression: HoldForm[1 + 1]. */
            if (g_inputform_depth > 0) {
                printf("HoldForm[");
                print_standard(e->data.function.args[0], 0);
                printf("]");
            } else {
                print_standard(e->data.function.args[0], parent_prec);
            }
        }
        /* TraditionalForm/StandardForm only pick a TYPESET notation (see
         * print_latex.c); as text they are invisible, like HoldForm -- except
         * under InputForm, which must keep them so the result re-reads. */
        else if ((head == SYM_TraditionalForm || head == SYM_StandardForm)
                 && e->data.function.arg_count == 1 && g_inputform_depth == 0) {
            print_standard(e->data.function.args[0], parent_prec);
        }
        else if (head == SYM_NumberForm && e->data.function.arg_count >= 1
                 && g_inputform_depth == 0) {
            /* Install the formatting context and render the wrapped expression;
             * numeric leaves route through numberform_render_number. */
            NumberFormCtx ctx;
            if (numberform_build_ctx(e, &ctx)) {
                const NumberFormCtx* saved = g_numberform_ctx;
                g_numberform_ctx = &ctx;
                print_standard(e->data.function.args[0], parent_prec);
                g_numberform_ctx = saved;
            } else {
                print_standard(e->data.function.head, 9500);
                printf("[]");
            }
        }
        else if (head == SYM_Row && e->data.function.arg_count >= 1
                 && e->data.function.args[0]->type == EXPR_FUNCTION
                 && e->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                 && e->data.function.args[0]->data.function.head->data.symbol.name == SYM_List
                 && g_inputform_depth == 0) {
            /* Row[{a, b, ...}] concatenates the OutputForm of each element (no
             * quotes on strings); Row[{...}, s] inserts the string s between. */
            Expr* list = e->data.function.args[0];
            /* Any separator expression, printed like the elements: Row[{a, b},
             * x] is axb, not ab (audit Row#1). */
            Expr* sep = e->data.function.arg_count >= 2 ? e->data.function.args[1] : NULL;
            bool saved_of = g_print_output_form;
            g_print_output_form = true;
            for (size_t i = 0; i < list->data.function.arg_count; i++) {
                if (i > 0 && sep) {
                    if (sep->type == EXPR_STRING) fputs(sep->data.string, MTH_POUT);
                    else print_standard(sep, 0);
                }
                print_standard(list->data.function.args[i], 0);
            }
            g_print_output_form = saved_of;
        }
        /* OutputForm[e] is a display wrapper: shown as e in OutputForm, never
         * as the literal OutputForm[...] (audit Print#1). */
        else if (head == SYM_OutputForm && e->data.function.arg_count == 1
                 && g_inputform_depth == 0) {
            bool saved_of = g_print_output_form;
            g_print_output_form = true;
            print_standard(e->data.function.args[0], parent_prec);
            g_print_output_form = saved_of;
        }
        /* StringForm["tmpl", e1, ...] displays the template with its `` / `n`
         * slots filled (audit Print#2). */
        else if (head == SYM_StringForm && e->data.function.arg_count >= 1
                 && e->data.function.args[0]->type == EXPR_STRING
                 && g_inputform_depth == 0) {
            char* txt = mth_stringform_fill(e->data.function.args[0]->data.string,
                                            e->data.function.args + 1,
                                            e->data.function.arg_count - 1);
            if (txt) { fputs(txt, MTH_POUT); free(txt); }
        }
        else if (head == SYM_Graphics && e->data.function.arg_count >= 1 && g_inputform_depth == 0) {
            printf("-Graphics-");
        }
        else if (head == SYM_Graphics3D && e->data.function.arg_count >= 1 && g_inputform_depth == 0) {
            printf("-Graphics3D-");
        }
        else if ((head == SYM_DirectedEdge || head == SYM_UndirectedEdge)
                 && e->data.function.arg_count == 2) {
            /* u -> v (directed) / u <-> v (undirected). Infix in both Standard
             * and Input forms so the literal round-trips through the parser. */
            print_standard(e->data.function.args[0], my_prec);
            printf("%s", head == SYM_DirectedEdge ? " -> " : " <-> ");
            print_standard(e->data.function.args[1], my_prec);
        }
        else if (head == SYM_Graph && g_inputform_depth == 0 && graph_is_valid(e)) {
            /* Terse summary in standard output (Graph[verts, edges] plus any
             * EdgeWeight / EdgeCapacity options); the option lists are not
             * shown here. InputForm/FullForm fall through to the literal
             * Graph[{...}, {...}, opts...] constructor (round-trippable). */
            unsigned long nv = (unsigned long)e->data.function.args[0]->data.function.arg_count;
            unsigned long ne = (unsigned long)e->data.function.args[1]->data.function.arg_count;
            printf("Graph[<%lu %s, %lu %s>]",
                   nv, nv == 1 ? "vertex" : "vertices",
                   ne, ne == 1 ? "edge" : "edges");
        }
        else if (g_inputform_depth == 0 && e->data.function.arg_count == 2
                 && head == hyp_sym_hypergraph() && hypergraph_is_valid(e)) {
            /* Hypergraph terse summary; InputForm/FullForm print the literal. */
            unsigned long nv = (unsigned long)e->data.function.args[0]->data.function.arg_count;
            unsigned long ne = (unsigned long)e->data.function.args[1]->data.function.arg_count;
            printf("Hypergraph[<%lu %s, %lu %s>]",
                   nv, nv == 1 ? "vertex" : "vertices",
                   ne, ne == 1 ? "hyperedge" : "hyperedges");
        }
        else if (head == SYM_Rational && print_rational_is_numeric(e)) {
            /* Only a Rational NUMBER is p/q: a held Rational[a, b] with
             * symbolic parts prints as written (audit InputForm#4). */
            print_standard(e->data.function.args[0], 5000);
            printf("/");
            print_standard(e->data.function.args[1], 5000);
        }
        else if (head == SYM_Complex && print_complex_is_numeric(e)
                 && g_inputform_depth > 0
                 && print_complex_rational_im_inputform(e)) {
            /* printed */
        }
        else if (head == SYM_Complex && print_complex_is_numeric(e)) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0);
            bool im_one = (im->type == EXPR_INTEGER && im->data.integer == 1);
            bool im_minus_one = (im->type == EXPR_INTEGER && im->data.integer == -1);

            bool im_neg = im_minus_one;
            Expr* im_abs = NULL;
            int64_t rn, rd;
            if (!im_neg) {
                if (im->type == EXPR_INTEGER && im->data.integer < 0) {
                    im_neg = true;
                    im_abs = print_neg_int(im->data.integer);
                } else if (im->type == EXPR_REAL && im->data.real < 0.0) {
                    im_neg = true;
                    im_abs = expr_new_real(-im->data.real);
#ifdef USE_MPFR
                } else if (im->type == EXPR_MPFR && mpfr_sgn(im->data.mpfr) < 0) {
                    im_neg = true;
                    im_abs = expr_new_mpfr_copy(im->data.mpfr);
                    mpfr_neg(im_abs->data.mpfr, im_abs->data.mpfr, MPFR_RNDN);
#endif
                } else if (is_rational(im, &rn, &rd) && rn < 0) {
                    im_neg = true;
                    im_abs = make_rational_neg(rn, rd);
                }
            }

            if (!re_zero) {
                print_standard(re, 3500);
                printf(im_neg ? " - " : " + ");
            } else if (im_neg) {
                printf("-");
            }
            if (im_one || im_minus_one) printf("I");
            else if (im_abs) { print_standard(im_abs, 4500); printf("*I"); expr_free(im_abs); }
            else { print_standard(im, 4500); printf("*I"); }
        }
        else if (head == SYM_Power && e->data.function.arg_count == 2) {
            Expr* exp = e->data.function.args[1];
            int64_t n, d;
            bool is_negative = false;
            Expr* pos_exp = NULL;
            if ((pos_exp = print_neg_int_exp(exp)) != NULL) {
                is_negative = true;
            } else if (is_rational(exp, &n, &d) && n < 0) {
                is_negative = true;
                pos_exp = make_rational_neg(n, d);
            } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                is_negative = true;
                pos_exp = expr_new_real(-exp->data.real);
            }

            if (is_negative && g_inputform_depth > 0
                && (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT)) {
                /* InputForm: x^(-1), x^(-2), as the reference writes them */
                print_standard(e->data.function.args[0], 6501);
                printf("^(");
                print_standard(exp, 0);
                printf(")");
                if (pos_exp) expr_free(pos_exp);
            } else if (is_negative) {
                /* The denominator is the RIGHT operand of a left-associative
                 * `/`, so it is printed one level above Divide (5001): a
                 * base that itself prints as a quotient (Power[x, -1],
                 * Rational, Divide -- all 5000) must be parenthesised, or
                 * Power[Power[x, -1], -1] prints "1/1/x", which re-reads as
                 * (1/1)/x = 1/x instead of x. */
                printf("1/");
                if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                    print_standard(e->data.function.args[0], 5001);
                } else {
                    Expr* args[2] = { expr_copy(e->data.function.args[0]), expr_copy(pos_exp) };
                    Expr* tmp = expr_new_function(expr_new_symbol(SYM_Power), args, 2);
                    print_standard(tmp, 5001);
                    expr_free(tmp);
                }
                if (pos_exp) expr_free(pos_exp);
            } else if (is_rational(exp, &n, &d) && n == 1 && d == 2) {
                printf("Sqrt[");
                print_standard(e->data.function.args[0], 0);
                printf("]");
            } else {
                /* Negative numeric atoms parse as unary-minus, which has lower
                 * precedence than Power. So `(-1)^(1/3)` round-trips as
                 * `-1^(1/3)` -> `-(1^(1/3))`. Wrap explicitly. */
                Expr* base = e->data.function.args[0];
                bool base_paren = false;
                int64_t bn, bd;
                if (base->type == EXPR_INTEGER && base->data.integer < 0) base_paren = true;
                else if (base->type == EXPR_REAL && base->data.real < 0.0) base_paren = true;
                else if (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) < 0) base_paren = true;
                else if (is_rational(base, &bn, &bd) && bn < 0) base_paren = true;
                if (base_paren) {
                    printf("(");
                    print_standard(base, 0);
                    printf(")");
                } else {
                    /* Context 6501 (one above Power's 6500) so a Power base is
                     * parenthesised: `^` is right-associative, so (a^2)^x must
                     * print as (a^2)^x, not a^2^x (= a^(2^x)). */
                    print_standard(base, 6501);
                }
                printf("^");
                print_standard(e->data.function.args[1], 6500);
            }
        }
        else if (head == SYM_Times && g_inputform_depth > 0) {
            print_times_inputform(e);
        }
        else if (head == SYM_Times) {
            size_t count = e->data.function.arg_count;

            /* Pull a leading negative numeric factor out as unary minus so
             * Times[-1, Infinity] prints "-Infinity", Times[-1, a, b] prints
             * "-a b", and Times[-3, x] prints "-3 x" (we keep the magnitude
             * as a factor when it isn't 1). Without this the numerator loop
             * below would emit "-1 Infinity". */
            int64_t lead_start = 0;
            bool flipped_sign = false;
            Expr* flipped_head = NULL; /* only used when lead arg != -1 */
            if (count > 0) {
                Expr* a0 = e->data.function.args[0];
                if (a0->type == EXPR_INTEGER && a0->data.integer < 0) {
                    flipped_sign = true;
                    if (a0->data.integer == -1 && count > 1) {
                        lead_start = 1; /* drop the -1, emit only "-" */
                    } else {
                        flipped_head = print_neg_int(a0->data.integer);
                    }
                } else if (a0->type == EXPR_REAL && a0->data.real < 0.0) {
                    flipped_sign = true;
                    if (a0->data.real == -1.0 && count > 1) {
                        lead_start = 1;
                    } else {
                        flipped_head = expr_new_real(-a0->data.real);
                    }
                } else if (a0->type == EXPR_FUNCTION &&
                           a0->data.function.head->type == EXPR_SYMBOL &&
                           a0->data.function.head->data.symbol.name == SYM_Rational &&
                           a0->data.function.arg_count == 2 &&
                           a0->data.function.args[0]->type == EXPR_INTEGER &&
                           a0->data.function.args[0]->data.integer < 0) {
                    flipped_sign = true;
                    Expr* num = a0->data.function.args[0];
                    Expr* den = a0->data.function.args[1];
                    Expr* rargs[2] = { print_neg_int(num->data.integer), expr_copy(den) };
                    flipped_head = expr_new_function(expr_new_symbol(SYM_Rational), rargs, 2);
                } else if (a0->type == EXPR_FUNCTION &&
                           a0->data.function.head->type == EXPR_SYMBOL &&
                           a0->data.function.head->data.symbol.name == SYM_Complex &&
                           a0->data.function.arg_count == 2 &&
                           a0->data.function.args[0]->type == EXPR_INTEGER &&
                           a0->data.function.args[0]->data.integer == 0 &&
                           a0->data.function.args[1]->type == EXPR_INTEGER &&
                           a0->data.function.args[1]->data.integer < 0) {
                    /* Pure-imaginary-negative lead factor: Times[Complex[0, -n], ...]
                     * pulls the sign out so it prints as "-n I ..." instead of
                     * "(-n I) ..." or, for n == 1, "-I ..." instead of "(-I) ...". */
                    flipped_sign = true;
                    Expr* im = a0->data.function.args[1];
                    if (im->data.integer == -1 && count > 1) {
                        /* Drop the imaginary unit's -1; the loop will print
                         * Complex[0, 1] = "I" as the new lead factor. */
                        Expr* cargs[2] = { expr_new_integer(0), expr_new_integer(1) };
                        flipped_head = expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
                    } else {
                        Expr* cargs[2] = { expr_new_integer(0), print_neg_int(im->data.integer) };
                        flipped_head = expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
                    }
                }
            }
            if (flipped_sign) printf("-");

            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t num_count = 0, den_count = 0;

            for (size_t i = lead_start; i < count; i++) {
                Expr* arg = (i == 0 && flipped_head) ? flipped_head : e->data.function.args[i];
                if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol.name == SYM_Power && arg->data.function.arg_count == 2) {
                    Expr* exp = arg->data.function.args[1];
                    bool is_neg = false;
                    int64_t n, d;
                    Expr* pos_exp = NULL;
                    if ((pos_exp = print_neg_int_exp(exp)) != NULL) {
                        is_neg = true;
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        is_neg = true;
                        pos_exp = make_rational_neg(n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        is_neg = true;
                        pos_exp = expr_new_real(-exp->data.real);
                    } else if ((pos_exp = print_neg_symbolic_exp(exp)) != NULL) {
                        is_neg = true;
                    }
                    if (is_neg) {
                        if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                            den_args[den_count++] = expr_copy(arg->data.function.args[0]);
                        } else {
                            Expr* p_args[2] = { expr_copy(arg->data.function.args[0]), pos_exp };
                            den_args[den_count++] = expr_new_function(expr_new_symbol(SYM_Power), p_args, 2);
                            pos_exp = NULL; 
                        }
                        if (pos_exp) expr_free(pos_exp);
                    } else {
                        num_args[num_count++] = expr_copy(arg);
                    }
                } else {
                    num_args[num_count++] = expr_copy(arg);
                }
            }
            
            fold_rational_coeff_into_den(num_args, &num_count, den_args, &den_count);

            if (den_count == 0) {
                for (size_t i = 0; i < num_count; i++) {
                    if (i > 0) printf(" ");
                    print_standard(num_args[i], 4500);
                }
            } else {
                Expr* num = NULL;
                if (num_count == 0) num = expr_new_integer(1);
                else if (num_count == 1) num = expr_copy(num_args[0]);
                else {
                    Expr** nc = malloc(sizeof(Expr*) * num_count);
                    for(size_t i=0; i<num_count; i++) nc[i] = expr_copy(num_args[i]);
                    num = expr_new_function(expr_new_symbol(SYM_Times), nc, num_count);
                    free(nc);
                }

                Expr* den = NULL;
                if (den_count == 1) den = expr_copy(den_args[0]);
                else {
                    Expr** dc = malloc(sizeof(Expr*) * den_count);
                    for(size_t i=0; i<den_count; i++) dc[i] = expr_copy(den_args[i]);
                    den = expr_new_function(expr_new_symbol(SYM_Times), dc, den_count);
                    free(dc);
                }
                
                print_standard(num, 5000);
                printf("/");
                /* 5001: the right operand of a left-associative `/` (see the
                 * Power branch) -- a quotient-valued denominator factor must
                 * be parenthesised, a/(1/x) not a/1/x. */
                print_standard(den, 5001);
                
                expr_free(num);
                expr_free(den);
            }
            
            for (size_t i = 0; i < num_count; i++) expr_free(num_args[i]);
            free(num_args);
            for (size_t i = 0; i < den_count; i++) expr_free(den_args[i]);
            free(den_args);
            if (flipped_head) expr_free(flipped_head);
        }

        else if (head == SYM_Inequality
                 && e->data.function.arg_count >= 3
                 && (e->data.function.arg_count & 1u) == 1) {
            /* Inequality[v0, op0, v1, op1, v2, ...]. Walk the alternating
             * (value, op-symbol) sequence, printing each binary operator
             * using the same glyphs as the corresponding binary heads. */
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* a = e->data.function.args[i];
                if ((i & 1u) == 0) {
                    print_standard(a, my_prec);
                } else {
                    const char* op = " ?? ";
                    if (a->type == EXPR_SYMBOL) {
                        const char* s = a->data.symbol.name;
                        if      (s == SYM_Equal)        op = " == ";
                        else if (s == SYM_Less)         op = " < ";
                        else if (s == SYM_Greater)      op = " > ";
                        else if (s == SYM_LessEqual)    op = " <= ";
                        else if (s == SYM_GreaterEqual) op = " >= ";
                        else if (s == SYM_Unequal)      op = " != ";
                    }
                    printf("%s", op);
                }
            }
        }
        else if (head == SYM_MessageName && e->data.function.arg_count == 2 &&
                 e->data.function.args[1]->type == EXPR_STRING) {
            /* MessageName[f, "tag"] -> f::tag  (the tag prints bare, not quoted). */
            print_standard(e->data.function.args[0], my_prec);
            printf("::%s", e->data.function.args[1]->data.string);
        }
#if UP_VALUES
        else if ((head == SYM_UpSet || head == SYM_UpSetDelayed) &&
                 e->data.function.arg_count >= 2) {
            /* lhs ^= rhs  /  lhs ^:= rhs */
            const char* op = (head == SYM_UpSetDelayed) ? " ^:= " : " ^= ";
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf("%s", op);
                print_standard(e->data.function.args[i], my_prec);
            }
        }
        else if ((head == SYM_TagSet || head == SYM_TagSetDelayed) &&
                 e->data.function.arg_count == 3) {
            /* TagSet[f,lhs,rhs] -> f /: lhs = rhs ; TagSetDelayed -> f /: lhs := rhs */
            print_standard(e->data.function.args[0], my_prec);
            printf(" /: ");
            print_standard(e->data.function.args[1], my_prec);
            printf("%s", (head == SYM_TagSetDelayed) ? " := " : " = ");
            print_standard(e->data.function.args[2], my_prec);
        }
        else if (head == SYM_TagUnset && e->data.function.arg_count == 2) {
            /* TagUnset[f, lhs] -> f /: lhs =. */
            print_standard(e->data.function.args[0], my_prec);
            printf(" /: ");
            print_standard(e->data.function.args[1], my_prec);
            printf(" =.");
        }
        else if (head == SYM_Definition && e->data.function.arg_count == 1 &&
                 e->data.function.args[0]->type == EXPR_SYMBOL) {
            /* Definition[sym] is an inert object that DISPLAYS its definitions
             * (FullForm still shows Definition[sym]). */
            print_definition_body(e->data.function.args[0]->data.symbol.name);
        }
#endif
        else if (head == SYM_Not && e->data.function.arg_count == 1) {
            /* Not[a] -> !a. The operand binds at least as tightly as a
             * comparison; a nested Not is parenthesised, !(!a), since `!!`
             * would read as one token. */
            printf("!");
            print_standard(e->data.function.args[0], my_prec + 1);
        }
        else if (head == SYM_Part && e->data.function.arg_count >= 2) {
            /* Part[a, i, j] -> a[[i, j]]: the subject binds as tightly as an
             * application, (a + b)[[1]]. */
            print_standard(e->data.function.args[0], 9500);
            printf("[[");
            for (size_t i = 1; i < e->data.function.arg_count; i++) {
                if (i > 1) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("]]");
        }
        else if (((head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_And || head == SYM_Or || head == SYM_Alternatives || head == SYM_NonCommutativeMultiply || head == SYM_Dot || head == SYM_StringExpression) && e->data.function.arg_count >= 2) ||
                 ((head == SYM_Set || head == SYM_SetDelayed || head == SYM_Rule || head == SYM_RuleDelayed || head == SYM_Condition) && e->data.function.arg_count == 2)) {
            const char* op = "";
            if (head == SYM_Equal) op = " == ";
            else if (head == SYM_Unequal) op = " != ";
            else if (head == SYM_Less) op = " < ";
            else if (head == SYM_Greater) op = " > ";
            else if (head == SYM_LessEqual) op = " <= ";
            else if (head == SYM_GreaterEqual) op = " >= ";
            else if (head == SYM_SameQ) op = " === ";
            else if (head == SYM_UnsameQ) op = " =!= ";
            else if (head == SYM_Set) op = " = ";
            else if (head == SYM_SetDelayed) op = " := ";
            else if (head == SYM_Rule) op = " -> ";
            else if (head == SYM_RuleDelayed) op = " :> ";
            else if (head == SYM_Condition) op = " /; ";
            else if (head == SYM_And) op = " && ";
            else if (head == SYM_Or) op = " || ";
            else if (head == SYM_Alternatives) op = " | ";
            else if (head == SYM_NonCommutativeMultiply) op = " ** ";
            else if (head == SYM_Dot) op = " . ";
            else if (head == SYM_StringExpression) op = "~~";   /* as Mathematica */

            /* Operands at the operator's own precedence are parenthesised
             * except on the side the operator associates to: the parser
             * builds a == b == c as one Equal, a -> b -> c as a -> (b -> c)
             * and a /; b /; c as (a /; b) /; c, so (a -> b) -> c and
             * Equal[Equal[a, b], c] need their parentheses to read back as
             * the same expression (audit Rule#1). */
            size_t n_ops = e->data.function.arg_count;
            bool right_assoc = (head == SYM_Set || head == SYM_SetDelayed ||
                                head == SYM_Rule || head == SYM_RuleDelayed);
            bool left_assoc = (head == SYM_Condition);
            /* a nested same-head operand of a Flat head (a held
             * And[And[a, b], c]) reads back as the same flat value */
            bool flat = (get_attributes(head) & ATTR_FLAT) != 0;
            for (size_t i = 0; i < n_ops; i++) {
                if (i > 0) printf("%s", op);
                int p = my_prec + 1;
                if (flat || (right_assoc && i == n_ops - 1) || (left_assoc && i == 0)) p = my_prec;
                print_standard(e->data.function.args[i], p);
            }
        }
        else if (head == SYM_Plus && g_inputform_depth == 0 &&
                 print_plus_spliced_series(e)) {
            /* printed: a SeriesData summand's terms were spliced in. */
        }
        else if (head == SYM_Plus) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* arg = e->data.function.args[i];
                bool is_negative = false;
                Expr* to_print = arg;
                Expr* t_copy = NULL;

                if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
                    is_negative = true;
                    t_copy = print_neg_int(arg->data.integer);
                    to_print = t_copy;
                } else if (arg->type == EXPR_BIGINT && mpz_sgn(arg->data.bigint) < 0) {
                    is_negative = true;
                    t_copy = negate_bigint(arg);
                    to_print = t_copy;
                } else if (arg->type == EXPR_REAL && arg->data.real < 0) {
                    is_negative = true;
                    t_copy = expr_new_real(-arg->data.real);
                    to_print = t_copy;
#ifdef USE_MPFR
                } else if (arg->type == EXPR_MPFR && mpfr_sgn(arg->data.mpfr) < 0) {
                    /* an arbitrary-precision term signs the operator, as a
                     * machine one does: a - 0.25`24, not a + -0.25`24 */
                    is_negative = true;
                    t_copy = expr_new_mpfr_copy(arg->data.mpfr);
                    mpfr_neg(t_copy->data.mpfr, t_copy->data.mpfr, MPFR_RNDN);
                    to_print = t_copy;
#endif
                } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
                    const char* h = arg->data.function.head->data.symbol.name;
                    if (h == SYM_Times && arg->data.function.arg_count > 0) {
                        Expr* f_arg = arg->data.function.args[0];
                        if (f_arg->type == EXPR_INTEGER && f_arg->data.integer < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = print_neg_int(f_arg->data.integer);

                            if (t_copy->data.function.args[0]->data.integer == 1 && t_copy->data.function.arg_count > 1) {
                                Expr** new_args = malloc(sizeof(Expr*) * (t_copy->data.function.arg_count - 1));
                                for (size_t k = 1; k < t_copy->data.function.arg_count; k++) {
                                    new_args[k-1] = expr_copy(t_copy->data.function.args[k]);
                                }
                                /* one factor left (-1*1 -> 1): print it bare, not as
                                 * the one-argument Times[1] InputForm prints empty */
                                Expr* new_times = (t_copy->data.function.arg_count == 2)
                                    ? new_args[0]
                                    : expr_new_function(expr_new_symbol(SYM_Times), new_args, t_copy->data.function.arg_count - 1);
                                expr_free(t_copy);
                                t_copy = new_times;
                                free(new_args);
                            }
                            to_print = t_copy;
                        } else if (f_arg->type == EXPR_BIGINT && mpz_sgn(f_arg->data.bigint) < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = negate_bigint(f_arg);

                            if (t_copy->data.function.args[0]->type == EXPR_INTEGER &&
                                t_copy->data.function.args[0]->data.integer == 1 &&
                                t_copy->data.function.arg_count > 1) {
                                Expr** new_args = malloc(sizeof(Expr*) * (t_copy->data.function.arg_count - 1));
                                for (size_t k = 1; k < t_copy->data.function.arg_count; k++) {
                                    new_args[k-1] = expr_copy(t_copy->data.function.args[k]);
                                }
                                /* one factor left (-1*1 -> 1): print it bare, not as
                                 * the one-argument Times[1] InputForm prints empty */
                                Expr* new_times = (t_copy->data.function.arg_count == 2)
                                    ? new_args[0]
                                    : expr_new_function(expr_new_symbol(SYM_Times), new_args, t_copy->data.function.arg_count - 1);
                                expr_free(t_copy);
                                t_copy = new_times;
                                free(new_args);
                            }
                            to_print = t_copy;
                        } else if (f_arg->type == EXPR_REAL && f_arg->data.real < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = expr_new_real(-f_arg->data.real);
                            to_print = t_copy;
#ifdef USE_MPFR
                        } else if (f_arg->type == EXPR_MPFR && mpfr_sgn(f_arg->data.mpfr) < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            Expr* pos = expr_new_mpfr_copy(f_arg->data.mpfr);
                            mpfr_neg(pos->data.mpfr, pos->data.mpfr, MPFR_RNDN);
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = pos;
                            to_print = t_copy;
#endif
                        } else if (f_arg->type == EXPR_FUNCTION && f_arg->data.function.head->type == EXPR_SYMBOL && f_arg->data.function.head->data.symbol.name == SYM_Rational) {
                            Expr* num = f_arg->data.function.args[0];
                            if (num->type == EXPR_INTEGER && num->data.integer < 0) {
                                is_negative = true;
                                t_copy = expr_unshare(expr_copy(arg));
                                /* Need to mutate the rational nested inside
                                 * t_copy.args[0]: unshare that level too. */
                                Expr* rat = expr_unshare(t_copy->data.function.args[0]);
                                t_copy->data.function.args[0] = rat;
                                expr_free(rat->data.function.args[0]);
                                rat->data.function.args[0] = print_neg_int(num->data.integer);
                                to_print = t_copy;
                            }
                        } else if (f_arg->type == EXPR_FUNCTION &&
                                   f_arg->data.function.head->type == EXPR_SYMBOL &&
                                   f_arg->data.function.head->data.symbol.name == SYM_Complex &&
                                   f_arg->data.function.arg_count == 2 &&
                                   f_arg->data.function.args[0]->type == EXPR_INTEGER &&
                                   f_arg->data.function.args[0]->data.integer == 0 &&
                                   f_arg->data.function.args[1]->type == EXPR_INTEGER &&
                                   f_arg->data.function.args[1]->data.integer < 0) {
                            /* Plus[..., Times[Complex[0, -n], rest...]]:
                             * print "- Times[Complex[0, n], rest...]" so the
                             * sign attaches to the operator instead of the
                             * factor. */
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            Expr* cplx = expr_unshare(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = cplx;
                            int64_t im_val = cplx->data.function.args[1]->data.integer;
                            expr_free(cplx->data.function.args[1]);
                            cplx->data.function.args[1] = expr_new_integer(-im_val);
                            to_print = t_copy;
                        }
                    } else if (h == SYM_Rational && arg->data.function.arg_count == 2) {
                        Expr* num = arg->data.function.args[0];
                        if (num->type == EXPR_INTEGER && num->data.integer < 0) {
                            is_negative = true;
                            t_copy = expr_unshare(expr_copy(arg));
                            expr_free(t_copy->data.function.args[0]);
                            t_copy->data.function.args[0] = print_neg_int(num->data.integer);
                            to_print = t_copy;
                        }
                    } else if (h == SYM_Complex && arg->data.function.arg_count == 2 &&
                               arg->data.function.args[0]->type == EXPR_INTEGER &&
                               arg->data.function.args[0]->data.integer == 0 &&
                               arg->data.function.args[1]->type == EXPR_INTEGER &&
                               arg->data.function.args[1]->data.integer < 0) {
                        /* Plus[..., Complex[0, -n]] -> "- Complex[0, n]" so the
                         * sign sits with the operator. */
                        is_negative = true;
                        int64_t im_val = arg->data.function.args[1]->data.integer;
                        Expr* cargs[2] = { expr_new_integer(0), expr_new_integer(-im_val) };
                        t_copy = expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
                        to_print = t_copy;
                    }
                }

                if (i > 0) {
                    if (is_negative) {
                        printf(" - ");
                    } else {
                        printf(" + ");
                    }
                } else if (g_inputform_depth > 0 && arg->type == EXPR_FUNCTION
                           && arg->data.function.head->type == EXPR_SYMBOL
                           && arg->data.function.head->data.symbol.name == SYM_Times) {
                    /* InputForm: a leading product keeps its own sign
                     * (-1/2*a - b, -(a*b) + c*d, as the reference) */
                    to_print = arg;
                } else {
                    if (is_negative) {
                        printf("-");
                    }
                }

                /* a held sum inside a sum keeps its parentheses:
                 * HoldForm[a - (b - c)], HoldForm[a + (b + c) + d] */
                bool nested_plus = to_print->type == EXPR_FUNCTION &&
                    to_print->data.function.head->type == EXPR_SYMBOL &&
                    to_print->data.function.head->data.symbol.name == SYM_Plus;
                print_standard(to_print, nested_plus ? 3501 : 3500);
                if (t_copy) expr_free(t_copy);
            }
        }
        else if (head == SYM_List) {
            printf("{");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("}");
        }
        else if (head == SYM_Association && assoc_prints_special(e)) {
            /* <|k1 -> v1, k2 -> v2|>. Each arg is a Rule[k,v] that
             * print_standard already renders as "k -> v". */
            printf("<|");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("|>");
        }
        else if (head == SYM_Blank && e->data.function.arg_count <= 1) {
            printf("_");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 9500);
        }
        else if (head == SYM_BlankSequence && e->data.function.arg_count <= 1) {
            printf("__");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 9500);
        }
        else if (head == SYM_BlankNullSequence && e->data.function.arg_count <= 1) {
            printf("___");
            if (e->data.function.arg_count == 1) print_standard(e->data.function.args[0], 9500);
        }
        else if (head == SYM_Slot && e->data.function.arg_count == 1) {
            /* Slot[n] -> #n (e.g. #1). Atom-level precedence already returned
             * by get_expr_prec via the fall-through, so no extra parens. */
            Expr* idx = e->data.function.args[0];
            printf("#");
            if (idx->type == EXPR_INTEGER) {
                printf("%" PRId64, idx->data.integer);
            } else if (idx->type == EXPR_STRING) {
                /* Named slot: Slot["a"] -> #a when "a" lexes back as a
                 * name ([A-Za-z$][A-Za-z0-9$`]*), else #"a b" (quoted). */
                const char* nm = idx->data.string;
                bool bare = nm[0] && (isalpha((unsigned char)nm[0]) || nm[0] == '$');
                for (const char* q = nm; bare && *q; q++)
                    if (!(isalnum((unsigned char)*q) || *q == '$' ||
                          (*q == '`' && (isalpha((unsigned char)q[1]) || q[1] == '$'))))
                        bare = false;
                if (bare) printf("%s", nm);
                else print_string_literal(nm);
            } else {
                print_standard(idx, 9500);
            }
        }
        else if (head == SYM_SlotSequence && e->data.function.arg_count == 1) {
            /* SlotSequence[n] -> ##n; SlotSequence[1] is conventionally "##". */
            Expr* idx = e->data.function.args[0];
            printf("##");
            if (idx->type == EXPR_INTEGER) {
                if (idx->data.integer != 1) printf("%" PRId64, idx->data.integer);
            } else {
                print_standard(idx, 9500);
            }
        }
        else if (head == SYM_Function && e->data.function.arg_count == 1) {
            /* Function[body] -> body &. The body normally has high precedence
             * (Plus/Times/Power etc.) so it prints without extra parens; if
             * something genuinely lower-precedence shows up, print_standard
             * inserts them at parent_prec=1000. */
            print_standard(e->data.function.args[0], 1000);
            printf(" &");
        }
        else if (head == SYM_InterpolatingFunction && e->data.function.arg_count >= 2 &&
                 g_inputform_depth == 0) {
            /* Only the domain element is printed explicitly; the remaining
             * data elements are abbreviated as <>, the reference's standard
             * output format. FullForm (handled separately) still reveals the
             * full structure, and InputForm prints the object in full so
             * ToExpression[ToString[InputForm[f]]] reads back the same function
             * (audit InterpolatingFunction#1). */
            printf("InterpolatingFunction[");
            print_standard(e->data.function.args[0], 0);
            printf(", <>]");
        }
        else if (head == SYM_InterpolatingFunction && e->data.function.arg_count >= 2 &&
                 e->data.function.args[1]->type == EXPR_NDARRAY &&
                 e->data.function.args[1]->data.ndarray.rank == 2 &&
                 e->data.function.args[1]->data.ndarray.dims[1] >= 2) {
            /* InputForm of a machine-data object: its private packed table
             * (rows {c1, ..., cm, v}) is written as the boxed table the
             * constructor reads, {{c, v}, ...} or {{{c1, ..., cm}, v}, ...},
             * since a held NDArray[...] call would not re-read as a buffer. */
            Expr* rows = ndarray_to_nested_list(e->data.function.args[1]);
            int64_t cols = e->data.function.args[1]->data.ndarray.dims[1];
            if (rows && cols > 2) {
                for (size_t i = 0; i < rows->data.function.arg_count; i++) {
                    Expr* row = rows->data.function.args[i];
                    size_t m = (size_t)cols - 1;
                    Expr** cs = malloc(sizeof(Expr*) * m);
                    for (size_t k = 0; k < m; k++) cs[k] = expr_copy(row->data.function.args[k]);
                    Expr* coords = expr_new_function(expr_new_symbol(SYM_List), cs, m);
                    free(cs);
                    Expr* pair = expr_new_function(expr_new_symbol(SYM_List),
                        (Expr*[]){ coords, expr_copy(row->data.function.args[m]) }, 2);
                    expr_free(row);
                    rows->data.function.args[i] = pair;
                }
            }
            printf("InterpolatingFunction[");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i) printf(", ");
                print_standard(i == 1 && rows ? rows : e->data.function.args[i], 0);
            }
            printf("]");
            if (rows) expr_free(rows);
        }
        else if ((head == SYM_PredictorFunction ||
                  head == SYM_DimensionReducerFunction ||
                  head == SYM_LearnedDistribution ||
                  head == SYM_ClassifierFunction) &&
                 e->data.function.arg_count >= 2 &&
                 e->data.function.args[0]->type == EXPR_STRING &&
                 g_inputform_depth == 0) {
            /* (InputForm prints the model in full, so ToString[m, InputForm]
             * re-reads as the same model.) */
            /* A fitted machine-learning model prints its METHOD and elides its
             * parameters, exactly as InterpolatingFunction above elides its data and
             * for the same reason. This is not cosmetic: a NearestNeighbors predictor
             * carries its entire training set as its parameter block, so the unabridged
             * form is unreadable at any real size and scrolls a useful answer off the
             * screen. The method name is the part a reader actually wants back.
             *
             * FullForm still reveals everything, which is what makes eliding safe --
             * the information is one keystroke away rather than gone. */
            print_standard(e->data.function.head, 0);
            printf("[");
            print_standard(e->data.function.args[0], 0);
            printf(", <>]");
        }
        else if (head == SYM_PatternTest && e->data.function.arg_count == 2) {
            /* PatternTest[p, test] -> p?test. `?` groups to the left, so the
             * pattern side prints at the operator's own level and the test at
             * one above it: a compound test such as a pure function prints
             * parenthesised, x_?(#1 > 0 & ), and the text re-reads to the
             * same tree. */
            print_standard(e->data.function.args[0], my_prec);
            printf("?");
            print_standard(e->data.function.args[1], my_prec + 1);
        }
        else if (head == SYM_Pattern && e->data.function.arg_count == 2 &&
                 e->data.function.args[0]->type == EXPR_SYMBOL &&
                 e->data.function.args[1]->type == EXPR_FUNCTION &&
                 e->data.function.args[1]->data.function.head->type == EXPR_SYMBOL &&
                 (e->data.function.args[1]->data.function.head->data.symbol.name == SYM_Blank ||
                  e->data.function.args[1]->data.function.head->data.symbol.name == SYM_BlankSequence ||
                  e->data.function.args[1]->data.function.head->data.symbol.name == SYM_BlankNullSequence) &&
                 e->data.function.args[1]->data.function.arg_count <= 1) {
            /* Pattern[name, Blank[...]] -> name_ / name_h / name__ / name___ */
            printf("%s", e->data.function.args[0]->data.symbol.name);
            print_standard(e->data.function.args[1], 9500);
        }
        else if (head == SYM_Pattern && e->data.function.arg_count == 2 &&
                 e->data.function.args[0]->type == EXPR_SYMBOL) {
            /* Pattern[x, p] -> x:p (the body binds tighter than `:`). */
            printf("%s:", e->data.function.args[0]->data.symbol.name);
            print_standard(e->data.function.args[1], my_prec + 1);
        }
        else if ((head == SYM_Repeated || head == SYM_RepeatedNull) &&
                 e->data.function.arg_count == 1) {
            /* Repeated[p] -> p.. ; RepeatedNull[p] -> p...
             * A bare blank is parenthesised, as WL prints it: `_..` would
             * re-read as Optional `_.` and a stray `.`, a syntax error, and
             * `x_...` as Repeated[x_.] -- so (_).. and (x_)... instead.  A
             * headed blank (`_h..`) needs no parentheses. */
            Expr* rp = e->data.function.args[0];
            Expr* rb = rp;
            if (rb->type == EXPR_FUNCTION && rb->data.function.head->type == EXPR_SYMBOL
                && rb->data.function.head->data.symbol.name == SYM_Pattern
                && rb->data.function.arg_count == 2)
                rb = rb->data.function.args[1];
            bool bare_blank = is_plain_blank(rb) && rb->data.function.arg_count == 0;
            if (bare_blank) printf("(");
            print_standard(rp, my_prec + 1);
            if (bare_blank) printf(")");
            printf(head == SYM_Repeated ? ".." : "...");
        }
        else {
            print_standard(e->data.function.head, 9500);
            printf("[");
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf(", ");
                print_standard(e->data.function.args[i], 0);
            }
            printf("]");
        }
    } else {
        print_standard(e->data.function.head, 9500);
        printf("[");
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (i > 0) printf(", ");
            print_standard(e->data.function.args[i], 0);
        }    
        printf("]");
    }

    if (need_parens) printf(")");
}

void expr_print(Expr* e) {
    bool saved_af = g_apply_formats;
    g_apply_formats = true;
    g_format_depth = 0;
    print_standard(e, 0);
    g_apply_formats = saved_af;
    printf("\n");
}

/* ---------------------------------------------------------------------------
 * Portable capture of the print sink into a heap string.
 *
 * expr_to_string & friends print into an in-memory FILE* and hand back the
 * bytes.  POSIX gives us open_memstream for exactly this; MinGW/UCRT (the
 * native Windows target) ships neither open_memstream nor fmemopen, so there we
 * route the sink through a temp FILE and read it back on finish.  String
 * rendering is a per-result operation, not an inner numeric loop, so the Windows
 * disk round-trip is acceptable -- the packed/NDArray hot paths never touch it.
 *
 * str_capture_open() returns a capture whose FILE* (->f) is pushed as the sink
 * (NULL on failure); str_capture_finish() closes it and returns a freshly
 * malloc'd, NUL-terminated string the caller frees (NULL on failure).
 *
 * The capture records live in a STATIC pool, not in the caller's frame.
 * open_memstream keeps the addresses of buf/len and writes through them on
 * every flush, so a capture abandoned by a non-local exit -- TimeConstrained's
 * timeout siglongjmp landing inside print_standard, or a FLINT exception --
 * would otherwise keep a FILE that writes into a DEAD stack frame.  Worse, the
 * skipped mth_out_pop left that FILE installed as the sink, so every later
 * print went into it: each write grew the buffer and stored the new pointer and
 * length into whatever frame now occupied that stack slot.  That was the
 * random stack corruption behind integrate_risch_transcendental_tests'
 * intermittent failures (2026-11-20): an early test's millisecond
 * TimeConstrained budgets unwound mid-print, and a later, unrelated case
 * recursed without end (risch_integrate_hypertangent_reduced overflowed the
 * stack), returned $Aborted, or the process aborted in exit() flushing the
 * orphan stream.  With the records here, mth_print_unwind can find and close
 * every abandoned capture and put the sink back (see mth_print_mark).
 *
 * The pool only grows; a slot's address never changes once handed out (the
 * FILE points into it), so it is an array of individually allocated slots.
 * ------------------------------------------------------------------------- */
typedef struct {
    FILE*  f;
#ifndef _WIN32
    char*  buf;
    size_t len;
#endif
} StrCapture;

static StrCapture** g_cap = NULL;   /* g_cap[0 .. g_ncap) are open, innermost last */
static int g_ncap = 0, g_cap_alloc = 0;

static StrCapture* str_capture_open(void) {
    if (g_ncap == g_cap_alloc) {
        int na = g_cap_alloc ? 2 * g_cap_alloc : 16;
        StrCapture** a = realloc(g_cap, (size_t)na * sizeof *a);
        if (!a) return NULL;
        for (int i = g_cap_alloc; i < na; i++) a[i] = NULL;
        g_cap = a;
        g_cap_alloc = na;
    }
    if (!g_cap[g_ncap] && !(g_cap[g_ncap] = malloc(sizeof(StrCapture)))) return NULL;
    StrCapture* c = g_cap[g_ncap];
    /* f is NULL until the stream exists and the slot is counted before it is
     * opened, so an unwind at any point sees either no stream or a live one */
    c->f = NULL;
    g_ncap++;
#ifdef _WIN32
    c->f = tmpfile();
#else
    c->buf = NULL;
    c->len = 0;
    c->f = open_memstream(&c->buf, &c->len);
#endif
    if (!c->f) { g_ncap--; return NULL; }
    return c;
}

static char* str_capture_finish(StrCapture* c) {
    g_ncap--;              /* released first: an unwind never closes it twice */
#ifdef _WIN32
    if (fflush(c->f) != 0) { fclose(c->f); return NULL; }
    long n = ftell(c->f);
    if (n < 0) { fclose(c->f); return NULL; }
    rewind(c->f);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(c->f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, c->f);
    buf[got] = '\0';
    fclose(c->f);
    return buf;
#else
    fclose(c->f);          /* flushes; finalises buf/len, NUL-terminated */
    return c->buf;
#endif
}

/* See print.h.  The mark is taken where a non-local exit can land
 * (tc_run_guarded_mode in core.c) and every capture opened after it, with the
 * sink and the formatting flags set inside it, is undone by the unwind. */
void mth_print_mark(MthPrintMark* m) {
    m->sink = g_print_sink;
    m->ncap = g_ncap;
    m->inputform_depth = g_inputform_depth;
    m->format_depth = g_format_depth;
    m->numberform_ctx = g_numberform_ctx;
    m->output_form = g_print_output_form;
    m->apply_formats = g_apply_formats;
    m->fullform_depth = g_fullform_depth;
}

void mth_print_unwind(const MthPrintMark* m) {
    /* the sink first: nothing may print into a stream about to be closed */
    g_print_sink = m->sink;
    while (g_ncap > m->ncap) {
        StrCapture* c = g_cap[--g_ncap];
        if (!c->f) continue;
        fclose(c->f);
#ifndef _WIN32
        free(c->buf);      /* the abandoned render's partial text */
#endif
        c->f = NULL;
    }
    g_inputform_depth = m->inputform_depth;
    g_format_depth = m->format_depth;
    g_numberform_ctx = m->numberform_ctx;
    g_print_output_form = m->output_form;
    g_apply_formats = m->apply_formats;
    g_fullform_depth = m->fullform_depth;
}

/* expr_to_string with display rules (Format) applied: what the front end shows
 * for a result. Caller frees. */
char* expr_to_string_display(Expr* e) {
    bool saved_af = g_apply_formats;
    int saved_fd = g_format_depth;
    g_apply_formats = true;
    g_format_depth = 0;
    char* s = expr_to_string(e);
    g_format_depth = saved_fd;
    g_apply_formats = saved_af;
    return s;
}

char* expr_to_string(Expr* e) {
    // Print to a memory buffer by diverting the print sink to it.
    StrCapture* c = str_capture_open();
    if (!c) return NULL;

    FILE* prev = mth_out_push(c->f);
    print_standard(e, 0);
    mth_out_pop(prev);

    return str_capture_finish(c);
}

/* OutputForm rendering into a heap string (caller frees): as expr_to_string,
 * except that strings render as their characters, recursively, exactly as
 * Print shows them -- ToString[{"a", 1}] is "{a, 1}". Linear syntax is kept
 * (Mathilda has no 2-D OutputForm), so ToString[1/2] is "1/2". */
char* expr_to_string_outputform(Expr* e) {
    StrCapture* c = str_capture_open();
    if (!c) return NULL;

    FILE* prev = mth_out_push(c->f);
    bool saved_of = g_print_output_form;
    bool saved_af = g_apply_formats;
    int saved_if = g_inputform_depth;
    int saved_fd = g_format_depth;
    g_print_output_form = true;
    g_apply_formats = true;
    g_inputform_depth = 0;
    g_format_depth = 0;
    print_standard(e, 0);
    g_format_depth = saved_fd;
    g_inputform_depth = saved_if;
    g_apply_formats = saved_af;
    g_print_output_form = saved_of;
    mth_out_pop(prev);

    return str_capture_finish(c);
}

/* InputForm text for ToString[expr, InputForm]: Format[..., InputForm]
 * display rules apply, so it is NOT for re-parsing by
 * the kernel itself -- internal callers use expr_to_string. */
char* expr_to_string_inputform_display(Expr* e) {
    Expr* inner = expr_copy(e);
    Expr* wrapper = expr_new_function(expr_new_symbol(SYM_InputForm), &inner, 1);
    bool saved_af = g_apply_formats;
    int saved_fd = g_format_depth;
    g_apply_formats = true;
    g_format_depth = 0;
    char* out = expr_to_string(wrapper);
    g_format_depth = saved_fd;
    g_apply_formats = saved_af;
    expr_free(wrapper);
    return out;
}

char* expr_to_string_fullform(Expr* e) {
    StrCapture* c = str_capture_open();
    if (!c) return NULL;

    FILE* prev = mth_out_push(c->f);
    expr_print_fullform(e);
    mth_out_pop(prev);

    return str_capture_finish(c);
}

/* As expr_to_string_fullform, with the user-facing FullForm spellings
 * (DirectedInfinity for the infinities): what ToString[e, FullForm] shows.
 * expr_to_string_fullform itself stays the kernel's structural text, which
 * tests and self-checks compare. Caller frees. */
char* expr_to_string_fullform_display(Expr* e) {
    StrCapture* c = str_capture_open();
    if (!c) return NULL;

    FILE* prev = mth_out_push(c->f);
    g_fullform_depth++;
    expr_print_fullform(e);
    g_fullform_depth--;
    mth_out_pop(prev);

    return str_capture_finish(c);
}

/* Render `e` (the result of a NumberForm NumberFormat function) to a heap
 * string in OutputForm: strings print raw, Row[...] concatenates. The active
 * NumberForm context is suspended so already-formatted mantissa/exponent
 * strings are not re-processed. Caller frees. */
char* numberform_format_result_to_string(const Expr* e) {
    StrCapture* c = str_capture_open();
    if (!c) return NULL;

    FILE* prev = mth_out_push(c->f);

    const NumberFormCtx* saved_ctx = g_numberform_ctx;
    bool saved_of = g_print_output_form;
    g_numberform_ctx = NULL;
    g_print_output_form = true;

    print_standard((Expr*)e, 0);

    g_print_output_form = saved_of;
    g_numberform_ctx = saved_ctx;
    mth_out_pop(prev);

    return str_capture_finish(c);
}

#include "symtab.h"


Expr* builtin_print(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* No display rules here: a script's Print writes the raw expression, as
     * wolframscript does (Format[k[x_]] := 11; Print[k[2]] shows k[2]), while
     * ToString[k[2]] is 11. */
    bool saved = g_print_output_form;
    bool saved_af = g_apply_formats;
    g_print_output_form = true;
    g_apply_formats = false;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        print_standard(res->data.function.args[i], 0);
    }
    g_apply_formats = saved_af;
    g_print_output_form = saved;
    printf("\n");
    /* B3: flush per Print so output appears incrementally down a pipe (stdout is
     * fully buffered when not a tty) and survives a timeout kill -- matching
     * wolframscript, which flushes on every Print. A tty is already line
     * buffered, so this is a no-op there. */
    fflush(stdout);
    return expr_new_symbol(SYM_Null);
}

/* Print one Echo line, ">> label value" (or ">> value" when label is NULL),
 * through Print itself so the rendering, flushing and output redirection are
 * Print's.  Both arguments are borrowed.  Shared by Echo and EchoTiming. */
void echo_print_line(Expr* label, Expr* value) {
    Expr* parts[4];
    size_t n = 0;
    parts[n++] = expr_new_string(">> ");
    if (label) {
        parts[n++] = expr_copy(label);
        parts[n++] = expr_new_string(" ");
    }
    parts[n++] = expr_copy(value);
    Expr* call = expr_new_function(expr_new_symbol("Print"), parts, n);
    Expr* r = builtin_print(call);
    expr_free(r);
    expr_free(call);
}

/* Echo[expr] prints ">> expr" and returns expr; Echo[expr, label] prints the
 * label first; Echo[expr, label, f] prints f[expr] in place of expr (the
 * returned value is still expr). */
Expr* builtin_echo(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return builtin_arg_error("Echo", argc, 1, 3);
    Expr* value = res->data.function.args[0];
    Expr* label = argc >= 2 ? res->data.function.args[1] : NULL;
    if (argc == 3) {
        Expr* shown = eval_and_free(expr_new_function(
            expr_copy(res->data.function.args[2]), (Expr*[]){ expr_copy(value) }, 1));
        echo_print_line(label, shown);
        expr_free(shown);
    } else {
        echo_print_line(label, value);
    }
    return expr_copy(value);
}

Expr* builtin_fullform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper
}

Expr* builtin_inputform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper
}

Expr* builtin_texform(Expr* res) {
    (void)res;
    return NULL; // Unevaluated wrapper; rendering happens in print_standard
}

/* TraditionalForm[expr] / StandardForm[expr]: inert display wrappers. They
 * choose the notation the NOTEBOOK typesets their subtree in (print_latex.c);
 * the plain-text printer renders straight through them. */
Expr* builtin_display_form(Expr* res) {
    (void)res;
    return NULL;
}

/*
 * Pretty-prints SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] as
 *     a0 + a1 (x-x0) + ... + a_{k-1} (x-x0)^((nmin+k-1)/den) + O[x-x0]^(nmax/den).
 * The rendering is built by constructing an equivalent Plus-expression of
 * Times/Power nodes and delegating to print_standard so that sign handling,
 * precedence, and fractional-exponent formatting are shared with ordinary
 * expressions. If the arguments are not in a shape we can render (for
 * example non-literal nmin/nmax/den or a non-List coefficient container),
 * we fall back to the generic head[args...] form.
 */

static int64_t abs_i64(int64_t v) { return v < 0 ? -v : v; }

/* If x0 is a syntactically-negative expansion point, sets *abs_out to a
 * freshly-owned positive form and returns true; otherwise returns false.
 * Recognises negative integers/bigints/reals/rationals and Times[] with a
 * negative leading numeric factor, so that compound centres such as
 * -Pi = Times[-1, Pi] and -1/E = Times[-1, Power[E, -1]] are normalised.
 * Without this, series_build_base would wrap such an x0 in Times[-1, x0]
 * and the Plus printer would emit a double minus ("x - -Pi"). */
static bool series_extract_negative(Expr* x0, Expr** abs_out) {
    int64_t rn, rd;
    if (x0->type == EXPR_INTEGER && x0->data.integer < 0) {
        *abs_out = print_neg_int(x0->data.integer); return true;
    }
    if (x0->type == EXPR_BIGINT && mpz_sgn(x0->data.bigint) < 0) {
        *abs_out = negate_bigint(x0); return true;
    }
    if (x0->type == EXPR_REAL && x0->data.real < 0.0) {
        *abs_out = expr_new_real(-x0->data.real); return true;
    }
    if (is_rational(x0, &rn, &rd) && rn < 0) {
        *abs_out = make_rational_neg(rn, rd); return true;
    }
    if (x0->type == EXPR_FUNCTION &&
        x0->data.function.head->type == EXPR_SYMBOL &&
        x0->data.function.head->data.symbol.name == SYM_Times &&
        x0->data.function.arg_count > 0) {
        Expr* f = x0->data.function.args[0];
        Expr* negf = NULL;
        if (f->type == EXPR_INTEGER && f->data.integer < 0)
            negf = print_neg_int(f->data.integer);
        else if (f->type == EXPR_BIGINT && mpz_sgn(f->data.bigint) < 0)
            negf = negate_bigint(f);
        else if (f->type == EXPR_REAL && f->data.real < 0.0)
            negf = expr_new_real(-f->data.real);
        else if (is_rational(f, &rn, &rd) && rn < 0)
            negf = make_rational_neg(rn, rd);
        if (!negf) return false;

        size_t n = x0->data.function.arg_count;
        /* A normalised leading +1 is dropped, mirroring the Plus printer. */
        if (negf->type == EXPR_INTEGER && negf->data.integer == 1 && n > 1) {
            expr_free(negf);
            if (n - 1 == 1) { *abs_out = expr_copy(x0->data.function.args[1]); return true; }
            Expr** na = malloc(sizeof(Expr*) * (n - 1));
            for (size_t k = 1; k < n; k++) na[k-1] = expr_copy(x0->data.function.args[k]);
            *abs_out = expr_new_function(expr_new_symbol(SYM_Times), na, n - 1);
            free(na);
            return true;
        }
        Expr* t = expr_unshare(expr_copy(x0));
        expr_free(t->data.function.args[0]);
        t->data.function.args[0] = negf;
        *abs_out = t;
        return true;
    }
    return false;
}

/* Builds (x - x0), or just x when x0 is a literal zero. The result is a
 * fresh expression owned by the caller. Negative centres are normalised
 * (e.g. x - (-2) becomes x + 2, x - (-Pi) becomes x + Pi) so that the Plus
 * printer produces a clean subtraction. */
static Expr* series_build_base(Expr* x, Expr* x0) {
    if (x0->type == EXPR_INTEGER && x0->data.integer == 0) return expr_copy(x);
    if (x0->type == EXPR_REAL && x0->data.real == 0.0) return expr_copy(x);

    Expr* absx0 = NULL;
    if (series_extract_negative(x0, &absx0)) {
        Expr* args[2] = { expr_copy(x), absx0 };
        return expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
    }

    Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(x0) };
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    Expr* plus_args[2] = { expr_copy(x), neg };
    return expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
}

/* Builds base^(num/den) as an Expr, reducing the exponent to lowest terms.
 * Returns NULL to signal "exponent is zero" (in which case the caller
 * should use the coefficient alone) and a fresh base-copy when the
 * reduced exponent is 1. */
static Expr* series_build_power(Expr* base, int64_t num, int64_t den) {
    if (num == 0) return NULL;
    int64_t g = gcd(abs_i64(num), abs_i64(den));
    if (g == 0) g = 1;
    int64_t n = num / g;
    int64_t d = den / g;
    if (d < 0) { d = -d; n = -n; }

    if (d == 1 && n == 1) return expr_copy(base);

    /* An expansion at +/-Infinity is a SeriesData in Power[x, -1], so the
     * term is (x^k)^n. For an INTEGER n that is exactly x^(k n) for every x,
     * so fold it: x + 1/x + O[1/x]^3, not (1/x)^-1 + 1/x + ..., matching
     * Mathematica. A fractional n is left as (x^k)^(n/d): (1/x)^(1/2) is not
     * x^(-1/2) on the negative reals, which a -Infinity expansion lives on. */
    if (d == 1 && base->type == EXPR_FUNCTION &&
        base->data.function.head->type == EXPR_SYMBOL &&
        base->data.function.head->data.symbol.name == SYM_Power &&
        base->data.function.arg_count == 2 &&
        base->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t k = base->data.function.args[1]->data.integer;
        if (k != 0 && abs_i64(k) <= INT32_MAX && abs_i64(n) <= INT32_MAX) {
            int64_t kn = k * n;
            if (kn == 1) return expr_copy(base->data.function.args[0]);
            Expr* pargs[2] = { expr_copy(base->data.function.args[0]),
                               expr_new_integer(kn) };
            return expr_new_function(expr_new_symbol(SYM_Power), pargs, 2);
        }
    }

    Expr* exp;
    if (d == 1) exp = expr_new_integer(n);
    else exp = make_rational(n, d);

    Expr* pargs[2] = { expr_copy(base), exp };
    return expr_new_function(expr_new_symbol(SYM_Power), pargs, 2);
}

/* Combines a coefficient with an optional power term into a single summand.
 * Returns NULL when the coefficient is a literal zero (i.e. the term is
 * omitted from the displayed sum). A NULL power corresponds to (x-x0)^0
 * and yields just the coefficient. */
static Expr* series_build_term(Expr* coef, Expr* power /* borrowed, may be NULL */) {
    if (coef->type == EXPR_INTEGER && coef->data.integer == 0) return NULL;
    if (coef->type == EXPR_REAL && coef->data.real == 0.0) return NULL;

    if (power == NULL) return expr_copy(coef);

    bool coef_is_one = (coef->type == EXPR_INTEGER && coef->data.integer == 1);
    if (coef_is_one) {
        /* A bare Plus would merge visually with the surrounding series
         * Plus (equal precedence -> no parens), so wrap it in a 1-arg
         * Times to force the Times printer to add parentheses. */
        if (power->type == EXPR_FUNCTION &&
            power->data.function.head->type == EXPR_SYMBOL &&
            power->data.function.head->data.symbol.name == SYM_Plus) {
            Expr* targs[1] = { expr_copy(power) };
            return expr_new_function(expr_new_symbol(SYM_Times), targs, 1);
        }
        return expr_copy(power);
    }

    /* A product coefficient (e.g. Times[-1/6, 2 C[1] + 13 C[2]]) is spliced into
     * one flat Times with the power: the Plus printer reads a summand's sign off
     * its leading Times factor, so a nested Times[Times[-1/6, ...], x^3] would
     * print as "+ -1/6 (...) x^3" instead of "- 1/6 (...) x^3". */
    if (coef->type == EXPR_FUNCTION &&
        coef->data.function.head->type == EXPR_SYMBOL &&
        coef->data.function.head->data.symbol.name == SYM_Times) {
        size_t n = coef->data.function.arg_count;
        Expr** fargs = malloc((n + 1) * sizeof(Expr*));
        if (!fargs) return NULL;
        for (size_t i = 0; i < n; i++) fargs[i] = expr_copy(coef->data.function.args[i]);
        fargs[n] = expr_copy(power);
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times), fargs, n + 1);
        free(fargs);
        return t;
    }

    Expr* args[2] = { expr_copy(coef), expr_copy(power) };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* Renders e as the generic head[arg1, arg2, ...] form. Used as a fallback
 * when SeriesData has unexpected argument shapes. */
static void print_series_generic(Expr* e, int parent_prec) {
    (void)parent_prec;
    print_standard(e->data.function.head, 9500);
    printf("[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i > 0) printf(", ");
        print_standard(e->data.function.args[i], 0);
    }
    printf("]");
}

/* Builds the human-readable equivalent of
 * SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] as a fresh
 * Plus-expression of Times/Power/O nodes:
 *     a0 + a1 (x-x0) + ... + O[x-x0]^(nmax/den).
 * Returns a freshly-owned expression (caller frees), or NULL when the
 * arguments are not in a renderable shape (non-literal nmin/nmax/den or a
 * non-List coefficient container) so the caller can fall back to the generic
 * head[args...] form. Shared by the plain printer (print_series_data) and the
 * LaTeX printer (src/print_latex.c) so the two renderings cannot drift. */
Expr* series_data_to_display_expr(Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 6) return NULL;

    Expr* x        = e->data.function.args[0];
    Expr* x0       = e->data.function.args[1];
    Expr* coefs    = e->data.function.args[2];
    Expr* nmin_e   = e->data.function.args[3];
    Expr* nmax_e   = e->data.function.args[4];
    Expr* den_e    = e->data.function.args[5];

    bool coefs_ok = (coefs->type == EXPR_FUNCTION &&
                     coefs->data.function.head->type == EXPR_SYMBOL &&
                     coefs->data.function.head->data.symbol.name == SYM_List);
    bool ints_ok  = (nmin_e->type == EXPR_INTEGER &&
                     nmax_e->type == EXPR_INTEGER &&
                     den_e->type  == EXPR_INTEGER &&
                     den_e->data.integer != 0);
    if (!coefs_ok || !ints_ok) return NULL;

    int64_t nmin = nmin_e->data.integer;
    int64_t nmax = nmax_e->data.integer;
    int64_t den  = den_e->data.integer;
    size_t  k    = coefs->data.function.arg_count;

    Expr* base = series_build_base(x, x0);

    /* Collect non-zero summands plus the trailing O-term. */
    Expr** terms = malloc(sizeof(Expr*) * (k + 1));
    size_t tcount = 0;

    for (size_t i = 0; i < k; i++) {
        int64_t num = nmin + (int64_t)i;
        Expr* power = series_build_power(base, num, den);
        Expr* term  = series_build_term(coefs->data.function.args[i], power);
        if (power) expr_free(power);
        if (term) terms[tcount++] = term;
    }

    /* Trailing O[base]^(nmax/den) term. Always emitted, even when all
     * coefficients are zero, because it conveys the truncation order. */
    {
        Expr* oargs[1] = { expr_copy(base) };
        Expr* o_call   = expr_new_function(expr_new_symbol(SYM_O), oargs, 1);
        int64_t g  = gcd(abs_i64(nmax), abs_i64(den));
        if (g == 0) g = 1;
        int64_t n  = nmax / g;
        int64_t d  = den  / g;
        if (d < 0) { d = -d; n = -n; }
        Expr* exp_expr = (d == 1) ? expr_new_integer(n) : make_rational(n, d);
        Expr* pargs[2] = { o_call, exp_expr };
        terms[tcount++] = expr_new_function(expr_new_symbol(SYM_Power), pargs, 2);
    }

    Expr* sum;
    if (tcount == 1) {
        sum = terms[0];
    } else {
        sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, tcount);
    }
    free(terms);

    expr_free(base);
    return sum;
}

static void print_series_data(Expr* e, int parent_prec) {
    Expr* sum = series_data_to_display_expr(e);
    if (!sum) {
        print_series_generic(e, parent_prec);
        return;
    }
    print_standard(sum, parent_prec);
    expr_free(sum);
}

static bool is_series_data6(const Expr* a) {
    return a->type == EXPR_FUNCTION &&
           a->data.function.head->type == EXPR_SYMBOL &&
           a->data.function.head->data.symbol.name == SYM_SeriesData &&
           a->data.function.arg_count == 6;
}

/* Plus[..., SeriesData[...], ...] -- e.g. Log[x] + (series in 1/x), which
 * Series[PolyGamma[0, x], {x, Infinity, n}] returns -- printed as one flat
 * sum: the series' display terms are spliced into the outer Plus so each
 * summand's sign is read off by the Plus printer ("Log[x] - 1/(2 x) ...",
 * not "Log[x] + -1/(2 x) ..."). The terms are spliced in place (no
 * reordering), so a nested multivariate series keeps its inner O[y]^n next
 * to its own terms, as before. Returns false (nothing printed) when no
 * summand is a renderable SeriesData, so the caller prints the Plus
 * normally. */
static bool print_plus_spliced_series(Expr* e) {
    size_t n = e->data.function.arg_count;
    bool any = false;
    for (size_t i = 0; i < n && !any; i++)
        any = is_series_data6(e->data.function.args[i]);
    if (!any) return false;

    /* Display sums of the series summands (NULL for a non-series or an
     * unrenderable one, which is then kept as-is). */
    Expr** disp = calloc(n, sizeof(Expr*));
    if (!disp) return false;
    size_t total = 0;
    bool spliced = false;
    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        if (is_series_data6(a)) disp[i] = series_data_to_display_expr(a);
        Expr* d = disp[i];
        if (d && d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL &&
            d->data.function.head->data.symbol.name == SYM_Plus) {
            total += d->data.function.arg_count;
        } else {
            total += 1;
        }
        if (d) spliced = true;
    }
    if (!spliced) { free(disp); return false; }

    Expr** terms = malloc(total * sizeof(Expr*));
    if (!terms) {
        for (size_t i = 0; i < n; i++) if (disp[i]) expr_free(disp[i]);
        free(disp);
        return false;
    }
    size_t nt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* d = disp[i];
        if (!d) { terms[nt++] = expr_copy(e->data.function.args[i]); continue; }
        bool is_sum = d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL &&
                      d->data.function.head->data.symbol.name == SYM_Plus;
        if (is_sum) {
            for (size_t j = 0; j < d->data.function.arg_count; j++)
                terms[nt++] = expr_copy(d->data.function.args[j]);
        } else {
            terms[nt++] = expr_copy(d);
        }
        expr_free(d);
    }
    free(disp);

    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    /* The caller's parentheses (if any) are already open: print at 0. A
     * series coefficient that is itself a SeriesData (multivariate) is
     * spliced again one level down, so the recursion is bounded by the
     * nesting depth of the input. */
    print_standard(sum, 0);
    expr_free(sum);
    return true;
}

/* ===================== TeXForm renderer =====================
 *
 * Renders an Expr as AMS-LaTeX-compatible TeX. TeXForm[expr] in
 * print_standard dispatches here. Precedence uses the same numeric
 * levels as get_expr_prec so parens are inserted consistently.
 * Single-character symbol names are emitted bare (italicised by TeX
 * default math mode); multi-character names are wrapped in \text{}.
 */

/* Emit a number in TeX form (integers, reals, bigints, mpfr). */
static void print_tex_number(Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: printf("%" PRId64, e->data.integer); break;
        case EXPR_REAL: {
            /* %g's digits; an exponent is \times 10^{e} (1e+15 is not TeX) */
            char buf[64];
            snprintf(buf, sizeof(buf), "%g", e->data.real == 0.0 ? 0.0 : e->data.real);
            char* ep = isfinite(e->data.real) ? strpbrk(buf, "eE") : NULL;
            int ev = 0;
            if (ep) { ev = atoi(ep + 1); *ep = '\0'; }
            printf("%s", buf);
            if (!strchr(buf, '.') && isfinite(e->data.real)) putchar('.');
            if (ep) printf(ev >= 0 && ev < 10 ? "\\times 10^%d" : "\\times 10^{%d}", ev);
            break;
        }
        case EXPR_BIGINT: {
            char* s = mpz_get_str(NULL, 10, e->data.bigint);
            printf("%s", s);
            free(s);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: print_mpfr_display(e->data.mpfr, true); break;
#endif
        default: expr_print_fullform(e); break;
    }
}

/* A string inside \text{...}: the characters TeX treats specially are
 * escaped, so the output compiles ("a_b%" is \text{a\_b\%}). */
static void print_tex_text(const char* s) {
    printf("\\text{");
    for (const char* c = s; *c; c++) {
        switch (*c) {
            case '_': case '%': case '&': case '#': case '$': case '{': case '}':
                printf("\\%c", *c); break;
            case '~':  printf("\\textasciitilde{}"); break;
            case '^':  printf("\\textasciicircum{}"); break;
            case '\\': printf("\\textbackslash{}"); break;
            default:   putchar(*c); break;
        }
    }
    printf("}");
}

/* The TeX command of a Greek-letter symbol name (UTF-8), or NULL. */
static const char* tex_greek(const char* name) {
    static const char* const tab[][2] = {
        {"\xce\xb1", "\\alpha"}, {"\xce\xb2", "\\beta"}, {"\xce\xb3", "\\gamma"},
        {"\xce\xb4", "\\delta"}, {"\xce\xb5", "\\epsilon"}, {"\xcf\xb5", "\\epsilon"},
        {"\xce\xb6", "\\zeta"}, {"\xce\xb7", "\\eta"}, {"\xce\xb8", "\\theta"},
        {"\xcf\x91", "\\vartheta"}, {"\xce\xb9", "\\iota"}, {"\xce\xba", "\\kappa"},
        {"\xce\xbb", "\\lambda"}, {"\xce\xbc", "\\mu"}, {"\xce\xbd", "\\nu"},
        {"\xce\xbe", "\\xi"}, {"\xcf\x80", "\\pi"}, {"\xcf\x96", "\\varpi"},
        {"\xcf\x81", "\\rho"}, {"\xcf\xb1", "\\varrho"}, {"\xcf\x83", "\\sigma"},
        {"\xcf\x82", "\\varsigma"}, {"\xcf\x84", "\\tau"}, {"\xcf\x85", "\\upsilon"},
        {"\xcf\x86", "\\phi"}, {"\xcf\x95", "\\phi"}, {"\xcf\x87", "\\chi"},
        {"\xcf\x88", "\\psi"}, {"\xcf\x89", "\\omega"},
        {"\xce\x93", "\\Gamma"}, {"\xce\x94", "\\Delta"}, {"\xce\x98", "\\Theta"},
        {"\xce\x9b", "\\Lambda"}, {"\xce\x9e", "\\Xi"}, {"\xce\xa0", "\\Pi"},
        {"\xce\xa3", "\\Sigma"}, {"\xce\xa5", "\\Upsilon"}, {"\xce\xa6", "\\Phi"},
        {"\xce\xa8", "\\Psi"}, {"\xce\xa9", "\\Omega"},
    };
    for (size_t i = 0; i < sizeof tab / sizeof tab[0]; i++)
        if (strcmp(name, tab[i][0]) == 0) return tab[i][1];
    return NULL;
}

static void print_tex_symbol(const char* raw_name) {
    const char* name = context_display_name(raw_name);
    const char* greek = tex_greek(name);
    if (greek) { printf("%s ", greek); return; }
    if (strcmp(name, "Pi") == 0)              { printf("\\pi "); return; }
    if (strcmp(name, "E") == 0)               { printf("e"); return; }
    if (strcmp(name, "I") == 0)               { printf("i"); return; }
    if (strcmp(name, "Infinity") == 0)        { printf("\\infty "); return; }
    if (strcmp(name, "ComplexInfinity") == 0) { printf("\\tilde{\\infty }"); return; }
    if (strcmp(name, "Indeterminate") == 0)   { printf("\\text{Indeterminate}"); return; }
    if (strcmp(name, "EulerGamma") == 0)      { printf("\\gamma "); return; }
    if (strcmp(name, "GoldenRatio") == 0)     { printf("\\phi "); return; }
    if (strcmp(name, "Catalan") == 0)         { printf("C"); return; }
    if (strcmp(name, "Degree") == 0)          { printf("{}^{\\circ }"); return; }
    if (strcmp(name, "True") == 0)            { printf("\\text{True}"); return; }
    if (strcmp(name, "False") == 0)           { printf("\\text{False}"); return; }
    if (strcmp(name, "Null") == 0)            { printf("\\text{Null}"); return; }
    size_t len = strlen(name);
    if (len == 1) printf("%s", name);
    else print_tex_text(name);
}

/* Looks up a trig/hyperbolic head and returns its TeX command and
 * whether it's an inverse variant. Returns false if the head is not a
 * standard trig/hyperbolic function. */
static bool tex_trig_lookup(const char* h, const char** tex, bool* inv) {
    static const struct { const char* head; const char* tex; bool inv; } table[] = {
        {"Sin",  "\\sin",  false}, {"Cos",  "\\cos",  false},
        {"Tan",  "\\tan",  false}, {"Cot",  "\\cot",  false},
        {"Sec",  "\\sec",  false}, {"Csc",  "\\csc",  false},
        {"Sinh", "\\sinh", false}, {"Cosh", "\\cosh", false},
        {"Tanh", "\\tanh", false}, {"Coth", "\\coth", false},
        {"Sech", "\\text{sech}", false}, {"Csch", "\\text{csch}", false},
        {"ArcSin", "\\sin", true}, {"ArcCos", "\\cos", true},
        {"ArcTan", "\\tan", true}, {"ArcCot", "\\cot", true},
        {"ArcSec", "\\sec", true}, {"ArcCsc", "\\csc", true},
        {"ArcSinh","\\sinh", true}, {"ArcCosh","\\cosh", true},
        {"ArcTanh","\\tanh", true}, {"ArcCoth","\\coth", true},
        {"ArcSech","\\text{sech}", true}, {"ArcCsch","\\text{csch}", true},
    };
    size_t n = sizeof(table) / sizeof(table[0]);
    for (size_t i = 0; i < n; i++) {
        if (strcmp(h, table[i].head) == 0) {
            *tex = table[i].tex;
            *inv = table[i].inv;
            return true;
        }
    }
    return false;
}

/* Maps a named head to a TeX symbol when the function is to be rendered
 * like f(x). Returns NULL if the head should use its default name. */
static const char* tex_function_name(const char* h) {
    if (h == SYM_Log)      return "\\log ";
    if (h == SYM_Exp)      return "\\exp ";
    if (h == SYM_Sqrt)     return NULL; /* handled via \sqrt{} */
    if (h == SYM_Abs)      return NULL; /* handled via |...|  */
    if (h == SYM_Floor)    return NULL; /* handled via \lfloor */
    if (h == SYM_Ceiling)  return NULL; /* handled via \lceil  */
    if (h == SYM_Gamma)    return "\\Gamma ";
    if (h == SYM_Re)       return "\\Re ";
    if (h == SYM_Zeta)     return "\\zeta ";
    if (h == SYM_Im)       return "\\Im ";
    if (h == SYM_Conjugate)return NULL; /* handled via \overline */
    return NULL;
}

/* A subscript: a digit or a one-letter symbol bare (J_n, x_1), anything
 * else braced (J_{10}, x_{n+1}). */
static void tex_subscript(Expr* sub) {
    bool bare = (sub->type == EXPR_INTEGER && sub->data.integer >= 0 && sub->data.integer <= 9)
             || (sub->type == EXPR_SYMBOL && strlen(context_display_name(sub->data.symbol.name)) == 1);
    if (!bare) printf("{");
    print_tex(sub, 0);
    if (!bare) printf("}");
}

/* Is e a matrix: a non-empty List of equally long, non-empty Lists none of
 * whose elements is a List? */
static bool tex_is_matrix(const Expr* e) {
    size_t n = e->data.function.arg_count, cols = 0;
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++) {
        const Expr* r = e->data.function.args[i];
        if (r->type != EXPR_FUNCTION || r->data.function.head->type != EXPR_SYMBOL ||
            r->data.function.head->data.symbol.name != SYM_List)
            return false;
        if (i == 0) cols = r->data.function.arg_count;
        if (cols == 0 || r->data.function.arg_count != cols) return false;
        for (size_t j = 0; j < cols; j++) {
            const Expr* x = r->data.function.args[j];
            if ((x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL &&
                 x->data.function.head->data.symbol.name == SYM_List) || x->type == EXPR_NDARRAY)
                return false;
        }
    }
    return true;
}

/* n when h is Derivative[n][f] with one positive machine integer n, else 0. */
static int64_t tex_derivative_order(const Expr* h) {
    if (h->type != EXPR_FUNCTION || h->data.function.arg_count != 1) return 0;
    const Expr* d = h->data.function.head;
    if (d->type != EXPR_FUNCTION || d->data.function.head->type != EXPR_SYMBOL ||
        d->data.function.head->data.symbol.name != SYM_Derivative || d->data.function.arg_count != 1)
        return 0;
    const Expr* n = d->data.function.args[0];
    return (n->type == EXPR_INTEGER && n->data.integer > 0) ? n->data.integer : 0;
}

/* Emits a comma-separated argument list wrapped in (...). */
static void print_tex_args_parens(Expr* e, size_t start) {
    printf("(");
    for (size_t i = start; i < e->data.function.arg_count; i++) {
        if (i > start) printf(",");
        print_tex(e->data.function.args[i], 0);
    }
    printf(")");
}

/* Detects whether `im` is a negative numeric, and if so produces a
 * freshly-owned Expr holding its absolute value (caller must free). */
static bool tex_extract_negative(Expr* im, Expr** abs_out) {
    int64_t rn, rd;
    if (im->type == EXPR_INTEGER && im->data.integer < 0) {
        *abs_out = print_neg_int(im->data.integer);
        return true;
    }
    if (im->type == EXPR_REAL && im->data.real < 0.0) {
        *abs_out = expr_new_real(-im->data.real);
        return true;
    }
    if (is_rational(im, &rn, &rd) && rn < 0) {
        *abs_out = make_rational_neg(rn, rd);
        return true;
    }
    return false;
}

/* Does the TeX rendering of `e` begin with a digit?  A nonnegative number
 * does, and so does a Power whose base is one (3^{x}).  Two such factors
 * side by side need a visible separator: math mode ignores spaces, so
 * "2 3^{x}" typesets as 23^{x} (audit TeXForm#1). */
static bool tex_leads_with_digit(const Expr* e) {
    for (;;) {
        if (e->type == EXPR_INTEGER) return e->data.integer >= 0;
        if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) >= 0;
        if (e->type == EXPR_REAL) return e->data.real >= 0.0;
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_Power
            && e->data.function.arg_count == 2) {
            e = e->data.function.args[0];
            continue;
        }
        return false;
    }
}

/* The separator between two adjacent Times factors in TeX: an explicit
 * space before a factor that leads with a digit, else a plain one. */
static void print_tex_factor_sep(const Expr* next) {
    printf(tex_leads_with_digit(next) ? "\\ " : " ");
}

static void print_tex(Expr* e, int parent_prec) {
    if (!e) { printf("\\text{Null}"); return; }

    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
    ) {
        print_tex_number(e);
        return;
    }
    if (e->type == EXPR_SYMBOL) { print_tex_symbol(e->data.symbol.name); return; }
    if (e->type == EXPR_STRING) { print_tex_text(e->data.string); return; }
    if (e->type == EXPR_NDARRAY) {
        if (is_packed_list(e)) {
            /* Renders as the List it is. Materialises rather than streaming,
             * unlike print_standard: TeX output is for display, so it is never
             * asked for on the 10^6-element values where streaming matters,
             * and print_tex's List branch handles nesting and separators. */
            Expr* nested = ndarray_to_nested_list(e);
            print_tex(nested, 0);
            expr_free(nested);
            return;
        }
        printf("\\text{NDArray}\\left[");
        Expr* nested = ndarray_to_nested_list(e);
        print_tex(nested, 0);
        expr_free(nested);
        printf("\\right]");
        return;
    }
    if (e->type == EXPR_SPARSE) {
        /* Typeset as the matrix / vector it represents. */
        Expr* dense = sp_normal(e->data.sparse);
        if (dense) { print_tex(dense, parent_prec); expr_free(dense); }
        else printf("\\text{SparseArray}");
        return;
    }
    if (e->type == EXPR_PVEC || e->type == EXPR_LINKEDLIST) {
        const DsOps* o = ds_ops(e->type);
        printf("\\text{%s}\\left[\\langle %" PRId64 "\\rangle \\right]", o->name, o->length(e));
        return;
    }
    if (e->type == EXPR_COMPILED) {
        const CompiledFunction* cf = e->data.compiled;
        size_t n = compiled_function_num_args(cf);
        const char* const* names = compiled_function_arg_names(cf);
        printf("\\text{CompiledFunction}\\left[\\{");
        for (size_t i = 0; i < n; i++) { if (i) printf(", "); print_tex_symbol(names[i]); }
        printf("\\}, ");
        print_tex((Expr*)compiled_function_body(cf), 0);
        printf("\\right]");
        return;
    }
    if (e->type != EXPR_FUNCTION) return;

    /* Display rules: Format[f[...], TeXForm] first, then a plain
     * Format[f[...]] (TeXForm is a display form). TeX is never fed back into
     * evaluation, so they always apply. */
    if (symtab_format_value_count && g_format_depth < FORMAT_DEPTH_LIMIT) {
        Expr* shown = format_value_rewrite(e, SYM_TeXForm);
        if (!shown) shown = format_value_rewrite(e, NULL);
        if (shown) {
            g_format_depth++;
            print_tex(shown, parent_prec);
            g_format_depth--;
            expr_free(shown);
            return;
        }
    }
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_SequenceForm) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* a = e->data.function.args[i];
            if (a->type == EXPR_STRING) print_tex_text(a->data.string);
            else print_tex(a, 0);
        }
        return;
    }
    /* HoldForm and the display wrappers TraditionalForm / StandardForm are
     * invisible in TeX, as in every other printer: TeXForm[HoldForm[1 + 1]]
     * is 1+1 (audit HoldForm#1, TraditionalForm#1). */
    if (e->data.function.head->type == EXPR_SYMBOL && e->data.function.arg_count == 1 &&
        (e->data.function.head->data.symbol.name == SYM_HoldForm ||
         e->data.function.head->data.symbol.name == SYM_TraditionalForm ||
         e->data.function.head->data.symbol.name == SYM_StandardForm)) {
        print_tex(e->data.function.args[0], parent_prec);
        return;
    }

    /* Pure functions in the language's own notation: Slot[n] is \#n, a named
     * slot \#\text{name}, SlotSequence[n] \#\#n, and Function[body] is body\&.
     * They were printed as calls -- \text{Root}(\text{Function}(\text{Slot}(1)^{3}
     * +2),1) -- which is how every Root and RootSum in an Integrate answer read.
     * A Function with named variables keeps the call form. */
    if (e->data.function.head->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* fh = e->data.function.head->data.symbol.name;
        Expr* a0 = e->data.function.args[0];
        if ((fh == SYM_Slot || fh == SYM_SlotSequence) &&
            (a0->type == EXPR_INTEGER || a0->type == EXPR_STRING)) {
            printf(fh == SYM_Slot ? "\\#" : "\\#\\#");
            if (a0->type == EXPR_INTEGER) printf("%" PRId64, a0->data.integer);
            else print_tex_text(a0->data.string);
            return;
        }
        if (fh == SYM_Function) {
            /* & binds loosest of all: parenthesised inside any operator */
            if (parent_prec > 0) printf("\\left(");
            print_tex(a0, 0);
            printf("\\&");
            if (parent_prec > 0) printf("\\right)");
            return;
        }
    }

    /* The big operators -- an unevaluated (held) Integrate / Sum / Product /
     * Limit -- typeset as the notation, not as a function call:
     * \int_0^1 f(x) \, dx, \sum _{i=1}^n a(i), \lim_{x\to 0} f(x).  Each binds
     * like a sum, so it is parenthesised inside a product or power. */
    if (e->data.function.head->type == EXPR_SYMBOL && e->data.function.arg_count == 2) {
        const char* bh = e->data.function.head->data.symbol.name;
        Expr* body = e->data.function.args[0];
        Expr* spec = e->data.function.args[1];
        bool is_list = spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
                       spec->data.function.head->data.symbol.name == SYM_List;
        size_t sn = is_list ? spec->data.function.arg_count : 0;
        bool rule = spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
                    spec->data.function.head->data.symbol.name == SYM_Rule &&
                    spec->data.function.arg_count == 2;
        int kind = 0;
        if (bh == SYM_Integrate && (spec->type == EXPR_SYMBOL || (is_list && sn == 3))) kind = 1;
        else if ((bh == SYM_Sum || bh == SYM_Product) && is_list && (sn == 2 || sn == 3)) kind = 2;
        else if (bh == SYM_Limit && rule) kind = 3;
        if (kind) {
            bool paren = parent_prec > 3500;
            if (paren) printf("\\left(");
            if (kind == 1) {
                Expr* var = is_list ? spec->data.function.args[0] : spec;
                printf("\\int");
                if (is_list) {
                    printf("_{"); print_tex(spec->data.function.args[1], 0);
                    printf("}^{"); print_tex(spec->data.function.args[2], 0); printf("}");
                }
                printf(" ");
                print_tex(body, 3500);
                printf(" \\, d");
                print_tex(var, 10000);
            } else if (kind == 2) {
                printf(bh == SYM_Sum ? "\\sum _{" : "\\prod _{");
                print_tex(spec->data.function.args[0], 0);
                printf("=");
                if (sn == 3) print_tex(spec->data.function.args[1], 0);
                else printf("1");
                printf("}^{");
                print_tex(spec->data.function.args[sn - 1], 0);
                printf("} ");
                print_tex(body, 3500);
            } else {
                printf("\\lim_{");
                print_tex(spec->data.function.args[0], 0);
                printf("\\to ");
                print_tex(spec->data.function.args[1], 0);
                printf("} ");
                print_tex(body, 3500);
            }
            if (paren) printf("\\right)");
            return;
        }
    }

    int my_prec = get_expr_prec(e);
    bool need_parens = (my_prec < parent_prec);
    if (need_parens) printf("\\left(");

    Expr* head_expr = e->data.function.head;
    if (head_expr->type == EXPR_SYMBOL) {
        const char* head = head_expr->data.symbol.name;
        size_t argc = e->data.function.arg_count;

        if (head == SYM_Rational && argc == 2) {
            printf("\\frac{");
            print_tex(e->data.function.args[0], 0);
            printf("}{");
            print_tex(e->data.function.args[1], 0);
            printf("}");
        }
        else if (head == SYM_Complex && argc == 2) {
            Expr* re = e->data.function.args[0];
            Expr* im = e->data.function.args[1];
            bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0);
            bool im_one = (im->type == EXPR_INTEGER && im->data.integer == 1);
            bool im_minus_one = (im->type == EXPR_INTEGER && im->data.integer == -1);
            Expr* im_abs = NULL;
            bool im_neg = im_minus_one;
            if (!im_neg) im_neg = tex_extract_negative(im, &im_abs);
            if (!re_zero) {
                print_tex(re, 3500);
                printf(im_neg ? "-" : "+");
            } else if (im_neg) {
                printf("-");
            }
            if (im_one || im_minus_one) printf("i");
            else if (im_abs) { print_tex(im_abs, 4500); printf(" i"); expr_free(im_abs); }
            else { print_tex(im, 4500); printf(" i"); }
        }
        else if (head == SYM_Sqrt && argc == 1) {
            printf("\\sqrt{");
            print_tex(e->data.function.args[0], 0);
            printf("}");
        }
        else if (head == SYM_Abs && argc == 1) {
            printf("\\left| ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right| ");
        }
        else if (head == SYM_Floor && argc == 1) {
            printf("\\left\\lfloor ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right\\rfloor ");
        }
        else if (head == SYM_Ceiling && argc == 1) {
            printf("\\left\\lceil ");
            print_tex(e->data.function.args[0], 0);
            printf("\\right\\rceil ");
        }
        else if (head == SYM_Conjugate && argc == 1) {
            printf("\\overline{");
            print_tex(e->data.function.args[0], 0);
            printf("}");
        }
        else if (head == SYM_Power && argc == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp = e->data.function.args[1];
            int64_t n, d;
            /* Negative exponent -> 1/base^|exp| with \frac */
            bool exp_neg = false;
            Expr* pos_exp = NULL;
            if ((pos_exp = print_neg_int_exp(exp)) != NULL) {
                exp_neg = true;
            } else if (is_rational(exp, &n, &d) && n < 0) {
                exp_neg = true; pos_exp = make_rational_neg(n, d);
            } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                exp_neg = true; pos_exp = expr_new_real(-exp->data.real);
            }
            if (exp_neg) {
                printf("\\frac{1}{");
                if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                    print_tex(base, 0);
                } else {
                    Expr* pargs[2] = { expr_copy(base), expr_copy(pos_exp) };
                    Expr* tmp = expr_new_function(expr_new_symbol(SYM_Power), pargs, 2);
                    print_tex(tmp, 0);
                    expr_free(tmp);
                }
                printf("}");
                expr_free(pos_exp);
            } else if (is_rational(exp, &n, &d) && n == 1 && d == 2) {
                printf("\\sqrt{");
                print_tex(base, 0);
                printf("}");
            } else if (is_rational(exp, &n, &d) && n == 1 && d > 2) {
                printf("\\sqrt[%" PRId64 "]{", d);
                print_tex(base, 0);
                printf("}");
            } else {
                /* Negative numeric atoms parse as unary-minus, so render
                 * `-1^{1/3}` as `(-1)^{1/3}` to keep the meaning unambiguous. */
                bool base_paren = false;
                int64_t bn, bd;
                if (base->type == EXPR_INTEGER && base->data.integer < 0) base_paren = true;
                else if (base->type == EXPR_REAL && base->data.real < 0.0) base_paren = true;
                else if (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) < 0) base_paren = true;
                else if (is_rational(base, &bn, &bd) && bn < 0) base_paren = true;
                if (base_paren) {
                    printf("\\left(");
                    print_tex(base, 0);
                    printf("\\right)");
                } else {
                    print_tex(base, 6700); /* force parens for compound bases */
                }
                printf("^{");
                print_tex(exp, 0);
                printf("}");
            }
        }
        else if (head == SYM_Times) {
            size_t count = argc;
            int64_t lead_start = 0;
            bool flipped_sign = false;
            Expr* flipped_head = NULL;
            if (count > 0) {
                Expr* a0 = e->data.function.args[0];
                if (a0->type == EXPR_INTEGER && a0->data.integer < 0) {
                    flipped_sign = true;
                    if (a0->data.integer == -1 && count > 1) lead_start = 1;
                    else flipped_head = print_neg_int(a0->data.integer);
                } else if (a0->type == EXPR_REAL && a0->data.real < 0.0) {
                    flipped_sign = true;
                    if (a0->data.real == -1.0 && count > 1) lead_start = 1;
                    else flipped_head = expr_new_real(-a0->data.real);
                } else {
                    int64_t rn, rd;
                    if (is_rational(a0, &rn, &rd) && rn < 0) {
                        flipped_sign = true;
                        flipped_head = make_rational_neg(rn, rd);
                    }
                }
            }
            if (flipped_sign) printf("-");

            /* Collect numerator / denominator via negative-exponent Powers,
             * matching print_standard so a/b renders as \frac{a}{b}. */
            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t num_count = 0, den_count = 0;
            for (size_t i = (size_t)lead_start; i < count; i++) {
                Expr* arg = (i == 0 && flipped_head) ? flipped_head : e->data.function.args[i];
                bool den_ok = false;
                if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL
                    && arg->data.function.head->data.symbol.name == SYM_Power
                    && arg->data.function.arg_count == 2) {
                    Expr* exp = arg->data.function.args[1];
                    int64_t n, d;
                    Expr* pos_exp = NULL;
                    bool neg = false;
                    if ((pos_exp = print_neg_int_exp(exp)) != NULL) {
                        neg = true;
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        neg = true; pos_exp = make_rational_neg(n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        neg = true; pos_exp = expr_new_real(-exp->data.real);
                    } else if ((pos_exp = print_neg_symbolic_exp(exp)) != NULL) {
                        neg = true;
                    }
                    if (neg) {
                        if (pos_exp->type == EXPR_INTEGER && pos_exp->data.integer == 1) {
                            den_args[den_count++] = expr_copy(arg->data.function.args[0]);
                            expr_free(pos_exp);
                        } else {
                            Expr* p_args[2] = { expr_copy(arg->data.function.args[0]), pos_exp };
                            den_args[den_count++] = expr_new_function(expr_new_symbol(SYM_Power), p_args, 2);
                        }
                        den_ok = true;
                    }
                }
                if (!den_ok) num_args[num_count++] = expr_copy(arg);
            }

            fold_rational_coeff_into_den(num_args, &num_count, den_args, &den_count);

            if (den_count == 0) {
                for (size_t i = 0; i < num_count; i++) {
                    if (i > 0) print_tex_factor_sep(num_args[i]);
                    print_tex(num_args[i], 4500);
                }
            } else {
                /* Context precedence inside a \frac slot.  0 ("never
                 * parenthesise") is right only while the factor is the slot's
                 * SOLE occupant: the brace group is then the grouping, and
                 * \frac{1}{1+x} needs no brackets.  With two or more factors
                 * they sit side by side in one brace group, joined by implicit
                 * multiplication, so each must be rendered at the `*`
                 * precedence -- otherwise a Plus among them loses its
                 * parentheses and the LaTeX states a DIFFERENT expression:
                 *   a/(b (c + d))            ->  \frac{a}{b c+d}   = (bc+d)
                 *   1/(Sqrt[10] (1-Sqrt[5])) ->  \frac{1}{\sqrt{10} 1-\sqrt{5}}
                 * Both were wrong answers in the typeset output, not merely
                 * ugly.  The den_count == 0 branch above already passes 4500
                 * for exactly this reason. */
                printf("\\frac{");
                if (num_count == 0) { printf("1"); }
                else {
                    int nprec = (num_count > 1) ? 4500 : 0;
                    for (size_t i = 0; i < num_count; i++) {
                        if (i > 0) print_tex_factor_sep(num_args[i]);
                        print_tex(num_args[i], nprec);
                    }
                }
                printf("}{");
                {
                    int dprec = (den_count > 1) ? 4500 : 0;
                    for (size_t i = 0; i < den_count; i++) {
                        if (i > 0) print_tex_factor_sep(den_args[i]);
                        print_tex(den_args[i], dprec);
                    }
                }
                printf("}");
            }
            for (size_t i = 0; i < num_count; i++) expr_free(num_args[i]);
            for (size_t i = 0; i < den_count; i++) expr_free(den_args[i]);
            free(num_args);
            free(den_args);
            if (flipped_head) expr_free(flipped_head);
        }
        else if (head == SYM_Plus) {
            for (size_t i = 0; i < argc; i++) {
                Expr* arg = e->data.function.args[i];
                bool is_neg = false;
                Expr* to_print = arg;
                Expr* owned = NULL;

                if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
                    is_neg = true; owned = print_neg_int(arg->data.integer); to_print = owned;
                } else if (arg->type == EXPR_BIGINT && mpz_sgn(arg->data.bigint) < 0) {
                    is_neg = true; owned = negate_bigint(arg); to_print = owned;
                } else if (arg->type == EXPR_REAL && arg->data.real < 0.0) {
                    is_neg = true; owned = expr_new_real(-arg->data.real); to_print = owned;
                } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL
                           && arg->data.function.head->data.symbol.name == SYM_Times
                           && arg->data.function.arg_count > 0) {
                    Expr* f0 = arg->data.function.args[0];
                    if (f0->type == EXPR_INTEGER && f0->data.integer < 0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = print_neg_int(f0->data.integer);
                        if (owned->data.function.args[0]->data.integer == 1 && owned->data.function.arg_count > 1) {
                            Expr** na = malloc(sizeof(Expr*) * (owned->data.function.arg_count - 1));
                            for (size_t k = 1; k < owned->data.function.arg_count; k++) na[k-1] = expr_copy(owned->data.function.args[k]);
                            Expr* t = expr_new_function(expr_new_symbol(SYM_Times), na, owned->data.function.arg_count - 1);
                            expr_free(owned); owned = t; free(na);
                        }
                        to_print = owned;
                    } else if (f0->type == EXPR_BIGINT && mpz_sgn(f0->data.bigint) < 0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = negate_bigint(f0);
                        if (owned->data.function.args[0]->type == EXPR_INTEGER &&
                            owned->data.function.args[0]->data.integer == 1 &&
                            owned->data.function.arg_count > 1) {
                            Expr** na = malloc(sizeof(Expr*) * (owned->data.function.arg_count - 1));
                            for (size_t k = 1; k < owned->data.function.arg_count; k++) na[k-1] = expr_copy(owned->data.function.args[k]);
                            Expr* t = expr_new_function(expr_new_symbol(SYM_Times), na, owned->data.function.arg_count - 1);
                            expr_free(owned); owned = t; free(na);
                        }
                        to_print = owned;
                    } else if (f0->type == EXPR_REAL && f0->data.real < 0.0) {
                        is_neg = true;
                        owned = expr_unshare(expr_copy(arg));
                        expr_free(owned->data.function.args[0]);
                        owned->data.function.args[0] = expr_new_real(-f0->data.real);
                        to_print = owned;
                    } else {
                        int64_t rn, rd;
                        if (is_rational(f0, &rn, &rd) && rn < 0) {
                            is_neg = true;
                            owned = expr_unshare(expr_copy(arg));
                            Expr* new_rat = make_rational_neg(rn, rd);
                            expr_free(owned->data.function.args[0]);
                            owned->data.function.args[0] = new_rat;
                            to_print = owned;
                        }
                    }
                } else {
                    int64_t rn, rd;
                    if (is_rational(arg, &rn, &rd) && rn < 0) {
                        is_neg = true; owned = make_rational_neg(rn, rd); to_print = owned;
                    }
                }

                if (i > 0) printf(is_neg ? "-" : "+");
                else if (is_neg) printf("-");

                print_tex(to_print, 3500);
                if (owned) expr_free(owned);
            }
        }
        else if (head == SYM_List && tex_is_matrix(e)) {
            /* A rectangular list of lists of non-lists is a matrix: an array
             * in parentheses, one row per line (audit TeXForm#7). */
            size_t cols = e->data.function.args[0]->data.function.arg_count;
            printf("\\left(\n\\begin{array}{");
            for (size_t j = 0; j < cols; j++) printf("c");
            printf("}\n");
            for (size_t i = 0; i < argc; i++) {
                Expr* row = e->data.function.args[i];
                printf(" ");
                for (size_t j = 0; j < cols; j++) {
                    if (j > 0) printf(" & ");
                    print_tex(row->data.function.args[j], 0);
                }
                printf(" \\\\\n");
            }
            printf("\\end{array}\n\\right)");
        }
        else if (head == SYM_List) {
            printf("\\{");
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf(",");
                print_tex(e->data.function.args[i], 0);
            }
            printf("\\}");
        }
        else if (head == SYM_Association && assoc_prints_special(e)) {
            printf("\\left\\langle\\!\\left|");
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf(",\\,");
                print_tex(e->data.function.args[i], 0);
            }
            printf("\\right|\\!\\right\\rangle");
        }
        else if (head == SYM_Inequality && argc >= 3 && (argc & 1u) == 1) {
            for (size_t i = 0; i < argc; i++) {
                Expr* a = e->data.function.args[i];
                if ((i & 1u) == 0) {
                    print_tex(a, my_prec);
                } else {
                    const char* op = "";
                    if (a->type == EXPR_SYMBOL) {
                        const char* s = a->data.symbol.name;
                        if      (s == SYM_Equal)        op = "=";
                        else if (s == SYM_Less)         op = "<";
                        else if (s == SYM_Greater)      op = ">";
                        else if (s == SYM_LessEqual)    op = "\\leq ";
                        else if (s == SYM_GreaterEqual) op = "\\geq ";
                        else if (s == SYM_Unequal)      op = "\\neq ";
                    }
                    printf("%s", op);
                }
            }
        }
        else if ((head == SYM_Equal || head == SYM_Unequal
                  || head == SYM_Less || head == SYM_Greater
                  || head == SYM_LessEqual || head == SYM_GreaterEqual
                  || head == SYM_SameQ || head == SYM_UnsameQ
                  || head == SYM_Rule || head == SYM_RuleDelayed
                  || head == SYM_Set || head == SYM_SetDelayed
                  || head == SYM_And || head == SYM_Or) && argc >= 2) {
            const char* op = "";
            if (head == SYM_Equal) op = "=";
            else if (head == SYM_Unequal) op = "\\neq ";
            else if (head == SYM_Less) op = "<";
            else if (head == SYM_Greater) op = ">";
            else if (head == SYM_LessEqual) op = "\\leq ";
            else if (head == SYM_GreaterEqual) op = "\\geq ";
            else if (head == SYM_SameQ) op = "\\equiv ";
            else if (head == SYM_UnsameQ) op = "\\not\\equiv ";
            else if (head == SYM_Rule) op = "\\to ";
            /* `:\to`, not `\to`: a delayed rule is a different object from an
             * immediate one, and one glyph for both makes `a -> b` and `a :> b`
             * indistinguishable once typeset. Mirrors SetDelayed's `:=`. */
            else if (head == SYM_RuleDelayed) op = ":\\to ";
            else if (head == SYM_Set) op = "=";
            else if (head == SYM_SetDelayed) op = ":=";
            else if (head == SYM_And) op = "\\land ";
            else if (head == SYM_Or)  op = "\\lor ";
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf("%s", op);
                print_tex(e->data.function.args[i], my_prec);
            }
        }
        else if (head == SYM_Not && argc == 1) {
            printf("\\neg ");
            print_tex(e->data.function.args[0], 9500);
        }
        /* The standard functions with a notation of their own (audit TeXForm#6). */
        else if (head == SYM_Factorial && argc == 1) {
            print_tex(e->data.function.args[0], 9500);
            printf("!");
        }
        else if (head == SYM_Binomial && argc == 2) {
            printf("\\binom{");
            print_tex(e->data.function.args[0], 0);
            printf("}{");
            print_tex(e->data.function.args[1], 0);
            printf("}");
        }
        else if ((head == SYM_BesselJ || head == SYM_BesselY || head == SYM_BesselI ||
                  head == SYM_BesselK) && argc == 2) {
            printf("%c_", head == SYM_BesselJ ? 'J' : head == SYM_BesselY ? 'Y'
                        : head == SYM_BesselI ? 'I' : 'K');
            tex_subscript(e->data.function.args[0]);
            printf("(");
            print_tex(e->data.function.args[1], 0);
            printf(")");
        }
        else if (strcmp(context_display_name(head), "Subscript") == 0 && argc >= 2) {
            print_tex(e->data.function.args[0], 9500);
            printf("_");
            if (argc == 2) tex_subscript(e->data.function.args[1]);
            else {
                printf("{");
                for (size_t i = 1; i < argc; i++) {
                    if (i > 1) printf(",");
                    print_tex(e->data.function.args[i], 0);
                }
                printf("}");
            }
        }
        else if (head == SYM_Dot && argc >= 2) {
            for (size_t i = 0; i < argc; i++) {
                if (i > 0) printf(".");
                print_tex(e->data.function.args[i], 9500);
            }
        }
        else if (strcmp(context_display_name(head), "C") == 0 && argc == 1) {
            /* The DSolve/Reduce/Integrate generated constant C[k] renders as
             * the subscripted lowercase c_k, matching Mathematica's TeXForm
             * (C[1] -> c_1). A single-character subscript stays bare (c_1,
             * c_n); anything longer is braced (c_{10}, c_{n+1}) so LaTeX
             * groups the whole subscript rather than only its first token. */
            Expr* sub = e->data.function.args[0];
            bool bare = (sub->type == EXPR_INTEGER
                         && sub->data.integer >= 0 && sub->data.integer <= 9)
                     || (sub->type == EXPR_SYMBOL
                         && strlen(context_display_name(sub->data.symbol.name)) == 1);
            printf("c_");
            if (!bare) printf("{");
            print_tex(sub, 0);
            if (!bare) printf("}");
        }
        else {
            const char* tex_fn = NULL;
            bool is_inv = false;
            if (tex_trig_lookup(head, &tex_fn, &is_inv) && argc == 1) {
                /* \sin(x) or \sin ^{-1}(x) */
                printf("%s", tex_fn);
                if (is_inv) printf(" ^{-1}");
                print_tex_args_parens(e, 0);
            } else {
                const char* named = tex_function_name(head);
                if (named) {
                    printf("%s", named);
                    print_tex_args_parens(e, 0);
                } else {
                    /* Generic: f(a, b). Head gets symbol treatment so single-
                     * char heads render italic, multi-char heads use \text{}. */
                    print_tex(head_expr, 9500);
                    print_tex_args_parens(e, 0);
                }
            }
        }
    } else if (tex_derivative_order(head_expr) > 0) {
        /* Derivative[n][f][x]: f(x) with n primes up to three, then f^{(n)}(x). */
        int64_t n = tex_derivative_order(head_expr);
        print_tex(head_expr->data.function.args[0], 9500);
        if (n <= 3) for (int64_t k = 0; k < n; k++) printf("'");
        else printf("^{(%" PRId64 ")}", n);
        print_tex_args_parens(e, 0);
    } else {
        /* Non-symbol head: render as (head)(args...). */
        print_tex(head_expr, 9500);
        print_tex_args_parens(e, 0);
    }

    if (need_parens) printf("\\right)");
}

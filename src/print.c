#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// Print

#include "print.h"
#include "expr.h"
#include "arithmetic.h"
#include "context.h"
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
static Expr* print_neg_int(int64_t v) {
    if (v != INT64_MIN) return expr_new_integer(-v);
    mpz_t z;
    mpz_init_set_si(z, v);
    mpz_neg(z, z);
    Expr* e = expr_new_bigint_from_mpz(z);
    mpz_clear(z);
    return e;
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

static int get_expr_prec(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 9500;
    if (e->data.function.head->type != EXPR_SYMBOL) return 9500;
    
    const char* head = e->data.function.head->data.symbol.name;

    if (head == SYM_Set || head == SYM_SetDelayed) return 500;
#if UP_VALUES
    if (head == SYM_UpSet || head == SYM_UpSetDelayed ||
        head == SYM_TagSet || head == SYM_TagSetDelayed || head == SYM_TagUnset) return 500;
#endif
    if (head == SYM_MessageName) return 8600;
    if (head == SYM_Rule || head == SYM_RuleDelayed) return 1500;
    if (head == SYM_DirectedEdge || head == SYM_UndirectedEdge) return 1500;
    if (head == SYM_Condition) return 1700;
    if (head == SYM_Alternatives) return 2300;
    if (head == SYM_StringExpression) return 2100;   /* a ~~ b, below `|` */
    if ((head == SYM_Repeated || head == SYM_RepeatedNull) &&
        e->data.function.arg_count == 1) return 2500;
    /* x:p (Pattern whose body is not a plain blank; x_ / x__ print atomically). */
    if (head == SYM_Pattern && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        !is_plain_blank(e->data.function.args[1])) return 2200;
    if (head == SYM_And) return 2800;
    if (head == SYM_Or) return 2700;   /* looser than And, so Or nested in And parenthesises */
    if (head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_Inequality) return 3200;
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
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) return 5000;
            int64_t n, d;
            if (is_rational(exp, &n, &d) && n < 0) return 5000;
            if (exp->type == EXPR_REAL && exp->data.real < 0.0) return 5000;
        }
        return 6500;
    }
    if (head == SYM_Rational) return 5000;
    if (head == SYM_Complex) {
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

/* The display rewrite of e under its symbol's FormatValues: the first plain
 * (non-MakeBoxes) rule whose pattern matches, its right-hand side evaluated.
 * NULL when no rule applies. */
static Expr* format_value_rewrite(Expr* e) {
    const Expr* k = e;
    while (k->type == EXPR_FUNCTION) k = k->data.function.head;
    if (k->type != EXPR_SYMBOL) return NULL;
    for (Rule* r = symtab_get_format_values(k->data.symbol.name); r; r = r->next) {
        const Expr* p = r->pattern;
        if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
            p->data.function.head->data.symbol.name == SYM_MakeBoxes) continue;
        MatchEnv* env = env_new();
        Expr* out = NULL;
        if (match(e, r->pattern, env)) {
            Expr* bound = replace_bindings(r->replacement, env);
            out = evaluate(bound);
            expr_free(bound);
        }
        env_free(env);
        if (out) return out;
    }
    return NULL;
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

static void print_function_fullform(Expr* e) {
    expr_print_fullform(e->data.function.head);
    printf("[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i > 0) printf(", ");
        expr_print_fullform(e->data.function.args[i]);
    }    
    printf("]");
}

/* InputForm/FullForm text of a machine real: the shortest decimal that reads
 * back to the same double, so the printed form re-parses to the same value
 * (OutputForm's %g keeps six digits and does not).  The parser makes a
 * literal of 17+ significant digits arbitrary precision, so a 17-digit
 * mantissa carries Mathematica's machine-precision mark "`".  Exponents use
 * the input syntax *^. */
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
    printf("%s%s%s", mant, digits >= 17 ? "`" : "", expo);
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
            char buf[64];
            /* -0.0 prints as 0.0: it is SameQ to 0.0, and printing the sign
             * made exact-zero products like 0 (-1.2) look negative */
            snprintf(buf, sizeof(buf), "%g", e->data.real == 0.0 ? 0.0 : e->data.real);
            if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL && strchr(buf, 'E') == NULL) printf("%s.0", buf);
            else printf("%s", buf);
            break;
        }
        case EXPR_SYMBOL: printf("%s", context_display_name(e->data.symbol.name)); break;
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
        case EXPR_MPFR: {
            /* Print with decimal-digit precision matching the value's
             * stored binary precision. `%.*Rg` uses the `g` style (drops
             * trailing zeros) so most common values print cleanly; we
             * append ".0" for integer-valued outputs for consistency
             * with EXPR_REAL formatting (e.g. 3. → 3.0). */
            long digits = (long)ceil((double)mpfr_get_prec(e->data.mpfr)
                                     / 3.3219280948873626);
            if (digits < 1) digits = 1;
            char* buf = NULL;
            int len = mpfr_asprintf(&buf, "%.*Rg", (int)digits, e->data.mpfr);
            if (len >= 0 && buf) {
                printf("%s", buf);
                if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL
                    && strchr(buf, 'E') == NULL && strchr(buf, 'n') == NULL
                    && strchr(buf, 'i') == NULL /* nan / inf */) {
                    printf(".0");
                }
                mpfr_free_str(buf);
            } else {
                printf("<mpfr print error>");
            }
            break;
        }
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
static void print_definition_body(const char* name) {
    Rule* ov = symtab_get_own_values(name);
    Rule* dv = symtab_get_down_values(name);
    Rule* uv = symtab_get_up_values(name);
    int any = 0;
    for (Rule* r = ov; r; r = r->next) {
        if (any++) printf("\n");
        print_standard(r->pattern, 0);
        printf(r->delayed ? " := " : " = ");
        print_standard(r->replacement, 0);
    }
    for (Rule* r = uv; r; r = r->next) {
        if (any++) printf("\n");
        if (upvalue_renders_as_upset(name, r->pattern)) {
            print_standard(r->pattern, 0);
            printf(r->delayed ? " ^:= " : " ^= ");
        } else {
            printf("%s /: ", name);
            print_standard(r->pattern, 0);
            printf(r->delayed ? " := " : " = ");
        }
        print_standard(r->replacement, 0);
    }
    for (Rule* r = dv; r; r = r->next) {
        if (any++) printf("\n");
        print_standard(r->pattern, 0);
        printf(r->delayed ? " := " : " = ");
        print_standard(r->replacement, 0);
    }
    for (Rule* r = symtab_get_sub_values(name); r; r = r->next) {   /* f[x_][y_] := ... */
        if (any++) printf("\n");
        print_standard(r->pattern, 0);
        printf(r->delayed ? " := " : " = ");
        print_standard(r->replacement, 0);
    }
    if (!any) printf("Null");
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

static void print_standard(Expr* e, int parent_prec) {
    if (!e) { printf("Null"); return; }
    /* Format[...] display rules. A rewrite is itself displayed (and so
     * formatted again: Format[a[x_]] := b[x], Format[b[x_]] := 9 shows 9); a
     * rule that keeps reproducing a match stops at FORMAT_DEPTH_LIMIT, as the
     * reference stops at $RecursionLimit with Format::forml. */
    if (g_apply_formats && symtab_format_value_count && g_inputform_depth == 0 &&
        (e->type == EXPR_SYMBOL || e->type == EXPR_FUNCTION)) {
        if (g_format_depth >= FORMAT_DEPTH_LIMIT) {
            if (g_format_depth == FORMAT_DEPTH_LIMIT)
                mth_message("Format", "forml",
                            "Recursion depth limit $RecursionLimit exceeded while formatting output.");
            g_format_depth++;   /* report once per print */
        } else {
            Expr* shown = format_value_rewrite(e);
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
        
        if (head == SYM_FullForm && e->data.function.arg_count == 1) {
            /* FullForm keeps string quotes even inside Print / ToString. */
            g_inputform_depth++;
            expr_print_fullform(e->data.function.args[0]);
            g_inputform_depth--;
        }
        else if (head == SYM_TeXForm && e->data.function.arg_count == 1) {
            print_tex(e->data.function.args[0], 0);
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
            print_standard(e->data.function.args[0], parent_prec);
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
            const char* sep = NULL;
            if (e->data.function.arg_count >= 2
                && e->data.function.args[1]->type == EXPR_STRING)
                sep = e->data.function.args[1]->data.string;
            bool saved_of = g_print_output_form;
            g_print_output_form = true;
            for (size_t i = 0; i < list->data.function.arg_count; i++) {
                if (i > 0 && sep) fputs(sep, MTH_POUT);
                print_standard(list->data.function.args[i], 0);
            }
            g_print_output_form = saved_of;
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
        else if (head == SYM_Rational && e->data.function.arg_count == 2) {
            print_standard(e->data.function.args[0], 5000);
            printf("/");
            print_standard(e->data.function.args[1], 5000);
        }
        else if (head == SYM_Complex && e->data.function.arg_count == 2) {
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
                    im_abs = make_rational(-rn, rd);
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
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                is_negative = true;
                pos_exp = print_neg_int(exp->data.integer);
            } else if (is_rational(exp, &n, &d) && n < 0) {
                is_negative = true;
                pos_exp = make_rational(-n, d);
            } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                is_negative = true;
                pos_exp = expr_new_real(-exp->data.real);
            }

            if (is_negative) {
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
                    if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                        is_neg = true;
                        pos_exp = print_neg_int(exp->data.integer);
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        is_neg = true;
                        pos_exp = make_rational(-n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        is_neg = true;
                        pos_exp = expr_new_real(-exp->data.real);
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
        else if ((head == SYM_Equal || head == SYM_Unequal || head == SYM_Less || head == SYM_Greater || head == SYM_LessEqual || head == SYM_GreaterEqual || head == SYM_SameQ || head == SYM_UnsameQ || head == SYM_Set || head == SYM_SetDelayed || head == SYM_Rule || head == SYM_RuleDelayed || head == SYM_Condition || head == SYM_And || head == SYM_Or || head == SYM_Alternatives || head == SYM_NonCommutativeMultiply || head == SYM_Dot || head == SYM_StringExpression) && e->data.function.arg_count >= 2) {
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

            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (i > 0) printf("%s", op);
                print_standard(e->data.function.args[i], my_prec);
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
                                Expr* new_times = expr_new_function(expr_new_symbol(SYM_Times), new_args, t_copy->data.function.arg_count - 1);
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
                                Expr* new_times = expr_new_function(expr_new_symbol(SYM_Times), new_args, t_copy->data.function.arg_count - 1);
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
                } else {
                    if (is_negative) {
                        printf("-");
                    }
                }

                print_standard(to_print, 3500);
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
        else if (head == SYM_InterpolatingFunction && e->data.function.arg_count >= 2) {
            /* Only the domain element is printed explicitly; the remaining
             * data elements are abbreviated as <>, matching Mathematica's
             * standard output format. FullForm (handled separately) still
             * reveals the full structure. */
            printf("InterpolatingFunction[");
            print_standard(e->data.function.args[0], 0);
            printf(", <>]");
        }
        else if ((head == SYM_PredictorFunction ||
                  head == SYM_DimensionReducerFunction ||
                  head == SYM_LearnedDistribution ||
                  head == SYM_ClassifierFunction) &&
                 e->data.function.arg_count >= 2 &&
                 e->data.function.args[0]->type == EXPR_STRING) {
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
 * str_capture_open() returns a FILE* to push as the sink (NULL on failure);
 * str_capture_finish() closes it and returns a freshly malloc'd, NUL-terminated
 * string the caller frees (NULL on failure).
 * ------------------------------------------------------------------------- */
typedef struct {
    FILE*  f;
#ifndef _WIN32
    char*  buf;
    size_t len;
#endif
} StrCapture;

static FILE* str_capture_open(StrCapture* c) {
#ifdef _WIN32
    c->f = tmpfile();
#else
    c->buf = NULL;
    c->len = 0;
    c->f = open_memstream(&c->buf, &c->len);
#endif
    return c->f;
}

static char* str_capture_finish(StrCapture* c) {
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
    StrCapture c;
    if (!str_capture_open(&c)) return NULL;

    FILE* prev = mth_out_push(c.f);
    print_standard(e, 0);
    mth_out_pop(prev);

    return str_capture_finish(&c);
}

/* OutputForm rendering into a heap string (caller frees): as expr_to_string,
 * except that strings render as their characters, recursively, exactly as
 * Print shows them -- ToString[{"a", 1}] is "{a, 1}". Linear syntax is kept
 * (Mathilda has no 2-D OutputForm), so ToString[1/2] is "1/2". */
char* expr_to_string_outputform(Expr* e) {
    StrCapture c;
    if (!str_capture_open(&c)) return NULL;

    FILE* prev = mth_out_push(c.f);
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

    return str_capture_finish(&c);
}

char* expr_to_string_fullform(Expr* e) {
    StrCapture c;
    if (!str_capture_open(&c)) return NULL;

    FILE* prev = mth_out_push(c.f);
    expr_print_fullform(e);
    mth_out_pop(prev);

    return str_capture_finish(&c);
}

/* Render `e` (the result of a NumberForm NumberFormat function) to a heap
 * string in OutputForm: strings print raw, Row[...] concatenates. The active
 * NumberForm context is suspended so already-formatted mantissa/exponent
 * strings are not re-processed. Caller frees. */
char* numberform_format_result_to_string(const Expr* e) {
    StrCapture c;
    if (!str_capture_open(&c)) return NULL;

    FILE* prev = mth_out_push(c.f);

    const NumberFormCtx* saved_ctx = g_numberform_ctx;
    bool saved_of = g_print_output_form;
    g_numberform_ctx = NULL;
    g_print_output_form = true;

    print_standard((Expr*)e, 0);

    g_print_output_form = saved_of;
    g_numberform_ctx = saved_ctx;
    mth_out_pop(prev);

    return str_capture_finish(&c);
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
        *abs_out = make_rational(-rn, rd); return true;
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
            negf = make_rational(-rn, rd);
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
            char buf[64];
            snprintf(buf, sizeof(buf), "%g", e->data.real == 0.0 ? 0.0 : e->data.real);
            if (strchr(buf, '.') == NULL && strchr(buf, 'e') == NULL && strchr(buf, 'E') == NULL) printf("%s.", buf);
            else printf("%s", buf);
            break;
        }
        case EXPR_BIGINT: {
            char* s = mpz_get_str(NULL, 10, e->data.bigint);
            printf("%s", s);
            free(s);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            long digits = (long)ceil((double)mpfr_get_prec(e->data.mpfr) / 3.3219280948873626);
            if (digits < 1) digits = 1;
            char* buf = NULL;
            int len = mpfr_asprintf(&buf, "%.*Rg", (int)digits, e->data.mpfr);
            if (len >= 0 && buf) { printf("%s", buf); mpfr_free_str(buf); }
            break;
        }
#endif
        default: expr_print_fullform(e); break;
    }
}

static void print_tex_symbol(const char* raw_name) {
    const char* name = context_display_name(raw_name);
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
    else printf("\\text{%s}", name);
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
    if (strcmp(h, "Log") == 0)      return "\\log ";
    if (strcmp(h, "Exp") == 0)      return "\\exp ";
    if (strcmp(h, "Sqrt") == 0)     return NULL; /* handled via \sqrt{} */
    if (strcmp(h, "Abs") == 0)      return NULL; /* handled via |...|  */
    if (strcmp(h, "Floor") == 0)    return NULL; /* handled via \lfloor */
    if (strcmp(h, "Ceiling") == 0)  return NULL; /* handled via \lceil  */
    if (strcmp(h, "Gamma") == 0)    return "\\Gamma ";
    if (strcmp(h, "Re") == 0)       return "\\Re ";
    if (strcmp(h, "Im") == 0)       return "\\Im ";
    if (strcmp(h, "Conjugate") == 0)return NULL; /* handled via \overline */
    return NULL;
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
        *abs_out = make_rational(-rn, rd);
        return true;
    }
    return false;
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
    if (e->type == EXPR_STRING) { printf("\\text{\"%s\"}", e->data.string); return; }
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
            if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                exp_neg = true; pos_exp = print_neg_int(exp->data.integer);
            } else if (is_rational(exp, &n, &d) && n < 0) {
                exp_neg = true; pos_exp = make_rational(-n, d);
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
                        flipped_head = make_rational(-rn, rd);
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
                    if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
                        neg = true; pos_exp = print_neg_int(exp->data.integer);
                    } else if (is_rational(exp, &n, &d) && n < 0) {
                        neg = true; pos_exp = make_rational(-n, d);
                    } else if (exp->type == EXPR_REAL && exp->data.real < 0.0) {
                        neg = true; pos_exp = expr_new_real(-exp->data.real);
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
                    if (i > 0) printf(" ");
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
                        if (i > 0) printf(" ");
                        print_tex(num_args[i], nprec);
                    }
                }
                printf("}{");
                {
                    int dprec = (den_count > 1) ? 4500 : 0;
                    for (size_t i = 0; i < den_count; i++) {
                        if (i > 0) printf(" ");
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
                            Expr* new_rat = make_rational(-rn, rd);
                            expr_free(owned->data.function.args[0]);
                            owned->data.function.args[0] = new_rat;
                            to_print = owned;
                        }
                    }
                } else {
                    int64_t rn, rd;
                    if (is_rational(arg, &rn, &rd) && rn < 0) {
                        is_neg = true; owned = make_rational(-rn, rd); to_print = owned;
                    }
                }

                if (i > 0) printf(is_neg ? "-" : "+");
                else if (is_neg) printf("-");

                print_tex(to_print, 3500);
                if (owned) expr_free(owned);
            }
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
    } else {
        /* Non-symbol head: render as (head)(args...). */
        print_tex(head_expr, 9500);
        print_tex_args_parens(e, 0);
    }

    if (need_parens) printf("\\right)");
}

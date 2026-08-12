/* Mathilda — NumberForm and a minimal Row.  See numberform.h for the design.
 *
 * NumberForm is a PRINT WRAPPER: builtin_numberform is inert (returns NULL) so
 * the head survives in the tree, and all display work happens here, driven by
 * print.c which installs an active NumberFormCtx and routes every numeric leaf
 * through numberform_render_number.
 *
 * The per-number pipeline (nf_format_parts):
 *   1. reject non-finite; special-case zero and exact integers.
 *   2. extract the value's `count` significant base-10 digits + decimal
 *      exponent (mpfr_get_str for MPFR, "%.*e" for machine reals).
 *   3. decide scientific vs decimal (ScientificNotationThreshold, or the
 *      caller's ExponentFunction), and the displayed exponent (ExponentStep).
 *   4. place the digits into integer/fractional strings; {n,f} re-rounds to f
 *      fractional digits, plain-n drops trailing zeros.
 *   5. apply DigitBlock grouping, then assemble sign / point / multiplier /
 *      NumberFormat.
 * The same nf_format_parts drives the measure pass, so alignment widths and the
 * printed output can never disagree.
 */
#include "numberform.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "ndarray.h"     /* is_packed_list / ndarray_to_nested_list */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#define NF_LOG2_10 3.3219280948873626

/* Rendering the result of a user NumberFormat function to an OutputForm string
 * (raw strings, Row concatenation, active NumberForm context suspended). Lives
 * in print.c, which owns print_standard and the output-form flag. */
extern char* numberform_format_result_to_string(const Expr* e);

/* reqsigz fires at most once per NumberForm print (reset in build_ctx). */
static bool g_nf_reqsigz_warned = false;

/* ------------------------------------------------------------------ */
/* small growable string buffer                                        */
/* ------------------------------------------------------------------ */
typedef struct { char* p; size_t len, cap; } SB;

static void sb_init(SB* b) { b->cap = 32; b->len = 0; b->p = malloc(b->cap); b->p[0] = 0; }
static void sb_ensure(SB* b, size_t extra) {
    if (b->len + extra + 1 > b->cap) {
        while (b->len + extra + 1 > b->cap) b->cap *= 2;
        b->p = realloc(b->p, b->cap);
    }
}
static void sb_putc(SB* b, char c) { sb_ensure(b, 1); b->p[b->len++] = c; b->p[b->len] = 0; }
static void sb_puts(SB* b, const char* s) {
    if (!s) return;
    size_t n = strlen(s); sb_ensure(b, n);
    memcpy(b->p + b->len, s, n); b->len += n; b->p[b->len] = 0;
}
static void sb_prepend_char(SB* b, char c) {
    sb_ensure(b, 1);
    memmove(b->p + 1, b->p, b->len + 1);
    b->p[0] = c; b->len++;
}
static void sb_free(SB* b) { free(b->p); b->p = NULL; }

static char* nf_strdup(const char* s) { return mathilda_strdup(s ? s : ""); }
static char* nf_ltoa(long v) { char t[32]; snprintf(t, sizeof t, "%ld", v); return nf_strdup(t); }

/* floor division for possibly-negative a, b >= 1 (ExponentStep). */
static long nf_floordiv(long a, long b) {
    if (b <= 0) return a;
    long q = a / b, r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) q--;
    return q;
}

/* ------------------------------------------------------------------ */
/* digit grouping (DigitBlock)                                          */
/* ------------------------------------------------------------------ */
/* Integer digits grouped from the RIGHT: "1000000000" / 3 -> "1,000,000,000". */
static char* nf_group_int(const char* d, long block, const char* sep) {
    long n = (long)strlen(d);
    if (block <= 0 || n <= block) return nf_strdup(d);
    SB b; sb_init(&b);
    for (long i = 0; i < n; i++) {
        if (i > 0 && ((n - i) % block) == 0) sb_puts(&b, sep);
        sb_putc(&b, d[i]);
    }
    return b.p;
}
/* Fractional digits grouped from the LEFT. */
static char* nf_group_frac(const char* d, long block, const char* sep) {
    long n = (long)strlen(d);
    if (block <= 0 || n <= block) return nf_strdup(d);
    SB b; sb_init(&b);
    for (long i = 0; i < n; i++) {
        if (i > 0 && (i % block) == 0) sb_puts(&b, sep);
        sb_putc(&b, d[i]);
    }
    return b.p;
}

/* ------------------------------------------------------------------ */
/* significant-digit extraction                                        */
/* ------------------------------------------------------------------ */
/* Fill *outdigits (malloc'd, free with free()) with the `count` most-significant
 * base-10 digits of |value(e)|, *exp10 with floor(log10|value|) after rounding,
 * *ndig with the digit count. `e` must be a nonzero finite Real or MPFR. */
static bool nf_digits_of(const Expr* e, long count, char** outdigits,
                         long* exp10, long* ndig) {
    if (count < 1) count = 1;
    if (e->type == EXPR_REAL) {
        double v = e->data.real; if (v < 0) v = -v;
        size_t bufsz = (size_t)count + 64;
        char* buf = malloc(bufsz);
        snprintf(buf, bufsz, "%.*e", (int)(count - 1), v);
        char* s = malloc((size_t)count + 1);
        long k = 0; const char* p = buf;
        while (*p && *p != 'e' && *p != 'E') { if (*p != '.' && k < count) s[k++] = *p; p++; }
        long ex = 0;
        if (*p) { p++; ex = strtol(p, NULL, 10); }
        s[k] = 0;
        free(buf);
        *outdigits = s; *ndig = k; *exp10 = ex;
        return true;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        mpfr_t x; mpfr_init2(x, mpfr_get_prec(e->data.mpfr));
        mpfr_abs(x, e->data.mpfr, MPFR_RNDN);
        mpfr_exp_t decexp = 0;
        char* ms = mpfr_get_str(NULL, &decexp, 10, (size_t)count, x, MPFR_RNDN);
        mpfr_clear(x);
        if (!ms) return false;
        const char* p = ms; if (*p == '-') p++;   /* defensive: abs already taken */
        char* s = malloc((size_t)count + 1);
        long k = 0; while (*p && k < count) s[k++] = *p++;
        s[k] = 0;
        mpfr_free_str(ms);
        *outdigits = s; *ndig = k; *exp10 = (long)decexp - 1;
        return true;
    }
#endif
    return false;
}

static long nf_capacity(const Expr* e) {
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        long c = (long)((double)mpfr_get_prec(e->data.mpfr) / NF_LOG2_10);
        return c < 1 ? 1 : c;
    }
#endif
    (void)e;
    return 17;   /* a machine double round-trips in <= 17 significant digits */
}

/* ------------------------------------------------------------------ */
/* placement + rounding                                                */
/* ------------------------------------------------------------------ */
/* Split `s` (ndig significant digits, most-significant first) into integer and
 * fractional digit strings, where s[0] has place value 10^eff_exp10.  When the
 * integer positions exceed the available digits (a requested precision lower
 * than the integer-digit count) the integer part is padded with trailing zeros. */
static void nf_place(const char* s, long ndig, long eff_exp10, SB* isb, SB* fsb) {
    if (eff_exp10 >= 0) {
        long intcount = eff_exp10 + 1;
        for (long i = 0; i < intcount; i++) sb_putc(isb, (i < ndig) ? s[i] : '0');
        for (long i = intcount; i < ndig; i++) sb_putc(fsb, s[i]);
    } else {
        sb_putc(isb, '0');
        for (long i = 0; i < (-eff_exp10 - 1); i++) sb_putc(fsb, '0');
        for (long i = 0; i < ndig; i++) sb_putc(fsb, s[i]);
    }
}

/* Force the fractional part to exactly f digits: pad with zeros, or round (half
 * up) with carry propagating through the integer part. */
static void nf_round_frac(SB* isb, SB* fsb, long f) {
    long fl = (long)fsb->len;
    if (fl <= f) { for (long i = fl; i < f; i++) sb_putc(fsb, '0'); return; }
    int roundup = (fsb->p[f] >= '5');
    fsb->len = (size_t)f; fsb->p[f] = 0;
    if (!roundup) return;
    long i = f - 1; int carry = 1;
    while (carry && i >= 0) {
        if (fsb->p[i] == '9') fsb->p[i] = '0';
        else { fsb->p[i]++; carry = 0; }
        i--;
    }
    if (carry) {
        long j = (long)isb->len - 1;
        while (carry && j >= 0) {
            if (isb->p[j] == '9') isb->p[j] = '0';
            else { isb->p[j]++; carry = 0; }
            j--;
        }
        if (carry) sb_prepend_char(isb, '1');
    }
}

/* ------------------------------------------------------------------ */
/* user-function options                                               */
/* ------------------------------------------------------------------ */
/* Evaluate ExponentFunction[exp10].  On success sets *isnull (Null return ->
 * force decimal) or *out (integer return -> displayed exponent). */
static bool nf_apply_exp_func(const NumberFormCtx* ctx, long exp10,
                              long* out, bool* isnull) {
    Expr* a = expr_new_integer(exp10);
    Expr* call = expr_new_function(expr_copy((Expr*)ctx->exp_func), &a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    bool ok = false; *isnull = false;
    if (r) {
        if (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_Null) { *isnull = true; ok = true; }
        else if (r->type == EXPR_INTEGER) { *out = r->data.integer; ok = true; }
    }
    expr_free(r);
    return ok;
}

/* Evaluate NumberFormat[mantissa, "10", exponent] and render to a string. */
static char* nf_call_number_format(const NumberFormCtx* ctx, const char* mant,
                                   const char* base, const char* expstr) {
    Expr* margs[3];
    margs[0] = expr_new_string(mant);
    margs[1] = expr_new_string(base);
    margs[2] = expr_new_string(expstr);
    Expr* call = expr_new_function(expr_copy((Expr*)ctx->num_format), margs, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    char* s = numberform_format_result_to_string(r);
    expr_free(r);
    return s;
}

/* ------------------------------------------------------------------ */
/* per-number formatting                                               */
/* ------------------------------------------------------------------ */
typedef struct {
    bool ok;
    const char* sign;   /* borrowed (ctx sign strings) */
    char* intstr;       /* owned */
    char* fracstr;      /* owned */
    bool has_point;
    char* sci;          /* owned or NULL */
    int signn, intdign, fracdign, scin;
} NFParts;

static void nf_parts_free(NFParts* p) {
    free(p->intstr); free(p->fracstr); free(p->sci);
    p->intstr = p->fracstr = p->sci = NULL;
}

static bool nf_is_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool nf_is_finite(const Expr* e) {
    if (e->type == EXPR_REAL) return isfinite(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_number_p(e->data.mpfr) != 0;
#endif
    return true;   /* integers/bigints always finite */
}

static bool nf_is_zero(const Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) == 0;
    return false;
}

static bool nf_is_negative(const Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real < 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_signbit(e->data.mpfr) && !mpfr_zero_p(e->data.mpfr);
#endif
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) < 0;
    return false;
}

/* Decompose `e` into printable parts under `ctx`. Returns false (parts left
 * clean) for a value the engine declines. */
static bool nf_format_parts(const Expr* e, const NumberFormCtx* ctx, NFParts* out) {
    memset(out, 0, sizeof *out);
    out->sign = "";

    if (!nf_is_number(e) || !nf_is_finite(e)) return false;

    bool neg = nf_is_negative(e);

    /* exact integers: shown in full, no rounding, no decimal point. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        mpz_t z; mpz_init(z); expr_to_mpz(e, z); mpz_abs(z, z);
        char* digits = mpz_get_str(NULL, 10, z);
        mpz_clear(z);
        out->intstr = nf_group_int(digits, ctx->block_left, ctx->sep_int);
        free(digits);
        out->fracstr = nf_strdup("");
        out->has_point = false;
        out->sign = neg ? ctx->sign_neg : ctx->sign_pos;
        out->signn = (int)strlen(out->sign);
        out->intdign = (int)strlen(out->intstr);
        out->ok = true;
        return true;
    }

    /* zero: 0. (or 0.000... under {n,f}), never scientific, never signed. */
    if (nf_is_zero(e)) {
        out->intstr = nf_strdup("0");
        if (ctx->have_f && ctx->f > 0) {
            SB f; sb_init(&f);
            for (long i = 0; i < ctx->f; i++) sb_putc(&f, '0');
            out->fracstr = f.p;
        } else {
            out->fracstr = nf_strdup("");
        }
        out->has_point = true;
        out->sign = ctx->sign_pos;
        out->signn = (int)strlen(out->sign);
        out->intdign = 1;
        out->fracdign = (int)strlen(out->fracstr);
        out->ok = true;
        return true;
    }

    long prec = (ctx->n >= 0) ? ctx->n : ctx->default_print_precision;
    if (prec < 1) prec = 1;
    bool explicit_n = (ctx->n >= 0);
    long capacity = nf_capacity(e);

    /* probe magnitude */
    char* pd = NULL; long exp10 = 0, pnd = 0;
    if (!nf_digits_of(e, capacity, &pd, &exp10, &pnd)) return false;
    free(pd);
    long intdig = (exp10 >= 0) ? (exp10 + 1) : 1;

    /* scientific decision + displayed exponent */
    bool sci; long displayed_exp = 0;
    if (ctx->exp_func) {
        long r = 0; bool isnull = false;
        if (nf_apply_exp_func(ctx, exp10, &r, &isnull)) {
            if (isnull) sci = false;
            else { sci = true; displayed_exp = r; }
        } else {
            sci = (exp10 < ctx->sci_lo || exp10 >= ctx->sci_hi);
            if (sci) displayed_exp = nf_floordiv(exp10, ctx->exp_step) * ctx->exp_step;
        }
    } else {
        sci = (exp10 < ctx->sci_lo || exp10 >= ctx->sci_hi);
        if (sci) displayed_exp = nf_floordiv(exp10, ctx->exp_step) * ctx->exp_step;
    }
    if (ctx->have_f) sci = false;   /* {n,f} is fixed-point */

    /* how many significant digits to render */
    long count;
    if (sci) count = prec;
    else if (explicit_n) count = prec;
    else count = (exp10 >= 0 && intdig > prec) ? intdig : prec;   /* default keeps all int digits */
    if (count > capacity) count = capacity;
    if (count < 1) count = 1;

    char* s = NULL; long ndig = 0, exp10b = 0;
    if (!nf_digits_of(e, count, &s, &exp10b, &ndig)) return false;

    /* rounding at `count` may bump the magnitude (9.99 -> 10); re-derive the
     * step-rounded exponent from the post-rounding magnitude. */
    if (sci && !ctx->exp_func)
        displayed_exp = nf_floordiv(exp10b, ctx->exp_step) * ctx->exp_step;
    long eff_exp10 = sci ? (exp10b - displayed_exp) : exp10b;

    if (!sci && eff_exp10 >= 0 && (eff_exp10 + 1) > ndig && !g_nf_reqsigz_warned) {
        g_nf_reqsigz_warned = true;
        fprintf(stderr, "NumberForm::reqsigz: Requested number precision is lower "
                        "than number of digits shown; padding with zeros.\n");
    }

    SB isb, fsb; sb_init(&isb); sb_init(&fsb);
    nf_place(s, ndig, eff_exp10, &isb, &fsb);
    free(s);

    if (ctx->have_f) {
        nf_round_frac(&isb, &fsb, ctx->f);
    } else {
        while (fsb.len > 0 && fsb.p[fsb.len - 1] == '0') { fsb.len--; fsb.p[fsb.len] = 0; }
    }

    char* intg = nf_group_int(isb.p, ctx->block_left, ctx->sep_int);
    char* fracg = nf_group_frac(fsb.p, ctx->block_right, ctx->sep_frac);
    sb_free(&isb); sb_free(&fsb);

    /* NumberFormat override: hand (mantissa, "10", exponent) to the user
     * function and print whatever it returns. */
    if (ctx->num_format) {
        SB mant; sb_init(&mant);
        sb_puts(&mant, intg); sb_puts(&mant, ctx->point); sb_puts(&mant, fracg);
        char* expstr = sci ? nf_ltoa(displayed_exp) : nf_strdup("");
        char* fmt = nf_call_number_format(ctx, mant.p, "10", expstr);
        const char* sgn = neg ? ctx->sign_neg : ctx->sign_pos;
        SB whole; sb_init(&whole);
        sb_puts(&whole, sgn); if (fmt) sb_puts(&whole, fmt);
        out->intstr = whole.p;
        out->fracstr = nf_strdup("");
        out->has_point = false;
        out->sign = "";
        out->signn = 0;
        out->intdign = (int)strlen(out->intstr);
        out->ok = true;
        free(fmt); free(expstr); sb_free(&mant);
        free(intg); free(fracg);
        return true;
    }

    char* sci_suffix = NULL;
    if (sci) {
        SB sfx; sb_init(&sfx);
        sb_puts(&sfx, ctx->mult);
        sb_puts(&sfx, "10^(");
        char eb[32]; snprintf(eb, sizeof eb, "%ld", displayed_exp); sb_puts(&sfx, eb);
        sb_puts(&sfx, ")");
        sci_suffix = sfx.p;
    }

    out->sign = neg ? ctx->sign_neg : ctx->sign_pos;
    out->intstr = intg;
    out->fracstr = fracg;
    out->has_point = true;
    out->sci = sci_suffix;
    out->signn = (int)strlen(out->sign);
    out->intdign = (int)strlen(intg);
    out->fracdign = (int)strlen(fracg);
    out->scin = sci_suffix ? (int)strlen(sci_suffix) : 0;
    out->ok = true;
    return true;
}

/* ------------------------------------------------------------------ */
/* measure pass + field width (padding)                                */
/* ------------------------------------------------------------------ */
static void nf_measure(const Expr* e, NumberFormCtx* ctx) {
    if (!e) return;
    if (is_packed_list(e)) {
        Expr* nested = ndarray_to_nested_list(e);
        nf_measure(nested, ctx);
        expr_free(nested);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nf_measure(e->data.function.args[i], ctx);
        return;
    }
    if (nf_is_number(e)) {
        NFParts p;
        if (nf_format_parts(e, ctx, &p) && p.ok) {
            if (p.signn > ctx->m_maxsignn)  ctx->m_maxsignn = p.signn;
            if (p.intdign > ctx->m_maxintdig) ctx->m_maxintdig = p.intdign;
            if (p.fracdign > ctx->m_maxfrac) ctx->m_maxfrac = p.fracdign;
            if (p.scin > ctx->m_maxscin)    ctx->m_maxscin = p.scin;
            if (p.has_point) ctx->m_anypoint = true;
            nf_parts_free(&p);
        }
    }
}

static void nf_compute_field(NumberFormCtx* ctx) {
    long prec = (ctx->n >= 0) ? ctx->n : ctx->default_print_precision;
    if (prec < 1) prec = 1;
    int frac_col = ctx->have_f ? (int)ctx->f : ctx->m_maxfrac;
    int int_col = ctx->m_maxintdig;
    int reserve = (int)prec - frac_col;
    if (reserve > int_col) int_col = reserve;
    ctx->field_sign = ctx->m_maxsignn;
    ctx->field_int = int_col;
    ctx->field_point = ctx->m_anypoint ? (int)strlen(ctx->point) : 0;
    ctx->field_frac = frac_col;
}

/* ------------------------------------------------------------------ */
/* emit                                                                */
/* ------------------------------------------------------------------ */
static void nf_emit_plain(const NFParts* p, const NumberFormCtx* ctx) {
    fputs(p->sign, stdout);
    fputs(p->intstr, stdout);
    if (p->has_point) fputs(ctx->point, stdout);
    fputs(p->fracstr, stdout);
    if (p->sci) fputs(p->sci, stdout);
}

static void nf_emit_padded(const NFParts* p, const NumberFormCtx* ctx) {
    int point_this = p->has_point ? (int)strlen(ctx->point) : 0;
    int field = ctx->field_sign + ctx->field_int + ctx->field_point
              + ctx->field_frac + ctx->m_maxscin;
    bool rpad_on = (ctx->pad_right && ctx->pad_right[0] != '\0');
    int efffrac = rpad_on ? ctx->field_frac : p->fracdign;
    int actual = p->signn + p->intdign + point_this + efffrac + p->scin;
    int left_pad = field - actual;
    if (left_pad < 0) left_pad = 0;

    if (ctx->sign_padding) {
        fputs(p->sign, stdout);
        for (int i = 0; i < left_pad; i++) fputs(ctx->pad_left, stdout);
    } else {
        for (int i = 0; i < left_pad; i++) fputs(ctx->pad_left, stdout);
        fputs(p->sign, stdout);
    }
    fputs(p->intstr, stdout);
    if (p->has_point) fputs(ctx->point, stdout);
    fputs(p->fracstr, stdout);
    if (rpad_on)
        for (int i = 0; i < ctx->field_frac - p->fracdign; i++) fputs(ctx->pad_right, stdout);
    if (p->sci) fputs(p->sci, stdout);
}

bool numberform_render_number(const Expr* e, const NumberFormCtx* ctx) {
    NFParts p;
    if (!nf_format_parts(e, ctx, &p) || !p.ok) { nf_parts_free(&p); return false; }
    if (ctx->has_padding) nf_emit_padded(&p, ctx);
    else                  nf_emit_plain(&p, ctx);
    nf_parts_free(&p);
    return true;
}

/* ------------------------------------------------------------------ */
/* option parsing                                                      */
/* ------------------------------------------------------------------ */
static const char* nf_as_str(const Expr* e, const char* fallback) {
    return (e && e->type == EXPR_STRING) ? e->data.string : fallback;
}
static bool nf_is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
            || e->data.function.head->data.symbol.name == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}
static bool nf_list2(const Expr* e, const Expr** a, const Expr** b) {
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count == 2) {
        *a = e->data.function.args[0]; *b = e->data.function.args[1];
        return true;
    }
    return false;
}
static bool nf_sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static void nf_apply_option_rule(NumberFormCtx* ctx, const Expr* rule) {
    if (!nf_is_rule(rule)) return;
    const Expr* key = rule->data.function.args[0];
    const Expr* val = rule->data.function.args[1];
    if (key->type != EXPR_SYMBOL) return;
    const char* k = key->data.symbol.name;
    const Expr *a, *b;

    if (k == SYM_NumberSigns) {
        if (nf_list2(val, &a, &b)) { ctx->sign_neg = nf_as_str(a, ctx->sign_neg);
                                     ctx->sign_pos = nf_as_str(b, ctx->sign_pos); }
    } else if (k == SYM_NumberPoint) {
        ctx->point = nf_as_str(val, ctx->point);
    } else if (k == SYM_NumberMultiplier) {
        ctx->mult = nf_as_str(val, ctx->mult);
    } else if (k == SYM_DigitBlock) {
        if (val->type == EXPR_INTEGER) {
            long v = val->data.integer < 0 ? 0 : val->data.integer;
            ctx->block_left = ctx->block_right = v;
        } else if (nf_sym_is(val, "Infinity")) {
            ctx->block_left = ctx->block_right = 0;
        } else if (nf_list2(val, &a, &b)) {
            if (a->type == EXPR_INTEGER) ctx->block_left = a->data.integer < 0 ? 0 : a->data.integer;
            if (b->type == EXPR_INTEGER) ctx->block_right = b->data.integer < 0 ? 0 : b->data.integer;
        }
    } else if (k == SYM_NumberSeparator) {
        if (val->type == EXPR_STRING) { ctx->sep_int = val->data.string; ctx->sep_frac = val->data.string; }
        else if (nf_list2(val, &a, &b)) { ctx->sep_int = nf_as_str(a, ctx->sep_int);
                                          ctx->sep_frac = nf_as_str(b, ctx->sep_frac); }
    } else if (k == SYM_NumberPadding) {
        if (val->type == EXPR_STRING) { ctx->pad_left = val->data.string; ctx->pad_right = val->data.string; }
        else if (nf_list2(val, &a, &b)) { ctx->pad_left = nf_as_str(a, ctx->pad_left);
                                          ctx->pad_right = nf_as_str(b, ctx->pad_right); }
    } else if (k == SYM_SignPadding) {
        if (nf_sym_is(val, "True")) ctx->sign_padding = true;
        else if (nf_sym_is(val, "False")) ctx->sign_padding = false;
    } else if (k == SYM_ExponentStep) {
        if (val->type == EXPR_INTEGER && val->data.integer >= 1) ctx->exp_step = val->data.integer;
    } else if (k == SYM_ExponentFunction) {
        ctx->exp_func = nf_sym_is(val, "Automatic") ? NULL : val;
    } else if (k == SYM_NumberFormat) {
        ctx->num_format = nf_sym_is(val, "Automatic") ? NULL : val;
    } else if (k == SYM_ScientificNotationThreshold) {
        if (nf_list2(val, &a, &b)) {
            if (a->type == EXPR_INTEGER) ctx->sci_lo = a->data.integer;
            if (b->type == EXPR_INTEGER) ctx->sci_hi = b->data.integer;
        }
    } else if (k == SYM_DefaultPrintPrecision) {
        if (val->type == EXPR_INTEGER && val->data.integer >= 1)
            ctx->default_print_precision = val->data.integer;
    }
}

static void nf_parse_spec(NumberFormCtx* ctx, const Expr* a1) {
    const Expr *a, *b;
    if (a1->type == EXPR_INTEGER) {
        ctx->n = a1->data.integer;
    } else if (nf_list2(a1, &a, &b) && a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
        ctx->n = a->data.integer;
        ctx->f = b->data.integer;
        ctx->have_f = true;
    }
}

bool numberform_build_ctx(const Expr* res, NumberFormCtx* ctx) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return false;

    g_nf_reqsigz_warned = false;

    ctx->n = -1; ctx->f = -1; ctx->have_f = false;
    ctx->sign_neg = "-"; ctx->sign_pos = "";
    ctx->point = "."; ctx->mult = "*";
    ctx->block_left = 0; ctx->block_right = 0;
    ctx->sep_int = ","; ctx->sep_frac = " ";
    ctx->pad_left = ""; ctx->pad_right = "";
    ctx->sign_padding = false;
    ctx->exp_step = 1;
    ctx->exp_func = NULL; ctx->num_format = NULL;
    ctx->sci_lo = -5; ctx->sci_hi = 6;
    ctx->default_print_precision = 6;
    ctx->has_padding = false;
    ctx->m_maxsignn = ctx->m_maxintdig = ctx->m_maxfrac = ctx->m_maxscin = 0;
    ctx->m_anypoint = false;
    ctx->field_sign = ctx->field_int = ctx->field_point = ctx->field_frac = 0;

    /* registered defaults (honours SetOptions[NumberForm, ...]) */
    Expr* defs = symtab_get_options("NumberForm");
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++)
            nf_apply_option_rule(ctx, defs->data.function.args[i]);
    }

    /* positional spec (first non-Rule argument) then explicit options */
    size_t argc = res->data.function.arg_count;
    size_t opt_start = 1;
    if (argc >= 2 && !nf_is_rule(res->data.function.args[1])) {
        nf_parse_spec(ctx, res->data.function.args[1]);
        opt_start = 2;
    }
    for (size_t i = opt_start; i < argc; i++)
        nf_apply_option_rule(ctx, res->data.function.args[i]);

    ctx->has_padding = (ctx->pad_left[0] != '\0') || (ctx->pad_right[0] != '\0');
    if (ctx->has_padding) {
        nf_measure(res->data.function.args[0], ctx);
        nf_compute_field(ctx);
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* builtins + registration                                             */
/* ------------------------------------------------------------------ */
Expr* builtin_numberform(Expr* res) { (void)res; return NULL; }  /* inert wrapper */
Expr* builtin_row(Expr* res)        { (void)res; return NULL; }  /* inert wrapper */

void numberform_init(void) {
    symtab_add_builtin("NumberForm", builtin_numberform);
    SymbolDef* d = symtab_get_def("NumberForm");
    if (d) d->attributes |= ATTR_NHOLDREST | ATTR_PROTECTED;

    symtab_add_builtin("Row", builtin_row);
    SymbolDef* dr = symtab_get_def("Row");
    if (dr) dr->attributes |= ATTR_PROTECTED;
}

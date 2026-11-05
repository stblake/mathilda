/*
 * print_fortran.c -- FortranForm / CForm renderers (see print_fortran.h).
 *
 * Two stages:
 *   1. build(): Expr -> layout tree.  A token is a run of text that is never
 *      split at an operator; a row is a sequence of children.  The rows mirror
 *      the reference's box nesting -- a sum, a product and a power are each a
 *      row, a parenthesised group is a row of "(", the group, ")", and a
 *      function call is ONE flat row of head, "(", arguments and commas, ")".
 *      The nesting depth of a row is what sets the indent of a continuation
 *      line broken inside it.
 *   2. layout(): flat text when it fits, else fixed-form Fortran lines.
 *
 * Text conventions (checked against the reference, one probe per rule):
 *   - Rational: alone, its machine value (2/3 -> 0.6666666666666666); as a
 *     coefficient, numerator in the numerator and the denominator as a real
 *     (a/100., (-2*x)/3., 1/(2.*x)), except -1/q, which prints as its value
 *     times the rest (-0.5*x/y).
 *   - A power with a negative exponent goes to the denominator inside a
 *     product (y/x**2) but keeps its exponent alone (x**(-2)); x^-1 is 1/x,
 *     x^(1/2) is Sqrt(x), x^(-1/2) is 1/Sqrt(x); rational exponents print as
 *     reals (x**1.5).
 *   - Reals: shortest round-trip digits, a trailing point when integral
 *     (100.), scientific "1.e6" from 10^6 up and below 10^-5.
 *   - Fortran: (re,im) complexes, .eq./.ne./.lt./.le./.gt./.ge., .and./.or./
 *     .not.; C: Complex(re,im), Power(b,e), == != < <= > >= && || ! %, a[i].
 *   - Every other call is head(arg,arg) (List(1,2), Part(a,1) in Fortran).
 */
#include "print_fortran.h"

#include <inttypes.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "context.h"
#include "print.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ */
/* Layout tree                                                         */
/* ------------------------------------------------------------------ */

enum { BRK_NONE = 0, BRK_TERM = 1, BRK_ARG = 2 };

/* Precedences of the printed forms (higher binds tighter). */
enum {
    P_LOWEST = 0,
    P_OR = 210,
    P_AND = 215,
    P_NOT = 230,
    P_REL = 290,
    P_PLUS = 310,
    P_UNARY = 395,    /* -x, a negative number */
    P_TIMES = 400,    /* a*b and a/b */
    P_POWER = 590,
    P_ATOM = 1000
};

typedef struct FNode {
    char* text;               /* token text (owned); NULL for a row */
    size_t len;               /* flat length, cached */
    struct FNode** kids;      /* row children */
    unsigned char* brk;       /* per child: break allowed before it */
    size_t n, cap;
} FNode;

static FNode* fn_alloc(void) {
    FNode* f = calloc(1, sizeof(FNode));
    if (!f) abort();
    return f;
}

static FNode* ftok(const char* s) {
    FNode* f = fn_alloc();
    size_t n = strlen(s);
    f->text = malloc(n + 1);
    if (!f->text) abort();
    memcpy(f->text, s, n + 1);
    f->len = n;
    return f;
}

static FNode* ftok_own(char* s) {
    FNode* f = fn_alloc();
    f->text = s;
    f->len = strlen(s);
    return f;
}

static FNode* frow(void) { return fn_alloc(); }

static void fnode_add(FNode* r, FNode* k, unsigned char brk) {
    if (r->n == r->cap) {
        size_t nc = r->cap ? 2 * r->cap : 4;
        r->kids = realloc(r->kids, nc * sizeof(FNode*));
        r->brk = realloc(r->brk, nc);
        if (!r->kids || !r->brk) abort();
        r->cap = nc;
    }
    r->kids[r->n] = k;
    r->brk[r->n] = brk;
    r->n++;
    r->len += k->len;
}

static void fadd_tok(FNode* r, const char* s) { fnode_add(r, ftok(s), BRK_NONE); }

static void ffree(FNode* f) {
    if (!f) return;
    for (size_t i = 0; i < f->n; i++) ffree(f->kids[i]);
    free(f->kids);
    free(f->brk);
    free(f->text);
    free(f);
}

static FNode* fparen(FNode* k) {
    FNode* r = frow();
    fadd_tok(r, "(");
    fnode_add(r, k, BRK_NONE);
    fadd_tok(r, ")");
    return r;
}

/* ------------------------------------------------------------------ */
/* Output buffer                                                       */
/* ------------------------------------------------------------------ */

typedef struct { char* s; size_t n, cap; } SBuf;

static void sb_putn(SBuf* b, const char* s, size_t n) {
    if (b->n + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 64;
        while (nc < b->n + n + 1) nc *= 2;
        b->s = realloc(b->s, nc);
        if (!b->s) abort();
        b->cap = nc;
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = '\0';
}

static void sb_puts(SBuf* b, const char* s) { sb_putn(b, s, strlen(s)); }

static char* sb_take(SBuf* b) {
    if (!b->s) { char* e = malloc(1); if (!e) abort(); e[0] = '\0'; return e; }
    return b->s;
}

/* ------------------------------------------------------------------ */
/* Numbers                                                             */
/* ------------------------------------------------------------------ */

static bool head_is(const Expr* e, const char* sym, size_t argc) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == argc &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

static bool is_int(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

static bool is_rational(const Expr* e) {
    return head_is(e, SYM_Rational, 2) && is_int(e->data.function.args[0]) &&
           is_int(e->data.function.args[1]);
}

static bool is_real_atom(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool is_complex(const Expr* e) {
    return head_is(e, SYM_Complex, 2);
}

static bool is_number(const Expr* e) {
    return is_int(e) || is_real_atom(e) || is_rational(e) || is_complex(e);
}

static int int_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (e->data.integer > 0) - (e->data.integer < 0);
    return mpz_sgn(e->data.bigint);
}

static int real_sign(const Expr* e) {
    if (e->type == EXPR_REAL) return (e->data.real > 0) - (e->data.real < 0);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_sgn(e->data.mpfr);
#endif
    return 0;
}

/* Sign of a real number; a complex number's is its real part's, or its
 * imaginary part's when the real part is zero (a - I b prints "a - (0,1)*b"). */
static int num_sign(const Expr* e) {
    if (is_int(e)) return int_sign(e);
    if (is_real_atom(e)) return real_sign(e);
    if (is_rational(e)) return int_sign(e->data.function.args[0]);
    if (is_complex(e)) {
        int s = num_sign(e->data.function.args[0]);
        return s ? s : num_sign(e->data.function.args[1]);
    }
    return 0;
}

static bool is_int_value(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

static bool is_rational_value(const Expr* e, int64_t p, int64_t q) {
    return is_rational(e) && is_int_value(e->data.function.args[0], p) &&
           is_int_value(e->data.function.args[1], q);
}

/* -e for a number, as a fresh (unevaluated) node. */
static Expr* num_negate(const Expr* e) {
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer != INT64_MIN) return expr_new_integer(-e->data.integer);
        mpz_t z; mpz_init(z); expr_to_mpz(e, z); mpz_neg(z, z);
        Expr* r = expr_new_bigint_from_mpz(z); mpz_clear(z);
        return r;
    }
    if (e->type == EXPR_BIGINT) {
        mpz_t z; mpz_init(z); mpz_neg(z, e->data.bigint);
        Expr* r = expr_new_bigint_from_mpz(z); mpz_clear(z);
        return r;
    }
    if (e->type == EXPR_REAL) return expr_new_real(-e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        Expr* r = expr_new_mpfr_bits(mpfr_get_prec(e->data.mpfr));
        mpfr_neg(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
        return r;
    }
#endif
    if (is_rational(e)) {
        Expr* a[2] = { num_negate(e->data.function.args[0]),
                       expr_copy(e->data.function.args[1]) };
        return expr_new_function(expr_new_symbol(SYM_Rational), a, 2);
    }
    if (is_complex(e)) {
        Expr* a[2] = { is_number(e->data.function.args[0]) ? num_negate(e->data.function.args[0])
                                                            : expr_copy(e->data.function.args[0]),
                       is_number(e->data.function.args[1]) ? num_negate(e->data.function.args[1])
                                                            : expr_copy(e->data.function.args[1]) };
        return expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
    }
    return expr_copy((Expr*)e);
}

/* The decimal text of sign * 0.D * 10^(X+1): D the significant digits (no
 * trailing zeros), X the decimal exponent of its first digit.  Scientific
 * from 10^6 up and below 10^-5 ("1.e6", "1.5e-7"), else positional with a
 * trailing point when integral ("100000.", "0.000015"). */
static char* fmt_decimal(bool neg, const char* D, long X) {
    SBuf b = {0};
    size_t nd = strlen(D);
    if (neg) sb_puts(&b, "-");
    if (nd == 0 || (nd == 1 && D[0] == '0')) {
        sb_puts(&b, "0.");
        return sb_take(&b);
    }
    if (X >= 6 || X < -5) {
        sb_putn(&b, D, 1);
        sb_puts(&b, ".");
        if (nd > 1) sb_putn(&b, D + 1, nd - 1);
        char eb[32];
        snprintf(eb, sizeof eb, "e%ld", X);
        sb_puts(&b, eb);
    } else if (X >= 0) {
        size_t ip = (size_t)X + 1;
        if (nd >= ip) {
            sb_putn(&b, D, ip);
            sb_puts(&b, ".");
            if (nd > ip) sb_putn(&b, D + ip, nd - ip);
        } else {
            sb_putn(&b, D, nd);
            for (size_t i = nd; i < ip; i++) sb_puts(&b, "0");
            sb_puts(&b, ".");
        }
    } else {
        sb_puts(&b, "0.");
        for (long i = 0; i < -X - 1; i++) sb_puts(&b, "0");
        sb_putn(&b, D, nd);
    }
    return sb_take(&b);
}

/* Shortest round-trip text of a machine real (15..17 significant digits). */
static char* fmt_double(double d) {
    if (d == 0.0 || d != d) return fmt_decimal(false, "0", 0);
    char buf[64];
    int p = 15;
    for (; p <= 17; p++) {
        snprintf(buf, sizeof buf, "%.*e", p - 1, d);
        if (strtod(buf, NULL) == d) break;
    }
    if (p > 17) snprintf(buf, sizeof buf, "%.16e", d);
    bool neg = buf[0] == '-';
    const char* s = buf + (neg ? 1 : 0);
    char D[40];
    size_t nd = 0;
    const char* q = s;
    for (; *q && *q != 'e'; q++)
        if (*q >= '0' && *q <= '9' && nd < sizeof D - 1) D[nd++] = *q;
    D[nd] = '\0';
    while (nd > 1 && D[nd - 1] == '0') D[--nd] = '\0';
    long X = (*q == 'e') ? strtol(q + 1, NULL, 10) : 0;
    return fmt_decimal(neg, D, X);
}

#ifdef USE_MPFR
static char* fmt_mpfr(const Expr* e) {
    mpfr_exp_t ex = 0;
    char* s = mpfr_get_str(NULL, &ex, 10, 0, e->data.mpfr, MPFR_RNDN);
    if (!s) return fmt_decimal(false, "0", 0);
    bool neg = s[0] == '-';
    const char* D0 = s + (neg ? 1 : 0);
    size_t nd = strlen(D0);
    char* D = malloc(nd + 1);
    if (!D) abort();
    memcpy(D, D0, nd + 1);
    while (nd > 1 && D[nd - 1] == '0') D[--nd] = '\0';
    bool zero = mpfr_zero_p(e->data.mpfr);
    char* r = fmt_decimal(neg, zero ? "0" : D, (long)ex - 1);
    free(D);
    mpfr_free_str(s);
    return r;
}
#endif

static char* fmt_int(const Expr* e) {
    if (e->type == EXPR_INTEGER) {
        char buf[32];
        snprintf(buf, sizeof buf, "%" PRId64, e->data.integer);
        char* r = malloc(strlen(buf) + 1);
        if (!r) abort();
        strcpy(r, buf);
        return r;
    }
    size_t n = mpz_sizeinbase(e->data.bigint, 10) + 2;
    char* r = malloc(n);
    if (!r) abort();
    mpz_get_str(r, 10, e->data.bigint);
    return r;
}

static double int_to_double(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    return mpz_get_d(e->data.bigint);
}

static double rational_to_double(const Expr* e) {
    mpq_t q;
    mpq_init(q);
    expr_to_mpz(e->data.function.args[0], mpq_numref(q));
    expr_to_mpz(e->data.function.args[1], mpq_denref(q));
    mpq_canonicalize(q);
    double d = mpq_get_d(q);
    mpq_clear(q);
    return d;
}

/* ------------------------------------------------------------------ */
/* Builder                                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    FormLang lang;
    const char* form;        /* SYM_FortranForm / SYM_CForm (Format rules) */
    int fmt_depth;
} Ctx;

#define FMT_DEPTH_LIMIT 256

static FNode* build(Ctx* c, const Expr* e, int* prec);

static FNode* build_wrapped(Ctx* c, const Expr* e, int need) {
    int p;
    FNode* f = build(c, e, &p);
    return p < need ? fparen(f) : f;
}

static FNode* wrap_node(FNode* f, int p, int need) {
    return p < need ? fparen(f) : f;
}

/* A real-number part of a complex number, or a rational shown as a real. */
static FNode* build_real_part(Ctx* c, const Expr* e) {
    int p;
    if (is_rational(e)) return ftok_own(fmt_double(rational_to_double(e)));
    return build(c, e, &p);
}

static FNode* build_number(Ctx* c, const Expr* e, int* prec) {
    *prec = num_sign(e) < 0 ? P_UNARY : P_ATOM;
    if (is_int(e)) return ftok_own(fmt_int(e));
    if (e->type == EXPR_REAL) return ftok_own(fmt_double(e->data.real));
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return ftok_own(fmt_mpfr(e));
#endif
    if (is_rational(e)) return ftok_own(fmt_double(rational_to_double(e)));
    /* Complex */
    *prec = P_ATOM;
    FNode* r = frow();
    if (c->lang == FORM_LANG_C) fadd_tok(r, "Complex");
    fadd_tok(r, "(");
    fnode_add(r, build_real_part(c, e->data.function.args[0]), BRK_NONE);
    fadd_tok(r, ",");
    fnode_add(r, build_real_part(c, e->data.function.args[1]), BRK_ARG);
    fadd_tok(r, ")");
    return r;
}

static FNode* build_sqrt(Ctx* c, const Expr* b) {
    FNode* r = frow();
    fadd_tok(r, "Sqrt");
    fadd_tok(r, "(");
    fnode_add(r, build_wrapped(c, b, P_LOWEST), BRK_NONE);
    fadd_tok(r, ")");
    return r;
}

/* An exponent: a rational prints as its machine value. */
static FNode* build_exponent(Ctx* c, const Expr* e, int* prec) {
    if (is_rational(e)) {
        *prec = num_sign(e) < 0 ? P_UNARY : P_ATOM;
        return ftok_own(fmt_double(rational_to_double(e)));
    }
    return build(c, e, prec);
}

static FNode* build_power(Ctx* c, const Expr* b, const Expr* e, int* prec) {
    if (is_rational_value(e, 1, 2)) { *prec = P_ATOM; return build_sqrt(c, b); }
    if (is_rational_value(e, -1, 2) || is_int_value(e, -1)) {
        *prec = P_TIMES;
        FNode* r = frow();
        fadd_tok(r, "1");
        fadd_tok(r, "/");
        if (is_int_value(e, -1)) fnode_add(r, build_wrapped(c, b, P_TIMES + 1), BRK_NONE);
        else fnode_add(r, build_sqrt(c, b), BRK_NONE);
        return r;
    }
    FNode* r = frow();
    if (c->lang == FORM_LANG_C) {
        *prec = P_ATOM;
        int pe;
        fadd_tok(r, "Power");
        fadd_tok(r, "(");
        fnode_add(r, build_wrapped(c, b, P_LOWEST), BRK_NONE);
        fadd_tok(r, ",");
        fnode_add(r, build_exponent(c, e, &pe), BRK_ARG);
        fadd_tok(r, ")");
        return r;
    }
    *prec = P_POWER;
    int pe;
    fnode_add(r, build_wrapped(c, b, P_POWER + 1), BRK_NONE);
    fadd_tok(r, "**");
    FNode* ex = build_exponent(c, e, &pe);
    fnode_add(r, wrap_node(ex, pe, P_POWER), BRK_NONE);
    return r;
}

/* Is e a power whose exponent is syntactically negative (a negative number,
 * or a product with a negative numeric coefficient)? */
static bool is_negexp_power(const Expr* f) {
    if (!head_is(f, SYM_Power, 2)) return false;
    const Expr* x = f->data.function.args[1];
    if (is_number(x) && !is_complex(x)) return num_sign(x) < 0;
    if (x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol.name == SYM_Times &&
        x->data.function.arg_count >= 2) {
        const Expr* k = x->data.function.args[0];
        return is_number(k) && !is_complex(k) && num_sign(k) < 0;
    }
    return false;
}

/* -x for an exponent x accepted by is_negexp_power, as a fresh node. */
static Expr* negate_exponent(const Expr* x) {
    if (is_number(x)) return num_negate(x);
    size_t n = x->data.function.arg_count;
    const Expr* k = x->data.function.args[0];
    if (is_int_value(k, -1)) {
        if (n == 2) return expr_copy(x->data.function.args[1]);
        Expr** a = malloc((n - 1) * sizeof(Expr*));
        if (!a) abort();
        for (size_t i = 1; i < n; i++) a[i - 1] = expr_copy(x->data.function.args[i]);
        Expr* r = expr_new_function(expr_new_symbol(SYM_Times), a, n - 1);
        free(a);
        return r;
    }
    Expr** a = malloc(n * sizeof(Expr*));
    if (!a) abort();
    a[0] = num_negate(k);
    for (size_t i = 1; i < n; i++) a[i] = expr_copy(x->data.function.args[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), a, n);
    free(a);
    return r;
}

typedef struct { FNode* node; int prec; bool coef; } Factor;
typedef struct { Factor* v; size_t n, cap; } FactorList;

static void fl_push(FactorList* l, FNode* node, int prec, bool coef) {
    if (l->n == l->cap) {
        l->cap = l->cap ? 2 * l->cap : 8;
        l->v = realloc(l->v, l->cap * sizeof(Factor));
        if (!l->v) abort();
    }
    l->v[l->n].node = node;
    l->v[l->n].prec = prec;
    l->v[l->n].coef = coef;
    l->n++;
}

/* The factors of l joined by '*' (each non-coefficient factor parenthesised
 * when it binds looser than a product), with the precedence of the result. */
static FNode* join_factors(FactorList* l, int* prec) {
    if (l->n == 0) { *prec = P_ATOM; return ftok("1"); }
    if (l->n == 1) {
        *prec = l->v[0].prec;
        return l->v[0].node;
    }
    FNode* r = frow();
    for (size_t i = 0; i < l->n; i++) {
        if (i) fadd_tok(r, "*");
        Factor* f = &l->v[i];
        fnode_add(r, f->coef ? f->node : wrap_node(f->node, f->prec, P_TIMES + 1), BRK_NONE);
    }
    *prec = P_TIMES;
    return r;
}

/* A product with an optional numeric coefficient (already sign-adjusted):
 * numerator factors / denominator factors. */
static FNode* build_times_core(Ctx* c, const Expr* coef, Expr* const* fs, size_t n, int* prec) {
    FactorList num = {0}, den = {0};
    if (coef) {
        if (is_rational(coef)) {
            const Expr* p = coef->data.function.args[0];
            if (!is_int_value(p, 1)) {
                int pp;
                FNode* pn = build_number(c, p, &pp);
                fl_push(&num, pn, pp, true);
            }
            fl_push(&den, ftok_own(fmt_double(int_to_double(coef->data.function.args[1]))),
                    P_ATOM, true);
        } else {
            int pp;
            FNode* cn = build_number(c, coef, &pp);
            fl_push(&num, cn, pp, true);
        }
    }
    for (size_t i = 0; i < n; i++) {
        const Expr* f = fs[i];
        int p;
        if (is_negexp_power(f)) {
            const Expr* b = f->data.function.args[0];
            Expr* x = negate_exponent(f->data.function.args[1]);
            FNode* d = is_int_value(x, 1) ? build(c, b, &p) : build_power(c, b, x, &p);
            expr_free(x);
            fl_push(&den, d, p, false);
        } else {
            FNode* m = build(c, f, &p);
            fl_push(&num, m, p, false);
        }
    }
    FNode* result;
    if (den.n == 0) {
        result = join_factors(&num, prec);
    } else {
        int pn, pd;
        bool num_single_coef = num.n == 1 && num.v[0].coef;
        FNode* nn = join_factors(&num, &pn);
        if (num.n >= 2) nn = fparen(nn);
        else if (!num_single_coef) nn = wrap_node(nn, pn, P_TIMES + 1);
        FNode* dn = join_factors(&den, &pd);
        if (den.n >= 2) dn = fparen(dn);
        else dn = wrap_node(dn, pd, P_TIMES + 1);
        result = frow();
        fnode_add(result, nn, BRK_NONE);
        fadd_tok(result, "/");
        fnode_add(result, dn, BRK_NONE);
        *prec = P_TIMES;
    }
    free(num.v);
    free(den.v);
    return result;
}

/* Times[args] (negated when neg). */
static FNode* build_times(Ctx* c, Expr* const* args, size_t n, bool neg, int* prec) {
    const Expr* coef = NULL;
    size_t start = 0;
    if (n > 0 && is_number(args[0])) { coef = args[0]; start = 1; }
    Expr* owned = NULL;
    if (neg) {
        owned = coef ? num_negate(coef) : expr_new_integer(-1);
        coef = owned;
    }
    if (coef && is_int_value(coef, 1)) coef = NULL;
    Expr* const* rest = args + start;
    size_t nr = n - start;
    FNode* r;
    if (coef && is_int_value(coef, -1) && nr == 1 && rest[0]->type == EXPR_SYMBOL &&
        rest[0]->data.symbol.name == SYM_Infinity) {
        r = ftok("DirectedInfinity(-1)");          /* -Infinity */
        *prec = P_ATOM;
    } else if (coef && is_int_value(coef, -1)) {
        if (nr == 0) {
            r = ftok("-1");
            *prec = P_UNARY;
        } else {
            int p;
            FNode* inner = nr == 1 ? build(c, rest[0], &p) : build_times_core(c, NULL, rest, nr, &p);
            r = frow();
            fadd_tok(r, "-");
            fnode_add(r, wrap_node(inner, p, P_TIMES + 1), BRK_NONE);
            *prec = P_UNARY;
        }
    } else if (coef && is_rational(coef) && is_int_value(coef->data.function.args[0], -1) && nr > 0) {
        int p;
        r = frow();
        fnode_add(r, ftok_own(fmt_double(rational_to_double(coef))), BRK_NONE);
        fadd_tok(r, "*");
        fnode_add(r, build_times_core(c, NULL, rest, nr, &p), BRK_NONE);
        *prec = P_TIMES;
    } else if (coef && nr == 0) {
        r = build_number(c, coef, prec);
    } else {
        r = build_times_core(c, coef, rest, nr, prec);
    }
    if (owned) expr_free(owned);
    return r;
}

static bool is_negative_term(const Expr* t) {
    if (is_number(t)) return num_sign(t) < 0;
    if (head_is(t, SYM_Times, t->type == EXPR_FUNCTION ? t->data.function.arg_count : 0) &&
        t->data.function.arg_count >= 1) {
        const Expr* k = t->data.function.args[0];
        return is_number(k) && num_sign(k) < 0;
    }
    return false;
}

static FNode* build_plus(Ctx* c, const Expr* e, int* prec) {
    FNode* r = frow();
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* t = e->data.function.args[i];
        int p;
        FNode* tn;
        if (i > 0 && is_negative_term(t)) {
            fadd_tok(r, " - ");
            if (is_number(t)) {
                Expr* m = num_negate(t);
                tn = build_number(c, m, &p);
                expr_free(m);
            } else {
                tn = build_times(c, t->data.function.args, t->data.function.arg_count, true, &p);
            }
        } else {
            if (i > 0) fadd_tok(r, " + ");
            tn = build(c, t, &p);
        }
        fnode_add(r, wrap_node(tn, p, P_PLUS + 1), i > 0 ? BRK_TERM : BRK_NONE);
    }
    *prec = P_PLUS;
    return r;
}

static FNode* build_infix(Ctx* c, const Expr* e, const char* op, int myprec, int* prec) {
    FNode* r = frow();
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (i) fadd_tok(r, op);
        fnode_add(r, build_wrapped(c, e->data.function.args[i], myprec + 1), BRK_NONE);
    }
    *prec = myprec;
    return r;
}

/* head(arg,arg,...) as one flat row; a break may precede every argument
 * after the first. */
static FNode* build_call(Ctx* c, const Expr* e, const char* open, const char* close,
                         const char* sep, int* prec) {
    FNode* r = frow();
    fnode_add(r, build_wrapped(c, e->data.function.head, P_ATOM), BRK_NONE);
    fadd_tok(r, open);
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (i) fadd_tok(r, sep);
        fnode_add(r, build_wrapped(c, e->data.function.args[i], P_LOWEST), i ? BRK_ARG : BRK_NONE);
    }
    fadd_tok(r, close);
    *prec = P_ATOM;
    return r;
}

static char* quote_string(const char* s) {
    SBuf b = {0};
    sb_puts(&b, "\"");
    for (const char* p = s; *p; p++) {
        if (*p == '"' || *p == '\\') sb_putn(&b, "\\", 1);
        sb_putn(&b, p, 1);
    }
    sb_puts(&b, "\"");
    return sb_take(&b);
}

static const char* rel_op(const char* h, FormLang lang, int* p) {
    bool f = lang == FORM_LANG_FORTRAN;
    *p = P_REL;
    if (h == SYM_Equal) return f ? ".eq." : " == ";
    if (h == SYM_Unequal) return f ? ".ne." : " != ";
    if (h == SYM_Less) return f ? ".lt." : " < ";
    if (h == SYM_LessEqual) return f ? ".le." : " <= ";
    if (h == SYM_Greater) return f ? ".gt." : " > ";
    if (h == SYM_GreaterEqual) return f ? ".ge." : " >= ";
    *p = P_AND;
    if (h == SYM_And) return f ? ".and." : " && ";
    *p = P_OR;
    if (h == SYM_Or) return f ? ".or." : " || ";
    return NULL;
}

static FNode* build(Ctx* c, const Expr* e, int* prec) {
    *prec = P_ATOM;
    if (!e) return ftok("Null");

    /* Format[f[...], FortranForm] := rhs (or CForm) */
    if ((e->type == EXPR_FUNCTION || e->type == EXPR_SYMBOL) && c->fmt_depth < FMT_DEPTH_LIMIT) {
        Expr* shown = print_form_format_rewrite((Expr*)e, c->form);
        if (shown) {
            c->fmt_depth++;
            FNode* r = build(c, shown, prec);
            c->fmt_depth--;
            expr_free(shown);
            return r;
        }
    }

    switch (e->type) {
    case EXPR_INTEGER:
    case EXPR_BIGINT:
    case EXPR_REAL:
#ifdef USE_MPFR
    case EXPR_MPFR:
#endif
        return build_number(c, e, prec);
    case EXPR_STRING:
        return ftok_own(quote_string(e->data.string));
    case EXPR_SYMBOL: {
        const char* nm = e->data.symbol.name;
        if (nm == SYM_Infinity) return ftok("DirectedInfinity(1)");
        if (nm == SYM_ComplexInfinity) return ftok("DirectedInfinity()");
        return ftok(context_display_name(nm));
    }
    case EXPR_FUNCTION:
        break;
    default: {
        char* s = expr_to_string_outputform((Expr*)e);
        return ftok_own(s ? s : (char*)calloc(1, 1));
    }
    }

    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    const char* h = hd->type == EXPR_SYMBOL ? hd->data.symbol.name : NULL;
    if (h) {
        if ((h == SYM_FortranForm || h == SYM_CForm || h == SYM_HoldForm) && n == 1)
            return build(c, e->data.function.args[0], prec);
        if (is_number(e)) return build_number(c, e, prec);
        if (h == SYM_Plus && n >= 2) return build_plus(c, e, prec);
        if (h == SYM_Times && n >= 1) return build_times(c, e->data.function.args, n, false, prec);
        if (h == SYM_Power && n == 2)
            return build_power(c, e->data.function.args[0], e->data.function.args[1], prec);
        if (h == SYM_SequenceForm) {
            FNode* r = frow();
            for (size_t i = 0; i < n; i++) {
                int p;
                fnode_add(r, build(c, e->data.function.args[i], &p), BRK_NONE);
            }
            *prec = P_ATOM;
            return r;
        }
        int rp;
        const char* op = rel_op(h, c->lang, &rp);
        if (op && n >= 2 && (rp != P_REL || n == 2)) return build_infix(c, e, op, rp, prec);
        if (h == SYM_Not && n == 1) {
            FNode* r = frow();
            fadd_tok(r, c->lang == FORM_LANG_FORTRAN ? ".not." : "!");
            fnode_add(r, build_wrapped(c, e->data.function.args[0], P_NOT + 1), BRK_NONE);
            *prec = P_NOT;
            return r;
        }
        if (c->lang == FORM_LANG_C) {
            if (h == SYM_Mod && n == 2) return build_infix(c, e, " % ", P_TIMES, prec);
            if (h == SYM_Part && n >= 2) {
                FNode* r = frow();
                fnode_add(r, build_wrapped(c, e->data.function.args[0], P_ATOM), BRK_NONE);
                for (size_t i = 1; i < n; i++) {
                    fadd_tok(r, "[");
                    fnode_add(r, build_wrapped(c, e->data.function.args[i], P_LOWEST), BRK_NONE);
                    fadd_tok(r, "]");
                }
                *prec = P_ATOM;
                return r;
            }
        }
    }
    return build_call(c, e, "(", ")", ",", prec);
}

/* ------------------------------------------------------------------ */
/* Layout                                                              */
/* ------------------------------------------------------------------ */

typedef struct {
    SBuf out;
    int width;
    size_t col;        /* current column */
    bool fresh;        /* nothing written since the last line start */
    FormLang lang;
} Lay;

static void flat(const FNode* f, SBuf* b) {
    if (f->text) { sb_putn(b, f->text, f->len); return; }
    for (size_t i = 0; i < f->n; i++) flat(f->kids[i], b);
}

/* Continuation line prefix: fixed-form Fortran "     -" and depth + 1
 * spaces; C depth + 2 spaces. */
static size_t cont_len(const Lay* L, int depth) {
    return L->lang == FORM_LANG_FORTRAN ? 6 + (size_t)depth + 1 : (size_t)depth + 2;
}

static void newline(Lay* L, int depth) {
    sb_puts(&L->out, L->lang == FORM_LANG_FORTRAN ? "\n     -" : "\n");
    size_t pad = L->lang == FORM_LANG_FORTRAN ? (size_t)depth + 1 : (size_t)depth + 2;
    for (size_t i = 0; i < pad; i++) sb_puts(&L->out, " ");
    L->col = cont_len(L, depth);
    L->fresh = true;
}

static void emit_text(Lay* L, const char* s, size_t n) {
    sb_putn(&L->out, s, n);
    L->col += n;
    if (n) L->fresh = false;
}

/* A token at nesting depth `depth`.  One that cannot fit even on a fresh
 * continuation line is split with a backslash. */
static void emit_token(Lay* L, const FNode* t, int depth) {
    const char* s = t->text;
    size_t n = t->len;
    size_t w = (size_t)L->width;
    if (L->col + n <= w || n + cont_len(L, depth) <= w) {
        emit_text(L, s, n);
        return;
    }
    while (L->col + n > w) {
        size_t room = (w > L->col + 1) ? w - L->col - 1 : 1;
        if (room >= n) break;
        emit_text(L, s, room);
        sb_puts(&L->out, "\\");
        newline(L, depth);
        s += room;
        n -= room;
    }
    emit_text(L, s, n);
}

static void emit(Lay* L, const FNode* f, int depth) {
    if (f->text) { emit_token(L, f, depth); return; }
    size_t w = (size_t)L->width;
    for (size_t i = 0; i < f->n; i++) {
        const FNode* k = f->kids[i];
        if (f->brk[i] != BRK_NONE && !L->fresh) {
            /* the child and the separator token that follows it */
            size_t need = k->len;
            if (i + 1 < f->n && f->kids[i + 1]->text) need += f->kids[i + 1]->len;
            if (L->col + need > w) {
                if (f->brk[i] == BRK_ARG || need + cont_len(L, depth) <= w)
                    newline(L, depth);
            }
        }
        if (k->text) emit_token(L, k, depth + 1);
        else emit(L, k, depth + 1);
    }
}

static char* layout(FNode* root, FormLang lang, int width) {
    SBuf b = {0};
    if (width <= 0 || root->len <= (size_t)width) {
        flat(root, &b);
        return sb_take(&b);
    }
    Lay L = { b, width, 0, false, lang };
    if (lang == FORM_LANG_FORTRAN) emit_text(&L, "        ", 8);   /* fixed-form statement column */
    L.fresh = true;
    emit(&L, root, 1);
    return sb_take(&L.out);
}

char* fortran_form_string(const Expr* e, FormLang lang, int page_width) {
    Ctx c = { lang, lang == FORM_LANG_C ? SYM_CForm : SYM_FortranForm, 0 };
    int p;
    FNode* root = build(&c, e, &p);
    char* s = layout(root, lang, page_width);
    ffree(root);
    return s;
}

char* fortran_form_string_seq(Expr* const* es, size_t n, FormLang lang, int page_width) {
    if (n == 1) return fortran_form_string(es[0], lang, page_width);
    Ctx c = { lang, lang == FORM_LANG_C ? SYM_CForm : SYM_FortranForm, 0 };
    FNode* root = frow();
    for (size_t i = 0; i < n; i++) {
        int p;
        fnode_add(root, build(&c, es[i], &p), BRK_NONE);
    }
    char* s = layout(root, lang, page_width);
    ffree(root);
    return s;
}

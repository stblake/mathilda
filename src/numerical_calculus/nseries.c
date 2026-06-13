/*
 * nseries.c — NSeries[f, {x, x0, n}, opts]
 *
 * Numerical Taylor/Laurent series expansion of `f` about x = x0, including the
 * terms (x - x0)^-n through (x - x0)^n, returned as a SeriesData object.
 *
 *   NSeries[f, {x, x0, n}]
 *
 * Method (Lyness & Sande, 1971; Bornemann, FoCM 2011). f is sampled at N
 * equispaced points on a circle of radius r centred at x0,
 *
 *     z_j = x0 + r e^(2 pi i j / N),   j = 0 .. N-1,
 *
 * and a discrete Fourier transform of the samples recovers the Laurent
 * coefficients via Cauchy's integral formula:
 *
 *     c_k = (1/N) sum_j f(z_j) e^(-2 pi i j k / N)        (DFT of the samples)
 *     a_e = c_(e mod N) * r^(-e)         for e = -n .. n  (coeff of (x-x0)^e)
 *
 * The upper-half DFT bins (k = N-m) supply the NEGATIVE-power coefficients, so
 * one transform yields both the principal part and the analytic part. This is
 * exact (no truncation) when f is analytic on an annulus containing the circle;
 * it fails if the disk centred at x0 contains a branch cut of f.
 *
 * N is chosen as a power of two with an oversampling margin,
 *     N = 2^(ceil(log2 n) + 2),
 * so the leading aliased term a_(k +/- N) r^(+/- N) is pushed below the
 * round-off floor. Because round-off breaks the conjugate symmetry of
 * real-coefficient functions, the result carries tiny spurious residuals —
 * Chop the result when needed.
 *
 * A direct O(N^2) DFT is used rather than an FFT: N is small (<= a few hundred)
 * and each sample requires a full symbolic evaluation of f, which dominates the
 * runtime by orders of magnitude. The same code path serves both the machine
 * (double _Complex) and arbitrary-precision (MPFR) computations; no
 * double-precision FFT library could serve the MPFR path anyway.
 *
 * Options (trailing Rule[...] in any order):
 *   Radius           -> r                contour radius (default 1)
 *   WorkingPrecision -> MachinePrecision | digits
 *
 * Memory: receives `res` owned by the evaluator. Returns a fresh Expr* (a
 * SeriesData) on success or NULL (unevaluated). Never frees `res`. All
 * temporary OwnValues are removed before returning, on every path.
 */

#include "nseries.h"

#include <complex.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#  include "numeric_complex.h"   /* mpfr_complex_exp (unused here, kept for parity) */
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    double radius;       /* contour radius (default 1)            */
    bool   prec_mpfr;    /* WorkingPrecision selects MPFR         */
    long   bits;         /* MPFR working precision in bits        */
} NsOpts;

static void ns_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NSeries::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers                                   *
 * ------------------------------------------------------------------ */

static bool ns_to_double_real(Expr* e, double* out) {
    if (!e) return false;
    int64_t rn, rd;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational(e, &rn, &rd)) { *out = (double)rn / (double)rd; return true; }
            return false;
        default: return false;
    }
}

static bool ns_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (ns_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (ns_to_double_real(re, &r) && ns_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* ns_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
static Expr* ns_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}
#endif

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors NResidue)                    *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;      /* interned                           */
    Rule*       saved_own; /* prior OwnValue chain (borrowed)    */
    uint32_t    saved_attrs;
    bool        valid;
} NsBind;

static void ns_bind_snapshot(NsBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void ns_bind_clear_temp(SymbolDef* def) {
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
}

static void ns_bind_set(NsBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    ns_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

static void ns_bind_restore(NsBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    ns_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Sampler context + callbacks                                        *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       f;        /* borrowed integrand (variable free)  */
    NsBind*     bind;     /* the active variable binding         */
    NumericSpec spec;     /* machine or MPFR                     */
    long        bits;     /* MPFR bits (0 => machine)            */
} NsCtx;

/* Evaluate f after binding the variable to `value` (consumed). Returns the
 * numericalized result (caller frees) or NULL. */
static Expr* ns_eval_at(NsCtx* c, Expr* value) {
    ns_bind_set(c->bind, value);   /* takes a copy internally */
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

static bool ns_sample_machine(NsCtx* c, double _Complex z, double _Complex* out) {
    double r = creal(z), i = cimag(z);
    Expr* zv = (i == 0.0) ? expr_new_real(r)
                          : make_complex(expr_new_real(r), expr_new_real(i));
    Expr* num = ns_eval_at(c, zv);
    if (!num) return false;
    bool ok = ns_to_complex(num, out);
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool ns_sample_mpfr(NsCtx* c, const mpfr_t z_re, const mpfr_t z_im,
                           mpfr_t out_re, mpfr_t out_im) {
    Expr* zv = mpfr_zero_p(z_im) ? expr_new_mpfr_copy(z_re)
                                 : make_complex(expr_new_mpfr_copy(z_re),
                                                expr_new_mpfr_copy(z_im));
    Expr* num = ns_eval_at(c, zv);
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    expr_free(num);
    return ok;
}
#endif

/* ------------------------------------------------------------------ *
 *  Sample count: N = 2^(ceil(log2 n) + 2), clamped so N >= 4, N > 2n   *
 * ------------------------------------------------------------------ */

static int ns_sample_count(int n) {
    int ne = n < 1 ? 1 : n;
    int p = 0;
    while ((1 << p) < ne) p++;            /* p = ceil(log2 ne)          */
    int N = 1 << (p + 2);                 /* 4x oversampling margin     */
    if (N < 4) N = 4;
    while (N <= 2 * n) N <<= 1;           /* hold the full -n..n window */
    return N;
}

/* ------------------------------------------------------------------ *
 *  SeriesData assembly                                                *
 *                                                                     *
 *  coefs[0 .. 2n] are the coefficients of (x-x0)^(-n) .. (x-x0)^n     *
 *  (the array is adopted; elements live on inside the result).        *
 * ------------------------------------------------------------------ */

static Expr* ns_build_seriesdata(const char* varname, Expr* x0_orig,
                                 Expr** coefs, int n) {
    size_t cnt = (size_t)(2 * n + 1);
    Expr* coefs_list = expr_new_function(expr_new_symbol("List"), coefs, cnt);
    Expr* args[6] = {
        expr_new_symbol(varname),
        expr_copy(x0_orig),
        coefs_list,
        expr_new_integer(-(int64_t)n),       /* nmin */
        expr_new_integer((int64_t)n + 1),    /* nmax (O-term exponent) */
        expr_new_integer(1)                  /* den (no fractional powers) */
    };
    return expr_new_function(expr_new_symbol("SeriesData"), args, 6);
}

/* ------------------------------------------------------------------ *
 *  Machine-precision DFT path                                         *
 * ------------------------------------------------------------------ */

static Expr* ns_compute_machine(NsCtx* ctx, const char* varname, Expr* x0_orig,
                                double _Complex z0, double r, int n) {
    int N = ns_sample_count(n);
    double _Complex* F = malloc(sizeof(double _Complex) * (size_t)N);
    if (!F) return NULL;

    for (int j = 0; j < N; j++) {
        double theta = 2.0 * M_PI * (double)j / (double)N;
        double _Complex z = z0 + r * (cos(theta) + I * sin(theta));
        if (!ns_sample_machine(ctx, z, &F[j])) { free(F); return NULL; }
    }

    size_t cnt = (size_t)(2 * n + 1);
    Expr** coefs = malloc(sizeof(Expr*) * cnt);
    if (!coefs) { free(F); return NULL; }

    for (int e = -n; e <= n; e++) {
        int k = ((e % N) + N) % N;
        double _Complex acc = 0.0;
        for (int j = 0; j < N; j++) {
            double ang = -2.0 * M_PI * (double)j * (double)k / (double)N;
            acc += F[j] * (cos(ang) + I * sin(ang));
        }
        acc /= (double)N;
        acc *= pow(r, -(double)e);        /* a_e = c_k * r^(-e) */
        coefs[e + n] = ns_from_complex_d(acc);
    }
    free(F);

    Expr* sd = ns_build_seriesdata(varname, x0_orig, coefs, n);
    free(coefs);
    return sd;
}

/* ------------------------------------------------------------------ *
 *  Arbitrary-precision (MPFR) DFT path                                *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR
static Expr* ns_compute_mpfr(NsCtx* ctx, const char* varname, Expr* x0_orig,
                             const mpfr_t z0r, const mpfr_t z0i,
                             double r, long bits, int n) {
    int N = ns_sample_count(n);

    /* Sample arrays: real and imaginary parts of f(z_j). */
    mpfr_t* Fr = malloc(sizeof(mpfr_t) * (size_t)N);
    mpfr_t* Fi = malloc(sizeof(mpfr_t) * (size_t)N);
    if (!Fr || !Fi) { free(Fr); free(Fi); return NULL; }

    mpfr_t pi2, theta, s, c, zr, zi, rr;
    mpfr_inits2(bits, pi2, theta, s, c, zr, zi, rr, (mpfr_ptr)0);
    mpfr_const_pi(pi2, MPFR_RNDN);
    mpfr_mul_2ui(pi2, pi2, 1, MPFR_RNDN);     /* 2*pi */
    mpfr_set_d(rr, r, MPFR_RNDN);

    bool ok = true;
    int filled = 0;
    for (int j = 0; j < N && ok; j++) {
        /* theta = 2*pi*j/N */
        mpfr_mul_si(theta, pi2, j, MPFR_RNDN);
        mpfr_div_si(theta, theta, N, MPFR_RNDN);
        mpfr_sin_cos(s, c, theta, MPFR_RNDN);
        /* z = z0 + r*(cos + i sin) */
        mpfr_fma(zr, rr, c, z0r, MPFR_RNDN);
        mpfr_fma(zi, rr, s, z0i, MPFR_RNDN);

        mpfr_init2(Fr[j], bits);
        mpfr_init2(Fi[j], bits);
        filled = j + 1;
        if (!ns_sample_mpfr(ctx, zr, zi, Fr[j], Fi[j])) ok = false;
    }

    Expr* sd = NULL;
    if (ok) {
        size_t cnt = (size_t)(2 * n + 1);
        Expr** coefs = malloc(sizeof(Expr*) * cnt);
        if (coefs) {
            mpfr_t accr, acci, ang, wc, ws, t1, t2, scale;
            mpfr_inits2(bits, accr, acci, ang, wc, ws, t1, t2, scale, (mpfr_ptr)0);

            for (int e = -n; e <= n; e++) {
                int k = ((e % N) + N) % N;
                mpfr_set_zero(accr, 1);
                mpfr_set_zero(acci, 1);
                for (int j = 0; j < N; j++) {
                    /* twiddle = cos(ang) + i sin(ang), ang = -2*pi*j*k/N */
                    mpfr_mul_si(ang, pi2, j, MPFR_RNDN);
                    mpfr_mul_si(ang, ang, k, MPFR_RNDN);
                    mpfr_div_si(ang, ang, N, MPFR_RNDN);
                    mpfr_neg(ang, ang, MPFR_RNDN);
                    mpfr_sin_cos(ws, wc, ang, MPFR_RNDN);
                    /* acc += F[j] * twiddle */
                    mpfr_mul(t1, Fr[j], wc, MPFR_RNDN);
                    mpfr_mul(t2, Fi[j], ws, MPFR_RNDN);
                    mpfr_sub(t1, t1, t2, MPFR_RNDN);
                    mpfr_add(accr, accr, t1, MPFR_RNDN);
                    mpfr_mul(t1, Fr[j], ws, MPFR_RNDN);
                    mpfr_mul(t2, Fi[j], wc, MPFR_RNDN);
                    mpfr_add(t1, t1, t2, MPFR_RNDN);
                    mpfr_add(acci, acci, t1, MPFR_RNDN);
                }
                /* c_k = acc / N ; a_e = c_k * r^(-e) */
                mpfr_div_si(accr, accr, N, MPFR_RNDN);
                mpfr_div_si(acci, acci, N, MPFR_RNDN);
                mpfr_pow_si(scale, rr, -e, MPFR_RNDN);
                mpfr_mul(accr, accr, scale, MPFR_RNDN);
                mpfr_mul(acci, acci, scale, MPFR_RNDN);
                coefs[e + n] = ns_from_complex_mpfr(accr, acci);
            }
            mpfr_clears(accr, acci, ang, wc, ws, t1, t2, scale, (mpfr_ptr)0);
            sd = ns_build_seriesdata(varname, x0_orig, coefs, n);
            free(coefs);
        }
    }

    for (int j = 0; j < filled; j++) { mpfr_clear(Fr[j]); mpfr_clear(Fi[j]); }
    free(Fr); free(Fi);
    mpfr_clears(pi2, theta, s, c, zr, zi, rr, (mpfr_ptr)0);
    return sd;
}
#endif

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */

static bool ns_is_known_option(const char* s) {
    return s == SYM_Radius || s == SYM_WorkingPrecision;
}

static bool ns_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && ns_is_known_option(lhs->data.symbol);
}

static bool ns_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!ns_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool ns_apply_option(Expr* rule, NsOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Radius) {
        double r;
        if (!ns_to_double_real(rhs, &r) || r <= 0.0) {
            ns_warn("badopt", "Radius must be a positive number");
            return false;
        }
        o->radius = r;
        return true;
    }
    if (name == SYM_WorkingPrecision) {
        if (!ns_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            ns_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nseries(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing options; exactly two positional args (f, spec) remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && ns_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!ns_is_option_arg(res->data.function.args[i])) {
            ns_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) return NULL;

    NsOpts opts;
    opts.radius = 1.0;        /* default radius 1 */
    opts.prec_mpfr = false;
    opts.bits = 0;
    for (size_t i = pos_end; i < argc; i++) {
        if (!ns_apply_option(res->data.function.args[i], &opts)) return NULL;
    }

    /* Parse the {x, x0, n} spec. */
    Expr* spec = res->data.function.args[1];
    if (spec->type != EXPR_FUNCTION
        || spec->data.function.head->type != EXPR_SYMBOL
        || spec->data.function.head->data.symbol != SYM_List
        || spec->data.function.arg_count != 3) {
        ns_warn("ivar", "second argument must be {x, x0, n}");
        return NULL;
    }
    Expr* var = spec->data.function.args[0];
    Expr* x0_orig = spec->data.function.args[1];
    Expr* n_expr = spec->data.function.args[2];
    if (var->type != EXPR_SYMBOL) {
        ns_warn("ivar", "the variable in {x, x0, n} must be a symbol");
        return NULL;
    }
    if (n_expr->type != EXPR_INTEGER || n_expr->data.integer < 0) {
        ns_warn("ivar", "the order n in {x, x0, n} must be a non-negative integer");
        return NULL;
    }
    int n = (int)n_expr->data.integer;

    /* Numeric spec. */
    NumericSpec spec_num;
#ifdef USE_MPFR
    if (opts.prec_mpfr) { spec_num.mode = NUMERIC_MODE_MPFR; spec_num.bits = opts.bits;
                          spec_num.preserve_inexact = false; }
    else                  spec_num = numeric_machine_spec();
#else
    spec_num = numeric_machine_spec();
#endif

    /* Numericalise x0 for sampling (keep x0_orig for the SeriesData base). */
    Expr* z0_raw = eval_and_free(expr_copy(x0_orig));
    Expr* z0 = z0_raw ? numericalize(z0_raw, spec_num) : NULL;
    expr_free(z0_raw);
    if (!z0) { ns_warn("nnum", "x0 is not numeric"); return NULL; }

    NsBind bind;
    ns_bind_snapshot(&bind, var->data.symbol);

    NsCtx ctx;
    ctx.f = res->data.function.args[0];
    ctx.bind = &bind;
    ctx.spec = spec_num;
    ctx.bits = opts.bits;

    Expr* result = NULL;

#ifdef USE_MPFR
    if (opts.prec_mpfr) {
        mpfr_t z0r, z0i;
        mpfr_init2(z0r, opts.bits); mpfr_init2(z0i, opts.bits);
        bool inexact;
        if (!get_approx_mpfr(z0, z0r, z0i, &inexact)) {
            ns_warn("nnum", "x0 is not numeric");
        } else {
            result = ns_compute_mpfr(&ctx, var->data.symbol, x0_orig,
                                     z0r, z0i, opts.radius, opts.bits, n);
            if (!result)
                ns_warn("nnum", "f could not be evaluated to a number on the contour");
        }
        mpfr_clears(z0r, z0i, (mpfr_ptr)0);
    } else
#endif
    {
        double _Complex z0c;
        if (!ns_to_complex(z0, &z0c)) {
            ns_warn("nnum", "x0 is not numeric");
        } else {
            result = ns_compute_machine(&ctx, var->data.symbol, x0_orig,
                                        z0c, opts.radius, n);
            if (!result)
                ns_warn("nnum", "f could not be evaluated to a number on the contour");
        }
    }

    ns_bind_restore(&bind);
    expr_free(z0);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nseries_init(void) {
    symtab_add_builtin("NSeries", builtin_nseries);
    /* Protected only. NOT ATTR_LISTABLE — generic threading would split the
     * {x, x0, n} spec across bogus calls. */
    symtab_get_def("NSeries")->attributes |= ATTR_PROTECTED;
}

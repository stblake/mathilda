/*
 * nderiv.c — ND[expr, x, x0] / ND[expr, {x, n}, x0, opts]
 *
 * Numerical approximation to the (n-th) derivative of `expr` w.r.t. `x` at
 * x = x0. Two methods, selected by Method (default EulerSum):
 *
 *   Method -> EulerSum   Richardson (Romberg/Neville) extrapolation of the
 *                        n-th *forward* finite difference taken along the
 *                        complex direction Scale:
 *
 *                          D(h) = (1/(s h)^n) sum_{k=0}^n (-1)^{n-k} C(n,k)
 *                                              f(x0 + k s h),   s = Scale,
 *                          h_i = 2^-i  (i = 0 .. Terms-1),
 *                          T(i,0) = D(h_i),
 *                          T(i,j) = T(i,j-1)
 *                                   + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
 *                          result = T(Terms-1, Terms-1).
 *
 *                        The forward (one-sided) stencil along `s` is what
 *                        gives directional/one-sided derivatives — e.g. the
 *                        left/right derivatives of Abs and the complex
 *                        direction Scale -> 1 + I. The error expansion of a
 *                        forward difference runs in *all* powers of h, hence
 *                        the (2^j - 1) Richardson denominator (not 4^j - 1).
 *                        Works for non-analytic f (only samples along `s`).
 *                        Requires integer order n >= 1.
 *
 *   Method -> NIntegrate Cauchy integral formula via the existing NResidue:
 *                          f^(n)(x0) = n! Res_{z=x0} f(z)/(z-x0)^(n+1)
 *                                    = Gamma(n+1) *
 *                                      NResidue[expr/(x-x0)^(n+1), {x, x0},
 *                                               Radius -> Scale].
 *                        Scale is the contour radius (default 1; NResidue's
 *                        own tiny 1/100 default would cause heavy 1/r^n
 *                        cancellation, so we always pass Radius -> Scale).
 *                        Gamma(n+1) (rather than n!) gives fractional /
 *                        complex order. Requires expr analytic near x0.
 *
 * Options: Method, Scale (default 1; step for EulerSum / contour radius for
 * NIntegrate; may be complex), Terms (default 7; EulerSum tableau depth),
 * WorkingPrecision (MachinePrecision | digits -> MPFR), PrecisionGoal,
 * MaxRecursion (NIntegrate only). Threads element-wise over a List in arg 1.
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated). Never frees `res`. Any temporary OwnValue
 * created for the EulerSum sampler is removed on every return path.
 */

#include "nderiv.h"
#include "nresidue.h"      /* NIntegrate reuses NResidue via the evaluator */

#include <complex.h>
#include <gmp.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#  include "numeric_complex.h"   /* mpfr_complex_div */
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* Hard cap on Terms: 2^(Terms-1) must stay an exact double for the machine
 * Richardson denominator, and any larger tableau is numerical nonsense. */
#define ND_MAX_TERMS 50

static void nd_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "ND::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers (mirrors nresidue.c)              *
 * ------------------------------------------------------------------ */

static bool nd_to_double_real(Expr* e, double* out) {
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

static bool nd_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (nd_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (nd_to_double_real(re, &r) && nd_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* nd_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
static Expr* nd_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}
#endif

/* Binomial coefficient C(n, k) as an exact integer (via GMP). */
static double nd_binom_d(int n, int k) {
    mpz_t b; mpz_init(b);
    mpz_bin_uiui(b, (unsigned long)n, (unsigned long)k);
    double r = mpz_get_d(b);
    mpz_clear(b);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors nresidue.c)                  *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} NdBind;

static void nd_bind_snapshot(NdBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void nd_bind_clear_temp(SymbolDef* def) {
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

static void nd_bind_set(NdBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    nd_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

static void nd_bind_restore(NdBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    nd_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Sampler: evaluate expr with the variable bound to a complex point  *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       f;     /* borrowed integrand (free in the variable)  */
    NdBind*     bind;
    NumericSpec spec;
} NdCtx;

/* Evaluate f after binding the variable to `value` (consumed). Returns the
 * numericalised result (caller frees) or NULL. */
static Expr* nd_eval_at(NdCtx* c, Expr* value) {
    nd_bind_set(c->bind, value);   /* takes a copy internally */
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

static bool nd_sample_machine(NdCtx* c, double _Complex z, double _Complex* out) {
    Expr* num = nd_eval_at(c, nd_from_complex_d(z));
    if (!num) return false;
    bool ok = nd_to_complex(num, out);
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nd_sample_mpfr(NdCtx* c, const mpfr_t z_re, const mpfr_t z_im,
                           mpfr_t out_re, mpfr_t out_im) {
    Expr* num = nd_eval_at(c, nd_from_complex_mpfr(z_re, z_im));
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    expr_free(num);
    return ok;
}
#endif

/* Numericalise an expression to a machine complex value (NULL expr => 1). */
static bool nd_numeric_complex_machine(Expr* e, double _Complex* out) {
    if (!e) { *out = 1.0 + 0.0 * I; return true; }
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
    expr_free(raw);
    bool ok = num && nd_to_complex(num, out);
    expr_free(num);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  EulerSum — machine precision                                       *
 * ------------------------------------------------------------------ */

static Expr* nd_eulersum_machine(Expr* expr, const char* var, int n,
                                 Expr* x0_expr, Expr* scale_expr, int terms) {
    NumericSpec spec = numeric_machine_spec();

    double _Complex x0, s;
    {
        Expr* raw = eval_and_free(expr_copy(x0_expr));
        Expr* num = raw ? numericalize(raw, spec) : NULL;
        expr_free(raw);
        bool ok = num && nd_to_complex(num, &x0);
        expr_free(num);
        if (!ok) { nd_warn("nnum", "the point x0 is not numeric"); return NULL; }
    }
    if (!nd_numeric_complex_machine(scale_expr, &s) || s == 0.0) {
        nd_warn("badscl", "Scale must evaluate to a nonzero number");
        return NULL;
    }

    NdBind bind; nd_bind_snapshot(&bind, var);
    NdCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = spec;

    double _Complex* T = malloc(sizeof(double _Complex) * (size_t)terms * terms);
    bool ok = (T != NULL);

    for (int i = 0; ok && i < terms; i++) {
        double _Complex step = s * ldexp(1.0, -i);   /* s * 2^-i */
        double _Complex acc = 0.0;
        for (int k = 0; k <= n; k++) {
            double _Complex fk;
            if (!nd_sample_machine(&ctx, x0 + (double)k * step, &fk)) { ok = false; break; }
            double sign = (((n - k) & 1) == 0) ? 1.0 : -1.0;
            acc += sign * nd_binom_d(n, k) * fk;
        }
        if (!ok) break;
        double _Complex sn = 1.0;
        for (int p = 0; p < n; p++) sn *= step;
        T[(size_t)i * terms + 0] = acc / sn;
    }

    double _Complex result = 0.0;
    if (ok) {
        for (int j = 1; j < terms; j++) {
            double denom = ldexp(1.0, j) - 1.0;       /* 2^j - 1 */
            for (int i = j; i < terms; i++) {
                double _Complex a = T[(size_t)i * terms + (j - 1)];
                double _Complex b = T[(size_t)(i - 1) * terms + (j - 1)];
                T[(size_t)i * terms + j] = a + (a - b) / denom;
            }
        }
        result = T[(size_t)(terms - 1) * terms + (terms - 1)];
    }

    free(T);
    nd_bind_restore(&bind);
    if (!ok) {
        nd_warn("nnum", "expr did not evaluate to a number at a sample point");
        return NULL;
    }
    return nd_from_complex_d(result);
}

/* ------------------------------------------------------------------ *
 *  EulerSum — MPFR precision                                          *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* (out_re, out_im) = (a_re + a_im i)(b_re + b_im i). Output may alias inputs. */
static void nd_cmul(mpfr_t out_re, mpfr_t out_im,
                    const mpfr_t a_re, const mpfr_t a_im,
                    const mpfr_t b_re, const mpfr_t b_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t t1, t2, t3, t4;
    mpfr_inits2(p, t1, t2, t3, t4, (mpfr_ptr)0);
    mpfr_mul(t1, a_re, b_re, MPFR_RNDN);
    mpfr_mul(t2, a_im, b_im, MPFR_RNDN);
    mpfr_mul(t3, a_re, b_im, MPFR_RNDN);
    mpfr_mul(t4, a_im, b_re, MPFR_RNDN);
    mpfr_sub(out_re, t1, t2, MPFR_RNDN);
    mpfr_add(out_im, t3, t4, MPFR_RNDN);
    mpfr_clears(t1, t2, t3, t4, (mpfr_ptr)0);
}

static Expr* nd_eulersum_mpfr(Expr* expr, const char* var, int n,
                              Expr* x0_expr, Expr* scale_expr,
                              int terms, long bits) {
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = bits;
    spec.preserve_inexact = false;
    mpfr_prec_t p = (mpfr_prec_t)bits;

    mpfr_t x0r, x0i, sr, si;
    mpfr_init2(x0r, p); mpfr_init2(x0i, p);
    mpfr_init2(sr, p);  mpfr_init2(si, p);

    bool prep_ok = true;
    {
        Expr* raw = eval_and_free(expr_copy(x0_expr));
        Expr* num = raw ? numericalize(raw, spec) : NULL;
        expr_free(raw);
        bool inexact;
        if (!num || !get_approx_mpfr(num, x0r, x0i, &inexact)) prep_ok = false;
        expr_free(num);
    }
    if (prep_ok) {
        if (scale_expr) {
            Expr* raw = eval_and_free(expr_copy(scale_expr));
            Expr* num = raw ? numericalize(raw, spec) : NULL;
            expr_free(raw);
            bool inexact;
            if (!num || !get_approx_mpfr(num, sr, si, &inexact)) prep_ok = false;
            expr_free(num);
        } else {
            mpfr_set_ui(sr, 1, MPFR_RNDN);
            mpfr_set_ui(si, 0, MPFR_RNDN);
        }
    }
    if (prep_ok && mpfr_zero_p(sr) && mpfr_zero_p(si)) prep_ok = false;
    if (!prep_ok) {
        mpfr_clears(x0r, x0i, sr, si, (mpfr_ptr)0);
        nd_warn("nnum", "the point x0 or Scale is not numeric/nonzero");
        return NULL;
    }

    /* Tableau storage. */
    size_t cells = (size_t)terms * terms;
    mpfr_t* Tr = malloc(sizeof(mpfr_t) * cells);
    mpfr_t* Ti = malloc(sizeof(mpfr_t) * cells);
    for (size_t c = 0; c < cells; c++) { mpfr_init2(Tr[c], p); mpfr_init2(Ti[c], p); }

    /* Scratch. */
    mpfr_t hr, str, sti, snr, sni, accr, acci, ptr, pti, fr, fi,
           binr, tr, ti, denom;
    mpfr_inits2(p, hr, str, sti, snr, sni, accr, acci, ptr, pti, fr, fi,
                binr, tr, ti, denom, (mpfr_ptr)0);

    NdBind bind; nd_bind_snapshot(&bind, var);
    NdCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = spec;

    bool ok = true;
    for (int i = 0; ok && i < terms; i++) {
        /* step = s * 2^-i */
        mpfr_set_ui(hr, 1, MPFR_RNDN);
        mpfr_div_2ui(hr, hr, (unsigned long)i, MPFR_RNDN);
        mpfr_mul(str, sr, hr, MPFR_RNDN);
        mpfr_mul(sti, si, hr, MPFR_RNDN);

        mpfr_set_ui(accr, 0, MPFR_RNDN);
        mpfr_set_ui(acci, 0, MPFR_RNDN);
        for (int k = 0; ok && k <= n; k++) {
            /* pt = x0 + k*step */
            mpfr_mul_ui(ptr, str, (unsigned long)k, MPFR_RNDN);
            mpfr_mul_ui(pti, sti, (unsigned long)k, MPFR_RNDN);
            mpfr_add(ptr, ptr, x0r, MPFR_RNDN);
            mpfr_add(pti, pti, x0i, MPFR_RNDN);
            if (!nd_sample_mpfr(&ctx, ptr, pti, fr, fi)) { ok = false; break; }
            /* coeff = (-1)^{n-k} C(n,k) */
            mpz_t b; mpz_init(b);
            mpz_bin_uiui(b, (unsigned long)n, (unsigned long)k);
            mpfr_set_z(binr, b, MPFR_RNDN);
            mpz_clear(b);
            if (((n - k) & 1) != 0) mpfr_neg(binr, binr, MPFR_RNDN);
            mpfr_fma(accr, binr, fr, accr, MPFR_RNDN);
            mpfr_fma(acci, binr, fi, acci, MPFR_RNDN);
        }
        if (!ok) break;
        /* sn = step^n */
        mpfr_set_ui(snr, 1, MPFR_RNDN);
        mpfr_set_ui(sni, 0, MPFR_RNDN);
        for (int q = 0; q < n; q++) nd_cmul(snr, sni, snr, sni, str, sti);
        /* T(i,0) = acc / sn */
        mpfr_complex_div(Tr[(size_t)i * terms + 0], Ti[(size_t)i * terms + 0],
                         accr, acci, snr, sni);
    }

    if (ok) {
        for (int j = 1; j < terms; j++) {
            /* denom = 2^j - 1 */
            mpfr_set_ui(denom, 1, MPFR_RNDN);
            mpfr_mul_2ui(denom, denom, (unsigned long)j, MPFR_RNDN);
            mpfr_sub_ui(denom, denom, 1, MPFR_RNDN);
            for (int i = j; i < terms; i++) {
                size_t cur = (size_t)i * terms + j;
                size_t a = (size_t)i * terms + (j - 1);
                size_t b = (size_t)(i - 1) * terms + (j - 1);
                /* T(i,j) = T(a) + (T(a) - T(b)) / denom */
                mpfr_sub(tr, Tr[a], Tr[b], MPFR_RNDN);
                mpfr_sub(ti, Ti[a], Ti[b], MPFR_RNDN);
                mpfr_div(tr, tr, denom, MPFR_RNDN);
                mpfr_div(ti, ti, denom, MPFR_RNDN);
                mpfr_add(Tr[cur], Tr[a], tr, MPFR_RNDN);
                mpfr_add(Ti[cur], Ti[a], ti, MPFR_RNDN);
            }
        }
    }

    Expr* result = NULL;
    if (ok) {
        size_t last = (size_t)(terms - 1) * terms + (terms - 1);
        result = nd_from_complex_mpfr(Tr[last], Ti[last]);
    }

    nd_bind_restore(&bind);
    for (size_t c = 0; c < cells; c++) { mpfr_clear(Tr[c]); mpfr_clear(Ti[c]); }
    free(Tr); free(Ti);
    mpfr_clears(x0r, x0i, sr, si, hr, str, sti, snr, sni, accr, acci,
                ptr, pti, fr, fi, binr, tr, ti, denom, (mpfr_ptr)0);

    if (!ok) {
        nd_warn("nnum", "expr did not evaluate to a number at a sample point");
        return NULL;
    }
    return result;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Small expression builders for the NIntegrate path                  *
 * ------------------------------------------------------------------ */

static Expr* mk1(const char* head, Expr* a) {
    Expr** v = malloc(sizeof(Expr*));
    v[0] = a;
    Expr* e = expr_new_function(expr_new_symbol(head), v, 1);
    free(v);
    return e;
}

static Expr* mk2(const char* head, Expr* a, Expr* b) {
    Expr** v = malloc(sizeof(Expr*) * 2);
    v[0] = a; v[1] = b;
    Expr* e = expr_new_function(expr_new_symbol(head), v, 2);
    free(v);
    return e;
}

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* method;        /* SYM_EulerSum (default) or SYM_NIntegrate */
    Expr*       scale;         /* borrowed Scale value (NULL => 1)         */
    int         terms;         /* EulerSum tableau depth (default 7)       */
    bool        prec_mpfr;     /* WorkingPrecision selects MPFR            */
    long        bits;          /* MPFR working precision in bits           */
    Expr*       wp_val;        /* borrowed WorkingPrecision (passthrough)  */
    Expr*       pg_val;        /* borrowed PrecisionGoal (passthrough)     */
    Expr*       mr_val;        /* borrowed MaxRecursion (passthrough)      */
} NdOpts;

static bool nd_is_known_option(const char* s) {
    return s == SYM_Method
        || s == SYM_Scale
        || s == SYM_Terms
        || s == SYM_WorkingPrecision
        || s == SYM_PrecisionGoal
        || s == SYM_MaxRecursion;
}

static bool nd_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nd_is_known_option(lhs->data.symbol.name);
}

static bool nd_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!nd_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool nd_apply_option(Expr* rule, NdOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol.name == SYM_EulerSum
                || rhs->data.symbol.name == SYM_NIntegrate)) {
            o->method = rhs->data.symbol.name; return true;
        }
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
            o->method = SYM_EulerSum; return true;
        }
        nd_warn("badmeth", "Method must be EulerSum or NIntegrate; using EulerSum");
        o->method = SYM_EulerSum;
        return true;
    }
    if (name == SYM_Scale) { o->scale = rhs; return true; }
    if (name == SYM_Terms) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 1) {
            o->terms = (int)rhs->data.integer;
            if (o->terms > ND_MAX_TERMS) o->terms = ND_MAX_TERMS;
            return true;
        }
        nd_warn("badopt", "Terms must be a positive integer");
        return false;
    }
    if (name == SYM_WorkingPrecision) {
        if (!nd_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            nd_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        o->wp_val = rhs;
        return true;
    }
    if (name == SYM_PrecisionGoal) { o->pg_val = rhs; return true; }
    if (name == SYM_MaxRecursion) {
        if (rhs->type != EXPR_INTEGER || rhs->data.integer < 0) {
            nd_warn("badopt", "MaxRecursion must be a non-negative integer");
            return false;
        }
        o->mr_val = rhs;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  NIntegrate — Gamma(n+1) * NResidue[expr/(x-x0)^(n+1), {x,x0}, ...]  *
 * ------------------------------------------------------------------ */

static Expr* nd_nintegrate(Expr* expr, const char* var, Expr* n_expr,
                           Expr* x0_expr, NdOpts* o) {
    /* integrand = expr * (x - x0)^(-(n+1)) */
    Expr* base = mk2("Plus", expr_new_symbol(var),
                     mk2("Times", expr_new_integer(-1), expr_copy(x0_expr)));
    Expr* expo = mk2("Times", expr_new_integer(-1),
                     mk2("Plus", expr_copy(n_expr), expr_new_integer(1)));
    Expr* integrand = mk2("Times", expr_copy(expr), mk2("Power", base, expo));

    Expr* speclist = mk2("List", expr_new_symbol(var), expr_copy(x0_expr));
    Expr* radius = o->scale
        ? mk2("Rule", expr_new_symbol(SYM_Radius), mk1("Abs", expr_copy(o->scale)))
        : mk2("Rule", expr_new_symbol(SYM_Radius), expr_new_integer(1));

    int extra = (o->wp_val ? 1 : 0) + (o->pg_val ? 1 : 0) + (o->mr_val ? 1 : 0);
    int cnt = 3 + extra;
    Expr** a = malloc(sizeof(Expr*) * (size_t)cnt);
    a[0] = integrand; a[1] = speclist; a[2] = radius;
    int idx = 3;
    if (o->wp_val)
        a[idx++] = mk2("Rule", expr_new_symbol(SYM_WorkingPrecision), expr_copy(o->wp_val));
    if (o->pg_val)
        a[idx++] = mk2("Rule", expr_new_symbol(SYM_PrecisionGoal), expr_copy(o->pg_val));
    if (o->mr_val)
        a[idx++] = mk2("Rule", expr_new_symbol(SYM_MaxRecursion), expr_copy(o->mr_val));
    Expr* nres = expr_new_function(expr_new_symbol(SYM_NResidue), a, cnt);
    free(a);

    Expr* gam = mk1("Gamma", mk2("Plus", expr_copy(n_expr), expr_new_integer(1)));
    Expr* full = mk2("Times", gam, nres);

    Expr* raw = eval_and_free(full);
    if (!raw) return NULL;

    NumericSpec spec;
#ifdef USE_MPFR
    if (o->prec_mpfr) { spec.mode = NUMERIC_MODE_MPFR; spec.bits = o->bits;
                        spec.preserve_inexact = false; }
    else              spec = numeric_machine_spec();
#else
    spec = numeric_machine_spec();
#endif
    Expr* num = numericalize(raw, spec);
    expr_free(raw);
    return num;
}

/* ------------------------------------------------------------------ *
 *  Manual list-threading over arg 0                                   *
 * ------------------------------------------------------------------ */

static Expr* nd_thread_over_list(Expr* res) {
    Expr* lst = res->data.function.args[0];
    size_t n = lst->data.function.arg_count;
    size_t argc = res->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr** ca = malloc(sizeof(Expr*) * argc);
        ca[0] = expr_copy(lst->data.function.args[i]);
        for (size_t j = 1; j < argc; j++) ca[j] = expr_copy(res->data.function.args[j]);
        items[i] = expr_new_function(expr_new_symbol(SYM_ND), ca, argc);
        free(ca);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return eval_and_free(out);
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nd(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;     /* expr, {x[,n]}, x0  (+ options) */

    /* Manual threading over a List in arg 0 (ND is intentionally NOT
     * ATTR_LISTABLE: that would split the {x, n} spec). */
    Expr* arg0 = res->data.function.args[0];
    if (arg0->type == EXPR_FUNCTION
        && arg0->data.function.head->type == EXPR_SYMBOL
        && arg0->data.function.head->data.symbol.name == SYM_List) {
        return nd_thread_over_list(res);
    }

    /* Peel trailing options; exactly three positional args must remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && nd_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_is_option_arg(res->data.function.args[i])) {
            nd_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 3) return NULL;

    NdOpts o;
    o.method = SYM_EulerSum;
    o.scale = NULL;
    o.terms = 7;
    o.prec_mpfr = false;
    o.bits = 0;
    o.wp_val = o.pg_val = o.mr_val = NULL;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_apply_option(res->data.function.args[i], &o)) return NULL;
    }

    /* Parse the {x, n} spec (or a bare symbol x => n = 1). */
    Expr* spec = res->data.function.args[1];
    const char* var = NULL;
    Expr* n_expr = NULL;
    Expr* n_owned = NULL;
    if (spec->type == EXPR_SYMBOL) {
        var = spec->data.symbol.name;
        n_owned = expr_new_integer(1);
        n_expr = n_owned;
    } else if (spec->type == EXPR_FUNCTION
               && spec->data.function.head->type == EXPR_SYMBOL
               && spec->data.function.head->data.symbol.name == SYM_List
               && spec->data.function.arg_count == 2
               && spec->data.function.args[0]->type == EXPR_SYMBOL) {
        var = spec->data.function.args[0]->data.symbol.name;
        n_expr = spec->data.function.args[1];
    } else {
        nd_warn("ivar", "the second argument must be x or {x, n} with x a symbol");
        return NULL;
    }

    Expr* expr = res->data.function.args[0];
    Expr* x0_expr = res->data.function.args[2];
    Expr* result = NULL;

    if (o.method == SYM_NIntegrate) {
        result = nd_nintegrate(expr, var, n_expr, x0_expr, &o);
    } else {
        /* EulerSum needs a non-negative integer order; n == 0 is just expr(x0). */
        if (n_expr->type != EXPR_INTEGER || n_expr->data.integer < 0) {
            nd_warn("ord", "Method -> EulerSum requires a non-negative integer "
                           "order; use Method -> NIntegrate for fractional order");
        } else {
            int n = (int)n_expr->data.integer;
#ifdef USE_MPFR
            if (o.prec_mpfr)
                result = nd_eulersum_mpfr(expr, var, n, x0_expr, o.scale, o.terms, o.bits);
            else
#endif
                result = nd_eulersum_machine(expr, var, n, x0_expr, o.scale, o.terms);
        }
    }

    expr_free(n_owned);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nd_init(void) {
    symtab_add_builtin("ND", builtin_nd);
    /* Protected only. List-threading over arg 0 is manual (see builtin_nd),
     * so ATTR_LISTABLE is deliberately NOT set — it would split the {x, n}
     * spec across bogus single-argument calls. */
    symtab_get_def("ND")->attributes |= ATTR_PROTECTED;
}

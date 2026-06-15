/*
 * hypergeopfq.c -- HypergeometricPFQ for Mathilda.
 *
 * HypergeometricPFQ[{a1,...,ap}, {b1,...,bq}, z] is the generalized
 * hypergeometric function pFq(a;b;z) with series
 *
 *     pFq(a;b;z) = Sum_{k>=0}  ( prod_i (a_i)_k / prod_j (b_j)_k ) z^k / k!
 *
 * where (x)_k is the Pochhammer (rising-factorial) symbol.  The evaluation
 * order in builtin_hypergeometric_pfq is:
 *
 *   1. validate (arity 3; first two args are Lists);
 *   2. thread over a List third argument;
 *   3. cancel parameters common to the upper and lower lists (generic case);
 *   4. z == 0 -> 1;
 *   5. terminate to an explicit polynomial when an upper parameter is a
 *      non-positive integer (exact, symbolic-z capable);
 *   6. numeric evaluation by direct series summation (machine / MPFR /
 *      complex) in the convergent regime;
 *   7. otherwise return NULL (stay unevaluated).
 *
 * The convergence regime is: p<=q (entire, all finite z); p==q+1 (|z|<1).
 * For p==q+1 with |z|>=1 and for p>q+1 (asymptotic only) the series is not
 * summed directly -- the call stays unevaluated.  Analytic continuation is a
 * deliberate follow-on (see docs/spec).
 *
 * Memory contract: builtin takes ownership of res but must not free it.
 */

#include "hypergeopfq.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include "arithmetic.h"   /* is_rational, is_complex, make_complex */
#include "numeric.h"      /* numeric_combined_bits, get_approx_mpfr, ... */
#ifdef USE_MPFR
#include "numeric_complex.h"
#include <mpfr.h>
#endif
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Hard caps on the number of series terms summed before giving up. */
#define HGPFQ_MACHINE_MAX_TERMS 200000
#define HGPFQ_MPFR_MAX_TERMS    1000000

/* ------------------------------------------------------------------ */
/*  Small structural helpers                                           */
/* ------------------------------------------------------------------ */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

/* True if e is a non-positive integer (..., -2, -1, 0); fills *neg with
 * its magnitude n>=0 (param == -n) when true. */
static bool is_nonpos_int(const Expr* e, long* neg) {
    if (e && e->type == EXPR_INTEGER && e->data.integer <= 0) {
        if (neg) *neg = (long)(-e->data.integer);
        return true;
    }
    return false;
}

/* Build HypergeometricPFQ[a, b, z] from owned argument copies and evaluate. */
static Expr* rebuild_eval(Expr* a, Expr* b, Expr* z) {
    Expr* args[3] = { a, b, z };
    Expr* call = expr_new_function(expr_new_symbol(SYM_HypergeometricPFQ), args, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Pochhammer product (x)_k = x (x+1) ... (x+k-1) as a symbolic Times,
 * or Integer 1 when k == 0.  `x` is copied. */
static Expr* poch_product(Expr* x, long k) {
    if (k <= 0) return expr_new_integer(1);
    Expr** fac = malloc(sizeof(Expr*) * (size_t)k);
    for (long i = 0; i < k; i++) {
        if (i == 0) {
            fac[i] = expr_copy(x);
        } else {
            Expr* pa[2] = { expr_copy(x), expr_new_integer(i) };
            fac[i] = expr_new_function(expr_new_symbol(SYM_Plus), pa, 2);
        }
    }
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times), fac, (size_t)k);
    free(fac);
    return t;
}

/* ------------------------------------------------------------------ */
/*  Parameter cancellation (generic): drop one element common to the   */
/*  upper and lower lists, unless it is a non-positive integer (which   */
/*  makes the series terminate and must be kept).                       */
/* ------------------------------------------------------------------ */

static Expr* try_cancel(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;
    for (size_t i = 0; i < p; i++) {
        Expr* ai = a->data.function.args[i];
        long dummy;
        if (is_nonpos_int(ai, &dummy)) continue;   /* keep terminating params */
        for (size_t j = 0; j < q; j++) {
            if (!expr_eq(ai, b->data.function.args[j])) continue;
            /* Build new lists with element i (upper) / j (lower) removed. */
            Expr** na = malloc(sizeof(Expr*) * (p - 1));
            Expr** nb = malloc(sizeof(Expr*) * (q - 1));
            size_t w = 0;
            for (size_t k = 0; k < p; k++)
                if (k != i) na[w++] = expr_copy(a->data.function.args[k]);
            w = 0;
            for (size_t k = 0; k < q; k++)
                if (k != j) nb[w++] = expr_copy(b->data.function.args[k]);
            Expr* nal = expr_new_function(expr_new_symbol(SYM_List), na, p - 1);
            Expr* nbl = expr_new_function(expr_new_symbol(SYM_List), nb, q - 1);
            free(na); free(nb);
            return rebuild_eval(nal, nbl, expr_copy(z));
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Termination to a polynomial                                        */
/* ------------------------------------------------------------------ */

/* If any upper parameter is a non-positive integer -n, the series truncates
 * at k = n (the smallest such n).  Returns the explicit Sum_{k=0}^{n} c_k z^k
 * evaluated, or NULL when termination does not apply or a lower non-positive
 * integer would force a 0 denominator inside the truncation range. */
static Expr* try_terminate(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;

    long n = -1;
    for (size_t i = 0; i < p; i++) {
        long m;
        if (is_nonpos_int(a->data.function.args[i], &m)) {
            if (n < 0 || m < n) n = m;
        }
    }
    if (n < 0) return NULL;   /* no terminating upper parameter */

    /* A lower non-positive integer -m gives (b_j)_k == 0 for k > m; if m < n
     * the denominator vanishes inside the range -> undefined.  Bail. */
    for (size_t j = 0; j < q; j++) {
        long m;
        if (is_nonpos_int(b->data.function.args[j], &m) && m < n) return NULL;
    }

    Expr** terms = malloc(sizeof(Expr*) * (size_t)(n + 1));
    for (long k = 0; k <= n; k++) {
        /* numerator: prod_i (a_i)_k */
        size_t nfac = p + 2;                 /* uppers + denom-power + z-power */
        Expr** mul = malloc(sizeof(Expr*) * nfac);
        size_t w = 0;
        for (size_t i = 0; i < p; i++)
            mul[w++] = poch_product(a->data.function.args[i], k);

        /* denominator: prod_j (b_j)_k * k!  -> Power[Times[...], -1] */
        size_t nden = q + 1;
        Expr** den = malloc(sizeof(Expr*) * nden);
        size_t dw = 0;
        for (size_t j = 0; j < q; j++)
            den[dw++] = poch_product(b->data.function.args[j], k);
        den[dw++] = expr_new_function(expr_new_symbol(SYM_Factorial),
                        (Expr*[]){ expr_new_integer(k) }, 1);
        Expr* denom = expr_new_function(expr_new_symbol(SYM_Times), den, nden);
        free(den);
        mul[w++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ denom, expr_new_integer(-1) }, 2);

        /* z^k */
        mul[w++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(z), expr_new_integer(k) }, 2);

        terms[k] = expr_new_function(expr_new_symbol(SYM_Times), mul, w);
        free(mul);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n + 1));
    free(terms);
    Expr* r = evaluate(sum);
    expr_free(sum);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Symbolic reductions to elementary functions                        */
/* ------------------------------------------------------------------ */

static bool is_int_val(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}
/* True if e == odd/2 (half-odd-integer); fills *twice with the numerator. */
static bool is_half_int(const Expr* e, int64_t* twice) {
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d) && d == 2) { if (twice) *twice = n; return true; }
    return false;
}

/* 1 - z as Plus[1, Times[-1, z]]. */
static Expr* one_minus(Expr* z) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(z) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_integer(1), neg }, 2);
}
/* Sqrt[z] as Power[z, 1/2]. */
static Expr* sqrt_of(Expr* z) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ expr_copy(z), make_rational(1, 2) }, 2);
}

/* Recognise a handful of closed forms that need no new function heads.
 * Returns the evaluated elementary expression, or NULL if no rule matches.
 * Valid as analytic continuations for all z, so safe to run before the
 * numeric series. */
static Expr* try_reduce(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;
    Expr** au = a->data.function.args;
    Expr** bl = b->data.function.args;
    Expr* out = NULL;

    if (p == 0 && q == 0) {
        /* 0F0(;;z) = E^z */
        out = expr_new_function(expr_new_symbol(SYM_Exp), (Expr*[]){ expr_copy(z) }, 1);
    } else if (p == 1 && q == 0) {
        /* 1F0(a;;z) = (1 - z)^(-a) */
        Expr* nega = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), expr_copy(au[0]) }, 2);
        out = expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ one_minus(z), nega }, 2);
    } else if (p == 0 && q == 1) {
        int64_t tw;
        if (is_half_int(bl[0], &tw) && tw == 1) {
            /* 0F1(;1/2;z) = Cosh[2 Sqrt[z]] */
            Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(2), sqrt_of(z) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Cosh), (Expr*[]){ arg }, 1);
        } else if (is_half_int(bl[0], &tw) && tw == 3) {
            /* 0F1(;3/2;z) = Sinh[2 Sqrt[z]] / (2 Sqrt[z]) */
            Expr* twosq = expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(2), sqrt_of(z) }, 2);
            Expr* sinh = expr_new_function(expr_new_symbol(SYM_Sinh),
                             (Expr*[]){ expr_copy(twosq) }, 1);
            Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ twosq, expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ sinh, inv }, 2);
        }
    } else if (p == 1 && q == 1) {
        /* 1F1(1;2;z) = (E^z - 1) / z */
        if (is_int_val(au[0], 1) && is_int_val(bl[0], 2)) {
            Expr* em1 = expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Exp),
                                           (Expr*[]){ expr_copy(z) }, 1),
                                       expr_new_integer(-1) }, 2);
            Expr* invz = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(z), expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ em1, invz }, 2);
        }
    } else if (p == 2 && q == 1) {
        /* 2F1(1,1;2;z) = -Log[1 - z] / z */
        if (is_int_val(au[0], 1) && is_int_val(au[1], 1) && is_int_val(bl[0], 2)) {
            Expr* lg = expr_new_function(expr_new_symbol(SYM_Log),
                           (Expr*[]){ one_minus(z) }, 1);
            Expr* invz = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(z), expr_new_integer(-1) }, 2);
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-1), lg, invz }, 3);
        }
    }

    if (!out) return NULL;
    Expr* r = evaluate(out);
    expr_free(out);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Numeric: machine double-complex                                    */
/* ------------------------------------------------------------------ */

/* Machine-precision numeric approximation of e into *out; sets *inexact if
 * any contributing leaf was inexact (Real / MPFR).  Returns false if e is
 * not a concrete number. */
static bool approx_machine(const Expr* e, double complex* out, bool* inexact) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; if (inexact) *inexact = true; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
                                   if (inexact) *inexact = true; return true; }
#endif
    int64_t nn, dd;
    if (is_rational((Expr*)e, &nn, &dd)) { *out = (double)nn / (double)dd; return true; }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        double complex r, i;
        if (approx_machine(re, &r, inexact) && approx_machine(im, &i, inexact)) {
            *out = creal(r) + creal(i) * I;
            return true;
        }
    }
    return false;
}

/* True if z is (near) a non-positive integer -- a lower-parameter pole. */
static bool near_nonpos_int(double complex z) {
    if (fabs(cimag(z)) > 1e-9) return false;
    double r = creal(z);
    return r <= 1e-9 && fabs(r - round(r)) < 1e-9;
}

static Expr* machine_sum(double complex zc,
                         const double complex* ac, size_t p,
                         const double complex* bc, size_t q,
                         bool all_real) {
    double complex term = 1.0, sum = 1.0;
    int settled = 0;
    for (long k = 0; k < HGPFQ_MACHINE_MAX_TERMS; k++) {
        double complex fac = zc / (double)(k + 1);
        for (size_t i = 0; i < p; i++) fac *= (ac[i] + (double)k);
        for (size_t j = 0; j < q; j++) {
            double complex d = bc[j] + (double)k;
            if (cabs(d) < 1e-300) return NULL;   /* denominator pole */
            fac /= d;
        }
        term *= fac;
        sum += term;
        if (cabs(term) <= cabs(sum) * 1e-16) {
            if (++settled >= 2) break;
        } else {
            settled = 0;
        }
        if (!isfinite(creal(sum)) || !isfinite(cimag(sum))) return NULL;
    }
    if (all_real && cimag(sum) == 0.0) return expr_new_real(creal(sum));
    return make_complex(expr_new_real(creal(sum)), expr_new_real(cimag(sum)));
}

/* ------------------------------------------------------------------ */
/*  Numeric: MPFR (real & complex via (re,im) pairs)                   */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
typedef struct { mpfr_t re, im; } cpx_t;

static void cpx_init(cpx_t* c, mpfr_prec_t bits) { mpfr_init2(c->re, bits); mpfr_init2(c->im, bits); }
static void cpx_clear(cpx_t* c) { mpfr_clear(c->re); mpfr_clear(c->im); }
static void cpx_set(cpx_t* d, const cpx_t* s) {
    mpfr_set(d->re, s->re, MPFR_RNDN); mpfr_set(d->im, s->im, MPFR_RNDN);
}
static void cpx_add(cpx_t* o, const cpx_t* x, const cpx_t* y) {
    mpfr_add(o->re, x->re, y->re, MPFR_RNDN);
    mpfr_add(o->im, x->im, y->im, MPFR_RNDN);
}
/* o = x * y  (o must not alias x or y) */
static void cpx_mul(cpx_t* o, const cpx_t* x, const cpx_t* y, mpfr_prec_t bits) {
    mpfr_t t1, t2;
    mpfr_init2(t1, bits); mpfr_init2(t2, bits);
    mpfr_mul(t1, x->re, y->re, MPFR_RNDN);
    mpfr_mul(t2, x->im, y->im, MPFR_RNDN);
    mpfr_sub(t1, t1, t2, MPFR_RNDN);          /* re = xr*yr - xi*yi */
    mpfr_mul(o->im, x->re, y->im, MPFR_RNDN);
    mpfr_mul(t2, x->im, y->re, MPFR_RNDN);
    mpfr_add(o->im, o->im, t2, MPFR_RNDN);    /* im = xr*yi + xi*yr */
    mpfr_set(o->re, t1, MPFR_RNDN);
    mpfr_clear(t1); mpfr_clear(t2);
}
/* o = x / y */
static void cpx_div(cpx_t* o, const cpx_t* x, const cpx_t* y) {
    mpfr_complex_div(o->re, o->im, x->re, x->im, y->re, y->im);
}
/* magnitude proxy max(|re|,|im|) into m */
static void cpx_absmax(mpfr_t m, const cpx_t* x) {
    mpfr_abs(m, x->re, MPFR_RNDN);
    mpfr_t t; mpfr_init2(t, mpfr_get_prec(m));
    mpfr_abs(t, x->im, MPFR_RNDN);
    if (mpfr_cmp(t, m) > 0) mpfr_set(m, t, MPFR_RNDN);
    mpfr_clear(t);
}

/* Extract e into a cpx pair at `bits` precision; false if not numeric. */
static bool cpx_from_expr(const Expr* e, cpx_t* c, mpfr_prec_t bits) {
    cpx_init(c, bits);
    bool inexact = false;
    if (!get_approx_mpfr(e, c->re, c->im, &inexact)) { cpx_clear(c); return false; }
    return true;
}

static Expr* hgpfq_mpfr_sum(const Expr* z, Expr* a, Expr* b,
                            size_t p, size_t q, long target_bits) {
    mpfr_prec_t work = (mpfr_prec_t)(target_bits + 64);

    cpx_t zc, *ac = NULL, *bc = NULL;
    if (!cpx_from_expr(z, &zc, work)) return NULL;
    if (p) ac = malloc(sizeof(cpx_t) * p);
    if (q) bc = malloc(sizeof(cpx_t) * q);
    bool ok = true;
    for (size_t i = 0; i < p && ok; i++)
        ok = cpx_from_expr(a->data.function.args[i], &ac[i], work);
    for (size_t j = 0; j < q && ok; j++)
        ok = cpx_from_expr(b->data.function.args[j], &bc[j], work);
    if (!ok) {
        cpx_clear(&zc);
        free(ac); free(bc);
        return NULL;
    }

    cpx_t term, sum, fac, tmp, num;
    cpx_init(&term, work); cpx_init(&sum, work);
    cpx_init(&fac, work);  cpx_init(&tmp, work); cpx_init(&num, work);
    mpfr_set_ui(term.re, 1, MPFR_RNDN); mpfr_set_ui(term.im, 0, MPFR_RNDN);
    mpfr_set_ui(sum.re, 1, MPFR_RNDN);  mpfr_set_ui(sum.im, 0, MPFR_RNDN);

    mpfr_t tmag, smag, thr;
    mpfr_init2(tmag, work); mpfr_init2(smag, work); mpfr_init2(thr, work);

    bool failed = false;
    int settled = 0;
    long k = 0;
    for (; k < HGPFQ_MPFR_MAX_TERMS; k++) {
        /* fac = z / (k+1) * prod(a_i + k) / prod(b_j + k) */
        cpx_set(&fac, &zc);
        mpfr_div_ui(fac.re, fac.re, (unsigned long)(k + 1), MPFR_RNDN);
        mpfr_div_ui(fac.im, fac.im, (unsigned long)(k + 1), MPFR_RNDN);
        for (size_t i = 0; i < p; i++) {
            cpx_set(&tmp, &ac[i]);
            mpfr_add_ui(tmp.re, tmp.re, (unsigned long)k, MPFR_RNDN);
            cpx_mul(&num, &fac, &tmp, work);
            cpx_set(&fac, &num);
        }
        for (size_t j = 0; j < q; j++) {
            cpx_set(&tmp, &bc[j]);
            mpfr_add_ui(tmp.re, tmp.re, (unsigned long)k, MPFR_RNDN);
            cpx_div(&num, &fac, &tmp);
            cpx_set(&fac, &num);
            if (mpfr_nan_p(fac.re) || mpfr_inf_p(fac.re) ||
                mpfr_nan_p(fac.im) || mpfr_inf_p(fac.im)) { failed = true; break; }
        }
        if (failed) break;
        cpx_mul(&num, &term, &fac, work);   /* term *= fac */
        cpx_set(&term, &num);
        cpx_add(&sum, &sum, &term);

        cpx_absmax(tmag, &term);
        cpx_absmax(smag, &sum);
        mpfr_mul_2si(thr, smag, -(long)work, MPFR_RNDN);   /* |sum| * 2^-work */
        if (mpfr_cmp(tmag, thr) <= 0) { if (++settled >= 2) break; }
        else settled = 0;
    }

    Expr* result = NULL;
    if (!failed && k < HGPFQ_MPFR_MAX_TERMS) {
        mpfr_t orr, oii;
        mpfr_init2(orr, (mpfr_prec_t)target_bits);
        mpfr_init2(oii, (mpfr_prec_t)target_bits);
        mpfr_set(orr, sum.re, MPFR_RNDN);
        mpfr_set(oii, sum.im, MPFR_RNDN);
        result = numeric_mpfr_make_complex(orr, oii);
        mpfr_clear(orr); mpfr_clear(oii);
    }

    mpfr_clear(tmag); mpfr_clear(smag); mpfr_clear(thr);
    cpx_clear(&term); cpx_clear(&sum); cpx_clear(&fac); cpx_clear(&tmp); cpx_clear(&num);
    cpx_clear(&zc);
    for (size_t i = 0; i < p; i++) cpx_clear(&ac[i]);
    for (size_t j = 0; j < q; j++) cpx_clear(&bc[j]);
    free(ac); free(bc);
    return result;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/*  Numeric dispatch                                                   */
/* ------------------------------------------------------------------ */

static Expr* try_numeric(Expr* a, Expr* b, Expr* z) {
    size_t p = a->data.function.arg_count;
    size_t q = b->data.function.arg_count;

    /* Machine extraction of z and all parameters: gating + machine path. */
    double complex zc; bool inexact = false;
    if (!approx_machine(z, &zc, &inexact)) return NULL;

    double complex* ac = p ? malloc(sizeof(double complex) * p) : NULL;
    double complex* bc = q ? malloc(sizeof(double complex) * q) : NULL;
    bool all_real = (cimag(zc) == 0.0);
    bool ok = true;
    for (size_t i = 0; i < p && ok; i++) {
        ok = approx_machine(a->data.function.args[i], &ac[i], &inexact);
        if (ok && cimag(ac[i]) != 0.0) all_real = false;
    }
    for (size_t j = 0; j < q && ok; j++) {
        ok = approx_machine(b->data.function.args[j], &bc[j], &inexact);
        if (ok && cimag(bc[j]) != 0.0) all_real = false;
        if (ok && near_nonpos_int(bc[j])) ok = false;   /* lower-parameter pole */
    }
    if (!ok || !inexact) { free(ac); free(bc); return NULL; }

    /* Convergence gate (radius of convergence by p vs q). */
    double az = cabs(zc);
    if (p > q + 1) { free(ac); free(bc); return NULL; }        /* divergent */
    if (p == q + 1 && az >= 1.0 - 1e-12) { free(ac); free(bc); return NULL; }

    Expr* result = NULL;
#ifdef USE_MPFR
    /* Precision contagion: the governing precision is the MINIMUM among the
     * inexact (Real/MPFR) leaves of z and every parameter. A machine-precision
     * argument forces the machine path even when others are high-precision
     * MPFR (matching Mathematica). Only when every inexact argument exceeds
     * machine precision do we take the arbitrary-precision MPFR sum. */
    {
        long gov = numeric_min_inexact_bits(z);
        for (size_t i = 0; i < p; i++) {
            long bi = numeric_min_inexact_bits(a->data.function.args[i]);
            if (bi > 0 && (gov == 0 || bi < gov)) gov = bi;
        }
        for (size_t j = 0; j < q; j++) {
            long bj = numeric_min_inexact_bits(b->data.function.args[j]);
            if (bj > 0 && (gov == 0 || bj < gov)) gov = bj;
        }
        if (gov > 53) {
            result = hgpfq_mpfr_sum(z, a, b, p, q, gov);
            free(ac); free(bc);
            return result;
        }
    }
#endif
    result = machine_sum(zc, ac, p, bc, q, all_real);
    free(ac); free(bc);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_hypergeometric_pfq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* z = res->data.function.args[2];
    if (!is_list(a) || !is_list(b)) return NULL;

    /* (2) Thread over a List third argument. */
    if (is_list(z)) {
        size_t n = z->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++)
            out[i] = rebuild_eval(expr_copy(a), expr_copy(b),
                                  expr_copy(z->data.function.args[i]));
        Expr* lst = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        return lst;
    }

    /* (3) Generic parameter cancellation. */
    Expr* c = try_cancel(a, b, z);
    if (c) return c;

    /* (4) z == 0 -> 1. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) return expr_new_integer(1);

    /* (5) Termination to a polynomial. */
    Expr* t = try_terminate(a, b, z);
    if (t) return t;

    /* (6) Reduction to elementary closed forms (analytic continuations). */
    Expr* red = try_reduce(a, b, z);
    if (red) return red;

    /* (7) Numeric evaluation. */
    Expr* n = try_numeric(a, b, z);
    if (n) return n;

    /* (8) Stay unevaluated. */
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Convenience heads: rewrite to HypergeometricPFQ                    */
/* ------------------------------------------------------------------ */

/* Wrap a scalar parameter into a one-element List (copying). */
static Expr* list1(Expr* x) {
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_copy(x) }, 1);
}

Expr* builtin_hypergeometric_0f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* b = res->data.function.args[0];
    Expr* z = res->data.function.args[1];
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    return rebuild_eval(empty, list1(b), expr_copy(z));
}

Expr* builtin_hypergeometric_1f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* z = res->data.function.args[2];
    return rebuild_eval(list1(a), list1(b), expr_copy(z));
}

Expr* builtin_hypergeometric_2f1(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* cc = res->data.function.args[2];
    Expr* z = res->data.function.args[3];
    Expr* up = expr_new_function(expr_new_symbol(SYM_List),
                   (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
    return rebuild_eval(up, list1(cc), expr_copy(z));
}

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void hypergeopfq_init(void) {
    symtab_add_builtin("HypergeometricPFQ", builtin_hypergeometric_pfq);
    symtab_get_def("HypergeometricPFQ")->attributes |=
        ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    symtab_add_builtin("Hypergeometric0F1", builtin_hypergeometric_0f1);
    symtab_get_def("Hypergeometric0F1")->attributes |=
        ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    symtab_add_builtin("Hypergeometric1F1", builtin_hypergeometric_1f1);
    symtab_get_def("Hypergeometric1F1")->attributes |=
        ATTR_NUMERICFUNCTION | ATTR_PROTECTED;

    symtab_add_builtin("Hypergeometric2F1", builtin_hypergeometric_2f1);
    symtab_get_def("Hypergeometric2F1")->attributes |=
        ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
}

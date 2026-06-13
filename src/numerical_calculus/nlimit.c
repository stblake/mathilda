/*
 * nlimit.c — NLimit[expr, z -> z0, opts]
 *
 * Numerically estimates  lim_{z -> z0} expr  by sampling expr on a geometric
 * sequence of points approaching z0 and accelerating the resulting sequence.
 *
 *   Sample points
 *   -------------
 *   Finite z0:    z_k = z0 - d * Scale * 2^-k          (k = 0 .. Terms-1)
 *                 d is the (unit) Direction vector; the points sit on the -d
 *                 side of z0, i.e. one moves *along* d to reach z0.  The
 *                 default Direction -> Automatic == -1 places the points at
 *                 z0 + Scale*2^-k, approaching "from larger values".
 *   Infinite z0:  z_k = u * Scale * 2^k                (k = 0 .. Terms-1)
 *                 u is the unit ray of the infinite limit point (the direction
 *                 of Infinity / I Infinity / DirectedInfinity[..]); the points
 *                 march outward to infinity along that ray.
 *
 *   Acceleration
 *   ------------
 *   Method -> EulerSum (default): Richardson / Romberg extrapolation of the
 *     sequence S_k, treated as a function of the geometric step.  Following the
 *     same convention as ND's "EulerSum" (see nderiv.c) the tableau uses the
 *     all-powers denominator 2^j - 1:
 *         T(i,0) = S_i,
 *         T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
 *         result = T(Terms-1, Terms-1).
 *
 *   Method -> SequenceLimit: Wynn's epsilon algorithm (iterated Shanks /
 *     Aitken).  With ε_{-1}=0, ε_0 = S_n,
 *         ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1/(ε_k^{(n+1)} - ε_k^{(n)});
 *     the even columns ε_{2d} are the limit estimates.  WynnDegree -> d selects
 *     column 2d and requires at least 2(d+1) terms for a convergence check.
 *
 *   Robustness
 *   ----------
 *   The last two extrapolates are compared; if they fail to agree to a loose
 *   tolerance relative to the result magnitude (or the result is non-finite),
 *   NLimit::noise is emitted and the form is returned unevaluated.  As with
 *   Mathematica, spurious tiny residuals are *not* recognised as zero — apply
 *   Chop when needed.
 *
 * Options: Method (EulerSum | SequenceLimit), WorkingPrecision (MachinePrecision
 * | digits -> MPFR), Direction (Automatic == -1, or a complex approach vector),
 * Scale (initial step / distance, default 1), Terms (default 7), WynnDegree
 * (default 1).
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.  Every temporary OwnValue
 * created for the sampler is removed on all return paths.
 */

#include "nlimit.h"

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
#  include "numeric_complex.h"   /* numeric_mpfr_make_complex, mpfr_complex_div */
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* Hard cap on Terms: 2^(Terms-1) must stay an exact double for the machine
 * Richardson denominator, and any larger tableau is numerical nonsense. */
#define NL_MAX_TERMS 50

/* Convergence / noise gate.  Measured against the *sample scale* (the largest
 * |S_k|), not the result, so that genuine near-zero limits are not mistaken for
 * noise.  The gate is deliberately lenient: its job is to reject clearly
 * divergent or non-convergent sequences (e.g. 1/x as x->0, whose extrapolant
 * escapes far beyond the sampled data), not to certify accuracy.  A bounded but
 * only roughly-resolved limit (e.g. Sin[x]/x as x->Infinity, ~0) is accepted —
 * matching NLimit, which returns such values rather than refusing.
 *
 *   - reject if the result is non-finite;
 *   - reject if |result| exceeds NL_DIVERGE_FACTOR * (max |S_k|)  (divergence);
 *   - reject if the last two extrapolates differ by more than
 *     NL_NOISE_RELTOL * (max |S_k|)                               (no settling). */
#define NL_DIVERGE_FACTOR 4.0
#define NL_NOISE_RELTOL   0.5
#define NL_NOISE_ABSTOL   1e-6

static bool nl_accept(double result_mag, double step, double maxsample) {
    if (!isfinite(result_mag) || !isfinite(step)) return false;
    if (maxsample < NL_NOISE_ABSTOL) maxsample = NL_NOISE_ABSTOL;
    if (result_mag > NL_DIVERGE_FACTOR * maxsample) return false;
    if (step > NL_NOISE_RELTOL * maxsample + NL_NOISE_ABSTOL) return false;
    return true;
}

static void nl_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NLimit::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers (mirrors nderiv.c / nresidue.c)   *
 * ------------------------------------------------------------------ */

static bool nl_to_double_real(Expr* e, double* out) {
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

static bool nl_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (nl_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (nl_to_double_real(re, &r) && nl_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* nl_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
static Expr* nl_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}

/* (out_re, out_im) = (a_re + a_im i)(b_re + b_im i). Output may alias inputs. */
static void nl_cmul(mpfr_t out_re, mpfr_t out_im,
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
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors nderiv.c)                    *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} NlBind;

static void nl_bind_snapshot(NlBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void nl_bind_clear_temp(SymbolDef* def) {
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

static void nl_bind_set(NlBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    nl_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

static void nl_bind_restore(NlBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    nl_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Sampler: evaluate expr with the variable bound to a complex point  *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       f;     /* borrowed expr (free in the variable) */
    NlBind*     bind;
    NumericSpec spec;
} NlCtx;

/* Evaluate f after binding the variable to `value` (consumed). Returns the
 * numericalised result (caller frees) or NULL. */
static Expr* nl_eval_at(NlCtx* c, Expr* value) {
    nl_bind_set(c->bind, value);   /* takes a copy internally */
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

static bool nl_sample_machine(NlCtx* c, double _Complex z, double _Complex* out) {
    Expr* num = nl_eval_at(c, nl_from_complex_d(z));
    if (!num) return false;
    bool ok = nl_to_complex(num, out);
    if (ok && (!isfinite(creal(*out)) || !isfinite(cimag(*out)))) ok = false;
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nl_sample_mpfr(NlCtx* c, const mpfr_t z_re, const mpfr_t z_im,
                           mpfr_t out_re, mpfr_t out_im) {
    Expr* num = nl_eval_at(c, nl_from_complex_mpfr(z_re, z_im));
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}
#endif

/* Numericalise a constant (var-free) expression to a machine complex value. */
static bool nl_const_complex_machine(Expr* e, double _Complex* out) {
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
    expr_free(raw);
    bool ok = num && nl_to_complex(num, out);
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nl_const_complex_mpfr(Expr* e, long bits, mpfr_t re, mpfr_t im) {
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = bits;
    spec.preserve_inexact = false;
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, spec) : NULL;
    expr_free(raw);
    bool inexact;
    bool ok = num && get_approx_mpfr(num, re, im, &inexact);
    expr_free(num);
    return ok;
}
#endif

/* ------------------------------------------------------------------ *
 *  Limit-point classification: finite vs. an infinite ray            *
 * ------------------------------------------------------------------ */

/* If z0 denotes an infinite limit point, return true and store in *ray an
 * owned expression for the (un-normalised) direction of the ray from the
 * origin.  Otherwise return false and leave *ray untouched.
 *
 * Recognised: Infinity (+1), ComplexInfinity / DirectedInfinity[] (+1),
 * DirectedInfinity[d] (d), and any Times[...] containing an Infinity or a
 * DirectedInfinity factor (e.g. -Infinity, I Infinity) whose ray is the
 * product of the remaining factors with the directed-infinity direction. */
static bool nl_infinite_ray(Expr* z0, Expr** ray) {
    if (!z0) return false;
    if (z0->type == EXPR_SYMBOL) {
        if (z0->data.symbol == SYM_Infinity
            || z0->data.symbol == SYM_ComplexInfinity) {
            *ray = expr_new_integer(1);
            return true;
        }
        return false;
    }
    if (z0->type != EXPR_FUNCTION) return false;
    Expr* head = z0->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;

    if (head->data.symbol == SYM_DirectedInfinity) {
        if (z0->data.function.arg_count == 1) {
            *ray = expr_copy(z0->data.function.args[0]);
        } else {
            *ray = expr_new_integer(1);
        }
        return true;
    }

    if (head->data.symbol == SYM_Times) {
        size_t n = z0->data.function.arg_count;
        bool found = false;
        Expr* prod = NULL;          /* accumulating Times of the remaining factors */
        for (size_t i = 0; i < n; i++) {
            Expr* f = z0->data.function.args[i];
            Expr* contrib = NULL;
            if (f->type == EXPR_SYMBOL
                && (f->data.symbol == SYM_Infinity
                    || f->data.symbol == SYM_ComplexInfinity)) {
                found = true;
                continue;           /* Infinity contributes unit magnitude */
            }
            if (f->type == EXPR_FUNCTION
                && f->data.function.head->type == EXPR_SYMBOL
                && f->data.function.head->data.symbol == SYM_DirectedInfinity) {
                found = true;
                contrib = (f->data.function.arg_count == 1)
                            ? expr_copy(f->data.function.args[0])
                            : expr_new_integer(1);
            } else {
                contrib = expr_copy(f);
            }
            if (!prod) {
                prod = contrib;
            } else {
                Expr** v = malloc(sizeof(Expr*) * 2);
                v[0] = prod; v[1] = contrib;
                prod = expr_new_function(expr_new_symbol("Times"), v, 2);
                free(v);
            }
        }
        if (!found) { expr_free(prod); return false; }
        *ray = prod ? prod : expr_new_integer(1);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Extrapolation — machine precision                                  *
 * ------------------------------------------------------------------ */

/* Richardson / Romberg (the "EulerSum" method).  Returns the accelerated value
 * and, via *step, the magnitude of the last extrapolation step (the difference
 * between the last two diagonal entries) for the convergence test. */
static bool nl_richardson_machine(const double _Complex* S, int terms,
                                  double _Complex* result, double* step) {
    double _Complex* T = malloc(sizeof(double _Complex) * (size_t)terms * terms);
    if (!T) return false;
    for (int i = 0; i < terms; i++) T[(size_t)i * terms + 0] = S[i];
    for (int j = 1; j < terms; j++) {
        double denom = ldexp(1.0, j) - 1.0;        /* 2^j - 1 */
        for (int i = j; i < terms; i++) {
            double _Complex a = T[(size_t)i * terms + (j - 1)];
            double _Complex b = T[(size_t)(i - 1) * terms + (j - 1)];
            T[(size_t)i * terms + j] = a + (a - b) / denom;
        }
    }
    double _Complex last = T[(size_t)(terms - 1) * terms + (terms - 1)];
    double _Complex prev = (terms >= 2)
        ? T[(size_t)(terms - 2) * terms + (terms - 2)] : last;
    *result = last;
    *step = cabs(last - prev);
    free(T);
    return true;
}

/* Wynn's epsilon algorithm (the "SequenceLimit" method). Columns are stored
 * shifted by one so that column c holds ε_{c-1}; column 0 is ε_{-1} = 0,
 * column 1 is ε_0 = S. The degree-d estimate lives in column 2d+1. */
static bool nl_wynn_machine(const double _Complex* S, int terms, int degree,
                            double _Complex* result, double* step) {
    int maxdeg = (terms - 1) / 2;
    if (maxdeg < 1) return false;          /* need >= 3 terms for one Shanks */
    if (degree > maxdeg) degree = maxdeg;
    if (degree < 1) degree = 1;

    int cols = terms + 1;
    size_t stride = (size_t)cols;
    double _Complex* eps = malloc(sizeof(double _Complex) * stride * stride);
    if (!eps) return false;
    for (size_t i = 0; i < stride * stride; i++) eps[i] = 0.0;
#define NL_E(c, n) eps[(size_t)(c) * stride + (size_t)(n)]
    for (int n = 0; n < terms; n++) NL_E(1, n) = S[n];   /* ε_0 */

    for (int c = 2; c <= terms; c++) {
        int len = terms + 1 - c;
        for (int n = 0; n < len; n++) {
            double _Complex d = NL_E(c - 1, n + 1) - NL_E(c - 1, n);
            double _Complex recip;
            if (cabs(d) <= 0.0) recip = 1e300;          /* singular bridge */
            else                recip = 1.0 / d;
            NL_E(c, n) = NL_E(c - 2, n + 1) + recip;
        }
    }

    /* Result is read from the deepest even column ε_{2d} (column 2d+1 in the
     * shifted layout).  Rather than blindly taking the bottom corner — whose
     * tiny Shanks denominators amplify roundoff — pick the entry that best
     * agrees with its neighbour, i.e. where the algorithm has converged. */
    int rc = 2 * degree + 1;               /* result column (ε_{2d}) */
    int rl = terms + 1 - rc;               /* its length */
    double _Complex res = NL_E(rc, rl - 1);
    double best_d = INFINITY;
    bool found = false;
    for (int n = 1; n < rl; n++) {
        double _Complex a = NL_E(rc, n), b = NL_E(rc, n - 1);
        if (!isfinite(creal(a)) || !isfinite(cimag(a))) continue;
        if (!isfinite(creal(b)) || !isfinite(cimag(b))) continue;
        double d = cabs(a - b);
        if (d < best_d) { best_d = d; res = a; found = true; }
    }
    if (!found) {                           /* single-entry column fallback */
        res = NL_E(rc, rl - 1);
        double _Complex prev = (rc >= 3)
            ? NL_E(rc - 2, terms + 1 - (rc - 2) - 1) : NL_E(1, terms - 1);
        best_d = cabs(res - prev);
    }

    *result = res;
    *step = best_d;
#undef NL_E
    free(eps);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Extrapolation — MPFR precision                                     *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* L1 magnitude |re| + |im| of (re, im) as a double, for the gauges. */
static double nl_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}

static Expr* nl_richardson_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                                long bits, double maxsample, bool* converged) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    size_t cells = (size_t)terms * terms;
    mpfr_t* Tr = malloc(sizeof(mpfr_t) * cells);
    mpfr_t* Ti = malloc(sizeof(mpfr_t) * cells);
    for (size_t c = 0; c < cells; c++) { mpfr_init2(Tr[c], p); mpfr_init2(Ti[c], p); }
    mpfr_t tr, ti, denom;
    mpfr_inits2(p, tr, ti, denom, (mpfr_ptr)0);

    for (int i = 0; i < terms; i++) {
        mpfr_set(Tr[(size_t)i * terms + 0], Sr[i], MPFR_RNDN);
        mpfr_set(Ti[(size_t)i * terms + 0], Si[i], MPFR_RNDN);
    }
    for (int j = 1; j < terms; j++) {
        mpfr_set_ui(denom, 1, MPFR_RNDN);
        mpfr_mul_2ui(denom, denom, (unsigned long)j, MPFR_RNDN);
        mpfr_sub_ui(denom, denom, 1, MPFR_RNDN);
        for (int i = j; i < terms; i++) {
            size_t cur = (size_t)i * terms + j;
            size_t a = (size_t)i * terms + (j - 1);
            size_t b = (size_t)(i - 1) * terms + (j - 1);
            mpfr_sub(tr, Tr[a], Tr[b], MPFR_RNDN);
            mpfr_sub(ti, Ti[a], Ti[b], MPFR_RNDN);
            mpfr_div(tr, tr, denom, MPFR_RNDN);
            mpfr_div(ti, ti, denom, MPFR_RNDN);
            mpfr_add(Tr[cur], Tr[a], tr, MPFR_RNDN);
            mpfr_add(Ti[cur], Ti[a], ti, MPFR_RNDN);
        }
    }

    size_t last = (size_t)(terms - 1) * terms + (terms - 1);
    size_t prev = (terms >= 2)
        ? (size_t)(terms - 2) * terms + (terms - 2) : last;
    mpfr_sub(tr, Tr[last], Tr[prev], MPFR_RNDN);
    mpfr_sub(ti, Ti[last], Ti[prev], MPFR_RNDN);
    double gd = nl_l1_d(tr, ti);
    double gm = nl_l1_d(Tr[last], Ti[last]);
    bool finite = mpfr_number_p(Tr[last]) && mpfr_number_p(Ti[last]);
    *converged = finite && nl_accept(gm, gd, maxsample);

    Expr* out = nl_from_complex_mpfr(Tr[last], Ti[last]);
    for (size_t c = 0; c < cells; c++) { mpfr_clear(Tr[c]); mpfr_clear(Ti[c]); }
    free(Tr); free(Ti);
    mpfr_clears(tr, ti, denom, (mpfr_ptr)0);
    return out;
}

static Expr* nl_wynn_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                          int degree, long bits, double maxsample,
                          bool* converged) {
    int maxdeg = (terms - 1) / 2;
    if (maxdeg < 1) { *converged = false; return NULL; }
    if (degree > maxdeg) degree = maxdeg;
    if (degree < 1) degree = 1;

    mpfr_prec_t p = (mpfr_prec_t)bits;
    int cols = terms + 1;
    size_t stride = (size_t)cols;
    size_t cells = stride * stride;
    mpfr_t* er = malloc(sizeof(mpfr_t) * cells);
    mpfr_t* ei = malloc(sizeof(mpfr_t) * cells);
    for (size_t i = 0; i < cells; i++) {
        mpfr_init2(er[i], p); mpfr_init2(ei[i], p);
        mpfr_set_ui(er[i], 0, MPFR_RNDN); mpfr_set_ui(ei[i], 0, MPFR_RNDN);
    }
    mpfr_t dr, di, rr, ri, one, zero;
    mpfr_inits2(p, dr, di, rr, ri, one, zero, (mpfr_ptr)0);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_ui(zero, 0, MPFR_RNDN);
#define NL_ER(c, n) er[(size_t)(c) * stride + (size_t)(n)]
#define NL_EI(c, n) ei[(size_t)(c) * stride + (size_t)(n)]
    for (int n = 0; n < terms; n++) {
        mpfr_set(NL_ER(1, n), Sr[n], MPFR_RNDN);
        mpfr_set(NL_EI(1, n), Si[n], MPFR_RNDN);
    }
    for (int c = 2; c <= terms; c++) {
        int len = terms + 1 - c;
        for (int n = 0; n < len; n++) {
            mpfr_sub(dr, NL_ER(c - 1, n + 1), NL_ER(c - 1, n), MPFR_RNDN);
            mpfr_sub(di, NL_EI(c - 1, n + 1), NL_EI(c - 1, n), MPFR_RNDN);
            if (mpfr_zero_p(dr) && mpfr_zero_p(di)) {
                mpfr_set_d(rr, 1e300, MPFR_RNDN);
                mpfr_set_ui(ri, 0, MPFR_RNDN);
            } else {
                mpfr_complex_div(rr, ri, one, zero, dr, di);
            }
            mpfr_add(NL_ER(c, n), NL_ER(c - 2, n + 1), rr, MPFR_RNDN);
            mpfr_add(NL_EI(c, n), NL_EI(c - 2, n + 1), ri, MPFR_RNDN);
        }
    }

    /* Pick the ε_{2d} entry that best agrees with its neighbour (see the
     * machine path for the rationale). */
    int rc = 2 * degree + 1;
    int rl = terms + 1 - rc;
    int best_n = rl - 1;
    double best_d = -1.0;
    for (int n = 1; n < rl; n++) {
        if (!mpfr_number_p(NL_ER(rc, n)) || !mpfr_number_p(NL_EI(rc, n))) continue;
        if (!mpfr_number_p(NL_ER(rc, n - 1)) || !mpfr_number_p(NL_EI(rc, n - 1))) continue;
        mpfr_sub(dr, NL_ER(rc, n), NL_ER(rc, n - 1), MPFR_RNDN);
        mpfr_sub(di, NL_EI(rc, n), NL_EI(rc, n - 1), MPFR_RNDN);
        double d = nl_l1_d(dr, di);
        if (best_d < 0.0 || d < best_d) { best_d = d; best_n = n; }
    }
    if (best_d < 0.0) {                      /* single-entry fallback */
        int pc = (rc >= 3) ? rc - 2 : 1;
        int pn = (rc >= 3) ? terms + 1 - pc - 1 : terms - 1;
        mpfr_sub(dr, NL_ER(rc, rl - 1), NL_ER(pc, pn), MPFR_RNDN);
        mpfr_sub(di, NL_EI(rc, rl - 1), NL_EI(pc, pn), MPFR_RNDN);
        best_d = nl_l1_d(dr, di);
        best_n = rl - 1;
    }
    double gd = best_d;
    double gm = nl_l1_d(NL_ER(rc, best_n), NL_EI(rc, best_n));
    bool finite = mpfr_number_p(NL_ER(rc, best_n)) && mpfr_number_p(NL_EI(rc, best_n));
    *converged = finite && nl_accept(gm, gd, maxsample);

    Expr* out = nl_from_complex_mpfr(NL_ER(rc, best_n), NL_EI(rc, best_n));
#undef NL_ER
#undef NL_EI
    for (size_t i = 0; i < cells; i++) { mpfr_clear(er[i]); mpfr_clear(ei[i]); }
    free(er); free(ei);
    mpfr_clears(dr, di, rr, ri, one, zero, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* method;    /* SYM_EulerSum (default) or SYM_SequenceLimit */
    Expr*       direction; /* borrowed Direction value (NULL => Automatic) */
    Expr*       scale;     /* borrowed Scale value (NULL => 1)             */
    int         terms;     /* sample count / tableau depth (default 7)     */
    int         wynn;      /* WynnDegree (default 1)                       */
    bool        prec_mpfr; /* WorkingPrecision selects MPFR                */
    long        bits;      /* MPFR working precision in bits               */
} NlOpts;

static bool nl_is_known_option(const char* s) {
    return s == SYM_Method
        || s == SYM_Direction
        || s == SYM_Scale
        || s == SYM_Terms
        || s == SYM_WynnDegree
        || s == SYM_WorkingPrecision;
}

static bool nl_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nl_is_known_option(lhs->data.symbol);
}

static bool nl_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!nl_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool nl_apply_option(Expr* rule, NlOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol == SYM_EulerSum
                || rhs->data.symbol == SYM_SequenceLimit)) {
            o->method = rhs->data.symbol; return true;
        }
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
            o->method = SYM_EulerSum; return true;
        }
        nl_warn("badmeth", "Method must be EulerSum or SequenceLimit; using EulerSum");
        o->method = SYM_EulerSum;
        return true;
    }
    if (name == SYM_Direction) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)
            o->direction = NULL;       /* Automatic == -1 */
        else
            o->direction = rhs;
        return true;
    }
    if (name == SYM_Scale) { o->scale = rhs; return true; }
    if (name == SYM_Terms) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 2) {
            o->terms = (int)rhs->data.integer;
            if (o->terms > NL_MAX_TERMS) o->terms = NL_MAX_TERMS;
            return true;
        }
        nl_warn("badopt", "Terms must be an integer >= 2");
        return false;
    }
    if (name == SYM_WynnDegree) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 1) {
            o->wynn = (int)rhs->data.integer;
            return true;
        }
        nl_warn("badopt", "WynnDegree must be a positive integer");
        return false;
    }
    if (name == SYM_WorkingPrecision) {
        if (!nl_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            nl_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Sequence construction + dispatch                                   *
 * ------------------------------------------------------------------ */

/* Machine path: build S[0..terms-1], extrapolate, return result Expr or NULL. */
static Expr* nl_run_machine(Expr* expr, const char* var, double _Complex z0,
                            double _Complex dir, double _Complex scale,
                            bool infinite, NlOpts* o) {
    int terms = o->terms;
    double _Complex* S = malloc(sizeof(double _Complex) * (size_t)terms);
    if (!S) return NULL;

    NlBind bind; nl_bind_snapshot(&bind, var);
    NlCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = numeric_machine_spec();

    bool ok = true;
    for (int k = 0; k < terms; k++) {
        double _Complex z;
        if (infinite) z = dir * scale * ldexp(1.0, k);     /* u * Scale * 2^k */
        else          z = z0 - dir * scale * ldexp(1.0, -k); /* z0 - d*Scale*2^-k */
        if (!nl_sample_machine(&ctx, z, &S[k])) { ok = false; break; }
    }
    nl_bind_restore(&bind);
    if (!ok) {
        free(S);
        nl_warn("notnum", "the expression is not numerical at a sample point");
        return NULL;
    }

    double maxsample = 0.0;
    for (int k = 0; k < terms; k++) {
        double v = cabs(S[k]);
        if (v > maxsample) maxsample = v;
    }

    double _Complex result;
    double step = 0.0;
    bool got;
    if (o->method == SYM_SequenceLimit)
        got = nl_wynn_machine(S, terms, o->wynn, &result, &step);
    else
        got = nl_richardson_machine(S, terms, &result, &step);
    free(S);
    if (!got) { nl_warn("ndterm", "not enough Terms for the chosen Method"); return NULL; }

    if (!nl_accept(cabs(result), step, maxsample)) {
        nl_warn("noise", "Cannot recognize a limiting value. This may be due to "
                         "noise from roundoff; try higher WorkingPrecision, fewer "
                         "Terms, or a different Scale.");
        return NULL;
    }
    return nl_from_complex_d(result);
}

#ifdef USE_MPFR
static Expr* nl_run_mpfr(Expr* expr, const char* var, Expr* z0_expr,
                         Expr* dir_expr, Expr* scale_expr, bool infinite,
                         Expr* ray_expr, NlOpts* o) {
    int terms = o->terms;
    long bits = o->bits;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; spec.preserve_inexact = false;

    /* Resolve z0 (finite) or ray (infinite), Direction and Scale to complex. */
    mpfr_t z0r, z0i, dr, di, sr, si, mag;
    mpfr_init2(z0r, p); mpfr_init2(z0i, p);
    mpfr_init2(dr, p);  mpfr_init2(di, p);
    mpfr_init2(sr, p);  mpfr_init2(si, p);
    mpfr_init2(mag, p);
    bool prep = true;

    if (infinite) {
        prep = nl_const_complex_mpfr(ray_expr, bits, dr, di);
        mpfr_set_ui(z0r, 0, MPFR_RNDN); mpfr_set_ui(z0i, 0, MPFR_RNDN);
    } else {
        prep = nl_const_complex_mpfr(z0_expr, bits, z0r, z0i);
        if (prep) {
            if (dir_expr) prep = nl_const_complex_mpfr(dir_expr, bits, dr, di);
            else { mpfr_set_si(dr, -1, MPFR_RNDN); mpfr_set_ui(di, 0, MPFR_RNDN); }
        }
    }
    /* Normalise the direction/ray to a unit vector. */
    if (prep) {
        mpfr_hypot(mag, dr, di, MPFR_RNDN);
        if (mpfr_zero_p(mag)) prep = false;
        else { mpfr_div(dr, dr, mag, MPFR_RNDN); mpfr_div(di, di, mag, MPFR_RNDN); }
    }
    if (prep) {
        if (scale_expr) prep = nl_const_complex_mpfr(scale_expr, bits, sr, si);
        else { mpfr_set_ui(sr, 1, MPFR_RNDN); mpfr_set_ui(si, 0, MPFR_RNDN); }
    }
    if (prep && mpfr_zero_p(sr) && mpfr_zero_p(si)) prep = false;
    if (!prep) {
        mpfr_clears(z0r, z0i, dr, di, sr, si, mag, (mpfr_ptr)0);
        nl_warn("notnum", "the limit point, Direction or Scale is not numerical");
        return NULL;
    }

    mpfr_t* Sr = malloc(sizeof(mpfr_t) * (size_t)terms);
    mpfr_t* Si = malloc(sizeof(mpfr_t) * (size_t)terms);
    for (int k = 0; k < terms; k++) { mpfr_init2(Sr[k], p); mpfr_init2(Si[k], p); }
    /* scratch for the point and the d*Scale product */
    mpfr_t dsr, dsi, zr, zi, hk;
    mpfr_inits2(p, dsr, dsi, zr, zi, hk, (mpfr_ptr)0);
    nl_cmul(dsr, dsi, dr, di, sr, si);          /* d * Scale (unit dir * scale) */

    NlBind bind; nl_bind_snapshot(&bind, var);
    NlCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = spec;

    bool ok = true;
    for (int k = 0; ok && k < terms; k++) {
        if (infinite) {
            mpfr_set_ui(hk, 1, MPFR_RNDN);
            mpfr_mul_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);   /* 2^k */
            mpfr_mul(zr, dsr, hk, MPFR_RNDN);
            mpfr_mul(zi, dsi, hk, MPFR_RNDN);
        } else {
            mpfr_set_ui(hk, 1, MPFR_RNDN);
            mpfr_div_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);   /* 2^-k */
            mpfr_mul(zr, dsr, hk, MPFR_RNDN);
            mpfr_mul(zi, dsi, hk, MPFR_RNDN);
            mpfr_sub(zr, z0r, zr, MPFR_RNDN);    /* z0 - d*Scale*2^-k */
            mpfr_sub(zi, z0i, zi, MPFR_RNDN);
        }
        if (!nl_sample_mpfr(&ctx, zr, zi, Sr[k], Si[k])) ok = false;
    }
    nl_bind_restore(&bind);

    double maxsample = 0.0;
    if (ok) {
        for (int k = 0; k < terms; k++) {
            double v = nl_l1_d(Sr[k], Si[k]);
            if (v > maxsample) maxsample = v;
        }
    }

    Expr* out = NULL;
    if (ok) {
        bool converged = false;
        if (o->method == SYM_SequenceLimit)
            out = nl_wynn_mpfr(Sr, Si, terms, o->wynn, bits, maxsample, &converged);
        else
            out = nl_richardson_mpfr(Sr, Si, terms, bits, maxsample, &converged);
        if (!converged) {
            expr_free(out); out = NULL;
            nl_warn("noise", "Cannot recognize a limiting value. This may be due "
                             "to noise from roundoff; try higher WorkingPrecision, "
                             "fewer Terms, or a different Scale.");
        }
    } else {
        nl_warn("notnum", "the expression is not numerical at a sample point");
    }

    for (int k = 0; k < terms; k++) { mpfr_clear(Sr[k]); mpfr_clear(Si[k]); }
    free(Sr); free(Si);
    mpfr_clears(z0r, z0i, dr, di, sr, si, mag, dsr, dsi, zr, zi, hk, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nlimit(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;     /* expr, z -> z0  (+ options) */

    /* Peel trailing options; exactly two positional args must remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && nl_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nl_is_option_arg(res->data.function.args[i])) {
            nl_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) return NULL;

    /* Parse the z -> z0 rule. */
    Expr* rule = res->data.function.args[1];
    if (rule->type != EXPR_FUNCTION
        || rule->data.function.head->type != EXPR_SYMBOL
        || (rule->data.function.head->data.symbol != SYM_Rule
            && rule->data.function.head->data.symbol != SYM_RuleDelayed)
        || rule->data.function.arg_count != 2
        || rule->data.function.args[0]->type != EXPR_SYMBOL) {
        return NULL;               /* not a var -> point spec: stay unevaluated */
    }
    const char* var = rule->data.function.args[0]->data.symbol;
    Expr* z0_expr = rule->data.function.args[1];
    Expr* expr = res->data.function.args[0];

    NlOpts o;
    o.method = SYM_EulerSum;
    o.direction = NULL;
    o.scale = NULL;
    o.terms = 7;
    o.wynn = 1;
    o.prec_mpfr = false;
    o.bits = 0;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nl_apply_option(res->data.function.args[i], &o)) return NULL;
    }

    /* Classify the limit point. */
    Expr* ray = NULL;
    bool infinite = nl_infinite_ray(z0_expr, &ray);

#ifdef USE_MPFR
    if (o.prec_mpfr) {
        Expr* out = nl_run_mpfr(expr, var, z0_expr, o.direction, o.scale,
                                infinite, ray, &o);
        expr_free(ray);
        return out;
    }
#endif

    /* Machine path: resolve z0/ray, Direction and Scale to complex up front. */
    double _Complex z0 = 0.0, dir, scale;
    bool prep;
    if (infinite) {
        prep = nl_const_complex_machine(ray, &dir);
    } else {
        prep = nl_const_complex_machine(z0_expr, &z0);
        if (prep) {
            if (o.direction) prep = nl_const_complex_machine(o.direction, &dir);
            else dir = -1.0 + 0.0 * I;          /* Automatic == -1 */
        }
    }
    if (prep) {
        double m = cabs(dir);
        if (m == 0.0) prep = false; else dir /= m;   /* unit direction/ray */
    }
    if (prep) {
        if (o.scale) prep = nl_const_complex_machine(o.scale, &scale);
        else scale = 1.0 + 0.0 * I;
        if (prep && scale == 0.0) prep = false;
    }
    expr_free(ray);
    if (!prep) {
        /* Non-numeric limit point / Direction / Scale: stay unevaluated unless
         * it was an explicitly bad option. */
        return NULL;
    }

    return nl_run_machine(expr, var, z0, dir, scale, infinite, &o);
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nlimit_init(void) {
    symtab_add_builtin("NLimit", builtin_nlimit);
    /* Protected only. NLimit is deliberately NOT ATTR_LISTABLE: threading would
     * split the z -> z0 spec across bogus single-argument calls. */
    symtab_get_def("NLimit")->attributes |= ATTR_PROTECTED;
}

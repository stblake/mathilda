/*
 * nsum.c — NSum[f, {i, imin, imax (, di)}, opts]   (see nsum.h)
 *
 * Strategy
 * --------
 * NSum holds its arguments, evaluates the iterator bounds, then Block-localises
 * the index and evaluates the summand once per term.  Terms are reindexed to
 * k = 0, 1, 2, … with the actual index value x_k = imin + k·di, so a step di is
 * handled uniformly and multidimensional sums fall out by making the summand of
 * the outer sum an inner NSum[...] (HoldAll + localisation lets a dependent
 * inner bound such as {k,1,n} see the bound outer index).
 *
 * Methods are layered: this file currently provides Direct (small finite sums)
 * and WynnEpsilon (partial-sum extrapolation, shared seqaccel kernels), machine
 * and MPFR, real and complex.  Euler–Maclaurin and Cohen–Villegas–Zagier are
 * added on top of the same term machinery.
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.  Every temporary index
 * binding is removed on all return paths.
 */

#include "nsum.h"

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
#  include "numeric_complex.h"
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "seqaccel.h"
#include "dequad.h"
#include "sym_names.h"
#include "symtab.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Defaults (Mathematica-compatible where applicable). */
#define NS_DEF_NSUM_TERMS   15
#define NS_DEF_WYNN          1
#define NS_MAX_SEQ          64    /* cap on the partial-sum sequence length     */
#define NS_DIVERGE_FACTOR  1e6    /* lenient accept gate (Shanks may sum series) */
#define NS_LADDER_MAX       16    /* geometric far-tail rungs: 4,8,…,4·2^15≈131k */
#define NS_EM_MAX_HEAD    4096    /* cap on the adaptive Euler–Maclaurin base    */

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */

static void ns_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NSum::%s: ", tag);
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

static double ns_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}
#endif

/* ------------------------------------------------------------------ *
 *  Block-style index binding (mirrors nlimit.c)                       *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;
    Rule*       saved_own;
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
    symtab_add_own_value(b->name, sym, value);   /* copies value internally */
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
 *  Summand evaluation context                                         *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       body;    /* borrowed summand (free in the index variable) */
    Expr*       imin;    /* borrowed evaluated lower bound (a number)      */
    Expr*       di;      /* borrowed evaluated step (a number, default 1)  */
    NsBind*     bind;
    NumericSpec spec;
} NsCtx;

/* The actual index value x_k = imin + k·di as a fresh evaluated number. */
static Expr* ns_index_value(NsCtx* c, long k) {
    Expr* kk = expr_new_integer(k);
    Expr* term;
    {
        Expr** v = malloc(sizeof(Expr*) * 2);
        v[0] = kk; v[1] = expr_copy(c->di);
        term = expr_new_function(expr_new_symbol(SYM_Times), v, 2);
        free(v);
    }
    Expr** v = malloc(sizeof(Expr*) * 2);
    v[0] = expr_copy(c->imin); v[1] = term;
    Expr* sumexpr = expr_new_function(expr_new_symbol(SYM_Plus), v, 2);
    free(v);
    return eval_and_free(sumexpr);
}

/* Evaluate `e` with the index bound to `value` (consumed). Returns the
 * numericalised result (caller frees) or NULL. */
static Expr* ns_eval_expr_at(NsCtx* c, Expr* e, Expr* value) {
    ns_bind_set(c->bind, value);
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(e));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

/* Evaluate the summand at the index value (consumed). */
static Expr* ns_eval_at(NsCtx* c, Expr* value) {
    return ns_eval_expr_at(c, c->body, value);
}

/* term_k as a machine complex. */
static bool ns_term_machine(NsCtx* c, long k, double _Complex* out) {
    Expr* num = ns_eval_at(c, ns_index_value(c, k));
    if (!num) return false;
    bool ok = ns_to_complex(num, out);
    if (ok && (!isfinite(creal(*out)) || !isfinite(cimag(*out)))) ok = false;
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
/* term_k as an MPFR complex (out_re, out_im pre-initialised by caller). */
static bool ns_term_mpfr(NsCtx* c, long k, mpfr_t out_re, mpfr_t out_im) {
    Expr* num = ns_eval_at(c, ns_index_value(c, k));
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}
#endif

/* ------------------------------------------------------------------ *
 *  Iterator-bound classification                                      *
 * ------------------------------------------------------------------ */

static bool ns_is_infinite(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_Infinity || e->data.symbol.name == SYM_ComplexInfinity;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_DirectedInfinity)
        return true;
    /* -Infinity etc.: Times[..., Infinity] */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (ns_is_infinite(e->data.function.args[i])) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Black-box summand detection                                        *
 *                                                                     *
 *  A summand that is itself a numeric routine (nested NSum/NProduct/…) *
 *  is only meaningful at integer indices and cannot be differentiated  *
 *  symbolically or sampled on a complex contour, so Euler–Maclaurin    *
 *  (which needs both) must never be applied to it — partial-sum         *
 *  extrapolation (Wynn / CVZ) is the only valid route.                 *
 * ------------------------------------------------------------------ */

static bool ns_head_is_numeric_blackbox(const char* s) {
    return s == SYM_NSum || s == SYM_NProduct || s == SYM_NIntegrate
        || s == SYM_NSeries || s == SYM_NLimit || s == SYM_NResidue
        || s == SYM_ND;
}

static bool ns_body_is_blackbox(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && ns_head_is_numeric_blackbox(e->data.function.head->data.symbol.name))
        return true;
    if (ns_body_is_blackbox(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ns_body_is_blackbox(e->data.function.args[i])) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* method;     /* SYM_Automatic / WynnEpsilon / EulerMaclaurin / ... */
    int   nsum_terms;       /* head terms summed explicitly                       */
    bool  nsum_terms_user;  /* user pinned NSumTerms (suppress adaptive base)     */
    int   extra_terms;      /* extrapolation-sequence length (-1 => auto)         */
    int   wynn;             /* WynnDegree                                         */
    int   levin_variant;    /* SeqaccelLevinVariant when method == SYM_Levin       */
    bool  verify;           /* VerifyConvergence                                  */
    bool  prec_mpfr;        /* WorkingPrecision selects MPFR                       */
    long  bits;             /* MPFR *internal* working precision in bits (guarded) */
    long  target_bits;      /* user-requested precision; result rounded back to it */
    double acc_goal;        /* AccuracyGoal digits (-1 => Infinity)               */
    double prec_goal;       /* PrecisionGoal digits (-1 => Automatic)            */
} NsOpts;

static bool ns_is_known_option(const char* s) {
    return s == SYM_Method || s == SYM_WorkingPrecision || s == SYM_NSumTerms
        || s == SYM_NSumExtraTerms || s == SYM_WynnDegree || s == SYM_VerifyConvergence
        || s == SYM_AccuracyGoal || s == SYM_PrecisionGoal
        || s == SYM_Compiled || s == SYM_EvaluationMonitor;
}

static bool ns_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && ns_is_known_option(lhs->data.symbol.name);
}

static bool ns_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
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

static bool ns_parse_goal(Expr* v, double* out) {
    if (v->type == EXPR_SYMBOL
        && (v->data.symbol.name == SYM_Infinity || v->data.symbol.name == SYM_Automatic)) {
        *out = -1.0; return true;
    }
    double d;
    if (ns_to_double_real(v, &d) && d > 0.0) { *out = d; return true; }
    return false;
}

static bool ns_apply_option(Expr* rule, NsOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_Method) {
        /* Accept a bare symbol, or Method -> {sym, subopts...} (subopts ignored
         * for now: the EulerMaclaurin NIntegrate sub-method tuning lands later). */
        Expr* m = rhs;
        if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL
            && m->data.function.head->data.symbol.name == SYM_List
            && m->data.function.arg_count >= 1)
            m = m->data.function.args[0];
        if (m->type == EXPR_STRING) {
            /* tolerate "WynnEpsilon"-style string method names */
            const char* s = m->data.string;
            if (!strcmp(s, "WynnEpsilon") || !strcmp(s, "SequenceLimit")) o->method = SYM_WynnEpsilon;
            else if (!strcmp(s, "EulerMaclaurin") || !strcmp(s, "Integrate")) o->method = SYM_EulerMaclaurin;
            else if (!strcmp(s, "AlternatingSigns")) o->method = SYM_AlternatingSigns;
            else if (!strcmp(s, "Levin") || !strcmp(s, "LevinU")) { o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_U; }
            else if (!strcmp(s, "LevinT")) { o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_T; }
            else if (!strcmp(s, "LevinV")) { o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_V; }
            else o->method = SYM_Automatic;
            return true;
        }
        if (m->type == EXPR_SYMBOL) {
            const char* sym = m->data.symbol.name;
            if (sym == SYM_WynnEpsilon || sym == SYM_SequenceLimit) o->method = SYM_WynnEpsilon;
            else if (sym == SYM_EulerMaclaurin || sym == SYM_Integrate
                     || sym == SYM_EulerSum) o->method = SYM_EulerMaclaurin;
            else if (sym == SYM_AlternatingSigns) o->method = SYM_AlternatingSigns;
            else if (sym == SYM_Levin) { o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_U; }
            else o->method = SYM_Automatic;
            return true;
        }
        o->method = SYM_Automatic;
        return true;
    }
    if (name == SYM_WorkingPrecision) {
        if (!ns_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            ns_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_NSumTerms) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 0) {
            o->nsum_terms = (int)rhs->data.integer; o->nsum_terms_user = true; return true;
        }
        ns_warn("badopt", "NSumTerms must be a non-negative integer"); return false;
    }
    if (name == SYM_NSumExtraTerms) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 1) {
            o->extra_terms = (int)rhs->data.integer; return true;
        }
        ns_warn("badopt", "NSumExtraTerms must be a positive integer"); return false;
    }
    if (name == SYM_WynnDegree) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 1) {
            o->wynn = (int)rhs->data.integer; return true;
        }
        ns_warn("badopt", "WynnDegree must be a positive integer"); return false;
    }
    if (name == SYM_VerifyConvergence) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_True)  { o->verify = true;  return true; }
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_False) { o->verify = false; return true; }
        ns_warn("badopt", "VerifyConvergence must be True or False"); return false;
    }
    if (name == SYM_AccuracyGoal)  return ns_parse_goal(rhs, &o->acc_goal);
    if (name == SYM_PrecisionGoal) return ns_parse_goal(rhs, &o->prec_goal);
    /* Compiled / EvaluationMonitor: accepted, ignored. */
    if (name == SYM_Compiled || name == SYM_EvaluationMonitor) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Accept gate                                                        *
 * ------------------------------------------------------------------ */

static bool ns_accept(double result_mag, double step, double maxsample) {
    if (!isfinite(result_mag) || !isfinite(step)) return false;
    if (maxsample < 1e-300) maxsample = 1e-300;
    if (result_mag > NS_DIVERGE_FACTOR * maxsample + 1.0) return false;
    if (step > 0.5 * maxsample + 1e-6) return false;
    return true;
}

/* Auto extrapolation-sequence length for the chosen precision. */
static int ns_seq_len(const NsOpts* o) {
    if (o->extra_terms > 0) {
        int L = o->extra_terms;
        if (L > NS_MAX_SEQ) L = NS_MAX_SEQ;
        return L;
    }
    int L = 12;
#ifdef USE_MPFR
    if (o->prec_mpfr) {
        /* TARGET precision: the guarded working precision would over-extend the
         * partial-sum sequence (and Wynn's bottom-corner denominators amplify). */
        long b = o->target_bits > 0 ? o->target_bits : o->bits;
        int want = (int)(b * 0.55) + 6;
        if (want > L) L = want;
    }
#endif
    if (L > NS_MAX_SEQ) L = NS_MAX_SEQ;
    if (L < 2 * o->wynn + 2) L = 2 * o->wynn + 2;
    return L;
}

/* ------------------------------------------------------------------ *
 *  Method: WynnEpsilon (partial-sum extrapolation)                    *
 * ------------------------------------------------------------------ */

/* Machine path. count<0 means infinite; otherwise terms beyond `count` are 0. */
static Expr* ns_wynn_machine(NsCtx* c, long count, const NsOpts* o) {
    int head = o->nsum_terms < 0 ? 0 : o->nsum_terms;
    int L = ns_seq_len(o);
    double _Complex* P = malloc(sizeof(double _Complex) * (size_t)L);
    if (!P) return NULL;

    double _Complex acc = 0.0;
    bool ok = true;
    long k = 0;
    for (; k < head; k++) {
        if (count >= 0 && k >= count) break;
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { ok = false; break; }
        acc += t;
    }
    P[0] = acc;
    for (int j = 1; ok && j < L; j++) {
        if (count >= 0 && k >= count) { P[j] = acc; k++; continue; }
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { ok = false; break; }
        acc += t;
        P[j] = acc;
        k++;
    }
    if (!ok) { free(P); ns_warn("nnum", "summand is not numerical at a term"); return NULL; }

    double maxsample = 0.0;
    for (int j = 0; j < L; j++) { double v = cabs(P[j]); if (v > maxsample) maxsample = v; }

    double _Complex result; double step = 0.0;
    bool got = seqaccel_wynn_machine(P, L, o->wynn, &result, &step);
    free(P);
    if (!got) { ns_warn("ndterm", "not enough terms for WynnEpsilon"); return NULL; }
    if (!ns_accept(cabs(result), step, maxsample)) {
        ns_warn("ncvg", "failed to converge; try more NSumExtraTerms or higher WorkingPrecision");
        return NULL;
    }
    return ns_from_complex_d(result);
}

#ifdef USE_MPFR
static Expr* ns_wynn_mpfr(NsCtx* c, long count, const NsOpts* o) {
    int head = o->nsum_terms < 0 ? 0 : o->nsum_terms;
    int L = ns_seq_len(o);
    long bits = o->bits;
    mpfr_prec_t p = (mpfr_prec_t)bits;

    mpfr_t* Pr = malloc(sizeof(mpfr_t) * (size_t)L);
    mpfr_t* Pi = malloc(sizeof(mpfr_t) * (size_t)L);
    for (int j = 0; j < L; j++) { mpfr_init2(Pr[j], p); mpfr_init2(Pi[j], p); }
    mpfr_t accr, acci, tr, ti;
    mpfr_inits2(p, accr, acci, tr, ti, (mpfr_ptr)0);
    mpfr_set_ui(accr, 0, MPFR_RNDN); mpfr_set_ui(acci, 0, MPFR_RNDN);

    bool ok = true;
    long k = 0;
    for (; k < head; k++) {
        if (count >= 0 && k >= count) break;
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(accr, accr, tr, MPFR_RNDN);
        mpfr_add(acci, acci, ti, MPFR_RNDN);
    }
    mpfr_set(Pr[0], accr, MPFR_RNDN); mpfr_set(Pi[0], acci, MPFR_RNDN);
    for (int j = 1; ok && j < L; j++) {
        if (count >= 0 && k >= count) {
            mpfr_set(Pr[j], accr, MPFR_RNDN); mpfr_set(Pi[j], acci, MPFR_RNDN); k++; continue;
        }
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(accr, accr, tr, MPFR_RNDN);
        mpfr_add(acci, acci, ti, MPFR_RNDN);
        mpfr_set(Pr[j], accr, MPFR_RNDN); mpfr_set(Pi[j], acci, MPFR_RNDN);
        k++;
    }

    Expr* out = NULL;
    if (ok) {
        double maxsample = 0.0;
        for (int j = 0; j < L; j++) { double v = ns_l1_d(Pr[j], Pi[j]); if (v > maxsample) maxsample = v; }
        mpfr_t rr, ri; mpfr_init2(rr, p); mpfr_init2(ri, p);
        double step = 0.0; bool finite = false;
        bool got = seqaccel_wynn_mpfr(Pr, Pi, L, o->wynn, bits, rr, ri, &step, &finite);
        if (got && finite && ns_accept(ns_l1_d(rr, ri), step, maxsample))
            out = ns_from_complex_mpfr(rr, ri);
        else
            ns_warn("ncvg", "failed to converge; try more NSumExtraTerms or higher WorkingPrecision");
        mpfr_clear(rr); mpfr_clear(ri);
    } else {
        ns_warn("nnum", "summand is not numerical at a term");
    }

    for (int j = 0; j < L; j++) { mpfr_clear(Pr[j]); mpfr_clear(Pi[j]); }
    free(Pr); free(Pi);
    mpfr_clears(accr, acci, tr, ti, (mpfr_ptr)0);
    return out;
}
#endif

/* ------------------------------------------------------------------ *
 *  Method: Levin (partial-sum transformation, u/t/v)                  *
 * ------------------------------------------------------------------ */

/* Machine path.  Structurally identical to ns_wynn_machine but drives the
 * shared Levin kernel with o->levin_variant. */
static Expr* ns_levin_machine(NsCtx* c, long count, const NsOpts* o) {
    int head = o->nsum_terms < 0 ? 0 : o->nsum_terms;
    int L = ns_seq_len(o);
    double _Complex* P = malloc(sizeof(double _Complex) * (size_t)L);
    if (!P) return NULL;

    double _Complex acc = 0.0;
    bool ok = true;
    long k = 0;
    for (; k < head; k++) {
        if (count >= 0 && k >= count) break;
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { ok = false; break; }
        acc += t;
    }
    P[0] = acc;
    for (int j = 1; ok && j < L; j++) {
        if (count >= 0 && k >= count) { P[j] = acc; k++; continue; }
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { ok = false; break; }
        acc += t;
        P[j] = acc;
        k++;
    }
    if (!ok) { free(P); ns_warn("nnum", "summand is not numerical at a term"); return NULL; }

    double maxsample = 0.0;
    for (int j = 0; j < L; j++) { double v = cabs(P[j]); if (v > maxsample) maxsample = v; }

    double _Complex result; double step = 0.0;
    bool got = seqaccel_levin_machine(P, L, o->levin_variant, 1.0, &result, &step);
    free(P);
    if (!got) { ns_warn("ndterm", "not enough terms for Levin"); return NULL; }
    if (!ns_accept(cabs(result), step, maxsample)) {
        ns_warn("ncvg", "failed to converge; try more NSumExtraTerms or higher WorkingPrecision");
        return NULL;
    }
    return ns_from_complex_d(result);
}

#ifdef USE_MPFR
static Expr* ns_levin_mpfr(NsCtx* c, long count, const NsOpts* o) {
    int head = o->nsum_terms < 0 ? 0 : o->nsum_terms;
    int L = ns_seq_len(o);
    long bits = o->bits;
    mpfr_prec_t p = (mpfr_prec_t)bits;

    mpfr_t* Pr = malloc(sizeof(mpfr_t) * (size_t)L);
    mpfr_t* Pi = malloc(sizeof(mpfr_t) * (size_t)L);
    for (int j = 0; j < L; j++) { mpfr_init2(Pr[j], p); mpfr_init2(Pi[j], p); }
    mpfr_t accr, acci, tr, ti;
    mpfr_inits2(p, accr, acci, tr, ti, (mpfr_ptr)0);
    mpfr_set_ui(accr, 0, MPFR_RNDN); mpfr_set_ui(acci, 0, MPFR_RNDN);

    bool ok = true;
    long k = 0;
    for (; k < head; k++) {
        if (count >= 0 && k >= count) break;
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(accr, accr, tr, MPFR_RNDN);
        mpfr_add(acci, acci, ti, MPFR_RNDN);
    }
    mpfr_set(Pr[0], accr, MPFR_RNDN); mpfr_set(Pi[0], acci, MPFR_RNDN);
    for (int j = 1; ok && j < L; j++) {
        if (count >= 0 && k >= count) {
            mpfr_set(Pr[j], accr, MPFR_RNDN); mpfr_set(Pi[j], acci, MPFR_RNDN); k++; continue;
        }
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(accr, accr, tr, MPFR_RNDN);
        mpfr_add(acci, acci, ti, MPFR_RNDN);
        mpfr_set(Pr[j], accr, MPFR_RNDN); mpfr_set(Pi[j], acci, MPFR_RNDN);
        k++;
    }

    Expr* out = NULL;
    if (ok) {
        double maxsample = 0.0;
        for (int j = 0; j < L; j++) { double v = ns_l1_d(Pr[j], Pi[j]); if (v > maxsample) maxsample = v; }
        mpfr_t rr, ri; mpfr_init2(rr, p); mpfr_init2(ri, p);
        double step = 0.0; bool finite = false;
        bool got = seqaccel_levin_mpfr(Pr, Pi, L, o->levin_variant, 1.0, bits, rr, ri, &step, &finite);
        if (got && finite && ns_accept(ns_l1_d(rr, ri), step, maxsample))
            out = ns_from_complex_mpfr(rr, ri);
        else
            ns_warn("ncvg", "failed to converge; try more NSumExtraTerms or higher WorkingPrecision");
        mpfr_clear(rr); mpfr_clear(ri);
    } else {
        ns_warn("nnum", "summand is not numerical at a term");
    }

    for (int j = 0; j < L; j++) { mpfr_clear(Pr[j]); mpfr_clear(Pi[j]); }
    free(Pr); free(Pi);
    mpfr_clears(accr, acci, tr, ti, (mpfr_ptr)0);
    return out;
}
#endif

/* ------------------------------------------------------------------ *
 *  Method: Direct (exact small finite sum)                            *
 * ------------------------------------------------------------------ */

static Expr* ns_direct_machine(NsCtx* c, long count) {
    double _Complex acc = 0.0;
    for (long k = 0; k < count; k++) {
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { ns_warn("nnum", "summand is not numerical at a term"); return NULL; }
        acc += t;
    }
    return ns_from_complex_d(acc);
}

#ifdef USE_MPFR
static Expr* ns_direct_mpfr(NsCtx* c, long count, long bits) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t accr, acci, tr, ti;
    mpfr_inits2(p, accr, acci, tr, ti, (mpfr_ptr)0);
    mpfr_set_ui(accr, 0, MPFR_RNDN); mpfr_set_ui(acci, 0, MPFR_RNDN);
    bool ok = true;
    for (long k = 0; k < count; k++) {
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(accr, accr, tr, MPFR_RNDN);
        mpfr_add(acci, acci, ti, MPFR_RNDN);
    }
    Expr* out = ok ? ns_from_complex_mpfr(accr, acci) : NULL;
    if (!ok) ns_warn("nnum", "summand is not numerical at a term");
    mpfr_clears(accr, acci, tr, ti, (mpfr_ptr)0);
    return out;
}
#endif

/* ------------------------------------------------------------------ *
 *  Method: Euler–Maclaurin                                            *
 *                                                                     *
 *  Σ_{x=N,N+di,…} f(x) ≈ (1/di)∫_N^∞ f dx + f(N)/2                    *
 *                        − Σ_{j≥1} B_{2j}/(2j)! · di^{2j-1} f^(2j-1)(N)*
 *  with N = imin + max(NSumTerms, settle)·di and the head terms summed         *
 *  exactly.  The tail integral uses exp-sinh quadrature (tolerance scaled to   *
 *  WorkingPrecision); the derivative corrections come from numerical contour   *
 *  (circle-DFT) Taylor coefficients and are truncated at the smallest term     *
 *  (the series is asymptotic).                                                 *
 * ------------------------------------------------------------------ */

/* Small held-expression builders. */
static Expr* ns_mk1(const char* h, Expr* a) {
    Expr** v = malloc(sizeof(Expr*)); v[0] = a;
    Expr* e = expr_new_function(expr_new_symbol(h), v, 1); free(v); return e;
}
static Expr* ns_mk2(const char* h, Expr* a, Expr* b) {
    Expr** v = malloc(sizeof(Expr*) * 2); v[0] = a; v[1] = b;
    Expr* e = expr_new_function(expr_new_symbol(h), v, 2); free(v); return e;
}

/* Numericalise a var-free constant expression (consumed) to machine complex. */
static bool ns_const_machine(Expr* e, double _Complex* out) {
    Expr* raw = eval_and_free(e);
    Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
    expr_free(raw);
    bool ok = num && ns_to_complex(num, out);
    expr_free(num);
    return ok;
}

/* Evaluate `e` at the index value `value` (consumed) to a machine complex. */
static bool ns_eval_complex_machine(NsCtx* c, Expr* e, Expr* value, double _Complex* out) {
    Expr* num = ns_eval_expr_at(c, e, value);
    if (!num) return false;
    bool ok = ns_to_complex(num, out);
    if (ok && (!isfinite(creal(*out)) || !isfinite(cimag(*out)))) ok = false;
    expr_free(num);
    return ok;
}

/* exp-sinh integrand sample: f at a continuous real x, machine. */
static bool ns_sample_x_machine(void* vc, double x, double _Complex* out) {
    NsCtx* c = vc;
    return ns_eval_complex_machine(c, c->body, expr_new_real(x), out);
}

/* ------------------------------------------------------------------ *
 *  Numerical contour (circle-DFT) derivatives for the EM corrections  *
 *                                                                     *
 *  The Euler–Maclaurin corrections need the odd derivatives f^(2j-1)  *
 *  at the base point N.  Iterated symbolic D makes the expression tree *
 *  explode for composite summands (Log[...], …), truncating the series *
 *  and capping accuracy.  Instead we recover the Taylor coefficients   *
 *      a_m = f^(m)(N)/m!                                               *
 *  of the summand (as a function of the continuous index) from one DFT *
 *  of N_s samples on a circle |z-N| = r — Cauchy's formula, exactly as *
 *  NSeries does (nseries.c).  This is stable to full precision and     *
 *  never balloons.  The EM term B_{2j}/(2j)!·di^(2j-1)·f^(2j-1)(N)      *
 *  folds the m! into the coefficient as B_{2j}/(2j)·di^(2j-1)·a_(2j-1).  *
 *  The radius is chosen adaptively so the sampling disk is             *
 *  singularity-free (the Fourier coefficients must decay).             *
 * ------------------------------------------------------------------ */

/* Number of correction terms (asymptotic series; the loop also stops at the
 * series minimum).  Grows with the target precision so high WorkingPrecision
 * reaches its EM floor; capped to bound the contour order. */
static int ns_em_maxj(const NsOpts* o) {
    /* Machine precision saturates at ~12 correction terms — more only inflates
     * the contour order (and sample count) for no gain.  Arbitrary precision
     * scales with the target digit count. */
    if (!o->prec_mpfr) return 12;
    double digits = numeric_bits_to_digits(o->target_bits);
    int J = (int)(digits * 0.8) + 6;
    if (J < 12) J = 12;
    if (J > 64) J = 64;
    return J;
}

/* Contour sample count: a power of two with 4x oversampling beyond order M. */
static int ns_contour_samples(int M) {
    int ne = M < 1 ? 1 : M;
    int p = 0; while ((1 << p) < ne) p++;
    int N = 1 << (p + 2);
    if (N < 8) N = 8;
    while (N <= 2 * M) N <<= 1;
    return N;
}

/* Machine Taylor coefficients a[0..M] about the real point Nd, radius r.
 * Returns false on a non-numeric sample or if the Fourier coefficients fail to
 * decay (the disk is not singularity-free). */
static bool ns_taylor_machine(NsCtx* c, double Nd, double r, int M, double _Complex* a) {
    int Ns = ns_contour_samples(M);
    double _Complex* F = malloc(sizeof(double _Complex) * (size_t)Ns);
    if (!F) return false;
    bool ok = true;
    for (int j = 0; j < Ns; j++) {
        double th = 2.0 * M_PI * (double)j / (double)Ns;
        double zr = Nd + r * cos(th), zi = r * sin(th);
        Expr* zv = (zi == 0.0) ? expr_new_real(zr)
                               : make_complex(expr_new_real(zr), expr_new_real(zi));
        if (!ns_eval_complex_machine(c, c->body, zv, &F[j])) { ok = false; break; }
    }
    if (ok) {
        double cmax = 0.0, cM = 0.0;
        for (int m = 0; m <= M; m++) {
            int k = ((m % Ns) + Ns) % Ns;
            double _Complex acc = 0.0;
            for (int j = 0; j < Ns; j++) {
                double ang = -2.0 * M_PI * (double)j * (double)k / (double)Ns;
                acc += F[j] * (cos(ang) + I * sin(ang));
            }
            acc /= (double)Ns;                  /* c_k = a_m · r^m              */
            double ck = cabs(acc);
            if (ck > cmax) cmax = ck;
            if (m == M) cM = ck;
            a[m] = acc * pow(r, -(double)m);     /* a_m = c_k · r^(-m)          */
        }
        if (!(cM <= 1e-4 * (cmax + 1e-300))) ok = false;   /* decay gate */
    }
    free(F);
    return ok;
}

/* Adaptive radius: start near N/2 and halve until the disk is singularity-free. */
static bool ns_taylor_machine_auto(NsCtx* c, double Nd, int M, double _Complex* a) {
    double r = 0.5 * fabs(Nd);
    if (r < 1.0) r = 1.0;
    for (int tries = 0; tries < 8; tries++) {
        if (ns_taylor_machine(c, Nd, r, M, a)) return true;
        r *= 0.5;
        if (r < 1e-3) break;
    }
    return false;
}

/* B_{2j}/(2j) as a real double (the EM coefficient with the m! folded into the
 * contour coefficient a_m = f^(m)/m!). */
static double ns_em_bcoeff_machine(int j) {
    Expr* b = ns_mk1("BernoulliB", expr_new_integer(2 * j));
    Expr* ratio = ns_mk2("Times", b,
                         ns_mk2("Power", expr_new_integer(2 * j), expr_new_integer(-1)));
    double _Complex v;
    return ns_const_machine(ratio, &v) ? creal(v) : 0.0;
}

/* B_{2j}/(2j)! as a real double (for the symbolic path, which supplies the raw
 * derivative value f^(2j-1)(N) rather than the Taylor coefficient). */
static double ns_em_coeff_machine(int j) {
    Expr* b = ns_mk1("BernoulliB", expr_new_integer(2 * j));
    Expr* f = ns_mk1("Factorial", expr_new_integer(2 * j));
    Expr* ratio = ns_mk2("Times", b, ns_mk2("Power", f, expr_new_integer(-1)));
    double _Complex v;
    return ns_const_machine(ratio, &v) ? creal(v) : 0.0;
}

/* Bounded node count (min(size, limit+1)) without walking a huge tree.  Iterated
 * symbolic derivatives of composite summands (Log[…], 1/Fibonacci[x]) blow up;
 * the hybrid correction loop uses this to switch from symbolic D to numerical
 * contour derivatives once the tree balloons. */
static int ns_expr_size(const Expr* e, int limit) {
    if (!e || limit <= 0) return 0;
    int n = 1;
    if (e->type == EXPR_FUNCTION) {
        n += ns_expr_size(e->data.function.head, limit - n);
        for (size_t i = 0; i < e->data.function.arg_count && n <= limit; i++)
            n += ns_expr_size(e->data.function.args[i], limit - n);
    }
    return n;
}
#define NS_EM_MAX_NODES 1200   /* switch to contour once a derivative exceeds this */

/* Effective head-term count: sum exactly through any peak (the decay onset at
 * `settle`) so the tail integral starts in the monotone region — unless the
 * user pinned NSumTerms. */
static int ns_em_head(const NsOpts* o, long settle) {
    int m = o->nsum_terms < 0 ? 0 : o->nsum_terms;
    if (!o->nsum_terms_user && settle > 0) {
        long want = settle + 4;
        if (want > NS_EM_MAX_HEAD) want = NS_EM_MAX_HEAD;
        if ((int)want > m) m = (int)want;
    }
    return m;
}

static Expr* ns_em_machine(NsCtx* c, const char* var, NsOpts* o, long settle) {
    int m = ns_em_head(o, settle);
    double _Complex H = 0.0, tprev = 0.0, tlast = 0.0;
    for (int k = 0; k < m; k++) {
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) return NULL;
        H += t;
        tprev = tlast; tlast = t;
    }
    /* A sign-alternating summand has an oscillatory continuous extension (e.g.
     * (-5)^x/x! = e^{x ln(-5)}/x!), so its values swing by many orders of
     * magnitude around the contour circle and the circle-DFT is ill-conditioned.
     * Detect it from the last two head terms' signs (no extra evaluation) and
     * keep to symbolic derivatives only — the original behaviour for such sums. */
    bool oscillatory = (m >= 2 && creal(tprev) * creal(tlast) < 0.0);

    Expr* Nval = ns_index_value(c, m);
    double Nd, ddi;
    if (!ns_to_double_real(Nval, &Nd) || !ns_to_double_real(c->di, &ddi) || ddi == 0.0) {
        expr_free(Nval); return NULL;
    }
    double _Complex fN;
    if (!ns_eval_complex_machine(c, c->body, expr_copy(Nval), &fN)) { expr_free(Nval); return NULL; }

    double _Complex Iv; double abserr;
    bool integ = dequad_halfline_machine(ns_sample_x_machine, c, Nd, 1e-12, 14, &Iv, &abserr);
    double relerr = abserr / (cabs(Iv) + 1e-300);
    if (!integ && relerr > 1e-6) { expr_free(Nval); return NULL; }   /* fall back */

    double _Complex tail = Iv / ddi + fN / 2.0;

    /* Hybrid corrections: symbolic D as the primary (robust for simple summands —
     * exactly the original behaviour); switch to numerical contour derivatives
     * for the order at which the derivative tree balloons (composite summands
     * like Log[1+1/n^2], where symbolic D would truncate and cap accuracy). */
    int J = ns_em_maxj(o);
    int M = 2 * J - 1;
    double _Complex* a = NULL;            /* contour coeffs, computed on demand */
    bool a_tried = false, a_ok = false;
    ns_bind_clear_temp(symtab_get_def(var));
    Expr* dcur = expr_copy(c->body);
    int cur = 0;
    bool ballooned = false;
    int napplied = 0;
    double prev_mag = INFINITY;
    double goal = (o->prec_goal > 0 ? o->prec_goal : NUMERIC_MACHINE_PRECISION_DIGITS - 2);
    double tol = pow(10.0, -goal) * (cabs(H + tail) + 1.0);
    double dipow = ddi;                          /* di^(2j-1): di^1, ·= di^2 each step */
    for (int j = 1; j <= J; j++) {
        int ord = 2 * j - 1;
        double _Complex corr; bool have_corr = false;
        if (!ballooned) {
            while (dcur && cur < ord) {
                dcur = eval_and_free(ns_mk2("D", dcur, expr_new_symbol(var)));
                cur++;
                if (dcur && ns_expr_size(dcur, NS_EM_MAX_NODES + 1) > NS_EM_MAX_NODES) {
                    ballooned = true; break;
                }
            }
            double _Complex dval;
            if (dcur && !ballooned
                && ns_eval_complex_machine(c, dcur, expr_copy(Nval), &dval)) {
                ns_bind_clear_temp(symtab_get_def(var));   /* free var for next D */
                corr = ns_em_coeff_machine(j) * dipow * dval;   /* B/(2j)!·di·f^(ord) */
                have_corr = true;
            } else if (!dcur) {
                break;                            /* derivative vanished: done */
            }
        }
        if (!have_corr) {                         /* symbolic failed/ballooned */
            if (oscillatory) {
                /* contour is unreliable for an oscillatory extension; if no
                 * correction was ever applied the bare integral is too poor —
                 * defer to Wynn (which sums alternating series well). */
                if (napplied == 0) { free(a); expr_free(dcur); expr_free(Nval); return NULL; }
                break;                            /* keep corrections summed so far */
            }
            if (!a_tried) {
                a_tried = true;
                a = malloc(sizeof(double _Complex) * (size_t)(M + 1));
                a_ok = a && ns_taylor_machine_auto(c, Nd, M, a);
            }
            if (!a_ok) break;                     /* keep corrections summed so far */
            corr = ns_em_bcoeff_machine(j) * dipow * a[ord];    /* B/(2j)·di·a_ord */
        }
        double mag = cabs(corr);
        if (j > 1 && mag > prev_mag) break;       /* asymptotic minimum */
        tail -= corr;
        napplied++;
        prev_mag = mag;
        if (mag < tol) break;
        dipow *= ddi * ddi;
    }
    free(a);
    expr_free(dcur);
    expr_free(Nval);
    /* Chop the sub-tolerance imaginary residual a contour DFT leaves on a
     * real-valued summand (genuinely complex sums keep |Im| ~ |Re|). */
    double _Complex r = H + tail;
    if (fabs(cimag(r)) < 1e-12 * (fabs(creal(r)) + 1.0)) r = creal(r);
    return ns_from_complex_d(r);
}

#ifdef USE_MPFR
/* f at a continuous real x, MPFR. */
static bool ns_sample_x_mpfr(void* vc, const mpfr_t x, mpfr_t out_re, mpfr_t out_im) {
    NsCtx* c = vc;
    Expr* num = ns_eval_expr_at(c, c->body, expr_new_mpfr_copy(x));
    if (!num) return false;
    bool inexact, ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}

/* Evaluate `e` at the index value (consumed) -> (re,im) MPFR. */
static bool ns_eval_complex_mpfr(NsCtx* c, Expr* e, Expr* value, mpfr_t re, mpfr_t im) {
    Expr* num = ns_eval_expr_at(c, e, value);
    if (!num) return false;
    bool inexact, ok = get_approx_mpfr(num, re, im, &inexact);
    if (ok && (!mpfr_number_p(re) || !mpfr_number_p(im))) ok = false;
    expr_free(num);
    return ok;
}

/* B_{2j}/(2j) as a real MPFR at `bits` (the EM coefficient with the m! folded
 * into the contour coefficient a_m = f^(m)/m!). */
static bool ns_em_bcoeff_mpfr(int j, long bits, mpfr_t out) {
    Expr* b = ns_mk1("BernoulliB", expr_new_integer(2 * j));
    Expr* ratio = ns_mk2("Times", b,
                         ns_mk2("Power", expr_new_integer(2 * j), expr_new_integer(-1)));
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; spec.preserve_inexact = false;
    Expr* raw = eval_and_free(ratio);
    Expr* num = raw ? numericalize(raw, spec) : NULL;
    expr_free(raw);
    bool ok = false, inexact;
    if (num) {
        mpfr_t im; mpfr_init2(im, (mpfr_prec_t)bits);
        ok = get_approx_mpfr(num, out, im, &inexact);
        mpfr_clear(im);
    }
    expr_free(num);
    return ok;
}

/* B_{2j}/(2j)! as a real MPFR (for the symbolic path: raw derivative value). */
static bool ns_em_coeff_mpfr(int j, long bits, mpfr_t out) {
    Expr* b = ns_mk1("BernoulliB", expr_new_integer(2 * j));
    Expr* f = ns_mk1("Factorial", expr_new_integer(2 * j));
    Expr* ratio = ns_mk2("Times", b, ns_mk2("Power", f, expr_new_integer(-1)));
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; spec.preserve_inexact = false;
    Expr* raw = eval_and_free(ratio);
    Expr* num = raw ? numericalize(raw, spec) : NULL;
    expr_free(raw);
    bool ok = false, inexact;
    if (num) {
        mpfr_t im; mpfr_init2(im, (mpfr_prec_t)bits);
        ok = get_approx_mpfr(num, out, im, &inexact);
        mpfr_clear(im);
    }
    expr_free(num);
    return ok;
}

/* MPFR Taylor coefficients a_re/a_im[0..M] about the base point Nbase (a number;
 * Nd is its double for the radius), radius r, at `bits`.  Decay-checked.  The
 * arrays must be pre-initialised by the caller. */
static bool ns_taylor_mpfr(NsCtx* c, const mpfr_t Nbase, double r, long bits,
                           int M, mpfr_t* a_re, mpfr_t* a_im) {
    int Ns = ns_contour_samples(M);
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t* Fr = malloc(sizeof(mpfr_t) * (size_t)Ns);
    mpfr_t* Fi = malloc(sizeof(mpfr_t) * (size_t)Ns);
    if (!Fr || !Fi) { free(Fr); free(Fi); return false; }

    mpfr_t pi2, theta, s, cth, zr, zi, rr;
    mpfr_inits2(p, pi2, theta, s, cth, zr, zi, rr, (mpfr_ptr)0);
    mpfr_const_pi(pi2, MPFR_RNDN); mpfr_mul_2ui(pi2, pi2, 1, MPFR_RNDN);
    mpfr_set_d(rr, r, MPFR_RNDN);

    bool ok = true; int filled = 0;
    for (int j = 0; j < Ns && ok; j++) {
        mpfr_mul_si(theta, pi2, j, MPFR_RNDN);
        mpfr_div_si(theta, theta, Ns, MPFR_RNDN);
        mpfr_sin_cos(s, cth, theta, MPFR_RNDN);
        mpfr_fma(zr, rr, cth, Nbase, MPFR_RNDN);    /* zr = N + r cos θ */
        mpfr_mul(zi, rr, s, MPFR_RNDN);             /* zi = r sin θ     */
        mpfr_init2(Fr[j], p); mpfr_init2(Fi[j], p); filled = j + 1;
        Expr* zv = mpfr_zero_p(zi) ? expr_new_mpfr_copy(zr)
                                   : make_complex(expr_new_mpfr_copy(zr),
                                                  expr_new_mpfr_copy(zi));
        if (!ns_eval_complex_mpfr(c, c->body, zv, Fr[j], Fi[j])) ok = false;
    }

    if (ok) {
        mpfr_t accr, acci, ang, wc, ws, t1, t2, scale, ckmag, cmax;
        mpfr_inits2(p, accr, acci, ang, wc, ws, t1, t2, scale, ckmag, cmax, (mpfr_ptr)0);
        mpfr_set_zero(cmax, 1);
        double cM_over = 0.0;
        for (int m = 0; m <= M; m++) {
            int k = ((m % Ns) + Ns) % Ns;
            mpfr_set_zero(accr, 1); mpfr_set_zero(acci, 1);
            for (int j = 0; j < Ns; j++) {
                mpfr_mul_si(ang, pi2, j, MPFR_RNDN);
                mpfr_mul_si(ang, ang, k, MPFR_RNDN);
                mpfr_div_si(ang, ang, Ns, MPFR_RNDN);
                mpfr_neg(ang, ang, MPFR_RNDN);
                mpfr_sin_cos(ws, wc, ang, MPFR_RNDN);
                mpfr_mul(t1, Fr[j], wc, MPFR_RNDN);
                mpfr_mul(t2, Fi[j], ws, MPFR_RNDN);
                mpfr_sub(t1, t1, t2, MPFR_RNDN); mpfr_add(accr, accr, t1, MPFR_RNDN);
                mpfr_mul(t1, Fr[j], ws, MPFR_RNDN);
                mpfr_mul(t2, Fi[j], wc, MPFR_RNDN);
                mpfr_add(t1, t1, t2, MPFR_RNDN); mpfr_add(acci, acci, t1, MPFR_RNDN);
            }
            mpfr_div_si(accr, accr, Ns, MPFR_RNDN);     /* c_k = a_m · r^m */
            mpfr_div_si(acci, acci, Ns, MPFR_RNDN);
            mpfr_hypot(ckmag, accr, acci, MPFR_RNDN);
            if (mpfr_cmp(ckmag, cmax) > 0) mpfr_set(cmax, ckmag, MPFR_RNDN);
            if (m == M) cM_over = mpfr_get_d(ckmag, MPFR_RNDN)
                                  / (mpfr_get_d(cmax, MPFR_RNDN) + 1e-300);
            mpfr_pow_si(scale, rr, -m, MPFR_RNDN);
            mpfr_mul(a_re[m], accr, scale, MPFR_RNDN);
            mpfr_mul(a_im[m], acci, scale, MPFR_RNDN);
        }
        if (!(cM_over <= 1e-4)) ok = false;             /* decay gate */
        mpfr_clears(accr, acci, ang, wc, ws, t1, t2, scale, ckmag, cmax, (mpfr_ptr)0);
    }

    for (int j = 0; j < filled; j++) { mpfr_clear(Fr[j]); mpfr_clear(Fi[j]); }
    free(Fr); free(Fi);
    mpfr_clears(pi2, theta, s, cth, zr, zi, rr, (mpfr_ptr)0);
    return ok;
}

static bool ns_taylor_mpfr_auto(NsCtx* c, const mpfr_t Nbase, double Nd, long bits,
                                int M, mpfr_t* a_re, mpfr_t* a_im) {
    double r = 0.5 * fabs(Nd);
    if (r < 1.0) r = 1.0;
    for (int tries = 0; tries < 8; tries++) {
        if (ns_taylor_mpfr(c, Nbase, r, bits, M, a_re, a_im)) return true;
        r *= 0.5;
        if (r < 1e-3) break;
    }
    return false;
}

static Expr* ns_em_mpfr(NsCtx* c, const char* var, NsOpts* o, long settle) {
    long bits = o->bits;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    int m = ns_em_head(o, settle);

    mpfr_t Hr, Hi, tr, ti, tail_r, tail_i, Iv_r, Iv_i, fN_r, fN_i;
    mpfr_t ddi_m, dipow, bcoeff, coeff, dvr, dvi, cr, ci, Na;
    mpfr_inits2(p, Hr, Hi, tr, ti, tail_r, tail_i, Iv_r, Iv_i, fN_r, fN_i,
                ddi_m, dipow, bcoeff, coeff, dvr, dvi, cr, ci, Na, (mpfr_ptr)0);
    mpfr_set_ui(Hr, 0, MPFR_RNDN); mpfr_set_ui(Hi, 0, MPFR_RNDN);

    bool ok = true;
    int sgn_prev = 0, sgn_last = 0;
    for (int k = 0; k < m; k++) {
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        mpfr_add(Hr, Hr, tr, MPFR_RNDN); mpfr_add(Hi, Hi, ti, MPFR_RNDN);
        sgn_prev = sgn_last; sgn_last = mpfr_sgn(tr);
    }
    /* Oscillatory (sign-alternating) summand: keep to symbolic derivatives — the
     * contour DFT is ill-conditioned on its continuous extension (see machine). */
    bool oscillatory = (m >= 2 && sgn_prev * sgn_last < 0);

    Expr* Nval = ok ? ns_index_value(c, m) : NULL;
    double ddi = 0.0, Nd = 0.0;
    if (ok && (!ns_to_double_real(c->di, &ddi) || ddi == 0.0)) ok = false;
    if (ok) {
        NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; spec.preserve_inexact = false;
        Expr* nn = numericalize(Nval, spec);
        bool inexact; mpfr_t dummy; mpfr_init2(dummy, p);
        ok = nn && get_approx_mpfr(nn, Na, dummy, &inexact);
        expr_free(nn); mpfr_clear(dummy);
        if (ok) Nd = mpfr_get_d(Na, MPFR_RNDN);
    }
    if (ok) mpfr_set_d(ddi_m, ddi, MPFR_RNDN);
    if (ok) ok = ns_eval_complex_mpfr(c, c->body, expr_copy(Nval), fN_r, fN_i);

    /* Tail-integral accuracy scaled to the TARGET precision (not the guarded
     * working precision): reltol must sit just ABOVE the achievable floor so the
     * double-exponential rule's convergence gate trips in ~12–16 levels.  A
     * reltol below the floor would never trip and would refine to a
     * catastrophic node count (the source of the Log-summand slowdown). */
    double digits = numeric_bits_to_digits(o->target_bits);
    double reltol = pow(10.0, -(digits - 2.0));
    if (reltol < 1e-290) reltol = 1e-290;
    int mlev = 14 + (int)(digits / 8.0);
    if (mlev > 22) mlev = 22;

    double abserr = INFINITY; bool integ = false;
    if (ok) {
        integ = dequad_halfline_mpfr(ns_sample_x_mpfr, c, Na, bits, reltol, mlev,
                                     Iv_r, Iv_i, &abserr);
        double relerr = abserr / (ns_l1_d(Iv_r, Iv_i) + 1e-300);
        if (!integ && relerr > 1e-6) ok = false;
    }

    Expr* out = NULL;
    if (ok) {
        mpfr_div(tail_r, Iv_r, ddi_m, MPFR_RNDN);
        mpfr_div(tail_i, Iv_i, ddi_m, MPFR_RNDN);
        mpfr_div_ui(tr, fN_r, 2, MPFR_RNDN); mpfr_add(tail_r, tail_r, tr, MPFR_RNDN);
        mpfr_div_ui(ti, fN_i, 2, MPFR_RNDN); mpfr_add(tail_i, tail_i, ti, MPFR_RNDN);

        /* Hybrid corrections: symbolic D primary (robust, original behaviour),
         * numerical contour derivatives as fallback once the derivative tree
         * balloons (composite summands). */
        int J = ns_em_maxj(o);
        int M = 2 * J - 1;
        mpfr_t* a_re = NULL; mpfr_t* a_im = NULL;   /* contour coeffs, on demand */
        bool a_tried = false, a_ok = false;
        mpfr_set(dipow, ddi_m, MPFR_RNDN);          /* di^1 */
        double prev_mag = INFINITY;
        double goal = (o->prec_goal > 0 ? o->prec_goal : digits - 2);
        double tol = pow(10.0, -goal)
                     * (ns_l1_d(Hr, Hi) + ns_l1_d(tail_r, tail_i) + 1.0);
        ns_bind_clear_temp(symtab_get_def(var));
        Expr* dcur = expr_copy(c->body);
        int cur = 0;
        bool ballooned = false, bail = false;
        int napplied = 0;
        for (int j = 1; j <= J; j++) {
            int ord = 2 * j - 1;
            bool have_corr = false;
            if (!ballooned) {
                while (dcur && cur < ord) {
                    dcur = eval_and_free(ns_mk2("D", dcur, expr_new_symbol(var)));
                    cur++;
                    if (dcur && ns_expr_size(dcur, NS_EM_MAX_NODES + 1) > NS_EM_MAX_NODES) {
                        ballooned = true; break;
                    }
                }
                if (dcur && !ballooned && ns_em_coeff_mpfr(j, bits, coeff)
                    && ns_eval_complex_mpfr(c, dcur, expr_copy(Nval), dvr, dvi)) {
                    ns_bind_clear_temp(symtab_get_def(var));
                    mpfr_mul(cr, coeff, dipow, MPFR_RNDN);   /* coeff·di^ord       */
                    mpfr_mul(ci, cr, dvi, MPFR_RNDN);        /* corr imag          */
                    mpfr_mul(cr, cr, dvr, MPFR_RNDN);        /* corr real          */
                    have_corr = true;
                } else if (!dcur) {
                    break;                          /* derivative vanished: done */
                }
            }
            if (!have_corr) {                       /* symbolic failed/ballooned */
                if (oscillatory) {
                    /* contour unreliable for an oscillatory extension; with no
                     * correction applied the bare integral is too poor — defer
                     * to Wynn. */
                    if (napplied == 0) bail = true;
                    break;
                }
                if (!a_tried) {
                    a_tried = true;
                    a_re = malloc(sizeof(mpfr_t) * (size_t)(M + 1));
                    a_im = malloc(sizeof(mpfr_t) * (size_t)(M + 1));
                    if (a_re && a_im) {
                        for (int i = 0; i <= M; i++) { mpfr_init2(a_re[i], p); mpfr_init2(a_im[i], p); }
                        a_ok = ns_taylor_mpfr_auto(c, Na, Nd, bits, M, a_re, a_im);
                    }
                }
                if (!a_ok) break;                   /* keep corrections so far */
                if (!ns_em_bcoeff_mpfr(j, bits, bcoeff)) break;
                mpfr_mul(cr, bcoeff, dipow, MPFR_RNDN);          /* bcoeff·di^ord  */
                mpfr_mul(ci, cr, a_im[ord], MPFR_RNDN);
                mpfr_mul(cr, cr, a_re[ord], MPFR_RNDN);
            }
            double mag = ns_l1_d(cr, ci);
            if (j > 1 && mag > prev_mag) break;     /* asymptotic minimum */
            mpfr_sub(tail_r, tail_r, cr, MPFR_RNDN);
            mpfr_sub(tail_i, tail_i, ci, MPFR_RNDN);
            napplied++;
            prev_mag = mag;
            if (mag < tol) break;
            mpfr_mul(dipow, dipow, ddi_m, MPFR_RNDN);
            mpfr_mul(dipow, dipow, ddi_m, MPFR_RNDN);
        }
        expr_free(dcur);
        if (a_re) { for (int i = 0; i <= M; i++) { mpfr_clear(a_re[i]); mpfr_clear(a_im[i]); } }
        free(a_re); free(a_im);

        if (!bail) {                                /* else NULL -> Wynn fallback */
            mpfr_add(tr, Hr, tail_r, MPFR_RNDN);
            mpfr_add(ti, Hi, tail_i, MPFR_RNDN);
            /* Chop the sub-tolerance imaginary residual on a real-valued summand. */
            {
                double thr = pow(10.0, -(digits - 2.0))
                             * (fabs(mpfr_get_d(tr, MPFR_RNDN)) + 1.0);
                if (fabs(mpfr_get_d(ti, MPFR_RNDN)) < thr) mpfr_set_zero(ti, 1);
            }
            out = ns_from_complex_mpfr(tr, ti);
        }
    }

    expr_free(Nval);
    mpfr_clears(Hr, Hi, tr, ti, tail_r, tail_i, Iv_r, Iv_i, fN_r, fN_i,
                ddi_m, dipow, bcoeff, coeff, dvr, dvi, cr, ci, Na, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Method: AlternatingSigns (Cohen–Villegas–Zagier, 2000)             *
 *                                                                     *
 *  For S = Σ_{k≥0} (-1)^k a_k the algorithm builds Chebyshev weights   *
 *  d_n = ((3+√8)^n + (3+√8)^{-n})/2 and a single pass over n terms      *
 *  yields ≈ n·log2(3+√8) ≈ 2.54 n bits.  Linear in the terms, so the   *
 *  a_k may be complex; the sign-normalised a_k = (-1)^k·term_k is used. *
 * ------------------------------------------------------------------ */

/* Number of CVZ terms for a target precision (bits). */
static int ns_cvz_terms(const NsOpts* o) {
    /* Scale the term count to the TARGET precision, not the guarded working
     * precision: CVZ on a non-completely-monotone (peaked/sawtooth) alternating
     * series is unstable with excess terms — its Chebyshev weights amplify the
     * deviation — so using guard-inflated bits here would lose accuracy. */
    double bits = o->prec_mpfr ? (double)o->target_bits : 53.0;
    int n = (int)ceil(bits / 2.54) + 6;
    if (n < 4) n = 4;
    return n;
}

static Expr* ns_cvz_machine(NsCtx* c, int n) {
    double d = pow(3.0 + sqrt(8.0), (double)n);
    d = (d + 1.0 / d) / 2.0;
    double b = -1.0, cc = -d;
    double _Complex s = 0.0;
    for (int k = 0; k < n; k++) {
        double _Complex term;
        if (!ns_term_machine(c, k, &term)) return NULL;
        double _Complex ak = ((k & 1) ? -1.0 : 1.0) * term;
        cc = b - cc;
        s += cc * ak;
        b = (double)(k + n) * (double)(k - n) * b / (((double)k + 0.5) * (double)(k + 1));
    }
    return ns_from_complex_d(s / d);
}

#ifdef USE_MPFR
static Expr* ns_cvz_mpfr(NsCtx* c, int n, long bits) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t d, b, cc, sr, si, tr, ti, akr, aki, tmp;
    mpfr_inits2(p, d, b, cc, sr, si, tr, ti, akr, aki, tmp, (mpfr_ptr)0);

    mpfr_set_ui(tmp, 8, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 3, MPFR_RNDN);          /* 3 + √8 */
    mpfr_pow_ui(d, tmp, (unsigned long)n, MPFR_RNDN);
    mpfr_ui_div(tmp, 1, d, MPFR_RNDN);
    mpfr_add(d, d, tmp, MPFR_RNDN);
    mpfr_div_2ui(d, d, 1, MPFR_RNDN);             /* d_n */
    mpfr_set_si(b, -1, MPFR_RNDN);
    mpfr_neg(cc, d, MPFR_RNDN);
    mpfr_set_ui(sr, 0, MPFR_RNDN); mpfr_set_ui(si, 0, MPFR_RNDN);

    bool ok = true;
    for (int k = 0; k < n; k++) {
        if (!ns_term_mpfr(c, k, tr, ti)) { ok = false; break; }
        if (k & 1) { mpfr_neg(akr, tr, MPFR_RNDN); mpfr_neg(aki, ti, MPFR_RNDN); }
        else       { mpfr_set(akr, tr, MPFR_RNDN); mpfr_set(aki, ti, MPFR_RNDN); }
        mpfr_sub(cc, b, cc, MPFR_RNDN);           /* c = b - c */
        mpfr_mul(tmp, cc, akr, MPFR_RNDN); mpfr_add(sr, sr, tmp, MPFR_RNDN);
        mpfr_mul(tmp, cc, aki, MPFR_RNDN); mpfr_add(si, si, tmp, MPFR_RNDN);
        long num = (long)(k + n) * (long)(k - n);  /* (k+n)(k-n) */
        mpfr_mul_si(b, b, num, MPFR_RNDN);
        mpfr_mul_2ui(b, b, 1, MPFR_RNDN);          /* /((k+1/2)(k+1)) = ·2/((2k+1)(k+1)) */
        mpfr_div_ui(b, b, (unsigned long)(2 * k + 1), MPFR_RNDN);
        mpfr_div_ui(b, b, (unsigned long)(k + 1), MPFR_RNDN);
    }

    Expr* out = NULL;
    if (ok) {
        mpfr_div(sr, sr, d, MPFR_RNDN);
        mpfr_div(si, si, d, MPFR_RNDN);
        out = ns_from_complex_mpfr(sr, si);
    }
    mpfr_clears(d, b, cc, sr, si, tr, ti, akr, aki, tmp, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Automatic method selection + convergence verification              *
 * ------------------------------------------------------------------ */

#define NS_PROBE_MAX 18

typedef struct { double _Complex t[NS_PROBE_MAX]; int n; bool ok; } NsProbe;

/* Full summand profile: head terms + a geometric far-tail magnitude ladder.
 * The ladder lets the classifier see structure (a peak, the decay onset, true
 * divergence) far beyond the first NS_PROBE_MAX terms — the head window alone
 * mishandles summands that peak or settle late (e.g. 1/(1+(k-20)^2)). */
typedef struct {
    NsProbe head;
    double  lmag[NS_LADDER_MAX];   /* |a_k| on the ladder k = 4,8,16,…    */
    long    lidx[NS_LADDER_MAX];
    int     nl;
    bool    ladder_ok;             /* every ladder sample was numeric     */
    bool    rising_far;            /* |a_k| still clearly growing far out  */
    long    settle;                /* decay-onset index (peak turned over), or -1 */
} NsProfile;

/* Sample the first few terms at machine precision (used to classify the sum). */
static NsProbe ns_probe(NsCtx* c, long count) {
    NsProbe pr; pr.n = 0; pr.ok = true;
    for (int k = 0; k < NS_PROBE_MAX; k++) {
        if (count >= 0 && k >= count) break;
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { pr.ok = false; break; }
        pr.t[pr.n++] = t;
    }
    return pr;
}

/* Build the head probe + far-tail magnitude ladder and derive the divergence
 * and decay-onset signals.  The ladder costs ~20 term evaluations and runs only
 * when verifying convergence or selecting the method automatically. */
static void ns_build_profile(NsCtx* c, long count, NsProfile* p) {
    p->head = ns_probe(c, count);
    p->nl = 0; p->ladder_ok = false; p->rising_far = false; p->settle = -1;
    /* The far-tail ladder samples deep indices.  For a nested numeric summand
     * with a dependent inner bound {k,1,n} that means a large inner computation
     * per sample for no benefit (the settle index is irrelevant there), so skip
     * the ladder for nested summands — but still allow Euler–Maclaurin, which is
     * valid for an independent inner bound {k,1,Infinity} and falls back on its
     * own otherwise. */
    if (ns_body_is_blackbox(c->body)) return;
    /* If the head probe already shows a monotone-decreasing tail the summand has
     * settled (no late peak, not divergent), so skip the far-tail ladder — this
     * avoids ~16 extra summand evaluations per profile on the common case (and,
     * via the shared evaluator's per-call allocations, keeps memory churn down).
     * A peaked/rising head (e.g. 1/(1+(k-20)^2), or a divergent 2^k) is NOT
     * monotone here, so the ladder still runs for the cases that need it. */
    {
        const NsProbe* h = &p->head;
        int from = h->n / 2;
        bool mono = (h->n - from >= 3);
        for (int kk = from + 1; mono && kk < h->n; kk++)
            if (cabs(h->t[kk]) > cabs(h->t[kk - 1]) * 1.0000001) mono = false;
        if (mono) return;
    }
    p->ladder_ok = true;
    long k = 4;
    for (int j = 0; j < NS_LADDER_MAX; j++, k *= 2) {
        if (count >= 0 && k >= count) break;
        double _Complex t;
        if (!ns_term_machine(c, k, &t)) { p->ladder_ok = false; break; }
        p->lidx[p->nl] = k; p->lmag[p->nl] = cabs(t); p->nl++;
    }
    if (p->nl >= 3) {
        double a = p->lmag[p->nl - 3], b = p->lmag[p->nl - 2], cc = p->lmag[p->nl - 1];
        p->rising_far = (cc > 1.1 * b) && (b >= a);   /* sustained far-tail growth */
        int peak = 0;
        for (int j = 1; j < p->nl; j++) if (p->lmag[j] > p->lmag[peak]) peak = j;
        if (peak < p->nl - 1) p->settle = p->lidx[peak + 1];   /* turned over -> known */
    }
}

/* Strictly real-and-sign-alternating terms? */
static bool ns_is_alternating(const NsProbe* pr) {
    if (pr->n < 4) return false;
    for (int k = 0; k < pr->n; k++) {
        double re = creal(pr->t[k]), im = cimag(pr->t[k]);
        if (fabs(im) > 1e-9 * (fabs(re) + 1e-300)) return false;
        if (re == 0.0) return false;
    }
    for (int k = 1; k < pr->n; k++)
        if (creal(pr->t[k]) * creal(pr->t[k - 1]) >= 0.0) return false;
    return true;
}

/* |term| non-increasing from index `from`? */
static bool ns_mag_decreasing(const NsProbe* pr, int from) {
    for (int k = from + 1; k < pr->n; k++)
        if (cabs(pr->t[k]) > cabs(pr->t[k - 1]) * 1.0000001) return false;
    return pr->n - from >= 3;
}

/* Clear divergence: the term magnitude is still growing far into the tail
 * (sampled on the geometric ladder).  Peaked-then-decaying and harmonic-like
 * sums turn over / decay and are NOT flagged — preserving NSum's deliberate
 * blindness to ratios → 1, while no longer mistaking a late peak for growth. */
static bool ns_diverges(const NsProfile* p) {
    return p->rising_far;
}

/* Resolve Automatic to a concrete method from the profile. */
static const char* ns_choose_method(const NsProfile* p) {
    const NsProbe* pr = &p->head;
    /* Strictly sign-alternating with magnitude monotone from the start ->
     * Cohen–Villegas–Zagier.  (A peaked alternating series like (-5)^i/i! is
     * NOT in this class: CVZ assumes |a_k| decreasing from k=0, so peaked or
     * sawtooth-magnitude alternating series fall through to Euler–Maclaurin /
     * Wynn instead.) */
    if (ns_is_alternating(pr) && ns_mag_decreasing(pr, 0)) return SYM_AlternatingSigns;
    if (!pr->ok || pr->n < 4) return SYM_WynnEpsilon;
    /* far tail located a peak then decays -> Euler–Maclaurin with adaptive base */
    if (p->ladder_ok && !p->rising_far && p->nl >= 2
        && p->lmag[p->nl - 1] < p->lmag[p->nl - 2])
        return SYM_EulerMaclaurin;
    if (ns_mag_decreasing(pr, pr->n / 2)) return SYM_EulerMaclaurin;   /* head tail monotone */
    return SYM_WynnEpsilon;
}

/* ------------------------------------------------------------------ *
 *  Single-variable driver                                             *
 * ------------------------------------------------------------------ */

/* Direct-evaluation cap for finite sums (above this, extrapolate the tail). */
static long ns_direct_cap(const NsOpts* o) {
    long cap = (long)(o->nsum_terms < 0 ? 0 : o->nsum_terms) + ns_seq_len(o) + 4;
    if (cap < 64) cap = 64;
    return cap;
}

/* Direct machine/MPFR sum dispatch over `count` terms. */
static Expr* ns_direct(NsCtx* c, long count, const NsOpts* o) {
    if (count < 0) count = 0;
#ifdef USE_MPFR
    return o->prec_mpfr ? ns_direct_mpfr(c, count, o->bits) : ns_direct_machine(c, count);
#else
    return ns_direct_machine(c, count);
#endif
}

/* Sum the infinite series Σ_{k≥0} f(imin + k·di) by method dispatch.  When
 * `do_verify`, a clearly divergent series returns the ComplexInfinity symbol. */
static Expr* ns_sum_infinite(NsCtx* c, const char* var, NsOpts* o, bool do_verify) {
    const char* method = o->method;
    bool is_auto = (o->method == SYM_Automatic);
    long settle = -1;
    if (do_verify || method == SYM_Automatic) {
        NsProfile prof;
        ns_build_profile(c, -1, &prof);
        if (do_verify && ns_diverges(&prof)) {
            ns_warn("div", "the sum does not appear to converge");
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        if (method == SYM_Automatic) method = ns_choose_method(&prof);
        settle = prof.settle;
    }

    Expr* out = NULL;
    bool used = false;
    if (method == SYM_AlternatingSigns) {
        int n = ns_cvz_terms(o);
#ifdef USE_MPFR
        out = o->prec_mpfr ? ns_cvz_mpfr(c, n, o->bits) : ns_cvz_machine(c, n);
#else
        out = ns_cvz_machine(c, n);
#endif
        used = (out != NULL);
    } else if (method == SYM_EulerMaclaurin) {
#ifdef USE_MPFR
        out = o->prec_mpfr ? ns_em_mpfr(c, var, o, settle) : ns_em_machine(c, var, o, settle);
#else
        out = ns_em_machine(c, var, o, settle);
#endif
        used = (out != NULL);                           /* else fall back to Wynn */
    } else if (method == SYM_Levin) {
#ifdef USE_MPFR
        out = o->prec_mpfr ? ns_levin_mpfr(c, -1, o) : ns_levin_machine(c, -1, o);
#else
        out = ns_levin_machine(c, -1, o);
#endif
        used = (out != NULL);                           /* else fall back to Wynn */
    }
    if (!used) {
#ifdef USE_MPFR
        out = o->prec_mpfr ? ns_wynn_mpfr(c, -1, o) : ns_wynn_machine(c, -1, o);
#else
        out = ns_wynn_machine(c, -1, o);
#endif
        /* Automatic last resort: when Wynn does not converge, give Levin's
         * u-transform a shot.  Purely additive — an existing Wynn result is
         * never replaced, so previously-converging Automatic sums are
         * unaffected. */
        if (!out && is_auto) {
#ifdef USE_MPFR
            out = o->prec_mpfr ? ns_levin_mpfr(c, -1, o) : ns_levin_machine(c, -1, o);
#else
            out = ns_levin_machine(c, -1, o);
#endif
        }
    }
    return out;
}

/* Run NSum over one iterator spec already split into (var, imin, imax, di),
 * with imin/imax/di evaluated to numbers and `infinite`/`count` describing the
 * range.  A large finite sum is computed as the difference of two infinite
 * tails, Σ_{imin}^∞ − Σ_{imax+di}^∞, when the summand decays. */
static Expr* ns_run_single(Expr* body, const char* var, Expr* imin, Expr* imax,
                           Expr* di, bool infinite, long count, NsOpts* o) {
    NsBind bind; ns_bind_snapshot(&bind, var);
    NsCtx ctx;
    ctx.body = body; ctx.imin = imin; ctx.di = di; ctx.bind = &bind;
#ifdef USE_MPFR
    if (o->prec_mpfr) { ctx.spec.mode = NUMERIC_MODE_MPFR; ctx.spec.bits = o->bits; ctx.spec.preserve_inexact = false; }
    else ctx.spec = numeric_machine_spec();
#else
    ctx.spec = numeric_machine_spec();
#endif

    Expr* out = NULL;
    if (infinite) {
        out = ns_sum_infinite(&ctx, var, o, o->verify);
    } else if (count <= ns_direct_cap(o)) {
        out = ns_direct(&ctx, count, o);
    } else {
        /* Large finite sum.  If the summand decays, use the difference of two
         * infinite tails; otherwise fall back to direct summation. */
        NsProfile prof;
        ns_build_profile(&ctx, -1, &prof);
        bool decays = prof.head.ok && !prof.rising_far
                      && ((prof.ladder_ok && prof.nl >= 2
                           && prof.lmag[prof.nl - 1] < prof.lmag[prof.nl - 2])
                          || ns_mag_decreasing(&prof.head, prof.head.n / 2));
        if (decays) {
            Expr* head_sum = ns_sum_infinite(&ctx, var, o, false);
            /* Second tail starts at imax + di. */
            Expr* imin2 = eval_and_free(ns_mk2("Plus", expr_copy(imax), expr_copy(di)));
            ctx.imin = imin2;
            Expr* tail_sum = head_sum ? ns_sum_infinite(&ctx, var, o, false) : NULL;
            ctx.imin = imin;
            expr_free(imin2);
            if (head_sum && tail_sum) {
                out = eval_and_free(ns_mk2("Subtract", head_sum, tail_sum));
            } else {
                expr_free(head_sum); expr_free(tail_sum);
            }
        } else if (count <= 2000000) {
            out = ns_direct(&ctx, count, o);
        } else {
            ns_warn("nsum", "large finite sum of non-decaying terms; supply a "
                            "closed form via Sum or fewer terms");
        }
    }

    ns_bind_restore(&bind);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Iterator-spec parsing                                              *
 * ------------------------------------------------------------------ */

/* A spec is a List {var, imax} | {var, imin, imax} | {var, imin, imax, di}. */
static bool ns_is_spec(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    return n >= 2 && n <= 4 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nsum(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && ns_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    if (pos_end < 2) return NULL;                       /* need f + >=1 spec */
    Expr* body = res->data.function.args[0];
    /* positional args [1 .. pos_end) must all be iterator specs. */
    for (size_t i = 1; i < pos_end; i++)
        if (!ns_is_spec(res->data.function.args[i])) return NULL;

    NsOpts o;
    o.method = SYM_Automatic;
    o.nsum_terms = NS_DEF_NSUM_TERMS;
    o.nsum_terms_user = false;
    o.extra_terms = -1;
    o.wynn = NS_DEF_WYNN;
    o.levin_variant = SEQACCEL_LEVIN_U;
    o.verify = true;
    o.prec_mpfr = false;
    o.bits = 0;
    o.target_bits = 0;
    o.acc_goal = -1.0;
    o.prec_goal = -1.0;
    for (size_t i = pos_end; i < argc; i++)
        if (!ns_apply_option(res->data.function.args[i], &o)) return NULL;

    /* Guard precision: the Euler–Maclaurin tail integral and contour derivatives
     * sample the summand at float arguments where a near-1 form (e.g.
     * Log[1+1/x^2]) cancels at fixed precision, and the exp-sinh weight then
     * amplifies that loss.  Run the internal MPFR work at ~2x the target so the
     * cancellation is absorbed, then round the result back to the target.  (At
     * the integer head points the summand stays exact, so this only matters off
     * the integers.)  Reaching the target needs ~2x bits of headroom. */
    if (o.prec_mpfr && o.bits > 0) {
        o.target_bits = o.bits;
        o.bits = 2 * o.target_bits + 32;
    }

    /* For now Automatic / EulerMaclaurin / AlternatingSigns all route through
     * WynnEpsilon (the dedicated methods are layered on next). */

    size_t nspecs = pos_end - 1;
    Expr* spec0 = res->data.function.args[1];

    /* Multidimensional: the summand of the outer sum is an inner NSum over the
     * remaining specs (HoldAll + index localisation makes dependent inner
     * bounds work). Build it once and treat as the single-variable body. */
    Expr* inner_body = NULL;   /* owned when built */
    Expr* eff_body = body;     /* borrowed unless inner_body is built */
    if (nspecs >= 2) {
        size_t ninner = (argc - pos_end) + 1 /*body*/ + (nspecs - 1) /*specs[1..]*/;
        Expr** v = malloc(sizeof(Expr*) * ninner);
        size_t w = 0;
        v[w++] = expr_copy(body);
        for (size_t i = 2; i < pos_end; i++) v[w++] = expr_copy(res->data.function.args[i]);
        for (size_t i = pos_end; i < argc; i++) v[w++] = expr_copy(res->data.function.args[i]);
        inner_body = expr_new_function(expr_new_symbol(SYM_NSum), v, w);
        free(v);
        eff_body = inner_body;
    }

    /* Parse spec0 into (var, imin, imax, di). */
    const char* var = spec0->data.function.args[0]->data.symbol.name;
    size_t sn = spec0->data.function.arg_count;
    Expr *imin_raw, *imax_raw, *di_raw;
    if (sn == 2) { imin_raw = NULL; imax_raw = spec0->data.function.args[1]; di_raw = NULL; }
    else if (sn == 3) { imin_raw = spec0->data.function.args[1]; imax_raw = spec0->data.function.args[2]; di_raw = NULL; }
    else { imin_raw = spec0->data.function.args[1]; imax_raw = spec0->data.function.args[2]; di_raw = spec0->data.function.args[3]; }

    Expr* imin = imin_raw ? eval_and_free(expr_copy(imin_raw)) : expr_new_integer(1);
    Expr* di   = di_raw   ? eval_and_free(expr_copy(di_raw))   : expr_new_integer(1);
    Expr* imax_e = eval_and_free(expr_copy(imax_raw));

    bool infinite = ns_is_infinite(imax_e);
    long count = -1;
    bool bad = false;
    if (!infinite) {
        /* count = floor((imax - imin)/di) + 1 */
        double dmin, dmax, ddi;
        if (ns_to_double_real(imin, &dmin) && ns_to_double_real(imax_e, &dmax)
            && ns_to_double_real(di, &ddi) && ddi != 0.0) {
            double n = floor((dmax - dmin) / ddi + 1e-9) + 1.0;
            if (n < 0) n = 0;
            if (n > 9.0e18) bad = true; else count = (long)n;
        } else {
            bad = true;     /* non-numeric finite bound: stay unevaluated */
        }
    }

    Expr* out = NULL;
    if (!bad)
        out = ns_run_single(eff_body, var, imin, imax_e, di, infinite, count, &o);

    /* Round the guarded internal result back to the requested precision. */
    if (out && o.prec_mpfr && o.target_bits > 0) {
        NumericSpec ts; ts.mode = NUMERIC_MODE_MPFR; ts.bits = o.target_bits;
        ts.preserve_inexact = false;
        Expr* rounded = numericalize(out, ts);
        if (rounded) { expr_free(out); out = rounded; }
    }

    expr_free(imin);
    expr_free(di);
    expr_free(imax_e);
    expr_free(inner_body);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nsum_init(void) {
    symtab_add_builtin("NSum", builtin_nsum);
    /* HoldAll: the summand and the iterator specs must not be pre-evaluated;
     * the index is Block-localised internally. Not Listable. */
    symtab_get_def("NSum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}

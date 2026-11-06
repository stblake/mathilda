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
#include "../modular.h"

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
#include "compile/autocompile.h"   /* machine fast path for the summand */
#include "core.h"         /* tc_register_unwind_counter */
#include "eval.h"
#include "nc_accuracy.h"  /* shared AccuracyGoal/PrecisionGoal handling */
#include "numeric.h"
#include "seqaccel.h"
#include "dequad.h"
#include "denint.h"
#include "gkadapt.h"
#include "sym_intern.h"   /* the scaled-Hankel placeholder heads */
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

/* > 0 while an internal second estimate runs whose failure is not the user's
 * (ns_wynn_nonsmooth): its messages are muted and kept out of Check[]. */
static int g_ns_quiet = 0;

static void ns_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(g_ns_quiet > 0, "NSum", tag, fmt, ap);
    va_end(ap);
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
    int         token;     /* binding frame (modular.h mth_bind_push) */
    bool        valid;
} NsBind;

static void ns_bind_snapshot(NsBind* b, const char* name) {
    b->name = name;
    /* Unwind-safe: a TimeConstrained timeout drains the frame. */
    b->token = mth_bind_push_name(name, MTH_BIND_ATTRS);
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
    mth_bind_pop(b->token);   /* frees the temporary value */
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
    /* Machine fast path for the summand, or NULL when the body is outside the
     * compilable subset.  Every acceleration method here — direct, Wynn, Levin,
     * Euler–Maclaurin, CVZ, the far-tail ladder and the profiler — draws its
     * terms through ns_term_machine, so compiling once at the top covers all of
     * them.  Only built for a machine-precision request: the MPFR path needs
     * more than 53 bits per term and must keep the interpreter. */
    AutoCompiled* ac;
    /* Compiled MPFR fast path for a high-WorkingPrecision sum: the compiled body
     * evaluates the summand at the working precision, replacing the per-term
     * tree-walk.  Same faithful-degrade contract as `ac` (NULL result ⇒ interpret
     * that term). */
    AutoCompiled* ac_prec;
    /* Complex-input variant, for the Euler–Maclaurin derivative step, which
     * evaluates the summand on a Cauchy contour (Complex[MPFR,MPFR] points). */
    AutoCompiled* ac_prec_z;
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
    /* Every summand sample runs through here — head terms, the far-tail ladder,
     * the Euler–Maclaurin quadrature/contour nodes, the continuity probe.  These
     * are internal trial evaluations and their per-sample messages are noise: a
     * summand like 1/Prime[x] sampled at a continuous x would otherwise emit one
     * Prime::intpp per node.  Wolfram's NSum is silent about such probes; mute
     * arithmetic warnings for the duration of the sample (the counter is nestable
     * so a genuinely user-facing message outside NSum is unaffected). */
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, c->spec) : NULL;
    arith_warnings_mute_pop();
    expr_free(raw);
    return num;
}

/* Evaluate the summand at the index value (consumed). */
static Expr* ns_eval_at(NsCtx* c, Expr* value) {
#ifdef USE_MPFR
    /* Compiled MPFR fast path: the index value is already a numeric Expr and the
     * compiled body returns the summand at the working precision, or NULL to
     * interpret this term. */
    if (c->ac_prec && c->spec.mode == NUMERIC_MODE_MPFR) {
        Expr* r = autocompiled_eval_mpfr(c->ac_prec, (const Expr* const*)&value);
        if (r) { expr_free(value); return r; }
    }
#endif
    return ns_eval_expr_at(c, c->body, value);
}

/* The summand at a real machine point via the compiled program, or false to
 * mean "use the interpreter for this one".  A failure is never an error: the
 * compiled program declines exactly where the interpreter would produce
 * something a double cannot carry (an exact value, a pole, a complex branch),
 * and the caller simply takes the slow path there. */
static bool ns_compiled_at(const NsCtx* c, double x, double _Complex* out) {
    if (!c->ac) return false;
    double _Complex v;
    if (!autocompiled_eval_complex(c->ac, &x, &v)) return false;
    if (!isfinite(creal(v)) || !isfinite(cimag(v))) return false;
    *out = v;
    return true;
}

/* The index value x_k = imin + k·di as a plain double, for the compiled path.
 * Returns false when either bound is not a real machine number (a complex or
 * symbolic step), which sends the term back to the interpreter. */
static bool ns_index_double(const NsCtx* c, long k, double* out) {
    double lo, step;
    if (!ns_to_double_real(c->imin, &lo) || !ns_to_double_real(c->di, &step)) return false;
    *out = lo + (double)k * step;
    return isfinite(*out);
}

/* term_k as a machine complex. */
static bool ns_term_machine(NsCtx* c, long k, double _Complex* out) {
    double idx;
    if (c->ac && ns_index_double(c, k, &idx) && ns_compiled_at(c, idx, out)) return true;
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

/* AccuracyGoal / PrecisionGoal parsing is shared (accepts MachinePrecision and
 * distinguishes Automatic from Infinity — see nc_accuracy.h). */
static bool ns_parse_goal(Expr* v, double* out) {
    return nc_parse_goal(v, out);
}

/* Working precision in decimal digits for the active mode. */
static double ns_wp_digits(const NsOpts* o) {
#ifdef USE_MPFR
    if (o->prec_mpfr) {
        long b = o->target_bits > 0 ? o->target_bits : o->bits;
        if (b > 0) return numeric_bits_to_digits(b);
    }
#endif
    return NUMERIC_MACHINE_PRECISION_DIGITS;
}

/* True unless both goals are Infinity (accuracy not used as a criterion). */
static bool ns_goal_active(const NsOpts* o) {
    return !(isinf(o->acc_goal) && isinf(o->prec_goal));
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

/* Machine path. count<0 means infinite; otherwise terms beyond `count` are 0.
 * With err_out non-NULL the extrapolation's own error estimate is stored there
 * and the goal warning is left to the caller (which may know the estimate is
 * not to be trusted, see ns_sum_infinite). */
static Expr* ns_wynn_machine(NsCtx* c, long count, const NsOpts* o, double* err_out) {
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
    if (err_out) *err_out = step;
    else if (ns_goal_active(o)) {
        double tol = nc_combined_tol(o->acc_goal, o->prec_goal, cabs(result), ns_wp_digits(o));
        if (step > tol) nc_warn_goal("NSum", step, tol);
    }
    return ns_from_complex_d(result);
}

#ifdef USE_MPFR
static Expr* ns_wynn_mpfr(NsCtx* c, long count, const NsOpts* o, double* err_out) {
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
        if (got && finite && ns_accept(ns_l1_d(rr, ri), step, maxsample)) {
            if (err_out) *err_out = step;
            else if (ns_goal_active(o)) {
                double tol = nc_combined_tol(o->acc_goal, o->prec_goal, ns_l1_d(rr, ri), ns_wp_digits(o));
                if (step > tol) nc_warn_goal("NSum", step, tol);
            }
            out = ns_from_complex_mpfr(rr, ri);
        } else
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
    if (ns_goal_active(o)) {
        double tol = nc_combined_tol(o->acc_goal, o->prec_goal, cabs(result), ns_wp_digits(o));
        if (step > tol) nc_warn_goal("NSum", step, tol);
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
        if (got && finite && ns_accept(ns_l1_d(rr, ri), step, maxsample)) {
            if (ns_goal_active(o)) {
                double tol = nc_combined_tol(o->acc_goal, o->prec_goal, ns_l1_d(rr, ri), ns_wp_digits(o));
                if (step > tol) nc_warn_goal("NSum", step, tol);
            }
            out = ns_from_complex_mpfr(rr, ri);
        } else
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

/* Evaluate `e` at the index value `value` (consumed) to a machine complex.
 *
 * The compiled program covers this too, not just the integer-indexed terms: the
 * Euler–Maclaurin correction samples the SAME summand at continuous real x
 * (through ns_sample_x_machine below, once per quadrature node), and that is
 * where most of the work in an EM sum actually goes.  It applies only when the
 * expression being sampled IS the summand — the derivative expressions this
 * function is also called with are different trees and are not compiled. */
static bool ns_eval_complex_machine(NsCtx* c, Expr* e, Expr* value, double _Complex* out) {
    double x;
    if (c->ac && e == c->body && ns_to_double_real(value, &x)
        && ns_compiled_at(c, x, out)) {
        expr_free(value);
        return true;
    }
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
 *  Tail integral for algebraically decaying summands (machine)        *
 *                                                                     *
 *  A summand f(x) ~ C x^-p with p just above 1 (1/(n^Pi+...)^(1/3),    *
 *  p = Pi/3) leaves ∫_X^∞ f dx ~ C X^(1-p)/(p-1) far from negligible   *
 *  at every X a double can hold: reaching 1e-12 relative needs         *
 *  X ~ 10^(12/(p-1)), i.e. 10^255 at p - 1 = 0.047, and the summand     *
 *  itself overflows (x^Pi) near 10^98.  The exp-sinh rule then loses    *
 *  its far tail to non-finite samples while its successive levels still *
 *  agree with each other, so it reports convergence on an integral that *
 *  is short by the whole unreached remainder (5e-5 relative here).      *
 *                                                                     *
 *  The tail is therefore split at a point X where the power law has     *
 *  been verified:                                                       *
 *      ∫_N^∞ f = ∫_{ln N}^{ln X} f(e^u) e^u du  +  X f(X) / (p(X) - 1),  *
 *  the first piece by adaptive Gauss–Kronrod on the logarithmic         *
 *  abscissa (the integrand there decays like e^{-(p-1)u}: smooth and    *
 *  short), the second analytically from the local exponent              *
 *      p(x) = -log2( f(2x) / f(x) ),                                    *
 *  which is exact for a pure power and converges to the true exponent   *
 *  as the lower-order terms x^-q die away.  X is the largest point of a *
 *  ladder N·10^(8j) at which f is a finite, normal, phase-stable power   *
 *  law; the same split taken one rung lower gives an independent        *
 *  estimate, and their difference is the reported error — so a          *
 *  remainder that is still moving (a slowly varying factor such as a     *
 *  log, or a correction x^-q with tiny q) shows up in the error rather   *
 *  than being silently absorbed.                                        *
 *                                                                     *
 *  Returns false when the summand is not an algebraic power law at      *
 *  large x (exponential decay underflows the ladder; oscillatory or     *
 *  logarithmic factors make p unstable; p <= 1 is not summable this way)*
 *  and the caller keeps the exp-sinh rule.                              *
 * ------------------------------------------------------------------ */

/* g(u) = f(e^u)·e^u — the tail integrand on the logarithmic abscissa. */
static bool ns_sample_logx_machine(void* vc, double u, double _Complex* out) {
    double x = exp(u);
    double _Complex f;
    if (!isfinite(x) || !ns_sample_x_machine(vc, x, &f)) return false;
    *out = f * x;
    return true;
}

/* A sample that carries full relative precision: finite and comfortably
 * normal.  An exact 0 is not "live" — it is either a genuine underflow or an
 * intermediate overflow (1/Infinity^(1/3)), and the two cannot be told apart. */
static bool ns_live_sample(double _Complex f) {
    double a = cabs(f);
    return isfinite(a) && a > 1e-280;
}

#define NS_ALG_LADDER_MAX 38          /* rungs N·10^(8j), j = 1..; < 10^300    */

/* Walk the ladder x_j = N·10^(8j).  Sets *far_mass to x·|f(x)| at the last
 * live rung (0 if none) — for any slowly varying decay, a lower estimate of
 * the integral mass beyond the reach of double-range quadrature.  Returns true
 * (with the split tail integral in *out, its error in *abserr) when the last
 * two live rungs show a settled, summable power law. */
static bool ns_tail_integral_algebraic_machine(NsCtx* c, double Nd,
                                               double _Complex* out, double* abserr,
                                               double* far_mass) {
    *far_mass = 0.0;
    if (!(Nd > 0.0)) return false;
    double X[2] = {0, 0}, P[2] = {0, 0};
    double _Complex F[2] = {0, 0};
    int n = 0;                         /* consecutive power-law rungs ending here */
    double x = Nd;
    for (int j = 0; j < NS_ALG_LADDER_MAX; j++) {
        x *= 1e8;
        if (!isfinite(2.0 * x)) break;
        double _Complex f1, f2;
        if (!ns_sample_x_machine(c, x, &f1) || !ns_live_sample(f1)) break;
        *far_mass = x * cabs(f1);
        /* Local exponent p(x) = -log2(f(2x)/f(x)): exact for a pure power,
         * whose ratio is real positive (an oscillating x^(-a-ib) or
         * sin(x)/x^2 is not). */
        bool pw = ns_sample_x_machine(c, 2.0 * x, &f2) && ns_live_sample(f2);
        double _Complex r = pw ? f2 / f1 : 0.0;
        pw = pw && creal(r) > 0.0 && fabs(cimag(r)) <= 1e-10 * creal(r);
        if (!pw) { n = 0; continue; }
        X[0] = X[1]; P[0] = P[1]; F[0] = F[1];
        X[1] = x; P[1] = -log2(creal(r)); F[1] = f1;
        n++;
    }
    if (n < 2) return false;
    /* A power law: the exponent has settled (tolerance far above the ~1e-15
     * roundoff of p, far below any visible drift) and is summable. */
    if (!(fabs(P[1] - P[0]) <= 1e-8 * fabs(P[1]))) return false;
    if (!(P[1] > 1.0 + 1e-6) || !(P[0] > 1.0 + 1e-6)) return false;
    int J = 1;

    double uN = log(Nd), u1 = log(X[J - 1]), u2 = log(X[J]);
    GkResult q1 = gk_integrate_machine(ns_sample_logx_machine, c, uN, u1,
                                       1e-300, 1e-14, 400, 0, false);
    if (q1.status == GK_NONNUMERIC) return false;
    GkResult q2 = gk_integrate_machine(ns_sample_logx_machine, c, u1, u2,
                                       1e-300, 1e-14, 400, 0, false);
    if (q2.status == GK_NONNUMERIC) return false;

    double _Complex R1 = X[J - 1] * F[J - 1] / (P[J - 1] - 1.0);
    double _Complex R2 = X[J] * F[J] / (P[J] - 1.0);
    *out = q1.value + q2.value + R2;
    *abserr = cabs((q2.value + R2) - R1) + q1.abs_err + q2.abs_err;
    return true;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Tail integral for logarithmically slow summands                    *
 *                                                                     *
 *  1/(x Log[x]^2) and 1/(x Log[x] Log[Log[x]]^2) are not power laws: *
 *  int_X^oo f = 1/Log[X] (resp. 1/Log[Log[X]]) is still 1e-3 (0.17)   *
 *  at the largest double, so neither the power-law split above nor     *
 *  exp-sinh can reach the remainder, and the MPFR exp-sinh rule refines *
 *  forever on an integrand that never ends.  After d logarithmic        *
 *  substitutions x = E_d(t) (E_0 = t, E_(i+1) = e^(E_i)),                *
 *      int_N^oo f(x) dx = int_(L_d N)^oo g_d(t) dt,                      *
 *      g_d(t) = f(E_d(t)) E_d(t) E_(d-1)(t) ... E_1(t),                  *
 *  the slowly varying factors become algebraic: g_1 = 1/t^2 for the first*
 *  example, g_2 = 1/t^2 for the second.  The power-law split is then     *
 *  applied in t (abscissa s = log t, integrand h(s) = g(e^s) e^s):        *
 *      S(T) = int_(log t_N)^(log T) h ds + T g(T)/(p(T) - 1),             *
 *  p(T) = -log2(g(2T)/g(T)), on the ladder T_j = t_N 4^j, until two       *
 *  successive splits agree.  E_d(t) is far outside the double range, so  *
 *  the summand is sampled at MPFR precision, whose exponent range         *
 *  reaches x = 2^(2^30): T up to ~7e8 for d = 1 and ~20 for d = 2.       *
 *  The last split difference is the error estimate.                     *
 * ------------------------------------------------------------------ */

typedef struct {
    NsCtx  c;        /* the summand context, spec forced to MPFR at `bits` */
    int    depth;
    long   bits;
} NsLogSub;

/* h(s) = f(E_d(e^s)) · E_d · ... · E_1 · e^s at MPFR precision. */
static bool ns_logsub_h_mpfr(void* vc, const mpfr_t s, mpfr_t out_re, mpfr_t out_im) {
    NsLogSub* L = vc;
    mpfr_prec_t p = (mpfr_prec_t)L->bits;
    mpfr_t x, J;
    mpfr_inits2(p, x, J, (mpfr_ptr)0);
    mpfr_exp(x, s, MPFR_RNDN);                    /* t = e^s            */
    mpfr_set(J, x, MPFR_RNDN);                    /* dt/ds = t          */
    for (int i = 0; i < L->depth; i++) {
        mpfr_exp(x, x, MPFR_RNDN);                /* E_(i+1) = e^(E_i)  */
        mpfr_mul(J, J, x, MPFR_RNDN);
    }
    bool ok = mpfr_number_p(x) && mpfr_number_p(J) && !mpfr_zero_p(x);
    if (ok) {
        Expr* num = ns_eval_at(&L->c, expr_new_mpfr_copy(x));
        bool inexact;
        ok = num && get_approx_mpfr(num, out_re, out_im, &inexact);
        expr_free(num);
        if (ok) {
            mpfr_mul(out_re, out_re, J, MPFR_RNDN);
            mpfr_mul(out_im, out_im, J, MPFR_RNDN);
            ok = mpfr_number_p(out_re) && mpfr_number_p(out_im);
        }
    }
    mpfr_clears(x, J, (mpfr_ptr)0);
    return ok;
}

/* Machine view of h(s) for the Gauss–Kronrod pieces. */
static bool ns_logsub_h_machine(void* vc, double s, double _Complex* out) {
    NsLogSub* L = vc;
    mpfr_prec_t p = (mpfr_prec_t)L->bits;
    mpfr_t sm, re, im;
    mpfr_inits2(p, sm, re, im, (mpfr_ptr)0);
    mpfr_set_d(sm, s, MPFR_RNDN);
    bool ok = ns_logsub_h_mpfr(vc, sm, re, im);
    if (ok) {
        *out = mpfr_get_d(re, MPFR_RNDN) + I * mpfr_get_d(im, MPFR_RNDN);
        ok = isfinite(creal(*out)) && isfinite(cimag(*out));
    }
    mpfr_clears(sm, re, im, (mpfr_ptr)0);
    return ok;
}

/* int_(log t_N)^oo h ds for one depth.  `Nbase` is the tail start N (> e^e for
 * depth 2 to keep L_d N real and positive).  `machine` integrates the pieces
 * in double (the samples are still MPFR); otherwise tanh-sinh at `bits`.
 * Writes the value and its error estimate; returns false if the ladder could
 * not be built (non-numeric samples, not a decaying power law in t). */
static bool ns_tail_logsub(const NsCtx* c0, const mpfr_t Nbase, int depth, long bits,
                           bool machine, double reltol,
                           mpfr_t out_re, mpfr_t out_im, double* err) {
    *err = INFINITY;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    NsLogSub L;
    L.c = *c0;
    L.c.spec.mode = NUMERIC_MODE_MPFR; L.c.spec.bits = bits;
    L.c.ac = NULL; L.c.ac_prec_z = NULL;
    L.c.ac_prec = NULL;
    {
        Expr* vsym = expr_new_symbol(c0->bind->name);
        const Expr* vs[1] = { vsym };
        L.c.ac_prec = autocompile_new_prec(c0->body, vs, 1, bits);
        expr_free(vsym);
    }
    L.depth = depth; L.bits = bits;

    mpfr_t tN, s0, s1, s2, h1r, h1i, h2r, h2i, acc_r, acc_i, pr, pi_, Sr, Si, Sp_r, Sp_i, tmp;
    mpfr_inits2(p, tN, s0, s1, s2, h1r, h1i, h2r, h2i, acc_r, acc_i, pr, pi_,
                Sr, Si, Sp_r, Sp_i, tmp, (mpfr_ptr)0);
    bool ok = true;
    mpfr_set(tN, Nbase, MPFR_RNDN);
    for (int i = 0; i < depth && ok; i++) {
        if (mpfr_cmp_ui(tN, 1) <= 0) ok = false;
        else mpfr_log(tN, tN, MPFR_RNDN);
    }
    if (ok && mpfr_sgn(tN) <= 0) ok = false;
    if (ok) mpfr_log(s0, tN, MPFR_RNDN);          /* s = log t */
    mpfr_set_ui(acc_r, 0, MPFR_RNDN); mpfr_set_ui(acc_i, 0, MPFR_RNDN);
    double qerr = 0.0, d_last = INFINITY, d_prev = INFINITY, hscale = 0.0;
    bool have_prev = false, converged = false;
    /* Ladder ratio q and power-probe ratio r: depth 2 has only t <= ~20 in
     * range, so it needs a finer ladder than depth 1 (t <= ~7e8). */
    const double ln4 = log(depth == 1 ? 4.0 : 1.5), ln2 = log(depth == 1 ? 2.0 : 1.1);
    for (int j = 1; ok && j <= 40; j++) {
        /* piece [s0 + (j-1) ln4, s0 + j ln4] */
        mpfr_set_d(tmp, (double)(j - 1) * ln4, MPFR_RNDN); mpfr_add(s1, s0, tmp, MPFR_RNDN);
        mpfr_set_d(tmp, (double)j * ln4, MPFR_RNDN);       mpfr_add(s2, s0, tmp, MPFR_RNDN);
        if (machine) {
            /* Absolute floor from the integrand's scale at the start: the
             * pieces shrink geometrically and a relative test alone would
             * chase roundoff on the small ones. */
            if (j == 1) {
                double _Complex h0;
                hscale = ns_logsub_h_machine(&L, mpfr_get_d(s0, MPFR_RNDN), &h0) ? cabs(h0) : 0.0;
            }
            GkResult q = gk_integrate_machine(ns_logsub_h_machine, &L,
                                              mpfr_get_d(s1, MPFR_RNDN), mpfr_get_d(s2, MPFR_RNDN),
                                              fmax(1e-300, 1e-16 * hscale), 1e-13, 50, 0, false);
            if (q.status == GK_NONNUMERIC) break;
            mpfr_set_d(pr, creal(q.value), MPFR_RNDN); mpfr_set_d(pi_, cimag(q.value), MPFR_RNDN);
            qerr += q.abs_err;
        } else {
            double ae = INFINITY;
            denint_tanhsinh_mpfr(ns_logsub_h_mpfr, &L, s1, s2, bits, 0.1 * reltol, 12,
                                 pr, pi_, &ae);
            if (!isfinite(ae) || !mpfr_number_p(pr) || !mpfr_number_p(pi_)) break;
            qerr += ae;
        }
        mpfr_add(acc_r, acc_r, pr, MPFR_RNDN); mpfr_add(acc_i, acc_i, pi_, MPFR_RNDN);
        /* remainder at T = e^(s2): h(log T)·ln r / ln(h(log T)/h(log rT)) */
        mpfr_set_d(tmp, ln2, MPFR_RNDN); mpfr_add(s1, s2, tmp, MPFR_RNDN);
        if (!ns_logsub_h_mpfr(&L, s2, h1r, h1i) || !ns_logsub_h_mpfr(&L, s1, h2r, h2i)) break;
        if (mpfr_sgn(h1r) == 0) break;
        double h1 = mpfr_get_d(h1r, MPFR_RNDN), h2 = mpfr_get_d(h2r, MPFR_RNDN);
        double h1im = mpfr_get_d(h1i, MPFR_RNDN), h2im = mpfr_get_d(h2i, MPFR_RNDN);
        /* a power law in t: same sign, real ratio in (0, 1) */
        if (!(h2 / h1 > 0.0 && h2 / h1 < 1.0)
            || fabs(h1im) > 1e-10 * fabs(h1) || fabs(h2im) > 1e-10 * fabs(h2)) {
            have_prev = false; continue;
        }
        mpfr_div(tmp, h1r, h2r, MPFR_RNDN);
        mpfr_log(tmp, tmp, MPFR_RNDN);              /* ln(h1/h2) = (p - 1) ln 2 */
        mpfr_t rem; mpfr_init2(rem, p);
        mpfr_set_d(rem, ln2, MPFR_RNDN);            /* ln r (exact enough: */
        mpfr_mul(rem, rem, h1r, MPFR_RNDN);         /*  r is the probe step) */
        mpfr_div(rem, rem, tmp, MPFR_RNDN);
        mpfr_add(Sr, acc_r, rem, MPFR_RNDN);
        mpfr_set(Si, acc_i, MPFR_RNDN);
        mpfr_clear(rem);
        if (have_prev) {
            mpfr_sub(tmp, Sr, Sp_r, MPFR_RNDN);
            d_prev = d_last;
            d_last = fabs(mpfr_get_d(tmp, MPFR_RNDN));
            double scale = fabs(mpfr_get_d(Sr, MPFR_RNDN)) + 1e-300;
            if (d_last <= reltol * scale && d_prev <= 10.0 * reltol * scale) {
                converged = true; break;
            }
            /* Not settling geometrically: the slowly varying factor survives
             * at this depth (a Log[Log[x]] seen from depth 1); give up on it. */
            if (j >= 3 && isfinite(d_prev) && d_last > 0.5 * d_prev) break;
        }
        mpfr_set(Sp_r, Sr, MPFR_RNDN); mpfr_set(Sp_i, Si, MPFR_RNDN);
        have_prev = true;
    }
    bool have = ok && isfinite(d_last);
    if (have) {
        mpfr_set(out_re, Sr, MPFR_RNDN);
        mpfr_set(out_im, Si, MPFR_RNDN);
        *err = d_last + qerr;
    }
    (void)converged;
    mpfr_clears(tN, s0, s1, s2, h1r, h1i, h2r, h2i, acc_r, acc_i, pr, pi_,
                Sr, Si, Sp_r, Sp_i, tmp, (mpfr_ptr)0);
    autocompiled_free(L.c.ac_prec);
    return have;
}

/* Best log-substituted tail over depths 1 and 2 (the smaller error wins). */
static bool ns_tail_logsub_best(const NsCtx* c, const mpfr_t Nbase, long bits, bool machine,
                                double reltol, mpfr_t out_re, mpfr_t out_im, double* err) {
    *err = INFINITY;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t r, i;
    mpfr_inits2(p, r, i, (mpfr_ptr)0);
    bool any = false;
    for (int d = 1; d <= 2; d++) {
        double e;
        if (ns_tail_logsub(c, Nbase, d, bits, machine, reltol, r, i, &e) && e < *err) {
            mpfr_set(out_re, r, MPFR_RNDN); mpfr_set(out_im, i, MPFR_RNDN);
            *err = e; any = true;
            if (e <= 10.0 * reltol * (fabs(mpfr_get_d(r, MPFR_RNDN)) + 1e-300)) break;
        }
    }
    mpfr_clears(r, i, (mpfr_ptr)0);
    return any;
}
#endif /* USE_MPFR */

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
    if (!o->nsum_terms_user) {
        long want = (settle > 0) ? settle + 4 : 0;
        /* Euler–Maclaurin's achievable accuracy is bounded by the asymptotic-
         * series minimum (~10^-2.4N in practice) at base N, so the base must
         * grow with the target precision: reaching `digits` places needs
         * N ~ 0.45*digits. The old base (settle+4) was precision-independent and
         * floored NSum at ~44 digits regardless of WorkingPrecision -- the tail
         * integral, correction count and guard bits already scale, only the base
         * did not. */
        long prec_base = (long)(0.45 * ns_wp_digits(o)) + 8;
        if (prec_base > want) want = prec_base;
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

    /* Tail integral: the split power-law rule when the summand decays
     * algebraically (exp-sinh cannot reach its remainder in double range),
     * else the exp-sinh half-line rule. */
    double _Complex Iv; double abserr;
    bool integ;
    double far_mass;
    if (ns_tail_integral_algebraic_machine(c, Nd, &Iv, &abserr, &far_mass)) {
        integ = abserr <= 1e-12 * cabs(Iv);
        far_mass = 0.0;                       /* the remainder is accounted for */
    } else {
        bool done = false;
#ifdef USE_MPFR
        /* Still live at the farthest double (a logarithmic factor): the tail
         * integral in log-substituted variables, sampled at MPFR range. */
        if (far_mass > 1e-15) {
            mpfr_t Nb, lr, li;
            mpfr_inits2(96, Nb, lr, li, (mpfr_ptr)0);
            mpfr_set_d(Nb, Nd, MPFR_RNDN);
            double le;
            if (ns_tail_logsub_best(c, Nb, 96, true, 1e-14, lr, li, &le)) {
                double _Complex v = mpfr_get_d(lr, MPFR_RNDN) + I * mpfr_get_d(li, MPFR_RNDN);
                if (le < far_mass) {
                    Iv = v; abserr = le; far_mass = 0.0;
                    integ = true; done = true;
                }
            }
            mpfr_clears(Nb, lr, li, (mpfr_ptr)0);
        }
#endif
        if (!done)
            integ = dequad_halfline_machine(ns_sample_x_machine, c, Nd, 1e-12, 14, &Iv, &abserr);
    }
    double relerr = abserr / (cabs(Iv) + 1e-300);
    if (!integ && relerr > 1e-6) { expr_free(Nval); return NULL; }   /* fall back */
    /* The exp-sinh level difference cannot see mass beyond its farthest finite
     * sample (every level loses it alike); while the summand is still live out
     * there (a log-damped 1/(x Log[x]^2)), that mass belongs in the error. */
    if (far_mass > abserr) abserr = far_mass;

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
    double tol = nc_combined_tol(o->acc_goal, o->prec_goal, cabs(H + tail),
                                 NUMERIC_MACHINE_PRECISION_DIGITS);
    bool met = false;
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
                met = true;                       /* derivative vanished: EM exact */
                break;
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
        if (mag < tol) { met = true; break; }
        dipow *= ddi * ddi;
    }
    /* Error estimate: the unconverged correction series plus the tail
     * integral's own estimate, which dominates for a slowly decaying summand. */
    {
        double est = abserr;
        if (!met && napplied > 0 && isfinite(prev_mag)) est += prev_mag;
        if (ns_goal_active(o) && est > tol) nc_warn_goal("NSum", est, tol);
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
/* f at a continuous real x, MPFR.  Routes through ns_eval_at (not
 * ns_eval_expr_at) so the Euler–Maclaurin tail integral also takes the compiled
 * MPFR fast path when the summand compiled. */
static bool ns_sample_x_mpfr(void* vc, const mpfr_t x, mpfr_t out_re, mpfr_t out_im) {
    NsCtx* c = vc;
    Expr* num = ns_eval_at(c, expr_new_mpfr_copy(x));
    if (!num) return false;
    bool inexact, ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}

/* Evaluate `e` at the index value (consumed) -> (re,im) MPFR. */
static bool ns_eval_complex_mpfr(NsCtx* c, Expr* e, Expr* value, mpfr_t re, mpfr_t im) {
    /* Compiled complex-input fast path — only for the summand itself (e ==
     * c->body); a symbolic derivative body (dcur) keeps the interpreter.  The
     * contour point `value` is a Complex[MPFR,MPFR] (or real); the compiled body
     * returns f[value] at the working precision, or NULL to interpret. */
    if (c->ac_prec_z && e == c->body && c->spec.mode == NUMERIC_MODE_MPFR) {
        Expr* r = autocompiled_eval_mpfr(c->ac_prec_z, (const Expr* const*)&value);
        if (r) {
            bool inexact, ok = get_approx_mpfr(r, re, im, &inexact);
            if (ok && (!mpfr_number_p(re) || !mpfr_number_p(im))) ok = false;
            expr_free(r);
            expr_free(value);
            return ok;
        }
    }
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
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits;
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
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits;
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
        NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits;
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
    double tail_err = 0.0;      /* certified-by-ladder tail error (log route) */
    bool logsub = false;
    if (ok) {
        /* A summand still carrying mass at x = 2^(2^29) has a logarithmically
         * slow tail (1/(x Log[x]^2)): exp-sinh would refine without end on it.
         * Integrate it in log-substituted variables instead. */
        mpfr_t xf, fr, fi;
        mpfr_inits2(p, xf, fr, fi, (mpfr_ptr)0);
        mpfr_set_ui_2exp(xf, 1, (mpfr_exp_t)1 << 29, MPFR_RNDN);
        double far = 0.0;
        if (ns_sample_x_mpfr(c, xf, fr, fi)) {
            mpfr_mul(fr, fr, xf, MPFR_RNDN); mpfr_mul(fi, fi, xf, MPFR_RNDN);
            far = ns_l1_d(fr, fi);
        }
        double scale = ns_l1_d(fN_r, fN_i) * fabs(Nd) + 1e-300;
        if (far > reltol * scale) {
            double le;
            if (ns_tail_logsub_best(c, Na, bits, false, 1e-3 * reltol, Iv_r, Iv_i, &le)) {
                logsub = true; integ = true; abserr = le; tail_err = le;
            }
        }
        mpfr_clears(xf, fr, fi, (mpfr_ptr)0);
    }
    if (ok && !logsub) {
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
        double tol = nc_combined_tol(o->acc_goal, o->prec_goal,
                                     ns_l1_d(Hr, Hi) + ns_l1_d(tail_r, tail_i),
                                     ns_wp_digits(o));
        bool met = false;
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
                    met = true;                     /* derivative vanished: EM exact */
                    break;
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
            if (mag < tol) { met = true; break; }
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
            double est = tail_err
                         + ((!met && napplied > 0 && isfinite(prev_mag)) ? prev_mag : 0.0);
            if (ns_goal_active(o) && est > tol) nc_warn_goal("NSum", est, tol);
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
     * deviation — so using guard-inflated bits here would lose accuracy.
     * The target's PRECISION, not its guard-padded width (numeric.h): the
     * guard bits are the rounding margin of the stored result, not digits
     * CVZ must converge to, and 16 bits more terms was enough to lose 11
     * digits of a peaked series at WorkingPrecision -> 30. */
    double bits = o->prec_mpfr ? (double)numeric_width_precision_bits(o->target_bits) : 53.0;
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
    bool    continuous;            /* summand numericalises off the integers */
} NsProfile;

/* Does the summand have a value at a non-integer index?  Euler–Maclaurin samples
 * the summand on a continuous real interval (its tail integral) and on a complex
 * contour (its derivative estimates), so a summand defined ONLY at integers —
 * 1/Prime[n], PartitionsP[n], anything that rejects a fractional argument —
 * cannot be summed by EM: every continuous sample fails, EM bails, and the work
 * (and the per-sample messages) is wasted.  We detect this behaviourally rather
 * than by a head blacklist: evaluate the summand (through the interpreter, not
 * the compiled fast path, so the answer reflects the true mathematical value) at
 * two generic non-integer index points and ask whether either is a finite
 * number.  Two points, not one, so a continuous summand with an accidental pole
 * at one probe is not mistaken for integer-only. */
static bool ns_summand_numeric_at(NsCtx* c, double x) {
    double _Complex v;
    Expr* num = ns_eval_expr_at(c, c->body, expr_new_real(x));
    if (!num) return false;
    bool ok = ns_to_complex(num, &v)
              && isfinite(creal(v)) && isfinite(cimag(v));
    expr_free(num);
    return ok;
}

static bool ns_summand_is_continuous(NsCtx* c) {
    double lo;
    if (!ns_to_double_real(c->imin, &lo)) return true;  /* can't probe → assume yes */
    return ns_summand_numeric_at(c, lo + 0.5) || ns_summand_numeric_at(c, lo + 0.3);
}

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
    p->continuous = true;
    /* The far-tail ladder samples deep indices.  For a nested numeric summand
     * with a dependent inner bound {k,1,n} that means a large inner computation
     * per sample for no benefit (the settle index is irrelevant there), so skip
     * the ladder for nested summands — but still allow Euler–Maclaurin, which is
     * valid for an independent inner bound {k,1,Infinity} and falls back on its
     * own otherwise. */
    if (ns_body_is_blackbox(c->body)) return;
    /* Continuity: does the summand extend off the integers?  Set before the
     * monotone early-out below, because the integer-only cases we must keep away
     * from Euler–Maclaurin (1/Prime[n], …) are exactly the monotone-decreasing
     * ones that would otherwise take it. */
    p->continuous = ns_summand_is_continuous(c);
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
    /* Euler–Maclaurin needs a summand that extends off the integers (it samples a
     * continuous tail integral and a complex derivative contour).  For an
     * integer-only summand — 1/Prime[n] and friends — that model does not exist,
     * so the only valid route is partial-sum extrapolation. */
    if (p->continuous) {
        /* far tail located a peak then decays -> Euler–Maclaurin with adaptive base */
        if (p->ladder_ok && !p->rising_far && p->nl >= 2
            && p->lmag[p->nl - 1] < p->lmag[p->nl - 2])
            return SYM_EulerMaclaurin;
        if (ns_mag_decreasing(pr, pr->n / 2)) return SYM_EulerMaclaurin;   /* head tail monotone */
    }
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

static double ns_expr_mag(const Expr* e);

/* Partial-sum extrapolation for a summand defined only at the integers
 * (DivisorSigma[0, n]/n^2, 1/Prime[n]^2, ...).  Such a summand is not a smooth
 * function of n, so its partial sums carry an irregular, arithmetic remainder
 * (d(n) jumps between 2 and ~n^eps): Wynn's epsilon algorithm has no model for
 * it, and the last-step difference it reports as its error describes the
 * transform's own convergence, not the distance to the sum.  On
 * Sum[DivisorSigma[0, n]/n^2] it reported 1.2e-7 for an answer 0.064 off.
 * The extrapolation is therefore run from two head lengths, N and 2N; both are
 * estimates of the same sum, so their difference bounds the error of the
 * better one honestly, and that (not the step) is what the accuracy-goal test
 * sees.  A summand with a genuine closed form never reaches here (see
 * ns_dirichlet); one that does gets the best value and a truthful warning. */
static Expr* ns_wynn_nonsmooth(NsCtx* c, const NsOpts* o) {
    double e1 = INFINITY, e2 = INFINITY;
#ifdef USE_MPFR
    Expr* w1 = o->prec_mpfr ? ns_wynn_mpfr(c, -1, o, &e1) : ns_wynn_machine(c, -1, o, &e1);
#else
    Expr* w1 = ns_wynn_machine(c, -1, o, &e1);
#endif
    if (!w1) return NULL;
    NsOpts o2 = *o;
    o2.nsum_terms = 2 * (o->nsum_terms > 0 ? o->nsum_terms : 1);
    g_ns_quiet++;                      /* the first estimate stands if this fails */
#ifdef USE_MPFR
    Expr* w2 = o->prec_mpfr ? ns_wynn_mpfr(c, -1, &o2, &e2) : ns_wynn_machine(c, -1, &o2, &e2);
#else
    Expr* w2 = ns_wynn_machine(c, -1, &o2, &e2);
#endif
    g_ns_quiet--;
    Expr* out;
    double err;
    if (w2) {
        Expr* d = eval_and_free(ns_mk2("Plus", expr_copy(w2),
                                       ns_mk2("Times", expr_new_integer(-1), expr_copy(w1))));
        double dm = ns_expr_mag(d);
        expr_free(d);
        err = (e2 > dm) ? e2 : dm;
        out = w2;
        expr_free(w1);
    } else {
        out = w1;
        err = e1;
    }
    if (ns_goal_active(o)) {
        double tol = nc_combined_tol(o->acc_goal, o->prec_goal, ns_expr_mag(out),
                                     ns_wp_digits(o));
        if (err > tol) nc_warn_goal("NSum", err, tol);
    }
    return out;
}

/* Sum the infinite series Σ_{k≥0} f(imin + k·di) by method dispatch.  When
 * `do_verify`, a clearly divergent series returns the ComplexInfinity symbol. */
static Expr* ns_sum_infinite(NsCtx* c, const char* var, NsOpts* o, bool do_verify) {
    const char* method = o->method;
    bool is_auto = (o->method == SYM_Automatic);
    bool continuous = true;
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
        continuous = prof.continuous;
        /* An integer-only summand cannot use Euler–Maclaurin, so its whole
         * accuracy rests on partial-sum extrapolation.  The prime-indexed tails
         * (1/Prime[n]^s) converge slowly and their partial sums are irregular, so
         * the default 15 head terms leave the extrapolator a large, bumpy tail
         * (P(2) came out 0.4485 vs 0.4522).  Summing more head terms first shrinks
         * that tail; the extra evaluations are cheap and a fast-converging integer
         * summand simply reaches full accuracy sooner.  (A longer extrapolation
         * *sequence* helps only slightly, and only the machine path is short
         * enough to matter — the MPFR sequence length already scales with the
         * bit budget and must not be clamped here.) */
        if (!prof.continuous) {
            if (!o->nsum_terms_user && o->nsum_terms < 100) o->nsum_terms = 100;
            if (!o->prec_mpfr && o->extra_terms < 24) o->extra_terms = 24;
        }
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
    if (!used && is_auto && !continuous) {
        out = ns_wynn_nonsmooth(c, o);
        used = (out != NULL);
    }
    if (!used) {
#ifdef USE_MPFR
        out = o->prec_mpfr ? ns_wynn_mpfr(c, -1, o, NULL) : ns_wynn_machine(c, -1, o, NULL);
#else
        out = ns_wynn_machine(c, -1, o, NULL);
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
    ctx.ac = NULL; ctx.ac_prec = NULL; ctx.ac_prec_z = NULL;
#ifdef USE_MPFR
    if (o->prec_mpfr) { ctx.spec.mode = NUMERIC_MODE_MPFR; ctx.spec.bits = o->bits; }
    else ctx.spec = numeric_machine_spec();
#else
    ctx.spec = numeric_machine_spec();
#endif

    /* Compile the summand once, as a function of the index.  Skipped for an
     * MPFR request, where every term needs more than a double can hold.  The
     * index is bound here through the symbol table, so an inner NSum's body
     * (which sees the outer index only as a symbol) is left to the interpreter
     * by the usual bail — no special case needed. */
    {
        Expr* vsym = expr_new_symbol(var);
        const Expr* vs[1] = { vsym };
        if (!o->prec_mpfr) {
            ctx.ac = autocompile_new(body, vs, 1);
        }
#ifdef USE_MPFR
        else if (o->bits > 0) {
            ctx.ac_prec   = autocompile_new_prec(body, vs, 1, (long)o->bits);
            /* Complex-input variant for the Euler–Maclaurin contour derivative. */
            ctx.ac_prec_z = autocompile_new_prec_z(body, vs, 1, (long)o->bits);
        }
#endif
        expr_free(vsym);
    }

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

    autocompiled_free(ctx.ac);
    autocompiled_free(ctx.ac_prec);
    autocompiled_free(ctx.ac_prec_z);
    ns_bind_restore(&bind);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Phase splitting: trigonometric and sign-alternating summands       *
 *                                                                     *
 *  A summand carrying a linear phase — Sin[a n], Cos[a n], (-1)^n,     *
 *  E^(I a n) — times a slowly decaying amplitude (Sin[n]/n,            *
 *  Cos[n]/Sqrt[n], (-1)^n Log[n]/n) defeats every method above: its    *
 *  continuous extension oscillates, so the Euler–Maclaurin tail        *
 *  integral and contour derivatives are wrong (an aliased frequency    *
 *  a > 2 Pi makes the correction series diverge outright), while a      *
 *  real summand mixes the two frequencies e^(+-i a n), which no single  *
 *  partial-sum model (Wynn, Levin) fits.  The numeric alternation       *
 *  detector also misses (-1)^n Log[n]/n, whose first term is 0.         *
 *                                                                     *
 *  The summand is therefore split symbolically,                         *
 *      f(n) = sum_j F_j(n) e^(lam_j n),   F_j free of any phase,        *
 *  by TrigToExp + Expand and stripping the factors b^(c n + d)          *
 *  (b free of n) from each term; terms are grouped by lam_j.  With the  *
 *  index n_k = n0 + k di and Lam = lam di, each group is                 *
 *      e^(lam n0) sum_k e^(Lam k) f(k),   f(x) = F(n0 + x di),           *
 *  where f is smooth (analytic for Re x >= 0 once n0 is past the head).  *
 *  e^(Lam k) depends on Lam only modulo 2 Pi i, so mu = Lam - 2 Pi i m   *
 *  is taken with Im mu in [0, 2 Pi), and:                                *
 *    - Im mu = 0 (no oscillation): the original terms join a "rest"     *
 *      summand summed by the ordinary methods above;                    *
 *    - mu = i Pi (alternating): Cohen–Villegas–Zagier on the group;     *
 *    - otherwise the Abel–Plana formula with the contour rotated onto    *
 *      the imaginary axis (g(x) = e^(mu x) f(x)):                        *
 *        sum_k g(k) = f(0)/2 + i int_0^oo [e^(i mu y) f(iy)              *
 *                     - e^(-i mu y - 2 Pi y) f(-iy)] / (1 - e^(-2 Pi y)) dy,*
 *      valid for 0 < Im mu < 2 Pi, Re mu <= 0 and f analytic and of      *
 *      sub-exponential growth in Re x >= 0.  Both kernels decay           *
 *      exponentially, at rates Im mu and 2 Pi - Im mu, so no slow         *
 *      amplitude (1/Sqrt[n], 1/Log[n]) matters.                           *
 *  The oscillatory part is computed from two base points n0 (head of     *
 *  m1 and of m2 terms summed directly); both are exact identities, so    *
 *  their difference is an honest error estimate that also exposes a     *
 *  summand that is not analytic where the formula needs it.             *
 * ------------------------------------------------------------------ */

static bool ns_free_of(const Expr* e, const char* var) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name != var;
    if (e->type != EXPR_FUNCTION) return true;
    if (!ns_free_of(e->data.function.head, var)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!ns_free_of(e->data.function.args[i], var)) return false;
    return true;
}

static bool ns_has_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* A var-dependent phase anywhere in `e`: b^ex with b free of var and ex not,
 * or a circular function of var. */
static bool ns_has_phase(const Expr* e, const char* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (ns_has_head(e, SYM_Power) && e->data.function.arg_count == 2
        && ns_free_of(e->data.function.args[0], var)
        && !ns_free_of(e->data.function.args[1], var))
        return true;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if ((h == SYM_Sin || h == SYM_Cos || h == SYM_Tan || h == SYM_Cot
             || h == SYM_Sec || h == SYM_Csc) && !ns_free_of(e, var))
            return true;
        /* J_nu(z), Y_nu(z) of a real argument oscillate like e^(+-i z): a hidden
         * linear phase when z is linear in the index (ns_bessel_to_hankel). */
        if ((h == SYM_BesselJ || h == SYM_BesselY) && e->data.function.arg_count == 2
            && !ns_free_of(e->data.function.args[1], var))
            return true;
    }
    if (ns_has_phase(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ns_has_phase(e->data.function.args[i], var)) return true;
    return false;
}

/* Split one evaluated product term into F(var)·e^(lam·var), F phase-free.
 * Returns F (new) and sets *lam (new, exact), or NULL when the term's phase
 * is not linear in var or is not a top-level factor. */
static Expr* ns_strip_phase(const Expr* term, const char* var, Expr** lam) {
    Expr* const* fs;
    Expr* one[1];
    size_t n;
    if (ns_has_head(term, SYM_Times)) {
        fs = term->data.function.args; n = term->data.function.arg_count;
    } else {
        one[0] = (Expr*)term; fs = one; n = 1;
    }
    Expr** keep = malloc(sizeof(Expr*) * n);
    if (!keep) return NULL;
    size_t nk = 0;
    Expr* L = expr_new_integer(0);
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr* f = fs[i];
        if (ns_has_head(f, SYM_Power) && f->data.function.arg_count == 2
            && ns_free_of(f->data.function.args[0], var)
            && !ns_free_of(f->data.function.args[1], var)) {
            Expr* b = f->data.function.args[0];
            Expr* ex = f->data.function.args[1];
            Expr* c = eval_and_free(ns_mk2("D", expr_copy(ex), expr_new_symbol(var)));
            if (!c || !ns_free_of(c, var)) { expr_free(c); ok = false; break; }
            L = eval_and_free(ns_mk2("Plus", L,
                                     ns_mk2("Times", c, ns_mk1("Log", expr_copy(b)))));
            /* b^(c var + d) = b^d · e^(c Log[b] var) for the principal power. */
            Expr* ex0 = eval_and_free(ns_mk2("ReplaceAll", expr_copy(ex),
                                             ns_mk2("Rule", expr_new_symbol(var),
                                                    expr_new_integer(0))));
            keep[nk++] = ns_mk2("Power", expr_copy(b), ex0);
        } else {
            keep[nk++] = expr_copy(f);
        }
    }
    if (!ok) {
        for (size_t i = 0; i < nk; i++) expr_free(keep[i]);
        free(keep); expr_free(L);
        return NULL;
    }
    Expr* s = expr_new_function(expr_new_symbol(SYM_Times), keep, nk);
    free(keep);
    s = eval_and_free(s);
    if (!s || ns_has_phase(s, var)) { expr_free(s); expr_free(L); return NULL; }
    *lam = L;
    return s;
}

typedef struct {
    Expr*          lam;      /* exact phase per unit of the index           */
    double _Complex lamd;
    Expr*          smooth;   /* F: sum of the phase-stripped terms           */
    Expr*          orig;     /* the group's original terms                   */
    Expr*          mu;       /* lam·di - 2 Pi i m, Im mu in [0, 2 Pi)         */
    bool           alternating;
    /* Conjugate pairing: a real summand's e^(+i a n) and e^(-i a n) groups
     * have F_k(x) = conj(F_j(conj x)), so group k sums to conj(group j) and
     * is not computed (it halves the work for Sin/Cos summands). */
    int            mirror_of;   /* j when this group is j's conjugate, else -1 */
    bool           has_mirror;
} NsPhaseGroup;

static Expr* ns_add_expr(Expr* acc, Expr* t) {
    if (!acc) return t;
    return ns_mk2("Plus", acc, t);
}

/* Context helper: a summation context for `body` from base `n0`. */
static void ns_ctx_setup(NsCtx* c, Expr* body, Expr* n0, Expr* di, NsBind* bind,
                         const NsOpts* o, const char* var, bool compile) {
    c->body = body; c->imin = n0; c->di = di; c->bind = bind;
    c->ac = NULL; c->ac_prec = NULL; c->ac_prec_z = NULL;
#ifdef USE_MPFR
    if (o->prec_mpfr) { c->spec.mode = NUMERIC_MODE_MPFR; c->spec.bits = o->bits; }
    else c->spec = numeric_machine_spec();
#else
    c->spec = numeric_machine_spec();
#endif
    if (compile) {
        Expr* vsym = expr_new_symbol(var);
        const Expr* vs[1] = { vsym };
        if (!o->prec_mpfr) c->ac = autocompile_new(body, vs, 1);
#ifdef USE_MPFR
        else if (o->bits > 0) c->ac_prec = autocompile_new_prec(body, vs, 1, (long)o->bits);
#endif
        expr_free(vsym);
    }
}

static void ns_ctx_teardown(NsCtx* c) {
    autocompiled_free(c->ac);
    autocompiled_free(c->ac_prec);
    c->ac = NULL; c->ac_prec = NULL;
}

/* Numericalise a var-free expression (consumed) at the context precision. */
static Expr* ns_numeric_const(Expr* e, const NsCtx* c) {
    Expr* raw = eval_and_free(e);
    Expr* num = raw ? numericalize(raw, c->spec) : NULL;
    expr_free(raw);
    return num;
}

/* L1 magnitude of a numeric Expr as a double (INFINITY if not numeric). */
static double ns_expr_mag(const Expr* e) {
    if (!e) return INFINITY;
#ifdef USE_MPFR
    {
        Expr *re, *im;
        if (e->type == EXPR_MPFR) return fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
        if (is_complex((Expr*)e, &re, &im)) {
            double r, i;
            if (ns_to_double_real(re, &r) && ns_to_double_real(im, &i)) return fabs(r) + fabs(i);
        }
    }
#endif
    double _Complex v;
    if (ns_to_complex((Expr*)e, &v)) return fabs(creal(v)) + fabs(cimag(v));
    return INFINITY;
}

/* One oscillatory group summed from base n0 by Abel–Plana:
 *   e^(lam n0) [ f(0)/2 + int_0^oo H(y) dy ],
 *   H(y) = I (e^(I mu y) F(n0 + I di y) - e^(-I mu y - 2 Pi y) F(n0 - I di y))
 *          / (1 - e^(-2 Pi y)).
 * *qerr receives the quadrature error estimate. */
/* F(n0 + sgn·I·di·var), evaluated (var stays symbolic: it is the new variable). */
static Expr* ns_shift_imag(const Expr* F, const char* var, Expr* n0, Expr* di, int sgn) {
    Expr* step = ns_mk2("Times", expr_new_integer(sgn),
                        ns_mk2("Times", expr_new_symbol(SYM_I),
                               ns_mk2("Times", expr_copy(di), expr_new_symbol(var))));
    return eval_and_free(ns_mk2("ReplaceAll", expr_copy((Expr*)F),
                                ns_mk2("Rule", expr_new_symbol(var),
                                       ns_mk2("Plus", expr_copy(n0), step))));
}

static Expr* ns_abel_plana_group(const NsPhaseGroup* g, const char* var, Expr* n0,
                                 Expr* di, NsBind* bind, const NsOpts* o, double* qerr) {
    *qerr = INFINITY;
    ns_bind_clear_temp(symtab_get_def(var));       /* symbolic: index unbound */
    Expr* shp = ns_shift_imag(g->smooth, var, n0, di, +1);
    Expr* shm = ns_shift_imag(g->smooth, var, n0, di, -1);
    if (!shp || !shm) { expr_free(shp); expr_free(shm); return NULL; }
    /* H(y), with y spelled as var. */
    Expr* imu_y = ns_mk2("Times", expr_new_symbol(SYM_I),
                         ns_mk2("Times", expr_copy(g->mu), expr_new_symbol(var)));
    Expr* e1 = ns_mk2("Power", expr_new_symbol(SYM_E), imu_y);
    Expr* e2 = ns_mk2("Power", expr_new_symbol(SYM_E),
                      ns_mk2("Plus",
                             ns_mk2("Times", expr_new_integer(-1), expr_copy(imu_y)),
                             ns_mk2("Times", ns_mk2("Times", expr_new_integer(-2),
                                                    expr_new_symbol(SYM_Pi)),
                                    expr_new_symbol(var))));
    Expr* num = ns_mk2("Plus", ns_mk2("Times", e1, shp),
                       ns_mk2("Times", expr_new_integer(-1), ns_mk2("Times", e2, shm)));
    Expr* den = ns_mk2("Plus", expr_new_integer(1),
                       ns_mk2("Times", expr_new_integer(-1),
                              ns_mk2("Power", expr_new_symbol(SYM_E),
                                     ns_mk2("Times", ns_mk2("Times", expr_new_integer(-2),
                                                            expr_new_symbol(SYM_Pi)),
                                            expr_new_symbol(var)))));
    Expr* H = eval_and_free(ns_mk2("Times", expr_new_symbol(SYM_I),
                                   ns_mk2("Times", num,
                                          ns_mk2("Power", den, expr_new_integer(-1)))));
    if (!H) return NULL;

    /* H needs only modest guard bits: its O(y) cancellation near y = 0 costs
     * a few bits, and the DE weights there are tiny (the factor-2 guard of the
     * Euler–Maclaurin path is for a different, larger cancellation). */
    NsOpts ho = *o;
    if (ho.prec_mpfr && ho.target_bits > 0) ho.bits = ho.target_bits + 48;
    NsCtx hc;
    ns_ctx_setup(&hc, H, n0, di, bind, &ho, var, true);
    Expr* integral = NULL;
#ifdef USE_MPFR
    if (o->prec_mpfr) {
        mpfr_prec_t p = (mpfr_prec_t)ho.bits;
        mpfr_t a, rr, ri;
        mpfr_inits2(p, a, rr, ri, (mpfr_ptr)0);
        mpfr_set_ui(a, 0, MPFR_RNDN);
        double digits = numeric_bits_to_digits(o->target_bits);
        double reltol = pow(10.0, -(0.6 * digits + 4.0));   /* DE: next level squares it */
        int mlev = 14 + (int)(digits / 8.0);
        if (mlev > 22) mlev = 22;
        double ae = INFINITY;
        dequad_halfline_mpfr(ns_sample_x_mpfr, &hc, a, ho.bits, reltol, mlev, rr, ri, &ae);
        if (isfinite(ae)) { integral = ns_from_complex_mpfr(rr, ri); *qerr = ae; }
        mpfr_clears(a, rr, ri, (mpfr_ptr)0);
    } else
#endif
    {
        /* Gauss–Kronrod near y = 0, where the two terms of H cancel to O(y)
         * against the O(y) denominator (its nodes stay clear of the endpoint),
         * then the exp-sinh half-line rule. */
        GkResult q1 = gk_integrate_machine(ns_sample_x_machine, &hc, 0.0, 2.0,
                                           1e-300, 1e-15, 200, 0, false);
        double _Complex v2 = 0.0; double e2r = 0.0;
        bool ok2 = (q1.status != GK_NONNUMERIC);
        if (ok2) {
            dequad_halfline_machine(ns_sample_x_machine, &hc, 2.0, 1e-14, 14, &v2, &e2r);
            ok2 = isfinite(creal(v2)) && isfinite(cimag(v2)) && isfinite(e2r);
        }
        if (ok2) {
            integral = ns_from_complex_d(q1.value + v2);
            *qerr = q1.abs_err + e2r;
        }
    }
    ns_ctx_teardown(&hc);
    expr_free(H);
    if (!integral) return NULL;

    /* e^(lam n0) (F(n0)/2 + integral) at the context precision. */
    NsCtx fc;
    ns_ctx_setup(&fc, g->smooth, n0, di, bind, o, var, false);
    Expr* f0 = ns_eval_at(&fc, expr_copy(n0));
    if (!f0) { expr_free(integral); return NULL; }
    Expr* r = ns_mk2("Times",
                     ns_mk2("Power", expr_new_symbol(SYM_E),
                            ns_mk2("Times", expr_copy(g->lam), expr_copy(n0))),
                     ns_mk2("Plus", ns_mk2("Times", ns_mk2("Power", expr_new_integer(2),
                                                           expr_new_integer(-1)), f0),
                            integral));
    Expr* out = ns_numeric_const(r, &fc);
    {
        /* Scale the quadrature error by |e^(lam n0)|. */
        double _Complex ln0;
        Expr* t = ns_mk2("Times", expr_copy(g->lam), expr_copy(n0));
        if (ns_const_machine(t, &ln0)) *qerr *= exp(creal(ln0));
    }
    return out;
}

/* One alternating group (mu = i Pi) from base n0 by Cohen–Villegas–Zagier. */
static Expr* ns_cvz_group(const NsPhaseGroup* g, const char* var, Expr* n0, Expr* di,
                          NsBind* bind, const NsOpts* o) {
    NsCtx cc;
    ns_ctx_setup(&cc, g->orig, n0, di, bind, o, var, true);
    int n = ns_cvz_terms(o);
    Expr* out;
#ifdef USE_MPFR
    out = o->prec_mpfr ? ns_cvz_mpfr(&cc, n, o->bits) : ns_cvz_machine(&cc, n);
#else
    out = ns_cvz_machine(&cc, n);
#endif
    ns_ctx_teardown(&cc);
    return out;
}

static Expr* ns_index_at(Expr* imin, Expr* di, long k) {
    return eval_and_free(ns_mk2("Plus", expr_copy(imin),
                                ns_mk2("Times", expr_new_integer(k), expr_copy(di))));
}

/* Sum of every oscillatory group from base index imin + m·di; NULL on failure. */
static Expr* ns_osc_groups_from(NsPhaseGroup* G, int ng, const char* var, Expr* imin,
                                Expr* di, long m, NsBind* bind, const NsOpts* o,
                                double* qerr) {
    Expr* n0 = ns_index_at(imin, di, m);
    Expr* acc = NULL;
    *qerr = 0.0;
    for (int j = 0; j < ng; j++) {
        if (!G[j].mu || G[j].mirror_of >= 0) continue;   /* rest / mirrored */
        ns_bind_clear_temp(symtab_get_def(var));
        double e = 0.0;
        Expr* v = G[j].alternating ? ns_cvz_group(&G[j], var, n0, di, bind, o)
                                   : ns_abel_plana_group(&G[j], var, n0, di, bind, o, &e);
        if (!v) { expr_free(acc); expr_free(n0); return NULL; }
        if (G[j].has_mirror) {
            v = eval_and_free(ns_mk2("Plus", v, ns_mk1("Conjugate", expr_copy(v))));
            e *= 2.0;
        }
        *qerr += e;
        acc = ns_add_expr(acc, v);
    }
    expr_free(n0);
    return acc ? eval_and_free(acc) : expr_new_integer(0);
}

static void ns_phase_groups_free(NsPhaseGroup* G, int ng) {
    for (int j = 0; j < ng; j++) {
        expr_free(G[j].lam); expr_free(G[j].smooth); expr_free(G[j].orig); expr_free(G[j].mu);
    }
    free(G);
}

/* Bessel functions as phase-carrying factors.
 *
 * J_nu(z) and Y_nu(z) of a real argument oscillate like cos, sin(z - phase):
 * Sum J_0(n) is Sin[n]/Sqrt[n]-like, but the phase is hidden inside the
 * function, so neither the symbolic split above nor any partial-sum model sees
 * it (the extrapolators returned 0.3348 for Sum_{n>=1} J_0(n) = 1/2).  The phase
 * is made explicit through the Hankel functions,
 *     J = (H1 + H2)/2,   Y = (H1 - H2)/(2 i),   H1,2(z) = e^(+-i z) S1,2(z),
 * where the scaled Hankel functions S1(z) = H1(z) e^(-i z) and S2(z) = H2(z) e^(i z)
 * are analytic and non-oscillatory for Re z > 0 (S ~ Sqrt[2/(Pi z)] e^(-+i(nu Pi/2
 * + Pi/4))).  With z = c n + d, c > 0, each half is a linear-phase group
 * e^(+-i c n) F(n), which the Abel-Plana rotation sums with exponentially
 * decaying kernels.  On the rotated contour z = c(n0 +- i y) + d the products
 * H1 e^(-i z) are evaluated without cancellation through the modified Bessel
 * function, H1_nu(w) = (2/Pi) i^(-nu-1) K_nu(-i w) and
 * H2_nu(w) = (2/Pi) i^(nu+1) K_nu(i w), valid for -Pi/2 < arg w <= Pi and
 * -Pi < arg w <= Pi/2 respectively -- the whole right half-plane.  (J + i Y
 * there would cancel e^(2 Im w) to leave e^(-Im w).)
 *
 * During the split the scaled functions are inert placeholder heads, so the
 * phase stripping sees only the explicit e^(+-i z) factors; they are expanded to
 * the K_nu form once the groups are formed. */
static const char* ns_hs_head(int kind) {
    static const char* h1 = NULL;
    static const char* h2 = NULL;
    if (!h1) { h1 = intern_symbol("NSum`HankelScaled1"); h2 = intern_symbol("NSum`HankelScaled2"); }
    return kind == 1 ? h1 : h2;
}

/* Copy of e with every BesselJ/BesselY whose argument depends on var replaced
 * by its Hankel decomposition.  *bad is set when such a factor is not of the
 * supported form (order depending on var, argument not real-linear in var with
 * a positive slope, or not in the right half-plane from n_check on). */
static Expr* ns_bessel_to_hankel(const Expr* e, const char* var, double n_check, bool* bad) {
    if (*bad) return NULL;
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (hd->type == EXPR_SYMBOL
        && (hd->data.symbol.name == SYM_BesselJ || hd->data.symbol.name == SYM_BesselY)
        && n == 2 && !ns_free_of(e->data.function.args[1], var)) {
        Expr* nu = e->data.function.args[0];
        Expr* z = e->data.function.args[1];
        double _Complex cd, zc;
        if (!ns_free_of(nu, var)) { *bad = true; return NULL; }
        Expr* c = eval_and_free(ns_mk2("D", expr_copy(z), expr_new_symbol(var)));
        bool lin = c && ns_free_of(c, var) && ns_const_machine(expr_copy(c), &cd)
                   && cimag(cd) == 0.0 && creal(cd) > 0.0;
        expr_free(c);
        if (lin) {
            Expr* zn = ns_mk2("ReplaceAll", expr_copy(z),
                              ns_mk2("Rule", expr_new_symbol(var), expr_new_real(n_check)));
            lin = ns_const_machine(zn, &zc) && cimag(zc) == 0.0 && creal(zc) > 0.0;
        }
        if (!lin) { *bad = true; return NULL; }
        bool is_y = (hd->data.symbol.name == SYM_BesselY);
        Expr* iz = ns_mk2("Times", expr_new_symbol(SYM_I), expr_copy(z));
        Expr* miz = ns_mk2("Times", expr_new_integer(-1), expr_copy(iz));
        Expr* t1 = ns_mk2("Times", ns_mk2("Power", expr_new_symbol(SYM_E), iz),
                          ns_mk2(ns_hs_head(1), expr_copy(nu), expr_copy(z)));
        Expr* t2 = ns_mk2("Times", ns_mk2("Power", expr_new_symbol(SYM_E), miz),
                          ns_mk2(ns_hs_head(2), expr_copy(nu), expr_copy(z)));
        if (is_y) t2 = ns_mk2("Times", expr_new_integer(-1), t2);
        /* J: 1/2 (t1 + t2);  Y: 1/(2 i) (t1 - t2) = -i/2 (t1 - t2) */
        Expr* pref = is_y ? ns_mk2("Times", ns_mk2("Rational", expr_new_integer(-1), expr_new_integer(2)),
                                   expr_new_symbol(SYM_I))
                          : ns_mk2("Rational", expr_new_integer(1), expr_new_integer(2));
        return ns_mk2("Times", pref, ns_mk2("Plus", t1, t2));
    }
    Expr* nh = ns_bessel_to_hankel(hd, var, n_check, bad);
    Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!v) { expr_free(nh); *bad = true; return NULL; }
    for (size_t i = 0; i < n; i++) {
        v[i] = ns_bessel_to_hankel(e->data.function.args[i], var, n_check, bad);
        if (*bad) {
            for (size_t j = 0; j < i; j++) expr_free(v[j]);
            free(v); expr_free(nh);
            return NULL;
        }
    }
    Expr* out = expr_new_function(nh, v, n);
    free(v);
    return out;
}

/* Copy of e with the placeholder heads replaced by their modified-Bessel form:
 *   S1_nu(z) = (2/Pi) I^(-nu-1) K_nu(-I z) E^(-I z),
 *   S2_nu(z) = (2/Pi) I^(nu+1)  K_nu(I z)  E^(I z).
 * (A fresh tree: nodes are shared by reference count, never edited in place.) */
static bool ns_has_hs(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type == EXPR_SYMBOL
        && (hd->data.symbol.name == ns_hs_head(1) || hd->data.symbol.name == ns_hs_head(2)))
        return true;
    if (ns_has_hs(hd)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ns_has_hs(e->data.function.args[i])) return true;
    return false;
}

static Expr* ns_hankel_expand(const Expr* e) {
    if (!ns_has_hs(e)) return expr_copy((Expr*)e);
    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    int kind = 0;
    if (hd->type == EXPR_SYMBOL && n == 2)
        kind = hd->data.symbol.name == ns_hs_head(1) ? 1
             : hd->data.symbol.name == ns_hs_head(2) ? 2 : 0;
    if (kind) {
        const Expr* nu = e->data.function.args[0];
        const Expr* z = e->data.function.args[1];
        int sg = (kind == 1) ? -1 : 1;
        Expr* isz = ns_mk2("Times", ns_mk2("Times", expr_new_integer(sg), expr_new_symbol(SYM_I)),
                           expr_copy((Expr*)z));
        Expr* ipow = ns_mk2("Power", expr_new_symbol(SYM_I),
                            ns_mk2("Times", expr_new_integer(sg),
                                   ns_mk2("Plus", expr_copy((Expr*)nu), expr_new_integer(1))));
        return ns_mk2("Times",
                      ns_mk2("Times", ns_mk2("Times", expr_new_integer(2),
                                             ns_mk2("Power", expr_new_symbol(SYM_Pi),
                                                    expr_new_integer(-1))), ipow),
                      ns_mk2("Times", ns_mk2(SYM_BesselK, expr_copy((Expr*)nu), expr_copy(isz)),
                             ns_mk2("Power", expr_new_symbol(SYM_E), isz)));
    }
    Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!v) return expr_copy((Expr*)e);
    for (size_t i = 0; i < n; i++) v[i] = ns_hankel_expand(e->data.function.args[i]);
    Expr* out = expr_new_function(ns_hankel_expand(hd), v, n);
    free(v);
    return out;
}

/* Attempt the phase split for an infinite Automatic sum.  Returns the sum, or
 * NULL when the summand has no oscillatory linear phase (or cannot be split),
 * in which case the caller proceeds with the ordinary method dispatch. */
static Expr* ns_phase_split(Expr* body, const char* var, Expr* imin, Expr* imax,
                            Expr* di, const NsOpts* o) {
    double dlo, ddi;
    if (!ns_to_double_real(imin, &dlo) || !ns_to_double_real(di, &ddi) || !(ddi > 0.0))
        return NULL;
    if (!ns_has_phase(body, var) || ns_body_is_blackbox(body)) return NULL;

    /* Bessel factors carry their phase implicitly: make it explicit first. */
    bool hk_bad = false;
    Expr* hk = ns_bessel_to_hankel(body, var, dlo + 16.0 * ddi, &hk_bad);
    if (hk_bad) return NULL;

    NsBind bind; ns_bind_snapshot(&bind, var);
    Expr* ex = eval_and_free(ns_mk1("Expand", ns_mk1("TrigToExp", hk)));
    NsPhaseGroup* G = NULL;
    int ng = 0, cap = 0, nosc = 0;
    bool ok = ex != NULL;
    Expr* rest = NULL;
    Expr* out = NULL;
    if (ok) {
        Expr* const* ts; Expr* one[1]; size_t nt;
        if (ns_has_head(ex, SYM_Plus)) { ts = ex->data.function.args; nt = ex->data.function.arg_count; }
        else { one[0] = ex; ts = one; nt = 1; }
        for (size_t i = 0; i < nt && ok; i++) {
            Expr* lam = NULL;
            Expr* s = ns_strip_phase(ts[i], var, &lam);
            double _Complex ld;
            if (!s || !ns_const_machine(expr_copy(lam), &ld)) {
                expr_free(s); expr_free(lam); ok = false; break;
            }
            int j = 0;
            for (; j < ng; j++)
                if (cabs(G[j].lamd - ld) <= 1e-12 * (1.0 + cabs(ld))) break;
            if (j == ng) {
                if (ng == cap) {
                    cap = cap ? 2 * cap : 4;
                    NsPhaseGroup* nG = realloc(G, sizeof(NsPhaseGroup) * (size_t)cap);
                    if (!nG) { expr_free(s); expr_free(lam); ok = false; break; }
                    G = nG;
                }
                G[ng].lam = lam; G[ng].lamd = ld; G[ng].smooth = s;
                G[ng].orig = expr_copy(ts[i]); G[ng].mu = NULL; G[ng].alternating = false;
                G[ng].mirror_of = -1; G[ng].has_mirror = false;
                ng++;
            } else {
                expr_free(lam);
                G[j].smooth = ns_mk2("Plus", G[j].smooth, s);
                G[j].orig = ns_mk2("Plus", G[j].orig, expr_copy(ts[i]));
            }
        }
    }
    for (int j = 0; ok && j < ng; j++) {
        if (ns_has_hs(G[j].smooth) || ns_has_hs(G[j].orig)) {
            Expr* s2 = ns_hankel_expand(G[j].smooth);
            Expr* o2 = ns_hankel_expand(G[j].orig);
            expr_free(G[j].smooth); expr_free(G[j].orig);
            G[j].smooth = s2; G[j].orig = o2;
        }
        G[j].smooth = eval_and_free(G[j].smooth);
        G[j].orig = eval_and_free(G[j].orig);
        if (!G[j].smooth || !G[j].orig) ok = false;
    }
    /* Classify each group by its per-step phase Lam = lam·di. */
    for (int j = 0; ok && j < ng; j++) {
        double _Complex L = G[j].lamd * ddi;
        if (creal(L) > 1e-12) { ok = false; break; }      /* growing: not ours */
        if (fabs(cimag(L)) < 1e-12) {                       /* no oscillation */
            rest = ns_add_expr(rest, expr_copy(G[j].orig));
            continue;
        }
        double twopi = 2.0 * M_PI;
        double m = floor(cimag(L) / twopi);
        double th = cimag(L) - m * twopi;
        if (th < 1e-9 || th > twopi - 1e-9) { ok = false; break; }   /* aliased 1 */
        G[j].mu = eval_and_free(ns_mk2("Plus", ns_mk2("Times", expr_copy(G[j].lam), expr_copy(di)),
                     ns_mk2("Times", expr_new_integer(-(int64_t)(2.0 * m)),
                            ns_mk2("Times", expr_new_symbol(SYM_Pi), expr_new_symbol(SYM_I)))));
        G[j].alternating = fabs(creal(L)) < 1e-12 && fabs(th - M_PI) < 1e-12;
        /* An undamped phase needs a decaying amplitude, else the series
         * diverges (Sin[n]); leave that to the ordinary dispatch. */
        if (fabs(creal(L)) < 1e-12) {
            NsCtx fc;
            ns_ctx_setup(&fc, G[j].smooth, imin, di, &bind, o, var, false);
            fc.spec = numeric_machine_spec();
            double _Complex a1, a2;
            bool d1 = ns_eval_complex_machine(&fc, G[j].smooth, expr_new_real(dlo + 1e4 * ddi), &a1);
            bool d2 = ns_eval_complex_machine(&fc, G[j].smooth, expr_new_real(dlo + 1e8 * ddi), &a2);
            if (!d1 || !d2 || !(cabs(a2) < cabs(a1))) { ok = false; break; }
        }
        nosc++;
    }
    if (ok && nosc == 0) ok = false;
    /* A plain alternating summand that the sign probe already recognises
     * (strictly alternating, magnitude decreasing from the first term) is
     * summed by the ordinary dispatch's single CVZ pass at a fraction of the
     * cost; the split is for the cases that probe misses. */
    if (ok && nosc == 1 && !rest) {
        bool only_alt = false;
        for (int j = 0; j < ng; j++) if (G[j].mu) only_alt = G[j].alternating;
        if (only_alt) {
            NsCtx pc;
            ns_ctx_setup(&pc, body, imin, di, &bind, o, var, false);
            pc.spec = numeric_machine_spec();
            NsProbe pr = ns_probe(&pc, -1);
            if (pr.ok && ns_is_alternating(&pr) && ns_mag_decreasing(&pr, 0)) ok = false;
        }
    }
    /* Pair conjugate groups (checked numerically at a generic complex point). */
    for (int j = 0; ok && j < ng; j++) {
        if (!G[j].mu || G[j].alternating || G[j].mirror_of >= 0 || G[j].has_mirror) continue;
        for (int k = j + 1; k < ng; k++) {
            if (!G[k].mu || G[k].alternating || G[k].mirror_of >= 0 || G[k].has_mirror) continue;
            if (!(cabs(G[k].lamd - conj(G[j].lamd)) <= 1e-12 * (1.0 + cabs(G[j].lamd)))) continue;
            NsCtx fc;
            ns_ctx_setup(&fc, G[j].smooth, imin, di, &bind, o, var, false);
            fc.spec = numeric_machine_spec();
            double zr = dlo + 2.371, zi = 0.613;
            double _Complex fj, fk;
            bool a1 = ns_eval_complex_machine(&fc, G[j].smooth,
                          make_complex(expr_new_real(zr), expr_new_real(zi)), &fj);
            bool a2 = ns_eval_complex_machine(&fc, G[k].smooth,
                          make_complex(expr_new_real(zr), expr_new_real(-zi)), &fk);
            if (a1 && a2 && cabs(fk - conj(fj)) <= 1e-12 * (cabs(fj) + 1e-300)) {
                G[k].mirror_of = j; G[j].has_mirror = true;
                break;
            }
        }
    }

    if (ok) {
        /* Head summed directly, then the oscillatory part from two bases. */
        const long m1 = 16, m2 = 32;
        NsCtx hc;
        ns_ctx_setup(&hc, body, imin, di, &bind, o, var, true);
        Expr* head = ns_direct(&hc, m1, o);
        ns_ctx_teardown(&hc);
        double qa = 0.0, qb = 0.0;
        Expr* oscA = head ? ns_osc_groups_from(G, ng, var, imin, di, m1, &bind, o, &qa) : NULL;
        Expr* oscB = NULL;
        if (oscA) {
            /* Base m2: the oscillatory terms m1..m2-1 directly, then the groups. */
            Expr* osc_orig = NULL;
            for (int j = 0; j < ng; j++)
                if (G[j].mu) osc_orig = ns_add_expr(osc_orig, expr_copy(G[j].orig));
            ns_bind_clear_temp(symtab_get_def(var));   /* symbolic: index unbound */
            osc_orig = eval_and_free(osc_orig);
            Expr* n1 = ns_index_at(imin, di, m1);
            NsCtx mc;
            ns_ctx_setup(&mc, osc_orig, n1, di, &bind, o, var, false);
            Expr* mid = ns_direct(&mc, m2 - m1, o);
            ns_ctx_teardown(&mc);
            expr_free(n1);
            Expr* tailB = mid ? ns_osc_groups_from(G, ng, var, imin, di, m2, &bind, o, &qb) : NULL;
            if (mid && tailB) oscB = eval_and_free(ns_mk2("Plus", mid, tailB));
            else { expr_free(mid); expr_free(tailB); }
            expr_free(osc_orig);
        }
        ns_bind_restore(&bind);

        Expr* restv = NULL;
        bool rest_ok = true;
        if (oscA && oscB && rest) {
            rest = eval_and_free(rest);
            Expr* n1 = ns_index_at(imin, di, m1);
            NsOpts ro = *o;
            restv = ns_run_single(rest, var, n1, imax, di, true, -1, &ro);
            expr_free(n1);
            rest_ok = restv != NULL;
        }
        if (head && oscA && oscB && rest_ok) {
            NsCtx sc;
            ns_ctx_setup(&sc, body, imin, di, NULL, o, var, false);
            Expr* dAB = ns_numeric_const(
                ns_mk2("Plus", expr_copy(oscA), ns_mk2("Times", expr_new_integer(-1),
                                                       expr_copy(oscB))), &sc);
            double err = ns_expr_mag(dAB);
            expr_free(dAB);
            Expr* tot = ns_mk2("Plus", head, oscA);
            if (restv) tot = ns_mk2("Plus", tot, restv);
            head = NULL; oscA = NULL; restv = NULL;
            out = ns_numeric_const(tot, &sc);
            /* A real summand gives a real sum: drop the roundoff imaginary part
             * left by the conjugate frequency pairs. */
            if (out) {
                bool real_summand = true;
                NsBind rb; ns_bind_snapshot(&rb, var);
                NsCtx pc;
                ns_ctx_setup(&pc, body, imin, di, &rb, o, var, false);
                pc.spec = numeric_machine_spec();
                for (long k = 0; k < 3 && real_summand; k++) {
                    double _Complex t;
                    Expr* num = ns_eval_expr_at(&pc, body, ns_index_value(&pc, k));
                    if (!num || !ns_to_complex(num, &t) || cimag(t) != 0.0) real_summand = false;
                    expr_free(num);
                }
                ns_bind_restore(&rb);
                if (real_summand) out = ns_numeric_const(ns_mk1("Re", out), &sc);
            }
            if (out) {
                err += qa;
                if (ns_goal_active(o)) {
                    double tol = nc_combined_tol(o->acc_goal, o->prec_goal, ns_expr_mag(out),
                                                 ns_wp_digits(o));
                    if (err > tol) nc_warn_goal("NSum", err, tol);
                }
            }
        }
        (void)qb;
        expr_free(head); expr_free(oscA); expr_free(oscB); expr_free(restv);
    } else {
        ns_bind_restore(&bind);
    }
    expr_free(rest);
    expr_free(ex);
    ns_phase_groups_free(G, ng);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Dirichlet series of arithmetic functions                           *
 * ------------------------------------------------------------------ *
 *
 * A summand g(n)/n^s with g an arithmetic function (DivisorSigma, EulerPhi,
 * MoebiusMu, LiouvilleLambda) is not a smooth function of n: d(n) is 2 at
 * every prime and large at highly composite n, phi(n)/n swings between 0 and 1.
 * Nothing that models the tail as smooth (Euler-Maclaurin's integral, the
 * partial-sum transforms) applies, and the tail of such a sum decays only like
 * a power of 1/N with an irregular remainder, so no affordable amount of direct
 * summation reaches machine accuracy either (Sum d(n)/n^2 truncated at 10^6
 * is still off by ~1e-5).  These series are exactly the Dirichlet generating
 * functions of the multiplicative functions, which Sum`Dirichlet closes in
 * zeta values -- d -> Zeta[s]^2, phi -> Zeta[s - 1]/Zeta[s], sigma_a ->
 * Zeta[s] Zeta[s - a], mu -> 1/Zeta[s], lambda -> Zeta[2 s]/Zeta[s], ... --
 * declining outside the half-plane of convergence.  The identity is exact, so
 * the sum is that closed form at the working precision, less the head terms
 * below the lower bound. */
static bool ns_has_arith_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_DivisorSigma || h == SYM_EulerPhi || h == SYM_MoebiusMu
            || h == SYM_LiouvilleLambda)
            return true;
    }
    if (ns_has_arith_head(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ns_has_arith_head(e->data.function.args[i])) return true;
    return false;
}

#define NS_DIRICHLET_MAX_HEAD 100000   /* head terms below imin summed exactly */

static Expr* ns_dirichlet(Expr* body, const char* var, Expr* imin, Expr* di,
                          const NsOpts* o) {
    if (!(di->type == EXPR_INTEGER && di->data.integer == 1)) return NULL;
    if (imin->type != EXPR_INTEGER || imin->data.integer < 1
        || imin->data.integer > NS_DIRICHLET_MAX_HEAD)
        return NULL;
    if (!ns_has_arith_head(body)) return NULL;

    NsBind bind; ns_bind_snapshot(&bind, var);        /* the index stays symbolic */
    Expr* args[4] = { expr_copy(body), expr_new_symbol(var),
                      expr_new_integer(1), expr_new_symbol(SYM_Infinity) };
    Expr* cf = eval_and_free(expr_new_function(expr_new_symbol("Sum`Dirichlet"), args, 4));
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (o->prec_mpfr) { spec.mode = NUMERIC_MODE_MPFR; spec.bits = o->bits; }
#else
    (void)o;
#endif
    Expr* out = NULL;
    Expr* num = cf ? numericalize(cf, spec) : NULL;
    expr_free(cf);
    if (num && ns_expr_mag(num) < INFINITY) {
        /* Subtract the terms below the lower bound, exactly. */
        Expr* head = expr_new_integer(0);
        for (int64_t k = 1; head && k < imin->data.integer; k++) {
            ns_bind_set(&bind, expr_new_integer(k));
            Expr* t = eval_and_free(expr_copy(body));
            if (!t) { expr_free(head); head = NULL; break; }
            head = eval_and_free(ns_mk2("Plus", head, t));
        }
        if (head) {
            Expr* diff = eval_and_free(ns_mk2("Plus", num,
                                              ns_mk2("Times", expr_new_integer(-1), head)));
            num = NULL;
            out = diff ? numericalize(diff, spec) : NULL;
            expr_free(diff);
        }
    }
    expr_free(num);
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
    o.acc_goal = NC_GOAL_MACHINE; /* AccuracyGoal -> MachinePrecision */
    o.prec_goal = NC_GOAL_AUTO;                    /* PrecisionGoal -> Automatic */
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
    /* Dirichlet series of arithmetic functions have exact zeta closed forms and
     * no smooth tail model (see ns_dirichlet). */
    if (!bad && infinite && o.method == SYM_Automatic && nspecs == 1)
        out = ns_dirichlet(eff_body, var, imin, di, &o);
    /* Linear-phase summands (Sin[n]/n, (-1)^n Log[n]/n) are split by frequency
     * before the ordinary dispatch, which cannot sum them (see ns_phase_split). */
    if (!out && !bad && infinite && o.method == SYM_Automatic && nspecs == 1)
        out = ns_phase_split(eff_body, var, imin, imax_e, di, &o);
    if (!bad && !out)
        out = ns_run_single(eff_body, var, imin, imax_e, di, infinite, count, &o);

    /* Round the guarded internal result back to the requested precision. */
#ifdef USE_MPFR
    if (out && o.prec_mpfr && o.target_bits > 0) {
        NumericSpec ts; ts.mode = NUMERIC_MODE_MPFR; ts.bits = o.target_bits;
        Expr* rounded = numericalize(out, ts);
        if (rounded) { expr_free(out); out = rounded; }
    }
#endif

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
    /* a timeout inside the muted second estimate must not leave NSum silent */
    tc_register_unwind_counter(&g_ns_quiet);
}

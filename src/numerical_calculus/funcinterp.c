/*
 * funcinterp.c
 *
 * FunctionInterpolation: sample an expression adaptively and return an
 * InterpolatingFunction that approximates it.
 *
 *   FunctionInterpolation[expr, {x, xmin, xmax}]
 *   FunctionInterpolation[expr, {x, xmin, xmax}, {y, ymin, ymax}, ...]
 *
 * HoldAll: `expr` is evaluated afresh at every sample point with the
 * variables Block-localised to numbers (iter_spec_shadow / iter_bind), so an
 * inner NIntegrate[..., {y, 0, 1}], a Module, or any other held construct
 * sees numeric x.
 *
 * Options (defaults in Options[FunctionInterpolation], so SetOptions works):
 *   InterpolationOrder     -> 3          order of the piecewise interpolant
 *   InterpolationPoints    -> 11         initial points per dimension (an
 *                                        integer, or a list per dimension)
 *   InterpolationPrecision -> Automatic  precision of the stored values (and
 *                                        of the sample coordinates); Automatic
 *                                        and MachinePrecision mean machine reals
 *   AccuracyGoal           -> Automatic  (Infinity: no absolute criterion)
 *   PrecisionGoal          -> Automatic  (half the InterpolationPrecision)
 *   MaxRecursion           -> 6          refinement passes
 *
 * Algorithm.  Each dimension starts from a uniform grid of InterpolationPoints
 * (at least order + 2) points; the samples form a tensor-product grid.  A
 * refinement pass estimates the interpolation error in every grid interval
 * without new samples: along each grid line, the (k+1)-th divided difference D
 * over the k+2 points centred on the interval, times the nodal polynomial of
 * the k+1 points nearest the interval's midpoint m, |D prod(m - t_i)|, is the
 * leading term of the order-k interpolation error at m.  An interval whose
 * estimate exceeds
 *     max(10^-AccuracyGoal, 10^-PrecisionGoal |f|)
 * on any grid line through it is bisected (the new grid line is sampled),
 * and the pass repeats until every interval meets the goals or MaxRecursion
 * passes have been made.  The estimator works in MPFR at the stored
 * precision plus guard bits when that precision exceeds machine precision,
 * so goals beyond 16 digits are measurable.
 *
 * The result is built by Interpolation[{{x1, f1}, ...}, InterpolationOrder ->
 * k] (n-D: {{{x1, y1}, f11}, ...}), so the object, its evaluation kernels,
 * derivatives and property queries are exactly those of interp.c.
 *
 * Surfaces: FunctionInterpolation returns an InterpolatingFunction object
 * (not a machine array) from held symbolic input, so it has no packed/NDArray
 * or Compile[] form -- see the EXEMPT entries of the audit tools.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "iter.h"
#include "sym_names.h"
#include "numeric.h"
#include "message.h"
#include "options.h"
#include "common.h"
#include "print.h"

#define FI_HEAD "FunctionInterpolation"
#define FI_MAXDIM 8
/* Hard cap on the number of grid samples (all dimensions together). */
#define FI_MAX_SAMPLES 4000000

/* ------------------------------------------------------------------ *
 *  State                                                              *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*   var;      /* borrowed iterator symbol                         */
    size_t  n;        /* number of grid coordinates                       */
    Expr**  cex;      /* owned: coordinates as computed (exact if bounds are) */
    Expr**  cnum;     /* owned: coordinates at the stored precision        */
    double* cd;       /* machine copies, for the estimator's geometry      */
} FiAxis;

typedef struct {
    size_t      m;            /* dimensionality */
    FiAxis      ax[FI_MAXDIM];
    Expr**      val;          /* owned: numericalized samples, row-major (last axis fastest) */
    size_t      nval;
    NumericSpec spec;         /* machine, or MPFR at the stored precision */
    int         order;        /* interpolation order k */
    double      acc_tol;      /* 10^-AccuracyGoal (0 for Infinity) */
    double      prec_tol;     /* 10^-PrecisionGoal (0 for Infinity) */
    Expr*       body;         /* borrowed held expression */
} FiState;

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool fi_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

static bool fi_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static Expr* fi_fn2(const char* head, Expr* a, Expr* b) {
    Expr* v[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), v, 2);
}

/* Machine value of a real number Expr (Integer, Real, MPFR, Rational). */
static bool fi_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real; return isfinite(*out);
#ifdef USE_MPFR
        case EXPR_MPFR:
            *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            return mpfr_number_p(e->data.mpfr) != 0;
#endif
        default: break;
    }
    Expr* v = numericalize(e, numeric_machine_spec());
    bool ok = v && v->type == EXPR_REAL && isfinite(v->data.real);
    if (ok) *out = v->data.real;
    expr_free(v);
    return ok;
}

/*
 * fi_numeric:
 *   A real numeric Expr at the stored precision: a machine Real in machine
 *   mode, an MPFR of exactly spec.bits in MPFR mode (an inexact value of a
 *   different precision is re-rounded, as SetPrecision would).  Returns NULL
 *   when `e` (borrowed) is not a real number.
 */
static Expr* fi_numeric(const Expr* e, NumericSpec spec) {
    Expr* v = numericalize(e, spec);
    if (!v) return NULL;
#ifdef USE_MPFR
    if (numeric_spec_is_mpfr(spec)) {
        mpfr_prec_t bits = (mpfr_prec_t)spec.bits;
        if (v->type == EXPR_MPFR) {
            if (!mpfr_number_p(v->data.mpfr)) { expr_free(v); return NULL; }
            if (mpfr_get_prec(v->data.mpfr) == bits) return v;
            Expr* r = expr_new_mpfr_bits(bits);
            mpfr_set(r->data.mpfr, v->data.mpfr, MPFR_RNDN);
            expr_free(v);
            return r;
        }
        if (v->type == EXPR_REAL) {
            if (!isfinite(v->data.real)) { expr_free(v); return NULL; }
            Expr* r = expr_new_mpfr_from_d(v->data.real, bits);
            expr_free(v);
            return r;
        }
        expr_free(v);
        return NULL;
    }
    if (v->type == EXPR_MPFR) {
        double d = mpfr_get_d(v->data.mpfr, MPFR_RNDN);
        expr_free(v);
        return isfinite(d) ? expr_new_real(d) : NULL;
    }
#endif
    if (v->type == EXPR_REAL) {
        if (!isfinite(v->data.real)) { expr_free(v); return NULL; }
        return v;
    }
    if (v->type == EXPR_INTEGER) {
        double d = (double)v->data.integer;
        expr_free(v);
        return expr_new_real(d);
    }
    expr_free(v);
    return NULL;
}

/* A goal option value: Automatic -> `dflt` digits, Infinity -> +inf,
 * a real number -> itself.  False on anything else. */
static bool fi_goal(const Expr* e, double dflt, double* out) {
    if (!e || fi_is_sym(e, SYM_Automatic)) { *out = dflt; return true; }
    if (fi_is_sym(e, SYM_Infinity) ||
        (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
         && e->data.function.head->data.symbol.name == SYM_DirectedInfinity)) {
        *out = INFINITY; return true;
    }
    double d;
    if (!fi_to_double(e, &d)) return false;
    *out = d;
    return true;
}

static bool fi_posint(const Expr* e, int64_t lo, int64_t hi, int64_t* out) {
    if (!e || e->type != EXPR_INTEGER) return false;
    if (e->data.integer < lo || e->data.integer > hi) return false;
    *out = e->data.integer;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Axis construction                                                  *
 * ------------------------------------------------------------------ */

static void fi_axis_free(FiAxis* a) {
    for (size_t i = 0; i < a->n; i++) { expr_free(a->cex[i]); expr_free(a->cnum[i]); }
    free(a->cex); free(a->cnum); free(a->cd);
    a->cex = a->cnum = NULL; a->cd = NULL; a->n = 0;
}

/* Set coordinate i of `a` from `cex` (consumed).  False if not real. */
static bool fi_axis_set(FiAxis* a, size_t i, Expr* cex, NumericSpec spec) {
    a->cex[i] = cex;
    a->cnum[i] = fi_numeric(cex, spec);
    if (!a->cnum[i]) return false;
    return fi_to_double(a->cnum[i], &a->cd[i]);
}

static bool fi_axis_alloc(FiAxis* a, size_t n) {
    a->n = n;
    a->cex = calloc(n, sizeof(Expr*));
    a->cnum = calloc(n, sizeof(Expr*));
    a->cd = calloc(n, sizeof(double));
    return a->cex && a->cnum && a->cd;
}

/* Uniform grid xmin + j (xmax - xmin)/(n - 1), j = 0..n-1, computed by the
 * evaluator (exact when the bounds are exact). */
static bool fi_axis_uniform(FiAxis* a, const Expr* lo, const Expr* hi, size_t n,
                            NumericSpec spec) {
    if (!fi_axis_alloc(a, n)) return false;
    for (size_t j = 0; j < n; j++) {
        Expr* c;
        if (j == 0) c = expr_copy((Expr*)lo);
        else if (j == n - 1) c = expr_copy((Expr*)hi);
        else {
            /* lo + j (hi - lo) / (n - 1) */
            Expr* diff = fi_fn2(SYM_Plus, expr_copy((Expr*)hi),
                                fi_fn2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)lo)));
            Expr* frac = fi_fn2(SYM_Rational, expr_new_integer((int64_t)j),
                                expr_new_integer((int64_t)(n - 1)));
            c = eval_and_free(fi_fn2(SYM_Plus, expr_copy((Expr*)lo),
                                     fi_fn2(SYM_Times, frac, diff)));
        }
        if (!fi_axis_set(a, j, c, spec)) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Sampling                                                           *
 * ------------------------------------------------------------------ */

/*
 * fi_sample:
 *   Evaluate the body at the grid point with per-axis indices idx[0..m).
 *   Returns the value at the stored precision, or NULL (after a message)
 *   when it is not a real number.
 */
static Expr* fi_sample(FiState* s, const size_t* idx) {
    IterShadow tok[FI_MAXDIM];
    for (size_t d = 0; d < s->m; d++) {
        tok[d] = iter_spec_shadow(s->ax[d].var);
        iter_bind(s->ax[d].var, s->ax[d].cnum[idx[d]]);
    }
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(s->body));
    for (size_t d = s->m; d-- > 0; ) iter_spec_restore(s->ax[d].var, tok[d]);
    Expr* v = raw ? fi_numeric(raw, s->spec) : NULL;
    if (!v) {
        char* pt = NULL;
        if (s->m == 1) pt = expr_to_string(s->ax[0].cnum[idx[0]]);
        char* rv = raw ? expr_to_string(raw) : NULL;
        mth_message(FI_HEAD, "nreal",
                    "The function value %s at %s is not a real number.",
                    rv ? rv : "$Failed", pt ? pt : "a sample point");
        free(pt); free(rv);
    }
    expr_free(raw);
    return v;
}

/* Fill every NULL entry of s->val by sampling. */
static bool fi_fill(FiState* s) {
    size_t idx[FI_MAXDIM];
    for (size_t lin = 0; lin < s->nval; lin++) {
        if (s->val[lin]) continue;
        size_t r = lin;
        for (size_t d = s->m; d-- > 0; ) { idx[d] = r % s->ax[d].n; r /= s->ax[d].n; }
        s->val[lin] = fi_sample(s, idx);
        if (!s->val[lin]) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Error estimate                                                     *
 * ------------------------------------------------------------------ */

/*
 * fi_line_marks:
 *   Mark (mark[j] = true) every interval j of one grid line whose estimated
 *   interpolation error exceeds the tolerance.  `t` are the line's machine
 *   coordinates (n of them), `v` the value Exprs reached by v[i*stride].
 */
static void fi_line_marks(const FiState* s, const double* t, Expr* const* v,
                          size_t stride, size_t n, bool* mark) {
    int k = s->order;
    size_t w = (size_t)k + 2;          /* window for the (k+1)-th difference */
    if (n < w) { for (size_t j = 0; j + 1 < n; j++) mark[j] = true; return; }
    if (s->acc_tol == 0.0 && s->prec_tol == 0.0) return;   /* no goals */

#ifdef USE_MPFR
    bool mp = numeric_spec_is_mpfr(s->spec);
    mpfr_prec_t wb = mp ? (mpfr_prec_t)(s->spec.bits + 32) : 53;
    mpfr_t dd[24];
    if (mp) for (size_t i = 0; i < w && i < 24; i++) mpfr_init2(dd[i], wb);
#endif
    double ddd[24];
    if (w > 24) w = 24;

    for (size_t j = 0; j + 1 < n; j++) {
        /* window s0 .. s0 + w - 1 centred on interval j */
        long s0 = (long)j - (long)(k / 2);
        if (s0 < 0) s0 = 0;
        if (s0 + (long)w > (long)n) s0 = (long)n - (long)w;
        double mid = 0.5 * (t[j] + t[j + 1]);

        /* |(k+1)-th divided difference| */
        double D;
#ifdef USE_MPFR
        if (mp) {
            for (size_t i = 0; i < w; i++)
                mpfr_set(dd[i], v[(size_t)(s0 + (long)i) * stride]->data.mpfr, MPFR_RNDN);
            for (size_t lev = 1; lev < w; lev++)
                for (size_t i = w - 1; i >= lev; i--) {
                    mpfr_sub(dd[i], dd[i], dd[i - 1], MPFR_RNDN);
                    mpfr_div_d(dd[i], dd[i], t[s0 + (long)i] - t[s0 + (long)(i - lev)], MPFR_RNDN);
                }
            D = fabs(mpfr_get_d(dd[w - 1], MPFR_RNDN));
        } else
#endif
        {
            for (size_t i = 0; i < w; i++) ddd[i] = v[(size_t)(s0 + (long)i) * stride]->data.real;
            for (size_t lev = 1; lev < w; lev++)
                for (size_t i = w - 1; i >= lev; i--)
                    ddd[i] = (ddd[i] - ddd[i - 1]) / (t[s0 + (long)i] - t[s0 + (long)(i - lev)]);
            D = fabs(ddd[w - 1]);
        }

        /* nodal polynomial over the k+1 window points nearest mid */
        long lo = s0, hi = s0 + (long)w - 1;
        if (fabs(mid - t[lo]) > fabs(t[hi] - mid)) lo++; else hi--;
        double P = 1.0;
        for (long i = lo; i <= hi; i++) P *= fabs(mid - t[i]);
        double est = D * P;

        double fa, fb;
        fi_to_double(v[j * stride], &fa);
        fi_to_double(v[(j + 1) * stride], &fb);
        double fmag = 0.5 * (fabs(fa) + fabs(fb));
        double tol = s->acc_tol;
        if (s->prec_tol * fmag > tol) tol = s->prec_tol * fmag;
        if (!(est <= tol)) mark[j] = true;
    }
#ifdef USE_MPFR
    if (mp) for (size_t i = 0; i < w; i++) mpfr_clear(dd[i]);
#endif
}

/* Per-axis interval marks over every grid line along that axis.  Returns
 * true when any interval is marked. */
static bool fi_marks(const FiState* s, bool** marks) {
    bool any = false;
    for (size_t d = 0; d < s->m; d++) {
        size_t n = s->ax[d].n;
        size_t stride = 1;
        for (size_t e = d + 1; e < s->m; e++) stride *= s->ax[e].n;
        size_t nlines = s->nval / n;
        memset(marks[d], 0, sizeof(bool) * n);
        for (size_t L = 0; L < nlines; L++) {
            /* base index: L enumerates the other axes; split L into the part
             * above axis d (outer) and below it (inner, < stride) */
            size_t inner = L % stride, outer = L / stride;
            size_t base = outer * n * stride + inner;
            fi_line_marks(s, s->ax[d].cd, s->val + base, stride, n, marks[d]);
        }
        for (size_t j = 0; j + 1 < n; j++) if (marks[d][j]) { any = true; break; }
    }
    return any;
}

/* ------------------------------------------------------------------ *
 *  Refinement                                                         *
 * ------------------------------------------------------------------ */

/*
 * fi_refine:
 *   Bisect every marked interval of every axis, carrying the existing
 *   samples over to the refined grid (new entries NULL, to be sampled).
 *   Returns false on allocation failure or a non-real midpoint.
 */
static bool fi_refine(FiState* s, bool** marks) {
    size_t* map[FI_MAXDIM] = { 0 };   /* old index -> new index, per axis */
    FiAxis nax[FI_MAXDIM];
    memset(nax, 0, sizeof nax);
    size_t newtotal = 1;
    bool ok = true;
    for (size_t d = 0; d < s->m && ok; d++) {
        FiAxis* a = &s->ax[d];
        size_t add = 0;
        for (size_t j = 0; j + 1 < a->n; j++) if (marks[d][j]) add++;
        nax[d].var = a->var;
        map[d] = malloc(sizeof(size_t) * a->n);
        if (!map[d] || !fi_axis_alloc(&nax[d], a->n + add)) { ok = false; break; }
        size_t q = 0;
        for (size_t j = 0; j < a->n && ok; j++) {
            map[d][j] = q;
            nax[d].cex[q] = expr_copy(a->cex[j]);
            nax[d].cnum[q] = expr_copy(a->cnum[j]);
            nax[d].cd[q] = a->cd[j];
            q++;
            if (j + 1 < a->n && marks[d][j]) {
                Expr* midx = eval_and_free(fi_fn2(SYM_Times,
                                 fi_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2)),
                                 fi_fn2(SYM_Plus, expr_copy(a->cex[j]), expr_copy(a->cex[j + 1]))));
                if (!fi_axis_set(&nax[d], q, midx, s->spec)) ok = false;
                q++;
            }
        }
        newtotal *= nax[d].n;
    }
    Expr** nval = ok ? calloc(newtotal, sizeof(Expr*)) : NULL;
    if (!nval) {
        for (size_t d = 0; d < s->m; d++) { fi_axis_free(&nax[d]); free(map[d]); }
        return false;
    }
    size_t idx[FI_MAXDIM];
    for (size_t lin = 0; lin < s->nval; lin++) {
        size_t r = lin;
        for (size_t d = s->m; d-- > 0; ) { idx[d] = r % s->ax[d].n; r /= s->ax[d].n; }
        size_t nl = 0;
        for (size_t d = 0; d < s->m; d++) nl = nl * nax[d].n + map[d][idx[d]];
        nval[nl] = s->val[lin];
    }
    free(s->val);
    s->val = nval;
    s->nval = newtotal;
    for (size_t d = 0; d < s->m; d++) {
        fi_axis_free(&s->ax[d]);
        s->ax[d] = nax[d];
        free(map[d]);
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Result                                                             *
 * ------------------------------------------------------------------ */

static Expr* fi_build(const FiState* s) {
    Expr** rows = malloc(sizeof(Expr*) * s->nval);
    if (!rows) return NULL;
    size_t idx[FI_MAXDIM];
    for (size_t lin = 0; lin < s->nval; lin++) {
        size_t r = lin;
        for (size_t d = s->m; d-- > 0; ) { idx[d] = r % s->ax[d].n; r /= s->ax[d].n; }
        Expr* coord;
        if (s->m == 1) coord = expr_copy(s->ax[0].cnum[idx[0]]);
        else {
            Expr* cs[FI_MAXDIM];
            for (size_t d = 0; d < s->m; d++) cs[d] = expr_copy(s->ax[d].cnum[idx[d]]);
            coord = expr_new_function(expr_new_symbol(SYM_List), cs, s->m);
        }
        rows[lin] = fi_fn2(SYM_List, coord, expr_copy(s->val[lin]));
    }
    Expr* data = expr_new_function(expr_new_symbol(SYM_List), rows, s->nval);
    free(rows);
    Expr* call = fi_fn2("Interpolation", data,
                        fi_fn2(SYM_Rule, expr_new_symbol(SYM_InterpolationOrder),
                               expr_new_integer(s->order)));
    return eval_and_free(call);
}

static void fi_state_free(FiState* s) {
    for (size_t d = 0; d < s->m; d++) fi_axis_free(&s->ax[d]);
    if (s->val) for (size_t i = 0; i < s->nval; i++) expr_free(s->val[i]);
    free(s->val);
    s->val = NULL;
}

/* ------------------------------------------------------------------ *
 *  Builtin                                                            *
 * ------------------------------------------------------------------ */

static Expr* fi_eval_opt(const Expr* e) {
    return e ? eval_and_free(expr_copy((Expr*)e)) : NULL;
}

static Expr* builtin_functioninterpolation(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return builtin_arg_error(FI_HEAD, argc, 2, 2 + FI_MAXDIM);

    const Expr *o_order = NULL, *o_points = NULL, *o_iprec = NULL,
               *o_acc = NULL, *o_prec = NULL, *o_maxrec = NULL;
    OptEntry ents[] = {
        { "InterpolationOrder",     &o_order,  NULL },
        { "InterpolationPoints",    &o_points, NULL },
        { "InterpolationPrecision", &o_iprec,  NULL },
        { "AccuracyGoal",           &o_acc,    NULL },
        { "PrecisionGoal",          &o_prec,   NULL },
        { "MaxRecursion",           &o_maxrec, NULL },
    };
    size_t npos = argc;
    if (!options_extract(res, FI_HEAD, ents, sizeof ents / sizeof ents[0], &npos)) {
        mth_message(FI_HEAD, "optx", "Unknown option in %s.", FI_HEAD);
        return NULL;
    }
    if (npos < 2 || npos - 1 > FI_MAXDIM) return NULL;

    /* Evaluated option values (the call is HoldAll). */
    Expr* v_order  = fi_eval_opt(o_order);
    Expr* v_points = fi_eval_opt(o_points);
    Expr* v_iprec  = fi_eval_opt(o_iprec);
    Expr* v_acc    = fi_eval_opt(o_acc);
    Expr* v_prec   = fi_eval_opt(o_prec);
    Expr* v_maxrec = fi_eval_opt(o_maxrec);
    Expr* lo[FI_MAXDIM] = { 0 };
    Expr* hi[FI_MAXDIM] = { 0 };
    Expr* result = NULL;
    FiState s;
    memset(&s, 0, sizeof s);
    bool* marks[FI_MAXDIM] = { 0 };

    /* InterpolationOrder */
    int64_t k = 3;
    if (v_order && !fi_is_sym(v_order, SYM_Automatic) && !fi_posint(v_order, 1, 20, &k)) {
        mth_message(FI_HEAD, "iord",
                    "Value of option InterpolationOrder should be an integer between 1 and 20.");
        goto done;
    }
    s.order = (int)k;

    /* InterpolationPrecision */
    s.spec = numeric_machine_spec();
    double ipd = NUMERIC_MACHINE_PRECISION_DIGITS;
    if (v_iprec && !fi_is_sym(v_iprec, SYM_Automatic) && !fi_is_sym(v_iprec, SYM_MachinePrecision)) {
        double p;
        if (!fi_to_double(v_iprec, &p) || !(p >= 1.0) || p > 1e5) {
            mth_message(FI_HEAD, "ipre",
                        "Value of option InterpolationPrecision should be a positive number.");
            goto done;
        }
#ifdef USE_MPFR
        s.spec.mode = NUMERIC_MODE_MPFR;
        s.spec.bits = numeric_digits_to_bits(p);
        ipd = p;
#endif
    }

    /* Goals: AccuracyGoal Automatic = Infinity, PrecisionGoal Automatic = ip/2 */
    double ag, pg;
    if (!fi_goal(v_acc, INFINITY, &ag) || !fi_goal(v_prec, ipd / 2.0, &pg)) {
        mth_message(FI_HEAD, "goal",
                    "Values of options AccuracyGoal and PrecisionGoal should be real numbers, "
                    "Infinity or Automatic.");
        goto done;
    }
    s.acc_tol  = isinf(ag) ? 0.0 : pow(10.0, -ag);
    s.prec_tol = isinf(pg) ? 0.0 : pow(10.0, -pg);

    /* MaxRecursion */
    int64_t maxrec = 6;
    if (v_maxrec && !fi_is_sym(v_maxrec, SYM_Automatic) && !fi_posint(v_maxrec, 0, 60, &maxrec)) {
        mth_message(FI_HEAD, "mxrec",
                    "Value of option MaxRecursion should be a non-negative integer.");
        goto done;
    }

    /* Ranges */
    s.m = npos - 1;
    s.body = res->data.function.args[0];
    for (size_t d = 0; d < s.m; d++) {
        Expr* sp = res->data.function.args[d + 1];
        if (!fi_is_list(sp) || sp->data.function.arg_count != 3
            || sp->data.function.args[0]->type != EXPR_SYMBOL)
            goto done;
        s.ax[d].var = sp->data.function.args[0];
        lo[d] = eval_and_free(expr_copy(sp->data.function.args[1]));
        hi[d] = eval_and_free(expr_copy(sp->data.function.args[2]));
        double a, b;
        if (!fi_to_double(lo[d], &a) || !fi_to_double(hi[d], &b) || !(a < b)) {
            mth_message(FI_HEAD, "range",
                        "The range specification at position %zu should be {var, min, max} "
                        "with real min < max.", d + 2);
            goto done;
        }
    }

    /* InterpolationPoints: an integer, or one per dimension */
    size_t s_m = s.m;
    for (size_t d = 0; d < s_m; d++) {
        int64_t np = 11;
        const Expr* pe = v_points;
        if (pe && fi_is_list(pe)) {
            pe = pe->data.function.arg_count == s_m ? pe->data.function.args[d] : NULL;
            if (!pe) {
                mth_message(FI_HEAD, "ipts",
                            "Value of option InterpolationPoints should be a positive integer "
                            "or a list of one per dimension.");
                goto done;
            }
        }
        if (pe && !fi_is_sym(pe, SYM_Automatic) && !fi_posint(pe, 2, 10000000, &np)) {
            mth_message(FI_HEAD, "ipts",
                        "Value of option InterpolationPoints should be a positive integer "
                        "or a list of one per dimension.");
            goto done;
        }
        if (np < k + 2) np = k + 2;
        if (!fi_axis_uniform(&s.ax[d], lo[d], hi[d], (size_t)np, s.spec)) goto done;
        marks[d] = NULL;
    }

    s.nval = 1;
    for (size_t d = 0; d < s.m; d++) s.nval *= s.ax[d].n;
    if (s.nval > FI_MAX_SAMPLES) goto done;
    s.val = calloc(s.nval, sizeof(Expr*));
    if (!s.val) goto done;
    if (!fi_fill(&s)) goto done;

    for (int64_t pass = 0; pass < maxrec; pass++) {
        for (size_t d = 0; d < s.m; d++) {
            free(marks[d]);
            marks[d] = calloc(s.ax[d].n, sizeof(bool));
            if (!marks[d]) goto done;
        }
        if (!fi_marks(&s, marks)) break;
        /* would the refined grid stay under the sample cap? */
        size_t tot = 1;
        for (size_t d = 0; d < s.m; d++) {
            size_t add = 0;
            for (size_t j = 0; j + 1 < s.ax[d].n; j++) if (marks[d][j]) add++;
            tot *= s.ax[d].n + add;
        }
        if (tot > FI_MAX_SAMPLES) break;
        if (!fi_refine(&s, marks)) goto done;
        if (!fi_fill(&s)) goto done;
    }

    result = fi_build(&s);

done:
    for (size_t d = 0; d < FI_MAXDIM; d++) { free(marks[d]); expr_free(lo[d]); expr_free(hi[d]); }
    fi_state_free(&s);
    expr_free(v_order); expr_free(v_points); expr_free(v_iprec);
    expr_free(v_acc); expr_free(v_prec); expr_free(v_maxrec);
    return result;
}

void funcinterp_init(void) {
    symtab_add_builtin(FI_HEAD, builtin_functioninterpolation);
    symtab_get_def(FI_HEAD)->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    Expr* rules[6] = {
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_AccuracyGoal), expr_new_symbol(SYM_Automatic)),
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_InterpolationOrder), expr_new_integer(3)),
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_InterpolationPoints), expr_new_integer(11)),
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_InterpolationPrecision), expr_new_symbol(SYM_Automatic)),
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_MaxRecursion), expr_new_integer(6)),
        fi_fn2(SYM_Rule, expr_new_symbol(SYM_PrecisionGoal), expr_new_symbol(SYM_Automatic)),
    };
    symtab_set_options(FI_HEAD, expr_new_function(expr_new_symbol(SYM_List), rules, 6));
    symtab_set_docstring(FI_HEAD,
        "FunctionInterpolation[expr, {x, xmin, xmax}] evaluates expr at numeric x, refining the "
        "sample grid adaptively, and returns an InterpolatingFunction approximating expr on "
        "[xmin, xmax].\n"
        "FunctionInterpolation[expr, {x, xmin, xmax}, {y, ymin, ymax}, ...] interpolates on a "
        "tensor-product grid.\n"
        "Options: InterpolationOrder (3), InterpolationPoints (11, initial points per dimension), "
        "InterpolationPrecision (Automatic = machine; the precision of the stored values), "
        "AccuracyGoal (Automatic = Infinity), PrecisionGoal (Automatic = half the "
        "InterpolationPrecision), MaxRecursion (6 refinement passes). An interval is bisected "
        "while its estimated error exceeds max(10^-AccuracyGoal, 10^-PrecisionGoal |f|).");
}

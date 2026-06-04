/*
 * interp.c
 *
 * InterpolatingFunction --- piecewise-polynomial interpolation of tabulated
 * data on a regular (tensor-product) grid, modelled on Mathematica's
 * InterpolatingFunction object.
 *
 *   InterpolatingFunction[domain, table]
 *   InterpolatingFunction[domain, table, ders]   (derivative-annotated form)
 *
 *     domain = {{x1min, x1max}, ..., {xmmin, xmmax}}   -- one interval per
 *              dimension; the number of intervals m is the dimensionality.
 *     table  = {{coord1, val1}, {coord2, val2}, ...}   -- data points.
 *              For m == 1 each coord is a scalar x; for m >= 2 each coord is
 *              a list {x1, ..., xm}.  The points must fill the Cartesian
 *              product of the per-dimension grids.
 *     ders   = {d1, ..., dm}   -- (optional) non-negative derivative orders.
 *              Absent means the value itself (all zero).
 *
 * The object is a normal form: it carries the domain, the data table and an
 * optional derivative annotation, and persists unevaluated until applied.
 * Application,
 *
 *   InterpolatingFunction[...][x1, ..., xm]
 *
 * is dispatched by the evaluator (eval.c) to interp_apply() below.
 *
 * Method
 * ------
 * Default order-3 (piecewise cubic) interpolation in every dimension,
 * matching Mathematica.  In each dimension a sliding window of (order + 1)
 * consecutive grid points, centred on the bracketing interval, is selected;
 * the value is the tensor product of the per-dimension Newton
 * divided-difference polynomials.  The order drops to min(3, n_k - 1) when a
 * dimension has fewer than four grid points.
 *
 * Derivatives
 * -----------
 * Derivative[d1, ..., dm][InterpolatingFunction[...]] is reduced (in eval.c,
 * via interp_make_derivative) to a fresh InterpolatingFunction carrying the
 * accumulated orders.  When applied, the mixed partial is obtained by
 * evaluating the d_k-th derivative of the Newton polynomial in dimension k.
 * Orders beyond the local polynomial degree yield 0.
 *
 * Outside the data range a warning (InterpolatingFunction::dmval) is issued
 * and the nearest window is used to extrapolate.  An exact argument tuple
 * that coincides with a grid node returns that node's value exactly (value
 * queries only); all other queries return a machine real.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "interp.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sym_names.h"

/* --- small numeric / structural helpers ------------------------------- */

/* True when `e` is List[...] (head is the List symbol). */
static bool interp_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

/*
 * node_to_double:
 *   Coerce a real-valued atom to a double.  Handles Integer, Real, BigInt,
 *   MPFR (when enabled) and exact Rational[n, d] (including big numerators /
 *   denominators).  Returns false for anything non-real, leaving *out
 *   untouched.
 */
static bool node_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }

    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        *out = (double)n / (double)d;
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den;
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        double dd = mpz_get_d(den);
        if (dd == 0.0) { mpz_clears(num, den, NULL); return false; }
        *out = mpz_get_d(num) / dd;
        mpz_clears(num, den, NULL);
        return true;
    }
    return false;
}

/* True for an exact (non-floating) real atom: Integer, BigInt or Rational. */
static bool is_exact_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return is_rational_like(e);
}

/* Extract a non-negative integer from `e` into *out. */
static bool node_to_order(const Expr* e, int* out) {
    if (e && e->type == EXPR_INTEGER && e->data.integer >= 0
        && e->data.integer < 1000000) {
        *out = (int)e->data.integer;
        return true;
    }
    return false;
}

/* --- interpolation kernel --------------------------------------------- */

/*
 * Locate the bracketing interval for `x`: the largest index i with
 * xs[i] <= x, clamped to [0, n-2].  xs is strictly increasing, n >= 2.
 */
static size_t bracket_interval(const double* xs, size_t n, double x) {
    if (x <= xs[0]) return 0;
    if (x >= xs[n - 1]) return n - 2;
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (xs[mid] <= x) lo = mid; else hi = mid;
    }
    return lo;
}

/*
 * newton_deriv_eval:
 *   Evaluate the d-th derivative at p of the Newton divided-difference
 *   interpolating polynomial through the w points (xs[k], ys[k]),
 *   k = 0..w-1 (window-local, zero-based).  d == 0 gives the value itself;
 *   d >= w (beyond the polynomial degree) gives 0.
 */
static double newton_deriv_eval(const double* xs, const double* ys,
                                size_t w, double p, size_t d) {
    double* c = malloc(sizeof(double) * w);
    for (size_t k = 0; k < w; k++) c[k] = ys[k];
    for (size_t j = 1; j < w; j++) {
        for (size_t k = w - 1; k >= j; k--) {
            c[k] = (c[k] - c[k - 1]) / (xs[k] - xs[k - j]);
            if (k == j) break; /* size_t underflow guard */
        }
    }

    /* Horner on the Newton form, carrying derivatives 0..d simultaneously.
     * After the sweep val[t] == P^(t)(p)/t!, so P^(d)(p) = d! * val[d]. */
    double* val = calloc(d + 1, sizeof(double));
    val[0] = c[w - 1];
    for (size_t ii = w - 1; ii-- > 0; ) {
        double dx = p - xs[ii];
        for (size_t t = d; t >= 1; t--) {
            val[t] = val[t] * dx + val[t - 1];
            if (t == 1) break;
        }
        val[0] = val[0] * dx + c[ii];
    }
    double fact = 1.0;
    for (size_t t = 2; t <= d; t++) fact *= (double)t;
    double result = fact * val[d];

    free(c);
    free(val);
    return result;
}

/* --- parsed object ----------------------------------------------------- */

typedef struct {
    size_t  m;          /* dimensionality */
    size_t* nk;         /* grid size per dimension, length m */
    double** grid;      /* grid[k] = sorted distinct abscissae, length nk[k] */
    size_t* stride;     /* row-major strides, length m */
    size_t  total;      /* prod(nk) */
    double* V;          /* flat value tensor, length total */
    double* dmin;       /* domain lower bounds, length m */
    double* dmax;       /* domain upper bounds, length m */
    bool*   has_range;  /* whether dmin/dmax[k] are valid, length m */
} IFun;

static void ifun_free(IFun* f) {
    if (f->grid) {
        for (size_t k = 0; k < f->m; k++) free(f->grid[k]);
        free(f->grid);
    }
    free(f->nk);
    free(f->stride);
    free(f->V);
    free(f->dmin);
    free(f->dmax);
    free(f->has_range);
    f->grid = NULL; f->nk = f->stride = NULL; f->V = f->dmin = f->dmax = NULL;
    f->has_range = NULL; f->m = f->total = 0;
}

/* Borrowed view of one data point's coordinates / value. */
static Expr* entry_value(Expr* entry) { return entry->data.function.args[1]; }

/* Coordinate k (0-based) of a data entry, given dimensionality m. */
static Expr* entry_coord(Expr* entry, size_t m, size_t k) {
    Expr* point = entry->data.function.args[0];
    if (m == 1) return point;
    if (!interp_is_list(point) || point->data.function.arg_count != m) return NULL;
    return point->data.function.args[k];
}

/* Insert `v` into a sorted-unique double array `g` of current length *len.
 * Keeps ascending order; ignores duplicates (within an exact compare). */
static void grid_insert(double* g, size_t* len, double v) {
    size_t i = 0;
    while (i < *len && g[i] < v) i++;
    if (i < *len && g[i] == v) return; /* duplicate */
    for (size_t j = *len; j > i; j--) g[j] = g[j - 1];
    g[i] = v;
    (*len)++;
}

/* Index of `v` in sorted array g[0..len-1], or len if absent (exact match). */
static size_t grid_index(const double* g, size_t len, double v) {
    size_t lo = 0, hi = len;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (g[mid] < v) lo = mid + 1; else hi = mid;
    }
    return (lo < len && g[lo] == v) ? lo : len;
}

/*
 * parse_domain_dim:
 *   Determine the dimensionality m from `domain` = {{...}, ...}.  Returns 0
 *   if domain is not a non-empty List.
 */
static size_t parse_domain_dim(Expr* domain) {
    if (!interp_is_list(domain)) return 0;
    return domain->data.function.arg_count;
}

/*
 * build_ifun:
 *   Parse `domain` and `table` into a tensor-product IFun.  Returns true on
 *   success (out fully populated), false on any malformation (no leaks).
 */
static bool build_ifun(Expr* domain, Expr* table, IFun* out) {
    out->m = 0; out->nk = out->stride = NULL; out->grid = NULL;
    out->V = out->dmin = out->dmax = NULL; out->has_range = NULL; out->total = 0;

    size_t m = parse_domain_dim(domain);
    if (m == 0) return false;
    if (!interp_is_list(table)) return false;
    size_t npts = table->data.function.arg_count;
    if (npts < 2) return false;

    /* Every entry must be a {coord, value} pair with real value and, for
     * m >= 2, an m-vector coordinate of reals. */
    for (size_t i = 0; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        if (!interp_is_list(e) || e->data.function.arg_count != 2) return false;
        double tmp;
        if (!node_to_double(entry_value(e), &tmp)) return false;
        for (size_t k = 0; k < m; k++) {
            Expr* ck = entry_coord(e, m, k);
            if (!ck || !node_to_double(ck, &tmp)) return false;
        }
    }

    size_t*  nk     = calloc(m, sizeof(size_t));
    double** grid   = calloc(m, sizeof(double*));
    size_t*  stride = calloc(m, sizeof(size_t));
    double*  dmin   = calloc(m, sizeof(double));
    double*  dmax   = calloc(m, sizeof(double));
    bool*    hasr   = calloc(m, sizeof(bool));
    if (!nk || !grid || !stride || !dmin || !dmax || !hasr) goto oom;

    /* Per-dimension grids: distinct abscissae over all data points. */
    for (size_t k = 0; k < m; k++) {
        grid[k] = malloc(sizeof(double) * npts);
        if (!grid[k]) goto oom;
        size_t len = 0;
        for (size_t i = 0; i < npts; i++) {
            double v;
            node_to_double(entry_coord(table->data.function.args[i], m, k), &v);
            grid_insert(grid[k], &len, v);
        }
        nk[k] = len;
        if (len < 2) goto fail;       /* need >= 2 nodes to interpolate */
    }

    /* Row-major strides and total tensor size. */
    size_t total = 1;
    for (size_t k = 0; k < m; k++) total *= nk[k];
    stride[m - 1] = 1;
    for (size_t k = m - 1; k-- > 0; ) stride[k] = stride[k + 1] * nk[k + 1];

    double* V = malloc(sizeof(double) * total);
    bool* present = calloc(total, sizeof(bool));
    if (!V || !present) { free(V); free(present); goto oom; }

    for (size_t i = 0; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        size_t flat = 0;
        for (size_t k = 0; k < m; k++) {
            double v;
            node_to_double(entry_coord(e, m, k), &v);
            flat += grid_index(grid[k], nk[k], v) * stride[k];
        }
        double val;
        node_to_double(entry_value(e), &val);
        V[flat] = val;
        present[flat] = true;
    }
    /* The data must fill the full Cartesian product. */
    for (size_t i = 0; i < total; i++) {
        if (!present[i]) { free(V); free(present); goto fail; }
    }
    free(present);

    /* Domain bounds per dimension (used only for the range warning). */
    for (size_t k = 0; k < m; k++) {
        Expr* iv = domain->data.function.args[k];
        if (interp_is_list(iv) && iv->data.function.arg_count == 2
            && node_to_double(iv->data.function.args[0], &dmin[k])
            && node_to_double(iv->data.function.args[1], &dmax[k])) {
            hasr[k] = true;
        }
    }

    out->m = m; out->nk = nk; out->grid = grid; out->stride = stride;
    out->total = total; out->V = V;
    out->dmin = dmin; out->dmax = dmax; out->has_range = hasr;
    return true;

oom:
fail:
    if (grid) { for (size_t k = 0; k < m; k++) free(grid[k]); free(grid); }
    free(nk); free(stride); free(dmin); free(dmax); free(hasr);
    return false;
}

/* --- recursive tensor-product evaluation ------------------------------- */

typedef struct {
    const IFun*   f;
    const double* p;    /* query coordinates, length m */
    const int*    ders; /* derivative orders, length m */
    const size_t* s;    /* window start per dimension */
    const size_t* w;    /* window size per dimension */
} EvalCtx;

static double eval_dim(const EvalCtx* ctx, size_t k, size_t base) {
    size_t wk = ctx->w[k], sk = ctx->s[k];
    double* tmp = malloc(sizeof(double) * wk);
    for (size_t t = 0; t < wk; t++) {
        size_t off = base + (sk + t) * ctx->f->stride[k];
        tmp[t] = (k + 1 == ctx->f->m) ? ctx->f->V[off]
                                      : eval_dim(ctx, k + 1, off);
    }
    double res = newton_deriv_eval(&ctx->f->grid[k][sk], tmp, wk,
                                   ctx->p[k], (size_t)ctx->ders[k]);
    free(tmp);
    return res;
}

/* --- public entry points ---------------------------------------------- */

Expr* interp_apply(Expr* ifun, Expr** call_args, size_t argc) {
    /* The object must be InterpolatingFunction[domain, table] or
     * InterpolatingFunction[domain, table, ders]. */
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc != 2 && obj_argc != 3) return NULL;
    Expr* domain = ifun->data.function.args[0];
    Expr* table  = ifun->data.function.args[1];

    size_t m = parse_domain_dim(domain);
    if (m == 0 || argc != m) return NULL;   /* arity must match dimensionality */

    /* All arguments must be real numbers; otherwise stay symbolic. */
    double* p = malloc(sizeof(double) * m);
    for (size_t k = 0; k < m; k++) {
        if (!node_to_double(call_args[k], &p[k])) { free(p); return NULL; }
    }

    /* Derivative orders (default all zero). */
    int* ders = calloc(m, sizeof(int));
    bool any_der = false;
    if (obj_argc == 3) {
        Expr* ds = ifun->data.function.args[2];
        if (!interp_is_list(ds) || ds->data.function.arg_count != m) {
            free(p); free(ders); return NULL;
        }
        for (size_t k = 0; k < m; k++) {
            if (!node_to_order(ds->data.function.args[k], &ders[k])) {
                free(p); free(ders); return NULL;
            }
            if (ders[k] > 0) any_der = true;
        }
    }

    IFun f;
    if (!build_ifun(domain, table, &f)) { free(p); free(ders); return NULL; }

    /* Exact node-coincident value query -> exact stored value. */
    if (!any_der) {
        bool all_exact = true;
        for (size_t k = 0; k < m; k++)
            if (!is_exact_real(call_args[k])) { all_exact = false; break; }
        if (all_exact) {
            size_t npts = table->data.function.arg_count;
            for (size_t i = 0; i < npts; i++) {
                Expr* e = table->data.function.args[i];
                bool hit = true;
                for (size_t k = 0; k < m; k++) {
                    if (!expr_eq(call_args[k], entry_coord(e, m, k))) { hit = false; break; }
                }
                if (hit) {
                    Expr* y = expr_copy(entry_value(e));
                    ifun_free(&f); free(p); free(ders);
                    return y;
                }
            }
        }
    }

    /* Range check: one dmval warning if any coordinate is out of domain. */
    for (size_t k = 0; k < m; k++) {
        if (f.has_range[k] && (p[k] < f.dmin[k] || p[k] > f.dmax[k])) {
            fprintf(stderr,
                    "InterpolatingFunction::dmval: Input value %g lies outside "
                    "the range of data in the interpolating function. "
                    "Extrapolation will be used.\n", p[k]);
            break;
        }
    }

    /* Per-dimension windows (order-3, centred, clamped). */
    size_t* s = malloc(sizeof(size_t) * m);
    size_t* w = malloc(sizeof(size_t) * m);
    for (size_t k = 0; k < m; k++) {
        size_t nkk = f.nk[k];
        size_t order = nkk - 1; if (order > 3) order = 3;
        size_t wk = order + 1;
        size_t i = bracket_interval(f.grid[k], nkk, p[k]);
        size_t shift = (wk >= 2) ? (wk / 2 - 1) : 0;
        size_t sk = (i < shift) ? 0 : i - shift;
        if (sk > nkk - wk) sk = nkk - wk;
        s[k] = sk; w[k] = wk;
    }

    EvalCtx ctx = { &f, p, ders, s, w };
    double value = eval_dim(&ctx, 0, 0);

    free(s); free(w);
    ifun_free(&f);
    free(p); free(ders);
    return expr_new_real(value);
}

/*
 * interp_make_derivative:
 *   Reduce Derivative[d1, ..., dm][InterpolatingFunction[...]] to a fresh
 *   InterpolatingFunction carrying the accumulated derivative orders.
 *   `deriv_head` is Derivative[d1, ..., dm]; `ifun` is the object.  Returns
 *   NULL (leaving the expression alone) if the orders do not match the
 *   object's dimensionality or are not non-negative integers.
 */
Expr* interp_make_derivative(Expr* deriv_head, Expr* ifun) {
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc != 2 && obj_argc != 3) return NULL;
    Expr* domain = ifun->data.function.args[0];
    Expr* table  = ifun->data.function.args[1];

    size_t m = parse_domain_dim(domain);
    if (m == 0) return NULL;
    if (deriv_head->data.function.arg_count != m) return NULL;

    /* New orders = existing orders (if any) + the requested increments. */
    int* orders = calloc(m, sizeof(int));
    if (!orders) return NULL;
    if (obj_argc == 3) {
        Expr* ds = ifun->data.function.args[2];
        if (!interp_is_list(ds) || ds->data.function.arg_count != m) { free(orders); return NULL; }
        for (size_t k = 0; k < m; k++)
            if (!node_to_order(ds->data.function.args[k], &orders[k])) { free(orders); return NULL; }
    }
    for (size_t k = 0; k < m; k++) {
        int inc;
        if (!node_to_order(deriv_head->data.function.args[k], &inc)) { free(orders); return NULL; }
        orders[k] += inc;
    }

    /* Build InterpolatingFunction[domain, table, {orders...}]. */
    Expr** der_args = malloc(sizeof(Expr*) * m);
    for (size_t k = 0; k < m; k++) der_args[k] = expr_new_integer(orders[k]);
    Expr* ders_list = expr_new_function(expr_new_symbol("List"), der_args, m);
    free(der_args);
    free(orders);

    Expr* args[3] = { expr_copy(domain), expr_copy(table), ders_list };
    return expr_new_function(expr_new_symbol("InterpolatingFunction"), args, 3);
}

/* --- registration ------------------------------------------------------ */

void interp_init(void) {
    SymbolDef* def = symtab_get_def("InterpolatingFunction");
    def->attributes |= ATTR_PROTECTED;
}

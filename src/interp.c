/*
 * interp.c
 *
 * InterpolatingFunction --- piecewise-polynomial interpolation of tabulated
 * data on a regular (tensor-product) grid, plus the Interpolation[] builder.
 * Modelled on Mathematica's InterpolatingFunction object.
 *
 *   InterpolatingFunction[domain, table]
 *   InterpolatingFunction[domain, table, ders]
 *   InterpolatingFunction[domain, table, ders, orders]
 *   InterpolatingFunction[domain, table, ders, orders, method]
 *
 *     domain = {{x1min, x1max}, ...}   -- one interval per dimension; the
 *              number of intervals m is the dimensionality.
 *     table  = {{coord, val}, ...}                     -- value-only data, or
 *              {{coord, val, grad, hess, ...}, ...}     -- derivative-supplied.
 *              coord is a scalar (1-D value-only) or an {x1,...,xm} list.
 *              grad = D[f,{vars,1}] (length-m vector), hess = D[f,{vars,2}]
 *              (m x m matrix), etc.
 *     ders   = {d1, ..., dm}   -- (optional) derivative-of-interpolant orders.
 *     orders = {o1, ..., om}   -- (optional) interpolation order per dimension.
 *     method = "Spline" | "Hermite"   -- (optional) interpolation method.
 *
 * Methods (all evaluate the ders-th mixed derivative so D[ifun[..],..] composes):
 *   default  : sliding-window Newton divided-difference (order min(3,n-1) or the
 *              requested InterpolationOrder), per dimension, tensor product.
 *   "Spline" : natural cubic spline (C2; second derivative 0 at the ends),
 *              tensor product over the full grid.
 *   "Hermite": tensor-product piecewise cubic Hermite with node slopes estimated
 *              by 3-point finite differences.
 *   supplied : derivative-annotated data is interpolated by tensor-product
 *              Hermite of per-dimension order k = max(K,1) where K is the highest
 *              supplied derivative order.  Mixed partials that are not supplied
 *              are filled by central finite differences across the grid.
 *
 * Precision: machine (double) by default; if the data/argument carry MPFR
 * arbitrary precision the MPFR kernels (interp_mpfr.c) are used instead and an
 * EXPR_MPFR is returned.
 *
 * Builtin ownership: interp_apply / the Interpolation builtin return a fresh
 * Expr* (or NULL to stay unevaluated); inputs are borrowed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "interp.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sym_names.h"
#ifdef USE_MPFR
#include "numeric.h"
#include "common.h"
#endif

/* Interpolation methods. */
enum { METHOD_DEFAULT = 0, METHOD_SPLINE = 1, METHOD_HERMITE = 2 };

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
 *   MPFR (when enabled) and exact Rational[n, d].  Returns false for anything
 *   non-real, leaving *out untouched.
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

/* True for a scalar real atom (Integer/Real/BigInt/MPFR/Rational). */
static bool is_scalar_real(const Expr* e) { double t; return node_to_double(e, &t); }

/* Extract a non-negative integer from `e` into *out. */
static bool node_to_order(const Expr* e, int* out) {
    if (e && e->type == EXPR_INTEGER && e->data.integer >= 0
        && e->data.integer < 1000000) {
        *out = (int)e->data.integer;
        return true;
    }
    return false;
}

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

/* --- 1-D kernels (double) --------------------------------------------- */

/*
 * newton_deriv_eval:
 *   d-th derivative at p of the Newton divided-difference polynomial through
 *   the w window points (xs[k], ys[k]).  d >= w gives 0.
 */
static double newton_deriv_eval(const double* xs, const double* ys,
                                size_t w, double p, size_t d) {
    double* c = malloc(sizeof(double) * w);
    for (size_t k = 0; k < w; k++) c[k] = ys[k];
    for (size_t j = 1; j < w; j++) {
        for (size_t k = w - 1; k >= j; k--) {
            c[k] = (c[k] - c[k - 1]) / (xs[k] - xs[k - j]);
            if (k == j) break;
        }
    }
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

/*
 * spline_eval:
 *   d-th derivative at p of the natural cubic spline through the n points
 *   (xs, ys).  Natural boundary: second derivative 0 at both ends.  For n == 2
 *   this reduces to the straight line.  d >= 4 gives 0.
 */
static double spline_eval(const double* xs, const double* ys,
                          size_t n, double p, size_t d) {
    double* M = calloc(n, sizeof(double));     /* second derivatives at nodes */
    if (n >= 3) {
        size_t ni = n - 2;                     /* interior unknowns M[1..n-2] */
        double* cc = malloc(sizeof(double) * ni);
        double* rr = malloc(sizeof(double) * ni);
        for (size_t j = 0; j < ni; j++) {
            size_t i = j + 1;
            double h0 = xs[i] - xs[i - 1], h1 = xs[i + 1] - xs[i];
            double diag = 2.0 * (h0 + h1), sub = h0, sup = h1;
            double d0 = (ys[i] - ys[i - 1]) / h0, d1 = (ys[i + 1] - ys[i]) / h1;
            double rhs = 6.0 * (d1 - d0);
            if (j == 0) { cc[j] = sup / diag; rr[j] = rhs / diag; }
            else {
                double denom = diag - sub * cc[j - 1];
                cc[j] = sup / denom;
                rr[j] = (rhs - sub * rr[j - 1]) / denom;
            }
        }
        M[n - 2] = rr[ni - 1];
        for (size_t j = ni - 1; j-- > 0; ) M[j + 1] = rr[j] - cc[j] * M[j + 2];
        free(cc); free(rr);
    }
    size_t i = bracket_interval(xs, n, p);
    double h = xs[i + 1] - xs[i];
    double A = xs[i + 1] - p, B = p - xs[i];
    double Mi = M[i], Mi1 = M[i + 1];
    double ci = ys[i] / h - Mi * h / 6.0, ci1 = ys[i + 1] / h - Mi1 * h / 6.0;
    double val;
    switch (d) {
        case 0: val = Mi * A * A * A / (6 * h) + Mi1 * B * B * B / (6 * h)
                      + ci * A + ci1 * B; break;
        case 1: val = -Mi * A * A / (2 * h) + Mi1 * B * B / (2 * h) - ci + ci1; break;
        case 2: val = Mi * A / h + Mi1 * B / h; break;
        case 3: val = (Mi1 - Mi) / h; break;
        default: val = 0.0; break;
    }
    free(M);
    return val;
}

/* --- Hermite tensor-product machinery (double) ------------------------ */

static double dfact(int n) { double r = 1; for (int i = 2; i <= n; i++) r *= i; return r; }
static double falling(int p, int g) { double r = 1; for (int i = 0; i < g; i++) r *= (p - i); return r; }
static size_t ipow_sz(size_t base, size_t e) { size_t r = 1; while (e--) r *= base; return r; }

static void decode_mi(size_t idx, size_t m, int k, int* a) {
    for (size_t d = m; d-- > 0; ) { a[d] = (int)(idx % (size_t)(k + 1)); idx /= (size_t)(k + 1); }
}
static int mi_sum(const int* a, size_t m) { int s = 0; for (size_t d = 0; d < m; d++) s += a[d]; return s; }

/*
 * build_basis:
 *   Coefficients of the 2(k+1) one-dimensional Hermite basis polynomials on
 *   [0,1].  basis[(e*(k+1)+j)*(2k+2) + p] is the t^p coefficient of phi_{e,j},
 *   the degree-(2k+1) polynomial with phi^{(i)}(e') = delta(e,e') delta(j,i)
 *   for e' in {0,1}, i in {0..k}.  Caller frees.
 */
static double* build_basis(int k) {
    int D = 2 * k + 2, K1 = k + 1;
    double* basis = calloc((size_t)2 * K1 * D, sizeof(double));
    for (int e = 0; e < 2; e++) {
        for (int j = 0; j < K1; j++) {
            double* c = basis + ((size_t)(e * K1 + j)) * D;
            /* low coeffs c_0..c_k fixed by the t=0 conditions: c_i = g0_i / i! */
            for (int i = 0; i < K1; i++) c[i] = (e == 0 && i == j) ? 1.0 / dfact(i) : 0.0;
            /* solve (K1 x K1) system for the high coeffs c_{k+1..2k+1} from the
             * t=1 conditions P^{(i)}(1) = g1_i. */
            double A[16][16], rhs[16];
            for (int i = 0; i < K1; i++) {
                double s = 0;
                for (int p = i; p <= k; p++) s += c[p] * falling(p, i);
                double g1 = (e == 1 && i == j) ? 1.0 : 0.0;
                rhs[i] = g1 - s;
                for (int q = 0; q < K1; q++) { int p = k + 1 + q; A[i][q] = (p >= i) ? falling(p, i) : 0.0; }
            }
            for (int col = 0; col < K1; col++) {
                int piv = col;
                for (int r = col + 1; r < K1; r++) if (fabs(A[r][col]) > fabs(A[piv][col])) piv = r;
                if (piv != col) {
                    for (int cc = 0; cc < K1; cc++) { double t = A[col][cc]; A[col][cc] = A[piv][cc]; A[piv][cc] = t; }
                    double t = rhs[col]; rhs[col] = rhs[piv]; rhs[piv] = t;
                }
                for (int r = 0; r < K1; r++) if (r != col) {
                    double f = A[r][col] / A[col][col];
                    for (int cc = col; cc < K1; cc++) A[r][cc] -= f * A[col][cc];
                    rhs[r] -= f * rhs[col];
                }
            }
            for (int q = 0; q < K1; q++) c[k + 1 + q] = rhs[q] / A[q][q];
        }
    }
    return basis;
}

/* g-th derivative at t of the polynomial with coefficients c[0..deg]. */
static double phi_eval(const double* c, int deg, int g, double t) {
    double s = 0;
    for (int p = deg; p >= g; p--) s = s * t + c[p] * falling(p, g);
    return s;
}

/* --- parsed grid ------------------------------------------------------- */

typedef struct {
    size_t  m;
    size_t* nk;
    double** grid;
    size_t* stride;
    size_t  total;
    double* dmin;
    double* dmax;
    bool*   has_range;
    Expr**  entryAt;     /* data entry sitting at each flat grid node */
    double* V;           /* flat value tensor (value-only kernels), or NULL */
} IFun;

static void ifun_free(IFun* f) {
    if (f->grid) { for (size_t k = 0; k < f->m; k++) free(f->grid[k]); free(f->grid); }
    free(f->nk); free(f->stride); free(f->V);
    free(f->dmin); free(f->dmax); free(f->has_range); free(f->entryAt);
    memset(f, 0, sizeof(*f));
}

static Expr* entry_value(Expr* entry) { return entry->data.function.args[1]; }

/* Coordinate k (0-based) of a data entry, given dimensionality m.  Accepts a
 * scalar coord (1-D value-only) or an {x1,...,xm} leading list. */
static Expr* entry_coord(Expr* entry, size_t m, size_t k) {
    Expr* point = entry->data.function.args[0];
    if (interp_is_list(point) && point->data.function.arg_count == m)
        return point->data.function.args[k];
    if (m == 1) return point;
    return NULL;
}

static void grid_insert(double* g, size_t* len, double v) {
    size_t i = 0;
    while (i < *len && g[i] < v) i++;
    if (i < *len && g[i] == v) return;
    for (size_t j = *len; j > i; j--) g[j] = g[j - 1];
    g[i] = v;
    (*len)++;
}

static size_t grid_index(const double* g, size_t len, double v) {
    size_t lo = 0, hi = len;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (g[mid] < v) lo = mid + 1; else hi = mid;
    }
    return (lo < len && g[lo] == v) ? lo : len;
}

static size_t parse_domain_dim(Expr* domain) {
    if (!interp_is_list(domain)) return 0;
    return domain->data.function.arg_count;
}

/*
 * build_grid:
 *   Parse `domain` and `table` into the tensor grid (abscissae, strides, node
 *   -> entry map, domain bounds).  Validates that coordinates are numeric and
 *   that the points fill the full Cartesian product.  Does NOT inspect values
 *   or derivative tensors.  Returns true on success (no leaks on failure).
 */
static bool build_grid(Expr* domain, Expr* table, size_t m, IFun* out) {
    memset(out, 0, sizeof(*out));
    if (!interp_is_list(table)) return false;
    size_t npts = table->data.function.arg_count;
    if (npts < 2) return false;

    for (size_t i = 0; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        if (!interp_is_list(e) || e->data.function.arg_count < 2) return false;
        double tmp;
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
    Expr**   eat    = NULL;
    if (!nk || !grid || !stride || !dmin || !dmax || !hasr) goto fail;

    for (size_t k = 0; k < m; k++) {
        grid[k] = malloc(sizeof(double) * npts);
        if (!grid[k]) goto fail;
        size_t len = 0;
        for (size_t i = 0; i < npts; i++) {
            double v;
            node_to_double(entry_coord(table->data.function.args[i], m, k), &v);
            grid_insert(grid[k], &len, v);
        }
        nk[k] = len;
        if (len < 2) goto fail;
    }

    size_t total = 1;
    for (size_t k = 0; k < m; k++) total *= nk[k];
    stride[m - 1] = 1;
    for (size_t k = m - 1; k-- > 0; ) stride[k] = stride[k + 1] * nk[k + 1];

    eat = calloc(total, sizeof(Expr*));
    if (!eat) goto fail;
    for (size_t i = 0; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        size_t flat = 0;
        for (size_t k = 0; k < m; k++) {
            double v;
            node_to_double(entry_coord(e, m, k), &v);
            flat += grid_index(grid[k], nk[k], v) * stride[k];
        }
        eat[flat] = e;
    }
    for (size_t i = 0; i < total; i++) if (!eat[i]) goto fail;  /* incomplete grid */

    for (size_t k = 0; k < m; k++) {
        Expr* iv = domain->data.function.args[k];
        if (interp_is_list(iv) && iv->data.function.arg_count == 2
            && node_to_double(iv->data.function.args[0], &dmin[k])
            && node_to_double(iv->data.function.args[1], &dmax[k]))
            hasr[k] = true;
    }

    out->m = m; out->nk = nk; out->grid = grid; out->stride = stride;
    out->total = total; out->dmin = dmin; out->dmax = dmax; out->has_range = hasr;
    out->entryAt = eat;
    return true;

fail:
    if (grid) { for (size_t k = 0; k < m; k++) free(grid[k]); free(grid); }
    free(nk); free(stride); free(dmin); free(dmax); free(hasr); free(eat);
    return false;
}

/* Fill the flat value tensor (node value) for the value-only kernels. */
static bool fill_values(IFun* f) {
    f->V = malloc(sizeof(double) * f->total);
    if (!f->V) return false;
    for (size_t i = 0; i < f->total; i++)
        if (!node_to_double(entry_value(f->entryAt[i]), &f->V[i])) return false;
    return true;
}

/* --- value-only tensor evaluation (default / spline) ------------------ */

typedef struct {
    const IFun*   f;
    const double* p;
    const int*    ders;
    const size_t* s;
    const size_t* w;
    int           kernel;   /* 0 = newton, 1 = spline */
} EvalCtx;

static double eval_dim(const EvalCtx* ctx, size_t k, size_t base) {
    size_t wk = ctx->w[k], sk = ctx->s[k];
    double* tmp = malloc(sizeof(double) * wk);
    for (size_t t = 0; t < wk; t++) {
        size_t off = base + (sk + t) * ctx->f->stride[k];
        tmp[t] = (k + 1 == ctx->f->m) ? ctx->f->V[off] : eval_dim(ctx, k + 1, off);
    }
    double res = ctx->kernel == 1
        ? spline_eval(&ctx->f->grid[k][sk], tmp, wk, ctx->p[k], (size_t)ctx->ders[k])
        : newton_deriv_eval(&ctx->f->grid[k][sk], tmp, wk, ctx->p[k], (size_t)ctx->ders[k]);
    free(tmp);
    return res;
}

/* --- supplied/estimated mixed-derivative tensor (Hermite) ------------- */

/* D^a f at the node held in `entry` (a has total order s = mi_sum(a)); read
 * from the supplied s-th derivative tensor (s==0 -> the value). */
static bool extract_supplied(Expr* entry, size_t m, const int* a, int s, double* out) {
    if ((size_t)(1 + s) >= entry->data.function.arg_count) return false;
    Expr* t = entry->data.function.args[1 + s];
    if (m == 1) {
        /* 1-D derivatives are scalars (D[f,x], D[f,{x,2}], ...).  Tolerate the
         * occasional singleton-list wrapping {df}. */
        while (interp_is_list(t) && t->data.function.arg_count == 1) t = t->data.function.args[0];
        return node_to_double(t, out);
    }
    /* m >= 2: navigate the rank-s derivative tensor with the index sequence of a
     * (dimension d repeated a_d times). */
    for (size_t d = 0; d < m; d++) {
        for (int r = 0; r < a[d]; r++) {
            if (!interp_is_list(t) || t->data.function.arg_count != m) return false;
            t = t->data.function.args[d];
        }
    }
    return node_to_double(t, out);
}

/* Estimate d/dx_d of the field F over the grid (3-point parabola interior,
 * one-sided at the ends), writing into G. */
static void fd_apply(const IFun* g, size_t d, const double* F, double* G) {
    size_t nd = g->nk[d], st = g->stride[d];
    const double* xs = g->grid[d];
    for (size_t n = 0; n < g->total; n++) {
        size_t id = (n / st) % nd;
        size_t base = n - id * st;
        double val;
        if (nd == 1) { val = 0; }
        else if (nd == 2) { val = (F[base + st] - F[base]) / (xs[1] - xs[0]); }
        else if (id == 0) {
            double h0 = xs[1] - xs[0], h1 = xs[2] - xs[1];
            double d0 = (F[base + st] - F[base]) / h0;
            double d1 = (F[base + 2 * st] - F[base + st]) / h1;
            val = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        } else if (id == nd - 1) {
            double h0 = xs[nd - 2] - xs[nd - 3], h1 = xs[nd - 1] - xs[nd - 2];
            double d0 = (F[base + (nd - 2) * st] - F[base + (nd - 3) * st]) / h0;
            double d1 = (F[base + (nd - 1) * st] - F[base + (nd - 2) * st]) / h1;
            val = ((2 * h1 + h0) * d1 - h1 * d0) / (h0 + h1);
        } else {
            double h0 = xs[id] - xs[id - 1], h1 = xs[id + 1] - xs[id];
            double d0 = (F[base + id * st] - F[base + (id - 1) * st]) / h0;
            double d1 = (F[base + (id + 1) * st] - F[base + id * st]) / h1;
            val = (h1 * d0 + h0 * d1) / (h0 + h1);
        }
        G[n] = val;
    }
}

/*
 * build_T:
 *   Full mixed-derivative tensor T[node*dcount + idx], idx enumerating the
 *   multi-indices a in {0..k}^m.  Supplied (sum(a) <= K) entries are read from
 *   the data; the rest are estimated by finite-differencing a supplied
 *   order-K derivative.  Returns NULL on malformed data.
 */
static double* build_T(const IFun* g, size_t m, int K, int k) {
    size_t dcount = ipow_sz((size_t)(k + 1), m);
    double* T = malloc(sizeof(double) * g->total * dcount);
    int* a = malloc(sizeof(int) * m);
    int* b = malloc(sizeof(int) * m);
    double* F = malloc(sizeof(double) * g->total);
    double* Gtmp = malloc(sizeof(double) * g->total);
    bool ok = T && a && b && F && Gtmp;
    for (size_t idx = 0; ok && idx < dcount; idx++) {
        decode_mi(idx, m, k, a);
        int s = mi_sum(a, m);
        if (s <= K) {
            for (size_t n = 0; n < g->total; n++) {
                double v;
                if (!extract_supplied(g->entryAt[n], m, a, s, &v)) { ok = false; break; }
                T[n * dcount + idx] = v;
            }
        } else {
            int rem = K;
            for (size_t d = 0; d < m; d++) { int bd = a[d] < rem ? a[d] : rem; b[d] = bd; rem -= bd; }
            for (size_t n = 0; n < g->total; n++) {
                double v;
                if (!extract_supplied(g->entryAt[n], m, b, K, &v)) { ok = false; break; }
                F[n] = v;
            }
            if (!ok) break;
            for (size_t d = 0; d < m; d++) {
                int e = a[d] - b[d];
                for (int r = 0; r < e; r++) { fd_apply(g, d, F, Gtmp); double* sw = F; F = Gtmp; Gtmp = sw; }
            }
            for (size_t n = 0; n < g->total; n++) T[n * dcount + idx] = F[n];
        }
    }
    free(a); free(b); free(F); free(Gtmp);
    if (!ok) { free(T); return NULL; }
    return T;
}

/*
 * hermite_tensor_eval:
 *   Evaluate the ders-th mixed derivative of the tensor-product Hermite
 *   interpolant at p, using the mixed-derivative tensor T and the 1-D basis.
 */
static double hermite_tensor_eval(const IFun* g, const double* T, size_t dcount,
                                  size_t m, int k, const double* basis, int D,
                                  const double* p, const int* ders) {
    int* lo = malloc(sizeof(int) * m);
    double* h = malloc(sizeof(double) * m);
    double* t = malloc(sizeof(double) * m);
    double* factor = malloc(sizeof(double) * m * 2 * (size_t)(k + 1));
    int* a = malloc(sizeof(int) * m);
    int* c = malloc(sizeof(int) * m);

    for (size_t d = 0; d < m; d++) {
        size_t i = bracket_interval(g->grid[d], g->nk[d], p[d]);
        lo[d] = (int)i;
        h[d] = g->grid[d][i + 1] - g->grid[d][i];
        t[d] = (p[d] - g->grid[d][i]) / h[d];
        for (int cd = 0; cd < 2; cd++)
            for (int ad = 0; ad <= k; ad++) {
                const double* coeffs = basis + ((size_t)(cd * (k + 1) + ad)) * D;
                double phi = phi_eval(coeffs, 2 * k + 1, ders[d], t[d]);
                factor[(d * 2 + cd) * (size_t)(k + 1) + ad] = phi * pow(h[d], (double)(ad - ders[d]));
            }
    }

    double result = 0;
    size_t ncorner = (size_t)1 << m;
    for (size_t cc = 0; cc < ncorner; cc++) {
        for (size_t d = 0; d < m; d++) c[d] = (int)((cc >> d) & 1);
        size_t node = 0;
        for (size_t d = 0; d < m; d++) node += (size_t)(lo[d] + c[d]) * g->stride[d];
        for (size_t idx = 0; idx < dcount; idx++) {
            double term = T[node * dcount + idx];
            if (term == 0.0) continue;
            decode_mi(idx, m, k, a);
            for (size_t d = 0; d < m; d++)
                term *= factor[(d * 2 + c[d]) * (size_t)(k + 1) + a[d]];
            result += term;
        }
    }
    free(lo); free(h); free(t); free(factor); free(a); free(c);
    return result;
}

/*
 * interp_eval_double:
 *   Evaluate the object at the (machine-precision) point.  Returns a fresh
 *   EXPR_REAL, or NULL if the point is non-numeric / the data is malformed.
 *   `Ksupplied` is the number of supplied derivative tensors per node (0 for
 *   value-only data).
 */
static Expr* interp_eval_double(Expr* domain, Expr* table, size_t m,
                                Expr** call_args, const int* ders,
                                const int* orders, int method, int Ksupplied) {
    double* p = malloc(sizeof(double) * m);
    for (size_t k = 0; k < m; k++)
        if (!node_to_double(call_args[k], &p[k])) { free(p); return NULL; }

    IFun f;
    if (!build_grid(domain, table, m, &f)) { free(p); return NULL; }

    for (size_t k = 0; k < m; k++)
        if (f.has_range[k] && (p[k] < f.dmin[k] || p[k] > f.dmax[k])) {
            fprintf(stderr,
                    "InterpolatingFunction::dmval: Input value %g lies outside "
                    "the range of data in the interpolating function. "
                    "Extrapolation will be used.\n", p[k]);
            break;
        }

    double value;
    bool ok = true;

    if (method == METHOD_HERMITE || Ksupplied >= 1) {
        int K = Ksupplied;
        int k = K > 1 ? K : 1;
        double* basis = build_basis(k);
        double* T = build_T(&f, m, K, k);
        if (!T) { ok = false; value = 0; }
        else {
            value = hermite_tensor_eval(&f, T, ipow_sz((size_t)(k + 1), m),
                                        m, k, basis, 2 * k + 2, p, ders);
            free(T);
        }
        free(basis);
    } else if (method == METHOD_SPLINE) {
        if (!fill_values(&f)) { ok = false; value = 0; }
        else {
            size_t* s = calloc(m, sizeof(size_t));
            size_t* w = malloc(sizeof(size_t) * m);
            for (size_t k = 0; k < m; k++) w[k] = f.nk[k];   /* full grid */
            EvalCtx ctx = { &f, p, ders, s, w, 1 };
            value = eval_dim(&ctx, 0, 0);
            free(s); free(w);
        }
    } else {
        if (!fill_values(&f)) { ok = false; value = 0; }
        else {
            size_t* s = malloc(sizeof(size_t) * m);
            size_t* w = malloc(sizeof(size_t) * m);
            for (size_t k = 0; k < m; k++) {
                size_t nkk = f.nk[k];
                size_t order;
                if (orders) { order = (size_t)orders[k]; if (order > nkk - 1) order = nkk - 1; }
                else { order = nkk - 1; if (order > 3) order = 3; }
                size_t wk = order + 1;
                size_t i = bracket_interval(f.grid[k], nkk, p[k]);
                size_t shift = (wk >= 2) ? (wk / 2 - 1) : 0;
                size_t sk = (i < shift) ? 0 : i - shift;
                if (sk > nkk - wk) sk = nkk - wk;
                s[k] = sk; w[k] = wk;
            }
            EvalCtx ctx = { &f, p, ders, s, w, 0 };
            value = eval_dim(&ctx, 0, 0);
            free(s); free(w);
        }
    }

    ifun_free(&f);
    free(p);
    return ok ? expr_new_real(value) : NULL;
}

#ifdef USE_MPFR
/* Implemented in interp_mpfr.c. */
Expr* interp_eval_mpfr(Expr* domain, Expr* table, size_t m,
                       Expr** call_args, const int* ders,
                       const int* orders, int method, int Ksupplied, long bits);
#endif

/* --- object inspection (precision + table form) ----------------------- */

/* Number of supplied derivative tensors per node (0 = value-only).  Returns -1
 * if entries are malformed or non-uniform in length. */
static int table_Ksupplied(Expr* table) {
    if (!interp_is_list(table) || table->data.function.arg_count < 1) return -1;
    size_t npts = table->data.function.arg_count;
    Expr* e0 = table->data.function.args[0];
    if (!interp_is_list(e0)) return -1;
    size_t L = e0->data.function.arg_count;
    for (size_t i = 1; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        if (!interp_is_list(e) || e->data.function.arg_count != L) return -1;
    }
    if (L < 2) return -1;
    return (int)(L - 2);
}

/* --- public entry points ---------------------------------------------- */

Expr* interp_apply(Expr* ifun, Expr** call_args, size_t argc) {
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc < 2 || obj_argc > 5) return NULL;
    Expr* domain = ifun->data.function.args[0];
    Expr* table  = ifun->data.function.args[1];

    size_t m = parse_domain_dim(domain);
    if (m == 0 || argc != m) return NULL;

    int* ders = calloc(m, sizeof(int));
    bool any_der = false;
    if (obj_argc >= 3) {
        Expr* ds = ifun->data.function.args[2];
        if (!interp_is_list(ds) || ds->data.function.arg_count != m) { free(ders); return NULL; }
        for (size_t k = 0; k < m; k++) {
            if (!node_to_order(ds->data.function.args[k], &ders[k])) { free(ders); return NULL; }
            if (ders[k] > 0) any_der = true;
        }
    }

    int* orders = NULL;
    if (obj_argc >= 4) {
        Expr* os = ifun->data.function.args[3];
        if (!interp_is_list(os) || os->data.function.arg_count != m) { free(ders); return NULL; }
        orders = calloc(m, sizeof(int));
        for (size_t k = 0; k < m; k++)
            if (!node_to_order(os->data.function.args[k], &orders[k])) { free(ders); free(orders); return NULL; }
    }

    int method = METHOD_DEFAULT;
    if (obj_argc == 5) {
        Expr* me = ifun->data.function.args[4];
        if (me->type == EXPR_STRING && me->data.string) {
            if (strcmp(me->data.string, "Spline") == 0) method = METHOD_SPLINE;
            else if (strcmp(me->data.string, "Hermite") == 0) method = METHOD_HERMITE;
        }
    }

    int Ksupplied = table_Ksupplied(table);
    if (Ksupplied < 0) { free(ders); free(orders); return NULL; }

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
                for (size_t k = 0; k < m; k++)
                    if (!expr_eq(call_args[k], entry_coord(e, m, k))) { hit = false; break; }
                if (hit) {
                    Expr* y = expr_copy(entry_value(e));
                    free(ders); free(orders);
                    return y;
                }
            }
        }
    }

    Expr* result;
#ifdef USE_MPFR
    /* Route to the MPFR kernels when the data or the argument carry
     * arbitrary precision. */
    bool want_mpfr = false;
    long bits = 53;
    {
        CommonInexactInfo di = common_scan_inexact(table);
        CommonInexactInfo ai = common_scan_inexact(domain);
        long mn = 0;
        bool any = false;
        if (di.has_inexact) { any = true; mn = di.min_bits; }
        if (ai.has_inexact) { any = true; mn = (!any || ai.min_bits < mn) ? ai.min_bits : mn; mn = (di.has_inexact && di.min_bits < mn) ? di.min_bits : mn; }
        for (size_t k = 0; k < m; k++) {
            CommonInexactInfo ci = common_scan_inexact(call_args[k]);
            if (ci.has_inexact) { if (!any || ci.min_bits < mn) mn = ci.min_bits; any = true; }
        }
        if ((numeric_expr_is_mpfr(table) || numeric_expr_is_mpfr(domain)) ) want_mpfr = true;
        for (size_t k = 0; k < m; k++) if (numeric_expr_is_mpfr(call_args[k])) want_mpfr = true;
        if (want_mpfr) bits = (mn > 53) ? mn : 53;
    }
    if (want_mpfr)
        result = interp_eval_mpfr(domain, table, m, call_args, ders, orders, method, Ksupplied, bits);
    else
#endif
        result = interp_eval_double(domain, table, m, call_args, ders, orders, method, Ksupplied);

    free(ders); free(orders);
    return result;
}

Expr* interp_make_derivative(Expr* deriv_head, Expr* ifun) {
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc < 2 || obj_argc > 5) return NULL;
    Expr* domain = ifun->data.function.args[0];
    Expr* table  = ifun->data.function.args[1];

    size_t m = parse_domain_dim(domain);
    if (m == 0) return NULL;
    if (deriv_head->data.function.arg_count != m) return NULL;

    int* ders = calloc(m, sizeof(int));
    if (!ders) return NULL;
    if (obj_argc >= 3) {
        Expr* ds = ifun->data.function.args[2];
        if (!interp_is_list(ds) || ds->data.function.arg_count != m) { free(ders); return NULL; }
        for (size_t k = 0; k < m; k++)
            if (!node_to_order(ds->data.function.args[k], &ders[k])) { free(ders); return NULL; }
    }
    for (size_t k = 0; k < m; k++) {
        int inc;
        if (!node_to_order(deriv_head->data.function.args[k], &inc)) { free(ders); return NULL; }
        ders[k] += inc;
    }

    Expr** der_args = malloc(sizeof(Expr*) * m);
    for (size_t k = 0; k < m; k++) der_args[k] = expr_new_integer(ders[k]);
    Expr* ders_list = expr_new_function(expr_new_symbol("List"), der_args, m);
    free(der_args);
    free(ders);

    size_t out_argc = (obj_argc >= 4) ? obj_argc : 3;
    Expr** args = malloc(sizeof(Expr*) * out_argc);
    args[0] = expr_copy(domain);
    args[1] = expr_copy(table);
    args[2] = ders_list;
    for (size_t k = 3; k < out_argc; k++) args[k] = expr_copy(ifun->data.function.args[k]);
    Expr* result = expr_new_function(expr_new_symbol("InterpolatingFunction"), args, out_argc);
    free(args);
    return result;
}

/* --- Interpolation builtin -------------------------------------------- */

/* List[a, b], taking ownership of a and b. */
static Expr* make_pair(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol("List"), args, 2);
}

/* List from an owned array of `n` exprs (array borrowed; caller frees). */
static Expr* make_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol("List"), items, n);
}

static void free_exprs(Expr** xs, size_t n) {
    for (size_t i = 0; i < n; i++) expr_free(xs[i]);
    free(xs);
}

static Expr* builtin_interpolation(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    Expr* data = res->data.function.args[0];
    if (!interp_is_list(data)) return NULL;
    size_t npts = data->data.function.arg_count;
    if (npts < 2) return NULL;

    /* options / immediate evaluation point */
    Expr* eval_point = NULL;
    int explicit_order = -1;
    int method = METHOD_DEFAULT;
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_FUNCTION && a->data.function.arg_count == 2
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol == SYM_Rule
                || a->data.function.head->data.symbol == SYM_RuleDelayed)) {
            Expr* lhs = a->data.function.args[0];
            Expr* rhs = a->data.function.args[1];
            if (lhs->type != EXPR_SYMBOL) return NULL;
            if (lhs->data.symbol == SYM_InterpolationOrder) {
                if (!node_to_order(rhs, &explicit_order)) return NULL;
            } else if (lhs->data.symbol == SYM_Method) {
                if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) { /* default */ }
                else if (rhs->type == EXPR_STRING && rhs->data.string
                         && strcmp(rhs->data.string, "Spline") == 0) method = METHOD_SPLINE;
                else if (rhs->type == EXPR_STRING && rhs->data.string
                         && strcmp(rhs->data.string, "Hermite") == 0) method = METHOD_HERMITE;
                else return NULL;
            } else if (lhs->data.symbol == SYM_PeriodicInterpolation) {
                if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_False)) return NULL;
            } else return NULL;
            continue;
        }
        if (eval_point) return NULL;
        eval_point = a;
    }

    /* determine the data form */
    Expr* e0 = data->data.function.args[0];
    size_t m;
    bool synth_x = false;     /* form 1 */
    int Ksupplied = 0;
    if (is_scalar_real(e0)) {
        m = 1; synth_x = true;
    } else if (interp_is_list(e0)) {
        size_t L = e0->data.function.arg_count;
        if (L < 2) return NULL;
        Expr* first = e0->data.function.args[0];
        if (L == 2 && is_scalar_real(first)) { m = 1; }
        else if (interp_is_list(first)) {
            m = first->data.function.arg_count;
            if (m < 1) return NULL;
            Ksupplied = (int)(L - 2);   /* L==2 -> value-only n-D; L>2 -> derivatives */
        } else return NULL;
    } else return NULL;

    /* build the {coord,...} table and the per-dimension domain bounds */
    Expr** entries = malloc(sizeof(Expr*) * npts);
    double* dmin = malloc(sizeof(double) * m);
    double* dmax = malloc(sizeof(double) * m);
    Expr**  dminE = malloc(sizeof(Expr*) * m);
    Expr**  dmaxE = malloc(sizeof(Expr*) * m);
    if (!entries || !dmin || !dmax || !dminE || !dmaxE) {
        free(entries); free(dmin); free(dmax); free(dminE); free(dmaxE); return NULL;
    }

    bool ok = true;
    size_t built = 0;
    for (size_t i = 0; i < npts && ok; i++) {
        Expr* e = data->data.function.args[i];
        double coord[64];
        Expr* coordE[64];
        if (m > 64) { ok = false; break; }

        if (synth_x) {
            if (!is_scalar_real(e)) { ok = false; break; }
            coord[0] = (double)(i + 1); coordE[0] = NULL;
        } else {
            if (!interp_is_list(e) || e->data.function.arg_count < 2) { ok = false; break; }
            if ((int)(e->data.function.arg_count - 2) != Ksupplied) { ok = false; break; }
            if (!is_scalar_real(e->data.function.args[1])) { ok = false; break; }  /* value */
            for (size_t k = 0; k < m; k++) {
                Expr* ck = entry_coord(e, m, k);
                if (!ck || !node_to_double(ck, &coord[k])) { ok = false; break; }
                coordE[k] = ck;
            }
            if (!ok) break;
        }

        for (size_t k = 0; k < m; k++) {
            if (i == 0) { dmin[k] = dmax[k] = coord[k]; dminE[k] = dmaxE[k] = coordE[k]; }
            else {
                if (coord[k] < dmin[k]) { dmin[k] = coord[k]; dminE[k] = coordE[k]; }
                if (coord[k] > dmax[k]) { dmax[k] = coord[k]; dmaxE[k] = coordE[k]; }
            }
        }

        if (synth_x) entries[i] = make_pair(expr_new_integer((int64_t)(i + 1)), expr_copy(e));
        else         entries[i] = expr_copy(e);
        built = i + 1;
    }

    if (!ok) {
        free_exprs(entries, built);
        free(dmin); free(dmax); free(dminE); free(dmaxE);
        return NULL;
    }

    Expr** dom_items = malloc(sizeof(Expr*) * m);
    for (size_t k = 0; k < m; k++) {
        Expr* lo = dminE[k] ? expr_copy(dminE[k]) : expr_new_integer((int64_t)(dmin[k] + 0.5));
        Expr* hi = dmaxE[k] ? expr_copy(dmaxE[k]) : expr_new_integer((int64_t)(dmax[k] + 0.5));
        dom_items[k] = make_pair(lo, hi);
    }
    Expr* domain = make_list(dom_items, m);
    free(dom_items);
    Expr* table = make_list(entries, npts);
    free(entries);
    free(dmin); free(dmax); free(dminE); free(dmaxE);

    /* assemble the object */
    Expr* object;
    bool need_full = (explicit_order >= 0) || (method != METHOD_DEFAULT);
    if (!need_full) {
        Expr* oargs[2] = { domain, table };
        object = expr_new_function(expr_new_symbol("InterpolatingFunction"), oargs, 2);
    } else {
        int ov = (explicit_order >= 0) ? explicit_order : 3;
        Expr** zeros = malloc(sizeof(Expr*) * m);
        Expr** ords  = malloc(sizeof(Expr*) * m);
        for (size_t k = 0; k < m; k++) { zeros[k] = expr_new_integer(0); ords[k] = expr_new_integer(ov); }
        Expr* ders_list  = make_list(zeros, m);
        Expr* order_list = make_list(ords, m);
        free(zeros); free(ords);
        if (method == METHOD_DEFAULT) {
            Expr* oargs[4] = { domain, table, ders_list, order_list };
            object = expr_new_function(expr_new_symbol("InterpolatingFunction"), oargs, 4);
        } else {
            Expr* mstr = expr_new_string(method == METHOD_SPLINE ? "Spline" : "Hermite");
            Expr* oargs[5] = { domain, table, ders_list, order_list, mstr };
            object = expr_new_function(expr_new_symbol("InterpolatingFunction"), oargs, 5);
        }
    }

    if (eval_point) {
        Expr* value;
        if (m == 1) {
            if (!is_scalar_real(eval_point)) { expr_free(object); return NULL; }
            value = interp_apply(object, &eval_point, 1);
        } else {
            if (!interp_is_list(eval_point) || eval_point->data.function.arg_count != m) {
                expr_free(object); return NULL;
            }
            value = interp_apply(object, eval_point->data.function.args, m);
        }
        expr_free(object);
        return value;
    }
    return object;
}

/* --- registration ------------------------------------------------------ */

void interp_init(void) {
    SymbolDef* def = symtab_get_def("InterpolatingFunction");
    def->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Interpolation", builtin_interpolation);
    SymbolDef* idef = symtab_get_def("Interpolation");
    idef->attributes |= ATTR_PROTECTED;
}

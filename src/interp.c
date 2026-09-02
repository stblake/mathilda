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
#include "ndarray.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "eval.h"         /* eval_and_free */
#include "zero_test.h"    /* zero_test_decide (Engine B consistency) */
#include "common.h"       /* builtin_arg_error */
#include "print.h"        /* expr_to_string (::poised message) */
#include <gmp.h>          /* Modulus path: exact modular linear solve */
#ifdef USE_MPFR
#include "numeric.h"
#endif

/* Interpolation methods. */
enum { METHOD_DEFAULT = 0, METHOD_SPLINE = 1, METHOD_HERMITE = 2 };

/* --- small numeric / structural helpers ------------------------------- */

/* True when `e` is List[...] (head is the List symbol). */
static bool interp_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
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
        && e->data.function.head->data.symbol.name == SYM_Rational
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

/* True for a scalar real or a (recursively) non-empty list of reals/arrays. */
static bool is_real_or_array(const Expr* e) {
    if (is_scalar_real(e)) return true;
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count > 0) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!is_real_or_array(e->data.function.args[i])) return false;
        return true;
    }
    return false;
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

/*
 * Locate the bracketing interval for `x`: the largest index i with
 * xs[i] <= x, clamped to [0, n-2].  xs is strictly increasing, n >= 2.
 */
/* Largest i in [0, n-2] with xs[i] <= x < xs[i+1].  Branchless exponential
 * search: the data-dependent comparison feeds an arithmetic step rather than a
 * conditional branch, so a batch of UNPREDICTABLE (e.g. random) query points
 * costs no misprediction penalty.  On a 2000-knot grid this took a random 10^5
 * evaluation from ~260 ns/point to a few tens of ns/point (a naive
 * if/else binary search mispredicts ~half its iterations). */
static size_t bracket_interval(const double* xs, size_t n, double x) {
    if (x <= xs[0]) return 0;
    if (x >= xs[n - 1]) return n - 2;
    size_t pos = 0;
    size_t step = 1;
    while ((step << 1) <= n - 1) step <<= 1;   /* largest power of two <= n-1 */
    for (; step > 0; step >>= 1) {
        size_t next = pos + step;
        size_t adv = (next < n && xs[next] <= x) ? 1u : 0u;
        pos += adv * step;                     /* branchless advance */
    }
    return pos;
}

/* --- 1-D kernels (double) --------------------------------------------- */

/*
 * newton_deriv_eval:
 *   d-th derivative at p of the Newton divided-difference polynomial through
 *   the w window points (xs[k], ys[k]).  d >= w gives 0.
 */
/* Core, with caller-provided scratch: `c` (>= w doubles) and `val` (>= d+1
 * doubles).  Split out so a vectorised loop allocates the scratch once instead
 * of two malloc/free per query point (the dominant cost of a batch evaluation).*/
static double newton_deriv_eval_buf(const double* xs, const double* ys, size_t w,
                                    double p, size_t d, double* c, double* val) {
    for (size_t k = 0; k < w; k++) c[k] = ys[k];
    for (size_t j = 1; j < w; j++) {
        for (size_t k = w - 1; k >= j; k--) {
            c[k] = (c[k] - c[k - 1]) / (xs[k] - xs[k - j]);
            if (k == j) break;
        }
    }
    for (size_t t = 0; t <= d; t++) val[t] = 0.0;
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
    return fact * val[d];
}

static double newton_deriv_eval(const double* xs, const double* ys,
                                size_t w, double p, size_t d) {
    double* c = malloc(sizeof(double) * w);
    double* val = malloc(sizeof(double) * (d + 1));
    double result = newton_deriv_eval_buf(xs, ys, w, p, d, c, val);
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
/* Natural cubic spline node second derivatives (M[0] = M[n-1] = 0), by the
 * Thomas algorithm, O(n).  Returns a malloc'd M[n] (caller frees).  Split out
 * from spline_eval so a vectorised evaluation can solve the system ONCE and
 * reuse M across every query point instead of re-solving per point. */
static double* spline_solve_M(const double* xs, const double* ys, size_t n) {
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
    return M;
}

/* d-th derivative at p of the natural cubic spline with precomputed node second
 * derivatives M (from spline_solve_M).  Binary-search interval, O(1) evaluate. */
static double spline_piece_eval(const double* xs, const double* ys, const double* M,
                                size_t n, double p, size_t d) {
    size_t i = bracket_interval(xs, n, p);
    double h = xs[i + 1] - xs[i];
    double A = xs[i + 1] - p, B = p - xs[i];
    double Mi = M[i], Mi1 = M[i + 1];
    double ci = ys[i] / h - Mi * h / 6.0, ci1 = ys[i + 1] / h - Mi1 * h / 6.0;
    switch (d) {
        case 0: return Mi * A * A * A / (6 * h) + Mi1 * B * B * B / (6 * h)
                       + ci * A + ci1 * B;
        case 1: return -Mi * A * A / (2 * h) + Mi1 * B * B / (2 * h) - ci + ci1;
        case 2: return Mi * A / h + Mi1 * B / h;
        case 3: return (Mi1 - Mi) / h;
        default: return 0.0;
    }
}

static double spline_eval(const double* xs, const double* ys,
                          size_t n, double p, size_t d) {
    double* M = spline_solve_M(xs, ys, n);
    double val = spline_piece_eval(xs, ys, M, n, p, d);
    free(M);
    return val;
}

/*
 * spline_eval_periodic:
 *   d-th derivative at p of the periodic (cyclic) cubic spline through the n
 *   distinct nodes (xs[0..n-1], ys[0..n-1]) with period P (so x_n = x_0 + P,
 *   y_n = y_0).  The C2 periodic conditions give a cyclic tridiagonal system for
 *   the node second derivatives, solved by Sherman-Morrison.  p is assumed
 *   already reduced into [x_0, x_0 + P).
 */
/* Periodic (cyclic) cubic spline node second derivatives, cyclic-tridiagonal via
 * Sherman-Morrison, O(n).  Returns malloc'd M[n] (caller frees).  Split from
 * spline_eval_periodic so the vectorised path solves it once. */
static double* spline_solve_M_periodic(const double* xs, const double* ys,
                                       size_t n, double P) {
    /* h[i] = x_{i+1} - x_i, with the wrap interval h[n-1] = (x_0 + P) - x_{n-1}. */
    double* h = malloc(sizeof(double) * n);
    for (size_t i = 0; i + 1 < n; i++) h[i] = xs[i + 1] - xs[i];
    h[n - 1] = (xs[0] + P) - xs[n - 1];

    double* M;
    if (n < 3) {
        /* Degenerate: treat as piecewise linear (M = 0). */
        M = calloc(n, sizeof(double));
    } else {
        /* Cyclic tridiagonal: a_i M_{i-1} + b_i M_i + c_i M_{i+1} = r_i (mod n). */
        double* a = malloc(sizeof(double) * n);
        double* b = malloc(sizeof(double) * n);
        double* c = malloc(sizeof(double) * n);
        double* r = malloc(sizeof(double) * n);
        for (size_t i = 0; i < n; i++) {
            double hm = h[(i + n - 1) % n], hi = h[i];
            a[i] = hm; b[i] = 2.0 * (hm + hi); c[i] = hi;
            double dnext = (ys[(i + 1) % n] - ys[i]) / hi;
            double dprev = (ys[i] - ys[(i + n - 1) % n]) / hm;
            r[i] = 6.0 * (dnext - dprev);
        }
        /* Sherman-Morrison: A = T + u v^T, where the corner entries a[0] (M_{n-1})
         * and c[n-1] (M_0) are moved into a rank-1 update. */
        double alpha = a[0], beta = c[n - 1];
        double gamma = -b[0];
        double* bb = malloc(sizeof(double) * n);
        for (size_t i = 0; i < n; i++) bb[i] = b[i];
        bb[0] = b[0] - gamma;
        bb[n - 1] = b[n - 1] - alpha * beta / gamma;
        /* Solve T y = r and T z = u (u = gamma e_0 + beta e_{n-1}). */
        double* u = calloc(n, sizeof(double));
        u[0] = gamma; u[n - 1] = beta;
        double* cp = malloc(sizeof(double) * n);
        double* y1 = malloc(sizeof(double) * n);
        double* z1 = malloc(sizeof(double) * n);
        /* Thomas on tridiag (a sub, bb diag, c super) for rhs r -> y1 and u -> z1. */
        cp[0] = c[0] / bb[0]; y1[0] = r[0] / bb[0]; z1[0] = u[0] / bb[0];
        for (size_t i = 1; i < n; i++) {
            double mden = bb[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / mden;
            y1[i] = (r[i] - a[i] * y1[i - 1]) / mden;
            z1[i] = (u[i] - a[i] * z1[i - 1]) / mden;
        }
        for (size_t i = n - 1; i-- > 0; ) { y1[i] -= cp[i] * y1[i + 1]; z1[i] -= cp[i] * z1[i + 1]; }
        /* x = y1 - ((v . y1)/(1 + v . z1)) z1, with v = e_0 + (alpha/gamma) e_{n-1}. */
        double vy = y1[0] + (alpha / gamma) * y1[n - 1];
        double vz = z1[0] + (alpha / gamma) * z1[n - 1];
        double fac = vy / (1.0 + vz);
        M = malloc(sizeof(double) * n);
        for (size_t i = 0; i < n; i++) M[i] = y1[i] - fac * z1[i];
        free(a); free(b); free(c); free(r); free(bb); free(u); free(cp); free(y1); free(z1);
    }
    free(h);
    return M;
}

/* d-th derivative at p (already reduced into [x_0, x_0 + P)) of the periodic
 * cubic spline with precomputed node second derivatives M. */
static double spline_piece_eval_periodic(const double* xs, const double* ys,
                                         const double* M, size_t n, double P,
                                         double p, size_t d) {
    /* Interval i (0..n-1); interval n-1 is the wrap [x_{n-1}, x_0+P]. */
    size_t i = (p >= xs[n - 1]) ? n - 1 : bracket_interval(xs, n, p);
    double xL = xs[i], xR = (i + 1 < n) ? xs[i + 1] : xs[0] + P;
    double yL = ys[i], yR = (i + 1 < n) ? ys[i + 1] : ys[0];
    double ML = M[i], MR = (i + 1 < n) ? M[i + 1] : M[0];
    double hh = xR - xL, A = xR - p, B = p - xL;
    double cL = yL / hh - ML * hh / 6.0, cR = yR / hh - MR * hh / 6.0;
    switch (d) {
        case 0: return ML * A * A * A / (6 * hh) + MR * B * B * B / (6 * hh) + cL * A + cR * B;
        case 1: return -ML * A * A / (2 * hh) + MR * B * B / (2 * hh) - cL + cR;
        case 2: return ML * A / hh + MR * B / hh;
        case 3: return (MR - ML) / hh;
        default: return 0.0;
    }
}

static double spline_eval_periodic(const double* xs, size_t n, const double* ys,
                                   double P, double p, size_t d) {
    double* M = spline_solve_M_periodic(xs, ys, n, P);
    double val = spline_piece_eval_periodic(xs, ys, M, n, P, p, d);
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
    bool*   periodic;    /* per-dimension periodicity, length m */
    double* period;      /* per-dimension period (grid span), 0 if aperiodic */
    /* Which component `V` currently holds, so that evaluating the same
     * interpolant at another point does not rebuild it -- see fill_values. */
    bool    v_valid;
    int     v_rank;
    int     v_cpath[16];
    /* Per-interval piecewise-polynomial coefficients (1-D scalar, aperiodic), a
     * la scipy's PPoly: (nk[0]-1) segments x ppoly_ncoef monomial coefficients
     * in t = x - xs[i], so a batch evaluation is bracket + Horner per point
     * instead of re-deriving the local polynomial. Built lazily by
     * interp_ppoly_ensure and cached by (method, order). */
    double* ppoly;
    bool    ppoly_valid;
    int     ppoly_method;
    int     ppoly_order;
    int     ppoly_ncoef;
} IFun;

static void ifun_free(IFun* f) {
    if (f->grid) { for (size_t k = 0; k < f->m; k++) free(f->grid[k]); free(f->grid); }
    free(f->nk); free(f->stride); free(f->V);
    free(f->dmin); free(f->dmax); free(f->has_range); free(f->entryAt);
    free(f->periodic); free(f->period);
    free(f->ppoly);
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

/* Ascending order for qsort. NaN cannot reach here: build_grid rejects a table
 * whose coordinates are not all node_to_double-convertible before this runs. */
static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
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
static bool build_grid(Expr* domain, Expr* table, size_t m,
                       const bool* periodic, IFun* out) {
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

    /* Gather then sort-and-unique, rather than inserting one coordinate at a
     * time into a sorted array. grid_insert scans from the front for the
     * insertion point, so a table whose coordinates ARRIVE IN ORDER -- which is
     * every table produced by Table[{x, f[x]}, {x, a, b, dx}] -- walked the
     * whole grid built so far on every one of its npts insertions: O(npts^2),
     * and the dominant cost of using a large interpolant at all. Building the
     * grid for a 160,000-node table took 4.96 s; it is now 21 ms.
     *
     * The result is identical: both produce the ascending distinct coordinates,
     * and grid_index still binary-searches it. */
    for (size_t k = 0; k < m; k++) {
        grid[k] = malloc(sizeof(double) * npts);
        if (!grid[k]) goto fail;
        for (size_t i = 0; i < npts; i++)
            node_to_double(entry_coord(table->data.function.args[i], m, k),
                           &grid[k][i]);
        qsort(grid[k], npts, sizeof(double), cmp_double);
        size_t len = 0;
        for (size_t i = 0; i < npts; i++)
            if (len == 0 || grid[k][i] != grid[k][len - 1]) grid[k][len++] = grid[k][i];
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

    bool*   per  = calloc(m, sizeof(bool));
    double* perd = calloc(m, sizeof(double));
    if (!per || !perd) { free(per); free(perd); goto fail; }
    for (size_t k = 0; k < m; k++) {
        per[k] = periodic && periodic[k];
        perd[k] = per[k] ? (grid[k][nk[k] - 1] - grid[k][0]) : 0.0;
    }

    out->m = m; out->nk = nk; out->grid = grid; out->stride = stride;
    out->total = total; out->dmin = dmin; out->dmax = dmax; out->has_range = hasr;
    out->entryAt = eat; out->periodic = per; out->period = perd;
    return true;

fail:
    if (grid) { for (size_t k = 0; k < m; k++) free(grid[k]); free(grid); }
    free(nk); free(stride); free(dmin); free(dmax); free(hasr); free(eat);
    return false;
}

/* --- parsed-grid cache (single entry) ---------------------------------
 *
 * Without caching, interp_apply rebuilt the entire tensor grid (abscissae,
 * strides, node->entry map, value tensor -- O(n) allocations with an O(n^2)
 * grid_insert term) from scratch on every evaluated point.  Cache the
 * most-recently-built IFun and reuse it when the next call carries a
 * structurally identical (domain, table, periodic).  A single entry exactly
 * serves the common build-once / evaluate-many pattern (e.g. Table[f[x],
 * {x, ...}]).  Keying on full structural identity (expr_eq) keeps it correct:
 * redefining the data yields a different table, misses, and rebuilds.
 * Machine-precision path only; the MPFR kernels are unaffected.  (The other
 * former per-point O(n) cost -- the evaluator re-evaluating the held-free
 * table arg -- is removed separately by the HoldAll attribute in interp_init.) */
static struct {
    bool   valid;
    Expr*  domain;       /* owned deep copy */
    Expr*  table;        /* owned deep copy; f.entryAt points into this */
    size_t m;
    bool   has_periodic; /* whether `periodic` below is meaningful */
    bool*  periodic;     /* owned copy, length m, or NULL */
    IFun   f;            /* built grid, referencing the cached `table` */
} g_grid_cache;

/* Release the cached entry (safe on a zeroed cache: frees NULLs / a zero IFun). */
static void grid_cache_free(void) {
    expr_free(g_grid_cache.domain);
    expr_free(g_grid_cache.table);
    free(g_grid_cache.periodic);
    ifun_free(&g_grid_cache.f);
    memset(&g_grid_cache, 0, sizeof(g_grid_cache));
}

/* Element-wise compare the cached periodicity against (periodic, periodic!=NULL). */
static bool periodic_match(const bool* periodic, size_t m) {
    bool has = (periodic != NULL);
    if (g_grid_cache.has_periodic != has) return false;
    if (!has) return true;
    for (size_t k = 0; k < m; k++)
        if (g_grid_cache.periodic[k] != periodic[k]) return false;
    return true;
}

/*
 * grid_cache_get:
 *   Return the parsed grid for (domain, table, m, periodic), building and
 *   caching it on a miss.  The returned IFun is borrowed (owned by the cache) --
 *   do NOT ifun_free it.  Returns NULL if the grid is malformed (cache stays
 *   empty), matching build_grid's failure contract.
 */
static IFun* grid_cache_get(Expr* domain, Expr* table, size_t m,
                            const bool* periodic) {
    /* Identity before equality. expr_copy is a refcount bump, so the cached
     * `table` is normally the very pointer the caller passes back on the next
     * evaluation of the same InterpolatingFunction -- and expr_eq on a table of
     * n entries is itself O(n), which made the cache HIT cost as much per point
     * as the work it was meant to save. Falls through to the structural test for
     * an equal table that is a different object. */
    if (g_grid_cache.valid && g_grid_cache.m == m
        && periodic_match(periodic, m)
        && ((g_grid_cache.domain == domain && g_grid_cache.table == table)
            || (expr_eq(g_grid_cache.domain, domain)
                && expr_eq(g_grid_cache.table, table))))
        return &g_grid_cache.f;

    grid_cache_free();
    g_grid_cache.domain = expr_copy(domain);
    g_grid_cache.table  = expr_copy(table);
    g_grid_cache.m      = m;
    g_grid_cache.has_periodic = (periodic != NULL);
    if (periodic) {
        g_grid_cache.periodic = malloc(sizeof(bool) * m);
        if (!g_grid_cache.periodic) { grid_cache_free(); return NULL; }
        memcpy(g_grid_cache.periodic, periodic, sizeof(bool) * m);
    }
    /* Build against the cached copies so f.entryAt points into stable memory. */
    if (!build_grid(g_grid_cache.domain, g_grid_cache.table, m,
                    g_grid_cache.periodic, &g_grid_cache.f)) {
        grid_cache_free();
        return NULL;
    }
    g_grid_cache.valid = true;
    return &g_grid_cache.f;
}

static bool value_component(Expr* entry, const int* cpath, int vrank, double* out);

/*
 * Fill the flat value tensor (component `cpath` of each node value).
 *
 * MEMOISED ON THE COMPONENT, because it is O(number of nodes) and was being run
 * once per EVALUATION POINT. `V` is a function of the table and the component
 * alone -- nothing here depends on where the interpolant is being evaluated --
 * yet eval_component_double called it unconditionally, so evaluating a 10^4-node
 * interpolant at 10^4 points walked 10^8 table entries, converting each Expr to
 * a double, and did 10^4 malloc/free pairs of 80 kB. That is what made the
 * per-point cost scale with the TABLE size instead of being constant: 151 us per
 * point at 10^4 nodes and 906 us at 5x10^4.
 *
 * `f` is owned by the grid cache and lives as long as the table it was built
 * from (grid_cache_get rebuilds it whenever either changes), so the memo is
 * valid for exactly as long as `f` is. Array-valued data evaluates one component
 * per point and so still refills per component -- correct, and no worse than
 * before; the scalar case, which is nearly all of them, refills once ever.
 */
static bool fill_values(IFun* f, const int* cpath, int vrank) {
    if (f->v_valid && f->v_rank == vrank && vrank <= 16) {
        bool same = true;
        for (int d = 0; d < vrank; d++)
            if (f->v_cpath[d] != cpath[d]) { same = false; break; }
        if (same) return true;
    }
    free(f->V);
    f->V = malloc(sizeof(double) * f->total);
    if (!f->V) return false;
    f->v_valid = false;
    for (size_t i = 0; i < f->total; i++)
        if (!value_component(f->entryAt[i], cpath, vrank, &f->V[i])) return false;
    if (vrank <= 16) {
        f->v_rank = vrank;
        for (int d = 0; d < vrank; d++) f->v_cpath[d] = cpath[d];
        f->v_valid = true;
    }
    return true;
}

/* --- value-only tensor evaluation (default / spline) ------------------ */

typedef struct {
    const IFun*   f;
    const double* p;
    const int*    ders;
    const int64_t* s;       /* signed window start (may be negative when periodic) */
    const size_t* w;
    int           kernel;   /* 0 = newton, 1 = spline */
} EvalCtx;

static double eval_dim(const EvalCtx* ctx, size_t k, size_t base) {
    const IFun* f = ctx->f;
    size_t wk = ctx->w[k];
    int64_t start = ctx->s[k];
    bool per = f->periodic && f->periodic[k];
    size_t ndist = per ? f->nk[k] - 1 : f->nk[k];
    double P = per ? f->period[k] : 0.0;

    double* xs = malloc(sizeof(double) * wk);
    double* tmp = malloc(sizeof(double) * wk);
    for (size_t t = 0; t < wk; t++) {
        int64_t j = start + (int64_t)t;
        size_t r; double xoff;
        if (per) {
            int64_t nd = (int64_t)ndist;
            int64_t rr = ((j % nd) + nd) % nd;
            xoff = (double)((j - rr) / nd) * P;
            r = (size_t)rr;
        } else { r = (size_t)j; xoff = 0.0; }
        xs[t] = f->grid[k][r] + xoff;
        size_t off = base + r * f->stride[k];
        tmp[t] = (k + 1 == f->m) ? f->V[off] : eval_dim(ctx, k + 1, off);
    }
    double res;
    if (ctx->kernel == 1)
        res = per ? spline_eval_periodic(xs, wk, tmp, P, ctx->p[k], (size_t)ctx->ders[k])
                  : spline_eval(xs, tmp, wk, ctx->p[k], (size_t)ctx->ders[k]);
    else
        res = newton_deriv_eval(xs, tmp, wk, ctx->p[k], (size_t)ctx->ders[k]);
    free(xs); free(tmp);
    return res;
}

/* --- supplied/estimated mixed-derivative tensor (Hermite) ------------- */

/* Navigate the component path `cpath` (length vrank) into a value/derivative
 * tensor `t` (its outermost vrank axes index the value-array components). */
static Expr* descend_cpath(Expr* t, const int* cpath, int vrank) {
    for (int d = 0; d < vrank; d++) {
        if (!interp_is_list(t) || (int)t->data.function.arg_count <= cpath[d]) return NULL;
        t = t->data.function.args[cpath[d]];
    }
    return t;
}

/* The scalar value of component `cpath` at node `entry` (vrank == 0 -> scalar
 * value). */
static bool value_component(Expr* entry, const int* cpath, int vrank, double* out) {
    Expr* t = descend_cpath(entry->data.function.args[1], cpath, vrank);
    return t && node_to_double(t, out);
}

/* D^a f_component at the node held in `entry` (a has total order s = mi_sum(a));
 * read from the supplied s-th derivative tensor (s==0 -> the value).  For
 * array-valued data, `cpath` (length vrank) selects the value-array component;
 * the derivative axes follow the component axes. */
static bool extract_supplied(Expr* entry, size_t m, const int* cpath, int vrank,
                             const int* a, int s, double* out) {
    if ((size_t)(1 + s) >= entry->data.function.arg_count) return false;
    Expr* t = descend_cpath(entry->data.function.args[1 + s], cpath, vrank);
    if (!t) return false;
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
 * one-sided at the ends; cyclic central differences on a periodic dimension),
 * writing into G. */
static void fd_apply(const IFun* g, size_t d, const double* F, double* G) {
    size_t nd = g->nk[d], st = g->stride[d];
    const double* xs = g->grid[d];
    bool per = g->periodic && g->periodic[d];
    if (per && nd >= 3) {
        size_t ndist = nd - 1;          /* node nd-1 ≡ node 0 */
        double P = g->period[d];
        for (size_t n = 0; n < g->total; n++) {
            size_t id = (n / st) % nd;
            size_t base = n - id * st;
            size_t idm = (id == nd - 1) ? 0 : id;       /* duplicate maps to node 0 */
            size_t L = (idm + ndist - 1) % ndist, R = (idm + 1) % ndist;
            double xc = xs[idm];
            double xl = xs[L] - (idm == 0 ? P : 0.0);
            double xr = xs[R] + (idm == ndist - 1 ? P : 0.0);
            double h0 = xc - xl, h1 = xr - xc;
            double d0 = (F[base + idm * st] - F[base + L * st]) / h0;
            double d1 = (F[base + R * st] - F[base + idm * st]) / h1;
            G[n] = (h1 * d0 + h0 * d1) / (h0 + h1);
        }
        return;
    }
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
static double* build_T(const IFun* g, size_t m, int K, int k,
                       const int* cpath, int vrank) {
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
                if (!extract_supplied(g->entryAt[n], m, cpath, vrank, a, s, &v)) { ok = false; break; }
                T[n * dcount + idx] = v;
            }
        } else {
            int rem = K;
            for (size_t d = 0; d < m; d++) { int bd = a[d] < rem ? a[d] : rem; b[d] = bd; rem -= bd; }
            for (size_t n = 0; n < g->total; n++) {
                double v;
                if (!extract_supplied(g->entryAt[n], m, cpath, vrank, b, K, &v)) { ok = false; break; }
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

/* --- value-array shape + result assembly ------------------------------ */

static Expr* make_list(Expr** items, size_t n);   /* defined below */

/* Shape of the value array `e` (from the first node): vshape[0..*vrank-1].
 * A scalar gives rank 0.  Caller provides space for up to 16 axes. */
static void value_shape(Expr* e, size_t* vshape, int* vrank) {
    int r = 0;
    while (interp_is_list(e) && r < 16) {
        vshape[r++] = e->data.function.arg_count;
        e = e->data.function.args[0];
    }
    *vrank = r;
}

/* Assemble a nested-List Expr of shape vshape from the row-major component
 * results `comps`.  Rank 0 yields a bare Real. */
static Expr* assemble_array(const double* comps, const size_t* vshape,
                            int vrank, int level, size_t* idx) {
    if (level == vrank) return expr_new_real(comps[(*idx)++]);
    size_t n = vshape[level];
    Expr** items = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) items[i] = assemble_array(comps, vshape, vrank, level + 1, idx);
    Expr* l = make_list(items, n);
    free(items);
    return l;
}

/* Evaluate a single scalar component (selected by `cpath`) of the interpolant
 * at `p`.  Returns false (no value) on malformed data. */
static bool eval_component_double(IFun* f, const int* cpath, int vrank, size_t m,
                                  const double* p, const int* ders,
                                  const int* orders, int method, int Ksupplied,
                                  double* out) {
    if (method == METHOD_HERMITE || Ksupplied >= 1) {
        int K = Ksupplied, k = K > 1 ? K : 1;
        double* basis = build_basis(k);
        double* T = build_T(f, m, K, k, cpath, vrank);
        if (!T) { free(basis); return false; }
        *out = hermite_tensor_eval(f, T, ipow_sz((size_t)(k + 1), m), m, k,
                                   basis, 2 * k + 2, p, ders);
        free(T); free(basis);
        return true;
    }
    if (!fill_values(f, cpath, vrank)) return false;
    int64_t* s = malloc(sizeof(int64_t) * m);
    size_t* w = malloc(sizeof(size_t) * m);
    int kernel = (method == METHOD_SPLINE) ? 1 : 0;
    for (size_t k = 0; k < m; k++) {
        size_t nkk = f->nk[k];
        bool per = f->periodic && f->periodic[k];
        size_t ndist = per ? nkk - 1 : nkk;            /* distinct cyclic nodes */
        if (kernel == 1) { s[k] = 0; w[k] = ndist; continue; }  /* spline: full grid */
        size_t order;
        if (orders) { order = (size_t)orders[k]; }
        else { order = nkk - 1; if (order > 3) order = 3; }
        if (order > ndist - 1) order = ndist - 1;
        size_t wk = order + 1;
        size_t i = bracket_interval(f->grid[k], nkk, p[k]);
        size_t shift = (wk >= 2) ? (wk / 2 - 1) : 0;
        if (per) {
            s[k] = (int64_t)i - (int64_t)shift;        /* signed; window wraps */
        } else {
            size_t sk = (i < shift) ? 0 : i - shift;
            if (sk > nkk - wk) sk = nkk - wk;
            s[k] = (int64_t)sk;
        }
        w[k] = wk;
    }
    EvalCtx ctx = { f, p, ders, s, w, kernel };
    *out = eval_dim(&ctx, 0, 0);
    free(s); free(w);
    return true;
}

/*
 * interp_eval_double:
 *   Evaluate the object at the (machine-precision) point.  Returns a fresh
 *   EXPR_REAL (or a nested List for array-valued data), or NULL if the point is
 *   non-numeric / the data is malformed.  `Ksupplied` is the number of supplied
 *   derivative tensors per node (0 for value-only data).
 */
static Expr* interp_eval_double(Expr* domain, Expr* table, size_t m,
                                Expr** call_args, const int* ders,
                                const int* orders, int method, int Ksupplied,
                                const bool* periodic) {
    double* p = malloc(sizeof(double) * m);
    for (size_t k = 0; k < m; k++)
        if (!node_to_double(call_args[k], &p[k])) { free(p); return NULL; }

    /* Borrowed: owned by the grid cache (do not ifun_free). */
    IFun* f = grid_cache_get(domain, table, m, periodic);
    if (!f) { free(p); return NULL; }

    /* Reduce periodic coordinates into [x0, x0 + P); warn on aperiodic dims. */
    for (size_t k = 0; k < m; k++) {
        if (f->periodic[k]) {
            double x0 = f->grid[k][0], P = f->period[k];
            double u = x0 + fmod(p[k] - x0, P);
            if (u < x0) u += P;
            p[k] = u;
        } else if (f->has_range[k] && (p[k] < f->dmin[k] || p[k] > f->dmax[k])) {
            fprintf(stderr,
                    "InterpolatingFunction::dmval: Input value %g lies outside "
                    "the range of data in the interpolating function. "
                    "Extrapolation will be used.\n", p[k]);
        }
    }

    size_t vshape[16];
    int vrank;
    value_shape(entry_value(f->entryAt[0]), vshape, &vrank);
    size_t vtotal = 1;
    for (int d = 0; d < vrank; d++) vtotal *= vshape[d];

    double* comps = malloc(sizeof(double) * vtotal);
    int* cpath = malloc(sizeof(int) * (vrank ? (size_t)vrank : 1));
    bool ok = true;
    for (size_t ci = 0; ci < vtotal && ok; ci++) {
        size_t rem = ci;
        for (int d = vrank; d-- > 0; ) { cpath[d] = (int)(rem % vshape[d]); rem /= vshape[d]; }
        ok = eval_component_double(f, cpath, vrank, m, p, ders, orders, method, Ksupplied, &comps[ci]);
    }

    free(p); free(cpath);
    if (!ok) { free(comps); return NULL; }
    size_t idx = 0;
    Expr* result = assemble_array(comps, vshape, vrank, 0, &idx);
    free(comps);
    return result;
}

#ifdef USE_MPFR
/* Implemented in interp_mpfr.c. */
Expr* interp_eval_mpfr(Expr* domain, Expr* table, size_t m,
                       Expr** call_args, const int* ders,
                       const int* orders, int method, int Ksupplied,
                       const bool* periodic, long bits);
#endif

/* --- object inspection (precision + table form) ----------------------- */

/*
 * A memo for the facts interp_apply must know about the TABLE before it can
 * evaluate at any point: how many derivative tensors each node carries, and
 * whether the data is arbitrary-precision. Each answer is a property of the
 * table alone, and each is computed by walking every node -- so asking them per
 * evaluation point made the cost of evaluating an interpolant scale with the
 * size of its own table. On a 5x10^4-node table these three walks
 * (table_Ksupplied, common_scan_inexact, numeric_expr_is_mpfr) were where the
 * profile spent essentially all of its samples.
 *
 * Keyed on the table's IDENTITY. Sound for the same reason the grid cache's
 * pointer fast path is: the memo holds a reference (expr_copy is a refcount
 * bump), so the address cannot be recycled underneath it while cached, and
 * expressions are immutable by convention once shared.
 */
static struct {
    Expr* table;              /* owned reference; identity is the key */
    int   Ksupplied;
#ifdef USE_MPFR
    bool  is_mpfr;
    CommonInexactInfo inexact;
#endif
} g_obj_cache;

static int table_Ksupplied(Expr* table);

static void obj_cache_load(Expr* table) {
    if (g_obj_cache.table == table) return;
    expr_free(g_obj_cache.table);
    memset(&g_obj_cache, 0, sizeof(g_obj_cache));
    g_obj_cache.table     = expr_copy(table);
    g_obj_cache.Ksupplied = table_Ksupplied(table);
#ifdef USE_MPFR
    g_obj_cache.is_mpfr   = numeric_expr_is_mpfr(table);
    g_obj_cache.inexact   = common_scan_inexact(table);
#endif
}

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

/* --- vectorised 1-D application (the scipy cs(array) equivalent) ------- */

/* True when `a` is a real array of query points: a rank-1 real NDArray, or a
 * non-empty List whose every element is a real scalar.  Sets *n to its length. */
static bool interp_arg_is_real_array(const Expr* a, size_t* n) {
    if (is_ndarray(a) && a->data.ndarray.rank == 1 && !ndt_is_complex(a->data.ndarray.dtype)) {
        *n = (size_t)a->data.ndarray.dims[0];
        return *n >= 1;
    }
    if (interp_is_list(a) && a->data.function.arg_count >= 1) {
        size_t k = a->data.function.arg_count;
        for (size_t i = 0; i < k; i++) {
            double t;
            if (!node_to_double(a->data.function.args[i], &t)) return false;
        }
        *n = k;
        return true;
    }
    return false;
}

/* Build (once, cached in `f`) the per-interval monomial coefficients of the 1-D
 * scalar interpolant: for each segment i, the local polynomial in t = x - xs[i]
 * with coefficient a[t] = P_i^(t)(xs[i]) / t!.  For the default method P_i is the
 * sliding-window Newton polynomial for that segment; for "Spline" it is the
 * natural cubic piece.  Cached by (method, order); returns false only on OOM.
 * This turns a batch evaluation into bracket + Horner per point, with the
 * divided-difference / tridiagonal work amortised over the whole grid. */
static bool interp_ppoly_ensure(IFun* f, const double* V, int method, int order) {
    if (f->ppoly_valid && f->ppoly_method == method && f->ppoly_order == order)
        return true;

    size_t nk = f->nk[0];
    if (nk < 2) return false;
    size_t nseg = nk - 1;
    const double* xs = f->grid[0];
    int ncoef = (method == METHOD_SPLINE) ? 4 : (order + 1);
    int deg = ncoef - 1;

    free(f->ppoly);
    f->ppoly = malloc(sizeof(double) * nseg * (size_t)ncoef);
    if (!f->ppoly) { f->ppoly_valid = false; return false; }

    if (method == METHOD_SPLINE) {
        double* M = spline_solve_M(xs, V, nk);
        for (size_t i = 0; i < nseg; i++) {
            double* a = f->ppoly + i * (size_t)ncoef;
            double fact = 1.0;
            for (int t = 0; t <= deg; t++) {
                a[t] = spline_piece_eval(xs, V, M, nk, xs[i], (size_t)t) / fact;
                fact *= (double)(t + 1);
            }
        }
        free(M);
    } else {
        size_t w = (size_t)order + 1;
        size_t shift = (w >= 2) ? (w / 2 - 1) : 0;
        double* c = malloc(sizeof(double) * w);
        double* vs = malloc(sizeof(double) * (size_t)(deg + 1));
        for (size_t i = 0; i < nseg; i++) {
            size_t sk = (i < shift) ? 0 : i - shift;
            if (sk > nk - w) sk = nk - w;
            double* a = f->ppoly + i * (size_t)ncoef;
            double fact = 1.0;
            for (int t = 0; t <= deg; t++) {
                a[t] = newton_deriv_eval_buf(xs + sk, V + sk, w, xs[i], (size_t)t, c, vs) / fact;
                fact *= (double)(t + 1);
            }
        }
        free(c); free(vs);
    }
    f->ppoly_valid = true;
    f->ppoly_method = method;
    f->ppoly_order = order;
    f->ppoly_ncoef = ncoef;
    return true;
}

/* der-th derivative at t of segment i's cached monomial polynomial (Horner). */
static double interp_ppoly_eval(const IFun* f, size_t i, double t, int der) {
    int ncoef = f->ppoly_ncoef, deg = ncoef - 1;
    const double* a = f->ppoly + i * (size_t)ncoef;
    double acc = 0.0;
    for (int k = deg; k >= der; k--)
        acc = acc * t + a[k] * falling(k, der);   /* falling(k,der)=k!/(k-der)! */
    return acc;
}

/*
 * Evaluate a 1-D, scalar-valued, value-only InterpolatingFunction at every point
 * of the real array `arg` in a single C loop over the grid double buffers --- the
 * counterpart of scipy's cs(array).  The spline system is solved ONCE and reused
 * across all points; the default (Newton) method uses the same per-point window
 * as the scalar path.  Returns a packed float64 result (inheriting `arg`'s
 * presentation for an NDArray, else a List of reals), or NULL if the object is
 * not in this simple form (array-valued, Hermite, supplied derivatives, or
 * periodic-with-default-method) so the caller can leave it unevaluated.
 */
static Expr* interp_vector_1d(Expr* domain, Expr* table, const Expr* arg, size_t n,
                              int der0, const int* orders, int method,
                              const bool* periodic, int Ksupplied) {
    if (Ksupplied >= 1 || method == METHOD_HERMITE) return NULL;

    IFun* f = grid_cache_get(domain, table, 1, periodic);
    if (!f) return NULL;

    size_t vshape[16]; int vrank;
    value_shape(entry_value(f->entryAt[0]), vshape, &vrank);
    if (vrank != 0) return NULL;                 /* scalar-valued only */
    if (!fill_values(f, NULL, 0)) return NULL;   /* populate f->V (memoised) */

    bool per = f->periodic && f->periodic[0];
    if (per && method != METHOD_SPLINE) return NULL;  /* periodic Newton: fall back */

    size_t nk = f->nk[0];
    const double* xs = f->grid[0];
    const double* V = f->V;
    double x0 = xs[0], P = per ? f->period[0] : 0.0;
    bool ranged = f->has_range && f->has_range[0];
    double dmin = ranged ? f->dmin[0] : xs[0];
    double dmax = ranged ? f->dmax[0] : xs[nk - 1];

    /* Query points: use the float64 buffer directly when possible (no copy). */
    const double* q;
    double* q_owned = NULL;
    if (is_ndarray(arg) && arg->data.ndarray.dtype == NDT_FLOAT64) {
        q = (const double*)arg->data.ndarray.data;
    } else if (is_ndarray(arg)) {
        const void* nd = arg->data.ndarray.data;
        NDType dt = arg->data.ndarray.dtype;
        q_owned = malloc(sizeof(double) * n);
        for (size_t i = 0; i < n; i++) { double re, im; ndt_get(nd, i, dt, &re, &im); q_owned[i] = re; }
        q = q_owned;
    } else {
        q_owned = malloc(sizeof(double) * n);
        for (size_t i = 0; i < n; i++)
            if (!node_to_double(arg->data.function.args[i], &q_owned[i])) { free(q_owned); return NULL; }
        q = q_owned;
    }

    double* out = malloc(sizeof(double) * n);
    size_t d = (size_t)der0;

    if (per) {
        /* Periodic cubic spline: solve M once, evaluate each reduced point. */
        size_t ndist = nk - 1;
        double* M = spline_solve_M_periodic(xs, V, ndist, P);
        for (size_t i = 0; i < n; i++) {
            double u = x0 + fmod(q[i] - x0, P); if (u < x0) u += P;
            out[i] = spline_piece_eval_periodic(xs, V, M, ndist, P, u, d);
        }
        free(M);
    } else {
        /* Aperiodic (default/Newton or "Spline"): precompute per-interval
         * coefficients ONCE, then bracket + Horner per point. */
        int order = (method == METHOD_SPLINE)
                        ? 3
                        : (orders ? orders[0] : (nk - 1 < 3 ? (int)(nk - 1) : 3));
        if (method != METHOD_SPLINE && order > (int)(nk - 1)) order = (int)(nk - 1);
        if (!interp_ppoly_ensure(f, V, method, order)) { free(q_owned); free(out); return NULL; }
        const double* pp = f->ppoly;
        int ncoef = f->ppoly_ncoef;
        size_t step0 = 1;
        while ((step0 << 1) <= nk - 1) step0 <<= 1;   /* bracket start, once */

        if (der0 == 0 && ncoef == 4) {
            /* Hot path: cubic value.  Inlined branchless bracket + Horner, no
             * per-point function calls -- this is the loop that has to keep pace
             * with scipy's cs(array). */
            for (size_t i = 0; i < n; i++) {
                double p = q[i];
                size_t pos = 0;
                for (size_t step = step0; step; step >>= 1) {
                    size_t next = pos + step;
                    size_t idx = next < nk ? next : nk - 1;
                    size_t adv = (size_t)((next < nk) & (xs[idx] <= p));
                    pos += adv * step;
                }
                if (pos > nk - 2) pos = nk - 2;
                const double* a = pp + pos * 4;
                double t = p - xs[pos];
                out[i] = ((a[3] * t + a[2]) * t + a[1]) * t + a[0];
            }
        } else {
            /* General order / derivative order. */
            for (size_t i = 0; i < n; i++) {
                size_t bi = bracket_interval(xs, nk, q[i]);
                out[i] = interp_ppoly_eval(f, bi, q[i] - xs[bi], (int)d);
            }
        }

        if (ranged) {
            for (size_t i = 0; i < n; i++)
                if (q[i] < dmin || q[i] > dmax) {
                    fprintf(stderr, "InterpolatingFunction::dmval: Input value %g lies "
                            "outside the range of data in the interpolating function. "
                            "Extrapolation will be used.\n", q[i]);
                    break;
                }
        }
    }
    free(q_owned);

    Expr* result;
    if (is_ndarray(arg)) {
        int64_t dims[1] = { (int64_t)n };
        result = expr_new_ndarray_like(arg, 1, dims, out, NDT_FLOAT64);  /* moves out */
    } else {
        Expr** cells = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) cells[i] = expr_new_real(out[i]);
        result = make_list(cells, n);
        free(cells);
        free(out);
    }
    return result;
}

/* --- public entry points ---------------------------------------------- */

Expr* interp_apply(Expr* ifun, Expr** call_args, size_t argc) {
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc < 2 || obj_argc > 6) return NULL;
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
    if (obj_argc >= 5) {
        Expr* me = ifun->data.function.args[4];
        if (me->type == EXPR_STRING && me->data.string) {
            if (strcmp(me->data.string, "Spline") == 0) method = METHOD_SPLINE;
            else if (strcmp(me->data.string, "Hermite") == 0) method = METHOD_HERMITE;
        }
    }

    /* Per-dimension periodicity (6th arg): a List of True/False. */
    bool* periodic = NULL;
    if (obj_argc == 6) {
        Expr* pe = ifun->data.function.args[5];
        if (!interp_is_list(pe) || pe->data.function.arg_count != m) { free(ders); free(orders); return NULL; }
        periodic = calloc(m, sizeof(bool));
        for (size_t k = 0; k < m; k++) {
            Expr* b = pe->data.function.args[k];
            periodic[k] = (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_True);
        }
    }

    obj_cache_load(table);
    int Ksupplied = g_obj_cache.Ksupplied;
    if (Ksupplied < 0) { free(ders); free(orders); free(periodic); return NULL; }

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
                    free(ders); free(orders); free(periodic);
                    return y;
                }
            }
        }
    }

    Expr* result = NULL;

    /* Vectorised application: ifn[{x1, x2, ...}] / ifn[packedArray] evaluates
     * every point in one C loop (the scipy cs(array) equivalent) instead of one
     * evaluator round-trip per point.  Checked BEFORE the per-argument MPFR
     * scan below, which is O(number of points) and would dominate a machine
     * batch: the machine fast path fires when the DATA is not arbitrary
     * precision (a packed float64 argument never is), so only the object's own
     * precision is consulted here, not the whole point array. */
    {
        size_t nvec = 0;
        bool data_mpfr = false;
#ifdef USE_MPFR
        data_mpfr = g_obj_cache.is_mpfr || numeric_expr_is_mpfr(domain);
#endif
        if (!data_mpfr && m == 1 && argc == 1 &&
            interp_arg_is_real_array(call_args[0], &nvec)) {
            result = interp_vector_1d(domain, table, call_args[0], nvec,
                                      ders[0], orders, method, periodic, Ksupplied);
            if (result) { free(ders); free(orders); free(periodic); return result; }
        }
    }

#ifdef USE_MPFR
    /* Route to the MPFR kernels when the data or the argument carry
     * arbitrary precision. */
    bool want_mpfr = false;
    long bits = 53;
    {
        /* `di` and the table's mpfr-ness come from the object memo: both walk
         * every node, and neither depends on the point being evaluated. The
         * domain is m intervals and the arguments are scalars, so those stay
         * inline. */
        CommonInexactInfo di = g_obj_cache.inexact;
        CommonInexactInfo ai = common_scan_inexact(domain);
        long mn = 0;
        bool any = false;
        if (di.has_inexact) { any = true; mn = di.min_bits; }
        if (ai.has_inexact) { any = true; mn = (!any || ai.min_bits < mn) ? ai.min_bits : mn; mn = (di.has_inexact && di.min_bits < mn) ? di.min_bits : mn; }
        for (size_t k = 0; k < m; k++) {
            CommonInexactInfo ci = common_scan_inexact(call_args[k]);
            if (ci.has_inexact) { if (!any || ci.min_bits < mn) mn = ci.min_bits; any = true; }
        }
        if (g_obj_cache.is_mpfr || numeric_expr_is_mpfr(domain)) want_mpfr = true;
        for (size_t k = 0; k < m; k++) if (numeric_expr_is_mpfr(call_args[k])) want_mpfr = true;
        if (want_mpfr) bits = (mn > 53) ? mn : 53;
    }
    if (want_mpfr)
        result = interp_eval_mpfr(domain, table, m, call_args, ders, orders, method, Ksupplied, periodic, bits);
    else
#endif
        result = interp_eval_double(domain, table, m, call_args, ders, orders,
                                    method, Ksupplied, periodic);

    free(ders); free(orders); free(periodic);
    return result;
}

Expr* interp_make_derivative(Expr* deriv_head, Expr* ifun) {
    if (ifun->type != EXPR_FUNCTION) return NULL;
    size_t obj_argc = ifun->data.function.arg_count;
    if (obj_argc < 2 || obj_argc > 6) return NULL;
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
    Expr* ders_list = expr_new_function(expr_new_symbol(SYM_List), der_args, m);
    free(der_args);
    free(ders);

    size_t out_argc = (obj_argc >= 4) ? obj_argc : 3;
    Expr** args = malloc(sizeof(Expr*) * out_argc);
    args[0] = expr_copy(domain);
    args[1] = expr_copy(table);
    args[2] = ders_list;
    for (size_t k = 3; k < out_argc; k++) args[k] = expr_copy(ifun->data.function.args[k]);
    Expr* result = expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), args, out_argc);
    free(args);
    return result;
}

/* --- Interpolation builtin -------------------------------------------- */

/* List[a, b], taking ownership of a and b. */
static Expr* make_pair(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_List), args, 2);
}

/* List from an owned array of `n` exprs (array borrowed; caller frees). */
static Expr* make_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

static void free_exprs(Expr** xs, size_t n) {
    for (size_t i = 0; i < n; i++) expr_free(xs[i]);
    free(xs);
}

static Expr* builtin_interpolation_impl(Expr* res) {
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
    Expr* periodic_spec = NULL;   /* PeriodicInterpolation rhs: True/False/List */
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_FUNCTION && a->data.function.arg_count == 2
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = a->data.function.args[0];
            Expr* rhs = a->data.function.args[1];
            if (lhs->type != EXPR_SYMBOL) return NULL;
            if (lhs->data.symbol.name == SYM_InterpolationOrder) {
                if (!node_to_order(rhs, &explicit_order)) return NULL;
            } else if (lhs->data.symbol.name == SYM_Method) {
                if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) { /* default */ }
                else if (rhs->type == EXPR_STRING && rhs->data.string
                         && strcmp(rhs->data.string, "Spline") == 0) method = METHOD_SPLINE;
                else if (rhs->type == EXPR_STRING && rhs->data.string
                         && strcmp(rhs->data.string, "Hermite") == 0) method = METHOD_HERMITE;
                else return NULL;
            } else if (lhs->data.symbol.name == SYM_PeriodicInterpolation) {
                periodic_spec = rhs;   /* validated against m below */
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

    /* Resolve PeriodicInterpolation against the dimensionality. */
    bool* periodic = NULL;
    bool any_periodic = false;
    if (periodic_spec) {
        bool all_true = (periodic_spec->type == EXPR_SYMBOL && periodic_spec->data.symbol.name == SYM_True);
        bool all_false = (periodic_spec->type == EXPR_SYMBOL && periodic_spec->data.symbol.name == SYM_False);
        if (!all_false) {
            periodic = calloc(m, sizeof(bool));
            if (all_true) { for (size_t k = 0; k < m; k++) periodic[k] = true; any_periodic = true; }
            else if (interp_is_list(periodic_spec) && periodic_spec->data.function.arg_count == m) {
                for (size_t k = 0; k < m; k++) {
                    Expr* b = periodic_spec->data.function.args[k];
                    periodic[k] = (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_True);
                    if (periodic[k]) any_periodic = true;
                }
            } else { free(periodic); return NULL; }   /* malformed periodic spec */
            if (!any_periodic) { free(periodic); periodic = NULL; }
        }
    }

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
            if (!is_real_or_array(e->data.function.args[1])) { ok = false; break; }  /* value */
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
        free(dmin); free(dmax); free(dminE); free(dmaxE); free(periodic);
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
    bool need_full = (explicit_order >= 0) || (method != METHOD_DEFAULT) || any_periodic;
    if (!need_full) {
        Expr* oargs[2] = { domain, table };
        object = expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), oargs, 2);
    } else {
        int ov = (explicit_order >= 0) ? explicit_order : 3;
        Expr** zeros = malloc(sizeof(Expr*) * m);
        Expr** ords  = malloc(sizeof(Expr*) * m);
        for (size_t k = 0; k < m; k++) { zeros[k] = expr_new_integer(0); ords[k] = expr_new_integer(ov); }
        Expr* ders_list  = make_list(zeros, m);
        Expr* order_list = make_list(ords, m);
        free(zeros); free(ords);
        if (!any_periodic && method == METHOD_DEFAULT) {
            Expr* oargs[4] = { domain, table, ders_list, order_list };
            object = expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), oargs, 4);
        } else {
            Expr* mslot = (method == METHOD_DEFAULT)
                ? expr_new_symbol(SYM_Automatic)
                : expr_new_string(method == METHOD_SPLINE ? "Spline" : "Hermite");
            if (!any_periodic) {
                Expr* oargs[5] = { domain, table, ders_list, order_list, mslot };
                object = expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), oargs, 5);
            } else {
                Expr** pf = malloc(sizeof(Expr*) * m);
                for (size_t k = 0; k < m; k++) pf[k] = expr_new_symbol(periodic[k] ? "True" : "False");
                Expr* per_list = make_list(pf, m);
                free(pf);
                Expr* oargs[6] = { domain, table, ders_list, order_list, mslot, per_list };
                object = expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), oargs, 6);
            }
        }
    }
    free(periodic);

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

/* Accept a packed / visible NDArray data table (a rank-2 {{x,v},...} matrix or a
 * rank-1 value vector) by delisting it once, then process it as a List --- so
 * Interpolation[NDArray[...]] / a packed data table works instead of returning
 * unevaluated, and (with Interpolation packed-aware) the gate does not
 * materialise it twice.  The object stores copies, so the temporary List is
 * freed here. */
static Expr* builtin_interpolation(Expr* res) {
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count >= 1 &&
        is_ndarray(res->data.function.args[0])) {
        Expr* orig = res->data.function.args[0];
        Expr* lst = ndarray_to_nested_list(orig);
        res->data.function.args[0] = lst;         /* borrow the delisted table */
        Expr* r = builtin_interpolation_impl(res);
        res->data.function.args[0] = orig;        /* restore; evaluator owns res */
        expr_free(lst);
        return r;
    }
    return builtin_interpolation_impl(res);
}

/* ===================================================================== *
 *  InterpolatingPolynomial
 *
 *  Exact single-polynomial interpolant (unlike InterpolatingFunction, which
 *  is a piecewise-spline object).  Two engines share one condition model:
 *
 *    Engine A  (1-D, all values present) -- Newton confluent divided
 *              differences -> nested Newton-Horner form.  Exact/symbolic via
 *              the CAS arithmetic heads; a double fast path for inexact data.
 *
 *    Engine B  (multivariate, OR 1-D with Automatic) -- minimal-total-degree
 *              solve over a graded monomial basis by exact Expr Gauss-Jordan
 *              (matrix numeric, RHS may be symbolic).  Success gate is
 *              CONSISTENCY (not full column rank); the degree loop is capped
 *              at d_sq = min d with C(d+m,m) >= N so noipf/poised fire at the
 *              Mathematica-correct degree.
 *
 *    Modulus->n -- the same graded solve done in Z/nZ (ip_solve_modular).
 * ===================================================================== */

/* Small consuming Expr-arithmetic helpers: each ADOPTS its arguments (they are
 * folded into the built node and freed with it) and returns a fresh evaluated
 * result -- the eval_and_free contract, as in central_moment.c. */
static Expr* ip_add(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                           (Expr*[]){ a, b }, 2));
}
static Expr* ip_mul(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                           (Expr*[]){ a, b }, 2));
}
static Expr* ip_sub(Expr* a, Expr* b) {  /* a - b */
    Expr* nb = expr_new_function(expr_new_symbol(SYM_Times),
                                 (Expr*[]){ expr_new_integer(-1), b }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                           (Expr*[]){ a, nb }, 2));
}
static Expr* ip_div(Expr* a, Expr* b) {  /* a / b */
    Expr* ib = expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ b, expr_new_integer(-1) }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                           (Expr*[]){ a, ib }, 2));
}
static Expr* ip_pow_i(Expr* base, int e) {  /* base^e, e >= 1 */
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                                           (Expr*[]){ base, expr_new_integer(e) }, 2));
}

/* True when e is the symbol Automatic. */
static bool ip_is_automatic(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Automatic;
}

/* Integer falling factorial e*(e-1)*...*(e-a+1): 0 if a>e, 1 if a==0. */
static int64_t ip_falling_i(int e, int a) {
    if (a > e) return 0;
    int64_t r = 1;
    for (int i = 0; i < a; i++) r *= (int64_t)(e - i);
    return r;
}
static int64_t ip_factorial_i(int j) {
    int64_t r = 1;
    for (int i = 2; i <= j; i++) r *= (int64_t)i;
    return r;
}

/* ---- condition model ---- */
/* Each known interpolation condition is D^alpha P(point) = value. Automatic
 * conditions are dropped at collection time, so `n` counts KNOWN conditions.
 * Points/values are OWNED copies; IPData owns and frees everything. */
typedef struct {
    size_t m;               /* number of variables */
    size_t n, cap;          /* known conditions */
    Expr**  pE;             /* n*m owned coord Exprs   (pE[i*m + k]) */
    double* pD;             /* n*m coord doubles       (pD[i*m + k]) */
    int*    alpha;          /* n*m multi-indices       (alpha[i*m + k]) */
    Expr**  val;            /* n owned value Exprs */
    bool    any_automatic;  /* any Automatic seen (routes 1-D to Engine B) */
    bool    coords_numeric; /* every coord node_to_double-convertible */
} IPData;

static void ipdata_free(IPData* D) {
    if (D->pE)  { for (size_t t = 0; t < D->n * D->m; t++) expr_free(D->pE[t]); free(D->pE); }
    if (D->val) { for (size_t t = 0; t < D->n; t++) expr_free(D->val[t]); free(D->val); }
    free(D->pD); free(D->alpha);
    memset(D, 0, sizeof(*D));
}

/* Append one condition. coordE[k] borrowed (copied here); value borrowed
 * (copied). alpha copied. */
static bool ip_push(IPData* D, Expr** coordE, const double* coordD,
                    const int* alpha, Expr* value) {
    if (D->n == D->cap) {
        size_t nc = D->cap ? D->cap * 2 : 8;
        Expr**  npE = realloc(D->pE,    nc * D->m * sizeof(Expr*));
        double* npD = realloc(D->pD,    nc * D->m * sizeof(double));
        int*    na  = realloc(D->alpha, nc * D->m * sizeof(int));
        Expr**  nv  = realloc(D->val,   nc * sizeof(Expr*));
        if (!npE || !npD || !na || !nv) { free(npE); free(npD); free(na); free(nv); return false; }
        D->pE = npE; D->pD = npD; D->alpha = na; D->val = nv; D->cap = nc;
    }
    size_t i = D->n;
    for (size_t k = 0; k < D->m; k++) {
        D->pE[i * D->m + k]    = expr_copy(coordE[k]);
        D->pD[i * D->m + k]    = coordD[k];
        D->alpha[i * D->m + k] = alpha[k];
    }
    D->val[i] = expr_copy(value);
    D->n++;
    return true;
}

/* Navigate a rank-`len` derivative tensor by the index sequence seq[0..len-1];
 * returns the borrowed leaf, or NULL if the shape does not match. */
static Expr* ip_tensor_at(Expr* t, const int* seq, int len) {
    for (int d = 0; d < len; d++) {
        if (!interp_is_list(t)) return NULL;
        if (seq[d] < 0 || (size_t)seq[d] >= t->data.function.arg_count) return NULL;
        t = t->data.function.args[seq[d]];
    }
    return t;
}

/* Enumerate weak compositions of `g` into `m` parts (first part slowest),
 * invoking cb(comp) for each. Used both for the order-g derivative tensor
 * multi-indices and for the graded monomial basis. */
typedef void (*ip_comp_cb)(const int* comp, void* ctx);
static void ip_comps_rec(int g, size_t m, size_t pos, int* cur, ip_comp_cb cb, void* ctx) {
    if (pos == m - 1) { cur[pos] = g; cb(cur, ctx); return; }
    for (int v = g; v >= 0; v--) { cur[pos] = v; ip_comps_rec(g - v, m, pos + 1, cur, cb, ctx); }
}

/* Collection context for one entry's tensor of order k. */
typedef struct { IPData* D; Expr** coordE; double* coordD; Expr* tensor; size_t m; bool ok; } IPTensorCtx;
static void ip_tensor_cb(const int* alpha, void* vctx) {
    IPTensorCtx* c = (IPTensorCtx*)vctx;
    if (!c->ok) return;
    /* sorted index sequence for alpha: axis 0 alpha_0 times, axis 1 ... */
    int seq[64]; int len = 0;
    for (size_t k = 0; k < c->m; k++) for (int r = 0; r < alpha[k]; r++) seq[len++] = (int)k;
    Expr* el = ip_tensor_at(c->tensor, seq, len);
    if (!el) return;                 /* malformed tensor: skip this component */
    if (ip_is_automatic(el)) { c->D->any_automatic = true; return; }
    if (!ip_push(c->D, c->coordE, c->coordD, alpha, el)) c->ok = false;
}

/* Parse `data`/`vars` into the condition list. Returns false on malformed data
 * (caller leaves the call unevaluated). */
static bool ip_collect(Expr* data, size_t m, IPData* D) {
    memset(D, 0, sizeof(*D));
    D->m = m;
    D->coords_numeric = true;
    if (!interp_is_list(data) || data->data.function.arg_count < 1) return false;
    size_t npts = data->data.function.arg_count;

    if (m == 1) {
        /* value-form if any top-level element is a scalar; else pairs-form. */
        bool all_lists = true;
        for (size_t i = 0; i < npts; i++)
            if (!interp_is_list(data->data.function.args[i])) { all_lists = false; break; }

        if (!all_lists) {
            /* value form: point = i+1; scalar -> value; list -> {f, df, ...}. */
            for (size_t i = 0; i < npts; i++) {
                Expr* e = data->data.function.args[i];
                double xd = (double)(i + 1);
                Expr* xE = expr_new_integer((int64_t)(i + 1));
                Expr* coordE[1] = { xE }; double coordD[1] = { xd };
                if (interp_is_list(e)) {
                    size_t L = e->data.function.arg_count;
                    if (L < 1) { expr_free(xE); return false; }
                    for (size_t k = 0; k < L; k++) {
                        Expr* v = e->data.function.args[k];
                        if (ip_is_automatic(v)) { D->any_automatic = true; continue; }
                        int a[1] = { (int)k };
                        if (!ip_push(D, coordE, coordD, a, v)) { expr_free(xE); return false; }
                    }
                } else {
                    int a[1] = { 0 };
                    if (!ip_push(D, coordE, coordD, a, e)) { expr_free(xE); return false; }
                }
                expr_free(xE);
            }
        } else {
            /* pairs form: entry = {x, f, df, ...}. x may be symbolic. */
            for (size_t i = 0; i < npts; i++) {
                Expr* e = data->data.function.args[i];
                if (e->data.function.arg_count < 2) return false;
                Expr* xE = e->data.function.args[0];
                double xd = 0.0;
                if (!node_to_double(xE, &xd)) D->coords_numeric = false;
                Expr* coordE[1] = { xE }; double coordD[1] = { xd };
                size_t L = e->data.function.arg_count;
                for (size_t k = 0; k + 1 < L; k++) {
                    Expr* v = e->data.function.args[1 + k];
                    if (ip_is_automatic(v)) { D->any_automatic = true; continue; }
                    int a[1] = { (int)k };
                    if (!ip_push(D, coordE, coordD, a, v)) return false;
                }
            }
        }
        return true;
    }

    /* m >= 2: entry = {{x1,...,xm}, f, t1, t2, ...} with tk the order-k tensor. */
    if (m > 64) return false;
    for (size_t i = 0; i < npts; i++) {
        Expr* e = data->data.function.args[i];
        if (!interp_is_list(e) || e->data.function.arg_count < 2) return false;
        Expr* coord = e->data.function.args[0];
        if (!interp_is_list(coord) || coord->data.function.arg_count != m) return false;
        Expr* coordE[64]; double coordD[64];
        for (size_t k = 0; k < m; k++) {
            coordE[k] = coord->data.function.args[k];
            if (!node_to_double(coordE[k], &coordD[k])) { D->coords_numeric = false; coordD[k] = 0.0; }
        }
        size_t L = e->data.function.arg_count;
        for (size_t k = 0; k + 1 < L; k++) {          /* k = derivative order */
            Expr* tk = e->data.function.args[1 + k];
            if (k == 0) {
                if (ip_is_automatic(tk)) { D->any_automatic = true; continue; }
                int a0[64]; for (size_t q = 0; q < m; q++) a0[q] = 0;
                if (!ip_push(D, coordE, coordD, a0, tk)) return false;
            } else {
                IPTensorCtx c = { D, coordE, coordD, tk, m, true };
                int cur[64];
                ip_comps_rec((int)k, m, 0, cur, ip_tensor_cb, &c);
                if (!c.ok) return false;
            }
        }
    }
    return true;
}

/* ---- Engine A: 1-D Newton (confluent) divided differences ---- */

/* Build the interpolant from Newton coefficients c[0..M-1] and nodes z[0..M-1]
 * (both borrowed): P = c0 + (x-z0)(c1 + (x-z1)(... )). Consumes nothing. */
static Expr* ip_newton_horner(Expr** c, Expr** z, size_t M, Expr* xvar) {
    Expr* P = expr_copy(c[M - 1]);
    for (size_t ii = M - 1; ii-- > 0; ) {
        Expr* xm = ip_sub(expr_copy(xvar), expr_copy(z[ii]));
        P = ip_add(expr_copy(c[ii]), ip_mul(xm, P));
    }
    return P;
}

/* Engine A. Returns the interpolating polynomial, or NULL only on a degenerate
 * table (which routes to the general engine). 1-D, no Automatic. */
static Expr* ip_newton_1d(IPData* D, Expr* xvar) {
    size_t N = D->n;
    if (N == 0) return NULL;

    /* Group conditions into distinct nodes (first-appearance order); collect
     * per-node derivative values indexed by order alpha[0]. */
    Expr** nodeX  = malloc(N * sizeof(Expr*));      /* distinct node coord (borrowed) */
    int*   nodeMx = calloc(N, sizeof(int));         /* max order + 1 at node */
    Expr** deriv  = calloc(N * N, sizeof(Expr*));   /* deriv[node*N + order] borrowed */
    if (!nodeX || !nodeMx || !deriv) { free(nodeX); free(nodeMx); free(deriv); return NULL; }
    size_t nnodes = 0;
    bool numeric_nodes = true, any_confluent = false;
    for (size_t i = 0; i < N; i++) {
        Expr* x = D->pE[i * 1 + 0];
        int   ord = D->alpha[i * 1 + 0];
        size_t ni = nnodes;
        for (size_t t = 0; t < nnodes; t++) if (expr_eq(nodeX[t], x)) { ni = t; break; }
        if (ni == nnodes) { nodeX[nnodes] = x; nnodes++; }
        if (ord + 1 > nodeMx[ni]) nodeMx[ni] = ord + 1;
        if (ord > 0) any_confluent = true;
        deriv[ni * N + ord] = D->val[i];
        double t; if (!node_to_double(x, &t)) numeric_nodes = false;
    }

    /* Double fast path: non-confluent, numeric nodes, at least one inexact
     * value -> compute divided differences in machine doubles. */
    bool any_inexact = false, all_real = true;
    for (size_t i = 0; i < N; i++) {
        Expr* v = D->val[i];
        if (v->type == EXPR_REAL) any_inexact = true;
        double tmp; if (!node_to_double(v, &tmp)) all_real = false;
    }
    if (!any_confluent && numeric_nodes && all_real && any_inexact && nnodes == N) {
        double* xs = malloc(N * sizeof(double));
        double* ys = malloc(N * sizeof(double));
        if (xs && ys) {
            for (size_t i = 0; i < N; i++) { node_to_double(nodeX[i], &xs[i]); node_to_double(deriv[i * N + 0], &ys[i]); }
            for (size_t j = 1; j < N; j++)
                for (size_t i = N; i-- > j; )
                    ys[i] = (ys[i] - ys[i - 1]) / (xs[i] - xs[i - j]);
            /* Drop trailing zero Newton coefficients: a zero top coefficient
             * means the data lie on a lower-degree polynomial, so the minimal
             * (Mathematica-matching) degree is used -- e.g. {1.,4.,9.,16.} is
             * degree 2, not a cubic with a dangling 0. (x-z_i) term. */
            size_t Meff = N;
            while (Meff > 1 && ys[Meff - 1] == 0.0) Meff--;
            Expr** cE = malloc(Meff * sizeof(Expr*));
            Expr** zE = malloc(Meff * sizeof(Expr*));
            for (size_t i = 0; i < Meff; i++) { cE[i] = expr_new_real(ys[i]); zE[i] = expr_new_real(xs[i]); }
            Expr* P = ip_newton_horner(cE, zE, Meff, xvar);
            for (size_t i = 0; i < Meff; i++) { expr_free(cE[i]); expr_free(zE[i]); }
            free(cE); free(zE); free(xs); free(ys);
            free(nodeX); free(nodeMx); free(deriv);
            return P;
        }
        free(xs); free(ys);
    }

    /* Exact / symbolic confluent path. Expand nodes with multiplicity into z[]. */
    size_t M = 0;
    for (size_t t = 0; t < nnodes; t++) M += (size_t)nodeMx[t];
    Expr** z    = malloc(M * sizeof(Expr*));   /* borrowed node Exprs */
    int*   znode = malloc(M * sizeof(int));    /* distinct-node index of each row */
    Expr** Q    = calloc(M * M, sizeof(Expr*));/* owned divided-difference table */
    if (!z || !znode || !Q) { free(z); free(znode); free(Q); free(nodeX); free(nodeMx); free(deriv); return NULL; }
    size_t row = 0;
    for (size_t t = 0; t < nnodes; t++)
        for (int r = 0; r < nodeMx[t]; r++) { z[row] = nodeX[t]; znode[row] = (int)t; row++; }

    for (size_t i = 0; i < M; i++) Q[i * M + 0] = expr_copy(deriv[znode[i] * N + 0]);
    for (size_t j = 1; j < M; j++) {
        for (size_t i = j; i < M; i++) {
            if (znode[i] == znode[i - j]) {                /* coincident run of length j+1 */
                Expr* dv = expr_copy(deriv[znode[i] * N + (int)j]);
                Q[i * M + j] = ip_div(dv, expr_new_integer(ip_factorial_i((int)j)));
            } else {
                Expr* num = ip_sub(expr_copy(Q[i * M + (j - 1)]), expr_copy(Q[(i - 1) * M + (j - 1)]));
                Expr* den = ip_sub(expr_copy(z[i]), expr_copy(z[i - j]));
                Q[i * M + j] = ip_div(num, den);
            }
        }
    }
    Expr** c = malloc(M * sizeof(Expr*));
    for (size_t k = 0; k < M; k++) c[k] = Q[k * M + k];   /* borrow the diagonal */
    Expr* P = ip_newton_horner(c, z, M, xvar);
    free(c);
    for (size_t t = 0; t < M * M; t++) expr_free(Q[t]);
    free(Q); free(z); free(znode);
    free(nodeX); free(nodeMx); free(deriv);
    return P;
}

/* ---- Engine B: minimal-total-degree graded solve ---- */

/* C(a, b) as int64 (small arguments). */
static int64_t ip_binom(int a, int b) {
    if (b < 0 || b > a) return 0;
    int64_t r = 1;
    for (int i = 0; i < b; i++) { r = r * (a - i) / (i + 1); }
    return r;
}

/* Monomial collection (graded). */
typedef struct { int* E; size_t m; size_t cnt; } IPMonos;
static void ip_mono_cb(const int* comp, void* vctx) {
    IPMonos* M = (IPMonos*)vctx;
    for (size_t k = 0; k < M->m; k++) M->E[M->cnt * M->m + k] = comp[k];
    M->cnt++;
}
/* All exponent vectors of total degree <= d, graded. Returns malloc'd E
 * (count*m) and sets *count. */
static int* ip_monomials(int d, size_t m, size_t* count) {
    int64_t total = 0;
    for (int g = 0; g <= d; g++) total += ip_binom(g + (int)m - 1, (int)m - 1);
    int* E = malloc((size_t)total * m * sizeof(int));
    if (!E) { *count = 0; return NULL; }
    IPMonos M = { E, m, 0 };
    int* cur = malloc(m * sizeof(int));
    for (int g = 0; g <= d; g++) ip_comps_rec(g, m, 0, cur, ip_mono_cb, &M);
    free(cur);
    *count = M.cnt;
    return E;
}

/* A[i][j] = D^{alpha_i}(x^{E_j})(p_i), as a fresh numeric Expr. 0^0 == 1. */
static Expr* ip_entry(IPData* D, size_t i, const int* Ej) {
    int64_t coeff = 1;
    for (size_t k = 0; k < D->m; k++) {
        coeff *= ip_falling_i(Ej[k], D->alpha[i * D->m + k]);
        if (coeff == 0) return expr_new_integer(0);
    }
    Expr* e = expr_new_integer(coeff);
    for (size_t k = 0; k < D->m; k++) {
        int exp = Ej[k] - D->alpha[i * D->m + k];
        if (exp > 0) e = ip_mul(e, ip_pow_i(expr_copy(D->pE[i * D->m + k]), exp));
    }
    return e;
}

/* Recursive Horner-by-variable over the terms `idxs` (indices into E/C),
 * factoring variable `k`. */
static Expr* ip_horner(const int* E, Expr** C, const int* idxs, size_t nid,
                       size_t m, size_t k, Expr** vars) {
    if (k == m) {
        Expr* s = expr_new_integer(0);
        for (size_t t = 0; t < nid; t++) s = ip_add(s, expr_copy(C[idxs[t]]));
        return s;
    }
    int pmax = 0;
    for (size_t t = 0; t < nid; t++) if (E[idxs[t] * m + k] > pmax) pmax = E[idxs[t] * m + k];
    Expr* result = NULL;
    int* grp = malloc(nid * sizeof(int));
    for (int p = pmax; p >= 0; p--) {
        size_t ng = 0;
        for (size_t t = 0; t < nid; t++) if (E[idxs[t] * m + k] == p) grp[ng++] = idxs[t];
        Expr* sub = ng ? ip_horner(E, C, grp, ng, m, k + 1, vars) : expr_new_integer(0);
        result = result ? ip_add(ip_mul(result, expr_copy(vars[k])), sub) : sub;
    }
    free(grp);
    return result;
}

/* Build the polynomial from coefficients C[0..L-1] over monomials E. */
static Expr* ip_build_poly(const int* E, Expr** C, size_t L, size_t m, Expr** vars) {
    int* idxs = malloc(L * sizeof(int));
    size_t nid = 0;
    for (size_t j = 0; j < L; j++)
        if (zero_test_decide(C[j]) != ZERO_TEST_TRUE) idxs[nid++] = (int)j;
    Expr* P = nid ? ip_horner(E, C, idxs, nid, m, 0, vars) : expr_new_integer(0);
    free(idxs);
    return P;
}

/* Solve the interpolation system at total degree d. On success returns the
 * polynomial; on an INCONSISTENT system returns NULL with *consistent=false so
 * the caller raises d. */
static Expr* ip_solve_degree(IPData* D, int d, Expr** vars, bool* consistent) {
    size_t m = D->m, N = D->n, L = 0;
    int* E = ip_monomials(d, m, &L);
    if (!E) { *consistent = false; return NULL; }

    Expr** A = malloc(N * L * sizeof(Expr*));
    Expr** f = malloc(N * sizeof(Expr*));
    int64_t* pivrow = malloc(L * sizeof(int64_t));
    for (size_t j = 0; j < L; j++) pivrow[j] = -1;
    for (size_t i = 0; i < N; i++) {
        f[i] = expr_copy(D->val[i]);
        for (size_t j = 0; j < L; j++) A[i * L + j] = ip_entry(D, i, &E[j * m]);
    }

    /* Gauss-Jordan; the matrix A is numeric so pivots are exact. */
    size_t prow = 0;
    for (size_t j = 0; j < L && prow < N; j++) {
        size_t piv = N;
        for (size_t r = prow; r < N; r++)
            if (zero_test_decide(A[r * L + j]) != ZERO_TEST_TRUE) { piv = r; break; }
        if (piv == N) continue;
        if (piv != prow) {
            for (size_t jj = 0; jj < L; jj++) { Expr* t = A[piv * L + jj]; A[piv * L + jj] = A[prow * L + jj]; A[prow * L + jj] = t; }
            Expr* tf = f[piv]; f[piv] = f[prow]; f[prow] = tf;
        }
        Expr* pivc = expr_copy(A[prow * L + j]);
        for (size_t jj = j; jj < L; jj++) A[prow * L + jj] = ip_div(A[prow * L + jj], expr_copy(pivc));
        f[prow] = ip_div(f[prow], expr_copy(pivc));
        expr_free(pivc);
        for (size_t r = 0; r < N; r++) {
            if (r == prow) continue;
            if (zero_test_decide(A[r * L + j]) == ZERO_TEST_TRUE) continue;
            Expr* facc = expr_copy(A[r * L + j]);
            for (size_t jj = j; jj < L; jj++)
                A[r * L + jj] = ip_sub(A[r * L + jj], ip_mul(expr_copy(facc), expr_copy(A[prow * L + jj])));
            f[r] = ip_sub(f[r], ip_mul(expr_copy(facc), expr_copy(f[prow])));
            expr_free(facc);
        }
        pivrow[j] = (int64_t)prow;
        prow++;
    }

    bool ok = true;
    for (size_t r = prow; r < N; r++)
        if (zero_test_decide(f[r]) != ZERO_TEST_TRUE) { ok = false; break; }

    Expr* result = NULL;
    if (ok) {
        Expr** C = malloc(L * sizeof(Expr*));
        for (size_t j = 0; j < L; j++)
            C[j] = (pivrow[j] >= 0) ? expr_copy(f[pivrow[j]]) : expr_new_integer(0);
        result = ip_build_poly(E, C, L, m, vars);
        for (size_t j = 0; j < L; j++) expr_free(C[j]);
        free(C);
    }

    for (size_t t = 0; t < N * L; t++) expr_free(A[t]);
    for (size_t i = 0; i < N; i++) expr_free(f[i]);
    free(A); free(f); free(pivrow); free(E);
    *consistent = ok;
    return result;
}

/* Emit the ::poised / ::noipf diagnostic and return NULL (unevaluated). */
static Expr* ip_fail(IPData* D, int d_sq) {
    if (D->m == 1) {
        fprintf(stderr, "InterpolatingPolynomial::noipf: Unable to find an "
                        "interpolating polynomial of total degree %d.\n", d_sq);
        return NULL;
    }
    /* Distinct abscissae in input order, for the message. */
    Expr** pts = malloc(D->n * sizeof(Expr*));
    size_t np = 0;
    for (size_t i = 0; i < D->n; i++) {
        Expr* coord;
        if (D->m == 1) coord = expr_copy(D->pE[i]);
        else {
            Expr** cs = malloc(D->m * sizeof(Expr*));
            for (size_t k = 0; k < D->m; k++) cs[k] = expr_copy(D->pE[i * D->m + k]);
            coord = make_list(cs, D->m); free(cs);
        }
        bool dup = false;
        for (size_t t = 0; t < np; t++) if (expr_eq(pts[t], coord)) { dup = true; break; }
        if (dup) { expr_free(coord); continue; }
        pts[np++] = coord;
    }
    Expr* lst = make_list(pts, np);   /* adopts pts[] elements (no refcount bump) */
    free(pts);
    char* s = expr_to_string(lst);
    expr_free(lst);
    fprintf(stderr, "InterpolatingPolynomial::poised: The interpolation points %s "
                    "are not poised, so an interpolating polynomial of total "
                    "degree %d could not be found.\n", s ? s : "{}", d_sq);
    free(s);
    return NULL;
}

/* Engine B driver: raise the degree to d_sq; return the first consistent
 * interpolant, else the poised/noipf diagnostic. */
static Expr* ip_solve_general(IPData* D, Expr** vars) {
    size_t N = D->n, m = D->m;
    int d_sq = 0;
    int64_t count = 0;
    while (count < (int64_t)N) { count += ip_binom(d_sq + (int)m - 1, (int)m - 1); if (count >= (int64_t)N) break; d_sq++; }
    for (int d = 0; d <= d_sq; d++) {
        bool consistent = false;
        Expr* r = ip_solve_degree(D, d, vars, &consistent);
        if (consistent) return r;   /* r may be a valid polynomial (incl. 0) */
    }
    return ip_fail(D, d_sq);
}

/* ---- driver + option parsing ---- */

static Expr* ip_solve_modular(IPData* D, Expr** vars, const mpz_t n);  /* below */

static Expr* ip_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* Separate positional args from option rules (only Modulus is recognised). */
    Expr* pos[2] = { NULL, NULL };
    size_t npos = 0;
    mpz_t modulus; mpz_init_set_ui(modulus, 0);
    for (size_t i = 0; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_FUNCTION && a->data.function.arg_count == 2
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.args[0]->type == EXPR_SYMBOL
            && a->data.function.args[0]->data.symbol.name == SYM_Modulus) {
            Expr* rhs = a->data.function.args[1];
            if (rhs->type == EXPR_INTEGER) mpz_set_si(modulus, rhs->data.integer);
            else if (rhs->type == EXPR_BIGINT) mpz_set(modulus, rhs->data.bigint);
            else { mpz_clear(modulus); return NULL; }
            continue;
        }
        if (npos < 2) pos[npos] = a;
        npos++;
    }
    if (npos != 2) { mpz_clear(modulus); return builtin_arg_error("InterpolatingPolynomial", argc, 2, 2); }

    Expr* data = pos[0];
    Expr* vars = pos[1];
    size_t m;
    Expr** varr;
    if (vars->type == EXPR_SYMBOL) { m = 1; varr = &vars; }
    else if (interp_is_list(vars) && vars->data.function.arg_count >= 1) {
        m = vars->data.function.arg_count;
        varr = vars->data.function.args;
        for (size_t k = 0; k < m; k++) if (varr[k]->type != EXPR_SYMBOL) { mpz_clear(modulus); return NULL; }
    } else { mpz_clear(modulus); return NULL; }

    IPData D;
    if (!ip_collect(data, m, &D)) { ipdata_free(&D); mpz_clear(modulus); return NULL; }

    Expr* result;
    if (mpz_sgn(modulus) > 0) {
        result = ip_solve_modular(&D, varr, modulus);
    } else if (m == 1 && !D.any_automatic) {
        result = ip_newton_1d(&D, varr[0]);
        if (!result) result = ip_solve_general(&D, varr);   /* degenerate fallback */
    } else {
        result = ip_solve_general(&D, varr);
    }
    ipdata_free(&D);
    mpz_clear(modulus);
    return result;
}

/* NDArray / packed intake: delist a packed or visible NDArray data table once,
 * then process it as a List (mirrors builtin_interpolation). */
static Expr* builtin_interpolatingpolynomial(Expr* res) {
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count >= 1
        && is_ndarray(res->data.function.args[0])) {
        Expr* orig = res->data.function.args[0];
        Expr* lst = ndarray_to_nested_list(orig);
        res->data.function.args[0] = lst;
        Expr* r = ip_impl(res);
        res->data.function.args[0] = orig;
        expr_free(lst);
        return r;
    }
    return ip_impl(res);
}

/* ---- Modulus->n path: the graded solve done exactly in Z/nZ ---- */

/* mpz residue in [0,n) -> Expr integer (or bigint). */
static Expr* ip_mpz_to_expr(const mpz_t z) {
    if (mpz_fits_slong_p(z)) return expr_new_integer((int64_t)mpz_get_si(z));
    return expr_new_bigint_from_mpz(z);
}

/* Reduce an exact value (Integer/BigInt/int64-Rational) mod n. Returns 1 on
 * success, 0 for an inexact/symbolic value or a non-invertible denominator. */
static int ip_mod_value(const Expr* v, const mpz_t n, mpz_t out) {
    if (v->type == EXPR_INTEGER) { mpz_set_si(out, (long)v->data.integer); mpz_mod(out, out, n); return 1; }
    if (v->type == EXPR_BIGINT)  { mpz_mod(out, v->data.bigint, n); return 1; }
    int64_t pn, pd;
    if (is_rational(v, &pn, &pd)) {
        mpz_t p, q, qi; mpz_init_set_si(p, (long)pn); mpz_init_set_si(q, (long)pd); mpz_init(qi);
        mpz_mod(q, q, n);
        int ok = mpz_invert(qi, q, n);
        if (ok) { mpz_mul(out, p, qi); mpz_mod(out, out, n); }
        mpz_clears(p, q, qi, NULL);
        return ok;
    }
    return 0;
}

static Expr* ip_solve_modular(IPData* D, Expr** vars, const mpz_t n) {
    size_t m = D->m, N = D->n;

    /* Coordinate and value residues (Modulus requires exact integer data). */
    mpz_t* pres = malloc(N * m * sizeof(mpz_t));
    mpz_t* bres = malloc(N * sizeof(mpz_t));
    bool bad = false;
    for (size_t t = 0; t < N * m; t++) mpz_init(pres[t]);
    for (size_t i = 0; i < N; i++) mpz_init(bres[i]);
    for (size_t i = 0; i < N && !bad; i++) {
        if (!ip_mod_value(D->val[i], n, bres[i])) bad = true;
        for (size_t k = 0; k < m; k++)
            if (!ip_mod_value(D->pE[i * m + k], n, pres[i * m + k])) bad = true;
    }
    if (bad) {
        for (size_t t = 0; t < N * m; t++) mpz_clear(pres[t]);
        for (size_t i = 0; i < N; i++) mpz_clear(bres[i]);
        free(pres); free(bres);
        return NULL;
    }

    int d_sq = 0; int64_t count = 0;
    while (count < (int64_t)N) { count += ip_binom(d_sq + (int)m - 1, (int)m - 1); if (count >= (int64_t)N) break; d_sq++; }

    Expr* result = NULL; bool solved = false;
    mpz_t inv, fac, prod, tpow; mpz_inits(inv, fac, prod, tpow, NULL);
    for (int d = 0; d <= d_sq && !solved; d++) {
        size_t L = 0;
        int* E = ip_monomials(d, m, &L);
        mpz_t* A = malloc(N * L * sizeof(mpz_t));
        mpz_t* b = malloc(N * sizeof(mpz_t));
        int64_t* pivrow = malloc(L * sizeof(int64_t));
        for (size_t t = 0; t < N * L; t++) mpz_init(A[t]);
        for (size_t i = 0; i < N; i++) { mpz_init(b[i]); mpz_set(b[i], bres[i]); }
        for (size_t j = 0; j < L; j++) pivrow[j] = -1;

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < L; j++) {
                int64_t coeff = 1;
                for (size_t k = 0; k < m; k++) { coeff *= ip_falling_i(E[j * m + k], D->alpha[i * m + k]); if (coeff == 0) break; }
                mpz_t* aij = &A[i * L + j];
                if (coeff == 0) { mpz_set_ui(*aij, 0); continue; }
                mpz_set_si(*aij, (long)coeff); mpz_mod(*aij, *aij, n);
                for (size_t k = 0; k < m; k++) {
                    int e = E[j * m + k] - D->alpha[i * m + k];
                    if (e > 0) { mpz_powm_ui(tpow, pres[i * m + k], (unsigned long)e, n); mpz_mul(*aij, *aij, tpow); mpz_mod(*aij, *aij, n); }
                }
            }
        }

        size_t prow = 0;
        for (size_t j = 0; j < L && prow < N; j++) {
            size_t piv = N;
            for (size_t r = prow; r < N; r++)
                if (mpz_sgn(A[r * L + j]) != 0 && mpz_invert(inv, A[r * L + j], n)) { piv = r; break; }
            if (piv == N) continue;
            if (piv != prow) {
                for (size_t jj = 0; jj < L; jj++) mpz_swap(A[piv * L + jj], A[prow * L + jj]);
                mpz_swap(b[piv], b[prow]);
            }
            mpz_invert(inv, A[prow * L + j], n);
            for (size_t jj = j; jj < L; jj++) { mpz_mul(A[prow * L + jj], A[prow * L + jj], inv); mpz_mod(A[prow * L + jj], A[prow * L + jj], n); }
            mpz_mul(b[prow], b[prow], inv); mpz_mod(b[prow], b[prow], n);
            for (size_t r = 0; r < N; r++) {
                if (r == prow || mpz_sgn(A[r * L + j]) == 0) continue;
                mpz_set(fac, A[r * L + j]);
                for (size_t jj = j; jj < L; jj++) { mpz_mul(prod, fac, A[prow * L + jj]); mpz_sub(A[r * L + jj], A[r * L + jj], prod); mpz_mod(A[r * L + jj], A[r * L + jj], n); }
                mpz_mul(prod, fac, b[prow]); mpz_sub(b[r], b[r], prod); mpz_mod(b[r], b[r], n);
            }
            pivrow[j] = (int64_t)prow; prow++;
        }

        bool ok = true;
        for (size_t r = prow; r < N; r++) if (mpz_sgn(b[r]) != 0) { ok = false; break; }
        if (ok) {
            Expr** C = malloc(L * sizeof(Expr*));
            for (size_t j = 0; j < L; j++)
                C[j] = (pivrow[j] >= 0) ? ip_mpz_to_expr(b[pivrow[j]]) : expr_new_integer(0);
            result = ip_build_poly(E, C, L, m, vars);
            for (size_t j = 0; j < L; j++) expr_free(C[j]);
            free(C);
            solved = true;
        }

        for (size_t t = 0; t < N * L; t++) mpz_clear(A[t]);
        for (size_t i = 0; i < N; i++) mpz_clear(b[i]);
        free(A); free(b); free(pivrow); free(E);
    }
    mpz_clears(inv, fac, prod, tpow, NULL);
    for (size_t t = 0; t < N * m; t++) mpz_clear(pres[t]);
    for (size_t i = 0; i < N; i++) mpz_clear(bres[i]);
    free(pres); free(bres);

    if (!solved) result = ip_fail(D, d_sq);
    return result;
}

/* --- registration ------------------------------------------------------ */

void interp_init(void) {
    SymbolDef* def = symtab_get_def("InterpolatingFunction");
    /* HoldAll: the object's `domain` and `table` arguments are constant,
     * already-evaluated numeric data (built by Interpolation[] or by
     * interp_make_derivative).  Without holding them, every application
     * `ifun[x]` makes the evaluator re-evaluate the entire N-point table on
     * each call -- an O(N) deep re-traversal that dominated per-point cost.
     * Holding them keeps the object a true normal form and makes evaluating it
     * O(1); the application's own argument `x` is still evaluated normally
     * (attributes of the object head, not InterpolatingFunction, govern that). */
    def->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    symtab_add_builtin("Interpolation", builtin_interpolation);
    SymbolDef* idef = symtab_get_def("Interpolation");
    idef->attributes |= ATTR_PROTECTED;
    /* Options[Interpolation] = {InterpolationOrder -> 3, Method -> Automatic,
     *                            PeriodicInterpolation -> False} — mirrors the
     * runtime option parser in builtin_interpolation_impl. */
    Expr* iopts = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_InterpolationOrder), expr_new_integer(3) }, 2),
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic) }, 2),
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_PeriodicInterpolation), expr_new_symbol(SYM_False) }, 2),
    }, 3);
    symtab_set_options("Interpolation", iopts);   /* takes ownership */

    symtab_add_builtin("InterpolatingPolynomial", builtin_interpolatingpolynomial);
    SymbolDef* pdef = symtab_get_def("InterpolatingPolynomial");
    pdef->attributes |= ATTR_PROTECTED;
    /* Options[InterpolatingPolynomial] = {Modulus -> 0}. */
    Expr* mrule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_new_symbol(SYM_Modulus), expr_new_integer(0) }, 2);
    Expr* popts = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ mrule }, 1);
    symtab_set_options("InterpolatingPolynomial", popts);
}

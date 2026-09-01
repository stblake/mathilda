#include "list_common.h"
#include "ndarray.h"     /* is_ndarray, ndarray_size, expr_new_ndarray_like, NDType */
#include "ndstruct.h"    /* ndstruct_delist_repack — packed/declined fallback */
#include "options.h"     /* options_extract / OptEntry */
#include "list_gradient.h"

/* ---- ListGradient --------------------------------------------------------
 *
 * A from-first-principles port of numpy.gradient: the numerical gradient of a
 * sampled array, by finite differences, per axis, at any rank.
 *
 *   ListGradient[f]                    unit spacing on every axis
 *   ListGradient[f, spacing]           uniform h, per-axis h, or coordinates
 *   ListGradient[f, spacing, opts...]  Method / DifferenceOrder / WindowLength / Axis
 *
 * DEFAULT semantics reproduce numpy.gradient exactly: interior points use a
 * second-order central difference and the two extreme endpoints a first-order
 * one-sided difference (numpy edge_order = 1). A rank-1 f gives one vector; a
 * rank-k f gives {g_1, ..., g_k}, one same-shape array per axis.
 *
 * ONE kernel underlies every Method / order / window / grid combination:
 * Fornberg's finite-difference-weight recurrence (B. Fornberg, "Generation of
 * finite difference formulas on arbitrarily spaced grids", Math. Comp. 51
 * (1988) 699-706). Given m grid nodes z[] and an evaluation point x0 it returns
 * the weights w[] of the first-derivative estimate f'(x0) ~= sum_k w[k] f(z[k]).
 * It is exact on non-uniform grids, so coordinate vectors need no special case.
 *
 * TWO instantiations of that one recurrence: fd_weights_double for the machine
 * buffer fast path (float64/float32 NDArray, packed or visible), and
 * fd_weights_expr which builds Expr weights the evaluator reduces — exact
 * Rationals for an exact grid, symbolic weights for a symbolic grid — so integer
 * and symbolic inputs answer exactly. The buffer path declines (returns NULL) for
 * anything it cannot represent in a float buffer, and the caller re-runs the
 * exact/symbolic List path via ndstruct_delist_repack.
 * -------------------------------------------------------------------------- */

typedef enum { LG_CENTERED, LG_FORWARD, LG_BACKWARD } LGMethod;

/* Largest finite-difference window we will build a stencil for. Well past any
 * sane accuracy order; guards against a WindowLength -> 10^6 blowing up the
 * per-(count,d) weight tables. */
#define LG_MAX_WINDOW 256

typedef struct {
    LGMethod method;
    int m;                          /* stencil size (window) */
    int rank;
    int64_t dims[NDARRAY_MAX_RANK];
    int axes[NDARRAY_MAX_RANK];     /* 0-based axes to differentiate, in order */
    int naxes;
    bool single_axis;               /* Axis -> <int>: caller wanted one axis */
    const Expr* spacing;            /* borrowed: NULL | scalar | List */
} LGParams;

/* ===================================================================== *
 *  Fornberg weights — double
 * ===================================================================== */

/* First-derivative weights at x0 for the `win` nodes nd[0..win-1]. Writes
 * w[0..win-1]. Standard Fornberg recurrence specialised to derivative order 1
 * (the c[.][0] and c[.][1] columns). */
static void fd_weights_double(const double* nd, int win, double x0, double* w) {
    double* c0 = calloc((size_t)win, sizeof(double)); /* derivative order 0 */
    double* c1a = calloc((size_t)win, sizeof(double)); /* derivative order 1 */
    double c1 = 1.0;
    double c4 = nd[0] - x0;
    c0[0] = 1.0;
    for (int i = 1; i < win; i++) {
        double c2 = 1.0;
        double c5 = c4;
        c4 = nd[i] - x0;
        for (int j = 0; j <= i - 1; j++) {
            double c3 = nd[i] - nd[j];
            c2 *= c3;
            if (j == i - 1) {
                c1a[i] = c1 * (c0[i - 1] - c5 * c1a[i - 1]) / c2;
                c0[i] = -c1 * c5 * c0[i - 1] / c2;
            }
            c1a[j] = (c4 * c1a[j] - c0[j]) / c3;
            c0[j] = c4 * c0[j] / c3;
        }
        c1 = c2;
    }
    for (int k = 0; k < win; k++) w[k] = c1a[k];
    free(c0);
    free(c1a);
}

/* ===================================================================== *
 *  Fornberg weights — Expr (exact / symbolic)
 * ===================================================================== */

static Expr* e_binop(const char* head, Expr* a, Expr* b) {
    Expr** args = malloc(2 * sizeof(Expr*));
    args[0] = expr_copy(a);
    args[1] = expr_copy(b);
    Expr* e = expr_new_function(expr_new_symbol(head), args, 2);
    free(args);
    return eval_and_free(e);
}
static Expr* e_sub(Expr* a, Expr* b) { return e_binop(SYM_Subtract, a, b); }
static Expr* e_mul(Expr* a, Expr* b) { return e_binop(SYM_Times, a, b); }
static Expr* e_div(Expr* a, Expr* b) { return e_binop(SYM_Divide, a, b); }
static Expr* e_neg(Expr* a) {
    Expr* z = expr_new_integer(-1);
    Expr* r = e_binop(SYM_Times, z, a);
    expr_free(z);
    return r;
}

/* Expr twin of fd_weights_double. nd[] and x0 are borrowed; w[0..win-1] are
 * returned as freshly-owned Exprs (the caller takes ownership). Exact-rational
 * for an integer/rational grid, symbolic for a symbolic grid. */
static void fd_weights_expr(Expr** nd, int win, Expr* x0, Expr** w) {
    Expr** c0 = malloc((size_t)win * sizeof(Expr*));
    Expr** c1a = malloc((size_t)win * sizeof(Expr*));
    for (int k = 0; k < win; k++) { c0[k] = expr_new_integer(0); c1a[k] = expr_new_integer(0); }
    Expr* c1 = expr_new_integer(1);
    Expr* c4 = e_sub(nd[0], x0);
    expr_free(c0[0]); c0[0] = expr_new_integer(1);
    for (int i = 1; i < win; i++) {
        Expr* c2 = expr_new_integer(1);
        Expr* c5 = expr_copy(c4);
        expr_free(c4); c4 = e_sub(nd[i], x0);
        for (int j = 0; j <= i - 1; j++) {
            Expr* c3 = e_sub(nd[i], nd[j]);
            Expr* nc2 = e_mul(c2, c3); expr_free(c2); c2 = nc2;
            if (j == i - 1) {
                /* c1a[i] = c1 * (c0[i-1] - c5*c1a[i-1]) / c2 */
                Expr* t1 = e_mul(c5, c1a[i - 1]);
                Expr* t2 = e_sub(c0[i - 1], t1); expr_free(t1);
                Expr* t3 = e_mul(c1, t2); expr_free(t2);
                Expr* t4 = e_div(t3, c2); expr_free(t3);
                expr_free(c1a[i]); c1a[i] = t4;
                /* c0[i] = -(c1*c5*c0[i-1]) / c2 */
                Expr* s1 = e_mul(c1, c5);
                Expr* s2 = e_mul(s1, c0[i - 1]); expr_free(s1);
                Expr* s3 = e_div(s2, c2); expr_free(s2);
                Expr* s4 = e_neg(s3); expr_free(s3);
                expr_free(c0[i]); c0[i] = s4;
            }
            /* c1a[j] = (c4*c1a[j] - c0[j]) / c3 */
            Expr* u1 = e_mul(c4, c1a[j]);
            Expr* u2 = e_sub(u1, c0[j]); expr_free(u1);
            Expr* u3 = e_div(u2, c3); expr_free(u2);
            expr_free(c1a[j]); c1a[j] = u3;
            /* c0[j] = c4*c0[j] / c3 */
            Expr* v1 = e_mul(c4, c0[j]);
            Expr* v2 = e_div(v1, c3); expr_free(v1);
            expr_free(c0[j]); c0[j] = v2;
            expr_free(c3);
        }
        expr_free(c1); c1 = c2;
        expr_free(c5);
    }
    for (int k = 0; k < win; k++) { w[k] = c1a[k]; expr_free(c0[k]); }
    free(c0);
    free(c1a);
    expr_free(c1);
    expr_free(c4);
}

/* ===================================================================== *
 *  Stencil selection
 * ===================================================================== */

/* For output index i along an axis of length n with window m, choose the node
 * range [*start, *start+count). numpy-matching boundaries: Centered uses a
 * symmetric stencil in the interior, shifting inward where it would run off an
 * edge, and drops the two EXTREME endpoints to a one-sided (m-1)-point stencil
 * (accuracy order p-1); Forward/Backward stay order p everywhere. */
static int lg_stencil(int64_t i, int64_t n, int m, LGMethod method, int64_t* start) {
    int count = m;
    if ((int64_t)count > n) count = (int)n; /* not enough points: use them all */
    int64_t s;
    if (method == LG_FORWARD) {
        s = i;
        if (s > n - count) s = n - count;
        if (s < 0) s = 0;
    } else if (method == LG_BACKWARD) {
        s = i - (count - 1);
        if (s < 0) s = 0;
        if (s > n - count) s = n - count;
    } else { /* LG_CENTERED */
        if (count >= 3 && (i == 0 || i == n - 1)) {
            count -= 1;                 /* endpoint reduction -> order p-1 */
            s = (i == 0) ? 0 : n - count;
        } else {
            s = i - count / 2;
            if (s < 0) s = 0;
            if (s > n - count) s = n - count;
        }
    }
    *start = s;
    return count;
}

/* ===================================================================== *
 *  Option / spacing parsing
 * ===================================================================== */

static bool lg_method_from_string(const Expr* v, LGMethod* out) {
    if (!v || v->type != EXPR_STRING) return false;
    const char* s = v->data.string;
    if (!strcmp(s, "Centered") || !strcmp(s, "Central")) { *out = LG_CENTERED; return true; }
    if (!strcmp(s, "Forward") || !strcmp(s, "Forwards")) { *out = LG_FORWARD; return true; }
    if (!strcmp(s, "Backward") || !strcmp(s, "Backwards")) { *out = LG_BACKWARD; return true; }
    return false;
}

/* Shape of a nested rectangular List (spine of first elements). */
static bool lg_list_shape(const Expr* f, int64_t* dims, int* rank) {
    int r = 0;
    const Expr* cur = f;
    while (cur->type == EXPR_FUNCTION
           && cur->data.function.head->type == EXPR_SYMBOL
           && cur->data.function.head->data.symbol.name == SYM_List) {
        if (r >= NDARRAY_MAX_RANK) return false;
        dims[r] = (int64_t)cur->data.function.arg_count;
        r++;
        if (cur->data.function.arg_count == 0) break;
        cur = cur->data.function.args[0];
    }
    *rank = r;
    return r >= 1;
}

static bool lg_parse(const Expr* res, LGParams* p, size_t* positional_argc) {
    if (res->type != EXPR_FUNCTION) return false;

    const Expr* v_method = NULL; bool g_method = false;
    const Expr* v_order  = NULL; bool g_order  = false;
    const Expr* v_window = NULL; bool g_window = false;
    const Expr* v_axis   = NULL; bool g_axis   = false;
    OptEntry ents[4] = {
        { "Method",          &v_method, &g_method },
        { "DifferenceOrder", &v_order,  &g_order  },
        { "WindowLength",    &v_window, &g_window },
        { "Axis",            &v_axis,   &g_axis   },
    };
    size_t argc = res->data.function.arg_count;
    if (!options_extract(res, "ListGradient", ents, 4, &argc)) return false;
    if (argc < 1 || argc > 2) return false;
    *positional_argc = argc;

    p->method = LG_CENTERED;
    if (v_method && !lg_method_from_string(v_method, &p->method)) return false;

    int64_t order = 2;
    if (v_order) {
        if (v_order->type != EXPR_INTEGER || v_order->data.integer < 1) return false;
        order = v_order->data.integer;
    }
    int64_t m;
    if (v_window && !(v_window->type == EXPR_SYMBOL
                      && v_window->data.symbol.name == SYM_Automatic)) {
        if (v_window->type != EXPR_INTEGER || v_window->data.integer < 2) return false;
        m = v_window->data.integer;
    } else {
        m = order + 1;
    }
    if (m < 2) m = 2;
    if (m > LG_MAX_WINDOW) return false;
    p->m = (int)m;

    Expr* f = res->data.function.args[0];
    if (is_ndarray(f)) {
        p->rank = f->data.ndarray.rank;
        for (int i = 0; i < p->rank; i++) p->dims[i] = f->data.ndarray.dims[i];
    } else if (is_listq(f)) {
        if (!lg_list_shape(f, p->dims, &p->rank)) return false;
    } else {
        return false;
    }
    if (p->rank < 1) return false;

    p->single_axis = false;
    if (!v_axis || (v_axis->type == EXPR_SYMBOL && v_axis->data.symbol.name == SYM_All)) {
        p->naxes = p->rank;
        for (int i = 0; i < p->rank; i++) p->axes[i] = i;
    } else if (v_axis->type == EXPR_INTEGER) {
        int64_t a = v_axis->data.integer;
        int ax = (a > 0) ? (int)(a - 1) : (int)(p->rank + a);
        if (ax < 0 || ax >= p->rank) return false;
        p->axes[0] = ax; p->naxes = 1; p->single_axis = true;
    } else if (is_listq((Expr*)v_axis)) {
        size_t k = v_axis->data.function.arg_count;
        if (k < 1 || k > (size_t)p->rank) return false;
        for (size_t i = 0; i < k; i++) {
            Expr* ae = v_axis->data.function.args[i];
            if (ae->type != EXPR_INTEGER) return false;
            int64_t a = ae->data.integer;
            int ax = (a > 0) ? (int)(a - 1) : (int)(p->rank + a);
            if (ax < 0 || ax >= p->rank) return false;
            p->axes[i] = ax;
        }
        p->naxes = (int)k;
    } else {
        return false;
    }

    p->spacing = (argc == 2) ? res->data.function.args[1] : NULL;
    return true;
}

/* The spacing spec for the ordinal-th computed axis (borrowed). See the header
 * of the source for the disambiguation rules. Returns false on a malformed
 * per-axis list length. */
static bool lg_axis_spec(const LGParams* p, int ordinal, const Expr** out) {
    const Expr* sp = p->spacing;
    if (sp == NULL) { *out = NULL; return true; }
    if (!is_listq((Expr*)sp)) { *out = sp; return true; } /* scalar / symbol / unit */
    size_t len = sp->data.function.arg_count;
    if (p->naxes == 1) {
        int64_t axislen = p->dims[p->axes[0]];
        if ((int64_t)len == axislen) { *out = sp; return true; }          /* coordinates */
        if (len == 1) { *out = sp->data.function.args[0]; return true; }  /* [spec] */
        return false;
    }
    if (len != (size_t)p->naxes) return false;
    *out = sp->data.function.args[ordinal];
    return true;
}

/* ===================================================================== *
 *  Buffer fast path (float64 / float32)
 * ===================================================================== */

/* Gradient of `src` (row-major double buffer, shape dims/rank) along axis ax,
 * into `dst` (same shape). `coords` NULL -> uniform spacing h; else a
 * length-dims[ax] coordinate array (h ignored). `ncomp` is 1 for a real buffer
 * and 2 for a complex buffer stored as interleaved (re, im) doubles: the
 * finite-difference weights are real, so a complex element is just `ncomp`
 * contiguous doubles carried through the same axpy. */
static bool lg_axis_double(const double* restrict src, double* restrict dst,
                           const int64_t* dims, int rank, int ax, LGMethod method,
                           int m, double h, const double* coords, int ncomp) {
    int64_t n = dims[ax];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < ax; i++) outer *= (size_t)dims[i];
    for (int i = ax + 1; i < rank; i++) inner *= (size_t)dims[i];
    size_t innerd = inner * (size_t)ncomp;   /* doubles per row (2x for complex) */
    int mm = m;
    if ((int64_t)mm > n) mm = (int)n;

    bool uniform = (coords == NULL);
    double hinv = uniform ? (1.0 / h) : 0.0;

    /* Uniform spacing: the weight vector depends only on (count, d=t-start), of
     * which there are a handful. Memoise into a small table indexed count*mm+d. */
    double** wt = NULL;
    /* Non-uniform: weights depend on the absolute index; precompute one vector
     * per t into a flat buffer reused across every line. */
    double* wflat = NULL;
    size_t* woff = NULL;

    if (uniform) {
        wt = calloc((size_t)(mm + 1) * mm, sizeof(double*));
        if (!wt) return false;
    } else {
        woff = malloc((size_t)n * sizeof(size_t));
        if (!woff) return false;
        size_t tot = 0;
        for (int64_t t = 0; t < n; t++) {
            int64_t st; int c = lg_stencil(t, n, mm, method, &st);
            woff[t] = tot; tot += (size_t)c;
        }
        wflat = malloc((tot ? tot : 1) * sizeof(double));
        if (!wflat) { free(woff); return false; }
        double nodes[LG_MAX_WINDOW] = {0};
        for (int64_t t = 0; t < n; t++) {
            int64_t st; int c = lg_stencil(t, n, mm, method, &st);
            for (int j = 0; j < c; j++) nodes[j] = coords[st + j];
            fd_weights_double(nodes, c, coords[t], wflat + woff[t]);
        }
    }

    /* Interior fast range [ilo, ihi): for uniform spacing every interior point
     * shares one stencil offset (start = t - id) and one weight vector, so the
     * whole run avoids the per-element lg_stencil() call and wt[] lookup that the
     * boundary/non-uniform path needs. This is what keeps the inner=1 shapes
     * (1-D, the strided axis) close to a raw finite-difference loop. Empty for
     * the non-uniform (coords) path, which keeps the general per-element loop. */
    int64_t ilo = 0, ihi = 0;
    int id = 0;
    double w_int[LG_MAX_WINDOW];
    if (uniform) {
        if (method == LG_FORWARD)       { id = 0;        ilo = 0;        ihi = n - mm + 1; }
        else if (method == LG_BACKWARD) { id = mm - 1;   ilo = mm - 1;   ihi = n; }
        else                            { id = mm / 2;   ilo = mm / 2;   ihi = n - mm + id + 1; }
        if (ilo < 0) ilo = 0;
        if (ihi > n) ihi = n;
        if (ihi < ilo) ihi = ilo;
        double nodes[LG_MAX_WINDOW] = {0};
        for (int j = 0; j < mm; j++) nodes[j] = (double)(j - id);
        fd_weights_double(nodes, mm, 0.0, w_int);
        for (int j = 0; j < mm; j++) w_int[j] *= hinv;
    }

    for (size_t o = 0; o < outer; o++) {
        size_t base = o * (size_t)n * innerd;
        int64_t t = 0;
        while (t < n) {
            if (t >= ilo && t < ihi) {
                /* Interior as a vectorizable tap-outer axpy over the contiguous
                 * block [ilo, ihi) x innerd: one straight `dst[i] (+)= w*src[i+off]`
                 * pass per stencil tap, which the compiler turns into SIMD FMA
                 * (restrict promises src/dst do not alias). This is what lets the
                 * inner=1 shapes (1-D and the strided axis) vectorize like the
                 * contiguous-axis shape, matching numpy's a*f[:-2]+b*f[1:-1]+... */
                double* od = dst + base + (size_t)ilo * innerd;
                const double* sb = src + base + (size_t)ilo * innerd;
                size_t len = (size_t)(ihi - ilo) * innerd;
                ptrdiff_t st0 = (ptrdiff_t)innerd;
                /* Fuse the common small windows into ONE pass (1 write, not mm):
                 * the default centered 3-point and the 2-point end schemes cover
                 * the overwhelming majority of calls, and a fused unrolled pass
                 * halves/thirds the memory-write traffic of the general
                 * tap-outer while still emitting SIMD FMA. */
                if (mm == 3) {
                    const double* s0 = sb + (ptrdiff_t)(0 - id) * st0;
                    const double* s1 = sb + (ptrdiff_t)(1 - id) * st0;
                    const double* s2 = sb + (ptrdiff_t)(2 - id) * st0;
                    double a = w_int[0], b = w_int[1], c = w_int[2];
                    for (size_t i = 0; i < len; i++) od[i] = a * s0[i] + b * s1[i] + c * s2[i];
                } else if (mm == 2) {
                    const double* s0 = sb + (ptrdiff_t)(0 - id) * st0;
                    const double* s1 = sb + (ptrdiff_t)(1 - id) * st0;
                    double a = w_int[0], b = w_int[1];
                    for (size_t i = 0; i < len; i++) od[i] = a * s0[i] + b * s1[i];
                } else {
                    const double* s0 = sb + (ptrdiff_t)(-id) * st0;
                    double w0 = w_int[0];
                    for (size_t i = 0; i < len; i++) od[i] = w0 * s0[i];
                    for (int k = 1; k < mm; k++) {
                        double wk = w_int[k];
                        const double* sk = sb + (ptrdiff_t)(k - id) * st0;
                        for (size_t i = 0; i < len; i++) od[i] += wk * sk[i];
                    }
                }
                t = ihi;
                continue;
            }
            int64_t st; int count = lg_stencil(t, n, mm, method, &st);
            const double* w;
            if (uniform) {
                int d = (int)(t - st);
                int idx = count * mm + d;
                if (!wt[idx]) {
                    double* cell = malloc((size_t)count * sizeof(double));
                    if (!cell) goto fail;
                    double nodes[LG_MAX_WINDOW] = {0};
                    for (int j = 0; j < count; j++) nodes[j] = (double)(j - d);
                    fd_weights_double(nodes, count, 0.0, cell);
                    for (int j = 0; j < count; j++) cell[j] *= hinv;
                    wt[idx] = cell;
                }
                w = wt[idx];
            } else {
                w = wflat + woff[t];
            }
            double* orow = dst + base + (size_t)t * innerd;
            const double* s0 = src + base + (size_t)st * innerd;
            double w0 = w[0];
            for (size_t in = 0; in < innerd; in++) orow[in] = w0 * s0[in];
            for (int k = 1; k < count; k++) {        /* k=0 folded into the init */
                double wk = w[k];
                const double* srow = s0 + (size_t)k * innerd;
                for (size_t in = 0; in < innerd; in++) orow[in] += wk * srow[in];
            }
            t++;
        }
    }

    if (wt) { for (int i = 0; i < (mm + 1) * mm; i++) free(wt[i]); free(wt); }
    free(wflat);
    free(woff);
    return true;

fail:
    if (wt) { for (int i = 0; i < (mm + 1) * mm; i++) free(wt[i]); free(wt); }
    free(wflat);
    free(woff);
    return false;
}

/* Resolve the spacing spec for computed-axis ordinal `i` (whose length is `n`)
 * into either uniform spacing `*h` or a freshly-allocated `*coords` array (the
 * caller frees it). Returns false to decline (symbolic/non-numeric/zero/bad
 * length), which sends the call to the exact/symbolic List path. */
static bool lg_resolve_axis_spacing(const LGParams* p, int i, int64_t n,
                                    double* h, double** coords) {
    *h = 1.0; *coords = NULL;
    const Expr* spec;
    if (!lg_axis_spec(p, i, &spec)) return false;
    if (spec == NULL) return true;
    if (is_listq((Expr*)spec)) {
        if (spec->data.function.arg_count != (size_t)n) return false;
        double* c = malloc((size_t)n * sizeof(double));
        if (!c) return false;
        for (int64_t j = 0; j < n; j++)
            if (!common_machine_real_value(spec->data.function.args[j], &c[j])) {
                free(c); return false;
            }
        *coords = c;
        return true;
    }
    if (!common_machine_real_value(spec, h) || *h == 0.0) return false;
    return true;
}

/* Narrow a double buffer of `n` components into a fresh float buffer. */
static float* lg_narrow_to_float(double* d, size_t n) {
    float* f = malloc(n * sizeof(float));
    if (!f) return NULL;
    for (size_t k = 0; k < n; k++) f[k] = (float)d[k];
    return f;
}

Expr* list_gradient_ndarray(Expr* res) {
    LGParams p; size_t pc;
    if (!lg_parse(res, &p, &pc)) return NULL;

    Expr* a = res->data.function.args[0];
    if (!is_ndarray(a)) return NULL;
    NDType dt = a->data.ndarray.dtype;
    /* Real and complex machine buffers compute here; int64/bool take the exact
     * List path (an integer gradient is Rational, a bool one is symbolic). The
     * finite-difference weights are real, so a complex buffer is handled by the
     * same kernel over its interleaved (re, im) doubles. */
    if (dt != NDT_FLOAT64 && dt != NDT_FLOAT32
        && dt != NDT_COMPLEX64 && dt != NDT_COMPLEX32) return NULL;
    int ncomp = ndt_components(dt);                 /* 1 real, 2 complex */
    bool is32 = (dt == NDT_FLOAT32 || dt == NDT_COMPLEX32);

    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    for (int i = 0; i < p.naxes; i++)
        if (dims[p.axes[i]] < 2) return NULL; /* gradient undefined on a length<2 axis */

    size_t total = ndarray_size(a);
    size_t ndbl = (total ? total : 1) * (size_t)ncomp;  /* doubles in the buffer */

    const double* src;
    double* srcbuf = NULL;
    if (!is32) {
        src = (const double*)a->data.ndarray.data;  /* float64 / complex64 already double */
    } else {
        srcbuf = malloc(ndbl * sizeof(double));
        if (!srcbuf) return NULL;
        const float* sf = (const float*)a->data.ndarray.data;
        for (size_t k = 0; k < total * (size_t)ncomp; k++) srcbuf[k] = (double)sf[k];
        src = srcbuf;
    }

    /* Single computed axis: one buffer, one array. */
    if (p.naxes == 1) {
        int ax = p.axes[0];
        double h; double* coords;
        if (!lg_resolve_axis_spacing(&p, 0, dims[ax], &h, &coords)) { free(srcbuf); return NULL; }
        double* dst = malloc(ndbl * sizeof(double));
        if (!dst) { free(coords); free(srcbuf); return NULL; }
        bool okk = lg_axis_double(src, dst, dims, rank, ax, p.method, p.m, h, coords, ncomp);
        free(coords); free(srcbuf);
        if (!okk) { free(dst); return NULL; }
        if (is32) {
            float* df = lg_narrow_to_float(dst, ndbl); free(dst);
            return df ? expr_new_ndarray_like(a, rank, dims, df, dt) : NULL;
        }
        return expr_new_ndarray_like(a, rank, dims, dst, dt);
    }

    /* Multiple axes on a PACKED input: build the stacked rank-(rank+1) output
     * directly, each axis gradient written into its own contiguous slice. This
     * is exactly the shape the transparency gate would otherwise produce by
     * copying naxes separate arrays into one — so building it in place skips that
     * copy (it was ~3x the kernel cost on a 1000x1000). A VISIBLE NDArray is not
     * stacked by the gate, so it keeps the List-of-arrays shape below. */
    if (is_packed_list(a) && rank + 1 <= NDARRAY_MAX_RANK) {
        double* big = malloc((size_t)p.naxes * ndbl * sizeof(double));
        if (!big) { free(srcbuf); return NULL; }
        bool okk = true;
        for (int i = 0; i < p.naxes && okk; i++) {
            int ax = p.axes[i];
            double h; double* coords;
            if (!lg_resolve_axis_spacing(&p, i, dims[ax], &h, &coords)) { okk = false; break; }
            okk = lg_axis_double(src, big + (size_t)i * ndbl, dims, rank, ax,
                                 p.method, p.m, h, coords, ncomp);
            free(coords);
        }
        free(srcbuf);
        if (!okk) { free(big); return NULL; }
        int64_t dims2[NDARRAY_MAX_RANK + 1];
        dims2[0] = p.naxes;
        for (int j = 0; j < rank; j++) dims2[j + 1] = dims[j];
        if (is32) {
            float* bigf = lg_narrow_to_float(big, (size_t)p.naxes * ndbl); free(big);
            return bigf ? expr_new_ndarray_like(a, rank + 1, dims2, bigf, dt) : NULL;
        }
        return expr_new_ndarray_like(a, rank + 1, dims2, big, dt);
    }

    /* Multiple axes, visible NDArray (or rank at the cap): a List of arrays. */
    Expr** outs = malloc((size_t)p.naxes * sizeof(Expr*));
    if (!outs) { free(srcbuf); return NULL; }
    int built = 0;
    bool ok = true;
    for (int i = 0; i < p.naxes && ok; i++) {
        int ax = p.axes[i];
        double h; double* coords;
        if (!lg_resolve_axis_spacing(&p, i, dims[ax], &h, &coords)) { ok = false; break; }
        double* dst = malloc(ndbl * sizeof(double));
        if (!dst) { free(coords); ok = false; break; }
        bool okk = lg_axis_double(src, dst, dims, rank, ax, p.method, p.m, h, coords, ncomp);
        free(coords);
        if (!okk) { free(dst); ok = false; break; }
        if (is32) {
            float* df = lg_narrow_to_float(dst, ndbl); free(dst);
            if (!df) { ok = false; break; }
            outs[built] = expr_new_ndarray_like(a, rank, dims, df, dt);
        } else {
            outs[built] = expr_new_ndarray_like(a, rank, dims, dst, dt);
        }
        built++;
    }
    free(srcbuf);
    if (!ok) {
        for (int i = 0; i < built; i++) expr_free(outs[i]);
        free(outs);
        return NULL;
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), outs, (size_t)p.naxes);
    free(outs);
    return result;
}

/* ===================================================================== *
 *  Exact / symbolic List path
 * ===================================================================== */

/* Gradient of the top level of `values` (a List) along its own axis 0, using
 * spec (NULL=unit / scalar h / coordinate List). Returns a fresh List. */
static Expr* lg_grad_along_top(const Expr* values, const Expr* spec, LGMethod method, int m) {
    int64_t n = (int64_t)values->data.function.arg_count;
    int mm = m;
    if ((int64_t)mm > n) mm = (int)n;

    bool uniform = (spec == NULL) || !is_listq((Expr*)spec);
    Expr* hinv = NULL; /* Power[h,-1] for uniform, non-unit spacing */
    if (uniform && spec != NULL &&
        !(spec->type == EXPR_INTEGER && spec->data.integer == 1)) {
        Expr* neg1 = expr_new_integer(-1);
        hinv = e_binop(SYM_Power, (Expr*)spec, neg1);
        expr_free(neg1);
    }
    Expr** X = uniform ? NULL : spec->data.function.args; /* coordinates, borrowed */

    /* Exact offset-weight cache for uniform spacing, keyed (count, d). */
    Expr*** cache = uniform ? calloc((size_t)(mm + 1) * mm, sizeof(Expr**)) : NULL;

    Expr** out = malloc((size_t)n * sizeof(Expr*));
    for (int64_t t = 0; t < n; t++) {
        int64_t st; int count = lg_stencil(t, n, mm, method, &st);
        Expr** w = malloc((size_t)count * sizeof(Expr*)); /* owned per-point weights */

        if (uniform) {
            int d = (int)(t - st);
            int idx = count * mm + d;
            if (!cache[idx]) {
                Expr** nodes = malloc((size_t)count * sizeof(Expr*));
                for (int j = 0; j < count; j++) nodes[j] = expr_new_integer(j - d);
                Expr* x0 = expr_new_integer(0);
                Expr** ow = malloc((size_t)count * sizeof(Expr*));
                fd_weights_expr(nodes, count, x0, ow);
                for (int j = 0; j < count; j++) expr_free(nodes[j]);
                free(nodes); expr_free(x0);
                cache[idx] = ow;
            }
            for (int k = 0; k < count; k++) {
                if (hinv) w[k] = e_mul(cache[idx][k], hinv);
                else w[k] = expr_copy(cache[idx][k]);
            }
        } else {
            Expr** nodes = malloc((size_t)count * sizeof(Expr*));
            for (int j = 0; j < count; j++) nodes[j] = X[st + j];
            fd_weights_expr(nodes, count, X[t], w);
            free(nodes);
        }

        Expr** terms = malloc((size_t)count * sizeof(Expr*));
        for (int k = 0; k < count; k++) {
            Expr** ta = malloc(2 * sizeof(Expr*));
            ta[0] = w[k]; /* ownership transfers into the Times node */
            ta[1] = expr_copy(values->data.function.args[st + k]);
            terms[k] = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
            free(ta);
        }
        free(w);
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)count);
        free(terms);
        out[t] = eval_and_free(sum);
    }

    if (cache) {
        for (int idx = 0; idx < (mm + 1) * mm; idx++) {
            if (!cache[idx]) continue;
            int count = idx / mm;
            for (int k = 0; k < count; k++) expr_free(cache[idx][k]);
            free(cache[idx]);
        }
        free(cache);
    }
    if (hinv) expr_free(hinv);

    Expr* r = expr_new_function(expr_copy(values->data.function.head), out, (size_t)n);
    free(out);
    return r;
}

/* Gradient of `values` along axis `ax` (0-based), preserving shape. Recurses
 * into the leading axes so the differencing always happens along the top level. */
static Expr* lg_grad_axis_list(const Expr* values, int ax, const Expr* spec,
                               LGMethod method, int m) {
    if (ax == 0) return lg_grad_along_top(values, spec, method, m);
    size_t nn = values->data.function.arg_count;
    Expr** out = malloc((nn ? nn : 1) * sizeof(Expr*));
    for (size_t i = 0; i < nn; i++)
        out[i] = lg_grad_axis_list(values->data.function.args[i], ax - 1, spec, method, m);
    Expr* r = expr_new_function(expr_copy(values->data.function.head), out, nn);
    free(out);
    return r;
}

/* ===================================================================== *
 *  Entry point
 * ===================================================================== */

Expr* builtin_list_gradient(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Native buffer path first (packed or visible NDArray); ndstruct_delist_repack
     * defines every form it declines (int/complex dtype, symbolic spacing), which
     * re-enters the exact/symbolic List path below with a materialised List. */
    if (is_ndarray(lst)) {
        Expr* fast = list_gradient_ndarray(res);
        if (fast) return fast;
        return ndstruct_delist_repack(res, lst);
    }

    LGParams p; size_t pc;
    if (!lg_parse(res, &p, &pc)) return NULL;
    if (!is_listq(lst)) return NULL;

    for (int i = 0; i < p.naxes; i++) {
        if (p.dims[p.axes[i]] < 2) {
            ndarray_warn_once("ListGradient: each differentiated dimension must have length >= 2.");
            return NULL;
        }
    }

    Expr** outs = malloc((size_t)p.naxes * sizeof(Expr*));
    bool ok = true;
    int built = 0;
    for (int i = 0; i < p.naxes && ok; i++) {
        const Expr* spec;
        if (!lg_axis_spec(&p, i, &spec)) { ok = false; break; }
        outs[built] = lg_grad_axis_list(lst, p.axes[i], spec, p.method, p.m);
        built++;
    }
    if (!ok) {
        for (int i = 0; i < built; i++) expr_free(outs[i]);
        free(outs);
        return NULL;
    }

    Expr* result;
    if (p.naxes == 1) {
        result = outs[0];
    } else {
        result = expr_new_function(expr_new_symbol(SYM_List), outs, (size_t)p.naxes);
    }
    free(outs);
    return result;
}

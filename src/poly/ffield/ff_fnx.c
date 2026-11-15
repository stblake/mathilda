/*
 * ff_fnx.c -- Q(v_1..v_r)[t]/(mu) as a FLINT generic ring (see ff_fnx.h).
 *
 * Only the core methods are written here; gr_method_tab_init fills every other
 * slot with FLINT's generic implementation (powers, dot products, polynomial
 * and matrix algorithms), which call back into these.
 */
#include "ff_fnx.h"

#ifdef USE_FLINT

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  element storage                                                    */
/* ------------------------------------------------------------------ */

static void fnx_init(gr_ptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++) fmpz_mpoly_q_init(FF_FNX_COEFF(x, j), D->mctx);
}

static void fnx_clear(gr_ptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++) fmpz_mpoly_q_clear(FF_FNX_COEFF(x, j), D->mctx);
}

static void fnx_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++) fmpz_mpoly_q_swap(FF_FNX_COEFF(x, j), FF_FNX_COEFF(y, j), D->mctx);
}

static void fnx_set_shallow(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) {
    memcpy(x, y, (size_t)FF_FNX_DATA(ctx)->n * sizeof(fmpz_mpoly_q_struct));
}

static int fnx_zero(gr_ptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++) fmpz_mpoly_q_zero(FF_FNX_COEFF(x, j), D->mctx);
    return GR_SUCCESS;
}

static int fnx_one(gr_ptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    fnx_zero(x, ctx);
    fmpz_mpoly_q_one(FF_FNX_COEFF(x, 0), D->mctx);
    return GR_SUCCESS;
}

static int fnx_randtest(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) {
    (void)state;
    return fnx_zero(x, ctx);
}

static int fnx_set(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    if (x == y) return GR_SUCCESS;
    for (slong j = 0; j < D->n; j++)
        fmpz_mpoly_q_set(FF_FNX_COEFF(x, j), FF_FNX_COEFF((gr_ptr)y, j), D->mctx);
    return GR_SUCCESS;
}

static int fnx_set_si(gr_ptr x, slong c, gr_ctx_t ctx) {
    fnx_zero(x, ctx);
    fmpz_mpoly_q_set_si(FF_FNX_COEFF(x, 0), c, FF_FNX_DATA(ctx)->mctx);
    return GR_SUCCESS;
}

static int fnx_set_ui(gr_ptr x, ulong c, gr_ctx_t ctx) {
    fmpz_t z; fmpz_init_set_ui(z, c);
    fnx_zero(x, ctx);
    fmpz_mpoly_q_set_fmpz(FF_FNX_COEFF(x, 0), z, FF_FNX_DATA(ctx)->mctx);
    fmpz_clear(z);
    return GR_SUCCESS;
}

static int fnx_set_fmpz(gr_ptr x, const fmpz_t c, gr_ctx_t ctx) {
    fnx_zero(x, ctx);
    fmpz_mpoly_q_set_fmpz(FF_FNX_COEFF(x, 0), c, FF_FNX_DATA(ctx)->mctx);
    return GR_SUCCESS;
}

static int fnx_set_fmpq(gr_ptr x, const fmpq_t c, gr_ctx_t ctx) {
    fnx_zero(x, ctx);
    fmpz_mpoly_q_set_fmpq(FF_FNX_COEFF(x, 0), c, FF_FNX_DATA(ctx)->mctx);
    return GR_SUCCESS;
}

/* ------------------------------------------------------------------ */
/*  predicates                                                         */
/* ------------------------------------------------------------------ */

static truth_t fnx_is_zero(gr_srcptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++)
        if (!fmpz_mpoly_q_is_zero(FF_FNX_COEFF((gr_ptr)x, j), D->mctx)) return T_FALSE;
    return T_TRUE;
}

static truth_t fnx_is_one(gr_srcptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    if (!fmpz_mpoly_q_is_one(FF_FNX_COEFF((gr_ptr)x, 0), D->mctx)) return T_FALSE;
    for (slong j = 1; j < D->n; j++)
        if (!fmpz_mpoly_q_is_zero(FF_FNX_COEFF((gr_ptr)x, j), D->mctx)) return T_FALSE;
    return T_TRUE;
}

static truth_t fnx_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++)
        if (!fmpz_mpoly_q_equal(FF_FNX_COEFF((gr_ptr)x, j), FF_FNX_COEFF((gr_ptr)y, j), D->mctx))
            return T_FALSE;
    return T_TRUE;
}

/* ------------------------------------------------------------------ */
/*  arithmetic                                                         */
/* ------------------------------------------------------------------ */

static int fnx_neg(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++)
        fmpz_mpoly_q_neg(FF_FNX_COEFF(x, j), FF_FNX_COEFF((gr_ptr)y, j), D->mctx);
    return GR_SUCCESS;
}

static int fnx_add(gr_ptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++)
        fmpz_mpoly_q_add(FF_FNX_COEFF(x, j), FF_FNX_COEFF((gr_ptr)y, j), FF_FNX_COEFF((gr_ptr)z, j), D->mctx);
    return GR_SUCCESS;
}

static int fnx_sub(gr_ptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j < D->n; j++)
        fmpz_mpoly_q_sub(FF_FNX_COEFF(x, j), FF_FNX_COEFF((gr_ptr)y, j), FF_FNX_COEFF((gr_ptr)z, j), D->mctx);
    return GR_SUCCESS;
}

/* acc[0..2n-2] reduced modulo the monic mu into acc[0..n-1]. */
static void fnx_reduce(fmpz_mpoly_q_struct* acc, slong len, ff_fnx_data* D) {
    fmpz_mpoly_q_t c, s; fmpz_mpoly_q_init(c, D->mctx); fmpz_mpoly_q_init(s, D->mctx);
    for (slong k = len - 1; k >= D->n; k--) {
        if (fmpz_mpoly_q_is_zero(acc + k, D->mctx)) continue;
        fmpz_mpoly_q_swap(c, acc + k, D->mctx);
        fmpz_mpoly_q_zero(acc + k, D->mctx);
        /* t^k = t^(k-n) t^n = -t^(k-n) Sum_{j<n} mu_j t^j */
        for (slong j = 0; j < D->n; j++) {
            if (fmpq_is_zero(D->mu + j)) continue;
            fmpz_mpoly_q_mul_fmpq(s, c, D->mu + j, D->mctx);
            fmpz_mpoly_q_sub(acc + (k - D->n + j), acc + (k - D->n + j), s, D->mctx);
        }
    }
    fmpz_mpoly_q_clear(c, D->mctx); fmpz_mpoly_q_clear(s, D->mctx);
}

static int fnx_mul(gr_ptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    slong n = D->n;
    if (n == 1) {
        fmpz_mpoly_q_mul(FF_FNX_COEFF(x, 0), FF_FNX_COEFF((gr_ptr)y, 0), FF_FNX_COEFF((gr_ptr)z, 0), D->mctx);
        return GR_SUCCESS;
    }
    slong len = 2 * n - 1;
    fmpz_mpoly_q_struct* acc = flint_malloc((size_t)len * sizeof(fmpz_mpoly_q_struct));
    fmpz_mpoly_q_t p; fmpz_mpoly_q_init(p, D->mctx);
    for (slong k = 0; k < len; k++) fmpz_mpoly_q_init(acc + k, D->mctx);
    for (slong i = 0; i < n; i++) {
        if (fmpz_mpoly_q_is_zero(FF_FNX_COEFF((gr_ptr)y, i), D->mctx)) continue;
        for (slong j = 0; j < n; j++) {
            if (fmpz_mpoly_q_is_zero(FF_FNX_COEFF((gr_ptr)z, j), D->mctx)) continue;
            fmpz_mpoly_q_mul(p, FF_FNX_COEFF((gr_ptr)y, i), FF_FNX_COEFF((gr_ptr)z, j), D->mctx);
            fmpz_mpoly_q_add(acc + i + j, acc + i + j, p, D->mctx);
        }
    }
    fnx_reduce(acc, len, D);
    for (slong j = 0; j < n; j++) fmpz_mpoly_q_swap(FF_FNX_COEFF(x, j), acc + j, D->mctx);
    for (slong k = 0; k < len; k++) fmpz_mpoly_q_clear(acc + k, D->mctx);
    fmpz_mpoly_q_clear(p, D->mctx);
    flint_free(acc);
    return GR_SUCCESS;
}

/* 1/y: solve M c = e_0 with M the matrix of multiplication by y on the basis
 * 1, t, .., t^(n-1) (column k = y t^k), by Gaussian elimination over Q(v).
 * M is invertible exactly when y != 0 (K is a field). */
static int fnx_inv(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    slong n = D->n;
    if (fnx_is_zero(y, ctx) == T_TRUE) return GR_DOMAIN;
    if (n == 1) {
        fmpz_mpoly_q_inv(FF_FNX_COEFF(x, 0), FF_FNX_COEFF((gr_ptr)y, 0), D->mctx);
        return GR_SUCCESS;
    }
    /* augmented matrix A = [M | e_0], n x (n + 1), row-major */
    slong w = n + 1;
    fmpz_mpoly_q_struct* A = flint_malloc((size_t)(n * w) * sizeof(fmpz_mpoly_q_struct));
    for (slong k = 0; k < n * w; k++) fmpz_mpoly_q_init(A + k, D->mctx);
    /* column k: y t^k reduced */
    fmpz_mpoly_q_struct* col = flint_malloc((size_t)(2 * n - 1) * sizeof(fmpz_mpoly_q_struct));
    for (slong k = 0; k < 2 * n - 1; k++) fmpz_mpoly_q_init(col + k, D->mctx);
    for (slong k = 0; k < n; k++) {
        for (slong i = 0; i < 2 * n - 1; i++) fmpz_mpoly_q_zero(col + i, D->mctx);
        for (slong i = 0; i < n; i++) fmpz_mpoly_q_set(col + i + k, FF_FNX_COEFF((gr_ptr)y, i), D->mctx);
        fnx_reduce(col, 2 * n - 1, D);
        for (slong i = 0; i < n; i++) fmpz_mpoly_q_swap(A + i * w + k, col + i, D->mctx);
    }
    fmpz_mpoly_q_one(A + 0 * w + n, D->mctx);
    int st = GR_SUCCESS;
    fmpz_mpoly_q_t f, s; fmpz_mpoly_q_init(f, D->mctx); fmpz_mpoly_q_init(s, D->mctx);
    for (slong c = 0; c < n && st == GR_SUCCESS; c++) {
        slong piv = -1;
        for (slong r = c; r < n; r++) if (!fmpz_mpoly_q_is_zero(A + r * w + c, D->mctx)) { piv = r; break; }
        if (piv < 0) { st = GR_DOMAIN; break; }       /* mu reducible: not a field */
        if (piv != c)
            for (slong k = 0; k < w; k++) fmpz_mpoly_q_swap(A + piv * w + k, A + c * w + k, D->mctx);
        fmpz_mpoly_q_inv(f, A + c * w + c, D->mctx);
        for (slong k = c; k < w; k++) fmpz_mpoly_q_mul(A + c * w + k, A + c * w + k, f, D->mctx);
        for (slong r = 0; r < n; r++) {
            if (r == c || fmpz_mpoly_q_is_zero(A + r * w + c, D->mctx)) continue;
            fmpz_mpoly_q_set(f, A + r * w + c, D->mctx);
            for (slong k = c; k < w; k++) {
                if (fmpz_mpoly_q_is_zero(A + c * w + k, D->mctx)) continue;
                fmpz_mpoly_q_mul(s, f, A + c * w + k, D->mctx);
                fmpz_mpoly_q_sub(A + r * w + k, A + r * w + k, s, D->mctx);
            }
        }
    }
    if (st == GR_SUCCESS)
        for (slong i = 0; i < n; i++) fmpz_mpoly_q_swap(FF_FNX_COEFF(x, i), A + i * w + n, D->mctx);
    fmpz_mpoly_q_clear(f, D->mctx); fmpz_mpoly_q_clear(s, D->mctx);
    for (slong k = 0; k < 2 * n - 1; k++) fmpz_mpoly_q_clear(col + k, D->mctx);
    for (slong k = 0; k < n * w; k++) fmpz_mpoly_q_clear(A + k, D->mctx);
    flint_free(col); flint_free(A);
    return st;
}

static int fnx_div(gr_ptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) {
    gr_ptr iz = gr_heap_init(ctx);
    int st = fnx_inv(iz, z, ctx);
    if (st == GR_SUCCESS) st = fnx_mul(x, y, iz, ctx);
    gr_heap_clear(iz, ctx);
    return st;
}

/* ------------------------------------------------------------------ */
/*  printing and the context                                           */
/* ------------------------------------------------------------------ */

static int fnx_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    int st = gr_stream_write(out, "[");
    for (slong j = 0; j < D->n; j++) {
        char* s = fmpz_mpoly_q_get_str_pretty(FF_FNX_COEFF((gr_ptr)x, j), NULL, D->mctx);
        st |= gr_stream_write(out, s);
        flint_free(s);
        if (j + 1 < D->n) st |= gr_stream_write(out, ", ");
    }
    st |= gr_stream_write(out, "]");
    return st;
}

static int fnx_ctx_write(gr_stream_t out, gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    int st = gr_stream_write(out, "Q(v_1..v_");
    st |= gr_stream_write_si(out, D->nvars);
    st |= gr_stream_write(out, ")[t]/(mu), deg mu = ");
    st |= gr_stream_write_si(out, D->n);
    return st;
}

static void fnx_ctx_clear(gr_ctx_t ctx) {
    ff_fnx_data* D = FF_FNX_DATA(ctx);
    for (slong j = 0; j <= D->n; j++) fmpq_clear(D->mu + j);
    flint_free(D->mu);
    fmpz_mpoly_ctx_clear(D->mctx);
    flint_free(D);
}

static truth_t fnx_true(gr_ctx_t ctx) { (void)ctx; return T_TRUE; }

static int fnx_methods_ready = 0;
static gr_static_method_table fnx_methods;

/* FLINT's method tables store every method as gr_funcptr; the casts are its
 * documented idiom, so GCC's -Wcast-function-type is silenced for the table. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
static gr_method_tab_input fnx_methods_input[] = {
    {GR_METHOD_CTX_WRITE,               (gr_funcptr) fnx_ctx_write},
    {GR_METHOD_CTX_CLEAR,               (gr_funcptr) fnx_ctx_clear},
    {GR_METHOD_CTX_IS_RING,             (gr_funcptr) fnx_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) fnx_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) fnx_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) fnx_true},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) fnx_true},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) fnx_true},
    {GR_METHOD_INIT,        (gr_funcptr) fnx_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) fnx_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) fnx_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) fnx_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) fnx_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) fnx_write},
    {GR_METHOD_ZERO,        (gr_funcptr) fnx_zero},
    {GR_METHOD_ONE,         (gr_funcptr) fnx_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) fnx_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) fnx_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) fnx_equal},
    {GR_METHOD_SET,         (gr_funcptr) fnx_set},
    {GR_METHOD_SET_SI,      (gr_funcptr) fnx_set_si},
    {GR_METHOD_SET_UI,      (gr_funcptr) fnx_set_ui},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) fnx_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) fnx_set_fmpq},
    {GR_METHOD_NEG,         (gr_funcptr) fnx_neg},
    {GR_METHOD_ADD,         (gr_funcptr) fnx_add},
    {GR_METHOD_SUB,         (gr_funcptr) fnx_sub},
    {GR_METHOD_MUL,         (gr_funcptr) fnx_mul},
    {GR_METHOD_INV,         (gr_funcptr) fnx_inv},
    {GR_METHOD_DIV,         (gr_funcptr) fnx_div},
    {0,                     (gr_funcptr) NULL},
};
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

int ff_fnx_ctx_init(gr_ctx_t ctx, slong nvars, const fmpq_poly_t mu) {
    slong n = mu ? fmpq_poly_degree(mu) : 1;
    if (n < 1) return GR_DOMAIN;
    if (mu) {
        fmpq_t lc; fmpq_init(lc);
        fmpq_poly_get_coeff_fmpq(lc, mu, n);
        int monic = fmpq_is_one(lc);
        fmpq_clear(lc);
        if (!monic) return GR_DOMAIN;
    }
    ff_fnx_data* D = flint_malloc(sizeof(ff_fnx_data));
    D->nvars = nvars;
    D->n = n;
    fmpz_mpoly_ctx_init(D->mctx, nvars > 0 ? nvars : 1, ORD_LEX);
    D->mu = flint_malloc((size_t)(n + 1) * sizeof(fmpq));
    for (slong j = 0; j <= n; j++) {
        fmpq_init(D->mu + j);
        if (mu) fmpq_poly_get_coeff_fmpq(D->mu + j, mu, j);
        else fmpq_set_si(D->mu + j, j == 1 ? 1 : 0, 1);
    }
    ctx->which_ring = GR_CTX_UNKNOWN_DOMAIN;
    ctx->sizeof_elem = (slong)((size_t)n * sizeof(fmpz_mpoly_q_struct));
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = D;
    if (!fnx_methods_ready) {
        gr_method_tab_init(fnx_methods, fnx_methods_input);
        fnx_methods_ready = 1;
    }
    ctx->methods = fnx_methods;
    return GR_SUCCESS;
}

#endif /* USE_FLINT */

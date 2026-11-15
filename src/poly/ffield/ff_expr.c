/*
 * ff_expr.c -- Expr <-> Q(p)[t]/(mu) and its polynomials in x (see ff_expr.h).
 *
 * The input is read into one fmpz_mpoly_q over all the symbols {x, t, p_1..p_r}
 * (lex order), so a sum of fractions is combined exactly by FLINT; the
 * numerator is then split by (deg_x, deg_t) into coefficients over the
 * parameters.  The denominator must be free of x and t.
 */
#include "ff_expr.h"

#ifdef USE_FLINT

#include <stdlib.h>
#include <string.h>

#include "sym_names.h"
#include "flint_bridge_internal.h"

/* variable index of a symbol among {x, t, p...} (x = 0, t = 1, p_i = 2 + i), -1 if none */
static slong ffe_var(const Expr* e, const ff_names* nm) {
    const char* s = e->data.symbol.name;
    if (nm->x && strcmp(s, nm->x) == 0) return 0;
    if (nm->t && strcmp(s, nm->t) == 0) return 1;
    for (slong i = 0; i < nm->np; i++) if (strcmp(s, nm->p[i]) == 0) return 2 + i;
    return -1;
}

static int ffe_pow(fmpz_mpoly_q_t out, const fmpz_mpoly_q_t b, slong k, const fmpz_mpoly_ctx_t ctx) {
    if (k < 0) {
        if (fmpz_mpoly_q_is_zero(b, ctx)) return 0;
        fmpz_mpoly_q_t ib; fmpz_mpoly_q_init(ib, ctx);
        fmpz_mpoly_q_inv(ib, b, ctx);
        int ok = ffe_pow(out, ib, -k, ctx);
        fmpz_mpoly_q_clear(ib, ctx);
        return ok;
    }
    /* (n/d)^k with n, d coprime stays reduced */
    fmpz_mpoly_q_t r; fmpz_mpoly_q_init(r, ctx);
    int ok = fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(r), fmpz_mpoly_q_numref(b), (ulong)k, ctx) &&
             fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(r), fmpz_mpoly_q_denref(b), (ulong)k, ctx);
    if (ok) fmpz_mpoly_q_swap(out, r, ctx);
    fmpz_mpoly_q_clear(r, ctx);
    return ok;
}

/* e -> out over the big context ctx (vars x, t, p...). */
static int ffe_read(fmpz_mpoly_q_t out, const Expr* e, const ff_names* nm, const fmpz_mpoly_ctx_t ctx) {
    switch (e->type) {
    case EXPR_INTEGER:
        fmpz_mpoly_q_set_si(out, (slong)e->data.integer, ctx);
        return 1;
    case EXPR_BIGINT: {
        fmpz_t z; fmpz_init(z); fmpz_set_mpz(z, e->data.bigint);
        fmpz_mpoly_q_set_fmpz(out, z, ctx);
        fmpz_clear(z);
        return 1;
    }
    case EXPR_SYMBOL: {
        slong v = ffe_var(e, nm);
        if (v < 0) return 0;
        fmpz_mpoly_q_gen(out, v, ctx);
        return 1;
    }
    case EXPR_FUNCTION: break;
    default: return 0;
    }
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return 0;
    const char* hn = h->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    if (hn == SYM_Rational) {
        fmpq_t q; fmpq_init(q);
        int ok = fb_expr_to_fmpq(e, q);
        if (ok) fmpz_mpoly_q_set_fmpq(out, q, ctx);
        fmpq_clear(q);
        return ok;
    }
    if (hn == SYM_Plus || hn == SYM_Times) {
        fmpz_mpoly_q_t s; fmpz_mpoly_q_init(s, ctx);
        if (hn == SYM_Plus) fmpz_mpoly_q_zero(out, ctx); else fmpz_mpoly_q_one(out, ctx);
        int ok = 1;
        for (size_t i = 0; i < n && ok; i++) {
            ok = ffe_read(s, a[i], nm, ctx);
            if (ok) {
                if (hn == SYM_Plus) fmpz_mpoly_q_add(out, out, s, ctx);
                else fmpz_mpoly_q_mul(out, out, s, ctx);
            }
        }
        fmpz_mpoly_q_clear(s, ctx);
        return ok;
    }
    if (hn == SYM_Power && n == 2 && a[1]->type == EXPR_INTEGER) {
        fmpz_mpoly_q_t b; fmpz_mpoly_q_init(b, ctx);
        int ok = ffe_read(b, a[0], nm, ctx) && ffe_pow(out, b, (slong)a[1]->data.integer, ctx);
        fmpz_mpoly_q_clear(b, ctx);
        return ok;
    }
    return 0;
}

/* the coefficient (over the parameters) of x^dx t^dt in the big numerator,
 * divided by the big denominator (free of x, t), into c (K's mctx). */
static void ffe_split(fmpz_mpoly_q_struct* acc, slong dxmax, slong dtmax,
                      const fmpz_mpoly_q_t big, const fmpz_mpoly_ctx_t bctx,
                      const ff_names* nm, ff_fnx_data* D) {
    slong nb = 2 + nm->np;
    ulong* eb = flint_malloc((size_t)nb * sizeof(ulong));
    ulong* ek = flint_malloc((size_t)(D->nvars > 0 ? D->nvars : 1) * sizeof(ulong));
    fmpz_t c; fmpz_init(c);
    fmpz_mpoly_struct* nums = flint_malloc((size_t)((dxmax + 1) * (dtmax + 1)) * sizeof(fmpz_mpoly_struct));
    for (slong k = 0; k < (dxmax + 1) * (dtmax + 1); k++) fmpz_mpoly_init(nums + k, D->mctx);
    const fmpz_mpoly_struct* N = fmpz_mpoly_q_numref(big);
    for (slong i = 0; i < fmpz_mpoly_length(N, bctx); i++) {
        fmpz_mpoly_get_term_exp_ui(eb, N, i, bctx);
        fmpz_mpoly_get_term_coeff_fmpz(c, N, i, bctx);
        for (slong j = 0; j < (D->nvars > 0 ? D->nvars : 1); j++) ek[j] = (j < nm->np) ? eb[2 + j] : 0;
        fmpz_mpoly_push_term_fmpz_ui(nums + eb[0] * (ulong)(dtmax + 1) + eb[1], c, ek, D->mctx);
    }
    /* the denominator in K's context */
    fmpz_mpoly_t den; fmpz_mpoly_init(den, D->mctx);
    const fmpz_mpoly_struct* Dn = fmpz_mpoly_q_denref(big);
    for (slong i = 0; i < fmpz_mpoly_length(Dn, bctx); i++) {
        fmpz_mpoly_get_term_exp_ui(eb, Dn, i, bctx);
        fmpz_mpoly_get_term_coeff_fmpz(c, Dn, i, bctx);
        for (slong j = 0; j < (D->nvars > 0 ? D->nvars : 1); j++) ek[j] = (j < nm->np) ? eb[2 + j] : 0;
        fmpz_mpoly_push_term_fmpz_ui(den, c, ek, D->mctx);
    }
    fmpz_mpoly_sort_terms(den, D->mctx); fmpz_mpoly_combine_like_terms(den, D->mctx);
    for (slong k = 0; k < (dxmax + 1) * (dtmax + 1); k++) {
        fmpz_mpoly_sort_terms(nums + k, D->mctx); fmpz_mpoly_combine_like_terms(nums + k, D->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(acc + k), nums + k, D->mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_denref(acc + k), den, D->mctx);
        fmpz_mpoly_q_canonicalise(acc + k, D->mctx);
        fmpz_mpoly_clear(nums + k, D->mctx);
    }
    fmpz_mpoly_clear(den, D->mctx);
    fmpz_clear(c);
    flint_free(nums); flint_free(eb); flint_free(ek);
}

/* t-polynomial coefficients acc[0..len-1] (over the parameters) -> element of K */
static int ffe_elem_from_tpoly(gr_ptr out, fmpz_mpoly_q_struct* acc, slong len, gr_ctx_t K) {
    ff_fnx_data* D = FF_FNX_DATA(K);
    int st = gr_zero(out, K);
    gr_ptr tj = gr_heap_init(K), tgen = gr_heap_init(K), term = gr_heap_init(K);
    st |= gr_one(tj, K);
    st |= gr_zero(tgen, K);
    if (D->n > 1) fmpz_mpoly_q_one(FF_FNX_COEFF(tgen, 1), D->mctx);
    else {
        /* n = 1 means mu = t: t = 0 in K (only t-free input is meaningful) */
    }
    for (slong j = 0; j < len && st == GR_SUCCESS; j++) {
        if (!fmpz_mpoly_q_is_zero(acc + j, D->mctx)) {
            for (slong i = 0; i < D->n; i++)
                fmpz_mpoly_q_mul(FF_FNX_COEFF(term, i), FF_FNX_COEFF(tj, i), acc + j, D->mctx);
            st |= gr_add(out, out, term, K);
        }
        if (j + 1 < len) st |= gr_mul(tj, tj, tgen, K);
    }
    gr_heap_clear(tj, K); gr_heap_clear(tgen, K); gr_heap_clear(term, K);
    return st == GR_SUCCESS;
}

static int ffe_read_big(fmpz_mpoly_q_t big, fmpz_mpoly_ctx_t bctx, const Expr* e, const ff_names* nm,
                        slong* dxmax, slong* dtmax) {
    if (!ffe_read(big, e, nm, bctx)) return 0;
    /* the denominator is free of x and t */
    if (fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(big), 0, bctx) > 0 ||
        fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(big), 1, bctx) > 0) return 0;
    *dxmax = fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(big), 0, bctx);
    *dtmax = fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(big), 1, bctx);
    if (*dxmax < 0) *dxmax = 0;
    if (*dtmax < 0) *dtmax = 0;
    return 1;
}

int ff_fnx_poly_from_expr(gr_poly_t P, const Expr* e, const ff_names* nm, gr_ctx_t K) {
    ff_fnx_data* D = FF_FNX_DATA(K);
    if (nm->np != D->nvars) return 0;
    fmpz_mpoly_ctx_t bctx; fmpz_mpoly_ctx_init(bctx, 2 + nm->np, ORD_LEX);
    fmpz_mpoly_q_t big; fmpz_mpoly_q_init(big, bctx);
    slong dx, dt;
    int ok = ffe_read_big(big, bctx, e, nm, &dx, &dt);
    if (ok && !nm->t && dt > 0) ok = 0;
    if (ok) ok = gr_poly_zero(P, K) == GR_SUCCESS;
    if (ok && !fmpz_mpoly_q_is_zero(big, bctx)) {
        slong cells = (dx + 1) * (dt + 1);
        fmpz_mpoly_q_struct* acc = flint_malloc((size_t)cells * sizeof(fmpz_mpoly_q_struct));
        for (slong k = 0; k < cells; k++) fmpz_mpoly_q_init(acc + k, D->mctx);
        ffe_split(acc, dx, dt, big, bctx, nm, D);
        gr_ptr c = gr_heap_init(K);
        for (slong d = 0; d <= dx && ok; d++) {
            ok = ffe_elem_from_tpoly(c, acc + d * (dt + 1), dt + 1, K);
            if (ok && gr_is_zero(c, K) != T_TRUE) ok = gr_poly_set_coeff_scalar(P, d, c, K) == GR_SUCCESS;
        }
        gr_heap_clear(c, K);
        for (slong k = 0; k < cells; k++) fmpz_mpoly_q_clear(acc + k, D->mctx);
        flint_free(acc);
    }
    fmpz_mpoly_q_clear(big, bctx);
    fmpz_mpoly_ctx_clear(bctx);
    return ok;
}

int ff_fnx_elem_from_expr(gr_ptr c, const Expr* e, const ff_names* nm, gr_ctx_t K) {
    ff_names n2 = *nm; n2.x = NULL;
    gr_poly_t P; gr_poly_init(P, K);
    int ok = ff_fnx_poly_from_expr(P, e, &n2, K);
    if (ok) ok = (P->length == 0 ? gr_zero(c, K) : gr_set(c, P->coeffs, K)) == GR_SUCCESS;
    gr_poly_clear(P, K);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  back to Expr                                                       */
/* ------------------------------------------------------------------ */

static Expr* ffe_mono(const char* v, slong k) {
    if (k == 1) return expr_new_symbol(v);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_symbol(v), expr_new_integer((int64_t)k) }, 2);
}

static Expr* ffe_fmpz_to_expr(const fmpz_t z) {
    fmpq_t q; fmpq_init(q);
    fmpz_set(fmpq_numref(q), z); fmpz_one(fmpq_denref(q));
    Expr* r = fb_fmpq_to_expr(q);
    fmpq_clear(q);
    return r;
}

static Expr* ffe_mpoly_to_expr(const fmpz_mpoly_t A, const ff_names* nm, const fmpz_mpoly_ctx_t ctx) {
    slong len = fmpz_mpoly_length(A, ctx);
    if (len == 0) return expr_new_integer(0);
    slong nv = nm->np > 0 ? nm->np : 1;
    ulong* ex = flint_malloc((size_t)nv * sizeof(ulong));
    Expr** terms = malloc((size_t)len * sizeof(Expr*));
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i < len; i++) {
        fmpz_mpoly_get_term_exp_ui(ex, A, i, ctx);
        fmpz_mpoly_get_term_coeff_fmpz(c, A, i, ctx);
        Expr** f = malloc((size_t)(nv + 1) * sizeof(Expr*));
        size_t nf = 0;
        f[nf++] = ffe_fmpz_to_expr(c);
        for (slong j = 0; j < nm->np; j++) if (ex[j] > 0) f[nf++] = ffe_mono(nm->p[j], (slong)ex[j]);
        terms[i] = nf == 1 ? f[0] : expr_new_function(expr_new_symbol(SYM_Times), f, nf);
        free(f);
    }
    fmpz_clear(c);
    flint_free(ex);
    Expr* r = len == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)len);
    free(terms);
    return r;
}

static Expr* ffe_q_to_expr(const fmpz_mpoly_q_t a, const ff_names* nm, const fmpz_mpoly_ctx_t ctx) {
    Expr* num = ffe_mpoly_to_expr(fmpz_mpoly_q_numref(a), nm, ctx);
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(a), ctx)) return num;
    Expr* den = ffe_mpoly_to_expr(fmpz_mpoly_q_denref(a), nm, ctx);
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2);
}

Expr* ff_fnx_elem_to_expr(gr_srcptr c, const ff_names* nm, gr_ctx_t K) {
    ff_fnx_data* D = FF_FNX_DATA(K);
    Expr** terms = malloc((size_t)D->n * sizeof(Expr*));
    size_t nt = 0;
    for (slong j = 0; j < D->n; j++) {
        const fmpz_mpoly_q_struct* a = FF_FNX_COEFF((gr_ptr)c, j);
        if (fmpz_mpoly_q_is_zero(a, D->mctx)) continue;
        Expr* q = ffe_q_to_expr(a, nm, D->mctx);
        terms[nt++] = j == 0 ? q : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ q, ffe_mono(nm->t, j) }, 2);
    }
    Expr* r = nt == 0 ? expr_new_integer(0) : nt == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return r;
}

Expr* ff_fnx_poly_to_expr(const gr_poly_t P, const ff_names* nm, gr_ctx_t K) {
    slong sz = gr_ctx_sizeof_elem(K);
    Expr** terms = malloc((size_t)(P->length > 0 ? P->length : 1) * sizeof(Expr*));
    size_t nt = 0;
    for (slong d = 0; d < P->length; d++) {
        gr_srcptr c = (gr_srcptr)((const char*)P->coeffs + d * sz);
        if (gr_is_zero(c, K) == T_TRUE) continue;
        Expr* e = ff_fnx_elem_to_expr(c, nm, K);
        /* each coefficient of x^d as ONE fraction (Together), the canonical form of
         * Trager's fieldFn (fnCollect): the package compares these forms */
        e = expr_new_function(expr_new_symbol("Together"), (Expr*[]){ e }, 1);
        terms[nt++] = d == 0 ? e : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ e, ffe_mono(nm->x, d) }, 2);
    }
    Expr* r = nt == 0 ? expr_new_integer(0) : nt == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return r;
}

#endif /* USE_FLINT */

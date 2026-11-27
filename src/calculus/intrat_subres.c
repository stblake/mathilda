/* intrat_subres.c — native subresultant PRS over Z[t][x] for the
 * Lazard-Rioboo-Trager log part of the rational-function integrator.
 *
 * Why this exists.  intrat_log_part_core needs, for each squarefree factor
 * Q_i(t) of the Rothstein-Trager resultant Res_x(D, A - t D'), the element of
 * the subresultant chain of D and A - t D' whose x-degree is i (Bronstein,
 * Symbolic Integration I, §2.5, IntRationalLogPart).  It used to take that
 * chain from the SubresultantPolynomialRemainders builtin, which is a plain
 * pseudo-remainder sequence (no content removal, no subresultant scaling) run
 * on Expr trees: its coefficients in t grow exponentially with the length of
 * the chain (Bronstein §1.5).  For Czichowski's example
 *   (7x^13 + 10x^8 + 4x^7 - 7x^6 - 4x^3 - 4x^2 + 3x + 3) /
 *   (x^14 - 2x^8 - 2x^7 - 2x^4 - 4x^3 - x^2 + 2x + 1)
 * the chain had a leaf count of 71924 and took 19 s, against 2 ms for the
 * resultant itself.
 *
 * This module computes the genuine subresultant PRS (algorithm SubResultant,
 * Bronstein p. 24, with the Brown-Collins gamma/beta scaling, every division
 * exact) on dense FLINT arrays: x is the main variable, each x-coefficient an
 * fmpz_poly in t.  The subresultants are determinants of Sylvester
 * sub-matrices, so their size grows only polynomially.  Same chain elements up
 * to a constant factor in Q(t), so the downstream primitive-part and degree
 * reads are unchanged; this is a pure speed change.
 *
 * Scope: exact rational coefficients in x and t only (the pure-rational
 * integrand).  Anything else -- a parameter, a radical, a Root object, the
 * transcendental-monomial callers -- returns NULL and the caller keeps the
 * generic Expr path.
 */
#include "intrat_subres.h"
#include "expand.h"
#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef USE_FLINT
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_mpoly.h>
#include "flint_bridge_internal.h"

/* Above these the generic path is no worse off and a pathological input
 * cannot allocate without bound. */
#define SUBRES_MAX_DEG_X 400
#define SUBRES_MAX_DEG_T 400

/* A polynomial in Z[t][x]: c[i] is the coefficient of x^i, deg = -1 for 0. */
typedef struct {
    int deg;
    int cap;
    fmpz_poly_struct* c;
} ZTX;

static void ztx_init(ZTX* p, int cap) {
    if (cap < 1) cap = 1;
    p->deg = -1;
    p->cap = cap;
    p->c = (fmpz_poly_struct*)malloc(sizeof(fmpz_poly_struct) * (size_t)cap);
    for (int i = 0; i < cap; i++) fmpz_poly_init(p->c + i);
}

static void ztx_clear(ZTX* p) {
    for (int i = 0; i < p->cap; i++) fmpz_poly_clear(p->c + i);
    free(p->c);
    p->c = NULL;
    p->cap = 0;
    p->deg = -1;
}

static void ztx_normalise(ZTX* p) {
    while (p->deg >= 0 && fmpz_poly_is_zero(p->c + p->deg)) p->deg--;
}

static void ztx_copy(ZTX* dst, const ZTX* src) {
    ztx_init(dst, src->deg + 1);
    for (int i = 0; i <= src->deg; i++) fmpz_poly_set(dst->c + i, src->c + i);
    dst->deg = src->deg;
}

/* ------------------------------------------------------------------ */
/* Expr <-> ZTX                                                         */
/* ------------------------------------------------------------------ */

static bool sym_is(const Expr* e, const Expr* s) {
    return e->type == EXPR_SYMBOL && s->type == EXPR_SYMBOL
        && (e->data.symbol.name == s->data.symbol.name
            || strcmp(e->data.symbol.name, s->data.symbol.name) == 0);
}

static bool head_is(const Expr* e, const char* h) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* An exact integer leaf into z. */
static bool int_leaf(const Expr* e, fmpz_t z) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(z, (slong)e->data.integer); return true; }
    if (e->type == EXPR_BIGINT) { fmpz_set_mpz(z, e->data.bigint); return true; }
    return false;
}

/* Multiply q by the exact rational number e; false if e is not one. */
static bool mul_rational(fmpq_t q, const Expr* e) {
    fmpz_t n, d;
    fmpz_init(n); fmpz_init(d);
    bool ok = false;
    if (int_leaf(e, n)) {
        fmpq_mul_fmpz(q, q, n);
        ok = true;
    } else if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2
               && int_leaf(e->data.function.args[0], n)
               && int_leaf(e->data.function.args[1], d) && !fmpz_is_zero(d)) {
        fmpq_t r;
        fmpq_init(r);
        fmpq_set_fmpz_frac(r, n, d);
        fmpq_mul(q, q, r);
        fmpq_clear(r);
        ok = true;
    }
    fmpz_clear(n); fmpz_clear(d);
    return ok;
}

/* One factor of a monomial: a number, x, t, x^k or t^k (k a positive
 * machine integer). */
static bool parse_factor(const Expr* f, const Expr* x, const Expr* t,
                         fmpq_t c, int* ex, int* et) {
    if (sym_is(f, x)) { (*ex)++; return true; }
    if (sym_is(f, t)) { (*et)++; return true; }
    if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
        && f->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t k = f->data.function.args[1]->data.integer;
        const Expr* b = f->data.function.args[0];
        if (k < 1 || k > SUBRES_MAX_DEG_X) return false;
        if (sym_is(b, x)) { *ex += (int)k; return true; }
        if (sym_is(b, t)) { *et += (int)k; return true; }
        return false;
    }
    return mul_rational(c, f);
}

/* One summand c x^i t^j: accumulate into the per-(i, j) rational table. */
typedef struct { int i, j; fmpq_t c; } Mono;

static bool parse_term(const Expr* term, const Expr* x, const Expr* t, Mono* m) {
    fmpq_one(m->c);
    m->i = m->j = 0;
    if (head_is(term, SYM_Times)) {
        for (size_t k = 0; k < term->data.function.arg_count; k++)
            if (!parse_factor(term->data.function.args[k], x, t, m->c, &m->i, &m->j))
                return false;
        return true;
    }
    return parse_factor(term, x, t, m->c, &m->i, &m->j);
}

/* p (any form) -> ZTX scaled by the lcm of its coefficient denominators. */
static bool ztx_from_expr(const Expr* p, const Expr* x, const Expr* t, ZTX* out) {
    Expr* e = expr_expand((Expr*)p);
    if (!e) return false;
    size_t n = head_is(e, SYM_Plus) ? e->data.function.arg_count : 1;
    Mono* ms = (Mono*)malloc(sizeof(Mono) * n);
    size_t nm = 0;
    bool ok = true;
    int maxi = -1;
    for (size_t k = 0; k < n && ok; k++) {
        const Expr* term = head_is(e, SYM_Plus) ? e->data.function.args[k] : e;
        fmpq_init(ms[nm].c);
        nm++;
        ok = parse_term(term, x, t, ms + nm - 1)
             && ms[nm - 1].i <= SUBRES_MAX_DEG_X && ms[nm - 1].j <= SUBRES_MAX_DEG_T;
        if (ok && !fmpq_is_zero(ms[nm - 1].c) && ms[nm - 1].i > maxi) maxi = ms[nm - 1].i;
    }
    if (ok) {
        fmpz_t L, v, cur;
        fmpz_init(L); fmpz_init(v); fmpz_init(cur);
        fmpz_one(L);
        for (size_t k = 0; k < nm; k++) fmpz_lcm(L, L, fmpq_denref(ms[k].c));
        ztx_init(out, maxi + 1);
        for (size_t k = 0; k < nm; k++) {
            if (fmpq_is_zero(ms[k].c)) continue;
            /* c * L is an integer */
            fmpz_divexact(v, L, fmpq_denref(ms[k].c));
            fmpz_mul(v, v, fmpq_numref(ms[k].c));
            fmpz_poly_get_coeff_fmpz(cur, out->c + ms[k].i, ms[k].j);
            fmpz_add(cur, cur, v);
            fmpz_poly_set_coeff_fmpz(out->c + ms[k].i, ms[k].j, cur);
        }
        out->deg = maxi;
        ztx_normalise(out);
        fmpz_clear(L); fmpz_clear(v); fmpz_clear(cur);
    }
    for (size_t k = 0; k < nm; k++) fmpq_clear(ms[k].c);
    free(ms);
    expr_free(e);
    return ok;
}

static Expr* expr_from_fmpz(const fmpz_t z) {
    if (fmpz_fits_si(z)) return expr_new_integer((int64_t)fmpz_get_si(z));
    mpz_t m;
    mpz_init(m);
    fmpz_get_mpz(m, z);
    return expr_new_bigint_take(m);
}

/* x^k as an Expr (1 for k = 0). */
static Expr* expr_pow(const Expr* v, slong k) {
    if (k == 1) return expr_copy((Expr*)v);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)v), expr_new_integer((int64_t)k) }, 2);
}

/* ZTX -> an expanded polynomial in x and t, put in normal form by the
 * evaluator once. */
static Expr* ztx_to_expr(const ZTX* p, const Expr* x, const Expr* t) {
    size_t n = 0;
    for (int i = 0; i <= p->deg; i++)
        for (slong j = 0; j < fmpz_poly_length(p->c + i); j++)
            if (!fmpz_is_zero(p->c[i].coeffs + j)) n++;
    if (n == 0) return expr_new_integer(0);
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    size_t k = 0;
    for (int i = 0; i <= p->deg; i++) {
        for (slong j = 0; j < fmpz_poly_length(p->c + i); j++) {
            const fmpz* cz = p->c[i].coeffs + j;
            if (fmpz_is_zero(cz)) continue;
            Expr* fs[3];
            size_t nf = 0;
            fs[nf++] = expr_from_fmpz(cz);
            if (j > 0) fs[nf++] = expr_pow(t, j);
            if (i > 0) fs[nf++] = expr_pow(x, i);
            terms[k++] = nf == 1 ? fs[0]
                       : expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
        }
    }
    Expr* sum = n == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
    free(terms);
    return eval_and_free(sum);
}

/* ------------------------------------------------------------------ */
/* Arithmetic in Z[t][x]                                               */
/* ------------------------------------------------------------------ */

/* R := prem(A, B) = lc(B)^(deg A - deg B + 1) A mod B.  B != 0. */
static void ztx_prem(ZTX* R, const ZTX* A, const ZTX* B) {
    int m = A->deg, n = B->deg;
    const fmpz_poly_struct* b = B->c + n;
    ztx_copy(R, A);
    int e = m - n + 1;
    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);
    while (R->deg >= n) {
        int k = R->deg - n;
        fmpz_poly_t lr;
        fmpz_poly_init(lr);
        fmpz_poly_set(lr, R->c + R->deg);
        for (int idx = 0; idx <= R->deg; idx++) fmpz_poly_mul(R->c + idx, R->c + idx, b);
        for (int idx = 0; idx <= n; idx++) {
            fmpz_poly_mul(tmp, lr, B->c + idx);
            fmpz_poly_sub(R->c + idx + k, R->c + idx + k, tmp);
        }
        fmpz_poly_clear(lr);
        /* the leading term cancels exactly */
        fmpz_poly_zero(R->c + R->deg);
        ztx_normalise(R);
        e--;
    }
    if (e > 0 && R->deg >= 0) {
        fmpz_poly_pow(tmp, b, (ulong)e);
        for (int idx = 0; idx <= R->deg; idx++) fmpz_poly_mul(R->c + idx, R->c + idx, tmp);
    }
    fmpz_poly_clear(tmp);
}

/* p := p / q coefficient-wise, exact; false if some division is inexact. */
static bool ztx_divexact(ZTX* p, const fmpz_poly_t q) {
    fmpz_poly_t quo;
    fmpz_poly_init(quo);
    bool ok = true;
    for (int i = 0; i <= p->deg && ok; i++) {
        if (fmpz_poly_is_zero(p->c + i)) continue;
        ok = fmpz_poly_divides(quo, p->c + i, q) != 0;
        if (ok) fmpz_poly_swap(p->c + i, quo);
    }
    fmpz_poly_clear(quo);
    return ok;
}

Expr* intrat_subres_chain(const Expr* A, const Expr* B, const Expr* x, const Expr* t) {
    if (!A || !B || !x || !t || x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    ZTX a, b;
    if (!ztx_from_expr(A, x, t, &a)) return NULL;
    if (!ztx_from_expr(B, x, t, &b)) { ztx_clear(&a); return NULL; }
    if (a.deg < b.deg) { ZTX s = a; a = b; b = s; }
    if (b.deg < 1) { ztx_clear(&a); ztx_clear(&b); return NULL; }

    int cap = a.deg + 3, k = 0;
    ZTX* R = (ZTX*)malloc(sizeof(ZTX) * (size_t)cap);
    R[k++] = a;
    R[k++] = b;

    /* SubResultant (Bronstein p. 24):  gamma_1 = -1, delta_1 = deg A - deg B,
     * beta_1 = (-1)^(delta_1 + 1);  R_{i+1} = prem(R_{i-1}, R_i) / beta_i;
     * gamma_{i+1} = (-r_i)^delta_i gamma_i^(1 - delta_i);
     * delta_{i+1} = deg R_i - deg R_{i+1};  beta_{i+1} = -r_i gamma_{i+1}^delta_{i+1}. */
    fmpz_poly_t gamma, beta, r, tmp, tmp2;
    fmpz_poly_init(gamma); fmpz_poly_init(beta); fmpz_poly_init(r);
    fmpz_poly_init(tmp); fmpz_poly_init(tmp2);
    fmpz_poly_set_si(gamma, -1);
    int delta = a.deg - b.deg;
    fmpz_poly_set_si(beta, (delta + 1) % 2 == 0 ? 1 : -1);
    bool ok = true;
    while (ok && R[k - 1].deg >= 0) {
        fmpz_poly_set(r, R[k - 1].c + R[k - 1].deg);
        ZTX nxt;
        ztx_prem(&nxt, R + k - 2, R + k - 1);
        if (nxt.deg < 0) { ztx_clear(&nxt); break; }
        if (!ztx_divexact(&nxt, beta)) { ztx_clear(&nxt); ok = false; break; }
        if (k >= cap) { cap *= 2; R = (ZTX*)realloc(R, sizeof(ZTX) * (size_t)cap); }
        R[k++] = nxt;
        /* gamma := (-r)^delta * gamma^(1 - delta) */
        fmpz_poly_neg(tmp, r);
        if (delta == 0) {
            /* gamma unchanged */
        } else if (delta == 1) {
            fmpz_poly_set(gamma, tmp);
        } else {
            fmpz_poly_pow(tmp2, tmp, (ulong)delta);
            fmpz_poly_pow(tmp, gamma, (ulong)(delta - 1));
            if (!fmpz_poly_divides(gamma, tmp2, tmp)) { ok = false; break; }
        }
        delta = R[k - 2].deg - R[k - 1].deg;
        /* beta := -r gamma^delta */
        fmpz_poly_pow(tmp, gamma, (ulong)delta);
        fmpz_poly_mul(beta, r, tmp);
        fmpz_poly_neg(beta, beta);
    }
    fmpz_poly_clear(gamma); fmpz_poly_clear(beta); fmpz_poly_clear(r);
    fmpz_poly_clear(tmp); fmpz_poly_clear(tmp2);

    Expr* out = NULL;
    if (ok) {
        Expr** els = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
        for (int i = 0; i < k; i++) els[i] = ztx_to_expr(R + i, x, t);
        out = expr_new_function(expr_new_symbol(SYM_List), els, (size_t)k);
        free(els);
    }
    for (int i = 0; i < k; i++) ztx_clear(R + i);
    free(R);
    return out;
}

/* ------------------------------------------------------------------ */
/* The same chain over Q[t, y1..ym][x]                                 */
/* ------------------------------------------------------------------ */

/* Parameters (Q(k)) and algebraic mode (a primitive element z, every
 * coefficient a polynomial in z reduced afterwards by the caller) put further
 * symbols beside t.  The subresultants are polynomials in the coefficients, so
 * they may be computed over Q[t, y][x] and reduced (modulo the minimal
 * polynomial of z) at the end: reduction is a ring homomorphism.  Each
 * x-coefficient is an fmpq_mpoly over every symbol of A and B except x. */
typedef struct { int deg, cap; fmpq_mpoly_struct* c; } MTX;

#define SUBRES_MV_MAX_TERMS 200000

static void mtx_init(MTX* p, int cap, const fmpq_mpoly_ctx_struct* ctx) {
    if (cap < 1) cap = 1;
    p->deg = -1; p->cap = cap;
    p->c = (fmpq_mpoly_struct*)malloc(sizeof(fmpq_mpoly_struct) * (size_t)cap);
    for (int i = 0; i < cap; i++) fmpq_mpoly_init(p->c + i, ctx);
}
static void mtx_clear(MTX* p, const fmpq_mpoly_ctx_struct* ctx) {
    for (int i = 0; i < p->cap; i++) fmpq_mpoly_clear(p->c + i, ctx);
    free(p->c); p->c = NULL; p->cap = 0; p->deg = -1;
}
static void mtx_normalise(MTX* p, const fmpq_mpoly_ctx_struct* ctx) {
    while (p->deg >= 0 && fmpq_mpoly_is_zero(p->c + p->deg, ctx)) p->deg--;
}
static void mtx_copy(MTX* dst, const MTX* src, const fmpq_mpoly_ctx_struct* ctx) {
    mtx_init(dst, src->deg + 1, ctx);
    for (int i = 0; i <= src->deg; i++) fmpq_mpoly_set(dst->c + i, src->c + i, ctx);
    dst->deg = src->deg;
}

/* P (an fmpq_mpoly over vs, x at index xi) split by the degree in x. */
static bool mtx_from_mpoly(const fmpq_mpoly_t P, slong xi, slong nv,
                           const fmpq_mpoly_ctx_struct* ctx, MTX* out) {
    slong len = fmpq_mpoly_length(P, ctx);
    ulong* e = (ulong*)malloc(sizeof(ulong) * (size_t)nv);
    int maxi = -1;
    for (slong k = 0; k < len; k++) {
        fmpq_mpoly_get_term_exp_ui(e, P, k, ctx);
        if ((int)e[xi] > maxi) maxi = (int)e[xi];
    }
    if (maxi > SUBRES_MAX_DEG_X) { free(e); return false; }
    mtx_init(out, maxi + 1, ctx);
    fmpq_t c, old; fmpq_init(c); fmpq_init(old);
    for (slong k = 0; k < len; k++) {
        fmpq_mpoly_get_term_exp_ui(e, P, k, ctx);
        fmpq_mpoly_get_term_coeff_fmpq(c, P, k, ctx);
        int i = (int)e[xi];
        e[xi] = 0;
        fmpq_mpoly_get_coeff_fmpq_ui(old, out->c + i, e, ctx);
        fmpq_add(old, old, c);
        fmpq_mpoly_set_coeff_fmpq_ui(out->c + i, old, e, ctx);
    }
    fmpq_clear(c); fmpq_clear(old);
    free(e);
    out->deg = maxi;
    mtx_normalise(out, ctx);
    return true;
}

static Expr* mtx_to_expr(const MTX* p, slong xi, slong nv, const VarSet* vs,
                         const fmpq_mpoly_ctx_struct* ctx) {
    size_t n = 0;
    for (int i = 0; i <= p->deg; i++) n += (size_t)fmpq_mpoly_length(p->c + i, ctx);
    if (n == 0) return expr_new_integer(0);
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    ulong* e = (ulong*)malloc(sizeof(ulong) * (size_t)nv);
    Expr** fs = (Expr**)malloc(sizeof(Expr*) * (size_t)(nv + 2));
    fmpq_t c; fmpq_init(c);
    size_t k = 0;
    for (int i = 0; i <= p->deg; i++) {
        for (slong j = 0; j < fmpq_mpoly_length(p->c + i, ctx); j++) {
            fmpq_mpoly_get_term_exp_ui(e, p->c + i, j, ctx);
            fmpq_mpoly_get_term_coeff_fmpq(c, p->c + i, j, ctx);
            e[xi] = (ulong)i;
            size_t nf = 0;
            fs[nf++] = fb_fmpq_to_expr(c);
            for (slong v = 0; v < nv; v++) {
                if (!e[v]) continue;
                Expr* sv = expr_new_symbol(vs->names[v]);
                if (e[v] == 1) fs[nf++] = sv;
                else {
                    fs[nf++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ sv, expr_new_integer((int64_t)e[v]) }, 2);
                }
            }
            terms[k++] = nf == 1 ? fs[0] : expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
        }
    }
    fmpq_clear(c);
    free(e); free(fs);
    Expr* sum = n == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
    free(terms);
    return eval_and_free(sum);
}

static void mtx_prem(MTX* R, const MTX* A, const MTX* B, const fmpq_mpoly_ctx_struct* ctx) {
    int m = A->deg, n = B->deg;
    const fmpq_mpoly_struct* b = B->c + n;
    mtx_copy(R, A, ctx);
    int e = m - n + 1;
    fmpq_mpoly_t tmp, lr;
    fmpq_mpoly_init(tmp, ctx); fmpq_mpoly_init(lr, ctx);
    while (R->deg >= n) {
        int k = R->deg - n;
        fmpq_mpoly_set(lr, R->c + R->deg, ctx);
        for (int idx = 0; idx <= R->deg; idx++) fmpq_mpoly_mul(R->c + idx, R->c + idx, b, ctx);
        for (int idx = 0; idx <= n; idx++) {
            fmpq_mpoly_mul(tmp, lr, B->c + idx, ctx);
            fmpq_mpoly_sub(R->c + idx + k, R->c + idx + k, tmp, ctx);
        }
        fmpq_mpoly_zero(R->c + R->deg, ctx);
        mtx_normalise(R, ctx);
        e--;
    }
    if (e > 0 && R->deg >= 0) {
        fmpq_mpoly_pow_ui(tmp, b, (ulong)e, ctx);
        for (int idx = 0; idx <= R->deg; idx++) fmpq_mpoly_mul(R->c + idx, R->c + idx, tmp, ctx);
    }
    fmpq_mpoly_clear(tmp, ctx); fmpq_mpoly_clear(lr, ctx);
}

static bool mtx_divexact(MTX* p, const fmpq_mpoly_t q, const fmpq_mpoly_ctx_struct* ctx) {
    fmpq_mpoly_t quo; fmpq_mpoly_init(quo, ctx);
    bool ok = true;
    for (int i = 0; i <= p->deg && ok; i++) {
        if (fmpq_mpoly_is_zero(p->c + i, ctx)) continue;
        ok = fmpq_mpoly_divides(quo, p->c + i, q, ctx) != 0;
        if (ok) fmpq_mpoly_swap(p->c + i, quo, ctx);
        if (ok && fmpq_mpoly_length(p->c + i, ctx) > SUBRES_MV_MAX_TERMS) ok = false;
    }
    fmpq_mpoly_clear(quo, ctx);
    return ok;
}

Expr* intrat_subres_chain_mv(const Expr* A, const Expr* B, const Expr* x, const Expr* t) {
    if (!A || !B || !x || !t || x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    Expr* Ae = expr_expand((Expr*)A);
    Expr* Be = expr_expand((Expr*)B);
    if (!Ae || !Be) { if (Ae) expr_free(Ae); if (Be) expr_free(Be); return NULL; }
    VarSet vs; memset(&vs, 0, sizeof vs);
    Expr* out = NULL;
    MTX* R = NULL;
    int k = 0;
    const fmpq_mpoly_ctx_struct* ctx = NULL;
    slong nv = 0, xi = -1;
    if (!fb_collect_vars(Ae, &vs) || !fb_collect_vars(Be, &vs) ||
        !fb_collect_vars(x, &vs) || !fb_collect_vars(t, &vs)) goto done;
    nv = (slong)vs.count;
    for (slong v = 0; v < nv; v++) if (strcmp(vs.names[v], x->data.symbol.name) == 0) xi = v;
    if (xi < 0 || nv > 16) goto done;
    ctx = fb_ctx_lex(nv);
    if (!ctx) goto done;
    {
        fmpq_mpoly_t pa, pb;
        fmpq_mpoly_init(pa, ctx); fmpq_mpoly_init(pb, ctx);
        MTX a, b;
        bool ok = fb_to_mpoly(Ae, pa, ctx, &vs) && fb_to_mpoly(Be, pb, ctx, &vs);
        bool ha = false, hb = false;
        if (ok) { ok = mtx_from_mpoly(pa, xi, nv, ctx, &a); ha = ok; }
        if (ok) { ok = mtx_from_mpoly(pb, xi, nv, ctx, &b); hb = ok; }
        fmpq_mpoly_clear(pa, ctx); fmpq_mpoly_clear(pb, ctx);
        if (!ok) {
            if (ha) mtx_clear(&a, ctx);
            if (hb) mtx_clear(&b, ctx);
            goto done;
        }
        if (a.deg < b.deg) { MTX s0 = a; a = b; b = s0; }
        if (b.deg < 1) { mtx_clear(&a, ctx); mtx_clear(&b, ctx); goto done; }
        int cap = a.deg + 3;
        R = (MTX*)malloc(sizeof(MTX) * (size_t)cap);
        R[k++] = a;
        R[k++] = b;
        fmpq_mpoly_t gamma, beta, r, tmp, tmp2;
        fmpq_mpoly_init(gamma, ctx); fmpq_mpoly_init(beta, ctx); fmpq_mpoly_init(r, ctx);
        fmpq_mpoly_init(tmp, ctx); fmpq_mpoly_init(tmp2, ctx);
        fmpq_mpoly_set_si(gamma, -1, ctx);
        int delta = a.deg - b.deg;
        fmpq_mpoly_set_si(beta, (delta + 1) % 2 == 0 ? 1 : -1, ctx);
        while (ok && R[k - 1].deg >= 0) {
            fmpq_mpoly_set(r, R[k - 1].c + R[k - 1].deg, ctx);
            MTX nxt;
            mtx_prem(&nxt, R + k - 2, R + k - 1, ctx);
            if (nxt.deg < 0) { mtx_clear(&nxt, ctx); break; }
            if (!mtx_divexact(&nxt, beta, ctx)) { mtx_clear(&nxt, ctx); ok = false; break; }
            if (k >= cap) { cap *= 2; R = (MTX*)realloc(R, sizeof(MTX) * (size_t)cap); }
            R[k++] = nxt;
            fmpq_mpoly_neg(tmp, r, ctx);
            if (delta == 0) {
                /* gamma unchanged */
            } else if (delta == 1) {
                fmpq_mpoly_set(gamma, tmp, ctx);
            } else {
                fmpq_mpoly_pow_ui(tmp2, tmp, (ulong)delta, ctx);
                fmpq_mpoly_pow_ui(tmp, gamma, (ulong)(delta - 1), ctx);
                if (!fmpq_mpoly_divides(gamma, tmp2, tmp, ctx)) { ok = false; break; }
            }
            delta = R[k - 2].deg - R[k - 1].deg;
            fmpq_mpoly_pow_ui(tmp, gamma, (ulong)delta, ctx);
            fmpq_mpoly_mul(beta, r, tmp, ctx);
            fmpq_mpoly_neg(beta, beta, ctx);
        }
        fmpq_mpoly_clear(gamma, ctx); fmpq_mpoly_clear(beta, ctx); fmpq_mpoly_clear(r, ctx);
        fmpq_mpoly_clear(tmp, ctx); fmpq_mpoly_clear(tmp2, ctx);
        if (ok) {
            Expr** els = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
            for (int i = 0; i < k; i++) els[i] = mtx_to_expr(R + i, xi, nv, &vs, ctx);
            out = expr_new_function(expr_new_symbol(SYM_List), els, (size_t)k);
            free(els);
        }
    }
done:
    if (R) {
        for (int i = 0; i < k; i++) mtx_clear(R + i, ctx);
        free(R);
    }
    fb_varset_free(&vs);
    expr_free(Ae); expr_free(Be);
    return out;
}

#else  /* !USE_FLINT */

Expr* intrat_subres_chain(const Expr* A, const Expr* B, const Expr* x, const Expr* t) {
    (void)A; (void)B; (void)x; (void)t;
    return NULL;
}
Expr* intrat_subres_chain_mv(const Expr* A, const Expr* B, const Expr* x, const Expr* t) {
    (void)A; (void)B; (void)x; (void)t;
    return NULL;
}

#endif

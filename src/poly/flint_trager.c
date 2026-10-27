/*
 * flint_trager.c
 * --------------
 * Factoring over an algebraic number field K, in any number of variables
 * (see flint_trager.h for the contract).
 *
 * THE METHOD (Trager 1976, "Algebraic factoring and rational function
 * integration").  K = Q(theta) with theta an algebraic integer of minimal
 * polynomial M(tau), found as the common field of every algebraic constant in
 * play (the extension generators and the input's own coefficients) by the same
 * primitive-element normalisation flint_field_gcd uses.  A polynomial A over K
 * is lifted to Q[x_1..x_k, tau] reduced mod M -- the "tau lift" -- and made
 * monic.  Then:
 *
 *   1. shift   A_s(x) = A(x_1 - s_1 tau, ..., x_k - s_k tau)   (s = 0 first,
 *              when A has algebraic coefficients at all);
 *   2. norm    N(x) = Res_tau(M(tau), A_s(x, tau)), a polynomial over Q;
 *   3. factor  N = prod G_j^E_j over Q (fmpq_mpoly_factor);
 *   4. lift    H_j = gcd_K(A_s, G_j) by the certified modular field GCD, and
 *              unshift.
 *
 * Every G_j is squarefree, hence so is H_j, and if H_j were reducible over K
 * each of its irreducible factors would have a norm that is a power of G_j, so
 * n deg(H_j) would be at least 2 deg(G_j).  Hence the test
 *
 *        [K:Q] * deg(H_j) == deg(G_j)
 *
 * proves H_j irreducible with norm exactly G_j -- and then no other factor of
 * A can divide G_j, so H_j's multiplicity in A is exactly E_j.  Every H_j that
 * passes is divided out of A (exactly, or the call declines); a G_j that fails
 * (two conjugate factors sharing a norm, or a factor whose coefficients lie in
 * a proper subfield) leaves its part of A for the next shift, on the cofactor
 * only.  So no squarefree pre-pass over K is needed at all: the review of the
 * univariate engine found 75% of its time in exactly that pre-pass.
 *
 * The answer is CERTIFIED at the end -- c f_1^e_1 ... f_k^e_k is recomputed in
 * Q[x, tau]/(M) and compared with the input -- so no unlucky shift, prime or
 * normalisation can produce a wrong factorisation; at worst the call declines.
 *
 * OUTPUT CONVENTIONS.  Each factor is written over the product basis of the
 * caller's own radicals (Sqrt[2] stays Sqrt[2]; Q(Sqrt[2], Sqrt[3]) is spelled
 * in Sqrt[2], Sqrt[3], Sqrt[6]), made primitive over Z in that basis, and
 * sign-normalised so that its leading term is positive, where the generators
 * rank above the variables (x^2 - 2 over Q(Sqrt[2]) is
 * -((Sqrt[2] - x)(Sqrt[2] + x))) and the variables rank in their sorted order.
 * I is a coefficient, not a generator: over Q(i) a factor is primitive over
 * Z[i] and its leading coefficient is the associate in the first quadrant
 * ((x - I y)(x + I y)); the Gaussian content goes in front.
 */
#include "flint_trager.h"

#include <stdlib.h>
#include <string.h>

#include "eval.h"
#include "sym_names.h"
#include "flint_qqbar.h"

#ifdef USE_FLINT

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>

#include "flint_bridge_internal.h"

/* How many shifts to try before declining.  Generically the first or second
 * works; the cap only bounds the pathological case. */
#define TRAGER_MAX_ATTEMPTS 12

/* ------------------------------------------------------------------ */
/*  Small Expr helpers                                                 */
/* ------------------------------------------------------------------ */

static const char* tr_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    return (h && h->type == EXPR_SYMBOL) ? h->data.symbol.name : NULL;
}

static int tr_is_exact_rational(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return 1;
    return tr_head(e) == SYM_Rational && e->data.function.arg_count == 2;
}

static int tr_is_int(const Expr* e, int64_t v) {
    return e->type == EXPR_INTEGER && e->data.integer == v;
}

static int tr_is_imag_unit(const Expr* e) {
    return tr_head(e) == SYM_Complex && e->data.function.arg_count == 2 &&
           tr_is_int(e->data.function.args[0], 0) && tr_is_int(e->data.function.args[1], 1);
}

int flint_trager_has_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    if (tr_head(e) == SYM_Complex && e->data.function.arg_count == 2)
        return tr_is_exact_rational(e->data.function.args[0]) &&
               tr_is_exact_rational(e->data.function.args[1]) &&
               !tr_is_int(e->data.function.args[1], 0);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (flint_trager_has_complex(e->data.function.args[i])) return 1;
    return 0;
}

/* Complex[a, b] -> a + b I, with I the atom Complex[0, 1] (left unevaluated),
 * so a Gaussian coefficient such as 2 I or 1 - I is a polynomial in the one
 * algebraic atom I instead of an atom of its own (review item F1: Complex[0, 2]
 * had been an independent generator, and Extension -> I left it alone). */
static Expr* tr_split_complex(const Expr* e) {
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_I)   /* unevaluated I */
        return expr_new_function(expr_new_symbol(SYM_Complex),
                                 (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    if (tr_head(e) == SYM_Complex && n == 2 && !tr_is_imag_unit(e) &&
        tr_is_exact_rational(e->data.function.args[0]) &&
        tr_is_exact_rational(e->data.function.args[1])) {
        const Expr* re = e->data.function.args[0];
        const Expr* im = e->data.function.args[1];
        Expr* iu = expr_new_function(expr_new_symbol(SYM_Complex),
            (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
        Expr* t = tr_is_int(im, 1) ? iu
            : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_copy((Expr*)im), iu }, 2);
        if (tr_is_int(re, 0)) return t;
        return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_copy((Expr*)re), t }, 2);
    }
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = tr_split_complex(e->data.function.args[i]);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}

/* ------------------------------------------------------------------ */
/*  The field context                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    const fmpq_mpoly_ctx_struct* ctx;
    slong nvars, ngens, tau_idx, n;
    fmpq_mpoly_t M;                 /* theta's minimal polynomial in tau */
    const Expr* mp;                 /* its integer coefficient List      */
    nf_t nf;                        /* K, for element inverses           */
} TCtx;

/* A reduced mod M, in place. */
static void tr_reduce(fmpq_mpoly_t A, TCtx* T) {
    fmpq_mpoly_t q, r;
    fmpq_mpoly_init(q, T->ctx); fmpq_mpoly_init(r, T->ctx);
    fmpq_mpoly_divrem(q, r, A, T->M, T->ctx);
    fmpq_mpoly_swap(A, r, T->ctx);
    fmpq_mpoly_clear(q, T->ctx); fmpq_mpoly_clear(r, T->ctx);
}

/* Total degree in the polynomial variables (tau excluded); -1 for zero. */
static slong tr_xdeg(const fmpq_mpoly_t A, TCtx* T) {
    slong len = fmpq_mpoly_length(A, T->ctx), best = -1;
    ulong* e = flint_malloc(sizeof(ulong) * (size_t)T->nvars);
    for (slong i = 0; i < len; i++) {
        fmpq_mpoly_get_term_exp_ui(e, A, i, T->ctx);
        slong d = 0;
        for (slong v = 0; v < T->ngens; v++) d += (slong)e[v];
        if (d > best) best = d;
    }
    flint_free(e);
    return best;
}

static int tr_has_tau(const fmpq_mpoly_t A, TCtx* T) {
    return fmpq_mpoly_degree_si(A, T->tau_idx, T->ctx) > 0;
}

/* The K-coefficient of A's LEX-leading monomial in the variables. */
static void tr_lead_kelem(fmpq_poly_t c, const fmpq_mpoly_t A, TCtx* T) {
    fmpq_poly_zero(c);
    slong len = fmpq_mpoly_length(A, T->ctx);
    if (len == 0) return;
    ulong* e0 = flint_malloc(sizeof(ulong) * (size_t)T->nvars);
    ulong* e  = flint_malloc(sizeof(ulong) * (size_t)T->nvars);
    fmpq_t q; fmpq_init(q);
    fmpq_mpoly_get_term_exp_ui(e0, A, 0, T->ctx);
    for (slong i = 0; i < len; i++) {
        fmpq_mpoly_get_term_exp_ui(e, A, i, T->ctx);
        int same = 1;
        for (slong v = 0; v < T->ngens; v++) if (e[v] != e0[v]) { same = 0; break; }
        if (!same) break;               /* tau is last in LEX: one run */
        fmpq_mpoly_get_term_coeff_fmpq(q, A, i, T->ctx);
        fmpq_poly_set_coeff_fmpq(c, (slong)e[T->tau_idx], q);
    }
    fmpq_clear(q); flint_free(e0); flint_free(e);
}

/* A K-element (polynomial in tau) as a constant of Q[x, tau]. */
static void tr_kelem_to_mpoly(fmpq_mpoly_t out, const fmpq_poly_t c, TCtx* T) {
    fmpq_mpoly_zero(out, T->ctx);
    ulong* e = flint_calloc((size_t)T->nvars, sizeof(ulong));
    fmpq_t q; fmpq_init(q);
    for (slong j = 0; j <= fmpq_poly_degree(c); j++) {
        fmpq_poly_get_coeff_fmpq(q, c, j);
        if (fmpq_is_zero(q)) continue;
        e[T->tau_idx] = (ulong)j;
        fmpq_mpoly_set_coeff_fmpq_ui(out, q, e, T->ctx);
    }
    fmpq_clear(q); flint_free(e);
}

/* A <- c A mod M. */
static void tr_scale(fmpq_mpoly_t A, const fmpq_poly_t c, TCtx* T) {
    fmpq_mpoly_t C; fmpq_mpoly_init(C, T->ctx);
    tr_kelem_to_mpoly(C, c, T);
    fmpq_mpoly_mul(A, A, C, T->ctx);
    tr_reduce(A, T);
    fmpq_mpoly_clear(C, T->ctx);
}

/* c <- 1/c in K. */
static void tr_kinv(fmpq_poly_t c, TCtx* T) {
    nf_elem_t a; nf_elem_init(a, T->nf);
    nf_elem_set_fmpq_poly(a, c, T->nf);
    nf_elem_inv(a, a, T->nf);
    nf_elem_get_fmpq_poly(c, a, T->nf);
    nf_elem_clear(a, T->nf);
}

/* c <- a b in K. */
static void tr_kmul(fmpq_poly_t c, const fmpq_poly_t a, const fmpq_poly_t b, TCtx* T) {
    nf_elem_t x, y; nf_elem_init(x, T->nf); nf_elem_init(y, T->nf);
    nf_elem_set_fmpq_poly(x, a, T->nf);
    nf_elem_set_fmpq_poly(y, b, T->nf);
    nf_elem_mul(x, x, y, T->nf);
    nf_elem_get_fmpq_poly(c, x, T->nf);
    nf_elem_clear(x, T->nf); nf_elem_clear(y, T->nf);
}

/* out <- A(x_i + sign s_i tau) mod M. */
static int tr_shift(fmpq_mpoly_t out, const fmpq_mpoly_t A, const slong* s, int sign, TCtx* T) {
    fmpq_mpoly_struct* C = flint_malloc(sizeof(fmpq_mpoly_struct) * (size_t)T->nvars);
    fmpq_mpoly_struct** Cp = flint_malloc(sizeof(fmpq_mpoly_struct*) * (size_t)T->nvars);
    fmpq_mpoly_t t; fmpq_mpoly_init(t, T->ctx);
    for (slong v = 0; v < T->nvars; v++) {
        fmpq_mpoly_init(C + v, T->ctx);
        fmpq_mpoly_gen(C + v, v, T->ctx);
        if (v < T->ngens && s[v] != 0) {
            fmpq_mpoly_gen(t, T->tau_idx, T->ctx);
            fmpq_mpoly_scalar_mul_si(t, t, sign * s[v], T->ctx);
            fmpq_mpoly_add(C + v, C + v, t, T->ctx);
        }
        Cp[v] = C + v;
    }
    int ok = fmpq_mpoly_compose_fmpq_mpoly(out, A, Cp, T->ctx, T->ctx);
    if (ok) tr_reduce(out, T);
    for (slong v = 0; v < T->nvars; v++) fmpq_mpoly_clear(C + v, T->ctx);
    fmpq_mpoly_clear(t, T->ctx);
    flint_free(C); flint_free(Cp);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  The factor list                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    fmpq_mpoly_struct* f;
    slong* e;
    size_t num, cap;
} TFac;

static void tfac_push(TFac* L, const fmpq_mpoly_t f, slong e, TCtx* T) {
    if (L->num == L->cap) {
        size_t nc = L->cap ? 2 * L->cap : 8;
        L->f = realloc(L->f, sizeof(fmpq_mpoly_struct) * nc);
        L->e = realloc(L->e, sizeof(slong) * nc);
        L->cap = nc;
    }
    fmpq_mpoly_init(L->f + L->num, T->ctx);
    fmpq_mpoly_set(L->f + L->num, f, T->ctx);
    L->e[L->num] = e;
    L->num++;
}

static void tfac_clear(TFac* L, TCtx* T) {
    for (size_t i = 0; i < L->num; i++) fmpq_mpoly_clear(L->f + i, T->ctx);
    free(L->f); free(L->e);
    memset(L, 0, sizeof *L);
}

/* The shift of attempt k >= 1: small nonzero integers, varied per variable
 * and per attempt.  Every variable is shifted, so a factor free of some
 * variable (y, or y - Sqrt[2] z) still gets a squarefree norm. */
static void tr_pick_shift(slong* s, slong ngens, int k) {
    for (slong v = 0; v < ngens; v++) {
        slong m = (slong)k + 1;
        slong a = 1 + (slong)(((ulong)k * 7u + (ulong)v * 5u + (ulong)(v * v)) % (ulong)m);
        s[v] = ((k + v) & 1) ? -a : a;
    }
}

/* Split the monic A (tau-lifted, reduced) into irreducible factors over K.
 * A is consumed (left as the unit 1).  1 on success, 0 to decline. */
static int tr_split(fmpq_mpoly_t A, TCtx* T, TFac* out) {
    slong* s = flint_calloc((size_t)(T->ngens ? T->ngens : 1), sizeof(slong));
    fmpq_mpoly_t As, N, H, Hs, Q;
    fmpq_mpoly_init(As, T->ctx); fmpq_mpoly_init(N, T->ctx);
    fmpq_mpoly_init(H, T->ctx); fmpq_mpoly_init(Hs, T->ctx); fmpq_mpoly_init(Q, T->ctx);
    int ok = 1;
    /* A with rational coefficients only (x^2 - 2 y^2 over Q(Sqrt[2])) has
     * the norm A^n unshifted, which can never pass: skip straight to a shift. */
    int attempt = tr_has_tau(A, T) ? 0 : 1;

    while (ok && tr_xdeg(A, T) > 0) {
        if (attempt > TRAGER_MAX_ATTEMPTS) { ok = 0; break; }
        if (attempt == 0) memset(s, 0, sizeof(slong) * (size_t)T->ngens);
        else tr_pick_shift(s, T->ngens, attempt);
        attempt++;

        if (!tr_shift(As, A, s, -1, T)) continue;
        if (!fmpq_mpoly_resultant(N, T->M, As, T->tau_idx, T->ctx)) continue;
        if (fmpq_mpoly_is_zero(N, T->ctx)) continue;

        fmpq_mpoly_factor_t fac;
        fmpq_mpoly_factor_init(fac, T->ctx);
        if (!fmpq_mpoly_factor(fac, N, T->ctx)) { fmpq_mpoly_factor_clear(fac, T->ctx); continue; }

        for (slong j = 0; j < fac->num && ok; j++) {
            fmpq_mpoly_struct* G = fac->poly + j;
            slong dG = tr_xdeg(G, T);
            if (dG <= 0) continue;
            if (!fb_field_mpoly_gcd(Hs, As, G, T->M, T->ctx, T->mp, T->ngens, T->tau_idx, T->n))
                continue;
            if (T->n * tr_xdeg(Hs, T) != dG) continue;      /* not irreducible yet */
            if (!tr_shift(H, Hs, s, +1, T)) continue;
            slong E = fmpz_get_si(fac->exp + j);
            for (slong k = 0; k < E && ok; k++) {
                if (!fb_field_mpoly_divides(Q, A, H, T->M, T->ctx)) ok = 0;
                else fmpq_mpoly_swap(A, Q, T->ctx);
            }
            if (ok) tfac_push(out, H, E, T);
        }
        fmpq_mpoly_factor_clear(fac, T->ctx);
    }

    fmpq_mpoly_clear(As, T->ctx); fmpq_mpoly_clear(N, T->ctx);
    fmpq_mpoly_clear(H, T->ctx); fmpq_mpoly_clear(Hs, T->ctx); fmpq_mpoly_clear(Q, T->ctx);
    flint_free(s);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Normalisation in the caller's basis                                */
/* ------------------------------------------------------------------ */

/* What each basis member is made of: the exponents of the non-I generators
 * (in rank order) and the exponent of I. */
typedef struct {
    slong n;
    slong ngk;                 /* non-I generators                    */
    slong* gk;                 /* n * ngk exponents                   */
    slong* ie;                 /* n exponents of I                    */
    slong* ip;                 /* the member times I (+-), or -1      */
    int gaussian_only;         /* K = Q(i), basis {1, I}              */
    int ready;                 /* a product basis exists              */
    fmpq_mat_t Pinv, P;        /* theta <-> basis coordinates         */
} TBasis;

static void tbasis_init(TBasis* B, const FgRender* R, const Expr* theta, slong n) {
    memset(B, 0, sizeof *B);
    B->n = n;
    B->ie = flint_calloc((size_t)n, sizeof(slong));
    B->ip = flint_malloc(sizeof(slong) * (size_t)n);
    for (slong m = 0; m < n; m++) B->ip[m] = -1;
    if (R->ready) {
        int iidx = -1;
        for (size_t s = 0; s < R->nsel; s++) if (tr_is_imag_unit(R->sel_atom[s])) iidx = (int)s;
        B->ngk = (slong)R->nsel - (iidx >= 0 ? 1 : 0);
        B->gk = flint_calloc((size_t)(n * (B->ngk ? B->ngk : 1)), sizeof(slong));
        for (slong m = 0; m < n; m++) {
            slong rest = m, g = 0;
            for (size_t s = 0; s < R->nsel; s++) {
                slong d = rest % R->sel_deg[s];
                rest /= R->sel_deg[s];
                if ((int)s == iidx) B->ie[m] = d;
                else B->gk[m * B->ngk + g++] = d;
            }
        }
        fmpq_mat_init(B->Pinv, n, n); fmpq_mat_init(B->P, n, n);
        fmpq_mat_set(B->Pinv, R->Pinv);
        if (fmpq_mat_inv(B->P, B->Pinv)) B->ready = 1;
        else { fmpq_mat_clear(B->Pinv); fmpq_mat_clear(B->P); }
        if (B->ready && iidx >= 0) {
            slong stride = 1;
            for (int s = 0; s < iidx; s++) stride *= R->sel_deg[s];
            for (slong m = 0; m < n; m++) B->ip[m] = B->ie[m] ? m - stride : m + stride;
        }
        if (!B->ready) {
            flint_free(B->gk); B->gk = NULL;
            memset(B->ie, 0, sizeof(slong) * (size_t)n);
        }
    }
    if (!B->ready) {
        /* The theta power basis: theta is the one generator (or I). */
        int is_i = tr_is_imag_unit(theta);
        B->ngk = is_i ? 0 : 1;
        B->gk = flint_calloc((size_t)n, sizeof(slong));
        for (slong m = 0; m < n; m++) {
            if (is_i) B->ie[m] = m; else B->gk[m] = m;
        }
    }
    B->gaussian_only = (n == 2 && B->ngk == 0);
}

static void tbasis_clear(TBasis* B) {
    if (B->ready) { fmpq_mat_clear(B->Pinv); fmpq_mat_clear(B->P); }
    flint_free(B->gk); flint_free(B->ie); flint_free(B->ip);
}

/* Coordinates of a reduced polynomial, grouped by monomial in the variables
 * (LEX-descending, so group 0 holds the leading monomial). */
typedef struct {
    slong nmon, n, nv;
    ulong* mons;               /* nmon * nv full exponent vectors (tau 0) */
    fmpq* co;                  /* nmon * n                                */
} TCoords;

static void tcoords_clear(TCoords* C) {
    if (C->co) { for (slong i = 0; i < C->nmon * C->n; i++) fmpq_clear(C->co + i); flint_free(C->co); }
    flint_free(C->mons);
    memset(C, 0, sizeof *C);
}

static void tcoords_from(TCoords* C, const fmpq_mpoly_t F, TCtx* T) {
    slong len = fmpq_mpoly_length(F, T->ctx);
    C->n = T->n; C->nv = T->nvars; C->nmon = 0;
    C->mons = flint_calloc((size_t)((len ? len : 1) * T->nvars), sizeof(ulong));
    C->co = flint_malloc(sizeof(fmpq) * (size_t)((len ? len : 1) * T->n));
    for (slong i = 0; i < (len ? len : 1) * T->n; i++) fmpq_init(C->co + i);
    ulong* e = flint_malloc(sizeof(ulong) * (size_t)T->nvars);
    fmpq_t q; fmpq_init(q);
    for (slong i = 0; i < len; i++) {
        fmpq_mpoly_get_term_exp_ui(e, F, i, T->ctx);
        int same = C->nmon > 0;
        if (same) {
            const ulong* cur = C->mons + (C->nmon - 1) * T->nvars;
            for (slong v = 0; v < T->ngens; v++) if (e[v] != cur[v]) { same = 0; break; }
        }
        if (!same) {
            ulong* cur = C->mons + C->nmon * T->nvars;
            memcpy(cur, e, sizeof(ulong) * (size_t)T->nvars);
            cur[T->tau_idx] = 0;
            C->nmon++;
        }
        fmpq_mpoly_get_term_coeff_fmpq(q, F, i, T->ctx);
        fmpq_set(C->co + (C->nmon - 1) * T->n + (slong)e[T->tau_idx], q);
    }
    fmpq_clear(q); flint_free(e);
}

static void tcoords_to(fmpq_mpoly_t F, const TCoords* C, TCtx* T) {
    fmpq_mpoly_zero(F, T->ctx);
    ulong* e = flint_malloc(sizeof(ulong) * (size_t)T->nvars);
    for (slong i = 0; i < C->nmon; i++) {
        memcpy(e, C->mons + i * T->nvars, sizeof(ulong) * (size_t)T->nvars);
        for (slong j = 0; j < T->n; j++) {
            if (fmpq_is_zero(C->co + i * T->n + j)) continue;
            e[T->tau_idx] = (ulong)j;
            fmpq_mpoly_set_coeff_fmpq_ui(F, C->co + i * T->n + j, e, T->ctx);
        }
    }
    flint_free(e);
}

/* Each coordinate vector v <- A v (A n x n). */
static void tcoords_apply(TCoords* C, const fmpq_mat_t A) {
    slong n = C->n;
    fmpq* w = flint_malloc(sizeof(fmpq) * (size_t)n);
    for (slong j = 0; j < n; j++) fmpq_init(w + j);
    fmpq_t t; fmpq_init(t);
    for (slong i = 0; i < C->nmon; i++) {
        fmpq* v = C->co + i * n;
        for (slong r = 0; r < n; r++) {
            fmpq_zero(w + r);
            for (slong k = 0; k < n; k++) {
                if (fmpq_is_zero(v + k)) continue;
                fmpq_mul(t, fmpq_mat_entry(A, r, k), v + k);
                fmpq_add(w + r, w + r, t);
            }
        }
        for (slong r = 0; r < n; r++) fmpq_swap(v + r, w + r);
    }
    fmpq_clear(t);
    for (slong j = 0; j < n; j++) fmpq_clear(w + j);
    flint_free(w);
}

/* q <- round(a / d), d > 0. */
static void tr_round_div(fmpz_t q, const fmpz_t a, const fmpz_t d) {
    fmpz_t t, d2; fmpz_init(t); fmpz_init(d2);
    fmpz_mul_2exp(t, a, 1); fmpz_add(t, t, d);
    fmpz_mul_2exp(d2, d, 1);
    fmpz_fdiv_q(q, t, d2);
    fmpz_clear(t); fmpz_clear(d2);
}

/* (ar + ai i) <- gcd over Z[i] of itself and (br + bi i); b clobbered. */
static void tr_gauss_gcd(fmpz_t ar, fmpz_t ai, fmpz_t br, fmpz_t bi) {
    fmpz_t nr, ni, N, qr, qi, t;
    fmpz_init(nr); fmpz_init(ni); fmpz_init(N); fmpz_init(qr); fmpz_init(qi); fmpz_init(t);
    while (!fmpz_is_zero(br) || !fmpz_is_zero(bi)) {
        fmpz_mul(nr, ar, br); fmpz_addmul(nr, ai, bi);
        fmpz_mul(ni, ai, br); fmpz_submul(ni, ar, bi);
        fmpz_mul(N, br, br); fmpz_addmul(N, bi, bi);
        tr_round_div(qr, nr, N); tr_round_div(qi, ni, N);
        fmpz_mul(t, qr, br); fmpz_submul(t, qi, bi); fmpz_sub(nr, ar, t);
        fmpz_mul(t, qr, bi); fmpz_addmul(t, qi, br); fmpz_sub(ni, ai, t);
        fmpz_swap(ar, br); fmpz_swap(ai, bi);
        fmpz_swap(br, nr); fmpz_swap(bi, ni);
    }
    fmpz_clear(nr); fmpz_clear(ni); fmpz_clear(N); fmpz_clear(qr); fmpz_clear(qi); fmpz_clear(t);
}

/* Make the basis coordinates of C integral and primitive over Z, then over
 * Z[i] when K = Q(i), then normalise the unit (see the file comment). */
/* Scale the coordinates to integers with no common factor. */
static void tr_z_primitive(TCoords* C) {
    slong tot = C->nmon * C->n;
    fmpz_t D, g; fmpz_init_set_ui(D, 1); fmpz_init(g);
    for (slong i = 0; i < tot; i++)
        if (!fmpq_is_zero(C->co + i)) fmpz_lcm(D, D, fmpq_denref(C->co + i));
    for (slong i = 0; i < tot; i++) fmpq_mul_fmpz(C->co + i, C->co + i, D);
    for (slong i = 0; i < tot; i++) fmpz_gcd(g, g, fmpq_numref(C->co + i));
    if (!fmpz_is_zero(g) && !fmpz_is_one(g))
        for (slong i = 0; i < tot; i++) fmpq_div_fmpz(C->co + i, C->co + i, g);
    fmpz_clear(D); fmpz_clear(g);
}

static void tr_normalise_coords(TCoords* C, const TBasis* B) {
    slong n = C->n, tot = C->nmon * n;
    tr_z_primitive(C);

    if (B->gaussian_only) {
        /* basis {1, I}: coordinate 0 is the real part, 1 the imaginary */
        fmpz_t gr, gi, br, bi, N, a, b;
        fmpz_init(gr); fmpz_init(gi); fmpz_init(br); fmpz_init(bi);
        fmpz_init(N); fmpz_init(a); fmpz_init(b);
        slong one = (B->ie[0] == 0) ? 0 : 1, im = 1 - one;
        for (slong i = 0; i < C->nmon; i++) {
            fmpz_set(br, fmpq_numref(C->co + i * n + one));
            fmpz_set(bi, fmpq_numref(C->co + i * n + im));
            tr_gauss_gcd(gr, gi, br, bi);
        }
        fmpz_mul(N, gr, gr); fmpz_addmul(N, gi, gi);
        if (!fmpz_is_zero(N) && !fmpz_is_one(N)) {
            for (slong i = 0; i < C->nmon; i++) {
                fmpq* re = C->co + i * n + one; fmpq* ip = C->co + i * n + im;
                /* (x + y i)(gr - gi i) / N */
                fmpz_mul(a, fmpq_numref(re), gr); fmpz_addmul(a, fmpq_numref(ip), gi);
                fmpz_mul(b, fmpq_numref(ip), gr); fmpz_submul(b, fmpq_numref(re), gi);
                fmpz_divexact(a, a, N); fmpz_divexact(b, b, N);
                fmpq_set_fmpz(re, a); fmpq_set_fmpz(ip, b);
            }
        }
        /* rotate by i until the leading coefficient is in the first quadrant */
        for (int k = 0; k < 4; k++) {
            const fmpq* r0 = C->co + one; const fmpq* i0 = C->co + im;
            if (fmpq_sgn(r0) > 0 && fmpq_sgn(i0) >= 0) break;
            for (slong i = 0; i < C->nmon; i++) {
                fmpq* re = C->co + i * n + one; fmpq* ip = C->co + i * n + im;
                fmpq_neg(ip, ip); fmpq_swap(re, ip);      /* (x + y i) i = -y + x i */
            }
        }
        fmpz_clear(gr); fmpz_clear(gi); fmpz_clear(br); fmpz_clear(bi);
        fmpz_clear(N); fmpz_clear(a); fmpz_clear(b);
    } else {
        /* the leading term: generators (in rank order) above the variables
         * (LEX, i.e. group order), I last */
        slong bi = -1, bm = -1;
        for (slong i = 0; i < C->nmon; i++)
            for (slong m = 0; m < n; m++) {
                if (fmpq_is_zero(C->co + i * n + m)) continue;
                if (bi < 0) { bi = i; bm = m; continue; }
                int c = 0;
                for (slong s = 0; s < B->ngk && c == 0; s++) {
                    slong x = B->gk[m * B->ngk + s], y = B->gk[bm * B->ngk + s];
                    if (x != y) c = x > y ? 1 : -1;
                }
                if (c == 0 && i != bi) c = i < bi ? 1 : -1;
                /* the real part of a coefficient leads its I part */
                if (c == 0 && B->ie[m] != B->ie[bm]) c = B->ie[m] < B->ie[bm] ? 1 : -1;
                if (c > 0) { bi = i; bm = m; }
            }
        /* With I in the field, the leading term's Gaussian coefficient
         * u + v I is made a positive rational by its conjugate: -I x +
         * Sqrt[2] y rather than x + I Sqrt[2] y, and Sqrt[2] - (1 + I) x
         * rather than (1 - I) Sqrt[2] - 2 x. */
        if (bi >= 0 && B->ip[bm] >= 0) {
            slong mr = B->ie[bm] ? B->ip[bm] : bm;          /* the real member */
            fmpq_t u, v, a, b, t;
            fmpq_init(u); fmpq_init(v); fmpq_init(a); fmpq_init(b); fmpq_init(t);
            fmpq_set(u, C->co + bi * n + mr);
            fmpq_set(v, C->co + bi * n + B->ip[mr]);
            if (!fmpq_is_zero(v)) {
                for (slong i = 0; i < C->nmon; i++)
                    for (slong m = 0; m < n; m++) {
                        if (B->ie[m] != 0 || B->ip[m] < 0) continue;
                        fmpq* x = C->co + i * n + m; fmpq* y = C->co + i * n + B->ip[m];
                        /* (x + y I)(u - v I) = (x u + y v) + (y u - x v) I */
                        fmpq_mul(a, x, u); fmpq_mul(t, y, v); fmpq_add(a, a, t);
                        fmpq_mul(b, y, u); fmpq_mul(t, x, v); fmpq_sub(b, b, t);
                        fmpq_set(x, a); fmpq_set(y, b);
                    }
                tr_z_primitive(C);
            }
            bm = mr;
            fmpq_clear(u); fmpq_clear(v); fmpq_clear(a); fmpq_clear(b); fmpq_clear(t);
        }
        if (bi >= 0 && fmpq_sgn(C->co + bi * n + bm) < 0)
            for (slong i = 0; i < tot; i++) fmpq_neg(C->co + i, C->co + i);
    }
}

/* ------------------------------------------------------------------ */
/*  Rendering                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    const VarSet* vs;
    const Expr* theta;
    const FgRender* R;
} TRender;

static Expr* tr_render_mpoly(const fmpq_mpoly_t F, TCtx* T, const TRender* TR) {
    Expr* e = fb_field_mpoly_to_expr(F, T->ctx, TR->vs, T->tau_idx, T->n, TR->theta);
    if (!e) return NULL;
    Expr* r = TR->R->ready ? fb_render_tree(e, TR->R) : fb_expand_algnum(e);
    expr_free(e);
    return r ? eval_and_free(r) : NULL;
}

/* ------------------------------------------------------------------ */
/*  The entry point                                                    */
/* ------------------------------------------------------------------ */

static int tr_cmp_exprs(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}

/* List[f, e] pairs by the canonical order of f. */
static int tr_cmp_pairs(const void* a, const void* b) {
    const Expr* x = *(Expr* const*)a;
    const Expr* y = *(Expr* const*)b;
    return expr_compare(x->data.function.args[0], y->data.function.args[0]);
}

static int tr_cmp_names(const void* a, const void* b) {
    return strcmp(*(char* const*)a, *(char* const*)b);
}

/* MATHILDA_NO_TRAGER=1 turns the engine off (every call declines), for A/B
 * comparison against the classical univariate path.  Read once. */
static int tr_disabled(void) {
    static int cached = -1;
    if (cached < 0) {
        const char* v = getenv("MATHILDA_NO_TRAGER");
        cached = (v && v[0] && v[0] != '0') ? 1 : 0;
    }
    return cached;
}

Expr* flint_factor_over_field(const Expr* poly, Expr* const* ext, size_t n_ext) {
    if (!poly || tr_disabled()) return NULL;
    Expr** atoms = NULL; size_t natoms = 0, acap = 0;
    Expr* P = tr_split_complex(poly);
    Expr* imgs = NULL;
    Expr* mp = NULL;
    Expr* Ps = NULL;
    Expr* Pt = NULL;
    Expr* out = NULL;
    VarSet vs; memset(&vs, 0, sizeof vs);
    int ok = 1;

    /* 1. The atoms: the generators first (they rank highest in the output),
     * then the input's own algebraic constants. */
    for (size_t i = 0; i < n_ext && ok; i++) {
        Expr* g = tr_split_complex(ext[i]);
        ok = fb_collect_atoms(g, &atoms, &natoms, &acap);
        expr_free(g);
    }
    if (ok) ok = fb_collect_atoms(P, &atoms, &natoms, &acap);
    /* The atoms in canonical order, whatever order the generators were
     * listed in (Extension -> {Sqrt[3], Sqrt[2]} answers as {Sqrt[2],
     * Sqrt[3]} does), and so the render basis is built from Sqrt[2] and
     * Sqrt[3] rather than from Sqrt[6] and Sqrt[2] (whose product Sqrt[12]
     * would make Sqrt[3] look like (Sqrt[6] Sqrt[2])/2). */
    if (ok && natoms > 1)
        qsort(atoms, natoms, sizeof(Expr*), tr_cmp_exprs);
    if (!ok || natoms == 0 || natoms > FG_RENDER_MAX_ATOMS) goto done;

    /* 2. The field: theta, its (monic, integral) minimal polynomial. */
    {
        const Expr* ci = fb_field_images(atoms, natoms);
        if (!ci) goto done;
        imgs = expr_copy((Expr*)ci);       /* the cache slot may be recycled */
    }
    const Expr* theta = NULL;
    for (size_t i = 0; i < natoms; i++) {
        const Expr* im = imgs->data.function.args[i];
        if (tr_head(im) != SYM_AlgebraicNumber || im->data.function.arg_count != 2) goto done;
        if (!theta) theta = im->data.function.args[0];
        else if (!expr_eq((Expr*)theta, im->data.function.args[0])) goto done;
    }
    mp = flint_qqbar_gen_minpoly_coeffs(theta);
    if (!mp || mp->type != EXPR_FUNCTION || mp->data.function.arg_count < 3) goto done;
    slong n = (slong)mp->data.function.arg_count - 1;
    if (!tr_is_int(mp->data.function.args[n], 1)) goto done;
    for (slong j = 0; j <= n; j++) {
        const Expr* c = mp->data.function.args[j];
        if (c->type != EXPR_INTEGER && c->type != EXPR_BIGINT) goto done;
    }

    /* 3. The tau lift and the variables (sorted, tau last). */
    Ps = fb_subst_atoms(P, atoms, imgs->data.function.args, natoms);
    if (!Ps) goto done;
    Pt = fb_field_subst_tau(Ps, tr_is_imag_unit(theta));
    if (!Pt || !fb_collect_vars(Pt, &vs)) goto done;
    {
        Expr* tau = expr_new_symbol(FIELD_TAU_NAME);
        int tok = fb_collect_vars(tau, &vs);
        expr_free(tau);
        if (!tok) goto done;
        size_t ti = 0;
        for (size_t i = 0; i < vs.count; i++) if (strcmp(vs.names[i], FIELD_TAU_NAME) == 0) ti = i;
        char* tn = vs.names[ti];
        vs.names[ti] = vs.names[vs.count - 1];
        vs.names[vs.count - 1] = tn;
        qsort(vs.names, vs.count - 1, sizeof(char*), tr_cmp_names);
    }
    if (vs.count < 2) goto done;                          /* a constant */

    {
    TCtx T; memset(&T, 0, sizeof T);
    T.nvars = (slong)vs.count; T.ngens = T.nvars - 1; T.tau_idx = T.nvars - 1; T.n = n;
    T.ctx = fb_ctx_lex(T.nvars);
    T.mp = mp;
    fmpq_mpoly_init(T.M, T.ctx);
    {
        fmpq_poly_t mq; fmpq_poly_init(mq);
        fmpq_t c; fmpq_init(c);
        for (slong j = 0; j <= n; j++) {
            fb_expr_to_fmpq(mp->data.function.args[j], c);
            fmpq_poly_set_coeff_fmpq(mq, j, c);
        }
        nf_init(T.nf, mq);
        fmpq_clear(c); fmpq_poly_clear(mq);
    }
    fmpq_mpoly_t A, A1, Prod, Fm;
    fmpq_poly_t lc, r, cK;
    fmpq_mpoly_init(A, T.ctx); fmpq_mpoly_init(A1, T.ctx);
    fmpq_mpoly_init(Prod, T.ctx); fmpq_mpoly_init(Fm, T.ctx);
    fmpq_poly_init(lc); fmpq_poly_init(r); fmpq_poly_init(cK);
    TFac L; memset(&L, 0, sizeof L);
    FgRender R; memset(&R, 0, sizeof R);
    TBasis B; memset(&B, 0, sizeof B);
    int have_render = 0, have_basis = 0;
    Expr** fexprs = NULL;

    ok = fb_to_mpoly(Pt, A, T.ctx, &vs) &&
         fb_field_build_minpoly(T.M, mp, T.ctx, T.nvars, T.tau_idx);
    if (ok) { tr_reduce(A, &T); ok = !fmpq_mpoly_is_zero(A, T.ctx) && tr_xdeg(A, &T) > 0; }

    /* 4. Monic, then split. */
    if (ok) {
        tr_lead_kelem(lc, A, &T);
        fmpq_poly_set(r, lc);
        tr_kinv(r, &T);
        fmpq_mpoly_set(A1, A, T.ctx);
        tr_scale(A1, r, &T);
        ok = tr_split(A1, &T, &L);
    }

    /* 5. Normalise each factor in the caller's basis; c = lc / prod lc(F)^e. */
    if (ok) {
        fb_render_init(&R, atoms, imgs->data.function.args, natoms, theta, n);
        have_render = 1;
        tbasis_init(&B, &R, theta, n);
        have_basis = 1;
        fmpq_poly_set(cK, lc);
        for (size_t i = 0; i < L.num; i++) {
            TCoords C; memset(&C, 0, sizeof C);
            tcoords_from(&C, L.f + i, &T);
            if (B.ready) tcoords_apply(&C, B.Pinv);
            tr_normalise_coords(&C, &B);
            if (B.ready) tcoords_apply(&C, B.P);
            tcoords_to(L.f + i, &C, &T);
            tcoords_clear(&C);
            tr_lead_kelem(r, L.f + i, &T);
            fmpq_poly_t rp; fmpq_poly_init(rp);
            fmpq_poly_one(rp);
            for (slong k = 0; k < L.e[i]; k++) tr_kmul(rp, rp, r, &T);
            tr_kinv(rp, &T);
            tr_kmul(cK, cK, rp, &T);
            fmpq_poly_clear(rp);
        }

        /* 6. Certify: c prod F^e == A in Q[x, tau]/(M). */
        tr_kelem_to_mpoly(Prod, cK, &T);
        for (size_t i = 0; i < L.num && ok; i++)
            for (slong k = 0; k < L.e[i]; k++) {
                fmpq_mpoly_mul(Prod, Prod, L.f + i, T.ctx);
                tr_reduce(Prod, &T);
            }
        ok = fmpq_mpoly_equal(Prod, A, T.ctx);
    }

    /* 7. Render. */
    if (ok) {
        TRender TR = { &vs, theta, &R };
        fexprs = calloc(L.num + 1, sizeof(Expr*));
        tr_kelem_to_mpoly(Fm, cK, &T);
        fexprs[0] = tr_render_mpoly(Fm, &T, &TR);
        ok = fexprs[0] != NULL;
        for (size_t i = 0; i < L.num && ok; i++) {
            Expr* f = tr_render_mpoly(L.f + i, &T, &TR);
            if (!f) { ok = 0; break; }
            fexprs[i + 1] = expr_new_function(expr_new_symbol(SYM_List),
                (Expr*[]){ f, expr_new_integer((int64_t)L.e[i]) }, 2);
        }
        if (ok) {
            /* the factors in canonical order (FactorList shows them so) */
            qsort(fexprs + 1, L.num, sizeof(Expr*), tr_cmp_pairs);
            out = expr_new_function(expr_new_symbol(SYM_List), fexprs, L.num + 1);
            free(fexprs);
            fexprs = NULL;
        }
    }
    if (fexprs) {
        for (size_t i = 0; i <= L.num; i++) if (fexprs[i]) expr_free(fexprs[i]);
        free(fexprs);
    }

    if (have_basis) tbasis_clear(&B);
    if (have_render) fb_render_clear(&R);
    tfac_clear(&L, &T);
    fmpq_poly_clear(lc); fmpq_poly_clear(r); fmpq_poly_clear(cK);
    fmpq_mpoly_clear(A, T.ctx); fmpq_mpoly_clear(A1, T.ctx);
    fmpq_mpoly_clear(Prod, T.ctx); fmpq_mpoly_clear(Fm, T.ctx);
    fmpq_mpoly_clear(T.M, T.ctx);
    nf_clear(T.nf);
    }

done:
    fb_varset_free(&vs);
    for (size_t i = 0; i < natoms; i++) expr_free(atoms[i]);
    free(atoms);
    expr_free(P); expr_free(Ps); expr_free(Pt); expr_free(mp); expr_free(imgs);
    return out;
}

/* {c, {f1, e1}, ...} -> c f1^e1 ... , evaluated. */
Expr* flint_factor_over_field_product(const Expr* poly, Expr* const* ext, size_t n_ext) {
    Expr* L = flint_factor_over_field(poly, ext, n_ext);
    if (!L) return NULL;
    size_t n = L->data.function.arg_count;
    Expr** t = malloc(sizeof(Expr*) * n);
    t[0] = expr_copy(L->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        const Expr* pr = L->data.function.args[i];
        t[i] = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(pr->data.function.args[0]), expr_copy(pr->data.function.args[1]) }, 2);
    }
    Expr* out = eval_and_free(n == 1 ? t[0] : expr_new_function(expr_new_symbol(SYM_Times), t, n));
    free(t);
    expr_free(L);
    return out;
}

#else  /* !USE_FLINT */

int flint_trager_has_complex(const Expr* e) { (void)e; return 0; }
Expr* flint_factor_over_field(const Expr* poly, Expr* const* ext, size_t n_ext) {
    (void)poly; (void)ext; (void)n_ext;
    return NULL;
}
Expr* flint_factor_over_field_product(const Expr* poly, Expr* const* ext, size_t n_ext) {
    (void)poly; (void)ext; (void)n_ext;
    return NULL;
}

#endif /* USE_FLINT */

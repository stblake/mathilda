/* nfrat.c -- Together / Cancel over the number field of the constant algebraic
 * coefficients, natively.  See nfrat.h for the contract.
 *
 * Representation.  K = Q(theta), theta an algebraic integer with monic integral
 * minimal polynomial M.  A polynomial over K in the variables x1..xn is an
 * fmpq_mpoly in (x1..xn, tau), LEX, tau last, with tau-degree < [K:Q] (reduced
 * modulo M(tau)).  A K constant is such a polynomial in tau alone.
 *
 * A rational function is kept factored:
 *
 *      c * prod_j Q_j^f_j  /  ( u * prod_i P_i^e_i )
 *
 * where c, u are K constants and every Q_j, P_i is a MONIC piece (the
 * coefficient of its lex-leading x-monomial is 1) interned in a per-call piece
 * table.  u collects the leading coefficients of the inverted sums and the
 * written constant denominators; it is never inverted, only multiplied, so a
 * numeric algebraic denominator survives as a denominator (the reference does
 * not rationalise).  A rational u is folded into c at once.
 *
 *   Times / Power  merge the piece lists (exponents add / scale).
 *   Plus           takes the least common denominator over the piece table --
 *                  identical pieces share an index, so 1/P + 1/P^2 is over P^2,
 *                  not P^3 -- and expands each term's numerator times its
 *                  cofactor into ONE new numerator piece.
 *   inversion      swaps the lists and the two constants.
 *
 * At the end every numerator piece is cancelled against every denominator piece
 * with the certified modular gcd over K (fb_field_mpoly_gcd, any number of
 * variables), the denominator is refined to a squarefree, gcd-free basis, and
 * the answer is written back in the caller's own radicals (FgRender).
 *
 * Nothing here calls the evaluator except on the final, already-reduced answer
 * (and on a handful of atom-sized constants).  The previous paths evaluated the
 * substituted tree or a huge expanded num/den with radical coefficients, which
 * re-entered Times' radical canonicaliser on every intermediate product. */

#include "nfrat.h"

#include "sym_names.h"
#include "eval.h"
#include "common.h"
#include "core.h"
#include "arithmetic.h"
#include "flint_qqbar.h"
#include "ratfun.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  The structural gate (no FLINT needed)                              */
/* ------------------------------------------------------------------ */

static bool nr_is_ratnum(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
           (e->type == EXPR_FUNCTION && head_is(e, SYM_Rational) && e->data.function.arg_count == 2);
}

/* Does `e` contain a symbol? */
static bool nr_has_symbol(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nr_has_symbol(e->data.function.args[i])) return true;
    return false;
}

/* A constant algebraic leaf of the rational structure: Complex[a, b] with
 * rational parts, a rational power of a constant, a constant Root object. */
static bool nr_is_alg_leaf(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2)
        return nr_is_ratnum(e->data.function.args[0]) && nr_is_ratnum(e->data.function.args[1]);
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        return head_is(x, SYM_Rational) && !nr_has_symbol(e->data.function.args[0]);
    }
    if (head_is(e, SYM_Root)) return !nr_has_symbol(e);
    return false;
}

static bool nr_gate_walk(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (nr_is_alg_leaf(e)) return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nr_gate_walk(e->data.function.args[i])) return true;
        return false;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER)
        return nr_gate_walk(e->data.function.args[0]);
    return false;
}

bool nfrat_has_algebraic(const Expr* e) { return e && nr_gate_walk(e); }

/* Root[f, k]^m in the rational structure with |m| >= deg f (a power the
 * reference keeps unreduced). */
static bool nr_root_power_big(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        const Expr* b = e->data.function.args[0];
        if (head_is(b, SYM_Root) && !nr_has_symbol(b)) {
            int64_t m = e->data.function.args[1]->data.integer;
            if (m < 0) m = -m;
            if (m >= 2) {
                Expr* mp = flint_qqbar_gen_minpoly_coeffs(b);
                int64_t d = (mp && mp->type == EXPR_FUNCTION) ? (int64_t)mp->data.function.arg_count - 1 : 0;
                if (mp) expr_free(mp);
                if (d > 0 && m >= d) return true;
            }
            return false;
        }
        return nr_root_power_big(b);
    }
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nr_root_power_big(e->data.function.args[i])) return true;
    }
    return false;
}

bool nfrat_has_unreduced_root_power(const Expr* e) { return e && nr_root_power_big(e); }

/* A negative power (integer or rational: 1/(1 + x), 1/Sqrt[2]) inside the
 * rational structure. */
bool nfrat_has_denominator(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER)
            return x->data.integer < 0 || nfrat_has_denominator(e->data.function.args[0]);
        if (head_is(x, SYM_Rational) && x->data.function.args[0]->type == EXPR_INTEGER)
            return x->data.function.args[0]->data.integer < 0;
        return false;
    }
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nfrat_has_denominator(e->data.function.args[i])) return true;
    }
    return false;
}

/* Distinct algebraic leaves, Complex counted as the one atom I. */
static void nr_count_walk(const Expr* e, const Expr** seen, int* n, int stop, bool* imag) {
    if (*n >= stop || e->type != EXPR_FUNCTION) return;
    if (nr_is_alg_leaf(e)) {
        if (head_is(e, SYM_Complex)) {
            if (!*imag) { *imag = true; (*n)++; }
            return;
        }
        /* b^(-p/q) and b^(p/q) are one generator */
        const Expr* key = e;
        if (head_is(e, SYM_Power)) key = e->data.function.args[0];
        for (int i = 0; i < *n; i++)
            if (seen[i] && expr_eq((Expr*)seen[i], (Expr*)key)) return;
        seen[(*n)++] = key;
        return;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nr_count_walk(e->data.function.args[i], seen, n, stop, imag);
}

int nfrat_count_atoms(const Expr* e, int stop) {
    if (!e || stop <= 0) return 0;
    const Expr** seen = calloc((size_t)stop + 1, sizeof(Expr*));
    if (!seen) return 0;
    int n = 0; bool imag = false;
    nr_count_walk(e, seen, &n, stop, &imag);
    free(seen);
    return n;
}

#ifdef USE_FLINT

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

#include "flint_bridge_internal.h"

/* Decline rather than expand past this many terms in one numerator. */
#define NR_TERM_BUDGET  400000
/* Decline above this field degree (the product-basis render is n x n). */
#define NR_DEGREE_CAP   64
/* Distinct polynomial variables. */
#define NR_MAX_VARS     16

/* ------------------------------------------------------------------ */
/*  Context                                                            */
/* ------------------------------------------------------------------ */

typedef struct { size_t i; slong e; } NrPE;

typedef struct {
    const fmpq_mpoly_ctx_struct* ctx;   /* cached by fb_ctx_lex, not ours */
    slong ngens, tau, n;
    fmpq_mpoly_t M;
    Expr* mp;                 /* owned: theta's minimal polynomial coefficients */
    Expr* theta;              /* owned copy                                    */
    nf_t nf;
    gr_ctx_t gr;              /* the same field as a gr ring (Kronecker division) */
    Expr** atoms; size_t na;  /* owned                                          */
    fmpq_mpoly_struct* aimg;  /* the atoms' images, tau polynomials             */
    const char* var[NR_MAX_VARS];
    FgRender R;
    /* tau^j mod M for j = n .. 2n - 2 (the tau-degrees a product of two
     * reduced polynomials reaches): the term-wise reduction */
    fmpq_poly_struct* tpow; slong ntpow;
    /* the piece table: monic, reduced, distinct */
    fmpq_mpoly_struct* pc;
    fmpq_mpoly_struct* wl;    /* written leading coefficient (1 if none)        */
    size_t np, pcap;
    /* known-coprime pairs (i < j) */
    size_t* cop; size_t ncop, ccop;
    int ready;
} Nr;

typedef struct {
    fmpq_mpoly_t c, u;
    NrPE* nu; size_t nnu, cnu;
    NrPE* de; size_t nde, cde;
} NrRf;

/* p := p mod M(tau).  M is monic of degree n in tau alone, so a term
 * c x^a tau^j with n <= j <= 2n - 2 is replaced by c x^a (tau^j mod M) from
 * the precomputed table: one pass over the terms, a sort and a merge.  The
 * general multivariate division by M was most of a large Together (each
 * product of reduced polynomials went through fmpq_mpoly_divrem). */
static void nr_reduce(Nr* N, fmpq_mpoly_t p) {
    if (fmpq_mpoly_is_zero(p, N->ctx)) return;
    slong dt = fmpq_mpoly_degree_si(p, N->tau, N->ctx);
    if (dt < N->n) return;
    if (N->tpow && dt - N->n < N->ntpow) {
        fmpq_mpoly_t out; fmpq_mpoly_init(out, N->ctx);
        ulong e[NR_MAX_VARS + 1];
        fmpq_t c, v; fmpq_init(c); fmpq_init(v);
        slong len = fmpq_mpoly_length(p, N->ctx);
        for (slong i = 0; i < len; i++) {
            fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
            fmpq_mpoly_get_term_coeff_fmpq(c, p, i, N->ctx);
            slong j = (slong)e[N->tau];
            if (j < N->n) { fmpq_mpoly_push_term_fmpq_ui(out, c, e, N->ctx); continue; }
            const fmpq_poly_struct* rj = N->tpow + (j - N->n);
            for (slong k = 0; k <= fmpq_poly_degree(rj); k++) {
                fmpq_poly_get_coeff_fmpq(v, rj, k);
                if (fmpq_is_zero(v)) continue;
                fmpq_mul(v, v, c);
                e[N->tau] = (ulong)k;
                fmpq_mpoly_push_term_fmpq_ui(out, v, e, N->ctx);
            }
        }
        fmpq_mpoly_sort_terms(out, N->ctx);
        fmpq_mpoly_combine_like_terms(out, N->ctx);
        fmpq_mpoly_swap(p, out, N->ctx);
        fmpq_mpoly_clear(out, N->ctx);
        fmpq_clear(c); fmpq_clear(v);
        return;
    }
    fmpq_mpoly_t q, r;
    fmpq_mpoly_init(q, N->ctx); fmpq_mpoly_init(r, N->ctx);
    fmpq_mpoly_divrem(q, r, p, N->M, N->ctx);
    fmpq_mpoly_swap(p, r, N->ctx);
    fmpq_mpoly_clear(q, N->ctx); fmpq_mpoly_clear(r, N->ctx);
}

/* dst = a * b reduced (dst may alias). */
static void nr_mul(Nr* N, fmpq_mpoly_t dst, const fmpq_mpoly_t a, const fmpq_mpoly_t b) {
    fmpq_mpoly_mul(dst, a, b, N->ctx);
    nr_reduce(N, dst);
}

/* p is a tau polynomial (no x): its lex-leading term has zero x-exponents. */
static bool nr_is_const(Nr* N, const fmpq_mpoly_t p) {
    if (fmpq_mpoly_is_zero(p, N->ctx)) return true;
    ulong e[NR_MAX_VARS + 1];
    fmpq_mpoly_get_term_exp_ui(e, p, 0, N->ctx);
    for (slong k = 0; k < N->ngens; k++) if (e[k]) return false;
    return true;
}

/* A constant with no tau: rational. */
static bool nr_is_rational_const(Nr* N, const fmpq_mpoly_t p) {
    if (fmpq_mpoly_is_zero(p, N->ctx)) return true;
    return fmpq_mpoly_is_fmpq(p, N->ctx);
}

/* tau polynomial <-> nf_elem */
static void nr_const_to_nf(Nr* N, nf_elem_t out, const fmpq_mpoly_t c) {
    fmpq_poly_t q; fmpq_poly_init(q);
    ulong e[NR_MAX_VARS + 1];
    fmpq_t v; fmpq_init(v);
    for (slong i = 0; i < fmpq_mpoly_length(c, N->ctx); i++) {
        fmpq_mpoly_get_term_exp_ui(e, c, i, N->ctx);
        fmpq_mpoly_get_term_coeff_fmpq(v, c, i, N->ctx);
        fmpq_poly_set_coeff_fmpq(q, (slong)e[N->tau], v);
    }
    nf_elem_set_fmpq_poly(out, q, N->nf);
    fmpq_clear(v); fmpq_poly_clear(q);
}

static void nr_nf_to_const(Nr* N, fmpq_mpoly_t out, const nf_elem_t a) {
    fmpq_poly_t q; fmpq_poly_init(q);
    nf_elem_get_fmpq_poly(q, a, N->nf);
    fmpq_mpoly_zero(out, N->ctx);
    ulong e[NR_MAX_VARS + 1];
    memset(e, 0, sizeof e);
    fmpq_t v; fmpq_init(v);
    for (slong j = 0; j <= fmpq_poly_degree(q); j++) {
        fmpq_poly_get_coeff_fmpq(v, q, j);
        if (fmpq_is_zero(v)) continue;
        e[N->tau] = (ulong)j;
        fmpq_mpoly_set_coeff_fmpq_ui(out, v, e, N->ctx);
    }
    fmpq_clear(v); fmpq_poly_clear(q);
}

/* out = 1/c over K (c a nonzero constant).  0 when c is zero. */
static int nr_const_inv(Nr* N, fmpq_mpoly_t out, const fmpq_mpoly_t c) {
    if (fmpq_mpoly_is_zero(c, N->ctx)) return 0;
    nf_elem_t a; nf_elem_init(a, N->nf);
    nr_const_to_nf(N, a, c);
    nf_elem_inv(a, a, N->nf);
    nr_nf_to_const(N, out, a);
    nf_elem_clear(a, N->nf);
    return 1;
}

/* The coefficient of the lex-leading x-monomial of p (a constant), and that
 * monomial's exponent vector. */
static void nr_lead_coeff(Nr* N, fmpq_mpoly_t lc, const fmpq_mpoly_t p) {
    fmpq_mpoly_zero(lc, N->ctx);
    ulong e0[NR_MAX_VARS + 1], e[NR_MAX_VARS + 1];
    fmpq_mpoly_get_term_exp_ui(e0, p, 0, N->ctx);
    fmpq_t v; fmpq_init(v);
    for (slong i = 0; i < fmpq_mpoly_length(p, N->ctx); i++) {
        fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
        bool same = true;
        for (slong k = 0; k < N->ngens && same; k++) same = (e[k] == e0[k]);
        if (!same) break;                    /* LEX, tau last: the group is a prefix */
        fmpq_mpoly_get_term_coeff_fmpq(v, p, i, N->ctx);
        ulong t[NR_MAX_VARS + 1];
        memset(t, 0, sizeof t);
        t[N->tau] = e[N->tau];
        fmpq_mpoly_set_coeff_fmpq_ui(lc, v, t, N->ctx);
    }
    fmpq_clear(v);
}

/* p := p / lc(p) (monic); lc returned.  p must be non-constant. */
static int nr_make_monic(Nr* N, fmpq_mpoly_t p, fmpq_mpoly_t lc) {
    nr_lead_coeff(N, lc, p);
    if (fmpq_mpoly_is_one(lc, N->ctx)) return 1;
    fmpq_mpoly_t inv; fmpq_mpoly_init(inv, N->ctx);
    int ok = 0;
    if (fmpq_mpoly_is_fmpq(lc, N->ctx)) {
        fmpq_t q; fmpq_init(q);
        fmpq_mpoly_get_fmpq(q, lc, N->ctx);
        if (!fmpq_is_zero(q)) {
            fmpq_inv(q, q);
            fmpq_mpoly_scalar_mul_fmpq(p, p, q, N->ctx);
            ok = 1;
        }
        fmpq_clear(q);
    } else if (nr_const_inv(N, inv, lc)) {
        nr_mul(N, p, p, inv);
        ok = 1;
    }
    fmpq_mpoly_clear(inv, N->ctx);
    return ok;
}

/* Intern a MONIC piece; `wl` (may be NULL) is its written leading coefficient. */
static size_t nr_piece(Nr* N, const fmpq_mpoly_t p, const fmpq_mpoly_t wl) {
    slong len = fmpq_mpoly_length(p, N->ctx);
    for (size_t i = 0; i < N->np; i++)
        if (fmpq_mpoly_length(N->pc + i, N->ctx) == len && fmpq_mpoly_equal(N->pc + i, p, N->ctx)) {
            if (wl && fmpq_mpoly_is_one(N->wl + i, N->ctx) && !fmpq_mpoly_is_fmpq(wl, N->ctx))
                fmpq_mpoly_set(N->wl + i, wl, N->ctx);
            return i;
        }
    if (N->np == N->pcap) {
        size_t nc = N->pcap ? 2 * N->pcap : 16;
        N->pc = realloc(N->pc, nc * sizeof(fmpq_mpoly_struct));
        N->wl = realloc(N->wl, nc * sizeof(fmpq_mpoly_struct));
        N->pcap = nc;
    }
    fmpq_mpoly_init(N->pc + N->np, N->ctx);
    fmpq_mpoly_set(N->pc + N->np, p, N->ctx);
    fmpq_mpoly_init(N->wl + N->np, N->ctx);
    if (wl && !fmpq_mpoly_is_fmpq(wl, N->ctx)) fmpq_mpoly_set(N->wl + N->np, wl, N->ctx);
    else fmpq_mpoly_one(N->wl + N->np, N->ctx);
    return N->np++;
}

static bool nr_coprime_known(Nr* N, size_t a, size_t b) {
    if (a > b) { size_t t = a; a = b; b = t; }
    for (size_t k = 0; k < N->ncop; k++)
        if (N->cop[2 * k] == a && N->cop[2 * k + 1] == b) return true;
    return false;
}

static void nr_coprime_note(Nr* N, size_t a, size_t b) {
    if (a > b) { size_t t = a; a = b; b = t; }
    if (N->ncop == N->ccop) {
        N->ccop = N->ccop ? 2 * N->ccop : 32;
        N->cop = realloc(N->cop, 2 * N->ccop * sizeof(size_t));
    }
    N->cop[2 * N->ncop] = a; N->cop[2 * N->ncop + 1] = b;
    N->ncop++;
}

/* ------------------------------------------------------------------ */
/*  Rational functions                                                 */
/* ------------------------------------------------------------------ */

static void rf_init(Nr* N, NrRf* r) {
    memset(r, 0, sizeof *r);
    fmpq_mpoly_init(r->c, N->ctx); fmpq_mpoly_one(r->c, N->ctx);
    fmpq_mpoly_init(r->u, N->ctx); fmpq_mpoly_one(r->u, N->ctx);
}

static void rf_clear(Nr* N, NrRf* r) {
    fmpq_mpoly_clear(r->c, N->ctx); fmpq_mpoly_clear(r->u, N->ctx);
    free(r->nu); free(r->de);
}


static void pe_add(NrPE** L, size_t* n, size_t* cap, size_t i, slong e) {
    if (e == 0) return;
    for (size_t k = 0; k < *n; k++)
        if ((*L)[k].i == i) {
            (*L)[k].e += e;
            if ((*L)[k].e == 0) { (*L)[k] = (*L)[*n - 1]; (*n)--; }
            return;
        }
    if (*n == *cap) { *cap = *cap ? 2 * *cap : 4; *L = realloc(*L, *cap * sizeof(NrPE)); }
    (*L)[*n].i = i; (*L)[*n].e = e; (*n)++;
}

static slong pe_get(const NrPE* L, size_t n, size_t i) {
    for (size_t k = 0; k < n; k++) if (L[k].i == i) return L[k].e;
    return 0;
}

/* A rational u is folded into c (rational inversion only). */
static int rf_fold_u(Nr* N, NrRf* r) {
    if (fmpq_mpoly_is_one(r->u, N->ctx)) return 1;
    if (fmpq_mpoly_is_zero(r->u, N->ctx)) return 0;          /* 1/0 */
    if (!fmpq_mpoly_is_fmpq(r->u, N->ctx)) return 1;
    fmpq_t q; fmpq_init(q);
    fmpq_mpoly_get_fmpq(q, r->u, N->ctx);
    fmpq_inv(q, q);
    fmpq_mpoly_scalar_mul_fmpq(r->c, r->c, q, N->ctx);
    fmpq_clear(q);
    fmpq_mpoly_one(r->u, N->ctx);
    return 1;
}

static void rf_set_zero(Nr* N, NrRf* r) {
    fmpq_mpoly_zero(r->c, N->ctx); fmpq_mpoly_one(r->u, N->ctx);
    r->nnu = r->nde = 0;
}

static bool rf_is_zero(Nr* N, const NrRf* r) { return fmpq_mpoly_is_zero(r->c, N->ctx); }

/* r := the polynomial p (consumed into r's numerator). */
static int rf_set_poly(Nr* N, NrRf* r, fmpq_mpoly_t p) {
    r->nnu = r->nde = 0;
    fmpq_mpoly_one(r->u, N->ctx);
    if (nr_is_const(N, p)) { fmpq_mpoly_swap(r->c, p, N->ctx); return 1; }
    if (!nr_make_monic(N, p, r->c)) return 0;
    size_t i = nr_piece(N, p, r->c);
    pe_add(&r->nu, &r->nnu, &r->cnu, i, 1);
    return 1;
}

/* dst := a * b */
static void rf_mul_into(Nr* N, NrRf* dst, const NrRf* b) {
    if (rf_is_zero(N, dst)) return;
    if (rf_is_zero(N, b)) { rf_set_zero(N, dst); return; }
    nr_mul(N, dst->c, dst->c, b->c);
    nr_mul(N, dst->u, dst->u, b->u);
    for (size_t k = 0; k < b->nnu; k++) pe_add(&dst->nu, &dst->nnu, &dst->cnu, b->nu[k].i, b->nu[k].e);
    for (size_t k = 0; k < b->nde; k++) pe_add(&dst->de, &dst->nde, &dst->cde, b->de[k].i, b->de[k].e);
    /* a piece on both sides cancels at once */
    for (size_t k = 0; k < dst->nnu; ) {
        slong d = pe_get(dst->de, dst->nde, dst->nu[k].i);
        if (d > 0) {
            slong m = d < dst->nu[k].e ? d : dst->nu[k].e;
            size_t idx = dst->nu[k].i;
            pe_add(&dst->de, &dst->nde, &dst->cde, idx, -m);
            size_t before = dst->nnu;
            pe_add(&dst->nu, &dst->nnu, &dst->cnu, idx, -m);
            if (dst->nnu < before) continue;      /* entry k replaced: recheck k */
        }
        k++;
    }
}

/* c^k for a constant, k >= 0 */
static void nr_const_pow(Nr* N, fmpq_mpoly_t c, ulong k) {
    fmpq_mpoly_t acc, b;
    fmpq_mpoly_init(acc, N->ctx); fmpq_mpoly_one(acc, N->ctx);
    fmpq_mpoly_init(b, N->ctx); fmpq_mpoly_set(b, c, N->ctx);
    while (k) {
        if (k & 1) nr_mul(N, acc, acc, b);
        k >>= 1;
        if (k) nr_mul(N, b, b, b);
    }
    fmpq_mpoly_swap(c, acc, N->ctx);
    fmpq_mpoly_clear(acc, N->ctx); fmpq_mpoly_clear(b, N->ctx);
}

/* r := r^k, any integer k.  0 on 1/0. */
static int rf_pow(Nr* N, NrRf* r, slong k) {
    if (k == 0) { rf_set_zero(N, r); fmpq_mpoly_one(r->c, N->ctx); return 1; }
    if (k < 0) {
        if (rf_is_zero(N, r)) return 0;
        fmpq_mpoly_swap(r->c, r->u, N->ctx);
        NrPE* t = r->nu; r->nu = r->de; r->de = t;
        size_t tn = r->nnu; r->nnu = r->nde; r->nde = tn;
        tn = r->cnu; r->cnu = r->cde; r->cde = tn;
        k = -k;
    }
    if (k != 1) {
        nr_const_pow(N, r->c, (ulong)k);
        nr_const_pow(N, r->u, (ulong)k);
        for (size_t i = 0; i < r->nnu; i++) r->nu[i].e *= k;
        for (size_t i = 0; i < r->nde; i++) r->de[i].e *= k;
    }
    return rf_fold_u(N, r);
}

/* acc := acc * P_i^e (expanded, reduced). */
static int nr_mul_piece_pow(Nr* N, fmpq_mpoly_t acc, size_t i, slong e) {
    if (e <= 0) return 1;
    double est = (double)fmpq_mpoly_length(acc, N->ctx);
    double pl = (double)fmpq_mpoly_length(N->pc + i, N->ctx);
    for (slong k = 0; k < e; k++) {
        est *= pl;
        if (est > NR_TERM_BUDGET * 50.0) break;
    }
    fmpq_mpoly_t pw; fmpq_mpoly_init(pw, N->ctx);
    int ok;
    if (e == 1) {
        nr_mul(N, acc, acc, N->pc + i);
        ok = 1;
    } else {
        ok = fmpq_mpoly_pow_ui(pw, N->pc + i, (ulong)e, N->ctx);
        if (ok) { nr_reduce(N, pw); nr_mul(N, acc, acc, pw); }
    }
    fmpq_mpoly_clear(pw, N->ctx);
    return ok && fmpq_mpoly_length(acc, N->ctx) <= NR_TERM_BUDGET;
}

/* dst := sum of the k terms (each consumed into the sum's structure). */
static int rf_sum(Nr* N, NrRf* dst, NrRf* T, size_t k) {
    /* the common denominator: max exponent per piece; the distinct u's */
    NrPE* L = NULL; size_t nl = 0, cl = 0;
    fmpq_mpoly_struct* U = NULL; size_t nu = 0;
    int ok = 1;
    size_t live = 0;
    for (size_t t = 0; t < k; t++) {
        if (rf_is_zero(N, T + t)) continue;
        live++;
        for (size_t j = 0; j < T[t].nde; j++) {
            slong have = pe_get(L, nl, T[t].de[j].i);
            if (T[t].de[j].e > have) pe_add(&L, &nl, &cl, T[t].de[j].i, T[t].de[j].e - have);
        }
        if (!fmpq_mpoly_is_one(T[t].u, N->ctx)) {
            bool seen = false;
            for (size_t j = 0; j < nu && !seen; j++) seen = fmpq_mpoly_equal(U + j, T[t].u, N->ctx);
            if (!seen) {
                U = realloc(U, (nu + 1) * sizeof(fmpq_mpoly_struct));
                fmpq_mpoly_init(U + nu, N->ctx);
                fmpq_mpoly_set(U + nu, T[t].u, N->ctx);
                nu++;
            }
        }
    }
    if (live == 0) { rf_set_zero(N, dst); free(L); return 1; }

    fmpq_mpoly_t S, acc;
    fmpq_mpoly_init(S, N->ctx); fmpq_mpoly_init(acc, N->ctx);
    /* A sum of simple fractions (every term over at most one piece): the
     * cofactor of the term over piece j is prod_{i != j} P_i^e_i, the product
     * of a prefix and a suffix -- 3 m multiplications for the whole sum
     * instead of m per term (a 17-term Gaussian Together made 272). */
    fmpq_mpoly_struct* pre = NULL;     /* pre[j] = prod_{i < j} P_i^e_i, j = 0..nl */
    fmpq_mpoly_struct* suf = NULL;     /* suf[j] = prod_{i >= j} P_i^e_i            */
    bool simple = nl >= 3;
    for (size_t t = 0; t < k && simple; t++)
        if (!rf_is_zero(N, T + t) && T[t].nde > 1) simple = false;
    if (simple) {
        pre = malloc((nl + 1) * sizeof(fmpq_mpoly_struct));
        suf = malloc((nl + 1) * sizeof(fmpq_mpoly_struct));
        for (size_t j = 0; j <= nl; j++) { fmpq_mpoly_init(pre + j, N->ctx); fmpq_mpoly_init(suf + j, N->ctx); }
        fmpq_mpoly_one(pre + 0, N->ctx);
        fmpq_mpoly_one(suf + nl, N->ctx);
        for (size_t j = 0; j < nl && ok; j++) {
            fmpq_mpoly_set(pre + j + 1, pre + j, N->ctx);
            ok = nr_mul_piece_pow(N, pre + j + 1, L[j].i, L[j].e);
        }
        for (size_t j = nl; j-- > 0 && ok; ) {
            fmpq_mpoly_set(suf + j, suf + j + 1, N->ctx);
            ok = nr_mul_piece_pow(N, suf + j, L[j].i, L[j].e);
        }
    }
    for (size_t t = 0; t < k && ok; t++) {
        if (rf_is_zero(N, T + t)) continue;
        fmpq_mpoly_set(acc, T[t].c, N->ctx);
        for (size_t j = 0; j < T[t].nnu && ok; j++)
            ok = nr_mul_piece_pow(N, acc, T[t].nu[j].i, T[t].nu[j].e);
        if (simple) {
            size_t jt = nl;                                   /* the term's piece in L */
            for (size_t j = 0; j < nl && T[t].nde == 1; j++)
                if (L[j].i == T[t].de[0].i) { jt = j; break; }
            if (jt == nl) nr_mul(N, acc, acc, pre + nl);      /* no denominator: all of L */
            else {
                nr_mul(N, acc, acc, pre + jt);
                nr_mul(N, acc, acc, suf + jt + 1);
                ok = nr_mul_piece_pow(N, acc, L[jt].i, L[jt].e - T[t].de[0].e);
            }
        } else {
            for (size_t j = 0; j < nl && ok; j++)
                ok = nr_mul_piece_pow(N, acc, L[j].i, L[j].e - pe_get(T[t].de, T[t].nde, L[j].i));
        }
        bool skipped = fmpq_mpoly_is_one(T[t].u, N->ctx);
        for (size_t j = 0; j < nu && ok; j++) {
            if (!skipped && fmpq_mpoly_equal(U + j, T[t].u, N->ctx)) { skipped = true; continue; }
            nr_mul(N, acc, acc, U + j);
        }
        if (ok) {
            fmpq_mpoly_add(S, S, acc, N->ctx);
            ok = fmpq_mpoly_length(S, N->ctx) <= NR_TERM_BUDGET;
        }
    }
    if (ok) {
        if (fmpq_mpoly_is_zero(S, N->ctx)) {
            rf_set_zero(N, dst);
        } else {
            ok = rf_set_poly(N, dst, S);
            if (ok) {
                for (size_t j = 0; j < nl; j++) pe_add(&dst->de, &dst->nde, &dst->cde, L[j].i, L[j].e);
                for (size_t j = 0; j < nu; j++) nr_mul(N, dst->u, dst->u, U + j);
                ok = rf_fold_u(N, dst);
            }
        }
    }
    fmpq_mpoly_clear(S, N->ctx); fmpq_mpoly_clear(acc, N->ctx);
    if (pre) {
        for (size_t j = 0; j <= nl; j++) { fmpq_mpoly_clear(pre + j, N->ctx); fmpq_mpoly_clear(suf + j, N->ctx); }
        free(pre); free(suf);
    }
    for (size_t j = 0; j < nu; j++) fmpq_mpoly_clear(U + j, N->ctx);
    free(U); free(L);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Expr -> rational function                                          */
/* ------------------------------------------------------------------ */

static long nr_atom_index(const Nr* N, const Expr* a) {
    for (size_t i = 0; i < N->na; i++) if (expr_eq(N->atoms[i], (Expr*)a)) return (long)i;
    return -1;
}

static long nr_var_index(const Nr* N, const char* s) {
    for (slong i = 0; i < N->ngens; i++) if (strcmp(N->var[i], s) == 0) return (long)i;
    return -1;
}

/* b^(p/q) with p < 0, written as (b^(-p/q))^-1: the positive radical. */
static Expr* nr_positive_radical(const Expr* e) {
    const Expr* x = e->data.function.args[1];
    const Expr* p = x->data.function.args[0];
    const Expr* q = x->data.function.args[1];
    if (p->type != EXPR_INTEGER || q->type != EXPR_INTEGER || p->data.integer >= 0) return NULL;
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy(e->data.function.args[0]),
                   make_rational(-p->data.integer, q->data.integer) }, 2);
}

static int nr_from_fmpq(Nr* N, NrRf* r, const fmpq_t q) {
    rf_set_zero(N, r);
    fmpq_mpoly_set_fmpq(r->c, q, N->ctx);
    return 1;
}

static int nr_eval(Nr* N, const Expr* e, NrRf* r);

static int nr_eval_atom(Nr* N, const Expr* e, NrRf* r) {
    long at = nr_atom_index(N, e);
    if (at < 0) return 0;
    rf_set_zero(N, r);
    fmpq_mpoly_set(r->c, N->aimg + at, N->ctx);
    return 1;
}

static int nr_eval(Nr* N, const Expr* e, NrRf* r) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: {
            fmpq_t q; fmpq_init(q);
            int ok = fb_expr_to_fmpq(e, q) && nr_from_fmpq(N, r, q);
            fmpq_clear(q);
            return ok;
        }
        case EXPR_SYMBOL: {
            long v = nr_var_index(N, e->data.symbol.name);
            if (v < 0) return 0;
            fmpq_mpoly_t p; fmpq_mpoly_init(p, N->ctx);
            fmpq_mpoly_gen(p, v, N->ctx);
            int ok = rf_set_poly(N, r, p);
            fmpq_mpoly_clear(p, N->ctx);
            return ok;
        }
        case EXPR_FUNCTION: break;
        default: return 0;
    }
    if (head_is(e, SYM_Rational)) {
        fmpq_t q; fmpq_init(q);
        int ok = fb_expr_to_fmpq(e, q) && nr_from_fmpq(N, r, q);
        fmpq_clear(q);
        return ok;
    }
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        /* a + b I */
        fmpq_t a, b; fmpq_init(a); fmpq_init(b);
        Expr* iu = make_complex(expr_new_integer(0), expr_new_integer(1));
        long at = nr_atom_index(N, iu);
        expr_free(iu);
        int ok = at >= 0 && fb_expr_to_fmpq(e->data.function.args[0], a) &&
                 fb_expr_to_fmpq(e->data.function.args[1], b);
        if (ok) {
            rf_set_zero(N, r);
            fmpq_mpoly_scalar_mul_fmpq(r->c, N->aimg + at, b, N->ctx);
            fmpq_mpoly_add_fmpq(r->c, r->c, a, N->ctx);
        }
        fmpq_clear(a); fmpq_clear(b);
        return ok;
    }
    bool plus = head_is(e, SYM_Plus);
    if (plus || head_is(e, SYM_Times)) {
        size_t k = e->data.function.arg_count;
        if (k == 0) return 0;
        if (!plus) {
            int ok = nr_eval(N, e->data.function.args[0], r);
            NrRf t; rf_init(N, &t);
            for (size_t i = 1; i < k && ok; i++) {
                ok = nr_eval(N, e->data.function.args[i], &t);
                if (ok) rf_mul_into(N, r, &t);
            }
            rf_clear(N, &t);
            return ok;
        }
        NrRf* T = malloc(k * sizeof(NrRf));
        int ok = T != NULL;
        size_t made = 0;
        for (size_t i = 0; i < k && ok; i++) {
            rf_init(N, T + i); made++;
            ok = nr_eval(N, e->data.function.args[i], T + i);
        }
        if (ok) ok = rf_sum(N, r, T, made);
        for (size_t i = 0; i < made; i++) rf_clear(N, T + i);
        free(T);
        return ok;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER) {
            if (!nr_eval(N, e->data.function.args[0], r)) return 0;
            return rf_pow(N, r, (slong)x->data.integer);
        }
        if (head_is(x, SYM_Rational)) {
            Expr* pos = nr_positive_radical(e);
            if (pos) {
                int ok = nr_eval_atom(N, pos, r) && rf_pow(N, r, -1);
                expr_free(pos);
                return ok;
            }
            return nr_eval_atom(N, e, r);
        }
        return 0;
    }
    if (head_is(e, SYM_Root)) return nr_eval_atom(N, e, r);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Setup: atoms, variables, field                                     */
/* ------------------------------------------------------------------ */

/* Root[f, k]^m with |m| below the degree of f: written in the field's power
 * basis unchanged. */
static bool nr_root_power_ok(const Expr* e) {
    const Expr* x = e->data.function.args[1];
    int64_t m = x->data.integer < 0 ? -x->data.integer : x->data.integer;
    if (m < 2) return true;
    Expr* mp = flint_qqbar_gen_minpoly_coeffs(e->data.function.args[0]);
    int64_t d = (mp && mp->type == EXPR_FUNCTION) ? (int64_t)mp->data.function.arg_count - 1 : 0;
    if (mp) expr_free(mp);
    return d > 0 && m < d;
}

static int nr_scan(const Expr* e, Nr* N, size_t* acap) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return 1;
        case EXPR_SYMBOL: {
            if (nr_var_index(N, e->data.symbol.name) >= 0) return 1;
            if (N->ngens == NR_MAX_VARS) return 0;
            N->var[N->ngens++] = e->data.symbol.name;
            return 1;
        }
        case EXPR_FUNCTION: break;
        default: return 0;
    }
    if (head_is(e, SYM_Rational)) return 1;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!nr_scan(e->data.function.args[i], N, acap)) return 0;
        return 1;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        /* Together / Cancel do not reduce a Root power modulo its polynomial
         * (the reference keeps Root[f, k]^5 as written); the field arithmetic
         * here would, so such an input is left to the generator paths. */
        if (head_is(e->data.function.args[0], SYM_Root) && !nr_root_power_ok(e)) return 0;
        return nr_scan(e->data.function.args[0], N, acap);
    }
    if (!nr_is_alg_leaf(e)) return 0;
    if (head_is(e, SYM_Power)) {
        Expr* pos = nr_positive_radical(e);
        if (pos) {
            int ok = fb_collect_atoms(pos, &N->atoms, &N->na, acap);
            expr_free(pos);
            return ok;
        }
    }
    return fb_collect_atoms(e, &N->atoms, &N->na, acap);
}

static int nr_cmp_name(const void* a, const void* b) {
    return strcmp(*(const char* const*)a, *(const char* const*)b);
}

static void nr_clear(Nr* N) {
    if (N->ready) {
        for (size_t i = 0; i < N->np; i++) {
            fmpq_mpoly_clear(N->pc + i, N->ctx);
            fmpq_mpoly_clear(N->wl + i, N->ctx);
        }
        for (size_t i = 0; i < N->na; i++) fmpq_mpoly_clear(N->aimg + i, N->ctx);
        fmpq_mpoly_clear(N->M, N->ctx);
        gr_ctx_clear(N->gr);
        nf_clear(N->nf);
        fb_render_clear(&N->R);
    }
    if (N->tpow) {
        for (slong j = 0; j < N->ntpow; j++) fmpq_poly_clear(N->tpow + j);
        free(N->tpow);
    }
    free(N->pc); free(N->wl); free(N->aimg); free(N->cop);
    for (size_t i = 0; i < N->na; i++) expr_free(N->atoms[i]);
    free(N->atoms);
    if (N->mp) expr_free(N->mp);
    if (N->theta) expr_free(N->theta);
    memset(N, 0, sizeof *N);
}

/* coordinates of an image (AlgebraicNumber[theta, {..}] or a rational) as a
 * tau polynomial */
static int nr_image_poly(Nr* N, const Expr* im, fmpq_mpoly_t out) {
    fmpq_mpoly_zero(out, N->ctx);
    fmpq_t q; fmpq_init(q);
    int ok = 1;
    if (fb_expr_to_fmpq(im, q)) {
        fmpq_mpoly_set_fmpq(out, q, N->ctx);
    } else if (im->type == EXPR_FUNCTION && head_is(im, SYM_AlgebraicNumber) &&
               im->data.function.arg_count == 2 && expr_eq(im->data.function.args[0], N->theta) &&
               head_is(im->data.function.args[1], SYM_List)) {
        const Expr* L = im->data.function.args[1];
        ulong ex[NR_MAX_VARS + 1]; memset(ex, 0, sizeof ex);
        for (size_t j = 0; j < L->data.function.arg_count && ok; j++) {
            if (!fb_expr_to_fmpq(L->data.function.args[j], q)) { ok = 0; break; }
            if (fmpq_is_zero(q)) continue;
            ex[N->tau] = (ulong)j;
            fmpq_mpoly_set_coeff_fmpq_ui(out, q, ex, N->ctx);
        }
        nr_reduce(N, out);
    } else ok = 0;
    fmpq_clear(q);
    return ok;
}

/* The generators of an explicit Extension -> alpha / {alpha, ...}: each a
 * constant algebraic number, joined to the atom set. */
static int nr_scan_ext(const Expr* ext, Nr* N, size_t* acap) {
    if (!ext) return 1;
    size_t n = 1;
    const Expr* const* it = &ext;
    if (head_is(ext, SYM_List)) {
        n = ext->data.function.arg_count;
        it = (const Expr* const*)ext->data.function.args;
    }
    for (size_t i = 0; i < n; i++) {
        const Expr* a = it[i];
        if (nr_is_ratnum(a)) continue;
        if (nr_has_symbol(a) || !nr_is_alg_leaf(a)) return 0;
        if (!fb_collect_atoms(a, &N->atoms, &N->na, acap)) return 0;
    }
    return 1;
}

static int nr_setup(Nr* N, const Expr* e, const Expr* ext) {
    memset(N, 0, sizeof *N);
    size_t acap = 0;
    if (!nr_scan(e, N, &acap) || !nr_scan_ext(ext, N, &acap) || N->na == 0 || N->ngens == 0)
        return 0;
    qsort(N->var, (size_t)N->ngens, sizeof(char*), nr_cmp_name);
    const Expr* imgs = fb_field_images(N->atoms, N->na);
    if (!imgs || imgs->type != EXPR_FUNCTION || imgs->data.function.arg_count != N->na) return 0;
    for (size_t i = 0; i < N->na && !N->theta; i++) {
        const Expr* im = imgs->data.function.args[i];
        if (im->type == EXPR_FUNCTION && head_is(im, SYM_AlgebraicNumber) && im->data.function.arg_count == 2)
            N->theta = expr_copy(im->data.function.args[0]);
    }
    if (!N->theta) return 0;
    N->mp = flint_qqbar_gen_minpoly_coeffs(N->theta);
    if (!N->mp || N->mp->type != EXPR_FUNCTION || N->mp->data.function.arg_count < 3) return 0;
    N->n = (slong)N->mp->data.function.arg_count - 1;
    if (N->n > NR_DEGREE_CAP) return 0;
    {
        const Expr* lc = N->mp->data.function.args[N->n];
        if (!(lc->type == EXPR_INTEGER && lc->data.integer == 1)) return 0;   /* not an algebraic integer */
    }
    N->tau = N->ngens;
    N->ctx = fb_ctx_lex(N->ngens + 1);
    if (!N->ctx) return 0;
    fmpq_mpoly_init(N->M, N->ctx);
    if (!fb_field_build_minpoly(N->M, N->mp, N->ctx, N->ngens + 1, N->tau)) {
        fmpq_mpoly_clear(N->M, N->ctx); return 0;
    }
    {
        fmpq_poly_t m; fmpq_poly_init(m);
        fmpq_t q; fmpq_init(q);
        for (slong j = 0; j <= N->n; j++) {
            fb_expr_to_fmpq(N->mp->data.function.args[j], q);
            fmpq_poly_set_coeff_fmpq(m, j, q);
        }
        nf_init(N->nf, m);
        gr_ctx_init_nf(N->gr, m);
        /* tau^j mod m, j = n .. 2n - 2 */
        N->ntpow = N->n - 1;
        N->tpow = N->ntpow > 0 ? malloc((size_t)N->ntpow * sizeof(fmpq_poly_struct)) : NULL;
        if (N->tpow) {
            fmpq_poly_t r, qq; fmpq_poly_init(r); fmpq_poly_init(qq);
            fmpq_poly_set_coeff_si(r, N->n - 1, 1);                  /* tau^(n-1) */
            for (slong j = 0; j < N->ntpow; j++) {
                fmpq_poly_shift_left(r, r, 1);
                fmpq_poly_rem(r, r, m);
                fmpq_poly_init(N->tpow + j);
                fmpq_poly_set(N->tpow + j, r);
            }
            fmpq_poly_clear(r); fmpq_poly_clear(qq);
        }
        fmpq_clear(q); fmpq_poly_clear(m);
    }
    N->aimg = malloc(N->na * sizeof(fmpq_mpoly_struct));
    for (size_t i = 0; i < N->na; i++) fmpq_mpoly_init(N->aimg + i, N->ctx);
    fb_render_init(&N->R, N->atoms, imgs->data.function.args, N->na, N->theta, N->n);
    N->ready = 1;
    if (!N->R.ready) return 0;                 /* the caller's atoms do not span K */
    for (size_t i = 0; i < N->na; i++)
        if (!nr_image_poly(N, imgs->data.function.args[i], N->aimg + i)) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Cancellation and the denominator basis                             */
/* ------------------------------------------------------------------ */

/* g = gcd(P_a, P_b) over K, monic.  1 when g is non-constant (g set), 0 when
 * coprime or undecided. */
static int nr_kron_divides(Nr* N, fmpq_mpoly_t q, const fmpq_mpoly_t A, const fmpq_mpoly_t G);

/* Total degree in the variables (tau excluded), stopping past `stop`. */
static slong nr_xdeg(Nr* N, const fmpq_mpoly_t p, slong stop) {
    ulong e[NR_MAX_VARS + 1];
    slong best = 0;
    for (slong i = 0; i < fmpq_mpoly_length(p, N->ctx) && best <= stop; i++) {
        fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
        slong d = 0;
        for (slong v = 0; v < N->ngens; v++) d += (slong)e[v];
        if (d > best) best = d;
    }
    return best;
}

static int nr_gcd(Nr* N, fmpq_mpoly_t g, size_t a, size_t b) {
    if (a == b) return 0;
    if (nr_coprime_known(N, a, b)) return 0;
    /* the exact cheap cases first (a Together of 17 linear Gaussian pieces
     * made 136 modular gcds): two distinct monic linear pieces are coprime,
     * and a linear piece shares a factor with another exactly when it
     * divides it */
    {
        slong da = nr_xdeg(N, N->pc + a, 1), db = nr_xdeg(N, N->pc + b, 1);
        if (da == 1 && db == 1) { nr_coprime_note(N, a, b); return 0; }
        if (da == 1 || db == 1) {
            size_t l = da == 1 ? a : b, o = da == 1 ? b : a;
            fmpq_mpoly_t q; fmpq_mpoly_init(q, N->ctx);
            int d = nr_kron_divides(N, q, N->pc + o, N->pc + l);
            fmpq_mpoly_clear(q, N->ctx);
            if (d == 1) { fmpq_mpoly_set(g, N->pc + l, N->ctx); return 1; }
            if (d == 0) { nr_coprime_note(N, a, b); return 0; }
        }
    }
    int got = fb_field_mpoly_gcd(g, N->pc + a, N->pc + b, N->M, N->ctx, N->mp,
                                 N->ngens, N->tau, N->n);
    if (!got || nr_is_const(N, g)) { nr_coprime_note(N, a, b); return 0; }
    return 1;
}

/* ---- exact division over K by Kronecker substitution ----
 *
 * x_i -> z^(prod_{j > i} D_j), D_j = deg_{x_j}(A) + 1, maps K[x] into K[z]
 * injectively on the polynomials of those degree bounds; the lex-leading
 * x-monomial goes to the leading z-power, so a monic G stays monic and the
 * division needs no inverse.  The quotient is read back and multiplied out
 * against A (the substitution can wrap a non-divisor onto a divisor), so a
 * success is certified.  The generic ideal division over {G, M} it replaces
 * was nearly all of a 3-variable Gaussian Together (59 s). */
#define NR_KRON_MAX 400000

static int nr_kron_to_grpoly(Nr* N, gr_poly_t out, const fmpq_mpoly_t p, const slong* stride) {
    slong len = fmpq_mpoly_length(p, N->ctx);
    ulong e[NR_MAX_VARS + 1];
    fmpq_t c; fmpq_init(c);
    fmpq_poly_t co; fmpq_poly_init(co);
    nf_elem_t el; nf_elem_init(el, N->nf);
    int ok = gr_poly_zero(out, N->gr) == GR_SUCCESS;
    slong i = 0;
    while (i < len && ok) {
        fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
        ulong g0[NR_MAX_VARS];
        memcpy(g0, e, (size_t)N->ngens * sizeof(ulong));
        slong z = 0;
        for (slong v = 0; v < N->ngens; v++) z += (slong)e[v] * stride[v];
        fmpq_poly_zero(co);
        while (i < len) {
            fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
            if (memcmp(g0, e, (size_t)N->ngens * sizeof(ulong)) != 0) break;
            fmpq_mpoly_get_term_coeff_fmpq(c, p, i, N->ctx);
            fmpq_poly_set_coeff_fmpq(co, (slong)e[N->tau], c);
            i++;
        }
        nf_elem_set_fmpq_poly(el, co, N->nf);
        ok = gr_poly_set_coeff_scalar(out, z, (gr_srcptr)el, N->gr) == GR_SUCCESS;
    }
    nf_elem_clear(el, N->nf); fmpq_poly_clear(co); fmpq_clear(c);
    return ok;
}

static int nr_kron_from_grpoly(Nr* N, fmpq_mpoly_t out, const gr_poly_t q, const slong* stride,
                               const slong* D) {
    fmpq_mpoly_zero(out, N->ctx);
    nf_elem_t el; nf_elem_init(el, N->nf);
    fmpq_poly_t co; fmpq_poly_init(co);
    fmpq_t c; fmpq_init(c);
    ulong e[NR_MAX_VARS + 1];
    int ok = 1;
    for (slong z = gr_poly_length(q, N->gr) - 1; z >= 0 && ok; z--) {
        if (gr_poly_get_coeff_scalar((gr_ptr)el, q, z, N->gr) != GR_SUCCESS) { ok = 0; break; }
        if (nf_elem_is_zero(el, N->nf)) continue;
        slong r = z;
        for (slong v = 0; v < N->ngens; v++) { e[v] = (ulong)(r / stride[v]); r %= stride[v]; }
        for (slong v = 0; v < N->ngens; v++) if ((slong)e[v] >= D[v]) ok = 0;
        nf_elem_get_fmpq_poly(co, el, N->nf);
        for (slong k = 0; k <= fmpq_poly_degree(co); k++) {
            fmpq_poly_get_coeff_fmpq(c, co, k);
            if (fmpq_is_zero(c)) continue;
            e[N->tau] = (ulong)k;
            fmpq_mpoly_push_term_fmpq_ui(out, c, e, N->ctx);
        }
    }
    fmpq_mpoly_sort_terms(out, N->ctx);
    fmpq_mpoly_combine_like_terms(out, N->ctx);
    fmpq_clear(c); fmpq_poly_clear(co); nf_elem_clear(el, N->nf);
    return ok;
}

/* q = A / G over K for a monic G; 1 when G divides A (q set), 0 when it does
 * not, -1 when the dense image would be too large (caller falls back). */
static int nr_kron_divides(Nr* N, fmpq_mpoly_t q, const fmpq_mpoly_t A, const fmpq_mpoly_t G) {
    slong D[NR_MAX_VARS], stride[NR_MAX_VARS];
    double size = 1;
    for (slong v = 0; v < N->ngens; v++) {
        slong da = fmpq_mpoly_degree_si(A, v, N->ctx), dg = fmpq_mpoly_degree_si(G, v, N->ctx);
        if (dg > da) return 0;
        D[v] = (da < 0 ? 0 : da) + 1;
        size *= (double)D[v];
    }
    /* the dense image costs len(A) len(G) field products: worth it only when
     * it is not much larger than the sparse operands (a sparse 3-variable
     * piece is better left to the ideal division) */
    if (size > NR_KRON_MAX || size > 32.0 * (double)(fmpq_mpoly_length(A, N->ctx) + 1)) return -1;
    slong s = 1;
    for (slong v = N->ngens - 1; v >= 0; v--) { stride[v] = s; s *= D[v]; }
    gr_poly_t a, g, qq, rr;
    gr_poly_init(a, N->gr); gr_poly_init(g, N->gr); gr_poly_init(qq, N->gr); gr_poly_init(rr, N->gr);
    int res = 0;
    if (nr_kron_to_grpoly(N, a, A, stride) && nr_kron_to_grpoly(N, g, G, stride) &&
        gr_poly_divrem(qq, rr, a, g, N->gr) == GR_SUCCESS && gr_poly_is_zero(rr, N->gr) == T_TRUE &&
        nr_kron_from_grpoly(N, q, qq, stride, D)) {
        /* certify: q G == A over K */
        fmpq_mpoly_t chk; fmpq_mpoly_init(chk, N->ctx);
        nr_mul(N, chk, q, G);
        res = fmpq_mpoly_equal(chk, A, N->ctx);
        fmpq_mpoly_clear(chk, N->ctx);
    }
    gr_poly_clear(a, N->gr); gr_poly_clear(g, N->gr); gr_poly_clear(qq, N->gr); gr_poly_clear(rr, N->gr);
    return res;
}

/* P_a = g * q: the quotient (monic, since both are), interned. */
static int nr_split(Nr* N, size_t a, const fmpq_mpoly_t g, size_t* gi, size_t* qi) {
    fmpq_mpoly_t q; fmpq_mpoly_init(q, N->ctx);
    int ok = nr_kron_divides(N, q, N->pc + a, g);
    if (ok < 0) ok = fb_field_mpoly_divides(q, N->pc + a, g, N->M, N->ctx);
    if (ok) {
        /* the cofactor keeps the written leading coefficient: wl P = (g)(wl q)
         * with g monic, so q can still be shown as written */
        fmpq_mpoly_t wl; fmpq_mpoly_init(wl, N->ctx);
        fmpq_mpoly_set(wl, N->wl + a, N->ctx);
        *gi = nr_piece(N, g, NULL);
        if (nr_is_const(N, q)) *qi = (size_t)-1;      /* g == P_a */
        else *qi = nr_piece(N, q, wl);
        fmpq_mpoly_clear(wl, N->ctx);
    }
    fmpq_mpoly_clear(q, N->ctx);
    return ok;
}

/* Replace piece a (exponent e) in list L by g^e q^e. */
static void pe_replace(NrPE** L, size_t* n, size_t* cap, size_t a, size_t gi, size_t qi) {
    slong e = pe_get(*L, *n, a);
    if (!e) return;
    pe_add(L, n, cap, a, -e);
    pe_add(L, n, cap, gi, e);
    if (qi != (size_t)-1) pe_add(L, n, cap, qi, e);
}

/* Cancel every numerator piece against every denominator piece. */
static int nr_cancel(Nr* N, NrRf* r) {
    fmpq_mpoly_t g; fmpq_mpoly_init(g, N->ctx);
    int ok = 1;
    bool changed = true;
    while (changed && ok) {
        changed = false;
        for (size_t i = 0; i < r->nnu && !changed && ok; i++)
            for (size_t j = 0; j < r->nde && !changed && ok; j++) {
                size_t a = r->nu[i].i, b = r->de[j].i;
                if (a == b) {
                    slong m = r->nu[i].e < r->de[j].e ? r->nu[i].e : r->de[j].e;
                    pe_add(&r->nu, &r->nnu, &r->cnu, a, -m);
                    pe_add(&r->de, &r->nde, &r->cde, a, -m);
                    changed = true;
                    break;
                }
                if (!nr_gcd(N, g, a, b)) continue;
                size_t ga, qa, gb, qb;
                if (!nr_split(N, a, g, &ga, &qa) || !nr_split(N, b, g, &gb, &qb)) { ok = 0; break; }
                pe_replace(&r->nu, &r->nnu, &r->cnu, a, ga, qa);
                pe_replace(&r->de, &r->nde, &r->cde, b, gb, qb);
                changed = true;
            }
    }
    fmpq_mpoly_clear(g, N->ctx);
    return ok;
}

/* The denominator over a squarefree, gcd-free basis. */
static int nr_den_basis(Nr* N, NrRf* r) {
    fmpq_mpoly_t g, d; fmpq_mpoly_init(g, N->ctx); fmpq_mpoly_init(d, N->ctx);
    int ok = 1;
    bool changed = true;
    /* squarefree: split each piece by its gcd with a partial derivative */
    while (changed && ok) {
        changed = false;
        for (size_t j = 0; j < r->nde && !changed && ok; j++) {
            size_t a = r->de[j].i;
            for (slong v = 0; v < N->ngens && !changed; v++) {
                if (fmpq_mpoly_degree_si(N->pc + a, v, N->ctx) < 2) continue;
                fmpq_mpoly_derivative(d, N->pc + a, v, N->ctx);
                if (fmpq_mpoly_is_zero(d, N->ctx)) continue;
                if (!fb_field_mpoly_gcd(g, N->pc + a, d, N->M, N->ctx, N->mp, N->ngens, N->tau, N->n))
                    continue;
                if (nr_is_const(N, g)) continue;
                size_t gi, qi;
                if (!nr_split(N, a, g, &gi, &qi)) { ok = 0; break; }
                pe_replace(&r->de, &r->nde, &r->cde, a, gi, qi);
                changed = true;
            }
        }
    }
    /* gcd-free: split any pair with a common factor */
    changed = true;
    while (changed && ok) {
        changed = false;
        for (size_t i = 0; i < r->nde && !changed && ok; i++)
            for (size_t j = i + 1; j < r->nde && !changed && ok; j++) {
                size_t a = r->de[i].i, b = r->de[j].i;
                if (!nr_gcd(N, g, a, b)) continue;
                size_t ga, qa, gb, qb;
                if (!nr_split(N, a, g, &ga, &qa) || !nr_split(N, b, g, &gb, &qb)) { ok = 0; break; }
                pe_replace(&r->de, &r->nde, &r->cde, a, ga, qa);
                pe_replace(&r->de, &r->nde, &r->cde, b, gb, qb);
                changed = true;
            }
    }
    fmpq_mpoly_clear(g, N->ctx); fmpq_mpoly_clear(d, N->ctx);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Presentation                                                       */
/* ------------------------------------------------------------------ */

/* The product-basis coordinates of a constant (tau polynomial): w = Pinv v. */
static void nr_coords(Nr* N, fmpq_mat_t w, const fmpq_mpoly_t c) {
    fmpq_mat_t v; fmpq_mat_init(v, N->n, 1);
    ulong e[NR_MAX_VARS + 1];
    fmpq_t q; fmpq_init(q);
    for (slong i = 0; i < fmpq_mpoly_length(c, N->ctx); i++) {
        fmpq_mpoly_get_term_exp_ui(e, c, i, N->ctx);
        fmpq_mpoly_get_term_coeff_fmpq(q, c, i, N->ctx);
        fmpq_set(fmpq_mat_entry(v, (slong)e[N->tau], 0), q);
    }
    fmpq_mat_mul(w, N->R.Pinv, v);
    fmpq_clear(q); fmpq_mat_clear(v);
}

/* A polynomial split into its x-monomial groups, each with product-basis
 * coordinates. */
typedef struct {
    slong ng;
    ulong* ex;          /* ng * ngens exponents */
    fmpq_mat_struct* w; /* ng columns (n x 1)   */
} NrGroups;

static void nr_groups_clear(Nr* N, NrGroups* G) {
    (void)N;
    for (slong i = 0; i < G->ng; i++) fmpq_mat_clear(G->w + i);
    free(G->w); free(G->ex);
    memset(G, 0, sizeof *G);
}

static void nr_groups(Nr* N, NrGroups* G, const fmpq_mpoly_t p) {
    memset(G, 0, sizeof *G);
    slong len = fmpq_mpoly_length(p, N->ctx);
    G->ex = malloc((size_t)(len ? len : 1) * (size_t)N->ngens * sizeof(ulong));
    G->w = malloc((size_t)(len ? len : 1) * sizeof(fmpq_mat_struct));
    ulong e[NR_MAX_VARS + 1];
    fmpq_mpoly_t c; fmpq_mpoly_init(c, N->ctx);
    fmpq_t q; fmpq_init(q);
    slong i = 0;
    while (i < len) {
        fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
        ulong* gx = G->ex + G->ng * N->ngens;
        memcpy(gx, e, (size_t)N->ngens * sizeof(ulong));
        fmpq_mpoly_zero(c, N->ctx);
        while (i < len) {
            fmpq_mpoly_get_term_exp_ui(e, p, i, N->ctx);
            if (memcmp(gx, e, (size_t)N->ngens * sizeof(ulong)) != 0) break;
            fmpq_mpoly_get_term_coeff_fmpq(q, p, i, N->ctx);
            ulong t[NR_MAX_VARS + 1]; memset(t, 0, sizeof t);
            t[N->tau] = e[N->tau];
            fmpq_mpoly_set_coeff_fmpq_ui(c, q, t, N->ctx);
            i++;
        }
        fmpq_mat_init(G->w + G->ng, N->n, 1);
        nr_coords(N, G->w + G->ng, c);
        G->ng++;
    }
    fmpq_clear(q); fmpq_mpoly_clear(c, N->ctx);
}

/* Number of nonzero coordinates (the rendered size). */
static slong nr_groups_size(Nr* N, const NrGroups* G) {
    slong s = 0;
    for (slong i = 0; i < G->ng; i++)
        for (slong m = 0; m < N->n; m++)
            if (!fmpq_is_zero(fmpq_mat_entry(G->w + i, m, 0))) s++;
    return s;
}

/* Nonzero coordinates on an irrational basis member (prod[m] != 1). */
static slong nr_groups_irr(Nr* N, const NrGroups* G) {
    slong s = 0;
    for (slong i = 0; i < G->ng; i++)
        for (slong m = 0; m < N->n; m++) {
            const Expr* pm = N->R.prod[m];
            if (pm->type == EXPR_INTEGER && pm->data.integer == 1) continue;
            if (!fmpq_is_zero(fmpq_mat_entry(G->w + i, m, 0))) s++;
        }
    return s;
}

/* s with s * (every coordinate) integral and jointly primitive, the sign
 * making the leading group's first nonzero coordinate positive. */
static void nr_groups_scale(Nr* N, const NrGroups* G, fmpq_t s) {
    fmpz_t L, g, t; fmpz_init(L); fmpz_init(g); fmpz_init(t);
    fmpz_one(L);
    for (slong i = 0; i < G->ng; i++)
        for (slong m = 0; m < N->n; m++) {
            const fmpq* c = fmpq_mat_entry(G->w + i, m, 0);
            if (!fmpq_is_zero(c)) fmpz_lcm(L, L, fmpq_denref(c));
        }
    for (slong i = 0; i < G->ng; i++)
        for (slong m = 0; m < N->n; m++) {
            const fmpq* c = fmpq_mat_entry(G->w + i, m, 0);
            if (fmpq_is_zero(c)) continue;
            fmpz_divexact(t, L, fmpq_denref(c));
            fmpz_mul(t, t, fmpq_numref(c));
            fmpz_gcd(g, g, t);
        }
    if (fmpz_is_zero(g)) fmpz_one(g);
    fmpq_set_fmpz_frac(s, L, g);
    if (G->ng > 0)
        for (slong m = 0; m < N->n; m++) {
            const fmpq* c = fmpq_mat_entry(G->w + 0, m, 0);
            if (fmpq_is_zero(c)) continue;
            if (fmpq_sgn(c) < 0) fmpq_neg(s, s);
            break;
        }
    fmpz_clear(L); fmpz_clear(g); fmpz_clear(t);
}

static Expr* nr_fmpq_expr(const fmpq_t q) { return fb_fmpq_to_expr(q); }

/* s * p as an (unevaluated) Plus over the product basis and x-monomials. */
static Expr* nr_groups_expr(Nr* N, const NrGroups* G, const fmpq_t s) {
    size_t cap = (size_t)(G->ng * N->n) + 1, nt = 0;
    Expr** terms = malloc(cap * sizeof(Expr*));
    fmpq_t c; fmpq_init(c);
    for (slong i = 0; i < G->ng; i++) {
        const ulong* gx = G->ex + i * N->ngens;
        for (slong m = 0; m < N->n; m++) {
            if (fmpq_is_zero(fmpq_mat_entry(G->w + i, m, 0))) continue;
            fmpq_mul(c, fmpq_mat_entry(G->w + i, m, 0), s);
            Expr* f[NR_MAX_VARS + 2];
            size_t nf = 0;
            f[nf++] = nr_fmpq_expr(c);
            const Expr* pm = N->R.prod[m];
            if (!(pm->type == EXPR_INTEGER && pm->data.integer == 1)) f[nf++] = expr_copy((Expr*)pm);
            for (slong v = 0; v < N->ngens; v++) {
                if (!gx[v]) continue;
                Expr* xv = expr_new_symbol(N->var[v]);
                f[nf++] = gx[v] == 1 ? xv
                    : expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ xv, expr_new_integer((int64_t)gx[v]) }, 2);
            }
            terms[nt++] = nf == 1 ? f[0] : expr_new_function(expr_new_symbol(SYM_Times), f, nf);
        }
    }
    fmpq_clear(c);
    Expr* r = nt == 0 ? expr_new_integer(0)
            : nt == 1 ? terms[0]
            : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return r;
}

/* A constant as an Expr over the product basis (unevaluated). */
static Expr* nr_const_expr(Nr* N, const fmpq_mpoly_t c) {
    NrGroups G; memset(&G, 0, sizeof G);
    G.ex = calloc((size_t)N->ngens + 1, sizeof(ulong));
    G.w = malloc(sizeof(fmpq_mat_struct));
    fmpq_mat_init(G.w, N->n, 1);
    nr_coords(N, G.w, c);
    G.ng = 1;
    fmpq_t one; fmpq_init(one); fmpq_one(one);
    Expr* r = nr_groups_expr(N, &G, one);
    fmpq_clear(one);
    nr_groups_clear(N, &G);
    return r;
}

/* Is c a rational multiple of I (its only product-basis coordinate on the
 * member I)?  A Gaussian unit is normalised away, as the reference does for
 * Gaussian coefficients (denominator pieces only; a numerator keeps its
 * written form): Together[1/(1 - I x) + 1/(1 + I x)] is
 * 2/((-I + x) (I + x)), while a non-unit (1 + (2 + I) x) keeps its form. */
static bool nr_is_gaussian_unit(Nr* N, const fmpq_mpoly_t c) {
    fmpq_mat_t w; fmpq_mat_init(w, N->n, 1);
    nr_coords(N, w, c);
    bool ok = false, other = false;
    for (slong m = 0; m < N->n; m++) {
        if (fmpq_is_zero(fmpq_mat_entry(w, m, 0))) continue;
        const Expr* pm = N->R.prod[m];
        if (pm->type == EXPR_FUNCTION && head_is(pm, SYM_Complex) &&
            pm->data.function.arg_count == 2 &&
            pm->data.function.args[0]->type == EXPR_INTEGER && pm->data.function.args[0]->data.integer == 0 &&
            pm->data.function.args[1]->type == EXPR_INTEGER && pm->data.function.args[1]->data.integer == 1)
            ok = true;
        else other = true;
    }
    fmpq_mat_clear(w);
    return ok && !other;
}

/* Does a basis member mention I? */
static bool nr_mentions_imag(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Complex)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nr_mentions_imag(e->data.function.args[i])) return true;
    return false;
}

/* Nonzero coordinates of the constant c on basis members that mention I. */
static slong nr_imag_coords(Nr* N, const fmpq_mpoly_t c, slong* total) {
    fmpq_mat_t w; fmpq_mat_init(w, N->n, 1);
    nr_coords(N, w, c);
    slong k = 0, t = 0;
    for (slong m = 0; m < N->n; m++) {
        if (fmpq_is_zero(fmpq_mat_entry(w, m, 0))) continue;
        t++;
        if (nr_mentions_imag(N->R.prod[m])) k++;
    }
    fmpq_mat_clear(w);
    if (total) *total = t;
    return k;
}

/* c := c / I when that removes every I from c's coordinates (c = I w with
 * w free of I).  I must be one of the field's atoms. */
static void nr_strip_imag_unit(Nr* N, fmpq_mpoly_t c) {
    long iu = -1;
    for (size_t i = 0; i < N->na && iu < 0; i++) {
        const Expr* a = N->atoms[i];
        if (head_is(a, SYM_Complex) && a->data.function.args[0]->type == EXPR_INTEGER &&
            a->data.function.args[0]->data.integer == 0 && a->data.function.args[1]->type == EXPR_INTEGER &&
            a->data.function.args[1]->data.integer == 1)
            iu = (long)i;
    }
    if (iu < 0) return;
    slong tc = 0, tw = 0;
    slong kc = nr_imag_coords(N, c, &tc);
    if (kc == 0 || kc != tc) return;               /* not I times an I-free constant */
    fmpq_mpoly_t w; fmpq_mpoly_init(w, N->ctx);
    nr_mul(N, w, c, N->aimg + iu);                 /* I c */
    fmpq_mpoly_neg(w, w, N->ctx);                  /* -I c = c / I */
    if (nr_imag_coords(N, w, &tw) == 0) fmpq_mpoly_swap(c, w, N->ctx);
    fmpq_mpoly_clear(w, N->ctx);
}

/* -p for a shown sum p (an unevaluated Plus of Times[c, ..] / number / atom
 * terms, as nr_groups_expr builds it): each coefficient negated. */
static Expr* nr_negate_sum(const Expr* p) {
    size_t n = p->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        const Expr* t = p->data.function.args[i];
        if (head_is(t, SYM_Times) && t->data.function.arg_count >= 1 &&
            nr_is_ratnum(t->data.function.args[0])) {
            size_t k = t->data.function.arg_count;
            Expr** f = malloc(k * sizeof(Expr*));
            f[0] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(t->data.function.args[0]) }, 2));
            for (size_t j = 1; j < k; j++) f[j] = expr_copy(t->data.function.args[j]);
            a[i] = expr_new_function(expr_new_symbol(SYM_Times), f, k);
            free(f);
        } else {
            a[i] = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)t) }, 2);
        }
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_Plus), a, n);
    free(a);
    return r;
}

/* The presented form of piece i: s * P_i * (wl_i or 1).  *scale receives the
 * total K-scalar applied (s * wl_i). */
static Expr* nr_piece_expr(Nr* N, size_t i, bool allow_written, bool den, fmpq_mpoly_t scale) {
    NrGroups Gm; nr_groups(N, &Gm, N->pc + i);
    fmpq_t s; fmpq_init(s);
    Expr* out = NULL;
    if (allow_written && !fmpq_mpoly_is_one(N->wl + i, N->ctx) && !(den && nr_is_gaussian_unit(N, N->wl + i))) {
        fmpq_mpoly_t pw; fmpq_mpoly_init(pw, N->ctx);
        nr_mul(N, pw, N->pc + i, N->wl + i);
        NrGroups Gw; nr_groups(N, &Gw, pw);
        /* the written form unless it is larger, or the same size with more
         * irrational coordinates (I (16 t - 100 t^4 ..) against the rational
         * 16 t - 100 t^4 ..) */
        slong sw = nr_groups_size(N, &Gw), sm = nr_groups_size(N, &Gm);
        if (sw < sm || (sw == sm && nr_groups_irr(N, &Gw) <= nr_groups_irr(N, &Gm))) {
            nr_groups_scale(N, &Gw, s);
            out = nr_groups_expr(N, &Gw, s);
            fmpq_mpoly_scalar_mul_fmpq(scale, N->wl + i, s, N->ctx);
        }
        nr_groups_clear(N, &Gw);
        fmpq_mpoly_clear(pw, N->ctx);
    }
    if (!out) {
        nr_groups_scale(N, &Gm, s);
        out = nr_groups_expr(N, &Gm, s);
        fmpq_mpoly_set_fmpq(scale, s, N->ctx);
    }
    fmpq_clear(s);
    nr_groups_clear(N, &Gm);
    return out;
}

static Expr* nr_pow_expr(Expr* b, slong e) {
    if (e == 1) return b;
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ b, expr_new_integer((int64_t)e) }, 2);
}

/* Write r (cancelled) as an Expr: factored (user level) or expanded num/den. */
static Expr* nr_present(Nr* N, NrRf* r, int factored) {
    if (rf_is_zero(N, r)) return expr_new_integer(0);
    size_t cap = r->nnu + r->nde + 4, na = 0;
    Expr** a = malloc(cap * sizeof(Expr*));
    fmpq_mpoly_t num_s, den_s, sc, t;
    fmpq_mpoly_init(num_s, N->ctx); fmpq_mpoly_init(den_s, N->ctx);
    fmpq_mpoly_init(sc, N->ctx); fmpq_mpoly_init(t, N->ctx);
    fmpq_mpoly_set(num_s, r->c, N->ctx);     /* value = num_s * shown_num / (den_s * shown_den) */
    fmpq_mpoly_set(den_s, r->u, N->ctx);
    int ok = 1;
    if (factored) {
        /* The numerator pieces as written (their written leading coefficient
         * kept: 1 + x + Sqrt[2] x over 1 + Sqrt[2], not the rationalised
         * -1 + Sqrt[2] + x) or monic.  Both are built; the written form is
         * kept unless only the monic one makes the scalar in front rational
         * (I (16 t - ..) would otherwise put an I into the numerator and a
         * compensating constant over the denominator). */
        Expr** na_w = malloc((r->nnu + 1) * sizeof(Expr*));
        Expr** na_m = malloc((r->nnu + 1) * sizeof(Expr*));
        fmpq_mpoly_t dw, dm, q;
        fmpq_mpoly_init(dw, N->ctx); fmpq_mpoly_init(dm, N->ctx); fmpq_mpoly_init(q, N->ctx);
        fmpq_mpoly_set(dw, den_s, N->ctx); fmpq_mpoly_set(dm, den_s, N->ctx);
        bool differ = false;
        for (size_t j = 0; j < r->nnu; j++) {
            for (int w = 0; w < 2; w++) {
                Expr* pe = nr_piece_expr(N, r->nu[j].i, w == 0, false, sc);
                fmpq_mpoly_set(t, sc, N->ctx); nr_const_pow(N, t, (ulong)r->nu[j].e);
                nr_mul(N, w == 0 ? dw : dm, w == 0 ? dw : dm, t);   /* shown = sc^e P^e */
                (w == 0 ? na_w : na_m)[j] = nr_pow_expr(pe, r->nu[j].e);
            }
            if (!expr_eq(na_w[j], na_m[j])) differ = true;
        }
        bool use_monic = false;
        /* With a written constant denominator u the written numerator is
         * kept: the monic one would absorb u (rationalise it), which the
         * reference never does. */
        if (differ && fmpq_mpoly_is_one(r->u, N->ctx)) {
            /* is num_s / dw rational?  num_s / dm? (num_s is still c here; the
             * denominator pieces below scale both alike) */
            bool rw = false, rm = false;
            if (nr_const_inv(N, q, dw)) { nr_mul(N, q, q, num_s); rw = nr_is_rational_const(N, q); }
            if (nr_const_inv(N, q, dm)) { nr_mul(N, q, q, num_s); rm = nr_is_rational_const(N, q); }
            use_monic = !rw && rm;
        }
        for (size_t j = 0; j < r->nnu; j++) {
            a[na++] = use_monic ? na_m[j] : na_w[j];
            expr_free(use_monic ? na_w[j] : na_m[j]);
        }
        fmpq_mpoly_set(den_s, use_monic ? dm : dw, N->ctx);
        fmpq_mpoly_clear(dw, N->ctx); fmpq_mpoly_clear(dm, N->ctx); fmpq_mpoly_clear(q, N->ctx);
        free(na_w); free(na_m);
        for (size_t j = 0; j < r->nde; j++) {
            Expr* pe = nr_piece_expr(N, r->de[j].i, true, true, sc);
            fmpq_mpoly_set(t, sc, N->ctx); nr_const_pow(N, t, (ulong)r->de[j].e);
            nr_mul(N, num_s, num_s, t);
            a[na++] = nr_pow_expr(pe, -r->de[j].e);
        }
    } else {
        for (int side = 0; side < 2 && ok; side++) {
            NrPE* L = side ? r->de : r->nu;
            size_t n = side ? r->nde : r->nnu;
            if (n == 0) continue;
            fmpq_mpoly_t P; fmpq_mpoly_init(P, N->ctx); fmpq_mpoly_one(P, N->ctx);
            for (size_t j = 0; j < n && ok; j++) ok = nr_mul_piece_pow(N, P, L[j].i, L[j].e);
            if (ok) {
                NrGroups G; nr_groups(N, &G, P);
                fmpq_t s; fmpq_init(s);
                nr_groups_scale(N, &G, s);
                Expr* pe = nr_groups_expr(N, &G, s);
                fmpq_mpoly_set_fmpq(sc, s, N->ctx);
                if (side) nr_mul(N, num_s, num_s, sc); else nr_mul(N, den_s, den_s, sc);
                a[na++] = side ? nr_pow_expr(pe, -1) : pe;
                fmpq_clear(s);
                nr_groups_clear(N, &G);
            }
            fmpq_mpoly_clear(P, N->ctx);
        }
    }
    /* The scalar sigma = num_s / den_s.  The only non-rational constant the
     * answer may show as a denominator is the one the input wrote, u (the
     * reference never rationalises it): sigma is shown as (sigma u) / u unless
     * it is rational or there is no written constant denominator.  The
     * inversions of the presentation's own scalings are internal. */
    if (ok) {
        if (!nr_const_inv(N, t, den_s)) ok = 0;
        else {
            nr_mul(N, t, t, num_s);                  /* sigma */
            /* a negative rational in front goes into the first numerator
             * sum of exponent 1, as the reference writes it:
             * (1 - Sqrt[2] x)/x, not -((-1 + Sqrt[2] x)/x) */
            if (nr_is_rational_const(N, t)) {
                fmpq_t q; fmpq_init(q);
                fmpq_mpoly_get_fmpq(q, t, N->ctx);
                size_t nnum = factored ? r->nnu : (r->nnu ? 1u : 0u);
                for (size_t j = 0; j < nnum && fmpq_sgn(q) < 0; j++) {
                    slong ej = factored ? r->nu[j].e : 1;
                    if (ej != 1 || !head_is(a[j], SYM_Plus)) continue;
                    Expr* ng = nr_negate_sum(a[j]);
                    expr_free(a[j]);
                    a[j] = ng;
                    fmpq_mpoly_neg(t, t, N->ctx);
                    break;
                }
                fmpq_clear(q);
            }
            if (nr_is_rational_const(N, t) || fmpq_mpoly_is_one(r->u, N->ctx)) {
                a[na++] = nr_const_expr(N, t);
            } else {
                /* u = q u0 with q its rational content: u0 is shown, q is
                 * folded into the numerator constant */
                fmpq_mpoly_t u0; fmpq_mpoly_init(u0, N->ctx);
                fmpq_mat_t w; fmpq_mat_init(w, N->n, 1);
                nr_coords(N, w, r->u);
                NrGroups G; memset(&G, 0, sizeof G);
                G.ex = calloc((size_t)N->ngens + 1, sizeof(ulong));
                G.w = w; G.ng = 1;
                fmpq_t s; fmpq_init(s);
                nr_groups_scale(N, &G, s);              /* s u0-coords integral, primitive */
                if (fmpq_sgn(s) < 0) fmpq_neg(s, s);     /* u keeps its written sign */
                G.w = NULL; G.ng = 0; free(G.ex);
                fmpq_mat_clear(w);
                fmpq_mpoly_scalar_mul_fmpq(u0, r->u, s, N->ctx);   /* u0 = s u */
                /* a Gaussian unit is not part of the written constant either:
                 * -577 I + 408 I Sqrt[2] is shown as -577 + 408 Sqrt[2] (times
                 * -I into the numerator), not with an I on both sides */
                nr_strip_imag_unit(N, u0);
                nr_mul(N, t, t, u0);                    /* sigma u0 */
                a[na++] = nr_const_expr(N, t);
                a[na++] = nr_pow_expr(nr_const_expr(N, u0), -1);
                fmpq_clear(s);
                fmpq_mpoly_clear(u0, N->ctx);
            }
        }
    }
    fmpq_mpoly_clear(num_s, N->ctx); fmpq_mpoly_clear(den_s, N->ctx);
    fmpq_mpoly_clear(sc, N->ctx); fmpq_mpoly_clear(t, N->ctx);
    if (!ok) {
        for (size_t i = 0; i < na; i++) expr_free(a[i]);
        free(a);
        return NULL;
    }
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), a, na));
    free(a);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Entry points                                                       */
/* ------------------------------------------------------------------ */

static Expr* nr_together_in(Nr* N, const Expr* e, int factored) {
    NrRf r; rf_init(N, &r);
    Expr* out = NULL;
    if (nr_eval(N, e, &r) && nr_cancel(N, &r) && (!factored || nr_den_basis(N, &r)))
        out = nr_present(N, &r, factored);
    rf_clear(N, &r);
    return out;
}

Expr* nfrat_together(const Expr* e, int factored, const Expr* ext) {
    if (!e || !nfrat_has_algebraic(e)) return NULL;
    Nr N;
    Expr* out = NULL;
    if (nr_setup(&N, e, ext)) out = nr_together_in(&N, e, factored);
    nr_clear(&N);
    return out;
}

/* Does `t` hold a denominator (a negative integer power) in its rational
 * structure? */
static bool nr_has_den(const Expr* t) {
    if (t->type != EXPR_FUNCTION) return false;
    if (head_is(t, SYM_Power) && t->data.function.arg_count == 2) {
        const Expr* x = t->data.function.args[1];
        if (x->type == EXPR_INTEGER) return x->data.integer < 0 || nr_has_den(t->data.function.args[0]);
        return false;
    }
    if (head_is(t, SYM_Plus) || head_is(t, SYM_Times)) {
        for (size_t i = 0; i < t->data.function.arg_count; i++)
            if (nr_has_den(t->data.function.args[i])) return true;
    }
    return false;
}

/* Cancel's structure, as in the reference: a sum is cancelled term by term; a
 * power of a sum is the power of the cancelled sum (Cancel[(1/(x+1) +
 * 1/(x+2))^2] stays as it is, Cancel[(1/(x+1) + (x^2-1)/(x-1))^2] is
 * (1 + x + 1/(1+x))^2); any other term with a denominator -- a product, a
 * single fraction -- is reduced to one fraction.  *changed is set when some
 * part was rewritten.  NULL on a decline. */
static Expr* nr_cancel_rec(Nr* N, const Expr* e, int factored, bool* changed) {
    if (head_is(e, SYM_Plus)) {
        size_t k = e->data.function.arg_count;
        Expr** a = malloc((k ? k : 1) * sizeof(Expr*));
        bool any = false;
        for (size_t i = 0; i < k; i++) {
            bool ch = false;
            a[i] = nr_cancel_rec(N, e->data.function.args[i], factored, &ch);
            if (!a[i]) {
                for (size_t j = 0; j < i; j++) expr_free(a[j]);
                free(a);
                return NULL;
            }
            any = any || ch;
        }
        Expr* out;
        if (!any) {
            for (size_t j = 0; j < k; j++) expr_free(a[j]);
            out = expr_copy((Expr*)e);
        } else {
            out = expr_new_function(expr_new_symbol(SYM_Plus), a, k);
            *changed = true;
        }
        free(a);
        return out;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER && head_is(e->data.function.args[0], SYM_Plus)) {
        bool ch = false;
        Expr* b = nr_cancel_rec(N, e->data.function.args[0], factored, &ch);
        if (!b) return NULL;
        if (!ch) { expr_free(b); return expr_copy((Expr*)e); }
        *changed = true;
        return expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ b, expr_copy(e->data.function.args[1]) }, 2);
    }
    if (!nr_has_den(e) || !nr_has_symbol(e)) return expr_copy((Expr*)e);
    Expr* r = nr_together_in(N, e, factored);
    if (r) *changed = true;
    return r;
}

Expr* nfrat_cancel(const Expr* e, int factored, const Expr* ext) {
    if (!e || !nfrat_has_algebraic(e)) return NULL;
    Nr N;
    Expr* out = NULL;
    if (nr_setup(&N, e, ext)) {
        bool changed = false;
        /* unchanged is an answer too: the input is already in Cancel's form */
        Expr* r = nr_cancel_rec(&N, e, factored, &changed);
        if (r) out = changed ? eval_and_free(r) : r;
    }
    nr_clear(&N);
    return out;
}

#else  /* !USE_FLINT */

Expr* nfrat_together(const Expr* e, int factored, const Expr* ext) { (void)e; (void)factored; (void)ext; return NULL; }
Expr* nfrat_cancel(const Expr* e, int factored, const Expr* ext) { (void)e; (void)factored; (void)ext; return NULL; }

#endif

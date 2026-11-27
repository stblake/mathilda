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
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"

#include <stdio.h>
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

/* ------------------------------------------------------------------ */
/*  Radicals of non-constant radicands                                 */
/* ------------------------------------------------------------------ */

/* Power[B, p/q] with a symbol in B. */
static bool nr_symrad(const Expr* e, const Expr** base, int64_t* p, int64_t* q) {
    if (e->type != EXPR_FUNCTION || !head_is(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    const Expr* x = e->data.function.args[1];
    if (!head_is(x, SYM_Rational) || x->data.function.arg_count != 2) return false;
    const Expr* xp = x->data.function.args[0];
    const Expr* xq = x->data.function.args[1];
    if (xp->type != EXPR_INTEGER || xq->type != EXPR_INTEGER) return false;
    if (!nr_has_symbol(e->data.function.args[0])) return false;
    if (base) *base = e->data.function.args[0];
    if (p) *p = xp->data.integer;
    if (q) *q = xq->data.integer;
    return true;
}

bool nfrat_has_symbolic_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (nr_symrad(e, NULL, NULL, NULL)) return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nfrat_has_symbolic_radical(e->data.function.args[i])) return true;
        return false;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER)
        return nfrat_has_symbolic_radical(e->data.function.args[0]);
    return false;
}

/* Together / Cancel read each radical of a non-constant radicand as a free
 * kernel, as the reference does: one generator g per radicand as written,
 * B^(p/q) -> g^(p Q/q) with Q the lcm of the radicand's indices (an integer
 * power of one principal radical is exact), and a bare-symbol radicand
 * eliminated, s -> g^Q (so 1/(x + Sqrt[x]) + 1/(x - Sqrt[x]) cancels the
 * common Sqrt[x]).  No relation g^Q = B is applied beyond that: the radicals
 * are otherwise independent atoms. */
#define NK_MAX 8
typedef struct {
    Expr*   base;      /* radicand as written (owned) */
    int64_t Q;
    int     elim;      /* bare-symbol radicand, eliminated */
    Expr*   sym;       /* fresh generator (owned) */
} NkGen;

typedef struct { NkGen g[NK_MAX]; int n; int bad; } NkTab;

static bool nk_reserved_symbol(const char* name) {
    SymbolDef* d = symtab_lookup(name);
    return d && (d->attributes & (ATTR_CONSTANT | ATTR_PROTECTED));
}

static void nk_collect(const Expr* e, NkTab* T) {
    if (T->bad || e->type != EXPR_FUNCTION) return;
    const Expr* b; int64_t p, q;
    if (nr_symrad(e, &b, &p, &q)) {
        if (q < 2) { T->bad = 1; return; }
        for (int i = 0; i < T->n; i++)
            if (expr_eq(T->g[i].base, (Expr*)b)) { T->g[i].Q = lcm(T->g[i].Q, q); return; }
        if (T->n == NK_MAX) { T->bad = 1; return; }
        NkGen* g = &T->g[T->n];
        g->base = expr_copy((Expr*)b);
        g->Q = q;
        g->elim = b->type == EXPR_SYMBOL && !nk_reserved_symbol(b->data.symbol.name);
        char name[48];
        snprintf(name, sizeof name, "~MthRadGen%02d", T->n + 1);   /* sorts after the variables */
        g->sym = expr_new_symbol(intern_symbol(name));
        T->n++;
        return;
    }
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times) ||
        (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
         e->data.function.args[1]->type == EXPR_INTEGER))
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nk_collect(e->data.function.args[i], T);
}

static Expr* nk_gen_power(const NkGen* g, int64_t k) {
    if (k == 1) return expr_copy(g->sym);
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(g->sym), expr_new_integer(k) }, 2);
}

static Expr* nk_subst(const Expr* e, const NkTab* T) {
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < T->n; i++)
            if (T->g[i].elim && strcmp(T->g[i].base->data.symbol.name, e->data.symbol.name) == 0)
                return nk_gen_power(&T->g[i], T->g[i].Q);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* b; int64_t p, q;
    if (nr_symrad(e, &b, &p, &q))
        for (int i = 0; i < T->n; i++)
            if (expr_eq(T->g[i].base, (Expr*)b)) return nk_gen_power(&T->g[i], p * (T->g[i].Q / q));
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = nk_subst(e->data.function.args[i], T);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

static Expr* nk_back(const Expr* e, const NkTab* T) {
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < T->n; i++)
            if (strcmp(T->g[i].sym->data.symbol.name, e->data.symbol.name) == 0)
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(T->g[i].base), make_rational(1, T->g[i].Q) }, 2);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = nk_back(e->data.function.args[i], T);
    Expr* r = expr_new_function(nk_back(e->data.function.head, T), a, n);
    free(a);
    return r;
}

static bool nr_field_degree_small(const Expr* e);   /* below, with FLINT */

Expr* nfrat_radical_kernels(const Expr* e, int factored, const Expr* ext, int cancel) {
    if (!e || !nfrat_has_symbolic_radical(e)) return NULL;
    /* building a large constant field costs more than the paths below */
    if (!nr_field_degree_small(e) || (ext && !nr_field_degree_small(ext))) return NULL;
    NkTab T; memset(&T, 0, sizeof T);
    nk_collect(e, &T);
    Expr* out = NULL;
    if (!T.bad && T.n > 0) {
        Expr* sub = nk_subst(e, &T);
        if (nfrat_has_algebraic(sub) && nfrat_has_denominator(sub) && !nfrat_has_unreduced_root_power(sub)) {
            Expr* r = cancel ? nfrat_cancel(sub, factored, ext) : nfrat_together(sub, factored, ext);
            if (r) {
                out = eval_and_free(nk_back(r, &T));
                expr_free(r);
            }
        }
        expr_free(sub);
    }
    for (int i = 0; i < T.n; i++) { expr_free(T.g[i].base); expr_free(T.g[i].sym); }
    return out;
}

#ifdef USE_FLINT

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_factor.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mpoly.h>
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
    fmpz_poly_struct* tpowz;  /* the same, integral (M monic integral), or NULL */
    /* the piece table: monic, reduced, distinct */
    fmpq_mpoly_struct* pc;
    fmpq_mpoly_struct* wl;    /* written leading coefficient (1 if none)        */
    size_t np, pcap;
    /* known-coprime pairs (i < j) */
    size_t* cop; size_t ncop, ccop;
    int ready;
    /* the radical function field's generator table (RnTab), or NULL */
    void* rad;
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
    if (N->tpowz && dt - N->n < N->ntpow) {
        /* M is monic and integral, so tau^j mod M is integral: the reduction
         * runs on p's integer part (p = content * zpoly) with no rational
         * rescaling per term, and the content is renormalised once. */
        const fmpz_mpoly_ctx_struct* zc = N->ctx->zctx;
        fmpz_mpoly_t out; fmpz_mpoly_init(out, zc);
        ulong e[NR_MAX_VARS + 1];
        fmpz_t v; fmpz_init(v);
        fmpz_mpoly_struct* zp = p->zpoly;
        slong len = fmpz_mpoly_length(zp, zc);
        fmpz_mpoly_fit_length(out, len, zc);
        for (slong i = 0; i < len; i++) {
            fmpz_mpoly_get_term_exp_ui(e, zp, i, zc);
            const fmpz* c = zp->coeffs + i;
            slong j = (slong)e[N->tau];
            if (j < N->n) { fmpz_mpoly_push_term_fmpz_ui(out, c, e, zc); continue; }
            const fmpz_poly_struct* rj = N->tpowz + (j - N->n);
            for (slong k = 0; k <= fmpz_poly_degree(rj); k++) {
                const fmpz* rk = rj->coeffs + k;
                if (fmpz_is_zero(rk)) continue;
                fmpz_mul(v, rk, c);
                e[N->tau] = (ulong)k;
                fmpz_mpoly_push_term_fmpz_ui(out, v, e, zc);
            }
        }
        fmpz_mpoly_sort_terms(out, zc);
        fmpz_mpoly_combine_like_terms(out, zc);
        fmpz_mpoly_swap(zp, out, zc);
        fmpz_mpoly_clear(out, zc);
        fmpz_clear(v);
        fmpq_mpoly_reduce(p, N->ctx);
        return;
    }
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
    /* The radical function field (nfrat_radical_normal_form): a radical of a
     * non-constant radicand is a generator over K(x); its radicand supplies
     * variables and atoms. */
    if (N->rad && head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        head_is(e->data.function.args[1], SYM_Rational) && nr_has_symbol(e->data.function.args[0]))
        return nr_scan(e->data.function.args[0], N, acap);
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
    if (N->tpowz) {
        for (slong j = 0; j < N->ntpow; j++) fmpz_poly_clear(N->tpowz + j);
        free(N->tpowz);
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

static int nr_setup_field(Nr* N);

/* allow_const: an expression with no variable is accepted (the canonical
 * form of a constant; Together / Cancel leave numbers to the evaluator). */
static int nr_setup_ex(Nr* N, const Expr* e, const Expr* ext, int allow_const) {
    memset(N, 0, sizeof *N);
    size_t acap = 0;
    if (!nr_scan(e, N, &acap) || !nr_scan_ext(ext, N, &acap) || N->na == 0 ||
        (N->ngens == 0 && !allow_const))
        return 0;
    return nr_setup_field(N);
}

/* The field and the polynomial context, once the atoms and the variables of
 * the input have been scanned into N. */
static int nr_setup_field(Nr* N) {
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
            N->tpowz = malloc((size_t)N->ntpow * sizeof(fmpz_poly_struct));
            for (slong j = 0; j < N->ntpow; j++) fmpz_poly_init(N->tpowz + j);
            for (slong j = 0; j < N->ntpow && N->tpowz; j++) {
                if (!fmpz_is_one(fmpq_poly_denref(N->tpow + j))) {
                    for (slong t = 0; t < N->ntpow; t++) fmpz_poly_clear(N->tpowz + t);
                    free(N->tpowz); N->tpowz = NULL;
                    break;
                }
                fmpq_poly_get_numerator(N->tpowz + j, N->tpow + j);
            }
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

static int nr_setup(Nr* N, const Expr* e, const Expr* ext) { return nr_setup_ex(N, e, ext, 0); }

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

/* ------------------------------------------------------------------ */
/*  The radical function field: Simplify's normal form                 */
/* ------------------------------------------------------------------ */
/*
 * An expression rational in variables x1..xm, in constant algebraic atoms
 * (the field K above) and in radicals y = B^(p/q) of NON-constant radicands
 * B in K(x) is an element of
 *
 *      L = K(x)[h_1, ..., h_r] / ( h_i^Q_i - R_i ),
 *
 * one generator per radicand class:
 *
 *   - a class is a radicand up to a POSITIVE RATIONAL factor.  B^(p/q) with
 *     B = c B_i, c > 0 rational, is c^(p/q) B_i^(p/q): valid for principal
 *     values because arg(c B_i) = arg(B_i).  No other identification is made
 *     (Sqrt[a] Sqrt[b] -> Sqrt[a b] and (-B)^(1/2) -> I B^(1/2) fail on part
 *     of the plane), so dependent radicands simply keep separate generators;
 *   - Q_i is the lcm of the indices of the class, g_i = B_i^(1/Q_i), and
 *     B_i^(p/q) = g_i^(p Q_i/q).  An integer power of ONE principal radical is
 *     always exact: B^(k/Q) = Exp[(k/Q) Log B] = (B^(1/Q))^k;
 *   - B_i = N/D over K[x] (D monic) is rescaled: h_i = D g_i is integral,
 *     h_i^Q_i = R_i = N D^(Q_i - 1).
 *
 * Arithmetic is exact in L.  An element is A / P_1^e_1 ... P_s^e_s with A a
 * polynomial over K in (x, h) reduced modulo the relations (deg_h_i A < Q_i)
 * and the P_j generator-free monic pieces of the per-call piece table, so a
 * sum is a merge of piece lists.  An inverse uses the norm: for each
 * generator (outermost first) the multiplication-by-A map on the free module
 * with basis 1, h, .., h^(Q-1) has characteristic polynomial
 * t^Q - e_1 t^(Q-1) + .. + (-1)^Q e_Q, whose coefficients follow from the
 * traces Tr(A^k) = Q * [h^0] A^k by Newton's identities (integer divisions
 * only), and Cayley-Hamilton gives A * adj = e_Q = Norm(A) with
 * adj = (-1)^(Q-1) (A^(Q-1) - e_1 A^(Q-2) + .. + (-1)^(Q-1) e_(Q-1)).  The
 * identity holds in any commutative ring, so the rewrite is exact whatever the
 * relations; a zero norm (a zero divisor: (x^2)^(1/2) - x) declines.
 *
 * At the end the content of A is cancelled against the pieces (certified gcds
 * over K), which makes the form canonical whenever L is a field.  Kummer
 * theory decides that (rn_kummer_field): h^Q - R is irreducible over K(x)
 * adjoined the other generators when R has a prime factor of multiplicity
 * prime to Q that divides no other radicand.  Then, and only then, a nonzero
 * normal form proves a nonzero expression; a zero normal form is a proof of
 * zero always.
 *
 * The normal form is read back as Sum_m c_m(x) y^m with c_m in K(x), each
 * c_m factored and written in the caller's own constant radicals, and a
 * radicand factor of c_m absorbed into the radical when that shortens it
 * (c y^m = (c B) y^(m - Q): (8 + x^3)^(2/3)/(8 + x^3) is (8 + x^3)^(-1/3)).
 */

#define RN_MAX_GEN    6
#define RN_MAX_DIM    64
#define RN_MAX_SPELL  32

typedef struct {
    Expr*   base;          /* the class's radicand as first written (owned)    */
    Expr*   key;           /* Expand[base] (owned)                              */
    int64_t Q;             /* g = base^(1/Q)                                    */
    const char* name;      /* interned generator symbol                         */
    slong   v;             /* its variable index                                */
    NrRf    B;             /* the radicand over K(x)                            */
    fmpq_mpoly_t R;        /* h^Q = R                                           */
    fmpq_mpoly_struct* Rp; /* Rp[k - 1] = R^k, cached                           */
    slong   nRp;
    int     live;          /* B, R initialised                                  */
    int     elim;          /* bare-symbol radicand s: s = t^Q with t a
                            * transcendental variable (exact on every branch,
                            * no relation, no dimension)                       */
} RnGen;

typedef struct {
    const Expr* spell;     /* a radicand as the input writes it (borrowed)      */
    int   gi;              /* its class                                         */
    Expr* ratio;           /* NULL, or the positive rational r: spell = r base  */
} RnSpell;

typedef struct {
    RnGen   g[RN_MAX_GEN];
    int     ng;
    RnSpell s[RN_MAX_SPELL];
    int     ns;
    Expr**  extra;         /* constants r^(p/q) the field must contain          */
    size_t  nextra, cextra;
} RnTab;

typedef struct {
    fmpq_mpoly_t A;        /* reduced numerator in (x, h, tau)                  */
    NrPE* de; size_t nde, cde;   /* generator-free monic pieces                 */
} RnE;

static void rn_tab_clear(Nr* N, RnTab* T) {
    for (int i = 0; i < T->ng; i++) {
        RnGen* g = &T->g[i];
        if (g->base) expr_free(g->base);
        if (g->key) expr_free(g->key);
        if (g->live && N && N->ctx) {
            rf_clear(N, &g->B);
            fmpq_mpoly_clear(g->R, N->ctx);
            for (slong k = 0; k < g->nRp; k++) fmpq_mpoly_clear(g->Rp + k, N->ctx);
        }
        free(g->Rp);
    }
    for (int i = 0; i < T->ns; i++) if (T->s[i].ratio) expr_free(T->s[i].ratio);
    for (size_t i = 0; i < T->nextra; i++) expr_free(T->extra[i]);
    free(T->extra);
    memset(T, 0, sizeof *T);
}

static Expr* rn_call1(const char* head, Expr* arg) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ arg }, 1));
}

static bool rn_is_positive_rational(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) > 0;
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        if (p->type == EXPR_INTEGER) return p->data.integer > 0;
        if (p->type == EXPR_BIGINT) return mpz_sgn(p->data.bigint) > 0;
    }
    return false;
}

static int rn_push_extra(RnTab* T, Expr* c) {
    if (T->nextra == T->cextra) {
        T->cextra = T->cextra ? 2 * T->cextra : 4;
        T->extra = realloc(T->extra, T->cextra * sizeof(Expr*));
    }
    T->extra[T->nextra++] = c;
    return 1;
}

/* r^(p/q) for the positive rational r of a merged spelling (evaluated). */
static Expr* rn_ratio_power(const Expr* r, int64_t p, int64_t q) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)r), make_rational(p, q) }, 2));
}

/* Record the radical b^(p/q) (b non-constant, non-nested). */
static int rn_register(RnTab* T, const Expr* b, int64_t p, int64_t q) {
    for (int i = 0; i < T->ns; i++)
        if (expr_eq((Expr*)T->s[i].spell, (Expr*)b)) {
            RnGen* g = &T->g[T->s[i].gi];
            g->Q = lcm(g->Q, q);
            if (T->s[i].ratio) rn_push_extra(T, rn_ratio_power(T->s[i].ratio, p, q));
            return 1;
        }
    if (T->ns == RN_MAX_SPELL) return 0;
    Expr* key = rn_call1(SYM_Expand, expr_copy((Expr*)b));
    if (!key) return 0;
    for (int i = 0; i < T->ng; i++)
        if (expr_eq(T->g[i].key, key)) {
            expr_free(key);
            T->s[T->ns++] = (RnSpell){ b, i, NULL };
            T->g[i].Q = lcm(T->g[i].Q, q);
            return 1;
        }
    /* the same radicand up to a positive rational factor shares the class */
    for (int i = 0; i < T->ng; i++) {
        Expr* r = rn_call1(SYM_Together, expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_copy(key),
                                 expr_new_function(expr_new_symbol(SYM_Power),
                                     (Expr*[]){ expr_copy(T->g[i].key), expr_new_integer(-1) }, 2) }, 2));
        if (r && rn_is_positive_rational(r)) {
            expr_free(key);
            T->s[T->ns++] = (RnSpell){ b, i, r };
            T->g[i].Q = lcm(T->g[i].Q, q);
            rn_push_extra(T, rn_ratio_power(r, p, q));
            return 1;
        }
        if (r) expr_free(r);
    }
    if (T->ng == RN_MAX_GEN) { expr_free(key); return 0; }
    RnGen* g = &T->g[T->ng];
    memset(g, 0, sizeof *g);
    g->base = expr_copy((Expr*)b);
    g->key = key;
    g->Q = q;
    g->elim = b->type == EXPR_SYMBOL && !nk_reserved_symbol(b->data.symbol.name);
    char name[48];
    snprintf(name, sizeof name, "~MthRadFieldGen%d", T->ng + 1);   /* sorts after the variables */
    g->name = intern_symbol(name);
    T->s[T->ns++] = (RnSpell){ b, T->ng, NULL };
    T->ng++;
    return 1;
}

/* Walk the rational structure: 1 when every leaf is a number, a symbol, a
 * constant algebraic atom or a radical of a radicand that is itself such a
 * rational function (no nested symbolic radical). */
static int rn_collect(const Expr* e, RnTab* T, int inside) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_SYMBOL: return 1;
        case EXPR_FUNCTION: break;
        default: return 0;
    }
    if (head_is(e, SYM_Rational)) return 1;
    if (head_is(e, SYM_Complex)) return nr_is_alg_leaf(e);
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!rn_collect(e->data.function.args[i], T, inside)) return 0;
        return 1;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER) return rn_collect(b, T, inside);
        if (!head_is(x, SYM_Rational)) return 0;
        if (!nr_has_symbol(b)) return nr_is_alg_leaf(e);
        if (inside) return 0;                          /* nested: algfield's business */
        const Expr* p = x->data.function.args[0];
        const Expr* q = x->data.function.args[1];
        if (p->type != EXPR_INTEGER || q->type != EXPR_INTEGER) return 0;
        if (!rn_collect(b, T, 1)) return 0;
        return rn_register(T, b, p->data.integer, q->data.integer);
    }
    if (head_is(e, SYM_Root)) return !nr_has_symbol(e);
    return 0;
}

/* ---- arithmetic in L ---- */

static void rn_init(Nr* N, RnE* a) { memset(a, 0, sizeof *a); fmpq_mpoly_init(a->A, N->ctx); }
static void rn_clear(Nr* N, RnE* a) { fmpq_mpoly_clear(a->A, N->ctx); free(a->de); }
static bool rn_is_zero(Nr* N, const RnE* a) { return fmpq_mpoly_is_zero(a->A, N->ctx); }

static void rn_set_const(Nr* N, RnE* a, const fmpq_mpoly_t c) {
    fmpq_mpoly_set(a->A, c, N->ctx);
    a->nde = 0;
}

/* R_i^k (k >= 1), cached. */
static const fmpq_mpoly_struct* rn_rpow(Nr* N, RnGen* g, slong k) {
    while (g->nRp < k) {
        g->Rp = realloc(g->Rp, (size_t)(g->nRp + 1) * sizeof(fmpq_mpoly_struct));
        fmpq_mpoly_init(g->Rp + g->nRp, N->ctx);
        if (g->nRp == 0) fmpq_mpoly_set(g->Rp, g->R, N->ctx);
        else nr_mul(N, g->Rp + g->nRp, g->Rp + g->nRp - 1, g->R);
        g->nRp++;
    }
    return g->Rp + (k - 1);
}

/* p := p reduced modulo the relations h_i^Q_i = R_i and M(tau).  The R_i are
 * generator-free, so one pass over the terms suffices. */
static void rn_reduce(Nr* N, fmpq_mpoly_t p) {
    RnTab* T = N->rad;
    bool need = false;
    for (int i = 0; i < T->ng && !need; i++)
        need = !T->g[i].elim && fmpq_mpoly_degree_si(p, T->g[i].v, N->ctx) >= T->g[i].Q;
    if (!need) { nr_reduce(N, p); return; }
    fmpq_mpoly_t out, acc, mono;
    fmpq_mpoly_init(out, N->ctx); fmpq_mpoly_init(acc, N->ctx); fmpq_mpoly_init(mono, N->ctx);
    ulong e[NR_MAX_VARS + 1];
    fmpq_t c; fmpq_init(c);
    slong len = fmpq_mpoly_length(p, N->ctx);
    for (slong t = 0; t < len; t++) {
        fmpq_mpoly_get_term_exp_ui(e, p, t, N->ctx);
        fmpq_mpoly_get_term_coeff_fmpq(c, p, t, N->ctx);
        slong kk[RN_MAX_GEN];
        bool any = false;
        for (int i = 0; i < T->ng; i++) {
            slong v = T->g[i].v;
            if (T->g[i].elim) { kk[i] = 0; continue; }
            kk[i] = (slong)e[v] / (slong)T->g[i].Q;
            e[v] = (ulong)((slong)e[v] % (slong)T->g[i].Q);
            if (kk[i]) any = true;
        }
        if (!any) { fmpq_mpoly_push_term_fmpq_ui(out, c, e, N->ctx); continue; }
        fmpq_mpoly_zero(mono, N->ctx);
        fmpq_mpoly_set_coeff_fmpq_ui(mono, c, e, N->ctx);
        for (int i = 0; i < T->ng; i++)
            if (kk[i]) fmpq_mpoly_mul(mono, mono, rn_rpow(N, &T->g[i], kk[i]), N->ctx);
        fmpq_mpoly_add(acc, acc, mono, N->ctx);
    }
    fmpq_mpoly_sort_terms(out, N->ctx);
    fmpq_mpoly_combine_like_terms(out, N->ctx);
    fmpq_mpoly_add(out, out, acc, N->ctx);
    nr_reduce(N, out);
    fmpq_mpoly_swap(p, out, N->ctx);
    fmpq_mpoly_clear(out, N->ctx); fmpq_mpoly_clear(acc, N->ctx); fmpq_mpoly_clear(mono, N->ctx);
    fmpq_clear(c);
}

static int rn_mulp(Nr* N, fmpq_mpoly_t dst, const fmpq_mpoly_t a, const fmpq_mpoly_t b) {
    fmpq_mpoly_mul(dst, a, b, N->ctx);
    rn_reduce(N, dst);
    return fmpq_mpoly_length(dst, N->ctx) <= NR_TERM_BUDGET;
}

/* dst := dst * b */
static int rn_mul(Nr* N, RnE* dst, const RnE* b) {
    if (!rn_mulp(N, dst->A, dst->A, b->A)) return 0;
    for (size_t k = 0; k < b->nde; k++) pe_add(&dst->de, &dst->nde, &dst->cde, b->de[k].i, b->de[k].e);
    return 1;
}

/* The [h_v^0] part of p. */
static void rn_coeff0(Nr* N, fmpq_mpoly_t out, const fmpq_mpoly_t p, slong v) {
    fmpq_mpoly_zero(out, N->ctx);
    ulong e[NR_MAX_VARS + 1];
    fmpq_t c; fmpq_init(c);
    for (slong t = 0; t < fmpq_mpoly_length(p, N->ctx); t++) {
        fmpq_mpoly_get_term_exp_ui(e, p, t, N->ctx);
        if (e[v]) continue;
        fmpq_mpoly_get_term_coeff_fmpq(c, p, t, N->ctx);
        fmpq_mpoly_push_term_fmpq_ui(out, c, e, N->ctx);
    }
    fmpq_mpoly_sort_terms(out, N->ctx);
    fmpq_mpoly_combine_like_terms(out, N->ctx);
    fmpq_clear(c);
}

/* A * adj = nrm with nrm free of every generator (see the header comment).
 * 0 on a zero norm (A a zero divisor) or past the budget. */
static int rn_norm(Nr* N, const fmpq_mpoly_t A, fmpq_mpoly_t adj, fmpq_mpoly_t nrm) {
    RnTab* T = N->rad;
    fmpq_mpoly_t cur; fmpq_mpoly_init(cur, N->ctx);
    fmpq_mpoly_set(cur, A, N->ctx);
    fmpq_mpoly_one(adj, N->ctx);
    int ok = 1;
    for (int i = T->ng - 1; i >= 0 && ok; i--) {
        slong v = T->g[i].v, Q = (slong)T->g[i].Q;
        if (T->g[i].elim || fmpq_mpoly_degree_si(cur, v, N->ctx) <= 0) continue;
        /* P[k] = cur^k, k = 0 .. Q; E[k] = e_k */
        fmpq_mpoly_struct* P = malloc((size_t)(Q + 1) * sizeof(fmpq_mpoly_struct));
        fmpq_mpoly_struct* E = malloc((size_t)(Q + 1) * sizeof(fmpq_mpoly_struct));
        fmpq_mpoly_struct* pt = malloc((size_t)(Q + 1) * sizeof(fmpq_mpoly_struct));
        for (slong k = 0; k <= Q; k++) {
            fmpq_mpoly_init(P + k, N->ctx); fmpq_mpoly_init(E + k, N->ctx); fmpq_mpoly_init(pt + k, N->ctx);
        }
        fmpq_mpoly_one(P + 0, N->ctx);
        fmpq_mpoly_set(P + 1, cur, N->ctx);
        for (slong k = 2; k <= Q && ok; k++) ok = rn_mulp(N, P + k, P + k - 1, cur);
        fmpq_mpoly_t t; fmpq_mpoly_init(t, N->ctx);
        fmpq_t qf; fmpq_init(qf);
        /* traces p_k = Q [h^0] cur^k */
        for (slong k = 1; k <= Q && ok; k++) {
            rn_coeff0(N, pt + k, P + k, v);
            fmpq_mpoly_scalar_mul_si(pt + k, pt + k, Q, N->ctx);
        }
        /* Newton: k e_k = sum_{j=1}^k (-1)^(j-1) e_(k-j) p_j */
        fmpq_mpoly_one(E + 0, N->ctx);
        for (slong k = 1; k <= Q && ok; k++) {
            fmpq_mpoly_zero(E + k, N->ctx);
            for (slong j = 1; j <= k && ok; j++) {
                ok = rn_mulp(N, t, E + k - j, pt + j);
                if (j % 2 == 1) fmpq_mpoly_add(E + k, E + k, t, N->ctx);
                else fmpq_mpoly_sub(E + k, E + k, t, N->ctx);
            }
            fmpq_set_si(qf, 1, (ulong)k);
            fmpq_mpoly_scalar_mul_fmpq(E + k, E + k, qf, N->ctx);
        }
        if (ok && fmpq_mpoly_is_zero(E + Q, N->ctx)) ok = 0;
        if (ok) {
            /* adj_i = (-1)^(Q-1) sum_{k=0}^{Q-1} (-1)^k e_k cur^(Q-1-k) */
            fmpq_mpoly_t ai; fmpq_mpoly_init(ai, N->ctx);
            for (slong k = 0; k < Q && ok; k++) {
                ok = rn_mulp(N, t, E + k, P + (Q - 1 - k));
                if (k % 2 == 0) fmpq_mpoly_add(ai, ai, t, N->ctx);
                else fmpq_mpoly_sub(ai, ai, t, N->ctx);
            }
            if ((Q - 1) % 2 == 1) fmpq_mpoly_neg(ai, ai, N->ctx);
            if (ok) ok = rn_mulp(N, adj, adj, ai);
            fmpq_mpoly_set(cur, E + Q, N->ctx);
            fmpq_mpoly_clear(ai, N->ctx);
        }
        fmpq_clear(qf);
        fmpq_mpoly_clear(t, N->ctx);
        for (slong k = 0; k <= Q; k++) {
            fmpq_mpoly_clear(P + k, N->ctx); fmpq_mpoly_clear(E + k, N->ctx); fmpq_mpoly_clear(pt + k, N->ctx);
        }
        free(P); free(E); free(pt);
    }
    if (ok) fmpq_mpoly_swap(nrm, cur, N->ctx);
    fmpq_mpoly_clear(cur, N->ctx);
    return ok && !fmpq_mpoly_is_zero(nrm, N->ctx);
}

/* a := 1 / a.  0 on a zero divisor. */
static int rn_inv(Nr* N, RnE* a) {
    if (rn_is_zero(N, a)) return 0;
    fmpq_mpoly_t adj, nrm, lc;
    fmpq_mpoly_init(adj, N->ctx); fmpq_mpoly_init(nrm, N->ctx); fmpq_mpoly_init(lc, N->ctx);
    int ok = rn_norm(N, a->A, adj, nrm);
    for (size_t j = 0; j < a->nde && ok; j++) ok = nr_mul_piece_pow(N, adj, a->de[j].i, a->de[j].e);
    if (ok) {
        a->nde = 0;
        if (nr_is_const(N, nrm)) {
            ok = nr_const_inv(N, lc, nrm);
            if (ok) nr_mul(N, adj, adj, lc);
        } else {
            ok = nr_make_monic(N, nrm, lc);
            if (ok) ok = nr_const_inv(N, lc, lc);
            if (ok) {
                nr_mul(N, adj, adj, lc);
                pe_add(&a->de, &a->nde, &a->cde, nr_piece(N, nrm, NULL), 1);
            }
        }
        fmpq_mpoly_swap(a->A, adj, N->ctx);
    }
    fmpq_mpoly_clear(adj, N->ctx); fmpq_mpoly_clear(nrm, N->ctx); fmpq_mpoly_clear(lc, N->ctx);
    return ok;
}

static int rn_pow(Nr* N, RnE* a, slong k) {
    if (k < 0) { if (!rn_inv(N, a)) return 0; k = -k; }
    /* (A / prod P^e)^k = A^k / prod P^(k e): the pieces only scale */
    for (size_t j = 0; j < a->nde; j++) a->de[j].e *= k;
    if (k == 0) a->nde = 0;
    fmpq_mpoly_t acc, b;
    fmpq_mpoly_init(acc, N->ctx); fmpq_mpoly_init(b, N->ctx);
    fmpq_mpoly_one(acc, N->ctx);
    fmpq_mpoly_set(b, a->A, N->ctx);
    int ok = 1;
    while (k && ok) {
        if (k & 1) ok = rn_mulp(N, acc, acc, b);
        k >>= 1;
        if (k && ok) ok = rn_mulp(N, b, b, b);
    }
    fmpq_mpoly_swap(a->A, acc, N->ctx);
    fmpq_mpoly_clear(acc, N->ctx); fmpq_mpoly_clear(b, N->ctx);
    return ok;
}

/* dst := sum of the k terms. */
static int rn_sum(Nr* N, RnE* dst, RnE* T, size_t k) {
    NrPE* L = NULL; size_t nl = 0, cl = 0;
    for (size_t t = 0; t < k; t++) {
        if (rn_is_zero(N, T + t)) continue;
        for (size_t j = 0; j < T[t].nde; j++) {
            slong have = pe_get(L, nl, T[t].de[j].i);
            if (T[t].de[j].e > have) pe_add(&L, &nl, &cl, T[t].de[j].i, T[t].de[j].e - have);
        }
    }
    fmpq_mpoly_t S, acc;
    fmpq_mpoly_init(S, N->ctx); fmpq_mpoly_init(acc, N->ctx);
    int ok = 1;
    for (size_t t = 0; t < k && ok; t++) {
        if (rn_is_zero(N, T + t)) continue;
        fmpq_mpoly_set(acc, T[t].A, N->ctx);
        for (size_t j = 0; j < nl && ok; j++)
            ok = nr_mul_piece_pow(N, acc, L[j].i, L[j].e - pe_get(T[t].de, T[t].nde, L[j].i));
        fmpq_mpoly_add(S, S, acc, N->ctx);
        ok = ok && fmpq_mpoly_length(S, N->ctx) <= NR_TERM_BUDGET;
    }
    if (ok) {
        fmpq_mpoly_swap(dst->A, S, N->ctx);
        free(dst->de);
        if (fmpq_mpoly_is_zero(dst->A, N->ctx)) { free(L); L = NULL; nl = cl = 0; }
        dst->de = L; dst->nde = nl; dst->cde = cl;
        L = NULL;
    }
    free(L);
    fmpq_mpoly_clear(S, N->ctx); fmpq_mpoly_clear(acc, N->ctx);
    return ok;
}

/* A constant (symbol-free) subexpression, through the field's own evaluator. */
static int rn_eval_const(Nr* N, const Expr* e, RnE* a) {
    NrRf r; rf_init(N, &r);
    int ok = nr_eval(N, e, &r) && r.nnu == 0 && r.nde == 0;
    if (ok) {
        fmpq_mpoly_t ui; fmpq_mpoly_init(ui, N->ctx);
        ok = nr_const_inv(N, ui, r.u);
        if (ok) { nr_mul(N, ui, ui, r.c); rn_set_const(N, a, ui); }
        fmpq_mpoly_clear(ui, N->ctx);
    }
    rf_clear(N, &r);
    return ok;
}

static int rn_spell_index(const RnTab* T, const Expr* b) {
    for (int i = 0; i < T->ns; i++) if (expr_eq((Expr*)T->s[i].spell, (Expr*)b)) return i;
    return -1;
}

static int rn_eval(Nr* N, const Expr* e, RnE* a) {
    RnTab* T = N->rad;
    if (!nr_has_symbol(e)) return rn_eval_const(N, e, a);
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < T->ng; i++)
            if (T->g[i].elim && strcmp(T->g[i].base->data.symbol.name, e->data.symbol.name) == 0) {
                fmpq_mpoly_gen(a->A, T->g[i].v, N->ctx);       /* s = t^Q */
                fmpq_mpoly_pow_ui(a->A, a->A, (ulong)T->g[i].Q, N->ctx);
                a->nde = 0;
                return 1;
            }
        long v = nr_var_index(N, e->data.symbol.name);
        if (v < 0) return 0;
        fmpq_mpoly_gen(a->A, v, N->ctx);
        a->nde = 0;
        return 1;
    }
    if (e->type != EXPR_FUNCTION) return 0;
    bool plus = head_is(e, SYM_Plus);
    if (plus || head_is(e, SYM_Times)) {
        size_t k = e->data.function.arg_count;
        if (k == 0) return 0;
        if (!plus) {
            int ok = rn_eval(N, e->data.function.args[0], a);
            RnE t; rn_init(N, &t);
            for (size_t i = 1; i < k && ok; i++) {
                ok = rn_eval(N, e->data.function.args[i], &t);
                if (ok) ok = rn_mul(N, a, &t);
            }
            rn_clear(N, &t);
            return ok;
        }
        RnE* Ts = malloc(k * sizeof(RnE));
        size_t made = 0;
        int ok = Ts != NULL;
        for (size_t i = 0; i < k && ok; i++) {
            rn_init(N, Ts + i); made++;
            ok = rn_eval(N, e->data.function.args[i], Ts + i);
        }
        if (ok) ok = rn_sum(N, a, Ts, made);
        for (size_t i = 0; i < made; i++) rn_clear(N, Ts + i);
        free(Ts);
        return ok;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER) {
            if (!rn_eval(N, b, a)) return 0;
            return rn_pow(N, a, (slong)x->data.integer);
        }
        if (!head_is(x, SYM_Rational)) return 0;
        int si = rn_spell_index(T, b);
        if (si < 0) return 0;
        RnGen* g = &T->g[T->s[si].gi];
        int64_t p = x->data.function.args[0]->data.integer;
        int64_t q = x->data.function.args[1]->data.integer;
        if (g->Q % q) return 0;
        /* g^k = (h / D)^k, D the radicand's denominator pieces */
        fmpq_mpoly_gen(a->A, g->v, N->ctx);
        a->nde = 0;
        for (size_t j = 0; j < g->B.nde; j++) pe_add(&a->de, &a->nde, &a->cde, g->B.de[j].i, g->B.de[j].e);
        if (!rn_pow(N, a, (slong)(p * (g->Q / q)))) return 0;
        if (T->s[si].ratio) {
            Expr* c = rn_ratio_power(T->s[si].ratio, p, q);
            RnE ce; rn_init(N, &ce);
            int ok = rn_eval_const(N, c, &ce) && rn_mul(N, a, &ce);
            rn_clear(N, &ce);
            expr_free(c);
            return ok;
        }
        return 1;
    }
    return 0;
}

/* Cancel the content of A against the pieces: canonical when L is a field. */
static int rn_cancel(Nr* N, RnE* a) {
    if (rn_is_zero(N, a)) { a->nde = 0; return 1; }
    fmpq_mpoly_t q, g; fmpq_mpoly_init(q, N->ctx); fmpq_mpoly_init(g, N->ctx);
    int ok = 1;
    size_t j = 0;
    while (j < a->nde && ok) {
        size_t pi = a->de[j].i;
        int d = nr_kron_divides(N, q, a->A, N->pc + pi);
        if (d < 0) d = fb_field_mpoly_divides(q, a->A, N->pc + pi, N->M, N->ctx);
        if (d == 1) {
            fmpq_mpoly_swap(a->A, q, N->ctx);
            size_t before = a->nde;
            pe_add(&a->de, &a->nde, &a->cde, pi, -1);
            if (a->nde < before) j = 0;              /* entry j replaced: rescan */
            continue;
        }
        if (fb_field_mpoly_gcd(g, a->A, N->pc + pi, N->M, N->ctx, N->mp, N->ngens, N->tau, N->n) &&
            !nr_is_const(N, g)) {
            size_t gi, qi;
            if (!nr_split(N, pi, g, &gi, &qi)) { ok = 0; break; }
            pe_replace(&a->de, &a->nde, &a->cde, pi, gi, qi);
            j = 0;
            continue;
        }
        j++;
    }
    fmpq_mpoly_clear(q, N->ctx); fmpq_mpoly_clear(g, N->ctx);
    return ok;
}

/* ---- setup ---- */

/* e with every eliminated radicand symbol s written t^Q. */
static Expr* rn_elim_subst(const Expr* e, const RnTab* T) {
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < T->ng; i++)
            if (T->g[i].elim && strcmp(T->g[i].base->data.symbol.name, e->data.symbol.name) == 0)
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol(T->g[i].name), expr_new_integer(T->g[i].Q) }, 2);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rn_elim_subst(e->data.function.args[i], T);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

/* e with every eliminated generator t written s^(1/Q) (the read-back). */
static Expr* rn_elim_back(const Expr* e, const RnTab* T) {
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < T->ng; i++)
            if (T->g[i].elim && strcmp(T->g[i].name, e->data.symbol.name) == 0)
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(T->g[i].base), make_rational(1, T->g[i].Q) }, 2);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rn_elim_back(e->data.function.args[i], T);
    Expr* r = expr_new_function(rn_elim_back(e->data.function.head, T), a, n);
    free(a);
    return r;
}

static int rn_setup(Nr* N, const Expr* e, RnTab* T) {
    memset(N, 0, sizeof *N);
    N->rad = T;
    size_t acap = 0;
    if (!nr_scan(e, N, &acap)) return 0;
    for (size_t i = 0; i < T->nextra; i++)
        if (!nr_scan(T->extra[i], N, &acap)) return 0;
    for (int i = 0; i < T->ng; i++) {
        if (N->ngens == NR_MAX_VARS) return 0;
        N->var[N->ngens++] = T->g[i].name;
    }
    if (N->na == 0) {
        /* K = Q: Q(I) holds the same rational functions (a field operation on
         * rational data stays rational), and the engine needs a nontrivial
         * primitive element. */
        Expr* iu = make_complex(expr_new_integer(0), expr_new_integer(1));
        int ok = fb_collect_atoms(iu, &N->atoms, &N->na, &acap);
        expr_free(iu);
        if (!ok) return 0;
    }
    if (!nr_setup_field(N)) return 0;
    for (int i = 0; i < T->ng; i++) {
        RnGen* g = &T->g[i];
        g->v = nr_var_index(N, g->name);
        if (g->v < 0) return 0;
        rf_init(N, &g->B);
        fmpq_mpoly_init(g->R, N->ctx);
        g->live = 1;
        Expr* bt = rn_elim_subst(g->base, T);
        int okb = nr_eval(N, bt, &g->B) && !rf_is_zero(N, &g->B);
        expr_free(bt);
        if (!okb) return 0;
        if (g->elim) continue;                 /* t is free: no relation */
        /* R = (c/u) prod nu^e * D^(Q-1), D = prod de^e */
        fmpq_mpoly_t num, den, ui;
        fmpq_mpoly_init(num, N->ctx); fmpq_mpoly_init(den, N->ctx); fmpq_mpoly_init(ui, N->ctx);
        int ok = nr_const_inv(N, ui, g->B.u);
        if (ok) { nr_mul(N, num, g->B.c, ui); fmpq_mpoly_one(den, N->ctx); }
        for (size_t j = 0; j < g->B.nnu && ok; j++) ok = nr_mul_piece_pow(N, num, g->B.nu[j].i, g->B.nu[j].e);
        for (size_t j = 0; j < g->B.nde && ok; j++) ok = nr_mul_piece_pow(N, den, g->B.de[j].i, g->B.de[j].e);
        for (int64_t k = 1; k < g->Q && ok; k++) ok = rn_mulp(N, num, num, den);
        if (ok) fmpq_mpoly_swap(g->R, num, N->ctx);
        fmpq_mpoly_clear(num, N->ctx); fmpq_mpoly_clear(den, N->ctx); fmpq_mpoly_clear(ui, N->ctx);
        if (!ok) return 0;
    }
    return 1;
}

/* Kummer: is L a field of degree prod Q_i over K(x)?  Each radicand needs a
 * piece of multiplicity prime to its index that is squarefree and coprime to
 * every other piece of every radicand (an Eisenstein prime for its generator,
 * unramified in the others). */
static bool rn_piece_squarefree(Nr* N, size_t a, fmpq_mpoly_t g, fmpq_mpoly_t d) {
    for (slong v = 0; v < N->ngens; v++) {
        if (fmpq_mpoly_degree_si(N->pc + a, v, N->ctx) < 1) continue;
        fmpq_mpoly_derivative(d, N->pc + a, v, N->ctx);
        if (fmpq_mpoly_is_zero(d, N->ctx)) continue;
        if (!fb_field_mpoly_gcd(g, N->pc + a, d, N->M, N->ctx, N->mp, N->ngens, N->tau, N->n)) return false;
        return nr_is_const(N, g);
    }
    return false;
}

static bool rn_kummer_field(Nr* N, RnTab* T) {
    fmpq_mpoly_t g, d; fmpq_mpoly_init(g, N->ctx); fmpq_mpoly_init(d, N->ctx);
    bool all = true;
    for (int i = 0; i < T->ng && all; i++) {
        if (T->g[i].elim) continue;            /* a free variable: no relation */
        bool found = false;
        for (int side = 0; side < 2 && !found; side++) {
            const NrPE* L = side ? T->g[i].B.de : T->g[i].B.nu;
            size_t n = side ? T->g[i].B.nde : T->g[i].B.nnu;
            for (size_t j = 0; j < n && !found; j++) {
                if (gcd(L[j].e < 0 ? -L[j].e : L[j].e, T->g[i].Q) != 1) continue;
                size_t a = L[j].i;
                if (!rn_piece_squarefree(N, a, g, d)) continue;
                bool coprime = true;
                for (int k = 0; k < T->ng && coprime; k++)
                    for (int s2 = 0; s2 < 2 && coprime; s2++) {
                        const NrPE* L2 = s2 ? T->g[k].B.de : T->g[k].B.nu;
                        size_t n2 = s2 ? T->g[k].B.nde : T->g[k].B.nnu;
                        for (size_t j2 = 0; j2 < n2 && coprime; j2++) {
                            size_t b = L2[j2].i;
                            if (k == i && s2 == side && j2 == j) continue;
                            if (b == a) { coprime = false; break; }
                            if (!fb_field_mpoly_gcd(g, N->pc + a, N->pc + b, N->M, N->ctx, N->mp,
                                                    N->ngens, N->tau, N->n) || !nr_is_const(N, g))
                                coprime = false;
                        }
                    }
                found = coprime;
            }
        }
        all = found;
    }
    fmpq_mpoly_clear(g, N->ctx); fmpq_mpoly_clear(d, N->ctx);
    return all;
}

/* ---- read back ---- */

static void rf_copy(Nr* N, NrRf* dst, const NrRf* src) {
    fmpq_mpoly_set(dst->c, src->c, N->ctx);
    fmpq_mpoly_set(dst->u, src->u, N->ctx);
    dst->nnu = dst->nde = 0;
    for (size_t k = 0; k < src->nnu; k++) pe_add(&dst->nu, &dst->nnu, &dst->cnu, src->nu[k].i, src->nu[k].e);
    for (size_t k = 0; k < src->nde; k++) pe_add(&dst->de, &dst->nde, &dst->cde, src->de[k].i, src->de[k].e);
}

/* Total x-degree of the pieces (numerator and denominator): the size the
 * absorption of a radicand factor is judged by. */
static slong rf_degree(Nr* N, const NrRf* r) {
    slong d = 0;
    for (size_t k = 0; k < r->nnu; k++) d += r->nu[k].e * nr_xdeg(N, N->pc + r->nu[k].i, WORD_MAX);
    for (size_t k = 0; k < r->nde; k++) d += r->de[k].e * nr_xdeg(N, N->pc + r->de[k].i, WORD_MAX);
    return d;
}

static Expr* rn_present(Nr* N, RnTab* T, const RnE* a) {
    if (rn_is_zero(N, a)) return expr_new_integer(0);
    /* the h-monomials of A and their coefficient polynomials */
    slong len = fmpq_mpoly_length(a->A, N->ctx);
    slong nm = 0;
    ulong (*mon)[RN_MAX_GEN] = malloc((size_t)(len ? len : 1) * sizeof *mon);
    fmpq_mpoly_struct* co = malloc((size_t)(len ? len : 1) * sizeof(fmpq_mpoly_struct));
    ulong e[NR_MAX_VARS + 1];
    fmpq_t c; fmpq_init(c);
    for (slong t = 0; t < len; t++) {
        fmpq_mpoly_get_term_exp_ui(e, a->A, t, N->ctx);
        fmpq_mpoly_get_term_coeff_fmpq(c, a->A, t, N->ctx);
        ulong m[RN_MAX_GEN] = { 0 };
        for (int i = 0; i < T->ng; i++) {
            if (T->g[i].elim) continue;        /* an ordinary variable of c_m */
            m[i] = e[T->g[i].v]; e[T->g[i].v] = 0;
        }
        slong k = 0;
        while (k < nm && memcmp(mon[k], m, sizeof m) != 0) k++;
        if (k == nm) { memcpy(mon[nm], m, sizeof m); fmpq_mpoly_init(co + nm, N->ctx); nm++; }
        fmpq_mpoly_push_term_fmpq_ui(co + k, c, e, N->ctx);
    }
    fmpq_clear(c);
    Expr** terms = malloc((size_t)(nm ? nm : 1) * sizeof(Expr*));
    int ok = 1;
    slong nt = 0;
    for (slong k = 0; k < nm && ok; k++) {
        fmpq_mpoly_sort_terms(co + k, N->ctx);
        fmpq_mpoly_combine_like_terms(co + k, N->ctx);
        /* c_m = a_m prod D_i^m_i / den */
        NrRf r; rf_init(N, &r);
        ok = rf_set_poly(N, &r, co + k);
        int64_t ex[RN_MAX_GEN];
        for (int i = 0; i < T->ng && ok; i++) {
            ex[i] = (int64_t)mon[k][i];
            for (size_t j = 0; j < T->g[i].B.nde; j++)
                pe_add(&r.nu, &r.nnu, &r.cnu, T->g[i].B.de[j].i, T->g[i].B.de[j].e * (slong)ex[i]);
        }
        for (size_t j = 0; j < a->nde && ok; j++) pe_add(&r.de, &r.nde, &r.cde, a->de[j].i, a->de[j].e);
        if (ok) ok = nr_cancel(N, &r);
        /* absorb a radicand factor: c y^m = (c B) y^(m - Q) or (c / B) y^(m + Q),
         * whichever shortens c (exact: y^Q = B for the same principal radical) */
        for (int i = 0; i < T->ng && ok; i++) {
            if (ex[i] == 0 || T->g[i].elim) continue;
            for (int dir = 1; dir >= -1; dir -= 2) {
                NrRf r2; rf_init(N, &r2);
                rf_copy(N, &r2, &r);
                NrRf b; rf_init(N, &b);
                rf_copy(N, &b, &T->g[i].B);
                int ok2 = dir > 0 || rf_pow(N, &b, -1);
                if (ok2) { rf_mul_into(N, &r2, &b); ok2 = nr_cancel(N, &r2); }
                bool take = ok2 && rf_degree(N, &r2) < rf_degree(N, &r);
                if (take) {
                    NrRf t = r; r = r2; r2 = t;
                    ex[i] -= dir * T->g[i].Q;
                }
                rf_clear(N, &b);
                rf_clear(N, &r2);
                if (take) break;
            }
        }
        if (ok) ok = nr_den_basis(N, &r);
        Expr* ce = ok ? nr_present(N, &r, 1) : NULL;
        rf_clear(N, &r);
        if (!ce) { ok = 0; break; }
        size_t nf = 1;
        Expr* f[RN_MAX_GEN + 1];
        f[0] = ce;
        for (int i = 0; i < T->ng; i++) {
            if (ex[i] == 0) continue;
            f[nf++] = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_copy(T->g[i].base), make_rational(ex[i], T->g[i].Q) }, 2);
        }
        terms[nt++] = nf == 1 ? f[0] : expr_new_function(expr_new_symbol(SYM_Times), f, nf);
    }
    for (slong k = 0; k < nm; k++) fmpq_mpoly_clear(co + k, N->ctx);
    free(co); free(mon);
    Expr* out = NULL;
    if (ok) {
        Expr* sum = nt == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt);
        out = eval_and_free(rn_elim_back(sum, T));
        expr_free(sum);
    }
    else for (slong k = 0; k < nt; k++) expr_free(terms[k]);
    free(terms);
    return out;
}

/* ---- an a-priori bound on [K:Q] ----
 *
 * The primitive element of K is a dense computation (a minimal polynomial of
 * a [K:Q]-square matrix); at degree 54 (2^(1/3), 3^(1/6), 37^(1/3)) it is
 * ~25 s, more than all the arithmetic after it.  The bound is read off the
 * atoms without building anything: a rational radical c^(p/q) lies in
 * Q(zeta_2q, p_i^(1/q_i)) for the primes p_i of c (q_i = q/gcd(q, e_i)), so
 * K lies in Q(zeta_2L, p^(1/L_p) ...) of degree <= phi(2L) prod L_p, times the
 * degree of a Root object and the index q of a nested radical c^(1/q) (whose
 * radicand c lies in the field of its own atoms, which are counted too). */
#define RN_MAX_FIELD_DEGREE 32

typedef struct {
    fmpz_t pr[16]; slong L[16]; int np;
    slong Lm1;        /* lcm of the root-of-unity orders / 2 */
    int   imag;       /* I is needed */
    double other;     /* product of the other atoms' degrees */
    int   over;       /* past the bound already */
    const Expr* seen[32]; int nseen;   /* the other atoms, each counted once */
} RnDeg;

/* First sight of the other atom e?  (Past 32 distinct ones: over.) */
static bool rn_deg_first(RnDeg* D, const Expr* e) {
    for (int i = 0; i < D->nseen; i++) if (expr_eq((Expr*)D->seen[i], (Expr*)e)) return false;
    if (D->nseen == 32) { D->over = 1; return false; }
    D->seen[D->nseen++] = e;
    return true;
}

static void rn_deg_prime(RnDeg* D, const fmpz_t p, slong q) {
    for (int i = 0; i < D->np; i++)
        if (fmpz_equal(D->pr[i], p)) { D->L[i] = (slong)lcm(D->L[i], q); return; }
    if (D->np == 16) { D->over = 1; return; }
    fmpz_init_set(D->pr[D->np], p);
    D->L[D->np++] = q;
}

static void rn_deg_rational_radical(RnDeg* D, const Expr* c, int64_t q) {
    fmpq_t v; fmpq_init(v);
    if (!fb_expr_to_fmpq(c, v)) { fmpq_clear(v); D->other *= (double)q; return; }
    if (fmpq_sgn(v) < 0) D->Lm1 = (slong)lcm(D->Lm1, q);
    for (int side = 0; side < 2 && !D->over; side++) {
        fmpz_t a; fmpz_init(a);
        fmpz_abs(a, side ? fmpq_denref(v) : fmpq_numref(v));
        if (fmpz_bits(a) > 128) { D->other *= (double)q; fmpz_clear(a); continue; }
        if (!fmpz_is_one(a)) {
            fmpz_factor_t F; fmpz_factor_init(F);
            fmpz_factor(F, a);
            for (slong i = 0; i < F->num && !D->over; i++) {
                slong qi = (slong)(q / gcd(q, (int64_t)F->exp[i]));
                if (qi > 1) rn_deg_prime(D, F->p + i, qi);
            }
            fmpz_factor_clear(F);
        }
        fmpz_clear(a);
    }
    fmpq_clear(v);
}

static void rn_deg_walk(RnDeg* D, const Expr* e) {
    if (D->over || e->type != EXPR_FUNCTION) return;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[1];
        if (!(b->type == EXPR_INTEGER && b->data.integer == 0)) D->imag = 1;
        return;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 && !nr_has_symbol(e) &&
        head_is(e->data.function.args[1], SYM_Rational)) {
        const Expr* qx = e->data.function.args[1]->data.function.args[1];
        int64_t q = qx->type == EXPR_INTEGER ? qx->data.integer : 64;
        const Expr* c = e->data.function.args[0];
        if (nr_is_ratnum(c)) rn_deg_rational_radical(D, c, q);
        else {
            /* a nested radical: c lies in the field of its own atoms, over
             * which c^(1/q) has degree at most q (b^(-p/q) is b^(p/q)'s
             * inverse: the base and index identify the atom) */
            if (rn_deg_first(D, c)) {
                D->other *= (double)q;
                rn_deg_walk(D, c);
            }
        }
        return;
    }
    if (head_is(e, SYM_Root)) {
        if (nr_has_symbol(e) || !rn_deg_first(D, e)) return;
        Expr* mp = flint_qqbar_gen_minpoly_coeffs(e);
        D->other *= (mp && mp->type == EXPR_FUNCTION) ? (double)(mp->data.function.arg_count - 1) : 1e9;
        if (mp) expr_free(mp);
        return;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) rn_deg_walk(D, e->data.function.args[i]);
}

/* Is [K:Q] provably at most RN_MAX_FIELD_DEGREE for the atoms of e and of
 * the extra constants? */
static bool rn_field_small(const Expr* e, const RnTab* T) {
    RnDeg D; memset(&D, 0, sizeof D);
    D.Lm1 = 1; D.other = 1.0;
    rn_deg_walk(&D, e);
    for (size_t i = 0; i < T->nextra; i++) rn_deg_walk(&D, T->extra[i]);
    double b = 1.0;
    for (int i = 0; i < D.np; i++) { b *= (double)D.L[i]; fmpz_clear(D.pr[i]); }
    ulong m = (ulong)(2 * D.Lm1);
    if (D.imag) m = (ulong)lcm((int64_t)m, 4);
    b *= (double)n_euler_phi(m);
    /* b is the degree of Q(zeta_2L, p^(1/L_p)) exactly when the rational
     * radicals are independent (the usual case: the 54 of 2^(1/3), 3^(1/6),
     * 37^(1/3) is real), while the factor of the nested radicals and Root
     * objects is loose (Sqrt[2 + Sqrt[6]] and Sqrt[4 + 2 Sqrt[6]] count twice),
     * so it gets more room. */
    return !D.over && b <= RN_MAX_FIELD_DEGREE && b * D.other <= 8 * RN_MAX_FIELD_DEGREE;
}

static bool nr_field_degree_small(const Expr* e) {
    RnTab T; memset(&T, 0, sizeof T);
    return rn_field_small(e, &T);
}

Expr* nfrat_radical_normal_form(const Expr* e, int* field) {
    if (field) *field = 0;
    if (!e || !nfrat_has_symbolic_radical(e)) return NULL;
    RnTab T; memset(&T, 0, sizeof T);
    Expr* out = NULL;
    Nr N; memset(&N, 0, sizeof N);
    int dim = 1;
    if (!rn_collect(e, &T, 0) || T.ng == 0) goto done;
    for (int i = 0; i < T.ng; i++) {
        if (T.g[i].elim) continue;
        if (T.g[i].Q > RN_MAX_DIM) goto done;
        dim *= (int)T.g[i].Q;
        if (dim > RN_MAX_DIM) goto done;
    }
    if (!rn_field_small(e, &T) || !rn_setup(&N, e, &T)) goto done;
    {
        RnE a; rn_init(&N, &a);
        if (rn_eval(&N, e, &a) && rn_cancel(&N, &a)) {
            out = rn_present(&N, &T, &a);
            if (out && field) *field = rn_kummer_field(&N, &T) ? 1 : 0;
        }
        rn_clear(&N, &a);
    }
done:
    rn_tab_clear(N.ctx ? &N : NULL, &T);
    nr_clear(&N);
    return out;
}

/* The canonical form over K(x1..xn) (nfrat_canonical): every constant
 * reduced in K, the written constant denominator u inverted over K, the
 * numerator and the denominator expanded, the denominator the product of the
 * monic pieces scaled to primitive integral product-basis coordinates (the
 * numerator takes the same rational scale).  Equal rational functions get
 * identical output, whatever the input's spelling. */
static Expr* nr_present_canon(Nr* N, NrRf* r) {
    if (rf_is_zero(N, r)) return expr_new_integer(0);
    fmpq_mpoly_t P, D, inv;
    fmpq_mpoly_init(P, N->ctx); fmpq_mpoly_init(D, N->ctx); fmpq_mpoly_init(inv, N->ctx);
    fmpq_mpoly_one(D, N->ctx);
    int ok = nr_const_inv(N, inv, r->u);
    if (ok) nr_mul(N, P, r->c, inv);
    for (size_t j = 0; j < r->nnu && ok; j++) ok = nr_mul_piece_pow(N, P, r->nu[j].i, r->nu[j].e);
    for (size_t j = 0; j < r->nde && ok; j++) ok = nr_mul_piece_pow(N, D, r->de[j].i, r->de[j].e);
    Expr* out = NULL;
    if (ok) {
        NrGroups GD, GP;
        nr_groups(N, &GD, D);
        nr_groups(N, &GP, P);
        fmpq_t s; fmpq_init(s);
        nr_groups_scale(N, &GD, s);
        Expr* num = nr_groups_expr(N, &GP, s);
        Expr* den = nr_groups_expr(N, &GD, s);
        out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2));
        fmpq_clear(s);
        nr_groups_clear(N, &GD); nr_groups_clear(N, &GP);
    }
    fmpq_mpoly_clear(P, N->ctx); fmpq_mpoly_clear(D, N->ctx); fmpq_mpoly_clear(inv, N->ctx);
    return out;
}

Expr* nfrat_canonical(const Expr* e) {
    if (!e || !nfrat_has_algebraic(e)) return NULL;

    Nr N;
    Expr* out = NULL;
    if (nr_setup_ex(&N, e, NULL, 1)) {
        NrRf r; rf_init(&N, &r);
        if (nr_eval(&N, e, &r) && nr_cancel(&N, &r)) out = nr_present_canon(&N, &r);
        rf_clear(&N, &r);
    }
    nr_clear(&N);
    return out;
}

#else  /* !USE_FLINT */

Expr* nfrat_together(const Expr* e, int factored, const Expr* ext) { (void)e; (void)factored; (void)ext; return NULL; }
Expr* nfrat_cancel(const Expr* e, int factored, const Expr* ext) { (void)e; (void)factored; (void)ext; return NULL; }
Expr* nfrat_radical_normal_form(const Expr* e, int* field) { (void)e; if (field) *field = 0; return NULL; }
static bool nr_field_degree_small(const Expr* e) { (void)e; return true; }
Expr* nfrat_canonical(const Expr* e) { (void)e; return NULL; }

#endif

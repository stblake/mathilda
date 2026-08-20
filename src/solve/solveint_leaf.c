/*
 * solveint_leaf.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"
#include "solveint_internal.h"

static void emit_solution(SearchState* st, int64_t leafval) {
    SICtx* c = st->ctx;
    if (st->nsol == st->cap) {
        st->cap = st->cap ? st->cap * 2 : 32;
        st->sols = (int64_t*)realloc(st->sols, sizeof(int64_t) * (size_t)st->cap * (size_t)c->n);
    }
    int64_t* row = st->sols + (size_t)st->nsol * (size_t)c->n;
    for (int i = 0; i < c->n; i++) row[i] = (i == st->leaf) ? leafval : st->val[i];
    st->nsol++;
}


/* Verify a full leaf assignment against the original conjunction. */
static bool verify_candidate(SearchState* st, int64_t leafval) {
    SICtx* c = st->ctx;
    int64_t vals[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) vals[i] = (i == st->leaf) ? leafval : st->val[i];
    return si_verify(c, vals);
}


/* Evaluate equation `eq`'s univariate-in-leaf coefficients a_0..a_d at the
 * current search assignment `st->val` (leaf left symbolic).  Reads the
 * precomputed monomials directly -- no MPoly allocation, the search hot
 * path.  Returns the leaf degree, -1 for identically zero, -2 if the leaf
 * degree exceeds the solver's reach. */
static int eval_leaf_coeffs(SearchState* st, const MPoly* eq, int leaf, mpz_t* a, mpz_t term) {
    SICtx* c = st->ctx;
    int n = c->n;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_set_ui(a[k], 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int el = ex[leaf];
        if (el > SI_LEAF_MAXDEG) return -2;
        mpz_set(term, eq->coefs[t]);
        for (int u = 0; u < n; u++) {
            if (u == leaf) continue;
            int e = ex[u];
            for (int k = 0; k < e; k++) mpz_mul_si(term, term, (long)st->val[u]);
        }
        mpz_add(a[el], a[el], term);
    }
    for (int k = SI_LEAF_MAXDEG; k >= 0; k--) if (mpz_sgn(a[k]) != 0) return k;
    return -1;
}


/* Build the per-equation int64 coefficient cache so the hot leaf loop never
 * calls mpz_get_si / mpz_fits_slong_p per node.  Idempotent. */
void si_build_leaf_cache(SearchState* st) {
    if (st->eqc_built) return;
    SICtx* c = st->ctx;
    for (int q = 0; q < c->neq; q++) {
        const MPoly* eq = c->eq[q];
        st->eqc[q] = (int64_t*)malloc(sizeof(int64_t) * (eq->n_terms ? eq->n_terms : 1));
        bool ok = true;
        for (size_t t = 0; t < eq->n_terms; t++) {
            if (mpz_fits_slong_p(eq->coefs[t])) st->eqc[q][t] = mpz_get_si(eq->coefs[t]);
            else { ok = false; break; }
        }
        st->eqc_ok[q] = ok;
    }
    st->eqc_built = true;
}

void si_free_leaf_cache(SearchState* st) {
    if (!st->eqc_built) return;
    for (int q = 0; q < st->ctx->neq; q++) { free(st->eqc[q]); st->eqc[q] = NULL; }
    st->eqc_built = false;
}


/* Integer roots of a2 x^2 + a1 x + a0 (degree d in {1,2}) inside [lo,hi], in
 * pure int64.  Sets *overflow and returns 0 when an intermediate would exceed
 * int64 (caller then falls back to the GMP path).  Returns the root count. */
static int si_uniroots_i64(int64_t a2, int64_t a1, int64_t a0, int d,
                           int64_t lo, int64_t hi, int64_t* roots, bool* overflow) {
    int nr = 0;
    *overflow = false;
    if (d == 1) {                                  /* a1 x + a0 = 0 */
        if (a1 != 0 && a0 % a1 == 0) {
            int64_t r = -a0 / a1;
            if (r >= lo && r <= hi) roots[nr++] = r;
        }
        return nr;
    }
    /* d == 2: discriminant a1^2 - 4 a2 a0. */
    int64_t a1sq, fac, disc, den;
    if (ci_mul_i64(a1, a1, &a1sq) || ci_mul_i64(a2, a0, &fac)
        || ci_mul_i64(fac, 4, &fac) || ci_sub_i64(a1sq, fac, &disc)
        || ci_mul_i64(a2, 2, &den)) { *overflow = true; return 0; }
    if (disc < 0) return 0;
    int64_t s = si_isqrt_i64(disc);
    int64_t s2;
    if (ci_mul_i64(s, s, &s2)) { *overflow = true; return 0; }
    if (s2 != disc) return 0;                      /* not a perfect square */
    for (int sign = -1; sign <= 1; sign += 2) {
        int64_t num = -a1 + sign * s;              /* |a1|,|s| already fit; sum fits */
        if (den != 0 && num % den == 0) {
            int64_t r = num / den;
            if (r >= lo && r <= hi) {
                bool dup = false;
                for (int i = 0; i < nr; i++) if (roots[i] == r) dup = true;
                if (!dup) roots[nr++] = r;
            }
        }
    }
    return nr;
}


/* Evaluate equation `eq`'s univariate-in-leaf coefficients at st->val entirely
 * in int64 and return the leaf roots in [lo,hi].  Return codes mirror the GMP
 * path: -1 = eq is identically zero in the leaf (no constraint), -2 = a nonzero
 * constant (infeasible), >=0 = root count.  Sets *need_mpz (caller falls back to
 * the exact GMP path) when the leaf degree exceeds 2 or any product/coefficient
 * would exceed int64.  This is the hot path for large ordered boxes. */
static int si_leaf_roots_i64(SearchState* st, const MPoly* eq, const int64_t* ci,
                             bool ci_ok, int leaf, int64_t lo, int64_t hi,
                             int64_t* roots, bool* need_mpz) {
    SICtx* c = st->ctx;
    int n = c->n;
    int64_t a0 = 0, a1 = 0, a2 = 0;
    *need_mpz = false;
    if (!ci_ok) { *need_mpz = true; return 0; }        /* a coefficient overflows int64 */
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int el = ex[leaf];
        if (el > 2) { *need_mpz = true; return 0; }
        int64_t val = ci[t];                            /* precomputed int64 coefficient */
        for (int u = 0; u < n; u++) {
            if (u == leaf) continue;
            for (int k = 0; k < ex[u]; k++)
                if (ci_mul_i64(val, st->val[u], &val)) { *need_mpz = true; return 0; }
        }
        int64_t* slot = (el == 0) ? &a0 : (el == 1) ? &a1 : &a2;
        if (ci_add_i64(*slot, val, slot)) { *need_mpz = true; return 0; }
    }
    int d = (a2 != 0) ? 2 : (a1 != 0) ? 1 : 0;
    if (d == 0) return (a0 != 0) ? -2 : -1;
    bool ovf;
    int nr = si_uniroots_i64(a2, a1, a0, d, lo, hi, roots, &ovf);
    if (ovf) { *need_mpz = true; return 0; }
    return nr;
}


/* Solve the leaf variable given a full search-var assignment in st->val. */
static void solve_leaf(SearchState* st) {
    SICtx* c = st->ctx;
    int leaf = st->leaf;
    int64_t lo = c->lo[leaf], hi = c->hi[leaf];
    if (++st->visits > st->max_visits) { st->overflow = true; return; }

    /* int64 fast path: no GMP allocation per node (the hot path for large
     * ordered boxes).  Falls back to the exact GMP body below whenever a leaf
     * degree exceeds 2 or an intermediate would exceed int64. */
    {
        int64_t inter[512]; int ninter = 0; bool have = false, feasible = true, need_mpz = false;
        for (int q = 0; q < c->neq && feasible; q++) {
            int64_t roots[8];
            int nr = si_leaf_roots_i64(st, c->eq[q], st->eqc[q], st->eqc_ok[q],
                                       st->leaf, lo, hi, roots, &need_mpz);
            if (need_mpz) break;
            if (nr == -1) continue;                  /* no constraint on leaf */
            if (nr == -2) { feasible = false; break; } /* nonzero constant: infeasible */
            if (!have) {
                for (int i = 0; i < nr && ninter < 512; i++) inter[ninter++] = roots[i];
                have = true;
            } else {
                int64_t keep[512]; int nk = 0;
                for (int i = 0; i < ninter; i++)
                    for (int j = 0; j < nr; j++)
                        if (inter[i] == roots[j]) { keep[nk++] = inter[i]; break; }
                memcpy(inter, keep, sizeof(int64_t) * (size_t)nk);
                ninter = nk;
            }
            if (ninter == 0) feasible = false;
        }
        if (!need_mpz) {
            if (feasible) {
                if (have) {
                    for (int i = 0; i < ninter; i++)
                        if (verify_candidate(st, inter[i])) emit_solution(st, inter[i]);
                } else {
                    if (hi - lo > 10000000LL) st->overflow = true;
                    else for (int64_t r = lo; r <= hi; r++)
                        if (verify_candidate(st, r)) emit_solution(st, r);
                }
            }
            return;
        }
    }

    int64_t inter[512]; int ninter = 0; bool have = false;
    mpz_t a[SI_LEAF_MAXDEG + 1], term;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_init(a[k]);
    mpz_init(term);

    bool feasible = true;
    for (int q = 0; q < c->neq && feasible; q++) {
        int d = eval_leaf_coeffs(st, c->eq[q], st->leaf, a, term);
        if (d == -1) continue;                       /* 0 == 0: no constraint */
        if (d == -2) { st->overflow = true; feasible = false; break; }
        if (d == 0) {                                /* constant == 0 ? */
            if (mpz_sgn(a[0]) != 0) feasible = false;
            continue;
        }
        int64_t roots[SI_LEAF_MAXDEG * 2];
        int nr = univariate_roots(a, d, lo, hi, roots);
        if (!have) {
            for (int i = 0; i < nr && ninter < 512; i++) inter[ninter++] = roots[i];
            have = true;
        } else {
            int64_t keep[512]; int nk = 0;
            for (int i = 0; i < ninter; i++)
                for (int j = 0; j < nr; j++)
                    if (inter[i] == roots[j]) { keep[nk++] = inter[i]; break; }
            memcpy(inter, keep, sizeof(int64_t) * (size_t)nk);
            ninter = nk;
        }
        if (ninter == 0) { feasible = false; }
    }

    if (feasible && !st->overflow) {
        if (have) {
            for (int i = 0; i < ninter; i++)
                if (verify_candidate(st, inter[i])) emit_solution(st, inter[i]);
        } else {
            /* No equation constrains the leaf at this branch: every value in
             * its range satisfying the constraints is a solution.  Only safe
             * for a genuinely small range; otherwise abort. */
            if (hi - lo > 10000000LL) { st->overflow = true; }
            else for (int64_t r = lo; r <= hi; r++)
                if (verify_candidate(st, r)) emit_solution(st, r);
        }
    }
    mpz_clear(term);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_clear(a[k]);
}


/* Effective [lo,hi] for search var at position, tightened by orderings
 * against already-assigned variables. */
static void effective_bounds(SearchState* st, int depth, int vi, int64_t* elo, int64_t* ehi) {
    SICtx* c = st->ctx;
    *elo = c->lo[vi]; *ehi = c->hi[vi];
    (void)depth;
    for (int k = 0; k < c->n_ord; k++) {
        int a = c->ord_a[k], b = c->ord_b[k]; int s = c->ord_strict[k] ? 1 : 0;
        /* var a < / <= var b */
        if (a == vi) {                               /* vi < b: hi <= val_b - s */
            bool assigned = false; int64_t vb = 0;
            for (int p = 0; p < depth; p++) if (st->order[p] == b) { assigned = true; vb = st->val[b]; }
            if (assigned && vb - s < *ehi) *ehi = vb - s;
        }
        if (b == vi) {                               /* a < vi: lo >= val_a + s */
            bool assigned = false; int64_t va = 0;
            for (int p = 0; p < depth; p++) if (st->order[p] == a) { assigned = true; va = st->val[a]; }
            if (assigned && va + s > *elo) *elo = va + s;
        }
    }
    /* Abs-orderings prune only the smaller side: |vi| < |val_b| once the larger
     * variable b is assigned, tightening vi to a symmetric window.  (The larger
     * side is a hole around 0, not a contiguous narrowing, so it is left to the
     * final verify.)  Realised only when b is enumerated before vi -- see the
     * abs-aware search-order builder. */
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k]; int s = c->abs_ord_strict[k] ? 1 : 0;
        if (a != vi) continue;
        bool assigned = false; int64_t vb = 0;
        for (int p = 0; p < depth; p++) if (st->order[p] == b) { assigned = true; vb = st->val[b]; }
        if (!assigned) continue;
        int64_t avb = vb < 0 ? -vb : vb;
        int64_t bnd = avb - s;                       /* |vi| <= bnd */
        if (bnd < 0) bnd = -1;
        if ( bnd < *ehi) *ehi =  bnd;
        if (-bnd > *elo) *elo = -bnd;
    }
}


static void resolve_peeled(SearchState* st, int pi);   /* multi-leaf resolver */
void search_rec(SearchState* st, int depth) {
    if (st->overflow) return;
    if (depth == st->n_search) {
        if (st->multileaf) resolve_peeled(st, 0); else solve_leaf(st);
        return;
    }

    int vi = st->order[depth];
    int64_t elo, ehi;
    effective_bounds(st, depth, vi, &elo, &ehi);

    for (int64_t r = elo; r <= ehi && !st->overflow; r++) {
        st->val[vi] = r;
        search_rec(st, depth + 1);
    }
}


/* Length of the longest <= / < ordering chain among the search vars `vars`.
 * A chain of length L cuts an otherwise-N^L box to ~N^L / L! ordered tuples,
 * so the search-space guard divides the raw estimate by L! -- otherwise an
 * ordered box like  0 < x <= y <= z < 1000  (raw 10^9, ordered ~1.7*10^8) is
 * wrongly declined.  The visit cap backstops any under-estimate, so this only
 * ever lets a genuinely-ordered box be attempted, never truncates a result. */
int si_longest_chain(const SICtx* c, const int* vars, int nv) {
    bool in[SI_MAX_VARS];
    int memo[SI_MAX_VARS];
    for (int i = 0; i < SI_MAX_VARS; i++) { in[i] = false; memo[i] = 1; }
    for (int i = 0; i < nv; i++) in[vars[i]] = true;
    /* memo[v] = longest chain starting at v; relax to a fixpoint (DAG, tiny). */
    for (int iter = 0; iter < nv; iter++) {
        for (int i = 0; i < nv; i++) {
            int v = vars[i], m = 1;
            for (int e = 0; e < c->n_ord; e++) {
                if (c->ord_a[e] == v && in[c->ord_b[e]]) {   /* edge v < / <= b */
                    int cand = 1 + memo[c->ord_b[e]];
                    if (cand > m) m = cand;
                }
            }
            for (int e = 0; e < c->n_abs_ord; e++) {         /* |v| < / <= |b| */
                if (c->abs_ord_a[e] == v && in[c->abs_ord_b[e]]) {
                    int cand = 1 + memo[c->abs_ord_b[e]];
                    if (cand > m) m = cand;
                }
            }
            if (m > memo[v]) memo[v] = m;
        }
    }
    int best = 1;
    for (int i = 0; i < nv; i++) if (memo[vars[i]] > best) best = memo[vars[i]];
    return best;
}


/* ================================================================== *
 *  Modular-sieved leaf search for a large (often non-separable) box.
 *
 *  When the ordering-pruned leaf box is still above SI_MAX_NODES the
 *  ordinary search declines.  This raises that ceiling for a single
 *  polynomial equation by pruning the innermost enumerated variable to
 *  the residues (mod a small M) for which the leaf equation can vanish
 *  -- a sound necessary condition (never drops a solution), so the box
 *  is still searched EXHAUSTIVELY: it returns the complete set / a proven
 *  {} when it finishes within the raised budget, and declines (never a
 *  partial list) otherwise.  Reuses the leaf/order/cache the caller has
 *  already built. */

/* P mod M at an integer residue assignment res[0..n). */
static int si_eval_mpoly_mod(const MPoly* p, const int* res, int M) {
    long long acc = 0;
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* ex = p->exps + t * (size_t)p->n_vars;
        long long term = (long long)mpz_fdiv_ui(p->coefs[t], (unsigned long)M);
        for (int v = 0; v < p->n_vars; v++) {
            int rv = res[v] % M; if (rv < 0) rv += M;
            for (int e = 0; e < ex[v]; e++) term = (term * rv) % M;
        }
        acc = (acc + term) % M;
    }
    return (int)acc;                                     /* in [0, M) */
}

static void si_modsieve_rec(SearchState* st, int depth, int M,
                            const unsigned char* feas, const long long* stride) {
    if (st->overflow) return;
    if (depth == st->n_search) { solve_leaf(st); return; }
    int vi = st->order[depth];
    int64_t elo, ehi; effective_bounds(st, depth, vi, &elo, &ehi);
    if (depth == st->n_search - 1) {                     /* innermost: prune by residue */
        long long base = 0;
        for (int i = 0; i < depth; i++) {
            int r = (int)(st->val[st->order[i]] % M); if (r < 0) r += M;
            base += (long long)r * stride[i];
        }
        long long s = stride[depth];
        for (int64_t v = elo; v <= ehi && !st->overflow; v++) {
            int r = (int)(v % M); if (r < 0) r += M;
            if (!feas[base + (long long)r * s]) continue;
            st->val[vi] = v;
            solve_leaf(st);
        }
        return;
    }
    for (int64_t v = elo; v <= ehi && !st->overflow; v++) {
        st->val[vi] = v;
        si_modsieve_rec(st, depth + 1, M, feas, stride);
    }
}

bool si_solve_box_modsieve(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int nsrch = st->n_search;
    if (nsrch < 1 || nsrch > 4) return false;            /* keep the feas table small */
    const MPoly* eq = c->eq[0];
    int leaf = st->leaf;

    /* M = product of small primes with M^nsrch <= ~4e6 (cache-resident table). */
    int Mmax = (int)powl(4.0e6L, 1.0L / (long double)nsrch);
    static const int primes[] = {2,3,5,7,11,13,17,19,23,29,31};
    int M = 1;
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); i++)
        if ((long long)M * primes[i] <= Mmax) M *= primes[i];
    if (M < 6) return false;

    long long tsz = 1; for (int i = 0; i < nsrch; i++) tsz *= M;
    unsigned char* feas = (unsigned char*)calloc((size_t)tsz, 1);
    if (!feas) return false;

    /* feas[idx] = does some leaf residue make P == 0 mod M, for the outer residue
     * tuple encoded by idx (digit i is the residue of st->order[i]). */
    int res[SI_MAX_VARS]; for (int i = 0; i < c->n; i++) res[i] = 0;
    int digit[SI_MAX_VARS]; for (int i = 0; i < nsrch; i++) digit[i] = 0;
    long long idx = 0;
    for (;;) {
        for (int i = 0; i < nsrch; i++) res[st->order[i]] = digit[i];
        bool ok = false;
        for (int L = 0; L < M && !ok; L++) { res[leaf] = L; if (si_eval_mpoly_mod(eq, res, M) == 0) ok = true; }
        feas[idx++] = ok ? 1 : 0;
        int d = 0; for (; d < nsrch; d++) { if (++digit[d] < M) break; digit[d] = 0; }
        if (d == nsrch) break;
    }

    long long stride[SI_MAX_VARS]; { long long s = 1; for (int i = 0; i < nsrch; i++) { stride[i] = s; s *= M; } }
    st->max_visits = 3000000000LL;                       /* raised ceiling; the sieve keeps real visits down */
    st->overflow = false;
    si_build_leaf_cache(st);
    si_modsieve_rec(st, 0, M, feas, stride);
    si_free_leaf_cache(st);
    free(feas);
    if (st->overflow) { st->nsol = 0; return false; }    /* did not exhaust -> decline (no partial list) */
    return true;                                         /* complete set (or proven {}) */
}


/* ------------------------------------------------------------------ *
 *  A2: multi-leaf staged elimination.                                 *
 *                                                                     *
 *  The single-leaf search enumerates every variable but one.  A system *
 *  like the Euler brick  x^2+y^2==a^2 && x^2+z^2==b^2 && y^2+z^2==c^2   *
 *  has THREE variables (a,b,c) each determined by a single equation --  *
 *  enumerating even two of them is hopeless.  Staged elimination peels  *
 *  every variable that appears in exactly one equation and is           *
 *  univariate-solvable there, resolving it per free-variable assignment *
 *  (an exact root, branching on multiplicity) instead of enumerating.   *
 *  Only the genuinely-coupled "free" variables are walked.              *
 * ------------------------------------------------------------------ */

/* Resolve peeled variable peel[pi] from its equation given the free-variable
 * assignment (and earlier-resolved peels) in st->val, recursing to the next
 * peel; at the end, verify the full assignment and emit.  Each peel_eq contains
 * only its own peeled variable plus free variables, so any resolution order is
 * valid. */
static void resolve_peeled(SearchState* st, int pi) {
    if (st->overflow) return;
    SICtx* c = st->ctx;
    if (pi == st->n_peel) {
        if (++st->visits > st->max_visits) { st->overflow = true; return; }
        if (si_verify(c, st->val)) emit_full(st, st->val);
        return;
    }
    int v = st->peel[pi], e = st->peel_eq[pi];
    int64_t lo = c->has_lo[v] ? c->lo[v] : -(1LL << 50);
    int64_t hi = c->has_hi[v] ? c->hi[v] :  (1LL << 50);

    /* int64 fast path (peeled var degree <= 2). */
    int64_t roots[8]; bool need_mpz = false;
    int nr = si_leaf_roots_i64(st, c->eq[e], st->eqc[e], st->eqc_ok[e], v, lo, hi, roots, &need_mpz);
    if (!need_mpz) {
        if (nr == -2) return;                       /* infeasible at this branch */
        if (nr == -1) { st->overflow = true; return; } /* v unconstrained here: decline */
        for (int i = 0; i < nr && !st->overflow; i++) { st->val[v] = roots[i]; resolve_peeled(st, pi + 1); }
        return;
    }

    /* exact GMP fallback (peeled var degree up to SI_LEAF_MAXDEG). */
    mpz_t a[SI_LEAF_MAXDEG + 1], term;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_init(a[k]);
    mpz_init(term);
    int d = eval_leaf_coeffs(st, c->eq[e], v, a, term);
    if (d == -2) st->overflow = true;               /* degree beyond the solver's reach */
    else if (d == -1) st->overflow = true;          /* v unconstrained here: decline */
    else if (d == 0) { /* constant: feasible iff 0, but then v is free -> decline */
        if (mpz_sgn(a[0]) != 0) { /* infeasible: no roots */ } else st->overflow = true;
    } else {
        int64_t r2[SI_LEAF_MAXDEG * 2];
        int nr2 = univariate_roots(a, d, lo, hi, r2);
        for (int i = 0; i < nr2 && !st->overflow; i++) { st->val[v] = r2[i]; resolve_peeled(st, pi + 1); }
    }
    mpz_clear(term);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_clear(a[k]);
}


/* Detect a staged-elimination structure and, if the free variables form a
 * tractable ordered box, run the multi-leaf search.  Returns true when it
 * settled the answer (candidates in st->sols); false to fall back. */
bool si_solve_multileaf(SICtx* c, SearchState* st) {
    int n = c->n;
    if (c->neq < 2) return false;                   /* single-leaf path handles neq == 1 */

    /* appear[v] = number of equations that contain v (deg > 0). */
    int appear[SI_MAX_VARS], only_eq[SI_MAX_VARS];
    for (int v = 0; v < n; v++) { appear[v] = 0; only_eq[v] = -1; }
    for (int q = 0; q < c->neq; q++)
        for (int v = 0; v < n; v++)
            if (mpoly_deg_var(c->eq[q], v) > 0) { appear[v]++; only_eq[v] = q; }

    /* A variable is a peel candidate iff it appears in exactly one equation and
     * is univariate-solvable there.  Each equation is claimed by at most one
     * peeled variable (widest domain / unbounded first); losers become free. */
    bool peeled[SI_MAX_VARS]; int claim[SI_MAX_VARS * 2];
    for (int v = 0; v < n; v++) peeled[v] = false;
    for (int q = 0; q < c->neq; q++) claim[q] = -1;
    for (int v = 0; v < n; v++) {
        if (appear[v] != 1) continue;
        int e = only_eq[v];
        int dg = mpoly_deg_var(c->eq[e], v);
        if (dg < 1 || dg > SI_LEAF_MAXDEG) continue;
        int cur = claim[e];
        if (cur < 0) { claim[e] = v; continue; }
        /* Prefer to peel the wider / unbounded variable (enumerate the narrow). */
        long double wv = (c->has_lo[v] && c->has_hi[v]) ? (long double)(c->hi[v] - c->lo[v]) : 1e30L;
        long double wc = (c->has_lo[cur] && c->has_hi[cur]) ? (long double)(c->hi[cur] - c->lo[cur]) : 1e30L;
        if (wv > wc) claim[e] = v;
    }
    st->n_peel = 0;
    for (int q = 0; q < c->neq; q++)
        if (claim[q] >= 0) { int v = claim[q]; peeled[v] = true; st->peel[st->n_peel] = v; st->peel_eq[st->n_peel] = q; st->n_peel++; }
    if (st->n_peel < 2) return false;               /* single leaf: ordinary path */

    /* Free variables: everything not peeled.  All must be finitely bounded. */
    int free_v[SI_MAX_VARS], nfree = 0;
    for (int v = 0; v < n; v++) if (!peeled[v]) {
        if (!(c->has_lo[v] && c->has_hi[v])) return false;
        free_v[nfree++] = v;
    }
    if (nfree == 0) return false;

    /* Ordered free-box estimate (raw product / longest-chain factorial). */
    long double est = 1.0L;
    for (int i = 0; i < nfree; i++) est *= (long double)(c->hi[free_v[i]] - c->lo[free_v[i]] + 1);
    int L = si_longest_chain(c, free_v, nfree);
    for (int i = 2; i <= L; i++) est /= (long double)i;
    if (est > (long double)SI_MAX_NODES) return false;

    /* Order free vars ascending domain (narrowest innermost). */
    st->n_search = 0;
    for (int i = 0; i < nfree; i++) st->order[st->n_search++] = free_v[i];
    for (int i = 0; i + 1 < st->n_search; i++) {
        int pk = i;
        for (int j = i + 1; j < st->n_search; j++)
            if ((c->hi[st->order[j]] - c->lo[st->order[j]]) < (c->hi[st->order[pk]] - c->lo[st->order[pk]])) pk = j;
        if (pk != i) { int t = st->order[i]; st->order[i] = st->order[pk]; st->order[pk] = t; }
    }
    st->multileaf = true;
    st->max_visits = SI_MAX_NODES;
    si_build_leaf_cache(st);
    for (int i = 0; i < n; i++) st->val[i] = 0;
    search_rec(st, 0);
    si_free_leaf_cache(st);
    return !st->overflow;
}

/*
 * solveint_reciprocal.c
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


/* Next lexicographic permutation of a[0..k) in place.  Starting from an
 * ascending-sorted array it enumerates every DISTINCT permutation (equal
 * elements never produce a duplicate).  Returns false after the last one. */
static bool si_next_permutation_i64(int64_t* a, int k) {
    if (k < 2) return false;
    int i = k - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return false;
    int j = k - 1;
    while (a[j] <= a[i]) j--;
    int64_t t = a[i]; a[i] = a[j]; a[j] = t;
    for (int lo = i + 1, hi = k - 1; lo < hi; lo++, hi--) { t = a[lo]; a[lo] = a[hi]; a[hi] = t; }
    return true;
}

/* The reciprocal equation is fully symmetric in its active variables (all a_i
 * equal -- verified by si_reciprocal_detect).  When the caller supplied no
 * ordering to pin a canonical representative, an ascending solution stands for
 * its whole permutation orbit; emit every distinct permutation so the returned
 * set matches the unordered answer.  Each permuted tuple is re-verified, so any
 * disequation / individual bound still filters correctly. */
static void si_recip_emit_perms(SICtx* c, SearchState* st, const int* order, int k,
                                const int64_t* full) {
    int64_t vals[SI_MAX_VARS], tmp[SI_MAX_VARS];
    for (int j = 0; j < k; j++) vals[j] = full[order[j]];    /* ascending */
    for (int i = 0; i < c->n; i++) tmp[i] = full[i];
    do {
        for (int j = 0; j < k; j++) tmp[order[j]] = vals[j];
        if (si_verify(c, tmp)) emit_full(st, tmp);
    } while (si_next_permutation_i64(vals, k));
}


/* --- Reciprocal (Egyptian-fraction) recursion. --- */
static void si_recip_rec(SICtx* c, SearchState* st, const int* order, int k,
                         int pos, const mpz_t p, const mpz_t q, int64_t prev,
                         int64_t* full, bool emit_perms) {
    if (st->overflow) return;
    if (mpz_sgn(p) <= 0) return;                 /* R must stay positive */
    int m = k - pos;                             /* remaining terms */
    int vi = order[pos];

    if (m == 1) {                                /* 1/v == p/q -> v = q/p */
        if (mpz_divisible_p(q, p)) {
            mpz_t vv; mpz_init(vv); mpz_divexact(vv, q, p);
            if (mpz_fits_slong_p(vv)) {
                int64_t v = mpz_get_si(vv);
                if (v >= prev) {
                    full[vi] = v;
                    if (emit_perms) si_recip_emit_perms(c, st, order, k, full);
                    else if (si_verify(c, full)) emit_full(st, full);
                }
            }
            mpz_clear(vv);
        }
        return;
    }

    /* 1/v < p/q  => v >= floor(q/p)+1;  1/v >= (p/q)/m => v <= floor(m q / p). */
    mpz_t lo_m, hi_m, tmp; mpz_init(lo_m); mpz_init(hi_m); mpz_init(tmp);
    mpz_fdiv_q(lo_m, q, p); mpz_add_ui(lo_m, lo_m, 1);
    mpz_mul_ui(tmp, q, (unsigned long)m); mpz_fdiv_q(hi_m, tmp, p);
    int64_t lo = mpz_fits_slong_p(lo_m) ? mpz_get_si(lo_m) : INT64_MAX;
    int64_t hi = mpz_fits_slong_p(hi_m) ? mpz_get_si(hi_m) : INT64_MAX;
    mpz_clear(lo_m); mpz_clear(hi_m); mpz_clear(tmp);
    if (prev > lo) lo = prev;

    for (int64_t v = lo; v <= hi && !st->overflow; v++) {
        if (++st->visits > st->max_visits) { st->overflow = true; return; }
        mpz_t np, nq, g; mpz_init(np); mpz_init(nq); mpz_init(g);
        mpz_mul_si(np, p, (long)v); mpz_sub(np, np, q);   /* p*v - q */
        mpz_mul_si(nq, q, (long)v);                        /* q*v */
        if (mpz_sgn(np) > 0) {
            mpz_gcd(g, np, nq); mpz_divexact(np, np, g); mpz_divexact(nq, nq, g);
            full[vi] = v;
            si_recip_rec(c, st, order, k, pos + 1, np, nq, v, full, emit_perms);
        }
        mpz_clear(np); mpz_clear(nq); mpz_clear(g);
    }
}


/* Detect the shape  c_full * prod x_i - a * sum_i prod_{j!=i} x_j == 0
 * (equivalently  a * sum 1/x_i == c_full), with all a_i equal.  Fills the
 * active-variable list, c_full and a. */
static bool si_reciprocal_detect(const MPoly* eq, int n, int* active, int* ka_out,
                                 mpz_t c_full, mpz_t a_out) {
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++) if (ex[v] > 1) return false;   /* multilinear only */
    }
    int ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) active[ka++] = v;
    if (ka < 2 || eq->n_terms != (size_t)(ka + 1)) return false;

    int* ex = (int*)calloc((size_t)n, sizeof(int));
    for (int j = 0; j < ka; j++) ex[active[j]] = 1;
    const mpz_t* cf = mpoly_get_coef(eq, ex);
    if (!cf) { free(ex); return false; }
    mpz_set(c_full, *cf);

    bool ok = true, first = true;
    mpz_t a0; mpz_init(a0);
    for (int j = 0; j < ka && ok; j++) {
        ex[active[j]] = 0;
        const mpz_t* cc = mpoly_get_coef(eq, ex);
        ex[active[j]] = 1;
        if (!cc) { ok = false; break; }
        if (first) { mpz_neg(a0, *cc); first = false; }
        else {
            mpz_t ai; mpz_init(ai); mpz_neg(ai, *cc);
            if (mpz_cmp(ai, a0) != 0) ok = false;
            mpz_clear(ai);
        }
    }
    free(ex);
    if (!ok) { mpz_clear(a0); return false; }
    if (mpz_sgn(c_full) < 0) { mpz_neg(c_full, c_full); mpz_neg(a0, a0); }
    if (mpz_sgn(c_full) <= 0 || mpz_sgn(a0) <= 0) { mpz_clear(a0); return false; }
    mpz_set(a_out, a0); mpz_clear(a0);
    *ka_out = ka;
    return true;
}

bool si_solve_reciprocal(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int active[SI_MAX_VARS], ka = 0;
    mpz_t c_full, a; mpz_init(c_full); mpz_init(a);
    if (!si_reciprocal_detect(c->eq[0], c->n, active, &ka, c_full, a)) {
        mpz_clear(c_full); mpz_clear(a); return false;
    }
    int order[SI_MAX_VARS];
    bool emit_perms = false;
    if (!si_build_total_order(c, active, ka, order)) {
        /* No ordering chain among the active vars.  If NO ordering constraint
         * touches them at all, the symmetric equation still has a well-defined
         * unordered solution set: enumerate the ascending canonical reps and
         * emit their permutations.  A PARTIAL order we cannot canonicalize
         * safely -> decline. */
        bool touches = false;
        for (int e = 0; e < c->n_ord && !touches; e++)
            for (int j = 0; j < ka; j++)
                if (c->ord_a[e] == active[j] || c->ord_b[e] == active[j]) { touches = true; break; }
        if (touches) { mpz_clear(c_full); mpz_clear(a); return false; }
        for (int j = 0; j < ka; j++) order[j] = active[j];   /* canonical order */
        emit_perms = true;
    }
    int64_t full[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) full[i] = 0;
    /* p/q = R = c_full / a. */
    st->max_visits = SI_MAX_NODES;
    si_recip_rec(c, st, order, ka, 0, c_full, a, 1, full, emit_perms);
    mpz_clear(c_full); mpz_clear(a);
    return true;
}

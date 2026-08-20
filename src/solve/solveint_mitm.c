/*
 * solveint_mitm.c
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


/* ------------------------------------------------------------------ *
 *  Meet-in-the-middle for a separable additive equation.              *
 *                                                                     *
 *  When the single equation is  sum_i g_i(x_i) == TARGET  with each    *
 *  g_i univariate (every term touches at most one variable), the       *
 *  ~N^(n-1) leaf search is replaced by splitting the variables into    *
 *  two groups, tabulating one group's partial sums and binary-         *
 *  searching the other -- ~N^ceil(n/2) work.  Correctness rests on the *
 *  same final verify_candidate against the original conjunction, so    *
 *  cross-group orderings / disequations need no special handling here. *
 * ------------------------------------------------------------------ */

typedef struct { int64_t sum; int64_t vals[SI_MAX_VARS]; } MitmEntry;

#define MITM_HASH_CAP 5000000LL
static int mitm_cmp(const void* pa, const void* pb) {
    int64_t a = ((const MitmEntry*)pa)->sum, b = ((const MitmEntry*)pb)->sum;
    return a < b ? -1 : (a > b ? 1 : 0);
}


/* Attempt the meet-in-the-middle path.  Returns true if it ran (results are
 * in st->sols); false to fall back to the general search. */
bool mitm_solve(SearchState* st) {
    SICtx* c = st->ctx;
    int n = c->n;
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0];

    /* Separability + constant term. */
    mpz_t c0; mpz_init(c0); mpz_set_ui(c0, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nv = 0;
        for (int v = 0; v < n; v++) if (ex[v] > 0) nv++;
        if (nv >= 2) { mpz_clear(c0); return false; }   /* not separable */
        if (nv == 0) mpz_add(c0, c0, eq->coefs[t]);
    }

    /* Fully bounded? */
    for (int i = 0; i < n; i++)
        if (!(c->has_lo[i] && c->has_hi[i]) || c->hi[i] < c->lo[i]) { mpz_clear(c0); return false; }

    /* TARGET = -c0 as int64. */
    mpz_neg(c0, c0);
    if (!mpz_fits_slong_p(c0)) { mpz_clear(c0); return false; }
    int64_t target = mpz_get_si(c0);
    mpz_clear(c0);

    /* Per-variable value tables g_i(v). */
    int64_t* gtab[SI_MAX_VARS];
    for (int i = 0; i < n; i++) gtab[i] = NULL;
    int64_t domain[SI_MAX_VARS];
    long double total_tab = 0.0L;
    bool ok = true;
    mpz_t gv, pw;
    mpz_init(gv); mpz_init(pw);
    for (int i = 0; i < n && ok; i++) {
        domain[i] = c->hi[i] - c->lo[i] + 1;
        total_tab += (long double)domain[i];
        if (total_tab > 100000000.0L) { ok = false; break; }   /* tables too big */
        gtab[i] = (int64_t*)malloc(sizeof(int64_t) * (size_t)domain[i]);
        for (int64_t d = 0; d < domain[i] && ok; d++) {
            int64_t v = c->lo[i] + d;
            mpz_set_ui(gv, 0);
            for (size_t t = 0; t < eq->n_terms; t++) {
                const int* ex = eq->exps + t * (size_t)n;
                if (ex[i] <= 0) continue;
                bool only_i = true;
                for (int u = 0; u < n; u++) if (u != i && ex[u] > 0) { only_i = false; break; }
                if (!only_i) continue;
                mpz_set(pw, eq->coefs[t]);
                for (int k = 0; k < ex[i]; k++) mpz_mul_si(pw, pw, (long)v);
                mpz_add(gv, gv, pw);
            }
            if (!mpz_fits_slong_p(gv)) { ok = false; break; }
            gtab[i][d] = mpz_get_si(gv);
        }
    }
    mpz_clear(gv); mpz_clear(pw);
    if (!ok) { for (int i = 0; i < n; i++) free(gtab[i]); return false; }

    /* Balance variables into two groups by product of domains (greedy). */
    int order[SI_MAX_VARS];
    for (int i = 0; i < n; i++) order[i] = i;
    for (int i = 0; i + 1 < n; i++) {           /* sort desc by domain */
        int pk = i;
        for (int j = i + 1; j < n; j++) if (domain[order[j]] > domain[order[pk]]) pk = j;
        if (pk != i) { int tmp = order[i]; order[i] = order[pk]; order[pk] = tmp; }
    }
    int gA[SI_MAX_VARS], gB[SI_MAX_VARS], nA = 0, nB = 0;
    long double pA = 1.0L, pB = 1.0L;
    for (int k = 0; k < n; k++) {
        int i = order[k];
        if (pA <= pB) { gA[nA++] = i; pA *= (long double)domain[i]; }
        else          { gB[nB++] = i; pB *= (long double)domain[i]; }
    }
    /* Hash the smaller group, iterate the larger. */
    int *hg, *ig, nh, ni;
    if (pA <= pB) { hg = gA; nh = nA; ig = gB; ni = nB; }
    else          { hg = gB; nh = nB; ig = gA; ni = nA; }
    long double p_hash = (pA <= pB) ? pA : pB;
    long double p_iter = (pA <= pB) ? pB : pA;
    if (p_hash > (long double)MITM_HASH_CAP || p_iter > (long double)SI_MAX_NODES) {
        for (int i = 0; i < n; i++) free(gtab[i]);
        return false;
    }

    /* Build the hash-group table via an odometer. */
    size_t hcap = (size_t)(p_hash) + 1, hcnt = 0;
    MitmEntry* H = (MitmEntry*)malloc(sizeof(MitmEntry) * hcap);
    int idx[SI_MAX_VARS];
    for (int j = 0; j < nh; j++) idx[j] = 0;
    for (;;) {
        int64_t sum = 0;
        MitmEntry* e = &H[hcnt];
        for (int j = 0; j < n; j++) e->vals[j] = 0;
        for (int j = 0; j < nh; j++) {
            int vi = hg[j]; int64_t v = c->lo[vi] + idx[j];
            e->vals[vi] = v; sum += gtab[vi][idx[j]];
        }
        e->sum = sum; hcnt++;
        int j = 0;                               /* odometer increment */
        for (; j < nh; j++) { if (++idx[j] < domain[hg[j]]) break; idx[j] = 0; }
        if (j == nh) break;
    }
    qsort(H, hcnt, sizeof(MitmEntry), mitm_cmp);

    /* Iterate the other group; binary-search complements. */
    int64_t full[SI_MAX_VARS];
    for (int j = 0; j < ni; j++) idx[j] = 0;
    if (ni == 0) {
        /* No iterate group: TARGET must be hit by hash entries alone. */
        for (size_t h = 0; h < hcnt; h++)
            if (H[h].sum == target && si_verify(c, H[h].vals)) emit_full(st, H[h].vals);
    } else for (;;) {
        int64_t sum = 0;
        for (int i = 0; i < n; i++) full[i] = 0;
        for (int j = 0; j < ni; j++) {
            int vi = ig[j]; int64_t v = c->lo[vi] + idx[j];
            full[vi] = v; sum += gtab[vi][idx[j]];
        }
        int64_t need = target - sum;
        /* lower_bound on need */
        size_t lo = 0, hi = hcnt;
        while (lo < hi) { size_t mid = (lo + hi) / 2; if (H[mid].sum < need) lo = mid + 1; else hi = mid; }
        for (size_t h = lo; h < hcnt && H[h].sum == need; h++) {
            for (int j = 0; j < nh; j++) full[hg[j]] = H[h].vals[hg[j]];
            if (si_verify(c, full)) emit_full(st, full);
        }
        int j = 0;
        for (; j < ni; j++) { if (++idx[j] < domain[ig[j]]) break; idx[j] = 0; }
        if (j == ni) break;
    }

    free(H);
    for (int i = 0; i < n; i++) free(gtab[i]);
    return true;
}

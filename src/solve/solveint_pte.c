/*
 * solveint_pte.c
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


/* --- Prouhet-Tarry-Escott  ->  {}  via Newton's identities. ---
 *
 * If two disjoint k-element variable groups G1, G2 satisfy the equal power sums
 *   Sum_{G1} v^j == Sum_{G2} v^j   for every j = 1..k,
 * then (Newton's identities) their elementary symmetric functions agree, so the
 * two groups are the SAME multiset.  With a strict ordering inside each group
 * this forces  G1[i] == G2[i]  elementwise, so a disequation like  a != d
 * between the corresponding smallest elements is contradictory: the system has
 * NO solution.  This settles the size-3, degree-1..3 case
 *   a+b+c==d+e+f && a^2+..==d^2+.. && a^3+..==d^3+.. && 0<a<b<c && 0<d<e<f && a!=d
 * as {} despite every variable being unbounded.  Returns the empty List when
 * proved empty, or NULL to decline (no forcing disequation -> parametric). */
static bool si_pset_eq(const bool* x, const bool* y, int n) {
    for (int i = 0; i < n; i++) if (x[i] != y[i]) return false;
    return true;
}

Expr* si_solve_power_sum_equal(SICtx* c) {
    int k = c->neq, n = c->n;
    if (k < 2 || n != 2 * k || k > SI_MAX_VARS) return NULL;

    bool Pset[SI_MAX_VARS], Mset[SI_MAX_VARS];
    bool have_sets = false;
    bool jseen[SI_MAX_VARS + 1];
    for (int j = 0; j <= k; j++) jseen[j] = false;

    for (int q = 0; q < k; q++) {
        const MPoly* eq = c->eq[q];
        bool p[SI_MAX_VARS], m[SI_MAX_VARS];
        for (int v = 0; v < n; v++) { p[v] = false; m[v] = false; }
        int np = 0, nm = 0, thisj = -1; bool ok = true;
        for (size_t t = 0; t < eq->n_terms && ok; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            int nv = -1, e = 0;
            for (int v = 0; v < n; v++) if (ex[v] > 0) { if (nv >= 0) { ok = false; break; } nv = v; e = ex[v]; }
            if (!ok) break;
            if (nv < 0) { ok = false; break; }             /* constant term: not PTE */
            if (thisj < 0) thisj = e; else if (e != thisj) { ok = false; break; }
            long cf = mpz_cmp_si(eq->coefs[t], 1) == 0 ? 1
                    : mpz_cmp_si(eq->coefs[t], -1) == 0 ? -1 : 0;
            if (cf == 1) { p[nv] = true; np++; }
            else if (cf == -1) { m[nv] = true; nm++; }
            else { ok = false; break; }
        }
        if (!ok || np != k || nm != k || thisj < 1 || thisj > k || jseen[thisj]) return NULL;
        jseen[thisj] = true;
        if (!have_sets) { for (int v = 0; v < n; v++) { Pset[v] = p[v]; Mset[v] = m[v]; } have_sets = true; }
        else if (!((si_pset_eq(p, Pset, n) && si_pset_eq(m, Mset, n))     /* same convention */
                 || (si_pset_eq(p, Mset, n) && si_pset_eq(m, Pset, n))))  /* or the sign flip */
            return NULL;
    }
    for (int j = 1; j <= k; j++) if (!jseen[j]) return NULL;              /* need degrees 1..k */

    int G1[SI_MAX_VARS], G2[SI_MAX_VARS], n1 = 0, n2 = 0;
    for (int v = 0; v < n; v++) { if (Pset[v]) G1[n1++] = v; if (Mset[v]) G2[n2++] = v; }
    if (n1 != k || n2 != k) return NULL;
    int ord1[SI_MAX_VARS], ord2[SI_MAX_VARS];
    if (!si_build_total_order(c, G1, k, ord1) || !si_build_total_order(c, G2, k, ord2)) return NULL;

    /* Newton forces ord1[i] == ord2[i]; a disequation between them is a proof of
     * emptiness. */
    for (int d = 0; d < c->n_neq; d++) {
        int a = c->neq_a[d], b = c->neq_b[d];
        for (int i = 0; i < k; i++)
            if ((a == ord1[i] && b == ord2[i]) || (a == ord2[i] && b == ord1[i]))
                return mk_list(NULL, 0);                                  /* {} */
    }
    return NULL;                                     /* no forcing disequation: decline */
}

/*
 * product_telescoping.c -- Product`Telescoping (Stage 1).
 *
 * Catches rational products whose anti-quotient is itself rational (no
 * Pochhammer/Gamma needed), giving the cleanest closed forms:
 *   Product[1 + 1/k, {k,1,n}]   -> 1 + n
 *   Product[k/(k+1), {k,1,n}]   -> 1/(1 + n)
 *   Product[1 - 1/k^2, {k,2,n}] -> (1 + n)/(2 n)
 *
 * Algorithm.  Write the factored rational as f(k) = prod_i (k - r_i)^{e_i}
 * (leading constant must be 1; numerator roots have e>0, denominator e<0).
 * We seek g(k) = prod (k - s)^{D_s} with g(k+1)/g(k) = f(k).  A generator
 * (k - s)^{D} contributes +D at root (s-1) and -D at root s, so the exponent of
 * the target at lattice position m satisfies E[m] = D[m+1] - D[m].  Per
 * integer-spaced chain this difference equation solves by prefix sums
 *   D[m] = sum_{j < m} E[j]
 * and closes (D returns to 0 outside the chain) iff sum_m E[m] == 0.  Then
 *   definite:   Product = Cancel[ g(imax+1) / g(imin) ]
 *   indefinite: Product[f, i] = g(i)   (up to the usual free constant).
 *
 * Returns NULL (fall through to Product`Rational) when f is non-rational, has a
 * symbolic exponent, has an irreducible/non-linear factor, has a leading
 * constant != 1, has a non-rational root, a chain that does not close, or a
 * chain span beyond CHAIN_SPAN_CAP (Product`Rational gives a compact Pochhammer
 * form in that case).
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

#define CHAIN_SPAN_CAP 1024

/* A merged (rational-root, integer-exponent) pair. */
typedef struct { Expr* val; int exp; } RootExp;

static bool is_num_rat(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational;
}

/* If evaluate(a - b) is a machine integer, store it and return true. */
static bool int_offset(const Expr* a, const Expr* b, int64_t* out) {
    Expr* negb = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)b) }, 2);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy((Expr*)a), negb }, 2);
    Expr* d = evaluate(sum);
    expr_free(sum);
    bool ok = (d->type == EXPR_INTEGER);
    if (ok) *out = d->data.integer;
    expr_free(d);
    return ok;
}

/* g-factor: (var - genval)^D, evaluated lazily by the caller's Times. */
static Expr* gfactor(Expr* var, const Expr* genval, int D) {
    Expr* negv = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)genval) }, 2);
    Expr* base = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(var), negv }, 2);
    Expr* be = evaluate(base);
    expr_free(base);
    Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ be, expr_new_integer(D) }, 2);  /* adopts be */
    return p;
}

static void free_roots(Expr** roots, int* mults, size_t n) {
    for (size_t i = 0; i < n; i++) expr_free(roots[i]);
    free(roots);
    free(mults);
}

/* Build g(var) for the whole rational from its merged root/exponent list.
 * Returns g (owned) on success, or NULL if any chain fails to close / is too
 * wide / has a non-rational root.  Appends nothing on failure. */
static Expr* build_g(RootExp* R, size_t nr, Expr* var) {
    /* All roots must be rational numbers for a rational anti-quotient. */
    for (size_t i = 0; i < nr; i++)
        if (!is_num_rat(R[i].val)) return NULL;

    /* Partition into integer-spaced chains via a union into representative
     * buckets.  chain[i] = index of the first root sharing i's chain. */
    size_t* chain = malloc(sizeof(size_t) * (nr ? nr : 1));
    for (size_t i = 0; i < nr; i++) {
        chain[i] = i;
        for (size_t j = 0; j < i; j++) {
            int64_t off;
            if (chain[j] == j && int_offset(R[i].val, R[j].val, &off)) {
                chain[i] = j; break;
            }
        }
    }

    size_t gcap = 8, gn = 0;
    Expr** gfac = malloc(sizeof(Expr*) * gcap);
    bool ok = true;

    for (size_t c = 0; c < nr && ok; c++) {
        if (chain[c] != c) continue;            /* not a chain representative */

        /* Gather this chain's members with lattice offsets from R[c]. */
        int64_t min_off = 0, max_off = 0;
        bool first = true;
        for (size_t i = 0; i < nr; i++) {
            if (chain[i] != c) continue;
            int64_t off = 0;
            if (i != c) {
                if (!int_offset(R[i].val, R[c].val, &off)) { ok = false; break; }
            }
            if (first) { min_off = max_off = off; first = false; }
            else { if (off < min_off) min_off = off; if (off > max_off) max_off = off; }
        }
        if (!ok) break;
        int64_t span = max_off - min_off;
        if (span < 0 || span > CHAIN_SPAN_CAP) { ok = false; break; }

        size_t L = (size_t)span + 1;
        int* E = calloc(L, sizeof(int));
        for (size_t i = 0; i < nr; i++) {
            if (chain[i] != c) continue;
            int64_t off = 0;
            if (i != c) (void)int_offset(R[i].val, R[c].val, &off);
            E[(size_t)(off - min_off)] += R[i].exp;
        }
        int total = 0;
        for (size_t m = 0; m < L; m++) total += E[m];
        if (total != 0) { free(E); ok = false; break; }

        /* D[m] = sum_{j<m} E[j]; generator at lattice m has value v0 + m. */
        int D = 0;
        for (size_t m = 0; m < L; m++) {
            /* D here is D[m] (prefix sum strictly below m). */
            if (D != 0) {
                /* genval = R[c].val + (min_off + m) */
                Expr* shift = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(R[c].val),
                               expr_new_integer(min_off + (int64_t)m) }, 2);
                Expr* gv = evaluate(shift);
                expr_free(shift);
                if (gn == gcap) { gcap *= 2; gfac = realloc(gfac, sizeof(Expr*) * gcap); }
                gfac[gn++] = gfactor(var, gv, D);
                expr_free(gv);
            }
            D += E[m];
        }
        free(E);
    }

    free(chain);
    if (!ok) {
        for (size_t i = 0; i < gn; i++) expr_free(gfac[i]);
        free(gfac);
        return NULL;
    }
    if (gn == 0) { free(gfac); return expr_new_integer(1); }
    Expr* times = expr_new_function(expr_new_symbol(SYM_Times), gfac, gn);
    free(gfac);
    Expr* g = evaluate(times);
    expr_free(times);
    return g;
}

Expr* builtin_product_telescoping(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol.name == SYM_Infinity)
        return NULL;
    if (prod_has_symbolic_power(f, var)) return NULL;

    Expr* tg  = prod_eval("Together",    (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = prod_eval("Numerator",   (Expr*[]){ expr_copy(tg) }, 1);
    Expr* den = prod_eval("Denominator", (Expr*[]){ expr_copy(tg) }, 1);
    expr_free(tg);

    Expr *lnum = NULL, *lden = NULL;
    Expr **rn = NULL, **rd = NULL;
    int *mn = NULL, *md = NULL;
    size_t nn = 0, nd = 0;
    bool aln = false, ald = false;
    bool okn = prod_linear_factors(num, var, &lnum, &rn, &mn, &nn, &aln);
    bool okd = okn && prod_linear_factors(den, var, &lden, &rd, &md, &nd, &ald);
    expr_free(num); expr_free(den);
    if (!okn || !okd || !aln || !ald) {
        if (okn) { free_roots(rn, mn, nn); expr_free(lnum); }
        if (okd) { free_roots(rd, md, nd); expr_free(lden); }
        return NULL;
    }

    /* Leading constant must be exactly 1 for a pure rational anti-quotient. */
    Expr* c = prod_div(lnum, lden);
    bool c_is_one = (c->type == EXPR_INTEGER && c->data.integer == 1);
    expr_free(c); expr_free(lnum); expr_free(lden);
    if (!c_is_one) { free_roots(rn, mn, nn); free_roots(rd, md, nd); return NULL; }

    /* Merge numerator (+mult) and denominator (-mult) roots by value. */
    size_t cap = nn + nd, rcount = 0;
    RootExp* R = malloc(sizeof(RootExp) * (cap ? cap : 1));
    for (size_t i = 0; i < nn; i++) { R[rcount].val = rn[i]; R[rcount].exp = mn[i]; rcount++; }
    for (size_t i = 0; i < nd; i++) { R[rcount].val = rd[i]; R[rcount].exp = -md[i]; rcount++; }
    /* (rn/rd Expr* ownership transferred into R; free only the arrays.) */
    free(rn); free(mn); free(rd); free(md);

    /* Coalesce equal-valued roots, summing exponents. */
    for (size_t i = 0; i < rcount; i++) {
        if (!R[i].val) continue;
        for (size_t j = i + 1; j < rcount; j++) {
            if (R[j].val && expr_eq(R[i].val, R[j].val)) {
                R[i].exp += R[j].exp;
                expr_free(R[j].val); R[j].val = NULL;
            }
        }
    }
    /* Compact, dropping zero-exponent (fully cancelled) roots. */
    size_t w = 0;
    for (size_t i = 0; i < rcount; i++) {
        if (R[i].val && R[i].exp != 0) R[w++] = R[i];
        else if (R[i].val) expr_free(R[i].val);
    }
    rcount = w;

    Expr* g = build_g(R, rcount, var);
    for (size_t i = 0; i < rcount; i++) expr_free(R[i].val);
    free(R);
    if (!g) return NULL;

    if (!definite) return g;   /* indefinite anti-quotient g(i) */

    /* Definite: Cancel[ g(imax+1) / g(imin) ]. */
    Expr* up = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy(imax), expr_new_integer(1) }, 2);
    Expr* upe = evaluate(up);
    expr_free(up);
    Expr* ghi = prod_subst(g, var, upe);
    Expr* glo = prod_subst(g, var, imin);
    expr_free(upe); expr_free(g);
    Expr* ratio = prod_div(ghi, glo);
    expr_free(ghi); expr_free(glo);
    Expr* out = prod_eval("Cancel", (Expr*[]){ ratio }, 1);  /* adopts ratio */
    return out;
}

void product_telescoping_init(void) {
    symtab_add_builtin("Product`Telescoping", builtin_product_telescoping);
    symtab_get_def("Product`Telescoping")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Telescoping",
        "Product`Telescoping[f, i, imin, imax] gives the closed form of a "
        "rational product whose anti-quotient is itself rational (Gamma-free), "
        "via integer-spaced root chains. Returns unevaluated to fall through "
        "when the anti-quotient is not rational.");
}

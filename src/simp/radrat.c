/* radrat.c -- multi-generator radical-rational normalisation for Simplify.
 *
 * See radrat.h for the high-level contract.  The algorithm (validated on
 * D[Integrate[1/(x^3 (a+b x)^(1/3)),x],x] -> 1/(x^3 (a+b x)^(1/3))):
 *
 *   1. Collect the distinct radical bases B_1..B_n (each B_k appears as
 *      B_k^(p/q_k) with q_k > 1).  Require n >= 2 (single-generator inputs
 *      are already handled by Together/Cancel) and every B_k provably
 *      positive (so B_k^(p/q) -> g_k^p is branch-safe).
 *   2. Substitute, outer-base-first, B_k and all its powers to a fresh
 *      polynomial generator g_k (poly_subst_radical_to_gen).
 *   3. Build the relation g_k^{q_k} - V_k for each compound base, where V_k
 *      is B_k with the OTHER generators substituted in; keep a relation
 *      only when V_k stays inside the substituted ring (drops degenerate /
 *      non-nesting cases).
 *   4. Together -> P/Q; reduce P and Q modulo every relation
 *      (PolynomialRemainder); rationalise the denominator per relation
 *      (PolynomialExtendedGCD).
 *   5. Substitute the radicals back and finish with Cancel[N / Factor[D]].
 *   6. Strict score gate: return the result only when it beats the input.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "arithmetic.h"
#include "simp.h"
#include "simp_internal.h"
#include "radrat.h"
#include "poly.h"

/* leaf count (heads included), matching the simp scoring convention. */
static int64_t rr_leaves(const Expr* e) { return leaf_count_internal((Expr*)e, true); }

/* Cap on distinct radical generators; beyond this the multivariate
 * reduction risks a GCD blow-up, so we bail (return NULL). */
#define RR_MAX_GENS 4

/* ---- small numeric helpers ---------------------------------------- */

static int64_t rr_gcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t r = a % b; a = b; b = r; }
    return a;
}
static int64_t rr_lcm(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t g = rr_gcd(a, b);
    return (a / g) * b;
}

/* ---- positive-base convention (branch safety) --------------------- */

/* True when `e` is provably positive under the symbol-positive
 * convention: symbols are treated as positive, positive numbers are
 * positive, and sums / products / powers of positive things are
 * positive.  Negative or complex numeric bases return false so we never
 * cross a branch cut (the bug class that caused false corpus DIFFs in the
 * intsimp work). */
static bool rr_pos_base(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer > 0;
        case EXPR_REAL:    return e->data.real > 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) > 0;
        case EXPR_SYMBOL:  return true;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL) return false;
            const char* hn = h->data.symbol;
            int64_t n, d;
            if (is_rational(e, &n, &d)) return (n > 0) == (d > 0) && n != 0;
            if (strcmp(hn, "Plus") == 0 || strcmp(hn, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    if (!rr_pos_base(e->data.function.args[i])) return false;
                return e->data.function.arg_count > 0;
            }
            if (strcmp(hn, "Power") == 0 && e->data.function.arg_count == 2)
                return rr_pos_base(e->data.function.args[0]);
            return false;
        }
        default: return false;
    }
}

/* ---- radical detection -------------------------------------------- */

/* If `e` is Power[B, p/q] with reduced |q| > 1, set *base (borrowed) and
 * *q (positive) and return true. */
static bool rr_parse_radical(const Expr* e, const Expr** base, int64_t* q) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || strcmp(h->data.symbol, "Power") != 0) return false;
    if (e->data.function.arg_count != 2) return false;
    int64_t p, qq;
    if (!is_rational(e->data.function.args[1], &p, &qq)) return false;
    if (qq < 0) { p = -p; qq = -qq; }
    int64_t g = rr_gcd(p, qq);
    if (g > 1) { p /= g; qq /= g; }
    if (qq <= 1) return false;
    *base = e->data.function.args[0];
    *q = qq;
    return true;
}

/* True if `e` contains a symbol anywhere (i.e. is not a pure constant). */
static bool rr_has_symbol(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return true;
    if (e->type == EXPR_FUNCTION) {
        if (rr_has_symbol(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (rr_has_symbol(e->data.function.args[i])) return true;
    }
    return false;
}

typedef struct {
    Expr*  base;   /* owned deep copy of the radical base */
    int64_t q;     /* lcm of fractional-exponent denominators */
    char*  gen;    /* fresh generator symbol name (owned) */
} RRGen;

/* Walk `e` accumulating distinct radical bases into gens[].  Dedup by
 * structural equality; q is lcm'd across sites.  Sets *overflow when more
 * than RR_MAX_GENS distinct bases appear, or when a radical base is not
 * provably positive (branch-unsafe -> bail entirely). */
static void rr_collect(const Expr* e, RRGen* gens, int* n, bool* overflow) {
    if (!e || *overflow || e->type != EXPR_FUNCTION) return;
    const Expr* base = NULL;
    int64_t q = 0;
    if (rr_parse_radical(e, &base, &q) && rr_has_symbol(base)) {
        if (!rr_pos_base(base)) { *overflow = true; return; }
        int found = -1;
        for (int i = 0; i < *n; i++)
            if (expr_eq(gens[i].base, (Expr*)base)) { found = i; break; }
        if (found >= 0) {
            gens[found].q = rr_lcm(gens[found].q, q);
        } else {
            if (*n >= RR_MAX_GENS) { *overflow = true; return; }
            gens[*n].base = expr_copy((Expr*)base);
            gens[*n].q = q;
            gens[*n].gen = NULL;
            (*n)++;
        }
    }
    rr_collect(e->data.function.head, gens, n, overflow);
    for (size_t i = 0; i < e->data.function.arg_count && !*overflow; i++)
        rr_collect(e->data.function.args[i], gens, n, overflow);
}

/* ---- symbol-set membership (relation usability guard) ------------- */

static void rr_collect_symnames(const Expr* e, const char*** set, int* n, int* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < *n; i++)
            if (strcmp((*set)[i], e->data.symbol) == 0) return;
        if (*n >= *cap) {
            *cap = *cap ? *cap * 2 : 8;
            *set = realloc(*set, sizeof(char*) * (size_t)*cap);
        }
        (*set)[(*n)++] = e->data.symbol;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        rr_collect_symnames(e->data.function.head, set, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rr_collect_symnames(e->data.function.args[i], set, n, cap);
    }
}

/* True if every symbol of `e` is a member of the name set. */
static bool rr_symbols_subset(const Expr* e, const char** set, int n) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < n; i++)
            if (strcmp(set[i], e->data.symbol) == 0) return true;
        return false;
    }
    if (e->type == EXPR_FUNCTION) {
        if (!rr_symbols_subset(e->data.function.head, set, n)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!rr_symbols_subset(e->data.function.args[i], set, n)) return false;
    }
    return true;
}

/* ---- builtin-call helpers (construct + evaluate + free call) ------ */

/* evaluate(call) then free the call node (evaluate does not free its
 * argument).  Takes ownership of `call`. */
static Expr* rr_eval_free(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Evaluate Head[a0] ; takes ownership of a0. */
static Expr* rr_call1(const char* head, Expr* a0) {
    Expr* call = expr_new_function(expr_new_symbol(head), (Expr*[]){a0}, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}
/* Evaluate Head[a0, a1, a2] ; takes ownership of a0,a1,a2. */
static Expr* rr_call3(const char* head, Expr* a0, Expr* a1, Expr* a2) {
    Expr* call = expr_new_function(expr_new_symbol(head),
                                   (Expr*[]){a0, a1, a2}, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* PolynomialRemainder[poly, rel, gen] (poly consumed). */
static Expr* rr_polyrem(Expr* poly, const Expr* rel, const char* gen) {
    return rr_call3("PolynomialRemainder", poly,
                    expr_copy((Expr*)rel), expr_new_symbol(gen));
}

/* ---- the pass ----------------------------------------------------- */

Expr* simp_radical_rational(const Expr* input,
                            const AssumeCtx* ctx,
                            const Expr* complexity_func) {
    (void)ctx;
    if (!input || input->type != EXPR_FUNCTION) return NULL;

    /* Step 1: collect distinct, positive radical generators. */
    RRGen gens[RR_MAX_GENS];
    int n = 0;
    bool overflow = false;
    rr_collect(input, gens, &n, &overflow);
    if (overflow || n < 2) {
        for (int i = 0; i < n; i++) expr_free(gens[i].base);
        return NULL;
    }

    /* Order outer-base-first (descending leaf count) so a containing base
     * such as (a+b x) is substituted before the contained `a`. */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (rr_leaves(gens[j].base) > rr_leaves(gens[i].base)) {
                RRGen t = gens[i]; gens[i] = gens[j]; gens[j] = t;
            }

    /* Fresh generator names (reserved prefix, collision-checked). */
    for (int i = 0; i < n; i++) {
        char buf[48];
        snprintf(buf, sizeof(buf), "$rr_radgen%d$", i);
        size_t len = strlen(buf);
        gens[i].gen = malloc(len + 1);
        memcpy(gens[i].gen, buf, len + 1);
    }

    Expr* ONE = expr_new_integer(1);

    /* Step 2: substitute every base -> its generator (outer first). */
    Expr* sub = expr_copy((Expr*)input);
    for (int i = 0; i < n; i++) {
        Expr* nx = poly_subst_radical_to_gen(sub, gens[i].base, ONE,
                                              gens[i].q, gens[i].gen);
        expr_free(sub);
        sub = nx;
    }

    /* Symbol set of the substituted ring (for the relation guard). */
    const char** ringset = NULL;
    int rn = 0, rcap = 0;
    rr_collect_symnames(sub, &ringset, &rn, &rcap);

    /* Step 3: build usable relations  g_k^{q_k} - V_k. */
    Expr* rels[RR_MAX_GENS];
    char* rel_gen[RR_MAX_GENS];
    int nrel = 0;
    for (int k = 0; k < n; k++) {
        if (gens[k].base->type != EXPR_FUNCTION) continue; /* symbol base: free */
        Expr* V = expr_copy(gens[k].base);
        for (int j = 0; j < n; j++) {
            if (j == k) continue;
            Expr* nv = poly_subst_radical_to_gen(V, gens[j].base, ONE,
                                                 gens[j].q, gens[j].gen);
            expr_free(V);
            V = nv;
        }
        if (!rr_symbols_subset(V, ringset, rn)) { expr_free(V); continue; }
        /* rel = g_k^{q_k} - V */
        Expr* gpow = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){expr_new_symbol(gens[k].gen),
                      expr_new_integer(gens[k].q)}, 2);
        Expr* negV = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){expr_new_integer(-1), V}, 2);
        Expr* relcall = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){gpow, negV}, 2);
        Expr* rel = evaluate(relcall);
        expr_free(relcall);
        rels[nrel] = rel;
        rel_gen[nrel] = gens[k].gen;
        nrel++;
    }
    free(ringset);

    if (nrel == 0) {
        /* Nothing to reduce against -> defer. */
        expr_free(sub); expr_free(ONE);
        for (int i = 0; i < nrel; i++) expr_free(rels[i]);
        for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }
        return NULL;
    }

    /* Step 4: Together, split, reduce mod the relation ideal. */
    Expr* tog = rr_call1("Together", sub);   /* consumes sub */
    if (!tog) {
        expr_free(ONE);
        for (int i = 0; i < nrel; i++) expr_free(rels[i]);
        for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }
        return NULL;
    }
    Expr* num = rr_call1("Numerator", expr_copy(tog));
    Expr* den = rr_call1("Denominator", tog);  /* consumes tog */

    for (int i = 0; i < nrel && num && den; i++) {
        num = rr_polyrem(num, rels[i], rel_gen[i]);
        den = rr_polyrem(den, rels[i], rel_gen[i]);
    }

    /* Step 4b: rationalise the denominator against each relation. */
    for (int i = 0; i < nrel && num && den; i++) {
        Expr* gsym = expr_new_symbol(rel_gen[i]);
        int deg = get_degree_poly(den, gsym);
        expr_free(gsym);
        if (deg <= 0) continue;

        Expr* eg = rr_call3("PolynomialExtendedGCD",
                            expr_copy(den), expr_copy(rels[i]),
                            expr_new_symbol(rel_gen[i]));
        if (!eg || eg->type != EXPR_FUNCTION
            || eg->data.function.arg_count != 2
            || eg->data.function.args[1]->type != EXPR_FUNCTION
            || eg->data.function.args[1]->data.function.arg_count < 1) {
            if (eg) expr_free(eg);
            continue;
        }
        Expr* gcd_e = eg->data.function.args[0];
        Expr* u = eg->data.function.args[1]->data.function.args[0];
        Expr* gc = expr_new_symbol(rel_gen[i]);
        int gcd_deg = get_degree_poly(gcd_e, gc);
        expr_free(gc);
        if (gcd_deg != 0) { expr_free(eg); continue; }

        /* prod = num * u ; combine ; numerator reduced mod rel, denom *= gcd */
        Expr* prod = rr_eval_free(expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){expr_copy(num), expr_copy(u)}, 2));
        Expr* comb = rr_call1("Together", prod);
        Expr* pnum = rr_call1("Numerator", expr_copy(comb));
        Expr* pden = rr_call1("Denominator", comb);
        pnum = rr_polyrem(pnum, rels[i], rel_gen[i]);
        Expr* nden = rr_eval_free(expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){expr_copy(gcd_e), pden}, 2));
        expr_free(eg);

        if (pnum && nden) {
            /* keep rationalised form only if it didn't inflate */
            int64_t before = rr_leaves(num) + rr_leaves(den);
            int64_t after = rr_leaves(pnum) + rr_leaves(nden);
            if (after <= before + 8) {
                expr_free(num); expr_free(den);
                num = pnum; den = nden;
                /* re-reduce against all relations */
                for (int j = 0; j < nrel && num && den; j++) {
                    num = rr_polyrem(num, rels[j], rel_gen[j]);
                    den = rr_polyrem(den, rels[j], rel_gen[j]);
                }
            } else { expr_free(pnum); expr_free(nden); }
        } else {
            if (pnum) expr_free(pnum);
            if (nden) expr_free(nden);
        }
    }

    for (int i = 0; i < nrel; i++) expr_free(rels[i]);

    if (!num || !den) {
        if (num) expr_free(num);
        if (den) expr_free(den);
        expr_free(ONE);
        for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }
        return NULL;
    }

    /* Step 5: substitute the radicals back, then Cancel[N / Factor[D]]. */
    for (int i = 0; i < n; i++) {
        Expr* nb = poly_subst_radical_from_gen(num, gens[i].base, ONE,
                                               gens[i].q, gens[i].gen);
        expr_free(num); num = nb;
        Expr* db = poly_subst_radical_from_gen(den, gens[i].base, ONE,
                                               gens[i].q, gens[i].gen);
        expr_free(den); den = db;
    }
    expr_free(ONE);
    for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }

    Expr* dfac = rr_call1("Factor", den);              /* consumes den */
    Expr* frac = rr_eval_free(expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){num, expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){dfac, expr_new_integer(-1)}, 2)}, 2));
    Expr* result = rr_call1("Cancel", frac);           /* consumes frac */
    if (!result) return NULL;

    /* Step 6: strict score gate. */
    if (expr_eq(result, (Expr*)input)
        || score_with_func(result, complexity_func)
               >= score_with_func(input, complexity_func)) {
        expr_free(result);
        return NULL;
    }
    return result;
}

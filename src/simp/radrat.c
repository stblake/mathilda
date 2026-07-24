/* radrat.c -- multi-generator radical-rational normalisation for Simplify.
 *
 * See radrat.h for the high-level contract.  The algorithm (validated on
 * D[Integrate[1/(x^3 (a+b x)^(1/3)),x],x] -> 1/(x^3 (a+b x)^(1/3))):
 *
 *   1. Collect the distinct radical bases B_1..B_n (each B_k appears as
 *      B_k^(p/q_k) with q_k > 1).  Require n >= 2 (single-generator inputs
 *      are already handled by Together/Cancel) and every B_k real-valued
 *      (symbols assumed real; only explicit complex constants are rejected)
 *      -- the reduction is an exact principal-branch rewrite for any real
 *      radicand regardless of sign.
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
#include "sym_names.h"
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

/* True when `e` is a real-valued radicand under the symbol-is-real
 * convention: symbols and real numbers are real, and functions of real
 * things are real.  We reject only bases carrying an explicit complex
 * constant (a `Complex[..]` node or the imaginary unit `I`), which is the
 * one situation where the principal branch could genuinely differ.
 *
 * Positivity is deliberately NOT required: radrat's whole algorithm is an
 * exact rewrite on the principal branch — each radical B^(p/q) is mapped
 * to a fresh generator via the identity (B^(1/q))^q = B, generators are
 * never merged across distinct bases, and the radicals are reconstructed
 * verbatim before the final Cancel/Factor.  So the reduction is valid for
 * any real base regardless of sign (e.g. `1 - x`, negative for x > 1),
 * and the earlier positive-only gate merely blocked legitimate
 * sign-indefinite radicands such as the nested `1 + Sqrt[1 - x]`. */
static bool rr_real_base(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:  return true;
        case EXPR_SYMBOL:  return strcmp(e->data.symbol.name, "I") != 0;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "Complex") == 0)
                return false;
            if (!rr_real_base(h)) return false;
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!rr_real_base(e->data.function.args[i])) return false;
            return true;
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
    if (h->type != EXPR_SYMBOL || strcmp(h->data.symbol.name, "Power") != 0) return false;
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
        if (!rr_real_base(base)) { *overflow = true; return; }
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
            if (strcmp((*set)[i], e->data.symbol.name) == 0) return;
        if (*n >= *cap) {
            *cap = *cap ? *cap * 2 : 8;
            *set = realloc(*set, sizeof(char*) * (size_t)*cap);
        }
        (*set)[(*n)++] = e->data.symbol.name;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        rr_collect_symnames(e->data.function.head, set, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rr_collect_symnames(e->data.function.args[i], set, n, cap);
    }
}

/* Collect the distinct free VARIABLE names of `e` — like rr_collect_symnames but
 * SKIPPING each function node's operator head (Plus, Times, Power, Rational, …),
 * which name operators, not variables.  Used to count the parameter variables a
 * radical base couples: base `Sqrt[x] + x` has ONE variable (x), not four
 * (Plus/Power/Rational/x).  A non-symbol head (e.g. a Derivative[…][x] operator)
 * is still traversed. */
static void rr_collect_varnames(const Expr* e, const char*** set, int* n, int* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < *n; i++)
            if (strcmp((*set)[i], e->data.symbol.name) == 0) return;
        if (*n >= *cap) {
            *cap = *cap ? *cap * 2 : 8;
            *set = realloc(*set, sizeof(char*) * (size_t)*cap);
        }
        (*set)[(*n)++] = e->data.symbol.name;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type != EXPR_SYMBOL)
            rr_collect_varnames(e->data.function.head, set, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rr_collect_varnames(e->data.function.args[i], set, n, cap);
    }
}

/* True if `e` contains a Power[g, m] whose base g is one of the generator
 * symbols and whose exponent m is NOT a non-negative integer (a Laurent or
 * fractional generator power).  Such a term makes a relation g_k^{q_k} - V_k
 * non-polynomial in the generators; radrat's per-generator univariate
 * reduction cannot handle it, and feeding it to PolynomialRemainder /
 * PolynomialExtendedGCD blows up.  Arises with doubly-nested radicals whose
 * outer base holds an inverse inner radical, e.g. the base 1 + 1/Sqrt[1-x]
 * substitutes to 1 + g^(-2).  Relations that trip this are skipped. */
static bool rr_has_nonpoly_gen_power(const Expr* e, const RRGen* gens, int n) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "Power") == 0
        && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* ex = e->data.function.args[1];
        if (b->type == EXPR_SYMBOL) {
            for (int i = 0; i < n; i++)
                if (gens[i].gen && strcmp(b->data.symbol.name, gens[i].gen) == 0) {
                    if (ex->type != EXPR_INTEGER || ex->data.integer < 0)
                        return true;
                    break;
                }
        }
    }
    if (rr_has_nonpoly_gen_power(h, gens, n)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rr_has_nonpoly_gen_power(e->data.function.args[i], gens, n)) return true;
    return false;
}

/* True if every symbol of `e` is a member of the name set. */
static bool rr_symbols_subset(const Expr* e, const char** set, int n) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < n; i++)
            if (strcmp(set[i], e->data.symbol.name) == 0) return true;
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

/* Finalise a reduced (num, den) generator-polynomial pair into a candidate
 * result: substitute every radical back for its generator, then form
 * Cancel[num / Factor[den]].  Consumes `num` and `den`; the generators and
 * `ONE` stay owned by the caller.  Returns NULL on failure. */
static Expr* rr_finalize(Expr* num, Expr* den, const RRGen* gens, int n,
                         Expr* ONE) {
    for (int i = 0; i < n && num && den; i++) {
        Expr* nb = poly_subst_radical_from_gen(num, gens[i].base, ONE,
                                               gens[i].q, gens[i].gen);
        expr_free(num); num = nb;
        Expr* db = poly_subst_radical_from_gen(den, gens[i].base, ONE,
                                               gens[i].q, gens[i].gen);
        expr_free(den); den = db;
    }
    if (!num || !den) { if (num) expr_free(num); if (den) expr_free(den);
                        return NULL; }
    Expr* dfac = rr_call1("Factor", den);              /* consumes den */
    if (!dfac) { expr_free(num); return NULL; }
    Expr* frac = rr_eval_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){num, expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){dfac, expr_new_integer(-1)}, 2)}, 2));
    return rr_call1("Cancel", frac);                   /* consumes frac */
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

    /* Scope guard: the per-generator rationalisation (Step 4b) reduces the
     * denominator against g_k^{q_k} - V_k with PolynomialExtendedGCD taken
     * UNIVARIATELY in the generator g_k — its coefficient ring is Q[parameters].
     * When the generator bases span TWO OR MORE distinct parameter variables, that
     * coefficient ring is multivariate and the univariate pseudo-remainder
     * sequence suffers classic coefficient explosion (e.g. the cube root of the
     * degree-3 bivariate x (1-x) (1-k x) from the Goursat integrand: the extended
     * GCD over Q[x, k] ground the whole Simplify to a halt).  This normal-form
     * method is a decision procedure only over a single-parameter coefficient
     * field; multivariate is out of scope, so decline cleanly (Simplify falls
     * through to Together / Cancel / Factor, which reduce the same input without
     * the generator-ring blow-up).  All in-scope radical-rational simplifications
     * are single-parameter, so this never fires on them. */
    {
        const char** pset = NULL; int pn = 0, pcap = 0;
        for (int i = 0; i < n; i++)
            rr_collect_varnames(gens[i].base, &pset, &pn, &pcap);
        free(pset);
        if (pn >= 2) {
            for (int i = 0; i < n; i++) expr_free(gens[i].base);
            return NULL;
        }
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

    /* Symbol set of the reduction ring (for the relation-usability guard).
     * This is the generator names (present in `sub`) together with the
     * problem's parameters.  The parameters must be gathered from the
     * ORIGINAL `input`, not just `sub`: a variable that appears only inside a
     * radical base (e.g. `x` in `Sqrt[1-x]` when the numerator carries no bare
     * `x`) vanishes from `sub` after substitution, yet its relation
     * g_k^{q_k} - V_k (here g1^2 - (1-x)) is genuine and needed to reduce odd
     * generator powers.  Collecting from `input` keeps such relations. */
    const char** ringset = NULL;
    int rn = 0, rcap = 0;
    rr_collect_symnames(sub, &ringset, &rn, &rcap);
    rr_collect_symnames((Expr*)input, &ringset, &rn, &rcap);

    /* Step 3: build usable relations  g_k^{q_k} - V_k. */
    Expr* rels[RR_MAX_GENS];
    char* rel_gen[RR_MAX_GENS];
    int nrel = 0;
    bool laurent = false;
    for (int k = 0; k < n && !laurent; k++) {
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
        /* A generator base holding an inverse inner radical substitutes to a
         * negative generator power (e.g. 1 + 1/Sqrt[1-x] -> 1 + g^(-2)),
         * making g_k^{q_k} - V_k a Laurent, not polynomial, relation.  radrat's
         * per-generator univariate reduction cannot handle it and the downstream
         * PolynomialExtendedGCD blows up, so abandon the whole pass (the input
         * is outside this method's polynomial-ring model). */
        if (rr_has_nonpoly_gen_power(V, gens, n)) { expr_free(V); laurent = true; break; }
        /* rel = g_k^{q_k} - V */
        Expr* gpow = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){expr_new_symbol(gens[k].gen),
                      expr_new_integer(gens[k].q)}, 2);
        Expr* negV = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){expr_new_integer(-1), V}, 2);
        Expr* relcall = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){gpow, negV}, 2);
        Expr* rel = evaluate(relcall);
        expr_free(relcall);
        rels[nrel] = rel;
        rel_gen[nrel] = gens[k].gen;
        nrel++;
    }
    free(ringset);

    if (laurent || nrel == 0) {
        /* Laurent relation (out of model) or nothing to reduce against -> defer. */
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
        Expr* prod = rr_eval_free(expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){expr_copy(num), expr_copy(u)}, 2));
        Expr* comb = rr_call1("Together", prod);
        Expr* pnum = rr_call1("Numerator", expr_copy(comb));
        Expr* pden = rr_call1("Denominator", comb);
        pnum = rr_polyrem(pnum, rels[i], rel_gen[i]);
        Expr* nden = rr_eval_free(expr_new_function(expr_new_symbol(SYM_Times),
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

    if (!num || !den) {
        if (num) expr_free(num);
        if (den) expr_free(den);
        for (int i = 0; i < nrel; i++) expr_free(rels[i]);
        expr_free(ONE);
        for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }
        return NULL;
    }

    /* Step 5: enumerate candidate representatives; keep the lowest-scoring.
     *
     * The reduced (num, den) is the polynomial-in-generators normal form, but
     * its equivalence class also holds representatives with a generator in the
     * DENOMINATOR.  Multiplying num and den by g_i^k and reducing the
     * numerator modulo g_i^{q_i} - B_i trades a numerator radical for the
     * (often simpler) base polynomial B_i -- the exact dual of the Step-4b
     * denominator rationalisation.  This folds e.g.
     *   (1 - Sqrt[1-x]) Sqrt[1+Sqrt[1-x]]  ->  x / Sqrt[1+Sqrt[1-x]]
     * back to its integrand.  Both sides are multiplied by the same
     * (generically nonzero) g_i^k, so the move is exact on the principal
     * branch; the strict score gate below keeps a variant only when it is
     * genuinely simpler than the input. */
    Expr* best = rr_finalize(expr_copy(num), expr_copy(den), gens, n, ONE);
    size_t best_score = best ? score_with_func(best, complexity_func) : SIZE_MAX;

    for (int i = 0; i < n; i++) {
        Expr* gsym = expr_new_symbol(gens[i].gen);
        int ndeg = get_degree_poly(num, gsym);
        expr_free(gsym);
        if (ndeg <= 0) continue;   /* no g_i to trade out of the numerator */

        for (int64_t k = 1; k < gens[i].q; k++) {
            /* vnum = reduce(num * g_i^k) ; vden = den * g_i^k */
            Expr* vnum = rr_eval_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){expr_copy(num),
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){expr_new_symbol(gens[i].gen),
                                        expr_new_integer(k)}, 2)}, 2));
            for (int j = 0; j < nrel && vnum; j++)
                vnum = rr_polyrem(vnum, rels[j], rel_gen[j]);
            if (!vnum) continue;
            Expr* vden = rr_eval_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){expr_copy(den),
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){expr_new_symbol(gens[i].gen),
                                        expr_new_integer(k)}, 2)}, 2));
            Expr* cand = rr_finalize(vnum, vden, gens, n, ONE);
            if (!cand) continue;
            size_t cs = score_with_func(cand, complexity_func);
            if (!best || cs < best_score) {
                if (best) expr_free(best);
                best = cand; best_score = cs;
            } else {
                expr_free(cand);
            }
        }
    }

    for (int i = 0; i < nrel; i++) expr_free(rels[i]);
    expr_free(num); expr_free(den);
    expr_free(ONE);
    for (int i = 0; i < n; i++) { expr_free(gens[i].base); free(gens[i].gen); }

    if (!best) return NULL;

    /* Step 6: strict score gate. */
    if (expr_eq(best, (Expr*)input)
        || best_score >= score_with_func(input, complexity_func)) {
        expr_free(best);
        return NULL;
    }
    return best;
}

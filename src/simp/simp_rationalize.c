#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* Denominator rationalisation: simp_rationalize_denom                     */
/* ----------------------------------------------------------------------- */

/*
 * simp_rationalize_denom rewrites Power[denom, -1] subtrees where
 * `denom` is a polynomial in radicals over a single integer base `c`,
 * by computing the inverse via extended-Euclidean in Q[α]/(α^n - c)
 * with α = c^(1/n) and n the LCM of the denominators of the rational
 * exponents that appear. Closes two phase-2 user cases:
 *
 *   Simplify[1/(2^(1/3) - 1)]              -> 4^(1/3) + 2^(1/3) + 1
 *   Simplify[1/(Sqrt[2] + 2^(1/3))]        -> ... in Q[2^(1/6)]/(α^6 - 2)
 *
 * Algorithm:
 *   1. Walk the tree top-down. At each Power[denom, k] with k a
 *      negative rational, attempt rationalisation.
 *   2. Find all Power[c, p_i/q_i] subexpressions of denom. Reject
 *      (return NULL) if they don't all share a single positive
 *      integer base c.
 *   3. Compute n = lcm(q_i). Substitute each Power[c, p_i/q_i] with
 *      Power[α, p_i * (n / q_i)] for a fresh variable α.
 *   4. Run PolynomialExtendedGCD on the substituted denom and α^n - c.
 *      Accept the result iff the gcd is a nonzero constant g; the
 *      first Bezout coefficient u(α) then satisfies
 *      u(α) * denom(α) ≡ g (mod α^n - c), so 1/denom = u(α)/g.
 *   5. Substitute α → Power[c, 1/n] and reduce. The result has no
 *      radicals in its denominator.
 *
 * Soundness: the rewrite is correctness-preserving for any positive
 * integer base c (so c^(1/n) is real-valued and the principal-branch
 * arithmetic is unambiguous). For c = 1 the radicals collapse to
 * constants and the substitution is trivial; we skip that case so
 * the round loop's other transforms can handle it.
 */

/* Forward declaration: simp_search drives this. */
Expr* simp_rationalize_denom(const Expr* e, const AssumeCtx* ctx);

/* Forward declaration: applied as a post-pass on the inverse so that
 * c^(p/n) forms with c = base^k canonicalise to a single power-of-base
 * representation. Without this post-pass, my candidate
 * "1 + 2^(1/3) + 2^(2/3)" subtracted from the user's
 * "1 + 2^(1/3) + 4^(1/3)" leaves 2^(2/3) - 4^(1/3) which is
 * mathematically zero but syntactically nonzero. PrimeRebase rewrites
 * 4^(1/3) -> 2^(2/3) and the cancellation completes. PrimeRebase
 * normally only runs at simp_dispatch on the input, not on candidates
 * produced inside simp_search, so we apply it explicitly here. */
Expr* transform_prime_rebase(const Expr* e);

/* GCD/LCM on int64_t (inputs assumed non-negative). */
static int64_t denom_gcd_i64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int64_t denom_lcm_i64(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t g = denom_gcd_i64(a, b);
    /* Overflow guard: skip when the product would exceed INT64_MAX. */
    if (a / g > INT64_MAX / b) return 0;
    return (a / g) * b;
}

/* Walk e collecting (base, q) pairs where each Power[base, p/q] uses
 * the SAME positive integer base. Populates *out_base (borrowed
 * pointer into e) and accumulates n = lcm of all q's into *out_n.
 * Returns false on:
 *   - any Power with non-rational exponent
 *   - any Power[c, p/q] with q == 1 (those are integer powers, not
 *     radicals — they pass through cleanly)
 *   - any Power whose base is not a positive integer
 *   - bases that disagree (multi-base extensions are out of scope here).
 *
 * Returns true with *out_base = NULL when e contains no radical
 * subexpressions at all (caller treats as "nothing to rationalise"). */
static bool denom_collect_radical_base(const Expr* e,
                                       const Expr** out_base,
                                       int64_t* out_n) {
    if (!e) return true;
    if (e->type != EXPR_FUNCTION) return true;

    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1) {
            /* Reject base != positive integer. */
            if (base->type != EXPR_INTEGER || base->data.integer < 2)
                return false;
            /* Reject base disagreement. */
            if (*out_base) {
                if (!expr_eq((Expr*)*out_base, (Expr*)base)) return false;
            } else {
                *out_base = base;
            }
            int64_t new_n = denom_lcm_i64(*out_n ? *out_n : 1, q);
            if (new_n == 0) return false;
            *out_n = new_n;
            /* Don't recurse into the base — it's an integer leaf. */
            return true;
        }
    }
    /* Recurse into children. */
    if (!denom_collect_radical_base(e->data.function.head, out_base, out_n))
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!denom_collect_radical_base(e->data.function.args[i],
                                        out_base, out_n))
            return false;
    }
    return true;
}

/* Substitute every Power[base, p/q] in e with Power[gen, p * (n/q)]
 * (an integer-exponent power of the fresh symbol). Returns a fresh
 * tree. */
static Expr* denom_subst_radical_to_gen(const Expr* e, const Expr* base,
                                        int64_t n, const char* gen) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* this_base = e->data.function.args[0];
        const Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1 &&
            expr_eq((Expr*)this_base, (Expr*)base)) {
            int64_t k = p * (n / q);
            Expr* pa[2] = { expr_new_symbol(gen), expr_new_integer(k) };
            Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = denom_subst_radical_to_gen(e->data.function.args[i],
                                                  base, n, gen);
    }
    Expr* new_head = denom_subst_radical_to_gen(e->data.function.head,
                                                 base, n, gen);
    Expr* result = expr_new_function(new_head, new_args, count);
    free(new_args);
    return eval_and_free(result);
}

/* Reverse: substitute Power[gen, k] -> Power[base, k * (1/n)] (so a
 * polynomial in `gen` becomes a polynomial in c^(1/n)). Returns a
 * fresh tree. */
static Expr* denom_subst_gen_to_radical(const Expr* e, const char* gen,
                                         const Expr* base, int64_t n) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, gen) == 0) {
        Expr* pa[2] = { expr_copy((Expr*)base), make_rational(1, n) };
        Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
        Expr* out = evaluate(pc);
        expr_free(pc);
        return out;
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(e->data.function.args[0]->data.symbol.name, gen) == 0) {
        const Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER) {
            int64_t k = exp->data.integer;
            /* gen^k = c^(k/n). Build directly. */
            Expr* pa[2] = { expr_copy((Expr*)base), make_rational(k, n) };
            Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = denom_subst_gen_to_radical(e->data.function.args[i],
                                                  gen, base, n);
    }
    Expr* new_head = denom_subst_gen_to_radical(e->data.function.head,
                                                 gen, base, n);
    Expr* result = expr_new_function(new_head, new_args, count);
    free(new_args);
    return eval_and_free(result);
}

/* Compute the inverse of `denom` (a polynomial in radicals over a
 * common positive integer base) via extended-Euclidean in
 * Q[α]/(α^n - c). Returns NULL when:
 *   - denom has no radicals, or has multiple bases.
 *   - the gcd of the substituted denom and α^n - c is not a nonzero
 *     constant (i.e. denom shares a factor with α^n - c — would mean
 *     the input expression is undefined when α takes that value).
 *
 * Soundness gate: only fires when c is a positive integer >= 2. */
static Expr* denom_compute_inverse(const Expr* denom) {
    const Expr* base = NULL;
    int64_t n = 0;
    /* Bail on a denominator that carries an explicit complex literal (I).  This
     * pass rationalises a REAL radical denominator via the conjugate over
     * Q[α]/(α^n - c); a complex denominator like x - (1 + I Sqrt[3])/2 substitutes
     * the radical for a generator but leaves I in the coefficients, so the
     * PolynomialExtendedGCD runs over Q(i) in {gen, x} and the generic
     * multivariate GCD blows up in exact_poly_div (Simplify/FullSimplify hang over
     * Q(i sqrt d)).  Clearing an i-times-radical denominator is not this pass's
     * job; decline so Together/Cancel (which handle it fine) do the work. */
    if (contains_explicit_complex(denom)) return NULL;
    if (!denom_collect_radical_base(denom, &base, &n)) return NULL;
    if (!base) return NULL;
    if (n < 2) return NULL;
    if (base->type != EXPR_INTEGER || base->data.integer < 2) return NULL;
    int64_t c = base->data.integer;

    /* Generate a fresh symbol name unlikely to clash with the user's
     * variables. The leading $ marks it as a system-generated name. */
    static int counter = 0;
    char gen[32];
    snprintf(gen, sizeof(gen), "$denomGen%d", counter++);

    /* p_in_a = denom rewritten as polynomial in `gen`. */
    Expr* p_in_a = denom_subst_radical_to_gen(denom, base, n, gen);
    if (!p_in_a) return NULL;

    /* Build the relation polynomial: gen^n - c. */
    Expr* gen_pow_args[2] = { expr_new_symbol(gen), expr_new_integer(n) };
    Expr* gen_pow = expr_new_function(expr_new_symbol(SYM_Power), gen_pow_args, 2);
    Expr* gen_pow_e = evaluate(gen_pow);
    expr_free(gen_pow);
    Expr* relation_args[2] = { gen_pow_e, expr_new_integer(-c) };
    Expr* relation_call = expr_new_function(expr_new_symbol(SYM_Plus),
                                            relation_args, 2);
    Expr* relation = evaluate(relation_call);
    expr_free(relation_call);

    /* Run PolynomialExtendedGCD[p_in_a, relation, gen]. */
    Expr* xgcd_args[3] = {
        p_in_a, relation, expr_new_symbol(gen)
    };
    Expr* xgcd_call = expr_new_function(expr_new_symbol(SYM_PolynomialExtendedGCD),
                                        xgcd_args, 3);
    Expr* xgcd_result = evaluate(xgcd_call);
    expr_free(xgcd_call);

    /* Validate result shape: List[gcd, List[u, v]]. */
    if (!xgcd_result ||
        xgcd_result->type != EXPR_FUNCTION ||
        xgcd_result->data.function.head->type != EXPR_SYMBOL ||
        xgcd_result->data.function.head->data.symbol.name != SYM_List ||
        xgcd_result->data.function.arg_count != 2) {
        if (xgcd_result) expr_free(xgcd_result);
        return NULL;
    }
    Expr* gcd_e = xgcd_result->data.function.args[0];
    Expr* coeffs = xgcd_result->data.function.args[1];
    if (!coeffs ||
        coeffs->type != EXPR_FUNCTION ||
        coeffs->data.function.head->type != EXPR_SYMBOL ||
        coeffs->data.function.head->data.symbol.name != SYM_List ||
        coeffs->data.function.arg_count < 1) {
        expr_free(xgcd_result);
        return NULL;
    }
    Expr* u_in_a = coeffs->data.function.args[0];

    /* Require gcd to be a nonzero rational constant (no `gen` in it).
     * If gcd has the gen, it shares a factor with α^n - c, meaning
     * denom vanishes at some root and the expression is undefined. */
    if (gcd_e->type != EXPR_INTEGER &&
        !is_rational_literal(gcd_e)) {
        expr_free(xgcd_result);
        return NULL;
    }
    if (gcd_e->type == EXPR_INTEGER && gcd_e->data.integer == 0) {
        expr_free(xgcd_result);
        return NULL;
    }

    /* Build inverse = u_in_a / gcd_e (still in `gen` form). */
    Expr* inv_gen_args[2] = { expr_copy(gcd_e), expr_new_integer(-1) };
    Expr* gcd_inv_call = expr_new_function(expr_new_symbol(SYM_Power),
                                           inv_gen_args, 2);
    Expr* gcd_inv = evaluate(gcd_inv_call);
    expr_free(gcd_inv_call);
    Expr* prod_args[2] = { expr_copy(u_in_a), gcd_inv };
    Expr* prod_call = expr_new_function(expr_new_symbol(SYM_Times),
                                        prod_args, 2);
    Expr* inv_in_gen = evaluate(prod_call);
    expr_free(prod_call);

    expr_free(xgcd_result);

    /* Substitute gen → c^(1/n) and re-evaluate. */
    Expr* inv_in_radical = denom_subst_gen_to_radical(inv_in_gen, gen, base, n);
    expr_free(inv_in_gen);
    return inv_in_radical;
}

/* Walk the expression tree, applying denominator rationalisation at
 * each Power[denom, k] subterm with k a negative rational. Returns
 * NULL when no rewrite fires. */
static Expr* simp_rationalize_denom_walk(const Expr* e,
                                         const AssumeCtx* ctx) {
    (void)ctx;
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return NULL;

    /* Recurse into children first so deeper rewrites bubble up. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = simp_rationalize_denom_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) {
                    new_args[j] = expr_copy(e->data.function.args[j]);
                }
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }
    Expr* current = NULL;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = evaluate(rebuilt);
        expr_free(rebuilt);
    }

    /* Check whether this node is a Power[denom, k] with k a negative
     * rational and denom a candidate for rationalisation. */
    const Expr* target = current ? current : e;
    if (target->type == EXPR_FUNCTION
        && target->data.function.head
        && target->data.function.head->type == EXPR_SYMBOL
        && target->data.function.head->data.symbol.name == SYM_Power
        && target->data.function.arg_count == 2) {
        const Expr* denom = target->data.function.args[0];
        const Expr* exp = target->data.function.args[1];
        bool is_neg_rational = false;
        int64_t k_p = 0, k_q = 1;
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
            is_neg_rational = true;
            k_p = exp->data.integer;
            k_q = 1;
        } else if (is_rational_literal(exp)) {
            int64_t pn, qn;
            if (is_rational(exp, &pn, &qn) && pn < 0 && qn == 1) {
                /* canonical Rational with denom 1 shouldn't occur, but
                 * be defensive. */
                is_neg_rational = true;
                k_p = pn;
                k_q = 1;
            }
        }
        /* For now, only handle k = -1. Power[denom, -k] for k > 1 can
         * be reduced to repeated multiplication of the rationalised
         * inverse, but isn't needed by the user's phase-2 cases. */
        if (is_neg_rational && k_p == -1 && k_q == 1) {
            Expr* inv = denom_compute_inverse(denom);
            if (inv) {
                if (current) expr_free(current);
                return inv;
            }
        }
    }
    return current;
}

Expr* simp_rationalize_denom(const Expr* e, const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = simp_rationalize_denom_walk(e, ctx);
    Expr* out = r ? r : expr_copy((Expr*)e);
    /* PrimeRebase post-pass on the FULL output. The user's RHS may use
     * a "compound" radical base (e.g. 4^(1/3)) while the inverse our
     * walker computes uses the PRIME base (2^(2/3)). Without rebasing
     * the full expression, the subtraction "2^(2/3) - 4^(1/3)" doesn't
     * cancel even though the values are equal. PrimeRebase normally
     * runs only at simp_dispatch on the input, not on candidates
     * produced inside simp_search, so we apply it here. */
    if (r) {
        Expr* rebased = transform_prime_rebase(out);
        if (rebased) {
            if (!expr_eq(rebased, out)) {
                expr_free(out);
                out = rebased;
            } else {
                expr_free(rebased);
            }
        }
    }
    if (dbg) simp_debug_log("RationalizeDenom", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}


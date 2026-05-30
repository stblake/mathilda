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
/* PrimeRebase: Power[c, e] -> Power[p, k*e] for c = p^k integer (k >= 2)  */
/* ----------------------------------------------------------------------- */

/* Soundness:  for positive integer c = p^k,
 *     c^e = (p^k)^e = p^(k*e)
 * holds for ALL complex e -- the (a^b)^c = a^(b*c) identity is sound
 * when a > 0, with no branch cut to worry about.
 *
 * Why this is needed:  Mathilda's canonical Power evaluator never rebases
 * composite integer bases ((2^2)^x stays as (2^2)^x), so factors like
 * 4^x and 2^x sit in different Times "base buckets" and the same-base
 * exponent combine in `times.c` cannot cancel them.  After rebasing all
 * such Power factors to a single canonical prime base, evaluate()
 * collapses the combined exponents in a single pass.
 *
 * Coverage in Simplify search:
 *     4^x * 2^(-x) * 2^(-x) - 1     ->  0  (top-level rebase)
 *     f[4^x * 2^(-x) * 2^(-x)] - f[1]  ->  0  (rebase inside f's arg)
 *     2^(2^(2 x) x) - 2^(x*4^x)     ->  0  (rebase inside an exponent)
 *
 * Branch-cut-sensitive shapes -- negative-integer-base split
 * ((-4)^x -> (-1)^x 4^x), constant-positive distribute
 * ((c1 c2)^e -> c1^e c2^e for ci > 0), and integer-exponent distribute
 * ((a b)^n -> a^n b^n for n integer) -- are handled by
 * `transform_power_distribute` (defined just below); see the comment
 * block there for the soundness argument and the dispatch wiring.
 */

/* If n is a perfect prime power p^k with k >= 2, set *p_out = p,
 * *k_out = k and return true.  Otherwise return false (n is < 4, prime,
 * or has at least two distinct prime factors).  Trial division up to
 * sqrt(n); since Mathilda's prime-base inputs are typically small literals
 * (4, 8, 9, 16, 25, 27, 32, 49, ...), this is microseconds-cheap. */
static bool prime_rebase_check(int64_t n, int64_t* p_out, int64_t* k_out) {
    if (n < 4) return false;
    int64_t p;
    if ((n & 1) == 0) {
        p = 2;
    } else {
        for (p = 3; ; p += 2) {
            if (p > n / p) return false;  /* n is prime */
            if (n % p == 0) break;
        }
    }
    int64_t k = 0;
    int64_t m = n;
    while (m > 1 && (m % p) == 0) {
        m /= p;
        k++;
    }
    if (m != 1 || k < 2) return false;
    *p_out = p;
    *k_out = k;
    return true;
}

/* Returns true if e contains any Power[c, _] with c an EXPR_INTEGER >= 4
 * that is a perfect prime power.  Cheap structural gate to avoid the full
 * walk + Expand reseed when nothing can fire (the common case). */
static bool has_rebaseable_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base && base->type == EXPR_INTEGER && base->data.integer >= 4) {
            int64_t p, k;
            if (prime_rebase_check(base->data.integer, &p, &k)) return true;
        }
    }
    if (has_rebaseable_power(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_rebaseable_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* Recursive structural copy that rebases every prime-power-base Power.
 * Walks children before checking the current node so nested rebases
 * (e.g. inside an exponent) bubble out.  Caller owns the returned tree. */
static Expr* prime_rebase_copy(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = prime_rebase_copy(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = prime_rebase_copy(e->data.function.args[i]);
        }
    }

    int64_t p_val, k_val;
    if (new_head && new_head->type == EXPR_SYMBOL
        && new_head->data.symbol == SYM_Power
        && n == 2 && new_args && new_args[0]
        && new_args[0]->type == EXPR_INTEGER
        && prime_rebase_check(new_args[0]->data.integer, &p_val, &k_val)) {
        /* Build Power[p, Times[k, new_args[1]]] -- new_args[1] is the
         * already-rebased exponent (so nested c'^e' inside the exponent
         * has already been rewritten before we get here). */
        Expr* p_expr = expr_new_integer(p_val);
        Expr* k_expr = expr_new_integer(k_val);
        Expr* times_args[2] = { k_expr, new_args[1] };
        Expr* times_call = expr_new_function(
            expr_new_symbol("Times"), times_args, 2);
        Expr* power_args[2] = { p_expr, times_call };
        Expr* result = expr_new_function(
            expr_new_symbol("Power"), power_args, 2);
        expr_free(new_args[0]);  /* the original integer c */
        free(new_args);
        expr_free(new_head);
        return result;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_prime_rebase_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    if (!has_rebaseable_power(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PrimeRebase", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* rebased = prime_rebase_copy(e);
    if (!rebased) rebased = expr_copy((Expr*)e);

    /* Re-evaluate so canonical Times same-base combine collapses the
     * rebased Power factors (e.g. 2^(2x) * 2^(-x) * 2^(-x) -> 2^0 = 1). */
    Expr* result = eval_and_free(rebased);
    if (!result) result = expr_copy((Expr*)e);

    if (dbg) simp_debug_log("PrimeRebase", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

Expr* transform_prime_rebase(const Expr* e) {
    return simp_memo_wrap(e, "$PrimeRebase", transform_prime_rebase_impl);
}

/* ----------------------------------------------------------------------- */
/* PowerOneify: combine `A * Power[A, e]` -> `Power[A, e+1]` in Times      */
/* ----------------------------------------------------------------------- */

/* Soundness: Power[a, 1] = a holds for ALL a, so
 *     A * Power[A, e] = Power[A, 1] * Power[A, e] = Power[A, e + 1]
 * is universally valid -- no branch cut, no positivity assumption.
 *
 * Why this is needed: Mathilda's Times-canonical-form same-base combine
 * groups factors by the base of any wrapping Power[base, exp].  A bare
 * factor A whose canonical form is itself a Power expression (e.g.
 * Power[x, -1] = 1/x) does NOT get re-bucketed as Power[A, 1] before
 * grouping, so Times[Power[x,-1], Power[Power[x,-1], Log[2]]] keeps the
 * two factors as distinct same-base bucket entries.  This blocks
 *     (1/x)^Log[2] / x - (1/x)^(1 + Log[2])  ->  0
 * from collapsing because the LHS Times can't combine the two (1/x)
 * factors into one.
 *
 * Implementation: walk the tree; in every Times node, look for any pair
 * (i, j) where args[j] = Power[B, e] and args[i] is structurally equal
 * to B.  Replace args[j] with Power[B, e + 1] and drop args[i] (which is
 * the bare A = Power[B, 1] absorbed into the exponent).  Then evaluate so
 * canonical Plus folds the new exponent and so any nested rebase cascades.
 *
 * Inert when no Times node contains an A-and-Power[A,_] pair (the common
 * case): one structural pass, microseconds.
 */

static Expr* power_oneify_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = power_oneify_walk(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = power_oneify_walk(e->data.function.args[i]);
        }
    }

    /* Only Times triggers the implicit-one combine. */
    bool combined = false;
    if (new_head && new_head->type == EXPR_SYMBOL
        && new_head->data.symbol == SYM_Times && n >= 2) {
        for (size_t i = 0; i < n; i++) {
            if (!new_args[i]) continue;
            Expr* a = new_args[i];
            for (size_t j = 0; j < n; j++) {
                if (i == j || !new_args[j]) continue;
                Expr* p = new_args[j];
                if (p->type != EXPR_FUNCTION) continue;
                Expr* ph = p->data.function.head;
                if (!(ph && ph->type == EXPR_SYMBOL
                      && ph->data.symbol == SYM_Power
                      && p->data.function.arg_count == 2)) continue;
                if (!expr_eq(p->data.function.args[0], a)) continue;
                /* args[i] = A, args[j] = Power[A, e].  Combine. */
                Expr* old_exp = p->data.function.args[1];
                Expr* one = expr_new_integer(1);
                Expr* plus_args[2] = { expr_copy(old_exp), one };
                Expr* new_exp = expr_new_function(
                    expr_new_symbol("Plus"), plus_args, 2);
                Expr* power_args[2] = {
                    expr_copy(p->data.function.args[0]), new_exp };
                Expr* new_power = expr_new_function(
                    expr_new_symbol("Power"), power_args, 2);
                expr_free(new_args[j]);
                new_args[j] = new_power;
                expr_free(new_args[i]);
                new_args[i] = NULL;
                combined = true;
                break;  /* args[i] consumed; move to next i */
            }
        }
    }

    if (combined) {
        /* Compact dropped slots. */
        size_t out_n = 0;
        for (size_t i = 0; i < n; i++) {
            if (new_args[i]) new_args[out_n++] = new_args[i];
        }
        Expr* times_call = expr_new_function(new_head, new_args, out_n);
        free(new_args);
        Expr* result = eval_and_free(times_call);
        if (!result) result = expr_copy((Expr*)e);
        return result;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_power_oneify_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* result = power_oneify_walk(e);
    if (!result) result = expr_copy((Expr*)e);
    if (dbg) simp_debug_log("PowerOneify", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

Expr* transform_power_oneify(const Expr* e) {
    return simp_memo_wrap(e, "$PowerOneify", transform_power_oneify_impl);
}

/* ----------------------------------------------------------------------- */
/* PowerDistribute: distribute Power over Times in three universally       */
/* sound shapes.  Walks the tree once; each shape independently rewrites   */
/* matching nodes in place, then evaluate() is called on the result.       */
/*                                                                         */
/*   (A) Power[neg_int, e]  ->  Power[-1, e] * Power[|neg_int|, e]         */
/*       Soundness: principal-branch identity for any e (real or           */
/*       symbolic).  Mathilda uses principal branches everywhere.           */
/*       Why: lets prime_rebase + Times same-base merge cancel mixed-sign  */
/*       products like (-4)^x * (-2)^(-x) * 2^(-x) -> 1.                   */
/*                                                                         */
/*   (B) Power[Times[c1,...,ck, u1,...,um], e]                             */
/*           ->  Times[c1^e, ..., ck^e, Power[Times[u1,...,um], e]]        */
/*       when each ci is a constant positive (literal positive numeric,    */
/*       or a recognised positive constant symbol like Pi/E/...).          */
/*       Soundness: (a b)^c = a^c b^c whenever a > 0, valid for any b      */
/*       and any c.  Why: collapses identities like                        */
/*           Exp[x] Exp[y] 2^x 2^y - (2 Exp[1])^(x+y) -> 0                 */
/*       which need (2 E)^(x+y) -> 2^(x+y) E^(x+y) to align with the LHS.  */
/*                                                                         */
/*   (C) Power[Times[u1,...,un], e]  ->  Power[u1, e] * ... * Power[un, e] */
/*       and Power[Power[u, p], e]   ->  Power[u, p*e]                     */
/*       when prov_int(ctx, e) (the exponent is provably integer, either   */
/*       as a literal or via an Element[_, Integers] assumption).          */
/*       Soundness: integer exponents distribute through products and      */
/*       through nested Power without branch cuts, for any complex base.   */
/*       Why: handles (y/x)^(-n) - x^n y^(-n) -> 0 under                   */
/*       Element[n, Integers].                                             */
/*                                                                         */
/* All three shapes are inert when no matching node exists.  After any     */
/* rewrite we evaluate() the result so the canonical Times same-base       */
/* merge collapses adjacent Power[u, ?] factors.                           */
/* ----------------------------------------------------------------------- */

static bool is_constant_positive_factor(const Expr* x) {
    if (!x) return false;
    if (numeric_sign(x) == 1) return true;
    if (x->type == EXPR_SYMBOL && is_positive_constant_symbol(x->data.symbol))
        return true;
    return false;
}

/* Cheap structural gate for shape (A). */
static bool has_distributable_power(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp_  = e->data.function.args[1];
        /* (A) negative integer base */
        if (base && base->type == EXPR_INTEGER && base->data.integer < -1)
            return true;
        /* (B)/(C) Power[Times[...], e] */
        if (base && base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Times
            && base->data.function.arg_count > 0) {
            /* (B) any constant-positive factor? */
            for (size_t i = 0; i < base->data.function.arg_count; i++) {
                if (is_constant_positive_factor(base->data.function.args[i]))
                    return true;
            }
            /* (C) integer exponent triggers full distribute */
            if (prov_int(ctx, exp_)) return true;
        }
        /* (C) Power[Power[u, p], e] with e integer */
        if (base && base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Power
            && base->data.function.arg_count == 2
            && prov_int(ctx, exp_)) {
            return true;
        }
    }
    if (has_distributable_power(head, ctx)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_distributable_power(e->data.function.args[i], ctx)) return true;
    }
    return false;
}

static Expr* power_distribute_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = power_distribute_walk(e->data.function.head, ctx);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = power_distribute_walk(e->data.function.args[i], ctx);
        }
    }

    /* Only Power[_, _] triggers any rewrite at this level. */
    if (!(new_head && new_head->type == EXPR_SYMBOL
          && new_head->data.symbol == SYM_Power
          && n == 2 && new_args && new_args[0] && new_args[1])) {
        Expr* result = expr_new_function(new_head, new_args, n);
        if (new_args) free(new_args);
        return result;
    }

    Expr* base = new_args[0];
    Expr* exp_ = new_args[1];

    /* (A) Power[neg_int, e] -> Power[-1, e] * Power[|neg_int|, e] */
    if (base->type == EXPR_INTEGER && base->data.integer < -1) {
        int64_t v = base->data.integer;
        Expr* p1_args[2] = { expr_new_integer(-1), expr_copy(exp_) };
        Expr* p1 = expr_new_function(expr_new_symbol("Power"), p1_args, 2);
        Expr* p2_args[2] = { expr_new_integer(-v), exp_ };
        Expr* p2 = expr_new_function(expr_new_symbol("Power"), p2_args, 2);
        Expr* tm_args[2] = { p1, p2 };
        Expr* tm = expr_new_function(expr_new_symbol("Times"), tm_args, 2);
        expr_free(base);             /* original neg-int base */
        free(new_args);
        expr_free(new_head);
        return tm;
    }

    /* (B) and (C) Power[Times[args...], e] */
    if (base->type == EXPR_FUNCTION && base->data.function.head
        && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Times
        && base->data.function.arg_count > 0) {
        size_t bn = base->data.function.arg_count;
        Expr** ba = base->data.function.args;
        bool exp_is_int = prov_int(ctx, exp_);

        /* (C) integer exponent: distribute over every factor. */
        if (exp_is_int) {
            Expr** powers = (Expr**)malloc(sizeof(Expr*) * bn);
            for (size_t i = 0; i < bn; i++) {
                Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                powers[i] = expr_new_function(
                    expr_new_symbol("Power"), pa, 2);
            }
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), powers, bn);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }

        /* (B) split off constant-positive factors. */
        size_t pos_count = 0;
        for (size_t i = 0; i < bn; i++) {
            if (is_constant_positive_factor(ba[i])) pos_count++;
        }
        if (pos_count > 0 && pos_count < bn) {
            /* Split: pos_count Power[ci, e] factors + Power[Times[rest], e] */
            Expr** out = (Expr**)malloc(sizeof(Expr*) * (pos_count + 1));
            size_t out_i = 0;
            size_t rest_n = bn - pos_count;
            Expr** rest = (Expr**)malloc(sizeof(Expr*) * rest_n);
            size_t rest_i = 0;
            for (size_t i = 0; i < bn; i++) {
                if (is_constant_positive_factor(ba[i])) {
                    Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                    out[out_i++] = expr_new_function(
                        expr_new_symbol("Power"), pa, 2);
                } else {
                    rest[rest_i++] = expr_copy(ba[i]);
                }
            }
            Expr* rest_times;
            if (rest_n == 1) {
                rest_times = rest[0];
                free(rest);
            } else {
                rest_times = expr_new_function(
                    expr_new_symbol("Times"), rest, rest_n);
            }
            Expr* pa[2] = { rest_times, expr_copy(exp_) };
            out[out_i++] = expr_new_function(
                expr_new_symbol("Power"), pa, 2);
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), out, pos_count + 1);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }
        if (pos_count == bn) {
            /* All factors constant-positive: full distribute. */
            Expr** powers = (Expr**)malloc(sizeof(Expr*) * bn);
            for (size_t i = 0; i < bn; i++) {
                Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                powers[i] = expr_new_function(
                    expr_new_symbol("Power"), pa, 2);
            }
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), powers, bn);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }
    }

    /* (C) Power[Power[u, p], e] -> Power[u, p*e] when e provably integer. */
    if (base->type == EXPR_FUNCTION && base->data.function.head
        && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Power
        && base->data.function.arg_count == 2
        && prov_int(ctx, exp_)) {
        Expr* u  = base->data.function.args[0];
        Expr* p  = base->data.function.args[1];
        Expr* tm_args[2] = { expr_copy(p), exp_ };
        Expr* prod = expr_new_function(expr_new_symbol("Times"), tm_args, 2);
        Expr* pa[2] = { expr_copy(u), prod };
        Expr* po = expr_new_function(expr_new_symbol("Power"), pa, 2);
        expr_free(base);
        free(new_args);
        expr_free(new_head);
        return po;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_power_distribute_impl(const Expr* e,
                                              const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    if (!has_distributable_power(e, ctx)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PowerDistribute", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* split = power_distribute_walk(e, ctx);
    if (!split) split = expr_copy((Expr*)e);

    /* Re-evaluate so canonical Times same-base merge collapses the
     * newly-introduced Power[u, e] factors against existing same-base
     * factors elsewhere in the surrounding Times. */
    Expr* result = eval_and_free(split);
    if (!result) result = expr_copy((Expr*)e);

    if (dbg) simp_debug_log("PowerDistribute", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

/* Note: ctx-dependent, so we cannot use simp_memo_wrap (which keys on
 * the input expression alone).  PowerDistribute is invoked at most a
 * few times per Simplify call from simp_dispatch and is structurally
 * cheap (single tree walk, no recursive search), so the caching loss
 * is negligible. */
Expr* transform_power_distribute(const Expr* e,
                                         const AssumeCtx* ctx) {
    return transform_power_distribute_impl(e, ctx);
}

/* ----------------------------------------------------------------------- */
/* RadicalCanon: split Power[Rational[a,b],q] and rationalise -p/q powers  */
/* ----------------------------------------------------------------------- */

/* Two related canonicalisations that Mathilda's standard Power evaluator
 * does NOT do, and which leave equivalent expressions in distinct shapes:
 *
 * (1)  Power[Rational[a, b], q]                       (a, b positive ints)
 *           ->  Power[a, q] * Power[b, -q]
 *      Soundness: a, b > 0, so (a/b)^q = a^q * b^(-q) for all complex q
 *      with no branch cut.  Worked example:
 *          Sqrt[1/2] = Power[Rational[1,2], 1/2]
 *                    -> Power[1, 1/2] * Power[2, -1/2]
 *                    -> 1/Sqrt[2]                    (after evaluate)
 *
 * (2)  Power[a, q]                          (a positive int >= 2,
 *                                            q rational, q < 0, denom > 1)
 *           ->  Power[a, r] * Rational[1, a^k]
 *      where k = ceil(-q) >= 1, r = q + k in [0, 1).
 *      Soundness: a > 0, so (a)^q = (a)^r / a^k.  Not folded back by
 *      canonical Times because Rational[1, a^k] occupies a distinct base
 *      bucket from Power[a, r].  Worked example:
 *          Power[2, -1/2]  ->  Power[2, 1/2] * Rational[1, 2]
 *                            =  Sqrt[2] / 2
 *
 * Why both: (1) alone pushes Sqrt[1/2] to 1/Sqrt[2] but leaves it in
 * the negative-exponent form; the surrounding additive context still
 * doesn't see it as same-shape with a Sqrt[2]/6 sibling.  (2) finishes
 * the job by rationalising the denominator so additive cancellations
 * fire.
 *
 * Coverage:
 *     Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2] + Cos[y]/(3 Sqrt[6]) +
 *     Sin[y]/(3 Sqrt[2]) - (Sqrt[6]/3) Sin[x+Pi/6] -
 *     (Sqrt[6]/9) Sin[y+Pi/6]                                  ->  0
 *
 * Inert when the input contains no Power[Rational, _] or no
 * Power[positive_int, negative_rational].
 */

/* Compute a^k with overflow-guard.  Returns false on overflow or
 * a, k out of range; the caller falls back to leaving the term alone. */
static bool radical_canon_pow_int(int64_t a, int64_t k, int64_t* out) {
    if (k < 0 || k > 62) return false;
    int64_t v = 1;
    for (int64_t i = 0; i < k; i++) {
        if (a != 0 && v > INT64_MAX / (a > 0 ? a : -a)) return false;
        v *= a;
    }
    *out = v;
    return true;
}

/* Match Rational[num, den] with den > 1, returning num and den.  The
 * canonical evaluator never produces Rational[_, 1] (it folds to a bare
 * integer) so we only see denominator > 1 here. */
static bool radical_canon_get_rational(const Expr* e,
                                        int64_t* num, int64_t* den) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Rational) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* na = e->data.function.args[0];
    Expr* da = e->data.function.args[1];
    if (na->type != EXPR_INTEGER || da->type != EXPR_INTEGER) return false;
    *num = na->data.integer;
    *den = da->data.integer;
    return true;
}

/* If e = Power[positive_integer_>=2, Rational[num<0, den>1]], rewrite to
 *     Times[Power[a, Rational[r_num, den]], Rational[1, a^k]]
 * Returns NULL when the rule does not apply. */
static Expr* radical_canon_rationalise_negexp(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return NULL;
    if (e->data.function.arg_count != 2) return NULL;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    if (base->type != EXPR_INTEGER || base->data.integer < 2) return NULL;
    int64_t num, den;
    if (!radical_canon_get_rational(exp, &num, &den)) return NULL;
    if (num >= 0 || den <= 1) return NULL;

    int64_t a = base->data.integer;
    int64_t neg_num = -num;
    int64_t k = neg_num / den;
    if (neg_num % den != 0) k++;
    int64_t r_num = num + k * den;  /* r_num in [0, den) */

    int64_t a_k;
    if (!radical_canon_pow_int(a, k, &a_k)) return NULL;
    if (a_k < 1) return NULL;

    /* Build Power[a, Rational[r_num, den]].  When r_num == 0, Power
     * evaluates to 1; when r_num == den (impossible here since
     * r_num < den), Power[a, 1] = a. */
    Expr* pow_args[2] = {
        expr_new_integer(a),
        expr_new_function(
            expr_new_symbol("Rational"),
            (Expr*[]){ expr_new_integer(r_num), expr_new_integer(den) }, 2)
    };
    Expr* pow_call = expr_new_function(
        expr_new_symbol("Power"), pow_args, 2);

    /* Rational[1, a^k].  When a^k == 1 (impossible for a >= 2, k >= 1),
     * the evaluator folds to the integer 1. */
    Expr* recip = expr_new_function(
        expr_new_symbol("Rational"),
        (Expr*[]){ expr_new_integer(1), expr_new_integer(a_k) }, 2);

    Expr* times_args[2] = { pow_call, recip };
    Expr* times_call = expr_new_function(
        expr_new_symbol("Times"), times_args, 2);
    Expr* out = evaluate(times_call);
    if (!out) out = expr_copy(times_call);  /* fall back to unevaluated form */
    expr_free(times_call);
    return out;
}

/* If e = Power[Rational[a, b], q] with positive integers a, b > 1,
 * rewrite to Times[Power[a, q], Power[b, -q]].  Returns NULL when not
 * applicable. */
static Expr* radical_canon_split_rational_base(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return NULL;
    if (e->data.function.arg_count != 2) return NULL;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    int64_t a, b;
    if (!radical_canon_get_rational(base, &a, &b)) return NULL;
    if (a < 1 || b < 2) return NULL;

    /* Power[a, q] * Power[b, -q].  a may be 1 (Sqrt[1/n] case): then
     * Power[1, q] evaluates to 1 and we're left with Power[b, -q]. */
    Expr* pow_a_args[2] = {
        expr_new_integer(a), expr_copy(exp)
    };
    Expr* pow_a = expr_new_function(
        expr_new_symbol("Power"), pow_a_args, 2);

    Expr* neg_exp_args[2] = {
        expr_new_integer(-1), expr_copy(exp)
    };
    Expr* neg_exp = expr_new_function(
        expr_new_symbol("Times"), neg_exp_args, 2);
    Expr* pow_b_args[2] = {
        expr_new_integer(b), neg_exp
    };
    Expr* pow_b = expr_new_function(
        expr_new_symbol("Power"), pow_b_args, 2);

    Expr* times_args[2] = { pow_a, pow_b };
    Expr* times_call = expr_new_function(
        expr_new_symbol("Times"), times_args, 2);
    Expr* out = eval_and_free(times_call);
    return out;
}

/* Recursive walker: applies both rules at each Power node and recurses
 * into children.  Returns a freshly allocated tree (caller owns). */
static Expr* radical_canon_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = radical_canon_walk(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = radical_canon_walk(e->data.function.args[i]);
        }
    }

    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);

    /* Try Rational-base split first (may produce a Power node that the
     * negative-exponent rule then handles). */
    Expr* split = radical_canon_split_rational_base(rebuilt);
    if (split) {
        expr_free(rebuilt);
        rebuilt = split;
    }
    Expr* rat = radical_canon_rationalise_negexp(rebuilt);
    if (rat) {
        expr_free(rebuilt);
        rebuilt = rat;
    }
    return rebuilt;
}

static Expr* transform_radical_canon_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* walked = radical_canon_walk(e);
    if (!walked) walked = expr_copy((Expr*)e);
    /* Re-evaluate so canonical Plus/Times fold the rewritten tree --
     * the walker only evaluates at each Power-level rewrite, leaving the
     * surrounding Plus/Times in unevaluated form (e.g. an arg pair like
     * 1/2 Sqrt[2] + -1 (1/2 Sqrt[2]) won't auto-cancel without this).
     * evaluate() does not consume its argument, so `walked` is still ours
     * to release afterwards (refcount-safe even when evaluate returns the
     * same shared node for an already-stable input). */
    Expr* result = evaluate(walked);
    if (!result) result = expr_copy((Expr*)e);
    expr_free(walked);
    if (dbg) simp_debug_log("RadicalCanon", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

Expr* transform_radical_canon(const Expr* e) {
    return simp_memo_wrap(e, "$RadicalCanon", transform_radical_canon_impl);
}


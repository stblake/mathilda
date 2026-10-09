#include "simp.h"
#include "simp_internal.h"
#include "simp_trigexp_zero.h"
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
#include "flint_qqbar.h"

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
/* builtin_simplify                                                        */
/* ----------------------------------------------------------------------- */

Expr* read_dollar_assumptions(void) {
    /* Read the OwnValue directly. We must NOT evaluate $Assumptions, because
     * once an assumption like Element[x, Reals] becomes the bound value, our
     * own Element evaluator would recurse on it (Element reads $Assumptions
     * to decide -> evaluator fires the OwnValue rule -> Element[x, Reals]
     * gets re-evaluated -> ...). The first OwnValue rule on a symbol is its
     * current value (newest first); we deep-copy its replacement. */
    Rule* r = symtab_get_own_values("$Assumptions");
    if (!r || !r->replacement) return expr_new_symbol(SYM_True);
    return expr_copy(r->replacement);
}

/* ----------------------------------------------------------------------- */
/* builtin_element -- Element[x, Domain]                                   */
/* ----------------------------------------------------------------------- */

static bool is_complex_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2;
}

bool is_rational_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2;
}

/* True iff `r` is exactly representable as a 64-bit integer. */
static bool real_is_integer(double r) {
    if (r != r) return false;                       /* NaN */
    if (r > 9.2233720368547758e18) return false;    /* > INT64_MAX */
    if (r < -9.2233720368547758e18) return false;
    long long i = (long long)r;
    return (double)i == r;
}

/* Element[x, dom] decision: 1 = True, 0 = False, -1 = undetermined. Exported
 * (declared in simp_internal.h) so Refine can decide Element statements under
 * its own assumption context. */
int element_decide(const Expr* x, const char* dom, const AssumeCtx* ctx) {
    if (!x || !dom) return -1;

    /* Direct fact lookup is always safe regardless of domain. */
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            if (fact_in_domain(ctx->facts[i], x, dom)) return 1;
        }
    }

    if (strcmp(dom, "Integers") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (x->type == EXPR_REAL) return real_is_integer(x->data.real) ? 1 : 0;
        if (is_rational_literal(x)) return 0;
        if (is_complex_literal(x)) return 0;
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Rationals") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (is_rational_literal(x)) return 1;
        if (x->type == EXPR_REAL) return 1;             /* every double is dyadic */
        if (is_complex_literal(x)) return 0;
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Algebraics") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (is_rational_literal(x)) return 1;
        if (x->type == EXPR_REAL) return 1;
        if (is_complex_literal(x)) return 1;            /* canonical Complex parts are rational */
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Reals") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT || x->type == EXPR_REAL) return 1;
        if (is_rational_literal(x)) return 1;
        if (is_complex_literal(x)) {
            /* canonical Complex always carries a non-zero imaginary part */
            Expr* im = x->data.function.args[1];
            if (im->type == EXPR_INTEGER && im->data.integer == 0) return 1;
            return 0;
        }
        if (prov_re(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Complexes") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT || x->type == EXPR_REAL) return 1;
        if (is_rational_literal(x)) return 1;
        if (is_complex_literal(x)) return 1;
        if (prov_re(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Booleans") == 0) {
        if (x->type == EXPR_SYMBOL) {
            if (x->data.symbol.name == SYM_True)  return 1;
            if (x->data.symbol.name == SYM_False) return 1;
        }
        return -1;
    }

    if (strcmp(dom, "Primes") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) {
            Expr* primeq = call_unary_copy("PrimeQ", x);
            int ans = -1;
            if (primeq && primeq->type == EXPR_SYMBOL) {
                if (primeq->data.symbol.name == SYM_True) ans = 1;
                if (primeq->data.symbol.name == SYM_False) ans = 0;
            }
            if (primeq) expr_free(primeq);
            return ans;
        }
        return -1;
    }

    if (strcmp(dom, "Composites") == 0) {
        if ((x->type == EXPR_INTEGER && x->data.integer >= 2) || x->type == EXPR_BIGINT) {
            Expr* primeq = call_unary_copy("PrimeQ", x);
            int ans = -1;
            if (primeq && primeq->type == EXPR_SYMBOL) {
                if (primeq->data.symbol.name == SYM_True) ans = 0;
                if (primeq->data.symbol.name == SYM_False) ans = 1;
            }
            if (primeq) expr_free(primeq);
            return ans;
        }
        return -1;
    }

    /* Sign domains as *queried* membership, decided by the sign provers (which
     * already fold numeric literals and assumption facts). Positive excludes 0;
     * a provably nonpositive x is therefore definitely not Positive, etc. */
    if (strcmp(dom, "Positive") == 0) {
        if (prov_pos(ctx, x)) return 1;
        if (prov_np(ctx, x))  return 0;   /* x <= 0 => not positive */
        return -1;
    }
    if (strcmp(dom, "Negative") == 0) {
        if (prov_neg(ctx, x)) return 1;
        if (prov_nn(ctx, x))  return 0;   /* x >= 0 => not negative */
        return -1;
    }
    if (strcmp(dom, "NonNegative") == 0) {
        if (prov_nn(ctx, x))  return 1;
        if (prov_neg(ctx, x)) return 0;   /* x < 0 => not nonnegative */
        return -1;
    }
    if (strcmp(dom, "NonPositive") == 0) {
        if (prov_np(ctx, x))  return 1;
        if (prov_pos(ctx, x)) return 0;   /* x > 0 => not nonpositive */
        return -1;
    }

    return -1;
}

Expr* builtin_element(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;

    Expr* x   = res->data.function.args[0];
    Expr* dom = res->data.function.args[1];

    /* Element[{x1, x2, ...}, dom] and Element[x1|x2|..., dom] are both
     * shorthand for the conjunction Element[x1, dom] && ... && Element[xN, dom].
     * We collapse to a single True/False only when every component decides;
     * otherwise we leave the original Element[{...}, dom] / Element[x|y, dom]
     * in place so downstream consumers (Simplify's AssumeCtx; see ctx_walk in
     * simp_assume.c) can still treat it as a multi-variable real/integer/...
     * assumption.
     *
     * Threading to a `List` of Element calls (the old behaviour) was
     * incorrect: a list of partial Element facts in the assumption position
     * would cause threading to split the joint assumption into per-variable
     * Simplify runs that each only see one of the facts. */
    if (x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        (x->data.function.head->data.symbol.name == SYM_List ||
         x->data.function.head->data.symbol.name == SYM_Alternatives)) {
        size_t n = x->data.function.arg_count;
        bool all_true = (n > 0);
        bool any_false = false;
        for (size_t i = 0; i < n; i++) {
            Expr* sub_args[2] = { expr_copy(x->data.function.args[i]), expr_copy(dom) };
            Expr* call = expr_new_function(expr_new_symbol(SYM_Element), sub_args, 2);
            Expr* sub = evaluate(call);
            expr_free(call);
            if (sub && sub->type == EXPR_SYMBOL) {
                if (sub->data.symbol.name == SYM_True) {
                    /* keep all_true */
                } else if (sub->data.symbol.name == SYM_False) {
                    any_false = true;
                    all_true = false;
                } else {
                    all_true = false;
                }
            } else {
                all_true = false;
            }
            if (sub) expr_free(sub);
            if (any_false) break;
        }
        if (any_false) return expr_new_symbol(SYM_False);
        if (all_true)  return expr_new_symbol(SYM_True);
        return NULL;
    }

    if (dom->type != EXPR_SYMBOL) return NULL;
    const char* d = dom->data.symbol.name;

    /* Build context from current $Assumptions. */
    Expr* dollar = read_dollar_assumptions();
    AssumeCtx* ctx = assume_ctx_from_expr(dollar);
    expr_free(dollar);

    int decision = element_decide(x, d, ctx);
    assume_ctx_free(ctx);

    if (decision == 1)  return expr_new_symbol(SYM_True);
    if (decision == 0)  return expr_new_symbol(SYM_False);
    return NULL;
}

static Expr* combine_with_and(Expr* a, Expr* b) {
    /* Both inputs owned and consumed. */
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_And), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ----------------------------------------------------------------------- */
/* Equation / inequality rebalancing                                       */
/*                                                                         */
/* For a binary relation `lhs OP rhs`, compute d = lhs - rhs as an         */
/* evaluated Plus, then rewrite as `pos OP neg` after dividing through by  */
/* the GCD of integer coefficients. Negative-coefficient terms move to    */
/* the opposite side. The result is correctness-preserving for both       */
/* equality and ordering relations (we never multiply or divide by a      */
/* negative quantity, only the positive integer GCD).                     */
/* ----------------------------------------------------------------------- */

bool simp_eq_head_sym(const Expr* e, const char* name) {
    return head_is(e, intern_symbol(name));
}

/* Decompose a Plus term into integer-coefficient * rest. Returns false
 * for terms whose leading numeric factor isn't an int64 (Real, BigInt,
 * Rational), which signals the caller to skip rebalancing -- mixing in
 * those would risk losing precision or introducing fractions. */
static bool simp_plus_term_int_coeff(const Expr* term, int64_t* coef,
                                     Expr** rest_out) {
    if (term->type == EXPR_INTEGER) {
        *coef = term->data.integer;
        *rest_out = expr_new_integer(1);
        return true;
    }
    if (term->type == EXPR_BIGINT || term->type == EXPR_REAL) return false;

    if (simp_eq_head_sym(term, "Times") &&
        term->data.function.arg_count >= 1) {
        const Expr* a0 = term->data.function.args[0];
        if (a0->type == EXPR_INTEGER) {
            *coef = a0->data.integer;
            size_t n = term->data.function.arg_count;
            if (n == 2) {
                *rest_out = expr_copy(term->data.function.args[1]);
            } else {
                Expr** args = (Expr**)calloc(n - 1, sizeof(Expr*));
                for (size_t i = 1; i < n; i++) {
                    args[i - 1] = expr_copy(term->data.function.args[i]);
                }
                *rest_out = expr_new_function(
                    expr_new_symbol(SYM_Times), args, n - 1);
                free(args);
            }
            return true;
        }
        if (a0->type == EXPR_BIGINT || a0->type == EXPR_REAL) return false;
        if (simp_eq_head_sym(a0, "Rational")) return false;
    }

    /* Generic term: implicit coefficient 1, rest = term. */
    *coef = 1;
    *rest_out = expr_copy((Expr*)term);
    return true;
}

/* Build `c * rest`, dropping a coefficient of 1 and Times wrappers when
 * rest = 1. Takes ownership of `rest`. */
static Expr* simp_make_term(int64_t c, Expr* rest) {
    if (rest->type == EXPR_INTEGER && rest->data.integer == 1) {
        expr_free(rest);
        return expr_new_integer(c);
    }
    if (c == 1) return rest;
    /* Flatten into existing Times; otherwise wrap. */
    if (simp_eq_head_sym(rest, "Times")) {
        size_t n = rest->data.function.arg_count;
        Expr** args = (Expr**)calloc(n + 1, sizeof(Expr*));
        args[0] = expr_new_integer(c);
        for (size_t i = 0; i < n; i++) {
            args[i + 1] = expr_copy(rest->data.function.args[i]);
        }
        Expr* out = expr_new_function(
            expr_new_symbol(SYM_Times), args, n + 1);
        free(args);
        expr_free(rest);
        return out;
    }
    Expr* args[2] = { expr_new_integer(c), rest };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* Returns NULL when no rebalanced form is produced (non-int64 coeffs,
 * fully symbolic d, or d = 0). The caller compares scores. */
static Expr* simp_try_rebalance_relation(const Expr* relation) {
    if (!relation || relation->type != EXPR_FUNCTION) return NULL;
    if (relation->data.function.arg_count != 2) return NULL;
    const Expr* h = relation->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return NULL;
    const char* hn = h->data.symbol.name;
    bool ok = (hn == SYM_Equal ||
               hn == SYM_Unequal ||
               hn == SYM_Less ||
               hn == SYM_LessEqual ||
               hn == SYM_Greater ||
               hn == SYM_GreaterEqual);
    if (!ok) return NULL;

    /* d = lhs - rhs, evaluated. */
    Expr* neg_args[2] = {
        expr_new_integer(-1),
        expr_copy(relation->data.function.args[1])
    };
    Expr* neg_rhs = expr_new_function(
        expr_new_symbol(SYM_Times), neg_args, 2);
    Expr* d_args[2] = {
        expr_copy(relation->data.function.args[0]),
        neg_rhs
    };
    Expr* d_call = expr_new_function(
        expr_new_symbol(SYM_Plus), d_args, 2);
    Expr* d_sum = eval_and_free(d_call);
    /* Expand so Times[2, Plus[...]] partitions term-by-term. The threaded
     * input may already have collected common factors via Collect, which
     * defeats coefficient-level rebalancing. */
    Expr* exp_args[1] = { d_sum };
    Expr* d_exp_call = expr_new_function(
        expr_new_symbol(SYM_Expand), exp_args, 1);
    Expr* d = eval_and_free(d_exp_call);

    Expr* d_singleton[1];
    Expr** terms;
    size_t n;
    if (simp_eq_head_sym(d, "Plus")) {
        n = d->data.function.arg_count;
        terms = d->data.function.args;
    } else {
        d_singleton[0] = d;
        terms = d_singleton;
        n = 1;
    }
    if (n == 0) { expr_free(d); return NULL; }

    /* Extract integer coefficients. Bail on non-int64. */
    int64_t* coefs = (int64_t*)calloc(n, sizeof(int64_t));
    Expr** rests = (Expr**)calloc(n, sizeof(Expr*));
    bool ok2 = true;
    for (size_t i = 0; i < n; i++) {
        if (!simp_plus_term_int_coeff(terms[i], &coefs[i], &rests[i])) {
            ok2 = false;
            for (size_t j = 0; j < i; j++) expr_free(rests[j]);
            break;
        }
    }
    if (!ok2) {
        free(coefs);
        free(rests);
        expr_free(d);
        return NULL;
    }

    /* GCD of |coefs|. */
    int64_t g = 0;
    for (size_t i = 0; i < n; i++) {
        int64_t c = coefs[i];
        if (c == INT64_MIN) { g = 1; break; }
        if (c < 0) c = -c;
        g = (g == 0) ? c : gcd(g, c);
    }
    if (g == 0) g = 1;

    /* Polarity: pick the first non-constant term's coefficient sign so the
     * leading variable term ends up positive after dividing through. This
     * turns `-2 x == 4` into `x == -2` rather than `0 == x + 2`. For strict
     * inequalities (Less, Greater) a negative divisor flips the operator;
     * the non-strict and equality forms are direction-symmetric. */
    int64_t divisor = g;
    bool flipped = false;
    for (size_t i = 0; i < n; i++) {
        bool is_const = (rests[i]->type == EXPR_INTEGER &&
                         rests[i]->data.integer == 1);
        if (!is_const) {
            if (coefs[i] < 0) { divisor = -g; flipped = true; }
            break;
        }
    }
    for (size_t i = 0; i < n; i++) coefs[i] /= divisor;

    const char* out_head = hn;
    if (flipped) {
        if      (hn == SYM_Less)         out_head = "Greater";
        else if (hn == SYM_Greater)      out_head = "Less";
        else if (hn == SYM_LessEqual)    out_head = "GreaterEqual";
        else if (hn == SYM_GreaterEqual) out_head = "LessEqual";
    }

    /* Build LHS from positive-coef variable terms, RHS from
     * negated-negative-coef variable terms plus the negated constant. */
    Expr** pos = (Expr**)calloc(n, sizeof(Expr*));
    Expr** neg = (Expr**)calloc(n, sizeof(Expr*));
    size_t pn = 0, nn = 0;
    int64_t const_sum = 0;       /* moves to RHS as -const_sum */
    bool const_overflow = false; /* on overflow, fall back to a Plus term */
    Expr** const_terms = (Expr**)calloc(n, sizeof(Expr*));
    size_t cn = 0;
    for (size_t i = 0; i < n; i++) {
        bool is_const = (rests[i]->type == EXPR_INTEGER &&
                         rests[i]->data.integer == 1);
        if (is_const) {
            int64_t c = coefs[i];
            /* Track sum but guard against int64 overflow. */
            int64_t sum;
            if (!const_overflow &&
                ((c > 0 && const_sum > INT64_MAX - c) ||
                 (c < 0 && const_sum < INT64_MIN - c))) {
                const_overflow = true;
            }
            if (!const_overflow) {
                sum = const_sum + c;
                const_sum = sum;
            }
            /* Always keep the term in case we hit overflow later. */
            const_terms[cn++] = simp_make_term(c, rests[i]);
        } else {
            if (coefs[i] > 0) {
                pos[pn++] = simp_make_term(coefs[i], rests[i]);
            } else if (coefs[i] < 0) {
                neg[nn++] = simp_make_term(-coefs[i], rests[i]);
            } else {
                expr_free(rests[i]);
            }
        }
    }

    Expr* new_lhs;
    if (pn == 0)      new_lhs = expr_new_integer(0);
    else if (pn == 1) new_lhs = pos[0];
    else              new_lhs = expr_new_function(
                          expr_new_symbol(SYM_Plus), pos, pn);

    /* RHS = (negated negative-coef vars) + (-const). */
    size_t total_rhs = nn + cn;
    Expr* new_rhs;
    if (total_rhs == 0) {
        new_rhs = expr_new_integer(0);
        for (size_t i = 0; i < cn; i++) expr_free(const_terms[i]);
    } else {
        Expr** rhs_terms = (Expr**)calloc(total_rhs, sizeof(Expr*));
        size_t rt = 0;
        for (size_t i = 0; i < nn; i++) rhs_terms[rt++] = neg[i];
        if (!const_overflow) {
            /* Single integer for the constant: -const_sum (zero is fine). */
            for (size_t i = 0; i < cn; i++) expr_free(const_terms[i]);
            if (const_sum != 0 || rt == 0) {
                /* Build -const_sum, watching INT64_MIN. */
                int64_t neg_const = (const_sum == INT64_MIN)
                                        ? INT64_MAX  /* impossible in practice */
                                        : -const_sum;
                rhs_terms[rt++] = expr_new_integer(neg_const);
            }
        } else {
            /* Overflow path: keep each constant term, negated. */
            for (size_t i = 0; i < cn; i++) {
                /* Negate the leading coefficient. */
                if (const_terms[i]->type == EXPR_INTEGER) {
                    /* Replace, don't mutate: the integer atom may be
                     * shared (M3 atom-sharing). */
                    int64_t v = -const_terms[i]->data.integer;
                    expr_free(const_terms[i]);
                    rhs_terms[rt++] = expr_new_integer(v);
                } else {
                    /* Wrap in Times[-1, ...]. */
                    Expr* args[2] = { expr_new_integer(-1), const_terms[i] };
                    rhs_terms[rt++] = expr_new_function(
                        expr_new_symbol(SYM_Times), args, 2);
                }
            }
        }
        if (rt == 0) {
            new_rhs = expr_new_integer(0);
            free(rhs_terms);
        } else if (rt == 1) {
            new_rhs = rhs_terms[0];
            free(rhs_terms);
        } else {
            new_rhs = expr_new_function(
                expr_new_symbol(SYM_Plus), rhs_terms, rt);
            free(rhs_terms);
        }
    }

    free(const_terms);
    free(pos);
    free(neg);
    free(coefs);
    free(rests);
    expr_free(d);

    /* Re-evaluate each side so canonical ordering / Plus flattening kicks in. */
    Expr* lhs_e = eval_and_free(new_lhs);
    Expr* rhs_e = eval_and_free(new_rhs);

    Expr* rel_args[2] = { lhs_e, rhs_e };
    Expr* out = expr_new_function(
        expr_new_symbol(out_head), rel_args, 2);
    Expr* out_eval = eval_and_free(out);
    return out_eval;
}

/* Thread Simplify over the elements of a List, deciding numeric contagion
 * independently per element. Numeric inexactness must not cross a List
 * boundary: an exact element that merely shares a list with an inexact
 * sibling stays exact (WL: `{1., 1} // Simplify` === `{1., 1}`). The blanket
 * rationalise-then-numericalise path (see builtin_simplify below) would
 * otherwise numericalise the whole result tree, turning every exact number
 * in the list inexact. Trailing options / assumptions are preserved on each
 * per-element call; nested lists (e.g. matrices) are handled by the
 * recursive re-entry through the evaluator. */
static Expr* simplify_thread_list(Expr* res) {
    Expr*  list = res->data.function.args[0];
    size_t n    = list->data.function.arg_count;
    size_t argc = res->data.function.arg_count;

    Expr** out = (n > 0) ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) {
        /* Build Simplify[el_i, opts...] and evaluate it. Re-entry lets each
         * element take the inexact-or-exact branch on its own merits. */
        Expr** call_args = (Expr**)malloc(argc * sizeof(Expr*));
        call_args[0] = expr_copy(list->data.function.args[i]);
        for (size_t j = 1; j < argc; j++)
            call_args[j] = expr_copy(res->data.function.args[j]);
        Expr* call = expr_new_function(expr_new_symbol(SYM_Simplify),
                                       call_args, argc);
        free(call_args);
        out[i] = eval_and_free(call);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    if (out) free(out);
    return result;
}

/* True iff some Power[base, Rational[p, q]] (q >= 2) in `e` has a *compound*
 * polynomial base -- base is a function whose head is not Rational/Complex
 * (e.g. Sqrt[6 + x^2]). That is a radicand the radical-fraction polish's
 * Factor call can actually expose. It deliberately EXCLUDES bare-symbol or
 * integer radicands like Sqrt[u] or Sqrt[6]: Factor cannot help there, and
 * over algebraically-dependent generators such as Sqrt[u] and u it would
 * hand poly_gcd_internal a degenerate pseudo-remainder sequence. Paired with
 * simp_has_rational_root as the polish gate below. */
static bool has_compound_radicand(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Power
        && argc == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION && exp->data.function.head
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol.name == SYM_Rational
            && exp->data.function.arg_count == 2) {
            const Expr* qq = exp->data.function.args[1];
            if (qq->type == EXPR_INTEGER && qq->data.integer >= 2
                && base->type == EXPR_FUNCTION && base->data.function.head
                && base->data.function.head->type == EXPR_SYMBOL
                && base->data.function.head->data.symbol.name != SYM_Rational
                && base->data.function.head->data.symbol.name != SYM_Complex) {
                return true;
            }
        }
    }
    for (size_t i = 0; i < argc; i++) {
        if (has_compound_radicand(e->data.function.args[i])) return true;
    }
    return false;
}

/* Parse a TimeConstraint option value into a per-subexpression budget in
 * seconds. Accepts a machine-real number (Integer/Real/BigInt/Rational),
 * Infinity / DirectedInfinity[1] (-> no limit), and a list {tLoc, ...} whose
 * first element is the per-subexpression budget (matching FullSimplify's tLoc).
 * Anything else, or a non-positive value, yields HUGE_VAL == "no limit" so the
 * default path stays inert. The WHOLE-CALL cap is parsed separately by
 * simp_parse_total_budget. */
static double simp_parse_time_budget(const Expr* e) {
    if (!e) return HUGE_VAL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity)
        return HUGE_VAL;
    if (e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_DirectedInfinity) return HUGE_VAL;   /* +Infinity */
        if (h == SYM_List && e->data.function.arg_count >= 1)
            return simp_parse_time_budget(e->data.function.args[0]);
    }
    double v;
    if (common_machine_real_value(e, &v) && v > 0.0) return v;
    return HUGE_VAL;
}

/* Parse a TimeConstraint option value into a WHOLE-CALL budget in seconds
 * (HUGE_VAL == no cap). This bounds the entire Simplify call, closing the gap
 * where a user-set budget was previously ignored on the SHAPE_RATIONAL
 * dispatch path and during bottom-up descent.
 *   scalar t          -> t          (the intuitive "bound the call to ~t")
 *   {tLoc, tTot, ...}  -> tTot       (power-user form; matches FullSimplify)
 *   {t}                -> t          (single-element list, treated like scalar)
 *   Infinity           -> HUGE_VAL   (no cap, the default)
 * Keeping the scalar form as a whole-call cap is the Mathematica-aligned
 * reading and what makes "guarantee termination" hold for a user-set budget. */
static double simp_parse_total_budget(const Expr* e) {
    if (!e) return HUGE_VAL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity)
        return HUGE_VAL;
    if (e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_DirectedInfinity) return HUGE_VAL;   /* +Infinity */
        if (h == SYM_List && e->data.function.arg_count >= 1) {
            size_t idx = e->data.function.arg_count >= 2 ? 1 : 0;
            return simp_parse_total_budget(e->data.function.args[idx]);
        }
    }
    double v;
    if (common_machine_real_value(e, &v) && v > 0.0) return v;
    return HUGE_VAL;
}

/* ----------------------------------------------------------------------- */
/* Inverse-trig rational-angle combination zero test (Machin-like)          */
/* ----------------------------------------------------------------------- */

typedef enum { IA_ATAN, IA_ASIN, IA_ACOS } IAKind;

/* Build e^{i * invtrig[a]} as an UNEVALUATED constant-algebraic Expr (so the
 * Complex[0,1] and Sqrt structure survives for to_qqbar):
 *   e^{i ArcTan a} = (1 + i a) / Sqrt[1 + a^2]
 *   e^{i ArcSin a} = Sqrt[1 - a^2] + i a
 *   e^{i ArcCos a} = a + i Sqrt[1 - a^2]
 * matching the principal branches (ArcTan,ArcSin in (-pi/2,pi/2]; ArcCos in
 * [0,pi]), so the exp product encodes the true signed angle sum. */
static Expr* ia_imag_unit(void) {
    return expr_new_function(expr_new_symbol(SYM_Complex),
             (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
}
static Expr* ia_sqrt_1pm_a2(const Expr* a, int plus) {
    Expr* a2 = expr_new_function(expr_new_symbol(SYM_Power),
                 (Expr*[]){ expr_copy((Expr*)a), expr_new_integer(2) }, 2);
    Expr* term = plus ? a2
        : expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), a2 }, 2);
    Expr* inside = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_new_integer(1), term }, 2);
    return expr_new_function(expr_new_symbol(SYM_Sqrt), (Expr*[]){ inside }, 1);
}
static Expr* ia_unit_expr(IAKind k, const Expr* a) {
    if (k == IA_ATAN) {
        Expr* num = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_new_integer(1),
                        expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ ia_imag_unit(), expr_copy((Expr*)a) }, 2) }, 2);
        Expr* invsq = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ ia_sqrt_1pm_a2(a, 1), expr_new_integer(-1) }, 2);
        return expr_new_function(expr_new_symbol(SYM_Times),
                 (Expr*[]){ num, invsq }, 2);
    }
    if (k == IA_ASIN) {
        return expr_new_function(expr_new_symbol(SYM_Plus),
                 (Expr*[]){ ia_sqrt_1pm_a2(a, 0),
                   expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ ia_imag_unit(), expr_copy((Expr*)a) }, 2) }, 2);
    }
    /* IA_ACOS */
    return expr_new_function(expr_new_symbol(SYM_Plus),
             (Expr*[]){ expr_copy((Expr*)a),
               expr_new_function(expr_new_symbol(SYM_Times),
                 (Expr*[]){ ia_imag_unit(), ia_sqrt_1pm_a2(a, 0) }, 2) }, 2);
}

/* Split `term` into a numeric coefficient and the single non-numeric atom.
 * Returns 0 if the term has no atom or more than one. coeff is owned. */
static int ia_term_coeff_atom(const Expr* term, Expr** coeff, const Expr** atom) {
    if (term->type == EXPR_FUNCTION && term->data.function.head->type == EXPR_SYMBOL &&
        term->data.function.head->data.symbol.name == SYM_Times) {
        size_t n = term->data.function.arg_count;
        Expr* c = expr_new_integer(1);
        const Expr* at = NULL;
        for (size_t i = 0; i < n; i++) {
            const Expr* f = term->data.function.args[i];
            if (f->type == EXPR_INTEGER || f->type == EXPR_BIGINT ||
                (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
                 f->data.function.head->data.symbol.name == SYM_Rational)) {
                Expr* nc = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ c, expr_copy((Expr*)f) }, 2);
                c = eval_and_free(nc);
            } else if (!at) {
                at = f;
            } else { expr_free(c); return 0; }   /* two non-numeric factors */
        }
        if (!at) { expr_free(c); return 0; }
        *coeff = c; *atom = at; return 1;
    }
    *coeff = expr_new_integer(1); *atom = term; return 1;
}

static int ia_is_integer(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

/* Returns 1 iff `e` is a Z-linear combination of ArcTan/ArcSin/ArcCos of
 * constant arguments plus a rational multiple of Pi that is identically 0.
 * EXACT certificate: e^{i e} == 1 (via the qqbar constant-algebraic engine),
 * which proves e == 0 (mod 2 Pi); a numeric screen then selects the branch
 * k = round(e / 2 Pi) and the identity holds only when k == 0.  The numeric
 * step merely disambiguates the integer multiple (a discrete choice with a
 * ~pi margin, far above double precision) -- the zero itself is certified
 * algebraically, never by sampling. */
static int simp_invtrig_combo_is_zero(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Plus) return 0;
    size_t nt = e->data.function.arg_count;

    /* W = product of e^{i * coeff * invtrig} times e^{i * piCoeff * Pi}. */
    Expr* W = expr_new_integer(1);
    Expr* piAccum = NULL;            /* sum of Pi coefficients (rational) */
    int have_invtrig = 0, ok = 1;

    for (size_t i = 0; i < nt && ok; i++) {
        Expr* coeff = NULL; const Expr* atom = NULL;
        if (!ia_term_coeff_atom(e->data.function.args[i], &coeff, &atom)) { ok = 0; break; }
        if (atom->type == EXPR_SYMBOL && atom->data.symbol.name == SYM_Pi) {
            piAccum = piAccum ? expr_new_function(expr_new_symbol(SYM_Plus),
                                  (Expr*[]){ piAccum, coeff }, 2)
                              : coeff;
            continue;
        }
        IAKind k;
        if (atom->type == EXPR_FUNCTION && atom->data.function.arg_count == 1 &&
            atom->data.function.head->type == EXPR_SYMBOL) {
            const char* h = atom->data.function.head->data.symbol.name;
            if      (h == SYM_ArcTan) k = IA_ATAN;
            else if (h == SYM_ArcSin) k = IA_ASIN;
            else if (h == SYM_ArcCos) k = IA_ACOS;
            else { expr_free(coeff); ok = 0; break; }
        } else { expr_free(coeff); ok = 0; break; }
        if (!ia_is_integer(coeff)) { expr_free(coeff); ok = 0; break; }  /* integer coeff only */
        have_invtrig = 1;
        Expr* unit = ia_unit_expr(k, atom->data.function.args[0]);
        Expr* powu = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ unit, coeff }, 2);     /* consumes coeff */
        W = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ W, powu }, 2);
    }

    if (ok && piAccum) {
        /* e^{i piCoeff Pi} = Power[E, Times[Complex[0, piCoeff], Pi]] */
        Expr* cplx = expr_new_function(expr_new_symbol(SYM_Complex),
                       (Expr*[]){ expr_new_integer(0), piAccum }, 2);  /* consumes piAccum */
        Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ cplx, expr_new_symbol(SYM_Pi) }, 2);
        Expr* epi = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_new_symbol(SYM_E), arg }, 2);
        W = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ W, epi }, 2);
        piAccum = NULL;
    } else if (piAccum) {
        expr_free(piAccum); piAccum = NULL;
    }

    if (!ok || !have_invtrig) { expr_free(W); return 0; }

    /* Exact: reduce W over the algebraic numbers; identity needs W == 1. */
    int proven = 0;
    Expr* Wval = flint_qqbar_canonical(W, QQBAR_METHOD_AUTOMATIC);
    expr_free(W);
    if (Wval) {
        int is_one = (Wval->type == EXPR_INTEGER && Wval->data.integer == 1);
        expr_free(Wval);
        if (is_one) {
            /* Branch selection: e == 2 Pi k; keep only k == 0.  N[e] has error
             * ~1e-15 << pi, so round(N[e]/2pi) is the exact k. */
            Expr* nv = eval_and_free(expr_new_function(expr_new_symbol("N"),
                         (Expr*[]){ expr_copy((Expr*)e) }, 1));
            double v;
            if (nv && common_machine_real_value(nv, &v)) {
                double twopi = 2.0 * acos(-1.0);
                double kf = v / twopi;
                long k = (long)(kf < 0 ? kf - 0.5 : kf + 0.5);   /* round */
                if (k == 0 && (v > -0.5 && v < 0.5)) proven = 1;
            }
            if (nv) expr_free(nv);
        }
    }
    return proven;
}

typedef enum { IH_ASINH, IH_ACOSH, IH_ATANH } IHKind;

/* Build e^{ArcHyp[a]} as an UNEVALUATED constant-algebraic Expr (real, for real
 * a in the function's domain):
 *   e^{ArcSinh a} = a + Sqrt[1 + a^2]
 *   e^{ArcCosh a} = a + Sqrt[a^2 - 1]        (a >= 1)
 *   e^{ArcTanh a} = Sqrt[1 + a] / Sqrt[1 - a]  (|a| < 1)
 * Since exp is strictly monotonic on R, a real combination c = sum ck ArcHyp[ak]
 * is 0 iff e^c = 1 -- NO branch screen needed (unlike the circular case). */
static Expr* ih_unit_expr(IHKind k, const Expr* a) {
    if (k == IH_ASINH)
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy((Expr*)a), ia_sqrt_1pm_a2(a, 1) }, 2));
    if (k == IH_ACOSH) {
        Expr* a2m1 = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy((Expr*)a), expr_new_integer(2) }, 2),
                       expr_new_integer(-1) }, 2);
        Expr* sq = expr_new_function(expr_new_symbol(SYM_Sqrt), (Expr*[]){ a2m1 }, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy((Expr*)a), sq }, 2));
    }
    /* IH_ATANH: Sqrt[1+a] * Sqrt[1-a]^(-1) */
    Expr* num = expr_new_function(expr_new_symbol(SYM_Sqrt),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_integer(1), expr_copy((Expr*)a) }, 2) }, 1);
    Expr* den = expr_new_function(expr_new_symbol(SYM_Sqrt),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_integer(1),
                      expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)a) }, 2) }, 2) }, 1);
    Expr* deninv = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ den, expr_new_integer(-1) }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ num, deninv }, 2));
}

/* Numeric domain screen: the argument `a` (a constant) lies in ArcHyp's real
 * domain (ArcCosh: a>=1; ArcTanh: -1<a<1; ArcSinh: all).  A discrete check with
 * margin -- it only confirms we are on the real branch (so e^c=1 <=> c=0); the
 * zero itself is certified algebraically. */
static int ih_arg_in_domain(IHKind k, const Expr* a) {
    if (k == IH_ASINH) return 1;
    Expr* nv = eval_and_free(expr_new_function(expr_new_symbol("N"),
                   (Expr*[]){ expr_copy((Expr*)a) }, 1));
    double v; int ok = 0;
    if (nv && common_machine_real_value(nv, &v)) {
        if (k == IH_ACOSH) ok = (v >= 1.0 - 1e-9);
        else /* IH_ATANH */ ok = (v > -1.0 + 1e-9 && v < 1.0 - 1e-9);
    }
    if (nv) expr_free(nv);
    return ok;
}

/* Returns 1 iff `e` is a Z-combination of ArcSinh/ArcCosh/ArcTanh of constant
 * (real, in-domain) arguments that is identically 0.  EXACT: e^e == 1 over the
 * algebraic numbers (exp injective on R, so no branch ambiguity). */
static int simp_invhyp_combo_is_zero(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Plus) return 0;
    size_t nt = e->data.function.arg_count;
    Expr* W = expr_new_integer(1);
    int have = 0, ok = 1;
    for (size_t i = 0; i < nt && ok; i++) {
        Expr* coeff = NULL; const Expr* atom = NULL;
        if (!ia_term_coeff_atom(e->data.function.args[i], &coeff, &atom)) { ok = 0; break; }
        IHKind k;
        if (atom->type == EXPR_FUNCTION && atom->data.function.arg_count == 1 &&
            atom->data.function.head->type == EXPR_SYMBOL) {
            const char* h = atom->data.function.head->data.symbol.name;
            if      (strcmp(h, "ArcSinh") == 0) k = IH_ASINH;
            else if (strcmp(h, "ArcCosh") == 0) k = IH_ACOSH;
            else if (strcmp(h, "ArcTanh") == 0) k = IH_ATANH;
            else { expr_free(coeff); ok = 0; break; }
        } else { expr_free(coeff); ok = 0; break; }
        if (!ia_is_integer(coeff) || !ih_arg_in_domain(k, atom->data.function.args[0])) {
            expr_free(coeff); ok = 0; break;
        }
        have = 1;
        Expr* unit = ih_unit_expr(k, atom->data.function.args[0]);
        Expr* powu = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ unit, coeff }, 2);
        W = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ W, powu }, 2);
    }
    if (!ok || !have) { expr_free(W); return 0; }
    int proven = 0;
    Expr* Wval = flint_qqbar_canonical(W, QQBAR_METHOD_AUTOMATIC);
    expr_free(W);
    if (Wval) {
        proven = (Wval->type == EXPR_INTEGER && Wval->data.integer == 1);
        expr_free(Wval);
    }
    return proven;
}

/* Dispatcher: a constant inverse-function combination is zero (circular via
 * the e^{i e}==1 certifier, hyperbolic via e^e==1). */
static int simp_invfn_combo_is_zero(const Expr* e) {
    return simp_invtrig_combo_is_zero(e) || simp_invhyp_combo_is_zero(e);
}

/* ----------------------------------------------------------------------- */
/* Phase 2: region-identity engine (derivative-constancy on a box)          */
/* ----------------------------------------------------------------------- */
/*
 * Proves f(x1..xn) == 0 on the assumption region R, for f a combination of
 * inverse trig/hyperbolic functions (the addition identities), by the theorem:
 * if R is connected, f is continuous on R, grad f == 0 on R, and f(p)=0 for one
 * p in R, then f == 0 on R.  Each hypothesis is discharged SOUNDLY:
 *   - R connected: require R to be a BOX (every assumption constrains a single
 *     variable) -- convex, hence connected;
 *   - grad f == 0: each partial D[f,xi] is proven 0 by the exact Simplify zero
 *     machinery (no sampling);
 *   - f continuous on R: each inverse-function subterm's argument is verified by
 *     Reduce/CAD entailment to stay in that function's real continuity domain
 *     throughout R, so no branch cut is crossed;
 *   - f(p)=0: a sample point p in R (FindInstance) yields a CONSTANT combination
 *     certified 0 exactly by simp_invfn_combo_is_zero.
 * Declines (0) on any failure; the result is never certified by sampling. */

static const char* region_inv_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1 ||
        !e->data.function.head || e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol.name;
    if (strcmp(h,"ArcTan")==0||strcmp(h,"ArcSin")==0||strcmp(h,"ArcCos")==0||
        strcmp(h,"ArcCot")==0||strcmp(h,"ArcSinh")==0||strcmp(h,"ArcCosh")==0||
        strcmp(h,"ArcTanh")==0||strcmp(h,"ArcCoth")==0) return h;
    return NULL;
}
static void region_vars(const Expr* e, const Expr*** v, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol.name)) return;
        for (size_t i=0;i<*n;i++) if ((*v)[i]->data.symbol.name==e->data.symbol.name) return;
        if (*n==*cap){*cap=*cap?*cap*2:8;*v=realloc((void*)*v,*cap*sizeof(Expr*));}
        (*v)[(*n)++]=e; return;
    }
    if (e->type==EXPR_FUNCTION)
        for (size_t i=0;i<e->data.function.arg_count;i++)
            region_vars(e->data.function.args[i],v,n,cap);
}
static size_t region_count_vars_in(const Expr* e, const Expr** vars, size_t nv) {
    const Expr** seen=NULL; size_t ns=0, cap=0;
    region_vars(e,&seen,&ns,&cap);
    size_t c=0;
    for (size_t i=0;i<nv;i++)
        for (size_t j=0;j<ns;j++)
            if (vars[i]->data.symbol.name==seen[j]->data.symbol.name){c++;break;}
    free((void*)seen);
    return c;
}
/* Is `e` a branch-sensitive node (an inverse trig/hyp, Log, Sqrt, or an
 * even-root Power) whose continuity on R must be verified? */
static bool region_branch_node(const Expr* e) {
    if (region_inv_head(e)) return true;
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if ((strcmp(h,"Log")==0 || strcmp(h,"Sqrt")==0) && n == 1) return true;
    if (strcmp(h,"Power")==0 && n == 2 && is_rational_literal(e->data.function.args[1])) {
        const Expr* q = e->data.function.args[1]->data.function.args[1];
        if (q->type == EXPR_INTEGER && q->data.integer != 1 && (q->data.integer % 2) == 0)
            return true;   /* even root: real only for nonneg base */
    }
    return false;
}
static void region_collect_branch(const Expr* e, const Expr*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (region_branch_node(e)) {
        if (*n==*cap){*cap=*cap?*cap*2:8;*out=realloc((void*)*out,*cap*sizeof(Expr*));}
        (*out)[(*n)++]=e;
    }
    if (e->type==EXPR_FUNCTION)
        for (size_t i=0;i<e->data.function.arg_count;i++)
            region_collect_branch(e->data.function.args[i],out,n,cap);
}
/* Continuity-domain predicate for a branch-sensitive node (owned), or NULL when
 * it is continuous for all real inputs with no finiteness issue. */
static Expr* region_domain_cond(const Expr* node) {
    const char* head = node->data.function.head->data.symbol.name;
    const Expr* g = node->data.function.args[0];
    if (strcmp(head,"ArcTanh")==0 || strcmp(head,"ArcSin")==0 || strcmp(head,"ArcCos")==0) {
        Expr* lo = expr_new_function(expr_new_symbol("Less"),
            (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)g) }, 2);
        Expr* hi = expr_new_function(expr_new_symbol("Less"),
            (Expr*[]){ expr_copy((Expr*)g), expr_new_integer(1) }, 2);
        return expr_new_function(expr_new_symbol("And"), (Expr*[]){ lo, hi }, 2);
    }
    if (strcmp(head,"ArcCosh")==0)
        return expr_new_function(expr_new_symbol("Greater"),
            (Expr*[]){ expr_copy((Expr*)g), expr_new_integer(1) }, 2);
    /* Log[g] and Sqrt[g] / even-root Power[g,_]: real-continuous for g > 0. */
    if (strcmp(head,"Log")==0 || strcmp(head,"Sqrt")==0 || strcmp(head,"Power")==0)
        return expr_new_function(expr_new_symbol("Greater"),
            (Expr*[]){ expr_copy((Expr*)g), expr_new_integer(0) }, 2);
    /* ArcTan/ArcSinh/ArcCot/ArcCoth: continuous for all real g except poles of g;
     * require Denominator[g] != 0 on R (trivial when g is polynomial). */
    Expr* den = eval_and_free(expr_new_function(expr_new_symbol("Denominator"),
                    (Expr*[]){ expr_copy((Expr*)g) }, 1));
    if (den && (den->type == EXPR_INTEGER || den->type == EXPR_BIGINT)) { expr_free(den); return NULL; }
    return expr_new_function(expr_new_symbol("Unequal"),
        (Expr*[]){ den, expr_new_integer(0) }, 2);
}
static Expr* region_assum_conj(const AssumeCtx* ctx) {
    if (!ctx || ctx->count==0) return expr_new_symbol("True");
    Expr** k = calloc(ctx->count, sizeof(Expr*)); size_t nk=0;
    for (size_t i=0;i<ctx->count;i++) k[nk++]=expr_copy(ctx->facts[i]);
    Expr* out = (nk==1)?k[0]:expr_new_function(expr_new_symbol("And"),k,nk);
    free(k); return out;
}

/* The engine fires only on a Plus whose TOP-LEVEL terms are each an inverse
 * function, a Log (Cap A may convert ArcSinh->Log), or a variable-free constant
 * -- the clean shape of an addition identity.  This both tightens the gate and,
 * crucially, keeps the recursive Simplify calls (gradient, sample value) from
 * re-entering: a gradient (rational), a constant (var-free), or a forward-of-
 * inverse residual like Cosh[ArcCosh x + ...] all fail this shape, so no
 * re-entrancy guard is needed (and none can leak across a TimeConstrained abort). */
static bool region_shape_ok(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Plus) return false;
    bool has_inv = false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* coeff = NULL; const Expr* atom = NULL;
        if (!ia_term_coeff_atom(e->data.function.args[i], &coeff, &atom)) return false;
        expr_free(coeff);
        if (region_inv_head(atom)) { has_inv = true; continue; }
        if (atom->type == EXPR_FUNCTION && atom->data.function.head->type == EXPR_SYMBOL &&
            strcmp(atom->data.function.head->data.symbol.name, "Log") == 0 &&
            atom->data.function.arg_count == 1) continue;
        const Expr** sv = NULL; size_t ns = 0, cap = 0;
        region_vars(atom, &sv, &ns, &cap);
        free((void*)sv);
        if (ns == 0) continue;            /* a variable-free constant term */
        return false;                      /* anything else (Cosh[..], poly, ...) */
    }
    return has_inv;
}

static int simp_region_identity_is_zero(const Expr* e, const AssumeCtx* ctx) {
    if (!ctx || ctx->count == 0 || ctx->inconsistent) return 0;
    if (!region_shape_ok(e)) return 0;

    const Expr** vars=NULL; size_t nv=0, vcap=0;
    region_vars(e,&vars,&nv,&vcap);
    if (nv == 0 || nv > 4) { free((void*)vars); return 0; }
    for (size_t i=0;i<nv;i++) if (!prov_re(ctx, vars[i])) { free((void*)vars); return 0; }

    /* Performance gate: the derivative of ArcSin/ArcCos/ArcSinh/ArcCosh is
     * algebraic, and its multi-variable zero test (Sqrt of a 2-var expression)
     * is not provable by the exact machinery and is costly -- so for >= 2
     * variables with such a head present, decline here rather than grind on the
     * gradient.  (Those addition identities are closed by the pointwise
     * hyperbolic recognizer, or genuinely need a 2-var multi-radical test.)
     * ArcTan/ArcTanh (rational derivative) and the 1-variable case are kept. */
    if (nv >= 2) {
        const Expr** br=NULL; size_t nb=0, bc=0;
        region_collect_branch(e,&br,&nb,&bc);
        int algebraic = 0;
        for (size_t i=0;i<nb;i++) {
            const char* h = region_inv_head(br[i]);
            if (h && (strcmp(h,"ArcSin")==0||strcmp(h,"ArcCos")==0||
                      strcmp(h,"ArcSinh")==0||strcmp(h,"ArcCosh")==0)) { algebraic = 1; break; }
        }
        free((void*)br);
        if (algebraic) { free((void*)vars); return 0; }
    }

    /* Region must be a BOX: every relational fact constrains <= 1 of our vars. */
    for (size_t i=0;i<ctx->count;i++)
        if (ctx->facts[i]->type==EXPR_FUNCTION &&
            region_count_vars_in(ctx->facts[i], vars, nv) >= 2) { free((void*)vars); return 0; }

    /* Gradient: each partial D[e,xi] must prove to literal 0 under R. */
    Expr* conj = region_assum_conj(ctx);
    int grad_ok = 1;
    for (size_t i=0;i<nv && grad_ok;i++) {
        Expr* d = expr_new_function(expr_new_symbol(SYM_D),
                    (Expr*[]){ expr_copy((Expr*)e), expr_copy((Expr*)vars[i]) }, 2);
        Expr* simp = eval_and_free(expr_new_function(expr_new_symbol("Simplify"),
                        (Expr*[]){ d, expr_copy(conj) }, 2));
        if (!(simp && simp->type==EXPR_INTEGER && simp->data.integer==0)) grad_ok = 0;
        if (simp) expr_free(simp);
    }
    if (!grad_ok) { expr_free(conj); free((void*)vars); return 0; }

    /* Continuity: each branch-sensitive argument (inverse fn, Log, Sqrt, even
     * root) stays in its real continuity domain throughout R, so e has no branch
     * cut on the (connected) box. */
    const Expr** br=NULL; size_t ni=0, icap=0;
    region_collect_branch(e,&br,&ni,&icap);
    int cont_ok = 1;
    for (size_t i=0;i<ni && cont_ok;i++) {
        Expr* cond = region_domain_cond(br[i]);
        if (cond) { if (assume_reduce_entails(ctx, cond) != 1) cont_ok = 0; expr_free(cond); }
    }
    free((void*)br);
    if (!cont_ok) { expr_free(conj); free((void*)vars); return 0; }

    /* Sample point in R via FindInstance; certify the constant value is 0. */
    Expr** vlist = malloc(nv*sizeof(Expr*));
    for (size_t i=0;i<nv;i++) vlist[i]=expr_copy((Expr*)vars[i]);
    free((void*)vars);
    Expr* vl = expr_new_function(expr_new_symbol("List"), vlist, nv);
    free(vlist);
    Expr* fi = eval_and_free(expr_new_function(expr_new_symbol("FindInstance"),
                   (Expr*[]){ conj, vl, expr_new_symbol("Reals") }, 3));
    int proven = 0;
    if (fi && fi->type==EXPR_FUNCTION && fi->data.function.head->type==EXPR_SYMBOL &&
        fi->data.function.head->data.symbol.name==SYM_List &&
        fi->data.function.arg_count >= 1) {
        Expr* ep = eval_and_free(expr_new_function(expr_new_symbol("ReplaceAll"),
                       (Expr*[]){ expr_copy((Expr*)e),
                                  expr_copy(fi->data.function.args[0]) }, 2));
        if (ep) {
            /* ep is a CONSTANT combination: 0 directly, via the exact inverse-fn
             * certifier, or (for Log-converted forms like ArcCosh[2]-Log[2+Sqrt3])
             * via a recursive Simplify of the constant.  All exact, no sampling. */
            if (ep->type==EXPR_INTEGER && ep->data.integer==0) proven = 1;
            else if (simp_invfn_combo_is_zero(ep)) proven = 1;
            else if (ep->type == EXPR_FUNCTION) {
                Expr* s = eval_and_free(expr_new_function(expr_new_symbol("Simplify"),
                              (Expr*[]){ expr_copy(ep) }, 1));
                if (s && s->type==EXPR_INTEGER && s->data.integer==0) proven = 1;
                if (s) expr_free(s);
            }
            expr_free(ep);
        }
    }
    if (fi) expr_free(fi);
    return proven;
}

/* Inverse-HYPERBOLIC addition identities on a region: ArcSinh[a]+ArcSinh[b] =
 * ArcSinh[c], and likewise ArcCosh / ArcTanh.  Unlike the circular Sin/Cos case
 * (2-to-1, needing a radical range condition Reduce cannot discharge), the
 * hyperbolic forward functions are injective enough that the identity is
 * POINTWISE once the INPUT arguments are in domain -- the target's domain
 * follows automatically from c = F[theta] (Cosh>=1, Tanh in (-1,1)), and
 * ArcF[F[theta]] = theta holds (no branch reflection: ArcSinh/ArcTanh are full
 * inverses on R; ArcCosh[Cosh[theta]]=theta for theta>=0, and theta =
 * ArcCosh[a]+ArcCosh[b] >= 0 automatically).  So NO connectedness, sample point,
 * or radical range check is needed: verify F[theta] = c (algebraic) and the two
 * input args lie in F's domain (a POLYNOMIAL condition Reduce handles). */
static int simp_invhyp_addition_is_zero(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Plus ||
        e->data.function.arg_count != 3) return 0;
    if (!ctx || ctx->count == 0 || ctx->inconsistent) return 0;

    const char* arc = NULL; const Expr* args[3]; int sg[3];
    for (size_t i = 0; i < 3; i++) {
        Expr* coeff = NULL; const Expr* atom = NULL;
        if (!ia_term_coeff_atom(e->data.function.args[i], &coeff, &atom)) return 0;
        int s = 0;
        if (coeff->type == EXPR_INTEGER && coeff->data.integer == 1) s = 1;
        else if (coeff->type == EXPR_INTEGER && coeff->data.integer == -1) s = -1;
        expr_free(coeff);
        if (!s) return 0;
        if (atom->type != EXPR_FUNCTION || atom->data.function.arg_count != 1 ||
            atom->data.function.head->type != EXPR_SYMBOL) return 0;
        const char* h = atom->data.function.head->data.symbol.name;
        if (strcmp(h,"ArcSinh") && strcmp(h,"ArcCosh") && strcmp(h,"ArcTanh")) return 0;
        if (!arc) arc = h; else if (strcmp(arc, h) != 0) return 0;
        args[i] = atom->data.function.args[0]; sg[i] = s;
    }
    /* Exactly one term has the minority sign; the other two sum to theta. */
    int sum = sg[0]+sg[1]+sg[2];
    int maj = (sum > 0) ? 1 : -1, oddi = -1, nodd = 0;
    for (int i = 0; i < 3; i++) if (sg[i] != maj) { oddi = i; nodd++; }
    if (nodd != 1) return 0;
    int ia = -1, ib = -1;
    for (int i = 0; i < 3; i++) if (i != oddi) { if (ia < 0) ia = i; else ib = i; }

    const char* Fh = (strcmp(arc,"ArcSinh")==0) ? "Sinh"
                   : (strcmp(arc,"ArcCosh")==0) ? "Cosh" : "Tanh";

    /* theta = ArcF[args[ia]] + ArcF[args[ib]]; step 1: F[theta] - c == 0. */
    Expr* th = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
        expr_new_function(expr_new_symbol(arc), (Expr*[]){ expr_copy((Expr*)args[ia]) }, 1),
        expr_new_function(expr_new_symbol(arc), (Expr*[]){ expr_copy((Expr*)args[ib]) }, 1) }, 2);
    /* TrigExpand forces the forward-of-inverse expansion F[ArcF a + ArcF b] ->
     * algebraic (e.g. Cosh[ArcCosh a + ArcCosh b] -> a b + Sqrt[a^2-1]Sqrt[b^2-1]),
     * which the plain pipeline does not always reach under an assumption set. */
    Expr* fexp = eval_and_free(expr_new_function(expr_new_symbol("TrigExpand"),
                     (Expr*[]){ expr_new_function(expr_new_symbol(Fh), (Expr*[]){ th }, 1) }, 1));
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
        fexp,
        expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)args[oddi]) }, 2) }, 2);
    /* Use the assumptions: Cap A may have rewritten the target's radical form
     * (Sqrt[a]Sqrt[b] -> Sqrt[a b]) under them, so matching needs them too. */
    Expr* s1 = eval_and_free(expr_new_function(expr_new_symbol("Simplify"),
                   (Expr*[]){ diff, region_assum_conj(ctx) }, 2));
    int step1 = (s1 && s1->type == EXPR_INTEGER && s1->data.integer == 0);
    if (s1) expr_free(s1);
    if (!step1) return 0;

    /* Domain of the two input args (a POLYNOMIAL condition): ArcCosh arg >= 1,
     * ArcTanh arg in (-1,1); ArcSinh any real (require reality of free vars). */
    for (int t = 0; t < 2; t++) {
        const Expr* a = (t == 0) ? args[ia] : args[ib];
        Expr* cond = NULL;
        if (strcmp(arc,"ArcCosh") == 0)
            cond = expr_new_function(expr_new_symbol("GreaterEqual"),
                       (Expr*[]){ expr_copy((Expr*)a), expr_new_integer(1) }, 2);
        else if (strcmp(arc,"ArcTanh") == 0) {
            Expr* lo = expr_new_function(expr_new_symbol("Less"),
                (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)a) }, 2);
            Expr* hi = expr_new_function(expr_new_symbol("Less"),
                (Expr*[]){ expr_copy((Expr*)a), expr_new_integer(1) }, 2);
            cond = expr_new_function(expr_new_symbol("And"), (Expr*[]){ lo, hi }, 2);
        }
        if (cond) { int ent = assume_reduce_entails(ctx, cond); expr_free(cond); if (ent != 1) return 0; }
        else { /* ArcSinh: require the argument's free symbols real. */
            const Expr** sv=NULL; size_t ns=0, cap=0;
            region_vars(a, &sv, &ns, &cap);
            int allr = (ns > 0);
            for (size_t i=0;i<ns;i++) if (!prov_re(ctx, sv[i])) { allr = 0; break; }
            free((void*)sv);
            if (!allr) return 0;
        }
    }
    return 1;
}

Expr* builtin_simplify(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* The simplification pipeline routes through Together/Cancel/Apart/
     * Factor and the polynomial GCD machinery, all of which need rational
     * coefficients. Rationalise on entry, run the exact pipeline, then
     * numericalise on the way out so callers still see inexact-in /
     * inexact-out semantics. */
    if (internal_args_contain_inexact(res)) {
        /* Numeric contagion must not cross List boundaries: thread over the
         * list so each element's inexactness is judged on its own, instead of
         * numericalising the whole result tree (see simplify_thread_list). */
        Expr* expr0 = res->data.function.args[0];
        if (expr0->type == EXPR_FUNCTION && head_is(expr0, SYM_List)) {
            return simplify_thread_list(res);
        }
        return internal_rationalize_then_numericalize(res, builtin_simplify);
    }

    Expr* expr = res->data.function.args[0];

    /* Parse remaining args: at most one positional assumption, plus
     * options Rule[Assumptions, X], Rule[ComplexityFunction, f], and
     * Rule[TransformationFunctions, spec]. */
    Expr* positional_assum = NULL;
    Expr* opt_assumptions  = NULL;
    Expr* opt_complexity   = NULL;
    Expr* opt_transform    = NULL;
    Expr* opt_timeconstraint = NULL;

    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            opt_assumptions = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "ComplexityFunction")) {
            opt_complexity = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "TransformationFunctions")) {
            opt_transform = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "TimeConstraint")) {
            /* Must precede the positional-assumption fallback: otherwise
             * TimeConstraint -> t is silently swallowed as an assumption. */
            opt_timeconstraint = a->data.function.args[1];
        } else if (positional_assum == NULL) {
            positional_assum = a;
        }
    }

    /* Per-subexpression wall-clock budget (HUGE_VAL == no limit, the default). */
    double time_budget = simp_parse_time_budget(opt_timeconstraint);
    /* Whole-call wall-clock cap (HUGE_VAL == no cap). Bounds the ENTIRE call,
     * so a user-set TimeConstraint is honoured on every path -- including the
     * SHAPE_RATIONAL dispatch route that bypasses simp_search's local window. */
    double time_total = simp_parse_total_budget(opt_timeconstraint);

    /* Resolve TransformationFunctions into (use_builtin, user_funcs[]).
     *   Automatic             -> built-in pipeline only (the default).
     *   {f1, ...}             -> ONLY f1, ... (built-in pipeline suppressed).
     *   {Automatic, f1, ...}  -> built-in pipeline plus f1, ...
     *   bare f                -> treat as the single-function list {f}.
     * The funcs are borrowed pointers into the option expression, which the
     * evaluator keeps alive until after this builtin returns. */
    bool   use_builtin = true;
    Expr** user_funcs  = NULL;
    size_t n_user_funcs = 0;
    if (opt_transform) {
        if (opt_transform->type == EXPR_SYMBOL &&
            opt_transform->data.symbol.name == SYM_Automatic) {
            /* Default: built-in only. */
        } else if (simp_eq_head_sym(opt_transform, "List")) {
            size_t m = opt_transform->data.function.arg_count;
            user_funcs = (Expr**)calloc(m ? m : 1, sizeof(Expr*));
            use_builtin = false; /* unless Automatic appears in the list */
            for (size_t i = 0; i < m; i++) {
                Expr* el = opt_transform->data.function.args[i];
                if (el->type == EXPR_SYMBOL && el->data.symbol.name == SYM_Automatic) {
                    use_builtin = true;
                } else {
                    user_funcs[n_user_funcs++] = el;
                }
            }
        } else {
            /* A bare function: {f}, built-ins suppressed. */
            user_funcs = (Expr**)calloc(1, sizeof(Expr*));
            user_funcs[0] = opt_transform;
            n_user_funcs = 1;
            use_builtin = false;
        }
    }

    /* ComplexityFunction -> Automatic is a synonym for the built-in
     * default. Treating it as NULL makes score_with_func use the fast
     * native simp_default_complexity path instead of evaluating
     * Automatic[candidate] (which would never reduce). */
    if (opt_complexity &&
        opt_complexity->type == EXPR_SYMBOL &&
        opt_complexity->data.symbol.name == SYM_Automatic) {
        opt_complexity = NULL;
    }

    /* Compute the effective assumption expression.
     *   - With Assumptions->X, X overrides the $Assumptions default.
     *   - Without, the positional assumption is appended to $Assumptions.
     * Then evaluate to canonicalise (e.g. And[True, x>0] -> x>0). */
    Expr* effective;
    if (opt_assumptions) {
        if (positional_assum) {
            effective = combine_with_and(expr_copy(positional_assum),
                                         expr_copy(opt_assumptions));
        } else {
            effective = evaluate(expr_copy(opt_assumptions));
        }
    } else {
        Expr* dollar = read_dollar_assumptions();
        if (positional_assum) {
            effective = combine_with_and(expr_copy(positional_assum), dollar);
        } else {
            effective = dollar;
        }
    }

    AssumeCtx* ctx = assume_ctx_from_expr(effective);
    expr_free(effective);

    /* If the input is a predicate that appears literally as one of our
     * assumed facts, it folds to True. This is a narrow win for simple
     * cases like Simplify[x > 0, x > 0]; it does not constitute a real
     * inequality reasoner (see Mathilda_spec.md for v1 gaps). */
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            if (expr_eq(expr, ctx->facts[i])) {
                assume_ctx_free(ctx);
                free(user_funcs);
                return expr_new_symbol(SYM_True);
            }
        }
    }

    /* Manual threading over Equal/Less/.../And/Or (List handled by
     * ATTR_LISTABLE on the Simplify symbol itself). For binary
     * relational heads we additionally try a rebalanced form
     * `pos OP neg` (after dividing through by the GCD of integer
     * coefficients) and pick the simpler of the two by SimplifyCount. */
    if (expr->type == EXPR_FUNCTION &&
        expr->data.function.head &&
        expr->data.function.head->type == EXPR_SYMBOL &&
        head_threads_over(expr->data.function.head->data.symbol.name)) {
        size_t n = expr->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** sub_args = (Expr**)calloc(argc, sizeof(Expr*));
            sub_args[0] = expr_copy(expr->data.function.args[i]);
            for (size_t k = 1; k < argc; k++) {
                sub_args[k] = expr_copy(res->data.function.args[k]);
            }
            Expr* call = expr_new_function(expr_new_symbol(SYM_Simplify), sub_args, argc);
            free(sub_args); /* expr_new_function copies the array contents */
            new_args[i] = evaluate(call);
            expr_free(call);
        }
        Expr* threaded = expr_new_function(expr_copy(expr->data.function.head), new_args, n);
        free(new_args);
        Expr* threaded_eval = eval_and_free(threaded);

        /* Rebalance candidate: only meaningful for a binary relation that
         * survived evaluation (Equal collapsed to True/False is not a
         * Function any more). */
        Expr* rebalanced = simp_try_rebalance_relation(threaded_eval);
        if (rebalanced && !expr_eq(rebalanced, threaded_eval)) {
            size_t s_threaded = score_with_func(threaded_eval, opt_complexity);
            size_t s_rebal    = score_with_func(rebalanced, opt_complexity);
            if (s_rebal < s_threaded) {
                expr_free(threaded_eval);
                threaded_eval = rebalanced;
            } else {
                expr_free(rebalanced);
            }
        } else if (rebalanced) {
            expr_free(rebalanced);
        }

        assume_ctx_free(ctx);
        /* The recursive per-component Simplify calls above re-parse the option
         * for each element, so the borrowed funcs list is not needed here. */
        free(user_funcs);
        return threaded_eval;
    }

    SimpMemo memo;
    simp_memo_init(&memo);

    FactorMemo* fmemo = factor_memo_new();
    factor_memo_push(fmemo);

    /* Top-level rational shortcut. simp_bottomup descends into every Plus /
     * Times child before dispatching at the top, and for a SHAPE_RATIONAL
     * input each child re-enters simp_dispatch -> simp_pipeline_rational.
     * Together / Cancel / Factor at the top combines all the children into
     * a single canonical num/den, so the per-child work is wasted: each
     * subnode's "best" form ends up subsumed by the top-level Together.
     *
     * Empirically, on multivariate rational inputs Simplify takes ~8 s vs
     * Cancel[Together[expr]] ~25 ms (~300x). Even when the search returns
     * the input unchanged, the cost is in the search itself. By dispatching
     * once at the top we cut directly to the pipeline that decides
     * acceptance against the input, bypassing the redundant per-subnode
     * traversal. The polish passes (lift_common_factor, PythagReduce,
     * canon_negate_pairs) still run on the result.
     *
     * Gated on SHAPE_RATIONAL: the classifier rejects inputs with trig,
     * log, abs, and non-integer powers, so we only take the shortcut when
     * the polynomial pipeline has full coverage. */
    /* The built-in transformation pipeline runs only when TransformationFunctions
     * permits it (Automatic, the default, or an explicit list that includes
     * Automatic). With an explicit user-only list we skip straight to applying
     * those functions to the input. */
    /* Arm the per-subexpression TimeConstraint for the simplification work
     * below. Save/restore the dynamically-scoped budget so nested Simplify
     * calls (e.g. a ComplexityFunction that itself simplifies, or the
     * per-element list threading) stack correctly and the outer budget is
     * always restored on return. All early returns above this point are
     * budget-neutral. */
    double saved_time_budget = simp_current_time_budget();
    simp_set_time_budget(time_budget);
    /* Arm the whole-call deadline too (save for restore; arm never relaxes an
     * outer deadline, so a nested Simplify can only tighten it). */
    double saved_call_deadline = simp_call_deadline();
    simp_arm_call_deadline(time_total);

    Expr* best = NULL;
    if (use_builtin) {
    if (simp_classify(expr) == SIMP_SHAPE_RATIONAL) {
        best = simp_dispatch(expr, ctx, opt_complexity);
    } else if ((best = transform_trigexp_vanish(expr)) != NULL) {
        /* Exp-kernel vanishing fast path. A rational function of trig/exp
         * kernels that is identically 0 — canonically a Risch antiderivative
         * diff-back D[G]-f — is proven 0 here, BEFORE the per-subnode
         * bottom-up descent that otherwise fires a full simp_search on every
         * internal node and hangs: >40 s on the multiple-angle Sec^n/Csc^n
         * forms, >90 s on symbolic-parameter I-laden forms (SIMPLIFY_GAPS.md
         * Families 1 & 3). TrigToExp + E^(k I x) -> t^k kernelization + exact
         * numerator zero-test — a genuine symbolic decision, no sampling.
         * Declines cheaply (NULL) on non-vanishing / non-single-kernel inputs,
         * so the general search below runs unchanged for everything else. */
    } else {
        /* Top-level algebraic-rational fast path. When the input is a
         * Plus over a multi-generator algebraic-number tower (e.g. the
         * output of D[Integrate[a x/(x^3+2), x], x] which is a sum of 3
         * fractions over {2^(1/3), Sqrt[3], Sqrt[radicand-with-α-inside]}),
         * Together[expr, Extension -> Automatic] is the one transform
         * that can collapse it back to (a x)/(x^3+2) in a single pass via
         * builtin_together's multi-gen single-α fallback (rat.c, Phase F).
         * simp_bottomup's per-subnode descent doesn't reach this combined-
         * over-common-denominator form on its own.  Strict leaf-count gate
         * ensures no regression on inputs where Together-with-Auto is a
         * no-op or worse. */
        Expr* alg_top = NULL;
        if (expr->type == EXPR_FUNCTION
            && expr->data.function.head
            && expr->data.function.head->type == EXPR_SYMBOL
            && expr->data.function.head->data.symbol.name == SYM_Plus
            && simp_has_rational_root(expr)
            && !contains_explicit_complex(expr)
            && !expr_has_nested_radical_radicand(expr)
            && !simp_contains_root_head(expr)) {
            /* Root[...] objects would be handed to Together[..., Extension ->
             * Automatic] / the poly engine as independent generators and blow
             * up (algebraically-dependent roots -> degenerate pseudo-remainder).
             * Skip this fast path; the qqbar coefficient pass does the algebra. */
            /* Two flavours of radical-fraction sum collapse here, both of
             * which simp_bottomup's per-subnode descent cannot reach on its
             * own (it simplifies each Plus child separately and never combines
             * them over a common denominator):
             *
             *   - Multi-generator algebraic-number towers (e.g. the sum of
             *     three fractions over {2^(1/3), Sqrt[3], ...} from
             *     D[Integrate[a x/(x^3+2), x], x]).  extension_autodetect
             *     returns n >= 2 and Together[expr, Extension -> Automatic]
             *     collapses it via builtin_together's Phase F multi-gen
             *     single-alpha fallback.
             *
             *   - A single polynomial-radicand radical (e.g. the cube-root
             *     tower (1 - b x)^(p/3) in D[Integrate[(1-b x)^(1/3)/x, x], x]).
             *     extension_autodetect returns NULL here -- the QAExt
             *     machinery only represents integer-base minimal polynomials --
             *     so plain Together does the work via its Phase E poly-radical
             *     reduction.
             *
             * The outer gate is simp_has_rational_root (a genuine Power[_,
             * Rational[_, q>=2]]) rather than has_non_integer_power so we never
             * feed a symbolic-exponent power (a^x) to Together, which can hang
             * on those.  The strict leaf-count gate below keeps the combined
             * fraction only when it is genuinely simpler, so there is no
             * regression when Together is a no-op or worse.
             *
             * Layer-0 prefilter shared with builtin_together/cancel: the
             * nested-radical-radicand condition above skips the expensive
             * primitive-element compositum construction inside
             * extension_autodetect. */
            QATower* qa_t = extension_autodetect(expr);
            /* Use Together[Extension -> Automatic] whenever ANY algebraic
             * tower is detected — including a single generator (n == 1).
             * For a single radical (Sqrt) plain Together happens to work
             * because Mathilda auto-reduces Sqrt powers, but for a single
             * cyclotomic generator (e.g. (-1)^(2/3)) plain Together treats
             * the root of unity as opaque and blows up super-polynomially;
             * the extension path reduces it modulo Φ_n instead.  The
             * leaf-count gate below still rejects any non-improving result,
             * so radical cases are unaffected. */
            bool use_extension = (qa_t != NULL);
            if (qa_t) qa_tower_free(qa_t);

            Expr* tog;
            if (use_extension) {
                tog = expr_new_function(
                    expr_new_symbol(SYM_Together),
                    (Expr*[]){
                        expr_copy(expr),
                        expr_new_function(expr_new_symbol(SYM_Rule),
                            (Expr*[]){expr_new_symbol(SYM_Extension),
                                      expr_new_symbol(SYM_Automatic)}, 2)
                    }, 2);
            } else {
                tog = expr_new_function(
                    expr_new_symbol(SYM_Together),
                    (Expr*[]){ expr_copy(expr) }, 1);
            }
            Expr* cand = evaluate(tog);
            expr_free(tog);
            if (cand && simp_default_complexity(cand)
                            < simp_default_complexity(expr)) {
                alg_top = cand;
            } else if (cand) {
                expr_free(cand);
            }
        }

        if (alg_top) {
            /* Use the algebraic collapse as the starting point.  Run
             * simp_bottomup on it to apply any further polish (rare for
             * already-canonical rational forms but harmless). */
            best = simp_bottomup(alg_top, ctx, opt_complexity, &memo, 0);
            size_t s_alg = score_with_func(alg_top, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_alg <= s_best) {
                expr_free(best);
                best = alg_top;
            } else {
                expr_free(alg_top);
            }
        } else {
            /* Top-level trig-rational fast path. Substitutes Sin/Cos/Sinh/
             * Cosh (and Tan/Cot/Sec/Csc/Tanh/etc. after preprocessing) plus
             * every opaque non-rational subtree (Log[...], Exp[...], etc.)
             * into fresh ground-field symbols so the algebraic core sees a
             * pure rational function; works in the quotient ring modulo the
             * trig/hyp ideals, then back-substitutes. Strict leaf-count gate
             * inside ensures it never regresses; on no-improvement or when
             * the input is out of budget it returns NULL and we fall through
             * to the normal bottom-up search. Doing this BEFORE simp_bottomup
             * means we bypass the per-subnode descent (which itself is
             * extremely slow on inputs like
             *   D[Integrate[Tan[x]^2 + Tan[x] + 1, x], x]
             * because every internal node fires a full simp_search). */
            Expr* tr = simp_trig_rational(expr, ctx, opt_complexity);
            if (tr) {
                best = tr;
            } else {
                best = simp_bottomup(expr, ctx, opt_complexity, &memo, 0);
            }
        }
    }

    /* Final-form polish: lift a shared algebraic generator out of a
     * top-level Plus (or out of a Plus child of a top-level Times -- the
     * numerator of a fraction with a non-integer-power denominator).
     * This catches:
     *   (8/105)(1+x^2)^(3/2) - (4/35)x^2(1+x^2)^(3/2) + (1/7)x^4(1+x^2)^(3/2)
     *     -> (1/105)(1+x^2)^(3/2)(8 - 12 x^2 + 15 x^4)
     *   (15 x^2 + 5 x^3)/(5+2x)^(3/2)
     *     -> (5 x^2 (3 + x))/(5+2x)^(3/2)
     * which Mathilda's polynomial Factor cannot reach because Variables[]
     * does not return non-integer-power generators. We apply it once at
     * the top level rather than as a seed in simp_search to avoid
     * destabilising the heuristic search on multi-variable trig inputs. */
    {
        Expr* lifted = simp_lift_common_factor(best);
        if (lifted && !expr_eq(lifted, best)) {
            size_t s_lift = score_with_func(lifted, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_lift <= s_best) {
                expr_free(best);
                best = lifted;
            } else {
                expr_free(lifted);
            }
        } else if (lifted) {
            expr_free(lifted);
        }
    }

    /* Pythagorean polish: PythagReduce already runs as a seed inside
     * simp_search, but its result enters update_best with a strict `<`
     * tiebreak, so structurally-collapsed forms that tie on
     * SimplifyCount (e.g. `-Sech[x]^2` vs `-1 + Tanh[x]^2`, both score
     * 7) lose to whatever arrived at the score plateau first. As a
     * polish, accept on `<=`: when the pythag rules turn the result
     * into a single Power-of-trig head with the same score or lower,
     * take it. Bypass when the Tanh/Coth/Tan/Cot rules cannot fire
     * (no relevant head present). */
    {
        Expr* reduced = transform_pythag_reduce(best);
        if (reduced && !expr_eq(reduced, best)) {
            size_t s_red  = score_with_func(reduced, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_red <= s_best) {
                expr_free(best);
                best = reduced;
            } else {
                expr_free(reduced);
            }
        } else if (reduced) {
            expr_free(reduced);
        }
    }

    /* Sign canonicalisation: flip pairs of negative-leading Plus factors
     * inside a top-level Times so each binomial leads with its
     * positive-coefficient term, e.g.
     *   ((-a + c) (-b + d))/(a b c d)  ->  ((a - c) (b - d))/(a b c d)
     * Value-preserving (flips occur in pairs so signs cancel). */
    {
        Expr* canon = canon_negate_pairs(best);
        if (canon) {
            expr_free(best);
            best = canon;
        }
    }

    /* Radical-fraction polish.  Simplify's Together/Cancel can leave a radical
     * Sqrt[p] in the NUMERATOR over an EXPANDED polynomial denominator, e.g.
     *   (Sqrt[6] Sqrt[6 + x^2]) / (6 x + x^3)
     * whose factored denominator x (6 + x^2) would cancel the radical back
     * into the denominator via the automatic Power[p, 1/2] Power[p, -1] ->
     * Power[p, -1/2] rule, giving WL's canonical Sqrt[6]/(x Sqrt[6 + x^2]).
     * The per-pass search never re-Factors its own Together output, so it
     * stops one step short.  Factoring the result exposes the shared radicand.
     * Gated on an actual radical (rational-root) present and a STRICT score
     * improvement, so non-radical results and no-improvement cases are
     * untouched, and Factor's own cost gates keep it bounded. */
    if (simp_has_rational_root(best) && has_compound_radicand(best)
            && !simp_contains_root_head(best)) {
        /* A Root[...] object is a constant algebraic number that Factor treats
         * as an independent polynomial generator; with several algebraically-
         * dependent roots the multivariate factoriser blows up exponentially
         * (the Sqrt[Tan[x]]-with-Root-coefficients hang). Decline here and let
         * the qqbar coefficient pass handle the algebra. */
        Expr* fac = expr_new_function(expr_new_symbol(SYM_Factor),
                                      (Expr*[]){ expr_copy(best) }, 1);
        Expr* factored = eval_and_free(fac);
        if (factored && !expr_eq(factored, best)) {
            size_t s_fac  = score_with_func(factored, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_fac < s_best) {
                expr_free(best);
                best = factored;
            } else {
                expr_free(factored);
            }
        } else if (factored) {
            expr_free(factored);
        }
    }

    /* Constant-algebraic fold.  A variable-free expression built from Sin/Cos/
     * Tan/... at rational multiples of Pi (together with radicals, roots of
     * unity and Root objects) is an exact algebraic number.  The shape router
     * sends a bare trig constant to simp_trig_rational, which never reaches
     * simp_search's qqbar pre-pass, so the fold would otherwise be missed:
     * Cos[Pi/7] Cos[2 Pi/7] Cos[3 Pi/7] - 1/8 -> 0, Tan[Pi/7] Tan[2 Pi/7]
     * Tan[3 Pi/7] - Sqrt[7] -> 0.  Taken when no worse (<=, matching the
     * pre-pass); flint_qqbar_reduce_coeffs is a no-op copy when nothing folds.
     * (contains_variable is unusable here: it counts the Cos/Sin head symbol as
     * a variable; is_constant_algebraic is the precise gate — 0 on any free
     * variable, 1 on trig-at-rational-Pi / radicals / Root.) */
    if (flint_qqbar_is_constant_algebraic(best)) {
        Expr* q = flint_qqbar_reduce_coeffs(best, QQBAR_METHOD_AUTOMATIC);
        if (q) q = eval_and_free(q);
        if (q && !expr_eq(q, best) &&
            score_with_func(q, opt_complexity)
                <= score_with_func(best, opt_complexity)) {
            expr_free(best);
            best = q;
        } else if (q) {
            expr_free(q);
        }
    }

    /* Inverse-trig rational-angle combination fold.  ArcTan/ArcSin/ArcCos of
     * constants are transcendental (not constant-algebraic, so the qqbar fold
     * above skips them), but a Z-combination that is identically a multiple of
     * Pi is decided exactly via e^{i e} == 1.  Closes the Machin-type identity
     * 4 ArcTan[1/5] - ArcTan[1/239] - Pi/4 -> 0 and ArcSin[3/5] + ArcSin[5/13]
     * - ArcSin[56/65] -> 0.  simp_invtrig_combo_is_zero returns 1 only on a
     * proven zero, so 0 (the simplest form) always wins. */
    if (best->type == EXPR_FUNCTION && simp_invfn_combo_is_zero(best)) {
        expr_free(best);
        best = expr_new_integer(0);
    }

    /* Region-identity fold: an inverse-function addition identity that holds on
     * the assumption region (a box) is decided by the derivative-constancy
     * engine -- grad == 0, branch-cut-free on the box (Reduce), value 0 at a
     * sample point.  Only under assumptions; declines otherwise. */
    /* Inverse-function identities on the assumption region: the cheap pointwise
     * hyperbolic-addition recognizer (ArcSinh/ArcCosh/ArcTanh) first, then the
     * general derivative-constancy region engine.  Both use region_shape_ok, so
     * their recursive Simplify calls cannot re-enter (no guard needed). */
    if (best->type == EXPR_FUNCTION && ctx && ctx->count > 0 &&
        (simp_invhyp_addition_is_zero(best, ctx) ||
         simp_region_identity_is_zero(best, ctx))) {
        expr_free(best);
        best = expr_new_integer(0);
    }
    } else {
        /* TransformationFunctions -> {f1, ...} without Automatic: the
         * built-in pipeline is suppressed, so the search starts from the
         * input itself. */
        best = expr_copy(expr);
    }

    /* Apply any user-supplied transformation functions, keeping the
     * lowest-complexity result. Runs after the built-in pipeline (when
     * enabled) so the user functions refine the best built-in form. */
    if (n_user_funcs > 0) {
        Expr* t = simp_apply_transformations(best, user_funcs, n_user_funcs,
                                             opt_complexity);
        expr_free(best);
        best = t;
    }
    free(user_funcs);

    factor_memo_pop();
    factor_memo_free(fmemo);

    simp_memo_free(&memo);
    assume_ctx_free(ctx);
    simp_set_time_budget(saved_time_budget);     /* restore dynamic scope */
    simp_set_call_deadline(saved_call_deadline); /* restore whole-call deadline */
    return best;
}

/* ----------------------------------------------------------------------- */
/* builtin_assuming -- desugar to Block[{$Assumptions = $A && a}, body]    */
/* ----------------------------------------------------------------------- */

Expr* builtin_assuming(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;

    Expr* assum = res->data.function.args[0];   /* already evaluated */
    Expr* body  = res->data.function.args[1];   /* held by HoldRest */

    /* Convert lists of assumptions to conjunctions, per Mathematica
     * semantics. */
    Expr* assum_norm;
    if (assum->type == EXPR_FUNCTION &&
        assum->data.function.head &&
        assum->data.function.head->type == EXPR_SYMBOL &&
        assum->data.function.head->data.symbol.name == SYM_List) {
        size_t n = assum->data.function.arg_count;
        Expr** copies = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) copies[i] = expr_copy(assum->data.function.args[i]);
        Expr* and_call = expr_new_function(expr_new_symbol(SYM_And), copies, n);
        free(copies);
        assum_norm = and_call;  /* not yet evaluated; Block will evaluate it */
    } else {
        assum_norm = expr_copy(assum);
    }

    /* Build $Assumptions && assum_norm */
    Expr* and_args[2] = { expr_new_symbol(SYM_DollarAssumptions), assum_norm };
    Expr* combined = expr_new_function(expr_new_symbol(SYM_And), and_args, 2);

    /* Build Set[$Assumptions, combined] -- represents
     * "$Assumptions = $Assumptions && a" inside the Block var list. */
    Expr* set_args[2] = { expr_new_symbol(SYM_DollarAssumptions), combined };
    Expr* set_call = expr_new_function(expr_new_symbol(SYM_Set), set_args, 2);

    /* Block[{Set[$Assumptions, ...]}, body] */
    Expr* var_list_args[1] = { set_call };
    Expr* var_list = expr_new_function(expr_new_symbol(SYM_List), var_list_args, 1);

    Expr* block_args[2] = { var_list, expr_copy(body) };
    Expr* block_call = expr_new_function(expr_new_symbol(SYM_Block), block_args, 2);

    Expr* result = evaluate(block_call);
    expr_free(block_call);
    return result;
}


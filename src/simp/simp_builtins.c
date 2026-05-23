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
/* builtin_simplify                                                        */
/* ----------------------------------------------------------------------- */

static Expr* read_dollar_assumptions(void) {
    /* Read the OwnValue directly. We must NOT evaluate $Assumptions, because
     * once an assumption like Element[x, Reals] becomes the bound value, our
     * own Element evaluator would recurse on it (Element reads $Assumptions
     * to decide -> evaluator fires the OwnValue rule -> Element[x, Reals]
     * gets re-evaluated -> ...). The first OwnValue rule on a symbol is its
     * current value (newest first); we deep-copy its replacement. */
    Rule* r = symtab_get_own_values("$Assumptions");
    if (!r || !r->replacement) return expr_new_symbol("True");
    return expr_copy(r->replacement);
}

/* ----------------------------------------------------------------------- */
/* builtin_element -- Element[x, Domain]                                   */
/* ----------------------------------------------------------------------- */

static bool is_complex_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2;
}

bool is_rational_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
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

/* Element[x, dom] decision: 1 = True, 0 = False, -1 = undetermined. */
static int element_decide(const Expr* x, const char* dom, const AssumeCtx* ctx) {
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
            if (x->data.symbol == SYM_True)  return 1;
            if (x->data.symbol == SYM_False) return 1;
        }
        return -1;
    }

    if (strcmp(dom, "Primes") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) {
            Expr* primeq = call_unary_copy("PrimeQ", x);
            int ans = -1;
            if (primeq && primeq->type == EXPR_SYMBOL) {
                if (primeq->data.symbol == SYM_True) ans = 1;
                if (primeq->data.symbol == SYM_False) ans = 0;
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
                if (primeq->data.symbol == SYM_True) ans = 0;
                if (primeq->data.symbol == SYM_False) ans = 1;
            }
            if (primeq) expr_free(primeq);
            return ans;
        }
        return -1;
    }

    return -1;
}

Expr* builtin_element(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;

    Expr* x   = res->data.function.args[0];
    Expr* dom = res->data.function.args[1];

    /* Element[{x1, x2, ...}, dom] is shorthand for the conjunction
     * Element[x1, dom] && ... && Element[xN, dom]. We collapse it to a
     * single True/False only when every component decides; otherwise we
     * leave the original Element[{...}, dom] in place so downstream
     * consumers (Simplify's AssumeCtx; see ctx_walk in simp_assume.c)
     * can still treat it as a multi-variable real/integer/... assumption.
     *
     * Threading to a `List` of Element calls (the old behaviour) was
     * incorrect: Simplify carries ATTR_LISTABLE on the assumption
     * argument, so a list of partial Element facts would cause Simplify
     * to thread itself, splitting the joint assumption into per-variable
     * runs that each only see one of the facts. */
    if (x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol == SYM_List) {
        size_t n = x->data.function.arg_count;
        bool all_true = (n > 0);
        bool any_false = false;
        for (size_t i = 0; i < n; i++) {
            Expr* sub_args[2] = { expr_copy(x->data.function.args[i]), expr_copy(dom) };
            Expr* call = expr_new_function(expr_new_symbol("Element"), sub_args, 2);
            Expr* sub = evaluate(call);
            expr_free(call);
            if (sub && sub->type == EXPR_SYMBOL) {
                if (sub->data.symbol == SYM_True) {
                    /* keep all_true */
                } else if (sub->data.symbol == SYM_False) {
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
        if (any_false) return expr_new_symbol("False");
        if (all_true)  return expr_new_symbol("True");
        return NULL;
    }

    if (dom->type != EXPR_SYMBOL) return NULL;
    const char* d = dom->data.symbol;

    /* Build context from current $Assumptions. */
    Expr* dollar = read_dollar_assumptions();
    AssumeCtx* ctx = assume_ctx_from_expr(dollar);
    expr_free(dollar);

    int decision = element_decide(x, d, ctx);
    assume_ctx_free(ctx);

    if (decision == 1)  return expr_new_symbol("True");
    if (decision == 0)  return expr_new_symbol("False");
    return NULL;
}

static Expr* combine_with_and(Expr* a, Expr* b) {
    /* Both inputs owned and consumed. */
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol("And"), args, 2);
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
                    expr_new_symbol("Times"), args, n - 1);
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
            expr_new_symbol("Times"), args, n + 1);
        free(args);
        expr_free(rest);
        return out;
    }
    Expr* args[2] = { expr_new_integer(c), rest };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Returns NULL when no rebalanced form is produced (non-int64 coeffs,
 * fully symbolic d, or d = 0). The caller compares scores. */
static Expr* simp_try_rebalance_relation(const Expr* relation) {
    if (!relation || relation->type != EXPR_FUNCTION) return NULL;
    if (relation->data.function.arg_count != 2) return NULL;
    const Expr* h = relation->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return NULL;
    const char* hn = h->data.symbol;
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
        expr_new_symbol("Times"), neg_args, 2);
    Expr* d_args[2] = {
        expr_copy(relation->data.function.args[0]),
        neg_rhs
    };
    Expr* d_call = expr_new_function(
        expr_new_symbol("Plus"), d_args, 2);
    Expr* d_sum = evaluate(d_call);
    /* Expand so Times[2, Plus[...]] partitions term-by-term. The threaded
     * input may already have collected common factors via Collect, which
     * defeats coefficient-level rebalancing. */
    Expr* exp_args[1] = { d_sum };
    Expr* d_exp_call = expr_new_function(
        expr_new_symbol("Expand"), exp_args, 1);
    Expr* d = evaluate(d_exp_call);

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
                          expr_new_symbol("Plus"), pos, pn);

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
                        expr_new_symbol("Times"), args, 2);
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
                expr_new_symbol("Plus"), rhs_terms, rt);
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
    Expr* lhs_e = evaluate(new_lhs);
    Expr* rhs_e = evaluate(new_rhs);

    Expr* rel_args[2] = { lhs_e, rhs_e };
    Expr* out = expr_new_function(
        expr_new_symbol(out_head), rel_args, 2);
    Expr* out_eval = evaluate(out);
    return out_eval;
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
        return internal_rationalize_then_numericalize(res, builtin_simplify);
    }

    Expr* expr = res->data.function.args[0];

    /* Parse remaining args: at most one positional assumption, plus
     * options Rule[Assumptions, X] and Rule[ComplexityFunction, f]. */
    Expr* positional_assum = NULL;
    Expr* opt_assumptions  = NULL;
    Expr* opt_complexity   = NULL;

    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            opt_assumptions = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "ComplexityFunction")) {
            opt_complexity = a->data.function.args[1];
        } else if (positional_assum == NULL) {
            positional_assum = a;
        }
    }

    /* ComplexityFunction -> Automatic is a synonym for the built-in
     * default. Treating it as NULL makes score_with_func use the fast
     * native simp_default_complexity path instead of evaluating
     * Automatic[candidate] (which would never reduce). */
    if (opt_complexity &&
        opt_complexity->type == EXPR_SYMBOL &&
        opt_complexity->data.symbol == SYM_Automatic) {
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
                return expr_new_symbol("True");
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
        head_threads_over(expr->data.function.head->data.symbol)) {
        size_t n = expr->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** sub_args = (Expr**)calloc(argc, sizeof(Expr*));
            sub_args[0] = expr_copy(expr->data.function.args[i]);
            for (size_t k = 1; k < argc; k++) {
                sub_args[k] = expr_copy(res->data.function.args[k]);
            }
            Expr* call = expr_new_function(expr_new_symbol("Simplify"), sub_args, argc);
            new_args[i] = evaluate(call);
            expr_free(call);
        }
        Expr* threaded = expr_new_function(expr_copy(expr->data.function.head), new_args, n);
        free(new_args);
        Expr* threaded_eval = evaluate(threaded);

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
    Expr* best;
    if (simp_classify(expr) == SIMP_SHAPE_RATIONAL) {
        best = simp_dispatch(expr, ctx, opt_complexity);
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
            && expr->data.function.head->data.symbol == SYM_Plus
            && has_non_integer_power(expr)
            && !contains_explicit_complex(expr)
            && !expr_has_nested_radical_radicand(expr)) {
            /* Layer-0 prefilter shared with builtin_together/cancel:
             * skip extension_autodetect entirely when the input has a
             * nested radical and no free polynomial variable.  Even
             * with builtin_together's own gate, this saves a redundant
             * autodetect call here (extension_autodetect on a nested-
             * radical input runs the primitive-element compositum
             * construction, which is the bulk of the cost). */
            QATower* qa_t = extension_autodetect(expr);
            if (qa_t && qa_t->n >= 2) {
                qa_tower_free(qa_t);
                Expr* tog = expr_new_function(
                    expr_new_symbol("Together"),
                    (Expr*[]){
                        expr_copy(expr),
                        expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){expr_new_symbol("Extension"),
                                      expr_new_symbol("Automatic")}, 2)
                    }, 2);
                Expr* cand = evaluate(tog);
                expr_free(tog);
                if (cand && simp_default_complexity(cand)
                                < simp_default_complexity(expr)) {
                    alg_top = cand;
                } else if (cand) {
                    expr_free(cand);
                }
            } else if (qa_t) {
                qa_tower_free(qa_t);
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
             *   D[Integrate`RischNorman[Tan[x]^2 + Tan[x] + 1, x], x]
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

    factor_memo_pop();
    factor_memo_free(fmemo);

    simp_memo_free(&memo);
    assume_ctx_free(ctx);
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
        assum->data.function.head->data.symbol == SYM_List) {
        size_t n = assum->data.function.arg_count;
        Expr** copies = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) copies[i] = expr_copy(assum->data.function.args[i]);
        Expr* and_call = expr_new_function(expr_new_symbol("And"), copies, n);
        free(copies);
        assum_norm = and_call;  /* not yet evaluated; Block will evaluate it */
    } else {
        assum_norm = expr_copy(assum);
    }

    /* Build $Assumptions && assum_norm */
    Expr* and_args[2] = { expr_new_symbol("$Assumptions"), assum_norm };
    Expr* combined = expr_new_function(expr_new_symbol("And"), and_args, 2);

    /* Build Set[$Assumptions, combined] -- represents
     * "$Assumptions = $Assumptions && a" inside the Block var list. */
    Expr* set_args[2] = { expr_new_symbol("$Assumptions"), combined };
    Expr* set_call = expr_new_function(expr_new_symbol("Set"), set_args, 2);

    /* Block[{Set[$Assumptions, ...]}, body] */
    Expr* var_list_args[1] = { set_call };
    Expr* var_list = expr_new_function(expr_new_symbol("List"), var_list_args, 1);

    Expr* block_args[2] = { var_list, expr_copy(body) };
    Expr* block_call = expr_new_function(expr_new_symbol("Block"), block_args, 2);

    Expr* result = evaluate(block_call);
    expr_free(block_call);
    return result;
}


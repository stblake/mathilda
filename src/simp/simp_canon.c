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
/* Common-factor lift across Plus terms                                    */
/* ----------------------------------------------------------------------- */
/*
 * lift_common: given a Plus expression whose terms share a multiplicative
 * factor (an algebraic generator like (1+x^2)^(3/2), a free symbol, an
 * integer, or a Power[base, n>=1] that splits into n copies of base),
 * factor the common piece outside the Plus.
 *
 * Why a dedicated transform? Mathilda's Factor / FactorTerms decompose
 * polynomials over K[x_1, ..., x_n] using Variables[] to discover the
 * generator set, and Variables[] does not return non-integer Power
 * expressions (e.g. Sqrt[x], (1+x^2)^(3/2)). So a Plus that obviously
 * shares (1+x^2)^(3/2) across all terms slips past Factor untouched.
 * This transform takes a structural multiset view: a non-numeric factor
 * is either an algebraic generator (Power with non-integer exponent,
 * or any Power exponent we can't reason about) treated as one opaque
 * token, or a Power[base, n] with n a small positive integer that we
 * split into n copies of base. The numeric coefficients merge via
 * rational GCD. Lifted result: gcd_coef * Times(common_tokens) *
 * Plus[t_i / lift_factor], with the division handed back to evaluate()
 * for cancellation.
 *
 * Cases this enables that Factor alone cannot:
 *   - Plus of c_i * (1+x^2)^(3/2) * x^k  -> (1+x^2)^(3/2) * Plus[c_i x^k]
 *   - Plus inside Times[Plus, Power[denom, neg]] (rational expressions
 *     with non-integer-power denominator): factor the numerator only.
 *
 * Returns NULL when no nontrivial lift is possible (single-term Plus,
 * coprime coefficients with no shared token). */

typedef struct {
    Expr** items;      /* aliased pointers into the term tree (no ownership) */
    size_t count;
    size_t cap;
} LiftTokList;

static void lift_tl_init(LiftTokList* t) {
    t->items = NULL; t->count = 0; t->cap = 0;
}
static void lift_tl_free(LiftTokList* t) { free(t->items); }
static void lift_tl_push(LiftTokList* t, Expr* e) {
    if (t->count == t->cap) {
        size_t nc = t->cap ? t->cap * 2 : 4;
        Expr** ni = (Expr**)realloc(t->items, sizeof(Expr*) * nc);
        if (!ni) { /* OOM: drop the push silently rather than abort. */
            return;
        }
        t->items = ni; t->cap = nc;
    }
    t->items[t->count++] = e;
}

/* Convert an Expr* to mpq_t. Recognises EXPR_INTEGER, EXPR_BIGINT, and
 * Rational[n, d]. Returns false for anything else. */
static bool lift_expr_to_mpq(const Expr* e, mpq_t out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpq_set_z(out, e->data.bigint);
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpq_set_si(out, (long)n, (unsigned long)d);
        mpq_canonicalize(out);
        return true;
    }
    return false;
}

/* Build a normalised numeric Expr from an mpq_t. Returns Integer when
 * the denominator is 1, otherwise a Rational[n, d]. */
static Expr* lift_mpq_to_expr(const mpq_t v) {
    if (mpz_cmp_ui(mpq_denref(v), 1) == 0) {
        if (mpz_fits_slong_p(mpq_numref(v))) {
            return expr_new_integer((int64_t)mpz_get_si(mpq_numref(v)));
        }
        return expr_new_bigint_from_mpz(mpq_numref(v));
    }
    Expr* num = mpz_fits_slong_p(mpq_numref(v))
                  ? expr_new_integer((int64_t)mpz_get_si(mpq_numref(v)))
                  : expr_new_bigint_from_mpz(mpq_numref(v));
    Expr* den = mpz_fits_slong_p(mpq_denref(v))
                  ? expr_new_integer((int64_t)mpz_get_si(mpq_denref(v)))
                  : expr_new_bigint_from_mpz(mpq_denref(v));
    Expr* args[2] = { num, den };
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

/* Decompose one Plus term into (mpq coefficient, list of token aliases).
 * Numeric leaves accumulate into *coef. Power[base, n] with 1<=n<=16 is
 * split into n copies of base. Power[base, exp] with any other exp shape
 * (rational, negative, symbolic) is treated as one opaque token.
 *
 * Recurses into nested Times: in practice Mathilda's Plus does not always
 * fully flatten Times children -- a literal-times-product subexpression
 * inside a Plus surfaces as Times[c, Times[a, b]] -- so we walk the
 * subtree rather than relying on a one-level-deep view. */
static void lift_decompose_term(Expr* term, mpq_t coef, LiftTokList* tokens) {
    if (!term) return;
    if (term->type == EXPR_FUNCTION
        && term->data.function.head
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < term->data.function.arg_count; i++) {
            lift_decompose_term(term->data.function.args[i], coef, tokens);
        }
        return;
    }
    mpq_t tmp; mpq_init(tmp);
    if (lift_expr_to_mpq(term, tmp)) {
        mpq_mul(coef, coef, tmp);
        mpq_clear(tmp);
        return;
    }
    mpq_clear(tmp);
    if (term->type == EXPR_FUNCTION
        && term->data.function.head
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Power
        && term->data.function.arg_count == 2) {
        Expr* exp = term->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer >= 1
            && exp->data.integer <= 16) {
            Expr* base = term->data.function.args[0];
            for (int64_t k = 0; k < exp->data.integer; k++) {
                lift_tl_push(tokens, base);
            }
            return;
        }
    }
    lift_tl_push(tokens, term);
}

/* Find common multiset of tokens. Greedy: walk the first list; for each
 * token, search for an unused match in every other list. If all match,
 * mark them used and add the token to the result. */
static void lift_find_common(LiftTokList* lists, size_t n_terms,
                             Expr*** out_common, size_t* out_count) {
    *out_common = NULL; *out_count = 0;
    if (n_terms == 0 || lists[0].count == 0) return;
    Expr** result = (Expr**)malloc(sizeof(Expr*) * lists[0].count);
    if (!result) return;
    size_t res_count = 0;

    char** used = (char**)malloc(sizeof(char*) * n_terms);
    if (!used) { free(result); return; }
    for (size_t i = 0; i < n_terms; i++) {
        used[i] = lists[i].count ? (char*)calloc(lists[i].count, 1) : NULL;
    }

    size_t* idx = (size_t*)malloc(sizeof(size_t) * n_terms);
    if (!idx) { free(result); for (size_t i = 0; i < n_terms; i++) free(used[i]); free(used); return; }

    for (size_t j = 0; j < lists[0].count; j++) {
        if (used[0][j]) continue;
        Expr* tok = lists[0].items[j];
        idx[0] = j;
        bool ok = true;
        for (size_t i = 1; i < n_terms; i++) {
            bool found = false;
            for (size_t k = 0; k < lists[i].count; k++) {
                if (!used[i][k] && expr_eq(tok, lists[i].items[k])) {
                    idx[i] = k; found = true; break;
                }
            }
            if (!found) { ok = false; break; }
        }
        if (ok) {
            for (size_t i = 0; i < n_terms; i++) used[i][idx[i]] = 1;
            result[res_count++] = tok;
        }
    }
    free(idx);
    for (size_t i = 0; i < n_terms; i++) free(used[i]);
    free(used);

    *out_common = result;
    *out_count = res_count;
}

/* GCD of an array of mpq values: gcd(numerators) / lcm(denominators).
 * Result is positive. n must be >= 1. */
static void lift_compute_mpq_gcd(mpq_t* values, size_t n, mpq_t out) {
    mpq_set(out, values[0]);
    if (mpq_sgn(out) < 0) mpq_neg(out, out);
    for (size_t i = 1; i < n; i++) {
        mpz_t b_num, g, lcm_d;
        mpz_inits(b_num, g, lcm_d, NULL);
        mpz_set(b_num, mpq_numref(values[i]));
        mpz_abs(b_num, b_num);
        mpz_gcd(g, mpq_numref(out), b_num);
        mpz_lcm(lcm_d, mpq_denref(out), mpq_denref(values[i]));
        mpz_set(mpq_numref(out), g);
        mpz_set(mpq_denref(out), lcm_d);
        mpq_canonicalize(out);
        mpz_clears(b_num, g, lcm_d, NULL);
    }
}

static Expr* lift_common_from_plus_impl(const Expr* plus_e) {
    if (!plus_e || plus_e->type != EXPR_FUNCTION
        || !plus_e->data.function.head
        || plus_e->data.function.head->type != EXPR_SYMBOL
        || plus_e->data.function.head->data.symbol != SYM_Plus) {
        return NULL;
    }
    size_t n = plus_e->data.function.arg_count;
    if (n < 2) return NULL;

    LiftTokList* lists = (LiftTokList*)malloc(sizeof(LiftTokList) * n);
    mpq_t* coefs = (mpq_t*)malloc(sizeof(mpq_t) * n);
    if (!lists || !coefs) { free(lists); free(coefs); return NULL; }
    for (size_t i = 0; i < n; i++) {
        lift_tl_init(&lists[i]);
        mpq_init(coefs[i]);
        mpq_set_ui(coefs[i], 1, 1);
        lift_decompose_term(plus_e->data.function.args[i], coefs[i], &lists[i]);
    }

    Expr** common = NULL;
    size_t common_count = 0;
    lift_find_common(lists, n, &common, &common_count);

    /* Restrict the firing condition to a real shared algebraic factor.
     * A coefficient-only GCD lift (e.g. Plus[a/3, b/9] -> (1/9)(3a + b))
     * doesn't reveal new structure that Together/Cancel haven't already
     * exposed, and feeding it through the round loop can blow up
     * downstream transform cost (the evaluator re-rationalises Sqrt
     * products, etc.). */
    if (common_count == 0) {
        for (size_t i = 0; i < n; i++) { lift_tl_free(&lists[i]); mpq_clear(coefs[i]); }
        free(lists); free(coefs); free(common);
        return NULL;
    }

    mpq_t cgcd; mpq_init(cgcd);
    lift_compute_mpq_gcd(coefs, n, cgcd);

    /* Build lifted multiplicative factor: cgcd * Times(common). */
    Expr* gcd_expr = lift_mpq_to_expr(cgcd);
    Expr* lift_factor;
    if (common_count == 0) {
        lift_factor = gcd_expr;
    } else {
        Expr** args = (Expr**)malloc(sizeof(Expr*) * (common_count + 1));
        args[0] = gcd_expr;
        for (size_t i = 0; i < common_count; i++) {
            args[i + 1] = expr_copy(common[i]);
        }
        lift_factor = expr_new_function(expr_new_symbol(SYM_Times), args, common_count + 1);
        free(args);
    }
    Expr* lift_factor_eval = eval_and_free(lift_factor);

    /* For each term, divide by lift_factor and let evaluate() cancel. */
    Expr** new_terms = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* inv_args[2] = { expr_copy(lift_factor_eval), expr_new_integer(-1) };
        Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
        Expr* mul_args[2] = { expr_copy(plus_e->data.function.args[i]), inv };
        Expr* div = expr_new_function(expr_new_symbol(SYM_Times), mul_args, 2);
        new_terms[i] = eval_and_free(div);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), new_terms, n);
    free(new_terms);
    Expr* sum_eval = eval_and_free(sum);

    Expr* res_args[2] = { lift_factor_eval, sum_eval };
    Expr* result = expr_new_function(expr_new_symbol(SYM_Times), res_args, 2);
    Expr* result_eval = eval_and_free(result);

    for (size_t i = 0; i < n; i++) { lift_tl_free(&lists[i]); mpq_clear(coefs[i]); }
    free(lists); free(coefs); free(common);
    mpq_clear(cgcd);

    /* If the lift was a no-op (the Plus stayed structurally identical
     * after evaluate normalised the round trip), report no improvement. */
    if (expr_eq(result_eval, plus_e)) {
        expr_free(result_eval);
        return NULL;
    }
    return result_eval;
}

/* Walker entry point. Tries the lift on the input directly (when it's a
 * Plus) or on a Plus child of a Times product (e.g. the numerator of a
 * Times[Plus, Power[denom, -negative]] fraction). Returns NULL when no
 * structural improvement is found. */
Expr* simp_lift_common_factor(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) {
        return NULL;
    }
    const char* sym = e->data.function.head->data.symbol;
    if (sym == SYM_Plus) {
        return lift_common_from_plus_impl(e);
    }
    if (sym == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* child = e->data.function.args[i];
            if (child && child->type == EXPR_FUNCTION
                && child->data.function.head
                && child->data.function.head->type == EXPR_SYMBOL
                && child->data.function.head->data.symbol == SYM_Plus) {
                Expr* lifted = lift_common_from_plus_impl(child);
                if (lifted) {
                    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * e->data.function.arg_count);
                    for (size_t j = 0; j < e->data.function.arg_count; j++) {
                        new_args[j] = (j == i)
                                          ? lifted
                                          : expr_copy(e->data.function.args[j]);
                    }
                    Expr* new_times = expr_new_function(expr_new_symbol(SYM_Times),
                                                        new_args, e->data.function.arg_count);
                    free(new_args);
                    Expr* res = eval_and_free(new_times);
                    if (expr_eq(res, e)) { expr_free(res); return NULL; }
                    return res;
                }
            }
        }
    }
    return NULL;
}

/* ----------------------------------------------------------------------- */
/* Sign canonicalization: paired negative-leading Plus factors            */
/* ----------------------------------------------------------------------- */
/*
 * Mathematica's Factor (and our facpoly) emit Plus subterms whose
 * canonically-sorted first argument is the smallest by name -- so for
 * `(c - a)(d - b)` the printed form is `((-a + c)(-b + d))` because the
 * Plus[Times[-1, a], c] sorts the negated `a` before the bare `c` (the
 * sort key strips the leading -1 coefficient).
 *
 * Mathematica's printed convention prefers each binomial to lead with a
 * positive coefficient. We achieve this post-hoc: when a Times has an
 * even number of "negatively-leading" Plus factors, flipping the sign
 * of each is value-preserving (each pair contributes (-1)*(-1) = 1) and
 * leaves the printed form leading with the positive coefficient.
 *
 * We do not attempt the odd-count case here: pulling an extra -1 onto
 * an outer numeric factor changes which token absorbs the sign and is
 * not always a canonical win for the score function. */

static bool plus_arg_is_negative_leading(const Expr* arg) {
    if (!arg) return false;
    if (arg->type == EXPR_INTEGER) return arg->data.integer < 0;
    if (arg->type == EXPR_BIGINT) return mpz_sgn(arg->data.bigint) < 0;
    if (is_rational_literal(arg)
        && arg->data.function.args[0]->type == EXPR_INTEGER) {
        return arg->data.function.args[0]->data.integer < 0;
    }
    if (arg->type == EXPR_FUNCTION
        && arg->data.function.head
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol == SYM_Times
        && arg->data.function.arg_count >= 1) {
        Expr* coef = arg->data.function.args[0];
        if (coef->type == EXPR_INTEGER) return coef->data.integer < 0;
        if (coef->type == EXPR_BIGINT) return mpz_sgn(coef->data.bigint) < 0;
        if (is_rational_literal(coef)
            && coef->data.function.args[0]->type == EXPR_INTEGER) {
            return coef->data.function.args[0]->data.integer < 0;
        }
    }
    return false;
}

static bool plus_is_negative_leading(const Expr* p) {
    if (!p || p->type != EXPR_FUNCTION
        || !p->data.function.head
        || p->data.function.head->type != EXPR_SYMBOL
        || p->data.function.head->data.symbol != SYM_Plus
        || p->data.function.arg_count < 1) {
        return false;
    }
    return plus_arg_is_negative_leading(p->data.function.args[0]);
}

/* Build a Plus equal to -p by negating every term and re-evaluating so
 * Mathilda re-canonicalises the argument order. */
static Expr* plus_negate(const Expr* p) {
    size_t n = p->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* neg_args[2] = { expr_new_integer(-1),
                              expr_copy(p->data.function.args[i]) };
        args[i] = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    }
    Expr* neg_plus = expr_new_function(expr_new_symbol(SYM_Plus), args, n);
    free(args);
    return eval_and_free(neg_plus);
}

/* Walk a Times. If two or more Plus children are negative-leading, flip
 * pairs of them. Only flips an even number; the odd remainder stays.
 * Returns NULL if no flip applies. */
Expr* canon_negate_pairs(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION
        || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol != SYM_Times) {
        return NULL;
    }
    size_t n = e->data.function.arg_count;
    if (n < 2) return NULL;

    size_t* neg_idx = (size_t*)malloc(sizeof(size_t) * n);
    if (!neg_idx) return NULL;
    size_t neg_count = 0;
    for (size_t i = 0; i < n; i++) {
        if (plus_is_negative_leading(e->data.function.args[i])) {
            neg_idx[neg_count++] = i;
        }
    }
    size_t flip_count = (neg_count / 2) * 2;
    if (flip_count == 0) {
        free(neg_idx);
        return NULL;
    }

    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = expr_copy(e->data.function.args[i]);
    }
    for (size_t j = 0; j < flip_count; j++) {
        size_t i = neg_idx[j];
        Expr* flipped = plus_negate(new_args[i]);
        expr_free(new_args[i]);
        new_args[i] = flipped;
    }
    free(neg_idx);
    Expr* new_times = expr_new_function(expr_copy(e->data.function.head),
                                        new_args, n);
    free(new_args);
    Expr* result = eval_and_free(new_times);
    if (expr_eq(result, e)) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Pythagorean reduction:
 *   1 - Cos[x]^2  -> Sin[x]^2
 *   1 - Sin[x]^2  -> Cos[x]^2
 *   Cosh[x]^2 - 1 -> Sinh[x]^2
 *   1 + Sinh[x]^2 -> Cosh[x]^2
 *
 * Each rule is a strict leaf-count reduction (the Plus collapses to a
 * single Power). The trailing `r___` inside the Plus lets the rule fire
 * when the matching pair sits among other terms (e.g.
 * `1 - Cos[x]^2 + 5` -> `5 + Sin[x]^2`). Idempotent on inputs that
 * don't match. */
/* Wrap a transform's `impl` with FactorMemo lookup + store.  When no
 * memo is active (i.e., we're not inside a Simplify call), the impl
 * runs directly with no overhead.  When active, identical inputs
 * return cached results; the memo key includes a $-prefixed pseudo-
 * head so it never collides with builtin keys (Factor[X], TrigFactor[X],
 * etc.) sharing the same memo. */
Expr* simp_memo_wrap(const Expr* e, const char* pseudo_head,
                            Expr* (*impl)(const Expr*)) {
    FactorMemo* memo = factor_memo_active();
    if (!memo) return impl(e);

    /* Note: we use raw-input keying here (not Together(Expand(.))
     * canonicalisation as in trig_memo_call).  Reason: the wrapped
     * transforms (PythagReduce, PythagSquareComplete, HalfAngle) use
     * pattern rules that look for specific surface structure --
     * `1 - Cos[x]^2`, `1 + 2 Sin Cos`, `Sin[x] / (1 + Cos[x])` etc.
     * Distributive Expand destroys those patterns (`a (-1 + Cos^2)`
     * becomes `-a + a Cos^2`, where the -1 disappears as a coefficient
     * adjustment), so the rules no longer fire on the canonical form.
     *
     * For the trig memos the canonical form is fine because those
     * transforms internally normalise via Together / TrigToExp before
     * pattern matching. */
    Expr* key_args[1] = { expr_copy((Expr*)e) };
    Expr* key = expr_new_function(expr_new_symbol(pseudo_head), key_args, 1);
    const Expr* hit = factor_memo_lookup(memo, key);
    if (hit) {
        Expr* cached = expr_copy((Expr*)hit);
        expr_free(key);
        return cached;
    }
    Expr* result = impl(e);
    if (result) factor_memo_store(memo, key, result);
    expr_free(key);
    return result;
}

/* Cheap structural check: does `e` contain any of Cos/Sin/Cosh/Sinh,
 * Tan/Cot/Tanh/Coth, or Sec/Csc/Sech/Csch as a function head?
 * PythagReduce / PythagCanon's rules can only fire on these heads, so
 * when the answer is no we can skip the ReplaceRepeated walk entirely.
 * Walks the tree once; cheaper by orders of magnitude than the
 * pattern-matching pass it gates.  Sec/Csc/Sech/Csch were added
 * alongside the Sec[x]^2 -> 1 + Tan[x]^2 substitution direction in
 * PythagCanon (and its three siblings); without listing them here,
 * inputs like (Sec[x]+1)(Sec[x]-1) - Tan[x]^2 would skip PythagCanon
 * via the gate and land in TrigFactor (700 ms+ on multi-variable
 * expansions). */
bool has_pythag_head(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;
        if (h == SYM_Cos || h == SYM_Sin ||
            h == SYM_Cosh || h == SYM_Sinh ||
            h == SYM_Tan || h == SYM_Cot ||
            h == SYM_Tanh || h == SYM_Coth ||
            h == SYM_Sec || h == SYM_Csc ||
            h == SYM_Sech || h == SYM_Csch) {
            return true;
        }
    }
    if (has_pythag_head(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_pythag_head(e->data.function.args[i])) return true;
    }
    return false;
}


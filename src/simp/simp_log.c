/*
 * simp_log.c -- Log-specific simplifications for Simplify.
 *
 * See simp_log.h for the high-level contract and the two-primitive
 * design (Pass A: Log[positive_rational] prime decomposition;
 * Pass B: Plus[ci Log[ai]] fuser). Wired into both
 * simp_pipeline_logexp and the simp_search seed phase so it fires
 * regardless of whether the input is classified LOGEXP or TRIG.
 */

#include "simp_log.h"
#include "expr.h"
#include "eval.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

/* assume_known_* / simp_default_complexity are exposed via simp.h. */

/* ----------------------------------------------------------------------- */
/* Small helpers                                                           */
/* ----------------------------------------------------------------------- */

static Expr* make_call1(const char* head_name, Expr* a0) {
    Expr* a[1] = { a0 };
    return expr_new_function(expr_new_symbol(head_name), a, 1);
}

static Expr* make_call2(const char* head_name, Expr* a0, Expr* a1) {
    Expr* a[2] = { a0, a1 };
    return expr_new_function(expr_new_symbol(head_name), a, 2);
}

/* Evaluate `owned` and return the result, freeing the input. */
static Expr* eval_take(Expr* owned) {
    Expr* r = evaluate(owned);
    expr_free(owned);
    return r;
}

/* Is e structurally Log[arg] (single arg)? */
static bool is_log1(const Expr* e) {
    return e &&
           e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_Log &&
           e->data.function.arg_count == 1;
}

/* Extract a positive rational p/q from `e` into freshly-init'd num/den.
 * Returns false (and leaves num/den uninitialised) when e isn't a positive
 * rational, an integer, or a bigint > 0. The caller owns num/den on success
 * and must mpz_clear them. */
static bool extract_pos_rational(const Expr* e, mpz_t num, mpz_t den) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer <= 0) return false;
        mpz_init_set_si(num, e->data.integer);
        mpz_init_set_ui(den, 1);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        if (mpz_sgn(e->data.bigint) <= 0) return false;
        mpz_init_set(num, e->data.bigint);
        mpz_init_set_ui(den, 1);
        return true;
    }
    /* Rational literal: Rational[p, q] with integer-typed children. */
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2) {
        Expr* p = e->data.function.args[0];
        Expr* q = e->data.function.args[1];
        if (!expr_is_integer_like(p) || !expr_is_integer_like(q)) return false;
        mpz_init(num); mpz_init(den);
        expr_to_mpz(p, num);
        expr_to_mpz(q, den);
        if (mpz_sgn(num) <= 0 || mpz_sgn(den) <= 0) {
            mpz_clear(num); mpz_clear(den);
            return false;
        }
        return true;
    }
    return false;
}

/* ----------------------------------------------------------------------- */
/* Pass A: Log[positive rational] prime decomposition                      */
/* ----------------------------------------------------------------------- */

/*
 * Build the decomposition Sum e_i * Log[p_i] for a single positive
 * integer n > 1 using FactorInteger. Returns NULL for n <= 1, or when the
 * factorisation has only one prime with exponent 1 (which is just
 * Log[p_i] itself -- no decomposition).
 *
 * The returned expression is evaluated (and so canonicalised).
 */
static Expr* decompose_log_integer(const mpz_t n) {
    if (mpz_cmp_ui(n, 1) <= 0) return NULL;

    Expr* n_expr = expr_new_bigint_from_mpz(n);
    Expr* call = make_call1("FactorInteger", n_expr);
    Expr* factors = eval_take(call);
    if (!factors ||
        factors->type != EXPR_FUNCTION ||
        !factors->data.function.head ||
        factors->data.function.head->type != EXPR_SYMBOL ||
        factors->data.function.head->data.symbol != SYM_List) {
        if (factors) expr_free(factors);
        return NULL;
    }

    size_t k = factors->data.function.arg_count;
    /* Single prime with exponent 1: no decomposition. */
    if (k == 1) {
        Expr* pair = factors->data.function.args[0];
        if (pair && pair->type == EXPR_FUNCTION &&
            pair->data.function.head &&
            pair->data.function.head->type == EXPR_SYMBOL &&
            pair->data.function.head->data.symbol == SYM_List &&
            pair->data.function.arg_count == 2) {
            Expr* e = pair->data.function.args[1];
            if (e && e->type == EXPR_INTEGER && e->data.integer == 1) {
                expr_free(factors);
                return NULL;
            }
        }
    }

    /* Build Plus[Times[e_i, Log[p_i]], ...] */
    Expr** terms = (Expr**)calloc(k ? k : 1, sizeof(Expr*));
    size_t emitted = 0;
    for (size_t i = 0; i < k; i++) {
        Expr* pair = factors->data.function.args[i];
        if (!pair || pair->type != EXPR_FUNCTION ||
            pair->data.function.arg_count != 2) continue;
        Expr* p = expr_copy(pair->data.function.args[0]);
        Expr* e = expr_copy(pair->data.function.args[1]);
        Expr* logp = make_call1("Log", p);
        Expr* term = make_call2("Times", e, logp);
        terms[emitted++] = term;
    }
    expr_free(factors);

    if (emitted == 0) {
        free(terms);
        return NULL;
    }
    Expr* plus = expr_new_function(expr_new_symbol(SYM_Plus), terms, emitted);
    free(terms);
    return eval_take(plus);
}

/*
 * Decompose Log[num/den] = Log[num] - Log[den] into prime logs.
 * Returns NULL when both num and den are 1 (Log[1] = 0; handled by the
 * Log builtin already) or when the decomposition is a no-op (single
 * Log[prime] form).
 */
static Expr* decompose_log_rational(const mpz_t num, const mpz_t den) {
    Expr* part_num = NULL;
    Expr* part_den = NULL;
    bool num_changed = false, den_changed = false;
    bool num_one = (mpz_cmp_ui(num, 1) == 0);
    bool den_one = (mpz_cmp_ui(den, 1) == 0);

    if (!num_one) {
        part_num = decompose_log_integer(num);
        num_changed = (part_num != NULL);
        if (!part_num) {
            /* Single prime: just Log[num]. */
            Expr* n_expr = expr_new_bigint_from_mpz(num);
            part_num = make_call1("Log", n_expr);
            part_num = eval_take(part_num);
        }
    }
    if (!den_one) {
        Expr* dec = decompose_log_integer(den);
        den_changed = (dec != NULL);
        if (!dec) {
            Expr* d_expr = expr_new_bigint_from_mpz(den);
            dec = make_call1("Log", d_expr);
            dec = eval_take(dec);
        }
        /* Negate: -Log[den] */
        part_den = eval_take(make_call2("Times",
                                        expr_new_integer(-1), dec));
    }

    /* Combine. */
    Expr* result;
    if (part_num && part_den) {
        result = eval_take(make_call2("Plus", part_num, part_den));
    } else if (part_num) {
        result = part_num;
    } else if (part_den) {
        result = part_den;
    } else {
        return NULL;
    }

    /* Bail when nothing decomposed: a single Log[prime] / Log[1/prime]
     * with no internal structure. The caller wants NULL for "no change". */
    if (!num_changed && !den_changed) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/*
 * Try to decompose a Log[..] expression. Returns NULL when the argument
 * isn't a positive rational or a Power of one. Otherwise returns the
 * decomposed sum.
 */
static Expr* try_decompose_one_log(const Expr* log_expr) {
    if (!is_log1(log_expr)) return NULL;
    Expr* arg = log_expr->data.function.args[0];

    /* Direct positive-rational arg. */
    mpz_t num, den;
    if (extract_pos_rational(arg, num, den)) {
        Expr* dec = decompose_log_rational(num, den);
        mpz_clear(num); mpz_clear(den);
        return dec;
    }

    /* Log[Power[positive_rational, expt]] -> expt * decomposed Log[r]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2) {
        Expr* base = arg->data.function.args[0];
        Expr* expt = arg->data.function.args[1];
        mpz_t bp, bq;
        if (extract_pos_rational(base, bp, bq)) {
            Expr* dec = decompose_log_rational(bp, bq);
            mpz_clear(bp); mpz_clear(bq);
            if (!dec) {
                /* Single-prime base: emit expt * Log[base]. */
                Expr* base_log = eval_take(make_call1("Log", expr_copy(base)));
                return eval_take(make_call2("Times",
                                            expr_copy(expt), base_log));
            }
            return eval_take(make_call2("Times", expr_copy(expt), dec));
        }
    }
    return NULL;
}

/*
 * Bottom-up walk applying try_decompose_one_log to every Log[..] subexpr.
 * Returns NULL when nothing changes; otherwise a freshly owned, evaluated
 * tree.
 */
static Expr* pass_a_decompose(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    /* Recurse into children. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = pass_a_decompose(e->data.function.args[i]);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n, sizeof(Expr*));
                for (size_t j = 0; j < i; j++)
                    new_args[j] = expr_copy(e->data.function.args[j]);
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current = NULL;
    const Expr* target;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = eval_take(rebuilt);
        target = current;
    } else {
        target = e;
    }

    /* Apply decomposition at this level if applicable. */
    Expr* top = try_decompose_one_log(target);
    if (top) {
        if (current) expr_free(current);
        return top;
    }
    return current;
}

/* ----------------------------------------------------------------------- */
/* Pass B: linear-combination-of-logs fuser                                */
/* ----------------------------------------------------------------------- */

/*
 * If `term` is one of:
 *     Log[a]                       -> coeff = 1, arg = a
 *     Times[c, Log[a]]             -> coeff = c, arg = a (c numeric or symbolic)
 *     Times[c1, c2, ..., Log[a]]   -> coeff = Times[c1, c2, ...], arg = a
 * write the (borrowed) coefficient and argument into *coeff_out / *arg_out
 * and return true. The caller does NOT own *coeff_out -- it's either a
 * freshly built (caller-owned) literal or a freshly built Times node;
 * in either case the caller must eventually free it.
 * Returns false when the term contains zero or multiple Log[..] factors,
 * or when any Log has the wrong arity.
 */
static bool extract_log_term(const Expr* term, Expr** coeff_out, Expr** arg_out) {
    if (is_log1(term)) {
        *coeff_out = expr_new_integer(1);
        *arg_out   = expr_copy(term->data.function.args[0]);
        return true;
    }
    if (!term || term->type != EXPR_FUNCTION ||
        !term->data.function.head ||
        term->data.function.head->type != EXPR_SYMBOL ||
        term->data.function.head->data.symbol != SYM_Times) {
        return false;
    }
    size_t n = term->data.function.arg_count;
    int log_idx = -1;
    for (size_t i = 0; i < n; i++) {
        if (is_log1(term->data.function.args[i])) {
            if (log_idx >= 0) return false; /* more than one Log */
            log_idx = (int)i;
        }
    }
    if (log_idx < 0) return false;

    /* Coefficient is the product of all other factors. */
    Expr** rest = (Expr**)calloc(n > 1 ? n - 1 : 1, sizeof(Expr*));
    size_t emitted = 0;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == log_idx) continue;
        rest[emitted++] = expr_copy(term->data.function.args[i]);
    }
    Expr* coeff;
    if (emitted == 0) {
        coeff = expr_new_integer(1);
        free(rest);
    } else if (emitted == 1) {
        coeff = rest[0];
        free(rest);
    } else {
        coeff = eval_take(expr_new_function(expr_new_symbol(SYM_Times),
                                            rest, emitted));
        free(rest);
    }
    *coeff_out = coeff;
    *arg_out   = expr_copy(term->data.function.args[log_idx]->data.function.args[0]);
    return true;
}

/* True iff `e` has a leaf-count complexity below `cap` -- proxy for
 * "collapsed to a near-atomic form" used to decide whether a fuse
 * is safe to take without provable positivity. */
static bool is_small_constant(const Expr* e, size_t cap) {
    return simp_default_complexity(e) <= cap;
}

/*
 * Try to fuse the log block of a single Plus[t1, ..., tn]. Returns NULL
 * when no useful fuse fires; otherwise a freshly-owned evaluated Plus
 * with the fused Log block in place of the original log terms.
 */
static Expr* try_fuse_plus(const Expr* plus_expr, const AssumeCtx* ctx) {
    if (!plus_expr || plus_expr->type != EXPR_FUNCTION ||
        !plus_expr->data.function.head ||
        plus_expr->data.function.head->type != EXPR_SYMBOL ||
        plus_expr->data.function.head->data.symbol != SYM_Plus) {
        return NULL;
    }
    size_t n = plus_expr->data.function.arg_count;
    if (n < 2) return NULL;

    Expr** coeffs = (Expr**)calloc(n, sizeof(Expr*));
    Expr** args   = (Expr**)calloc(n, sizeof(Expr*));
    Expr** rest   = (Expr**)calloc(n, sizeof(Expr*));
    size_t n_log = 0, n_rest = 0;
    bool all_args_positive = true;

    for (size_t i = 0; i < n; i++) {
        Expr* t = plus_expr->data.function.args[i];
        Expr* c = NULL; Expr* a = NULL;
        if (extract_log_term(t, &c, &a)) {
            coeffs[n_log] = c;
            args[n_log]   = a;
            if (!assume_known_positive(ctx, a)) all_args_positive = false;
            n_log++;
        } else {
            rest[n_rest++] = expr_copy(t);
        }
    }

    if (n_log < 2) {
        for (size_t i = 0; i < n_log; i++) {
            expr_free(coeffs[i]); expr_free(args[i]);
        }
        for (size_t i = 0; i < n_rest; i++) expr_free(rest[i]);
        free(coeffs); free(args); free(rest);
        return NULL;
    }

    /* Build product: Product[args[i] ^ coeffs[i]]. */
    Expr** factors = (Expr**)calloc(n_log, sizeof(Expr*));
    for (size_t i = 0; i < n_log; i++) {
        factors[i] = make_call2("Power", args[i], coeffs[i]);
    }
    /* args[] and coeffs[] are consumed by the Power factors above. */
    free(coeffs); free(args);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), factors, n_log);
    free(factors);
    Expr* prod_eval = eval_take(prod);

    /* Surface algebraic cancellations inside the fused argument so the
     * gate downstream can recognise a collapsed form.
     *
     *   - Together folds (x^2 - y^2) / ((x-y)(x+y)) -> 1 (algebraic
     *     test cases 6, 8).
     *   - Expand opens (Sec-Tan)(Sec+Tan) into Sec^2 - Tan^2 so the
     *     Pythagorean Simplify pass can fire on it.
     *   - A one-shot Simplify of the inner argument lets identities
     *     like Sec^2 - Tan^2 -> 1 collapse the fused Log to Log[1]=0.
     *     Recursion into Pass B is blocked via simp_log_depth so this
     *     never re-enters the same Plus. */
    Expr* tg = eval_take(make_call1("Together", expr_copy(prod_eval)));
    if (tg && simp_default_complexity(tg) < simp_default_complexity(prod_eval)) {
        expr_free(prod_eval);
        prod_eval = tg;
    } else if (tg) {
        expr_free(tg);
    }
    Expr* ex = eval_take(make_call1("Expand", expr_copy(prod_eval)));
    if (ex && simp_default_complexity(ex) < simp_default_complexity(prod_eval)) {
        expr_free(prod_eval);
        prod_eval = ex;
    } else if (ex) {
        expr_free(ex);
    }
    /* Recursive Simplify on the argument only -- bounded by simp_log_depth
     * (defined in the driver) and by the strict structural decrease from
     * the original Plus to its log-block argument. */
    extern int simp_log_depth; /* see driver below */
    if (simp_log_depth <= 1) {
        Expr* si = eval_take(make_call1("Simplify", expr_copy(prod_eval)));
        if (si &&
            simp_default_complexity(si) < simp_default_complexity(prod_eval)) {
            expr_free(prod_eval);
            prod_eval = si;
        } else if (si) {
            expr_free(si);
        }
    }

    /* Build new Plus[Log[prod_eval], rest...] and evaluate. */
    Expr* fused_log = eval_take(make_call1("Log", prod_eval));

    Expr** all = (Expr**)calloc(n_rest + 1, sizeof(Expr*));
    all[0] = fused_log;
    for (size_t i = 0; i < n_rest; i++) all[i + 1] = rest[i];
    free(rest);

    Expr* new_plus;
    if (n_rest + 1 == 1) {
        new_plus = all[0];
        free(all);
    } else {
        new_plus = expr_new_function(expr_new_symbol(SYM_Plus), all, n_rest + 1);
        free(all);
        new_plus = eval_take(new_plus);
    }

    /* Decide whether to take the fusion.
     *
     *   - Provably positive args: accept any strict improvement.
     *   - Otherwise: accept only when the result collapsed to a small
     *     constant (leaf count <= 4) -- e.g. 0, Log[2], -2 x. This is
     *     the safe gate that prevents speculative fusion from being
     *     "wrong but smaller" on inputs whose log args might be
     *     negative or complex. */
    size_t old_score = simp_default_complexity(plus_expr);
    size_t new_score = simp_default_complexity(new_plus);
    bool take;
    if (all_args_positive) {
        take = new_score < old_score;
    } else {
        take = new_score < old_score && is_small_constant(new_plus, 4);
    }
    if (!take) {
        expr_free(new_plus);
        return NULL;
    }
    return new_plus;
}

/*
 * Bottom-up walk applying try_fuse_plus to every Plus subexpr.
 */
static Expr* pass_b_fuse(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = pass_b_fuse(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n, sizeof(Expr*));
                for (size_t j = 0; j < i; j++)
                    new_args[j] = expr_copy(e->data.function.args[j]);
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current = NULL;
    const Expr* target;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = eval_take(rebuilt);
        target = current;
    } else {
        target = e;
    }

    Expr* top = try_fuse_plus(target, ctx);
    if (top) {
        if (current) expr_free(current);
        return top;
    }
    return current;
}

/* ----------------------------------------------------------------------- */
/* Driver -- iterate the two passes to a fixed point.                      */
/* ----------------------------------------------------------------------- */

/* Recursion guard. Pass B may call Simplify on the fused argument, which
 * routes back through Simplify -> simp_dispatch -> simp_log_apply. The
 * argument is strictly smaller than the parent Plus so termination is
 * guaranteed structurally, but we cap depth at 2 to avoid unbounded work
 * on pathological inputs. Single-threaded by design (matches the rest of
 * the simp module). */
int simp_log_depth = 0;

/* Cheap pre-check: skip the walker when the input has no Log subexpr. */
static bool has_log(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Log) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_log(e->data.function.args[i])) return true;
    }
    return false;
}

Expr* simp_log_apply(const Expr* e, const AssumeCtx* ctx) {
    if (!e) return NULL;
    if (!has_log(e)) return NULL;
    if (simp_log_depth >= 2) return NULL;

    simp_log_depth++;
    Expr* cur = expr_copy((Expr*)e);
    bool any_change = false;

    /* Bounded fixed-point loop. Each pass either makes progress or
     * leaves cur unchanged. Typical inputs converge in 1-2 iterations;
     * 6 is a safety cap, not a typical iteration count. */
    for (int iter = 0; iter < 6; iter++) {
        bool changed = false;

        Expr* a = pass_a_decompose(cur);
        if (a) {
            if (!expr_eq(a, cur)) {
                expr_free(cur);
                cur = a;
                changed = true;
                any_change = true;
            } else {
                expr_free(a);
            }
        }

        Expr* b = pass_b_fuse(cur, ctx);
        if (b) {
            if (!expr_eq(b, cur)) {
                expr_free(cur);
                cur = b;
                changed = true;
                any_change = true;
            } else {
                expr_free(b);
            }
        }
        if (!changed) break;
    }

    simp_log_depth--;
    if (!any_change || expr_eq(cur, e)) {
        expr_free(cur);
        return NULL;
    }
    return cur;
}

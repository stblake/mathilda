/* intsimp.c — simplification helpers for the rational-function
 * integrator (`Integrate`` package).  See intsimp.h for the public
 * surface and the rationale for the split out of intrat.c.
 *
 * Two layers of helper live here:
 *   - Algorithmic simplification used inside the integrator
 *     (Sqrt under positive-symbol assumption, sign tests,
 *      canonic-zero test, radical detection / Simplify guard).
 *   - Output cleanup applied to the resulting integral
 *     (Log pairing, Log[c·p] -> Log[p], Plus distribution,
 *      ArcTan / ArcTanh sign normalisation).
 *
 * Memory contract follows the standard Mathilda convention: every
 * helper returns a freshly-allocated Expr* the caller owns; none of
 * them free their input arguments unless explicitly noted in
 * intsimp.h.
 */

#include "intsimp.h"
#include "intrat_internal.h"
#include "expr.h"
#include "eval.h"
#include "internal.h"
#include "expand.h"
#include "poly.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ====================================================================
 * Radical detection + Simplify guard.
 * ====================================================================
 *
 * Used by the LRT pipeline's primitivePart step: when the residue r
 * carries an unsimplified radical (e.g. an Sqrt[…] that survived
 * extEuc), the downstream PolynomialRemainder over an algebraic
 * extension fails to terminate.  Pre-applying Simplify[] reduces the
 * radical-coefficient terms in the underlying algebraic field so the
 * pseudo-division loop converges.
 */

bool intsimp_has_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol == SYM_Rational) {
            return true;
        }
    }
    if (intsimp_has_radical(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (intsimp_has_radical(e->data.function.args[i])) return true;
    }
    return false;
}

Expr* intsimp_simplify_if_radical(Expr* e) {
    if (!e || !intsimp_has_radical(e)) return e;
    Expr* call = expr_new_function(
        expr_new_symbol("Simplify"), (Expr*[]){ e }, 1);
    Expr* simplified = evaluate(call);
    expr_free(call);
    return simplified;
}

/* ====================================================================
 * Sign tests.
 * ====================================================================
 *
 * intsimp_sign_pos_assumption: returns +1 / -1 / 0 for the sign of
 * `e` under the assumption that every free symbol denotes a positive
 * real.  Used by the LogToReal dispatcher so a quadratic factor with
 * a symbolic-but-provably-signed discriminant (e.g. -4 a^2 b^2 with
 * a, b parametric) can take the appropriate real / complex branch
 * instead of falling through to NaiveLogPart's held RootSum form.
 *
 * intsimp_numeric_sign: complement that decides the sign by
 * numerical evaluation through `N[e]`.  Handles `Sqrt[5] - 5` and
 * similar mixed-sign Plus aggregates that the symbolic walker leaves
 * at sign-unknown.  Returns +1 / -1 only when N produces a
 * definite-sign Real / Integer; 0 otherwise.
 */

int intsimp_sign_pos_assumption(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t n = e->data.function.args[0]->data.integer;
        int64_t d = e->data.function.args[1]->data.integer;
        int sign = 1;
        if (n < 0) { n = -n; sign = -sign; }
        if (d < 0) { d = -d; sign = -sign; }
        if (n == 0) return 0;
        return sign;
    }
    if (e->type == EXPR_SYMBOL) {
        /* Free symbols treated as positive reals. */
        return 1;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp = e->data.function.args[1];
        int b_sign = intsimp_sign_pos_assumption(base);
        if (exp->type == EXPR_INTEGER) {
            int64_t k = exp->data.integer;
            if (k % 2 == 0 && b_sign != 0) return 1;
            return b_sign;
        }
        if (b_sign > 0) return 1;
        return 0;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        int sign = 1;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int s = intsimp_sign_pos_assumption(e->data.function.args[i]);
            if (s == 0) return 0;
            sign *= s;
        }
        return sign;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Plus) {
        int sign = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int s = intsimp_sign_pos_assumption(e->data.function.args[i]);
            if (s == 0) return 0;
            if (sign == 0) sign = s;
            else if (sign != s) return 0;
        }
        return sign;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Sqrt
        && e->data.function.arg_count == 1) {
        int s = intsimp_sign_pos_assumption(e->data.function.args[0]);
        if (s > 0) return 1;
        return 0;
    }
    return 0;
}

int intsimp_numeric_sign(Expr* e) {
    if (!e) return 0;
    Expr* call = expr_new_function(expr_new_symbol("N"),
        (Expr*[]){expr_copy(e)}, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    int sign = 0;
    if (r) {
        if (r->type == EXPR_REAL) {
            double v = r->data.real;
            /* 1e-12 dead zone keeps Sqrt[2] - 1.4142... = 1.36e-13
             * style round-off out of the wrong sign bin. */
            if (v >  1e-12) sign =  1;
            else if (v < -1e-12) sign = -1;
        } else if (r->type == EXPR_INTEGER) {
            if (r->data.integer >  0) sign =  1;
            else if (r->data.integer <  0) sign = -1;
        }
        expr_free(r);
    }
    return sign;
}

/* ====================================================================
 * Canonic-zero test.
 * ====================================================================
 * zeroQ[e] = TrueQ[Cancel[Together[e]] === 0].
 */

bool intsimp_zero_q(Expr* e) {
    Expr* canon = intrat_canonic(e);
    bool ok = is_zero_poly(canon);
    expr_free(canon);
    return ok;
}

/* ====================================================================
 * Sqrt simplification under the positive-symbol assumption.
 * ====================================================================
 *
 * Returns the simplification of `Sqrt[e]` under the assumption that
 * every free symbol denotes a positive real, or a fresh `Sqrt[e]` if
 * no factors can be extracted.  The caller owns the returned tree.
 * Walks `e` factor-by-factor:
 *   integer n           -> peel out the largest perfect-square divisor
 *   Rational p/q        -> recurse on each part
 *   Power[base, k]      -> base^(k/2) when k is even and base > 0,
 *                          base^((k-1)/2) Sqrt[base] for odd positive k
 *   Times[f1, f2, ...]  -> product of recursive results
 *   anything else       -> Sqrt[anything] (irreducible)
 */

Expr* intsimp_pos_sqrt_factor(Expr* e) {
    if (e->type == EXPR_INTEGER) {
        int64_t v = e->data.integer;
        if (v < 0) {
            /* Sqrt of a negative integer: leave inside Sqrt. */
            return expr_new_function(expr_new_symbol("Sqrt"),
                (Expr*[]){expr_copy(e)}, 1);
        }
        if (v == 0) return expr_new_integer(0);
        if (v == 1) return expr_new_integer(1);
        int64_t s = 1, r = v;
        for (int64_t p = 2; p * p <= r; p++) {
            while (r % (p * p) == 0) {
                s *= p;
                r /= p * p;
            }
        }
        if (s == 1) {
            return expr_new_function(expr_new_symbol("Sqrt"),
                (Expr*[]){expr_new_integer(r)}, 1);
        }
        if (r == 1) return expr_new_integer(s);
        Expr* sqrt_part = expr_new_function(expr_new_symbol("Sqrt"),
            (Expr*[]){expr_new_integer(r)}, 1);
        return eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(s), sqrt_part}, 2));
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        Expr* num = e->data.function.args[0];
        Expr* den = e->data.function.args[1];
        Expr* sn = intsimp_pos_sqrt_factor(num);
        Expr* sd = intsimp_pos_sqrt_factor(den);
        return eval_and_free(internal_divide((Expr*[]){sn, sd}, 2));
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp_e = e->data.function.args[1];
        if (exp_e->type == EXPR_INTEGER) {
            int64_t k = exp_e->data.integer;
            int b_sign = intsimp_sign_pos_assumption(base);
            if (k % 2 == 0 && b_sign > 0) {
                return expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){expr_copy(base),
                              expr_new_integer(k / 2)}, 2);
            }
            if (b_sign > 0) {
                /* Odd k, positive base: base^((k-1)/2) * Sqrt[base]. */
                int64_t kh = (k - 1) / 2;  /* may be negative for k < 0 */
                Expr* p1;
                if (kh == 0) p1 = expr_new_integer(1);
                else p1 = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){expr_copy(base), expr_new_integer(kh)}, 2);
                Expr* p2 = expr_new_function(expr_new_symbol("Sqrt"),
                    (Expr*[]){expr_copy(base)}, 1);
                return eval_and_free(internal_times(
                    (Expr*[]){p1, p2}, 2));
            }
        }
        return expr_new_function(expr_new_symbol("Sqrt"),
            (Expr*[]){expr_copy(e)}, 1);
    }
    if (e->type == EXPR_SYMBOL) {
        /* Sqrt[symbol] is irreducible — keep wrapped. */
        return expr_new_function(expr_new_symbol("Sqrt"),
            (Expr*[]){expr_copy(e)}, 1);
    }
    /* Fallback: no structural simplification. */
    return expr_new_function(expr_new_symbol("Sqrt"),
        (Expr*[]){expr_copy(e)}, 1);
}

Expr* intsimp_pos_sqrt(Expr* e) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        size_t n = e->data.function.arg_count;
        Expr** out = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            out[i] = intsimp_pos_sqrt_factor(e->data.function.args[i]);
        }
        Expr* prod = internal_times(out, n);
        free(out);
        return eval_and_free(prod);
    }
    return intsimp_pos_sqrt_factor(e);
}

/* ====================================================================
 * Output cleanup — Log pairing.
 * ====================================================================
 *
 * Combines `c Log[A] + c Log[B] -> c Log[A B]` and
 * `c Log[A] - c Log[B] -> c Log[A/B]` into single logarithms, then
 * recognises the `c Log[A] - c Log[B] -> 2 c ArcTanh[(B-A)/(B+A)]`
 * pattern when the ArcTanh argument simplifies to a rational
 * function in x.  Direct port of IntegrateRational.m:1722-1761
 * (LogToArcTanh) and :1896-1958 (LogToArcTan).
 *
 * Implemented as direct C transformations on Plus[...] of Log[...]
 * terms (per the plan note) rather than pattern-rewriting at the
 * Mathilda rule-engine layer.  The transformations are
 * differentiation-equivalent — they only beautify the output — so
 * the universal correctness check stays green either way.
 */

/* Walk a Times tree (handling unflattened nested Times produced by e.g.
 * builtin_subtract building Plus[a, Times[-1, b]] without re-evaluating)
 * and append each non-Times leaf into `factors`.  Returns true on success;
 * `*log_arg_out` is set to the inner argument of the (unique) Log factor
 * encountered, or left NULL if none was found. Returns false if more than
 * one Log factor is found. */
static bool collect_times_factors(Expr* term,
                                  Expr*** factors, size_t* fc, size_t* fcap,
                                  Expr** log_arg_out) {
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Times) {
        size_t n = term->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            if (!collect_times_factors(term->data.function.args[i],
                                       factors, fc, fcap, log_arg_out)) {
                return false;
            }
        }
        return true;
    }
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Log
        && term->data.function.arg_count == 1) {
        if (*log_arg_out) return false; /* multiple Log factors */
        *log_arg_out = term->data.function.args[0];
        return true;
    }
    if (*fc == *fcap) {
        *fcap = (*fcap) ? (*fcap) * 2 : 4;
        *factors = (Expr**)realloc(*factors, sizeof(Expr*) * (*fcap));
    }
    (*factors)[(*fc)++] = term;
    return true;
}

/* Decompose a term of a Plus head into (coeff, log_arg) when it
 * matches `coeff * Log[log_arg]`.  Returns true on success. */
static bool decompose_log_term(Expr* term, Expr** coeff_out, Expr** log_arg_out) {
    /* Plain Log[x]. */
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Log
        && term->data.function.arg_count == 1) {
        *coeff_out = expr_new_integer(1);
        *log_arg_out = expr_copy(term->data.function.args[0]);
        return true;
    }
    /* Times-tree (possibly unflattened) containing one Log factor; the
     * rest become the coefficient. */
    if (term->type != EXPR_FUNCTION
        || term->data.function.head->type != EXPR_SYMBOL
        || term->data.function.head->data.symbol != SYM_Times) return false;

    Expr** factors = NULL;
    size_t fc = 0, fcap = 0;
    Expr* log_arg = NULL;
    if (!collect_times_factors(term, &factors, &fc, &fcap, &log_arg)) {
        free(factors);
        return false;
    }
    if (!log_arg) { free(factors); return false; }

    /* Build coeff = product of all collected non-Log factors. */
    Expr* coeff;
    if (fc == 0) {
        coeff = expr_new_integer(1);
    } else if (fc == 1) {
        coeff = expr_copy(factors[0]);
    } else {
        Expr** copies = (Expr**)malloc(sizeof(Expr*) * fc);
        for (size_t i = 0; i < fc; i++) copies[i] = expr_copy(factors[i]);
        coeff = eval_and_free(internal_times(copies, fc));
        free(copies);
    }
    free(factors);
    *coeff_out = coeff;
    *log_arg_out = expr_copy(log_arg);
    return true;
}

Expr* intsimp_log_to_arctanh(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (!(e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Plus)) {
        return expr_copy(e);
    }

    size_t n = e->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    bool* used = (bool*)calloc(n, sizeof(bool));
    Expr** coeffs = (Expr**)calloc(n, sizeof(Expr*));
    Expr** logargs = (Expr**)calloc(n, sizeof(Expr*));
    bool* is_log = (bool*)calloc(n, sizeof(bool));

    for (size_t i = 0; i < n; i++) {
        terms[i] = e->data.function.args[i];
        is_log[i] = decompose_log_term(terms[i], &coeffs[i], &logargs[i]);
    }

    Expr** out = (Expr**)malloc(sizeof(Expr*) * (n * 2));
    size_t out_n = 0;

    for (size_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (!is_log[i]) {
            out[out_n++] = expr_copy(terms[i]);
            used[i] = true;
            continue;
        }
        /* Try to pair with another Log term j > i. */
        bool merged = false;
        for (size_t j = i + 1; j < n; j++) {
            if (used[j] || !is_log[j]) continue;
            if (!intrat_freeq_test(coeffs[i], x) || !intrat_freeq_test(coeffs[j], x)) continue;

            /* Same coefficient -> Log[A*B]. */
            Expr* delta = internal_subtract(
                (Expr*[]){expr_copy(coeffs[i]), expr_copy(coeffs[j])}, 2);
            if (intsimp_zero_q(delta)) {
                expr_free(delta);
                Expr* prod = expr_expand(internal_times(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2));
                Expr* logp = expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){prod}, 1);
                Expr* combined = internal_times(
                    (Expr*[]){expr_copy(coeffs[i]), logp}, 2);
                out[out_n++] = eval_and_free(combined);
                used[i] = used[j] = true; merged = true; break;
            }
            expr_free(delta);

            /* Opposite coefficient: try ArcTanh first, then Log[A/B]. */
            Expr* sumcoef = internal_plus(
                (Expr*[]){expr_copy(coeffs[i]), expr_copy(coeffs[j])}, 2);
            if (intsimp_zero_q(sumcoef)) {
                expr_free(sumcoef);
                /* Mathematica's log2ArcTanhRule emits
                 *   (c2 - c1) ArcTanh[(A + B) / (B - A)]
                 * for c1 Log[A] + c2 Log[B] with c1 + c2 == 0.
                 * Equivalent to 2 c2 ArcTanh[…].  We require the
                 * ArcTanh argument's Denominator to be free of x; if
                 * not, fall back to the simpler Log[A/B] rewrite. */
                Expr* sumAB_raw = internal_plus(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2);
                Expr* sumAB = expr_expand(sumAB_raw); expr_free(sumAB_raw);
                Expr* diffBA_raw = internal_subtract(
                    (Expr*[]){expr_copy(logargs[j]), expr_copy(logargs[i])}, 2);
                Expr* diffBA = expr_expand(diffBA_raw); expr_free(diffBA_raw);

                Expr* arg_raw = internal_divide(
                    (Expr*[]){expr_copy(sumAB), expr_copy(diffBA)}, 2);
                Expr* arg_can = intrat_canonic(arg_raw); expr_free(arg_raw);
                Expr* arg_den = intrat_denominator(arg_can);
                Expr* arg_den_eval = eval_and_free(arg_den);
                bool atanh_ok = intrat_freeq_test(arg_den_eval, x)
                              && !is_zero_poly(arg_can);
                expr_free(arg_den_eval);
                if (atanh_ok) {
                    /* (c2 - c1) ArcTanh[(A + B) / (B - A)] */
                    Expr* coef_diff = internal_subtract(
                        (Expr*[]){expr_copy(coeffs[j]), expr_copy(coeffs[i])}, 2);
                    coef_diff = eval_and_free(coef_diff);
                    Expr* atanh = expr_new_function(expr_new_symbol("ArcTanh"),
                        (Expr*[]){arg_can}, 1);
                    Expr* term = internal_times(
                        (Expr*[]){coef_diff, atanh}, 2);
                    out[out_n++] = eval_and_free(term);
                    expr_free(sumAB); expr_free(diffBA);
                    used[i] = used[j] = true; merged = true; break;
                }
                expr_free(arg_can);

                /* Log[A/B] with denominator free of x. */
                Expr* AoverB_raw = internal_divide(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2);
                Expr* AoverB = intrat_canonic(AoverB_raw); expr_free(AoverB_raw);
                Expr* AoverB_den = intrat_denominator(AoverB);
                Expr* AoverB_den_eval = eval_and_free(AoverB_den);
                bool divlog_ok = intrat_freeq_test(AoverB_den_eval, x);
                expr_free(AoverB_den_eval);
                if (divlog_ok) {
                    Expr* logp = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){AoverB}, 1);
                    Expr* term = internal_times(
                        (Expr*[]){expr_copy(coeffs[i]), logp}, 2);
                    out[out_n++] = eval_and_free(term);
                    expr_free(sumAB); expr_free(diffBA);
                    used[i] = used[j] = true; merged = true; break;
                }
                expr_free(AoverB);
                expr_free(sumAB); expr_free(diffBA);
                continue;
            }
            expr_free(sumcoef);
        }
        if (!merged) {
            out[out_n++] = expr_copy(terms[i]);
            used[i] = true;
        }
    }

    for (size_t i = 0; i < n; i++) {
        if (coeffs[i])  expr_free(coeffs[i]);
        if (logargs[i]) expr_free(logargs[i]);
    }
    free(coeffs); free(logargs); free(is_log); free(used); free(terms);

    Expr* result;
    if (out_n == 0) { free(out); result = expr_new_integer(0); }
    else if (out_n == 1) { result = out[0]; free(out); }
    else { result = eval_and_free(internal_plus(out, out_n)); free(out); }
    return result;
}

/* ====================================================================
 * Output cleanup — Log[c · p] -> Log[p] when c is free of x.
 * ====================================================================
 *
 * Walks the input top-down and rewrites every Log subexpression
 * whose argument is a Times-headed product with at least one
 * constant (FreeQ[…, x]) factor.  The additive `Log[c]` term is
 * harmless to the antiderivative — it is a constant of integration —
 * so dropping it merely beautifies the result and keeps the
 * downstream `intsimp_log_to_arctanh` pairing rule able to fire (it
 * requires `Log[Cancel[A/B]]` denominator to be free of x, which a
 * stray `Log[1/2·…]` factor inside would block).
 *
 * Direct mirror of IntegrateRational.m:1746 / :1750.
 */

Expr* intsimp_strip_log_constants(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    /* Recurse first so nested Logs get cleaned bottom-up. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = intsimp_strip_log_constants(e->data.function.args[i], x);
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* rebuilt = expr_new_function(head, new_args, n);
    free(new_args);

    /* Only Log heads need post-rewrite. */
    if (rebuilt->type != EXPR_FUNCTION
        || rebuilt->data.function.head->type != EXPR_SYMBOL
        || rebuilt->data.function.head->data.symbol != SYM_Log
        || rebuilt->data.function.arg_count != 1) {
        return rebuilt;
    }

    Expr* arg = rebuilt->data.function.args[0];
    if (arg->type != EXPR_FUNCTION
        || arg->data.function.head->type != EXPR_SYMBOL
        || arg->data.function.head->data.symbol != SYM_Times) {
        return rebuilt;
    }

    /* Partition arg's Times factors into (const, x-dependent). */
    size_t af = arg->data.function.arg_count;
    Expr** keep = (Expr**)malloc(sizeof(Expr*) * af);
    size_t nkeep = 0;
    bool any_dropped = false;
    for (size_t i = 0; i < af; i++) {
        Expr* fac = arg->data.function.args[i];
        if (intrat_freeq_test(fac, x)) {
            any_dropped = true;
        } else {
            keep[nkeep++] = expr_copy(fac);
        }
    }
    if (!any_dropped || nkeep == 0) {
        for (size_t i = 0; i < nkeep; i++) expr_free(keep[i]);
        free(keep);
        return rebuilt;
    }

    Expr* new_inner;
    if (nkeep == 1) { new_inner = keep[0]; }
    else            { new_inner = eval_and_free(internal_times(keep, nkeep)); }
    free(keep);

    Expr* new_log = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){new_inner}, 1);
    expr_free(rebuilt);
    return eval_and_free(new_log);
}

/* ====================================================================
 * Output cleanup — distribute scalar Times across Plus.
 * ====================================================================
 *
 * After per-summand `c_k * piece_int_k` scaling, the result frequently
 * shows up as
 *   Times[c_k, Plus[term1, term2, …]]
 * because Times has no auto-distribution attribute.  A single
 * `expr_expand` pass turns it into `c_k term1 + c_k term2`, which is
 * what the `Collect[expanded, _Log|_ArcTan|_ArcTanh, …]` step in
 * IntegrateRational.m:99 produces.  Mathematica's Distribute / Expand
 * has the same effect; we already have expr_expand which threads
 * Times-over-Plus polynomial-style, which is exactly what we want.
 * Calling this twice (before and after intsimp_log_to_arctanh) is
 * cheap and lets the log-pairing rule see fully-distributed sums.
 */
Expr* intsimp_distribute_plus(Expr* e) {
    if (!e) return NULL;
    return expr_expand(e);
}

/* ====================================================================
 * Output cleanup — ArcTan / ArcTanh sign normalisation.
 * ====================================================================
 *
 * Pull a leading minus out of the argument so the printed form is
 * canonical.  Walks the input expression top-down, rewriting just the
 * relevant heads.
 */
Expr* intsimp_normalize_inverse_trig_signs(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);

    /* ArcTan[-(...)] -> -ArcTan[...], same for ArcTanh.  We look for
     * an argument that starts with a unary minus or a Times whose
     * first factor is a literal -1. */
    if (e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol == SYM_ArcTan
            || e->data.function.head->data.symbol == SYM_ArcTanh)
        && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        bool negative = false;
        Expr* stripped = NULL;
        if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
            negative = true;
            stripped = expr_new_integer(-arg->data.integer);
        } else if (arg->type == EXPR_FUNCTION
            && arg->data.function.head->type == EXPR_SYMBOL
            && arg->data.function.head->data.symbol == SYM_Times
            && arg->data.function.arg_count >= 1
            && arg->data.function.args[0]->type == EXPR_INTEGER
            && arg->data.function.args[0]->data.integer == -1) {
            negative = true;
            size_t n = arg->data.function.arg_count;
            if (n == 2) stripped = expr_copy(arg->data.function.args[1]);
            else {
                Expr** rest = (Expr**)malloc(sizeof(Expr*) * (n - 1));
                for (size_t i = 1; i < n; i++) rest[i - 1] = expr_copy(arg->data.function.args[i]);
                stripped = eval_and_free(internal_times(rest, n - 1));
                free(rest);
            }
        }
        if (negative && stripped) {
            Expr* inner = expr_new_function(
                expr_copy(e->data.function.head),
                (Expr*[]){stripped}, 1);
            Expr* neg = internal_times(
                (Expr*[]){expr_new_integer(-1), inner}, 2);
            return eval_and_free(neg);
        }
    }

    /* Recurse: rebuild the function with normalised children. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = intsimp_normalize_inverse_trig_signs(e->data.function.args[i]);
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* result = expr_new_function(head, new_args, n);
    free(new_args);
    /* Run a single evaluator pass so Plus / Times canonicalise. */
    return eval_and_free(result);
}

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
#include "arithmetic.h"
#include "core.h"
#include "sym_names.h"
#include "attr.h"
#include "symtab.h"
#include "numeric.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

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
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol.name == SYM_Rational) {
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
        expr_new_symbol(SYM_Simplify), (Expr*[]){ e }, 1);
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

/* True when `e` is an opaque parameter: an application of an UNDEFINED head
 * (not protected, no builtin, no rules) that is not itself a numeric quantity --
 * DSolve's generated constants C[2], a user's c[2] or f[a].  It denotes an
 * arbitrary value exactly like a bare symbol, so the positive-parameter sign
 * convention applies to it too.  Without this, the discriminant 16 C[2] of the
 * log part of Integrate[1/(C[2] - u^2), u] read as sign-unknown and the whole
 * ArcTanh contribution was lost (the integral came back 0). */
static bool intsimp_is_opaque_parameter(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (get_attributes(h) & ATTR_PROTECTED) return false;
    SymbolDef* def = symtab_lookup(h);
    if (def && (def->builtin_func || def->own_values || def->down_values)) return false;
    return !expr_is_numeric_quantity(e);
}

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
        && e->data.function.head->data.symbol.name == SYM_Rational
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
    if (intsimp_is_opaque_parameter(e)) {
        /* An opaque parameter such as DSolve's C[2] or a user's c[2]
         * is a free parameter exactly like a bare symbol. */
        return 1;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
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
        && e->data.function.head->data.symbol.name == SYM_Times) {
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
        && e->data.function.head->data.symbol.name == SYM_Plus) {
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
        && e->data.function.head->data.symbol.name == SYM_Sqrt
        && e->data.function.arg_count == 1) {
        int s = intsimp_sign_pos_assumption(e->data.function.args[0]);
        if (s > 0) return 1;
        return 0;
    }
    return 0;
}

int intsimp_numeric_sign(Expr* e) {
    if (!e) return 0;
    Expr* call = expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){expr_copy(e)}, 1);
    /* N[e] can land on a pole (1/0); this is a sign/branch probe whose value is
     * discarded but for its sign, so mute the Power::infy -- it is noise, and it
     * must not be seen by an enclosing Check[].  Covers the six intrat.c sign
     * fallbacks that funnel through here. */
    arith_warnings_mute_push();
    Expr* r = evaluate(call);
    arith_warnings_mute_pop();
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
 *   a + b Sqrt[c]       -> Sqrt[(a+d)/2] + sgn(b) Sqrt[(a-d)/2] when
 *                          d = Sqrt[a^2 - b^2 c] is rational (denesting)
 *   anything else       -> Sqrt[anything] (irreducible)
 */

/* Exact rational (Integer, BigInt or Rational) — no inexact numbers. */
static bool intsimp_is_exact_rational(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational;
}

/* Match `b * c^(±1/2)` with b an exact rational and c a positive integer,
 * returning b' and c such that the term equals b' * Sqrt[c]
 * (c^(-1/2) = Sqrt[c]/c).  *b_out is owned by the caller. */
static bool intsimp_match_rational_surd(const Expr* t, Expr** b_out,
                                        const Expr** c_out) {
    const Expr* coef = NULL;
    const Expr* pw = t;
    if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL
        && t->data.function.head->data.symbol.name == SYM_Times
        && t->data.function.arg_count == 2
        && intsimp_is_exact_rational(t->data.function.args[0])) {
        coef = t->data.function.args[0];
        pw = t->data.function.args[1];
    }
    if (!(pw->type == EXPR_FUNCTION && pw->data.function.head->type == EXPR_SYMBOL
          && pw->data.function.head->data.symbol.name == SYM_Power
          && pw->data.function.arg_count == 2))
        return false;
    const Expr* c = pw->data.function.args[0];
    const Expr* ex = pw->data.function.args[1];
    if (!(c->type == EXPR_INTEGER && c->data.integer > 1)) return false;
    if (!(ex->type == EXPR_FUNCTION && ex->data.function.head->type == EXPR_SYMBOL
          && ex->data.function.head->data.symbol.name == SYM_Rational
          && ex->data.function.arg_count == 2
          && ex->data.function.args[0]->type == EXPR_INTEGER
          && ex->data.function.args[1]->type == EXPR_INTEGER
          && ex->data.function.args[1]->data.integer == 2))
        return false;
    int64_t num = ex->data.function.args[0]->data.integer;
    if (num != 1 && num != -1) return false;
    Expr* b = coef ? expr_copy((Expr*)coef) : expr_new_integer(1);
    if (num == -1)
        b = eval_and_free(internal_divide((Expr*[]){b, expr_copy((Expr*)c)}, 2));
    *b_out = b;
    *c_out = c;
    return true;
}

/* Denest Sqrt[a + b Sqrt[c]] for exact rationals a > 0, b and integer c > 1:
 * when a^2 - b^2 c = d^2 with d rational (necessarily 0 <= d < a),
 *     Sqrt[a + b Sqrt[c]] = Sqrt[(a+d)/2] + sgn(b) Sqrt[(a-d)/2].
 * Both radicands are non-negative rationals, so the right side is the
 * positive real root.  Returns NULL when `e` is not of that shape or does
 * not denest.  The residue roots of a biquadratic Rothstein-Trager resultant
 * (e.g. 256 t^4 + 96 t^2 + 1 for (x^2+x+1)/(x^4+1), issue #133) land here:
 * Sqrt[3/4 - 1/Sqrt[2]] = 1/Sqrt[2] - 1/2. */
static Expr* intsimp_denest_quadratic_surd(const Expr* e) {
    if (!(e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
          && e->data.function.head->data.symbol.name == SYM_Plus
          && e->data.function.arg_count == 2))
        return NULL;
    const Expr* a = NULL;
    const Expr* sterm = NULL;
    for (int i = 0; i < 2; i++) {
        if (intsimp_is_exact_rational(e->data.function.args[i])) {
            a = e->data.function.args[i];
            sterm = e->data.function.args[1 - i];
        }
    }
    if (!a || intsimp_sign_pos_assumption((Expr*)a) <= 0) return NULL;
    Expr* b = NULL;
    const Expr* c = NULL;
    if (!intsimp_match_rational_surd(sterm, &b, &c)) return NULL;

    /* d^2 = a^2 - b^2 c (an exact rational). */
    Expr* d2 = eval_and_free(internal_plus((Expr*[]){
        internal_power((Expr*[]){expr_copy((Expr*)a), expr_new_integer(2)}, 2),
        internal_times((Expr*[]){expr_new_integer(-1),
            internal_power((Expr*[]){expr_copy(b), expr_new_integer(2)}, 2),
            expr_copy((Expr*)c)}, 3)}, 2));
    Expr* d = NULL;
    if (intsimp_is_exact_rational(d2) && intsimp_sign_pos_assumption(d2) >= 0) {
        Expr* sd = intsimp_pos_sqrt_factor(d2);
        if (intsimp_is_exact_rational(sd)) d = sd; else expr_free(sd);
    }
    expr_free(d2);
    if (!d) { expr_free(b); return NULL; }

    int bsign = intsimp_sign_pos_assumption(b);
    expr_free(b);
    Expr* half_p = eval_and_free(internal_divide((Expr*[]){
        internal_plus((Expr*[]){expr_copy((Expr*)a), expr_copy(d)}, 2),
        expr_new_integer(2)}, 2));
    Expr* half_m = eval_and_free(internal_divide((Expr*[]){
        internal_subtract((Expr*[]){expr_copy((Expr*)a), d}, 2),
        expr_new_integer(2)}, 2));
    Expr* r1 = intsimp_pos_sqrt_factor(half_p);
    Expr* r2 = intsimp_pos_sqrt_factor(half_m);
    expr_free(half_p); expr_free(half_m);
    if (bsign < 0)
        r2 = eval_and_free(internal_times((Expr*[]){expr_new_integer(-1), r2}, 2));
    return eval_and_free(internal_plus((Expr*[]){r1, r2}, 2));
}

Expr* intsimp_pos_sqrt_factor(Expr* e) {
    if (e->type == EXPR_INTEGER) {
        int64_t v = e->data.integer;
        if (v < 0) {
            /* Sqrt of a negative integer: leave inside Sqrt. */
            return expr_new_function(expr_new_symbol(SYM_Sqrt),
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
            return expr_new_function(expr_new_symbol(SYM_Sqrt),
                (Expr*[]){expr_new_integer(r)}, 1);
        }
        if (r == 1) return expr_new_integer(s);
        Expr* sqrt_part = expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_new_integer(r)}, 1);
        return eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(s), sqrt_part}, 2));
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        Expr* num = e->data.function.args[0];
        Expr* den = e->data.function.args[1];
        Expr* sn = intsimp_pos_sqrt_factor(num);
        Expr* sd = intsimp_pos_sqrt_factor(den);
        return eval_and_free(internal_divide((Expr*[]){sn, sd}, 2));
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp_e = e->data.function.args[1];
        if (exp_e->type == EXPR_INTEGER) {
            int64_t k = exp_e->data.integer;
            int b_sign = intsimp_sign_pos_assumption(base);
            if (k % 2 == 0 && b_sign > 0) {
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_copy(base),
                              expr_new_integer(k / 2)}, 2);
            }
            if (b_sign > 0) {
                /* Odd k, positive base: base^((k-1)/2) * Sqrt[base]. */
                int64_t kh = (k - 1) / 2;  /* may be negative for k < 0 */
                Expr* p1;
                if (kh == 0) p1 = expr_new_integer(1);
                else p1 = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_copy(base), expr_new_integer(kh)}, 2);
                Expr* p2 = expr_new_function(expr_new_symbol(SYM_Sqrt),
                    (Expr*[]){expr_copy(base)}, 1);
                return eval_and_free(internal_times(
                    (Expr*[]){p1, p2}, 2));
            }
        }
        return expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(e)}, 1);
    }
    if (e->type == EXPR_SYMBOL) {
        /* Sqrt[symbol] is irreducible — keep wrapped. */
        return expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(e)}, 1);
    }
    /* Quadratic surd a + b Sqrt[c]: denest when possible. */
    Expr* dn = intsimp_denest_quadratic_surd(e);
    if (dn) return dn;
    /* Fallback: no structural simplification. */
    return expr_new_function(expr_new_symbol(SYM_Sqrt),
        (Expr*[]){expr_copy(e)}, 1);
}

Expr* intsimp_pos_sqrt(Expr* e) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times) {
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
        && term->data.function.head->data.symbol.name == SYM_Times) {
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
        && term->data.function.head->data.symbol.name == SYM_Log
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
        && term->data.function.head->data.symbol.name == SYM_Log
        && term->data.function.arg_count == 1) {
        *coeff_out = expr_new_integer(1);
        *log_arg_out = expr_copy(term->data.function.args[0]);
        return true;
    }
    /* Times-tree (possibly unflattened) containing one Log factor; the
     * rest become the coefficient. */
    if (term->type != EXPR_FUNCTION
        || term->data.function.head->type != EXPR_SYMBOL
        || term->data.function.head->data.symbol.name != SYM_Times) return false;

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
        && e->data.function.head->data.symbol.name == SYM_Plus)) {
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
                Expr* prod_raw = internal_times(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2);
                Expr* prod = expr_expand(prod_raw); expr_free(prod_raw);
                Expr* logp = expr_new_function(expr_new_symbol(SYM_Log),
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
                    Expr* atanh = expr_new_function(expr_new_symbol(SYM_ArcTanh),
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
                    Expr* logp = expr_new_function(expr_new_symbol(SYM_Log),
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
        || rebuilt->data.function.head->data.symbol.name != SYM_Log
        || rebuilt->data.function.arg_count != 1) {
        return rebuilt;
    }

    Expr* arg = rebuilt->data.function.args[0];
    if (arg->type != EXPR_FUNCTION
        || arg->data.function.head->type != EXPR_SYMBOL
        || arg->data.function.head->data.symbol.name != SYM_Times) {
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

    Expr* new_log = expr_new_function(expr_new_symbol(SYM_Log),
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
/* Sign of the leading coefficient of the polynomial p under the
 * lexicographic order, variables in canonical order: +1, -1, or 0 when p is
 * not a polynomial with a rational leading coefficient. */
static int it_lead_sign_rec(Expr* p, Expr** vars, size_t nv) {
    if (nv == 0) {
        int64_t n, d;
        if (p->type == EXPR_BIGINT) return mpz_sgn(p->data.bigint);
        if (is_rational(p, &n, &d)) return ((n > 0) == (d > 0)) ? (n != 0) : -1;
        return 0;
    }
    int deg = get_degree_poly(p, vars[0]);
    Expr* c = get_coeff(p, vars[0], deg);
    int sg = it_lead_sign_rec(c, vars + 1, nv - 1);
    expr_free(c);
    return sg;
}

static int it_poly_lead_sign(Expr* arg) {
    size_t vc = 0, vcap = 8;
    Expr** vars = malloc(sizeof(Expr*) * vcap);
    collect_variables(arg, &vars, &vc, &vcap);
    int sg = 0;
    if (vc > 0 && is_polynomial(arg, vars, vc)) {
        bool plain = true;
        for (size_t i = 0; i < vc; i++) if (vars[i]->type != EXPR_SYMBOL) plain = false;
        if (plain) {
            qsort(vars, vc, sizeof(Expr*), compare_expr_ptrs);
            Expr* ex = expr_expand(arg);
            sg = it_lead_sign_rec(ex, vars, vc);
            expr_free(ex);
        }
    }
    for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
    free(vars);
    return sg;
}

Expr* intsimp_normalize_inverse_trig_signs(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);

    /* ArcTan[-(...)] -> -ArcTan[...], same for ArcTanh.  We look for
     * an argument that starts with a unary minus or a Times whose
     * first factor is a literal -1. */
    if (e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_ArcTan
            || e->data.function.head->data.symbol.name == SYM_ArcTanh)
        && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        bool negative = false;
        Expr* stripped = NULL;
        if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
            negative = true;
            stripped = expr_new_integer(-arg->data.integer);
        } else if (arg->type == EXPR_FUNCTION
            && arg->data.function.head->type == EXPR_SYMBOL
            && arg->data.function.head->data.symbol.name == SYM_Times
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
        /* A polynomial argument whose lexicographic leading coefficient is
         * negative (variables in canonical order): ArcTan[1 - 2 x^2] is
         * -ArcTan[-1 + 2 x^2].  Odd functions, so this is exact; it keeps
         * the printed form independent of which factoriser produced the
         * argument's orientation. */
        if (!negative && arg->type == EXPR_FUNCTION
            && arg->data.function.head->type == EXPR_SYMBOL
            && arg->data.function.head->data.symbol.name == SYM_Plus
            && it_poly_lead_sign(arg) < 0) {
            negative = true;
            stripped = eval_and_free(internal_times(
                (Expr*[]){ expr_new_integer(-1), expr_copy(arg) }, 2));
            Expr* ex = expr_expand(stripped);
            expr_free(stripped);
            stripped = ex;
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

/* ====================================================================
 * Final antiderivative normalisation — intsimp_finalize.
 * ====================================================================
 *
 * Applied once to every closed antiderivative (all Integrate methods).
 * Two independent clean-ups that the per-method pipelines leave on the
 * table:
 *
 *  (1) x-free reciprocal / product powers.  Substitution-based radical
 *      integrators emit constants like (1/a)^(1/3) that the core Power
 *      evaluator deliberately will NOT flatten (branch-cut safety: `a`
 *      is not provably positive).  Integration output carries implicit
 *      domain assumptions, so a PowerExpand restricted to *x-free*
 *      subexpressions is the correct, Mathematica-matching cleanup; it
 *      turns (1/a)^(2/3) -> a^(-2/3), which then collects against the
 *      surrounding a^k into the familiar a^(7/3) denominators.
 *
 *  (2) algebraic recombination.  The linear/quadratic-radical methods
 *      back-substitute u -> (a x + b)^(1/n) into an *un-combined* sum of
 *      Hermite terms whose denominators are polynomials in the compound
 *      radicand (a x + b).  A single Cancel[Together[.]] over the
 *      algebraic part collapses e.g.
 *        2 a^3 - 4 a^2 (a+b x) + 2 a (a+b x)^2  ->  2 a b^2 x^2
 *      because the integer-power radicand expands in x and the radical
 *      (a x + b)^(2/3) rides along as an opaque generator.  Gated to fire
 *      only when it strictly shrinks the algebraic part, so it never
 *      enlarges or reshapes a result that was already tidy.
 *
 *  (3) the inverse-trig arguments are then distributed (Expand) and
 *      sign-normalised so ArcTan[-(p+q)] -> -ArcTan[p+q].
 *
 * Memory: returns a fresh owned Expr*; does NOT free `r`.
 */

/* Cancel[Together[e]], consuming `e`. */
static Expr* it_cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_and_free(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_and_free(internal_cancel((Expr*[]){ t }, 1));
}

/* Positive under the integration convention (free symbols denote positive
 * reals).  PowerExpand is only branch-safe on such bases: rewriting
 * (1/a)^(1/3) -> a^(-1/3) is sound when a > 0, but the same rewrite on a
 * negative or complex numeric constant (e.g. (-3)^(1/2)) crosses a branch
 * cut and silently changes the value — which would corrupt the constant
 * coefficient of an antiderivative and break differentiate-back. */
static bool it_pos_base(const Expr* c) {
    if (!c) return false;
    if (c->type == EXPR_SYMBOL) return true;       /* assume positive */
    if (c->type == EXPR_INTEGER) return c->data.integer > 0;
    if (c->type == EXPR_REAL)    return c->data.real > 0.0;
    {
        int64_t n, d;
        if (is_rational(c, &n, &d)) return ((n > 0) == (d > 0)) && n != 0;
    }
    if (c->type == EXPR_FUNCTION
        && c->data.function.head->type == EXPR_SYMBOL) {
        const char* h = c->data.function.head->data.symbol.name;
        if (h == SYM_Times) {
            for (size_t i = 0; i < c->data.function.arg_count; i++)
                if (!it_pos_base(c->data.function.args[i])) return false;
            return true;
        }
        if (h == SYM_Power && c->data.function.arg_count == 2)
            return it_pos_base(c->data.function.args[0]);  /* b>0 => b^e>0 */
    }
    return false;
}

/* Rewrite every x-free Power[c, p/q] (q != 1, c positive under the
 * symbol-positive convention) node via PowerExpand so (1/c)^r -> c^(-r)
 * and (c1 c2)^r -> c1^r c2^r.  Only x-free, positive-base nodes are
 * touched; x-dependent radicals such as (a x + b)^(2/3) and negative /
 * complex numeric constants are left intact.  Returns a fresh copy;
 * `e` is untouched. */
static Expr* it_normalize_xfree_powers(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        int64_t pn, pd;
        if (is_rational(e->data.function.args[1], &pn, &pd) && pd != 1
            && it_pos_base(e->data.function.args[0])
            && intrat_freeq_test(e, x)) {
            Expr* pe = expr_new_function(expr_new_symbol(SYM_PowerExpand),
                                         (Expr*[]){ expr_copy(e) }, 1);
            return eval_and_free(pe);
        }
    }

    size_t n = e->data.function.arg_count;
    Expr** na = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        na[i] = it_normalize_xfree_powers(e->data.function.args[i], x);
    Expr* head = it_normalize_xfree_powers(e->data.function.head, x);
    Expr* r = expr_new_function(head, na, n);
    free(na);
    return r;
}

/* True if `e` contains a Power[base, p/q] (q != 1) whose base depends on
 * x — i.e. a genuine algebraic radical of the integration variable. */
static bool it_has_x_radical(const Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        int64_t pn, pd;
        if (is_rational(e->data.function.args[1], &pn, &pd) && pd != 1
            && !intrat_freeq_test(e->data.function.args[0], x))
            return true;
    }
    if (it_has_x_radical(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (it_has_x_radical(e->data.function.args[i], x)) return true;
    return false;
}

static bool it_is_invtrig_head(const char* s) {
    return s == SYM_ArcTan  || s == SYM_ArcTanh
        || s == SYM_ArcSin  || s == SYM_ArcCos  || s == SYM_ArcCot
        || s == SYM_ArcSinh || s == SYM_ArcCosh || s == SYM_ArcCoth;
}

/* True if `e` contains a Log or inverse-trig head anywhere. */
static bool it_term_is_transcendental(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* s = e->data.function.head->data.symbol.name;
        if (s == SYM_Log || it_is_invtrig_head(s)) return true;
    }
    if (it_term_is_transcendental(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (it_term_is_transcendental(e->data.function.args[i])) return true;
    return false;
}

/* Leading numeric sign of a single Times/atom term is negative. */
static bool it_summand_negative(const Expr* t) {
    if (!t) return false;
    if (t->type == EXPR_INTEGER) return t->data.integer < 0;
    if (t->type == EXPR_REAL)    return t->data.real < 0.0;
    if (t->type == EXPR_FUNCTION
        && t->data.function.head->type == EXPR_SYMBOL) {
        const char* s = t->data.function.head->data.symbol.name;
        if (s == SYM_Rational && t->data.function.arg_count == 2
            && t->data.function.args[0]->type == EXPR_INTEGER)
            return t->data.function.args[0]->data.integer < 0;
        if (s == SYM_Times && t->data.function.arg_count >= 1)
            return it_summand_negative(t->data.function.args[0]);
    }
    return false;
}

/* True if `e` has an overall negative leading sign: a Plus all of whose
 * summands are negative, or a single negative term. */
static bool it_all_negative(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Plus) {
        if (e->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!it_summand_negative(e->data.function.args[i])) return false;
        return true;
    }
    return it_summand_negative(e);
}

/* Distribute (Expand) and sign-normalise inverse-trig arguments:
 * ArcTan[(-p - q)/d] -> -ArcTan[(p + q)/d] (distributed). */
static Expr* it_tidy_invtrig_args(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;
    Expr** na = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        na[i] = it_tidy_invtrig_args(e->data.function.args[i]);
    Expr* head = expr_copy(e->data.function.head);
    Expr* rebuilt = expr_new_function(head, na, n);
    free(na);

    if (rebuilt->data.function.head->type == EXPR_SYMBOL
        && it_is_invtrig_head(rebuilt->data.function.head->data.symbol.name)
        && rebuilt->data.function.arg_count == 1) {
        Expr* ax = expr_expand(rebuilt->data.function.args[0]);
        Expr* ae = eval_and_free(ax);
        Expr* node;
        if (it_all_negative(ae)) {
            /* Distribute the negation so the argument becomes a genuinely
             * positive sum.  A bare Times[-1, Plus[...]] would be re-pulled
             * by ArcTan's own oddness, silently undoing the flip. */
            Expr* neg = internal_times(
                (Expr*[]){ expr_new_integer(-1), ae }, 2);    /* consumes ae */
            Expr* pos = expr_expand(neg);
            expr_free(neg);
            pos = eval_and_free(pos);
            Expr* inner = expr_new_function(
                expr_copy(rebuilt->data.function.head),
                (Expr*[]){ pos }, 1);
            node = eval_and_free(internal_times(
                (Expr*[]){ expr_new_integer(-1), inner }, 2));
        } else {
            node = expr_new_function(
                expr_copy(rebuilt->data.function.head),
                (Expr*[]){ ae }, 1);
        }
        expr_free(rebuilt);
        return node;
    }
    return rebuilt;
}

/* Recombine an algebraic (radical) sum into a single reduced fraction.
 * Cancel[Together] alone leaves the denominator as a polynomial in the
 * compound radicand (e.g. 6 a^4 - 12 a^3 (a+b x) + 6 a^2 (a+b x)^2);
 * splitting into Numerator / Denominator and Expand-ing each collapses
 * those integer-power compounds in x (-> 6 a^2 b^2 x^2) while the
 * fractional-power radicals ride through untouched.  A final Cancel
 * clears the shared polynomial factor.  Consumes nothing; returns a
 * fresh owned Expr* (or NULL on failure). */
static Expr* it_recombine_algebraic(Expr* alg) {
    Expr* cand = it_cancel_together(expr_copy(alg));
    if (!cand) return NULL;

    Expr* num = eval_and_free(internal_numerator(
        (Expr*[]){ expr_copy(cand) }, 1));
    Expr* den = eval_and_free(internal_denominator(
        (Expr*[]){ expr_copy(cand) }, 1));
    if (!num || !den) {
        if (num) expr_free(num);
        if (den) expr_free(den);
        return cand;
    }

    Expr* num_e = expr_expand(num);  expr_free(num);
    Expr* den_e = expr_expand(den);  expr_free(den);
    Expr* frac  = eval_and_free(internal_divide((Expr*[]){ num_e, den_e }, 2));
    Expr* out   = eval_and_free(internal_cancel((Expr*[]){ frac }, 1));
    expr_free(cand);
    return out;
}

/* ====================================================================
 * Output cleanup — real-branch repair of ArcTanh.
 * ====================================================================
 *
 * A substitution integrator can back-substitute into an antiderivative
 * carrying ArcTanh[g(x)] whose argument g is real with |g| > 1 on the
 * integrand's real domain.  ArcTanh's real interval is (-1, 1); for
 * |g| > 1 the value sits on the branch cut (1, Inf) U (-Inf, -1) and
 * gains a constant +-I Pi/2, so the antiderivative of a *real* integrand
 * comes out complex-valued (off by a constant imaginary part).  The
 * companion function ArcCoth has the SAME derivative — d/dz ArcTanh[z] =
 * d/dz ArcCoth[z] = 1/(1 - z^2) — and is real precisely where |z| > 1.
 * Rewriting ArcTanh[g] -> ArcCoth[g] there therefore selects the real
 * branch without changing the derivative (differentiate-back is
 * unaffected), matching the real closed form Mathematica returns.  E.g.
 * Integrate[Sqrt[x + Sqrt[x]], x] would otherwise carry
 * ArcTanh[Sqrt[1 + Sqrt[x]]/x^(1/4)] with the argument > 1 for all x > 0.
 */

/* Evaluate N[e]; on a real machine number store it in *out and return true.
 * Returns false for complex, infinite (a probe that hit a pole ->
 * ComplexInfinity), or still-symbolic results.  Borrows `e`. */
static bool it_num_real(const Expr* e, double* out) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_N),
                                   (Expr*[]){ expr_copy((Expr*)e) }, 1);
    /* N[e] of a probe point may hit a pole -> ComplexInfinity (the docstring
     * above expects this and returns false); mute the Power::infy it emits so a
     * branch probe does not leak or trip an enclosing Check[]. */
    arith_warnings_mute_push();
    Expr* r = eval_and_free(call);
    arith_warnings_mute_pop();
    bool ok = false;
    if (r) {
        if (r->type == EXPR_REAL)         { *out = r->data.real;            ok = true; }
        else if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
        else {
            int64_t nn, dd;
            if (is_rational(r, &nn, &dd) && dd != 0) {
                *out = (double)nn / (double)dd; ok = true;
            }
        }
        expr_free(r);
    }
    return ok;
}

/* True iff g(x) is real with |g| > 1 across a spread of positive-real sample
 * points — the branch condition under which ArcCoth[g] is the real form and
 * ArcTanh[g] sits on its cut.  Probes that hit a singularity / complex value
 * are skipped; a definite real |g| <= 1 or a nonzero imaginary part at any
 * probe vetoes the rewrite (conservative: the swap is derivative-preserving,
 * so a false negative merely leaves the pre-existing form). */
/* The free symbols of g other than x and the numeric constants, collected
 * once (borrowed pointers into g). */
static void it_free_params(const Expr* e, const Expr* x, const Expr*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (expr_eq(e, x) || expr_is_numeric_quantity(e)) return;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*out)[i], e)) return;
        if (*n == *cap) { *cap = *cap ? 2 * *cap : 4; *out = realloc(*out, *cap * sizeof(Expr*)); }
        (*out)[(*n)++] = e;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        it_free_params(e->data.function.args[i], x, out, n, cap);
}

/* True iff g(x) is real with |g| > 1 wherever it is real across a spread of
 * sample points -- the branch condition under which ArcCoth[g] = ArcTanh[1/g]
 * is the real form and ArcTanh[g] sits on its cut.  The probes run over
 * (0, 1) as well as beyond, so an antiderivative real only on a bounded
 * interval (Sqrt[1 - x^2] inside g) is still decided; a probe where g is not
 * real is outside that interval and is skipped.  Other free symbols (the
 * parameters: a in a/Sqrt[a^2 - y^2]) take fixed positive sample values.  A
 * definite real |g| <= 1 at any probe vetoes the rewrite (conservative: the
 * swap is derivative-preserving, so a false negative merely leaves the
 * pre-existing form). */
static bool it_arg_real_outside_unit(const Expr* g, Expr* x) {
    static const double probes[] = { 0.05, 0.15, 0.3, 0.45, 0.7, 0.9, 1.3, 2.0, 3.5, 8.0, 15.0 };
    static const double pvals[] = { 1.37, 2.11, 1.73, 2.9 };
    const size_t np = sizeof(probes) / sizeof(probes[0]);
    const Expr** ps = NULL;
    size_t npar = 0, pcap = 0;
    it_free_params(g, x, &ps, &npar, &pcap);
    if (npar > 4) { free(ps); return false; }
    Expr* gp = expr_copy((Expr*)g);
    for (size_t k = 0; k < npar; k++) {
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy((Expr*)ps[k]), expr_new_real(pvals[k]) }, 2);
        gp = eval_and_free(internal_replace_all((Expr*[]){ gp, rule }, 2));
    }
    free(ps);
    int valid = 0;
    bool veto = false;
    for (size_t i = 0; i < np && !veto && gp; i++) {
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy(x), expr_new_real(probes[i]) }, 2);
        arith_warnings_mute_push();          /* a probe may sit on a pole */
        Expr* sub = eval_and_free(internal_replace_all(
            (Expr*[]){ expr_copy(gp), rule }, 2));
        arith_warnings_mute_pop();
        if (!sub) continue;
        Expr* im_call = expr_new_function(expr_new_symbol(SYM_Im),
            (Expr*[]){ expr_copy(sub) }, 1);
        Expr* ab_call = expr_new_function(expr_new_symbol(SYM_Abs),
            (Expr*[]){ sub }, 1);                 /* consumes sub */
        double im, ab;
        bool im_ok = it_num_real(im_call, &im);
        bool ab_ok = it_num_real(ab_call, &ab);
        expr_free(im_call);
        expr_free(ab_call);
        if (!im_ok || !ab_ok) continue;           /* singular / complex-infinite probe */
        if (fabs(im) > 1e-7)   continue;          /* g not real here: outside the domain */
        if (ab <= 1.0 + 1e-7)  veto = true;       /* |g| <= 1: ArcTanh is the real branch */
        else valid++;
    }
    if (gp) expr_free(gp);
    return !veto && valid >= 2;
}

/* Walk `e`, rewriting every x-dependent ArcTanh[g] whose argument is real
 * with |g| > 1 on the real domain into the derivative-identical, real-valued
 * ArcCoth[g].  Returns a fresh owned tree; `e` is untouched. */
static Expr* it_realify_arctanh(const Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr** na = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        na[i] = it_realify_arctanh(e->data.function.args[i], x);
    Expr* head = it_realify_arctanh(e->data.function.head, x);
    Expr* rebuilt = expr_new_function(head, na, n);
    free(na);

    if (rebuilt->data.function.head->type == EXPR_SYMBOL
        && rebuilt->data.function.head->data.symbol.name == SYM_ArcTanh
        && rebuilt->data.function.arg_count == 1) {
        Expr* g = rebuilt->data.function.args[0];
        if (!intrat_freeq_test(g, x) && it_arg_real_outside_unit(g, x)) {
            /* ArcCoth[g] == ArcTanh[1/g]; the reciprocal form when it is
             * at most two leaves larger (a/Sqrt[a^2 - y^2] -> Sqrt[a^2 - y^2]/a). */
            Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_copy(g), expr_new_integer(-1) }, 2));
            Expr* out;
            if (inv && leaf_count_internal(inv, true) <= leaf_count_internal(g, true) + 2) {
                out = expr_new_function(expr_new_symbol(SYM_ArcTanh), (Expr*[]){ inv }, 1);
            } else {
                if (inv) expr_free(inv);
                out = expr_new_function(expr_new_symbol(SYM_ArcCoth),
                    (Expr*[]){ expr_copy(g) }, 1);
            }
            expr_free(rebuilt);
            return out;
        }
    }
    return rebuilt;
}

/* Largest algebraic part we will attempt to Cancel[Together] — a guard
 * against the multivariate-GCD blow-up on high-generator results. */
#define IT_RECOMBINE_GATE 800

/* True if `e` contains an unevaluated Integrate[...] head anywhere — a
 * partially-closed result we must not re-evaluate (it would re-enter the
 * integrator and can blow the recursion limit). */
static bool it_contains_integrate(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Integrate) return true;
    if (it_contains_integrate(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (it_contains_integrate(e->data.function.args[i])) return true;
    return false;
}

Expr* intsimp_finalize(Expr* r, Expr* x) {
    if (!r) return NULL;

    /* This pass is scoped to radical-bearing antiderivatives — the only
     * class where the per-method pipelines leave the deficiency (the
     * substitution integrators back-substitute an un-combined sum and emit
     * (1/a)^(k/q) constants).  Pure-rational / Risch / trig output is passed
     * straight through: re-touching it perturbs forms the corpus verifier
     * pins (and reciprocal-power flattening on a generic symbol would not be
     * branch-sound without an assumption the verifier does not share).  A
     * partially-closed result (nested Integrate) is likewise left as-is — re-
     * evaluating it would re-enter the integrator and risk the recursion
     * limit.  Both checks are single tree-walks; this function takes
     * ownership of `r` and returns it (or a derived tree). */
    if (it_contains_integrate(r) || !it_has_x_radical(r, x)) return r;

    /* (1) x-free power normalisation, then one evaluator pass so the
     *     freed reciprocal powers collect against neighbouring a^k. */
    Expr* p = it_normalize_xfree_powers(r, x);
    expr_free(r);
    p = eval_and_free(p);

    /* (2) split into algebraic vs transcendental summands. */
    bool is_plus = (p->type == EXPR_FUNCTION
        && p->data.function.head->type == EXPR_SYMBOL
        && p->data.function.head->data.symbol.name == SYM_Plus);
    size_t n = is_plus ? p->data.function.arg_count : 1;

    Expr** alg = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** trn = (Expr**)malloc(sizeof(Expr*) * n);
    size_t na2 = 0, nt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* term = is_plus ? p->data.function.args[i] : p;
        if (it_term_is_transcendental(term)) trn[nt++] = expr_copy(term);
        else                                 alg[na2++] = expr_copy(term);
    }

    Expr* alg_sum = NULL;
    if (na2 == 1) alg_sum = alg[0];
    else if (na2 > 1) alg_sum = eval_and_free(internal_plus(alg, na2));
    free(alg);

    Expr* trn_sum = NULL;
    if (nt == 1) trn_sum = trn[0];
    else if (nt > 1) trn_sum = eval_and_free(internal_plus(trn, nt));
    free(trn);

    /* (3) recombine the algebraic part when it strictly shrinks. */
    if (alg_sum && leaf_count_internal(alg_sum, true) <= IT_RECOMBINE_GATE) {
        Expr* cand = it_recombine_algebraic(alg_sum);
        if (cand && leaf_count_internal(cand, true)
                    < leaf_count_internal(alg_sum, true)) {
            expr_free(alg_sum);
            alg_sum = cand;
        } else if (cand) {
            expr_free(cand);
        }
    }

    /* (4) select the real branch of ArcTanh (ArcTanh[g] -> ArcCoth[g] where
     *     g is real with |g| > 1 on the real domain), then distribute +
     *     sign-normalise inverse-trig arguments. */
    if (trn_sum) {
        Expr* rb = it_realify_arctanh(trn_sum, x);
        expr_free(trn_sum);
        trn_sum = rb;
        Expr* t2 = it_tidy_invtrig_args(trn_sum);
        expr_free(trn_sum);
        trn_sum = t2;
    }

    expr_free(p);

    /* (5) reassemble + final evaluator pass. */
    Expr* out;
    if (alg_sum && trn_sum)
        out = eval_and_free(internal_plus((Expr*[]){ alg_sum, trn_sum }, 2));
    else if (alg_sum) out = alg_sum;
    else if (trn_sum) out = trn_sum;
    else out = expr_new_integer(0);
    return out;
}

/* ====================================================================
 * Algebraic-constant reduction — intsimp_reduce_algebraic_constants.
 * ====================================================================
 *
 * The Rothstein-Trager / LogToReal pipeline does its arithmetic on
 * algebraic constants (Sqrt[2], 2^(1/4), ...) without normalising them:
 * like the reference, the evaluator never rationalises a denominator, so
 * Integrate[1/(x^2 + 3 + 2 Sqrt[2]), x] came out as
 *   (4 ArcTan[(3 x + 2 Sqrt[2] x)/(7 + 5 Sqrt[2])])/(12 + 8 Sqrt[2])
 *     + (4 Sqrt[2] ArcTan[(3 x + 2 Sqrt[2] x)/(7 + 5 Sqrt[2])])/(12 + 8 Sqrt[2])
 * (leaf count 82) for (-1 + Sqrt[2]) ArcTan[(-1 + Sqrt[2]) x] (14).
 *
 * This pass rewrites a closed rational-function antiderivative so that
 *   - the argument of every Log / ArcTan / ArcTanh kernel is a polynomial
 *     in x (or a ratio of polynomials, denominator made monic), each
 *     coefficient an exact algebraic constant in RootReduce normal form;
 *   - terms sharing one kernel are merged, and the summed coefficient is
 *     put in RootReduce normal form.
 * Every step is an exact identity on numbers (RootReduce / ToRadicals).  A
 * coefficient is replaced only when its normal form is radical (no Root[]
 * object; see ar_reduce_const for the ToRadicals fallback past degree 2) and
 * no larger; the
 * rewritten antiderivative is kept only if its leaf count does not grow.
 * Results with no algebraic constant, or with inexact numbers, are
 * returned unchanged at the cost of one tree walk.
 */

static bool ar_is_exact_rational(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION && head_is(e, SYM_Rational);
}

/* Contains a Power[rational, non-integer rational] (Sqrt[2], 2^(1/4), ...). */
static bool ar_contains_algebraic(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2
        && ar_is_exact_rational(e->data.function.args[0])
        && e->data.function.args[1]->type == EXPR_FUNCTION
        && head_is(e->data.function.args[1], SYM_Rational))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ar_contains_algebraic(e->data.function.args[i])) return true;
    return false;
}

static bool ar_contains_inexact_or_root(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Root)) return true;
    if (ar_contains_inexact_or_root(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ar_contains_inexact_or_root(e->data.function.args[i])) return true;
    return false;
}

static int64_t ar_leaves(Expr* e) { return leaf_count_internal(e, false); }

static bool ar_polyq(Expr* p /* borrowed */, Expr* x) {
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
        (Expr*[]){ expr_copy(p), expr_copy(x) }, 2));
    bool ok = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
    if (r) expr_free(r);
    return ok;
}

static Expr* ar_call1(const char* head, Expr* a /* owned */) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1));
}

static Expr* ar_call2(const char* head, Expr* a /* owned */, Expr* b /* owned */) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2));
}

/* ToRadicals fallbacks allowed per antiderivative (see ar_reduce_const). */
#define AR_RADICAL_BUDGET 32

/* Denest every Sqrt[a + b Sqrt[c]] (and its reciprocal) in `e` that
 * intsimp_denest_quadratic_surd can open.  Borrows e, returns a new tree. */
static Expr* ar_denest_walk(Expr* e /* borrowed */) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* ex = e->data.function.args[1];
        bool half = ex->type == EXPR_FUNCTION && head_is(ex, SYM_Rational)
            && ex->data.function.args[0]->type == EXPR_INTEGER
            && ex->data.function.args[1]->type == EXPR_INTEGER
            && ex->data.function.args[1]->data.integer == 2
            && (ex->data.function.args[0]->data.integer == 1
                || ex->data.function.args[0]->data.integer == -1);
        if (half) {
            Expr* base = ar_denest_walk(e->data.function.args[0]);
            Expr* be = expr_expand(base);
            expr_free(base);
            Expr* d = intsimp_denest_quadratic_surd(be);
            if (d) {
                expr_free(be);
                if (ex->data.function.args[0]->data.integer == 1) return d;
                return eval_and_free(internal_power((Expr*[]){ d, expr_new_integer(-1) }, 2));
            }
            return eval_and_free(internal_power((Expr*[]){ be, expr_copy(ex) }, 2));
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = ar_denest_walk(e->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return eval_and_free(out);
}

/* Normal form of an x-free algebraic constant, when radical and no larger.
 * RootReduce (sub-millisecond) is canonical for a quadratic field:
 * (4 + 4 Sqrt[2])/(12 + 8 Sqrt[2]) -> -1 + Sqrt[2].  Past degree 2 it
 * answers a Root[] object, never acceptable in an antiderivative; then the
 * Root is converted back with ToRadicals and its quadratic surds denested
 * (about 1 ms, at most *budget times per antiderivative):
 *   (4 Sqrt[3] - 4 Sqrt[2])/(20 - 8 Sqrt[6]) -> Root[1 - 10 #1^2 + #1^4 &, 4]
 *     -> Sqrt[5 + 2 Sqrt[6]] -> Sqrt[2] + Sqrt[3]
 *   Sqrt[16 + 16 Sqrt[2]]/(4 + 4 Sqrt[2]) -> Sqrt[-1 + Sqrt[2]].
 * Simplify finds the same forms but costs ~10 ms per constant, which
 * doubled the integrator's time on biquadratic answers. */
static Expr* ar_reduce_const(Expr* c /* borrowed */, int* budget) {
    if (!ar_contains_algebraic(c)) return expr_copy(c);
    Expr* r = ar_call1("RootReduce", expr_copy(c));
    if (r && head_is(r, SYM_Root) && *budget > 0) {
        (*budget)--;
        Expr* rad = ar_call1("ToRadicals", r);
        r = rad ? ar_denest_walk(rad) : NULL;
        if (rad) expr_free(rad);
    }
    if (!r || ar_contains_inexact_or_root(r) || ar_leaves(r) > ar_leaves(c)) {
        if (r) expr_free(r);
        return expr_copy(c);
    }
    return r;
}

/* p a polynomial in x: rebuild Sum[reduce(c_k) x^k].  Returns NULL when p
 * is not polynomial in x. */
static Expr* ar_poly_coeffs(Expr* p /* borrowed */, Expr* x, int* budget) {
    Expr* pe = expr_expand(p);
    if (!ar_polyq(pe, x)) { expr_free(pe); return NULL; }
    Expr* cl = ar_call2("CoefficientList", pe, expr_copy(x));
    if (!cl || !head_is(cl, SYM_List)) { if (cl) expr_free(cl); return NULL; }
    size_t n = cl->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t k = 0; k < n; k++) {
        Expr* c = ar_reduce_const(cl->data.function.args[k], budget);
        terms[k] = internal_times((Expr*[]){ c,
            internal_power((Expr*[]){ expr_copy(x), expr_new_integer((int64_t)k) }, 2) }, 2);
    }
    expr_free(cl);
    Expr* out = eval_and_free(internal_plus(terms, n));
    free(terms);
    return out;
}

/* A kernel argument (rational function of x): normalise its coefficients.
 * Falls back to a copy of `a` when that does not shrink it. */
static Expr* ar_tidy_arg(Expr* a /* borrowed */, Expr* x, int* budget) {
    if (!ar_contains_algebraic(a)) return expr_copy(a);
    if (intrat_freeq_test(a, x)) return ar_reduce_const(a, budget);
    Expr* tg = ar_call1("Together", expr_copy(a));
    Expr* num = tg ? ar_call1("Numerator", expr_copy(tg)) : NULL;
    Expr* den = tg ? ar_call1("Denominator", expr_copy(tg)) : NULL;
    if (tg) expr_free(tg);
    Expr* out = NULL;
    if (num && den && intrat_freeq_test(den, x)) {
        Expr* q = eval_and_free(internal_divide((Expr*[]){ expr_copy(num), expr_copy(den) }, 2));
        out = ar_poly_coeffs(q, x, budget);
        expr_free(q);
    } else if (num && den && ar_polyq(den, x)) {
        /* make the denominator monic in x, then normalise both sides */
        int dd = get_degree_poly(den, x);
        Expr* lc = dd > 0 ? get_coeff(den, x, dd) : NULL;
        if (lc && intrat_freeq_test(lc, x)) {
            Expr* nq = eval_and_free(internal_divide((Expr*[]){ expr_copy(num), expr_copy(lc) }, 2));
            Expr* dq = eval_and_free(internal_divide((Expr*[]){ expr_copy(den), expr_copy(lc) }, 2));
            Expr* n2 = ar_poly_coeffs(nq, x, budget);
            Expr* d2 = ar_poly_coeffs(dq, x, budget);
            expr_free(nq); expr_free(dq);
            if (n2 && d2) out = eval_and_free(internal_divide((Expr*[]){ n2, d2 }, 2));
            else { if (n2) expr_free(n2); if (d2) expr_free(d2); }
        }
        if (lc) expr_free(lc);
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    if (!out || ar_leaves(out) > ar_leaves(a)) {
        if (out) expr_free(out);
        return expr_copy(a);
    }
    return out;
}

static bool ar_is_kernel(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && (head_is(e, SYM_Log) || head_is(e, SYM_ArcTan) || head_is(e, SYM_ArcTanh));
}

/* Rebuild every kernel with tidied arguments (a kernel inside a kernel's
 * argument is not expected in a rational antiderivative and is left). */
static Expr* ar_tidy_kernels(Expr* e /* borrowed */, Expr* x, int* budget) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    bool kern = ar_is_kernel(e);
    for (size_t i = 0; i < n; i++)
        args[i] = kern ? ar_tidy_arg(e->data.function.args[i], x, budget)
                       : ar_tidy_kernels(e->data.function.args[i], x, budget);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return eval_and_free(out);
}

/* Split a summand into (x-free coefficient) * (rest); both owned. */
static void ar_split_term(Expr* t /* borrowed */, Expr* x, Expr** coef, Expr** rest) {
    if (intrat_freeq_test(t, x)) { *coef = expr_copy(t); *rest = expr_new_integer(1); return; }
    if (!(t->type == EXPR_FUNCTION && head_is(t, SYM_Times))) {
        *coef = expr_new_integer(1); *rest = expr_copy(t); return;
    }
    size_t n = t->data.function.arg_count;
    Expr** cf = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** rf = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nc = 0, nr = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = t->data.function.args[i];
        if (intrat_freeq_test(f, x)) cf[nc++] = expr_copy(f);
        else rf[nr++] = expr_copy(f);
    }
    *coef = nc ? eval_and_free(internal_times(cf, nc)) : expr_new_integer(1);
    *rest = nr ? eval_and_free(internal_times(rf, nr)) : expr_new_integer(1);
    free(cf); free(rf);
}

Expr* intsimp_reduce_algebraic_constants(Expr* e /* borrowed */, Expr* x) {
    if (!e || !ar_contains_algebraic(e) || ar_contains_inexact_or_root(e)
        || it_contains_integrate(e))
        return expr_copy(e);
    int budget = AR_RADICAL_BUDGET;
    Expr* k = ar_tidy_kernels(e, x, &budget);
    if (!k) return expr_copy(e);

    /* Merge summands that share their x-dependent part. */
    Expr* out = NULL;
    if (k->type == EXPR_FUNCTION && head_is(k, SYM_Plus)) {
        size_t n = k->data.function.arg_count;
        Expr** rests = (Expr**)malloc(sizeof(Expr*) * n);
        Expr** coefs = (Expr**)malloc(sizeof(Expr*) * n);
        size_t ng = 0;
        for (size_t i = 0; i < n; i++) {
            Expr *c = NULL, *r = NULL;
            ar_split_term(k->data.function.args[i], x, &c, &r);
            size_t g = 0;
            while (g < ng && !expr_eq(rests[g], r)) g++;
            if (g < ng) {
                coefs[g] = eval_and_free(internal_plus((Expr*[]){ coefs[g], c }, 2));
                expr_free(r);
            } else {
                rests[ng] = r; coefs[ng] = c; ng++;
            }
        }
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * (ng ? ng : 1));
        for (size_t g = 0; g < ng; g++) {
            Expr* c = ar_reduce_const(coefs[g], &budget);
            expr_free(coefs[g]);
            terms[g] = internal_times((Expr*[]){ c, rests[g] }, 2);
        }
        out = eval_and_free(internal_plus(terms, ng));
        free(terms); free(rests); free(coefs);
        expr_free(k);
    } else {
        Expr *c = NULL, *r = NULL;
        ar_split_term(k, x, &c, &r);
        Expr* c2 = ar_reduce_const(c, &budget);
        expr_free(c);
        out = eval_and_free(internal_times((Expr*[]){ c2, r }, 2));
        expr_free(k);
    }
    if (!out || ar_leaves(out) > ar_leaves(e)) {
        if (out) expr_free(out);
        return expr_copy(e);
    }
    return out;
}

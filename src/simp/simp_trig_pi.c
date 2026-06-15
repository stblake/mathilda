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
/* Trig at rational multiples of Pi: shortest-numerator canonicalization   */
/* ----------------------------------------------------------------------- */

/* For inputs like Cos[4 Pi/9] - Sin[Pi/18] the difference is exactly 0
 * by the complement identity Cos[Pi/2 - x] == Sin[x] (here Pi/2 - 4Pi/9
 * = Pi/18).  builtin_cos / builtin_sin handle a whitelist of "nice"
 * denominators (1,2,3,4,5,6,10,12) and otherwise leave the call as is,
 * so the two terms never become structurally equal and additive
 * cancellation in simp_search never fires.
 *
 * This transform rewrites every Sin/Cos/Tan/Cot/Sec/Csc of a rational
 * multiple of Pi into a unique form: pick the (Sin vs Cos / Tan vs Cot
 * / Sec vs Csc) representation whose reduced fraction has the smaller
 * numerator.  After the rewrite, Cos[4 Pi/9] and Sin[Pi/18] both land
 * at Sin[Pi/18] and the surrounding Plus collapses to 0.  Cos[5 Pi/9]
 * lands at -Sin[Pi/18] (since 5/9 > 1/2 picks up a sign via Cos[Pi -
 * x] = -Cos[x] before the complement swap), which lets the Morrie's-
 * law product 1/8 (Cos[4Pi/9] + Cos[5Pi/9]) collapse to 0 too.
 *
 * Idempotent: re-applying to the canonical form returns it unchanged,
 * so it is safe to seed without bounding the round count.
 *
 * Why a Simplify-only transform (rather than wiring into builtin_cos
 * / builtin_sin): the user-facing default print for Cos[4 Pi/9] should
 * stay as-written, both to match Mathematica and to keep regression
 * tests stable.  Only Simplify needs the unified representation. */

static int64_t trig_pi_i64_gcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static void trig_pi_reduce_frac(int64_t* n, int64_t* d) {
    if (*d == 0) return;
    int64_t g = trig_pi_i64_gcd(*n, *d);
    if (g > 1) { *n /= g; *d /= g; }
    if (*d < 0) { *d = -*d; *n = -*n; }
}

/* Returns true if `arg` is structurally (n/d) * Pi or Pi (n=d=1) and
 * sets n, d.  Mirrors trig.c's extract_pi_multiplier; duplicated here
 * so simp.c does not have to expose that file-static helper. */
static bool trig_pi_extract(const Expr* arg, int64_t* n_out, int64_t* d_out) {
    if (arg->type == EXPR_SYMBOL && arg->data.symbol == SYM_Pi) {
        *n_out = 1; *d_out = 1; return true;
    }
    if (arg->type != EXPR_FUNCTION || !arg->data.function.head ||
        arg->data.function.head->type != EXPR_SYMBOL ||
        arg->data.function.head->data.symbol != SYM_Times ||
        arg->data.function.arg_count != 2) {
        return false;
    }
    Expr* a0 = arg->data.function.args[0];
    Expr* a1 = arg->data.function.args[1];
    if (!(a1->type == EXPR_SYMBOL && a1->data.symbol == SYM_Pi)) return false;
    return is_rational(a0, n_out, d_out);
}

/* Build (out_n / out_d) * Pi as an Expr*. */
static Expr* trig_pi_make_arg(int64_t n, int64_t d) {
    if (n == 0) return expr_new_integer(0);
    if (n == 1 && d == 1) return expr_new_symbol(SYM_Pi);
    Expr* coeff;
    if (d == 1) {
        coeff = expr_new_integer(n);
    } else {
        coeff = expr_new_function(expr_new_symbol(SYM_Rational),
            (Expr*[]){ expr_new_integer(n), expr_new_integer(d) }, 2);
    }
    Expr* args[2] = { coeff, expr_new_symbol(SYM_Pi) };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), args, 2));
}

/* Compute the canonical (head, n, d, sign) for `head[(n/d) Pi]`.  See
 * the file-header comment block above for the rule set.
 *
 * `head_in` is one of the SYM_* trig pointers (Sin/Cos/Tan/Cot/Sec/Csc).
 * On success returns true and fills the outputs; on failure (head not
 * recognised, integer overflow, etc.) returns false and outputs are
 * untouched. */
static bool trig_pi_canon_one(const char* head_in, int64_t n, int64_t d,
                              const char** head_out, int64_t* n_out,
                              int64_t* d_out, int* sign_out) {
    if (d <= 0) return false;
    /* Guard against multiplication overflow when computing 2*d, alt_d. */
    if (d > (INT64_MAX / 4)) return false;
    if (n > (INT64_MAX / 4) || n < -(INT64_MAX / 4)) return false;

    int sign = 1;

    bool is_tan_family = (head_in == SYM_Tan || head_in == SYM_Cot);
    bool is_sin_like = (head_in == SYM_Sin || head_in == SYM_Csc);
    bool is_cos_like = (head_in == SYM_Cos || head_in == SYM_Sec);
    if (!is_tan_family && !is_sin_like && !is_cos_like) return false;

    if (is_tan_family) {
        /* Period Pi.  Reduce to [0, d). */
        n = n % d;
        if (n < 0) n += d;
        /* Tan[Pi - x] = -Tan[x], Cot[Pi - x] = -Cot[x]. */
        if (2 * n > d) {
            sign = -sign;
            n = d - n;
        }
        /* n in [0, d/2]. */
    } else {
        int64_t two_d = 2 * d;
        n = n % two_d;
        if (n < 0) n += two_d;
        if (is_sin_like) {
            /* Sin[Pi + x] = -Sin[x], Csc[Pi + x] = -Csc[x]. */
            if (n > d) { sign = -sign; n -= d; }
            /* Now n in [0, d].  Sin[Pi - x] = Sin[x], Csc same. */
            if (2 * n > d) n = d - n;
        } else { /* cos-like */
            /* Cos[2 Pi - x] = Cos[x], Sec same. */
            if (n > d) n = two_d - n;
            /* Now n in [0, d].  Cos[Pi - x] = -Cos[x], Sec same. */
            if (2 * n > d) { sign = -sign; n = d - n; }
        }
        /* n in [0, d/2]. */
    }

    trig_pi_reduce_frac(&n, &d);

    /* n == 0: argument collapsed to 0 (Sin[0]=0, Tan[0]=0, Cos[0]=1, ...).
     * builtin_sin / builtin_cos / etc. already handle the EXPR_INTEGER 0
     * case, so just emit head[0] -- the surrounding evaluate() pass will
     * resolve it. */
    if (n == 0) {
        *head_out = head_in;
        *n_out = 0;
        *d_out = 1;
        *sign_out = sign;
        return true;
    }

    /* Complement: head[(n/d) Pi] == alt_head[(d - 2n)/(2d) Pi].  The
     * sign relation is the identity Sin = Cos ∘ (Pi/2 - .) etc., which
     * is sign-preserving under the post-quadrant reduction we did above. */
    int64_t alt_n = d - 2 * n;
    int64_t alt_d = 2 * d;
    trig_pi_reduce_frac(&alt_n, &alt_d);

    const char* alt_head;
    if      (head_in == SYM_Sin) alt_head = "Cos";
    else if (head_in == SYM_Cos) alt_head = "Sin";
    else if (head_in == SYM_Tan) alt_head = "Cot";
    else if (head_in == SYM_Cot) alt_head = "Tan";
    else if (head_in == SYM_Sec) alt_head = "Csc";
    else                          alt_head = "Sec"; /* SYM_Csc */

    /* Pick the rep with the smaller reduced numerator.  Tie-break: keep
     * the input head so we never thrash on already-canonical forms. */
    if (alt_n < n) {
        *head_out = alt_head;
        *n_out = alt_n;
        *d_out = alt_d;
    } else {
        *head_out = head_in;
        *n_out = n;
        *d_out = d;
    }
    *sign_out = sign;
    return true;
}

/* Walk `e`, applying trig_pi_canon_one at every Sin/Cos/Tan/Cot/Sec/Csc
 * call whose argument is a rational multiple of Pi.  Returns a freshly
 * owned, evaluated tree (always non-NULL).  Idempotent on the canonical
 * form. */
static Expr* simp_trig_pi_canon_walk(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Apply at this node when shape matches. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.arg_count == 1) {
        const char* h = e->data.function.head->data.symbol;
        bool is_trig = (h == SYM_Sin || h == SYM_Cos ||
                        h == SYM_Tan || h == SYM_Cot ||
                        h == SYM_Sec || h == SYM_Csc);
        if (is_trig) {
            int64_t n, d;
            if (trig_pi_extract(e->data.function.args[0], &n, &d)) {
                const char* out_h;
                int64_t out_n, out_d;
                int out_sign;
                if (trig_pi_canon_one(h, n, d, &out_h, &out_n, &out_d, &out_sign)) {
                    Expr* arg = trig_pi_make_arg(out_n, out_d);
                    Expr* call = expr_new_function(
                        expr_new_symbol(out_h),
                        (Expr*[]){ arg }, 1);
                    Expr* call_eval = eval_and_free(call);
                    if (out_sign == -1) {
                        Expr* neg = expr_new_function(
                            expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(-1), call_eval }, 2);
                        return eval_and_free(neg);
                    }
                    return call_eval;
                }
            }
        }
    }

    /* Recurse into children. */
    size_t n_args = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n_args);
    bool any_changed = false;
    for (size_t i = 0; i < n_args; i++) {
        Expr* nc = simp_trig_pi_canon_walk(e->data.function.args[i]);
        if (!expr_eq(nc, e->data.function.args[i])) any_changed = true;
        new_args[i] = nc;
    }
    /* Walk the head too in case of e.g. Hold[Sin][...]. */
    Expr* new_head = simp_trig_pi_canon_walk(e->data.function.head);
    if (!expr_eq(new_head, e->data.function.head)) any_changed = true;

    if (!any_changed) {
        for (size_t i = 0; i < n_args; i++) expr_free(new_args[i]);
        free(new_args);
        expr_free(new_head);
        return expr_copy((Expr*)e);
    }
    Expr* res = expr_new_function(new_head, new_args, n_args);
    free(new_args);
    return eval_and_free(res);
}

/* Idempotent (re-application is a structural fixed point).  Inert on
 * inputs without a Sin/Cos/Tan/Cot/Sec/Csc of a rational multiple of
 * Pi -- the walker descends but never triggers a rewrite. */
Expr* simp_trig_pi_canon(const Expr* e) {
    return simp_trig_pi_canon_walk(e);
}

static Expr* transform_pythag_reduce_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ 1 - Cos[x_]^2 + r___  :> Sin[x]^2 + r, "
            "  1 - Sin[x_]^2 + r___  :> Cos[x]^2 + r, "
            "  -1 + Cos[x_]^2 + r___ :> -Sin[x]^2 + r, "
            "  -1 + Sin[x_]^2 + r___ :> -Cos[x]^2 + r, "
            "  -1 + Cosh[x_]^2 + r___ :> Sinh[x]^2 + r, "
            "  1 + Sinh[x_]^2 + r___ :> Cosh[x]^2 + r, "
            "  1 - Cosh[x_]^2 + r___ :> -Sinh[x]^2 + r, "
            "  -1 - Sinh[x_]^2 + r___ :> -Cosh[x]^2 + r, "
            /* Reciprocal-pair identities. tanh^2 + sech^2 == 1, so
             *   1 - Tanh^2 -> Sech^2  and  -1 + Tanh^2 -> -Sech^2.
             * coth^2 - csch^2 == 1, so
             *   -1 + Coth^2 -> Csch^2 and  1 - Coth^2 -> -Csch^2.
             * tan^2 + 1 == sec^2, cot^2 + 1 == csc^2 (real-valued
             * Pythagorean trig).  These resolve a tied-score plateau
             * where the simp_search round loop's strict `<` tiebreak
             * would otherwise prefer the bare Plus form (e.g. score 7 =
             * score 7 for `-1 + Tanh^2` vs `-Sech^2`); fired here, the
             * structural collapse to a single Power head wins outright. */
            "  1 - Tanh[x_]^2 + r___  :> Sech[x]^2 + r, "
            "  -1 + Tanh[x_]^2 + r___ :> -Sech[x]^2 + r, "
            "  -1 + Coth[x_]^2 + r___ :> Csch[x]^2 + r, "
            "  1 - Coth[x_]^2 + r___  :> -Csch[x]^2 + r, "
            "  1 + Tan[x_]^2 + r___   :> Sec[x]^2 + r, "
            "  -1 - Tan[x_]^2 + r___  :> -Sec[x]^2 + r, "
            "  1 + Cot[x_]^2 + r___   :> Csc[x]^2 + r, "
            "  -1 - Cot[x_]^2 + r___  :> -Csc[x]^2 + r, "
            /* Reverse-direction reciprocal-pair identities.  Sec^2 - 1 ==
             * Tan^2, Csc^2 - 1 == Cot^2 (real-valued Pythagorean trig);
             * Sech^2 + (-1) doesn't hold (Sech^2 = 1 - Tanh^2, so 1 -
             * Sech^2 == Tanh^2); 1 + Csch^2 == Coth^2.  Without these
             * rules, a Plus shape like `-1 + Sec[x]^2 - Tan[x]^2` (the
             * post-Expand form of `(Sec+1)(Sec-1) - Tan^2`) would only
             * collapse via the TAN -> SEC direction `-1 - Tan^2 -> -Sec^2`,
             * which works for the simple shape but leaves the SEC term
             * "stranded" in coefficient-bearing forms PythagCanon must
             * separately rewrite.  Including both directions here lets
             * the cheap PythagReduce rule fire on either presentation. */
            "  -1 + Sec[x_]^2 + r___  :> Tan[x]^2 + r, "
            "  1 - Sec[x_]^2 + r___   :> -Tan[x]^2 + r, "
            "  -1 + Csc[x_]^2 + r___  :> Cot[x]^2 + r, "
            "  1 - Csc[x_]^2 + r___   :> -Cot[x]^2 + r, "
            "  1 - Sech[x_]^2 + r___  :> Tanh[x]^2 + r, "
            "  -1 + Sech[x_]^2 + r___ :> -Tanh[x]^2 + r, "
            "  1 + Csch[x_]^2 + r___  :> Coth[x]^2 + r, "
            "  -1 - Csch[x_]^2 + r___ :> -Coth[x]^2 + r }");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Fast-skip: every PythagReduce rule LHS has a Cos/Sin/Cosh/Sinh
     * pattern.  If the input contains none of those heads, the
     * ReplaceRepeated walk would visit every node and try every rule
     * and find nothing -- which on huge sum-of-exponentials inputs
     * costs 50-120 ms per call.  Skip the rewrite and return a copy
     * unchanged. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagReduce", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol(SYM_ReplaceRepeated), args, 2);
    Expr* out = eval_and_free(call);
    if (dbg) simp_debug_log("PythagReduce", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* PythagReduce sees the highest call volume of any simp transform
 * (~200 calls in the Tan double-angle case, with ~20 % unique inputs).
 * The memo dedupes the rest. */
Expr* transform_pythag_reduce(const Expr* e) {
    return simp_memo_wrap(e, "$PythagReduce", transform_pythag_reduce_impl);
}

/* PythagCanon: substitution-based Pythagorean canonicalizer.
 *
 * The bare PythagReduce rules look for `1 +/- Cos[x_]^2 + r___` style
 * shapes -- they fire only when the unit constant and the squared trig
 * appear additively, with unit coefficient. After a generic Expand of an
 * input like
 *     18 (Cos[x]+1)(Cos[x]-1)(Cos[y]^2-1)^2 (x-1) + 18 (x-1) Sin[x]^2 Sin[y]^4
 * we get a polynomial in Cos[x], Cos[y] whose individual Cos^2 / Cos^4
 * factors carry coefficients other than 1, so PythagReduce misses every
 * one of them -- and the cancellation against the Sin^2 Sin^4 term is
 * never recognised by the round loop.
 *
 * This transform is coefficient-blind. It substitutes every even-power
 *     Cos[x_]^(2k) -> (1 - Sin[x_]^2)^k
 * (and the reverse Sin -> Cos, plus the hyperbolic counterparts) via a
 * single ReplaceRepeated pass, then Expands. For each of the four
 * directions it scores the result against the input and keeps the best
 * strict win. Idempotent on inputs without an even Cos^k / Sin^k power;
 * inert (returns a structural copy of the input) when no direction beats
 * the input score.
 *
 * Why all four directions: the choice of "all-Sin" vs "all-Cos" depends
 * on which side already has more mass. A user input that is mostly
 * Sin^2 + small Cos^2 minus 1 wants Cos -> 1 - Sin; the reverse leaves
 * the small Cos^2 term and bloats the rest. Trying both keeps us
 * agnostic. The hyperbolic pair is exactly analogous via
 * Cosh^2 - Sinh^2 = 1. */
static Expr* transform_pythag_canon_impl(const Expr* e) {
    static Expr* rules_to_sin = NULL;
    static Expr* rules_to_cos = NULL;
    static Expr* rules_to_sinh = NULL;
    static Expr* rules_to_cosh = NULL;
    /* Reciprocal-pair canonicalisation directions.  Same shape as the
     * Sin/Cos rules above, but for the Pythagorean identities
     *     Sec^2 = 1 + Tan^2,   Csc^2 = 1 + Cot^2
     *     Sech^2 = 1 - Tanh^2, Csch^2 = -1 + Coth^2.
     * These are needed for inputs like `(Sec[x]+1)(Sec[x]-1) - Tan[x]^2`
     * whose post-Expand form `-1 + Sec[x]^2 - Tan[x]^2` collapses when
     * Sec[x]^2 -> 1 + Tan[x]^2 (or Tan[x]^2 -> -1 + Sec[x]^2) is
     * substituted globally.  Without these directions the round loop
     * would have to ride out a TrigFactor pass (Sec -> 1/Cos rewrite +
     * polynomial reformulation), which on the multi-variable Sec/Tan
     * test case takes 700ms+ for a single call.  All eight rules guard
     * with `n >= 2 && EvenQ[n]` to keep odd-power forms intact. */
    static Expr* rules_to_tan = NULL;
    static Expr* rules_to_sec = NULL;
    static Expr* rules_to_cot = NULL;
    static Expr* rules_to_csc = NULL;
    static Expr* rules_to_tanh = NULL;
    static Expr* rules_to_sech = NULL;
    static Expr* rules_to_coth = NULL;
    static Expr* rules_to_csch = NULL;
    if (!rules_to_sin) {
        rules_to_sin = parse_expression(
            "{ Cos[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Sin[x]^2)^(n/2) }");
        rules_to_cos = parse_expression(
            "{ Sin[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Cos[x]^2)^(n/2) }");
        rules_to_sinh = parse_expression(
            "{ Cosh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Sinh[x]^2)^(n/2) }");
        rules_to_cosh = parse_expression(
            "{ Sinh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Cosh[x]^2)^(n/2) }");
        rules_to_tan = parse_expression(
            "{ Sec[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Tan[x]^2)^(n/2) }");
        rules_to_sec = parse_expression(
            "{ Tan[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Sec[x]^2)^(n/2) }");
        rules_to_cot = parse_expression(
            "{ Csc[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Cot[x]^2)^(n/2) }");
        rules_to_csc = parse_expression(
            "{ Cot[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Csc[x]^2)^(n/2) }");
        rules_to_tanh = parse_expression(
            "{ Sech[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Tanh[x]^2)^(n/2) }");
        rules_to_sech = parse_expression(
            "{ Tanh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Sech[x]^2)^(n/2) }");
        rules_to_coth = parse_expression(
            "{ Csch[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Coth[x]^2)^(n/2) }");
        rules_to_csch = parse_expression(
            "{ Coth[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Csch[x]^2)^(n/2) }");
    }
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Same fast-skip as PythagReduce: every rule LHS targets a
     * Cos/Sin/Cosh/Sinh head. Without one nothing can fire. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagCanon", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* best = expr_copy((Expr*)e);
    size_t best_score = score_with_func(best, NULL);

    /* Expand first: the substitution rule matches Power[Cos[x], 2k]
     * literally, which only appears after distributing factored forms
     * like (Cos+1)(Cos-1) -> Cos^2 - 1. Without the pre-Expand the
     * rules cannot fire on the user's typical input shape. */
    Expr* pre_args[1] = { expr_copy((Expr*)e) };
    Expr* pre_call = expr_new_function(
        expr_new_symbol(SYM_Expand), pre_args, 1);
    Expr* pre_expanded = eval_and_free(pre_call);
    if (!pre_expanded) {
        if (dbg) simp_debug_log("PythagCanon", e, best,
                                simp_debug_elapsed_ms(t0));
        return best;
    }

    /* Twelve directions: the four base Cos/Sin/Cosh/Sinh substitutions
     * plus the reciprocal-pair pairs (Sec<->Tan, Csc<->Cot, Sech<->Tanh,
     * Csch<->Coth).  Per-direction cost is one ReplaceRepeated walk +
     * one Expand; rules whose LHS head isn't present in `e` no-op out
     * almost instantly because the matcher rejects on head mismatch.
     * The keep-best-strict-win selection in the body below means firing
     * a no-op direction is only an ms-scale speed cost, never a
     * correctness cost. */
    Expr* rule_sets[12] = {
        rules_to_sin,  rules_to_cos,  rules_to_sinh, rules_to_cosh,
        rules_to_tan,  rules_to_sec,  rules_to_cot,  rules_to_csc,
        rules_to_tanh, rules_to_sech, rules_to_coth, rules_to_csch
    };
    for (int i = 0; i < 12; i++) {
        if (!rule_sets[i]) continue;
        Expr* ra_args[2] = { expr_copy(pre_expanded),
                              expr_copy(rule_sets[i]) };
        Expr* ra_call = expr_new_function(
            expr_new_symbol(SYM_ReplaceRepeated), ra_args, 2);
        Expr* substituted = eval_and_free(ra_call);
        if (!substituted) continue;
        if (expr_eq(substituted, pre_expanded)) {
            expr_free(substituted);
            continue;
        }
        Expr* exp_args[1] = { substituted };
        Expr* exp_call = expr_new_function(
            expr_new_symbol(SYM_Expand), exp_args, 1);
        Expr* expanded = eval_and_free(exp_call);
        if (!expanded) continue;
        size_t s = score_with_func(expanded, NULL);
        if (s < best_score) {
            expr_free(best);
            best = expanded;
            best_score = s;
        } else {
            expr_free(expanded);
        }
    }
    expr_free(pre_expanded);
    if (dbg) simp_debug_log("PythagCanon", e, best,
                            simp_debug_elapsed_ms(t0));
    return best;
}

Expr* transform_pythag_canon(const Expr* e) {
    return simp_memo_wrap(e, "$PythagCanon", transform_pythag_canon_impl);
}


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
/* Trig/exp roundtrip composite                                            */
/* ----------------------------------------------------------------------- */

Expr* transform_trig_roundtrip(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Two-level memo lookup.
     *
     * Level 1 (cheap): keyed on the raw input `e`.  Lets repeated
     * calls on the *same* expression short-circuit the entire
     * pipeline including the expensive TrigToExp stage.  This is
     * the common case during candidate-set iteration -- the same
     * sub-expression flows through many transforms, and most of
     * the time the TrigRoundtrip result is unchanged.
     *
     * Level 2 (canonical): keyed on TrigToExp(e).  Catches
     * equivalent forms (e.g., `Cos[x]^2 Sec[2x]` and
     * `1/4 Sec[2x] (2 + 2 Cos[2x])`) which collapse to the same
     * exponential expression.  Pays the TrigToExp cost (which we'd
     * incur for stage 1 anyway), but saves the rest of the pipeline.
     *
     * On a miss at both levels, the result is stored under BOTH
     * the raw and canonical keys, so future identical AND
     * equivalent calls hit Level 1 / Level 2 respectively. */
    FactorMemo* memo = factor_memo_active();
    Expr* raw_key = NULL;
    if (memo) {
        Expr* raw_args[1] = { expr_copy((Expr*)e) };
        raw_key = expr_new_function(expr_new_symbol("TrigRoundtrip"),
                                    raw_args, 1);
        const Expr* hit = factor_memo_lookup(memo, raw_key);
        if (hit) {
            Expr* cached = expr_copy((Expr*)hit);
            expr_free(raw_key);
            if (dbg) simp_debug_log("TrigRoundtrip", e, cached,
                                    simp_debug_elapsed_ms(t0));
            return cached;
        }
    }

    /* Stage 1 of the pipeline: convert trig atoms to exponential form. */
    Expr* a = call_unary_copy("TrigToExp", e);

    /* Level 2 lookup keyed on TrigToExp(input). */
    Expr* canon_key = NULL;
    if (memo) {
        Expr* canon_args[1] = { expr_copy(a) };
        canon_key = expr_new_function(expr_new_symbol("TrigRoundtrip"),
                                      canon_args, 1);
        const Expr* hit = factor_memo_lookup(memo, canon_key);
        if (hit) {
            Expr* cached = expr_copy((Expr*)hit);
            /* Promote to Level 1 for next time the same `e` arrives. */
            if (raw_key) {
                factor_memo_store(memo, raw_key, cached);
                expr_free(raw_key);
            }
            expr_free(canon_key);
            expr_free(a);
            if (dbg) simp_debug_log("TrigRoundtrip", e, cached,
                                    simp_debug_elapsed_ms(t0));
            return cached;
        }
    }
    /* Explosion guard: TrigToExp is structurally expanding -- a
     * single `Cos[x] Cos[y]` (complexity 7) maps to a sum of four
     * exponentials (complexity 77, 11x growth).  Together / Cancel /
     * ExpToTrig on that intermediate is expensive AND the final
     * result tends to use Cosh / Sinh of imaginary arguments rather
     * than Cos / Sin, leaving us with a complex-coefficient form
     * that's worse for the simp candidate-set search than the
     * input.
     *
     * If TrigToExp expanded the input by more than 5x, abort the
     * round-trip: skip the slow Together / Cancel / ExpToTrig stages
     * and return the input unchanged.  Other transforms in the
     * candidate set still see the input form.
     *
     * Verified safe on the user-reference case (Sin[x]^3 + Sin[3x] -
     * 3 Sin[x] expands by ~3x at TrigToExp stage but still benefits
     * from the round-trip).  Triggers on inputs like Cos[x] Cos[y]
     * where TrigToExp blows up 11x. */
    size_t in_score = simp_default_complexity(e);
    size_t exp_score = simp_default_complexity(a);
    if (dbg) {
        fprintf(stderr, "  TrigRoundtrip complexity: in=%zu exp=%zu ratio=%.2f\n",
                in_score, exp_score,
                in_score > 0 ? (double)exp_score / in_score : 0.0);
    }
    Expr* d;
    if (in_score > 0 && exp_score > 5 * in_score) {
        expr_free(a);
        d = expr_copy((Expr*)e);
    } else {
        Expr* b = call_unary_owned("Together", a);
        Expr* c = call_unary_owned("Cancel", b);
        d = call_unary_owned("ExpToTrig", c);
    }

    /* Store under both keys so future identical or canonically-equivalent
     * calls hit the appropriate level. */
    if (raw_key) {
        factor_memo_store(memo, raw_key, d);
        expr_free(raw_key);
    }
    if (canon_key) {
        factor_memo_store(memo, canon_key, d);
        expr_free(canon_key);
    }

    if (dbg) simp_debug_log("TrigRoundtrip", e, d, simp_debug_elapsed_ms(t0));
    return d;
}

/* Roots-of-unity simplification.
 *
 * Recognises every (-1)^(p/q) and E^(I p Pi / q) atom in the input,
 * lifts the expression to a univariate polynomial in
 *   omega = (-1)^(1/Q),  Q = LCM of denominators,
 * reduces modulo the cyclotomic polynomial Phi_{2Q}(omega) (the minimal
 * polynomial of omega = e^(I Pi / Q) over Q), and substitutes back. The
 * reduction is exact: omega is a primitive (2Q)-th root of unity, so
 * Phi_{2Q}(omega) = 0, and any polynomial p(omega) is identically zero
 * iff Phi_{2Q}(x) divides p(x). The substitute -> reduce -> substitute
 * round-trip preserves correctness for any polynomial in omega regardless
 * of the choice of free coefficients.
 *
 * Handles e.g.
 *   1 - (-1)^(1/3) + (-1)^(2/3)                 -> 0
 *   1 - (-1)^(1/5) + (-1)^(2/5) - ... + (-1)^(4/5) -> 0
 *   3 + 2 E^(-2 I Pi/3) + 2 E^(2 I Pi/3)        -> 1
 *
 * Implemented as a small Mathematica-syntax helper installed lazily
 * into the symbol table on first call. The cyclotomic polynomial is
 * computed on-the-fly by recursive division: Phi_n(x) = (x^n - 1) /
 * Prod_{d | n, d < n} Phi_d(x). Cache pressure is light because the
 * recursion is bounded by the LCM 2Q (typically < 30 for hand-written
 * inputs) and PolynomialQuotient memoises subresults via the term
 * structure of x^n - 1. */
void simp_install_roots_of_unity_helpers(void) {
    static bool installed = false;
    if (installed) return;
    /* Definitions are added as DownValues on internal `$ru*` symbols so
     * they don't shadow anything user-visible. parse_expression returns
     * a SetDelayed Expr*; evaluate runs the assignment and returns Null
     * (we free that). */
    const char* defs[] = {
        "$ruCyclotomic[1, x_] := x - 1",
        "$ruCyclotomic[n_Integer, x_] := Module["
        "  {d, num = x^n - 1, denom = 1},"
        "  Do[If[Mod[n, d] == 0, denom = denom * $ruCyclotomic[d, x]], {d, 1, n - 1}];"
        "  PolynomialQuotient[num, denom, x]]",
        /* Main simplifier: collect denominators, lift to polynomial in
         * $ru, reduce mod Phi_{2Q}($ru), substitute back. The mod 2Q
         * normalisation on the exponent handles negative-exponent forms
         * like E^(-I Pi p / q) without leaving x^(-k) terms that
         * PolynomialRemainder would reject. */
        "$ruSimplify[expr_] := Module["
        "  {denoms, Q, polyForm, phiPoly, reduced},"
        "  denoms = Union[Join["
        "    Cases[expr, Power[-1, Rational[_, q_]] :> q, {0, Infinity}],"
        "    Cases[expr, Power[E, Times[Complex[0, Rational[_, q_]], Pi]] :> q, {0, Infinity}]]];"
        "  If[denoms === {}, expr,"
        "    Q = Apply[LCM, denoms];"
        "    polyForm = expr /. {"
        "      Power[-1, Rational[a_, b_]] :> $ru^Mod[a Q/b, 2 Q],"
        "      Power[E, Times[Complex[0, Rational[a_, b_]], Pi]] :> $ru^Mod[a Q/b, 2 Q]};"
        "    phiPoly = $ruCyclotomic[2 Q, $ru];"
        "    reduced = PolynomialRemainder[polyForm, phiPoly, $ru];"
        "    reduced /. $ru -> Power[-1, 1/Q]]]"
    };
    for (size_t i = 0; i < sizeof(defs)/sizeof(defs[0]); i++) {
        Expr* parsed = parse_expression(defs[i]);
        if (!parsed) continue;
        Expr* r = eval_and_free(parsed);
        if (r) expr_free(r);
    }
    installed = true;
}

Expr* simp_roots_of_unity(const Expr* e) {
    simp_install_roots_of_unity_helpers();
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* args[1] = { expr_copy((Expr*)e) };
    Expr* call = expr_new_function(
        expr_new_symbol("$ruSimplify"), args, 1);
    Expr* out = eval_and_free(call);
    if (dbg) simp_debug_log("RootsOfUnity", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* Pythagorean perfect-square completion: 1 +/- 2 Sin[x] Cos[x]
 * = (Sin[x] +/- Cos[x])^2. Lets Simplify reach factored forms like
 * (Sin + Cos)^4 from a Factor result of (1 + 2 Sin Cos)^2. We keep
 * this as its own transform (separate from TrigFactor) because
 * TrigFactor's identity rule list also contains the linear-combination
 * rule a Sin[x] + b Cos[x] -> Sqrt[a^2+b^2] Sin[x + ArcTan[a, b]],
 * which would re-rewrite (Sin + Cos) into a single trig and obscure the
 * factored form. As a standalone seed the rewrite produces a candidate
 * Simplify can score directly. */
/* Forward declaration -- definition lives near transform_pythag_reduce. */
Expr* simp_memo_wrap(const Expr* e, const char* pseudo_head,
                            Expr* (*impl)(const Expr*));
bool has_pythag_head(const Expr* e);
bool has_non_integer_power(const Expr* e);
bool is_rational_literal(const Expr* e);

static Expr* transform_pythag_square_complete_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ 1 + 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] + Cos[x])^2 + r, "
            "  1 - 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] - Cos[x])^2 + r, "
            "  -1 + Cosh[x_]^2 + Sinh[x_]^2 + 2 Sinh[x_] Cosh[x_] + r___ "
            "      :> (Sinh[x] + Cosh[x])^2 - 1 + r, "
            "  -1 + Cosh[x_]^2 + Sinh[x_]^2 - 2 Sinh[x_] Cosh[x_] + r___ "
            "      :> (Sinh[x] - Cosh[x])^2 - 1 + r }");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Same fast-skip as PythagReduce: every rule LHS contains a
     * Cos/Sin/Cosh/Sinh head, so on inputs without any of those
     * the ReplaceRepeated walk finds nothing. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagSquareComplete", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceRepeated"), args, 2);
    Expr* out = eval_and_free(call);
    if (dbg) simp_debug_log("PythagSquareComplete", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

Expr* transform_pythag_square_complete(const Expr* e) {
    return simp_memo_wrap(e, "$PythagSquareComplete",
                          transform_pythag_square_complete_impl);
}

/* Half-angle tangent identity, applied to both circular and hyperbolic
 * functions. Folds the Weierstrass forms
 *
 *   Sin[x] / (1 + Cos[x])              -> Tan[x/2]
 *   Sin[x]^a (1 + Cos[x])^(-a)         -> Tan[x/2]^a
 *   Sin[x] / (c (1 + Cos[x]))          -> Tan[x/2] / c        (FreeQ[c, x])
 *   Sin[x]^a (c (1 + Cos[x]))^(-a)     -> Tan[x/2]^a / c^a    (FreeQ[c, x])
 *   (1 - Cos[x]) / Sin[x]              -> Tan[x/2]
 *   (1 - Cos[x])^a Sin[x]^(-a)         -> Tan[x/2]^a
 *
 * and the analogous Sinh/Cosh -> Tanh[x/2] family (with the sign
 * difference (Cosh[x] - 1)/Sinh[x] == Tanh[x/2]). Each rule has a
 * trailing `r___` BlankNullSequence inside the Times so the rule fires
 * on subterms inside larger products (e.g. (1/2) Sin[x] / (1 + Cos[x])
 * still rewrites). The conditional-pattern guards (a + b === 0,
 * FreeQ[c, x]) are what keeps the rules general -- there are no
 * specific-numeric variants.
 *
 * Output complexity is uniformly less than or equal to the input on
 * every shape that fires, so simp_search's leaf-count tiebreak takes
 * the rewritten form. */
static Expr* transform_halfangle_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ "
            /* Trig: Sin / (1 + Cos) -> Tan[x/2] */
            "  Sin[x_] Power[1 + Cos[x_], -1] r___ :> Tan[x/2] r, "
            "  Sin[x_]^a_ (1 + Cos[x_])^b_ r___ "
            "    /; a + b === 0 :> Tan[x/2]^a r, "
            /* Trig: Sin / (c (1 + Cos)) -> Tan[x/2]/c */
            "  Sin[x_] Power[c_ + c_ Cos[x_], -1] r___ "
            "    /; FreeQ[c, x] :> Tan[x/2] c^(-1) r, "
            "  Sin[x_]^a_ (c_ + c_ Cos[x_])^b_ r___ "
            "    /; a + b === 0 && FreeQ[c, x] :> Tan[x/2]^a c^b r, "
            /* Trig: (1 - Cos) / Sin -> Tan[x/2] */
            "  (1 - Cos[x_]) Power[Sin[x_], -1] r___ :> Tan[x/2] r, "
            "  (1 - Cos[x_])^a_ Sin[x_]^b_ r___ "
            "    /; a + b === 0 :> Tan[x/2]^a r, "
            /* Hyperbolic: Sinh / (1 + Cosh) -> Tanh[x/2] */
            "  Sinh[x_] Power[1 + Cosh[x_], -1] r___ :> Tanh[x/2] r, "
            "  Sinh[x_]^a_ (1 + Cosh[x_])^b_ r___ "
            "    /; a + b === 0 :> Tanh[x/2]^a r, "
            /* Hyperbolic: Sinh / (c (1 + Cosh)) -> Tanh[x/2]/c */
            "  Sinh[x_] Power[c_ + c_ Cosh[x_], -1] r___ "
            "    /; FreeQ[c, x] :> Tanh[x/2] c^(-1) r, "
            "  Sinh[x_]^a_ (c_ + c_ Cosh[x_])^b_ r___ "
            "    /; a + b === 0 && FreeQ[c, x] :> Tanh[x/2]^a c^b r, "
            /* Hyperbolic: (Cosh - 1) / Sinh -> Tanh[x/2] */
            "  (-1 + Cosh[x_]) Power[Sinh[x_], -1] r___ :> Tanh[x/2] r, "
            "  (-1 + Cosh[x_])^a_ Sinh[x_]^b_ r___ "
            "    /; a + b === 0 :> Tanh[x/2]^a r "
            "}");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Every HalfAngle rule LHS uses Sin/Cos or Sinh/Cosh.  Skip the
     * ReplaceRepeated walk on inputs without any of those heads. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("HalfAngle", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceRepeated"), args, 2);
    Expr* out = eval_and_free(call);
    if (dbg) simp_debug_log("HalfAngle", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

Expr* transform_halfangle(const Expr* e) {
    return simp_memo_wrap(e, "$HalfAngle", transform_halfangle_impl);
}


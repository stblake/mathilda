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
/* Heuristic search                                                        */
/* ----------------------------------------------------------------------- */

/* Forward declarations: defined below alongside the rest of the
 * simp_factorial cluster. transform_can_fire (a few hundred lines down)
 * needs to ask "does this input contain a Factorial atom?" before
 * firing the FactorialRules seed, and simp_search uses simp_eq_head_sym
 * to count factorials in candidates -- both must be visible here. */
bool contains_factorial(const Expr* e);
bool simp_eq_head_sym(const Expr* e, const char* name);

static const char* SIMP_TRANSFORMS[] = {
    "Together",
    "Cancel",
    "Expand",
    "ExpandNumerator",
    "ExpandDenominator",
    "Factor",
    "FactorSquareFree",
    "FactorTerms",
    "Apart",
    "TrigExpand",
    "TrigFactor",
    /* TrigReduce shrinks angle-addition forms (Sin[a]Cos[b]+Cos[a]Sin[b]
     * -> Sin[a+b]) and pulls integer powers / products of single-arg
     * trig calls into single trig calls of compound arguments. The
     * round-loop's score-gate keeps the reduction only when its leaf
     * count is no greater than the parent's, which is the typical case
     * for genuine angle-addition shapes. */
    "TrigReduce",
    /* TrigToExp surfaces the exp form directly so that hyperbolic
     * combinations whose exp form is strictly simpler (e.g. Sinh[x] +
     * Cosh[x] -> E^x) win the score tiebreak. The full trig roundtrip
     * (TrigToExp -> Together -> Cancel -> ExpToTrig) ends with ExpToTrig
     * and converts E^x back to Cosh[x] + Sinh[x], hiding the simpler
     * intermediate; offering TrigToExp as its own seed avoids that. For
     * pure trig inputs, TrigToExp yields a complex-coefficient exp form
     * with strictly higher complexity than the original, so the original
     * still wins -- no regression on Sin/Cos identities. */
    "TrigToExp"
};
static const size_t SIMP_TRANSFORM_COUNT =
    sizeof(SIMP_TRANSFORMS) / sizeof(SIMP_TRANSFORMS[0]);

/* Returns true if the expression contains any Power with a non-integer
 * exponent (e.g. Sqrt forms, Rational exponents, symbolic exponents).
 * Mathilda's Factor / FactorSquareFree call its trial-division loop in
 * factor_roots which can stall on multivariate inputs that include such
 * Power atoms, so we skip those transforms when this returns true. */
bool has_non_integer_power(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type != EXPR_INTEGER && exp->type != EXPR_BIGINT) return true;
    }
    if (has_non_integer_power(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_non_integer_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* Centralised cheap-precondition gate. Returns false only when the
 * predicate proves the named transform cannot possibly fire on `e` (and
 * `ctx`, where applicable). Conservative: returns true if uncertain.
 *
 * Names cover both the SIMP_TRANSFORMS table entries and the open-coded
 * seed/round transforms ("AbsRules", "LogExpRules", "AssumptionRules",
 * "TrigRoundtrip", "CollectPerVariable") so every call site shares one
 * dispatch table.
 *
 * `ctx` may be NULL; transforms that don't consult it ignore the parameter. */
bool transform_can_fire(const char* name, const Expr* e,
                               const AssumeCtx* ctx) {
    /* Polynomial machinery: Factor / FactorSquareFree / FactorTerms / TrigFactor
     * stall on non-integer Power exponents (Sqrt, Rational exponents, ...). */
    if (strcmp(name, "Factor") == 0 ||
        strcmp(name, "FactorSquareFree") == 0 ||
        strcmp(name, "FactorTerms") == 0 ||
        strcmp(name, "TrigFactor") == 0) {
        if (has_non_integer_power(e)) return false;
    }
    /* Polynomial / partial-fraction transforms are no-ops on numeric-only
     * inputs (no symbol leaves). */
    if (strcmp(name, "Factor") == 0 ||
        strcmp(name, "FactorSquareFree") == 0 ||
        strcmp(name, "FactorTerms") == 0 ||
        strcmp(name, "Apart") == 0) {
        if (!contains_variable(e)) return false;
    }
    /* Trig family: skip when there's no trig or hyperbolic head anywhere.
     * The roundtrip composite is gated identically. */
    if (strcmp(name, "TrigExpand") == 0 ||
        strcmp(name, "TrigFactor") == 0 ||
        strcmp(name, "TrigToExp") == 0 ||
        strcmp(name, "TrigRoundtrip") == 0 ||
        strcmp(name, "TrigReduce") == 0) {
        if (!contains_trig_or_hyperbolic(e)) return false;
    }
    /* TrigReduce additionally needs Plus or Times anywhere in the input;
     * a single trig call with no enclosing arithmetic produces no
     * product-to-sum work and the rule list is a no-op. Without this
     * gate every Sin/Cos seed would seed a TrigReduce candidate and
     * waste a memo slot per leaf. */
    if (strcmp(name, "TrigReduce") == 0) {
        if (!contains_plus_or_times(e)) return false;
    }
    /* Abs rules. */
    if (strcmp(name, "AbsRules") == 0) {
        if (!contains_abs(e)) return false;
    }
    /* Sqrt[_^2] rules. */
    if (strcmp(name, "SqrtSquareRules") == 0) {
        if (!contains_sqrt_of_square(e)) return false;
    }
    /* LogExp identity cascade: requires Log or Power somewhere (the
     * positivity-aware rewrites — Log[a*b] -> Log[a]+Log[b] etc. — are
     * gated on ctx facts internally; the unconditional Exp-distribute
     * rewrite Power[E, Plus[...]] -> Product Power[E, ti] fires
     * regardless and is needed for `Exp[k(Log a + Log b)] -> a^k b^k`
     * with integer k where no positivity assumption is required). */
    if (strcmp(name, "LogExpRules") == 0) {
        if (!contains_log(e) && !contains_power(e)) return false;
    }
    /* simp_log rewriter: nothing fires without a Log somewhere. */
    if (strcmp(name, "SimpLogRules") == 0) {
        if (!contains_log(e)) return false;
    }
    /* Assumption rewriter: nothing fires without facts. */
    if (strcmp(name, "AssumptionRules") == 0) {
        if (!ctx_has_facts(ctx)) return false;
    }
    /* Per-variable Collect: meaningful only when there are at least two
     * distinct variables to choose between. */
    if (strcmp(name, "CollectPerVariable") == 0) {
        if (expr_variables_count_capped(e, 2) < 2) return false;
    }
    /* Factorial rewriter: nothing fires unless the input mentions
     * Factorial somewhere in the tree. */
    if (strcmp(name, "FactorialRules") == 0) {
        if (!contains_factorial(e)) return false;
    }
    return true;
}

/* Score a candidate; if it beats the running best, replace best. The
 * candidate `c` is *borrowed* (caller still owns the source); `best` is
 * a slot the caller manages. */
static void update_best(Expr** best, size_t* best_score, const Expr* c,
                        const Expr* complexity_func) {
    size_t s = score_with_func(c, complexity_func);
    if (s < *best_score) {
        expr_free(*best);
        *best = expr_copy((Expr*)c);
        *best_score = s;
    }
}

/* True when `best` is the literal integer 0 -- the canonical
 * "input was an identity" outcome. Once Simplify reaches 0 there
 * is nothing left to improve, but the round loop would otherwise
 * keep grinding through the remaining transforms (Factor, TrigFactor,
 * etc.) on the still-large seeds. Early-exiting on this hit
 * trims the (Sin[x]+Cos[x])^4 - (1+Sin[2x])^2 case from ~220 ms to a
 * few ms by stopping as soon as TrigReduce produces 0.
 *
 * Restricted to literal 0 (rather than "any atomic with score 1") so
 * a custom complexity_func that ranks structurally-larger forms below
 * an atomic still gets full search coverage. The literal-zero case
 * is unambiguous: no transform we apply can produce a structurally
 * simpler form than 0. */
static bool simp_best_is_zero(const Expr* best) {
    return best && best->type == EXPR_INTEGER && best->data.integer == 0;
}

/* Apply Collect[expr, v] for each free variable v of expr, scoring the
 * results against `best` and adding novel ones to `next`. Collect can
 * surface a more compact form by grouping like powers (e.g. it can
 * recover x*(a+b) from a*x + b*x), and which variable to collect by is
 * not knowable up front -- Mathematica's Simplify likewise tries each
 * variable. We rely on Variables[] to enumerate the candidates. */
static void try_collect_per_variable(const Expr* seed, size_t parent_score,
                                     CandSet* next,
                                     Expr** best, size_t* best_score,
                                     const Expr* complexity_func) {
    Expr* vars = call_unary_copy("Variables", seed);
    if (!vars) return;
    if (vars->type != EXPR_FUNCTION ||
        !vars->data.function.head ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List) {
        expr_free(vars);
        return;
    }
    size_t nv = vars->data.function.arg_count;
    bool dbg = simp_debug_enabled();
    for (size_t i = 0; i < nv; i++) {
        Expr* v = vars->data.function.args[i];
        Expr* args[2] = { expr_copy((Expr*)seed), expr_copy(v) };
        Expr* call = expr_new_function(expr_new_symbol("Collect"), args, 2);
        clock_t t0 = dbg ? clock() : 0;
        Expr* r = evaluate(call);
        expr_free(call);
        if (dbg) {
            char* vname = expr_to_string(v);
            char buf[64];
            snprintf(buf, sizeof(buf), "Collect[_,%s]",
                     vname ? vname : "?");
            simp_debug_log(buf, seed, r, simp_debug_elapsed_ms(t0));
            free(vname);
        }
        if (!r) continue;
        update_best(best, best_score, r, complexity_func);
        /* Branch-pruning: only propagate if the candidate is no worse
         * than its parent. Strictly worse forms cannot lead to a better
         * result through any further unary transform we run. */
        if (expr_eq(r, seed)) {
            expr_free(r);
        } else if (score_with_func(r, complexity_func) > parent_score) {
            expr_free(r);
        } else {
            cs_add_or_free(next, r);
        }
    }
    expr_free(vars);
}

/* ----------------------------------------------------------------------- */
/* Shape classifier + pipeline dispatch                                    */
/* ----------------------------------------------------------------------- */

/*
 * simp_classify performs one O(n) walk over the input and assigns it to a
 * shape bucket. simp_dispatch then routes to a specialised pipeline that
 * runs only the transforms relevant for that shape. The classifier is
 * conservative: borderline inputs fall through to SIMP_SHAPE_GENERAL,
 * which calls simp_search (the original full search), so misclassification
 * cannot change behaviour, only performance.
 *
 * Priority order (first matching wins):
 *   1. TRIG       -- any Sin/Cos/.../Sinh/.../ArcXxx head present
 *   2. LOGEXP     -- Log present (after trig is excluded)
 *   3. POLYNOMIAL -- no trig/log/abs, every Power has integer exponent,
 *                    PolynomialQ over Variables[input] is True
 *   4. RATIONAL   -- no trig/log/abs, Together[input] is poly/poly
 *   5. GENERAL    -- everything else
 *
 * Inputs containing Abs route to GENERAL: Abs rewrites are best handled
 * by the existing search machinery (the bottom-up walker hits Abs heads
 * directly via apply_abs_rules), and a specialised Abs pipeline buys
 * little.
 */
/* SimpShape is declared in simp_internal.h (simp_builtins.c uses it
 * for the SIMP_SHAPE_RATIONAL fast path). */

/* Helper: PolynomialQ[e, Variables[e]] using the existing builtins. */
static bool simp_is_polynomial_in_own_vars(const Expr* e) {
    Expr* vars = call_unary_copy("Variables", e);
    if (!vars) return false;
    if (vars->type != EXPR_FUNCTION ||
        !vars->data.function.head ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List) {
        expr_free(vars);
        return false;
    }
    /* Zero-variable input: a numeric literal. PolynomialQ returns True
     * trivially, but a numeric leaf doesn't benefit from the polynomial
     * pipeline (no Factor / Collect target), so report false to fall
     * through to GENERAL. */
    if (vars->data.function.arg_count == 0) {
        expr_free(vars);
        return false;
    }
    Expr* args[2] = { expr_copy((Expr*)e), vars };
    Expr* pq = expr_new_function(expr_new_symbol("PolynomialQ"), args, 2);
    Expr* r = evaluate(pq);
    expr_free(pq);
    bool ok = (r && r->type == EXPR_SYMBOL &&
               r->data.symbol == SYM_True);
    if (r) expr_free(r);
    return ok;
}

/* Helper: build Together[e], extract Numerator/Denominator, check both
 * are polynomial in their own variables. Returns false (and frees nothing
 * external) if any step fails.
 *
 * Numeric collapse: when Together[e] reduces to a numeric leaf (or
 * Rational[p,q]), Numerator/Denominator yield zero-variable
 * expressions which simp_is_polynomial_in_own_vars rejects (its
 * empty-Variables clause is intentional for the polynomial classifier
 * but spurious here). Treat that case as trivially rational so the
 * shortcut still fires -- e.g. ((x-y)/(x+y) - (x+y)/(x-y)) /
 * (1 - (x^2-x*y-y^2)/(x^2-y^2)) collapses to -4 and routes to the
 * rational pipeline instead of falling through to GENERAL/bottomup. */
static bool simp_is_rational(const Expr* e) {
    Expr* tg = call_unary_copy("Together", e);
    if (!tg) return false;
    bool tg_is_numeric_leaf =
        tg->type == EXPR_INTEGER ||
        tg->type == EXPR_BIGINT ||
        tg->type == EXPR_REAL ||
        is_rational_literal(tg);
    if (tg_is_numeric_leaf) {
        expr_free(tg);
        return true;
    }
    Expr* num = call_unary_copy("Numerator", tg);
    Expr* den = call_unary_copy("Denominator", tg);
    bool ok = false;
    if (num && den &&
        !has_non_integer_power(num) &&
        !has_non_integer_power(den) &&
        simp_is_polynomial_in_own_vars(num) &&
        simp_is_polynomial_in_own_vars(den)) {
        ok = true;
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(tg);
    return ok;
}

SimpShape simp_classify(const Expr* e) {
    if (contains_trig_or_hyperbolic(e)) return SIMP_SHAPE_TRIG;
    if (contains_abs(e))                return SIMP_SHAPE_GENERAL;
    /* Factorial-bearing inputs MUST flow through the general pipeline so
     * the FactorialRules seed in simp_search can fire. The rational and
     * polynomial pipelines short-circuit past it, leaving Factorial[arg]
     * structurally opaque to Together/Cancel and missing the shift-
     * normalization that closes (n+3)!/(n+1)! -> (n+2)(n+3) etc. */
    if (contains_factorial(e))          return SIMP_SHAPE_GENERAL;
    if (contains_log(e))                return SIMP_SHAPE_LOGEXP;
    /* No trig, no abs, no log. Decide poly / rational / general. */
    if (!has_non_integer_power(e) && simp_is_polynomial_in_own_vars(e)) {
        return SIMP_SHAPE_POLYNOMIAL;
    }
    if (simp_is_rational(e)) return SIMP_SHAPE_RATIONAL;
    return SIMP_SHAPE_GENERAL;
}

/* Forward declaration: simp_search is the GENERAL pipeline. */
Expr* simp_search(const Expr* original_input, const AssumeCtx* ctx,
                         const Expr* complexity_func);

/* Forward declaration for simp_factorial (defined alongside simp_bottomup
 * further down). contains_factorial was forward-declared above. */
Expr* simp_factorial(const Expr* e);

/* ----------------------------------------------------------------------- */
/* Specialised pipelines                                                   */
/* ----------------------------------------------------------------------- */

/*
 * Polynomial pipeline. Skips every transform that targets trig, log,
 * abs, rational, or assumption-driven structure -- none of which can
 * fire on a SHAPE_POLYNOMIAL input. The candidate set is:
 *   - the input itself
 *   - Expand[input]
 *   - Factor[input], FactorTerms[input]   (gated by has_non_integer_power)
 *   - Collect[input, v]   for each variable v
 *   - Collect[Expand[input], v]   for each variable v
 *
 * The Expand-then-Collect path is the one that recovers c + x*(a+b)
 * from a*x + b*x + c (test_simplify_collect_by_variable). Factor wins
 * on cases like (x-1)*(x+1)*(x^2+1)+1 -> x^4 because it emits the
 * already-expanded form when no nontrivial factorisation exists, but
 * it ties on score and the Expand candidate wins by SimplifyCount.
 *
 * No round loop: every winning move on a polynomial is reachable in
 * one application of the listed transforms; iterating doesn't help.
 */
static Expr* simp_pipeline_polynomial(const Expr* input,
                                      const AssumeCtx* ctx,
                                      const Expr* complexity_func) {
    (void)ctx;
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    Expr* expanded = traced_call_unary("Expand", input);
    if (expanded) update_best(&best, &bs, expanded, complexity_func);

    if (!has_non_integer_power(input)) {
        Expr* factored = traced_call_unary("Factor", input);
        if (factored) update_best(&best, &bs, factored, complexity_func);
        if (factored) expr_free(factored);
        Expr* fterms = traced_call_unary("FactorTerms", input);
        if (fterms) update_best(&best, &bs, fterms, complexity_func);
        if (fterms) expr_free(fterms);
    }

    /* Per-variable Collect on input AND on the expanded form. The
     * Expand-then-Collect path catches cases the input-Collect misses,
     * because Collect groups like powers of v and a Plus-of-Times input
     * already in factored shape doesn't expose them. */
    if (transform_can_fire("CollectPerVariable", input, NULL)) {
        CandSet next; cs_init(&next);
        try_collect_per_variable(input, bs, &next, &best, &bs, complexity_func);
        if (expanded) {
            try_collect_per_variable(expanded, bs, &next, &best, &bs, complexity_func);
        }
        cs_free(&next);
    }

    if (expanded) expr_free(expanded);
    return best;
}

/*
 * Rational pipeline. For inputs whose Together-form is poly/poly with
 * no trig/log/abs heads. The candidate set covers:
 *   - Together  (canonical fraction form)
 *   - Cancel    (gcd reduction)
 *   - ExpandNumerator / ExpandDenominator
 *   - Apart     (partial-fraction decomposition)
 *   - Factor on the Cancel'd form
 *   - Per-variable Collect on the Cancel'd form
 *   - AssumptionRules force-take, when ctx has facts (covers cases like
 *     (1-a^2)/b^2 with a^2+b^2==1)
 *
 * No round loop. Every meaningful win on a rational expression is one
 * of these transforms applied directly.
 */
static Expr* simp_pipeline_rational(const Expr* input,
                                    const AssumeCtx* ctx,
                                    const Expr* complexity_func) {
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    /* AssumptionRules force-take, mirroring the seed-phase semantics in
     * simp_search: assumption-driven rewrites are correctness-preserving
     * under the assumption set and qualitatively "more simplified", so we
     * accept them at equal or lower complexity. */
    Expr* assum_seed = NULL;
    if (transform_can_fire("AssumptionRules", input, ctx)) {
        Expr* ar = apply_assumption_rules(input, ctx);
        if (ar && !expr_eq(ar, input)) {
            size_t s = score_with_func(ar, complexity_func);
            if (s <= bs) {
                expr_free(best);
                best = expr_copy(ar);
                bs = s;
            }
            assum_seed = ar;  /* keep for downstream Together/Cancel */
        } else if (ar) {
            expr_free(ar);
        }
    }

    const Expr* seed = assum_seed ? assum_seed : input;

    Expr* tg = traced_call_unary("Together", seed);
    if (tg) update_best(&best, &bs, tg, complexity_func);

    Expr* cn_src = tg ? tg : (Expr*)seed;
    Expr* cn = traced_call_unary("Cancel", cn_src);
    if (cn) update_best(&best, &bs, cn, complexity_func);

    Expr* en = traced_call_unary("ExpandNumerator", cn ? cn : (Expr*)seed);
    if (en) update_best(&best, &bs, en, complexity_func);

    Expr* ed = traced_call_unary("ExpandDenominator", cn ? cn : (Expr*)seed);
    if (ed) update_best(&best, &bs, ed, complexity_func);

    Expr* ap = traced_call_unary("Apart", seed);
    if (ap) update_best(&best, &bs, ap, complexity_func);

    /* Factor on the most-canonical form available.  Prefer the Cancel'd
     * form when it exists, but fall through to seed/input when Cancel
     * returned NULL (no change to apply) -- otherwise Simplify would
     * never try Factor on inputs like `x/(x^3 + a b + a x + b x^2)`
     * which are already cancelled.  Gated against non-integer powers,
     * which shouldn't appear on a SHAPE_RATIONAL input but is
     * defensive.
     *
     * Push a NULL factor_memo around the Factor call so it uses its
     * normal (separate num/den variable lists) behaviour rather than
     * the inside-Simplify combined-scope path.  The combined-scope
     * mode refuses to factor parametric denominators like
     * `x^3 + a b + a x + b x^2` (which factors over Z[a, b][x] as
     * (a + x^2)(b + x)); the separate-scope mode handles them
     * correctly. */
    {
        const Expr* fc_src = cn ? cn : (tg ? tg : seed);
        if (!has_non_integer_power(fc_src)) {
            factor_memo_push(NULL);
            Expr* fc = traced_call_unary("Factor", fc_src);
            factor_memo_pop();
            if (fc) {
                update_best(&best, &bs, fc, complexity_func);
                expr_free(fc);
            }
        }
    }

    /* Per-variable Collect on the most-canonical form (Cancel'd if
     * available, otherwise seed/input -- same fall-through reason as
     * the Factor branch above). */
    {
        const Expr* col_src = cn ? cn : (tg ? tg : seed);
        if (transform_can_fire("CollectPerVariable", col_src, NULL)) {
            CandSet next; cs_init(&next);
            try_collect_per_variable(col_src, bs, &next, &best, &bs, complexity_func);
            cs_free(&next);
        }
    }

    if (tg) expr_free(tg);
    if (cn) expr_free(cn);
    if (en) expr_free(en);
    if (ed) expr_free(ed);
    if (ap) expr_free(ap);
    if (assum_seed) expr_free(assum_seed);
    return best;
}

/*
 * Log/Exp pipeline. For inputs containing Log (and no trig/hyperbolic).
 * The cascade in apply_logexp_rules implements the positivity-aware
 * Log/Power identities (Log[a*b] -> Log[a]+Log[b], (a*b)^c -> a^c b^c,
 * Log[x^p] -> p Log[x], (x^p)^q -> x^(p*q)). It is force-take in
 * simp_search (correctness-preserving under positivity assumptions) and
 * we replicate that here.
 *
 * Candidate set:
 *   - LogExpRules cascade (force-take when it changes the form)
 *   - AssumptionRules (e.g. Log[Exp[x]] -> x via Log[E^x] -> x*Log[E])
 *   - Together, Cancel, Expand on the cascade output
 *
 * No round loop; no trig transforms.
 */
static Expr* simp_pipeline_logexp(const Expr* input,
                                  const AssumeCtx* ctx,
                                  const Expr* complexity_func) {
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    /* SimpLogRules + LogExpRules + AssumptionRules iterated to fixed
     * point.  Order matters: SimpLogRules can produce a fused Log whose
     * argument is then collapsible by LogExpRules (e.g. fusion yields
     * Log[E^(-x)] which LogExpRules turns into -x under realness), and
     * the AssumptionRules pass can in turn expose a fresh
     * Log[positive_rational] for SimpLogRules. Two rounds suffice for
     * all observed test cases; cap at 3 to keep the cost bounded. */
    for (int rd = 0; rd < 3; rd++) {
        bool changed = false;
        if (transform_can_fire("SimpLogRules", best, ctx)) {
            clock_t t0 = simp_debug_enabled() ? clock() : 0;
            Expr* sl = simp_log_apply(best, ctx);
            if (simp_debug_enabled()) {
                simp_debug_log("SimpLogRules", best, sl,
                               simp_debug_elapsed_ms(t0));
            }
            if (sl && !expr_eq(sl, best)) {
                expr_free(best);
                best = expr_copy(sl);
                bs = score_with_func(best, complexity_func);
                changed = true;
            }
            if (sl) expr_free(sl);
        }
        if (transform_can_fire("LogExpRules", best, ctx)) {
            Expr* lr = apply_logexp_rules(best, ctx);
            if (lr && !expr_eq(lr, best)) {
                expr_free(best);
                best = expr_copy(lr);
                bs = score_with_func(best, complexity_func);
                changed = true;
            }
            if (lr) expr_free(lr);
        }
        if (!changed) break;
    }

    /* Assumption-driven rewrites on the (possibly rewritten) best.
     * Force-take semantics: an assumption-aware rewrite that actually
     * changed the form is correctness-preserving under the assumption
     * set and counts as "more simplified" even when the leaf-count
     * tiebreak is even (e.g. Log[a^p] -> p Log[a] both score 4). This
     * matches the seed-phase behaviour in simp_search. */
    if (transform_can_fire("AssumptionRules", best, ctx)) {
        Expr* ar = apply_assumption_rules(best, ctx);
        if (ar && !expr_eq(ar, best)) {
            size_t s = score_with_func(ar, complexity_func);
            if (s <= bs) {
                expr_free(best);
                best = expr_copy(ar);
                bs = s;
            }
        }
        if (ar) expr_free(ar);
    }

    /* Standard cleanup. */
    Expr* tg = traced_call_unary("Together", best);
    if (tg) update_best(&best, &bs, tg, complexity_func);
    Expr* cn = traced_call_unary("Cancel", best);
    if (cn) update_best(&best, &bs, cn, complexity_func);
    Expr* ex = traced_call_unary("Expand", best);
    if (ex) update_best(&best, &bs, ex, complexity_func);

    if (tg) expr_free(tg);
    if (cn) expr_free(cn);
    if (ex) expr_free(ex);
    return best;
}

/* ----------------------------------------------------------------------- */
/* Additive-subexpression splitter                                         */
/* ----------------------------------------------------------------------- */

/*
 * When the input is a Plus whose addends partition into >=2 disjoint
 * connected components in the free-symbol-sharing graph (e.g.
 * Cos[x] + Sin[x] + Cos[y] + Sin[y] splits into the {x} group and the
 * {y} group), simplify each component independently and return the
 * sum.
 *
 * Why: simp_search's trig transforms (Together, TrigExpand, Cancel)
 * cost super-linearly in the number of addends. When the addends
 * decompose over disjoint user-symbol sets, those pieces cannot
 * interact under any of those transforms -- simplifying each component
 * independently preserves the result while keeping each call cheap.
 * The concrete trigger was the 6-term sum
 *   Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2] + Cos[y]/Sqrt[6]/3 + Sin[y]/Sqrt[2]/3
 *     - Sqrt[6] Sin[x+Pi/6]/3 - Sqrt[6] Sin[y+Pi/6]/9
 * which hits $RecursionLimit and runs for minutes when simplified as a
 * whole; each three-term piece simplifies to 0 in ~1s on its own.
 *
 * Free-symbol semantics: an addend's "free symbols" are symbol leaves
 * appearing in its argument tree, EXCLUDING numeric constants
 * (Pi, E, EulerGamma, ...) and symbols that appear ONLY as function
 * heads. So Sin[x] contributes {x}, Cos[x+y] contributes {x, y}, but
 * Sin and Cos themselves are not counted -- otherwise every all-trig
 * addend would glom into a single component and the split would never
 * fire on inputs like the case above.
 */

typedef struct {
    const char** names;     /* borrowed pointers into the input expr */
    size_t count;
    size_t cap;
} SplitSymSet;

static void split_symset_init(SplitSymSet* s) {
    s->names = NULL; s->count = 0; s->cap = 0;
}
static void split_symset_free(SplitSymSet* s) {
    free((void*)s->names);
    s->names = NULL; s->count = 0; s->cap = 0;
}
static bool split_symset_contains(const SplitSymSet* s, const char* name) {
    for (size_t i = 0; i < s->count; i++) {
        if (strcmp(s->names[i], name) == 0) return true;
    }
    return false;
}
static void split_symset_add(SplitSymSet* s, const char* name) {
    if (split_symset_contains(s, name)) return;
    if (s->count == s->cap) {
        size_t nc = s->cap == 0 ? 4 : s->cap * 2;
        const char** nn = (const char**)realloc((void*)s->names,
                                                nc * sizeof(char*));
        if (!nn) return;
        s->names = nn;
        s->cap = nc;
    }
    s->names[s->count++] = name;
}
static bool split_symset_intersects(const SplitSymSet* a,
                                    const SplitSymSet* b) {
    for (size_t i = 0; i < a->count; i++) {
        if (split_symset_contains(b, a->names[i])) return true;
    }
    return false;
}

/* Walk into args only (skip the head). Symbol leaves that aren't
 * numeric constants are added to `out`. */
static void split_collect_addend_symbols(const Expr* e, SplitSymSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol)) return;
        split_symset_add(out, e->data.symbol);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        split_collect_addend_symbols(e->data.function.args[i], out);
    }
}

/* Tiny union-find (path-compressed). */
static int split_uf_find(int* parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}
static void split_uf_union(int* parent, int x, int y) {
    int rx = split_uf_find(parent, x);
    int ry = split_uf_find(parent, y);
    if (rx != ry) parent[rx] = ry;
}

/* Forward declaration: the splitter recurses through simp_dispatch. */
Expr* simp_dispatch(const Expr* input, const AssumeCtx* ctx,
                           const Expr* complexity_func);

/* Returns NULL when the input doesn't decompose, has fewer than 4
 * addends, or every component is a singleton (no win). On a successful
 * split, returns a freshly allocated, evaluated sum of per-component
 * simplifications. */
static Expr* simp_split_additive(const Expr* input, const AssumeCtx* ctx,
                                 const Expr* complexity_func) {
    if (!input || input->type != EXPR_FUNCTION) return NULL;
    if (!input->data.function.head ||
        input->data.function.head->type != EXPR_SYMBOL ||
        input->data.function.head->data.symbol != SYM_Plus)
        return NULL;
    size_t n = input->data.function.arg_count;
    if (n < 4) return NULL;

    SplitSymSet* sets = (SplitSymSet*)calloc(n, sizeof(SplitSymSet));
    int* parent = (int*)malloc(n * sizeof(int));
    if (!sets || !parent) {
        free(sets); free(parent);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        split_symset_init(&sets[i]);
        split_collect_addend_symbols(input->data.function.args[i], &sets[i]);
        parent[i] = (int)i;
    }
    /* Pairwise union by symbol intersection. */
    for (size_t i = 0; i < n; i++) {
        if (sets[i].count == 0) continue;
        for (size_t j = i + 1; j < n; j++) {
            if (sets[j].count == 0) continue;
            if (split_symset_intersects(&sets[i], &sets[j])) {
                split_uf_union(parent, (int)i, (int)j);
            }
        }
    }
    /* Glom each constant addend (empty symset) into the first variable
     * component so it gets simplified alongside variable terms. Without
     * this, splitting `1 + Plus[..variable terms..]` into `{1}` and the
     * rest can mask outer/inner constant cancellations: simplifying the
     * variable component independently may introduce a new constant
     * (e.g. Factor turning -1 - x^3/6 - x^4/6 into -1/6*(6 + x^3 + x^4))
     * that no longer cancels with the split-out `1`, leaving Simplify
     * non-idempotent. The glom is into a single non-empty component, not
     * duplicated, so the variable structure of the split is preserved. */
    {
        int first_non_empty = -1;
        for (size_t i = 0; i < n; i++) {
            if (sets[i].count != 0) { first_non_empty = (int)i; break; }
        }
        if (first_non_empty >= 0) {
            for (size_t i = 0; i < n; i++) {
                if (sets[i].count == 0) {
                    split_uf_union(parent, (int)i, first_non_empty);
                }
            }
        }
    }
    int* root_per = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_per[i] = split_uf_find(parent, (int)i);
    int* size_by_root = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) size_by_root[root_per[i]]++;

    int comp_count = 0, max_size = 0;
    for (size_t i = 0; i < n; i++) {
        if (size_by_root[i] == 0) continue;
        comp_count++;
        if (size_by_root[i] > max_size) max_size = size_by_root[i];
    }
    /* Useful split: 2+ components AND at least one component has 2+
     * addends (otherwise we'd just be re-dispatching every atom for no
     * win, and bottomup has already simplified each addend in
     * isolation). */
    if (comp_count < 2 || max_size < 2) {
        for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
        free(sets); free(parent); free(root_per); free(size_by_root);
        return NULL;
    }

    /* Compact root -> 0..comp_count-1 index. */
    int* root_to_idx = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_to_idx[i] = -1;
    int next_idx = 0;
    for (size_t i = 0; i < n; i++) {
        int r = root_per[i];
        if (root_to_idx[r] == -1) root_to_idx[r] = next_idx++;
    }
    /* Bucket addend pointers per component. */
    Expr*** comp_addends = (Expr***)calloc(comp_count, sizeof(Expr**));
    int* comp_alloc = (int*)calloc(comp_count, sizeof(int));
    int* comp_n = (int*)calloc(comp_count, sizeof(int));
    for (size_t i = 0; i < n; i++) {
        int idx = root_to_idx[root_per[i]];
        if (comp_n[idx] == comp_alloc[idx]) {
            int nc = comp_alloc[idx] == 0 ? 4 : comp_alloc[idx] * 2;
            comp_addends[idx] = (Expr**)realloc(comp_addends[idx],
                                                nc * sizeof(Expr*));
            comp_alloc[idx] = nc;
        }
        comp_addends[idx][comp_n[idx]++] =
            expr_copy(input->data.function.args[i]);
    }

    /* Simplify each component, place results into a Plus. */
    Expr** results = (Expr**)calloc(comp_count, sizeof(Expr*));
    for (int c = 0; c < comp_count; c++) {
        Expr* sub;
        if (comp_n[c] == 1) {
            sub = comp_addends[c][0];
        } else {
            Expr* p = expr_new_function(expr_new_symbol("Plus"),
                                        comp_addends[c], comp_n[c]);
            sub = evaluate(p);
            expr_free(p);
        }
        results[c] = simp_dispatch(sub, ctx, complexity_func);
        expr_free(sub);
    }
    Expr* sum_raw = expr_new_function(expr_new_symbol("Plus"),
                                      results, comp_count);
    Expr* sum = evaluate(sum_raw);
    expr_free(sum_raw);

    for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
    free(sets); free(parent); free(root_per); free(size_by_root);
    free(root_to_idx);
    for (int c = 0; c < comp_count; c++) free(comp_addends[c]);
    free(comp_addends); free(comp_alloc); free(comp_n);
    free(results);

    return sum;
}

/* simp_split_multiplicative: the multiplicative analog of
 * simp_split_additive. Decomposes a Times node whose factors fall into
 * 2+ variable-disjoint components (with at least one component holding
 * 2+ factors). Each component's sub-Times is dispatched in isolation,
 * then the per-component results are multiplied.
 *
 * Same correctness argument as the additive splitter: every transform
 * in SIMP_TRANSFORMS (Together, Cancel, Expand, Factor, TrigExpand,
 * TrigToExp, TrigFactor, TrigRoundtrip, Collect, Pythag*, half-angle,
 * Radicals) acts within a single variable's algebraic/trigonometric
 * structure -- variable-disjoint factors cannot interact under any of
 * them, so component-wise simplification preserves the answer.
 *
 * Without this pass, a stray independent factor inflates simp_search's
 * effective free-symbol budget and the heuristic gives up. With it,
 *   Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y])
 *     -> Tan[z] * simp_dispatch[Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y])]
 *     -> Tan[z] Tan[x+y]
 */
static Expr* simp_split_multiplicative(const Expr* input,
                                       const AssumeCtx* ctx,
                                       const Expr* complexity_func) {
    if (!input || input->type != EXPR_FUNCTION) return NULL;
    if (!input->data.function.head ||
        input->data.function.head->type != EXPR_SYMBOL ||
        input->data.function.head->data.symbol != SYM_Times)
        return NULL;
    size_t n = input->data.function.arg_count;
    if (n < 3) return NULL;

    SplitSymSet* sets = (SplitSymSet*)calloc(n, sizeof(SplitSymSet));
    int* parent = (int*)malloc(n * sizeof(int));
    if (!sets || !parent) {
        free(sets); free(parent);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        split_symset_init(&sets[i]);
        split_collect_addend_symbols(input->data.function.args[i], &sets[i]);
        parent[i] = (int)i;
    }
    /* Pairwise union by symbol intersection. Empty sets (constants) stay
     * in their own singleton component -- no point dragging them in. */
    for (size_t i = 0; i < n; i++) {
        if (sets[i].count == 0) continue;
        for (size_t j = i + 1; j < n; j++) {
            if (sets[j].count == 0) continue;
            if (split_symset_intersects(&sets[i], &sets[j])) {
                split_uf_union(parent, (int)i, (int)j);
            }
        }
    }
    int* root_per = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_per[i] = split_uf_find(parent, (int)i);
    int* size_by_root = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) size_by_root[root_per[i]]++;

    int comp_count = 0, max_size = 0;
    for (size_t i = 0; i < n; i++) {
        if (size_by_root[i] == 0) continue;
        comp_count++;
        if (size_by_root[i] > max_size) max_size = size_by_root[i];
    }
    /* Useful split: 2+ components AND at least one component holds 2+
     * factors (singleton-only splits would just re-dispatch each atom
     * for no win). */
    if (comp_count < 2 || max_size < 2) {
        for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
        free(sets); free(parent); free(root_per); free(size_by_root);
        return NULL;
    }

    int* root_to_idx = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_to_idx[i] = -1;
    int next_idx = 0;
    for (size_t i = 0; i < n; i++) {
        int r = root_per[i];
        if (root_to_idx[r] == -1) root_to_idx[r] = next_idx++;
    }
    Expr*** comp_factors = (Expr***)calloc(comp_count, sizeof(Expr**));
    int* comp_alloc = (int*)calloc(comp_count, sizeof(int));
    int* comp_n = (int*)calloc(comp_count, sizeof(int));
    for (size_t i = 0; i < n; i++) {
        int idx = root_to_idx[root_per[i]];
        if (comp_n[idx] == comp_alloc[idx]) {
            int nc = comp_alloc[idx] == 0 ? 4 : comp_alloc[idx] * 2;
            comp_factors[idx] = (Expr**)realloc(comp_factors[idx],
                                                nc * sizeof(Expr*));
            comp_alloc[idx] = nc;
        }
        comp_factors[idx][comp_n[idx]++] =
            expr_copy(input->data.function.args[i]);
    }

    Expr** results = (Expr**)calloc(comp_count, sizeof(Expr*));
    for (int c = 0; c < comp_count; c++) {
        Expr* sub;
        if (comp_n[c] == 1) {
            sub = comp_factors[c][0];
        } else {
            Expr* p = expr_new_function(expr_new_symbol("Times"),
                                        comp_factors[c], comp_n[c]);
            sub = evaluate(p);
            expr_free(p);
        }
        results[c] = simp_dispatch(sub, ctx, complexity_func);
        expr_free(sub);
    }
    Expr* prod_raw = expr_new_function(expr_new_symbol("Times"),
                                       results, comp_count);
    Expr* prod = evaluate(prod_raw);
    expr_free(prod_raw);

    for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
    free(sets); free(parent); free(root_per); free(size_by_root);
    free(root_to_idx);
    for (int c = 0; c < comp_count; c++) free(comp_factors[c]);
    free(comp_factors); free(comp_alloc); free(comp_n);
    free(results);

    return prod;
}

/* simp_dispatch is the public entry point. It runs the shape classifier
 * and forwards to a specialised pipeline. SHAPE_TRIG and SHAPE_GENERAL
 * fall through to simp_search; trig is the heuristic search's strongest
 * domain, and general inputs need the full machinery. */
Expr* simp_dispatch(const Expr* input, const AssumeCtx* ctx,
                           const Expr* complexity_func) {
    /* Universal Power-identity seeds: PrimeRebase and PowerOneify are
     * cheap, inert-by-default rewrites that benefit every pipeline
     * (POLYNOMIAL, RATIONAL, LOGEXP, TRIG, GENERAL alike).  Both are
     * idempotent on their fixed point, so the strict-win recursion
     * below terminates after at most one rebase + one oneify cascade.
     * Score-gated: only adopt the rewrite when it strictly beats the
     * input; otherwise fall through to the normal classifier route. */
    {
        Expr* alt = transform_prime_rebase(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    {
        Expr* alt = transform_power_oneify(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    /* PowerDistribute: principal-branch power-distribution rewrites that
     * are universally sound -- (-c)^e split, constant-positive base
     * distribute, and integer-exponent distribute (via ctx).  Same
     * recursion shape as the rebase / oneify gates above: only adopt
     * the rewrite when it strictly beats the input on score. */
    {
        Expr* alt = transform_power_distribute(input, ctx);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    /* RadicalCanon: split Power[Rational[a,b], q] into Power[a,q] *
     * Power[b,-q] and rationalise Power[positive_int, negative_rational]
     * into Power[..., positive_residue] / a^k.  Forces equivalent radical
     * shapes (Sqrt[1/2] vs 1/Sqrt[2] vs Sqrt[2]/2) into a single canonical
     * representation so additive cancellations like
     *     -Sqrt[1/2]/3 + Sqrt[2]/6  ->  0
     * fire instead of being trapped behind shape divergence.  This is
     * NOT a strict-win gate: even when the score is the same, the rewrite
     * is correctness-preserving and always yields the more useful form
     * for downstream additive cancellation; the score check below uses
     * <= rather than <. */
    {
        Expr* alt = transform_radical_canon(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt <= s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }

    /* Try to decompose a Plus into disjoint-variable components and
     * simplify each independently. Each connected component's sub-Plus
     * cannot itself decompose further, so the recursion through
     * simp_dispatch terminates in one level. */
    Expr* split = simp_split_additive(input, ctx, complexity_func);
    if (split) return split;
    /* Same idea for Times: lift variable-disjoint factors out of the
     * search space. */
    Expr* tsplit = simp_split_multiplicative(input, ctx, complexity_func);
    if (tsplit) return tsplit;

    SimpShape shape = simp_classify(input);
    switch (shape) {
        case SIMP_SHAPE_POLYNOMIAL:
            return simp_pipeline_polynomial(input, ctx, complexity_func);
        case SIMP_SHAPE_RATIONAL:
            return simp_pipeline_rational(input, ctx, complexity_func);
        case SIMP_SHAPE_LOGEXP:
            return simp_pipeline_logexp(input, ctx, complexity_func);
        case SIMP_SHAPE_TRIG: {
            /* Fast algebraic normal-form pre-check for rational
             * functions of Sin/Cos/Sinh/Cosh (and Tan/Cot/Sec/Csc/Tanh/
             * etc. after preprocessing). Strict leaf-count gate inside:
             * returns NULL when the algorithm does not apply or does
             * not strictly improve the input, so we fall through to
             * simp_search unchanged in those cases. */
            Expr* tr = simp_trig_rational(input, ctx, complexity_func);
            if (tr) return tr;
            return simp_search(input, ctx, complexity_func);
        }
        case SIMP_SHAPE_GENERAL:
        default:
            return simp_search(input, ctx, complexity_func);
    }
}

Expr* simp_search(const Expr* original_input, const AssumeCtx* ctx,
                         const Expr* complexity_func) {
    /* Phase 0: pre-apply the Abs structural rules. These (idempotent,
     * force-take) rewrites canonicalise Abs[Times[...]] -> product of Abs,
     * Abs[Abs[x]] -> Abs[x], Abs[E^z] -> E^Re[z], etc. Doing it here
     * (rather than as a regular seed) means the rest of the search starts
     * from the rewritten form -- otherwise the original (often smaller in
     * leaf count) re-enters the candidate set and the leaf-count tiebreak
     * brings us right back. */
    Expr* abs_pre = transform_can_fire("AbsRules", original_input, ctx)
                        ? apply_abs_rules(original_input, ctx)
                        : NULL;
    const Expr* input;
    if (abs_pre && !expr_eq(abs_pre, original_input)) {
        if (simp_debug_enabled()) {
            simp_debug_log("AbsRules", original_input, abs_pre, 0.0);
        }
        input = abs_pre;
    } else {
        if (abs_pre) { expr_free(abs_pre); abs_pre = NULL; }
        input = original_input;
    }

    /* Phase 0b: pre-apply Sqrt[_^2] rewrites. Same rationale as the Abs
     * pre-apply above: turn Sqrt[e^2] into e / -e / Abs[e] before the
     * candidate search starts, so the rewritten form anchors the leaf-
     * count tiebreak. Built on top of the (already canonicalised) `input`
     * so the two phases compose. */
    Expr* sqsq_pre = transform_can_fire("SqrtSquareRules", input, ctx)
                         ? apply_sqrt_of_square_rules(input, ctx)
                         : NULL;
    if (sqsq_pre && !expr_eq(sqsq_pre, input)) {
        if (simp_debug_enabled()) {
            simp_debug_log("SqrtSquareRules", input, sqsq_pre, 0.0);
        }
        input = sqsq_pre;
    } else {
        if (sqsq_pre) { expr_free(sqsq_pre); sqsq_pre = NULL; }
    }

    Expr* best = expr_copy((Expr*)input);
    size_t best_score = score_with_func(best, complexity_func);

    CandSet seeds;
    cs_init(&seeds);
    cs_add_or_free(&seeds, expr_copy((Expr*)input));

    /* Assumption-derived alternatives. An assumption-aware rewrite that
     * actually changed the form is correctness-preserving under the
     * assumption set and is by definition more "simplified" than the
     * input, even if the leaf-count tiebreak is even (e.g.
     * Log[x^p] -> p Log[x] both score 4). Force-take it as the new best
     * so long as it isn't strictly worse. */
    if (transform_can_fire("AssumptionRules", input, ctx)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = apply_assumption_rules(input, ctx);
        if (dbg) simp_debug_log("AssumptionRules", input, alt, simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            size_t s = score_with_func(alt, complexity_func);
            if (s <= best_score) {
                expr_free(best);
                best = expr_copy(alt);
                best_score = s;
            }
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Logexp rewriter seed. The Log/Power identities (Log[a b] -> Log[a] +
     * Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), ...) typically
     * INCREASE leaf count, so they cannot win the standard complexity
     * tiebreak. We force them as the new best regardless of score: the
     * cascade is correctness-preserving under positivity assumptions and
     * the user's intent (per the documented rule cascade) is that they
     * fire whenever conditions are met. Downstream transforms (Cancel,
     * Together, ...) cannot recombine an expanded log/power, so the
     * expanded form persists through the rest of the search. */
    if (transform_can_fire("LogExpRules", input, ctx)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = apply_logexp_rules(input, ctx);
        if (dbg) simp_debug_log("LogExpRules", input, alt, simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            expr_free(best);
            best = expr_copy(alt);
            best_score = score_with_func(best, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* SimpLogRules seed. Force-take when it changes the form: Pass A
     * (prime decomposition of Log[rational]) and Pass B (linear
     * combination fuser) are correctness-preserving under positivity
     * and frequently expose dramatic cancellations on inputs with both
     * trig and log heads (e.g. Log[Sin] - Log[Tan] - Log[Cos] in
     * SIMP_SHAPE_TRIG inputs that bypass simp_pipeline_logexp). */
    if (contains_log(input)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = simp_log_apply(input, ctx);
        if (dbg) simp_debug_log("SimpLogRules", input, alt,
                                simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            expr_free(best);
            best = expr_copy(alt);
            best_score = score_with_func(best, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Roundtrip seed. Score-gate the seed propagation: TrigRoundtrip
     * runs Together inside, which on pure-real Sinh/Cosh exp forms can
     * factor out an asymmetric E^(-kx) that introduces a high-frequency
     * E^(2kx) term in the cofactor (e.g.
     * Together[-1/2 + 1/4 E^(-2x) + 1/4 E^(2x)]
     *   ->  1/4 E^(-2x) (1 - 2 E^(2x) + E^(4x))).
     * ExpToTrig then turns that into a 16-term Cosh/Sinh product that
     * subsequent transforms (Together, Cancel, Factor) all consume
     * 100s of ms on. update_best still picks the result if it's a win;
     * we only refuse to propagate dramatic blow-ups as a seed for
     * further exploration. */
    if (transform_can_fire("TrigRoundtrip", input, NULL)) {
        Expr* alt = transform_trig_roundtrip(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            size_t alt_score = score_with_func(alt, complexity_func);
            size_t input_score = score_with_func(input, complexity_func);
            if (alt_score <= 2 * input_score + 8) {
                cs_add_or_free(&seeds, alt);
            } else {
                expr_free(alt);
            }
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* ExpToTrig seed. TrigRoundtrip is gated on contains_trig_or_hyperbolic
     * which rejects pure-Exp inputs (e.g.
     * `-1 + (-E^-x + E^x)^2/(E^-x + E^x)^2` has no Cosh/Sinh head and so
     * misses the entire trig pipeline). ExpToTrig converts E^... into
     * Cosh/Sinh form, after which PythagReduce / Together / Cancel can
     * collapse the result (the example above lands at `-Sech[x]^2`).
     * Score-gate the seed propagation the same way TrigRoundtrip does:
     * keep wins as the new best, but only forward seeds that haven't
     * blown up the structure. Skip when the input already has trig or
     * hyperbolic heads -- TrigRoundtrip handles that case. */
    if (contains_exp_form(input) && !contains_trig_or_hyperbolic(input)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = call_unary_copy("ExpToTrig", input);
        if (dbg) simp_debug_log("ExpToTrigSeed", input, alt,
                                simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            size_t alt_score = score_with_func(alt, complexity_func);
            size_t input_score = score_with_func(input, complexity_func);
            if (alt_score <= 2 * input_score + 8) {
                cs_add_or_free(&seeds, alt);
            } else {
                expr_free(alt);
            }
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean square-completion seed. Idempotent on inputs without the
     * 1 +/- 2 Sin Cos shape, so always cheap to try. */
    {
        Expr* alt = transform_pythag_square_complete(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean reduction seed (1 - Cos^2 -> Sin^2, etc.). Strict
     * leaf-count win when it fires; inert otherwise. */
    {
        Expr* alt = transform_pythag_reduce(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean canonicalizer seed: substitute Cos^(2k) -> (1 - Sin^2)^k
     * (or the reverse) globally and Expand. Catches difference-of-squares
     * trig products that, after Expand, leave coefficient-bearing Cos^2
     * factors PythagReduce cannot match. The transform itself only
     * returns a strict-score win, so passing it through update_best is
     * safe: it never propagates a worse seed. */
    {
        Expr* alt = transform_pythag_canon(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* TrigReduce seed.  The product-to-sum / power-reduction pass
     * collapses inputs like `(Sin[x]+Cos[x])^4 - (1+Sin[2x])^2` directly
     * to 0 (after internal Expand + angle-addition collapse) in a few
     * ms.  The same call appears later in the round loop's
     * SIMP_TRANSFORMS table at index 11, which means without this seed
     * the search runs Factor (~28 ms), TrigFactor (~14 ms),
     * FactorSquareFree (~9 ms), etc. on the input *before* TrigReduce
     * gets a turn, even though TrigReduce alone is enough.  Promoting
     * it to the seed phase + relying on simp_best_is_zero to short-
     * circuit afterwards turns ~220 ms into ~5 ms on this shape.  Inert
     * (returns the input unchanged via the trig_memo no-op cycle) when
     * the input has no trig-product structure to reduce. */
    if (transform_can_fire("TrigReduce", input, NULL)) {
        Expr* alt = traced_call_unary("TrigReduce", input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                update_best(&best, &best_score, alt, complexity_func);
                size_t alt_score = score_with_func(alt, complexity_func);
                size_t input_score = score_with_func(input, complexity_func);
                /* Same blow-up guard TrigRoundtrip uses: keep the wins,
                 * drop catastrophically larger candidates from the seed
                 * set so they don't seed a wave of expensive Factor/
                 * TrigFactor work in subsequent rounds. */
                if (alt_score <= 2 * input_score + 8) {
                    cs_add_or_free(&seeds, alt);
                } else {
                    expr_free(alt);
                }
            } else {
                expr_free(alt);
            }
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* Trig-at-rational-Pi canonicalization seed.  Picks a unique
     * Sin-vs-Cos / Tan-vs-Cot / Sec-vs-Csc representation per pair,
     * which lets pairs like Cos[4 Pi/9] - Sin[Pi/18] (both lower to
     * Sin[Pi/18]) and 1/8 (Cos[4 Pi/9] + Cos[5 Pi/9]) (lowers to
     * 1/8 (Sin[Pi/18] + (-Sin[Pi/18])) = 0) cancel additively in
     * the surrounding Plus.  Idempotent and inert when the input
     * has no rational-Pi trig calls. */
    {
        Expr* alt = simp_trig_pi_canon(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* TanAddition seed: search for triples (a, b, c=a+b) among the trig
     * args appearing in the input and rewrite Tan[c]/Cot[c]/etc. via
     * the angle-addition formula in terms of Tan[a], Tan[b].  Catches
     * shapes like Tan[2] Tan[3] - 1 + (Tan[2]+Tan[3])/Tan[5] -> 0 (since
     * 5 = 2+3) which would otherwise go un-recognised because TrigExpand
     * only fires on literal-Plus arguments.  Inert when fewer than three
     * distinct trig args appear or when no pair sums to a third arg. */
    {
        Expr* alt = transform_tan_addition(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* Half-angle tangent / Tanh seed. Idempotent on inputs without
     * the Sin/(1+Cos) (resp. Sinh/(1+Cosh)) shape. */
    {
        Expr* alt = transform_halfangle(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Radical product seed: collapses Sqrt[a]*Sqrt[b] -> Sqrt[a*b] for
     * positive integer a, b (and similarly for higher rational
     * exponents).  Inert on inputs without Power[+integer, Rational]
     * factors. */
    {
        Expr* alt = simp_radicals(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* The common-factor lift is applied as a final-form polish after
     * simp_bottomup completes (see builtin_simplify), not as a search
     * seed. Wiring it into the round loop or the seed phase changes the
     * structural shape of intermediate forms in ways that interact badly
     * with TrigExpand / TrigReduce on multi-variable trig inputs. */

    /* Sqrt-of-Sqrt denesting seed (phase-1 radical denesting). Fires the
     * half-sum identity Sqrt[A + b Sqrt[C]] -> Sqrt[(A+s)/2] +/- Sqrt[(A-s)/2]
     * where s^2 = A^2 - b^2 C is detected via integer/rational
     * perfect-square or polynomial FactorSquareFree. Inert (returns the
     * input copy) when no Sqrt[Plus[...]] subtree denests cleanly under
     * the active assumption set. See the simp_denest_sqrt section above
     * for branch-soundness preconditions. */
    {
        Expr* alt = simp_denest_sqrt(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Denominator-rationalisation seed (phase-2 radical simplification).
     * Walks the tree looking for Power[denom, -1] where denom is a
     * polynomial in radicals over a single positive integer base c,
     * and computes the inverse via extended-Euclidean in
     * Q[α]/(α^n - c). Closes user cases like 1/(2^(1/3) - 1) ->
     * 4^(1/3) + 2^(1/3) + 1. Inert when no such Power[denom, -1]
     * subtree exists or when denom involves multiple distinct integer
     * bases (out of scope for the current single-extension
     * implementation). */
    {
        Expr* alt = simp_rationalize_denom(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Cube-root denesting / sum-of-conjugates seed (phase-3). Pattern A:
     * (a + b sqrt(c))^(1/3) -> p + q sqrt(c) via small-grid search.
     * Pattern B: (a+b sqrt(c))^(1/3) + (a-b sqrt(c))^(1/3) -> rational
     * via Cardano discriminant + rational-root test. Inert otherwise. */
    {
        Expr* alt = simp_cuberoot(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Algebraic-extension seed: substitute each Sqrt[u_i] by a fresh
     * generator g_i, reduce in K(vars)[g_1,...,g_n]/(g_i^2 - u_i), and
     * rationalise the denominator by successive sigma-conjugation.
     * Collapses identities ordinary Together / Cancel cannot see, e.g.
     *   (x/Sqrt[x^2+1] + 1) / ((Sqrt[x^2+1] + x)^2 + 1)  ->  1/(2 + 2 x^2)
     * Inert when there is no surd, when any surd has q != 2, or when
     * the input contains an explicit complex literal. */
    {
        Expr* alt = simp_algebraic(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Roots-of-unity seed. Reduces sums of (-1)^(p/q) and E^(I p Pi/q)
     * via the minimal polynomial Phi_{2Q}(x); see
     * simp_roots_of_unity above. Idempotent and inert when the input
     * has no root-of-unity atoms. */
    {
        Expr* alt = simp_roots_of_unity(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Factorial seed. Decomposes every Factorial[arg] as arg = sym + c,
     * groups by sym, and rewrites each call into Factorial[base] times
     * an explicit Pochhammer product so common Factorial[base] factors
     * cancel/absorb under Together / Cancel / Expand. A re-fold pass
     * then absorbs (base+1)(base+2)...(base+k) cofactors back into
     * Factorial[base+k] (with up to two missing-factor gap-fills). The
     * separate (2v)!/(2^v v!) -> Factorial2[2v-1] check piggybacks on
     * the same walk. Inert when the input has no Factorial atoms.
     *
     * Force-take semantics (mirroring LogExpRules / AssumptionRules): a
     * shape-changing factorial rewrite is by definition "more simplified"
     * even when the SimplifyCount tiebreak ties. The default complexity
     * measure cannot distinguish (n-2)!/n! (count 11) from 1/(n^2-n)
     * (count 11) but the latter is unambiguously the canonical form. We
     * force-take whenever the number of Factorial atoms strictly drops
     * (or stays the same and the leaf count is no worse). */
    if (transform_can_fire("FactorialRules", input, NULL)) {
        Expr* alt = simp_factorial(input);
        if (alt && !expr_eq(alt, input)) {
            size_t s = score_with_func(alt, complexity_func);
            /* Count Factorial atoms in input and alt. */
            size_t fc_in = 0, fc_alt = 0;
            {
                /* Inline tiny tree walk -- avoid threading a counter into
                 * the contains_factorial signature. */
                CandSet stack;
                cs_init(&stack);
                cs_add_or_free(&stack, expr_copy((Expr*)input));
                while (stack.count) {
                    Expr* top = stack.items[stack.count - 1];
                    stack.count--;
                    if (top && top->type == EXPR_FUNCTION) {
                        if (simp_eq_head_sym(top, "Factorial")) fc_in++;
                        cs_add_or_free(&stack, expr_copy(top->data.function.head));
                        for (size_t i = 0; i < top->data.function.arg_count; i++) {
                            cs_add_or_free(&stack,
                                expr_copy(top->data.function.args[i]));
                        }
                    }
                    expr_free(top);
                }
                cs_free(&stack);
            }
            {
                CandSet stack;
                cs_init(&stack);
                cs_add_or_free(&stack, expr_copy(alt));
                while (stack.count) {
                    Expr* top = stack.items[stack.count - 1];
                    stack.count--;
                    if (top && top->type == EXPR_FUNCTION) {
                        if (simp_eq_head_sym(top, "Factorial")) fc_alt++;
                        cs_add_or_free(&stack, expr_copy(top->data.function.head));
                        for (size_t i = 0; i < top->data.function.arg_count; i++) {
                            cs_add_or_free(&stack,
                                expr_copy(top->data.function.args[i]));
                        }
                    }
                    expr_free(top);
                }
                cs_free(&stack);
            }
            bool force_take = (fc_alt < fc_in) ||
                              (fc_alt == fc_in && s <= best_score);
            if (force_take) {
                expr_free(best);
                best = expr_copy(alt);
                best_score = s;
            } else {
                update_best(&best, &best_score, alt, complexity_func);
            }
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Per-variable Collect seed. parent_score = score(input). */
    if (transform_can_fire("CollectPerVariable", input, NULL)) {
        try_collect_per_variable(input, best_score, &seeds, &best, &best_score,
                                 complexity_func);
    }

    /*
     * Round loop. Branch-pruning rule: a candidate produced by a
     * transform on `seed` is propagated to `next` only when its
     * complexity is no greater than `seed`'s. Strictly worse forms
     * cannot lead to a better best through any further unary transform
     * we apply (they'd just feed transforms more work and grow the
     * candidate set), so we drop them. They may still beat the running
     * best on this very step (update_best already handles that), but
     * they won't seed further exploration.
     *
     * Assumption-aware rewrites and the logexp cascade keep their
     * "force-win" behaviour for the best slot (they're correctness-
     * preserving under the assumption set and qualitatively "more
     * simplified"), but the same parent-score gate applies to whether
     * they propagate as a new seed.
     */
    for (int round = 0; round < SIMP_ROUNDS; round++) {
        CandSet next;
        cs_init(&next);
        for (size_t i = 0; i < seeds.count; i++) {
            const Expr* seed = seeds.items[i];
            size_t parent_score = score_with_func(seed, complexity_func);

            for (size_t t = 0; t < SIMP_TRANSFORM_COUNT; t++) {
                if (!transform_can_fire(SIMP_TRANSFORMS[t], seed, ctx)) continue;
                Expr* r = traced_call_unary(SIMP_TRANSFORMS[t], seed);
                if (!r) continue;
                update_best(&best, &best_score, r, complexity_func);
                /* Chain a PythagReduce pass on the transform output so a
                 * candidate produced in the final round (e.g. FactorSquareFree
                 * surfacing Cos[x]^2 - 1 inside a product) still gets the
                 * Pythagorean rewrite applied -- otherwise the reduction
                 * would only fire when this candidate became a seed in a
                 * subsequent round, which doesn't happen at SIMP_ROUNDS-1. */
                {
                    Expr* pr = transform_pythag_reduce(r);
                    if (pr) {
                        if (!expr_eq(pr, r)) {
                            update_best(&best, &best_score, pr, complexity_func);
                        }
                        expr_free(pr);
                    }
                }
                /* Chain a simp_radicals pass on the transform output.
                 * Together / Cancel after TrigExpand can surface
                 * intermediate forms with adjacent radical factors
                 * (Sqrt[a] Sqrt[b]) -- e.g.
                 * Together[4/3 Sin/Sqrt[2] - 2/9 Sqrt[3] Sqrt[6] Sin + ...]
                 *   -> (-12 Sqrt[2] Sqrt[3] Sin + 12 Sqrt[6] Sin)/(9 Sqrt[2] Sqrt[6])
                 * whose combination collapses the numerator to 0. Without
                 * this chain, simp_radicals would only fire when r becomes
                 * a seed next round -- which doesn't happen at the final
                 * round. */
                {
                    Expr* rd = simp_radicals(r);
                    if (rd) {
                        if (!expr_eq(rd, r)) {
                            update_best(&best, &best_score, rd, complexity_func);
                        }
                        expr_free(rd);
                    }
                }
                /* Chain a trig-Pi canonicalization pass on the transform
                 * output.  TrigReduce on Cos[Pi/9] Cos[2Pi/9] Cos[3Pi/9]
                 * Cos[4Pi/9] - 1/16 surfaces the intermediate
                 * 1/8 (Cos[4 Pi/9] + Cos[5 Pi/9]); only after canonicalising
                 * Cos[5 Pi/9] -> -Sin[Pi/18] (equivalently -Cos[4 Pi/9])
                 * does the Plus collapse to 0.  Without this chain the
                 * canonicalizer would have to wait for r to become a seed
                 * next round, which does not happen at SIMP_ROUNDS-1. */
                {
                    Expr* tc = simp_trig_pi_canon(r);
                    if (tc) {
                        if (!expr_eq(tc, r)) {
                            update_best(&best, &best_score, tc, complexity_func);
                        }
                        expr_free(tc);
                    }
                }
                if (expr_eq(r, seed)) {
                    expr_free(r);
                } else if (score_with_func(r, complexity_func) > parent_score) {
                    /* TrigExpand expands Sin[a+b]/Cos[a+b] into
                     * Sin[a]Cos[b]+Cos[a]Sin[b] etc., which usually grows
                     * the leaf count but surfaces radical products
                     * (Sqrt[a]Sqrt[b]) and rationalisable coefficients that
                     * subsequent transforms (Together, simp_radicals,
                     * Collect) can collapse -- sometimes all the way to 0.
                     * The default score gate would drop those candidates;
                     * loosen it for TrigExpand using the same bound
                     * TrigRoundtrip uses at the seed phase (2*input + 8) so
                     * pathological blow-ups (Sin[a+b+c+d]) still get
                     * pruned. */
                    if (strcmp(SIMP_TRANSFORMS[t], "TrigExpand") == 0 &&
                        score_with_func(r, complexity_func) <=
                            2 * parent_score + 8) {
                        cs_add_or_free(&next, r);
                    } else {
                        expr_free(r);
                    }
                } else {
                    cs_add_or_free(&next, r);
                }
            }
            /* Also try the trig roundtrip on each candidate. */
            if (transform_can_fire("TrigRoundtrip", seed, NULL)) {
                Expr* tr = transform_trig_roundtrip(seed);
                if (tr) {
                    update_best(&best, &best_score, tr, complexity_func);
                    if (expr_eq(tr, seed)) {
                        expr_free(tr);
                    } else if (score_with_func(tr, complexity_func) > parent_score) {
                        expr_free(tr);
                    } else {
                        cs_add_or_free(&next, tr);
                    }
                }
            }
            /* Pythagorean square completion on each candidate. The Factor
             * transform run earlier may have produced (1 + 2 Sin Cos)^2;
             * this round step lets the rule fire on that intermediate. */
            {
                Expr* psc = transform_pythag_square_complete(seed);
                if (psc) {
                    update_best(&best, &best_score, psc, complexity_func);
                    if (expr_eq(psc, seed)) {
                        expr_free(psc);
                    } else if (score_with_func(psc, complexity_func) > parent_score) {
                        expr_free(psc);
                    } else {
                        cs_add_or_free(&next, psc);
                    }
                }
            }
            /* Pythagorean reduction on each candidate (1 - Cos^2 -> Sin^2,
             * Cosh^2 - 1 -> Sinh^2, etc.). Strict leaf-count win when
             * matched. */
            {
                Expr* pr = transform_pythag_reduce(seed);
                if (pr) {
                    update_best(&best, &best_score, pr, complexity_func);
                    if (expr_eq(pr, seed)) {
                        expr_free(pr);
                    } else if (score_with_func(pr, complexity_func) > parent_score) {
                        expr_free(pr);
                    } else {
                        cs_add_or_free(&next, pr);
                    }
                }
            }
            /* Half-angle tangent / Tanh on each candidate. Lets Together
             * /Cancel intermediates (which can surface (1+Cos[x]) or
             * (1+Cosh[x]) factors after partial cancellation) feed into
             * the rule. */
            {
                Expr* ha = transform_halfangle(seed);
                if (ha) {
                    update_best(&best, &best_score, ha, complexity_func);
                    if (expr_eq(ha, seed)) {
                        expr_free(ha);
                    } else if (score_with_func(ha, complexity_func) > parent_score) {
                        expr_free(ha);
                    } else {
                        cs_add_or_free(&next, ha);
                    }
                }
            }
            /* Radical product combine on each candidate. Together /
             * Cancel can surface fresh Sqrt[a]*Sqrt[b] products in their
             * output; this lets the combine fire on the intermediate. */
            {
                Expr* rd = simp_radicals(seed);
                if (rd) {
                    update_best(&best, &best_score, rd, complexity_func);
                    if (expr_eq(rd, seed)) {
                        expr_free(rd);
                    } else if (score_with_func(rd, complexity_func) > parent_score) {
                        expr_free(rd);
                    } else {
                        cs_add_or_free(&next, rd);
                    }
                }
            }
            /* Common-factor lift not applied per-candidate -- see comment
             * in seed phase above and the final-form polish in
             * builtin_simplify. */
            /* Abs rewrites on each candidate. Same force-take semantics
             * as the seed phase (idempotent rules, structural answer). */
            if (transform_can_fire("AbsRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t0 = dbg ? clock() : 0;
                Expr* abr = apply_abs_rules(seed, ctx);
                if (dbg) simp_debug_log("AbsRules", seed, abr, simp_debug_elapsed_ms(t0));
                if (abr) {
                    if (!expr_eq(abr, seed)) {
                        size_t s = score_with_func(abr, complexity_func);
                        expr_free(best);
                        best = expr_copy(abr);
                        best_score = s;
                        if (s <= parent_score) {
                            cs_add_or_free(&next, abr);
                        } else {
                            expr_free(abr);
                        }
                    } else {
                        expr_free(abr);
                    }
                }
            }
            /* Sqrt[_^2] rewrites on each candidate. Same force-take. New
             * Sqrt[X^2] subexpressions can surface mid-search via Expand
             * or Together, so re-running per candidate (not just at seed
             * time) catches the late-arriving ones. */
            if (transform_can_fire("SqrtSquareRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t0 = dbg ? clock() : 0;
                Expr* sqr = apply_sqrt_of_square_rules(seed, ctx);
                if (dbg) simp_debug_log("SqrtSquareRules", seed, sqr, simp_debug_elapsed_ms(t0));
                if (sqr) {
                    if (!expr_eq(sqr, seed)) {
                        size_t s = score_with_func(sqr, complexity_func);
                        if (s < best_score) {
                            expr_free(best);
                            best = expr_copy(sqr);
                            best_score = s;
                        }
                        if (s <= parent_score) {
                            cs_add_or_free(&next, sqr);
                        } else {
                            expr_free(sqr);
                        }
                    } else {
                        expr_free(sqr);
                    }
                }
            }
            /* And per-variable Collect on each candidate. */
            if (transform_can_fire("CollectPerVariable", seed, NULL)) {
                try_collect_per_variable(seed, parent_score, &next, &best, &best_score,
                                         complexity_func);
            }
            /* And the assumption rewriter on each candidate. Bias as in
             * the seed phase: prefer assumption-driven forms at equal
             * complexity for `best`; gate seeding on parent_score. */
            if (transform_can_fire("AssumptionRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t0 = dbg ? clock() : 0;
                Expr* ar = apply_assumption_rules(seed, ctx);
                if (dbg) simp_debug_log("AssumptionRules", seed, ar, simp_debug_elapsed_ms(t0));
                if (ar) {
                    if (!expr_eq(ar, seed)) {
                        size_t s = score_with_func(ar, complexity_func);
                        if (s <= best_score) {
                            expr_free(best);
                            best = expr_copy(ar);
                            best_score = s;
                        }
                        if (s <= parent_score) {
                            cs_add_or_free(&next, ar);
                        } else {
                            expr_free(ar);
                        }
                    } else {
                        expr_free(ar);
                    }
                }
            }
            if (transform_can_fire("LogExpRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t1 = dbg ? clock() : 0;
                Expr* lr = apply_logexp_rules(seed, ctx);
                if (dbg) simp_debug_log("LogExpRules", seed, lr, simp_debug_elapsed_ms(t1));
                if (lr) {
                    if (!expr_eq(lr, seed)) {
                        /* Force-win for `best` (correctness-preserving
                         * under positivity assumptions). Still gate the
                         * seeding step. */
                        size_t s = score_with_func(lr, complexity_func);
                        expr_free(best);
                        best = expr_copy(lr);
                        best_score = s;
                        if (s <= parent_score) {
                            cs_add_or_free(&next, lr);
                        } else {
                            expr_free(lr);
                        }
                    } else {
                        expr_free(lr);
                    }
                }
            }
        }
        cs_free(&seeds);
        seeds = next;
        if (seeds.count == 0) break;
        if (simp_best_is_zero(best)) break;
    }
search_done:
    cs_free(&seeds);
    if (abs_pre) expr_free(abs_pre);
    if (sqsq_pre) expr_free(sqsq_pre);
    return best;
}


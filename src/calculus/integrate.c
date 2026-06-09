/* integrate.c
 *
 * `Integrate[f, x]` System` dispatcher.  Cascades through three
 * stages and supports an explicit `Method -> "..."` option:
 *
 *   1. Integrate`BronsteinRational   — polynomial / rational integrands
 *   2. Integrate`RischNorman         — parallel-Risch (Bronstein pmint)
 *   3. Integrate`CRCTable            — CRC integral table (lazy-loaded)
 *
 * Method values: "Automatic" (default, full cascade), "BronsteinRational",
 * "RischNorman", "CRCTable" (strict passthrough, no fallback).
 *
 * The CRC table is large and most sessions never need it, so its
 * .m file is Get-loaded on first invocation of try_crctable() rather
 * than at startup.  See LAZY_LOAD_CRC below.
 */

#include "integrate.h"
#include "integrate_interp.h"
#include "integrate_unknown.h"
#include "integrate_derivdivides.h"
#include "integrate_linrad.h"
#include "integrate_quadrad.h"
#include "integrate_linratiorad.h"
#include "integrate_jeffrey.h"
#include "intrat.h"
#include "intrischnorman.h"
#include "common.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* True iff `e` is the symbol `True`.  The PolynomialQ / rationalQ
 * predicates we call below return either True or False. */
static bool is_true_symbol(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* Test whether `f` is a polynomial in `x`.  Calls the existing
 * PolynomialQ builtin so we get the same definition the rest of
 * Mathilda uses. */
static bool is_polynomial_in(Expr* f, Expr* x) {
    Expr* args[2] = { expr_copy(f), expr_copy(x) };
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = is_true_symbol(val);
    expr_free(val);
    return ok;
}

/* Test whether `f` is a rational function in `x`: Together[f] must
 * have a non-trivial denominator that is itself polynomial in `x`,
 * with a polynomial numerator. */
static bool is_rational_in(Expr* f, Expr* x) {
    Expr* together = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(together);
    expr_free(together);
    if (!combined) return false;

    Expr* num = internal_numerator((Expr*[]){expr_copy(combined)}, 1);
    Expr* den = internal_denominator((Expr*[]){expr_copy(combined)}, 1);
    expr_free(combined);
    Expr* num_v = evaluate(num);
    Expr* den_v = evaluate(den);
    expr_free(num); expr_free(den);

    bool ok = false;
    if (num_v && den_v) {
        bool num_is_poly = is_polynomial_in(num_v, x);
        bool den_is_poly = is_polynomial_in(den_v, x);
        bool den_is_unit = (den_v->type == EXPR_INTEGER && den_v->data.integer == 1);
        ok = num_is_poly && den_is_poly && !den_is_unit;
    }
    if (num_v) expr_free(num_v);
    if (den_v) expr_free(den_v);
    return ok;
}

/* True iff `result` is the unevaluated call `head_name[...]`. */
static bool result_is_unresolved(const Expr* result, const char* head_name) {
    return head_is(result, intern_symbol(head_name));
}

/* True iff any subexpression of `e` is an unevaluated call `head_name[...]`.
 * Needed for the CRC pipeline: the inner table rules can fire partially
 * (e.g. extracting a `-1` via `IntegrateTable[a_ f_, x_] := a IntegrateTable[f, x]`)
 * and leave a stray `IntegrateTable[...]` nested inside a Times / Plus.
 * The top-level head_is check would miss those, so we walk the tree. */
static bool result_contains_head(const Expr* e, const char* head_name) {
    if (!e) return false;
    if (head_is(e, intern_symbol(head_name))) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (result_contains_head(e->data.function.head, head_name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (result_contains_head(e->data.function.args[i], head_name)) return true;
    }
    return false;
}

/* Helper: build and evaluate `head_name[f, x]`, freeing the call
 * expression.  Returns whatever evaluate() produces. */
static Expr* call_stage(const char* head_name, Expr* f, Expr* x) {
    Expr* call = expr_new_function(
        expr_new_symbol(head_name),
        (Expr*[]){ expr_copy(f), expr_copy(x) }, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

/* Stage 0: undefined-function integrator (Roach 1992, §1.7).  Handles
 * integrands rational in unknown functions u[x] and their derivatives,
 * e.g. Integrate[x f'[x] + f[x], x] -> x f[x].  Cheaply gated: returns
 * NULL immediately unless the integrand contains an undefined-function
 * derivative, so ordinary integrands skip it. */
static Expr* try_undefined(Expr* f, Expr* x) {
    return integrate_unknown_try(f, x);
}

/* Stage 1: rational-function integrator.  Returns the antiderivative on
 * success, NULL when the integrand is non-rational or pmint gives up. */
static Expr* try_rational(Expr* f, Expr* x) {
    if (!is_polynomial_in(f, x) && !is_rational_in(f, x)) return NULL;
    Expr* result = call_stage("Integrate`BronsteinRational", f, x);
    if (!result) return NULL;
    if (result_is_unresolved(result, "Integrate`BronsteinRational")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage 1b: derivative-divides substitution.  Recognises integrands of the
 * shape c h(u(x)) u'(x) and reduces to Integrate[h[u], u].  Runs the fast,
 * branch-correct direct-quotient strategy first and then the more thorough
 * Eliminate/Solve branch-search (which closes radical substitutions such as
 * u = Sqrt[Tan[x]]); the latter's Eliminate diagnostics are muted while the
 * integrator drives it (see integrate_derivdivides.c). */
static Expr* try_derivdivides(Expr* f, Expr* x) {
    return integrate_derivdivides_full(f, x);
}

/* Stage 1c: linear-radical substitution.  Recognises a rational function of x
 * and radicals (a x + b)^(m/n) of one shared linear argument, rationalises via
 * u = (a x + b)^(1/n), integrates the rational result, and verifies. */
static Expr* try_linrad(Expr* f, Expr* x) {
    return integrate_linrad_try(f, x);
}

/* Stage 1d: quadratic-radical (Euler) substitution.  Recognises a rational
 * function of x and square roots (a x^2 + b x + c)^(m/2) of one shared quadratic
 * argument, applies a single real-valued Euler substitution, integrates the
 * rational result, and back-substitutes. */
static Expr* try_quadrad(Expr* f, Expr* x) {
    return integrate_quadrad_try(f, x);
}

/* Stage 1e: linear-ratio-radical (Möbius) substitution.  Recognises a rational
 * function of x and radicals ((a x + b)/(c x + d))^(m/n) of one shared
 * linear-fractional argument, rationalises via u = ((a x + b)/(c x + d))^(1/n),
 * integrates the rational result, and back-substitutes. */
static Expr* try_linratiorad(Expr* f, Expr* x) {
    return integrate_linratiorad_try(f, x);
}

/* Stage 1f: continuous Weierstrass substitution (Jeffrey & Rich 1994).
 * Recognises a rational function of the trig kernels Sin/Cos/Tan/Cot/Sec/Csc[x]
 * (or the hyperbolic kernels) carrying a kernel in a denominator, substitutes
 * u = Tan[x/2] (Tanh[x/2] for hyperbolic), integrates the resulting rational
 * function of u, back-substitutes, and -- for the trig case -- adds the
 * K Floor[(x - b)/p] secular term that removes the spurious discontinuities at
 * the poles of Tan[x/2].  Runs ahead of Risch-Norman so genuine rational-trig
 * integrands get the real, continuous antiderivative rather than pmint's
 * complex-logarithm form. */
static Expr* try_weierstrass(Expr* f, Expr* x) {
    return integrate_jeffrey_try(f, x);
}

/* Stage 2: Risch-Norman heuristic (Bronstein pmint). */
static Expr* try_risch(Expr* f, Expr* x) {
    Expr* result = call_stage("Integrate`RischNorman", f, x);
    if (!result) return NULL;
    if (result_is_unresolved(result, "Integrate`RischNorman")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage 3: CRC integral table.  Loaded from disk on first invocation
 * (lazy) so unaffected sessions pay nothing for it. */
#define MAX_CRC_DEPTH 256
static int crc_depth = 0;
static bool crc_load_attempted = false;
static bool crc_load_succeeded  = false;

/* Try Get["src/internal/CRCMathTablesIntegrals.m"]; tolerate the same
 * relative-then-fallback pattern test_integrals.c uses so it works
 * from any CWD with a sane layout. */
static void crc_lazy_load(void) {
    if (crc_load_attempted) return;
    crc_load_attempted = true;

    const char* paths[] = {
        "Get[\"src/internal/CRCMathTablesIntegrals.m\"]",
        "Get[\"../src/internal/CRCMathTablesIntegrals.m\"]",
        "Get[\"../../src/internal/CRCMathTablesIntegrals.m\"]",
        NULL
    };
    for (const char** p = paths; *p; p++) {
        Expr* parsed = parse_expression(*p);
        if (!parsed) continue;
        Expr* res = evaluate(parsed);
        expr_free(parsed);
        bool failed = res && res->type == EXPR_SYMBOL
                          && strcmp(res->data.symbol, "$Failed") == 0;
        if (res) expr_free(res);
        if (!failed) { crc_load_succeeded = true; return; }
    }
    fprintf(stderr,
        "Integrate`CRCTable::nofile: cannot locate "
        "src/internal/CRCMathTablesIntegrals.m on disk.\n");
}

static Expr* try_crctable(Expr* f, Expr* x) {
    crc_lazy_load();
    if (!crc_load_succeeded) return NULL;

    if (crc_depth >= MAX_CRC_DEPTH) {
        fprintf(stderr,
            "Integrate`CRCTable::depth: rule recursion exceeded %d levels; "
            "the table has a divergent rule for this integrand.\n",
            MAX_CRC_DEPTH);
        return NULL;
    }
    crc_depth++;
    Expr* result = call_stage("Integrate`CRCTable", f, x);
    crc_depth--;
    if (!result) return NULL;

    /* "Failed" sentinels: either the public head is unresolved (no
     * rule matched the CRCTable wrapper) or the internal IntegrateTable
     * head leaked through (table lookup found no matching formula).
     * The IntegrateTable scan is deep — a partial-match rule like
     * `IntegrateTable[a_ f_, x_] /; FreeQ[a, x] := a IntegrateTable[f, x]`
     * can leave `IntegrateTable[...]` nested under Times/Plus when the
     * inner call has no matching formula; that must also count as
     * unresolved so Integrate falls back to its own unevaluated form. */
    if (result_is_unresolved(result, "Integrate`CRCTable") ||
        result_contains_head(result, "IntegrateTable")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Method-option parsing.  Mirrors the canonical SYM_Method / SYM_Rule
 * idiom (see src/list.c:1480-1491 and src/facint.c:1276-1310). */
typedef enum {
    METHOD_AUTOMATIC = 0,
    METHOD_RATIONAL,
    METHOD_DERIVATIVE_DIVIDES,
    METHOD_LINEAR_RADICALS,
    METHOD_QUADRATIC_RADICALS,
    METHOD_LINEAR_RATIO_RADICALS,
    METHOD_WEIERSTRASS,
    METHOD_RISCH,
    METHOD_CRCTABLE,
    METHOD_UNDEFINED,
    METHOD_INVALID
} IntegrateMethod;

static IntegrateMethod parse_method_option(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return METHOD_INVALID;
    if (opt->data.function.head->type != EXPR_SYMBOL) return METHOD_INVALID;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return METHOD_INVALID;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_Method)
        return METHOD_INVALID;

    /* Accept either a string ("Automatic") or the symbol Automatic. */
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)
        return METHOD_AUTOMATIC;
    if (rhs->type != EXPR_STRING) return METHOD_INVALID;
    if (strcmp(rhs->data.string, "Automatic")   == 0) return METHOD_AUTOMATIC;
    if (strcmp(rhs->data.string, "BronsteinRational") == 0) return METHOD_RATIONAL;
    if (strcmp(rhs->data.string, "DerivativeDivides") == 0) return METHOD_DERIVATIVE_DIVIDES;
    if (strcmp(rhs->data.string, "LinearRadicals") == 0) return METHOD_LINEAR_RADICALS;
    if (strcmp(rhs->data.string, "QuadraticRadicals") == 0) return METHOD_QUADRATIC_RADICALS;
    if (strcmp(rhs->data.string, "LinearRatioRadicals") == 0) return METHOD_LINEAR_RATIO_RADICALS;
    if (strcmp(rhs->data.string, "Weierstrass") == 0) return METHOD_WEIERSTRASS;
    if (strcmp(rhs->data.string, "RischNorman") == 0) return METHOD_RISCH;
    if (strcmp(rhs->data.string, "CRCTable")    == 0) return METHOD_CRCTABLE;
    if (strcmp(rhs->data.string, "Undefined")   == 0) return METHOD_UNDEFINED;
    return METHOD_INVALID;
}

Expr* builtin_integrate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    /* The integration variable must be a (single) symbol so we can
     * compute derivatives, partial fractions, ... in it. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Applied InterpolatingFunction objects integrate to a fresh antiderivative
     * InterpolatingFunction, mirroring how D differentiates them.  This is
     * handled before the inexact-rationalisation scan below: the object embeds
     * Real sample data that must NOT be force-rationalised. */
    {
        Expr* interp = integrate_interp(f, x);
        if (interp) return interp;
    }

    /* Parse the optional Method -> "..." option. */
    IntegrateMethod method = METHOD_AUTOMATIC;
    if (argc == 3) {
        method = parse_method_option(res->data.function.args[2]);
        if (method == METHOD_INVALID) {
            static uint64_t last_warned_hash = 0;
            uint64_t h = expr_hash(res);
            if (h != last_warned_hash) {
                fprintf(stderr,
                    "Integrate::method: Method option value is not one of "
                    "\"Automatic\", \"BronsteinRational\", \"DerivativeDivides\", "
                    "\"LinearRadicals\", \"QuadraticRadicals\", "
                    "\"LinearRatioRadicals\", \"Weierstrass\", \"RischNorman\", "
                    "\"CRCTable\".\n");
                last_warned_hash = h;
            }
            return NULL;
        }
    }

    /* Inexact integrands: route through the shared preprocessor in
     * common.c -- if `f` contains any inexact leaf, force-rationalise
     * it (with bit-exact ½-ulp fallback so transcendental floats like
     * N[Pi] still become rationals), integrate exactly, then
     * numericalise the result so the user observes inexact-in /
     * inexact-out semantics.  Same contract as Solve.
     *
     * The scan also captures the *minimum* precision (in bits) across
     * every inexact leaf, which is then used both as the rationalisation
     * tolerance and as the precision of the final numericalised result. */
    CommonInexactInfo inexact = common_scan_inexact(f);
    Expr* coerced = NULL;
    Expr* effective_f = f;
    if (inexact.has_inexact) {
        coerced = common_rationalize_input(f, inexact.min_bits);
        if (!coerced) return NULL;
        effective_f = coerced;
    }

    Expr* result = NULL;
    switch (method) {
        case METHOD_AUTOMATIC:
            result = try_undefined(effective_f, x);
            if (!result) result = try_rational(effective_f, x);
            if (!result) result = try_linrad(effective_f, x);
            if (!result) result = try_quadrad(effective_f, x);
            if (!result) result = try_linratiorad(effective_f, x);
            /* Weierstrass before derivative-divides: it is a domain-specific,
             * deterministic algorithm for rational trig/hyperbolic integrands
             * that is guaranteed to close (and verified by construction), so it
             * runs ahead of the more expensive Eliminate/Solve substitution
             * search and ahead of Risch-Norman's complex-logarithm forms. */
            if (!result) result = try_weierstrass(effective_f, x);
            if (!result) result = try_derivdivides(effective_f, x);
            if (!result) result = try_risch(effective_f, x);
            if (!result) result = try_crctable(effective_f, x);
            break;
        case METHOD_RATIONAL:
            result = try_rational(effective_f, x);
            break;
        case METHOD_DERIVATIVE_DIVIDES:
            /* Explicit method runs the thorough Eliminate/Solve search in
             * addition to the direct quotient (the Automatic cascade uses the
             * direct strategy only). */
            result = integrate_derivdivides_full(effective_f, x);
            break;
        case METHOD_LINEAR_RADICALS:
            result = try_linrad(effective_f, x);
            break;
        case METHOD_QUADRATIC_RADICALS:
            result = try_quadrad(effective_f, x);
            break;
        case METHOD_LINEAR_RATIO_RADICALS:
            result = try_linratiorad(effective_f, x);
            break;
        case METHOD_WEIERSTRASS:
            result = integrate_jeffrey_full(effective_f, x);
            break;
        case METHOD_RISCH:
            result = try_risch(effective_f, x);
            break;
        case METHOD_CRCTABLE:
            result = try_crctable(effective_f, x);
            break;
        case METHOD_UNDEFINED:
            result = try_undefined(effective_f, x);
            break;
        case METHOD_INVALID:
            break;  /* unreachable: handled above */
    }

    if (coerced) expr_free(coerced);

    if (inexact.has_inexact && result) {
        Expr* numeric = common_numericalize_result(result, inexact.min_bits);
        expr_free(result);
        result = numeric;
    }

    return result;
}

void integrate_init(void) {
    symtab_add_builtin("Integrate", builtin_integrate);
    symtab_get_def("Integrate")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring("Integrate",
        "Integrate[f, x] gives the indefinite integral of f with respect to x.\n"
        "Integrate[f, x, Method -> \"<name>\"] dispatches directly to a single\n"
        "subroutine, bypassing the default cascade.  Accepted method names:\n"
        "  \"Automatic\"          — try BronsteinRational, then RischNorman, then CRCTable (default)\n"
        "  \"BronsteinRational\"  — Integrate`BronsteinRational (polynomial / rational)\n"
        "  \"DerivativeDivides\"  — Integrate`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)\n"
        "  \"LinearRadicals\"     — Integrate`LinearRadicals (rationalise radicals of a x + b)\n"
        "  \"QuadraticRadicals\"  — Integrate`QuadraticRadicals (Euler substitution for Sqrt[a x^2 + b x + c])\n"
        "  \"LinearRatioRadicals\" — Integrate`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d))\n"
        "  \"Weierstrass\"        — Integrate`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution)\n"
        "  \"RischNorman\"        — Integrate`RischNorman (Bronstein pmint heuristic)\n"
        "  \"CRCTable\"           — Integrate`CRCTable (lazy-loaded CRC integral table)\n"
        "  \"Undefined\"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)\n"
        "Named methods are strict: failure returns unevaluated, with no fallback.\n"
        "The CRCTable rules are loaded from disk on first use only.\n"
        "An applied 1-D InterpolatingFunction integrates to its antiderivative\n"
        "InterpolatingFunction (mirroring D).");

    /* Initialise the Integrate` package: HermiteReduce, IntegratePolynomial,
     * helpers, and the explicit `Integrate`BronsteinRational` entry. */
    intrat_init();

    /* Undefined-function integrator (Roach §1.7): Integrate`Undefined. */
    integrate_unknown_init();

    /* Derivative-divides substitution: Integrate`DerivativeDivides. */
    integrate_derivdivides_init();

    /* Linear-radical rationalising substitution: Integrate`LinearRadicals. */
    integrate_linrad_init();

    /* Quadratic-radical (Euler) substitution: Integrate`QuadraticRadicals. */
    integrate_quadrad_init();

    /* Linear-ratio-radical (Möbius) substitution: Integrate`LinearRatioRadicals. */
    integrate_linratiorad_init();

    /* Continuous Weierstrass substitution (Jeffrey & Rich 1994):
     * Integrate`Weierstrass. */
    integrate_jeffrey_init();

    /* Initialise the parallel-Risch / Risch-Norman heuristic
     * (Bronstein's pmint).  Provides `Integrate`RischNorman[f, x]`,
     * the fall-through for transcendental integrands. */
    intrischnorman_init();

    /* Pre-register Integrate`CRCTable so ?Integrate`CRCTable shows a
     * docstring even before the lazy load fires.  The actual rule
     * definition is added when CRCMathTablesIntegrals.m is Get-loaded
     * on first call — that DownValue install must NOT be blocked by
     * Protected, so we deliberately leave attributes empty here.  The
     * lazy load itself promotes Integrate`CRCTable to Protected via
     * SetAttributes inside CRCMathTablesIntegrals.m. */
    symtab_get_def("Integrate`CRCTable");
    symtab_set_docstring("Integrate`CRCTable",
        "Integrate`CRCTable[f, x] looks up f in the CRC Standard Mathematical\n"
        "Tables (31st ed.) integral table. Returns the antiderivative on a hit,\n"
        "or the unevaluated form if no rule applies. The table file\n"
        "(src/internal/CRCMathTablesIntegrals.m) is loaded from disk on first\n"
        "use only.");
}

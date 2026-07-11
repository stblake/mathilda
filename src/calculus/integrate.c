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
#include "integrate_chebychev.h"
#include "integrate_goursat.h"
#include "integrate_jeffrey.h"
#include "integrate_newton_leibniz.h"
#include "integrate_symmetry.h"
#include "integrate_beta.h"
#include "integrate_line.h"
#include "integrate_residue.h"
#include "integrate_diffunderint.h"
#include "integrate_ramanujan.h"
#include "intrat.h"
#include "intrischnorman.h"
#include "integrate_risch_transcendental.h"
#include "intsimp.h"
#include "common.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "loadmodule.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "series.h"

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

/* True iff `f` contains `x` as a subexpression (structural FreeQ negation). */
static bool depends_on_var(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return true;
    if (f->type == EXPR_FUNCTION) {
        if (depends_on_var(f->data.function.head, x)) return true;
        for (size_t i = 0; i < f->data.function.arg_count; i++)
            if (depends_on_var(f->data.function.args[i], x)) return true;
    }
    return false;
}

/* True iff `f` contains a Power[b, e] whose base involves `x` and whose
 * exponent is not a concrete number -- a symbolic-exponent power such as x^k,
 * x^(k-1) or (1-x)^l.  Such an integrand is provably not a rational function of
 * x, yet feeding it to the Together / PolynomialGCD rationality probe (or the
 * derivative-divides quotient normalisation) drives a pseudo-remainder-sequence
 * blow-up over the symbolic powers as independent generators.  Callers use it
 * to decline those Together-backed paths structurally, up front. */
static bool has_symbolic_power_in(const Expr* f, const Expr* x) {
    if (!f || f->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)f, SYM_Power) && f->data.function.arg_count == 2) {
        const Expr* b = f->data.function.args[0];
        const Expr* e = f->data.function.args[1];
        if (depends_on_var(b, x) && !expr_is_numeric_like(e)) return true;
    }
    if (has_symbolic_power_in(f->data.function.head, x)) return true;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (has_symbolic_power_in(f->data.function.args[i], x)) return true;
    return false;
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
    /* A symbolic-exponent power of x (x^k, (1-x)^l, ...) is never rational in
     * x; decline structurally so the Together-based rationality probe below is
     * never handed an integrand that sends PolynomialGCD into a blow-up. */
    if (has_symbolic_power_in(f, x)) return NULL;
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

/* Stage 1e2: Chebychev binomial differential.  Recognises an integrand
 * x^p (a x^r + b)^q with p, q, r rational and a, b free of x, and -- when one
 * of q, (p+1)/r, q+(p+1)/r is an integer (Chebychev's theorem) -- applies the
 * matching rationalising substitution (x = u^N, u^s = a x^r + b, or
 * u = x^r then t^s = (a u + b)/u), integrates the rational result, and
 * back-substitutes.  Recognition is a single structural scan, so it is cheap
 * enough to run ahead of the Eliminate/Solve search in derivative-divides.
 * Non-elementary binomials return NULL, so the cascade falls through to the
 * later methods (which may one day carry special-function representations). */
static Expr* try_chebychev(Expr* f, Expr* x) {
    return integrate_chebychev_try(f, x);
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

/* Stage 1e3: Goursat's pseudo-elliptic algorithm and its cube-/fourth-root
 * generalisations.  Recognises F(x)/R(x)^p (p in {1/2,1/3,2/3,1/4,3/4}); when
 * a Mobius automorphism of the roots of R makes the integrand pseudo-elliptic,
 * descends to genus-0 curves and integrates in closed form.  Deterministic and
 * correct by construction, so it runs ahead of the Eliminate/Solve search and
 * Risch-Norman.  Non-elementary (obstructed) integrands return NULL. */
static Expr* try_goursat(Expr* f, Expr* x) {
    return integrate_goursat_try(f, x);
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

/* Stage 2b: recursive transcendental Risch integrator.
 * A decision procedure over a differential transcendental tower, with
 * rational / logarithmic / exponential / special-function cases, each
 * correct by construction (it fires only behind an exact certificate, so
 * it emits no wrong closed forms and needs no differentiation check).
 * Distinct from the parallel-Risch (pmint) heuristic Integrate`RischNorman
 * run just above it. */
static Expr* try_rischtranscendental(Expr* f, Expr* x) {
    Expr* result = call_stage("Integrate`RischTranscendental", f, x);
    if (!result) return NULL;
    if (result_is_unresolved(result, "Integrate`RischTranscendental")) {
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

/* Load the CRC integral table lazily. Path resolution (mathilda_resolve_
 * internal) is CWD-independent and matches what LoadModule uses, so this works
 * from a relocated/installed binary as well as from the repo root or tests/. */
static void crc_lazy_load(void) {
    if (crc_load_attempted) return;
    crc_load_attempted = true;

    char path[2048];
    if (!mathilda_resolve_internal("CRCMathTablesIntegrals.m",
                                   path, sizeof(path))) {
        fprintf(stderr,
            "Integrate`CRCTable::nofile: cannot locate "
            "src/internal/CRCMathTablesIntegrals.m on disk.\n");
        return;
    }

    int opened = 0;
    Expr* res = mathilda_run_file(path, &opened);
    bool failed = res && res->type == EXPR_SYMBOL
                      && strcmp(res->data.symbol, "$Failed") == 0;
    if (res) expr_free(res);
    if (opened && !failed) crc_load_succeeded = true;
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
    METHOD_CHEBYCHEV,
    METHOD_GOURSAT,
    METHOD_WEIERSTRASS,
    METHOD_RISCH,
    METHOD_RISCH_TRANSCENDENTAL,
    METHOD_CRCTABLE,
    METHOD_UNDEFINED,
    METHOD_NEWTON_LEIBNIZ,   /* definite-only: selects the real-axis FTC mechanism */
    METHOD_LINE_INTEGRAL,    /* definite-only: selects the complex contour mechanism */
    METHOD_RESIDUE,          /* definite-only: selects the residue-theorem mechanism */
    METHOD_DIFF_UNDER_INT,   /* definite-only: differentiation under the integral sign */
    METHOD_RAMANUJAN,        /* definite-only: Mellin / Ramanujan Master Theorem (half-line) */
    METHOD_SYMMETRY,         /* definite-only: origin-symmetry (odd -> 0, even -> 2 half) */
    METHOD_BETA,             /* definite-only: Euler-Beta on [0,1] */
    METHOD_TRIG_POWER,       /* definite-only: Sin^m Cos^n over [0,Pi/2]/[0,Pi]/[0,2Pi] */
    METHOD_SINPOW_MONO,      /* definite-only: Sin[r x]^k / x^m on [0,Inf) */
    METHOD_OSC_POWER,        /* definite-only: Cos/Sin[b x^n] on [0,Inf) */
    METHOD_RATIONAL_LOG,     /* definite-only: R(x) Log[x]^n on [0,Inf) */
    METHOD_INVALID
} IntegrateMethod;

/* Map a method-name string to its enum, or METHOD_INVALID if unrecognised. */
static IntegrateMethod method_from_string(const char* s) {
    if (strcmp(s, "Automatic")   == 0) return METHOD_AUTOMATIC;
    if (strcmp(s, "BronsteinRational") == 0) return METHOD_RATIONAL;
    if (strcmp(s, "DerivativeDivides") == 0) return METHOD_DERIVATIVE_DIVIDES;
    if (strcmp(s, "LinearRadicals") == 0) return METHOD_LINEAR_RADICALS;
    if (strcmp(s, "QuadraticRadicals") == 0) return METHOD_QUADRATIC_RADICALS;
    if (strcmp(s, "LinearRatioRadicals") == 0) return METHOD_LINEAR_RATIO_RADICALS;
    if (strcmp(s, "ChebychevAlgebraic") == 0) return METHOD_CHEBYCHEV;
    if (strcmp(s, "GoursatAlgebraic") == 0) return METHOD_GOURSAT;
    if (strcmp(s, "Weierstrass") == 0) return METHOD_WEIERSTRASS;
    if (strcmp(s, "RischNorman") == 0) return METHOD_RISCH;
    if (strcmp(s, "RischTranscendental") == 0) return METHOD_RISCH_TRANSCENDENTAL;
    if (strcmp(s, "CRCTable")    == 0) return METHOD_CRCTABLE;
    if (strcmp(s, "Undefined")   == 0) return METHOD_UNDEFINED;
    if (strcmp(s, "NewtonLeibniz") == 0) return METHOD_NEWTON_LEIBNIZ;
    if (strcmp(s, "LineIntegral") == 0) return METHOD_LINE_INTEGRAL;
    if (strcmp(s, "Residue") == 0 || strcmp(s, "ContourResidue") == 0)
        return METHOD_RESIDUE;
    if (strcmp(s, "DiffUnderInt") == 0 ||
        strcmp(s, "DifferentiationUnderIntegral") == 0)
        return METHOD_DIFF_UNDER_INT;
    if (strcmp(s, "Symmetry") == 0) return METHOD_SYMMETRY;
    if (strcmp(s, "Beta") == 0) return METHOD_BETA;
    if (strcmp(s, "TrigPower") == 0) return METHOD_TRIG_POWER;
    if (strcmp(s, "SinPowerMonomial") == 0) return METHOD_SINPOW_MONO;
    if (strcmp(s, "OscillatoryPower") == 0) return METHOD_OSC_POWER;
    if (strcmp(s, "RationalLog") == 0) return METHOD_RATIONAL_LOG;
    if (strcmp(s, "RamanujanMasterTheorem") == 0 || strcmp(s, "Mellin") == 0)
        return METHOD_RAMANUJAN;
    return METHOD_INVALID;
}

/* Parse the `Method -> ...` option.  The RHS may be:
 *   - the symbol `Automatic` or a bare method-name string;
 *   - a list `{"<method>", subopt -> val, ...}` carrying method sub-options
 *     (mirrors the FactorInteger list-form in src/facint.c:1173).  The only
 *     recognised sub-option is `"Substitution" -> u`, valid only for
 *     "DerivativeDivides": on success `*out_sub` is set to an owned copy of `u`
 *     (the caller must free it).
 * `*out_sub` is always initialised (to NULL) before any other work. */
static IntegrateMethod parse_method_option(Expr* opt, Expr** out_sub) {
    *out_sub = NULL;
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

    /* List form: {"<method>", subopt -> val, ...}. */
    if (rhs->type == EXPR_FUNCTION &&
        rhs->data.function.head->type == EXPR_SYMBOL &&
        rhs->data.function.head->data.symbol == SYM_List &&
        rhs->data.function.arg_count >= 1) {
        Expr* mname = rhs->data.function.args[0];
        if (mname->type != EXPR_STRING) return METHOD_INVALID;
        IntegrateMethod m = method_from_string(mname->data.string);
        if (m == METHOD_INVALID) return METHOD_INVALID;
        for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
            Expr* so = rhs->data.function.args[i];
            if (so->type != EXPR_FUNCTION ||
                so->data.function.head->type != EXPR_SYMBOL) goto bad_subopt;
            const char* sh = so->data.function.head->data.symbol;
            if ((sh != SYM_Rule && sh != SYM_RuleDelayed) ||
                so->data.function.arg_count != 2) goto bad_subopt;
            Expr* skey = so->data.function.args[0];
            Expr* sval = so->data.function.args[1];
            if (skey->type == EXPR_STRING &&
                strcmp(skey->data.string, "Substitution") == 0 &&
                m == METHOD_DERIVATIVE_DIVIDES) {
                if (*out_sub) expr_free(*out_sub);
                *out_sub = expr_copy(sval);
                continue;
            }
            goto bad_subopt;
        }
        return m;
    bad_subopt:
        if (*out_sub) { expr_free(*out_sub); *out_sub = NULL; }
        return METHOD_INVALID;
    }

    if (rhs->type != EXPR_STRING) return METHOD_INVALID;
    return method_from_string(rhs->data.string);
}

/* Parse a definite integral's trailing `Method -> "..."` option.  `*name` is set
 * to a borrowed method-name string to pass through to the inner indefinite
 * Integrate (NULL for Automatic and for the definite mechanisms NewtonLeibniz /
 * LineIntegral / Residue, which select the mechanism itself).  `*mech` is set to
 * the chosen definite mechanism so the dispatcher knows whether to engage the
 * residue method and/or fall through to Newton-Leibniz.  Returns false on a
 * malformed / unrecognised option. */
static bool definite_parse_method(Expr* opt, const char** name,
                                  IntegrateMethod* mech) {
    *name = NULL;
    *mech = METHOD_AUTOMATIC;
    if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2)
        return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_Method) return false;
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) return true;
    if (rhs->type != EXPR_STRING) return false;
    IntegrateMethod m = method_from_string(rhs->data.string);
    if (m == METHOD_INVALID) return false;
    *mech = m;
    /* NewtonLeibniz / LineIntegral / Residue / DiffUnderInt name the definite
     * mechanism itself; the actual mechanism is chosen from the spec type, so
     * they pass NULL through to the inner indefinite Integrate. */
    if (m == METHOD_AUTOMATIC || m == METHOD_NEWTON_LEIBNIZ ||
        m == METHOD_LINE_INTEGRAL || m == METHOD_RESIDUE ||
        m == METHOD_DIFF_UNDER_INT || m == METHOD_RAMANUJAN ||
        m == METHOD_SYMMETRY || m == METHOD_BETA ||
        m == METHOD_TRIG_POWER || m == METHOD_SINPOW_MONO ||
        m == METHOD_OSC_POWER || m == METHOD_RATIONAL_LOG) return true;
    *name = rhs->data.string;   /* borrowed: valid while `res` is alive */
    return true;
}

/* True iff `opt` is a Rule/RuleDelayed whose LHS is the symbol `sym`. */
static bool option_lhs_is(Expr* opt, const char* sym) {
    if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
    Expr* lhs = opt->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && lhs->data.symbol == sym;
}

/* Definite / iterated integration Integrate[f, {x,a,b}, {y,c,d}, ..., opts].
 * Reduces innermost-first (the last spec is the inner integral) so an inner
 * bound may depend on an outer variable.  Returns NULL (unevaluated) if any
 * stage cannot be evaluated. */
static Expr* integrate_definite(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr* f = res->data.function.args[0];

    /* Leading run of range / contour specs, then at most one trailing Method
     * option.  A contour spec `{x, z0, ..., zn}` (polyline, or any non-real
     * endpoint) is a superset of the real 3-element range spec. */
    size_t nspecs = 0;
    while (1 + nspecs < argc &&
           integrate_line_is_contour_spec(res->data.function.args[1 + nspecs]))
        nspecs++;
    if (nspecs == 0) return NULL;

    const char* method = NULL;
    IntegrateMethod mech = METHOD_AUTOMATIC;
    Expr* assumptions = NULL;   /* borrowed from res: the Assumptions option value */
    bool principal_value = false;
    size_t tail = 1 + nspecs;
    /* Trailing options, in any order: `Method -> ...` and/or `Assumptions -> ...`.
     * Any other trailing argument (unknown option, stray expression) is rejected
     * with the Integrate::method diagnostic, leaving the call unevaluated. */
    for (size_t t = tail; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (option_lhs_is(opt, SYM_Assumptions)) {
            assumptions = opt->data.function.args[1];   /* borrowed */
            continue;
        }
        /* PrincipalValue -> True / False (not an interned system symbol; match
         * by name).  Selects the Cauchy principal value for interior poles. */
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol == SYM_Rule ||
             opt->data.function.head->data.symbol == SYM_RuleDelayed) &&
            opt->data.function.args[0]->type == EXPR_SYMBOL &&
            strcmp(opt->data.function.args[0]->data.symbol, "PrincipalValue") == 0) {
            Expr* rhs = opt->data.function.args[1];
            principal_value = (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_True);
            continue;
        }
        if (option_lhs_is(opt, SYM_Method)) {
            if (definite_parse_method(opt, &method, &mech)) continue;
        }
        static uint64_t last_warned_hash = 0;
        uint64_t h = expr_hash(res);
        if (h != last_warned_hash) {
            fprintf(stderr,
                "Integrate::method: Method option value is not a "
                "recognised integration method name.\n");
            last_warned_hash = h;
        }
        return NULL;
    }

    /* Fold from the innermost (last) spec outward. */
    Expr* cur = expr_copy(f);
    for (size_t k = nspecs; k >= 1; k--) {
        Expr* spec = res->data.function.args[1 + (k - 1)];
        Expr* r;
        if (integrate_line_spec_is_complex(spec)) {
            /* Complex line / contour integral (non-real endpoint or polyline). */
            r = integrate_line_from_spec(cur, spec, method);
        } else {
            /* Real-axis definite integral.  Under Automatic (or an explicit
             * Method -> "Residue") try the residue-theorem method first: for
             * improper / periodic forms it yields cleaner, NIntegrate-verified
             * closed forms.  On NULL it falls through to Newton-Leibniz, except
             * when the user pinned "Residue" (strict: no fallback). */
            Expr* x = spec->data.function.args[0];
            Expr* a = spec->data.function.args[1];
            Expr* b = spec->data.function.args[2];
            r = NULL;
            bool diverges = false;
            if (mech == METHOD_AUTOMATIC || mech == METHOD_RESIDUE)
                r = integrate_residue_try(cur, x, a, b, assumptions, &diverges);
            /* The residue method conclusively found a pole on the integration
             * contour: the integral does not converge.  Emit Integrate::idiv
             * and stop -- do NOT fall through to Newton-Leibniz (which would
             * re-derive the same divergence and emit a duplicate warning). */
            if (!r && diverges) {
                integrate_emit_idiv(cur, a, b);
                expr_free(cur);
                return NULL;
            }
            /* Origin-symmetry reduction, after residue (which keeps its clean
             * closed forms for the families it owns) but before the general FTC:
             * an odd integrand over a symmetric interval is 0 with no
             * antiderivative, and an even one halves the interval -- often
             * unlocking the half-line Ramanujan/Mellin path for a non-rational
             * even integrand that residue does not own.  It only claims a value
             * once the half converges, so a divergent principal value is never
             * reported as 0, and NULL always falls through unchanged. */
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_SYMMETRY))
                r = integrate_symmetry_try(cur, x, a, b, assumptions);
            /* Newton-Leibniz (FTC) unless the user pinned Residue, the
             * parameter-differentiation mechanism, the Ramanujan/Mellin
             * mechanism, the symmetry mechanism, or a Beta mechanism. */
            if (!r && mech != METHOD_RESIDUE && mech != METHOD_DIFF_UNDER_INT &&
                mech != METHOD_RAMANUJAN && mech != METHOD_SYMMETRY &&
                mech != METHOD_BETA && mech != METHOD_TRIG_POWER &&
                mech != METHOD_SINPOW_MONO && mech != METHOD_OSC_POWER &&
                mech != METHOD_RATIONAL_LOG)
                r = integrate_newton_leibniz_try_pv(cur, x, a, b, method,
                                                    principal_value);
            /* Euler-Beta reductions: x^(k-1)(1-x)^(l-1) on [0,1] and
             * Sin^m Cos^n over a canonical trig interval -- non-elementary
             * antiderivatives (incomplete Beta) that FTC cannot reach.  After
             * Newton-Leibniz (which owns the integer-power cases via an
             * elementary antiderivative), before Ramanujan. */
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_BETA))
                r = integrate_beta_try(cur, x, a, b, assumptions);
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_TRIG_POWER))
                r = integrate_trigpower_try(cur, x, a, b, assumptions);
            /* Mellin / Ramanujan Master Theorem: half-line ∫₀^∞ x^(s-1) f(x) dx
             * of a transcendental f (Gaussian moments, Gamma/Bessel/trig
             * transforms) that residue and FTC do not close.  Under Automatic it
             * runs before DiffUnderInt; the pinned mechanism has no fallback. */
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_RAMANUJAN ||
                       mech == METHOD_OSC_POWER))
                r = integrate_ramanujan_try(cur, x, a, b, assumptions);
            /* Sin[r x]^k / x^m half-line (ssp) and R(x) Log[x]^n (log*rat):
             * under Automatic these are pre-passes inside Ramanujan; the pinned
             * methods route directly here. */
            if (!r && mech == METHOD_SINPOW_MONO)
                r = integrate_sinpowmono_try(cur, x, a, b, assumptions);
            if (!r && mech == METHOD_RATIONAL_LOG)
                r = integrate_ratlogpow_try(cur, x, a, b, assumptions);
            /* Differentiation under the integral sign: last resort under
             * Automatic (parameter-dependent improper/periodic integrals that
             * residue and FTC cannot close), or the pinned mechanism. */
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_DIFF_UNDER_INT))
                r = integrate_diffunderint_try(cur, x, a, b, assumptions);
        }
        expr_free(cur);
        if (!r) return NULL;
        cur = r;
    }
    return cur;
}

Expr* builtin_integrate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f = res->data.function.args[0];

    /* Thread over a list integrand: Integrate[{f1, ..., fn}, spec...] gives
     * {Integrate[f1, spec...], ..., Integrate[fn, spec...]}.  Integrate is
     * deliberately NOT Listable (that would also thread over a `{x, a, b}`
     * range spec), so the integrand-only threading is handled explicitly here.
     * We build the List of sub-Integrate calls and let the evaluator recurse
     * into each element (which stays symbolic if a piece cannot be done). */
    if (f->type == EXPR_FUNCTION &&
        f->data.function.head->type == EXPR_SYMBOL &&
        f->data.function.head->data.symbol == SYM_List) {
        size_t n = f->data.function.arg_count;
        Expr** elems = malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** sub = malloc(argc * sizeof(Expr*));
            sub[0] = expr_copy(f->data.function.args[i]);
            for (size_t j = 1; j < argc; j++)
                sub[j] = expr_copy(res->data.function.args[j]);
            elems[i] = expr_new_function(expr_new_symbol(SYM_Integrate), sub, argc);
            free(sub);
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), elems, n);
        free(elems);
        return out;
    }

    /* Definite form: the second argument is a `{x, a, b}` range spec or a
     * `{x, z0, ..., zn}` complex line/contour spec.  Handles the iterated
     * multi-spec form too. */
    if (integrate_line_is_contour_spec(res->data.function.args[1]))
        return integrate_definite(res);

    if (argc > 3) return NULL;

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

    /* SeriesData integrates term-by-term to a new SeriesData. Handle this
     * before the inexact-rationalisation scan and method cascade, none of
     * which understand the SeriesData head. */
    if (f->type == EXPR_FUNCTION &&
        f->data.function.head->type == EXPR_SYMBOL &&
        f->data.function.head->data.symbol == SYM_SeriesData &&
        f->data.function.arg_count == 6) {
        Expr* r = series_integrate(f, x);
        if (r) return r;   /* residue/unsupported -> fall through, stays unevaluated */
    }

    /* Parse the optional Method -> "..." option.  `method_sub` captures the
     * user-pinned substitution from the
     * `{"DerivativeDivides", "Substitution" -> u}` list form (NULL otherwise);
     * it is freed below alongside `coerced`. */
    IntegrateMethod method = METHOD_AUTOMATIC;
    Expr* method_sub = NULL;
    /* An `Assumptions -> ...` option is accepted on the indefinite form for
     * surface compatibility (it constrains parameter domains); the indefinite
     * cascade does not yet consume it, so it is simply skipped here rather than
     * mis-parsed as a Method value. */
    if (argc == 3 && !option_lhs_is(res->data.function.args[2], SYM_Assumptions)) {
        method = parse_method_option(res->data.function.args[2], &method_sub);
        if (method == METHOD_INVALID) {
            static uint64_t last_warned_hash = 0;
            uint64_t h = expr_hash(res);
            if (h != last_warned_hash) {
                fprintf(stderr,
                    "Integrate::method: Method option value is not one of "
                    "\"Automatic\", \"BronsteinRational\", \"DerivativeDivides\", "
                    "\"LinearRadicals\", \"QuadraticRadicals\", "
                    "\"LinearRatioRadicals\", \"ChebychevAlgebraic\", "
                    "\"GoursatAlgebraic\", "
                    "\"Weierstrass\", \"RischNorman\", \"RischTranscendental\", "
                    "\"CRCTable\", \"NewtonLeibniz\", \"LineIntegral\".\n");
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
            /* Chebychev binomial differentials: a fast, deterministic
             * rationalising substitution that closes (correct by construction),
             * so it runs ahead of the Eliminate/Solve search and Risch-Norman. */
            if (!result) result = try_chebychev(effective_f, x);
            /* Goursat pseudo-elliptic (and cube-/fourth-root) reductions:
             * deterministic, correct-by-construction descents to genus-0
             * curves, so they run ahead of the Eliminate/Solve search and
             * Risch-Norman. */
            if (!result) result = try_goursat(effective_f, x);
            /* Weierstrass before derivative-divides: it is a domain-specific,
             * deterministic algorithm for rational trig/hyperbolic integrands
             * that is guaranteed to close (and verified by construction), so it
             * runs ahead of the more expensive Eliminate/Solve substitution
             * search and ahead of Risch-Norman's complex-logarithm forms. */
            if (!result) result = try_weierstrass(effective_f, x);
            if (!result) result = try_derivdivides(effective_f, x);
            if (!result) result = try_risch(effective_f, x);
            /* Recursive transcendental Risch: runs after the pmint heuristic
             * and is correct by construction, so it only adds closed forms
             * the earlier stages missed (logarithmic polynomials, Gaussians
             * -> Erf, exp/x -> ExpIntegralEi, 1/Log -> LogIntegral,
             * Log[1+ax]/x -> PolyLog), never changing an existing answer. */
            if (!result) result = try_rischtranscendental(effective_f, x);
            if (!result) result = try_crctable(effective_f, x);
            break;
        case METHOD_RATIONAL:
            result = try_rational(effective_f, x);
            break;
        case METHOD_DERIVATIVE_DIVIDES:
            /* Explicit method runs the thorough Eliminate/Solve search in
             * addition to the direct quotient (the Automatic cascade uses the
             * direct strategy only).  A user-pinned `"Substitution" -> u`
             * restricts the search to that single kernel. */
            result = method_sub
                   ? integrate_derivdivides_with_sub(effective_f, x, method_sub)
                   : integrate_derivdivides_full(effective_f, x);
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
        case METHOD_CHEBYCHEV:
            result = try_chebychev(effective_f, x);
            break;
        case METHOD_GOURSAT:
            result = try_goursat(effective_f, x);
            break;
        case METHOD_WEIERSTRASS:
            result = integrate_jeffrey_full(effective_f, x);
            break;
        case METHOD_RISCH:
            result = try_risch(effective_f, x);
            break;
        case METHOD_RISCH_TRANSCENDENTAL:
            result = try_rischtranscendental(effective_f, x);
            break;
        case METHOD_CRCTABLE:
            result = try_crctable(effective_f, x);
            break;
        case METHOD_UNDEFINED:
            result = try_undefined(effective_f, x);
            break;
        case METHOD_NEWTON_LEIBNIZ:
        case METHOD_LINE_INTEGRAL:
        case METHOD_RESIDUE:
        case METHOD_DIFF_UNDER_INT:
        case METHOD_RAMANUJAN:
        case METHOD_SYMMETRY:
        case METHOD_BETA:
        case METHOD_TRIG_POWER:
        case METHOD_SINPOW_MONO:
        case METHOD_OSC_POWER:
        case METHOD_RATIONAL_LOG:
            /* Definite-only mechanisms.  Meaningless on the indefinite form
             * Integrate[f, x, Method -> "NewtonLeibniz" / "LineIntegral" /
             * "Residue" / "DiffUnderInt" / "RamanujanMasterTheorem"]; leave
             * unevaluated. */
            break;
        case METHOD_INVALID:
            break;  /* unreachable: handled above */
    }

    /* Single post-method normalisation chokepoint: flatten x-free
     * reciprocal/product powers across every method's output, and — for
     * radical-bearing antiderivatives — recombine the algebraic part and
     * tidy inverse-trig arguments.  Runs on the symbolic result before the
     * inexact path numericalises it.  See intsimp_finalize. */
    if (result) result = intsimp_finalize(result, x);   /* takes ownership */

    if (coerced) expr_free(coerced);
    if (method_sub) expr_free(method_sub);

    if (inexact.has_inexact && result) {
        Expr* numeric = common_numericalize_result(result, inexact.min_bits);
        expr_free(result);
        result = numeric;
    }

    return result;
}

void integrate_init(void) {
    symtab_add_builtin("Integrate", builtin_integrate);
    /* NOT Listable: Listable would thread over a definite-integral range
     * spec `{x, a, b}` element-wise (producing garbage like
     * {Integrate[f,x], Integrate[f,a], Integrate[f,b]}).  The `{x,a,b}` form
     * is recognised explicitly at the top of builtin_integrate and routed to
     * the definite (Newton-Leibniz) path. */
    symtab_get_def("Integrate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate",
        "Integrate[f, x] gives the indefinite integral of f with respect to x.\n"
        "Integrate[f, {x, xmin, xmax}] gives the definite integral by the\n"
        "fundamental theorem of calculus (Method -> \"NewtonLeibniz\").\n"
        "Integrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] gives the iterated\n"
        "multiple integral (innermost/last spec integrated first; inner bounds\n"
        "may depend on outer variables).  See also Integrate`SingularPoints.\n"
        "Integrate[f, x, Method -> \"<name>\"] dispatches directly to a single\n"
        "subroutine, bypassing the default cascade.  Accepted method names:\n"
        "  \"Automatic\"          — try BronsteinRational, then RischNorman, then CRCTable (default)\n"
        "  \"BronsteinRational\"  — Integrate`BronsteinRational (polynomial / rational)\n"
        "  \"DerivativeDivides\"  — Integrate`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)\n"
        "  \"LinearRadicals\"     — Integrate`LinearRadicals (rationalise radicals of a x + b)\n"
        "  \"QuadraticRadicals\"  — Integrate`QuadraticRadicals (Euler substitution for Sqrt[a x^2 + b x + c])\n"
        "  \"LinearRatioRadicals\" — Integrate`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d))\n"
        "  \"ChebychevAlgebraic\" — Integrate`ChebychevAlgebraic (binomial x^p (a x^r + b)^q via Chebychev's theorem)\n"
        "  \"GoursatAlgebraic\"   — Integrate`GoursatAlgebraic (pseudo-elliptic F/R^p, p in {1/2,1/3,2/3,1/4,3/4}, via Mobius eigendescent)\n"
        "  \"Weierstrass\"        — Integrate`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution)\n"
        "  \"RischNorman\"        — Integrate`RischNorman (Bronstein pmint heuristic)\n"
        "  \"RischTranscendental\"       — Integrate`RischTranscendental (recursive transcendental Risch; correct by construction)\n"
        "  \"CRCTable\"           — Integrate`CRCTable (lazy-loaded CRC integral table)\n"
        "  \"Undefined\"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)\n"
        "  \"NewtonLeibniz\"       — real definite integrals via F(b)-F(a) (implicit for the {x,a,b} form)\n"
        "  \"LineIntegral\"        — complex contour integrals (implicit for the {x,z0,...,zn} form)\n"
        "  \"Residue\"             — improper/periodic real definite integrals by the residue theorem\n"
        "                          (rational/Fourier on (-Inf,Inf), rational-in-Sin/Cos over a period,\n"
        "                          principal values, even half-lines); tried before NewtonLeibniz under Automatic\n"
        "  \"DiffUnderInt\"         — parameter-dependent definite integrals by differentiation under the\n"
        "  (\"DifferentiationUnderIntegral\") integral sign (Feynman's trick): Integrate`DiffUnderInt;\n"
        "                          Laplace/Fourier, sinc, and even-rational half-line families;\n"
        "                          tried after Residue and NewtonLeibniz in the definite cascade\n"
        "  \"RamanujanMasterTheorem\" — half-line Int_0^Inf x^(s-1) f(x) dx by the Mellin transform /\n"
        "  (\"Mellin\")              Ramanujan Master Theorem: Integrate`RamanujanMasterTheorem;\n"
        "                          exp/Gaussian/algebraic/Cos/Sin/ArcTan/Log/BesselJ/pFq/PolyLog\n"
        "                          kernels (monomial x^k substitution; Erf, incomplete Gamma, BesselJ^2\n"
        "                          reduced to pFq); also the exp-geometric kernel 1/(E^(c x)+g)\n"
        "                          (Bose-Einstein / Fermi-Dirac -> Gamma*PolyLog), a Frullani pre-pass\n"
        "                          (f(a x)-f(b x))/x -> (f(0)-f(Inf)) Log[b/a], and a Log[x]^k weight;\n"
        "                          strip-gated, yielding a ConditionalExpression when\n"
        "                          Assumptions do not prove convergence; after NewtonLeibniz under Automatic\n"
        "Method -> {\"DerivativeDivides\", \"Substitution\" -> u} pins the kernel u(x),\n"
        "trialing only that substitution.\n"
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

    /* Chebychev binomial differential: Integrate`ChebychevAlgebraic. */
    integrate_chebychev_init();

    /* Goursat pseudo-elliptic / cube-/fourth-root: Integrate`GoursatAlgebraic. */
    integrate_goursat_init();

    /* Continuous Weierstrass substitution (Jeffrey & Rich 1994):
     * Integrate`Weierstrass. */
    integrate_jeffrey_init();

    /* Definite integration by the fundamental theorem of calculus:
     * Integrate`NewtonLeibniz and the pole detector Integrate`SingularPoints. */
    integrate_newton_leibniz_init();

    /* Origin-symmetry reduction: Integrate`Symmetry (odd -> 0, even -> 2 half).
     * Tried first in the definite dispatch. */
    integrate_symmetry_init();

    /* Euler-Beta reductions: Integrate`Beta ([0,1]) and Integrate`TrigPower
     * (Sin^m Cos^n over a canonical trig interval). */
    integrate_beta_init();

    /* Complex line / contour integration: Integrate`LineIntegral and the
     * on-path singularity detector Integrate`PathSingularPoints. */
    integrate_line_init();

    /* Definite integration by the residue theorem: Integrate`ContourResidue.
     * Engaged before Newton-Leibniz for a single real spec under Automatic. */
    integrate_residue_init();

    /* Definite integration by differentiation under the integral sign
     * (Feynman's trick): Integrate`DiffUnderInt.  Last resort in the definite
     * cascade, after residue and Newton-Leibniz. */
    integrate_diffunderint_init();
    integrate_ramanujan_init();

    /* Initialise the parallel-Risch / Risch-Norman heuristic
     * (Bronstein's pmint).  Provides `Integrate`RischNorman[f, x]`,
     * the fall-through for transcendental integrands. */
    intrischnorman_init();

    /* Recursive transcendental Risch integrator:
     * Integrate`RischTranscendental.  Correct by construction; inserted into the
     * Automatic cascade after the parallel-Risch RischNorman. */
    integrate_risch_transcendental_init();

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

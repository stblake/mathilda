/*
 * irrpolyq.c -- IrreduciblePolynomialQ[poly, opts].
 *
 * Always returns True or False on a structurally valid call.
 * Wrong arg count emits `IrreduciblePolynomialQ::argx` and returns NULL.
 * Malformed options emit `IrreduciblePolynomialQ::nonopt` and return NULL.
 *
 * Algorithm:
 *   1. Resolve the factoring field from (GaussianIntegers, Extension)
 *      and a Complex-coefficient sniff of poly.  Precedence:
 *        Extension -> All           : absolute irreducibility (see below).
 *        Extension -> α | {α_i}     : Q(α) / compositum.
 *        Extension -> Automatic     : extension_autodetect on poly.
 *        GaussianIntegers -> True   : Q(i).
 *        complex coef in poly       : Q(i).
 *        otherwise                  : Q.
 *
 *   2. Reject non-polynomial inputs (Sin[x], free symbol-only, ...).
 *      Pure constants (no polynomial variable) -> False (not irreducible).
 *
 *   3. Factor poly using the resolved field, then count non-constant
 *      factors with multiplicity:
 *        0     -> False (constant)
 *        1     -> True
 *        >= 2  -> False
 *
 *   Extension -> All:
 *     Univariate degree-1 polynomials are absolutely irreducible; everything
 *     else of degree >= 2 univariate splits into linear factors over C, so
 *     False.  For multivariate inputs, we approximate absolute irreducibility
 *     by factoring over Q(i) (catches conjugate-pair factorisations like
 *     x^2 + y^2 = (x+iy)(x-iy)).  This is incomplete -- it does not detect
 *     reducibility over Q(sqrt(d)) for d > 0 -- but covers the headline
 *     Wolfram examples and never produces a false "irreducible" for the
 *     Q(i) cases that motivate the option in practice.
 *
 *   IrreduciblePolynomialQ is registered with ATTR_LISTABLE so List inputs
 *   thread element-wise via the evaluator before this builtin runs; the
 *   code here only ever sees a scalar first argument.
 */

#include "irrpolyq.h"
#include "facpoly.h"
#include "qafactor.h"
#include "poly.h"
#include "eval.h"
#include "expand.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "internal.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===================================================================== */
/* Small utilities                                                       */
/* ===================================================================== */

static bool is_sym_eq(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

static bool is_rule_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol == SYM_Rule ||
            e->data.function.head->data.symbol == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

/* True iff `e` (or any descendant) has the `Complex` head.  Used to flip
 * into Gaussian mode when the user's input carries an explicit complex
 * coefficient (e.g. `x^2 + 2 I x - 1`).  Cheap recursive walk; the input
 * is already evaluated. */
static bool contains_complex_head(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_complex_head(e->data.function.args[i])) return true;
    }
    return false;
}

/* Rewrite every `Complex[a, b]` leaf in `e` into `a + b * I` so that the
 * algebraic-factoring path (which lifts polynomials over Q(I) by
 * substituting the surface symbol `I` -> α) can see the imaginary part
 * as a coefficient of `I` rather than an atomic Complex literal.
 *
 * CRITICAL: do NOT call `evaluate` on the result.  Mathilda's evaluator
 * canonicalises `Times[b, I]` straight back to `Complex[0, b]`, undoing
 * the lift.  The qa-factoring path tolerates un-evaluated Plus/Times
 * shapes -- its first step is `expr_subst(I -> alpha_internal)` followed
 * by `expr_expand`, and once `I` has been substituted away the Complex-
 * canonicalisation has nothing to latch onto.  Caller owns the returned
 * Expr. */
static Expr* expand_complex_to_i(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex &&
        e->data.function.arg_count == 2) {
        Expr* re = expand_complex_to_i(e->data.function.args[0]);
        Expr* im = expand_complex_to_i(e->data.function.args[1]);
        Expr** tm = malloc(sizeof(Expr*) * 2);
        tm[0] = im;
        tm[1] = expr_new_symbol("I");
        Expr* mul = expr_new_function(expr_new_symbol("Times"), tm, 2);
        free(tm);
        Expr** pl = malloc(sizeof(Expr*) * 2);
        pl[0] = re;
        pl[1] = mul;
        return expr_new_function(expr_new_symbol("Plus"), pl, 2);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = expand_complex_to_i(e->data.function.args[i]);
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* out = expr_new_function(head, new_args, n);
    free(new_args);
    return out;
}

/* ===================================================================== */
/* Diagnostics                                                           */
/* ===================================================================== */

static Expr* irrpolyq_emit_argx(size_t argc) {
    fprintf(stderr,
            "IrreduciblePolynomialQ::argx: IrreduciblePolynomialQ called with "
            "%zu arguments; 1 argument is expected.\n", argc);
    return NULL;
}

static Expr* irrpolyq_emit_nonopt(Expr* bad, size_t pos, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "IrreduciblePolynomialQ::nonopt: Options expected (instead of %s) "
            "beyond position %zu in %s. An option must be a rule or a list "
            "of rules.\n",
            bad_str ? bad_str : "?", pos, call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* ===================================================================== */
/* Count non-constant factors in a Factor result, with multiplicity.     */
/*                                                                       */
/* Examples (vars = {x}):                                                */
/*   x^2 + 1                       -> 1   (irreducible)                  */
/*   Times[-1+x, 1+x]              -> 2   (reducible)                    */
/*   Power[-1+x, 2]                -> 2   (one factor, mult 2)           */
/*   Times[2, -1+x, 1+x]           -> 2   (the 2 is a constant unit)     */
/*   5                             -> 0   (constant -- not irreducible)  */
/* ===================================================================== */

/* True when `e` is a constant in `vars` (i.e. a unit of Q[vars] when
 * counted with the convention above).  Numeric heads (Integer, BigInt,
 * Real, Rational, Complex) count as constants, as does any symbolic
 * expression that mentions no variable from `vars`. */
static bool factor_is_constant(const Expr* e, Expr** vars, size_t v_count) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
        e->type == EXPR_REAL || e->type == EXPR_STRING) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational || h == SYM_Complex) return true;
    }
    for (size_t i = 0; i < v_count; i++) {
        if (contains_any_symbol_from((Expr*)e, vars[i])) return false;
    }
    return true;
}

static int count_nonconstant_factors(const Expr* f, Expr** vars, size_t v_count) {
    if (!f) return 0;
    if (factor_is_constant(f, vars, v_count)) return 0;

    if (f->type == EXPR_FUNCTION &&
        f->data.function.head->type == EXPR_SYMBOL) {
        const char* h = f->data.function.head->data.symbol;
        if (h == SYM_Times) {
            int total = 0;
            for (size_t i = 0; i < f->data.function.arg_count; i++) {
                total += count_nonconstant_factors(f->data.function.args[i],
                                                   vars, v_count);
            }
            return total;
        }
        if (h == SYM_Power && f->data.function.arg_count == 2) {
            Expr* base = f->data.function.args[0];
            Expr* exp  = f->data.function.args[1];
            if (exp->type == EXPR_INTEGER && exp->data.integer > 0 &&
                !factor_is_constant(base, vars, v_count)) {
                int sub = count_nonconstant_factors(base, vars, v_count);
                /* sub == 0 should not happen for !factor_is_constant base,
                 * but be defensive: treat the Power as a single opaque
                 * non-constant factor. */
                if (sub == 0) sub = 1;
                return sub * (int)exp->data.integer;
            }
        }
    }
    /* Any non-constant atom or unrecognised function with a variable in
     * it is one factor (Plus, Power with non-int exp, opaque heads). */
    return 1;
}

/* ===================================================================== */
/* Factor wrappers                                                       */
/* ===================================================================== */

/* True iff `e` is an "atomic algebraic constant" of the form:
 *   Sqrt[int]                            (Mathilda canonicalises to Power[int, 1/2])
 *   Power[int, Rational[p, q]] with q != 1, |int| >= 2.
 * These shapes are exactly what Mathilda's multivariate Factor implicitly
 * promotes into the auto-extension generator, even with `Extension -> None`. */
static bool is_atomic_algebraic_const(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (!(base->type == EXPR_INTEGER || base->type == EXPR_BIGINT)) return false;
        if (exp->type != EXPR_FUNCTION) return false;
        if (exp->data.function.head->type != EXPR_SYMBOL) return false;
        if (exp->data.function.head->data.symbol != SYM_Rational) return false;
        if (exp->data.function.arg_count != 2) return false;
        Expr* den = exp->data.function.args[1];
        return (den->type == EXPR_INTEGER && den->data.integer > 1);
    }
    if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        return (a->type == EXPR_INTEGER || a->type == EXPR_BIGINT);
    }
    return false;
}

/* Replace every atomic-algebraic-constant subterm of `e` (Sqrt[int],
 * Power[int, p/q]) with a fresh placeholder symbol, ensuring distinct
 * algebraic constants get distinct placeholders.  The result is an
 * un-evaluated Expr where the algebraic constants now look like free
 * symbols to downstream polynomial code -- which is what we want when
 * the user has asked for *no* extension (Mathematica's Extension -> None
 * default) but Mathilda's multivariate Factor would otherwise auto-extend
 * and report `(x + Sqrt[3] y)^2` as a factorisation of
 * `x^2 + 2 Sqrt[3] x y + 3 y^2`.
 *
 * Counter `*next_id` is shared across the walk so we never collide
 * two distinct algebraic constants on the same placeholder.  Distinct-
 * but-equal occurrences (multiple Sqrt[3] in one polynomial) all map to
 * the same placeholder, which preserves the factorisation count. */
static Expr* freeze_alg_constants_walk(const Expr* e,
                                       Expr*** subterms_ptr,
                                       char*** placeholders_ptr,
                                       size_t* count, size_t* cap,
                                       int* next_id) {
    if (!e) return NULL;
    if (is_atomic_algebraic_const(e)) {
        /* See if we already have a placeholder for an equal subterm. */
        for (size_t i = 0; i < *count; i++) {
            if (expr_eq((Expr*)e, (*subterms_ptr)[i])) {
                return expr_new_symbol((*placeholders_ptr)[i]);
            }
        }
        /* Allocate a new placeholder. */
        if (*count == *cap) {
            *cap = (*cap == 0) ? 4 : (*cap * 2);
            *subterms_ptr     = realloc(*subterms_ptr,     sizeof(Expr*) * (*cap));
            *placeholders_ptr = realloc(*placeholders_ptr, sizeof(char*) * (*cap));
        }
        char buf[32];
        snprintf(buf, sizeof buf, "$irrpoly$alg$%d$", (*next_id)++);
        char* dup = malloc(strlen(buf) + 1);
        strcpy(dup, buf);
        (*subterms_ptr)[*count]     = expr_copy((Expr*)e);
        (*placeholders_ptr)[*count] = dup;
        (*count)++;
        return expr_new_symbol(dup);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = freeze_alg_constants_walk(e->data.function.args[i],
                                                subterms_ptr, placeholders_ptr,
                                                count, cap, next_id);
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* out = expr_new_function(head, new_args, n);
    free(new_args);
    return out;
}

/* Top-level wrapper for `freeze_alg_constants_walk`.  Frees the table of
 * (subterm, placeholder) pairs after the walk -- the substituted polynomial
 * keeps its own copies of the placeholder symbol names, so the table is
 * disposable. */
static Expr* freeze_algebraic_constants(const Expr* e) {
    Expr** subterms = NULL;
    char** placeholders = NULL;
    size_t count = 0, cap = 0;
    int next_id = 0;
    Expr* out = freeze_alg_constants_walk(e, &subterms, &placeholders,
                                          &count, &cap, &next_id);
    for (size_t i = 0; i < count; i++) {
        expr_free(subterms[i]);
        free(placeholders[i]);
    }
    free(subterms);
    free(placeholders);
    return out;
}

/* Call `Factor[poly]` -- plain rational-coefficient factoring.  We freeze
 * any Sqrt[int] / Power[int, p/q] atoms first so Mathilda's multivariate
 * Factor doesn't silently apply Extension -> Automatic and report a
 * factorisation in terms of the algebraic generator.  See
 * `freeze_algebraic_constants` for the rationale. */
static Expr* factor_over_q(const Expr* poly) {
    Expr* frozen = freeze_algebraic_constants(poly);
    Expr** args = malloc(sizeof(Expr*) * 1);
    args[0] = frozen;  /* hands ownership over */
    Expr* call = expr_new_function(expr_new_symbol("Factor"), args, 1);
    free(args);
    return eval_and_free(call);
}

/* Call `Factor[poly, Extension -> alpha]`.  `alpha` is borrowed (copied
 * here).  Returns a fresh Expr owned by the caller. */
static Expr* factor_over_extension(const Expr* poly, const Expr* alpha) {
    Expr** rule_args = malloc(sizeof(Expr*) * 2);
    rule_args[0] = expr_new_symbol("Extension");
    rule_args[1] = expr_copy((Expr*)alpha);
    Expr* rule = expr_new_function(expr_new_symbol("Rule"), rule_args, 2);
    free(rule_args);

    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy((Expr*)poly);
    args[1] = rule;
    Expr* call = expr_new_function(expr_new_symbol("Factor"), args, 2);
    free(args);
    return eval_and_free(call);
}

/* ===================================================================== */
/* Auto-collect polynomial variables (bare symbols only).                */
/*                                                                       */
/* Returns NULL on success with `*out_count == 0` for a constant input.  */
/* Returns a malloc'd Expr** array of borrowed pointers into `e` on      */
/* success otherwise.  Caller must `free()` the array; entries are not   */
/* owned (collect_variables makes copies; we keep those copies in the    */
/* array and free them via the `for + free` loop in the caller).         */
/*                                                                       */
/* When `e` references a non-symbol "variable" (e.g. Sin[x], Sqrt[2]),   */
/* those entries are dropped from the returned list.  If after filtering */
/* the list is empty but the input is non-numeric, sets `*out_nonpoly`   */
/* so the caller can short-circuit to False.                             */
/* ===================================================================== */
static Expr** auto_collect_vars(const Expr* e, size_t* out_count,
                                bool* out_nonpoly) {
    *out_nonpoly = false;
    size_t cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * cap);
    size_t vc = 0;
    collect_variables((Expr*)e, &vars, &vc, &cap);

    size_t kept = 0;
    bool dropped_any = false;
    for (size_t i = 0; i < vc; i++) {
        if (vars[i]->type == EXPR_SYMBOL) {
            vars[kept++] = vars[i];
        } else {
            expr_free(vars[i]);
            dropped_any = true;
        }
    }
    vc = kept;

    /* If we dropped non-symbol "variables" and nothing usable remains,
     * the input was a non-polynomial expression (e.g. Sin[x]).  Signal
     * that so the caller returns False without proceeding to Factor. */
    if (vc == 0 && dropped_any) *out_nonpoly = true;
    *out_count = vc;
    return vars;
}

/* ===================================================================== */
/* Multivariate Hilbert-irreducibility specialisation probe              */
/*                                                                       */
/* Mathilda's Factor[poly, Extension -> α] currently routes only the     */
/* single-polynomial-variable case through the qa-factoring path         */
/* (facpoly_factor_builtin.inc:100); multivariate inputs silently fall   */
/* back to plain Q-factoring, so factorisations that *require* the       */
/* extension to surface are missed.  Canonical example:                  */
/*                                                                       */
/*   x^4 - 3 y^2 = (x^2 - Sqrt[3] y) (x^2 + Sqrt[3] y)                   */
/*                                                                       */
/* stays unfactored, and a naive count of non-constant factors reports   */
/* the polynomial as irreducible.                                        */
/*                                                                       */
/* The probe fixes a "primary" indeterminate, specialises every other    */
/* variable to each integer c in a small probe set, factors the          */
/* resulting univariate over the extension (which the existing path     */
/* handles correctly), and counts non-constant factors.  By Hilbert's    */
/* irreducibility theorem, an irreducible multivariate p stays           */
/* irreducible under almost every integer specialisation; conversely, a  */
/* reducible p factors under almost every specialisation.  Requiring     */
/* unanimous "reducible" verdicts across the probe set minimises the     */
/* false-False risk (where an irreducible p coincidentally factors at    */
/* every probed c) at the cost of missing reducible polynomials whose    */
/* factors degenerate at every probed c.                                 */
/* ===================================================================== */

/* Recursive structural substitution: replace every sub-expression
 * structurally equal to `target` with a fresh copy of `repl`.  Returns
 * a freshly-allocated tree; caller owns. */
static Expr* irrpolyq_subst(const Expr* e, const Expr* target, const Expr* repl) {
    if (!e) return NULL;
    if (expr_eq((Expr*)e, (Expr*)target)) return expr_copy((Expr*)repl);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    Expr* new_head = irrpolyq_subst(e->data.function.head, target, repl);
    size_t n = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = irrpolyq_subst(e->data.function.args[i], target, repl);
    }
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Returns 1 iff every valid specialisation factors into >= 2 non-constant
 * factors over Q(α); 0 otherwise (including no valid specialisation).
 * `alpha_expr` is borrowed.  Picks the variable of maximum degree as the
 * surviving univariate indeterminate so the probed polynomial carries the
 * most signal. */
static int irr_multivariate_specialize_probe(const Expr* poly,
                                             Expr** vars, size_t v_count,
                                             const Expr* alpha_expr) {
    if (v_count < 2 || !alpha_expr) return 0;

    int best_deg = -1;
    size_t primary = 0;
    for (size_t i = 0; i < v_count; i++) {
        int d = get_degree_poly((Expr*)poly, vars[i]);
        if (d > best_deg) { best_deg = d; primary = i; }
    }
    if (best_deg < 2) return 0;

    static const int64_t cs[] = {2, 3, 5};
    const size_t n_specs = sizeof(cs) / sizeof(cs[0]);

    int agree_reducible = 0;
    int valid = 0;
    Expr* one_var[1] = { vars[primary] };

    for (size_t k = 0; k < n_specs; k++) {
        int64_t c = cs[k];

        Expr* cur = expr_copy((Expr*)poly);
        for (size_t i = 0; i < v_count; i++) {
            if (i == primary) continue;
            Expr* val = expr_new_integer(c);
            Expr* next = irrpolyq_subst(cur, vars[i], val);
            expr_free(val);
            expr_free(cur);
            cur = next;
        }
        Expr* spec_poly = evaluate(cur);
        expr_free(cur);

        /* A specialisation that drops the degree in the primary variable
         * doesn't reflect the original polynomial; skip it. */
        if (get_degree_poly(spec_poly, vars[primary]) != best_deg) {
            expr_free(spec_poly);
            continue;
        }

        Expr* factored = factor_over_extension(spec_poly, alpha_expr);
        expr_free(spec_poly);
        if (!factored) continue;
        int nonconst = count_nonconstant_factors(factored, one_var, 1);
        expr_free(factored);

        valid++;
        if (nonconst >= 2) agree_reducible++;
    }

    return (valid > 0 && agree_reducible == valid) ? 1 : 0;
}

/* ===================================================================== */
/* Top-level dispatcher                                                  */
/* ===================================================================== */

typedef enum {
    EXT_MODE_NONE = 0,        /* Q (no extension) */
    EXT_MODE_AUTOMATIC,       /* Extension -> Automatic */
    EXT_MODE_ALL,             /* Extension -> All (absolute) */
    EXT_MODE_EXPLICIT         /* Extension -> α  or  Extension -> {α_i} */
} ExtMode;

static bool irr_dispatch(Expr* poly, ExtMode mode, Expr* alpha_expr,
                         bool gaussian) {
    /* Non-polynomial inputs (Sin[x], etc.) and pure constants both
     * resolve to False -- a constant is not an irreducible polynomial. */
    bool nonpoly = false;
    size_t v_count = 0;
    Expr** vars = auto_collect_vars(poly, &v_count, &nonpoly);
    if (nonpoly) {
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        return false;
    }
    if (v_count == 0) {
        /* No polynomial variable -- constant or non-polynomial atom.
         * False per Mathematica: constants aren't irreducible. */
        free(vars);
        return false;
    }

    if (!is_polynomial(poly, vars, v_count)) {
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        return false;
    }

    /* Extension -> All: absolute irreducibility.
     *   Univariate: only degree-1 polynomials are absolutely irreducible.
     *   Multivariate: approximate via Factor[poly, Extension -> I].  This
     *   catches conjugate-pair factorisations (x^2+y^2 = (x+iy)(x-iy))
     *   without claiming completeness over arbitrary algebraic extensions. */
    Expr* factored = NULL;
    if (mode == EXT_MODE_ALL) {
        if (v_count == 1) {
            int deg = get_degree_poly(poly, vars[0]);
            for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
            free(vars);
            return (deg == 1);
        }
        Expr* i_sym = expr_new_symbol("I");
        factored = factor_over_extension(poly, i_sym);
        expr_free(i_sym);
    } else if (mode == EXT_MODE_EXPLICIT && alpha_expr) {
        factored = factor_over_extension(poly, alpha_expr);
    } else if (mode == EXT_MODE_AUTOMATIC) {
        Expr* auto_sym = expr_new_symbol("Automatic");
        factored = factor_over_extension(poly, auto_sym);
        expr_free(auto_sym);
    } else if (gaussian) {
        /* Lift Complex[a, b] leaves into a + b*I so the algebraic-
         * factoring path can substitute the surface `I` -> α; without
         * this, qa_expr_to_qaupoly sees an atomic Complex literal as
         * an opaque coefficient and falls back to leaving poly unchanged.
         * The lifted form is intentionally un-evaluated (see
         * `expand_complex_to_i`'s docstring for why). */
        Expr* lifted = expand_complex_to_i(poly);
        Expr* i_sym  = expr_new_symbol("I");
        if (v_count == 1) {
            /* Direct call into qa_factor_with_extension bypasses the
             * public Factor entry point and, with it, the Together /
             * Numerator preprocessing that would re-canonicalise our
             * un-evaluated `Times[b, I]` back into `Complex[0, b]`. */
            factored = qa_factor_with_extension(lifted, i_sym, vars[0]);
        }
        if (!factored) {
            /* Multivariate Q(i) factoring is not yet supported by the
             * qa machinery; the best-effort fallback factors over Q.
             * This is incomplete -- it will miss factorisations that
             * require I, e.g. x^2 + y^2 = (x + i y)(x - i y) -- but
             * preserves the contract that this builtin always returns
             * a boolean rather than leaving the call unevaluated. */
            factored = factor_over_q(poly);
        }
        expr_free(lifted);
        expr_free(i_sym);
    } else {
        factored = factor_over_q(poly);
    }

    int nonconst = count_nonconstant_factors(factored, vars, v_count);
    bool result = (nonconst == 1);

    if (factored) expr_free(factored);

    /* Multivariate gap: Factor[poly, Extension -> α] only applies the
     * extension to univariate inputs, so multivariate cases that need
     * the extension to expose a factorisation slip through the cheap
     * path.  When the cheap path returned "irreducible" on a multivariate
     * input that names an algebraic extension, run a Hilbert-style
     * specialisation probe: it can only flip the verdict from True to
     * False, never the other way (the existing factor result is
     * authoritative for direct decompositions like x*y). */
    if (result && v_count > 1) {
        const Expr* probe_alpha = NULL;
        Expr* probe_alpha_owned = NULL;
        if (mode == EXT_MODE_EXPLICIT && alpha_expr) {
            probe_alpha = alpha_expr;
        } else if (mode == EXT_MODE_ALL) {
            probe_alpha_owned = expr_new_symbol("I");
            probe_alpha = probe_alpha_owned;
        } else if (gaussian && mode == EXT_MODE_NONE) {
            probe_alpha_owned = expr_new_symbol("I");
            probe_alpha = probe_alpha_owned;
        }
        if (probe_alpha &&
            irr_multivariate_specialize_probe(poly, vars, v_count, probe_alpha)) {
            result = false;
        }
        if (probe_alpha_owned) expr_free(probe_alpha_owned);
    }

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);
    return result;
}

/* ===================================================================== */
/* Entry point                                                           */
/* ===================================================================== */

Expr* builtin_irreduciblepolynomialq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return irrpolyq_emit_argx(0);

    Expr* poly = res->data.function.args[0];

    /* Parse options at positions [1..argc).  Recognised:
     *   GaussianIntegers -> True|False|Automatic
     *   Extension        -> None | All | Automatic | α | {α_i}
     * Anything else (or unknown option name) marks `last_bad`; the
     * diagnostic reports the LAST offending arg (Mathematica's chosen
     * phrasing) only AFTER the scan completes. */
    int gaussian_setting = 0;   /* 0 = Automatic, 1 = True, -1 = False */
    ExtMode ext_mode = EXT_MODE_NONE;
    Expr* alpha_expr = NULL;    /* borrowed when EXT_MODE_EXPLICIT */
    Expr* last_bad = NULL;

    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!is_rule_head(opt) || opt->data.function.args[0]->type != EXPR_SYMBOL) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];

        if (name == SYM_GaussianIntegers || strcmp(name, "GaussianIntegers") == 0) {
            if (is_sym_eq(val, "True"))       gaussian_setting = 1;
            else if (is_sym_eq(val, "False")) gaussian_setting = -1;
            else if (is_sym_eq(val, "Automatic")) gaussian_setting = 0;
            else last_bad = opt;
        } else if (name == SYM_Extension || strcmp(name, "Extension") == 0) {
            if (is_sym_eq(val, "None")) {
                ext_mode = EXT_MODE_NONE;
                alpha_expr = NULL;
            } else if (is_sym_eq(val, "All")) {
                ext_mode = EXT_MODE_ALL;
                alpha_expr = NULL;
            } else if (is_sym_eq(val, "Automatic")) {
                ext_mode = EXT_MODE_AUTOMATIC;
                alpha_expr = NULL;
            } else {
                ext_mode = EXT_MODE_EXPLICIT;
                alpha_expr = val;  /* borrowed pointer into res */
            }
        } else {
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        return irrpolyq_emit_nonopt(last_bad, 1, res);
    }

    /* Resolve GaussianIntegers default and complex-coefficient sniff.
     * GaussianIntegers -> True wins outright; -> False forces off; the
     * default (Automatic) flips on iff the polynomial mentions Complex. */
    bool gaussian;
    if (gaussian_setting == 1) gaussian = true;
    else if (gaussian_setting == -1) gaussian = false;
    else gaussian = contains_complex_head(poly);

    /* Extension dominates GaussianIntegers when both are set non-default:
     * Mathematica routes irreducibility testing through the explicit
     * extension first, only falling back to the Gaussian-rationals path
     * when no Extension option is given.  This matches the "If any
     * coefficients in poly are complex numbers, irreducibility testing
     * is done over the Gaussian rationals" rule -- a no-op when a richer
     * extension already covers I (e.g. Extension -> All in C, or an
     * algebraic-tower extension that already includes I). */
    bool result = irr_dispatch(poly, ext_mode, alpha_expr, gaussian);
    return expr_new_symbol(result ? "True" : "False");
}

/* ===================================================================== */
/* Init                                                                  */
/* ===================================================================== */

void irrpolyq_init(void) {
    symtab_add_builtin("IrreduciblePolynomialQ", builtin_irreduciblepolynomialq);
    symtab_get_def("IrreduciblePolynomialQ")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
    /* Docstring lives in info.c alongside the other *Q predicates. */
}

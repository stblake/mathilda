/*
 * nsolve.c — NSolve[expr, vars, dom, prec, opts]
 *
 * Numerical equation solver.  NSolve returns approximate solutions of an
 * equation or system of equations as a list of replacement-rule lists:
 *
 *     {}                          no solutions
 *     {{x -> r1}, {x -> r2}, ...} one rule list per solution
 *     {{}}                        universal solution (every point satisfies)
 *
 * Strategy (two specialists, matching the Wolfram Language's "Symbolic" idea):
 *
 *   1. Univariate polynomial equations  ->  NRoots.
 *      When the input reduces to a single polynomial equation lhs == rhs in a
 *      single variable, NSolve calls NRoots (the state-of-the-art Aberth /
 *      companion-matrix / Jenkins–Traub engine) and repackages its disjunction
 *      var==r1 || var==r2 || ...  as the rule-list form.  This covers integer,
 *      real, and complex coefficients, multiple roots (repeated by
 *      multiplicity), machine and arbitrary working precision, and the Reals
 *      domain (by discarding the complex roots).
 *
 *   2. Everything else  ->  Solve, then numericalise.
 *      Linear systems, radical and inverse-function equations, etc. are solved
 *      symbolically by Solve and the exact result is rounded to the requested
 *      working precision.  This is the "Symbolic" method.  Inputs Solve cannot
 *      handle (e.g. genuine nonlinear polynomial systems) leave NSolve
 *      unevaluated.
 *
 * Options:  MaxRoots, Method, WorkingPrecision, VerifySolutions, RandomSeeding,
 *           PrecisionGoal, MaxIterations.  (Method and the verification/seeding
 *           options are accepted for compatibility; the polynomial engine is
 *           always the NRoots default.)
 *
 * Positional grammar:  NSolve[expr [, vars [, dom [, prec]]], opts...].
 *   dom  in {Reals, Complexes, Integers}; default Complexes.
 *   prec a number giving the working precision in decimal digits.
 *
 * Memory contract (builtin): takes ownership of `res`; returns a fresh Expr* on
 * success (the evaluator frees `res`) or NULL to leave NSolve unevaluated.
 */

#include "nroots.h"          /* builtin_nroots */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numeric.h"         /* numeric_digits_to_bits, NumericSpec */
#include "common.h"          /* common_numericalize_result, rationalize_input */
#include "arithmetic.h"      /* is_number, is_rational */
#include "poly/poly.h"       /* is_polynomial */
#include "nsolve_system.h"   /* nsolve_polynomial_system */
#include "findroot.h"        /* builtin_findroot */

#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Largest univariate polynomial degree NSolve will hand to the (simultaneous)
 * NRoots engine.  Beyond this, computing every root is infeasible and would
 * hang, so NSolve leaves the input unevaluated with an NSolve::deg message. */
#define NSOLVE_MAX_POLY_DEGREE 10000

/* ------------------------------------------------------------------ *
 *  Small Expr helpers
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* A value counts as real unless it is wrapped in a Complex[_, _] head. */
static bool value_is_real(const Expr* v) {
    return !is_head(v, SYM_Complex);
}

/* True iff `e` is a domain specifier symbol (Reals / Complexes / Integers). */
static bool is_domain_symbol(const Expr* e) {
    return e && e->type == EXPR_SYMBOL
        && (e->data.symbol.name == SYM_Reals
            || e->data.symbol.name == SYM_Complexes
            || e->data.symbol.name == SYM_Integers);
}

/* Read an exact numeric Expr (Integer / Real / BigInt / Rational / MPFR) as a
 * double.  Returns false for non-numeric arguments. */
static bool expr_to_double(const Expr* e, double* out) {
    int64_t n, d;
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational((Expr*)e, &n, &d) && d != 0) {
                *out = (double)n / (double)d; return true;
            }
            return false;
        default: return false;
    }
}

/* True iff `e` is a valid variable specification: a non-numeric symbol, or a
 * non-empty List of such symbols.  Domain symbols are excluded so a trailing
 * Reals/Complexes is never mistaken for the variable list. */
static bool is_var_spec(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return !is_domain_symbol(e);
    if (is_head(e, SYM_List)) {
        size_t n = e->data.function.arg_count;
        if (n == 0) return false;
        for (size_t i = 0; i < n; i++)
            if (e->data.function.args[i]->type != EXPR_SYMBOL) return false;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Variable collection (for the NSolve[expr] auto-variable form)
 * ------------------------------------------------------------------ */

static bool is_constant_symbol(const char* s) {
    return s == SYM_Pi || s == SYM_E || s == SYM_I || s == SYM_Degree
        || s == SYM_EulerGamma || s == SYM_GoldenRatio || s == SYM_Catalan
        || s == SYM_Glaisher || s == SYM_Khinchin || s == SYM_Indeterminate
        || s == SYM_Infinity || s == SYM_ComplexInfinity
        || s == SYM_True || s == SYM_False;
}

/* Recursively gather the distinct non-constant symbols appearing anywhere in
 * `e` (descending through every function's arguments but not its head), in
 * order of first appearance.  Used only when the caller omits the variable
 * specification. */
static void collect_symbols(const Expr* e, Expr*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_constant_symbol(e->data.symbol.name)) return;
        for (size_t i = 0; i < *n; i++)
            if ((*out)[i]->data.symbol.name == e->data.symbol.name) return;
        if (*n == *cap) { *cap *= 2; *out = realloc(*out, sizeof(Expr*) * (*cap)); }
        (*out)[(*n)++] = expr_new_symbol(e->data.symbol.name);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_symbols(e->data.function.args[i], out, n, cap);
    }
}

/* Build a List of the auto-collected variables; NULL if none were found. */
static Expr* auto_variable_list(const Expr* expr) {
    size_t n = 0, cap = 8;
    Expr** vars = malloc(sizeof(Expr*) * cap);
    collect_symbols(expr, &vars, &n, &cap);
    if (n == 0) { free(vars); return NULL; }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), vars, n);
    free(vars);
    return list;
}

/* ------------------------------------------------------------------ *
 *  Options
 * ------------------------------------------------------------------ */

/* System-solver method selection. */
typedef enum {
    NSM_AUTOMATIC = 0,   /* eigenvalue (endomorphism) engine               */
    NSM_HOMOTOPY,        /* routed to the eigenvalue engine for now        */
    NSM_SYMBOLIC         /* skip the numeric system engine; use Solve      */
} NSolveMethod;

typedef struct {
    long         max_roots;    /* -1 = all (Automatic / Infinity)              */
    double       prec_digits;  /* -1 = machine (Automatic / MachinePrecision)  */
    NSolveMethod method;       /* polynomial-system method                     */
    int          verify;       /* -1 auto (=on), 0 off, 1 on                   */
    unsigned long seed;        /* RandomSeeding for the generic form           */
} NSolveOpts;

static bool nsolve_is_known_option(const char* s) {
    return s == SYM_MaxRoots || s == SYM_Method || s == SYM_WorkingPrecision
        || s == SYM_VerifySolutions || s == SYM_RandomSeeding
        || s == SYM_PrecisionGoal || s == SYM_MaxIterations;
}

static bool nsolve_is_option_arg(const Expr* e) {
    if (!is_head(e, SYM_Rule) && !is_head(e, SYM_RuleDelayed)) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nsolve_is_known_option(lhs->data.symbol.name);
}

/* Apply one option rule.  Unknown values are tolerated (left at default). */
static void nsolve_apply_option(const Expr* rule, NSolveOpts* o) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_MaxRoots) {
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol.name == SYM_Automatic || rhs->data.symbol.name == SYM_Infinity)) {
            o->max_roots = -1;
        } else if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 0) {
            o->max_roots = (long)rhs->data.integer;
        }
        return;
    }
    if (name == SYM_WorkingPrecision || name == SYM_PrecisionGoal) {
        double v;
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol.name == SYM_Automatic
                || rhs->data.symbol.name == SYM_MachinePrecision)) {
            o->prec_digits = -1.0;
        } else if (expr_to_double(rhs, &v) && v > 0.0) {
            o->prec_digits = v;
        }
        return;
    }
    if (name == SYM_Method) {
        if (rhs->type == EXPR_STRING) {
            const char* s = rhs->data.string;
            if (strcmp(s, "Symbolic") == 0)            o->method = NSM_SYMBOLIC;
            else if (strcmp(s, "Homotopy") == 0)       o->method = NSM_HOMOTOPY;
            else /* EndomorphismMatrix / Monodromy / ... */ o->method = NSM_AUTOMATIC;
        }
        return;
    }
    if (name == SYM_VerifySolutions) {
        if (rhs->type == EXPR_SYMBOL) {
            if (rhs->data.symbol.name == SYM_True)       o->verify = 1;
            else if (rhs->data.symbol.name == SYM_False) o->verify = 0;
            else                                     o->verify = -1; /* Automatic */
        }
        return;
    }
    if (name == SYM_RandomSeeding) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 0)
            o->seed = (unsigned long)rhs->data.integer;
        return;
    }
    /* MaxIterations: accepted, no effect on the polynomial engine. */
}

/* ------------------------------------------------------------------ *
 *  NRoots output  ->  rule-list form
 * ------------------------------------------------------------------ */

/* Wrap a root value as the inner solution list  {var -> val}. */
static Expr* make_rule_list(Expr* var, Expr* val_owned) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(var), val_owned }, 2);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
}

/* Convert NRoots' disjunction (var==r1 || ... , or a bare var==r) into the
 * NSolve list-of-rule-lists, applying the Reals filter and MaxRoots cap.
 * `nr` is borrowed.  Returns a freshly owned List. */
static Expr* rules_from_nroots(Expr* nr, Expr* var,
                               bool reals_only, long max_roots) {
    /* Gather the individual  var == val  equations. */
    Expr** eqs;
    size_t neq;
    Expr* single[1];
    if (is_head(nr, SYM_Or)) {
        eqs = nr->data.function.args;
        neq = nr->data.function.arg_count;
    } else {
        single[0] = nr;
        eqs = single;
        neq = 1;
    }

    Expr** sols = malloc(sizeof(Expr*) * (neq ? neq : 1));
    size_t ns = 0;
    for (size_t i = 0; i < neq; i++) {
        Expr* eq = eqs[i];
        if (!is_head(eq, SYM_Equal) || eq->data.function.arg_count != 2) continue;
        Expr* val = eq->data.function.args[1];
        if (reals_only && !value_is_real(val)) continue;
        if (max_roots >= 0 && (long)ns >= max_roots) break;
        sols[ns++] = make_rule_list(var, expr_copy(val));
    }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), sols, ns);
    free(sols);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Univariate polynomial path
 * ------------------------------------------------------------------ */

/* Residual lhs - rhs of an equation (or the expression itself for a bare
 * polynomial), expanded.  Returns a freshly owned Expr*. */
static Expr* equation_residual(Expr* eqn) {
    Expr* sub;
    if (is_head(eqn, SYM_Equal) && eqn->data.function.arg_count == 2) {
        Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(eqn->data.function.args[1]) }, 2);
        sub = expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){ expr_copy(eqn->data.function.args[0]), neg }, 2);
    } else {
        sub = expr_copy(eqn);
    }
    return eval_and_free(sub);
}

/* Build the equation to hand to NRoots: the input equation unchanged, or
 * Equal[expr, 0] for a bare polynomial expression.  Freshly owned. */
static Expr* as_equation(Expr* expr) {
    if (is_head(expr, SYM_Equal) && expr->data.function.arg_count == 2)
        return expr_copy(expr);
    return expr_new_function(expr_new_symbol(SYM_Equal),
               (Expr*[]){ expr_copy(expr), expr_new_integer(0) }, 2);
}

/* True if `e` mentions the variable symbol `var` anywhere. */
static bool expr_mentions(const Expr* e, const Expr* var) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return var->type == EXPR_SYMBOL && e->data.symbol.name == var->data.symbol.name;
    if (e->type == EXPR_FUNCTION) {
        if (expr_mentions(e->data.function.head, var)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_mentions(e->data.function.args[i], var)) return true;
    }
    return false;
}

/* Cheap structural check (no dense-polynomial materialisation): is there a
 * Power[b, n] with a literal integer exponent n > cap whose base b mentions
 * `var`?  Used to bail on absurd-degree inputs (x^1000000) before the polynomial
 * machinery — which would otherwise allocate a million-entry coefficient array
 * just to discover the degree. */
static bool has_huge_power(const Expr* e, const Expr* var, long cap) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_head(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* ex = e->data.function.args[1];
        bool too_big = false;
        if (ex->type == EXPR_BIGINT) too_big = true;          /* > any long cap */
        else if (ex->type == EXPR_INTEGER && (long)ex->data.integer > cap)
            too_big = true;
        if (too_big && expr_mentions(e->data.function.args[0], var))
            return true;
    }
    if (has_huge_power(e->data.function.head, var, cap)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_huge_power(e->data.function.args[i], var, cap)) return true;
    return false;
}

/* Try the polynomial path.  Returns the solution list on success, or NULL to
 * fall through to the Solve specialist.  `expr` and `var` are borrowed. */
static Expr* nsolve_polynomial(Expr* expr, Expr* var,
                               bool reals_only, const NSolveOpts* o) {
    /* Confirm the residual is a polynomial in `var`. */
    Expr* residual = equation_residual(expr);
    bool is_poly = is_polynomial(residual, &var, 1);
    expr_free(residual);
    if (!is_poly) return NULL;

    /* Assemble  NRoots[eqn, var, PrecisionGoal -> digits]  and call the engine
     * directly (NRoots never frees its argument). */
    Expr* eqn = as_equation(expr);
    Expr* nrexpr;
    if (o->prec_digits > 0.0) {
        Expr* pg = expr_new_function(expr_new_symbol(SYM_Rule),
                       (Expr*[]){ expr_new_symbol(SYM_PrecisionGoal),
                                  expr_new_real(o->prec_digits) }, 2);
        nrexpr = expr_new_function(expr_new_symbol(SYM_NRoots),
                     (Expr*[]){ eqn, expr_copy(var), pg }, 3);
    } else {
        nrexpr = expr_new_function(expr_new_symbol(SYM_NRoots),
                     (Expr*[]){ eqn, expr_copy(var) }, 2);
    }

    Expr* nr = builtin_nroots(nrexpr);
    expr_free(nrexpr);
    if (!nr) return NULL;

    /* Constant equations collapse to True/False inside NRoots. */
    if (nr->type == EXPR_SYMBOL && nr->data.symbol.name == SYM_True) {
        expr_free(nr);
        Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ empty }, 1);
    }
    if (nr->type == EXPR_SYMBOL && nr->data.symbol.name == SYM_False) {
        expr_free(nr);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    Expr* out = rules_from_nroots(nr, var, reals_only, o->max_roots);
    expr_free(nr);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Solve fallback (linear systems, radical / inverse equations, ...)
 * ------------------------------------------------------------------ */

/* Truncate a result List to at most `max_roots` solutions (no-op if -1 or the
 * result is not a plain List of solutions). `lst` is consumed; returns owned. */
static Expr* cap_solutions(Expr* lst, long max_roots) {
    if (max_roots < 0 || !is_head(lst, SYM_List)) return lst;
    size_t n = lst->data.function.arg_count;
    if ((long)n <= max_roots) return lst;
    Expr** keep = malloc(sizeof(Expr*) * (size_t)(max_roots ? max_roots : 1));
    for (long i = 0; i < max_roots; i++)
        keep[i] = expr_copy(lst->data.function.args[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), keep, (size_t)max_roots);
    free(keep);
    expr_free(lst);
    return out;
}

/* If `e` numericalizes to a concrete number, set *mag = |e| and return true. */
static bool nsolve_number_abs(Expr* e, double* mag) {
    Expr* nexpr = expr_new_function(expr_new_symbol(SYM_N),
                      (Expr*[]){ expr_copy(e) }, 1);
    Expr* n = eval_and_free(nexpr);
    double re, im;
    bool ok = false;
    if (expr_to_double(n, &re)) { *mag = fabs(re); ok = true; }
    else if (is_head(n, SYM_Complex) && n->data.function.arg_count == 2 &&
             expr_to_double(n->data.function.args[0], &re) &&
             expr_to_double(n->data.function.args[1], &im)) {
        *mag = hypot(re, im); ok = true;
    }
    expr_free(n);
    return ok;
}

/* True if a single solution `sol` ({var -> val, ...}) is provably extraneous for
 * the equation / system `expr`: substituting it into some equation's residual
 * numericalizes to a value clearly away from zero, or (under `reals_only`) a
 * bound value is a manifestly non-real number.  Solutions that do not
 * numericalize (symbolic, or ConditionalExpression families with free
 * parameters) are kept (return false) — only provably wrong roots are dropped.
 * This filters the extraneous roots that radical substitution (t = x^(1/q),
 * x = t^q) introduces, which Solve cannot always discharge symbolically. */
static bool solution_is_extraneous(Expr* expr, Expr* sol, bool reals_only) {
    if (!is_head(sol, SYM_List)) return false;     /* not a rule list -> keep */

    /* Reals filter: drop if any bound value is a manifestly non-real number. */
    if (reals_only) {
        for (size_t i = 0; i < sol->data.function.arg_count; i++) {
            Expr* rule = sol->data.function.args[i];
            if (is_head(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                Expr* val = rule->data.function.args[1];
                if (is_head(val, SYM_Complex) && val->data.function.arg_count == 2) {
                    double im;
                    if (expr_to_double(val->data.function.args[1], &im) &&
                        fabs(im) > 1e-9)
                        return true;
                }
            }
        }
    }

    /* Residual check: substitute the solution into each equation's residual. */
    Expr** eqs; size_t neq; Expr* single[1];
    if (is_head(expr, SYM_And) || is_head(expr, SYM_List)) {
        eqs = expr->data.function.args;
        neq = expr->data.function.arg_count;
    } else {
        single[0] = expr; eqs = single; neq = 1;
    }
    for (size_t i = 0; i < neq; i++) {
        Expr* resid  = equation_residual(eqs[i]);  /* lhs - rhs (owned) */
        Expr* subbed = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                           (Expr*[]){ resid, expr_copy(sol) }, 2);
        Expr* ev = eval_and_free(subbed);
        double mag;
        bool isnum = nsolve_number_abs(ev, &mag);
        expr_free(ev);
        if (isnum && mag > 1e-6) return true;      /* provably wrong */
    }
    return false;
}

/* Drop provably-extraneous / non-real solutions from a numericalized result
 * List.  `lst` is consumed; returns a freshly owned List. */
static Expr* filter_solutions(Expr* lst, Expr* expr, bool reals_only) {
    if (!is_head(lst, SYM_List)) return lst;
    size_t n = lst->data.function.arg_count;
    Expr** keep = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nk = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = lst->data.function.args[i];
        if (!solution_is_extraneous(expr, sol, reals_only))
            keep[nk++] = expr_copy(sol);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), keep, nk);
    free(keep);
    expr_free(lst);
    return out;
}

/* Solve `expr` for `varlist` over `dom` symbolically, then round the bindings
 * to the requested working precision.  Returns NULL (leaving NSolve
 * unevaluated) when Solve cannot reduce the input. `expr`, `varlist`, `dom`
 * are borrowed. */
static Expr* nsolve_via_solve(Expr* expr, Expr* varlist,
                              Expr* dom, bool want_machine, long bits,
                              long max_roots) {
    Expr* solve_call;
    if (dom) {
        solve_call = expr_new_function(expr_new_symbol(SYM_Solve),
                         (Expr*[]){ expr_copy(expr), expr_copy(varlist),
                                    expr_copy(dom) }, 3);
    } else {
        solve_call = expr_new_function(expr_new_symbol(SYM_Solve),
                         (Expr*[]){ expr_copy(expr), expr_copy(varlist) }, 2);
    }
    Expr* sol = eval_and_free(solve_call);

    /* Solve returned unevaluated -> NSolve cannot handle this input. */
    if (is_head(sol, SYM_Solve)) { expr_free(sol); return NULL; }

    Expr* numeric = common_numericalize_result(sol, want_machine ? 53 : bits);
    expr_free(sol);

    /* Drop roots that radical/inverse substitution introduced but that do not
     * satisfy the original equation, and (under Reals) manifestly complex roots.
     * Solve can leave such extraneous roots when it cannot discharge the
     * substitution branch symbolically. */
    bool reals_only = dom && dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals;
    numeric = filter_solutions(numeric, expr, reals_only);

    return cap_solutions(numeric, max_roots);
}

/* ------------------------------------------------------------------ *
 *  Multivariate polynomial-system path
 * ------------------------------------------------------------------ */

/* True iff `e` is a relational head other than Equal (skipped by the
 * polynomial-system collector). */
static bool is_skipped_relation(const Expr* e) {
    return is_head(e, SYM_Unequal) || is_head(e, SYM_Greater)
        || is_head(e, SYM_GreaterEqual) || is_head(e, SYM_Less)
        || is_head(e, SYM_LessEqual) || is_head(e, SYM_Inequality);
}

/* Recursively collect residual polynomials (lhs - rhs, or a bare expression)
 * from an And / List / Equal system.  Inequalities are skipped. */
static void collect_residuals_rec(Expr* e, Expr*** out, size_t* n, size_t* cap) {
    if (is_head(e, SYM_And) || is_head(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_residuals_rec(e->data.function.args[i], out, n, cap);
        return;
    }
    if (is_skipped_relation(e)) return;

    Expr* resid;
    if (is_head(e, SYM_Equal) && e->data.function.arg_count == 2) {
        Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(e->data.function.args[1]) }, 2);
        resid = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(e->data.function.args[0]), neg }, 2);
    } else {
        resid = expr_copy(e);  /* bare expression == 0 */
    }
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 4; *out = realloc(*out, sizeof(Expr*) * (*cap)); }
    (*out)[(*n)++] = resid;
}

/* Try the numeric polynomial-system engine.  Builds residuals, rationalises
 * inexact coefficients, checks each is polynomial in `vars`, then calls
 * nsolve_polynomial_system.  Returns the solution list on success, NULL to
 * fall through to the Solve specialist.  `expr` and `varlist` are borrowed. */
static Expr* nsolve_system_path(Expr* expr, Expr* varlist, bool reals_only,
                                bool want_machine, long bits,
                                const NSolveOpts* o) {
    Expr** vars = varlist->data.function.args;
    int nvar = (int)varlist->data.function.arg_count;

    Expr** resid = NULL; size_t nres = 0, cap = 0;
    collect_residuals_rec(expr, &resid, &nres, &cap);
    if (nres == 0) { free(resid); return NULL; }

    /* Rationalise inexact coefficients and confirm each residual is a
     * polynomial in the variables. */
    long rbits = want_machine ? 53 : bits;
    bool all_poly = true;
    Expr** polys = (Expr**)malloc(sizeof(Expr*) * nres);
    for (size_t i = 0; i < nres; i++) {
        polys[i] = common_rationalize_input(resid[i], rbits);
        if (!is_polynomial(polys[i], vars, (size_t)nvar)) { all_poly = false; }
    }
    for (size_t i = 0; i < nres; i++) expr_free(resid[i]);
    free(resid);

    Expr* out = NULL;
    if (all_poly) {
        int verify = (o->verify != 0);      /* Automatic (-1) and True (1) both verify */
        if (o->method != NSM_SYMBOLIC) {
            /* Default: eigenvalue / multiplication-matrix engine. */
            out = nsolve_polynomial_system(polys, (int)nres, vars, nvar,
                                           NSYS_ENDOMORPHISM, reals_only,
                                           want_machine, bits, o->max_roots,
                                           verify, o->seed);
        }
        /* Method -> "Symbolic", or eigenvalue engine bailed: elimination. */
        if (!out) {
            out = nsolve_system_eliminate(polys, (int)nres, vars, nvar,
                                          reals_only, want_machine, bits,
                                          o->max_roots, verify);
        }
    }
    for (size_t i = 0; i < nres; i++) expr_free(polys[i]);
    free(polys);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Transcendental grid-seeding (univariate, last resort after Solve)
 * ------------------------------------------------------------------ */

/* Read a numeric Expr as a complex double; false for symbolic input. */
static bool to_cdouble(const Expr* e, double* re, double* im) {
    int64_t n, d;
    *im = 0.0;
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *re = (double)e->data.integer; return true;
        case EXPR_REAL:    *re = e->data.real;            return true;
        case EXPR_BIGINT:  *re = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational((Expr*)e, &n, &d) && d != 0) { *re = (double)n / (double)d; return true; }
            if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2) {
                double a, b, t;
                if (to_cdouble(e->data.function.args[0], &a, &t)
                    && to_cdouble(e->data.function.args[1], &b, &t)) { *re = a; *im = b; return true; }
            }
            return false;
        default: return false;
    }
}

/* |residual(var -> val)| at machine precision. */
static double residual_mag_at(Expr* residual, Expr* var, Expr* val) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(var), expr_copy(val) }, 2);
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
    Expr* ra = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                   (Expr*[]){ expr_copy(residual), rl }, 2);
    Expr* sub = eval_and_free(ra);
    NumericSpec spec = numeric_machine_spec();
    Expr* num = eval_and_free(numericalize(sub, spec));
    expr_free(sub);
    double re, im;
    double m = to_cdouble(num, &re, &im) ? hypot(re, im) : INFINITY;
    expr_free(num);
    return m;
}

/* Newton-solve `residual == 0` from the seed (re0, im0) via FindRoot; returns
 * the owned value Expr, or NULL if FindRoot fails. */
static Expr* findroot_at(Expr* residual, Expr* var, double re0, double im0,
                         bool want_machine, long bits) {
    Expr* x0 = (im0 == 0.0)
        ? expr_new_real(re0)
        : expr_new_function(expr_new_symbol(SYM_Complex),
              (Expr*[]){ expr_new_real(re0), expr_new_real(im0) }, 2);
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_copy(var), x0 }, 2);
    Expr* call;
    if (!want_machine) {
        Expr* wp = expr_new_function(expr_new_symbol(SYM_Rule),
                       (Expr*[]){ expr_new_symbol(SYM_WorkingPrecision),
                                  expr_new_real(numeric_bits_to_digits(bits)) }, 2);
        call = expr_new_function(expr_new_symbol(SYM_FindRoot),
                   (Expr*[]){ expr_copy(residual), spec, wp }, 3);
    } else {
        call = expr_new_function(expr_new_symbol(SYM_FindRoot),
                   (Expr*[]){ expr_copy(residual), spec }, 2);
    }
    Expr* fr = builtin_findroot(call);
    expr_free(call);
    if (!fr) return NULL;
    Expr* val = NULL;
    if (is_head(fr, SYM_List) && fr->data.function.arg_count >= 1) {
        Expr* rule = fr->data.function.args[0];
        if (is_head(rule, SYM_Rule) && rule->data.function.arg_count == 2)
            val = expr_copy(rule->data.function.args[1]);
    }
    expr_free(fr);
    return val;
}

/* Seed FindRoot from a grid and collect distinct verified roots of a single
 * non-polynomial equation.  Inherently incomplete (a finite sample). */
static Expr* nsolve_seed_transcendental(Expr* expr, Expr* var, bool reals_only,
                                        bool want_machine, long bits,
                                        long max_roots) {
    Expr* residual = equation_residual(expr);
    if (!residual) return NULL;

    double vtol = want_machine ? 1e-6 : pow(10.0, -0.4 * numeric_bits_to_digits(bits));
    if (vtol < 1e-12) vtol = 1e-12;
    long cap = (max_roots >= 0) ? max_roots : 16;

    double* gr = NULL; double* gi = NULL;   /* found-root grid for dedup */
    Expr**  vals = NULL; int nv = 0;

    /* Real seeds across a wide interval; for the complexes, a coarse grid. */
    for (int phase = 0; phase < 2 && nv < cap; phase++) {
        if (phase == 1 && reals_only) break;
        for (double s = -30.0; s <= 30.0 && nv < cap; s += 1.5) {
            int npt = (phase == 0) ? 1 : 2;
            for (int p = 0; p < npt && nv < cap; p++) {
                double im0 = (phase == 0) ? 0.0 : ((p == 0) ? 2.0 : -2.0);
                Expr* val = findroot_at(residual, var, s, im0, want_machine, bits);
                if (!val) continue;
                double vre, vim;
                bool num = to_cdouble(val, &vre, &vim);
                if (!num || residual_mag_at(residual, var, val) > vtol
                    || (reals_only && fabs(vim) > 1e-6)) { expr_free(val); continue; }
                bool dup = false;
                for (int k = 0; k < nv; k++)
                    if (hypot(vre - gr[k], vim - gi[k]) < 1e-5) { dup = true; break; }
                if (dup) { expr_free(val); continue; }
                gr = realloc(gr, sizeof(double) * (size_t)(nv + 1));
                gi = realloc(gi, sizeof(double) * (size_t)(nv + 1));
                vals = realloc(vals, sizeof(Expr*) * (size_t)(nv + 1));
                gr[nv] = vre; gi[nv] = vim; vals[nv] = val; nv++;
            }
        }
    }
    expr_free(residual);
    free(gr); free(gi);

    if (nv == 0) { free(vals); return NULL; }
    Expr** sols = (Expr**)malloc(sizeof(Expr*) * (size_t)nv);
    for (int i = 0; i < nv; i++) sols[i] = make_rule_list(var, vals[i]);
    free(vals);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), sols, (size_t)nv);
    free(sols);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Builtin entry
 * ------------------------------------------------------------------ */

Expr* builtin_nsolve(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && nsolve_is_option_arg(res->data.function.args[pos_end - 1]))
        pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nsolve_is_option_arg(res->data.function.args[i])) return NULL;
    }
    if (pos_end < 1 || pos_end > 4) return NULL;

    NSolveOpts opts = { -1, -1.0, NSM_AUTOMATIC, -1, 1234UL };
    for (size_t i = pos_end; i < argc; i++)
        nsolve_apply_option(res->data.function.args[i], &opts);

    Expr* expr = res->data.function.args[0];

    /* Tautology / contradiction short-circuits (the equation system is
     * evaluated before NSolve sees it). */
    if (expr->type == EXPR_SYMBOL && expr->data.symbol.name == SYM_True) {
        Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ empty }, 1);
    }
    if (expr->type == EXPR_SYMBOL && expr->data.symbol.name == SYM_False)
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    /* Parse the optional positional vars / domain / precision arguments. */
    size_t idx = 1;
    Expr* vars = NULL;   /* borrowed from res, or NULL = auto */
    Expr* dom  = NULL;   /* borrowed from res, or NULL = Complexes */

    if (idx < pos_end && is_var_spec(res->data.function.args[idx])) {
        vars = res->data.function.args[idx];
        idx++;
    }
    if (idx < pos_end && is_domain_symbol(res->data.function.args[idx])) {
        dom = res->data.function.args[idx];
        idx++;
    }
    {
        double v;
        if (idx < pos_end && expr_to_double(res->data.function.args[idx], &v)) {
            if (v > 0.0 && opts.prec_digits < 0.0) opts.prec_digits = v;
            idx++;
        }
    }
    if (idx != pos_end) return NULL;  /* leftover positional args */

    bool reals_only = dom && dom->data.symbol.name == SYM_Reals;

    /* Working precision -> machine flag + bit width. */
    bool want_machine = (opts.prec_digits < 0.0);
    long bits = 53;
    if (!want_machine) {
        bits = numeric_digits_to_bits(opts.prec_digits);
        if (bits <= 53) { bits = 53; want_machine = true; }
    }

    /* Resolve the variable list (owned). */
    Expr* varlist;
    if (vars) {
        varlist = (vars->type == EXPR_SYMBOL)
            ? expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_copy(vars) }, 1)
            : expr_copy(vars);
    } else {
        varlist = auto_variable_list(expr);
        if (!varlist) return NULL;  /* no variables found */
    }
    size_t nvars = varlist->data.function.arg_count;

    /* Univariate polynomial -> NRoots (the explicit requirement).  The Integers
     * domain is left to Solve, which can restrict to concrete integers. */
    bool integers_dom = dom && dom->data.symbol.name == SYM_Integers;
    Expr* out = NULL;
    if (nvars == 1 && !integers_dom) {
        Expr* var = varlist->data.function.args[0];
        /* Degree guard: the NRoots engine (Aberth) computes every root, and
         * Solve cannot factor a degree-millions polynomial either, so an
         * enormous literal exponent (e.g. x^1000000 - 2 x + 3) would hang both
         * paths regardless of MaxRoots.  Bail cheaply, leaving NSolve
         * unevaluated, before either path is entered. */
        if (has_huge_power(expr, var, NSOLVE_MAX_POLY_DEGREE)) {
            fprintf(stderr,
                    "NSolve::deg: the polynomial degree exceeds the supported "
                    "limit (%d); leaving unevaluated.\n", NSOLVE_MAX_POLY_DEGREE);
            expr_free(varlist);
            return NULL;
        }
        out = nsolve_polynomial(expr, var, reals_only, &opts);
    }

    /* Multivariate polynomial system -> numeric system engine (eigenvalue /
     * multiplication-matrix method).  Returns NULL for non-polynomial,
     * positive-dimensional, or inexact-nonlinear systems -> Solve fallback. */
    if (!out && nvars >= 2 && !integers_dom) {
        out = nsolve_system_path(expr, varlist, reals_only, want_machine, bits, &opts);
    }

    /* Fallback: Solve symbolically, then numericalise. */
    if (!out) {
        out = nsolve_via_solve(expr, varlist, dom, want_machine, bits, opts.max_roots);
    }

    /* Last resort for a univariate non-polynomial equation Solve could not
     * reduce: seed FindRoot from a grid (best-effort, finite sample). */
    if (!out && nvars == 1 && !integers_dom) {
        Expr* var = varlist->data.function.args[0];
        out = nsolve_seed_transcendental(expr, var, reals_only, want_machine,
                                         bits, opts.max_roots);
    }

    expr_free(varlist);
    return out;  /* evaluator frees res on non-NULL return */
}

/* ------------------------------------------------------------------ *
 *  Registration
 * ------------------------------------------------------------------ */
void nsolve_init(void) {
    symtab_add_builtin("NSolve", builtin_nsolve);
    /* Protected only.  The equation system is evaluated normally before NSolve
     * runs (numeric (in)equalities fold, variables stay symbolic). */
    symtab_get_def("NSolve")->attributes |= ATTR_PROTECTED;
}

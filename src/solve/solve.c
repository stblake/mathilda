/*
 * solve.c
 *
 * The `Solve` router: classifies the input equation system, parses
 * options, and dispatches to a specialist solver.  The only specialist
 * wired up in this initial cut is Solve`SolvePolynomialEquality
 * (src/solvepoly.c) for a single polynomial equality in one variable.
 *
 * `Solve` does not hold its arguments -- the evaluator delivers
 * `expr` and `vars` already evaluated, matching Mathematica's
 * attribute set ({Protected}).  When `vars` has been substituted to
 * a non-symbol (typically because the user previously assigned
 * `x = 5` and then called `Solve[..., x]`), the router emits
 * `Solve::ivar` and returns unevaluated.
 */

#include "solve.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "print.h"
#include "solveinv.h"
#include "solvelinsys.h"
#include "solvemod.h"
#include "solveint.h"
#include "solvenlsys.h"
#include "solvepoly.h"
#include "reduce_zerodim.h"
#include "solverad.h"
#include "solvetrig.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"
#include "zero_test.h"

/* ------------------------------------------------------------------ *
 *  Option parsing.                                                    *
 * ------------------------------------------------------------------ */

typedef struct {
    SolvePolyOpts poly;
    SolveInvOpts  inv;
    Expr* dom;             /* borrowed; default = NULL ( = Complexes) */
    Expr* modulus;         /* borrowed; rhs of Modulus -> p, else NULL  */
    bool  verify_on;       /* VerifySolutions -> True                   */
} SolveOpts;

/* Recognised Solve option-name symbols. */
static bool is_known_option_name(const char* s) {
    return s == SYM_Cubics
        || s == SYM_Quartics
        || s == SYM_GeneratedParameters
        || s == SYM_VerifySolutions
        || s == SYM_Assumptions
        || s == SYM_InverseFunctions
        || s == SYM_Method
        || s == SYM_Modulus;
}

/* True iff `e` is Rule[opt, _] or RuleDelayed[opt, _] for a recognised
 * Solve option name.  Used to peel options off the end of the arg list. */
static bool is_option_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return is_known_option_name(lhs->data.symbol.name);
}

/* Returns true iff `e` is the symbol True. */
static bool is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}

/* True iff `e` is the symbol False. */
static bool is_false(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False;
}

/* Apply a single option rule to `opts`.  Unknown values do not abort
 * (they are silently ignored for now) -- only unknown option *names*
 * are rejected, by is_known_option_name. */
static void apply_option(const Expr* rule, SolveOpts* opts) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    if (name == SYM_Cubics)   { opts->poly.cubics_radical = is_true(rhs); return; }
    if (name == SYM_Quartics) { opts->poly.quartics_radical = is_true(rhs); return; }
    if (name == SYM_InverseFunctions) {
        /* InverseFunctions -> False disables the specialist; any other
         * value (True / Automatic / unrecognised) leaves the default
         * `enabled = true` in place. */
        if (is_false(rhs)) opts->inv.enabled = false;
        return;
    }
    if (name == SYM_GeneratedParameters) {
        /* Bare symbol form only -- the Function form is reserved. */
        if (rhs && rhs->type == EXPR_SYMBOL) {
            opts->inv.param_head = rhs->data.symbol.name;
        }
        return;
    }
    if (name == SYM_Modulus) {
        /* Store the modulus rhs (borrowed).  builtin_solve dispatches to
         * the modular pre-pass when this is set; the value is validated
         * there. */
        opts->modulus = (Expr*)rhs;
        return;
    }
    if (name == SYM_VerifySolutions) {
        /* VerifySolutions -> True enables the post-dispatch
         * PossibleZeroQ back-substitution filter (below).  False /
         * Automatic keep the default per-specialist behaviour. */
        opts->verify_on = is_true(rhs);
        return;
    }
    /* Assumptions / Method: parsed but not yet wired into the
     * polynomial specialist. */
}

/* Warn once per distinct unevaluated form that the second argument
 * is not a valid variable specification (a symbol or a list of
 * symbols).  Mirrors Mathematica's `Solve::ivar`. */
static void warn_ivar(const Expr* vars) {
    static uint64_t last_warned_hash = 0;
    if (!vars) return;
    uint64_t h = expr_hash(vars);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    char* shown = expr_to_string((Expr*)vars);
    fprintf(stderr,
        "Solve::ivar: %s is not a valid variable.\n",
        shown ? shown : "?");
    free(shown);
}

/* Warn once per distinct unevaluated form about an unrecognised
 * option.  Mirrors the integrate.c:254-262 idiom. */
static void warn_bad_option(const Expr* res, const Expr* opt) {
    static uint64_t last_warned_hash = 0;
    uint64_t h = expr_hash(res);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    const Expr* lhs = (opt && opt->type == EXPR_FUNCTION
                       && opt->data.function.arg_count == 2)
        ? opt->data.function.args[0] : NULL;
    const char* name = (lhs && lhs->type == EXPR_SYMBOL)
        ? lhs->data.symbol.name : "?";
    fprintf(stderr,
        "Solve::optx: Unknown option %s in Solve.\n",
        name);
}

/* ------------------------------------------------------------------ *
 *  Argument classification.                                           *
 * ------------------------------------------------------------------ */

/* Returns the single variable contained in `vars` (a symbol or a
 * length-1 List of one symbol) and writes it to *var_out.  Returns
 * false if `vars` is not a supported shape (e.g. multivariate). */
static bool classify_single_var(Expr* vars, Expr** var_out) {
    if (vars->type == EXPR_SYMBOL) { *var_out = vars; return true; }
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List
        && vars->data.function.arg_count == 1
        && vars->data.function.args[0]->type == EXPR_SYMBOL) {
        *var_out = vars->data.function.args[0];
        return true;
    }
    return false;
}

/* True iff `expr` is the kind of compound that the linear-system
 * specialist accepts: an `And` of equations, a `List` of equations,
 * or a single `Equal` that the caller has marked as multivariate.
 * The detailed shape check (each conjunct is `Equal[_, _]`) is done
 * inside solvelinsys_solve_linear_system itself. */
static bool is_conjunction_of_equations(const Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return false;
    if (expr->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = expr->data.function.head->data.symbol.name;
    return h == SYM_And || h == SYM_List;
}

/* True for anything Mathematica would reject as a "solve variable":
 * a bare numeric atom (Integer / BigInt / Real / MPFR), a string, or
 * a packaged numeric head (`Rational[_, _]`, `Complex[_, _]`).  Bare
 * compound expressions like `Dt[y]`, `f[a, b]`, or `x^2` are allowed
 * as solve variables -- the dispatch later substitutes them through
 * a fresh internal symbol. */
static bool is_numeric_literal(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
        case EXPR_STRING:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION:
            if (e->data.function.head->type == EXPR_SYMBOL) {
                const char* h = e->data.function.head->data.symbol.name;
                if (h == SYM_Complex || h == SYM_Rational) return true;
            }
            return false;
        default:
            return false;
    }
}

/* Mathematica-style validation of the `vars` argument.  Accepts a
 * single non-numeric expression (symbol or compound), or a non-empty
 * List of such expressions.  Everything else triggers `Solve::ivar`
 * and leaves Solve unevaluated. */
static bool is_valid_solve_vars(const Expr* vars) {
    if (!vars) return false;
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        if (vars->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < vars->data.function.arg_count; i++) {
            const Expr* v = vars->data.function.args[i];
            if (!v || is_numeric_literal(v)) return false;
        }
        return true;
    }
    return !is_numeric_literal(vars);
}

/* ------------------------------------------------------------------ *
 *  Compound-variable substitution.                                    *
 *                                                                     *
 *  Mathematica lets `Solve[..., g]` take a non-symbol generalised      *
 *  variable -- typical examples are `Dt[y]`, `f[a, b]`, or `x^2`.     *
 *  Mathilda's specialists (polynomial, linear-system, inverse-       *
 *  function, radicals, trig) all expect a bare symbol, so the router  *
 *  rewrites compound vars by substituting each one with a fresh       *
 *  internal symbol (`Solve$var$N`) before dispatch, then reverses     *
 *  the substitution on the result so the user sees `{{g -> ...}}`.   *
 *                                                                     *
 *  Cap of 32 substitutions per call is more than enough -- a real    *
 *  Solve call rarely has more than a handful of distinct variables.   *
 * ------------------------------------------------------------------ */

#define SOLVE_MAX_VAR_SUBS 32

typedef struct {
    Expr* original;     /* borrowed from caller (lives as long as `res`) */
    const char* fresh;  /* interned symbol name                          */
} SolveVarSub;

/* Per-process counter for generating fresh internal symbol names.
 * Monotonic so distinct Solve calls don't collide -- the symbols are
 * never visible to the user (they only exist between the dispatch
 * and unsubst). */
static const char* gen_fresh_var_name(void) {
    static int counter = 0;
    char buf[64];
    snprintf(buf, sizeof(buf), "Solve$var$%d", ++counter);
    return intern_symbol(buf);
}

/* Replace every subexpression structurally equal to `from` with a
 * copy of `to`.  Used to substitute a compound variable with its
 * fresh internal symbol throughout the equation. */
static Expr* subst_expr(Expr* e, Expr* from, Expr* to) {
    if (!e) return NULL;
    if (expr_eq(e, from)) return expr_copy(to);
    if (e->type == EXPR_FUNCTION) {
        Expr* new_head = subst_expr(e->data.function.head, from, to);
        size_t n = e->data.function.arg_count;
        Expr** new_args = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = subst_expr(e->data.function.args[i], from, to);
        }
        Expr* result = expr_new_function(new_head, new_args, n);
        if (new_args) free(new_args);
        return result;
    }
    return expr_copy(e);
}

/* Reverse pass: replace every fresh-symbol leaf in `e` with the
 * original compound it stood for.  Each Rule LHS produced by the
 * dispatch carries the fresh symbol; this puts the user's `Dt[y]`
 * (etc.) back so the result reads `{{Dt[y] -> ...}}`. */
static Expr* unsubst_compound_vars(
    Expr* e,
    const SolveVarSub* subs,
    size_t n_subs)
{
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n_subs; i++) {
            if (e->data.symbol.name == subs[i].fresh) {
                return expr_copy(subs[i].original);
            }
        }
        return expr_copy(e);
    }
    if (e->type == EXPR_FUNCTION) {
        Expr* new_head = unsubst_compound_vars(e->data.function.head, subs, n_subs);
        size_t n = e->data.function.arg_count;
        Expr** new_args = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = unsubst_compound_vars(e->data.function.args[i], subs, n_subs);
        }
        Expr* result = expr_new_function(new_head, new_args, n);
        if (new_args) free(new_args);
        return result;
    }
    return expr_copy(e);
}

/* Pre-pass that wraps the dispatch.  Walks `vars_in`, allocates a
 * fresh symbol for every non-symbol entry, substitutes that entry
 * through `*expr_inout`, and builds `*vars_out_owned` -- a freshly
 * owned variable specification where every element is a symbol.
 *
 * Ownership:
 *   - `*expr_inout` enters owned and stays owned (possibly replaced).
 *   - `*vars_out_owned` is freshly allocated; caller must free.
 *   - `subs[i].original` borrows from the caller's `vars_in`; valid
 *     as long as `res` is.
 *
 * Returns the number of fresh-symbol substitutions installed (0
 * means the user's vars were already symbol-only and no substitution
 * happened -- `*vars_out_owned` is then just a copy of `vars_in`). */
static size_t collect_and_subst_compound_vars(
    Expr* vars_in,
    Expr** expr_inout,
    SolveVarSub* subs,
    Expr** vars_out_owned)
{
    size_t n_subs = 0;

    bool is_list = (vars_in->type == EXPR_FUNCTION
        && vars_in->data.function.head->type == EXPR_SYMBOL
        && vars_in->data.function.head->data.symbol.name == SYM_List);

    if (!is_list) {
        if (vars_in->type == EXPR_SYMBOL) {
            *vars_out_owned = expr_copy(vars_in);
            return 0;
        }
        /* Bare compound variable. */
        const char* fresh = gen_fresh_var_name();
        Expr* fresh_sym = expr_new_symbol(fresh);
        Expr* new_expr = subst_expr(*expr_inout, vars_in, fresh_sym);
        expr_free(*expr_inout);
        *expr_inout = new_expr;
        subs[0].original = vars_in;
        subs[0].fresh = fresh;
        *vars_out_owned = fresh_sym;
        return 1;
    }

    size_t n = vars_in->data.function.arg_count;
    Expr** new_args = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* v = vars_in->data.function.args[i];
        if (v->type == EXPR_SYMBOL || n_subs >= SOLVE_MAX_VAR_SUBS) {
            new_args[i] = expr_copy(v);
            continue;
        }
        const char* fresh = gen_fresh_var_name();
        Expr* fresh_sym = expr_new_symbol(fresh);
        Expr* new_expr = subst_expr(*expr_inout, v, fresh_sym);
        expr_free(*expr_inout);
        *expr_inout = new_expr;
        subs[n_subs].original = v;
        subs[n_subs].fresh = fresh;
        n_subs++;
        new_args[i] = fresh_sym;
    }
    *vars_out_owned = expr_new_function(expr_new_symbol(SYM_List), new_args, n);
    if (new_args) free(new_args);
    return n_subs;
}

/* True iff `vars` is a List of at least two symbols.  Used to route
 * single-Equal-but-multi-variable inputs through the linear-system
 * specialist. */
static bool is_multi_var_list(const Expr* vars) {
    return vars
        && vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List
        && vars->data.function.arg_count >= 2;
}

/* Wrap a single Expr* (borrowed) in a freshly allocated `List[expr]`.
 * The caller takes ownership of the returned wrapper. */
static Expr* wrap_in_list(Expr* expr) {
    return expr_new_function(expr_new_symbol(SYM_List),
                             (Expr*[]){ expr_copy(expr) }, 1);
}

/* Recognise Equal[Abs[u], 0] / Equal[0, Abs[u]] and rewrite as
 * Equal[u, 0] so the polynomial dispatch can solve u directly.
 * Mirrors Maxima's easy-cases `mabs` shortcut.  Returns a freshly
 * owned Expr* on rewrite, NULL otherwise. */
static Expr* try_abs_zero_rewrite(const Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return NULL;
    if (expr->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (expr->data.function.head->data.symbol.name != SYM_Equal) return NULL;
    if (expr->data.function.arg_count != 2) return NULL;
    const Expr* lhs = expr->data.function.args[0];
    const Expr* rhs = expr->data.function.args[1];
    const Expr* abs_side = NULL;
    const Expr* zero_side = NULL;
    if (lhs->type == EXPR_FUNCTION
        && lhs->data.function.head->type == EXPR_SYMBOL
        && lhs->data.function.head->data.symbol.name == SYM_Abs
        && lhs->data.function.arg_count == 1) {
        abs_side = lhs; zero_side = rhs;
    } else if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol.name == SYM_Abs
        && rhs->data.function.arg_count == 1) {
        abs_side = rhs; zero_side = lhs;
    }
    if (!abs_side) return NULL;
    if (zero_side->type != EXPR_INTEGER || zero_side->data.integer != 0)
        return NULL;
    Expr* u = abs_side->data.function.args[0];
    return expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy(u), expr_new_integer(0) }, 2);
}

/* ------------------------------------------------------------------ *
 *  Builtin entry.                                                     *
 * ------------------------------------------------------------------ */

/* True unless some Equal in `eqn` back-substitutes, under the solution
 * rule-list `sol`, to a residual that PossibleZeroQ proves non-zero
 * (zero_test_decide == ZERO_TEST_FALSE).  And / List systems are
 * descended into.  Undecidable residuals (Root[], free parameters,
 * ConditionalExpression) are ZERO_TEST_UNKNOWN and therefore kept --
 * matching Mathematica's "verify only when you can" policy. */
static bool solution_verifies(const Expr* eqn, Expr* sol) {
    if (!eqn || eqn->type != EXPR_FUNCTION
        || eqn->data.function.head->type != EXPR_SYMBOL) {
        return true;
    }
    const char* h = eqn->data.function.head->data.symbol.name;
    if (h == SYM_And || h == SYM_List) {
        for (size_t i = 0; i < eqn->data.function.arg_count; i++)
            if (!solution_verifies(eqn->data.function.args[i], sol)) return false;
        return true;
    }
    if (h == SYM_Equal && eqn->data.function.arg_count == 2) {
        Expr* lhs = eqn->data.function.args[0];
        Expr* rhs = eqn->data.function.args[1];
        Expr* residual = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){
                expr_copy(lhs),
                expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2)
            }, 2));
        Expr* subbed = eval_and_free(internal_replace_all(
            (Expr*[]){ residual, expr_copy(sol) }, 2));
        ZeroTestResult zt = zero_test_decide(subbed);
        expr_free(subbed);
        return zt != ZERO_TEST_FALSE;
    }
    return true;  /* non-equation operand (e.g. True): trivially satisfied */
}

/* VerifySolutions -> True: drop every solution rule-list in `out` whose
 * back-substitution into `equation` is decidably non-zero.  `out` is
 * List[List[Rule...], ...]; takes ownership and returns the filtered
 * List.  `equation` is borrowed (the substituted / rationalised form
 * seen by the dispatch). */
static Expr* verify_solutions_filter(Expr* out, const Expr* equation) {
    if (!out || out->type != EXPR_FUNCTION
        || out->data.function.head->type != EXPR_SYMBOL
        || out->data.function.head->data.symbol.name != SYM_List) {
        return out;
    }
    size_t n = out->data.function.arg_count;
    Expr** kept = n ? (Expr**)malloc(sizeof(Expr*) * n) : NULL;
    size_t nkept = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = out->data.function.args[i];
        if (solution_verifies(equation, sol)) {
            kept[nkept++] = sol;
        } else {
            expr_free(sol);
        }
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), kept, nkept);
    free(kept);
    /* Detach the moved args so freeing the old wrapper does not touch
     * them (kept args now belong to `result`; dropped ones are freed). */
    out->data.function.arg_count = 0;
    expr_free(out);
    return result;
}

/* Read a *concrete* numeric leaf (Integer / Real / BigInt / Rational / MPFR)
 * as a double.  Returns false for symbolic / non-numeric forms.  Used only
 * to inspect the imaginary part of a numericalized Complex[re, im]. */
static bool solve_concrete_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        double n, d;
        if (solve_concrete_double(e->data.function.args[0], &n)
            && solve_concrete_double(e->data.function.args[1], &d) && d != 0.0) {
            *out = n / d;
            return true;
        }
    }
    return false;
}

/* True iff `val` numericalizes to a Complex number with a concrete,
 * non-negligible imaginary part -- i.e. `val` is *provably* non-real.  A
 * Root[] object (companion-matrix + Sturm + Newton) numericalizes to a real
 * MPFR/Real for a real root and to Complex[re, im] for a complex one, so this
 * distinguishes them; likewise for closed-form radical values and any
 * arithmetic expression over them (e.g. a system's x-component).  Values that
 * do NOT reduce to a concrete Complex -- ordinary reals, or symbolic /
 * parametric residues like Sqrt[a] and un-numericalizable Root objects -- are
 * not provably non-real and return false, so they are kept.  `val` is
 * borrowed. */
static bool value_is_provably_nonreal(Expr* val) {
    Expr* nv = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy(val) }, 1));
    bool nonreal = false;
    if (nv && nv->type == EXPR_FUNCTION
        && nv->data.function.head->type == EXPR_SYMBOL
        && nv->data.function.head->data.symbol.name == SYM_Complex
        && nv->data.function.arg_count == 2) {
        double im;
        if (solve_concrete_double(nv->data.function.args[1], &im)
            && fabs(im) > 1e-9) {
            nonreal = true;
        }
    }
    expr_free(nv);
    return nonreal;
}

/* Reals-domain reality filter.  Drops every solution rule-list in `out` that
 * binds any variable to a *provably* non-real value.  `out` is
 * List[List[Rule...], ...] -- one Rule per variable, several for a system --
 * so the whole tuple is dropped when any component is complex.  Takes
 * ownership and returns the filtered List; mirrors verify_solutions_filter's
 * alloc / detach / free shape.  Conservative by construction (see
 * value_is_provably_nonreal): symbolic and parametric real answers survive. */
static Expr* filter_reals_solutions(Expr* out) {
    if (!out || out->type != EXPR_FUNCTION
        || out->data.function.head->type != EXPR_SYMBOL
        || out->data.function.head->data.symbol.name != SYM_List) {
        return out;
    }
    size_t n = out->data.function.arg_count;
    Expr** kept = n ? (Expr**)malloc(sizeof(Expr*) * n) : NULL;
    size_t nkept = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = out->data.function.args[i];
        bool nonreal = false;
        if (sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol.name == SYM_List) {
            for (size_t j = 0; j < sol->data.function.arg_count && !nonreal; j++) {
                Expr* rule = sol->data.function.args[j];
                if (rule->type == EXPR_FUNCTION
                    && rule->data.function.head->type == EXPR_SYMBOL
                    && (rule->data.function.head->data.symbol.name == SYM_Rule
                        || rule->data.function.head->data.symbol.name
                               == SYM_RuleDelayed)
                    && rule->data.function.arg_count == 2) {
                    if (value_is_provably_nonreal(rule->data.function.args[1]))
                        nonreal = true;
                }
            }
        }
        if (nonreal) {
            expr_free(sol);
        } else {
            kept[nkept++] = sol;
        }
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), kept, nkept);
    free(kept);
    out->data.function.arg_count = 0;
    expr_free(out);
    return result;
}

/* True iff `val` is provably NOT an integer: a concrete numeric value that
 * does not sit on an integer.  Symbolic / parametric values (free variables,
 * ConditionalExpression, un-numericalizable Root objects) are not provably
 * non-integer and return false, so they survive -- matching the conservative
 * "only drop what you can decide" policy of the reals filter.  `val` is
 * borrowed. */
static bool value_is_provably_non_integer(Expr* val) {
    if (!val) return false;
    if (val->type == EXPR_INTEGER || val->type == EXPR_BIGINT) return false;
    Expr* nv = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy(val) }, 1));
    bool nonint = false;
    if (nv && nv->type == EXPR_FUNCTION
        && nv->data.function.head->type == EXPR_SYMBOL
        && nv->data.function.head->data.symbol.name == SYM_Complex
        && nv->data.function.arg_count == 2) {
        double im;
        if (solve_concrete_double(nv->data.function.args[1], &im) && fabs(im) > 1e-9)
            nonint = true;                     /* non-real -> non-integer */
    } else {
        double d;
        if (solve_concrete_double(nv, &d)) {
            double r = (d < 0) ? -floor(-d + 0.5) : floor(d + 0.5);
            if (fabs(d - r) > 1e-9) nonint = true;
        }
    }
    expr_free(nv);
    return nonint;
}

/* Integers-domain restriction: drop every solution tuple that binds any
 * variable to a provably-non-integer value.  Mirrors the alloc / detach /
 * free shape of filter_reals_solutions.  Takes ownership of `out`. */
static Expr* filter_integers_solutions(Expr* out) {
    if (!out || out->type != EXPR_FUNCTION
        || out->data.function.head->type != EXPR_SYMBOL
        || out->data.function.head->data.symbol.name != SYM_List) {
        return out;
    }
    size_t n = out->data.function.arg_count;
    Expr** kept = n ? (Expr**)malloc(sizeof(Expr*) * n) : NULL;
    size_t nkept = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = out->data.function.args[i];
        bool bad = false;
        if (sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol.name == SYM_List) {
            for (size_t j = 0; j < sol->data.function.arg_count && !bad; j++) {
                Expr* rule = sol->data.function.args[j];
                if (rule->type == EXPR_FUNCTION
                    && rule->data.function.head->type == EXPR_SYMBOL
                    && (rule->data.function.head->data.symbol.name == SYM_Rule
                        || rule->data.function.head->data.symbol.name == SYM_RuleDelayed)
                    && rule->data.function.arg_count == 2) {
                    if (value_is_provably_non_integer(rule->data.function.args[1]))
                        bad = true;
                }
            }
        }
        if (bad) expr_free(sol); else kept[nkept++] = sol;
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), kept, nkept);
    free(kept);
    out->data.function.arg_count = 0;
    expr_free(out);
    return result;
}

/* Does the interned symbol name `sym` occur anywhere in `e`? */
static bool expr_mentions_symbol(const Expr* e, const char* sym) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == sym;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_mentions_symbol(e->data.function.head, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_mentions_symbol(e->data.function.args[i], sym)) return true;
    return false;
}

/* True iff `eq` is a lone `Equal` (not a conjunction / list / system) that
 * mentions at least two of the solve variables -- i.e. a positive-dimensional
 * curve/surface.  Its integer points are the exclusive province of the
 * `solveint` pre-pass; once that has declined, the Complexes-oriented
 * parametric dispatch (try_single_eq_multivar / the poly specialist over
 * Integers) can only fabricate a spurious `{}` -- e.g. y^2 == x^3 - 2 or
 * y == x^2 -- because the closed-form root is not integer-valued for the
 * symbolic parameter.  This predicate lets the caller convert that `{}` into
 * an unevaluated Solve (safe: an unbounded curve has no finite integer
 * enumeration we could offer here anyway). */
static bool is_single_multivar_equation(const Expr* eq, const Expr* vars) {
    if (!eq || eq->type != EXPR_FUNCTION
        || eq->data.function.head->type != EXPR_SYMBOL
        || eq->data.function.head->data.symbol.name != SYM_Equal
        || eq->data.function.arg_count != 2) return false;
    if (!vars || vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) return false;
    int mentioned = 0;
    for (size_t i = 0; i < vars->data.function.arg_count; i++) {
        Expr* v = vars->data.function.args[i];
        if (v->type == EXPR_SYMBOL
            && expr_mentions_symbol(eq, v->data.symbol.name)) mentioned++;
    }
    return mentioned >= 2;
}

/* True iff `e` mentions any of the solve variables in `vars` (a bare symbol
 * or a List of symbols).  Borrowed. */
static bool mentions_any_solve_var(const Expr* e, const Expr* vars) {
    if (!e || !vars) return false;
    if (vars->type == EXPR_SYMBOL)
        return expr_mentions_symbol(e, vars->data.symbol.name);
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < vars->data.function.arg_count; i++) {
            Expr* v = vars->data.function.args[i];
            if (v && v->type == EXPR_SYMBOL
                && expr_mentions_symbol(e, v->data.symbol.name)) return true;
        }
    }
    return false;
}

/* True iff the Solve result `out` (a List of solution branches, each a List
 * of Rule[var, rhs]) binds some variable to an expression that still mentions
 * a solve variable -- the signature of a positive-dimensional / parametric
 * answer (the `Solve::svars` case, e.g. {{y -> (4 + 2 x)/5, z -> ...}} for an
 * underdetermined linear system, or {{y -> Sqrt[x^3 - 2]}} for a lone curve).
 *
 * Over the Integers this distinguishes a genuine empty set -- every concrete
 * branch failed the integer test, a real proof of no solution (e.g. x + y ==
 * 2 && x - y == 1 -> {3/2, 1/2} -> {}) -- from a parametric family the integer
 * filter merely could not represent, which must NOT collapse to a spurious
 * `{}` (an underdetermined system has integer points; leaving Solve
 * unevaluated is the honest answer until the HNF path expresses them). */
static bool solution_set_is_parametric(const Expr* out, const Expr* vars) {
    if (!out || out->type != EXPR_FUNCTION
        || out->data.function.head->type != EXPR_SYMBOL
        || out->data.function.head->data.symbol.name != SYM_List) return false;
    for (size_t b = 0; b < out->data.function.arg_count; b++) {
        Expr* branch = out->data.function.args[b];
        if (!branch || branch->type != EXPR_FUNCTION
            || branch->data.function.head->type != EXPR_SYMBOL
            || branch->data.function.head->data.symbol.name != SYM_List) continue;
        for (size_t r = 0; r < branch->data.function.arg_count; r++) {
            Expr* rule = branch->data.function.args[r];
            if (rule && rule->type == EXPR_FUNCTION
                && rule->data.function.head->type == EXPR_SYMBOL
                && (rule->data.function.head->data.symbol.name == SYM_Rule
                    || rule->data.function.head->data.symbol.name == SYM_RuleDelayed)
                && rule->data.function.arg_count == 2
                && mentions_any_solve_var(rule->data.function.args[1], vars))
                return true;
        }
    }
    return false;
}

/* Single equation in >= 2 variables: solve for the earliest-listed
 * variable the equation is polynomial in, treating the other variables
 * as symbolic parameters.  Mathematica's Solve returns explicit rules
 * for these -- {{x -> 1/y}} for x y == 1, {{x -> -Sqrt[1-y^2]}, {x ->
 * Sqrt[1-y^2]}} for x^2 + y^2 == 1 -- NOT a Reduce-style case split, so
 * this stays firmly on the explicit-rules side of the Solve/Reduce line.
 * Returns the solution List, or NULL to defer to the nonlinear-system
 * specialist (which emits Solve::nsdim for genuinely positive-
 * dimensional systems it cannot express as rules for one variable).
 * All arguments are borrowed. */
static Expr* try_single_eq_multivar(Expr* equation, Expr* vars_list,
                                    Expr* dom, const SolvePolyOpts* polyopts) {
    if (!equation || equation->type != EXPR_FUNCTION
        || equation->data.function.head->type != EXPR_SYMBOL
        || equation->data.function.head->data.symbol.name != SYM_Equal
        || equation->data.function.arg_count != 2) {
        return NULL;
    }
    if (!vars_list || vars_list->type != EXPR_FUNCTION
        || vars_list->data.function.head->type != EXPR_SYMBOL
        || vars_list->data.function.head->data.symbol.name != SYM_List
        || vars_list->data.function.arg_count < 2) {
        return NULL;
    }
    for (size_t i = 0; i < vars_list->data.function.arg_count; i++) {
        Expr* v = vars_list->data.function.args[i];
        if (v->type != EXPR_SYMBOL) continue;
        if (!expr_mentions_symbol(equation, v->data.symbol.name)) continue;
        Expr* sol = solvepoly_solve_polynomial_equality(
            equation, v, dom, polyopts);
        if (sol) return sol;
    }
    return NULL;
}

Expr* builtin_solve(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Walk trailing args, peeling options.  Position of first option
     * = end of positional args. */
    size_t pos_end = argc;
    while (pos_end > 0) {
        Expr* a = res->data.function.args[pos_end - 1];
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol.name;
            if (is_known_option_name(name)) {
                pos_end--;
                continue;
            }
            /* Trailing Rule[] that is not a recognised option name is
             * a syntax error in Solve -- bail. */
            if (!is_option_arg(a)) {
                /* shape is right (Rule[sym, _]) but name is wrong */
                warn_bad_option(res, a);
                return NULL;
            }
        }
        break;
    }

    /* Positional args: expr [, vars [, dom]] */
    if (pos_end < 2 || pos_end > 3) {
        if (pos_end < 2) return NULL;
        /* pos_end > 3: too many positional args before options.  */
        return NULL;
    }

    /* Parse options. */
    SolveOpts opts = {
        { false, false },                /* poly: cubics, quartics */
        { true, intern_symbol("C") },    /* inv: enabled, param_head */
        NULL,                            /* dom */
        NULL,                            /* modulus */
        false                            /* verify_on */
    };
    for (size_t i = pos_end; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_option_arg(a)) apply_option(a, &opts);
    }

    /* Solve does not hold its args: both `expr` and `vars` arrive
     * already evaluated.  Take an owned copy of `expr` so the
     * downstream rationalisation/abs-zero rewrites can free-and-
     * replace it; `vars` and `dom` stay borrowed from `res`. */
    Expr* vars = res->data.function.args[1];
    Expr* dom  = (pos_end >= 3) ? res->data.function.args[2] : NULL;

    /* Mathematica-compatible bad-variable handling: emit `Solve::ivar`
     * and return unevaluated if `vars` is not a symbol or list of
     * symbols.  This is the path taken when a previously-assigned
     * OwnValue substitutes the variable to a non-symbol (most often a
     * number). */
    if (!is_valid_solve_vars(vars)) {
        warn_ivar(vars);
        return NULL;
    }

    Expr* expr = expr_copy(res->data.function.args[0]);

    /* Compound-variable pre-pass: every non-symbol entry in `vars` is
     * substituted with a fresh internal symbol throughout `expr` so
     * the dispatch specialists -- which only understand symbol
     * variables -- can run.  The substitution is reversed on the
     * result so the user sees `{{Dt[y] -> ...}}` (etc.) verbatim. */
    SolveVarSub subs[SOLVE_MAX_VAR_SUBS];
    Expr* vars_subst = NULL;
    size_t n_subs = collect_and_subst_compound_vars(
        vars, &expr, subs, &vars_subst);
    vars = vars_subst;  /* dispatch sees the symbol-only spec */

    /* Approximate-number preprocessing: if the equation system contains
     * any inexact numeric leaf, force-rationalise it so the downstream
     * specialists (polynomial / linear-system) -- which assume exact
     * arithmetic -- can run.  The result is numericalised back at the
     * tail so the user observes inexact-in / inexact-out semantics
     * consistent with Integrate and the exact-symbolic builtins
     * (Apart, Cancel, Together, Factor, ...).
     *
     * The scan also captures the *minimum* precision (in bits) across
     * every inexact leaf: that precision is then used both as the
     * rationalisation tolerance and as the precision of the final
     * numericalised result, so a 100-bit MPFR input flows back out at
     * 100 bits, a mixed Real + MPFR input drops to the lower 53 bits,
     * etc.  The vars argument (always a symbol or a list of symbols)
     * is never touched. */
    CommonInexactInfo inexact = common_scan_inexact(expr);
    if (inexact.has_inexact) {
        Expr* rationalised = common_rationalize_input(expr, inexact.min_bits);
        expr_free(expr);
        expr = rationalised;
    }

    /* True/False short-circuits regardless of the var shape.
     *   True  -> {{}}   (tautology: full-dimensional solution set)
     *   False -> {}     (contradiction: no solutions)             */
    Expr* out = NULL;
    Expr* verify_eq = NULL;   /* owned copy of the equation for VerifySolutions */
    bool integer_prepass_used = false;  /* solveint returned a (possibly empty) answer */
    bool lone_multivar_int_eq = false;  /* Integers + lone Equal in >= 2 vars */
    if (expr->type == EXPR_SYMBOL && expr->data.symbol.name == SYM_True) {
        Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        out = expr_new_function(expr_new_symbol(SYM_List),
                                (Expr*[]){ empty }, 1);
        expr_free(expr);
        expr_free(vars_subst);
        return out;
    }
    if (expr->type == EXPR_SYMBOL && expr->data.symbol.name == SYM_False) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        expr_free(expr);
        expr_free(vars_subst);
        return out;
    }

    /* Easy-case: Abs[u] == 0  -->  u == 0.  Lets the polynomial
     * dispatch solve u directly without seeing the non-polynomial Abs
     * head.  Mirrors Maxima's easy-cases `mabs` branch.  Only the bare
     * `Abs[u] == 0` shape is rewritten -- products and powers of Abs
     * are intentionally left to the standard dispatch. */
    {
        Expr* rewritten = try_abs_zero_rewrite(expr);
        if (rewritten) {
            expr_free(expr);
            expr = rewritten;
        }
    }

    /* Modulus pre-pass: Solve[poly == 0, x, Modulus -> p] solves over
     * the finite ring Z/pZ, a different value domain from the ordinary
     * characteristic-0 dispatch (whose radical / Cardano / Root[]
     * machinery is meaningless mod p).  It must therefore run before
     * dispatch and, on refusal (systems, multivariable, non-polynomial,
     * or an out-of-range modulus), return NULL so Solve stays
     * unevaluated -- never silently discarding the Modulus option.  The
     * result flows through the shared unsubst / numericalise tail (which
     * is a no-op for the integer residues modular solving returns). */
    if (opts.modulus) {
        out = solvemod_solve_modular(expr, vars, dom, opts.modulus);
        expr_free(expr);
        goto solve_finish;
    }

    /* Integer-domain (Diophantine) pre-pass: Solve[eqns && constraints,
     * vars, Integers].  Runs before the ordinary dispatch because the
     * default machinery cannot separate inequality constraints from
     * equations (an And of Greater[] conjuncts is otherwise mis-routed to
     * the linear/nonlinear system solvers, which refuse it).  Returns a
     * finite solution list, the empty List for a provably empty set, or
     * NULL to decline -- in which case we fall through to the ordinary
     * dispatch (so a bare Solve[x^2 == 4, x, Integers] still works via the
     * polynomial specialist + the integer-restriction filter below). */
    if (dom && dom->type == EXPR_SYMBOL
        && dom->data.symbol.name == SYM_Integers) {
        Expr* iout = solveint_solve_integer(expr, vars, dom);
        if (iout) {
            integer_prepass_used = true;
            expr_free(expr);
            out = iout;
            goto solve_finish;
        }
        /* solveint declined.  Remember whether this was a lone multivariable
         * equation so the shared tail can suppress a spurious `{}` fabricated
         * by the Complexes-oriented parametric dispatch (see the guard at
         * solve_finish). */
        lone_multivar_int_eq = is_single_multivar_equation(expr, vars);
    }

    /* Equations-with-constraints pre-pass (Complexes / Reals): Solve[eqns &&
     * ineqs, vars] where the equational subsystem is zero-dimensional.  Like
     * the Integers pre-pass above, the ordinary system specialists refuse an
     * And that mixes Equal with Less/Greater/Unequal, so separate the two:
     * solve the zero-dimensional equations and keep only the branches that
     * satisfy the side relations (and, over the Reals, are real), decided
     * exactly with the algebraic-number oracle.  Declines (NULL) for a
     * positive-dimensional or non-polynomial system, or when there is no side
     * constraint at all, in which case we fall through to the ordinary
     * dispatch. */
    {
        Expr* cout = reduce_zerodim_solve(expr, vars, dom, &opts.poly);
        if (cout) { expr_free(expr); out = cout; goto solve_finish; }
    }

    /* VerifySolutions -> True: snapshot the (substituted / rationalised)
     * equation now, before dispatch frees `expr`, so the post-dispatch
     * PossibleZeroQ filter can back-substitute against it. */
    if (opts.verify_on) verify_eq = expr_copy(expr);

    /* Dispatch.
     *
     *   Multi-var single Equal  ->  linear-system specialist.
     *   And/List of Equals      ->  linear-system specialist.
     *   Single var single Equal ->  polynomial-equality specialist.
     *
     * The linear-system specialist canonicalises each equation to
     * `lhs - rhs` itself and returns NULL when the input is non-affine
     * in the vars, in which case we leave Solve unevaluated. */
    bool conj = is_conjunction_of_equations(expr);
    bool multi_var = is_multi_var_list(vars);

    if (conj || multi_var) {
        /* The linear-system specialist wants `vars` as a List of
         * symbols.  When the caller passed a bare symbol, wrap it. */
        Expr* vars_list = NULL;
        if (vars->type == EXPR_SYMBOL) {
            vars_list = wrap_in_list(vars);
        } else {
            vars_list = expr_copy(vars);
        }
        out = solvelinsys_solve_linear_system(expr, vars_list, dom);
        /* A single non-affine equation in several variables is
         * underdetermined but still expressible as explicit rules for
         * one variable (x y == 1 -> {{x -> 1/y}}); try that before the
         * Gröbner specialist, which would refuse it as positive-
         * dimensional.  Only fires for a lone Equal, never a conjunction. */
        if (!out) {
            out = try_single_eq_multivar(expr, vars_list, dom, &opts.poly);
        }
        /* Remaining non-affine systems fall through to the nonlinear
         * specialist, which solves zero-dimensional polynomial systems
         * via a lex Gröbner basis and triangular back-substitution.
         * `expr` and `vars_list` are still owned here and borrowed by
         * the call. */
        if (!out) {
            out = solvenlsys_solve_nonlinear_system(expr, vars_list, dom,
                                                    &opts.poly);
        }
        expr_free(vars_list);
        expr_free(expr);
    } else {
        /* Single-variable path. */
        Expr* var = NULL;
        if (!classify_single_var(vars, &var)) {
            expr_free(expr);
            expr_free(vars_subst);
            expr_free(verify_eq);
            return NULL;
        }

        if (expr->type == EXPR_FUNCTION
            && expr->data.function.head->type == EXPR_SYMBOL
            && expr->data.function.head->data.symbol.name == SYM_Equal
            && expr->data.function.arg_count == 2) {
            out = solvepoly_solve_polynomial_equality(expr, var, dom, &opts.poly);
            /* Polynomial specialist returns NULL when the equation is
             * not a polynomial in `var` -- typically because it carries
             * a transcendental head over var (Sin, Log, ...) or radical
             * subterms (Sqrt, x^(p/q), nested radicals).  Try the
             * inverse-function specialist first (cheap if no peelable
             * head is present), then the radicals specialist. */
            if (!out && opts.inv.enabled
                && solveinv_looks_invertible(expr, var)) {
                out = solveinv_solve_inverse_equality(
                    expr, var, dom, &opts.inv);
            }
            /* Trig canonicalisation pre-pass: handles multi-trig
             * equations that the inverse-function isolator can't peel
             * because more than one trig head over var is present. */
            if (!out && opts.inv.enabled
                && solvetrig_has_trig(expr, var)) {
                out = solvetrig_solve_trig_equality(
                    expr, var, dom, &opts.inv);
            }
            /* Polynomial in a single transcendental kernel g(x) -- e.g.
             * E^(2x)-3E^x+2 (u=E^x) or Log[x]^2-3Log[x]+2 (u=Log[x]):
             * substitute u=g(x), solve the polynomial, unwind each root
             * through the inverse-function peel.  Runs after the single-
             * peel and multi-trig passes have declined. */
            if (!out && opts.inv.enabled) {
                out = solvetrig_solve_poly_in_kernel(
                    expr, var, dom, &opts.inv);
            }
            if (!out) {
                out = solverad_solve_radicals_equality(expr, var, dom);
            }
        }
        expr_free(expr);
    }

solve_finish:
    /* VerifySolutions -> True: drop solutions that fail back-substitution.
     * Runs before the unsubst / numericalise tail, while `out` and
     * `verify_eq` are both in the specialists' substituted / exact form. */
    if (verify_eq) {
        if (out) out = verify_solutions_filter(out, verify_eq);
        expr_free(verify_eq);
        verify_eq = NULL;
    }

    /* Over the Integers, snapshot whether the raw dispatch answer was a
     * positive-dimensional parametric family BEFORE the reals/integers filters
     * can empty it.  A subsequent collapse to `{}` is then not a proof of no
     * solutions but the Complexes/Reals parametric answer failing the integer
     * test for its free parameter -- see the guard below.  Only meaningful
     * when the integer pre-pass declined (a genuine solveint `{}` is a proof). */
    bool parametric_int_family =
        (out && dom && dom->type == EXPR_SYMBOL
         && dom->data.symbol.name == SYM_Integers && !integer_prepass_used)
        ? solution_set_is_parametric(out, vars) : false;

    /* Reals-domain reality filter.  A radical / polynomial / system
     * specialist can hand back Root[] objects (or Root tuples) that are
     * complex -- e.g. the two extraneous branches of Sqrt[x] + 3 x^(1/3) == 5,
     * or the complex roots of an irreducible cubic/quintic.  Over Reals (and
     * its subsets Integers / Rationals) those are not solutions.  Running the
     * drop here, at the shared post-dispatch funnel, covers every specialist
     * uniformly (single-variable and systems alike) and complements the
     * radical solver's own satisfaction check, which handles the default
     * Complexes domain.  Conservative: only *provably* non-real values are
     * dropped, so symbolic / parametric real answers survive. */
    if (out && dom && dom->type == EXPR_SYMBOL
        && (dom->data.symbol.name == SYM_Reals
            || dom->data.symbol.name == SYM_Integers
            || dom->data.symbol.name == SYM_Rationals)) {
        out = filter_reals_solutions(out);
    }

    /* Integers-domain: additionally drop any binding to a provably-non-integer
     * value.  The dedicated integer solver already emits only integers, so this
     * mainly backstops the fall-through polynomial path (e.g. Solve[2 x == 3,
     * x, Integers] -> {}, Solve[x^2 == 2, x, Integers] -> {}). */
    if (out && dom && dom->type == EXPR_SYMBOL
        && dom->data.symbol.name == SYM_Integers) {
        out = filter_integers_solutions(out);
    }

    /* Spurious-`{}` guard (correctness).  When solveint DECLINED and the
     * Complexes-oriented dispatch produced a positive-dimensional parametric
     * answer over the Integers, any empty result reaching here is not a proof
     * of no solutions -- it is the closed-form root / free-parameter family
     * failing the integer test.  Two shapes reach it: a lone multivariable
     * equation (y^2 == x^3 - 2, y == x^2, ...) and an underdetermined linear
     * *system* ({x + 2 y + 3 z == 10, x - y + z == 2} -> parametric in x,
     * which the integer filter emptied to a silent wrong `{}`).  Leave Solve
     * unevaluated instead -- an unproven `{}` would be a silent wrong answer,
     * and these families do have integer points.  (A legitimate solveint `{}`
     * took the `goto solve_finish` above with integer_prepass_used set, so it
     * is never touched here; a fully-determined concrete set that filters to
     * `{}` is a real proof and is not flagged parametric.) */
    if (out && !integer_prepass_used
        && (lone_multivar_int_eq || parametric_int_family)
        && out->type == EXPR_FUNCTION
        && out->data.function.head->type == EXPR_SYMBOL
        && out->data.function.head->data.symbol.name == SYM_List
        && out->data.function.arg_count == 0) {
        expr_free(out);
        expr_free(vars_subst);
        return NULL;
    }

    /* Unsubst pass: if we substituted any compound variables with
     * fresh symbols on the way in, restore the user's original
     * expression in every Rule LHS (and anywhere else the fresh
     * symbol leaked through). */
    if (n_subs > 0 && out) {
        Expr* restored = unsubst_compound_vars(out, subs, n_subs);
        expr_free(out);
        out = restored;
    }
    expr_free(vars_subst);

    /* If we rationalised the input, round-trip the bindings back to
     * floating-point at the original (minimum) precision.  The
     * traversal recurses through List / Rule so {{x -> 1/2}} comes out
     * as {{x -> 0.5}}; the unknown LHS symbols are left alone.  When
     * the original inputs were MPFR at > 53 bits, the result also
     * carries MPFR precision. */
    if (inexact.has_inexact && out) {
        Expr* numeric = common_numericalize_result(out, inexact.min_bits);
        expr_free(out);
        out = numeric;
    }

    return out;  /* evaluator frees res on non-NULL return */
}

/* ------------------------------------------------------------------ *
 *  Init.                                                              *
 * ------------------------------------------------------------------ */

void solve_init(void) {
    symtab_add_builtin("Solve", builtin_solve);
    SymbolDef* def = symtab_get_def("Solve");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve",
        "Solve[expr, vars]\n"
        "\tAttempts to solve the equation or system expr for the\n"
        "\tvariables vars.\n"
        "Solve[expr, vars, dom]\n"
        "\tSolves over the domain dom.  Default Complexes; Reals filters\n"
        "\tdown to real roots via per-degree discriminant and sign tests;\n"
        "\tIntegers finds all integer solutions.  For a single univariate\n"
        "\tpolynomial it filters the roots to concrete integers; a\n"
        "\tpolynomial equation or system with constraints invokes a\n"
        "\tDiophantine engine -- linear via Hermite Normal Form, Pell and\n"
        "\tgeneralised/negative Pell by continued fractions, binary and\n"
        "\tternary quadratic forms, Mordell curves, Thue equations\n"
        "\t(Tzanakis-de Weger), sum-of-three-cubes (Booker; the mod-9\n"
        "\timpossibility globally), exponential Diophantine (Catalan,\n"
        "\tRamanujan-Nagell), and additive power-sum searches\n"
        "\t(meet-in-the-middle).  An empty result {} is always a proof of\n"
        "\tno solution; an out-of-reach input is left unevaluated, never\n"
        "\tguessed.\n"
        "\n"
        "Options:\n"
        "    Cubics              -> False     (radical form for cubics)\n"
        "    Quartics            -> False     (radical form for quartics)\n"
        "    InverseFunctions    -> Automatic (use inverse-function peel)\n"
        "    GeneratedParameters -> C         (head for parameters C[k])\n"
        "    VerifySolutions     -> Automatic (True: drop non-verifying)\n"
        "    Modulus             -> 0         (solve over Z/pZ when p>0)\n"
        "\n"
        "Solves single polynomial equalities, radical equations, linear\n"
        "systems, zero-dimensional nonlinear polynomial systems (via a\n"
        "lexicographic Groebner basis and triangular back-substitution;\n"
        "positive-dimensional systems emit Solve::nsdim and stay\n"
        "unevaluated), a single non-affine equation in several variables\n"
        "(solved for the earliest variable it is polynomial in, e.g.\n"
        "x y == 1 -> {{x -> 1/y}}), equations that are a polynomial in one\n"
        "transcendental kernel g(x) (u = g(x); e.g. E^(2x)-3E^x+2 and\n"
        "Log[x]^2-3Log[x]+2), and -- via the inverse-function specialist --\n"
        "single-variable equations whose outermost dependence is an elementary\n"
        "invertible head (Log, Exp, Sin/Cos/Tan/Cot/Sec/Csc, their\n"
        "hyperbolic counterparts, the inverse trig/hyperbolic forms,\n"
        "and Power[g, n] for integer n >= 2).  Multi-branch heads\n"
        "introduce an integer parameter C[k] wrapped in\n"
        "ConditionalExpression[..., Element[C[k], Integers]].  Emits\n"
        "Solve::ifun the first time inverse functions are used.");

    symtab_set_docstring("Cubics",
        "Cubics is an option for Solve that controls whether cubic\n"
        "\tequations are solved via explicit radical formulas\n"
        "\t(Cubics -> True) or returned as held Root[] objects\n"
        "\t(default Cubics -> False).");
    symtab_set_docstring("Quartics",
        "Quartics is an option for Solve that controls whether quartic\n"
        "\tequations are solved via explicit radical formulas\n"
        "\t(Quartics -> True) or returned as held Root[] objects\n"
        "\t(default Quartics -> False).");
    symtab_set_docstring("GeneratedParameters",
        "GeneratedParameters is an option for Solve specifying the\n"
        "\thead used for fresh integer-parameter symbols introduced by\n"
        "\tthe inverse-function specialist.  Default: C, giving\n"
        "\tC[1], C[2], ...  Only the bare-symbol form is honoured;\n"
        "\tthe Function form is reserved.");
    symtab_set_docstring("InverseFunctions",
        "InverseFunctions is an option for Solve that enables the\n"
        "\tinverse-function specialist for elementary invertible heads\n"
        "\t(Log, Exp, Sin, Cos, Tan, ArcSin, ArcCos, Sinh, ..., and\n"
        "\tinteger Power).  Default: Automatic (enabled).  Setting it\n"
        "\tto False disables the specialist; equations that can only\n"
        "\tbe solved through inversion then return unevaluated.");
    symtab_set_docstring("VerifySolutions",
        "VerifySolutions is an option for Solve.  With VerifySolutions ->\n"
        "\tTrue every returned solution is back-substituted into the\n"
        "\tequation(s) and dropped when PossibleZeroQ proves the residual\n"
        "\tnon-zero; solutions that verify or are undecidable (Root[],\n"
        "\tfree parameters, ConditionalExpression) are kept.  Default:\n"
        "\tAutomatic (per-specialist verification, e.g. radicals).");
    symtab_set_docstring("Modulus",
        "Modulus is an option for Solve.  Solve[poly == 0, x, Modulus -> p]\n"
        "\tsolves a single-variable polynomial equation over the finite\n"
        "\tring Z/pZ by residue enumeration, returning {{x -> r}, ...}\n"
        "\twith r ascending in [0, p).  Supported for 2 <= p <= 100000;\n"
        "\tsystems, multivariable specs, non-polynomial equations, or an\n"
        "\tout-of-range modulus leave Solve unevaluated.");

    solvepoly_init();
    solvelinsys_init();
    solvenlsys_init();
    solverad_init();
    solveinv_init();
    solvetrig_init();
    solvemod_init();
}

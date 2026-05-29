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

#include "attr.h"
#include "common.h"
#include "expr.h"
#include "print.h"
#include "solveinv.h"
#include "solvelinsys.h"
#include "solvepoly.h"
#include "solverad.h"
#include "solvetrig.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Option parsing.                                                    *
 * ------------------------------------------------------------------ */

typedef struct {
    SolvePolyOpts poly;
    SolveInvOpts  inv;
    Expr* dom;             /* borrowed; default = NULL ( = Complexes) */
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
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return is_known_option_name(lhs->data.symbol);
}

/* Returns true iff `e` is the symbol True. */
static bool is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* True iff `e` is the symbol False. */
static bool is_false(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_False;
}

/* Apply a single option rule to `opts`.  Unknown values do not abort
 * (they are silently ignored for now) -- only unknown option *names*
 * are rejected, by is_known_option_name. */
static void apply_option(const Expr* rule, SolveOpts* opts) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;
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
            opts->inv.param_head = rhs->data.symbol;
        }
        return;
    }
    /* VerifySolutions / Assumptions / Method / Modulus: parsed but
     * not yet wired into the polynomial specialist. */
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
        ? lhs->data.symbol : "?";
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
        && vars->data.function.head->data.symbol == SYM_List
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
    const char* h = expr->data.function.head->data.symbol;
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
                const char* h = e->data.function.head->data.symbol;
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
        && vars->data.function.head->data.symbol == SYM_List) {
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
            if (e->data.symbol == subs[i].fresh) {
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
        && vars_in->data.function.head->data.symbol == SYM_List);

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
    *vars_out_owned = expr_new_function(expr_new_symbol("List"), new_args, n);
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
        && vars->data.function.head->data.symbol == SYM_List
        && vars->data.function.arg_count >= 2;
}

/* Wrap a single Expr* (borrowed) in a freshly allocated `List[expr]`.
 * The caller takes ownership of the returned wrapper. */
static Expr* wrap_in_list(Expr* expr) {
    return expr_new_function(expr_new_symbol("List"),
                             (Expr*[]){ expr_copy(expr) }, 1);
}

/* Recognise Equal[Abs[u], 0] / Equal[0, Abs[u]] and rewrite as
 * Equal[u, 0] so the polynomial dispatch can solve u directly.
 * Mirrors Maxima's easy-cases `mabs` shortcut.  Returns a freshly
 * owned Expr* on rewrite, NULL otherwise. */
static Expr* try_abs_zero_rewrite(const Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return NULL;
    if (expr->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (expr->data.function.head->data.symbol != SYM_Equal) return NULL;
    if (expr->data.function.arg_count != 2) return NULL;
    const Expr* lhs = expr->data.function.args[0];
    const Expr* rhs = expr->data.function.args[1];
    const Expr* abs_side = NULL;
    const Expr* zero_side = NULL;
    if (lhs->type == EXPR_FUNCTION
        && lhs->data.function.head->type == EXPR_SYMBOL
        && lhs->data.function.head->data.symbol == SYM_Abs
        && lhs->data.function.arg_count == 1) {
        abs_side = lhs; zero_side = rhs;
    } else if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol == SYM_Abs
        && rhs->data.function.arg_count == 1) {
        abs_side = rhs; zero_side = lhs;
    }
    if (!abs_side) return NULL;
    if (zero_side->type != EXPR_INTEGER || zero_side->data.integer != 0)
        return NULL;
    Expr* u = abs_side->data.function.args[0];
    return expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ expr_copy(u), expr_new_integer(0) }, 2);
}

/* ------------------------------------------------------------------ *
 *  Builtin entry.                                                     *
 * ------------------------------------------------------------------ */

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
            && (a->data.function.head->data.symbol == SYM_Rule
                || a->data.function.head->data.symbol == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol;
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
        NULL
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
    if (expr->type == EXPR_SYMBOL && expr->data.symbol == SYM_True) {
        Expr* empty = expr_new_function(expr_new_symbol("List"), NULL, 0);
        out = expr_new_function(expr_new_symbol("List"),
                                (Expr*[]){ empty }, 1);
        expr_free(expr);
        expr_free(vars_subst);
        return out;
    }
    if (expr->type == EXPR_SYMBOL && expr->data.symbol == SYM_False) {
        out = expr_new_function(expr_new_symbol("List"), NULL, 0);
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
        expr_free(vars_list);
        expr_free(expr);
    } else {
        /* Single-variable path. */
        Expr* var = NULL;
        if (!classify_single_var(vars, &var)) {
            expr_free(expr);
            expr_free(vars_subst);
            return NULL;
        }

        if (expr->type == EXPR_FUNCTION
            && expr->data.function.head->type == EXPR_SYMBOL
            && expr->data.function.head->data.symbol == SYM_Equal
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
            if (!out) {
                out = solverad_solve_radicals_equality(expr, var, dom);
            }
        }
        expr_free(expr);
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
        "\tIntegers further restricts the output to provably concrete\n"
        "\tinteger solutions (Integer / BigInt only -- Rationals, Sqrt[],\n"
        "\tand held Root[] objects are dropped).\n"
        "\n"
        "Options:\n"
        "    Cubics              -> False     (radical form for cubics)\n"
        "    Quartics            -> False     (radical form for quartics)\n"
        "    InverseFunctions    -> Automatic (use inverse-function peel)\n"
        "    GeneratedParameters -> C         (head for parameters C[k])\n"
        "    VerifySolutions     -> Automatic (reserved)\n"
        "\n"
        "Solves single polynomial equalities, radical equations, linear\n"
        "systems, and -- via the inverse-function specialist -- single-\n"
        "variable equations whose outermost dependence is an elementary\n"
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
        "VerifySolutions is an option for Solve that decides whether to\n"
        "\tverify each returned solution by back-substitution.\n"
        "\tDefault: Automatic.  Reserved.");

    solvepoly_init();
    solvelinsys_init();
    solverad_init();
    solveinv_init();
    solvetrig_init();
}

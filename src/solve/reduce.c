/*
 * reduce.c
 *
 * `Reduce` -- front-end and dispatch skeleton (REDUCE_PLAN.md, Phase 0).
 *
 * Mirrors builtin_solve's argument handling: positional `expr [, vars [, dom]]`,
 * a Solve::ivar-style bad-variable diagnostic, and a True/False short-circuit.
 * The input logical combination is normalised into the internal DNF layer
 * (reduce_form.h); Phase 0 returns True/False when the statement decides and
 * leaves everything else unevaluated (NULL).  The per-domain solving engines
 * (equational, sign-diagram, Fourier-Motzkin, CAD, integer, quantifier) are
 * wired on top of this skeleton in later phases.
 */
#include "reduce.h"
#include "reduce_form.h"
#include "reduce_opts.h"
#include "reduce_eq.h"
#include "reduce_univar.h"
#include "reduce_fm.h"
#include "reduce_int.h"
#include "reduce_sys.h"
#include "reduce_cad.h"
#include "reduce_realfn.h"
#include "reduce_realdiag.h"
#include "reduce_qe.h"
#include "reduce_zerodim.h"
#include "reduce_companions.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "attr.h"
#include "expr.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* True iff `var` occurs anywhere in `e` (interned-name pointer compare). */
static bool reduce_contains_var(const Expr* e, const Expr* var) {
    if (!e || !var || var->type != EXPR_SYMBOL) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == var->data.symbol.name;
    if (e->type == EXPR_FUNCTION) {
        if (reduce_contains_var(e->data.function.head, var)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (reduce_contains_var(e->data.function.args[i], var)) return true;
    }
    return false;
}

/* True iff `e` contains a Log[g(var)] or an exponential b^g(var) (var-free
 * base b, so Exp = E^g and general a^t) somewhere.  This is exactly the
 * Log/Exp surface this feature targets; it deliberately excludes forward-trig
 * (Sin/Cos/...) and inverse-trig, which solveinv_looks_invertible also matches
 * but which the Reals sign-diagram engines already own. */
static bool reduce_has_log_or_exp(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Log && e->data.function.arg_count == 1
            && reduce_contains_var(e->data.function.args[0], var)) return true;
        if (hn == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            const Expr* exp_ = e->data.function.args[1];
            if (!reduce_contains_var(base, var) && reduce_contains_var(exp_, var))
                return true;   /* b^g(var): Exp (base E) or a general base */
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_has_log_or_exp(e->data.function.args[i], var)) return true;
    return false;
}

/* True iff the residual `resid` (lhs - rhs) is a genuine Log/Exp transcendental
 * in `var`: it carries a Log[g(var)] or exponential b^g(var) AND is not a
 * positive-integer-degree polynomial in `var`.  The degree test excludes mixed
 * forms such as t^2 + E^t (no closed form) and, together with the Log/Exp
 * restriction, keeps ordinary polynomials -- which have neither a Log/Exp head
 * nor a zero degree -- on the polynomial engine (preserving, e.g., the
 * leading-coefficient case split of Reduce[a x^2 + b x + c == 0, x]). */
static bool reduce_is_transcendental_resid(const Expr* resid, const Expr* var) {
    if (!reduce_has_log_or_exp(resid, var)) return false;
    Expr* e = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy((Expr*)resid), expr_copy((Expr*)var) }, 2));
    bool poly_pos_deg = (e->type == EXPR_INTEGER && e->data.integer >= 1);
    expr_free(e);
    return !poly_pos_deg;
}

/* ------------------------------------------------------------------ *
 *  Option parsing (mirrors solve.c's peeler)                          *
 * ------------------------------------------------------------------ */

/* Recognised Reduce option-name symbols (see Options[Reduce]). */
static bool is_reduce_option_name(const char* s) {
    return s == SYM_Backsubstitution
        || s == SYM_Cubics
        || s == SYM_GeneratedParameters
        || s == SYM_Method
        || s == SYM_Modulus
        || s == SYM_Quartics
        || s == SYM_WorkingPrecision;
}

/* True iff `e` is Rule[opt,_] / RuleDelayed[opt,_] for a recognised name. */
static bool reduce_is_option_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return is_reduce_option_name(lhs->data.symbol.name);
}

/* Apply one option rule to `opts`.  Values borrowed from `res`. */
static void reduce_apply_option(const Expr* rule, ReduceOpts* opts) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    if (name == SYM_Cubics)   { opts->poly.cubics_radical   = is_sym(rhs, SYM_True); return; }
    if (name == SYM_Quartics) { opts->poly.quartics_radical = is_sym(rhs, SYM_True); return; }
    if (name == SYM_GeneratedParameters) {
        if (rhs && rhs->type == EXPR_SYMBOL) opts->param_head = rhs->data.symbol.name;
        return;
    }
    if (name == SYM_Modulus) {
        /* Any non-zero value routes to the modular pre-pass; 0 (the default)
         * stays characteristic 0.  A symbolic or out-of-range modulus makes
         * the pre-pass decline, so Reduce stays unevaluated rather than
         * silently solving over the complexes. */
        bool zero_int = (rhs && rhs->type == EXPR_INTEGER && rhs->data.integer == 0);
        if (rhs && !zero_int) opts->modulus = (Expr*)rhs;
        return;
    }
    if (name == SYM_Backsubstitution) { opts->backsub = is_sym(rhs, SYM_True); return; }
    if (name == SYM_WorkingPrecision) { opts->working_precision = (Expr*)rhs; return; }
    if (name == SYM_Method)           { opts->method = (Expr*)rhs; return; }
}

/* Warn once per distinct form about an unrecognised trailing option. */
static void warn_reduce_bad_option(const Expr* res, const Expr* opt) {
    static uint64_t last_warned_hash = 0;
    uint64_t h = expr_hash((Expr*)res);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    const Expr* lhs = (opt && opt->type == EXPR_FUNCTION
                       && opt->data.function.arg_count == 2)
        ? opt->data.function.args[0] : NULL;
    const char* name = (lhs && lhs->type == EXPR_SYMBOL)
        ? lhs->data.symbol.name : "?";
    fprintf(stderr, "Reduce::optx: Unknown option %s in Reduce.\n", name);
}

/* Warn once per distinct form that the variable spec is invalid.  Mirrors
 * solve.c's warn_ivar, with the Reduce::ivar tag. */
static void warn_reduce_ivar(const Expr* vars) {
    static uint64_t last_warned_hash = 0;
    if (!vars) return;
    uint64_t h = expr_hash(vars);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    char* shown = expr_to_string((Expr*)vars);
    fprintf(stderr, "Reduce::ivar: %s is not a valid variable.\n",
            shown ? shown : "?");
    free(shown);
}

/* Phase 0 accepts a single symbol or a non-empty List of symbols.  (Compound
 * generalised variables, as Solve allows, are a later refinement.) */
static bool reduce_valid_vars(const Expr* vars) {
    if (!vars) return false;
    if (vars->type == EXPR_SYMBOL) return true;
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        if (vars->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < vars->data.function.arg_count; i++)
            if (vars->data.function.args[i]->type != EXPR_SYMBOL) return false;
        return true;
    }
    return false;
}

/* Flatten `vars` (symbol or List of symbols) into a borrowed pointer array. */
static Expr** collect_vars(Expr* vars, int* nv_out) {
    if (vars->type == EXPR_SYMBOL) {
        Expr** v = malloc(sizeof(Expr*));
        v[0] = vars;
        *nv_out = 1;
        return v;
    }
    int nv = (int)vars->data.function.arg_count;
    Expr** v = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++) v[i] = vars->data.function.args[i];
    *nv_out = nv;
    return v;
}

/* ------------------------------------------------------------------ *
 *  builtin                                                            *
 * ------------------------------------------------------------------ */

Expr* builtin_reduce(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing option Rules; the first non-option arg from the end marks
     * the end of the positional args (expr, vars, [dom]).  A trailing
     * Rule[sym,_] whose name is not a known Reduce option is a syntax error. */
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
            if (is_reduce_option_name(name)) { pos_end--; continue; }
            warn_reduce_bad_option(res, a);
            return NULL;
        }
        break;
    }
    if (pos_end < 2 || pos_end > 3) return NULL;

    ReduceOpts opts;
    reduce_opts_default(&opts);
    for (size_t i = pos_end; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (reduce_is_option_arg(a)) reduce_apply_option(a, &opts);
    }

    Expr* expr = res->data.function.args[0];
    Expr* vars = res->data.function.args[1];
    Expr* dom  = (pos_end >= 3) ? res->data.function.args[2] : NULL;

    if (!reduce_valid_vars(vars)) { warn_reduce_ivar(vars); return NULL; }

    /* A list of statements in the `expr` slot is their conjunction, exactly as
     * in Mathematica: Reduce[{e1, e2, ...}, vars] == Reduce[e1 && e2 && ..., vars].
     * Rewrite the List into an And of (copies of) its elements so the rest of
     * the pipeline sees a single logical statement.  `owned_list` is a new owned
     * temporary (like `owned_pre` below) and must be freed on every return path
     * reachable from here; expr_copy is a refcount bump, so freeing it only
     * drops this reference.  Mirrors reduce_int.c's (n==1 ? parts[0] : And) idiom. */
    Expr* owned_list = NULL;
    if (is_head(expr, SYM_List) && expr->data.function.arg_count >= 1) {
        size_t nl = expr->data.function.arg_count;
        Expr** parts = malloc(nl * sizeof(Expr*));
        for (size_t i = 0; i < nl; i++)
            parts[i] = expr_copy(expr->data.function.args[i]);
        owned_list = (nl == 1) ? parts[0]
                   : expr_new_function(expr_new_symbol(SYM_And), parts, nl);
        free(parts);   /* expr_new_function memcpy's the array and adopts the elements */
        expr = owned_list;
    }

    /* Reduce does not hold its args: `expr` arrives already evaluated, so a
     * decidable statement is frequently already True/False here. */
    if (is_sym(expr, SYM_True))  { expr_free(owned_list); return expr_new_symbol(SYM_True); }
    if (is_sym(expr, SYM_False)) { expr_free(owned_list); return expr_new_symbol(SYM_False); }

    /* Phase 7: a top-level Exists / ForAll makes this a quantifier-elimination
     * problem -- eliminate the bound variables and reduce over the remaining
     * free ones.  A NULL / Reals domain proceeds; another explicit domain
     * declines (leaving the input unevaluated). */
    if (is_head(expr, SYM_Exists) || is_head(expr, SYM_ForAll)) {
        Expr* r = reduce_qe_dispatch(expr, dom);   /* borrows expr */
        expr_free(owned_list);
        return r;
    }

    /* Modulus -> p (p != 0): residue enumeration over Z/pZ overrides the
     * domain.  Reuses Solve's modular engine and reformats the result; a
     * symbolic / out-of-range modulus makes it decline (unevaluated). */
    if (opts.modulus) {
        Expr* r = reduce_modular(expr, vars, &opts);   /* borrows expr */
        expr_free(owned_list);
        return r;
    }

    int nv = 0;
    Expr** vlist = collect_vars(vars, &nv);

    /* The equation as it stands BEFORE the realfn rewrite below reassigns
     * `expr`.  A single invertible transcendental equation (Log/Exp/inverse-
     * trig over the sole variable) is force-routed to the Reals sign-diagram
     * engines, which invert real-domain identities (Log[x^2]==2Log[-x] -> x<0)
     * but decline the genuinely invertible ones; those get a Solve-based
     * fallback in the Reals branch, and this borrowed pointer feeds it the
     * un-rewritten residual.  Owned by `res` (or `owned_list`), alive until
     * return. */
    const Expr* orig_eq = expr;

    /* Phase 9: a univariate statement built from Abs / real radicals / Log /
     * bounded-domain inverse-trig / Floor / Mod is a real-domain object.  Route
     * it to the Reals, and rewrite Abs (sign-split), Mod->Floor and isolated
     * integer-part relations away so the sign-diagram engines can consume it. */
    bool force_reals = false;
    Expr* owned_pre = NULL;
    bool reals_dom = (dom == NULL
                      || (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals));
    if (nv == 1 && reals_dom && reduce_stmt_has_realfn(expr, vlist[0])) {
        force_reals = true;
        bool changed = false;
        Expr* pre = reduce_realfn_preprocess(expr, vlist[0], &changed);
        if (pre) { owned_pre = pre; expr = pre; }
    } else if (nv >= 2 && reals_dom && (reduce_stmt_has_piecewise(expr, vlist, nv)
                                        || reduce_stmt_has_radical(expr, vlist, nv))) {
        /* Multivariate piecewise (Max/Min/Abs/Piecewise/...) and/or real radicals
         * (Sqrt[u]): case-split the selectors into polynomial branches and
         * rationalize the radicals into polynomial constraints, then let
         * Fourier-Motzkin / CAD solve over the Reals. */
        force_reals = true;
        bool changed = false;
        Expr* pre = reduce_piecewise_preprocess(expr, vlist, nv, &changed);
        if (pre) { owned_pre = pre; expr = pre; }
    }

    /* Normalise the input logical combination into DNF, then run the
     * constant-atom simplifier. */
    bool ok = true;
    RForm* f = reduce_form_from_expr(expr, vlist, nv, &ok);
    if (!ok) { rform_free(f); free(vlist); expr_free(owned_pre); expr_free(owned_list); return NULL; }

    rform_simplify(f, vlist, nv);

    Expr* out = NULL;
    if (f->is_true)      out = expr_new_symbol(SYM_True);
    else if (f->n == 0)  out = expr_new_symbol(SYM_False);
    else {
        /* Ordering inequalities (< <= > >=, canonicalised to R_LT/R_LE) are only
         * meaningful over an ordered field, so a statement that contains one and
         * was given no explicit domain defaults to the Reals -- matching
         * Mathematica.  Equations and Unequal (!=) remain over the default
         * Complexes. */
        bool default_reals = force_reals;
        if (dom == NULL && !default_reals) {
            for (int i = 0; i < f->n && !default_reals; i++)
                for (int k = 0; k < f->c[i]->n; k++)
                    if (f->c[i]->a[k].rel == R_LT || f->c[i]->a[k].rel == R_LE) {
                        default_reals = true;
                        break;
                    }
        }

        /* Phase 1: a single univariate polynomial equation over Complexes gets
         * its complete (parametric) solution set.  Everything else still needs a
         * per-domain engine not yet wired (Phases 2+) -- leave it unevaluated. */
        bool complexes = !default_reals && ((dom == NULL)
            || (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Complexes));
        bool reals = default_reals || (dom && dom->type == EXPR_SYMBOL
                      && dom->data.symbol.name == SYM_Reals);
        bool integers = (dom && dom->type == EXPR_SYMBOL
                         && (dom->data.symbol.name == SYM_Integers
                             || dom->data.symbol.name == SYM_Rationals));
        if (integers) {
            /* Phase 5: Integers / Rationals via the Diophantine engine. */
            out = reduce_integers(expr, vars, f, vlist, nv, dom, &opts);
        } else if (complexes) {
            /* Over Complexes only equations are meaningful. */
            bool all_eq = true;
            for (int i = 0; i < f->n && all_eq; i++)
                for (int k = 0; k < f->c[i]->n; k++)
                    if (f->c[i]->a[k].rel != R_EQ) { all_eq = false; break; }
            if (nv == 1 && f->n == 1 && f->c[0]->n == 1 && all_eq) {
                const Expr* apoly = f->c[0]->a[0].poly;
                /* A single transcendental equation (Log/Exp/inverse-trig over
                 * the variable) is degree 0 as a polynomial in var, so the
                 * polynomial engine below would only echo `poly == 0`.  Route
                 * it back through Solve -- which combines multi-log residuals
                 * and inverts Exp kernels -- and render the rule-list as a
                 * logical formula.  Falls through when Solve declines. */
                if (reduce_is_transcendental_resid(apoly, vlist[0])) {
                    out = reduce_eq_transcendental(apoly, vlist[0], dom, &opts);
                }
                if (!out) {
                    /* Phase 1: single univariate polynomial equation. */
                    bool ok2 = true;
                    RForm* sol = reduce_eq_univariate(apoly, vlist[0],
                                                      vlist, nv, &ok2, &opts);
                    if (ok2) {
                        rform_simplify(sol, vlist, nv);
                        out = rform_to_expr(sol, vlist, nv);
                    }
                    rform_free(sol);
                }
            } else if (all_eq) {
                /* Phase 4: parametric linear system (declines if non-linear).
                 * Backsubstitution is accepted/echoed at the front-end; the
                 * current linear engine always emits the fully-solved
                 * (grafted) form, which is Reduce's Backsubstitution -> False
                 * default and also what -> True requests, so there is no
                 * behavioural fork to thread here.  A nonlinear (but
                 * zero-dimensional) system falls through to reduce_zerodim. */
                out = reduce_eq_system(f, vlist, nv);
                if (!out) out = reduce_zerodim(f, vlist, nv, false, &opts);
            } else {
                /* Equations plus disequations (!=) over the complexes: the
                 * zero-dimensional engine solves the equations and filters the
                 * != side relations exactly. */
                out = reduce_zerodim(f, vlist, nv, false, &opts);
            }
        } else if (reals && nv == 1) {
            /* Phase 2: any univariate combination of polynomial equations and
             * inequalities over the reals -> sign diagram.  Phase 9: when an atom
             * is a real radical / pole / bounded-domain transcendental (so the
             * exact polynomial engine declines), the general sign diagram takes
             * over. */
            out = reduce_univar(f, vlist[0], vlist, nv, &opts);
            if (!out) out = reduce_univar_general(f, vlist[0], vlist, nv, &opts);
            if (!out && is_head(orig_eq, SYM_Equal)
                && orig_eq->data.function.arg_count == 2) {
                /* The sign-diagram engines declined this univariate equation.
                 * If it is an invertible transcendental (e.g. Log[2t]-Log[x+t]
                 * ==C[1]), re-enter Solve on the original (pre-realfn-rewrite)
                 * residual and render its rule-list as a formula.  A NULL domain
                 * here is the realfn-forced default, not an explicit Reals
                 * request, so solve over Complexes to match how Mathematica
                 * inverts a bare equation. */
                Expr* resid = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(orig_eq->data.function.args[0]),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1),
                                       expr_copy(orig_eq->data.function.args[1]) },
                                   2) }, 2));
                if (reduce_is_transcendental_resid(resid, vlist[0])) {
                    const Expr* solve_dom =
                        (dom && dom->type == EXPR_SYMBOL
                         && dom->data.symbol.name == SYM_Reals) ? dom : NULL;
                    out = reduce_eq_transcendental(resid, vlist[0], solve_dom, &opts);
                }
                expr_free(resid);
            }
        } else if (reals && nv >= 2) {
            /* Phase 3: a multivariate LINEAR system over the reals ->
             * Fourier-Motzkin (declines to NULL if it is not linear).  Phase 6:
             * anything nonlinear falls through to the CAD engine.  A
             * zero-dimensional nonlinear system whose fibres are irrational
             * (which CAD declines) is finally caught by reduce_zerodim: solve
             * the equations exactly and filter the inequalities. */
            out = reduce_fm(f, vlist, nv);
            if (!out) out = reduce_cad(f, vlist, nv);
            if (!out) out = reduce_zerodim(f, vlist, nv, true, &opts);
        }
    }

    rform_free(f);
    free(vlist);
    expr_free(owned_pre);
    expr_free(owned_list);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void reduce_init(void) {
    symtab_add_builtin("Reduce", builtin_reduce);
    SymbolDef* def = symtab_get_def("Reduce");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Reduce",
        "Reduce[expr, vars]\n"
        "\tReduces the statement expr -- a logical combination (&&, ||, !,\n"
        "\tImplies, Xor) of equations (==, !=) and inequalities (<, <=, >,\n"
        "\t>=) -- to a complete, quantifier-free description of its solution\n"
        "\tset for the variables vars.  The default domain is Complexes, or\n"
        "\tReals when expr contains an ordering inequality (ordering is\n"
        "\tundefined over the complexes), so e.g. Reduce[-5 < 3x+7 <= 22, x]\n"
        "\tis solved over the Reals.\n"
        "Reduce[expr, vars, dom]\n"
        "\tReduces over the domain dom: Complexes, Reals, Integers, or\n"
        "\tRationals.\n"
        "A list {e1, e2, ...} in the expr slot is taken as the conjunction\n"
        "e1 && e2 && ... , so Reduce[{x + y == 3, x - y == 1}, {x, y}] is the\n"
        "same as Reduce[x + y == 3 && x - y == 1, {x, y}].\n"
        "\n"
        "Where Solve returns the generic solution of a set of equations as a\n"
        "list of replacement rules and silently drops the degenerate cases,\n"
        "Reduce returns an And/Or tree of equations and inequalities that\n"
        "describes the WHOLE solution set -- every parametric and boundary\n"
        "case kept -- and it solves inequalities over the reals.  A statement\n"
        "that decides returns True or False; an out-of-reach input is left\n"
        "unevaluated, never guessed.\n"
        "\n"
        "Handled so far:\n"
        "  - Complexes: univariate polynomial equations carrying the full\n"
        "    leading-coefficient case tree -- Reduce[a x == b, x] ->\n"
        "    (a != 0 && x == b/a) || (a == 0 && b == 0) -- parametric linear\n"
        "    systems by symbolic Gaussian elimination, and zero-dimensional\n"
        "    nonlinear systems solved exactly (finitely many points).\n"
        "  - Reals, one variable: any Boolean combination of polynomial and\n"
        "    rational-function equations and inequalities, solved as a union\n"
        "    of intervals and points on an exact real-algebraic sign diagram\n"
        "    (denominator roots are breakpoints and are excluded as poles);\n"
        "    plus statements built from Abs, real radicals, Log, bounded\n"
        "    inverse-trig, Floor/Ceiling/Round/Mod/IntegerPart, Min/Max, and\n"
        "    the piecewise heads (Piecewise, Sign, UnitStep, Ramp, Clip,\n"
        "    HeavisideTheta, Boole, ...).\n"
        "  - Reals, several variables: linear systems by Fourier-Motzkin\n"
        "    elimination and nonlinear systems (conics and beyond) by\n"
        "    Cylindrical Algebraic Decomposition (McCallum projection),\n"
        "    including multivariate Abs/Min/Max/piecewise selectors and\n"
        "    square-root radical rationalization.  A zero-dimensional system\n"
        "    (equations pinning finitely many points, e.g. three circles with\n"
        "    sign constraints) is solved exactly and its branches filtered by\n"
        "    the inequalities via the algebraic-number oracle.\n"
        "  - Integers / Rationals: the Solve Diophantine engine, reformatted\n"
        "    as an Or of Ands with Element[C[k], dom] for a free parameter.\n"
        "\n"
        "Options (Options[Reduce]):\n"
        "  Backsubstitution -> False     the linear-system output is fully\n"
        "                                solved (grafted); accepted/echoed.\n"
        "  Cubics -> False               emit radicals (True) or Root[]\n"
        "  Quartics -> False             (False) for irreducible cubic /\n"
        "                                quartic equations.\n"
        "  GeneratedParameters -> C      head of the free parameters C[k]\n"
        "                                for parametric Integers/Rationals\n"
        "                                solutions.\n"
        "  Method -> Automatic           reserved (Automatic is the only\n"
        "                                method).\n"
        "  Modulus -> 0                  a nonzero p solves the equations\n"
        "                                over Z/pZ by residue enumeration,\n"
        "                                overriding the domain.\n"
        "  WorkingPrecision -> Infinity  numeric-fallback tolerance for\n"
        "                                transcendental sign decisions;\n"
        "                                Infinity keeps the exact-first path.\n"
        "\n"
        "Reduce is sound over complete: an undecidable sign, an unsupported\n"
        "construct, or a positive-dimensional nonlinear system over Complexes\n"
        "(no CAD there yet) leaves the input unevaluated rather than risk a\n"
        "wrong formula.");

    reduce_qe_init();
    reduce_companions_init();
}

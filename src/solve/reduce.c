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
#include "reduce_eq.h"
#include "reduce_univar.h"
#include "reduce_fm.h"
#include "reduce_int.h"
#include "reduce_sys.h"
#include "reduce_cad.h"
#include "reduce_realfn.h"
#include "reduce_realdiag.h"

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
    /* Phase 0: positional expr, vars, [dom].  Options are a later phase. */
    if (argc < 2 || argc > 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* vars = res->data.function.args[1];
    Expr* dom  = (argc >= 3) ? res->data.function.args[2] : NULL;

    if (!reduce_valid_vars(vars)) { warn_reduce_ivar(vars); return NULL; }

    /* Reduce does not hold its args: `expr` arrives already evaluated, so a
     * decidable statement is frequently already True/False here. */
    if (is_sym(expr, SYM_True))  return expr_new_symbol(SYM_True);
    if (is_sym(expr, SYM_False)) return expr_new_symbol(SYM_False);

    int nv = 0;
    Expr** vlist = collect_vars(vars, &nv);

    /* Phase 9: a univariate statement built from Abs / real radicals / Log /
     * bounded-domain inverse-trig / Floor / Mod is a real-domain object.  Route
     * it to the Reals, and rewrite Abs (sign-split), Mod->Floor and isolated
     * integer-part relations away so the sign-diagram engines can consume it. */
    bool force_reals = false;
    Expr* owned_pre = NULL;
    if (nv == 1 && (dom == NULL
                    || (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals))
        && reduce_stmt_has_realfn(expr, vlist[0])) {
        force_reals = true;
        bool changed = false;
        Expr* pre = reduce_realfn_preprocess(expr, vlist[0], &changed);
        if (pre) { owned_pre = pre; expr = pre; }
    }

    /* Normalise the input logical combination into DNF, then run the
     * constant-atom simplifier. */
    bool ok = true;
    RForm* f = reduce_form_from_expr(expr, vlist, nv, &ok);
    if (!ok) { rform_free(f); free(vlist); expr_free(owned_pre); return NULL; }

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
            out = reduce_integers(expr, vars, f, vlist, nv, dom);
        } else if (complexes) {
            /* Over Complexes only equations are meaningful. */
            bool all_eq = true;
            for (int i = 0; i < f->n && all_eq; i++)
                for (int k = 0; k < f->c[i]->n; k++)
                    if (f->c[i]->a[k].rel != R_EQ) { all_eq = false; break; }
            if (all_eq) {
                if (nv == 1 && f->n == 1 && f->c[0]->n == 1) {
                    /* Phase 1: single univariate polynomial equation. */
                    bool ok2 = true;
                    RForm* sol = reduce_eq_univariate(f->c[0]->a[0].poly, vlist[0],
                                                      vlist, nv, &ok2);
                    if (ok2) {
                        rform_simplify(sol, vlist, nv);
                        out = rform_to_expr(sol, vlist, nv);
                    }
                    rform_free(sol);
                } else {
                    /* Phase 4: parametric linear system (declines if non-linear). */
                    out = reduce_eq_system(f, vlist, nv);
                }
            }
        } else if (reals && nv == 1) {
            /* Phase 2: any univariate combination of polynomial equations and
             * inequalities over the reals -> sign diagram.  Phase 9: when an atom
             * is a real radical / pole / bounded-domain transcendental (so the
             * exact polynomial engine declines), the general sign diagram takes
             * over. */
            out = reduce_univar(f, vlist[0], vlist, nv);
            if (!out) out = reduce_univar_general(f, vlist[0], vlist, nv);
        } else if (reals && nv >= 2) {
            /* Phase 3: a multivariate LINEAR system over the reals ->
             * Fourier-Motzkin (declines to NULL if it is not linear).  Phase 6:
             * anything nonlinear falls through to the CAD engine. */
            out = reduce_fm(f, vlist, nv);
            if (!out) out = reduce_cad(f, vlist, nv);
        }
    }

    rform_free(f);
    free(vlist);
    expr_free(owned_pre);
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
        "Reduce[expr, vars] reduces the statement expr by solving the equations "
        "and inequalities in it for vars, over the default domain Complexes "
        "(or Reals when the statement contains an ordering inequality), and "
        "returns a complete logical (And/Or) combination of equations and "
        "inequalities describing the whole solution set. Reduce[expr, vars, dom] "
        "works over the domain dom (Complexes, Reals, Integers, or Rationals). "
        "Unlike Solve, Reduce keeps every parametric and degenerate case.");
}

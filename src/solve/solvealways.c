/*
 * solvealways.c
 *
 * Implementation of `SolveAlways[eqns, vars]`.
 *
 * For each equation `lhs == rhs` we form the polynomial `p = lhs - rhs`,
 * treat `p` as a polynomial in `vars` via `CoefficientList[p, vars]`,
 * and require every coefficient to vanish.  The remaining symbols (those
 * appearing in `eqns` but not in `vars`) become the "parameters"; the
 * collected coefficient equations are then passed to `Solve` with the
 * parameters as unknowns.
 *
 * Scope (v1):
 *   - Equations may be a single `Equal[lhs, rhs]`, a `List[Equal, ...]`,
 *     or an `And[Equal, ...]`.
 *   - Variables may be a single symbol or a `List` of symbols.
 *   - Inequations (`Unequal`), disjunctions (`Or`), radicals, and
 *     `Series` strip are NOT handled here; they are deferred.
 */

#include "solvealways.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly/poly.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/*  Tiny utilities                                                     */
/* ------------------------------------------------------------------ */

static bool head_is(const Expr* e, const char* sym) {
    return e
        && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static bool is_list(const Expr* e)  { return head_is(e, SYM_List); }
static bool is_and(const Expr* e)   { return head_is(e, SYM_And); }
static bool is_equal(const Expr* e) { return head_is(e, SYM_Equal); }

static bool is_syntactic_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    return false;
}

/* ------------------------------------------------------------------ */
/*  Diagnostics (one-shot dedupe via hash; mirrors solve.c:warn_ivar)  */
/* ------------------------------------------------------------------ */

static void emit_argt(size_t argc) {
    fprintf(stderr,
        "SolveAlways::argt: SolveAlways called with %zu argument(s); "
        "2 expected.\n", argc);
}

static void emit_eqf(void) {
    fprintf(stderr,
        "SolveAlways::eqf: equations must be Equal[lhs, rhs] (==) or a "
        "list/And of such.\n");
}

static void emit_ivar(const Expr* vars) {
    static uint64_t last_warned_hash = 0;
    if (!vars) return;
    uint64_t h = expr_hash((Expr*)vars);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    fprintf(stderr,
        "SolveAlways::ivar: variable list must be a non-empty symbol or "
        "List of symbols.\n");
}

/* ------------------------------------------------------------------ */
/*  Step 2 -- normalise `vars`                                         */
/* ------------------------------------------------------------------ */

/* On success returns true and either sets *wrap_owned = NULL (borrowing
 * directly into a pre-existing `List[...]` argument) or hands the caller
 * a freshly allocated single-element `List[sym]` to free on exit.
 * `*items` always points at an Expr** of length `*n`. */
static bool normalise_vars(Expr* vars_arg, Expr*** items, size_t* n,
                           Expr** wrap_owned) {
    *wrap_owned = NULL;
    if (is_list(vars_arg)) {
        *n = vars_arg->data.function.arg_count;
        *items = vars_arg->data.function.args;
        if (*n == 0) return false;
        return true;
    }
    if (vars_arg && vars_arg->type == EXPR_SYMBOL) {
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        if (!wrapped) return false;
        wrapped[0] = expr_copy(vars_arg);
        Expr* wrap = expr_new_function(expr_new_symbol(SYM_List), wrapped, 1);
        free(wrapped);
        *wrap_owned = wrap;
        *items = wrap->data.function.args;
        *n = 1;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Step 3 -- normalise `eqns`                                         */
/* ------------------------------------------------------------------ */

typedef enum {
    EQNS_OK,        /* *items / *n populated; caller must validate Equal */
    EQNS_TAUT,      /* sentinel True / empty -- result is {{}}            */
    EQNS_UNSAT,     /* sentinel False        -- result is {}               */
    EQNS_BADSHAPE   /* not Equal / List / And / True / False               */
} EqnsKind;

/* Borrowed-only: *items points into eqns_arg's children (or to a single
 * one-slot scratch pointer the caller owns); no allocation, no copy. */
static EqnsKind normalise_eqns(Expr* eqns_arg, Expr** scratch_slot,
                               Expr*** items, size_t* n) {
    /* Pre-evaluated sentinels (SolveAlways is not HoldAll). */
    if (eqns_arg && eqns_arg->type == EXPR_SYMBOL) {
        if (eqns_arg->data.symbol.name == SYM_True)  return EQNS_TAUT;
        if (eqns_arg->data.symbol.name == SYM_False) return EQNS_UNSAT;
    }
    if (is_equal(eqns_arg)) {
        *scratch_slot = eqns_arg;
        *items = scratch_slot;
        *n = 1;
        return EQNS_OK;
    }
    if (is_list(eqns_arg) || is_and(eqns_arg)) {
        *n = eqns_arg->data.function.arg_count;
        *items = eqns_arg->data.function.args;
        if (*n == 0) return EQNS_TAUT;
        /* Validate each child shape. */
        for (size_t i = 0; i < *n; i++) {
            Expr* x = (*items)[i];
            if (x && x->type == EXPR_SYMBOL) {
                if (x->data.symbol.name == SYM_True)  continue;
                if (x->data.symbol.name == SYM_False) return EQNS_UNSAT;
            }
            if (!is_equal(x)) return EQNS_BADSHAPE;
        }
        return EQNS_OK;
    }
    return EQNS_BADSHAPE;
}

/* ------------------------------------------------------------------ */
/*  Step 4 -- form lhs - rhs and expand                                */
/* ------------------------------------------------------------------ */

/* True iff `ex` is a negative number (integer / bigint / rational), or a Times whose
 * leading factor is one — the exponent shape of a denominator Power[base, ex]. */
static bool sa_is_negative_number(const Expr* ex) {
    if (!ex) return false;
    if (ex->type == EXPR_INTEGER) return ex->data.integer < 0;
    if (ex->type == EXPR_BIGINT)  return mpz_sgn(ex->data.bigint) < 0;
    if (head_is(ex, SYM_Rational) && ex->data.function.arg_count == 2
        && ex->data.function.args[0]->type == EXPR_INTEGER)
        return ex->data.function.args[0]->data.integer < 0;
    if (head_is(ex, SYM_Times) && ex->data.function.arg_count >= 1)
        return sa_is_negative_number(ex->data.function.args[0]);
    return false;
}

/* True iff `e` contains any of the `vars` as a subexpression. */
static bool sa_involves_any_var(const Expr* e, Expr* const* vars, size_t nvars) {
    for (size_t i = 0; i < nvars; i++)
        if (expr_eq((Expr*)e, vars[i])) return true;
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (sa_involves_any_var(e->data.function.head, vars, nvars)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sa_involves_any_var(e->data.function.args[i], vars, nvars)) return true;
    return false;
}

/* True iff `e` contains a denominator that involves one of the vars — a
 * Power[base, neg] with `base` depending on a var.  Such an equation is a genuine
 * rational function of the vars whose denominator must be cleared before
 * coefficient extraction. */
static bool sa_has_var_denominator(const Expr* e, Expr* const* vars, size_t nvars) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2
        && sa_is_negative_number(e->data.function.args[1])
        && sa_involves_any_var(e->data.function.args[0], vars, nvars))
        return true;
    if (sa_has_var_denominator(e->data.function.head, vars, nvars)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sa_has_var_denominator(e->data.function.args[i], vars, nvars)) return true;
    return false;
}

/* Returns an owned, expanded `lhs - rhs`.  The caller must `expr_free`.
 *
 * When `lhs - rhs` is a rational function of the vars (has a var-denominator), the
 * identity `lhs == rhs` for ALL vars holds iff the NUMERATOR of the combined
 * fraction vanishes — so we clear the denominator with Together (robust / FLINT-fast,
 * incl. Q(i)) + Numerator before extracting coefficients.  Feeding the raw rational
 * function to CoefficientList instead yields a spurious no-solution (e.g.
 * SolveAlways[a/(x-1)+b/(x-2) == (3x-5)/((x-1)(x-2)), x] returned {} rather than
 * {a->2, b->1}).  A diff with no var-denominator takes the byte-identical
 * Expand-only path (zero extra cost on the polynomial hot path). */
static Expr* normalise_to_polynomial(const Expr* eq, Expr* const* vars, size_t nvars) {
    Expr* lhs = expr_copy(eq->data.function.args[0]);
    Expr* rhs = expr_copy(eq->data.function.args[1]);
    /* internal_subtract consumes its args. */
    Expr* diff = internal_subtract((Expr*[]){ lhs, rhs }, 2);
    if (sa_has_var_denominator(diff, vars, nvars)) {
        Expr* tog = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
            (Expr*[]){ diff }, 1));                       /* adopts diff */
        Expr* num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Numerator),
            (Expr*[]){ tog }, 1));                        /* adopts tog */
        Expr* expanded = internal_expand((Expr*[]){ num }, 1);
        return eval_and_free(expanded);
    }
    /* Then Expand to canonicalise into a sum of monomials. */
    Expr* expanded = internal_expand((Expr*[]){ diff }, 1);
    /* One evaluate pass lets the evaluator simplify any leftover heads
     * (e.g. Times[-1, Plus[...]] residue). */
    return eval_and_free(expanded);
}

/* ------------------------------------------------------------------ */
/*  Step 5 -- multivariate coefficient extraction                      */
/* ------------------------------------------------------------------ */

static Expr* call_coefficient_list(Expr* poly, Expr* vars_list) {
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(poly);
    args[1] = expr_copy(vars_list);
    Expr* call = expr_new_function(expr_new_symbol(SYM_CoefficientList),
                                   args, 2);
    free(args);
    return eval_and_free(call);
}

/* ------------------------------------------------------------------ */
/*  Step 6 -- flatten nested coefficient lists                          */
/* ------------------------------------------------------------------ */

static void push_owned(Expr*** arr, size_t* n, size_t* cap, Expr* val) {
    if (*n == *cap) {
        *cap = (*cap == 0) ? 8 : (*cap * 2);
        *arr = (Expr**)realloc(*arr, sizeof(Expr*) * (*cap));
    }
    (*arr)[(*n)++] = val;
}

static void flatten_coeffs(const Expr* nested,
                           Expr*** out, size_t* n, size_t* cap) {
    if (is_list(nested)) {
        for (size_t i = 0; i < nested->data.function.arg_count; i++) {
            flatten_coeffs(nested->data.function.args[i], out, n, cap);
        }
        return;
    }
    if (is_syntactic_zero(nested)) return;
    push_owned(out, n, cap, expr_copy((Expr*)nested));
}

/* ------------------------------------------------------------------ */
/*  Step 7 -- compute parameters = symbols(eqns) \ vars                */
/* ------------------------------------------------------------------ */

static bool list_contains_eq(Expr* const* list, size_t n, const Expr* x) {
    for (size_t i = 0; i < n; i++) {
        if (expr_eq(list[i], x)) return true;
    }
    return false;
}

/* Walks every (lhs, rhs) of every equation through `collect_variables`,
 * then drops any symbol that matches a `vars` entry.  Returns owned. */
static void compute_params(Expr* const* eqs, size_t n_eq,
                           Expr* const* var_items, size_t n_vars,
                           Expr*** out, size_t* out_n) {
    /* `collect_variables` doubles `cap` on growth, so a 0 start would
     * stay 0 -- callers always pre-allocate a non-zero buffer. */
    size_t count = 0, cap = 16;
    Expr** buf = (Expr**)malloc(sizeof(Expr*) * cap);
    for (size_t i = 0; i < n_eq; i++) {
        const Expr* eq = eqs[i];
        if (!is_equal(eq)) continue;          /* sentinel True survived */
        collect_variables(eq->data.function.args[0], &buf, &count, &cap);
        collect_variables(eq->data.function.args[1], &buf, &count, &cap);
    }
    /* Filter against vars. */
    size_t kept = 0;
    for (size_t i = 0; i < count; i++) {
        if (list_contains_eq(var_items, n_vars, buf[i])) {
            expr_free(buf[i]);
            buf[i] = NULL;
        } else {
            buf[kept++] = buf[i];
        }
    }
    *out   = buf;
    *out_n = kept;
}

/* ------------------------------------------------------------------ */
/*  Step 9 -- build and evaluate Solve[coeffEqs, params]               */
/* ------------------------------------------------------------------ */

/* Consumes `coeff_eqs[i]` and `params[i]`: they are handed into the new
 * `Solve[...]` tree.  Returns the evaluated result (owned). */
static Expr* call_solve(Expr** coeff_eqs, size_t n_eq,
                        Expr** params, size_t n_params) {
    /* Build List[ Equal[c_i, 0] ... ]. */
    Expr** eq_nodes = (Expr**)malloc(sizeof(Expr*) * n_eq);
    for (size_t i = 0; i < n_eq; i++) {
        Expr** ea = (Expr**)malloc(sizeof(Expr*) * 2);
        ea[0] = coeff_eqs[i];
        ea[1] = expr_new_integer(0);
        eq_nodes[i] = expr_new_function(expr_new_symbol(SYM_Equal), ea, 2);
        free(ea);
    }
    Expr* eq_list = expr_new_function(expr_new_symbol(SYM_List),
                                      eq_nodes, n_eq);
    free(eq_nodes);

    /* Build List[ params... ]. */
    Expr* param_list = expr_new_function(expr_new_symbol(SYM_List),
                                         params, n_params);

    /* Build Solve[eq_list, param_list]. */
    Expr** sa = (Expr**)malloc(sizeof(Expr*) * 2);
    sa[0] = eq_list;
    sa[1] = param_list;
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve), sa, 2);
    free(sa);

    return eval_and_free(call);
}

/* ------------------------------------------------------------------ */
/*  Public entry                                                       */
/* ------------------------------------------------------------------ */

Expr* builtin_solvealways(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) {
        emit_argt(argc);
        return NULL;
    }
    Expr* eqns_arg = res->data.function.args[0];
    Expr* vars_arg = res->data.function.args[1];

    /* ----- Vars ----- */
    Expr* wrap_vars = NULL;
    Expr** var_items = NULL;
    size_t n_vars = 0;
    if (!normalise_vars(vars_arg, &var_items, &n_vars, &wrap_vars)) {
        emit_ivar(vars_arg);
        return NULL;
    }

    /* ----- Eqns ----- */
    Expr* eqns_scratch = NULL;
    Expr** eq_items = NULL;
    size_t n_eq = 0;
    EqnsKind kind = normalise_eqns(eqns_arg, &eqns_scratch,
                                   &eq_items, &n_eq);
    if (kind == EQNS_BADSHAPE) {
        emit_eqf();
        if (wrap_vars) expr_free(wrap_vars);
        return NULL;
    }
    if (kind == EQNS_UNSAT) {
        if (wrap_vars) expr_free(wrap_vars);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    /* EQNS_TAUT (True / empty list) falls through: the eqns loop below
     * just produces no coefficients, and the final "no constraints"
     * branch decides whether to report {} (no params) or {{}}
     * (free params). */

    /* Build a synthetic vars-as-List view for CoefficientList.  If we
     * already own a wrapper (single-symbol case) reuse it; otherwise
     * borrow the user's `vars_arg` (it is already a List). */
    Expr* vars_list_view = wrap_vars ? wrap_vars : vars_arg;

    /* ----- Collect all flattened coefficients across every eq ----- */
    Expr** all_coeffs = NULL;
    size_t n_coeffs = 0, cap_coeffs = 0;
    for (size_t i = 0; i < n_eq; i++) {
        Expr* eq = eq_items[i];
        /* Skip sentinel True survivors from list/And filtering. */
        if (eq && eq->type == EXPR_SYMBOL && eq->data.symbol.name == SYM_True)
            continue;
        Expr* poly  = normalise_to_polynomial(eq, var_items, n_vars);
        Expr* clist = call_coefficient_list(poly, vars_list_view);
        expr_free(poly);
        if (!clist) {
            /* CoefficientList declined: treat poly as a single leaf. */
            Expr* poly2 = normalise_to_polynomial(eq, var_items, n_vars);
            if (!is_syntactic_zero(poly2)) {
                push_owned(&all_coeffs, &n_coeffs, &cap_coeffs, poly2);
            } else {
                expr_free(poly2);
            }
            continue;
        }
        flatten_coeffs(clist, &all_coeffs, &n_coeffs, &cap_coeffs);
        expr_free(clist);
    }

    /* ----- Parameters ----- */
    Expr** params = NULL;
    size_t n_params = 0;
    compute_params(eq_items, n_eq, var_items, n_vars,
                   &params, &n_params);

    if (n_params == 0) {
        /* Per spec: zero-parameter case always returns {} -- there is
         * nothing to report regardless of whether the system is a
         * tautology or a contradiction. */
        for (size_t i = 0; i < n_coeffs; i++) expr_free(all_coeffs[i]);
        free(all_coeffs);
        free(params);
        if (wrap_vars) expr_free(wrap_vars);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    if (n_coeffs == 0) {
        /* No constraints survive but at least one parameter is free:
         * vacuously satisfied for any parameter value -> {{}}. */
        for (size_t i = 0; i < n_params; i++) expr_free(params[i]);
        free(params);
        free(all_coeffs);
        if (wrap_vars) expr_free(wrap_vars);
        Expr* inner = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        Expr** outer = (Expr**)malloc(sizeof(Expr*));
        outer[0] = inner;
        Expr* r = expr_new_function(expr_new_symbol(SYM_List), outer, 1);
        free(outer);
        return r;
    }

    /* ----- Dispatch to Solve ----- */
    Expr* result = call_solve(all_coeffs, n_coeffs, params, n_params);
    free(all_coeffs);
    free(params);
    if (wrap_vars) expr_free(wrap_vars);
    return result;
}

void solvealways_init(void) {
    symtab_add_builtin("SolveAlways", builtin_solvealways);
    SymbolDef* def = symtab_get_def("SolveAlways");
    if (def) def->attributes |= ATTR_PROTECTED;
}

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

#include "sort.h"   /* expr_orderless_sort */
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
#include "message.h"   /* mth_message: Quiet/Check funnel */

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
    mth_message("SolveAlways", "argt", "SolveAlways called with %zu argument(s); 2 expected.", argc);
}

static void emit_eqf(void) {
    mth_message("SolveAlways", "eqf", "equations must be Equal[lhs, rhs] (==) or a list/And of such.");
}

static void emit_ivar(const Expr* vars) {
    static uint64_t last_warned_hash = 0;
    if (!vars) return;
    uint64_t h = expr_hash((Expr*)vars);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    mth_message("SolveAlways", "ivar", "variable list must be a non-empty symbol or List of symbols.");
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
    /* Canonical order, so the rules come out in the order the reference
     * gives them ({a -> 2, b -> 1}), whatever the order of the terms. */
    expr_orderless_sort(buf, kept);
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
/*  Transcendental identities (SolveAlways#1)                          */
/* ------------------------------------------------------------------ */

/* True iff `e` has a var-dependent subterm that is not built from the vars
 * by Plus, Times and integer powers: Sin[x], Log[x], E^x, Sqrt[x], ...
 * collect_variables enrols such a term as one more "parameter", so
 * SolveAlways[a Sin[x] + b == 0, x] became Solve[..., {a, b, Sin[x]}]. */
static bool sa_has_kernel(const Expr* e, Expr* const* vars, size_t nv) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!sa_involves_any_var(e, vars, nv)) return false;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (sa_has_kernel(e->data.function.args[i], vars, nv)) return true;
        return false;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER)
        return sa_has_kernel(e->data.function.args[0], vars, nv);
    return true;
}

typedef struct { Expr** k; Expr** sym; size_t n, cap; } SaKernels;

/* The generator standing for kernel `e` (borrowed), allocating it on first
 * sight. */
static Expr* sa_kernel_symbol(SaKernels* K, const Expr* e) {
    for (size_t i = 0; i < K->n; i++)
        if (expr_eq((Expr*)e, K->k[i])) return expr_copy(K->sym[i]);
    if (K->n == K->cap) {
        K->cap = K->cap ? 2 * K->cap : 8;
        K->k = (Expr**)realloc(K->k, sizeof(Expr*) * K->cap);
        K->sym = (Expr**)realloc(K->sym, sizeof(Expr*) * K->cap);
    }
    char name[64];
    snprintf(name, sizeof name, "SolveAlways`Private`k$%u", (unsigned)(K->n + 1));
    K->k[K->n] = expr_copy((Expr*)e);
    K->sym[K->n] = expr_new_symbol(name);
    return expr_copy(K->sym[K->n++]);
}

/* Is `v` one of the vars? */
static bool sa_is_var(const Expr* v, Expr* const* vars, size_t nv) {
    for (size_t i = 0; i < nv; i++) if (expr_eq((Expr*)v, vars[i])) return true;
    return false;
}

/* Replace every kernel of the expanded exp-log polynomial `e` by a generator.
 * The admitted kernels are E^(c v) (c a nonzero number, v a var) and Log[v]:
 * after TrigToExp + Expand the functions  v^j E^(c v) Log[v]^l  with distinct
 * (j, c, l) are linearly independent over the constants, so an identity holds
 * iff every coefficient vanishes.  Any other kernel sets *bad. */
static Expr* sa_subst_kernels(const Expr* e, Expr* const* vars, size_t nv,
                              SaKernels* K, bool* bad) {
    if (*bad) return NULL;
    if (!e || e->type != EXPR_FUNCTION || !sa_involves_any_var(e, vars, nv))
        return expr_copy((Expr*)e);
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == SYM_E) {
        const Expr* g = e->data.function.args[1];
        bool ok = sa_is_var(g, vars, nv);
        if (!ok && head_is(g, SYM_Times) && g->data.function.arg_count == 2
            && sa_is_var(g->data.function.args[1], vars, nv)) {
            Expr* nq = eval_and_free(expr_new_function(expr_new_symbol(SYM_NumericQ),
                (Expr*[]){ expr_copy(g->data.function.args[0]) }, 1));
            ok = nq->type == EXPR_SYMBOL && nq->data.symbol.name == SYM_True;
            expr_free(nq);
        }
        if (!ok) { *bad = true; return NULL; }
        return sa_kernel_symbol(K, e);
    }
    if (head_is(e, SYM_Log) && e->data.function.arg_count == 1
        && sa_is_var(e->data.function.args[0], vars, nv))
        return sa_kernel_symbol(K, e);
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)
        || (head_is(e, SYM_Power) && e->data.function.arg_count == 2
            && e->data.function.args[1]->type == EXPR_INTEGER)) {
        size_t n = e->data.function.arg_count;
        Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            a[i] = sa_subst_kernels(e->data.function.args[i], vars, nv, K, bad);
            if (*bad) { for (size_t j = 0; j < i; j++) expr_free(a[j]); free(a); return NULL; }
        }
        Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
        free(a);
        return r;
    }
    *bad = true;
    return NULL;
}

/* SolveAlways over exp-log identities: rewrite each lhs - rhs with TrigToExp,
 * expand, name the kernels, and require every coefficient of the resulting
 * polynomial in vars + kernels to vanish.  Returns the Solve result, or NULL
 * (with *declined) when a kernel outside the admitted family appears -- the
 * call then stays unevaluated rather than solving for a function of x. */
static Expr* sa_solve_transcendental(Expr** eqs, size_t n_eq, Expr* const* vars,
                                     size_t nv, bool* declined) {
    *declined = false;
    SaKernels K = { NULL, NULL, 0, 0 };
    Expr** polys = (Expr**)malloc(sizeof(Expr*) * (n_eq ? n_eq : 1));
    size_t np = 0;
    bool bad = false;
    for (size_t i = 0; i < n_eq && !bad; i++) {
        Expr* eq = eqs[i];
        if (!is_equal(eq)) continue;
        Expr* diff = internal_subtract((Expr*[]){ expr_copy(eq->data.function.args[0]),
                                                  expr_copy(eq->data.function.args[1]) }, 2);
        Expr* te = eval_and_free(expr_new_function(expr_new_symbol(SYM_TrigToExp),
                                                   (Expr*[]){ diff }, 1));
        if (sa_has_var_denominator(te, vars, nv)) {
            Expr* tog = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                                                        (Expr*[]){ te }, 1));
            te = eval_and_free(expr_new_function(expr_new_symbol(SYM_Numerator),
                                                 (Expr*[]){ tog }, 1));
        }
        Expr* ex = eval_and_free(internal_expand((Expr*[]){ te }, 1));
        Expr* sub = sa_subst_kernels(ex, vars, nv, &K, &bad);
        expr_free(ex);
        if (!bad) polys[np++] = eval_and_free(sub);
    }
    Expr* result = NULL;
    if (!bad) {
        /* generators: vars, then kernels */
        size_t ng = nv + K.n;
        Expr** gv = (Expr**)malloc(sizeof(Expr*) * (ng ? ng : 1));
        for (size_t i = 0; i < nv; i++) gv[i] = expr_copy(vars[i]);
        for (size_t i = 0; i < K.n; i++) gv[nv + i] = expr_copy(K.sym[i]);
        Expr* glist = expr_new_function(expr_new_symbol(SYM_List), gv, ng);
        free(gv);
        Expr** coeffs = NULL;
        size_t nc = 0, cc = 0;
        for (size_t i = 0; i < np && !bad; i++) {
            if (sa_has_kernel(polys[i], vars, nv)) { bad = true; break; }
            Expr* cl = call_coefficient_list(polys[i], glist);
            if (!cl || !is_list(cl)) { expr_free(cl); bad = true; break; }
            flatten_coeffs(cl, &coeffs, &nc, &cc);
            expr_free(cl);
        }
        /* parameters: the symbols of the generator polynomials less the generators */
        size_t count = 0, cap = 16;
        Expr** buf = (Expr**)malloc(sizeof(Expr*) * cap);
        for (size_t i = 0; i < np; i++) collect_variables(polys[i], &buf, &count, &cap);
        size_t kept = 0;
        for (size_t i = 0; i < count; i++) {
            if (list_contains_eq(glist->data.function.args, ng, buf[i])
                || sa_involves_any_var(buf[i], vars, nv)) expr_free(buf[i]);
            else buf[kept++] = buf[i];
        }
        expr_orderless_sort(buf, kept);
        expr_free(glist);
        if (bad || kept == 0) {
            /* No parameters: the identity either holds (every coefficient is
             * certifiably zero: {{}}), certainly fails (a coefficient is a
             * nonzero number: {}), or is undecided (decline).  Sin[x]^2 +
             * Cos[x]^2 == 1 is {{}}, Sin[x] == Sin[x] + 1 is {}. */
            int verdict = 1;                     /* 1 holds, 0 fails, -1 undecided */
            for (size_t i = 0; !bad && i < nc && verdict != 0; i++) {
                Expr* z = eval_and_free(expr_new_function(expr_new_symbol(SYM_PossibleZeroQ),
                    (Expr*[]){ expr_copy(coeffs[i]) }, 1));
                bool zero = z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_True;
                expr_free(z);
                if (zero) continue;
                Expr* nq = eval_and_free(expr_new_function(expr_new_symbol(SYM_NumericQ),
                    (Expr*[]){ expr_copy(coeffs[i]) }, 1));
                bool num = nq->type == EXPR_SYMBOL && nq->data.symbol.name == SYM_True;
                expr_free(nq);
                verdict = num ? 0 : -1;
            }
            for (size_t i = 0; i < nc; i++) expr_free(coeffs[i]);
            free(coeffs);
            for (size_t i = 0; i < kept; i++) expr_free(buf[i]);
            free(buf);
            if (!bad && verdict == 0)
                result = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            else if (!bad && verdict == 1) {
                Expr* inner = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
                result = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
            } else if (!bad)
                bad = true;                      /* undecided: stay unevaluated */
        } else if (nc == 0) {
            for (size_t i = 0; i < kept; i++) expr_free(buf[i]);
            free(buf);
            free(coeffs);
            Expr* inner = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            result = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
        } else {
            result = call_solve(coeffs, nc, buf, kept);
            free(coeffs);
            free(buf);
        }
    }
    for (size_t i = 0; i < np; i++) expr_free(polys[i]);
    free(polys);
    for (size_t i = 0; i < K.n; i++) { expr_free(K.k[i]); expr_free(K.sym[i]); }
    free(K.k); free(K.sym);
    if (bad) *declined = true;
    return result;
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

    /* ----- Transcendental (exp-log-trig) identities ----- */
    {
        bool any_kernel = false;
        for (size_t i = 0; i < n_eq && !any_kernel; i++)
            if (is_equal(eq_items[i])
                && (sa_has_kernel(eq_items[i]->data.function.args[0], var_items, n_vars)
                    || sa_has_kernel(eq_items[i]->data.function.args[1], var_items, n_vars)))
                any_kernel = true;
        if (any_kernel) {
            bool declined = false;
            Expr* r = sa_solve_transcendental(eq_items, n_eq, var_items, n_vars, &declined);
            if (wrap_vars) expr_free(wrap_vars);
            return declined ? NULL : r;
        }
    }

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

/*
 * reduce_zerodim.c
 *
 * Zero-dimensional polynomial-system engine for `Reduce` and `Solve`.
 * See reduce_zerodim.h for the algorithm and contract.
 */
#include "reduce_zerodim.h"

#include "eval.h"
#include "sym_names.h"
#include "reduce_real_util.h"   /* rru_sign_of */
#include "flint_qqbar.h"        /* flint_qqbar_is_real / _equal */

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* True if `e` (or any subexpression) is one of the solve variables -- used to
 * reject a parametric (positive-dimensional) Solve answer. */
static bool contains_solve_var(const Expr* e, Expr** vars, int nv) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < nv; i++)
            if (vars[i]->type == EXPR_SYMBOL
                && vars[i]->data.symbol.name == e->data.symbol.name)
                return true;
        return false;
    }
    if (e->type == EXPR_FUNCTION) {
        if (contains_solve_var(e->data.function.head, vars, nv)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_solve_var(e->data.function.args[i], vars, nv)) return true;
    }
    return false;
}

/* Index of `sym` in vars[0..nv-1], or -1. */
static int var_index(const Expr* sym, Expr** vars, int nv) {
    if (!sym || sym->type != EXPR_SYMBOL) return -1;
    for (int i = 0; i < nv; i++)
        if (vars[i]->type == EXPR_SYMBOL
            && vars[i]->data.symbol.name == sym->data.symbol.name)
            return i;
    return -1;
}

/* Substitute vars[i] -> val[i] into `poly` and evaluate to a constant. */
static Expr* subst_branch(const Expr* poly, Expr** vars, Expr** val, int nv) {
    Expr** rules = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++)
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy(vars[i]), expr_copy(val[i]) }, 2);
    Expr* rlist = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nv);
    free(rules);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rlist }, 2));
}

/* ------------------------------------------------------------------ *
 *  Branch accumulator: an owned array of solution branches, each a    *
 *  List[Rule[var, val], ...] over all nv variables.                   *
 * ------------------------------------------------------------------ */

typedef struct { Expr** b; int n, cap; } Branches;

static void branches_push(Branches* B, Expr* branch /*owned*/) {
    /* Structural dedup: a solution reached via two conjuncts appears once. */
    for (int i = 0; i < B->n; i++)
        if (expr_eq(B->b[i], branch)) { expr_free(branch); return; }
    if (B->n == B->cap) {
        B->cap = B->cap ? B->cap * 2 : 4;
        B->b = realloc(B->b, (size_t)B->cap * sizeof(Expr*));
    }
    B->b[B->n++] = branch;
}

static void branches_free(Branches* B) {
    for (int i = 0; i < B->n; i++) expr_free(B->b[i]);
    free(B->b);
}

/* ------------------------------------------------------------------ *
 *  Per-conjunct solver                                                *
 * ------------------------------------------------------------------ */

/* Decide whether the side relation `a` (rel over its poly) holds at branch
 * value V = poly(branch).  Returns 1 (holds), 0 (fails), -1 (undecidable). */
static int side_holds(RRel rel, const Expr* V) {
    switch (rel) {
        case R_LT: { int s = rru_sign_of(V); return (s == -2) ? -1 : (s < 0); }
        case R_LE: { int s = rru_sign_of(V); return (s == -2) ? -1 : (s <= 0); }
        case R_NE: {
            Expr* zero = expr_new_integer(0);
            int e = flint_qqbar_equal(V, zero);   /* 1 eq, 0 ne, -1 undecided */
            expr_free(zero);
            if (e == -1) return -1;
            return e == 0;                        /* holds iff not equal */
        }
        default: return -1;                       /* R_EQ/R_ELEM never here */
    }
}

/* Solve the zero-dimensional equational subsystem of `conj` over the complexes,
 * filter by its side relations (and realness when `reals`), and append the
 * surviving branches to `B`.  Returns 1 (handled -- branches appended, possibly
 * none) or 0 (DECLINE). */
static int zerodim_conj(const RConj* conj, Expr** vars, int nv, bool reals,
                        const ReduceOpts* opts, Branches* B) {
    if (conj->is_false) return 1;                 /* no solutions, handled */

    /* Partition atoms into equations and side relations. */
    Expr** eqs = malloc((size_t)(conj->n ? conj->n : 1) * sizeof(Expr*));
    int neq = 0;
    const RAtom** side = malloc((size_t)(conj->n ? conj->n : 1) * sizeof(RAtom*));
    int nside = 0;
    int rc = 1;                                   /* 1 = handled, 0 = decline */

    for (int k = 0; k < conj->n; k++) {
        const RAtom* a = &conj->a[k];
        if (a->rel == R_ELEM || a->nonconst_denom) { rc = 0; goto cleanup; }
        if (a->rel == R_EQ) {
            eqs[neq++] = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy(a->poly), expr_new_integer(0) }, 2);
        } else {
            side[nside++] = a;                    /* R_NE / R_LT / R_LE */
        }
    }
    if (neq == 0) { rc = 0; goto cleanup; }       /* not pinned -> not our job */

    /* Build and evaluate Solve[And(eqs), List(vars)] over the complexes. */
    Expr* eqconj = (neq == 1) ? eqs[0]
        : expr_new_function(expr_new_symbol(SYM_And), eqs, (size_t)neq);
    Expr** vcopy = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++) vcopy[i] = expr_copy(vars[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vcopy, (size_t)nv);
    free(vcopy);
    Expr* base[2] = { eqconj, vlist };            /* ownership passes on */
    neq = 0;                                       /* eqs[] now owned by eqconj */
    Expr* call = reduce_opts_build_solve(base, 2, opts);
    Expr* sols = eval_and_free(call);

    if (!is_head(sols, SYM_List)) { expr_free(sols); rc = 0; goto cleanup; }

    /* Each element is a solution branch: List[Rule[var, val], ...]. */
    Expr** val = malloc((size_t)nv * sizeof(Expr*));
    for (size_t t = 0; t < sols->data.function.arg_count && rc; t++) {
        Expr* row = sols->data.function.args[t];
        for (int i = 0; i < nv; i++) val[i] = NULL;

        /* Read the branch into val[] (indexed by variable), rejecting a
         * parametric / underdetermined / conditional answer. */
        bool clean = is_head(row, SYM_List);
        for (size_t j = 0; clean && j < row->data.function.arg_count; j++) {
            Expr* rule = row->data.function.args[j];
            if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) { clean = false; break; }
            int vi = var_index(rule->data.function.args[0], vars, nv);
            Expr* rv = rule->data.function.args[1];
            if (vi < 0 || val[vi] != NULL
                || is_head(rv, SYM_ConditionalExpression)
                || contains_solve_var(rv, vars, nv)) { clean = false; break; }
            val[vi] = rv;                         /* borrowed from sols */
        }
        if (clean)
            for (int i = 0; i < nv; i++) if (!val[i]) { clean = false; break; }
        if (!clean) { rc = 0; break; }            /* positive-dimensional -> decline */

        /* Filter: realness (Reals domain) then each side relation. */
        bool keep = true;
        if (reals) {
            for (int i = 0; i < nv && keep; i++) {
                int r = flint_qqbar_is_real(val[i]);
                if (r == -1) { rc = 0; break; }   /* undecidable -> decline */
                if (r == 0) keep = false;         /* non-real -> drop branch */
            }
        }
        for (int s = 0; s < nside && keep && rc; s++) {
            Expr* V = subst_branch(side[s]->poly, vars, val, nv);
            int h = side_holds(side[s]->rel, V);
            expr_free(V);
            if (h == -1) { rc = 0; break; }       /* undecidable -> decline */
            if (h == 0) keep = false;             /* violated -> drop branch */
        }
        if (!rc) break;
        if (!keep) continue;

        /* Surviving branch: emit List[Rule[var_i, val_i], ...] in var order. */
        Expr** rules = malloc((size_t)nv * sizeof(Expr*));
        for (int i = 0; i < nv; i++)
            rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_copy(vars[i]), expr_copy(val[i]) }, 2);
        branches_push(B, expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nv));
        free(rules);
    }
    free(val);
    expr_free(sols);

cleanup:
    for (int i = 0; i < neq; i++) expr_free(eqs[i]);
    free(eqs);
    free((void*)side);
    return rc;
}

/* Fill `B` from every conjunct of `F`; returns false (and frees `B`) on any
 * decline. */
static bool collect_branches(const RForm* F, Expr** vars, int nv, bool reals,
                             const ReduceOpts* opts, Branches* B) {
    memset(B, 0, sizeof(*B));
    for (int i = 0; i < F->n; i++)
        if (!zerodim_conj(F->c[i], vars, nv, reals, opts, B)) {
            branches_free(B);
            return false;
        }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Reduce formatter: Or of And(var == val)                            *
 * ------------------------------------------------------------------ */

/* One branch (List[Rule[var,val]]) -> And(Equal[var,val], ...). */
static Expr* branch_to_and(const Expr* branch) {
    size_t n = branch->data.function.arg_count;
    Expr** parts = malloc(n * sizeof(Expr*));
    for (size_t j = 0; j < n; j++) {
        Expr* rule = branch->data.function.args[j];
        parts[j] = expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(rule->data.function.args[0]),
                       expr_copy(rule->data.function.args[1]) }, 2);
    }
    Expr* out = (n == 1) ? parts[0]
        : expr_new_function(expr_new_symbol(SYM_And), parts, n);
    free(parts);
    return out;
}

Expr* reduce_zerodim(const RForm* F, Expr** vars, int nv, bool reals,
                     const ReduceOpts* opts) {
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    Branches B;
    if (!collect_branches(F, vars, nv, reals, opts, &B)) return NULL;

    Expr* out;
    if (B.n == 0) out = expr_new_symbol(SYM_False);
    else {
        Expr** disj = malloc((size_t)B.n * sizeof(Expr*));
        for (int i = 0; i < B.n; i++) disj[i] = branch_to_and(B.b[i]);
        out = (B.n == 1) ? disj[0]
            : expr_new_function(expr_new_symbol(SYM_Or), disj, (size_t)B.n);
        free(disj);
    }
    branches_free(&B);
    return eval_and_free(out);                    /* flatten And/Or */
}

/* ------------------------------------------------------------------ *
 *  Solve formatter: List of solution rule-lists                       *
 * ------------------------------------------------------------------ */

/* Collect vars (symbol or List of symbols) into a borrowed pointer array. */
static Expr** collect_vars(const Expr* vars, int* nv_out) {
    if (vars->type == EXPR_SYMBOL) {
        Expr** v = malloc(sizeof(Expr*)); v[0] = (Expr*)vars; *nv_out = 1; return v;
    }
    int nv = (int)vars->data.function.arg_count;
    Expr** v = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++) v[i] = vars->data.function.args[i];
    *nv_out = nv;
    return v;
}

/* Does `e` contain an ordering inequality head (Less/LessEqual/Greater/
 * GreaterEqual/Inequality)?  Such a relation is only meaningful over an ordered
 * field, so its presence forces the Reals when no domain was named. */
static bool has_ordering(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Less || n == SYM_LessEqual || n == SYM_Greater
            || n == SYM_GreaterEqual || n == SYM_Inequality) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_ordering(e->data.function.args[i])) return true;
    return false;
}

/* Does `e` (an And / bare relation) carry at least one inequality or
 * disequation conjunct -- the constraint that the plain system specialists
 * refuse?  Without one, this engine yields to Solve's own dispatch. */
static bool has_side_constraint(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Less || n == SYM_LessEqual || n == SYM_Greater
            || n == SYM_GreaterEqual || n == SYM_Inequality || n == SYM_Unequal)
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_side_constraint(e->data.function.args[i])) return true;
    return false;
}

Expr* reduce_zerodim_solve(const Expr* expr, const Expr* vars, const Expr* dom,
                           const SolvePolyOpts* poly) {
    /* Only step in for the mixed equation + constraint case; a pure equation
     * system is left to Solve's own specialists (byte-for-byte unchanged). */
    if (!has_side_constraint(expr)) return NULL;

    /* Domain: Reals / Complexes explicit; else Reals iff an ordering
     * inequality is present (undefined over the complexes otherwise). */
    bool reals;
    if (dom && dom->type == EXPR_SYMBOL) {
        if (dom->data.symbol.name == SYM_Reals) reals = true;
        else if (dom->data.symbol.name == SYM_Complexes) reals = false;
        else return NULL;                          /* Integers/Rationals: not here */
    } else {
        reals = has_ordering(expr);
    }

    int nv = 0;
    Expr** vlist = collect_vars(vars, &nv);
    if (nv == 0) { free(vlist); return NULL; }

    bool ok = true;
    RForm* F = reduce_form_from_expr(expr, vlist, nv, &ok);
    if (!ok) { rform_free(F); free(vlist); return NULL; }
    rform_simplify(F, vlist, nv);

    ReduceOpts opts;
    reduce_opts_default(&opts);
    if (poly) opts.poly = *poly;

    Expr* out = NULL;
    if (F->is_true) {
        /* A tautology among constraints alone is not a discrete solution set;
         * decline rather than fabricate {{}}. */
        out = NULL;
    } else if (F->n == 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);   /* {} */
    } else {
        Branches B;
        if (collect_branches(F, vlist, nv, reals, &opts, &B)) {
            Expr** rows = B.n ? malloc((size_t)B.n * sizeof(Expr*)) : NULL;
            for (int i = 0; i < B.n; i++) rows[i] = expr_copy(B.b[i]);
            out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)B.n);
            free(rows);
            branches_free(&B);
        }
    }

    rform_free(F);
    free(vlist);
    return out;
}

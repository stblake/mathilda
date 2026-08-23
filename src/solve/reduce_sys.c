/*
 * reduce_sys.c
 *
 * Parametric linear-system engine for `Reduce` (REDUCE_PLAN.md, Phase 4).
 * See reduce_sys.h.  Coefficients are rational-function Exprs; arithmetic and
 * simplification reuse the evaluator (Plus/Times/Together/Coefficient/Solve).
 */
#include "reduce_sys.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

#define SYS_MAX_ROWS 64
#define SYS_MAX_VARS 24

/* ------------------------------------------------------------------ *
 *  Rational-function scalar helpers                                   *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}
static bool is_num(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: return is_head(e, SYM_Rational);
        default: return false;
    }
}
static bool is_lit_zero(const Expr* e) { return e->type == EXPR_INTEGER && e->data.integer == 0; }

/* r := Together[e], evaluated (single canonical fraction).  Consumes `e`. */
static Expr* rf(Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
        (Expr*[]){ e }, 1));
}
static Expr* rf_add(const Expr* a, const Expr* b) {
    return rf(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
}
static Expr* rf_mul(const Expr* a, const Expr* b) {
    return rf(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
}
static Expr* rf_neg(const Expr* a) {
    return rf(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)a) }, 2));
}
static Expr* rf_div(const Expr* a, const Expr* b) {
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)b), expr_new_integer(-1) }, 2));
    Expr* r = rf_mul(a, inv);
    expr_free(inv);
    return r;
}
/* Provably zero: Together collapses it to the integer 0. */
static bool rf_is_zero(const Expr* e) {
    Expr* t = rf(expr_copy((Expr*)e));
    bool z = is_lit_zero(t);
    expr_free(t);
    return z;
}
static bool rf_is_nonzero_const(const Expr* e) {
    return is_num(e) && !is_lit_zero(e);
}

/* ------------------------------------------------------------------ *
 *  Solution representation (DNF of cases)                             *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr** cond; int ncond;              /* param relation Exprs (owned) */
    Expr** var;  Expr** val; int nass;   /* assignments var == val (owned) */
} LCase;
typedef struct { LCase* c; int n, cap; } LSol;

static void lcase_free(LCase* c) {
    for (int i = 0; i < c->ncond; i++) expr_free(c->cond[i]);
    for (int i = 0; i < c->nass; i++) { expr_free(c->var[i]); expr_free(c->val[i]); }
    free(c->cond); free(c->var); free(c->val);
}
static LSol* lsol_new(void) { return calloc(1, sizeof(LSol)); }
static void lsol_free(LSol* s) {
    if (!s) return;
    for (int i = 0; i < s->n; i++) lcase_free(&s->c[i]);
    free(s->c); free(s);
}
static void lsol_push(LSol* s, LCase c) {
    if (s->n == s->cap) { s->cap = s->cap ? s->cap * 2 : 4;
                          s->c = realloc(s->c, (size_t)s->cap * sizeof(LCase)); }
    s->c[s->n++] = c;
}
static LSol* lsol_true(void) {
    LSol* s = lsol_new(); LCase c; memset(&c, 0, sizeof c); lsol_push(s, c); return s;
}
static LSol* lsol_false(void) { return lsol_new(); }   /* no cases */

static void lcase_add_cond(LCase* c, Expr* rel /*owned*/) {
    c->cond = realloc(c->cond, (size_t)(c->ncond + 1) * sizeof(Expr*));
    c->cond[c->ncond++] = rel;
}
static void lcase_add_assign(LCase* c, Expr* var /*owned*/, Expr* val /*owned*/) {
    c->var = realloc(c->var, (size_t)(c->nass + 1) * sizeof(Expr*));
    c->val = realloc(c->val, (size_t)(c->nass + 1) * sizeof(Expr*));
    c->var[c->nass] = var; c->val[c->nass] = val; c->nass++;
}
/* OR: move all cases of b into a; frees b. */
static LSol* lsol_or(LSol* a, LSol* b) {
    for (int i = 0; i < b->n; i++) lsol_push(a, b->c[i]);
    b->n = 0; lsol_free(b);
    return a;
}
/* AND a param relation into every case (copied). */
static void lsol_and_cond(LSol* s, const Expr* rel) {
    for (int i = 0; i < s->n; i++) lcase_add_cond(&s->c[i], expr_copy((Expr*)rel));
}

/* ------------------------------------------------------------------ *
 *  Row helpers                                                        *
 * ------------------------------------------------------------------ */

/* A row is an Expr*[nv+1]: coeffs of vars[0..nv-1] then the constant term;
 * the equation is  sum coeff*var + const == 0. */
static void row_free(Expr** r, int nv) {
    for (int i = 0; i <= nv; i++) expr_free(r[i]);
    free(r);
}
static void rows_free(Expr*** R, int nR, int nv) {
    for (int i = 0; i < nR; i++) row_free(R[i], nv);
    free(R);
}

/* Substitute param -> val throughout every row entry (then Together). */
static Expr*** rows_subst(Expr*** R, int nR, int nv, const Expr* param, const Expr* val) {
    Expr*** d = malloc((size_t)(nR ? nR : 1) * sizeof(Expr**));
    for (int i = 0; i < nR; i++) {
        d[i] = malloc((size_t)(nv + 1) * sizeof(Expr*));
        for (int k = 0; k <= nv; k++) {
            Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_copy((Expr*)param), expr_copy((Expr*)val) }, 2);
            Expr* ra = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                (Expr*[]){ expr_copy(R[i][k]), rule }, 2);
            d[i][k] = rf(ra);
        }
    }
    return d;
}

/* Collect distinct free symbols of `e` that are not solve variables. */
static void collect_params(const Expr* e, Expr** vars, int nv, Expr*** out, int* n) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        for (int j = 0; j < nv; j++)
            if (vars[j]->type == EXPR_SYMBOL
                && vars[j]->data.symbol.name == e->data.symbol.name) return;
        for (int i = 0; i < *n; i++) if (expr_eq((*out)[i], e)) return;
        *out = realloc(*out, (size_t)(*n + 1) * sizeof(Expr*));
        (*out)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        collect_params(e->data.function.head, vars, nv, out, n);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_params(e->data.function.args[i], vars, nv, out, n);
    }
}

/* ------------------------------------------------------------------ *
 *  The recursive solver                                               *
 * ------------------------------------------------------------------ */

static LSol* lin_solve(Expr*** R, int nR, bool* remaining, int nv, Expr** vars, bool* ok);

/* Back-substitute a case's assignments into `xexpr`, then record var == value. */
static void graft(LSol* sub, const Expr* var, const Expr* xexpr) {
    for (int i = 0; i < sub->n; i++) {
        LCase* c = &sub->c[i];
        Expr* xval;
        if (c->nass == 0) {
            xval = rf(expr_copy((Expr*)xexpr));
        } else {
            Expr** rules = malloc((size_t)c->nass * sizeof(Expr*));
            for (int k = 0; k < c->nass; k++)
                rules[k] = expr_new_function(expr_new_symbol(SYM_Rule),
                    (Expr*[]){ expr_copy(c->var[k]), expr_copy(c->val[k]) }, 2);
            Expr* rlist = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)c->nass);
            free(rules);
            xval = rf(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                (Expr*[]){ expr_copy((Expr*)xexpr), rlist }, 2));
        }
        lcase_add_assign(c, expr_copy((Expr*)var), xval);
    }
}

/* Branch B: the pivot `p` vanishes.  Solve p==0 for a parameter, substitute, and
 * recurse; OR over multiple roots.  Returns lsol_false() on an inconsistent
 * branch; sets *ok=false when p==0 cannot be solved rationally. */
static LSol* solve_pzero(Expr*** R, int nR, bool* remaining, int nv, Expr** vars,
                         const Expr* p, bool* ok) {
    Expr** params = NULL; int nparams = 0;
    collect_params(p, vars, nv, &params, &nparams);

    for (int pi = 0; pi < nparams; pi++) {
        Expr* call = expr_new_function(expr_new_symbol(SYM_Solve),
            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Equal),
                          (Expr*[]){ expr_copy((Expr*)p), expr_new_integer(0) }, 2),
                       expr_copy(params[pi]) }, 2);
        Expr* sols = eval_and_free(call);
        if (!is_head(sols, SYM_List) || sols->data.function.arg_count == 0) {
            expr_free(sols); continue;                 /* try another parameter */
        }
        bool clean = true;
        for (size_t t = 0; t < sols->data.function.arg_count && clean; t++) {
            Expr* row = sols->data.function.args[t];
            if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) clean = false;
            else {
                Expr* rule = row->data.function.args[0];
                if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) clean = false;
            }
        }
        if (!clean) { expr_free(sols); continue; }

        LSol* acc = lsol_false();
        for (size_t t = 0; t < sols->data.function.arg_count; t++) {
            Expr* rule = sols->data.function.args[t]->data.function.args[0];
            Expr* val = rule->data.function.args[1];
            Expr*** R2 = rows_subst(R, nR, nv, params[pi], val);
            bool ok2 = true;
            LSol* sub = lin_solve(R2, nR, remaining, nv, vars, &ok2);
            rows_free(R2, nR, nv);
            if (!ok2) {
                lsol_free(sub); lsol_free(acc);
                for (int i = 0; i < nparams; i++) expr_free(params[i]);
                free(params); expr_free(sols); *ok = false; return lsol_false();
            }
            Expr* cond = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy(params[pi]), expr_copy(val) }, 2);
            lsol_and_cond(sub, cond);
            expr_free(cond);
            acc = lsol_or(acc, sub);
        }
        expr_free(sols);
        for (int i = 0; i < nparams; i++) expr_free(params[i]);
        free(params);
        return acc;
    }

    for (int i = 0; i < nparams; i++) expr_free(params[i]);
    free(params);
    *ok = false; return lsol_false();
}

static LSol* lin_solve(Expr*** R, int nR, bool* remaining, int nv, Expr** vars, bool* ok) {
    if (!*ok) return lsol_false();

    /* Find a pivot (row r, remaining var j) with a nonzero coefficient, preferring
     * a nonzero constant (which needs no case split). */
    int pr = -1, pj = -1; bool p_const = false;
    for (int r = 0; r < nR && !p_const; r++) {
        for (int j = 0; j < nv; j++) {
            if (!remaining[j]) continue;
            if (rf_is_nonzero_const(R[r][j])) { pr = r; pj = j; p_const = true; break; }
        }
    }
    if (pr < 0) {
        for (int r = 0; r < nR && pr < 0; r++) {
            for (int j = 0; j < nv; j++) {
                if (!remaining[j]) continue;
                if (!rf_is_zero(R[r][j])) { pr = r; pj = j; break; }
            }
        }
    }

    if (pr < 0) {
        /* No remaining variable has a nonzero coefficient: they are all free.
         * The rows are pure parameter conditions (const == 0). */
        LSol* s = lsol_true();
        LCase* c = &s->c[0];
        for (int r = 0; r < nR; r++) {
            Expr* k = R[r][nv];
            if (rf_is_nonzero_const(k)) { lsol_free(s); return lsol_false(); }  /* inconsistent */
            if (rf_is_zero(k)) continue;
            lcase_add_cond(c, expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy(k), expr_new_integer(0) }, 2));
        }
        return s;
    }

    Expr* p = R[pr][pj];
    /* x_pj = -( sum_{i != pj} R[pr][i] * vars[i] + const ) / p */
    Expr* acc = expr_copy(R[pr][nv]);
    for (int i = 0; i < nv; i++) {
        if (i == pj || rf_is_zero(R[pr][i])) continue;
        Expr* term = rf_mul(R[pr][i], vars[i]);
        Expr* na = rf_add(acc, term);
        expr_free(acc); expr_free(term); acc = na;
    }
    Expr* negacc = rf_neg(acc); expr_free(acc);
    Expr* xexpr = rf_div(negacc, p); expr_free(negacc);

    /* Eliminate x_pj from the other rows using the pivot row. */
    int mR = nR - 1;
    Expr*** R2 = malloc((size_t)(mR ? mR : 1) * sizeof(Expr**));
    int w = 0;
    for (int r = 0; r < nR; r++) {
        if (r == pr) continue;
        Expr* fct = rf_div(R[r][pj], p);        /* factor = R[r][pj]/p */
        Expr** nr = malloc((size_t)(nv + 1) * sizeof(Expr*));
        for (int i = 0; i <= nv; i++) {
            Expr* prod = rf_mul(fct, R[pr][i]);
            /* Subtract takes ownership of prod; rf() then frees it. */
            nr[i] = rf(expr_new_function(expr_new_symbol(SYM_Subtract),
                (Expr*[]){ expr_copy(R[r][i]), prod }, 2));
        }
        expr_free(fct);
        R2[w++] = nr;
    }

    bool rem2[SYS_MAX_VARS];
    for (int i = 0; i < nv; i++) rem2[i] = remaining[i];
    rem2[pj] = false;

    LSol* result;
    if (p_const) {
        LSol* sub = lin_solve(R2, mR, rem2, nv, vars, ok);
        graft(sub, vars[pj], xexpr);
        result = sub;
    } else {
        /* Branch A: p != 0. */
        bool okA = *ok;
        LSol* subA = lin_solve(R2, mR, rem2, nv, vars, &okA);
        graft(subA, vars[pj], xexpr);
        Expr* ne = expr_new_function(expr_new_symbol(SYM_Unequal),
            (Expr*[]){ expr_copy(p), expr_new_integer(0) }, 2);
        lsol_and_cond(subA, ne); expr_free(ne);
        /* Branch B: p == 0 (substitute into the ORIGINAL rows and recurse). */
        bool okB = *ok;
        LSol* subB = solve_pzero(R, nR, remaining, nv, vars, p, &okB);
        if (!okA || !okB) { *ok = false; lsol_free(subA); lsol_free(subB); result = lsol_false(); }
        else result = lsol_or(subA, subB);
    }

    expr_free(xexpr);
    rows_free(R2, mR, nv);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Coefficient extraction / emission                                  *
 * ------------------------------------------------------------------ */

/* Extract the linear row of one EQ atom's poly; NULL (via *ok=false) if the poly
 * is not linear in the solve variables. */
static Expr** atom_row(const Expr* poly, Expr** vars, int nv, bool* ok) {
    Expr** row = malloc((size_t)(nv + 1) * sizeof(Expr*));
    /* constant term: poly with all vars -> 0 */
    Expr** rr = malloc((size_t)nv * sizeof(Expr*));
    for (int j = 0; j < nv; j++)
        rr[j] = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy(vars[j]), expr_new_integer(0) }, 2);
    Expr* rlist = expr_new_function(expr_new_symbol(SYM_List), rr, (size_t)nv);
    free(rr);
    row[nv] = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rlist }, 2));
    for (int j = 0; j < nv; j++)
        row[j] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
            (Expr*[]){ expr_copy((Expr*)poly), expr_copy(vars[j]) }, 2));
    /* linearity: poly - (const + sum coeff*var) == 0 */
    Expr* recon = expr_copy(row[nv]);
    for (int j = 0; j < nv; j++) {
        Expr* term = rf_mul(row[j], vars[j]);
        Expr* na = rf_add(recon, term);
        expr_free(recon); expr_free(term); recon = na;
    }
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)poly), recon }, 2));
    if (!rf_is_zero(diff)) *ok = false;
    expr_free(diff);
    if (!*ok) { row_free(row, nv); return NULL; }
    return row;
}

static Expr* lcase_to_expr(const LCase* c) {
    int total = c->ncond + c->nass;
    if (total == 0) return expr_new_symbol(SYM_True);
    Expr** parts = malloc((size_t)total * sizeof(Expr*));
    int p = 0;
    for (int i = 0; i < c->ncond; i++) parts[p++] = expr_copy(c->cond[i]);
    for (int i = 0; i < c->nass; i++)
        parts[p++] = expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(c->var[i]), expr_copy(c->val[i]) }, 2);
    Expr* out = (total == 1) ? parts[0]
        : expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)total);
    free(parts);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

static LSol* solve_conjunction(const RConj* conj, Expr** vars, int nv, bool* ok) {
    if (conj->n > SYS_MAX_ROWS || nv > SYS_MAX_VARS) { *ok = false; return lsol_false(); }
    Expr*** R = malloc((size_t)(conj->n ? conj->n : 1) * sizeof(Expr**));
    int nR = 0;
    for (int k = 0; k < conj->n; k++) {
        const RAtom* at = &conj->a[k];
        if (at->rel != R_EQ || at->nonconst_denom) { *ok = false; break; }
        Expr** row = atom_row(at->poly, vars, nv, ok);
        if (!*ok) break;
        R[nR++] = row;
    }
    if (!*ok) { rows_free(R, nR, nv); return lsol_false(); }

    bool remaining[SYS_MAX_VARS];
    for (int i = 0; i < nv; i++) remaining[i] = true;
    LSol* sol = lin_solve(R, nR, remaining, nv, vars, ok);
    rows_free(R, nR, nv);
    return sol;
}

Expr* reduce_eq_system(const RForm* F, Expr** vars, int nv) {
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    bool ok = true;
    LSol* all = lsol_false();
    for (int i = 0; i < F->n && ok; i++) {
        LSol* sub = solve_conjunction(F->c[i], vars, nv, &ok);
        if (!ok) { lsol_free(sub); break; }
        all = lsol_or(all, sub);
    }
    if (!ok) { lsol_free(all); return NULL; }

    Expr* out;
    if (all->n == 0) out = expr_new_symbol(SYM_False);
    else {
        Expr** disj = malloc((size_t)all->n * sizeof(Expr*));
        for (int i = 0; i < all->n; i++) disj[i] = lcase_to_expr(&all->c[i]);
        out = (all->n == 1) ? disj[0]
            : expr_new_function(expr_new_symbol(SYM_Or), disj, (size_t)all->n);
        free(disj);
    }
    lsol_free(all);
    return eval_and_free(out);
}
